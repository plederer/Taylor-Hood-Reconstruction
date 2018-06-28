#include "reconstruction.hpp"

void ReconstructionVertex::CalcReconstructionVertex(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfsigma)
{
  gfsigma->GetVector()=0.0;

  shared_ptr<MeshAccess> ma  = gfu->GetMeshAccess();
  shared_ptr<FESpace> fespace_s = gfsigma->GetFESpace();
  shared_ptr<FESpace> fespace_u = gfu->GetFESpace();

  static mutex add_indirect;
  
  ParallelForRange(nv, [&](IntRange r)
		   {
		     LocalHeap lh(100000, "local_reconstruction");
		     for (auto i : r) 
		       {
			 HeapReset hr(lh);
      
			 int ldofs = globaldofs[i].Size();			 
			 FlatVector<> localrhs(ldofs+meanvaldofs, lh);
			 localrhs=0.0;
      
			 for (auto j : IntRange(patch_elements[i].Size()))
			   {
			     HeapReset hr(lh);
			     int el = patch_elements[i][j];
			     ElementId ei(VOL,el);

			     const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
	  
			     const FiniteElement & felu = fespace_u -> GetFE(ei,lh);
			     const CompoundFiniteElement & cfelu = dynamic_cast<const CompoundFiniteElement&>(felu);
			     const ScalarFiniteElement<2> & h1fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelu[0]);

			     const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
			     const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
			     const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelur[1]);	

			     int nd_vr = l2fel.GetNDof();

			     Array<int> dnumsu(felu.GetNDof(), lh), dnumss(felur.GetNDof(), lh);

			     fespace_u->GetDofNrs(ei, dnumsu);
			     fespace_s->GetDofNrs(ei, dnumss);

			     IntRange udofs = cfelu.GetRange(0);
			     IntRange vdofs = cfelu.GetRange(1);

			     FlatVector<> eluv(felu.GetNDof(), lh), shape(nd_vr, lh), elf(nd_vr,lh);

			     gfu->GetElementVector(dnumsu, eluv);	  
			     /*
			     IntegrationRule ir(eltrans.GetElementType(), felu.Order()*2);
			     
			     elf = 0.0;
			     for(auto k : IntRange(ir.GetNIP()))
			       {
				 MappedIntegrationPoint<2,2> mip(ir[k], eltrans);

				 Vec<2> gradu;
				 Vec<2> gradv;
				 DiffOpGradient<2>::Apply(h1fel, mip, eluv.Range(udofs), gradu, lh);
				 DiffOpGradient<2>::Apply(h1fel, mip, eluv.Range(vdofs), gradv, lh);
				 Vec<1> divu = gradu[0]+gradv[1];			     

				 double fac = ir[k].Weight() * mip.GetMeasure();
				 l2fel.CalcShape(ir[k], shape);
				 elf += (fac * divu(0))*shape;
			       }
			     */

			     SIMD_IntegrationRule ir(eltrans.GetElementType(), felu.Order()*2);

			     elf = 0.0;
			     SIMD_MappedIntegrationRule<2,2> mir(ir, eltrans,lh);

			     FlatMatrix<SIMD<double>> gradu(2,mir.Size(),lh);
			     FlatMatrix<SIMD<double>> gradv(2,mir.Size(),lh);
			     
			     h1fel.EvaluateGrad(mir, eluv.Range(udofs), gradu);
			     h1fel.EvaluateGrad(mir, eluv.Range(vdofs), gradv);

			     FlatVector<SIMD<double>> divu(mir.Size(),lh);
			     divu = gradu.Row(0) +  gradv.Row(1);

			     for (size_t i =0; i< mir.Size(); i++)
			       divu(i)*=mir[i].GetWeight();
			     
			     l2fel.AddTrans(ir, divu, elf);
			     
			     
			     FlatArray<int> dnumss_l2part = dnumss.Range(cfelur.GetRange(1));

			     for(auto k : dnumss_l2part.Range())
			       localrhs(globaldofs[i].Pos(dnumss_l2part[k])) = elf(k);

			     /*
			     for(auto k : IntRange(dnumss_l2part.Size()))
			       for (auto l : IntRange(globaldofs[i].Size()))
				 if(dnumss_l2part[k] == globaldofs[i][l])
				    localrhs[l] = elf[k];
			     */
			   }

			 Vector<> vertexpatchrhs(ldofs+meanvaldofs);
			 vertexpatchrhs = 0.0;
			 vertexpatchrhs.Range(0,ldofs) = extrosc[i] * localrhs.Range(0,ldofs);

			 localrhs = patchmatinv[i] * vertexpatchrhs;
			 {
			   lock_guard<mutex> guard(add_indirect);
			   gfsigma->GetVector().AddIndirect(globaldofs[i], localrhs.Range(0,ldofs));
			 }
		       }
		   });
}		   

void ReconstructionVertex::CalcReconstructionVertexTrans(shared_ptr<LinearForm> fsigma, shared_ptr<LinearForm> fu)
{
  fu->GetVector()=0.0;
  shared_ptr<MeshAccess> ma  = fu->GetMeshAccess();
  shared_ptr<FESpace> fespace_s = fsigma->GetFESpace();
  shared_ptr<FESpace> fespace_u = fu->GetFESpace();

  static mutex add_indirect;
  
  ParallelForRange(nv, [&](IntRange r)
		   {
		     LocalHeap lh(100000, "local_reconstruction");
		     for (auto i : r) 
		       {
			 HeapReset hr(lh);
			 int ldofs = globaldofs[i].Size();      
			 Vector<> localrhs(ldofs);
			 Vector<> vertexpatchrhs(ldofs+meanvaldofs);
			 vertexpatchrhs = 0.0;
      
			 Vector<> mean_sol(ldofs+meanvaldofs);
			 Vector<> sol(ldofs);

			 fsigma->GetElementVector(globaldofs[i], vertexpatchrhs.Range(0,ldofs));

			 mean_sol = patchmatinv[i]* vertexpatchrhs;
			 sol = Trans(extrosc[i]) * mean_sol.Range(0,ldofs);   

			 for (auto j : IntRange(patch_elements[i].Size()))
			   {
			     HeapReset hr(lh);
			     int el = patch_elements[i][j];
			     ElementId ei(VOL,el);	
 
			     const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

			     const FiniteElement & felu = fespace_u -> GetFE(ei,lh);
			     const CompoundFiniteElement & cfelu = dynamic_cast<const CompoundFiniteElement&>(felu);
			     const ScalarFiniteElement<2> & h1fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelu[0]);

			     const FiniteElement & felur = fespace_s -> GetFE(ei,lh);	  
			     const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);	  	  
			     const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelur[1]);	

			     int nd_vr = l2fel.GetNDof();

			     Array<int> dnumsu(felu.GetNDof(), lh), dnumss(felur.GetNDof(), lh);

			     fespace_u->GetDofNrs(ei, dnumsu);
			     fespace_s->GetDofNrs(ei, dnumss);

			     IntRange udofs = cfelu.GetRange(0);
			     IntRange vdofs = cfelu.GetRange(1);

			     FlatVector<>  shape(nd_vr, lh), elu(felu.GetNDof(),lh), elfsig(nd_vr, lh);
	  
			     elu = 0.0;
			     elfsig = 0.0;
	  
			     FlatArray<int> dnumss_l2part = dnumss.Range(cfelur.GetRange(1));

			     for(auto k : IntRange(dnumss_l2part.Size()))
			       for (auto l : IntRange(globaldofs[i].Size()))
				 if((dnumss_l2part[k] == globaldofs[i][l]))
				   elfsig[k] = sol[l];			     	    

			     SIMD_IntegrationRule ir(eltrans.GetElementType(), felu.Order()*2);
			     SIMD_MappedIntegrationRule<2,2> mir(ir, eltrans,lh);

			     FlatVector<SIMD<double>> fac(mir.Size(), lh);
			     l2fel.Evaluate(ir,elfsig,fac);

			     for (size_t i =0; i< mir.Size(); i++)
			       fac(i)*=mir[i].GetWeight();

			     FlatMatrix<SIMD<double>> u1(2,mir.Size(),lh);
			     FlatMatrix<SIMD<double>> u2(2,mir.Size(),lh);

			     u1 = 0.0;
			     u2 = 0.0;
			     u1.Row(0) = fac;
			     u2.Row(1) = fac;
			     			     
			     h1fel.AddGradTrans(mir, u1, elu.Range(udofs));
			     h1fel.AddGradTrans(mir, u2, elu.Range(vdofs));			     
			     			     
			     {
			       lock_guard<mutex> guard(add_indirect);
			       fu->GetVector().AddIndirect(dnumsu, elu);
			     }
			   }
		       }
		   });
}


void ReconstructionVertex::CalcVertexPatches(shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini)
{  
  shared_ptr<MeshAccess> ma  = fespace_u->GetMeshAccess();
  
  BitArray fineedges(ma->GetNEdges());
  fineedges.Clear();

  IntRange numberdofs = dynamic_pointer_cast<CompoundFESpace>(fespace_s)->GetRange(2); 

  int l2order = (*(dynamic_pointer_cast<CompoundFESpace>(fespace_s)).get())[1]->GetOrder();

  if ((l2order>1) && (mini == false))
    {
      meanvalorder =l2order-1;
      L2HighOrderFE<ET_TRIG> meanvalfe(meanvalorder-1);
      meanvaldofs = meanvalfe.GetNDof();      
    }
  else
    {
      meanvalorder = 0;
      meanvaldofs = 0;
    }
  
  nv = ma->GetNV();
  ne = ma->GetNE();
  nedges = ma->GetNEdges();

  //Check for edges on fine mesh
  for (auto i : IntRange(ma->GetNE()))    
    for(int e : ma->GetElEdges(ElementId(VOL,i)))
      fineedges.Set(e);      
    

  TableCreator<int> creator(nv); 
  TableCreator<int> creatoredges(nv); 
  
  for( ; !creator.Done(); creator++)
    for (auto i : Range(ne))	  
      for (int j: ma->GetElVertices(ElementId(VOL,i)))
	creator.Add(j,i);
  
  BitArray diredge(ma->GetNEdges());
  diredge.Clear();

  auto fespace_u0 = (*dynamic_pointer_cast<CompoundFESpace>(fespace_u))[0];
  for (auto i : IntRange(ma->GetNSE()))
    if (fespace_u0->IsDirichletBoundary(ma->GetElIndex(ElementId(BND,i))))
      diredge.Set( ma->GetElEdges(ElementId(BND,i))[0]);      

  for (auto i : IntRange(fineedges.Size()))
    {
      if(!fineedges[i])
	diredge.Set(i);
    }
  
  for( ; !creatoredges.Done(); creatoredges++)
    for (auto i : IntRange(nedges))
      if(!diredge[i])
	for (auto v : ma->GetEdgePNums(i))
	  creatoredges.Add(v, i);	
  
  patch_elements = creator.MoveTable();
  patch_edges = creatoredges.MoveTable();

  const SparseMatrix<double> & mata = dynamic_cast<SparseMatrix<double> &>(bfa->GetMatrix());

  TableCreator<int> h1tol2creator(h1fespace->GetNDof());
  LocalHeap lh(1000000, "local_reconstruction");

  for( ; !h1tol2creator.Done(); h1tol2creator++)
    {
      for (auto i : IntRange(ne))
	{
	  Array<int> elementdofsh1;
	  h1fespace->GetDofNrs(ElementId(VOL,i), elementdofsh1);      

	  Array<int> elementdofs;
	  HeapReset hr(lh);

	  ElementId ei(VOL,i);
	  const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
	  const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
	  const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelur[1]);
	  
	  IntRange l2dofs = cfelur.GetRange(1);
	  fespace_s->GetDofNrs(ElementId(VOL,i), elementdofs);

	  int locall2dofs = l2fel.GetNDof();
	 
	  if (mini == true)
	    if (l2order == 2)
	      {		
		locall2dofs -= 3;
	      }
	    else
	      throw Exception("High order mini elements not implemented yet");
	  
	  for (auto j : IntRange(locall2dofs))
	    {
	      
	      h1tol2creator.Add(elementdofsh1[j],elementdofs[l2dofs.First()+j]);	      
	    }	 
	}
    }
  Table<int> h1tol2table = h1tol2creator.MoveTable();
  
  globaldofs.SetSize(nv);
  patchmatinv.SetSize(nv);
  extrosc.SetSize(nv);
  
  ParallelForRange(nv, [&](IntRange r)		   
		   {
		     LocalHeap slh=lh.Split();
		     for (auto i : r)
		       {
			 HeapReset hr(slh);
      
			 Array<int> & globallocaldofs = globaldofs[i];
            
			 for (auto j : IntRange(patch_elements[i].Size()))
			   {
			     Array<int> dofs;
			     fespace_s->GetInnerDofNrs(patch_elements[i][j], dofs);
			     globallocaldofs+=dofs;
			   }
         
			 for (auto j : IntRange(patch_edges[i].Size()))
			   {
			     Array<int> dofs;
			     fespace_s->GetEdgeDofNrs(patch_edges[i][j], dofs);
			     globallocaldofs+=dofs;
			   }

			 globallocaldofs+=numberdofs;

			 int ldofs = globallocaldofs.Size();
			 QuickSort(globallocaldofs);
      
			 Matrix<> localmata(ldofs+meanvaldofs);
			 Matrix<> bubbleprojmat(ldofs);
			 bubbleprojmat = Identity(ldofs);
			 localmata=0.0;

			 for (int j : IntRange(ldofs))
			   {	 
			     auto d = mata.GetRowIndices(globallocaldofs[j]);
			     auto v = mata.GetRowValues(globallocaldofs[j]);
			     int k=0;
			     int l=0;
	  
			     while( k < globallocaldofs.Size())
			       {
				 if(d[l] == globallocaldofs[k])
				   {
				     localmata(j,k) = v[l];
				     k++;
				     l++;
				   }
				 else if(d[l]>globallocaldofs[k])
				   k++;
				 else
				   l++;	    
			       }
			   }    

			 //Calculate bubble Projection and meanvals
			 Vec<2> Vertex0(0,0);
			 ma->GetPoint(i,Vertex0);

			 IntegrationRule bubbleweightpoints;
			 for (auto l : IntRange(l2order+1))
			   {
			     for (auto k : IntRange(l2order+1-l))
			       {
				 IntegrationPoint ip(double(l)/(l2order), double(k)/(l2order));
				 bubbleweightpoints.AddIntegrationPoint(ip);
			       }	  
			   }
            
			 for (auto j : patch_elements[i].Range())
			   {	      
			     int id = patch_elements[i][j];
			     ElementId ei(VOL,id);
	      
			     const ElementTransformation & eltrans = ma->GetTrafo (ei, slh);
			     const FiniteElement & felur = fespace_s -> GetFE(ei,slh);
			     const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
			     const HDivFiniteElement<2> & hdivfel = dynamic_cast<const HDivFiniteElement<2> &>(cfelur[0]);
			     const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelur[1]);
	  
			     L2HighOrderFE<ET_TRIG> meanvalfe(meanvalorder-1);	     
			     int nd = hdivfel.GetNDof();
			     int nd_vr = l2fel.GetNDof();
	  
			     FlatMatrix<> localmeanval(meanvaldofs, nd, slh);
			     FlatMatrixFixWidth<2> shape (hdivfel.GetNDof(), slh);
	      
			     IntegrationRule ir(eltrans.GetElementType(), hdivfel.Order()*2);

			     localmeanval = 0.0;

			     FlatVector<> shapemeanval(meanvaldofs, slh);
			     shapemeanval = 0.0;

			     if(meanvalorder !=0)
			       {
				 for (auto l : IntRange(ir.GetNIP()))
				   {
				     const MappedIntegrationPoint<2,2> mip(ir[l],eltrans);
				     double det = mip.GetJacobiDet();
				     IntegrationPoint physip(mip(0),mip(1));
		  
				     hdivfel.CalcMappedShape(mip, shape);
				     meanvalfe.CalcShape(physip, shapemeanval);
		  
				     Vec<2> xycoord(mip(1)-Vertex0(1),-1*mip(0)-Vertex0(0));
		  
				     localmeanval += fabs(det) *ir[l].Weight()*shapemeanval*Trans(shape* xycoord);
				   }
				 Array<int>  dnumss(felur.GetNDof(), slh);
				 fespace_s->GetDofNrs(ei, dnumss);	  	      
				 FlatArray<int> dnums_sigma = dnumss.Range(cfelur.GetRange(0));
	  
				 for(auto k : dnums_sigma.Range())
				   for(auto l : globaldofs[i].Range())		
				     if(dnums_sigma[k] == globaldofs[i][l])
				       {
					 if(meanvalorder !=0)
					   {
					     localmata.Rows(ldofs,ldofs+meanvaldofs).Col(l) += localmeanval.Col(k);
					     localmata.Cols(ldofs,ldofs+meanvaldofs).Row(l) += localmeanval.Col(k);
					   }
				       }
			       }

			     //Bubble Projection
			     Array<int> vnums;
			     vnums = ma->GetElVertices(ei);
			     int vertexdof =0;
			     for (int k = 0; k < 3; k++)
			       {
				 if(vnums[k] == i)
				   {
				     vertexdof = k;
				   }
			       }
	  
			     FlatMatrix<> basistransformation(nd_vr,nd_vr, slh);
			     FlatMatrix<> invbasistransformation(nd_vr,nd_vr, slh);
			     FlatMatrix<> shape_hats(3, nd_vr, slh);
	  
			     l2fel.CalcShape(bubbleweightpoints, basistransformation);
			     ScalarFE<ET_TRIG,1> p1fe;

			     p1fe.CalcShape(bubbleweightpoints, shape_hats);
	  
			     auto help = Trans(basistransformation);

			     CalcInverse(basistransformation, invbasistransformation);
	  
			     FlatMatrix<> localbubbleproj(nd_vr,nd_vr,slh);
			     FlatMatrix<> diag_shape_hat(nd_vr,nd_vr,slh);
			     diag_shape_hat = 0.0;
			     diag_shape_hat.Diag() = shape_hats.Row(vertexdof);

			     localbubbleproj = 0.0;
			     localbubbleproj = Trans(invbasistransformation) * diag_shape_hat * help;
	  
			     Array<int>  dnumss(felur.GetNDof(), slh);
			     fespace_s->GetDofNrs(ei, dnumss);	  	      
			     FlatArray<int> dnums_l2part = dnumss.Range(cfelur.GetRange(1));
	  
			     for(auto k : dnums_l2part.Range())
			       for(auto l : globaldofs[i].Range())		
				 if(dnums_l2part[k] == globaldofs[i][l])
				   {
				     for(auto  z : dnums_l2part.Range())
				       for(auto r : globaldofs[i].Range())
					 if(dnums_l2part[z] == globaldofs[i][r])
					   {
					     bubbleprojmat(r,l) = localbubbleproj(k,z);
					     bubbleprojmat(l,r) = localbubbleproj(z,k);
					   }
				   }
			   }
	

			 Matrix<> extrosclocal(ldofs);
			 Vector<> smooth(ldofs);
      
			 extrosclocal = Identity(ldofs);
     
			 //add all H1-dofs to elementdofsh1 that they do not appear twice 
			 Array<int> elementdofsh1;      
			 for (auto l : patch_elements[i].Range())
			   {	      
			     int id = patch_elements[i][l];
			     Array<int> elementdofselement;
			     h1fespace->GetDofNrs(ElementId(VOL,id), elementdofselement);
	  
			     for(auto j : elementdofselement.Range())
			       if(!(elementdofsh1.Contains(elementdofselement[j])))
				 elementdofsh1.Append(elementdofselement[j]);
			   }

			 //loop over all coupling dofs
			 //note that also the dofs on the boundary of the vertex are included
			 //but those dofs will be set to zero in the evaluation of the operator later
			 for(auto j : elementdofsh1.Range())
			   {
			     //if dofnr is no inner dofs
			     Array<int> smoothdofs;
			     smooth=0.0;
			     //dofs that couple
			     for (auto k : h1tol2table[elementdofsh1[j]].Range())
			       {
				 //check in globallocaldofs
				 for(auto s : globallocaldofs.Range())
				   if(globallocaldofs[s] == h1tol2table[elementdofsh1[j]][k])
				     {
				       smooth(s)=1;
				       smoothdofs.Append(s);
				     }
			       }
			     for (auto k : smoothdofs.Range())
				 for (auto s : smoothdofs.Range())
				     extrosclocal(smoothdofs[k],smoothdofs[s]) -= 1.0/h1tol2table[elementdofsh1[j]].Size();
			   }
      
			 CalcInverse(localmata);
			 
			 patchmatinv[i] = std::move(localmata);

			 extrosc[i].SetSize(ldofs);
			 extrosc[i] = (bubbleprojmat*extrosclocal | Lapack) ;
		       }		   
		   });
}
