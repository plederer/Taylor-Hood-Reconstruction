#include "reconstruction.hpp"

void ReconstructionVertex::CalcReconstructionVertex(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfsigma)
{
  static Timer timer("Reconstruction::CalcReconstructionVertex");
  RegionTimer reg (timer);
  cout<<"CalcReconstructionVertex"<<endl;
  gfsigma->GetVector()=0.0;

  shared_ptr<MeshAccess> ma  = gfu->GetMeshAccess();
  shared_ptr<FESpace> fespace_s = gfsigma->GetFESpace();
  shared_ptr<FESpace> fespace_u = gfu->GetFESpace();
 #pragma omp parallel for
  for (int i = 0; i < nv; ++i)
    {
      int ldofs = globaldofs[i].Size();
      Vector<> localrhs(ldofs+meanvaldofs);
      localrhs=0.0;
      
      LocalHeap lh(100000, "local_reconstruction");
      for (int j = 0; j < patch_elements[i].Size() ; j++)
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

	  IntegrationRule ir(eltrans.GetElementType(), felu.Order()*2);
	  
	   elf = 0.0;
	   for(int k = 0; k< ir.GetNIP(); k++)
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
	   	  
	  FlatArray<int> dnumss_l2part = dnumss.Range(cfelur.GetRange(1));
	  
	  for(int k=0; k< dnumss_l2part.Size(); k++)
	    for (int l = 0; l<globaldofs[i].Size(); ++l)
	      if(dnumss_l2part[k] == globaldofs[i][l])
		  localrhs[l] = elf[k];
	}

      Vector<> vertexpatchrhs(ldofs+meanvaldofs);
      vertexpatchrhs = 0.0;
      vertexpatchrhs.Range(0,ldofs) = extrosc[i] * localrhs.Range(0,ldofs);
      
      localrhs = patchmatinv[i] * vertexpatchrhs;
 #pragma omp critical    
      gfsigma->GetVector().AddIndirect(globaldofs[i], localrhs.Range(0,ldofs));
    }
}

void ReconstructionVertex::CalcReconstructionVertexTrans(shared_ptr<LinearForm> fsigma, shared_ptr<LinearForm> fu)
{
  static Timer timer("Reconstruction::CalcReconstructionVertexTrans");
  RegionTimer reg (timer);
  fu->GetVector()=0.0;
  shared_ptr<MeshAccess> ma  = fu->GetMeshAccess();
  shared_ptr<FESpace> fespace_s = fsigma->GetFESpace();
  shared_ptr<FESpace> fespace_u = fu->GetFESpace();
  
#pragma omp parallel for  
  for(int i = 0; i < nv; i++)
    {
      int ldofs = globaldofs[i].Size();      
      Vector<> localrhs(ldofs);
      Vector<> vertexpatchrhs(ldofs+meanvaldofs);
      vertexpatchrhs = 0.0;
      
      Vector<> mean_sol(ldofs+meanvaldofs);
      Vector<> sol(ldofs);

      fsigma->GetElementVector(globaldofs[i], vertexpatchrhs.Range(0,ldofs));
      
      //sol = 0.0;

      mean_sol = patchmatinv[i]* vertexpatchrhs;
      sol = Trans(extrosc[i]) * mean_sol.Range(0,ldofs);
     
      LocalHeap lh(100000, "local_reconstruction");

      for (int j = 0; j < patch_elements[i].Size() ; j++)
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

	  FlatVector<>  shape(nd_vr, lh), elu(felu.GetNDof(),lh), elui(felu.GetNDof(),lh), elfsig(nd_vr, lh);
	  
	  elu = 0.0;
	  elfsig = 0.0;
	  
	  IntegrationRule ir(eltrans.GetElementType(), felu.Order()*2);
	  
	  FlatArray<int> dnumss_l2part = dnumss.Range(cfelur.GetRange(1));
	  
	  for(int k=0; k<dnumss_l2part.Size(); k++)
	     for (int l = 0; l<globaldofs[i].Size(); ++l)
	       if((dnumss_l2part[k] == globaldofs[i][l]))
		 elfsig[k] = sol[l];
	  	  
	  elui = 0; 
	  for(int k = 0; k< ir.GetNIP(); k++)
	    {
	      MappedIntegrationPoint<2,2> mip(ir[k], eltrans);
	      
	      l2fel.CalcShape(ir[k], shape);

	      double fac = InnerProduct(elfsig, shape);
	      fac *= ir[k].Weight() * mip.GetMeasure();

	      Vec<2> u1 = { fac, 0 }, u2={0,fac};
	      
	      FlatVector<double> eluu=  elui.Range(udofs);
	      FlatVector<double> eluv=  elui.Range(vdofs);

	      DiffOpGradient<2>::ApplyTrans(h1fel, mip, u1, eluu, lh);
	      DiffOpGradient<2>::ApplyTrans(h1fel, mip, u2, eluv, lh);
	      elu += elui;
	    }
#pragma omp critical
	  fu->GetVector().AddIndirect(dnumsu, elu);
	}      
   }
}


void ReconstructionVertex::CalcVertexPatches(shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini)
{
  static Timer timer("Reconstruction::CalcVertexPatches");
  static Timer timertable("Reconstruction::CreateTables");
  static Timer timerh1("Reconstruction::H1toL2");
  static Timer timerlocala("Reconstruction::locala");
  static Timer timerbubble("Reconstruction::Bubbleproj");
  static Timer timermeanvals("Reconstruction::meanvals");
  static Timer timersmooth("Reconstruction::smooth");
  static Timer timersmoothtest("Reconstruction::smoothtest");
  static Timer timerinv("Reconstruction::inv");
  static Timer timermult("Reconstruction::mult");
  RegionTimer reg (timer);
  cout<<"CalcVertexPatches"<<endl;
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
  for (int i = 0; i < ma->GetNE(); ++i)
    {
      Array<int> enums;
      enums = ma->GetElEdges(ElementId(VOL,i));
      for(int e : enums)
	fineedges.Set(e);      
    }

  timertable.Start();
  TableCreator<int> creator(nv); 
  TableCreator<int> creatoredges(nv); 

  
  for( ; !creator.Done(); creator++)
    {
      for (int i = 0; i < ne; ++i)
	{
	  Array<int> vnums;
	  vnums = ma->GetElVertices(i);      
	  for (int j: vnums)
	    creator.Add(j,i);

	}      
    }
  
  BitArray diredge(ma->GetNEdges());
  diredge.Clear();
  
  for (int i = 0; i < ma->GetNSE(); ++i)
    {
      Array<int> enums;
      
      enums = ma->GetElEdges(ElementId(BND,i));
      if ((*dynamic_pointer_cast<CompoundFESpace>(fespace_u))[0]->IsDirichletBoundary(ma->GetElIndex(ElementId(BND,i))))
	diredge.Set(enums[0]);      
    }

  for (int i = 0; i < fineedges.Size(); ++i)
    {
      if(!fineedges[i])
	diredge.Set(i);
    }
  
  for( ; !creatoredges.Done(); creatoredges++)
    {
      for (int i = 0; i < nedges; ++i)
	{
	  if(!diredge[i])
	    {
	      int v1,v2;
	      Array<int> vnums;
	      vnums = ma->GetEdgePNums(i);
	      
	      //ma->GetEdgePNums(i,v1,v2);
	      
	      //creatoredges.Add(v1,i);
	      //creatoredges.Add(v2,i);
	      creatoredges.Add(vnums[0],i);
	      creatoredges.Add(vnums[1],i);

	    }
	}      
    }
  
  patch_elements = creator.MoveTable();
  patch_edges = creatoredges.MoveTable();

  timertable.Stop();
  const SparseMatrix<double> & mata = dynamic_cast<SparseMatrix<double> &>(bfa->GetMatrix());

  timerh1.Start();
  TableCreator<int> h1tol2creator(h1fespace->GetNDof());
  LocalHeap lh(1000000, "local_reconstruction");

  for( ; !h1tol2creator.Done(); h1tol2creator++)
    {
      for (int i = 0; i < ne; ++i)
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
	  
	  for (int j = 0; j < locall2dofs ; j++ )
	    {
	      
	      h1tol2creator.Add(elementdofsh1[j],elementdofs[l2dofs.First()+j]);	      
	    }	 
	}
    }
  Table<int> h1tol2table = h1tol2creator.MoveTable();
  timerh1.Stop();
  
  globaldofs.SetSize(nv);
  patchmatinv.SetSize(nv);
  extrosc.SetSize(nv);
  
#pragma omp parallel for
    //ParallelFor(Range(nv),[&](int i)
  for (int i = 0; i < nv; ++i)
    {
      //LocalHeap lh=lhglobal.Split();
      HeapReset hr(lh);
      
      Array<int> & globallocaldofs = globaldofs[i];
            
      for (int j = 0; j < patch_elements[i].Size() ; j++)
	{
	  Array<int> dofs;
	  fespace_s->GetInnerDofNrs(patch_elements[i][j], dofs);
	  globallocaldofs+=dofs;
	}
         
      for (int j = 0; j < patch_edges[i].Size() ; j++)
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

      timerlocala.Start();
      for (int j = 0; j < ldofs; ++j)
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
      
      timerlocala.Stop();

      //Calculate bubble Projection and meanvals
      Vec<2> Vertex0(0,0);
      ma->GetPoint(i,Vertex0);

      IntegrationRule bubbleweightpoints;
      for (int l = 0; l < l2order+1; ++l)
	{
	  for (int k = 0; k < l2order+1-l; ++k)
	    {
	      IntegrationPoint ip(double(l)/(l2order), double(k)/(l2order));
	      bubbleweightpoints.AddIntegrationPoint(ip);
	    }	  
	}
            
      for (int j = 0; j < patch_elements[i].Size(); ++j)
	{
	      
	  int id = patch_elements[i][j];
	  ElementId ei(VOL,id);
	      
	  const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
	  const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
	  const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
	  const HDivFiniteElement<2> & hdivfel = dynamic_cast<const HDivFiniteElement<2> &>(cfelur[0]);
	  const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelur[1]);
	  
	  L2HighOrderFE<ET_TRIG> meanvalfe(meanvalorder-1);	     
	  int nd = hdivfel.GetNDof();
	  int nd_vr = l2fel.GetNDof();
	  
	  FlatMatrix<> localmeanval(meanvaldofs, nd, lh);
	  FlatMatrixFixWidth<2> shape (hdivfel.GetNDof(), lh);
	      
	  IntegrationRule ir(eltrans.GetElementType(), hdivfel.Order()*2);

	  localmeanval = 0.0;

	  FlatVector<> shapemeanval(meanvaldofs, lh);
	  shapemeanval = 0.0;
	  timermeanvals.Start();
	  if(meanvalorder !=0)
	    {
	      for (int l = 0; l < ir.GetNIP(); l++)
		{
		  const MappedIntegrationPoint<2,2> mip(ir[l],eltrans);
		  double det = mip.GetJacobiDet();
		  IntegrationPoint physip(mip(0),mip(1));
		  
		  hdivfel.CalcMappedShape(mip, shape);
		  meanvalfe.CalcShape(physip, shapemeanval);
		  
		  Vec<2> xycoord(mip(1)-Vertex0(1),-1*mip(0)-Vertex0(0));
		  
		  localmeanval += fabs(det) *ir[l].Weight()*shapemeanval*Trans(shape* xycoord);
		}
	      Array<int>  dnumss(felur.GetNDof(), lh);
	      fespace_s->GetDofNrs(ei, dnumss);	  	      
	      FlatArray<int> dnums_sigma = dnumss.Range(cfelur.GetRange(0));
	  
	      for(int k = 0; k < dnums_sigma.Size(); k++)
		for(int l = 0; l<globaldofs[i].Size(); l++)		
		  if(dnums_sigma[k] == globaldofs[i][l])
		    {
		      if(meanvalorder !=0)
			{
			  localmata.Rows(ldofs,ldofs+meanvaldofs).Col(l) += localmeanval.Col(k);
			  localmata.Cols(ldofs,ldofs+meanvaldofs).Row(l) += localmeanval.Col(k);
			}
		    }
	    }
	  timermeanvals.Stop();

	  timerbubble.Start();
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
	  
	  FlatMatrix<> basistransformation(nd_vr,nd_vr, lh);
	  FlatMatrix<> invbasistransformation(nd_vr,nd_vr, lh);
	  FlatMatrix<> shape_hats(3, nd_vr, lh);
	  
	  l2fel.CalcShape(bubbleweightpoints, basistransformation);
	  ScalarFE<ET_TRIG,1> p1fe;

	  p1fe.CalcShape(bubbleweightpoints, shape_hats);
	  
	  auto help = Trans(basistransformation);

	  CalcInverse(basistransformation, invbasistransformation);
	  
	  FlatMatrix<> localbubbleproj(nd_vr,nd_vr,lh);
	  FlatMatrix<> diag_shape_hat(nd_vr,nd_vr,lh);
	  diag_shape_hat = 0.0;
	  diag_shape_hat.Diag() = shape_hats.Row(vertexdof);

	  localbubbleproj = 0.0;
	  localbubbleproj = Trans(invbasistransformation) * diag_shape_hat * help;
	  
	  Array<int>  dnumss(felur.GetNDof(), lh);
	  fespace_s->GetDofNrs(ei, dnumss);	  	      
	  FlatArray<int> dnums_l2part = dnumss.Range(cfelur.GetRange(1));
	  
	  for(int k = 0; k < dnums_l2part.Size(); k++)
	      for(int l = 0; l<globaldofs[i].Size(); l++)		
		if(dnums_l2part[k] == globaldofs[i][l])
		{
		  for(int z = 0; z < dnums_l2part.Size(); z++)
		    for(int r = 0; r<globaldofs[i].Size(); r++)
		      if(dnums_l2part[z] == globaldofs[i][r])
			{
			  bubbleprojmat(r,l) = localbubbleproj(k,z);
			  bubbleprojmat(l,r) = localbubbleproj(z,k);
			}
		}
	  timerbubble.Stop();
	}
	

      Matrix<> extrosclocal(ldofs);
      Vector<> smooth(ldofs);
      
      extrosclocal = Identity(ldofs);
     
      //add all H1-dofs to elementdofsh1 that they do not appear twice 
      Array<int> elementdofsh1;      
      for (int l = 0; l < patch_elements[i].Size(); ++l)
	{	      
	  int id = patch_elements[i][l];
	  Array<int> elementdofselement;
	  h1fespace->GetDofNrs(ElementId(VOL,id), elementdofselement);
	  
	  for(int j = 0; j <elementdofselement.Size(); j++)
	    if(!(elementdofsh1.Contains(elementdofselement[j])))
	      elementdofsh1.Append(elementdofselement[j]);
	}

      //loop over all coupling dofs
      //note that also the dofs on the boundary of the vertex are included
      //but those dofs will be set to zero in the evaluation of the operator later


 
      for(int j = 0; j < elementdofsh1.Size(); j++)
	{
	  //if dofnr is no inner dofs
	  Array<int> smoothdofs;
	  smooth=0.0;
	  //smoothmat=0.0;
	  // if(h1tol2table[elementdofsh1[j]].Size()>0)
	  //{
	      //dofs that couple
	      for (int k = 0; k < h1tol2table[elementdofsh1[j]].Size(); ++k)
		{
		  //check in globallocaldofs
		  for(int s = 0; s<globallocaldofs.Size(); s++)
		    if(globallocaldofs[s] == h1tol2table[elementdofsh1[j]][k])
		      {
			smooth(s)=1;
			smoothdofs.Append(s);
		      }
		}
	      timersmooth.Start();
	      //extrosclocal-=1.0/h1tol2table[elementdofsh1[j]].Size()*(smooth*Trans(smooth));
	      for (int k =0; k < smoothdofs.Size(); ++k)
		{
		  for (int s =0; s < smoothdofs.Size(); ++s)
		    {
		      extrosclocal(smoothdofs[k],smoothdofs[s]) -= 1.0/h1tol2table[elementdofsh1[j]].Size();
		    }
		  
		}
	      timersmooth.Stop();
	}

     
      
      
      timerinv.Start();
      CalcInverse(localmata);
      timerinv.Stop();

      //patchmatinv[i].SetSize(ldofs+meanvaldofs);
      patchmatinv[i] = std::move(localmata);

      extrosc[i].SetSize(ldofs);
      timermult.Start();
      extrosc[i] = (bubbleprojmat*extrosclocal | Lapack) ;
      timermult.Stop();
      //timerinv.AddFlops(bubbleprojmat.Width()*bubbleprojmat.Width()*bubbleprojmat.Width());

    }
}
