#include "reconstruction.hpp"
#include <python_ngstd.hpp>

void ReconstructionElement::CalcReconstructionElement(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfsigma)
//void ReconstructionElement::CalcReconstructionElement(PyGF gfu,PyGF gfsigma)
{
  gfsigma->GetVector()=0.0;
  shared_ptr<MeshAccess> ma  = gfu->GetMeshAccess();
  shared_ptr<FESpace> fespace_s = gfsigma->GetFESpace();
  shared_ptr<FESpace> fespace_u = gfu->GetFESpace();
 
  IntRange numberdofs = dynamic_pointer_cast<CompoundFESpace>(fespace_s)->GetRange(2); 

  LocalHeap clh(100000, "local_reconstruction");

#pragma omp parallel for
  for (int i = 0; i < ne; ++i)
    {
      LocalHeap lh = clh.Split();
      int ldofs = globaldofs[i].Size();
	
      //Calc local rhs:
      Vector<> localrhs(ldofs);
      Vector<> elementpatchrhs(ldofs+meanvaldofs);

      elementpatchrhs = 0.0;
      localrhs=0.0;

      ElementId ei(VOL,i);
      const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

      	  
      const FiniteElement & felu = fespace_u -> GetFE(ei,lh);
      
      const CompoundFiniteElement & cfelu = dynamic_cast<const CompoundFiniteElement&>(felu);
      
      const ScalarFiniteElement<2> & h1fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelu[0]);
            
      const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
      const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
      const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelur[1]);	
      int nd_u = h1fel.GetNDof();
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
	  elf += (fac * divu(0)) * shape;
	}
	  
      FlatArray<int> dnumss_l2part = dnumss.Range(cfelur.GetRange(1));
      for(int k=0; k< dnumss_l2part.Size(); k++)
	for (int l = 0; l<globaldofs[i].Size(); ++l)
	  if(dnumss_l2part[k] == globaldofs[i][l])
	    localrhs[l] = elf[k];
      
      elementpatchrhs.Range(0,ldofs)= extrosc[i] * localrhs;
      //elementpatchrhs.Range(0,ldofs)=  localrhs;

      Vector<> sol(ldofs+meanvaldofs);
      
      sol = 0.0;
      sol = patchmatinv[i] * elementpatchrhs;
#pragma omp critical
      gfsigma->GetVector().AddIndirect(globaldofs[i], sol.Range(0,ldofs));
    }
}

void ReconstructionElement::CalcReconstructionElementTrans(shared_ptr<LinearForm> fsigma, shared_ptr<LinearForm> fu)
//void ReconstructionElement::CalcReconstructionElementTrans(PyLF fsigma, PyLF fu)
{
  fu->GetVector()=0.0;
  shared_ptr<MeshAccess> ma  = fu->GetMeshAccess();
  shared_ptr<FESpace> fespace_s = fsigma->GetFESpace();
  shared_ptr<FESpace> fespace_u = fu->GetFESpace();

  LocalHeap clh(1000000, "local_reconstruction");
  
#pragma omp parallel for
  for (int i = 0; i < ne; ++i)   
    {
      LocalHeap lh = clh.Split();
      int ldofs = globaldofs[i].Size();
      
      Vector<> localsol(ldofs);
      localsol=0.0;
      Vector<> elementpatchrhs(ldofs);
      elementpatchrhs = 0.0;
      
      fsigma->GetElementVector(globaldofs[i], elementpatchrhs);
      
      //Vector<> sol(ldofs);

      Vector<> mean_elementpatchrhs(ldofs+meanvaldofs);
      mean_elementpatchrhs = 0.0;
      //cout<<"mean_elementpatchrhs = "<<mean_elementpatchrhs<<endl;
      //cout<<"elementpatchrhs = "<<elementpatchrhs<<endl;

      mean_elementpatchrhs.Range(0,ldofs) = elementpatchrhs;
      //cout<<"mean_elementpatchrhs = "<<mean_elementpatchrhs<<endl;
      //int s;
      //cin>>s;
      Vector<> mean_sol(ldofs+meanvaldofs);
      mean_sol=0.0;
            
      mean_sol = patchmatinv[i] * mean_elementpatchrhs;
      
      
      localsol = extrosc[i] * mean_sol.Range(0,ldofs);
      
      ElementId ei(VOL,i);
      const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
      	  
      const FiniteElement & felu = fespace_u -> GetFE(ei,lh);
      const CompoundFiniteElement & cfelu = dynamic_cast<const CompoundFiniteElement&>(felu);
      const ScalarFiniteElement<2> & h1fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelu[0]);
            
      const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
      const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
      const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelur[1]);	
      int nd_u = h1fel.GetNDof();
      int nd_vr = l2fel.GetNDof();

      Array<int> dnumsu(felu.GetNDof(), lh), dnumss(felur.GetNDof(), lh);

      fespace_u->GetDofNrs(ei, dnumsu);
      fespace_s->GetDofNrs(ei, dnumss);


      IntRange udofs = cfelu.GetRange(0);
      IntRange vdofs = cfelu.GetRange(1);

      FlatVector<> elu(felu.GetNDof(), lh), shape(nd_vr, lh), elf(nd_vr,lh), elui(felu.GetNDof(), lh);
      
      elu = 0.0;      
      elf = 0.0;

      FlatArray<int> dnumss_l2part = dnumss.Range(cfelur.GetRange(1));
      for(int k=0; k< dnumss_l2part.Size(); k++)
	for (int l = 0; l<globaldofs[i].Size(); ++l)
	  if(dnumss_l2part[k] == globaldofs[i][l])
	    elf[k] = localsol[l];    
      
      IntegrationRule ir(eltrans.GetElementType(), felu.Order()*2);

      elui = 0.0; //for p component
      for(int k = 0; k< ir.GetNIP(); k++)
	{
	  MappedIntegrationPoint<2,2> mip(ir[k], eltrans);
	      
	  l2fel.CalcShape(ir[k], shape);

	  double fac = InnerProduct(shape, elf);
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


void ReconstructionElement::CalcElementPatches(shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini)
//void ReconstructionElement::CalcElementPatches(PyFES  fespace_u, PyFES fespace_s, PyBF bfa, PyFES h1fespace)
{
  shared_ptr<MeshAccess> ma  = fespace_u->GetMeshAccess();

  IntRange numberdofs = dynamic_pointer_cast<CompoundFESpace>(fespace_s)->GetRange(2); 
  int s;
  /*meanvalorder = (*(dynamic_pointer_cast<CompoundFESpace>(fespace_u)).get())[0]->GetOrder();

  if (mini)
    cout<<"Mini element not considered for element reconstruction"<<endl;
  
  if (meanvalorder>2)
    {
      meanvalorder -=2;
      L2HighOrderFE<ET_TRIG> meanvalfe(meanvalorder-1);
      meanvaldofs = meanvalfe.GetNDof();      
    }
  else
    {
      meanvalorder = 0;
      meanvaldofs = 0;
    }
  */
  
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

  BitArray fineedges(ma->GetNEdges());
  fineedges.Clear();
  //Check for edges on fine mesh
  for (int i = 0; i < ma->GetNE(); ++i)
    {
      Array<int> enums;
      enums = ma->GetElEdges(ElementId(VOL,i));
      for(int e : enums)
	fineedges.Set(e);      
    }


  TableCreator<int> creator(ne); 
  TableCreator<int> creatoredges(ne); 
  
  patchmatinv.SetSize(ne);
  extrosc.SetSize(ne);
  globaldofs.SetSize(ne);
  
  for( ; !creator.Done(); creator++)
    {
      for (int i = 0; i < ne; ++i)
	{
	  Array<int> elnums;	  
	  Array<int> vnums;
	  vnums = ma->GetElVertices(ElementId(VOL,i));
	  
	  for (int j=0;j<3;j++)
	    {
	      //add elements
	      Array<int> venums;
	      ma->GetVertexElements(vnums[j],venums);
	      for (int k : venums)
		if (!elnums.Contains(k))
		  elnums.Append(k);
	    }
	  for(int k=0;k<elnums.Size(); k++)
	    {
	      if((elnums[k] != i))
		creator.Add(i,elnums[k]);
	    }
	}
    }
  

  Table<int> patch_elements = creator.MoveTable();

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
      for (int i = 0; i < ne; ++i)
	{
	  Array<int> ednums;
	  Array<int> vnums;
	  vnums = ma->GetElVertices(ElementId(VOL,i));
	  
	  for (int k = 0; k < nedges; ++k)
	    {
	      if(!diredge[k])
		{
		  int v1,v2;
		  Array<int> pnums;
		  pnums = ma->GetEdgePNums(k);
		  v1 = pnums[0];
		  v2 = pnums[1];
		  if((vnums.Contains(v1)) || (vnums.Contains(v2)))
		    {
		      if (!ednums.Contains(k))
			ednums.Append(k);
		    }		      		  
		}
	    }
	  creatoredges.Add(i, ednums);
	}      
    }
  Table<int> patch_edges = creatoredges.MoveTable();

  const SparseMatrix<double> & mata = dynamic_cast<SparseMatrix<double> &>(bfa->GetMatrix());

  TableCreator<int> h1tol2creator(h1fespace->GetNDof());
  LocalHeap lh(1000000, "local_reconstruction");

  for( ; !h1tol2creator.Done(); h1tol2creator++)
    {
      for (int i = 0; i < ne; ++i)
	{
	  Array<int> elementdofsh1;
	  h1fespace->GetDofNrs(ElementId(VOL,i), elementdofsh1);      

	  //Array<int> innerdofs;
	  //h1fespace->GetInnerDofNrs(i, innerdofs);

	  Array<int> elementdofs;
	  HeapReset hr(lh);
	  //ElementId ei(VOL,patch_elements[i][k]);
	  ElementId ei(VOL,i);
	  const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
	  const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
	  const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelur[1]);
	  
	  IntRange l2dofs = cfelur.GetRange(1);
	  fespace_s->GetDofNrs(ElementId(VOL,i), elementdofs);
	  
	  for (int j = 0; j < l2fel.GetNDof(); j++ )
	    {
	      //just vertex and edge dofs couple!!!
	      //if(!innerdofs.Contains(elementdofsh1[j]))
	      h1tol2creator.Add(elementdofsh1[j],elementdofs[l2dofs.First()+j]);	      
	    }	 
	}
    }
  Table<int> h1tol2table = h1tol2creator.MoveTable();
  
  
  for (int i = 0; i < ne; ++i)
    {
      
      HeapReset hr(lh);
      //cout<<"patch: "<<i<<endl;
      Array<int> & globallocaldofs = globaldofs[i];
            
      for (int j = 0; j < patch_elements[i].Size() ; j++)
	{
	  Array<int> dofs;
	  fespace_s->GetInnerDofNrs(patch_elements[i][j], dofs);
	  globallocaldofs+=dofs;
	}
      
      //element i
      Array<int> dofs;
      fespace_s->GetInnerDofNrs(i, dofs);
      globallocaldofs+=dofs;
     
      for (int j = 0; j < patch_edges[i].Size() ; j++)
	{
	  Array<int> dofs;
	  fespace_s->GetEdgeDofNrs(patch_edges[i][j], dofs);
	  globallocaldofs+=dofs;
	}

      globallocaldofs+=numberdofs;

      int  ldofs = globallocaldofs.Size();
      QuickSort(globallocaldofs);
      
      Matrix<> localmata(ldofs+meanvaldofs);
      localmata=0.0;
      for (int j = 0; j < ldofs; ++j)
	{
	  for (int k = 0; k < ldofs; ++k)
	    {
	      localmata(j,k) = mata(globallocaldofs[j],globallocaldofs[k]);
	    }
	}
    
      if(meanvalorder!=0)
	{
	  //cout<<"i = "<<i<<endl;
	  //cout<<"patch_elements[i] = "<<patch_elements[i]<<endl;
	  //cin>>s;

	  Array<int> all_patch_elements(patch_elements[i].Size()+1);
	  all_patch_elements[0]=i;	  
	  for(int j = 1; j < patch_elements[i].Size()+1; ++j)
	    all_patch_elements[j]=patch_elements[i][j-1];
	  //cout<<"all_patch_elements = "<<all_patch_elements<<endl;
	  //cout<<"patch_elements[i] = "<<patch_elements[i]<<endl;
	  //cin>>s;

	  
	  Array<int> vnums;
	  vnums = ma->GetElVertices(ElementId(VOL,i));
	  
	  Vec<2> Vertex0(0,0);
	  ma->GetPoint(vnums[0],Vertex0);

	  
	  for (int j = 0; j < all_patch_elements.Size(); ++j)
	    {
	      int id = all_patch_elements[j];
	      ElementId ei(VOL,id);

	      const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
	            
	      const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
	      const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
	      const HDivFiniteElement<2> & hdivfel = dynamic_cast<const HDivFiniteElement<2> &>(cfelur[0]);

	      L2HighOrderFE<ET_TRIG> meanvalfe(meanvalorder-1);
	      
	      
	      int nd = hdivfel.GetNDof();
	      FlatMatrix<> localmeanval(meanvaldofs, nd, lh);
	      FlatMatrixFixWidth<2> shape (hdivfel.GetNDof(), lh);
	      
	      IntegrationRule ir(eltrans.GetElementType(), hdivfel.Order()*2);

	      localmeanval = 0.0;
	      //cout<<"localmeanval = "<<localmeanval<<endl;
	      //cin>>s;	     	      

	      FlatVector<> shapemeanval(meanvaldofs, lh);
	      shapemeanval = 0.0;
	      
	      for (int l = 0; l < ir.GetNIP(); l++)
		{
		  const MappedIntegrationPoint<2,2> mip(ir[l],eltrans);

		  double det = mip.GetJacobiDet();
		  //FlatMatrix<> jacobian(mip.GetJacobian());
		  hdivfel.CalcMappedShape(mip, shape);
		  //IntegrationPoint physip(mip(0)-cent(0),mip(1)-cent(1));
		  IntegrationPoint physip(mip(0),mip(1));
		  //cout<<"shapemeanval = "<<shapemeanval<<endl;

		  meanvalfe.CalcShape(physip, shapemeanval);
		  //cout<<"physip = "<<physip<<endl;
		  //cout<<"mip = "<<mip<<endl;
		  //cout<<"cent = "<<cent<<endl;


		  //cout<<"shapemeanval = "<<shapemeanval<<endl;
		  
		  //hdivfel.CalcShape(ir[l], shape);
		  //@Joachim: bei dem schon den schwerpunkt abziehen???
		  //Vec<2> xycoord(mip(1)-cent(1),-1*(mip(0)-cent(0)));
		  Vec<2> xycoord(mip(1)-Vertex0(1),-1*mip(0)-Vertex0(0));
		  //Vec<2> xycoord(mip(1),-1*mip(0));
		  //Vec<2> xycoord(mip(1),-1*mip(0));
		
		  localmeanval += fabs(det) *ir[l].Weight()*shapemeanval*Trans(shape* xycoord);
		  
		  //cout<<"localmeanval.Height() = "<<localmeanval.Height()<<endl;
		  //cout<<"localmeanval.Width() = "<<localmeanval.Width()<<endl;

		  //cin>>s;

		  //cout<<"localmeanval = "<<localmeanval<<endl;
		  
		  //cin>>s;
		}
	      //cout<<"BBB"<<endl; 
	      Array<int>  dnumss(felur.GetNDof(), lh);

	      fespace_s->GetDofNrs(ElementId(VOL,id), dnumss);
	      
	      FlatArray<int> dnums_sigma = dnumss.Range(cfelur.GetRange(0));
	      for(int k = 0; k < dnums_sigma.Size(); k++)
		for(int l = 0; l<globaldofs[i].Size(); l++)
		  if(dnums_sigma[k] == globaldofs[i][l])
		    {
		      localmata.Rows(ldofs,ldofs+meanvaldofs).Col(l) += localmeanval.Col(k);
		      localmata.Cols(ldofs,ldofs+meanvaldofs).Row(l) += localmeanval.Col(k);		      
		    }	     
	    }
	}

      //cin>>s;
      

      Matrix<> extrosclocal(ldofs);
      Vector<> smooth(ldofs);
      
      extrosclocal = Identity(ldofs);

      //Array<int> elementdofsl2;
      //fespace_s->GetDofNrs(i, elementdofsl2);

      Array<int> elementdofsh1;
      h1fespace->GetDofNrs(ElementId(VOL,i), elementdofsh1);

      for(int j = 0; j < elementdofsh1.Size(); j++)
	{
	  //if dofnr is no inner dofs	 
	  smooth=0.0;
	  //smoothmat=0.0;
	  //if(h1tol2table[elementdofsh1[j]].Size()>0)
	  //{
	      //dofs that couple
	      for (int k = 0; k < h1tol2table[elementdofsh1[j]].Size(); ++k)
		{
		  //check in globallocaldofs
		  for(int s = 0; s<globallocaldofs.Size(); s++)
		    if(globallocaldofs[s] == h1tol2table[elementdofsh1[j]][k])
		      smooth(s)=1;		    		    
		}	      
	      extrosclocal-=1.0/h1tol2table[elementdofsh1[j]].Size()*smooth*Trans(smooth);
	      //  }
	}

      CalcInverse(localmata);
      
      patchmatinv[i].SetSize(ldofs+meanvaldofs);
      patchmatinv[i] = localmata;

      extrosc[i].SetSize(ldofs);
      extrosc[i] = extrosclocal;                   
    }
}

