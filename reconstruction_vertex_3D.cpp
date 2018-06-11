//#define DEBUG
//#include <python_ngstd.hpp>
//#include <solve.hpp>
//#include <fem.hpp>  


//using namespace ngsolve;
//using namespace ngcomp;
//using namespace ngfem;
#include "reconstruction.hpp"

void ReconstructionVertex3D::CalcReconstructionVertex(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfsigma)
//void ReconstructionVertex3D::CalcReconstructionVertex(PyGF gfu, PyGF gfsigma)
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
      Vector<> BubbleProj(ldofs);
      BubbleProj = 0.0;
      localrhs=0.0;
      
      LocalHeap lh(1000000, "local_reconstruction");
      for (int j = 0; j < patch_elements[i].Size() ; j++)
	{	  
	  HeapReset hr(lh);
	  int el = patch_elements[i][j];	  

	  ElementId ei(VOL,el);

	  const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
	  
	  const FiniteElement & felu = fespace_u -> GetFE(ei,lh);
	  const CompoundFiniteElement & cfelu = dynamic_cast<const CompoundFiniteElement&>(felu);
	  const ScalarFiniteElement<3> & h1fel = dynamic_cast<const ScalarFiniteElement<3> &> (cfelu[0]);

	  const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
	  const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
	  const ScalarFiniteElement<3> & l2fel = dynamic_cast<const ScalarFiniteElement<3> &> (cfelur[1]);	

	  int nd_vr = l2fel.GetNDof();

	  Array<int> dnumsu(felu.GetNDof(), lh), dnumss(felur.GetNDof(), lh);

	  fespace_u->GetDofNrs(ei, dnumsu);
	  fespace_s->GetDofNrs(ei, dnumss);

	  IntRange udofs = cfelu.GetRange(0);
	  IntRange vdofs = cfelu.GetRange(1);
	  IntRange wdofs = cfelu.GetRange(2);

	  FlatVector<> eluv(felu.GetNDof(), lh), shape(nd_vr, lh), elf(nd_vr,lh);

	  gfu->GetElementVector(dnumsu, eluv);	  
	  
	  IntegrationRule ir(eltrans.GetElementType(), felu.Order()*2);
	  FlatVector<> BubbleProjEl(nd_vr, lh);
	  BubbleProjEl = 0.0;
	  Array<int> vnums;
	  vnums = ma->GetElVertices(ei);
	  //ma->GetElVertices(ei, vnums);
	  int vertexdof =0;
	  for (int k = 0; k < 4; k++)
	    {
	      if(vnums[k] == i)
	  	{
	  	  BubbleProjEl[k] = 1;
	  	}
	    }
	  
	   elf = 0.0;
	   for(int k = 0; k< ir.GetNIP(); k++)
	    {
	      MappedIntegrationPoint<3,3> mip(ir[k], eltrans);

	      Vec<3> gradu;
	      Vec<3> gradv;
	      Vec<3> gradw;
	      DiffOpGradient<3>::Apply(h1fel, mip, eluv.Range(udofs), gradu, lh);
	      DiffOpGradient<3>::Apply(h1fel, mip, eluv.Range(vdofs), gradv, lh);
	      DiffOpGradient<3>::Apply(h1fel, mip, eluv.Range(wdofs), gradw, lh);
	      Vec<1> divu = gradu[0] + gradv[1] + gradw[2];			     

	      double fac = ir[k].Weight() * mip.GetMeasure();
	      l2fel.CalcShape(ir[k], shape);
	      elf += (fac * divu(0))*shape;
	    }
	   
	  FlatArray<int> dnumss_l2part = dnumss.Range(cfelur.GetRange(1));
	  
	  for(int k=0; k< dnumss_l2part.Size(); k++)
	    for (int l = 0; l<globaldofs[i].Size(); ++l)
	      if(dnumss_l2part[k] == globaldofs[i][l])
		{
		  localrhs[l] = elf[k];
		  BubbleProj[l] = BubbleProjEl[k];
		}
	}

      Vector<> vertexpatchrhs(ldofs+meanvaldofs);
      vertexpatchrhs = 0.0;
      vertexpatchrhs.Range(0,ldofs) = extrosc[i] * localrhs.Range(0,ldofs);
      //cout<<"Norm(vertexpatchrhs) = "<<L2Norm(vertexpatchrhs)<<endl;
      for (int l = 0; l<globaldofs[i].Size(); ++l)
      	vertexpatchrhs[l] *= BubbleProj[l];

      localrhs = patchmatinv[i] * vertexpatchrhs;
 #pragma omp critical    
      gfsigma->GetVector().AddIndirect(globaldofs[i], localrhs.Range(0,ldofs));
    }
}

void ReconstructionVertex3D::CalcReconstructionVertexTrans(shared_ptr<LinearForm> fsigma, shared_ptr<LinearForm> fu, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace)
//void ReconstructionVertex3D::CalcReconstructionVertexTrans(PyLF fsigma, PyLF fu, PyBF bfa, PyFES h1fespace)
{
  static Timer timer("Reconstruction::CalcReconstructionVertexTrans");
  static Timer timerh1("Reconstruction::H1toL2");
  static Timer timerinv("Reconstruction::inverse");
  static Timer timerdofs("Reconstruction::dofs");
  static Timer timersmooth("Reconstruction::smooth");
  static Timer timerapplytrans("Reconstruction::applytrans");
  RegionTimer reg (timer);
  //cout<<"CalcReconstructionVertexTrans"<<endl;
  fu->GetVector()=0.0;
  shared_ptr<MeshAccess> ma  = fu->GetMeshAccess();
  shared_ptr<FESpace> fespace_s = fsigma->GetFESpace();
  shared_ptr<FESpace> fespace_u = fu->GetFESpace();

  LocalHeap lhglobal(1000000000, "local_reconstruction");

  IntRange numberdofs = dynamic_pointer_cast<CompoundFESpace>(fespace_s)->GetRange(2);

  TableCreator<int> h1tol2creator(h1fespace->GetNDof());

  timerh1.Start();
  for( ; !h1tol2creator.Done(); h1tol2creator++)
    {
      for (int i = 0; i < ne; ++i)
	{
	  Array<int> elementdofsh1;
	  h1fespace->GetDofNrs(ElementId(VOL,i), elementdofsh1);      

	  Array<int> elementdofs;
	  HeapReset hr(lhglobal);

	  ElementId ei(VOL,i);
	  const FiniteElement & felur = fespace_s -> GetFE(ei,lhglobal);
	  const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
	  const ScalarFiniteElement<3> & l2fel = dynamic_cast<const ScalarFiniteElement<3> &> (cfelur[1]);
	  
	  IntRange l2dofs = cfelur.GetRange(1);
	  fespace_s->GetDofNrs(ElementId(VOL,i), elementdofs);
	  
	  for (int j = 0; j < l2fel.GetNDof(); j++ )
	    {
	      h1tol2creator.Add(elementdofsh1[j],elementdofs[l2dofs.First()+j]);	      
	    }	 
	}
    }
  Table<int> h1tol2table = h1tol2creator.MoveTable();
  timerh1.Stop();
  //#pragma omp parallel for  
  //for(int i = 0; i < nv; i++)
  ParallelFor(Range(nv),[&](int i)
    {
      LocalHeap lh=lhglobal.Split();
      HeapReset hr(lh);
      const SparseMatrix<double> & mata = dynamic_cast<SparseMatrix<double> &>(bfa->GetMatrix());
      //Array<int> & globallocaldofs = globaldofs[i];
      Array<int> globallocaldofs;
      timerdofs.Start();     
      for (int j = 0; j < patch_elements[i].Size() ; j++)
	{
	  Array<int> dofs;
	  fespace_s->GetInnerDofNrs(patch_elements[i][j], dofs);
	  globallocaldofs+=dofs;
	}
         
      for (int j = 0; j < patch_faces[i].Size() ; j++)
	{
	  Array<int> dofs;
	  fespace_s->GetFaceDofNrs(patch_faces[i][j], dofs);
	  globallocaldofs+=dofs;
	}

      globallocaldofs+=numberdofs;

      int ldofs = globallocaldofs.Size();
      QuickSort(globallocaldofs);
      
      FlatMatrix<> localmata(ldofs,lh);
      FlatMatrix<> bubbleprojmat(ldofs,lh);
      bubbleprojmat = Identity(ldofs);
      localmata=0.0;

      for (int j = 0; j < ldofs; ++j)
	{
	  //for (int k = 0; k < ldofs; ++k)
	  //  {
	  //    localmata(j,k) = 1;//mata(globallocaldofs[j],globallocaldofs[k]);
	  //  }
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
      timerdofs.Stop();
      timerinv.Start();
      CalcInverse(localmata);
      timerinv.Stop();
      
  FlatMatrix<> extrosclocal(ldofs,lh);
  //FlatVector<> smooth(ldofs,lh);
  timersmooth.Start(); 
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

      for(int j = 0; j < elementdofsh1.Size(); j++)
	{
	  Array<int> smoothdofs;
	  for (int k = 0; k < h1tol2table[elementdofsh1[j]].Size(); ++k)
	    {
	      //check in globallocaldofs
	      for(int s = 0; s<globallocaldofs.Size(); s++)
		if(globallocaldofs[s] == h1tol2table[elementdofsh1[j]][k])
		  {
		    //smooth(s)=1;
		    smoothdofs.Append(s);
		  }
	    }
	  //extrosclocal-=1.0/h1tol2table[elementdofsh1[j]].Size()*(smooth*Trans(smooth));
	  for (int k =0; k < smoothdofs.Size(); ++k)
	    {
	      for (int s =0; s < smoothdofs.Size(); ++s)
		{
		  extrosclocal(smoothdofs[k],smoothdofs[s]) -= 1.0/h1tol2table[elementdofsh1[j]].Size();
		}
	      
	    }
	}
      timersmooth.Stop();
      
      //##################################################
      //int ldofs = globallocaldofs;// globaldofs[i].Size();      
      FlatVector<> localrhs(ldofs,lh);
      FlatVector<> vertexpatchrhs(ldofs,lh);
      vertexpatchrhs = 0.0;
      
      FlatVector<> mean_sol(ldofs,lh);
      FlatVector<> sol(ldofs,lh);

      //fsigma->GetElementVector(globaldofs[i], vertexpatchrhs.Range(0,ldofs));
      fsigma->GetElementVector(globallocaldofs, vertexpatchrhs);
      
      //sol = 0.0;

      //mean_sol = patchmatinv[i]* vertexpatchrhs;
      mean_sol = localmata* vertexpatchrhs;
      //sol = Trans(extrosc[i]) * mean_sol.Range(0,ldofs);
      sol = extrosclocal* mean_sol;
     
      timerapplytrans.Start();

      for (int j = 0; j < patch_elements[i].Size() ; j++)
	{
	  HeapReset hr(lh);
	  int el = patch_elements[i][j];
	  ElementId ei(VOL,el);	
 
	  const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

	  const FiniteElement & felu = fespace_u -> GetFE(ei,lh);
	  const CompoundFiniteElement & cfelu = dynamic_cast<const CompoundFiniteElement&>(felu);
	  const ScalarFiniteElement<3> & h1fel = dynamic_cast<const ScalarFiniteElement<3> &> (cfelu[0]);

	  const FiniteElement & felur = fespace_s -> GetFE(ei,lh);	  
	  const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);	  	  
	  const ScalarFiniteElement<3> & l2fel = dynamic_cast<const ScalarFiniteElement<3> &> (cfelur[1]);	

	  int nd_vr = l2fel.GetNDof();

	  Array<int> dnumsu(felu.GetNDof(), lh), dnumss(felur.GetNDof(), lh);

	  fespace_u->GetDofNrs(ei, dnumsu);
	  fespace_s->GetDofNrs(ei, dnumss);

	  IntRange udofs = cfelu.GetRange(0);
	  IntRange vdofs = cfelu.GetRange(1);
	  IntRange wdofs = cfelu.GetRange(2);

	  FlatVector<>  shape(nd_vr, lh), elu(felu.GetNDof(),lh), elui(felu.GetNDof(),lh), elfsig(nd_vr, lh);
	  
	  elu = 0.0;
	  elfsig = 0.0;
	  
	  IntegrationRule ir(eltrans.GetElementType(), felu.Order()*2);
	  
	  FlatArray<int> dnumss_l2part = dnumss.Range(cfelur.GetRange(1));
	  FlatVector<> BubbleProjEl(nd_vr, lh);
	  BubbleProjEl = 0.0;
	  Array<int> vnums;
	  vnums = ma->GetElVertices(ei);
	  int vertexdof =0;
	  for (int k = 0; k < 4; k++)
	    {
	      if(vnums[k] == i)
	  	{
	  	  BubbleProjEl[k] = 1;
	  	}
	    }
	  
	  for(int k=0; k<dnumss_l2part.Size(); k++)
	     for (int l = 0; l<globallocaldofs.Size(); ++l)
	       if((dnumss_l2part[k] == globallocaldofs[l]))
		 elfsig[k] = sol[l]* BubbleProjEl[k];
	  	  
	  elui = 0; 
	  for(int k = 0; k< ir.GetNIP(); k++)
	    {
	      MappedIntegrationPoint<3,3> mip(ir[k], eltrans);
	      
	      l2fel.CalcShape(ir[k], shape);

	      double fac = InnerProduct(elfsig, shape);
	      fac *= ir[k].Weight() * mip.GetMeasure();

	      Vec<3> u1 = { fac, 0,0 }, u2={0,fac,0}, u3={0,0,fac};
	      
	      FlatVector<double> eluu=  elui.Range(udofs);
	      FlatVector<double> eluv=  elui.Range(vdofs);
	      FlatVector<double> eluw=  elui.Range(wdofs);

	      DiffOpGradient<3>::ApplyTrans(h1fel, mip, u1, eluu, lh);
	      DiffOpGradient<3>::ApplyTrans(h1fel, mip, u2, eluv, lh);
	      DiffOpGradient<3>::ApplyTrans(h1fel, mip, u3, eluw, lh);
	      elu += elui;
	    }

#pragma omp critical
	  fu->GetVector().AddIndirect(dnumsu, elu);
	}
      timerapplytrans.Stop();
    });
}


void ReconstructionVertex3D::CalcVertexPatches(shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini)
//void ReconstructionVertex3D::CalcVertexPatches(PyFES fespace_u, PyFES fespace_s, PyBF bfa, PyFES h1fespace)
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
  RegionTimer reg(timer);
  cout<<"CalcVertexPatches"<<endl;
  shared_ptr<MeshAccess> ma  = fespace_u->GetMeshAccess();
  
  BitArray finefaces(ma->GetNFaces());
  finefaces.Clear();

  IntRange numberdofs = dynamic_pointer_cast<CompoundFESpace>(fespace_s)->GetRange(2); 

  int l2order = (*(dynamic_pointer_cast<CompoundFESpace>(fespace_s)).get())[1]->GetOrder();

  meanvalorder = 0;
  meanvaldofs = 0;

  if ((l2order>1))
    {
      throw Exception("High order reconstruction not implemented!");
      //meanvalorder =l2order-1;
      //L2HighOrderFE<ET_TRIG> meanvalfe(meanvalorder-1);
      //meanvaldofs = meanvalfe.GetNDof();      
    }
  
  nv = ma->GetNV();
  ne = ma->GetNE();
  
  nfaces = ma->GetNFaces();

  for (int i = 0; i < ne; ++i)
    {      
      Array<int> fnums;
      fnums = ma->GetElFaces(ElementId(VOL,i));
      for(int f : fnums)
	finefaces.Set(f);      
    }

  timertable.Start();
  TableCreator<int> creator(nv); 
  TableCreator<int> creatorfaces(nv); 

  
  for( ; !creator.Done(); creator++)
    {
      for (int i = 0; i < ne; ++i)
	{
	  Array<int> vnums;
	  vnums = ma->GetElVertices(ElementId(VOL,i));      
	  for (int j: vnums)
	    creator.Add(j,i);

	}      
    }
  
  //BitArray dirfaces(ma->GetNFaces());
  BitArray dirfaces(nfaces);
  dirfaces.Clear();

  for (int i = 0; i < ma->GetNSE(); ++i)
    {
      
      //dirfaces.Set(ma->GetSElFace(i));

      //facenr = ma->GetSElFace(i);
      //Array<int> enums;
      //     ma->GetSElFace(i,enums);
      //if ((*dynamic_pointer_cast<CompoundFESpace>(fespace_u))[0]->IsDirichletBoundary(ma->GetSElIndex(i)))
      if ((*dynamic_pointer_cast<CompoundFESpace>(fespace_u))[0]->IsDirichletBoundary(ma->GetElIndex(ElementId(BND,i))))
	dirfaces.Set(ma->GetSElFace(i));     
    }


  for (int i = 0; i < finefaces.Size(); ++i)
    {
      if(!finefaces[i])
	dirfaces.Set(i);
    }

  //cout<<"dirfaces = "<<dirfaces<<endl;

    
  for( ; !creatorfaces.Done(); creatorfaces++)
    {
      for (int i = 0; i < nfaces; ++i)
	{
	  if(!dirfaces[i])
	    {
	      //int v1,v2,v3;
	      Array<int> vnums;
	      vnums=ma->GetFacePNums(i);
	      
	      creatorfaces.Add(vnums[0],i);
	      creatorfaces.Add(vnums[1],i);
	      creatorfaces.Add(vnums[2],i);
	    }
	}      
    }
  patch_elements = creator.MoveTable();
  patch_faces = creatorfaces.MoveTable();

  //cout<<"patch_elements = "<<patch_elements<<endl;
  //cout<<"patch_faces = "<<patch_faces<<endl;

  timertable.Stop();
  /*const SparseMatrix<double> & mata = dynamic_cast<SparseMatrix<double> &>(bfa->GetMatrix());

  timerh1.Start();
  TableCreator<int> h1tol2creator(h1fespace->GetNDof());
  LocalHeap lh(1000000, "local_reconstruction");

  for( ; !h1tol2creator.Done(); h1tol2creator++)
    {
      for (int i = 0; i < ne; ++i)
	{
	  Array<int> elementdofsh1;
	  h1fespace->GetDofNrs(i, elementdofsh1);      

	  Array<int> elementdofs;
	  HeapReset hr(lh);

	  ElementId ei(VOL,i);
	  const FiniteElement & felur = fespace_s -> GetFE(ei,lh);
	  const CompoundFiniteElement & cfelur = dynamic_cast<const CompoundFiniteElement&>(felur);
	  const ScalarFiniteElement<3> & l2fel = dynamic_cast<const ScalarFiniteElement<3> &> (cfelur[1]);
	  
	  IntRange l2dofs = cfelur.GetRange(1);
	  fespace_s->GetDofNrs(i, elementdofs);
	  
	  for (int j = 0; j < l2fel.GetNDof(); j++ )
	    {
	      h1tol2creator.Add(elementdofsh1[j],elementdofs[l2dofs.First()+j]);	      
	    }	 
	}
	}*/
  // Table<int> h1tol2table = h1tol2creator.MoveTable();
  //timerh1.Stop();
  
  //globaldofs.SetSize(nv);
  //patchmatinv.SetSize(nv);
  //extrosc.SetSize(nv);
  /* 
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
         
      for (int j = 0; j < patch_faces[i].Size() ; j++)
	{
	  Array<int> dofs;
	  fespace_s->GetFaceDofNrs(patch_faces[i][j], dofs);
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
	  //for (int k = 0; k < ldofs; ++k)
	  //  {
	  //    localmata(j,k) = 1;//mata(globallocaldofs[j],globallocaldofs[k]);
	  //  }
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
      //Vec<2> Vertex0(0,0);
      //ma->GetPoint(i,Vertex0);

      IntegrationRule bubbleweightpoints;
      for (int l = 0; l < l2order+1; ++l)
	{
	  for (int k = 0; k < l2order+1-l; ++k)
	    {
	      for (int j = 0; j < l2order+1-l-k; ++j)
		{
		  IntegrationPoint ip(double(l)/(l2order), double(k)/(l2order), double(j)/l2order);
		  bubbleweightpoints.AddIntegrationPoint(ip);
		}
	    }	  
	}
      
      ////for (int j = 0; j < patch_elements[i].Size(); ++j)
      //{	  	  
	//shapemeanval = 0.0;
	  
	  //timermeanvals.Start();
	  //f(meanvalorder !=0)
	  // {
	  //   for (int l = 0; l < ir.GetNIP(); l++)
	  //	{
	  //	  const MappedIntegrationPoint<2,2> mip(ir[l],eltrans);
	  //	  double det = mip.GetJacobiDet();
	  //	  IntegrationPoint physip(mip(0),mip(1));
	  //	  
	  //	  hdivfel.CalcMappedShape(mip, shape);
	  //	  meanvalfe.CalcShape(physip, shapemeanval);
	  //	  
	  //	  Vec<2> xycoord(mip(1)-Vertex0(1),-1*mip(0)-Vertex0(0));
	  //	  
	  //	  localmeanval += fabs(det) *ir[l].Weight()*shapemeanval*Trans(shape* xycoord);
	  //	}
	  //   Array<int>  dnumss(felur.GetNDof(), lh);
	  //   fespace_s->GetDofNrs(id, dnumss);	  	      
	  //   FlatArray<int> dnums_sigma = dnumss.Range(cfelur.GetRange(0));
	  //
	  //   for(int k = 0; k < dnums_sigma.Size(); k++)
	  //	for(int l = 0; l<globaldofs[i].Size(); l++)		
	  //	  if(dnums_sigma[k] == globaldofs[i][l])
	  //	    {
	  //	      if(meanvalorder !=0)
	  //		{
	  //		  localmata.Rows(ldofs,ldofs+meanvaldofs).Col(l) += localmeanval.Col(k);
	  //		  localmata.Cols(ldofs,ldofs+meanvaldofs).Row(l) += localmeanval.Col(k);
	  //		}
	  //	    }
	  // }
	  //timermeanvals.Stop();
	  /*
	  timerbubble.Start();
	  //Bubble Projection
	  Array<int> vnums;
	  ma->GetElVertices(ei, vnums);
	  int vertexdof =0;
	  for (int k = 0; k < 4; k++)
	    {
	      if(vnums[k] == i)
		{
		  vertexdof = k;
		}
	    }
	  
	  FlatMatrix<> basistransformation(nd_vr,nd_vr, lh);
	  FlatMatrix<> invbasistransformation(nd_vr,nd_vr, lh);
	  FlatMatrix<> shape_hats(4, nd_vr, lh);
	  
	  l2fel.CalcShape(bubbleweightpoints, basistransformation);
	  ScalarFE<ET_TET,1> p1fe;

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
	  fespace_s->GetDofNrs(id, dnumss);	  	      
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
	  h1fespace->GetDofNrs(id, elementdofselement);
	  
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
      //extrosc[i] = (bubbleprojmat*extrosclocal | Lapack) ;
      extrosc[i] = extrosclocal;
      timermult.Stop();
      //timerinv.AddFlops(bubbleprojmat.Width()*bubbleprojmat.Width()*bubbleprojmat.Width());
      
    }*/
  
}
