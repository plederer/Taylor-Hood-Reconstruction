#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>  

using namespace ngsolve;
using namespace ngcomp;
using namespace ngfem;

class ReconstructionElement
{
 private:
  Array<Matrix<double>> patchmatinv;
  Array<Matrix<double>> extrosc;
  Array<Array<int>> globaldofs;

  int ne;
  int nv;
  int nedges;
  int meanvalorder;
  int meanvaldofs;
  
 public:
  ReconstructionElement()
    {      
      ne = 0;
      nv = 0;
      nedges = 0;
      meanvalorder =0;
      meanvaldofs = 0.0;
    }

  void CalcElementPatches(shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini);  
  void CalcReconstructionElement(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfsigma);
  void CalcReconstructionElementTrans(shared_ptr<LinearForm> fsigma, shared_ptr<LinearForm> fu);
};

class ReconstructionVertex
{
 private:
  Array<Matrix<double>> patchmatinv;
  Array<Matrix<double>> extrosc;
  Array<Array<int>> globaldofs;

  Table<int> patch_elements;
  Table<int> patch_edges;
  
  int ne;
  int nv;
  int nedges;
  int meanvalorder;
  int meanvaldofs;
  
 public:
  ReconstructionVertex()
    {      
      ne = 0;
      nv = 0;
      nedges = 0;
      meanvalorder =0;
      meanvaldofs = 0.0;      
    }

  void CalcVertexPatches(shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini);
  void CalcReconstructionVertex(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfsigma);
  void CalcReconstructionVertexTrans(shared_ptr<LinearForm> fsigma, shared_ptr<LinearForm> fu);

};

class ReconstructionVertex3D
{
 private:
  Array<Matrix<double>> patchmatinv;
  Array<Matrix<double>> extrosc;
  Array<Array<int>> globaldofs;

  Table<int> patch_elements;
  Table<int> patch_faces;
  
  int ne;
  int nv;
  int nfaces;
  int meanvalorder;
  int meanvaldofs;
  
 public:
  ReconstructionVertex3D()
    {      
      ne = 0;
      nv = 0;
      nfaces = 0;
      meanvalorder =0;
      meanvaldofs = 0.0;      
    }

  void CalcVertexPatches(shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini);  
  void CalcReconstructionVertex(shared_ptr<GridFunction> gfu, shared_ptr<GridFunction> gfsigma);
  void CalcReconstructionVertexTrans(shared_ptr<LinearForm> fsigma, shared_ptr<LinearForm> fu, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace);
};

#endif
