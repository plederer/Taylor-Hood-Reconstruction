#include "reconstruction.hpp"
#include <python_ngstd.hpp>

//L2HightOrderFESpace with dofs on the edges/vertices due to the H1-FiniteElement
//GetFE is used to get the right dofs 
class L2HighOrderFESpaceH1Dofs : public L2HighOrderFESpace
{
public:    
  L2HighOrderFESpaceH1Dofs(shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false) : L2HighOrderFESpace(ama,flags,parseflags){};
  
  string GetClassName() const
  {
    return "L2HighOrderFESpace with H1 dof coupling";
  }

  virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const
  {
    int elnr = ei.Nr();
    Ngs_Element ngel = ma->GetElement(ei);
    if (!DefinedOn (ElementId (VOL, elnr)))
      {
	throw Exception ("BAD ELEMENT");
      }
    try
      {
	if(ma->GetDimension()==2)
	  {
	    H1HighOrderFE<ET_TRIG> * hofe =  new (lh) H1HighOrderFE<ET_TRIG> (order_inner[elnr][0]);
	    hofe -> SetVertexNumbers (ngel.Vertices());
	    return *hofe;
	  }
	else
	  {
	    H1HighOrderFE<ET_TET> * hofe =  new (lh) H1HighOrderFE<ET_TET> (order_inner[elnr][0]);
	    hofe -> SetVertexNumbers (ngel.Vertices());
	    return *hofe;
	  }

      }
    catch (Exception & e)
      {
	e.Append ("in L2HoFESpaceH1Element::GetElement\n");
	throw;
      }
  }
};

static RegisterFESpace<L2HighOrderFESpaceH1Dofs> initl2hoh1space("l2hoh1dofs");

void ExportRec(py::module &m) 
{

  py::class_<ReconstructionElement>(m, "ReconstructionElement")
    .def(py::init<>())
    .def("Setup", [] (ReconstructionElement & self, shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini) {
	self.CalcElementPatches(fespace_u,fespace_s,bfa, h1fespace, mini);
      }, "Setup for reconstruction", py::arg("fespace_u"), py::arg("fespace_s"), py::arg("bfa"), py::arg("h1fespace"), py::arg("mini")=false)

    
    .def("Calc", &ReconstructionElement::CalcReconstructionElement)
    .def("CalcTrans", &ReconstructionElement::CalcReconstructionElementTrans)
    ;
  

  py::class_<ReconstructionVertex>(m, "ReconstructionVertex")
    .def(py::init<>())
    .def("Setup", [] (ReconstructionVertex & self, shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini) {
	self.CalcVertexPatches(fespace_u,fespace_s,bfa, h1fespace, mini);
      }, "Setup for reconstruction", py::arg("fespace_u"), py::arg("fespace_s"), py::arg("bfa"), py::arg("h1fespace"), py::arg("mini")=false)

    
    .def("Calc", &ReconstructionVertex::CalcReconstructionVertex)
    .def("CalcTrans", &ReconstructionVertex::CalcReconstructionVertexTrans)
    ;

  py::class_<ReconstructionVertex3D>(m, "ReconstructionVertex3D")
    .def(py::init<>())
    .def("Setup", [] (ReconstructionVertex3D & self, shared_ptr<FESpace> fespace_u, shared_ptr<FESpace> fespace_s, shared_ptr<BilinearForm> bfa, shared_ptr<FESpace> h1fespace, bool mini) {
	self.CalcVertexPatches(fespace_u,fespace_s,bfa, h1fespace, mini);
	 }, "Setup for reconstruction", py::arg("fespace_u"), py::arg("fespace_s"), py::arg("bfa"), py::arg("h1fespace"), py::arg("mini")=false)

    .def("Calc", &ReconstructionVertex3D::CalcReconstructionVertex)
    .def("CalcTrans", &ReconstructionVertex3D::CalcReconstructionVertexTrans)
    ;

}


PYBIND11_MODULE(reclib,m) 
{
  ExportRec(m);

}
