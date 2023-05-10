#include <py_auto_interface.h>

#include "py_macros.h"

#include <ppath_struct.h>

namespace Python {
void py_ppath(py::module_& m) {
  py::class_<GridPos>(m, "GridPos")
      .def(py::init([]() { return std::make_unique<GridPos>(); }))
      .PythonInterfaceCopyValue(GridPos)
      .PythonInterfaceWorkspaceVariableConversion(GridPos)
      .def(py::init([](Index i, std::array<Numeric, 2> x) { return std::make_unique<GridPos>(GridPos{i, x}); }))
      .PythonInterfaceFileIO(GridPos)
      .PythonInterfaceBasicRepresentation(GridPos)
      .def_readwrite("idx", &GridPos::idx)
      .def_readwrite("fd", &GridPos::fd)
      .def(py::pickle(
          [](const GridPos& self) {
            return py::make_tuple(self.idx, self.fd);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")
            return std::make_unique<GridPos>(GridPos{t[0].cast<Index>(), t[1].cast<std::array<Numeric, 2>>()});
          }))
      .PythonInterfaceWorkspaceDocumentation(GridPos);

  py::class_<ArrayOfGridPos>(m, "ArrayOfGridPos")
      .PythonInterfaceArrayDefault(GridPos)
      .PythonInterfaceBasicRepresentation(ArrayOfGridPos);
  py::implicitly_convertible<std::vector<GridPos>, ArrayOfGridPos>();

  py::class_<Ppath>(m, "Ppath")
      .def(py::init([]() { return std::make_unique<Ppath>(); }))
      .PythonInterfaceCopyValue(Ppath)
      .PythonInterfaceWorkspaceVariableConversion(Ppath)
      .def("__repr__", [](Ppath&) { return "Ppath"; })
      .PythonInterfaceFileIO(Ppath)
      .def_readwrite("dim", &Ppath::dim)
      .def_readwrite("np", &Ppath::np)
      .def_readwrite("constant", &Ppath::constant)
      .def_readwrite("background", &Ppath::background)
      .def_readwrite("start_pos", &Ppath::start_pos)
      .def_readwrite("start_los", &Ppath::start_los)
      .def_readwrite("start_lstep", &Ppath::start_lstep)
      .def_readwrite("pos", &Ppath::pos)
      .def_readwrite("los", &Ppath::los)
      .def_readwrite("r", &Ppath::r)
      .def_readwrite("lstep", &Ppath::lstep)
      .def_readwrite("end_pos", &Ppath::end_pos)
      .def_readwrite("end_los", &Ppath::end_los)
      .def_readwrite("end_lstep", &Ppath::end_lstep)
      .def_readwrite("nreal", &Ppath::nreal)
      .def_readwrite("ngroup", &Ppath::ngroup)
      .def_readwrite("gp_p", &Ppath::gp_p)
      .def_readwrite("gp_lat", &Ppath::gp_lat)
      .def_readwrite("gp_lon", &Ppath::gp_lon)
      .def(py::pickle(
          [](const Ppath& self) {
            return py::make_tuple(self.dim,
                                  self.np,
                                  self.constant,
                                  self.background,
                                  self.start_pos,
                                  self.start_los,
                                  self.start_lstep,
                                  self.pos,
                                  self.los,
                                  self.r,
                                  self.lstep,
                                  self.end_pos,
                                  self.end_los,
                                  self.end_lstep,
                                  self.nreal,
                                  self.ngroup,
                                  self.gp_p,
                                  self.gp_lat,
                                  self.gp_lon);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 19, "Invalid state!")
            return new Ppath{t[0].cast<Index>(),
                             t[1].cast<Index>(),
                             t[2].cast<Numeric>(),
                             t[3].cast<PpathBackground>(),
                             t[4].cast<Vector>(),
                             t[5].cast<Vector>(),
                             t[6].cast<Numeric>(),
                             t[7].cast<Matrix>(),
                             t[8].cast<Matrix>(),
                             t[9].cast<Vector>(),
                             t[10].cast<Vector>(),
                             t[11].cast<Vector>(),
                             t[12].cast<Vector>(),
                             t[13].cast<Numeric>(),
                             t[14].cast<Vector>(),
                             t[15].cast<Vector>(),
                             t[16].cast<ArrayOfGridPos>(),
                             t[17].cast<ArrayOfGridPos>(),
                             t[18].cast<ArrayOfGridPos>()};
          }))
      .PythonInterfaceWorkspaceDocumentation(Ppath);

  py::class_<ArrayOfPpath>(m, "ArrayOfPpath")
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfPpath)
      .PythonInterfaceFileIO(ArrayOfPpath)
      .def("__repr__", [](ArrayOfPpath&){return "ArrayOfPpath";})
      .PythonInterfaceArrayDefault(Ppath)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfPpath);
  py::implicitly_convertible<std::vector<Ppath>, ArrayOfPpath>();
}
}  // namespace Python
