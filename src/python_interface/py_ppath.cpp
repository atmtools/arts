#include <py_auto_interface.h>

#include "py_macros.h"

#include <ppath_struct.h>

namespace Python {
void py_ppath(py::module_& m) {
  py::class_<GridPos>(m, "GridPos")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(GridPos)
      .def(py::init<Index, std::array<Numeric, 2>>())
      .PythonInterfaceFileIO(GridPos)
      .PythonInterfaceBasicRepresentation(GridPos)
      .def_readwrite("idx", &GridPos::idx)
      .def_readwrite("fd", &GridPos::fd);

  py::class_<ArrayOfGridPos>(m, "ArrayOfGridPos")
      .PythonInterfaceArrayDefault(GridPos)
      .PythonInterfaceBasicRepresentation(ArrayOfGridPos);
  py::implicitly_convertible<std::vector<GridPos>, ArrayOfGridPos>();

  py::class_<Ppath>(m, "Ppath")
      .def(py::init<>())
      .def(py::init<Index,
                    Index,
                    Numeric,
                    String,
                    Vector,
                    Vector,
                    Numeric,
                    Matrix,
                    Matrix,
                    Vector,
                    Vector,
                    Vector,
                    Vector,
                    Numeric,
                    Vector,
                    Vector,
                    ArrayOfGridPos,
                    ArrayOfGridPos,
                    ArrayOfGridPos>())
      .PythonInterfaceWorkspaceVariableConversion(Ppath)
      .def("__repr__", [](Ppath&){return "Ppath";})
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
      .def_readwrite("gp_lon", &Ppath::gp_lon);

  py::class_<ArrayOfPpath>(m, "ArrayOfPpath")
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfPpath)
      .PythonInterfaceFileIO(ArrayOfPpath)
      .def("__repr__", [](ArrayOfPpath&){return "ArrayOfPpath";})
      .PythonInterfaceArrayDefault(Ppath);
  py::implicitly_convertible<std::vector<Ppath>, ArrayOfPpath>();
}
}  // namespace Python