#include <py_auto_interface.h>

#include "py_macros.h"


namespace Python {
void py_star(py::module_& m) {
  py::class_<Star>(m, "Star")
      .def(py::init([]() { return new Star{}; }))
      .PythonInterfaceCopyValue(Star)
//      .PythonInterfaceWorkspaceVariableConversion(Star)
      .def("__repr__", [](Star&) { return "Star"; })
      .PythonInterfaceFileIO(Star)
      .def_readwrite("description", &Star::description)
      .def_readwrite("spectrum", &Star::spectrum)
      .def_readwrite("radius", &Star::radius)
      .def_readwrite("distance", &Star::distance)
      .def_readwrite("latitude", &Star::latitude)
      .def_readwrite("longitude", &Star::longitude)
      .def(py::pickle(
          [](const Star& self) {
            return py::make_tuple(self.description,
                                  self.spectrum,
                                  self.radius,
                                  self.distance,
                                  self.latitude,
                                  self.longitude);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 6, "Invalid state!")
            return new Star{t[0].cast<String>(),
                             t[1].cast<Matrix>(),
                             t[2].cast<Numeric>(),
                             t[3].cast<Numeric>(),
                             t[4].cast<Numeric>(),
                             t[5].cast<Numeric>()};
          }));

  py::class_<ArrayOfStar>(m, "ArrayOfStar")
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfStar)
      .PythonInterfaceFileIO(ArrayOfStar)
      .def("__repr__", [](ArrayOfStar&) { return "ArrayOfStar"; })
      .PythonInterfaceArrayDefault(Star);
  py::implicitly_convertible<std::vector<Star>, ArrayOfStar>();
}
}  // namespace Python