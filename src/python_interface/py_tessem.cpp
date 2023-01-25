#include <py_auto_interface.h>

#include "matpackI.h"
#include "py_macros.h"

namespace Python {
void py_tessem(py::module_& m) {
  py::class_<TessemNN>(m, "TessemNN")
      .def(py::init([]() { return std::make_unique<TessemNN>(); }))
      .PythonInterfaceCopyValue(TessemNN)
      .PythonInterfaceWorkspaceVariableConversion(TessemNN)
      .PythonInterfaceFileIO(TessemNN)
      .def("__repr__", [](TessemNN&) { return "TessemNN"; })
      .def_readwrite("nb_inputs", &TessemNN::nb_inputs)
      .def_readwrite("nb_outputs", &TessemNN::nb_outputs)
      .def_readwrite("nb_cache", &TessemNN::nb_cache)
      .def_readwrite("b1", &TessemNN::b1)
      .def_readwrite("b2", &TessemNN::b2)
      .def_readwrite("w1", &TessemNN::w1)
      .def_readwrite("w2", &TessemNN::w2)
      .def_readwrite("x_min", &TessemNN::x_min)
      .def_readwrite("x_max", &TessemNN::x_max)
      .def_readwrite("y_min", &TessemNN::y_min)
      .def_readwrite("y_max", &TessemNN::y_max)
      .def(py::pickle(
          [](const TessemNN& self) {
            return py::make_tuple(self.nb_inputs,
                                  self.nb_outputs,
                                  self.nb_cache,
                                  self.b1,
                                  self.b2,
                                  self.w1,
                                  self.w2,
                                  self.x_min,
                                  self.x_max,
                                  self.y_min,
                                  self.y_max);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 11, "Invalid state!")
            return std::make_unique<TessemNN>(TessemNN{t[0].cast<Index>(),
                                t[1].cast<Index>(),
                                t[2].cast<Index>(),
                                t[3].cast<Vector>(),
                                t[4].cast<Vector>(),
                                t[5].cast<Matrix>(),
                                t[6].cast<Matrix>(),
                                t[7].cast<Vector>(),
                                t[8].cast<Vector>(),
                                t[9].cast<Vector>(),
                                t[10].cast<Vector>()});
          }))
      .PythonInterfaceWorkspaceDocumentation(TessemNN);
}
}  // namespace Python