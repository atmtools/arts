#include <auto_md.h>
#include <xml_io.h>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_tessem(py::module_& m) {
  py::class_<TessemNN>(m, "TessemNN")
      .def(py::init<>())
      .PythonInterfaceFileIO(TessemNN)
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
      .def_readwrite("y_max", &TessemNN::y_max);
}
}  // namespace Python