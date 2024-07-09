#include <python_interface.h>

#include "hpy_arts.h"

namespace Python {
void py_tessem(py::module_& m) try {
  py::class_<TessemNN> tess(m, "TessemNN");
  workspace_group_interface(tess);
  tess.def("__repr__", [](TessemNN&) { return "TessemNN"; })
      .def_rw("nb_inputs", &TessemNN::nb_inputs, ":class:`int`")
      .def_rw("nb_outputs", &TessemNN::nb_outputs, ":class:`int`")
      .def_rw("nb_cache", &TessemNN::nb_cache, ":class:`int`")
      .def_rw("b1", &TessemNN::b1, ":class:`~pyarts.arts.Vector`")
      .def_rw("b2", &TessemNN::b2, ":class:`~pyarts.arts.Vector`")
      .def_rw("w1", &TessemNN::w1, ":class:`~pyarts.arts.Matrix`")
      .def_rw("w2", &TessemNN::w2, ":class:`~pyarts.arts.Matrix`")
      .def_rw("x_min", &TessemNN::x_min, ":class:`~pyarts.arts.Vector`")
      .def_rw("x_max", &TessemNN::x_max, ":class:`~pyarts.arts.Vector`")
      .def_rw("y_min", &TessemNN::y_min, ":class:`~pyarts.arts.Vector`")
      .def_rw("y_max", &TessemNN::y_max, ":class:`~pyarts.arts.Vector`")
      .def("__getstate__",
           [](const TessemNN& self) {
             return std::tuple<Index,
                               Index,
                               Index,
                               Vector,
                               Vector,
                               Matrix,
                               Matrix,
                               Vector,
                               Vector,
                               Vector,
                               Vector>{self.nb_inputs,
                                       self.nb_outputs,
                                       self.nb_cache,
                                       self.b1,
                                       self.b2,
                                       self.w1,
                                       self.w2,
                                       self.x_min,
                                       self.x_max,
                                       self.y_min,
                                       self.y_max};
           })
      .def("__setstate__",
           [](TessemNN* self,
              const std::tuple<Index,
                               Index,
                               Index,
                               Vector,
                               Vector,
                               Matrix,
                               Matrix,
                               Vector,
                               Vector,
                               Vector,
                               Vector>& state) {
             new (self) TessemNN{std::get<0>(state),
                                 std::get<1>(state),
                                 std::get<2>(state),
                                 std::get<3>(state),
                                 std::get<4>(state),
                                 std::get<5>(state),
                                 std::get<6>(state),
                                 std::get<7>(state),
                                 std::get<8>(state),
                                 std::get<9>(state),
                                 std::get<10>(state)};
           });
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize tessem\n", e.what()));
}
}  // namespace Python