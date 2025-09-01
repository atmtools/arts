#include <nanobind/stl/bind_vector.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_vector.h"
#include "mystring.h"
#include "sun.h"

namespace Python {
void py_star(py::module_& m) try {
  py::class_<Sun> suns(m, "Sun");
  generic_interface(suns);
  suns.def_rw("description",
              &Sun::description,
              "Sun description\n\n.. :class:`~pyarts3.arts.String`")
      .def_rw(
          "spectrum",
          &Sun::spectrum,
          "Sun spectrum, monochromatic radiance spectrum at the surface of the sun\n\n.. :class:`~pyarts3.arts.Matrix`")
      .def_rw("radius", &Sun::radius, "Sun radius\n\n.. :class:`float`")
      .def_rw("distance", &Sun::distance, "Sun distance\n\n.. :class:`float`")
      .def_rw("latitude", &Sun::latitude, "Sun latitude\n\n.. :class:`float`")
      .def_rw(
          "longitude", &Sun::longitude, "Sun longitude\n\n.. :class:`float`")
      .def("__getstate__",
           [](const Sun& self) {
             return std::
                 tuple<String, Matrix, Numeric, Numeric, Numeric, Numeric>{
                     self.description,
                     self.spectrum,
                     self.radius,
                     self.distance,
                     self.latitude,
                     self.longitude};
           })
      .def("__setstate__",
           [](Sun* self,
              const std::
                  tuple<String, Matrix, Numeric, Numeric, Numeric, Numeric>&
                      state) {
             new (self) Sun{std::get<0>(state),
                            std::get<1>(state),
                            std::get<2>(state),
                            std::get<3>(state),
                            std::get<4>(state),
                            std::get<5>(state)};
           });

  auto a1 = py::bind_vector<ArrayOfSun, py::rv_policy::reference_internal>(
      m, "ArrayOfSun");
  generic_interface(a1);
  vector_interface(a1);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize star\n{}", e.what()));
}
}  // namespace Python
