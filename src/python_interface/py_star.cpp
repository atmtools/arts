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
              ":class:`~pyarts3.arts.String` Sun description")
      .def_rw(
          "spectrum",
          &Sun::spectrum,
          ":class:`~pyarts3.arts.Matrix` Sun spectrum, monochrmatic radiance spectrum at the surface of the sun")
      .def_rw("radius", &Sun::radius, ":class:`float` Sun radius")
      .def_rw("distance", &Sun::distance, ":class:`float` Sun distance")
      .def_rw("latitude", &Sun::latitude, ":class:`float` Sun latitude")
      .def_rw("longitude", &Sun::longitude, ":class:`float` Sun longitude")
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
