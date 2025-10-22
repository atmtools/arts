#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_vector.h"
#include "mystring.h"
#include "rtepack.h"
#include "sun.h"

namespace Python {
void py_star(py::module_& m) try {
  py::class_<Sun> suns(m, "Sun");
  generic_interface(suns);
  suns.def_rw("description",
              &Sun::description,
              "Sun description\n\n.. :class:`~pyarts3.arts.String`")
      .def_prop_rw(
          "spectrum",
          [](Sun& self) -> Matrix& { return self.spectrum; },
          [](Sun& self, const std::variant<Matrix, Vector, StokvecVector>& value) {
            if (std::holds_alternative<Matrix>(value)) {
              self.spectrum = std::get<Matrix>(value);
            } else if (std::holds_alternative<Vector>(value)) {
              const auto& v = std::get<Vector>(value);
              self.spectrum.resize(v.size(), 4);
              self.spectrum = 0.0;
              self.spectrum[joker, 0] = v;
            } else if (std::holds_alternative<StokvecVector>(value)) {
              const auto& sv = std::get<StokvecVector>(value);
              self.spectrum.resize(sv.size(), 4);
              for (Size i = 0; i < sv.size(); ++i) {
                self.spectrum[i, 0] = sv[i].I();
                self.spectrum[i, 1] = sv[i].Q();
                self.spectrum[i, 2] = sv[i].U();
                self.spectrum[i, 3] = sv[i].V();
              }
            }
          },
          "Sun spectrum, monochromatic radiance spectrum at the surface of the sun\n\nAccepts :class:`~pyarts3.arts.Matrix`, :class:`~pyarts3.arts.Vector` (sets first Stokes), or :class:`~pyarts3.arts.StokvecVector`\n\n.. :class:`~pyarts3.arts.Matrix`")
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
             new (self) Sun{.description = std::get<0>(state),
                            .spectrum    = std::get<1>(state),
                            .radius      = std::get<2>(state),
                            .distance    = std::get<3>(state),
                            .latitude    = std::get<4>(state),
                            .longitude   = std::get<5>(state)};
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
