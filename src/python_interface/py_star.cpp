#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_vector.h"
#include "mystring.h"
#include "rtepack.h"
#include "sun.h"
#include <sun_methods.h>

namespace Python {
void py_star(py::module_& m) try {
  py::class_<Sun> suns(m, "Sun");
  generic_interface(suns);
  suns.def_rw("description", &Sun::description, "Sun description\n\n.. :class:`~pyarts3.arts.String`")
      .def_prop_rw(
          "spectrum",
          [](Sun& self) -> Matrix& { return self.spectrum; },
          [](Sun& self, const std::variant<Matrix, Vector, StokvecVector>& value) {
            if (std::holds_alternative<Matrix>(value)) {
              self.spectrum = std::get<Matrix>(value);
            } else if (std::holds_alternative<Vector>(value)) {
              const auto& v = std::get<Vector>(value);
              self.spectrum.resize(v.size(), 4);
              self.spectrum           = 0.0;
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
      .def_rw("longitude", &Sun::longitude, "Sun longitude\n\n.. :class:`float`");

  auto a1 = py::bind_vector<ArrayOfSun, py::rv_policy::reference_internal>(m, "ArrayOfSun");
  generic_interface(a1);
  vector_interface(a1);

  auto sun = m.def_submodule("sun");
  sun.doc() = "Contains helper functions to deal with suns";

  sun.def("geometric_los",
          &sun_geometric_los,
          "sun"_a,
          "pos"_a,
          "surf_field"_a,
          R"(Geometric line-of-sight from an observer towards a sun (no refraction)

The sun is placed at geodetic position
``[sun.distance - surf_field["h"](lat, lon), sun.latitude, sun.longitude]``
relative to the observer's latitude and longitude, and the line-of-sight is
the straight-line direction from the observer to that point.

Parameters
----------
sun : ~pyarts3.arts.Sun
  The sun
pos : ~pyarts3.arts.Vector3
  Observer position [alt, lat, lon]
surf_field : ~pyarts3.arts.SurfaceField
  Surface field providing the surface altitude and reference ellipsoid

Returns
-------
los : ~pyarts3.arts.Vector2
  Line-of-sight [zenith, azimuth] in degrees, local ENU at the observer
)");

  sun.def("refractive_los",
          &sun_refractive_los,
          "ws"_a,
          "sun"_a,
          "pos"_a,
          "surf_field"_a,
          "ray_path_observer_agenda"_a,
          py::kw_only(),
          "angle_cut"_a = 0.0,
          "refinement"_a = 1,
          R"(Refraction-aware line-of-sight from an observer towards a sun

Runs the iterative sun-path search (as *sun_pathFromObserverAgenda* with
``just_hit = 1``) and returns the observer line-of-sight of the resulting
path.  With a geometric *ray_path_observer_agenda* this reduces to the
geometric line-of-sight.

Parameters
----------
ws : ~pyarts3.arts.Workspace
  The workspace
sun : ~pyarts3.arts.Sun
  The sun
pos : ~pyarts3.arts.Vector3
  Observer position [alt, lat, lon]
surf_field : ~pyarts3.arts.SurfaceField
  Surface field providing the surface altitude and reference ellipsoid
ray_path_observer_agenda : ~pyarts3.arts.Agenda
  Agenda for tracing observer paths (refraction-aware)
angle_cut : float
  The angle delta-cutoff in the iterative solver [0.0, ...]
refinement : int
  The refinement of the search algorithm (twice the power of this is the resolution)

Returns
-------
los : ~pyarts3.arts.Vector2
  Line-of-sight [zenith, azimuth] in degrees, local ENU at the observer
)");
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize star\n{}", e.what()));
}
}  // namespace Python
