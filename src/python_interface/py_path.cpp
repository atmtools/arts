#include <nanobind/stl/bind_vector.h>
#include <path_point.h>
#include <python_interface.h>

#include <algorithm>

#include "hpy_arts.h"
#include "hpy_vector.h"

namespace Python {
void py_path(py::module_& m) try {
  py::class_<PropagationPathPoint> pppp(m, "PropagationPathPoint");
  generic_interface(pppp);

  pppp.def_rw("pos_type",
              &PropagationPathPoint::pos_type,
              "Path position type\n\n.. :class:`PathPositionType`")
      .def_rw("los_type",
              &PropagationPathPoint::los_type,
              "Path line-of-sight type\n\n.. :class:`PathPositionType`")
      .def_rw("pos",
              &PropagationPathPoint::pos,
              "Path position\n\n.. :class:`~pyarts3.arts.Vector3`")
      .def_rw("los",
              &PropagationPathPoint::los,
              "Path line-of-sight\n\n.. :class:`~pyarts3.arts.Vector2`")
      .def_rw("nreal",
              &PropagationPathPoint::nreal,
              "Path real refractive index\n\n.. :class:`Numeric`")
      .def_rw("ngroup",
              &PropagationPathPoint::ngroup,
              "Path group refractive index\n\n.. :class:`Numeric`");

  auto a1 = py::bind_vector<ArrayOfPropagationPathPoint,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfPropagationPathPoint");
  a1.def(
      "distances",
      [](const ArrayOfPropagationPathPoint& path, const Vector2& ellipsoid) {
        return path::distance(path, ellipsoid);
      },
      "ellipsoid"_a,
      "Calculate distances between neighboring points in the path");
  a1.def(
      "altitudes",
      [](const ArrayOfPropagationPathPoint& path) {
        return Vector{std::from_range,
                      path | stdv::transform([](const PropagationPathPoint& p) {
                        return p.altitude();
                      })};
      },
      "Get all the altitudes in the path");
  a1.def(
      "latitudes",
      [](const ArrayOfPropagationPathPoint& path) {
        return Vector{std::from_range,
                      path | stdv::transform([](const PropagationPathPoint& p) {
                        return p.latitude();
                      })};
      },
      "Get all the latitudes in the path");
  a1.def(
      "longitudes",
      [](const ArrayOfPropagationPathPoint& path) {
        return Vector{std::from_range,
                      path | stdv::transform([](const PropagationPathPoint& p) {
                        return p.longitude();
                      })};
      },
      "Get all the longitudes in the path");
  a1.def(
      "zeniths",
      [](const ArrayOfPropagationPathPoint& path) {
        return Vector{std::from_range,
                      path | stdv::transform([](const PropagationPathPoint& p) {
                        return p.zenith();
                      })};
      },
      "Get all the zeniths in the path");
  a1.def(
      "azimuths",
      [](const ArrayOfPropagationPathPoint& path) {
        return Vector{std::from_range,
                      path | stdv::transform([](const PropagationPathPoint& p) {
                        return p.azimuth();
                      })};
      },
      "Get all the azimuths in the path");
  generic_interface(a1);
  vector_interface(a1);

  auto a2 = py::bind_vector<ArrayOfArrayOfPropagationPathPoint,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfPropagationPathPoint");
  generic_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfArrayOfArrayOfPropagationPathPoint,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfArrayOfPropagationPathPoint");
  generic_interface(a3);
  vector_interface(a3);

  auto path  = m.def_submodule("path");
  path.doc() = "Submodule for path-related functionality";

  path.def("mirror", &path::mirror, "los"_a, "Mirror the line of sight");

  path.def(
      "lowest_point",
      [](const ArrayOfPropagationPathPoint& ray_path) -> PropagationPathPoint {
        if (ray_path.empty()) return {};

        return *stdr::min_element(
            ray_path,
            [](Vector3 pos1, Vector3 pos2) { return pos1[0] < pos2[0]; },
            &PropagationPathPoint::pos);
      },
      "ray_path"_a,
      "Get the lowest point of a ray path");

  path.def(
      "distance",
      [](const Vector3& pos1, const Vector3& pos2, const Vector2& ell) {
        return path::distance(pos1, pos2, ell);
      },
      "pos1"_a,
      "pos2"_a,
      "ell"_a,
      "Calculate the distance between two positions");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize path\n{}", e.what()));
}
}  // namespace Python
