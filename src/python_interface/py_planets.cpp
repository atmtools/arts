#include <python_interface.h>

#include "core/functional/functional_gravity.h"
#include "core/util/planet_data.h"
#include "hpy_arts.h"
#include "hpy_numpy.h"

namespace Python {

void py_planets(py::module_& m) try {
  py::class_<EllipsoidGravity> eg(m, "EllipsoidGravity");
  eg.doc() = "Ellipsoid gravity operator - allows calculation of gravity on ellipsoidal bodies";
  generic_interface(eg);
  eg.def(py::init<>())
      .def_rw(
          "GM",
          &EllipsoidGravity::GM,
          "Standard gravitational parameter in m^3/s^2.\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("a",
              &EllipsoidGravity::a,
              "Semi-major axis in m.\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("e",
              &EllipsoidGravity::e,
              "Eccentricity.\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def(
          "__call__",
          [](const EllipsoidGravity& self,
             py::object h,
             py::object lat,
             py::object lon) {
            return vectorize(
                [&self](Numeric h, Numeric lat, Numeric lon) {
                  return self(h, lat, lon);
                },
                h,
                lat,
                lon);
          },
          "h"_a,
          "lat"_a,
          "lon"_a,
          "Calculate gravity at a given position.")
      .def_static("Mercury",
                  &EllipsoidGravity::Mercury,
                  "Mercury ellipsoid gravity operator")
      .def_static(
          "Venus", &EllipsoidGravity::Venus, "Venus ellipsoid gravity operator")
      .def_static(
          "Earth", &EllipsoidGravity::Earth, "Earth ellipsoid gravity operator")
      .def_static(
          "Moon", &EllipsoidGravity::Moon, "Moon ellipsoid gravity operator")
      .def_static(
          "Mars", &EllipsoidGravity::Mars, "Mars ellipsoid gravity operator")
      .def_static("Jupiter",
                  &EllipsoidGravity::Jupiter,
                  "Jupiter ellipsoid gravity operator")
      .def_static("Saturn",
                  &EllipsoidGravity::Saturn,
                  "Saturn ellipsoid gravity operator");
  auto planets = m.def_submodule("planets", "Planetary constants");

  auto earth              = planets.def_submodule("Earth", "Earth constants");
  earth.attr("ellipsoid") = Vector2{Body::Earth::a, Body::Earth::b};
  earth.attr("GM")        = Body::Earth::GM;
  earth.def(
      "gravity",
      &EllipsoidGravity::Earth,
      "Get Earth's gravity object.\n\n.. :class:`~pyarts3.arts.EllipsoidGravity`");

  auto jupiter = planets.def_submodule("Jupiter", "Jupiter constants");
  jupiter.attr("ellipsoid") = Vector2{Body::Jupiter::a, Body::Jupiter::b};
  jupiter.attr("GM")        = Body::Jupiter::GM;
  jupiter.def(
      "gravity",
      &EllipsoidGravity::Jupiter,
      "Get Jupiter's gravity object.\n\n.. :class:`~pyarts3.arts.EllipsoidGravity`");

  auto mars              = planets.def_submodule("Mars", "Mars constants");
  mars.attr("ellipsoid") = Vector2{Body::Mars::a, Body::Mars::b};
  mars.attr("GM")        = Body::Mars::GM;
  mars.def(
      "gravity",
      &EllipsoidGravity::Mars,
      "Get Mars' gravity object.\n\n.. :class:`~pyarts3.arts.EllipsoidGravity`");

  auto moon              = planets.def_submodule("Moon", "Moon constants");
  moon.attr("ellipsoid") = Vector2{Body::Moon::a, Body::Moon::b};
  moon.attr("GM")        = Body::Moon::GM;
  moon.def(
      "gravity",
      &EllipsoidGravity::Moon,
      "Get Moon's gravity object.\n\n.. :class:`~pyarts3.arts.EllipsoidGravity`");

  auto mercury = planets.def_submodule("Mercury", "Mercury constants");
  mercury.attr("ellipsoid") = Vector2{Body::Mercury::a, Body::Mercury::b};
  mercury.attr("GM")        = Body::Mercury::GM;
  mercury.def(
      "gravity",
      &EllipsoidGravity::Mercury,
      "Get Mercury's gravity object.\n\n.. :class:`~pyarts3.arts.EllipsoidGravity`");

  auto venus              = planets.def_submodule("Venus", "Venus constants");
  venus.attr("ellipsoid") = Vector2{Body::Venus::a, Body::Venus::b};
  venus.attr("GM")        = Body::Venus::GM;
  venus.def(
      "gravity",
      &EllipsoidGravity::Venus,
      "Get Venus' gravity object.\n\n.. :class:`~pyarts3.arts.EllipsoidGravity`");

  auto saturn = planets.def_submodule("Saturn", "Saturn constants");
  saturn.attr("ellipsoid") = Vector2{Body::Saturn::a, Body::Saturn::b};
  saturn.attr("GM")        = Body::Saturn::GM;
  saturn.def(
      "gravity",
      &EllipsoidGravity::Saturn,
      "Get Saturn's gravity object.\n\n.. :class:`~pyarts3.arts.EllipsoidGravity`");

} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize planets\n{}", e.what()));
}

}  // namespace Python