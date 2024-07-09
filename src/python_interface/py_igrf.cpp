#include <igrf13.h>
#include <python_interface.h>

namespace Python {
void py_igrf(py::module_& m) try {
  m.def("igrf",
        &IGRF::igrf,
        "pos"_a,
        "ell"_a = Vector2{6378137.0, 6356752.314245},
        "t"_a   = Time{},
        R"--(Compute the magnetic field according to IGRF

The coefficients are builtin and available for the years 2000-2020.
No secular variation is taken into account beyond an interpolation
between the available years.

Parameters
----------
pos : Vector3
    The position in [alt, lat, lon] (geodetic)
ell : Vector2, optional
    The ellipsoid (a, b), default is WGS84 [6378137.0, 6356752.314245]
t : Time, optional
    A time stamp, default is the current time
)--");
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize IGRF\n", e.what()));
}
}  // namespace Python
