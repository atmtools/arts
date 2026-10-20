#include <igrf.h>
#include <python_interface.h>

namespace Python {
void py_igrf(py::module_& m) try {
  m.def("igrf",
        &IGRF::igrf,
        "pos"_a,
        "ell"_a = Vector2{6378137.0, 6356752.314245},
        "t"_a   = Time{},
        R"--(Compute the magnetic field according to IGRF

The IGRF coefficients are compiled from ``mag-field/igrf.txt`` in arts-cat-data.
They are interpolated linearly in time between the epochs of the model, and the
secular variation of the model is used after its last epoch.  The first and the
last valid decimal year, where 2001.0 is the start of 2001, are
``pyarts3.arts.igrf_valid_range``.  Other times are an error.

Parameters
----------
pos : Vector3
    The position in [alt, lat, lon] (geodetic)
ell : Vector2, optional
    The ellipsoid (a, b), default is WGS84 [6378137.0, 6356752.314245]
t : Time, optional
    A time stamp, default is the current time

Returns
-------
mag : Vector3
    The magnetic field [T] in geodetic [east, north, up]
)--");

  m.attr("igrf_valid_range") = py::make_tuple(IGRF::Data::valid_range.first, IGRF::Data::valid_range.second);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize IGRF\n{}", e.what()));
}
}  // namespace Python
