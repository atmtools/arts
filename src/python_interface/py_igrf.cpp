#include <igrf13.h>
#include <python_interface.h>

namespace Python {
void py_igrf(py::module_ &m) try {
  m.def(
      "igrf",
      [](const Tensor3 &z_field,
         const Vector &lat_grid,
         const Vector &lon_grid,
         const Time &t,
         const Vector &ell) {
        auto [u, v, w] = IGRF::compute(z_field, lat_grid, lon_grid, t, ell);
        return std::array{std::move(u), std::move(v), std::move(w)};
      },
      py::arg("z_field"),
      py::arg("lat_grid"),
      py::arg("lon_grid"),
      py::arg("t") = Time{},
      py::arg("ell") = Vector{6378137, 0.081819190842621},
      R"--(Compute the magnetic field according to IGRF

The coefficients are builtin and available for the years 2000-2020.
No secular variation is taken into account beyond an interpolation
between the available years.

Parameters
----------
z_field : Tensor3
    The altitude of the field points.
lat_grid : Vector
    The latitude grid.
lon_grid : Vector
    The longitude grid.
t : Time, optional
    The time of the magnetic field, default is the current time.
ell : Vector, optional
    The ellipsoid model. The default is [6378137, 0.081819190842621].

Returns
-------
u : Tensor3
    The u component of the magnetic field.
v : Tensor3
    The v component of the magnetic field.
w : Tensor3
    The w component of the magnetic field.
)--");
} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize IGRF\n", e.what()));
}
}  // namespace Python
