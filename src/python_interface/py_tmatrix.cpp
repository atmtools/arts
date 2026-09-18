#include <hpy_arts.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <tmatrix.h>

namespace Python {
namespace py = nanobind;
using namespace nanobind::literals;

void py_tmatrix(py::module_& m) {
  auto tm = m.def_submodule("tmatrix", R"(Original Mishchenko T-matrix solver.

Use available() to check whether this solver is available. All angles are in
degrees. Supply radius and wavelength in the same length unit (metres for SI
results). Results remain valid after subsequent calls.

See :doc:`user.tmatrix` for interface conventions, :doc:`concept.tmatrix` for
physical definitions, and :doc:`dev.tmatrix` for build and implementation details.
)");
  tm.def("available", &tmatrix::available, "Whether the Fortran backend is enabled.");
  tm.def("extended_precision",
         &tmatrix::extended_precision,
         "Whether the selected backend uses extended-precision internals.");
  py::class_<tmatrix::FixedResult>(tm, "FixedResult")
      .def_ro("order", &tmatrix::FixedResult::order, "Converged multipole order.\n\n.. :class:`int`")
      .def_ro("scattering",
              &tmatrix::FixedResult::scattering,
              "Orientation-averaged scattering cross section, in length squared.\n\n.. :class:`float`")
      .def_ro("extinction",
              &tmatrix::FixedResult::extinction,
              "Orientation-averaged extinction cross section, in length squared.\n\n.. :class:`float`")
      .def_prop_ro(
          "amplitude",
          [](const tmatrix::FixedResult& x) { return ComplexMatrix{x.amplitude}; },
          "Complex 2 by 2 amplitude matrix, in length units.\n\n.. :class:`~pyarts3.arts.ComplexMatrix`")
      .def_ro(
          "phase",
          &tmatrix::FixedResult::phase,
          "4 by 4 phase matrix in the original Stokes convention, in length squared.\n\n.. :class:`~pyarts3.arts.Muelmat`")
      .doc() = "Fixed form results of T-matrix evaluation for one geometry";
  py::class_<tmatrix::RandomResult>(tm, "RandomResult")
      .def_ro("effective_radius",
              &tmatrix::RandomResult::effective_radius,
              "Effective radius in the input length unit and radius convention.\n\n.. :class:`float`")
      .def_ro("effective_variance",
              &tmatrix::RandomResult::effective_variance,
              "Dimensionless effective variance.\n\n.. :class:`float`")
      .def_ro("extinction",
              &tmatrix::RandomResult::extinction,
              "Distribution-averaged extinction cross section, in length squared.\n\n.. :class:`float`")
      .def_ro("scattering",
              &tmatrix::RandomResult::scattering,
              "Distribution-averaged scattering cross section, in length squared.\n\n.. :class:`float`")
      .def_ro("albedo", &tmatrix::RandomResult::albedo, "Single-scattering albedo.\n\n.. :class:`float`")
      .def_ro("asymmetry", &tmatrix::RandomResult::asymmetry, "Mean cosine of scattering angle.\n\n.. :class:`float`")
      .def_ro(
          "phase",
          &tmatrix::RandomResult::phase,
          "MuelmatVector at equally spaced angles from 0 to 180 degrees in the scattering-plane basis. Dimensionless: F11 integrates to 4*pi over solid angle. Multiply by scattering/(4*pi) for differential cross sections.\n\n.. :class:`~pyarts3.arts.MuelmatVector`")
      .doc() = "Random form results of T-matrix evaluation for one distribution";
  tm.def(
      "fixed",
      &tmatrix::fixed,
      "Compute a T-matrix and evaluate one geometry. radius_ratio=1 selects equal-volume radius; any other positive value selects equal-surface-area radius for spheroids/cylinders. shape=-1 selects a spheroid, -2 a cylinder. aspect_ratio is horizontal/rotational axis for spheroids, diameter/length for cylinders. alpha and beta specify the particle Euler angles.",
      "radius"_a,
      "wavelength"_a,
      "aspect_ratio"_a,
      "refractive_real"_a,
      "refractive_imag"_a,
      "theta_incident"_a,
      "theta_scattered"_a,
      "phi_incident"_a,
      "phi_scattered"_a,
      "alpha"_a,
      "beta"_a,
      "accuracy"_a     = 0.001,
      "radius_ratio"_a = 1.,
      "shape"_a        = -1,
      py::call_guard<py::gil_scoped_release>());
  tm.def(
      "fixed_batch",
      [](Numeric       radius,
         Numeric       wavelength,
         Numeric       aspect_ratio,
         Numeric       refractive_real,
         Numeric       refractive_imag,
         const Matrix& geometries,
         Numeric       accuracy,
         Numeric       radius_ratio,
         int           shape) {
        return tmatrix::fixed_batch(radius,
                                    wavelength,
                                    aspect_ratio,
                                    refractive_real,
                                    refractive_imag,
                                    geometries,
                                    accuracy,
                                    radius_ratio,
                                    shape);
      },
      R"(Compute one T-matrix and evaluate multiple geometries.

Each row of geometries contains theta_incident, theta_scattered,
phi_incident, phi_scattered, alpha, beta, in degrees. Returns a list of
FixedResult objects in row order. Parameters and units match fixed().
The solver lock covers the entire batch; results own their data.
)",
      "radius"_a,
      "wavelength"_a,
      "aspect_ratio"_a,
      "refractive_real"_a,
      "refractive_imag"_a,
      "geometries"_a,
      "accuracy"_a     = 0.001,
      "radius_ratio"_a = 1.,
      "shape"_a        = -1,
      py::call_guard<py::gil_scoped_release>());
  tm.def("random",
         &tmatrix::random,
         R"(Compute a randomly oriented particle size distribution.

The defaults approximate a monodisperse particle using a narrow radius interval.
For the original TMD distribution definitions: 1 modified gamma, 2 lognormal,
3 Hansen-Travis power law, 4 gamma, 5 modified power law. b and gamma have the
meaning documented in 3rdparty/tmatrix/tmd.lp.f. size_quadrature is NKMAX;
surface_quadrature is NDGS. A call evaluates one distribution (NPNAX=1).
For distribution=3, radius and b are the effective radius and variance, and
Fortran computes the integration bounds. Other distributions use the given
lower/upper radius ratios. Radius and shape conventions match fixed().
)",
         "radius"_a,
         "wavelength"_a,
         "aspect_ratio"_a,
         "refractive_real"_a,
         "refractive_imag"_a,
         "angles"_a             = 19,
         "accuracy"_a           = 0.001,
         "radius_ratio"_a       = 1.,
         "shape"_a              = -1,
         "distribution"_a       = 4,
         "b"_a                  = 0.1,
         "gamma"_a              = 1.,
         "size_quadrature"_a    = -1,
         "surface_quadrature"_a = 2,
         "lower_radius_ratio"_a = 0.9999999,
         "upper_radius_ratio"_a = 1.0000001,
         py::call_guard<py::gil_scoped_release>());
}
}  // namespace Python
