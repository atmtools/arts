#include <hpy_arts.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <py_auto_options.h>
#include <tmatrix.h>

namespace Python {
namespace py = nanobind;
using namespace nanobind::literals;

void py_tmatrix(py::module_& m) {
  auto tm = m.def_submodule("tmatrix", R"(The T-matrix solver of Mishchenko, ported to C++.

Use available() to check whether this solver is available. All angles are in
degrees. Supply radius and wavelength in the same length unit (metres for SI
results). Results remain valid after subsequent calls.

See :doc:`user.tmatrix` for interface conventions, :doc:`concept.tmatrix` for
physical definitions, and :doc:`dev.tmatrix` for build and implementation details.
)");
  tm.def("available", &tmatrix::available, "Whether the T-matrix solver is built (ENABLE_TMATRIX).");
  py::class_<tmatrix::Geometry>(tm, "Geometry")
      .def(
          "__init__",
          [](tmatrix::Geometry* g,
             Numeric            theta_incident,
             Numeric            theta_scattered,
             Numeric            phi_incident,
             Numeric            phi_scattered,
             Numeric            alpha,
             Numeric            beta) {
            new (g) tmatrix::Geometry{.theta_incident  = theta_incident,
                                      .theta_scattered = theta_scattered,
                                      .phi_incident    = phi_incident,
                                      .phi_scattered   = phi_scattered,
                                      .alpha           = alpha,
                                      .beta            = beta};
          },
          "theta_incident"_a  = 0.0,
          "theta_scattered"_a = 0.0,
          "phi_incident"_a    = 0.0,
          "phi_scattered"_a   = 0.0,
          "alpha"_a           = 0.0,
          "beta"_a            = 0.0)
      .def_rw("theta_incident", &tmatrix::Geometry::theta_incident, "Incident zenith angle [deg].\n\n.. :class:`float`")
      .def_rw(
          "theta_scattered", &tmatrix::Geometry::theta_scattered, "Scattered zenith angle [deg].\n\n.. :class:`float`")
      .def_rw("phi_incident", &tmatrix::Geometry::phi_incident, "Incident azimuth angle [deg].\n\n.. :class:`float`")
      .def_rw("phi_scattered", &tmatrix::Geometry::phi_scattered, "Scattered azimuth angle [deg].\n\n.. :class:`float`")
      .def_rw("alpha", &tmatrix::Geometry::alpha, "First Euler angle of the particle axis [deg].\n\n.. :class:`float`")
      .def_rw("beta", &tmatrix::Geometry::beta, "Second Euler angle of the particle axis [deg].\n\n.. :class:`float`")
      .doc() = "A fixed-orientation geometry: incident and scattered directions and the particle's Euler angles";
  py::class_<tmatrix::FixedResult>(tm, "FixedResult")
      .def_ro("order", &tmatrix::FixedResult::order, "Converged multipole order.\n\n.. :class:`int`")
      .def_ro("scattering",
              &tmatrix::FixedResult::scattering,
              "Orientation-averaged scattering cross section, in length squared.\n\n.. :class:`float`")
      .def_ro("extinction",
              &tmatrix::FixedResult::extinction,
              "Orientation-averaged extinction cross section, in length squared.\n\n.. :class:`float`")
      .def_ro("amplitude",
              &tmatrix::FixedResult::amplitude,
              "The Jones matrix [[S_vv, S_vh], [S_hv, S_hh]] between ARTS's meridional bases of the two rays, in "
              "length units.\n\n.. :class:`~pyarts3.arts.Jones`")
      .def_ro("phase",
              &tmatrix::FixedResult::phase,
              "The Mueller matrix of amplitude in ARTS's Stokes convention, in length squared.\n\n.. "
              ":class:`~pyarts3.arts.Muelmat`")
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
          "The scattering matrix [F11, F12, F22, F33, F34, F44] at equally spaced angles from 0 to 180 degrees in the scattering-plane basis, in ARTS's Stokes convention. Dimensionless: F11 integrates to 4*pi over solid angle. Multiply by scattering/(4*pi) for differential cross sections.\n\n.. :class:`~pyarts3.arts.CompactPlanarMuelmatVector`")
      .def_ro(
          "legendre",
          &tmatrix::RandomResult::legendre,
          "The same scattering matrix as its exact Legendre series: the coefficients a_l, l = 0..l_max, on Y_l0 = sqrt((2 l + 1) / 4 pi) P_l(cos(Theta)) of each element, as TRO spectral data hold them.\n\n.. :class:`~pyarts3.arts.CompactPlanarMuelmatVector`")
      .doc() = "Random form results of T-matrix evaluation for one distribution";
  tm.def(
      "fixed",
      &tmatrix::fixed,
      "Compute a T-matrix and evaluate one geometry. equivalent_radius is an EquivalentRadius: the radius is that of the sphere of the same volume or of the same surface area. shape is a ParticleShape: a spheroid, a cylinder, or a Chebyshev particle r_0 (1 + aspect_ratio cos(degree theta)). aspect_ratio is horizontal/rotational axis for spheroids, diameter/length for cylinders, and the deformation, of magnitude below 1, for Chebyshev particles. degree is at least 1 for Chebyshev particles and 0 for the other shapes. coefficients, c_0, ..., c_N of r_0 (1 + c_0 + sum_n c_n cos(n theta)), are the shape of a generalized Chebyshev particle, e.g. a distorted water drop, and are empty for the other shapes; its aspect_ratio is not used. alpha and beta specify the particle Euler angles. surface_quadrature is the number of surface quadrature points per degree, which highly aspherical particles may need more of. The convergence fails, with an error, once the T-matrix is no longer finite.",
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
      "accuracy"_a           = 0.001,
      "equivalent_radius"_a  = EquivalentRadius::volume,
      "shape"_a              = ParticleShape::spheroid,
      "degree"_a             = 0,
      "coefficients"_a       = Vector{},
      "surface_quadrature"_a = 2,
      py::call_guard<py::gil_scoped_release>());
  tm.def(
      "fixed_batch",
      [](Numeric                               radius,
         Numeric                               wavelength,
         Numeric                               aspect_ratio,
         Numeric                               refractive_real,
         Numeric                               refractive_imag,
         const std::vector<tmatrix::Geometry>& geometries,
         Numeric                               accuracy,
         EquivalentRadius                      equivalent_radius,
         ParticleShape                         shape,
         Index                                 degree,
         const Vector&                         coefficients,
         Index                                 surface_quadrature) {
        return tmatrix::fixed_batch(radius,
                                    wavelength,
                                    aspect_ratio,
                                    refractive_real,
                                    refractive_imag,
                                    geometries,
                                    accuracy,
                                    equivalent_radius,
                                    shape,
                                    degree,
                                    coefficients,
                                    surface_quadrature);
      },
      R"(Compute one T-matrix and evaluate multiple geometries.

geometries is a list of :class:`Geometry`. Returns a list of FixedResult
objects in the same order. Parameters and units match fixed(). One T-matrix
serves all geometries; results own their data.
)",
      "radius"_a,
      "wavelength"_a,
      "aspect_ratio"_a,
      "refractive_real"_a,
      "refractive_imag"_a,
      "geometries"_a,
      "accuracy"_a           = 0.001,
      "equivalent_radius"_a  = EquivalentRadius::volume,
      "shape"_a              = ParticleShape::spheroid,
      "degree"_a             = 0,
      "coefficients"_a       = Vector{},
      "surface_quadrature"_a = 2,
      py::call_guard<py::gil_scoped_release>());
  py::class_<tmatrix::AzimuthallyRandomResult>(tm, "AzimuthallyRandomResult")
      .def_ro("order", &tmatrix::AzimuthallyRandomResult::order, "Converged multipole order.\n\n.. :class:`int`")
      .def_ro("phase",
              &tmatrix::AzimuthallyRandomResult::phase,
              "The phase matrix [za_inc, delta_aa, za_scat] between ARTS's meridional bases, in ARTS's Stokes "
              "convention, in length squared.\n\n.. :class:`~pyarts3.arts.MuelmatTensor3`")
      .def_ro("extinction",
              &tmatrix::AzimuthallyRandomResult::extinction,
              "K11, K12 and K34 of the extinction matrix [za_inc], in length squared.\n\n.. "
              ":class:`~pyarts3.arts.CompactPlanarPropmatVector`")
      .def_ro("absorption",
              &tmatrix::AzimuthallyRandomResult::absorption,
              "I and Q of the absorption vector [za_inc], in length squared.\n\n.. "
              ":class:`~pyarts3.arts.CompactPlanarStokvecVector`")
      .doc() = "Results of azimuthally random orientation on the grids of the call";
  tm.def("azimuthally_random",
         &tmatrix::azimuthally_random,
         R"(Compute particles of a tilt distribution with random azimuth (ARO).

The phase matrix, extinction matrix and absorption vector on ARTS's grids:
za_inc and za_scat the zenith angles of the propagation directions and
delta_aa the scattered minus the incident azimuth, in degrees, as ARTS's
scattering data hold them.  tilt_angles are the angles beta between the
particle's symmetry axis and the vertical, in degrees, with tilt_weights,
which need not be normalized; a single tilt fixes beta.  The averages over
the random azimuth and the integral of the absorption vector are exact for the
truncated T-matrix.  The other arguments match fixed().
)",
         "radius"_a,
         "wavelength"_a,
         "aspect_ratio"_a,
         "refractive_real"_a,
         "refractive_imag"_a,
         "za_inc"_a,
         "delta_aa"_a,
         "za_scat"_a,
         "tilt_angles"_a,
         "tilt_weights"_a,
         "accuracy"_a           = 0.001,
         "equivalent_radius"_a  = EquivalentRadius::volume,
         "shape"_a              = ParticleShape::spheroid,
         "degree"_a             = 0,
         "coefficients"_a       = Vector{},
         "surface_quadrature"_a = 2,
         py::call_guard<py::gil_scoped_release>());
  tm.def("random",
         &tmatrix::random,
         R"(Compute a randomly oriented particle size distribution.

The defaults approximate a monodisperse particle using a narrow radius interval.
distribution is a :class:`~pyarts3.arts.TMatrixSizeDistribution`, with radius,
b and gamma as its parameters (see :doc:`concept.tmatrix`). size_quadrature + 2
is the number of size quadrature points (one for -1); surface_quadrature is the
number of surface quadrature points per degree. For the power law, radius and b
are the effective radius and variance, which determine the integration bounds. Other
distributions use the given lower/upper radius ratios. angles is the number of
equally spaced scattering angles of phase, 0 for none; legendre is always
given. Radius and shape conventions match fixed(). The Wigner tables must be
initialized (WignerInit).
)",
         "radius"_a,
         "wavelength"_a,
         "aspect_ratio"_a,
         "refractive_real"_a,
         "refractive_imag"_a,
         "angles"_a             = 0,
         "accuracy"_a           = 0.001,
         "equivalent_radius"_a  = EquivalentRadius::volume,
         "shape"_a              = ParticleShape::spheroid,
         "degree"_a             = 0,
         "coefficients"_a       = Vector{},
         "distribution"_a       = TMatrixSizeDistribution::gamma,
         "b"_a                  = 0.1,
         "gamma"_a              = 1.,
         "size_quadrature"_a    = -1,
         "surface_quadrature"_a = 2,
         "lower_radius_ratio"_a = 0.9999999,
         "upper_radius_ratio"_a = 1.0000001,
         py::call_guard<py::gil_scoped_release>());
}
}  // namespace Python
