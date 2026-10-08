#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <rt3.h>
#include <rt3_arts.h>

#include "hpy_arts.h"

NB_MAKE_OPAQUE(std::vector<rt3::scattering_set>);

namespace Python {
void py_rt3(py::module_& m) try {
  auto rt = m.def_submodule("rt3", R"(Evans' RT3 (polradtran) polarized doubling-adding solver.

A plane-parallel solver for randomly oriented particles with a solar beam and
thermal sources, for every Fourier azimuth mode and the Stokes components
[I], [I, Q], [I, Q, U] or [I, Q, U, V], kept as an external reference for
other solvers.  It has no workspace layer.  Use available() to check whether
the optional Fortran backend is built (ENABLE_RT3=ON); otherwise solve()
raises.

Conventions: a right-handed frame with z up; the direct beam propagates
downward toward azimuth 0, and phi is the azimuth of the propagation direction
of a ray (VDISORT's phi for phi0 = 0).  Stokes basis [I, Q, U, V] with
:math:`h = k x z / |k x z|`, :math:`v = h x k`, :math:`Q = I_v - I_h` and :math:`U = 2 Re(E_v E_h*)`, the
same basis in both hemispheres.  Streams are :math:`\mu = |cos(zenith)|`, ascending
quadrature nodes followed by the zero-weight ``extra_mu``.  Layers and levels
are top-down.  Radiances are in W m-2 Hz-1 sr-1, fluxes in W m-2 Hz-1.
RT3Result.up and RT3Result.down are Fourier coefficients: the radiance is
sum_m c_m cos(m phi) for I, Q and sum_m c_m sin(m phi) for U, V
(azimuth_radiance() sums them).  Lengths and extinctions must use reciprocal
units.

See :doc:`dev.rt3` for all conventions, limitations and the benchmarks.
)");

  rt.def("available", &rt3::available, "Whether the Fortran backend is enabled (ENABLE_RT3=ON).");

  py::enum_<rt3::quadrature_type>(rt, "QuadratureType", "RT3 quadrature rules on one hemisphere")
      .value("gauss",
             rt3::quadrature_type::gauss,
             "RT3 'G' ('E' with extra angles): positive half of a 2*nmu-point Gauss-Legendre rule on [-1, 1]")
      .value("double_gauss", rt3::quadrature_type::double_gauss, "RT3 'D': nmu-point Gauss-Legendre rule on [0, 1]")
      .value("lobatto",
             rt3::quadrature_type::lobatto,
             "RT3 'L': positive half of a 2*nmu-point Lobatto rule on [-1, 1]; includes mu = 1");

  py::class_<rt3::quadrature>(rt, "Quadrature")
      .def_ro("mu", &rt3::quadrature::mu, "Ascending nodes in (0, 1]\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_ro("weights",
              &rt3::quadrature::weights,
              "Weights for the integral over mu in [0, 1], summing to 1 (no 2 pi azimuth factor)\n\n.. "
              ":class:`~pyarts3.arts.Vector`")
      .doc() = "The streams of one hemisphere";

  rt.def("get_quadrature",
         &rt3::get_quadrature,
         "nmu"_a,
         "type"_a = rt3::quadrature_type::gauss,
         "RT3's streams, nmu >= 1 nodes per hemisphere: the positive half of ARTS's 2 nmu-point rule, as "
         "RT3 uses them.",
         py::call_guard<py::gil_scoped_release>());

  rt.def("max_legendre_degree",
         &rt3::max_legendre_degree,
         "nmu"_a,
         "type"_a = rt3::quadrature_type::gauss,
         R"(RT3's NLEGLIM: the highest Legendre degree it keeps for nmu quadrature nodes.

gauss 4 nmu - 3, double_gauss 2 nmu - 3, lobatto 4 nmu - 5, at least 1.  RT3
truncates longer series (after delta-M scaling); solve() raises instead when
that would drop a non-zero coefficient.)");

  py::class_<rt3::scattering_set>(rt, "ScatteringSet")
      .def(
          "__init__",
          [](rt3::scattering_set* s, Numeric extinction, Numeric scattering, const Matrix& legendre) {
            new (s) rt3::scattering_set{.extinction = extinction, .scattering = scattering, .legendre = legendre};
          },
          "extinction"_a,
          "scattering"_a,
          "legendre"_a)
      .def_rw("extinction",
              &rt3::scattering_set::extinction,
              "Particle extinction coefficient per unit length\n\n.. :class:`float`")
      .def_rw("scattering",
              &rt3::scattering_set::scattering,
              "Particle scattering coefficient per unit length\n\n.. :class:`float`")
      .def_rw("legendre",
              &rt3::scattering_set::legendre,
              "[nleg + 1, 6] Legendre coefficients in cos(Theta) of (F11, F12, F33, F34, F22, F44), legendre[0, 0] "
              "= 1\n\n.. :class:`~pyarts3.arts.Matrix`")
      .doc() = R"(Single-scattering properties of a homogeneous particle population.

The scattering-plane phase matrix [[F11, F12, 0, 0], [F12, F22, 0, 0],
[0, 0, F33, F34], [0, 0, -F34, F44]] with Q = I_par - I_perp, each element a
plain Legendre series F_c(Theta) = sum_l legendre[l, c] P_l(cos(Theta)) in
the column order of RT3's scattering files, c = 0: F11, 1: F12, 2: F33,
3: F34, 4: F22, 5: F44.  The coefficients include the factor 2 l + 1
(Henyey-Greenstein is (2 l + 1) g^l) and the phase function is normalised to
1 over 4 pi, legendre[0, 0] = 1.  Rayleigh scattering is
[[1, -1/2, 0, 0, 1, 0], [0, 0, 3/2, 0, 0, 3/2], [1/2, 1/2, 0, 0, 1/2, 0]].)";

  auto aoss =
      py::bind_vector<std::vector<rt3::scattering_set>, py::rv_policy::reference_internal>(rt, "ArrayOfScatteringSet");
  aoss.doc() = "A list of ScatteringSet";

  py::class_<rt3::lambertian_surface>(rt, "LambertianSurface")
      .def(
          "__init__",
          [](rt3::lambertian_surface* s, Numeric albedo) { new (s) rt3::lambertian_surface{.albedo = albedo}; },
          "albedo"_a = 0.0)
      .def_rw("albedo", &rt3::lambertian_surface::albedo, "Albedo A\n\n.. :class:`float`")
      .doc() =
      "RT3 'L': reflection 2 A mu_j w_j of the m = 0 mode, I to I only; emission (1 - A) B with thermal set; "
      "direct-beam reflection A F / pi.  Energy is conserved on the streams only for double_gauss quadrature.";

  py::class_<rt3::fresnel_surface>(rt, "FresnelSurface")
      .def(
          "__init__",
          [](rt3::fresnel_surface* s, Complex n) { new (s) rt3::fresnel_surface{.refractive_index = n}; },
          "refractive_index"_a)
      .def_rw("refractive_index",
              &rt3::fresnel_surface::refractive_index,
              "Complex refractive index of the surface (medium above has index 1)\n\n.. :class:`complex`")
      .doc() =
      "RT3 'F': specular Fresnel reflection with :math:`R1 = (|r_v|^2 + |r_h|^2) / 2, R2 = (|r_v|^2 - |r_h|^2) / 2, "
      "R3 = Re(r_v r_h*), R4 = Im(r_v r_h*)`; emission :math:`[(1 - R1) B, -R2 B, 0, 0]`, also when thermal is false.  "
      "Not allowed with a direct beam.";

  const rt3::problem d{};
  py::class_<rt3::problem>(rt, "Problem")
      .def(
          "__init__",
          [](rt3::problem*                           p,
             Index                                   nstokes,
             Index                                   nmu,
             rt3::quadrature_type                    quad,
             const Vector&                           extra_mu,
             Index                                   aziorder,
             Numeric                                 max_delta_tau,
             bool                                    delta_m,
             Numeric                                 direct_flux,
             Numeric                                 direct_mu,
             bool                                    thermal,
             Numeric                                 frequency,
             const Vector&                           height,
             const Vector&                           temperature,
             const Vector&                           gas_extinction,
             const std::vector<rt3::scattering_set>& scattering_sets,
             const ArrayOfIndex&                     layer_scattering_index,
             Numeric                                 sky_temperature,
             Numeric                                 surface_temperature,
             const rt3::surface&                     ground) {
            new (p) rt3::problem{.nstokes                = nstokes,
                                 .nmu                    = nmu,
                                 .quad                   = quad,
                                 .extra_mu               = extra_mu,
                                 .aziorder               = aziorder,
                                 .max_delta_tau          = max_delta_tau,
                                 .delta_m                = delta_m,
                                 .direct_flux            = direct_flux,
                                 .direct_mu              = direct_mu,
                                 .thermal                = thermal,
                                 .frequency              = frequency,
                                 .height                 = height,
                                 .temperature            = temperature,
                                 .gas_extinction         = gas_extinction,
                                 .scattering_sets        = scattering_sets,
                                 .layer_scattering_index = layer_scattering_index,
                                 .sky_temperature        = sky_temperature,
                                 .surface_temperature    = surface_temperature,
                                 .ground                 = ground};
          },
          "nstokes"_a                = d.nstokes,
          "nmu"_a                    = d.nmu,
          "quad"_a                   = d.quad,
          "extra_mu"_a               = d.extra_mu,
          "aziorder"_a               = d.aziorder,
          "max_delta_tau"_a          = d.max_delta_tau,
          "delta_m"_a                = d.delta_m,
          "direct_flux"_a            = d.direct_flux,
          "direct_mu"_a              = d.direct_mu,
          "thermal"_a                = d.thermal,
          "frequency"_a              = d.frequency,
          "height"_a                 = d.height,
          "temperature"_a            = d.temperature,
          "gas_extinction"_a         = d.gas_extinction,
          "scattering_sets"_a        = d.scattering_sets,
          "layer_scattering_index"_a = d.layer_scattering_index,
          "sky_temperature"_a        = d.sky_temperature,
          "surface_temperature"_a    = d.surface_temperature,
          "ground"_a                 = d.ground)
      .def_rw("nstokes", &rt3::problem::nstokes, "1 to 4 for [I] ... [I, Q, U, V]\n\n.. :class:`int`")
      .def_rw("nmu", &rt3::problem::nmu, "Quadrature nodes per hemisphere\n\n.. :class:`int`")
      .def_rw("quad", &rt3::problem::quad, "Quadrature rule\n\n.. :class:`~pyarts3.arts.rt3.QuadratureType`")
      .def_rw("extra_mu",
              &rt3::problem::extra_mu,
              "Zero-weight output angles appended after the quadrature streams, each in (0, 1]; gauss only\n\n.. "
              ":class:`~pyarts3.arts.Vector`")
      .def_rw("aziorder", &rt3::problem::aziorder, "Highest Fourier azimuth mode, >= 0\n\n.. :class:`int`")
      .def_rw("max_delta_tau",
              &rt3::problem::max_delta_tau,
              "Maximum vertical optical thickness of the initial doubling sublayer, > 0\n\n.. :class:`float`")
      .def_rw("delta_m",
              &rt3::problem::delta_m,
              "Delta-M scaling of every scattering set, M = 2 nmu_total\n\n.. :class:`bool`")
      .def_rw("direct_flux",
              &rt3::problem::direct_flux,
              "Direct-beam flux on the horizontal at the top [W m-2 Hz-1], >= 0; 0 switches the beam off\n\n.. "
              ":class:`float`")
      .def_rw("direct_mu",
              &rt3::problem::direct_mu,
              "Cosine of the zenith angle of the direct beam, in (0, 1]\n\n.. :class:`float`")
      .def_rw("thermal",
              &rt3::problem::thermal,
              "Thermal emission of the layers and of a Lambertian surface (the sky and a Fresnel surface always "
              "emit)\n\n.. :class:`bool`")
      .def_rw("frequency", &rt3::problem::frequency, "Frequency [Hz]\n\n.. :class:`float`")
      .def_rw("height",
              &rt3::problem::height,
              "[nlay + 1] layer interfaces, top-down; only ``|differences|`` are used, in the reciprocal of the "
              "extinction unit\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("temperature",
              &rt3::problem::temperature,
              "[nlay + 1] interface temperatures [K], top-down; > 0 with thermal, unused otherwise\n\n.. "
              ":class:`~pyarts3.arts.Vector`")
      .def_rw("gas_extinction",
              &rt3::problem::gas_extinction,
              "[nlay] scalar, unpolarized gas extinction per unit length, >= 0\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("scattering_sets",
              &rt3::problem::scattering_sets,
              "Scattering sets, at most 200\n\n.. :class:`~pyarts3.arts.rt3.ArrayOfScatteringSet`")
      .def_rw("layer_scattering_index",
              &rt3::problem::layer_scattering_index,
              "[nlay] index into scattering_sets per layer, or < 0 for a gas-only layer\n\n.. "
              ":class:`~pyarts3.arts.ArrayOfIndex`")
      .def_rw("sky_temperature",
              &rt3::problem::sky_temperature,
              "Temperature [K] of the isotropic, unpolarized blackbody incident at the top\n\n.. :class:`float`")
      .def_rw("surface_temperature", &rt3::problem::surface_temperature, "Surface temperature [K]\n\n.. :class:`float`")
      .def_prop_rw(
          "ground",
          [](const rt3::problem& p) { return p.ground; },
          [](rt3::problem& p, const rt3::surface& g) { p.ground = g; },
          "The surface (a copy is returned; assign to change it)\n\n.. :class:`~pyarts3.arts.rt3.LambertianSurface` "
          "| :class:`~pyarts3.arts.rt3.FresnelSurface`")
      .doc() = "An RT3 problem; see :doc:`dev.rt3` for the conventions";

  py::class_<rt3::result>(rt, "RT3Result")
      .def_ro("mu", &rt3::result::mu, "[nmu_total] streams\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_ro("weights",
              &rt3::result::weights,
              "[nmu_total] quadrature weights, 0 for the extra angles\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_ro("up",
              &rt3::result::up,
              "[nlay + 1 level (0 = top), aziorder + 1 mode, nmu_total, nstokes] Fourier coefficients of the "
              "upward-propagating radiance in W m-2 Hz-1 sr-1\n\n.. :class:`~pyarts3.arts.Tensor4`")
      .def_ro("down",
              &rt3::result::down,
              "[nlay + 1 level (0 = top), aziorder + 1 mode, nmu_total, nstokes] Fourier coefficients of the "
              "downward-propagating diffuse radiance in W m-2 Hz-1 sr-1\n\n.. :class:`~pyarts3.arts.Tensor4`")
      .def_ro("up_flux",
              &rt3::result::up_flux,
              "[nlay + 1, nstokes] upward flux 2 pi sum_i w_i mu_i c_0 in W m-2 Hz-1\n\n.. "
              ":class:`~pyarts3.arts.Matrix`")
      .def_ro("down_flux",
              &rt3::result::down_flux,
              "[nlay + 1, nstokes] downward flux in W m-2 Hz-1, including the direct beam\n\n.. "
              ":class:`~pyarts3.arts.Matrix`")
      .doc() = "The RT3 solution at every level";

  rt.def("solve",
         &rt3::solve,
         "problem"_a,
         R"(Run RT3.

Validates the shapes and every precondition on which the Fortran code would
stop the process or overrun a buffer, rejects a Legendre series that RT3
would silently truncate, then calls RADTRAN.  Calls are serialised by one
lock (separate from RT4's).)",
         py::call_guard<py::gil_scoped_release>());

  rt.def("azimuth_radiance",
         &rt3::azimuth_radiance,
         "coefficients"_a,
         "phi"_a,
         R"(The radiance at azimuths phi [rad] from Fourier coefficients.

coefficients is [nlevel, nmode, nmu, nstokes], e.g. RT3Result.up or RT3Result.down.
Returns [nlevel, len(phi), nmu, nstokes] with sum_m c_m cos(m phi) for I, Q
and sum_m c_m sin(m phi) for U, V, as rt3.f's OUTPUT_FILE does.)");

  rt.def("scattering_optics",
         &rt3::scattering_optics,
         "scattering_species"_a,
         "atm_point"_a,
         "frequency"_a,
         "degree"_a,
         "normalisation_tolerance"_a = 1e-3,
         R"(The RT3 scattering set of ARTS scattering species at one atmospheric point.

The species give their totally randomly oriented (TRO) Legendre series up to
``degree`` themselves (``get_bulk_scattering_properties_tro_spectral``);
gridded particle data must be converted to a Legendre series first
(``ParticleHabit.to_tro_spectral_with_report``).  ARTS's elements
[F11, F12, F22, F33, F34, F44] are reordered to RT3's columns
(F11, F12, F33, F34, F22, F44) without sign changes, and the series is
normalised so that ``legendre[0, 0] == 1``.  ``extinction`` is K11 and
``scattering`` is K11 - a1; the phase-function integral must match the latter
to ``normalisation_tolerance`` times the extinction.

For the same physical sphere, ARTS's Mie F34 is -1 times that of Evans and
Stephens (1991), so V computed from ARTS data is -1 times V in their
convention.  See :doc:`dev.rt3`.
)");

  const rt3::path_settings ds{};
  py::class_<rt3::path_settings>(rt, "PathSettings")
      .def(
          "__init__",
          [](rt3::path_settings*  s,
             Index                nstokes,
             Index                nmu,
             rt3::quadrature_type quad,
             const Vector&        extra_mu,
             Index                aziorder,
             Numeric              max_delta_tau,
             bool                 delta_m,
             Index                legendre_degree,
             Numeric              normalisation_tolerance) {
            new (s) rt3::path_settings{.nstokes                 = nstokes,
                                       .nmu                     = nmu,
                                       .quad                    = quad,
                                       .extra_mu                = extra_mu,
                                       .aziorder                = aziorder,
                                       .max_delta_tau           = max_delta_tau,
                                       .delta_m                 = delta_m,
                                       .legendre_degree         = legendre_degree,
                                       .normalisation_tolerance = normalisation_tolerance};
          },
          "nstokes"_a                 = ds.nstokes,
          "nmu"_a                     = ds.nmu,
          "quad"_a                    = ds.quad,
          "extra_mu"_a                = ds.extra_mu,
          "aziorder"_a                = ds.aziorder,
          "max_delta_tau"_a           = ds.max_delta_tau,
          "delta_m"_a                 = ds.delta_m,
          "legendre_degree"_a         = ds.legendre_degree,
          "normalisation_tolerance"_a = ds.normalisation_tolerance)
      .def_rw("nstokes", &rt3::path_settings::nstokes, "1 to 4\n\n.. :class:`int`")
      .def_rw("nmu", &rt3::path_settings::nmu, "Quadrature nodes per hemisphere\n\n.. :class:`int`")
      .def_rw("quad", &rt3::path_settings::quad, "Quadrature rule\n\n.. :class:`~pyarts3.arts.rt3.QuadratureType`")
      .def_rw("extra_mu",
              &rt3::path_settings::extra_mu,
              "Zero-weight output angles (gauss only)\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("aziorder", &rt3::path_settings::aziorder, "Highest Fourier azimuth mode\n\n.. :class:`int`")
      .def_rw("max_delta_tau",
              &rt3::path_settings::max_delta_tau,
              "Maximum vertical optical thickness of the initial doubling sublayer\n\n.. :class:`float`")
      .def_rw("delta_m", &rt3::path_settings::delta_m, "RT3's delta-M scaling\n\n.. :class:`bool`")
      .def_rw("legendre_degree",
              &rt3::path_settings::legendre_degree,
              "Degree of the Legendre series; negative selects RT3's maximum (at least 2 nmu_total with "
              "delta_m)\n\n.. :class:`int`")
      .def_rw("normalisation_tolerance",
              &rt3::path_settings::normalisation_tolerance,
              "Allowed mismatch of phase-function integral and scattering coefficient, relative to the "
              "extinction; infinity checks nothing\n\n.. :class:`float`")
      .doc() = "Solver settings of problem_from_path";

  rt.def("problem_from_path",
         &rt3::problem_from_path,
         "ray_path"_a,
         "atm_path"_a,
         "spectral_propmat_path"_a,
         "freq_grid"_a,
         "freq_index"_a,
         "scattering_species"_a,
         "settings"_a,
         "ground"_a,
         "surface_temperature"_a,
         "sky_temperature"_a,
         R"(An RT3 problem from an ARTS propagation path, with the DISORT path conventions.

The path conventions are those of :func:`pyarts3.arts.rt4.problem_from_path`.
A layer's scattering set has the mean extinction and scattering of its two
levels' :func:`scattering_optics` and their scattering-weighted mean series;
a layer with zero mean extinction is gas-only.  The problem has thermal
emission and no beam: set ``direct_flux``, ``direct_mu`` and ``thermal`` on
the result for other sources.  See :doc:`dev.rt3`.
)");
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize rt3\n{}", e.what()));
}
}  // namespace Python
