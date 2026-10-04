#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <rt4.h>
#include <rt4_arts.h>

#include "hpy_arts.h"

NB_MAKE_OPAQUE(std::vector<rt4::layer_optics>);

namespace Python {
void py_rt4(py::module_& m) try {
  auto rt = m.def_submodule("rt4", R"(Evans' RT4 (polradtran) polarized doubling-adding solver.

A thermal-only, plane-parallel solver for azimuthally symmetric media and the
Stokes components [I] or [I, Q], kept as an external reference for other
solvers.  It has no workspace layer.  Use available() to check whether the
optional Fortran backend is built (ENABLE_RT4=ON); otherwise get_quadrature()
and solve() raise.

Conventions: Stokes basis [I, Q] with Q = I_v - I_h in the meridional plane,
the same basis in both hemispheres.  Hemisphere index ``down`` (0) is
radiation propagating downward, ``up`` (1) propagating upward.  Streams are
mu = |cos(zenith)|, ascending quadrature nodes followed by the zero-weight
``extra_mu``.  Layers and levels are top-down.  Radiances are in
W m-2 Hz-1 sr-1.  Lengths and extinctions must use reciprocal units.

See :doc:`dev.rt4` for all conventions, limitations and the mapping to VDISORT.
)");

  rt.def("available", &rt4::available, "Whether the Fortran backend is enabled (ENABLE_RT4=ON).");
  rt.attr("down") = rt4::down;
  rt.attr("up")   = rt4::up;

  py::enum_<rt4::quadrature_type>(rt, "QuadratureType", "RT4 quadrature rules on one hemisphere")
      .value("double_gauss", rt4::quadrature_type::double_gauss, "RT4 'D': nmu-point Gauss-Legendre rule on [0, 1]")
      .value("gauss",
             rt4::quadrature_type::gauss,
             "RT4 'G': positive half of a 2*nmu-point Gauss-Legendre rule on [-1, 1]")
      .value("lobatto",
             rt4::quadrature_type::lobatto,
             "RT4 'L': positive half of a 2*nmu-point Lobatto rule on [-1, 1]; includes mu = 1");

  py::class_<rt4::quadrature>(rt, "Quadrature")
      .def_ro("mu", &rt4::quadrature::mu, "Ascending nodes in (0, 1]\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_ro("weights",
              &rt4::quadrature::weights,
              "Weights for the integral over mu in [0, 1], summing to 1 (no 2 pi azimuth factor)\n\n.. "
              ":class:`~pyarts3.arts.Vector`")
      .doc() = "The streams of one hemisphere";

  rt.def("get_quadrature",
         &rt4::get_quadrature,
         "nmu"_a,
         "type"_a = rt4::quadrature_type::double_gauss,
         "RT4's own quadrature with nmu >= 1 nodes per hemisphere.",
         py::call_guard<py::gil_scoped_release>());

  py::class_<rt4::layer_optics>(rt, "LayerOptics")
      .def(
          "__init__",
          [](rt4::layer_optics* o, const Tensor4& extinction, const Tensor3& absorption, const Tensor6& phase) {
            new (o) rt4::layer_optics{.extinction = extinction, .absorption = absorption, .phase = phase};
          },
          "extinction"_a,
          "absorption"_a,
          "phase"_a)
      .def_rw("extinction",
              &rt4::layer_optics::extinction,
              "Extinction matrix K per unit length, [2 hemisphere, nmu_total, nstokes row, nstokes col]\n\n.. "
              ":class:`~pyarts3.arts.Tensor4`")
      .def_rw("absorption",
              &rt4::layer_optics::absorption,
              "Absorption vector a per unit length, [2 hemisphere, nmu_total, nstokes]; RT4 multiplies it by "
              "the Planck function\n\n.. :class:`~pyarts3.arts.Tensor3`")
      .def_rw("phase",
              &rt4::layer_optics::phase,
              "Azimuthal mean (1 / 2 pi) int Z dphi per unit length and steradian, including the number density "
              "and no quadrature weights, [2 out hemisphere, 2 in hemisphere, nmu_total out, nmu_total in, "
              "nstokes out, nstokes in]\n\n.. :class:`~pyarts3.arts.Tensor6`")
      .doc() = R"(Particle optics of one homogeneous layer on the solver streams.

The medium must be mirror symmetric between the hemispheres
(``extinction[down] == extinction[up]``, ``absorption[down] == absorption[up]``,
``phase[down, down] == phase[up, up]``, ``phase[down, up] == phase[up, down]``).
Energy conservation, K11(h, mu_j) = a1(h, mu_j) + 2 pi sum_i w_i
[Z(up <- h)(1, i; 1, j) + Z(down <- h)(1, i; 1, j)], is not enforced.)";

  auto aolo =
      py::bind_vector<std::vector<rt4::layer_optics>, py::rv_policy::reference_internal>(rt, "ArrayOfLayerOptics");
  aolo.doc() = "A list of LayerOptics";

  py::class_<rt4::lambertian_surface>(rt, "LambertianSurface")
      .def(
          "__init__",
          [](rt4::lambertian_surface* s, Numeric albedo) { new (s) rt4::lambertian_surface{.albedo = albedo}; },
          "albedo"_a = 0.0)
      .def_rw("albedo", &rt4::lambertian_surface::albedo, "Albedo A\n\n.. :class:`float`")
      .doc() =
      "RT4 'L': reflection 2 A mu_j w_j into every stream, I to I only; emission [(1 - A) B, 0].  Energy is "
      "conserved on the streams only for double_gauss quadrature.";

  py::class_<rt4::fresnel_surface>(rt, "FresnelSurface")
      .def(
          "__init__",
          [](rt4::fresnel_surface* s, Complex n) { new (s) rt4::fresnel_surface{.refractive_index = n}; },
          "refractive_index"_a)
      .def_rw("refractive_index",
              &rt4::fresnel_surface::refractive_index,
              "Complex refractive index of the surface (medium above has index 1)\n\n.. :class:`complex`")
      .doc() =
      "RT4 'F': specular Fresnel reflection R = [[R1, R2], [R2, R1]], R1 = (|r_v|^2 + |r_h|^2) / 2, "
      "R2 = (|r_v|^2 - |r_h|^2) / 2; emission [(1 - R1) B, -R2 B].";

  py::class_<rt4::specular_surface>(rt, "SpecularSurface")
      .def(
          "__init__",
          [](rt4::specular_surface* s, const Matrix& reflectivity) {
            new (s) rt4::specular_surface{.reflectivity = reflectivity};
          },
          "reflectivity"_a)
      .def_rw("reflectivity",
              &rt4::specular_surface::reflectivity,
              "[nstokes, nstokes] R(out, in)\n\n.. :class:`~pyarts3.arts.Matrix`")
      .doc() = "RT4 'S': reflectivity applied specularly to every stream; emission [(1 - R(I, I)) B, -R(Q, I) B].";

  py::class_<rt4::discrete_surface>(rt, "DiscreteSurface")
      .def(
          "__init__",
          [](rt4::discrete_surface* s, const Tensor4& reflection, const Matrix& emission) {
            new (s) rt4::discrete_surface{.reflection = reflection, .emission = emission};
          },
          "reflection"_a,
          "emission"_a)
      .def_rw("reflection",
              &rt4::discrete_surface::reflection,
              "[nmu_total out (up), nmu_total in (down), nstokes out, nstokes in], the discrete operator "
              "including quadrature factors\n\n.. :class:`~pyarts3.arts.Tensor4`")
      .def_rw("emission",
              &rt4::discrete_surface::emission,
              "[nmu_total, nstokes] upwelling emission in W m-2 Hz-1 sr-1\n\n.. :class:`~pyarts3.arts.Matrix`")
      .doc() = R"(RT4 'A': I_up(i) = sum_j reflection(i, j) I_down(j) + emission(i).

A Lambertian surface is reflection 2 A mu_j w_j.  Problem.surface_temperature
is not used for this surface.)";

  const rt4::problem d{};
  py::class_<rt4::problem>(rt, "Problem")
      .def(
          "__init__",
          [](rt4::problem*                         p,
             Index                                 nstokes,
             Index                                 nmu,
             rt4::quadrature_type                  quad,
             const Vector&                         extra_mu,
             Numeric                               max_delta_tau,
             Numeric                               frequency,
             const Vector&                         height,
             const Vector&                         temperature,
             const Vector&                         gas_extinction,
             const std::vector<rt4::layer_optics>& optics,
             const ArrayOfIndex&                   layer_optics_index,
             Numeric                               sky_temperature,
             Numeric                               surface_temperature,
             const rt4::surface&                   ground) {
            new (p) rt4::problem{.nstokes             = nstokes,
                                 .nmu                 = nmu,
                                 .quad                = quad,
                                 .extra_mu            = extra_mu,
                                 .max_delta_tau       = max_delta_tau,
                                 .frequency           = frequency,
                                 .height              = height,
                                 .temperature         = temperature,
                                 .gas_extinction      = gas_extinction,
                                 .optics              = optics,
                                 .layer_optics_index  = layer_optics_index,
                                 .sky_temperature     = sky_temperature,
                                 .surface_temperature = surface_temperature,
                                 .ground              = ground};
          },
          "nstokes"_a             = d.nstokes,
          "nmu"_a                 = d.nmu,
          "quad"_a                = d.quad,
          "extra_mu"_a            = d.extra_mu,
          "max_delta_tau"_a       = d.max_delta_tau,
          "frequency"_a           = d.frequency,
          "height"_a              = d.height,
          "temperature"_a         = d.temperature,
          "gas_extinction"_a      = d.gas_extinction,
          "optics"_a              = d.optics,
          "layer_optics_index"_a  = d.layer_optics_index,
          "sky_temperature"_a     = d.sky_temperature,
          "surface_temperature"_a = d.surface_temperature,
          "ground"_a              = d.ground)
      .def_rw("nstokes", &rt4::problem::nstokes, "1 for [I], 2 for [I, Q]\n\n.. :class:`int`")
      .def_rw("nmu", &rt4::problem::nmu, "Quadrature nodes per hemisphere\n\n.. :class:`int`")
      .def_rw("quad", &rt4::problem::quad, "Quadrature rule\n\n.. :class:`~pyarts3.arts.rt4.QuadratureType`")
      .def_rw("extra_mu",
              &rt4::problem::extra_mu,
              "Zero-weight output angles appended after the quadrature streams, each in (0, 1]\n\n.. "
              ":class:`~pyarts3.arts.Vector`")
      .def_rw("max_delta_tau",
              &rt4::problem::max_delta_tau,
              "Maximum vertical optical thickness of the initial doubling sublayer, > 0\n\n.. :class:`float`")
      .def_rw("frequency", &rt4::problem::frequency, "Frequency [Hz]\n\n.. :class:`float`")
      .def_rw("height",
              &rt4::problem::height,
              "[nlay + 1] layer interfaces, top-down; only |differences| are used, in the reciprocal of the "
              "extinction unit\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("temperature",
              &rt4::problem::temperature,
              "[nlay + 1] interface temperatures [K], top-down, > 0; the Planck function is linear within each "
              "layer\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("gas_extinction",
              &rt4::problem::gas_extinction,
              "[nlay] scalar, unpolarized gas extinction per unit length, >= 0\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw(
          "optics", &rt4::problem::optics, "Particle optics sets\n\n.. :class:`~pyarts3.arts.rt4.ArrayOfLayerOptics`")
      .def_rw("layer_optics_index",
              &rt4::problem::layer_optics_index,
              "[nlay] index into optics per layer, or < 0 for a gas-only layer\n\n.. "
              ":class:`~pyarts3.arts.ArrayOfIndex`")
      .def_rw("sky_temperature",
              &rt4::problem::sky_temperature,
              "Temperature [K] of the isotropic, unpolarized blackbody incident at the top\n\n.. :class:`float`")
      .def_rw("surface_temperature",
              &rt4::problem::surface_temperature,
              "Surface temperature [K] for the Lambertian, Fresnel and specular surfaces\n\n.. :class:`float`")
      .def_prop_rw(
          "ground",
          [](const rt4::problem& p) { return p.ground; },
          [](rt4::problem& p, const rt4::surface& g) { p.ground = g; },
          "The surface (a copy is returned; assign to change it)\n\n.. :class:`~pyarts3.arts.rt4.LambertianSurface` "
          "| :class:`~pyarts3.arts.rt4.FresnelSurface` | :class:`~pyarts3.arts.rt4.SpecularSurface` | "
          ":class:`~pyarts3.arts.rt4.DiscreteSurface`")
      .doc() = "An RT4 problem; see :doc:`dev.rt4` for the conventions";

  py::class_<rt4::result>(rt, "RT4Result")
      .def_ro("mu", &rt4::result::mu, "[nmu_total] streams\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_ro("weights",
              &rt4::result::weights,
              "[nmu_total] quadrature weights, 0 for the extra angles\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_ro("up",
              &rt4::result::up,
              "[nlay + 1 level (0 = top), nmu_total, nstokes] upward-propagating radiance in W m-2 Hz-1 sr-1\n\n.. "
              ":class:`~pyarts3.arts.Tensor3`")
      .def_ro("down",
              &rt4::result::down,
              "[nlay + 1 level (0 = top), nmu_total, nstokes] downward-propagating radiance in W m-2 Hz-1 "
              "sr-1\n\n.. :class:`~pyarts3.arts.Tensor3`")
      .doc() = "The RT4 solution at every level";

  rt.def("solve",
         &rt4::solve,
         "problem"_a,
         R"(Run RT4.

Validates the shapes, every precondition on which the Fortran code would stop
the process, and the mirror symmetry RT4 requires, then calls RADTRANO.
Calls are serialised by one global lock.)",
         py::call_guard<py::gil_scoped_release>());

  rt.def("scattering_optics",
         &rt4::scattering_optics,
         "scattering_species"_a,
         "atm_point"_a,
         "frequency"_a,
         "mu"_a,
         "nstokes"_a       = 2,
         "azimuth_count"_a = 64,
         R"(The particle optics of ARTS scattering species at one atmospheric point on RT4's streams.

Uses ARTS's laboratory-frame (ARO gridded) bulk scattering properties, so
azimuthally randomly oriented species work as well as totally randomly
oriented ones.  RT4's stream ``(down, mu)`` is ARTS's propagation zenith
angle ``180 - acos(mu)`` deg, ``(up, mu)`` is ``acos(mu)``.

Parameters
----------
scattering_species : ~pyarts3.arts.ArrayOfScatteringSpecies
atm_point : ~pyarts3.arts.AtmPoint
frequency : float
    [Hz].
mu : ~pyarts3.arts.Vector
    The stream cosines of one hemisphere, each in (0, 1]: the quadrature
    nodes followed by the extra angles.
nstokes : int
    1 or 2.
azimuth_count : int
    Even number N of azimuth differences of the midpoint rule for the
    azimuthal mean, at (k + 1/2) 360 / N deg; exact for N > L when the
    scattering matrix is a regular Legendre series of degree L.

Returns
-------
LayerOptics
    Extinction ([[K11, K12], [K12, K11]]), absorption ([a1, a2]) and the
    azimuthal mean of the [I, Q] block of the phase matrix, per metre and
    steradian.  ARTS's GasScatterer and HenyeyGreensteinScatterer
    interpolate their laboratory-frame data on a 1 deg grid of scattering
    angles (an error of up to 5.7e-5 sigma / (4 pi) for Rayleigh).
)");

  const rt4::path_settings ds{};
  py::class_<rt4::path_settings>(rt, "PathSettings")
      .def(
          "__init__",
          [](rt4::path_settings*  s,
             Index                nstokes,
             Index                nmu,
             rt4::quadrature_type quad,
             const Vector&        extra_mu,
             Numeric              max_delta_tau,
             Index                azimuth_count) {
            new (s) rt4::path_settings{.nstokes       = nstokes,
                                       .nmu           = nmu,
                                       .quad          = quad,
                                       .extra_mu      = extra_mu,
                                       .max_delta_tau = max_delta_tau,
                                       .azimuth_count = azimuth_count};
          },
          "nstokes"_a       = ds.nstokes,
          "nmu"_a           = ds.nmu,
          "quad"_a          = ds.quad,
          "extra_mu"_a      = ds.extra_mu,
          "max_delta_tau"_a = ds.max_delta_tau,
          "azimuth_count"_a = ds.azimuth_count)
      .def_rw("nstokes", &rt4::path_settings::nstokes, "1 for [I], 2 for [I, Q]\n\n.. :class:`int`")
      .def_rw("nmu", &rt4::path_settings::nmu, "Quadrature nodes per hemisphere\n\n.. :class:`int`")
      .def_rw("quad", &rt4::path_settings::quad, "Quadrature rule\n\n.. :class:`~pyarts3.arts.rt4.QuadratureType`")
      .def_rw("extra_mu",
              &rt4::path_settings::extra_mu,
              "Zero-weight output angles, each in (0, 1]\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("max_delta_tau",
              &rt4::path_settings::max_delta_tau,
              "Maximum vertical optical thickness of the initial doubling sublayer\n\n.. :class:`float`")
      .def_rw("azimuth_count",
              &rt4::path_settings::azimuth_count,
              "Azimuth differences of the azimuthal mean, even\n\n.. :class:`int`")
      .doc() = "Solver settings of problem_from_path";

  rt.def("problem_from_path",
         &rt4::problem_from_path,
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
         R"(An RT4 problem from an ARTS propagation path, with the DISORT path conventions.

``ray_path``, ``atm_path`` and ``spectral_propmat_path`` have one entry per
level, top first, with strictly decreasing altitudes (the heights, in m).
The level temperatures are ``atm_path``'s.  ``spectral_propmat_path`` is the
unpolarized gas propagation matrix (no particles) per metre; a layer's gas
extinction is the mean of A at its two levels.  The frequency is
``freq_grid[freq_index]``.  Each layer gets the mean of its two levels'
:func:`scattering_optics`, or is gas-only when that is all zero.  The
streams are RT4's quadrature for ``settings``.  See :doc:`dev.rt4`.
)");
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize rt4\n{}", e.what()));
}
}  // namespace Python
