#include <disort-brdf.h>
#include <disort.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <pydocs.h>
#include <python_interface.h>
#include <vdisort-brdf.h>
#include <vdisort.h>
#include <vdisort_arts.h>
#include <vdisort_settings.h>

#include <concepts>
#include <limits>
#include <optional>
#include <ranges>

#include "configtypes.h"
#include "debug.h"
#include "hpy_arts.h"
#include "operators.h"
#include "sorting.h"

namespace Python {

void py_disort(py::module_& m) try {
  auto disort_nm  = m.def_submodule("disort");
  disort_nm.doc() = "DISORT solver internal types";

  disort_nm.def(
      "delta_m_plus",
      [](const Matrix& phase_moments, const Index nleg) {
        auto scaling = disort::delta_m_plus(phase_moments, nleg);
        return py::make_tuple(std::move(scaling.fraction), std::move(scaling.moments));
      },
      "phase_moments"_a,
      "nleg"_a,
      "Construct DISORT 4 delta-M-plus fractions and removed-peak moments");

  // VDISORT PYTHON INTERFACE BEGIN: polarized solver namespace and constants.
  auto vdisort_nm                     = m.def_submodule("vdisort");
  vdisort_nm.doc()                    = "VDISORT polarized solver internal types";
  vdisort_nm.attr("stokes_dimension") = vdisort::stokes_dimension;
  vdisort_nm.attr("cosine_mode")      = vdisort::cosine_mode;
  vdisort_nm.attr("sine_mode")        = vdisort::sine_mode;
  // VDISORT PYTHON INTERFACE END

  py::class_<disort::surface_reflection> dsr(disort_nm, "SurfaceReflection");
  dsr.def(py::init<>(), "No reflection")
      .def(py::init<Index, Index, Index, Numeric>(),
           "NBRDF"_a,
           "n_out"_a,
           "N"_a,
           "beam_mu"_a = std::numeric_limits<Numeric>::quiet_NaN(),
           "Zero reflection of NBRDF modes from N streams into n_out directions, for the beam cosine beam_mu (NaN: "
           "none)")
      .def(py::init<Tensor3, Matrix, Numeric>(),
           "brdf"_a,
           "beam_brdf"_a,
           "beam_mu"_a = std::numeric_limits<Numeric>::quiet_NaN(),
           "The reflection brdf [NBRDF, n_out, N] and beam_brdf [NBRDF, n_out], for the beam cosine beam_mu (NaN: "
           "none)")
      .def("like",
           &disort::surface_reflection::like,
           "n_out"_a,
           "Zero reflection of this one's modes from its streams into n_out directions, for its beam cosine")
      .def(
          "set",
          [](disort::surface_reflection& x, const Tensor3& brdf, const Matrix& beam_brdf, Numeric beam_mu) {
            x.set(brdf, beam_brdf, beam_mu);
          },
          "brdf"_a,
          "beam_brdf"_a,
          "beam_mu"_a,
          "Set the reflection in place, for the beam cosine beam_mu (NaN: none); the shapes must be this one's")
      .def_prop_ro(
          "brdf",
          [](const disort::surface_reflection& x) { return Tensor3{x.brdf()}; },
          "[NBRDF, n_out, N] Fourier modes R_m(mu_out, mu_j) of the BRDF from the downward streams, without the "
          "solver's quadrature factor (a copy)\n\n.. :class:`~pyarts3.arts.Tensor3`")
      .def_prop_ro(
          "beam_brdf",
          [](const disort::surface_reflection& x) { return Matrix{x.beam_brdf()}; },
          "[NBRDF, n_out] the same from the direct beam (a copy)\n\n.. :class:`~pyarts3.arts.Matrix`")
      .def_prop_ro("beam_mu",
                   &disort::surface_reflection::beam_mu,
                   "The beam cosine of beam_brdf, NaN when computed without a beam\n\n.. :class:`float`")
      .def_prop_ro("modes", &disort::surface_reflection::modes, "The number of Fourier modes, NBRDF\n\n.. :class:`int`")
      .def_prop_ro("directions",
                   &disort::surface_reflection::directions,
                   "The number of upward directions reflected into, n_out\n\n.. :class:`int`")
      .def_prop_ro("streams",
                   &disort::surface_reflection::streams,
                   "The number of downward streams reflected from, N\n\n.. :class:`int`");
  dsr.doc() = R"(A surface's reflection into upward directions, in DISORT's Fourier convention.

Its shapes are those it is constructed with, and its data is set in place
with ``set``.  On the solver's upward streams it is the solver's surface; at
user directions it is the reflection the user-angle radiances need there,
for which ``like`` gives one with the solver surface's modes, streams and
beam cosine.)";

  disort_nm.def(
      "lambertian_fourier_modes",
      [](Numeric albedo, Index number_of_modes, const Vector& mu_out, const Vector& mu_in, Numeric mu0) {
        return disort::brdf::lambertian_fourier_modes(albedo, number_of_modes, mu_out, mu_in, mu0);
      },
      "albedo"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a = std::numeric_limits<Numeric>::quiet_NaN(),
      R"(The exact scalar Lambertian surface: R_0 = A, higher modes zero.

``mu_out`` are the cosines of the upward directions reflected into (the
solver's upward streams, or user directions), ``mu_in`` those of the solver's
downward streams reflected from (their magnitudes are used), and ``mu0`` the
cosine of the direct beam, NaN without a beam.)");
  disort_nm.def(
      "fourier_modes",
      [](const disort::brdf::RawFunction& raw,
         Index                            number_of_modes,
         const Vector&                    mu_out,
         const Vector&                    mu_in,
         Numeric                          mu0,
         Index                            azimuth_quadrature_points) {
        return disort::brdf::fourier_modes(raw, number_of_modes, mu_out, mu_in, mu0, azimuth_quadrature_points);
      },
      "raw"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      R"(The cosine Fourier modes of a scalar BRDF ``raw(mu_out, mu_in, relative_azimuth)`` per steradian.

``mu_out`` are the cosines of the upward directions reflected into (the
solver's upward streams, or user directions), ``mu_in`` those of the solver's
downward streams reflected from (their magnitudes are used), and ``mu0`` the
cosine of the direct beam, NaN without a beam.)");
  disort_nm.def("combine_fourier_modes",
                &disort::brdf::combine_fourier_modes,
                "first"_a,
                "first_weight"_a,
                "second"_a,
                "second_weight"_a,
                "The weighted sum of two scalar surface reflections on the same directions and for the same beam");
  disort_nm.def(
      "cox_munk_lambertian_fourier_modes",
      [](Numeric       cox_munk_fraction,
         Numeric       lambertian_albedo,
         Numeric       wind_speed,
         Complex       refractive_index,
         bool          shadowing,
         Index         number_of_modes,
         const Vector& mu_out,
         const Vector& mu_in,
         Numeric       mu0,
         Index         azimuth_quadrature_points) {
        return disort::brdf::cox_munk_lambertian_fourier_modes(cox_munk_fraction,
                                                               lambertian_albedo,
                                                               wind_speed,
                                                               refractive_index,
                                                               shadowing,
                                                               number_of_modes,
                                                               mu_out,
                                                               mu_in,
                                                               mu0,
                                                               azimuth_quadrature_points);
      },
      "cox_munk_fraction"_a,
      "lambertian_albedo"_a,
      "wind_speed"_a,
      "refractive_index"_a,
      "shadowing"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      R"(A scalar Cox-Munk/Lambertian surface mixture.

``mu_out`` are the cosines of the upward directions reflected into (the
solver's upward streams, or user directions), ``mu_in`` those of the solver's
downward streams reflected from (their magnitudes are used), and ``mu0`` the
cosine of the direct beam, NaN without a beam.)");

  // VDISORT PYTHON INTERFACE BEGIN: a surface's reflection on given directions
  py::class_<vdisort::surface_reflection> vsr(vdisort_nm, "SurfaceReflection");
  vsr.def(py::init<>(), "No reflection")
      .def(py::init<Index, Index, Index, Numeric, bool>(),
           "NBRDF"_a,
           "n_out"_a,
           "N"_a,
           "beam_mu"_a  = std::numeric_limits<Numeric>::quiet_NaN(),
           "specular"_a = false,
           "Zero reflection of NBRDF modes from N streams into n_out directions, for the beam cosine beam_mu (NaN: "
           "none), with a (zero) specular part if specular")
      .def(py::init<rtepack::muelmat_tensor4, rtepack::muelmat_tensor3, Numeric, rtepack::muelmat_vector>(),
           "brdf"_a,
           "beam_brdf"_a,
           "beam_mu"_a  = std::numeric_limits<Numeric>::quiet_NaN(),
           "specular"_a = rtepack::muelmat_vector{},
           "The reflection brdf [2, NBRDF, n_out, N] and beam_brdf [2, NBRDF, n_out], for the beam cosine beam_mu "
           "(NaN: none), and the specular part [n_out] or none")
      .def("like",
           &vdisort::surface_reflection::like,
           "n_out"_a,
           "Zero reflection of this one's modes from its streams into n_out directions, for its beam cosine, with a "
           "(zero) specular part if it has one")
      .def(
          "set",
          [](vdisort::surface_reflection&    x,
             const rtepack::muelmat_tensor4& brdf,
             const rtepack::muelmat_tensor3& beam_brdf,
             Numeric                         beam_mu,
             const rtepack::muelmat_vector&  specular) { x.set(brdf, beam_brdf, beam_mu, specular); },
          "brdf"_a,
          "beam_brdf"_a,
          "beam_mu"_a,
          "specular"_a = rtepack::muelmat_vector{},
          "Set the reflection in place, for the beam cosine beam_mu (NaN: none); the shapes must be this one's")
      .def_prop_ro(
          "brdf",
          [](const vdisort::surface_reflection& x) { return rtepack::muelmat_tensor4{x.brdf()}; },
          "[2, NBRDF, n_out, N] combined cosine and sine Fourier modes of the polarized BRDF from the downward "
          "streams, before the pi w mu quadrature factor (a copy)\n\n.. :class:`~pyarts3.arts.MuelmatTensor4`")
      .def_prop_ro(
          "beam_brdf",
          [](const vdisort::surface_reflection& x) { return rtepack::muelmat_tensor3{x.beam_brdf()}; },
          "[2, NBRDF, n_out] the same from the direct beam (a copy)\n\n.. :class:`~pyarts3.arts.MuelmatTensor3`")
      .def_prop_ro("beam_mu",
                   &vdisort::surface_reflection::beam_mu,
                   "The beam cosine of beam_brdf, NaN when computed without a beam\n\n.. :class:`float`")
      .def_prop_ro(
          "specular",
          [](const vdisort::surface_reflection& x) { return rtepack::muelmat_vector{x.specular()}; },
          "[n_out] the specular (mirror) part R(mu_out), or empty (a copy)\n\n.. :class:`~pyarts3.arts.MuelmatVector`")
      .def_prop_ro("modes",
                   &vdisort::surface_reflection::modes,
                   "The number of Fourier modes of the BRDF, NBRDF\n\n.. :class:`int`")
      .def_prop_ro("directions",
                   &vdisort::surface_reflection::directions,
                   "The number of upward directions reflected into, n_out\n\n.. :class:`int`")
      .def_prop_ro("streams",
                   &vdisort::surface_reflection::streams,
                   "The number of downward streams reflected from, N\n\n.. :class:`int`")
      .def_prop_ro("has_specular",
                   &vdisort::surface_reflection::has_specular,
                   "Whether there is a specular part\n\n.. :class:`bool`");
  vsr.doc() = R"(A surface's reflection into upward directions, in VDISORT's combined representation.

Its shapes, and whether there is a specular part, are those it is
constructed with, and its data is set in place with ``set``.  On the
solver's upward streams it is the solver's surface; at user directions it is
the reflection the user-angle radiances need there, for which ``like`` gives
one with the solver surface's modes, streams, beam cosine and specular part.)";

  py::class_<vdisort::delta_m_transport_data> delta_m_transport(vdisort_nm, "DeltaMTransportData");
  delta_m_transport
      .def_rw("tau",
              &vdisort::delta_m_transport_data::tau,
              "Delta-M optical-depth grid\n\n.. :class:`~pyarts3.arts.AscendingGrid`")
      .def_rw("omega",
              &vdisort::delta_m_transport_data::omega,
              "Delta-M single-scattering albedo\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("phase_matrix",
              &vdisort::delta_m_transport_data::phase_matrix,
              "Delta-M diffuse phase matrices\n\n.. :class:`~pyarts3.arts.MuelmatTensor5`")
      .def_rw("beam_phase_matrix",
              &vdisort::delta_m_transport_data::beam_phase_matrix,
              "Delta-M beam phase matrices\n\n.. :class:`~pyarts3.arts.MuelmatTensor4`")
      .def_rw("source_coordinate_scale",
              &vdisort::delta_m_transport_data::source_coordinate_scale,
              "Scale of the affine physical-source coordinate map\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("source_coordinate_offset",
              &vdisort::delta_m_transport_data::source_coordinate_offset,
              "Offset of the affine physical-source coordinate map\n\n.. :class:`~pyarts3.arts.Vector`");
  delta_m_transport.doc() = "Solver-ready result of an explicitly specified polarized delta-M transform";

  vdisort_nm.def("combine_phase_matrices",
                 &vdisort::combine_phase_matrices,
                 "cosine"_a,
                 "sine"_a,
                 "Convert ordinary cosine/sine phase coefficients to the combined VDISORT representation.");
  vdisort_nm.def("combine_beam_phase_matrices",
                 &vdisort::combine_beam_phase_matrices,
                 "cosine"_a,
                 "sine"_a,
                 "Convert ordinary cosine/sine beam phase coefficients to the combined VDISORT representation.");
  vdisort_nm.def("delta_m_preprocess",
                 &vdisort::delta_m_preprocess,
                 "physical_tau"_a,
                 "physical_omega"_a,
                 "fraction"_a,
                 "original_phase_matrix"_a,
                 "removed_phase_matrix"_a,
                 "original_beam_phase_matrix"_a = vdisort::beam_phase_matrix_data{},
                 "removed_beam_phase_matrix"_a  = vdisort::beam_phase_matrix_data{},
                 "Apply a caller-defined polarized delta-M split and return solver-ready transport inputs");
  vdisort_nm.def(
      "cox_munk_reflection",
      [](const Numeric outgoing_mu,
         const Numeric incoming_mu,
         const Numeric relative_azimuth,
         const Numeric wind_speed,
         const Complex refractive_index,
         const bool    shadowing) {
        return vdisort::brdf::CoxMunk{wind_speed, refractive_index, shadowing}(
            outgoing_mu, incoming_mu, relative_azimuth);
      },
      "outgoing_mu"_a,
      "incoming_mu"_a,
      "relative_azimuth"_a,
      "wind_speed"_a       = 5.0,
      "refractive_index"_a = Complex{1.34, 0.0},
      "shadowing"_a        = true,
      "Evaluate the raw polarized Cox-Munk BPrDF in the ARTS Stokes basis");
  vdisort_nm.def(
      "cox_munk_fourier_modes",
      [](Numeric       wind_speed,
         Complex       refractive_index,
         bool          shadowing,
         Index         number_of_modes,
         const Vector& mu_out,
         const Vector& mu_in,
         Numeric       mu0,
         Index         azimuth_quadrature_points) {
        return vdisort::brdf::cox_munk_fourier_modes(
            wind_speed, refractive_index, shadowing, number_of_modes, mu_out, mu_in, mu0, azimuth_quadrature_points);
      },
      "wind_speed"_a,
      "refractive_index"_a,
      "shadowing"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      R"(A polarized Cox-Munk ocean.

``mu_out`` are the cosines of the upward directions reflected into (the
solver's upward streams, or user directions), ``mu_in`` those of the solver's
downward streams reflected from (their magnitudes are used), and ``mu0`` the
cosine of the direct beam, NaN without a beam.)");
  vdisort_nm.def(
      "fresnel_reflection",
      [](const Numeric incident_mu, const Complex refractive_index) {
        return vdisort::brdf::Fresnel{refractive_index}(incident_mu);
      },
      "incident_mu"_a,
      "refractive_index"_a = Complex{1.5, 0.0},
      "Evaluate the polarized Fresnel reflection matrix of a flat dielectric interface");
  vdisort_nm.def(
      "fresnel_fourier_modes",
      [](Complex refractive_index, Index number_of_modes, const Vector& mu_out, const Vector& mu_in, Numeric mu0) {
        return vdisort::brdf::fresnel_fourier_modes(refractive_index, number_of_modes, mu_out, mu_in, mu0);
      },
      "refractive_index"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a = std::numeric_limits<Numeric>::quiet_NaN(),
      R"(An ideal flat Fresnel surface: no reflection kernel, the Fresnel matrix as the specular part.

It refuses a beam.  ``mu_out`` are the cosines of the upward directions reflected into (the
solver's upward streams, or user directions), ``mu_in`` those of the solver's
downward streams reflected from (their magnitudes are used), and ``mu0`` the
cosine of the direct beam, NaN without a beam.)");
  vdisort_nm.def(
      "lambertian_fourier_modes",
      [](Numeric albedo, Index number_of_modes, const Vector& mu_out, const Vector& mu_in, Numeric mu0) {
        return vdisort::brdf::lambertian_fourier_modes(albedo, number_of_modes, mu_out, mu_in, mu0);
      },
      "albedo"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a = std::numeric_limits<Numeric>::quiet_NaN(),
      R"(The exact, fully depolarizing Lambertian surface.

``mu_out`` are the cosines of the upward directions reflected into (the
solver's upward streams, or user directions), ``mu_in`` those of the solver's
downward streams reflected from (their magnitudes are used), and ``mu0`` the
cosine of the direct beam, NaN without a beam.)");
  vdisort_nm.def(
      "fourier_modes",
      [](const vdisort::brdf::RawFunction& raw,
         Index                             number_of_modes,
         const Vector&                     mu_out,
         const Vector&                     mu_in,
         Numeric                           mu0,
         Index                             azimuth_quadrature_points) {
        return vdisort::brdf::fourier_modes(raw, number_of_modes, mu_out, mu_in, mu0, azimuth_quadrature_points);
      },
      "raw"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      R"(The combined Fourier modes of a polarized BPrDF ``raw(mu_out, mu_in, relative_azimuth)`` per steradian.

``mu_out`` are the cosines of the upward directions reflected into (the
solver's upward streams, or user directions), ``mu_in`` those of the solver's
downward streams reflected from (their magnitudes are used), and ``mu0`` the
cosine of the direct beam, NaN without a beam.)");
  vdisort_nm.def(
      "depolarizing_fourier_modes",
      [](const vdisort::brdf::ScalarRawFunction& raw,
         Index                                   number_of_modes,
         const Vector&                           mu_out,
         const Vector&                           mu_in,
         Numeric                                 mu0,
         Index                                   azimuth_quadrature_points) {
        return vdisort::brdf::depolarizing_fourier_modes(
            raw, number_of_modes, mu_out, mu_in, mu0, azimuth_quadrature_points);
      },
      "raw"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      R"(A scalar BRDF ``raw(mu_out, mu_in, relative_azimuth)`` as a fully depolarizing surface.

``mu_out`` are the cosines of the upward directions reflected into (the
solver's upward streams, or user directions), ``mu_in`` those of the solver's
downward streams reflected from (their magnitudes are used), and ``mu0`` the
cosine of the direct beam, NaN without a beam.)");
  vdisort_nm.def(
      "hapke_fourier_modes",
      [](Numeric       opposition_amplitude,
         Numeric       opposition_width,
         Numeric       single_scattering_albedo,
         Index         number_of_modes,
         const Vector& mu_out,
         const Vector& mu_in,
         Numeric       mu0,
         Index         azimuth_quadrature_points) {
        return vdisort::brdf::hapke_fourier_modes(opposition_amplitude,
                                                  opposition_width,
                                                  single_scattering_albedo,
                                                  number_of_modes,
                                                  mu_out,
                                                  mu_in,
                                                  mu0,
                                                  azimuth_quadrature_points);
      },
      "opposition_amplitude"_a,
      "opposition_width"_a,
      "single_scattering_albedo"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      "A fully depolarizing Hapke surface; arguments after number_of_modes as for fourier_modes");
  vdisort_nm.def(
      "rpv_fourier_modes",
      [](Numeric       rho0,
         Numeric       kappa,
         Numeric       asymmetry,
         Numeric       hotspot,
         Index         number_of_modes,
         const Vector& mu_out,
         const Vector& mu_in,
         Numeric       mu0,
         Index         azimuth_quadrature_points) {
        return vdisort::brdf::rpv_fourier_modes(
            rho0, kappa, asymmetry, hotspot, number_of_modes, mu_out, mu_in, mu0, azimuth_quadrature_points);
      },
      "rho0"_a,
      "kappa"_a,
      "asymmetry"_a,
      "hotspot"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      "A fully depolarizing RPV surface; arguments after number_of_modes as for fourier_modes");
  vdisort_nm.def(
      "ross_li_fourier_modes",
      [](Numeric       isotropic,
         Numeric       volumetric,
         Numeric       geometric,
         Numeric       hotspot_angle,
         Index         number_of_modes,
         const Vector& mu_out,
         const Vector& mu_in,
         Numeric       mu0,
         Index         azimuth_quadrature_points) {
        return vdisort::brdf::ross_li_fourier_modes(isotropic,
                                                    volumetric,
                                                    geometric,
                                                    hotspot_angle,
                                                    number_of_modes,
                                                    mu_out,
                                                    mu_in,
                                                    mu0,
                                                    azimuth_quadrature_points);
      },
      "isotropic"_a,
      "volumetric"_a,
      "geometric"_a,
      "hotspot_angle"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      "A fully depolarizing Ross-Li surface; arguments after number_of_modes as for fourier_modes");
  vdisort_nm.def("combine_fourier_modes",
                 &vdisort::brdf::combine_fourier_modes,
                 "first"_a,
                 "first_weight"_a,
                 "second"_a,
                 "second_weight"_a,
                 "The weighted sum of two polarized surface reflections on the same directions and for the same beam");
  vdisort_nm.def(
      "fresnel_lambertian_fourier_modes",
      [](Numeric       fresnel_fraction,
         Numeric       lambertian_albedo,
         Complex       refractive_index,
         Index         number_of_modes,
         const Vector& mu_out,
         const Vector& mu_in,
         Numeric       mu0) {
        return vdisort::brdf::fresnel_lambertian_fourier_modes(
            fresnel_fraction, lambertian_albedo, refractive_index, number_of_modes, mu_out, mu_in, mu0);
      },
      "fresnel_fraction"_a,
      "lambertian_albedo"_a,
      "refractive_index"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a = std::numeric_limits<Numeric>::quiet_NaN(),
      "An ideal Fresnel/depolarizing-Lambertian mixture; it refuses a beam");
  vdisort_nm.def(
      "cox_munk_lambertian_fourier_modes",
      [](Numeric       cox_munk_fraction,
         Numeric       lambertian_albedo,
         Numeric       wind_speed,
         Complex       refractive_index,
         bool          shadowing,
         Index         number_of_modes,
         const Vector& mu_out,
         const Vector& mu_in,
         Numeric       mu0,
         Index         azimuth_quadrature_points) {
        return vdisort::brdf::cox_munk_lambertian_fourier_modes(cox_munk_fraction,
                                                                lambertian_albedo,
                                                                wind_speed,
                                                                refractive_index,
                                                                shadowing,
                                                                number_of_modes,
                                                                mu_out,
                                                                mu_in,
                                                                mu0,
                                                                azimuth_quadrature_points);
      },
      "cox_munk_fraction"_a,
      "lambertian_albedo"_a,
      "wind_speed"_a,
      "refractive_index"_a,
      "shadowing"_a,
      "number_of_modes"_a,
      "mu_out"_a,
      "mu_in"_a,
      "mu0"_a                       = std::numeric_limits<Numeric>::quiet_NaN(),
      "azimuth_quadrature_points"_a = 100,
      "A rough Fresnel Cox-Munk/depolarizing-Lambertian mixture; arguments as for fourier_modes");
  // VDISORT PYTHON INTERFACE END

  py::class_<disort::coupling_result> coupling_result(disort_nm, "CouplingResult");
  generic_interface(coupling_result);
  coupling_result
      .def_rw("iterations",
              &disort::coupling_result::iterations,
              "Number of fixed-point iterations\n\n.. :class:`~pyarts3.arts.Index`")
      .def_rw("max_relative_change",
              &disort::coupling_result::max_relative_change,
              "Maximum relative interface update in the last iteration\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("converged",
              &disort::coupling_result::converged,
              "Whether the interface exchange converged\n\n.. :class:`bool`");
  coupling_result.doc() = "The result of a DISORT interface coupling";

  disort_nm.def("couple",
                &disort::couple,
                "atmosphere"_a,
                "subsurface"_a,
                "tolerance"_a      = 1e-6,
                "max_iterations"_a = 16,
                "relaxation"_a     = 1.0,
                "Iteratively exchange DISORT interface boundary conditions.");

  py::class_<disort::main_data> x(m, "cppdisort");
  x.doc() = unwrap_stars(R"(A DISORT object.

This offers a low level interface to the DISORT solver.  See *DisortSettings*
for a higher level interface.  Especially, see the workspace variables for the
type as the workspace methods that operate on them explains the interface on a
higher level.

The implementation is based on the Pythonic-DISORT implementation, which is
a from scratch reimplementation of DISORT in Python.  The interface here is
mostly mimicking the Pythonic-DISORT interface, with some exceptions to
improve performance and usability.

The two main differences are that we use a custom Legendre-Gauss quadrature
implementation, that we use the BandMatrix LAPACK solver for the left-hand
side of the linear system, and that we use a pure real eigenvalue solver
for the matrix decomposition that's been ported and optimized in C++.

Unlike Pythonic-DISORT, ``b_pos`` is not added to the surface reflection: it
is the radiance beneath the surface, which the surface transmits with
``1 - R``, ``R`` the reflection of ``surface`` on the streams.  For
a Planck radiance ``B(T)`` in mode 0 this is Kirchhoff's emission
``(1 - albedo(mu)) B(T)``, as original DISORT's ``BTEMP``.  Without a surface
``b_pos`` is the upward radiance itself.

.. warning::

    The DISORT implementation is still being tested.  Initial results look
    promising, but please report any issues you find.  Initial tests show
    that the implementation is about 6x faster than CDISORT.  We do not
    have numbers on the performance compared to Pythonic-DISORT because
    Pythonic-DISORT is not optimized for speed.

.. warning::

    The internals of this implementation calls LAPACK routines.  Please
    ensure that your LAPACK installation is either single threaded or uses
    OpenMP.  Mixing multiple threading implementations will lead to
    significant slowdowns (or a complete stall of the program).

The relevant references are:

- Pythonic-DISORT: :cite:t:`Ho2024`
- Original DISORT: :cite:t:`Stamnes88`
- Legendre-Gauss quadrature: :cite:t:`Bogaert2014`
- BandMatrix solver: :cite:t:`Barrett1994`
- Real eigenvalue solver (original sources, the executed code is ported to C++): :cite:t:`buras2011`, :cite:t:`Dongarra1984`, :cite:t:`Parlett1969`, :cite:t:`Mitchell1967`
)");
  x.def(
      "__init__",
      [](disort::main_data*             n,
         const AscendingGrid&           tau_arr,
         const Vector&                  omega_arr,
         const Index                    NQuad,
         const Matrix&                  Leg_coeffs_all,
         Numeric                        mu0,
         Numeric                        I0,
         Numeric                        phi0,
         const std::optional<Index>     NLeg_,
         const std::optional<Index>     NFourier_,
         const std::optional<Matrix>&   b_pos,
         const std::optional<Matrix>&   b_neg,
         const std::optional<Vector>&   f_arr,
         const std::optional<disort::surface_reflection>& surface,
         const std::optional<Matrix>&   s_poly_coeffs,
         const std::optional<Matrix>&   delta_m_peak) {
        const Index NFourier = NFourier_.value_or(NQuad);
        const Index NLeg     = NLeg_.value_or(NQuad);
        const Index NLayers  = tau_arr.size();

        new (n) disort::main_data(NQuad,
                                  NLeg,
                                  NFourier,
                                  tau_arr,
                                  omega_arr,
                                  Leg_coeffs_all,
                                  b_pos.value_or(Matrix(NFourier, NQuad / 2, 0.0)),
                                  b_neg.value_or(Matrix(NFourier, NQuad / 2, 0.0)),
                                  f_arr.value_or(Vector(NLayers, 0.0)),
                                  s_poly_coeffs.value_or(Matrix(NLayers, 0, 0.0)),
                                  surface.value_or(disort::surface_reflection{}),
                                  mu0,
                                  I0,
                                  phi0,
                                  delta_m_peak.value_or(Matrix{}));
      },
      "Run disort, mostly mimicying the 0.7 Pythonic-DISORT interface.\n",
      "tau_arr"_a,
      "omega_arr"_a,
      "NQuad"_a,
      "Leg_coeffs_all"_a,
      "mu0"_a,
      "I0"_a,
      "phi0"_a,
      "NLeg"_a.none()          = py::none(),
      "NFourier"_a.none()      = py::none(),
      "b_pos"_a.none()         = py::none(),
      "b_neg"_a.none()         = py::none(),
      "f_arr"_a.none()         = py::none(),
      "surface"_a.none()       = py::none(),
      "s_poly_coeffs"_a.none() = py::none(),
      "delta_m_peak"_a.none()  = py::none());
  x.def(
       "u",
       [](disort::main_data& dis, const AscendingGrid& tau, const Vector& phi) {
         Tensor3 out(tau.size(), phi.size(), dis.mu().size());
         dis.ungridded_u(out, tau, phi);
         return out;
       },
       "tau"_a,
       "phi"_a,
       "Compute the intensity")
      .def(
          "u_user",
          [](disort::main_data&                               dis,
             const Vector&                                    mu,
             const AscendingGrid&                             tau,
             const Vector&                                    phi,
             const std::optional<disort::surface_reflection>& surface) {
            Tensor3             out(mu.size(), tau.size(), phi.size());
            disort::user_u_data data;
            for (Size t = 0; t < tau.size(); ++t)
              for (Size p = 0; p < phi.size(); ++p) {
                dis.u_user(data, tau[t], phi[p], mu, surface.value_or(disort::surface_reflection{}));
                out[joker, t, p] = data.intensities;
              }
            return out;
          },
          "mu"_a,
          "tau"_a,
          "phi"_a,
          "surface"_a.none() = py::none(),
          "Compute intensity at user polar-angle cosines using DISORT source reconstruction and formal ray integration")
      .def(
          "u_user_corr",
          [](disort::main_data&                               dis,
             const Vector&                                    mu,
             const AscendingGrid&                             tau,
             const Vector&                                    phi,
             const std::optional<disort::surface_reflection>& surface) {
            Tensor3             out(mu.size(), tau.size(), phi.size());
            disort::user_u_data data;
            disort::tms_data    tms;
            Vector              ims;
            for (Size t = 0; t < tau.size(); ++t)
              for (Size p = 0; p < phi.size(); ++p) {
                dis.u_user_corr(data, ims, tms, tau[t], phi[p], mu, surface.value_or(disort::surface_reflection{}));
                out[joker, t, p] = data.intensities;
              }
            return out;
          },
          "mu"_a,
          "tau"_a,
          "phi"_a,
          "surface"_a.none() = py::none(),
          "Compute IMS/TMS-corrected intensity at user polar-angle cosines")
      .def(
          "flux",
          [](disort::main_data& dis, const AscendingGrid& tau) {
            Matrix out(4, tau.size());
            dis.ungridded_flux(out[0], out[1], out[2], out[3], tau);
            return out;
          },
          "tau"_a,
          "Compute upward, downward-diffuse, downward-direct flux and DFDT")
      .def(
          "pydisort_u",
          [](disort::main_data& dis, Vector tau_, const Vector& phi) {
            std::vector<Index> sorting(tau_.size());
            stdr::iota(sorting, 0);
            stdr::sort(stdv::zip(sorting, tau_), {}, [](const auto& x) { return std::get<1>(x); });

            AscendingGrid tau{std::move(tau_)};
            Tensor3       res(tau.size(), phi.size(), dis.mu().size());
            dis.ungridded_u(res, tau, phi);

            Tensor3 out(dis.mu().size(), tau.size(), phi.size());
            for (Size i = 0; i < tau.size(); i++) { out[joker, i, joker] = transpose(res[sorting[i]]); }
            return out;
          },
          "tau"_a,
          "phi"_a,
          "Compute the intensity");
  generic_interface(x);

  // VDISORT PYTHON INTERFACE BEGIN: low-level polarized counterpart of
  // cppdisort.  Radiation fields retain the scalar axes and append Stokes.
  py::class_<vdisort::main_data> vx(m, "cppvdisort");
  vx.doc() = unwrap_stars(R"(A low-level polarized VDISORT object.

The calling style mirrors :class:`~pyarts3.arts.cppdisort`, but scalar phase coefficients,
boundary values, sources, and beam intensity are replaced by their polarized
counterparts.  Stokes components are ordered ``[I, Q, U, V]``.

``phase_matrix`` has shape ``[2, NFourier, NLayers, NQuad, NQuad, 4, 4]``.
The leading dimension contains the combined cosine and sine equations.
``b_pos`` and ``b_neg`` have shape ``[2, NFourier, NQuad/2, 4]`` and
``s_poly_coeffs`` has shape ``[NLayers, Ncoeffs, 4]``.  ``b_pos`` is the
radiance beneath the surface, which the surface transmits with ``1 - R``
(``R`` the reflection of ``surface``, specular part included):
for the Planck radiance ``[B, 0, 0, 0]`` in the cosine system of mode 0 this
is Kirchhoff's polarized emission ``B ([1, 0, 0, 0] - R [1, 0, 0, 0])``.
Without a surface ``b_pos`` is the upward radiance itself.  Its coefficients are
the Stokes source function ``B = [B_I, B_Q, B_U, B_V]``; VDISORT applies the
layer factor ``1 - omega`` internally.  An unpolarized source is therefore
supplied as ``[B, 0, 0, 0]``.
)");
  vx.def(
      "__init__",
      [](vdisort::main_data*                                   n,
         const AscendingGrid&                                  tau_arr,
         const Vector&                                         omega_arr,
         const Index                                           NQuad,
         const vdisort::phase_matrix_data&                     phase_matrix,
         Numeric                                               mu0,
         const rtepack::stokvec&                               beam_stokes,
         Numeric                                               phi0,
         const std::optional<Index>                            NFourier_,
         const std::optional<rtepack::stokvec_tensor3>&        b_pos,
         const std::optional<rtepack::stokvec_tensor3>&        b_neg,
         const std::optional<vdisort::surface_reflection>&     surface,
         const std::optional<rtepack::stokvec_matrix>&         s_poly_coeffs,
         const std::optional<vdisort::beam_phase_matrix_data>& beam_phase_matrix,
         const std::optional<Vector>&                          source_coordinate_scale,
         const std::optional<Vector>&                          source_coordinate_offset) {
        const Index NFourier = NFourier_.value_or(phase_matrix.shape()[1]);
        const Index NLayers  = tau_arr.size();

        new (n) vdisort::main_data(NQuad,
                                   NFourier,
                                   tau_arr,
                                   omega_arr,
                                   phase_matrix,
                                   b_pos.value_or(rtepack::stokvec_tensor3(2, NFourier, NQuad / 2)),
                                   b_neg.value_or(rtepack::stokvec_tensor3(2, NFourier, NQuad / 2)),
                                   s_poly_coeffs.value_or(rtepack::stokvec_matrix(NLayers, 0)),
                                   surface.value_or(vdisort::surface_reflection{}),
                                   mu0,
                                   beam_stokes,
                                   phi0,
                                   beam_phase_matrix.value_or(vdisort::beam_phase_matrix_data{}),
                                   source_coordinate_scale.value_or(Vector{}),
                                   source_coordinate_offset.value_or(Vector{}));
      },
      "Run polarized VDISORT with an interface parallel to cppdisort.\n",
      "tau_arr"_a,
      "omega_arr"_a,
      "NQuad"_a,
      "phase_matrix"_a,
      "mu0"_a,
      "beam_stokes"_a,
      "phi0"_a,
      "NFourier"_a.none()                 = py::none(),
      "b_pos"_a.none()                    = py::none(),
      "b_neg"_a.none()                    = py::none(),
      "surface"_a.none()                  = py::none(),
      "s_poly_coeffs"_a.none()            = py::none(),
      "beam_phase_matrix"_a.none()        = py::none(),
      "source_coordinate_scale"_a.none()  = py::none(),
      "source_coordinate_offset"_a.none() = py::none());
  vx.def_prop_ro(
        "tau",
        [](const vdisort::main_data& dis) { return dis.tau(); },
        "Optical depth at the bottom of each layer\n\n.. :class:`~pyarts3.arts.AscendingGrid`")
      .def_prop_ro(
          "omega",
          [](const vdisort::main_data& dis) { return dis.omega(); },
          "Single-scattering albedo of each layer\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_prop_ro(
          "mu",
          [](const vdisort::main_data& dis) { return dis.mu(); },
          "Stream cosines: NQuad / 2 ascending upward streams, then the downward ones (-mu)\n\n"
          ".. :class:`~pyarts3.arts.Vector`")
      .def_prop_ro(
          "weights",
          [](const vdisort::main_data& dis) { return dis.weights(); },
          "Quadrature weights of the upward streams\n\n.. :class:`~pyarts3.arts.Vector`");
  vx.def(
        "u",
        [](vdisort::main_data& dis, const AscendingGrid& tau, const Vector& phi) {
          Tensor4 out(tau.size(), phi.size(), dis.mu().size(), vdisort::stokes_dimension);
          dis.ungridded_u(out, tau, phi);
          return out;
        },
        "tau"_a,
        "phi"_a,
        "Compute the Stokes radiance with shape [tau, phi, stream, stokes]")
      .def(
          "flux",
          [](vdisort::main_data& dis, const AscendingGrid& tau) {
            Matrix out(4, tau.size());
            dis.ungridded_flux(out[0], out[1], out[2], out[3], tau);
            return out;
          },
          "tau"_a,
          "Compute Stokes-I upward, downward-diffuse, downward-direct flux and DFDT")
      .def("has_complex_eigensolutions",
           &vdisort::main_data::has_complex_eigensolutions,
           "tolerance"_a = 1.0e-12,
           "Return whether any retained transport eigenvalue is significantly complex")
      .def(
          "pydisort_u",
          [](vdisort::main_data& dis, Vector tau_, const Vector& phi) {
            std::vector<Index> sorting(tau_.size());
            stdr::iota(sorting, 0);
            stdr::sort(stdv::zip(sorting, tau_), {}, [](const auto& x) { return std::get<1>(x); });

            AscendingGrid tau{std::move(tau_)};
            Tensor4       res(tau.size(), phi.size(), dis.mu().size(), vdisort::stokes_dimension);
            dis.ungridded_u(res, tau, phi);

            Tensor4 out(dis.mu().size(), tau.size(), phi.size(), vdisort::stokes_dimension);
            for (Size i = 0; i < tau.size(); ++i)
              for (Size p = 0; p < phi.size(); ++p)
                for (Size stream = 0; stream < dis.mu().size(); ++stream)
                  out[stream, i, p, joker] = res[sorting[i], p, stream, joker];
            return out;
          },
          "tau"_a,
          "phi"_a,
          "Compute Stokes radiance with shape [stream, tau, phi, stokes]");
  generic_interface(vx);
  // VDISORT PYTHON INTERFACE END

  py::class_<DisortSettings> disort_settings(m, "DisortSettings");
  generic_interface(disort_settings);
  disort_settings.def_rw(
      "quadrature_dimension", &DisortSettings::quadrature_dimension, ".. :class:`~pyarts3.arts.Index`");
  disort_settings.def_rw("legendre_polynomial_dimension",
                         &DisortSettings::legendre_polynomial_dimension,
                         ".. :class:`~pyarts3.arts.Index`");
  disort_settings.def_rw(
      "fourier_mode_dimension", &DisortSettings::fourier_mode_dimension, ".. :class:`~pyarts3.arts.Index`");
  disort_settings.def_rw("freq_grid", &DisortSettings::freq_grid, ".. :class:`~pyarts3.arts.AscendingGrid`");
  disort_settings.def_rw("alt_grid", &DisortSettings::alt_grid, ".. :class:`~pyarts3.arts.DescendingGrid`");
  disort_settings.def_rw(
      "solar_azimuth_angle", &DisortSettings::solar_azimuth_angle, ".. :class:`~pyarts3.arts.Vector`");
  disort_settings.def_rw("solar_zenith_angle", &DisortSettings::solar_zenith_angle, ".. :class:`~pyarts3.arts.Vector`");
  disort_settings.def_rw("solar_source", &DisortSettings::solar_source, ".. :class:`~pyarts3.arts.Vector`");
  disort_settings.def_rw("surface_brdf",
                         &DisortSettings::surface_brdf,
                         "[nfreq, nbrdf, N, N] Fourier modes of the surface BRDF between the streams\n\n"
                         ".. :class:`~pyarts3.arts.Tensor4`");
  disort_settings.def_rw("surface_beam_brdf",
                         &DisortSettings::surface_beam_brdf,
                         "[nfreq, nbrdf, N] the same from the direct beam\n\n.. :class:`~pyarts3.arts.Tensor3`");
  disort_settings.def_rw("surface_beam_mu",
                         &DisortSettings::surface_beam_mu,
                         "[nfreq] beam cosine of surface_beam_brdf, NaN where computed without a beam\n\n"
                         ".. :class:`~pyarts3.arts.Vector`");
  disort_settings.def_rw(
      "optical_thicknesses", &DisortSettings::optical_thicknesses, ".. :class:`~pyarts3.arts.Matrix`");
  disort_settings.def_rw(
      "single_scattering_albedo", &DisortSettings::single_scattering_albedo, ".. :class:`~pyarts3.arts.Matrix`");
  disort_settings.def_rw(
      "fractional_scattering", &DisortSettings::fractional_scattering, ".. :class:`~pyarts3.arts.Matrix`");
  disort_settings.def_rw(
      "delta_m_peak_moments", &DisortSettings::delta_m_peak_moments, ".. :class:`~pyarts3.arts.Tensor3`");
  disort_settings.def_rw("source_polynomial", &DisortSettings::source_polynomial, ".. :class:`~pyarts3.arts.Tensor3`");
  disort_settings.def_rw(
      "legendre_coefficients", &DisortSettings::legendre_coefficients, ".. :class:`~pyarts3.arts.Tensor3`");
  disort_settings.def_rw(
      "upward_boundary_condition", &DisortSettings::upward_boundary_condition, ".. :class:`~pyarts3.arts.Tensor3`");
  disort_settings.def_rw(
      "downward_boundary_condition", &DisortSettings::downward_boundary_condition, ".. :class:`~pyarts3.arts.Tensor3`");

  py::class_<FluxComponentsProfile> df(m, "FluxComponentsProfile");
  generic_interface(df);
  df.def_rw("freq_grid",
            &FluxComponentsProfile::freq_grid,
            "Frequency grid of the fluxes\n\n.. :class:`~pyarts3.arts.AscendingGrid`");
  df.def_rw("alt_grid",
            &FluxComponentsProfile::alt_grid,
            "Altitude grid of the fluxes (level values)\n\n.. :class:`~pyarts3.arts.DescendingGrid`");
  df.def_rw("up", &FluxComponentsProfile::up, "Upwelling flux (layer values)\n\n.. :class:`~pyarts3.arts.Matrix`");
  df.def_rw("down_diffuse",
            &FluxComponentsProfile::down_diffuse,
            "Downward diffuse flux (layer values)\n\n.. :class:`~pyarts3.arts.Matrix`");
  df.def_rw("down_direct",
            &FluxComponentsProfile::down_direct,
            "Downward direct flux (layer values)\n\n.. :class:`~pyarts3.arts.Matrix`");
  df.def_rw(
      "dfdt",
      &FluxComponentsProfile::dfdt,
      "Derivative of net upward flux with respect to downward optical depth [W/(m^2 Hz)], "
      "at each layer's lower boundary (alt_grid[1:]), using that layer's optical properties. "
      "Not a temperature tendency. See pyarts3.recipe.heating_rates.from_disort.\n\n.. :class:`~pyarts3.arts.Matrix`");

  py::class_<DisortRadiance> dr(m, "DisortRadiance");
  generic_interface(dr);
  dr.def_rw("freq_grid",
            &DisortRadiance::freq_grid,
            "Frequency grid of the fluxes\n\n.. :class:`~pyarts3.arts.AscendingGrid`");
  dr.def_rw("alt_grid",
            &DisortRadiance::alt_grid,
            "Altitude grid of the fluxes (level values)\n\n.. :class:`~pyarts3.arts.DescendingGrid`");
  dr.def_rw("zen_grid", &DisortRadiance::zen_grid, "Zenith grid\n\n.. :class:`~pyarts3.arts.ZenGrid`");
  dr.def_rw("azi_grid", &DisortRadiance::azi_grid, "Azimuth grid\n\n.. :class:`~pyarts3.arts.AziGrid`");
  dr.def_rw("data", &DisortRadiance::data, "Radiance field (layer values)\n\n.. :class:`~pyarts3.arts.Tensor4`");

  // VDISORT workspace groups (vdisort_settings.h), the counterparts of the Disort groups above
  py::class_<VDisortSettings> vdisort_settings(m, "VDisortSettings");
  generic_interface(vdisort_settings);
  vdisort_settings.def_rw(
      "quadrature_dimension", &VDisortSettings::quadrature_dimension, ".. :class:`~pyarts3.arts.Index`");
  vdisort_settings.def_rw("legendre_polynomial_dimension",
                          &VDisortSettings::legendre_polynomial_dimension,
                          ".. :class:`~pyarts3.arts.Index`");
  vdisort_settings.def_rw(
      "fourier_mode_dimension", &VDisortSettings::fourier_mode_dimension, ".. :class:`~pyarts3.arts.Index`");
  vdisort_settings.def_rw("freq_grid", &VDisortSettings::freq_grid, ".. :class:`~pyarts3.arts.AscendingGrid`");
  vdisort_settings.def_rw("alt_grid", &VDisortSettings::alt_grid, ".. :class:`~pyarts3.arts.DescendingGrid`");
  vdisort_settings.def_rw(
      "solar_azimuth_angle", &VDisortSettings::solar_azimuth_angle, ".. :class:`~pyarts3.arts.Vector`");
  vdisort_settings.def_rw(
      "solar_zenith_angle", &VDisortSettings::solar_zenith_angle, ".. :class:`~pyarts3.arts.Vector`");
  vdisort_settings.def_rw(
      "solar_source", &VDisortSettings::solar_source, ".. :class:`~pyarts3.arts.StokvecVector`");
  vdisort_settings.def_rw("surface_brdf",
                          &VDisortSettings::surface_brdf,
                          "[nfreq, 2, nbrdf, N, N] combined Fourier modes of the surface BRDF between the streams\n\n"
                          ".. :class:`~pyarts3.arts.MuelmatTensor5`");
  vdisort_settings.def_rw("surface_beam_brdf",
                          &VDisortSettings::surface_beam_brdf,
                          "[nfreq, 2, nbrdf, N] the same from the direct beam\n\n.. :class:`~pyarts3.arts.MuelmatTensor4`");
  vdisort_settings.def_rw("surface_beam_mu",
                          &VDisortSettings::surface_beam_mu,
                          "[nfreq] beam cosine of surface_beam_brdf, NaN where computed without a beam\n\n"
                          ".. :class:`~pyarts3.arts.Vector`");
  vdisort_settings.def_rw("surface_specular",
                          &VDisortSettings::surface_specular,
                          "[nfreq, N] specular part on the streams, or [nfreq, 0]\n\n.. :class:`~pyarts3.arts.MuelmatMatrix`");
  vdisort_settings.def_rw(
      "optical_thicknesses", &VDisortSettings::optical_thicknesses, ".. :class:`~pyarts3.arts.Matrix`");
  vdisort_settings.def_rw(
      "single_scattering_albedo", &VDisortSettings::single_scattering_albedo, ".. :class:`~pyarts3.arts.Matrix`");
  vdisort_settings.def_rw(
      "source_polynomial", &VDisortSettings::source_polynomial, ".. :class:`~pyarts3.arts.StokvecTensor3`");
  vdisort_settings.def_rw("legendre_coefficients",
                          &VDisortSettings::legendre_coefficients,
                          "Legendre moments of [F11, F12, F22, F33, F34, F44]\n\n.. :class:`~pyarts3.arts.Tensor4`");
  vdisort_settings.def_rw("upward_boundary_condition",
                          &VDisortSettings::upward_boundary_condition,
                          ".. :class:`~pyarts3.arts.StokvecTensor4`");
  vdisort_settings.def_rw("downward_boundary_condition",
                          &VDisortSettings::downward_boundary_condition,
                          ".. :class:`~pyarts3.arts.StokvecTensor4`");

  py::class_<VDisortRadiance> vdr(m, "VDisortRadiance");
  generic_interface(vdr);
  vdr.def_rw("freq_grid",
             &VDisortRadiance::freq_grid,
             "Frequency grid of the radiances\n\n.. :class:`~pyarts3.arts.AscendingGrid`");
  vdr.def_rw("alt_grid",
             &VDisortRadiance::alt_grid,
             "Altitude grid of the radiances (level values)\n\n.. :class:`~pyarts3.arts.DescendingGrid`");
  vdr.def_rw("zen_grid", &VDisortRadiance::zen_grid, "Zenith grid\n\n.. :class:`~pyarts3.arts.ZenGrid`");
  vdr.def_rw("azi_grid", &VDisortRadiance::azi_grid, "Azimuth grid\n\n.. :class:`~pyarts3.arts.AziGrid`");
  vdr.def_rw("data",
             &VDisortRadiance::data,
             "Stokes radiance field (layer values)\n\n.. :class:`~pyarts3.arts.StokvecTensor4`");

  // VDISORT PYTHON INTERFACE BEGIN: inputs from ARTS-native data (vdisort_arts.h)
  py::class_<vdisort::fourier_optics> fo(vdisort_nm, "FourierOptics");
  fo.def_ro("extinction", &vdisort::fourier_optics::extinction, "K11 of the particles per metre\n\n.. :class:`float`")
      .def_ro("scattering",
              &vdisort::fourier_optics::scattering,
              "K11 - a1, the extinction minus the absorption, per metre\n\n.. :class:`float`")
      .def_ro("cosine",
              &vdisort::fourier_optics::cosine,
              "[nfourier, mu_out, mu_in] ordinary cosine coefficients C^m of the normalised phase matrix\n\n.. "
              ":class:`~pyarts3.arts.MuelmatTensor3`")
      .def_ro("sine",
              &vdisort::fourier_optics::sine,
              "[nfourier, mu_out, mu_in] ordinary sine coefficients S^m of the normalised phase matrix\n\n.. "
              ":class:`~pyarts3.arts.MuelmatTensor3`");
  fo.doc() = "Normalised phase-matrix Fourier coefficients of scattering species, with their extinction";

  vdisort_nm.def("scattering_optics",
                 &vdisort::scattering_optics,
                 "scattering_species"_a,
                 "atm_point"_a,
                 "frequency"_a,
                 "mu_out"_a,
                 "mu_in"_a,
                 "nfourier"_a,
                 "normalisation_tolerance"_a = 1e-3,
                 R"(The phase-matrix Fourier coefficients of ARTS scattering species at one atmospheric point.

``mu_out`` and ``mu_in`` are signed direction cosines (> 0 upward), e.g.
VDISORT's streams and, for the beam column, ``[-mu0]``.  Returns the
ordinary coefficients ``C^m, S^m = (1 / 2 pi) int P(mu_o, 0; mu_i, phi)
{cos, sin}(m phi) dphi`` (no 2 - delta_m0) of the laboratory-frame phase
matrix ``P = 4 pi Z / sigma``, normalised to 1 over 4 pi, built by vector
geometry in VDISORT's meridional basis (Q = I_v - I_h, U = 2 Re(E_v E_h*)).
The modes, extinction, absorption and phase integral ``sigma = int Z11 dOmega``
all come from the species' Fourier modes at exactly these zenith angles
(``get_bulk_scattering_properties_aro_fourier``).  VDISORT's optical depth and
albedo are scalars, so to ``normalisation_tolerance`` times K11 the extinction,
absorption and ``sigma`` must not depend on the incidence angle nor polarize,
and ``sigma`` must match K11 - a1.  ARTS's own
laboratory-frame phase matrix for the propagation directions (za, aa) is
this one with ``mu = cos(za)`` and ``phi = -aa``.  Pass the results to
:func:`combine_phase_matrices` and :func:`combine_beam_phase_matrices`.
)");

  py::class_<vdisort::lambertian_surface>(vdisort_nm, "LambertianSurface")
      .def(
          "__init__",
          [](vdisort::lambertian_surface* s, Numeric albedo) { new (s) vdisort::lambertian_surface{.albedo = albedo}; },
          "albedo"_a = 0.0)
      .def_rw("albedo", &vdisort::lambertian_surface::albedo, "Albedo A\n\n.. :class:`float`")
      .doc() = "Depolarizing Lambertian surface for main_data_from_path: emission [(1 - A) B, 0, 0, 0]";

  py::class_<vdisort::fresnel_surface>(vdisort_nm, "FresnelSurface")
      .def(
          "__init__",
          [](vdisort::fresnel_surface* s, Complex n) { new (s) vdisort::fresnel_surface{.refractive_index = n}; },
          "refractive_index"_a)
      .def_rw("refractive_index",
              &vdisort::fresnel_surface::refractive_index,
              "Complex refractive index (medium above has index 1)\n\n.. :class:`complex`")
      .doc() =
      "Flat Fresnel surface for main_data_from_path: emission B ([1, 0, 0, 0] - R[:, 0]); reflection only "
      "between quadrature streams";

  const vdisort::path_settings vps{};
  py::class_<vdisort::path_settings>(vdisort_nm, "PathSettings")
      .def(
          "__init__",
          [](vdisort::path_settings* s,
             Index                   nquad,
             Index                   nfourier,
             Numeric                 normalisation_tolerance,
             bool                    thermal,
             Numeric                 beam_flux,
             Numeric                 beam_mu,
             Numeric                 beam_azimuth) {
            new (s) vdisort::path_settings{.nquad                   = nquad,
                                           .nfourier                = nfourier,
                                           .normalisation_tolerance = normalisation_tolerance,
                                           .thermal                 = thermal,
                                           .beam_flux               = beam_flux,
                                           .beam_mu                 = beam_mu,
                                           .beam_azimuth            = beam_azimuth};
          },
          "nquad"_a                   = vps.nquad,
          "nfourier"_a                = vps.nfourier,
          "normalisation_tolerance"_a = vps.normalisation_tolerance,
          "thermal"_a                 = vps.thermal,
          "beam_flux"_a               = vps.beam_flux,
          "beam_mu"_a                 = vps.beam_mu,
          "beam_azimuth"_a            = vps.beam_azimuth)
      .def_rw("nquad", &vdisort::path_settings::nquad, "Number of streams, even\n\n.. :class:`int`")
      .def_rw("nfourier", &vdisort::path_settings::nfourier, "Number of Fourier modes\n\n.. :class:`int`")
      .def_rw("normalisation_tolerance",
              &vdisort::path_settings::normalisation_tolerance,
              "Allowed mismatch of phase-function integral and scattering coefficient, relative to the "
              "extinction; infinity checks nothing\n\n.. :class:`float`")
      .def_rw("thermal",
              &vdisort::path_settings::thermal,
              "Thermal emission of the layers and the surface\n\n.. :class:`bool`")
      .def_rw("beam_flux",
              &vdisort::path_settings::beam_flux,
              "Beam flux on the horizontal at the top [W m-2 Hz-1], 0 for none\n\n.. :class:`float`")
      .def_rw("beam_mu", &vdisort::path_settings::beam_mu, "Cosine of the beam zenith angle\n\n.. :class:`float`")
      .def_rw("beam_azimuth",
              &vdisort::path_settings::beam_azimuth,
              "VDISORT's beam azimuth phi0 [rad]\n\n.. :class:`float`")
      .doc() = "Solver settings of main_data_from_path";

  vdisort_nm.def("main_data_from_path",
                 &vdisort::main_data_from_path,
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
                 R"(A solved VDISORT problem (:class:`~pyarts3.arts.cppvdisort`) from an ARTS propagation path.

The path conventions are those of :func:`pyarts3.arts.rt4.problem_from_path`
(one entry per level, top first, unpolarized gas propagation matrix).  A
layer has the mean extinction and scattering of its two levels'
:func:`scattering_optics` on VDISORT's streams and their scattering-weighted
mean Fourier coefficients; ``tau`` is the cumulative (gas + particle)
optical depth and ``omega`` the scattering over the total extinction.  With
``settings.thermal`` the Planck function at the level temperatures is linear
in optical depth within each layer and the surface emits; the sky is a
blackbody at ``sky_temperature``.  A beam has the Stokes irradiance
``[beam_flux / beam_mu, 0, 0, 0]`` normal to it.
)");
  // VDISORT PYTHON INTERFACE END
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize disort\n{}", e.what()));
}
}  // namespace Python
