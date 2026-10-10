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

#include <concepts>
#include <optional>
#include <ranges>

#include "configtypes.h"
#include "debug.h"
#include "hpy_arts.h"
#include "operators.h"
#include "sorting.h"

NB_MAKE_OPAQUE(std::vector<vdisort::BDRF>);

namespace Python {
using DisortBDRFOperator = CustomOperator<Matrix, const Vector&, const Vector&>;
using bdrf_func          = DisortBDRFOperator::func_t;

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

  py::class_<DisortBDRFOperator> bdrfop(m, "DisortBDRFOperator");
  bdrfop.doc() = "A BDRF operator for DISORT";
  bdrfop
      .def("__init__",
           [](DisortBDRFOperator* op, DisortBDRFOperator::func_t f) {
             new (op) DisortBDRFOperator([f](const Vector& x, const Vector& y) {
               py::gil_scoped_acquire gil{};
               return f(x, y);
             });
           })
      .def("__call__", [](DisortBDRFOperator& f, const Vector& x, const Vector& y) { return f.f(x, y); }, "x"_a, "y"_a);
  generic_interface(bdrfop);  // FIXME OLE
  py::implicitly_convertible<DisortBDRFOperator::func_t, DisortBDRFOperator>();

  py::class_<DisortBDRF> disbdrf(m, "DisortBDRF");
  disbdrf
      .def(
          "__init__",
          [](DisortBDRF* b, const DisortBDRFOperator& f) {
            new (b)
                DisortBDRF(DisortBDRF::func_t{[f](MatrixView mat, const ConstVectorView& a, const ConstVectorView& b) {
                  const Matrix out = f(Vector{a}, Vector{b});
                  if (out.shape() != mat.shape()) {
                    throw std::runtime_error(
                        std::format("BDRF function returned wrong shape\n{:B,} vs {:B,}", out.shape(), mat.shape()));
                  }
                  mat = out;
                }});
          },
          py::keep_alive<0, 1>())
      .def("__call__", [](const DisortBDRF& bdrf, const Vector& a, const Vector& b) {
        Matrix out(a.size(), b.size());
        bdrf(out, a, b);
        return out;
      });
  generic_interface(disbdrf);
  py::implicitly_convertible<bdrf_func, DisortBDRF>();

  py::class_<MatrixOfDisortBDRF> mat_disbdrf(m, "MatrixOfDisortBDRF");
  generic_interface(mat_disbdrf);

  auto vecs  = py::bind_vector<std::vector<DisortBDRF>, py::rv_policy::reference_internal>(disort_nm, "ArrayOfBDRF");
  vecs.doc() = "An array of BDRF functions";
  generic_interface(vecs);

  disort_nm.def("lambertian_fourier_modes",
                &disort::brdf::lambertian_fourier_modes,
                "albedo"_a,
                "number_of_modes"_a,
                "Construct exact scalar Lambertian Fourier modes");
  disort_nm.def("combine_fourier_modes",
                &disort::brdf::combine_fourier_modes,
                "first"_a,
                "first_weight"_a,
                "second"_a,
                "second_weight"_a,
                "Form a weighted sum of two scalar Fourier-mode surface models");
  disort_nm.def("cox_munk_lambertian_fourier_modes",
                &disort::brdf::cox_munk_lambertian_fourier_modes,
                "cox_munk_fraction"_a,
                "lambertian_albedo"_a,
                "wind_speed"_a,
                "refractive_index"_a,
                "shadowing"_a,
                "number_of_modes"_a,
                "azimuth_quadrature_points"_a = 100,
                "Construct a scalar Cox-Munk/Lambertian Fourier-mode mixture");

  // VDISORT PYTHON INTERFACE BEGIN: a Fourier BRDF mode has cosine and sine
  // Mueller-matrix callbacks.  Each callback returns [4*n_out, 4*n_in].
  const auto polarized_bdrf_callback = [](const DisortBDRFOperator& f) {
    return vdisort::BDRF::func_t{
        [f](rtepack::muelmat_matrix_view mat, const ConstVectorView& mu_out, const ConstVectorView& mu_in) {
          const Matrix     out = f(Vector{mu_out}, Vector{mu_in});
          const std::array expected{4 * mat.nrows(), 4 * mat.ncols()};
          if (out.shape() != expected) {
            throw std::runtime_error(
                std::format("Polarized BDRF function returned wrong shape\n{:B,} vs {:B,}", out.shape(), expected));
          }
          for (Index i = 0; i < mat.nrows(); ++i)
            for (Index j = 0; j < mat.ncols(); ++j)
              for (Index so = 0; so < vdisort::stokes_dimension; ++so)
                for (Index si = 0; si < vdisort::stokes_dimension; ++si)
                  mat[i, j][so, si] = out[4 * i + so, 4 * j + si];
        }};
  };

  py::class_<vdisort::BDRF> vdisbdrf(m, "VDisortBDRF");
  vdisbdrf.doc() = unwrap_stars(R"(A polarized VDISORT BDRF Fourier mode.

The cosine and optional sine callables receive ``(mu_out, mu_in)`` and return
a matrix of shape ``(4 * len(mu_out), 4 * len(mu_in))``.  Each 4-by-4 block is
the Mueller reflection matrix for one outgoing/incident stream pair.
)");
  vdisbdrf
      .def(
          "__init__",
          [polarized_bdrf_callback](vdisort::BDRF* b, const DisortBDRFOperator& cosine) {
            new (b) vdisort::BDRF{
                .cosine      = polarized_bdrf_callback(cosine),
                .sine        = vdisort::BDRF::func_t{[](rtepack::muelmat_matrix_view mat,
                                                        const ConstVectorView&,
                                                        const ConstVectorView&) { mat = rtepack::muelmat{0.0}; }},
                .beam_cosine = {},
                .beam_sine   = {}};
          },
          "cosine"_a,
          py::keep_alive<0, 1>())
      .def(
          "__init__",
          [polarized_bdrf_callback](
              vdisort::BDRF* b, const DisortBDRFOperator& cosine, const DisortBDRFOperator& sine) {
            new (b) vdisort::BDRF{.cosine      = polarized_bdrf_callback(cosine),
                                  .sine        = polarized_bdrf_callback(sine),
                                  .beam_cosine = {},
                                  .beam_sine   = {}};
          },
          "cosine"_a,
          "sine"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>())
      .def(
          "__call__",
          [](const vdisort::BDRF& bdrf, const Index alpha, const Vector& mu_out, const Vector& mu_in) {
            rtepack::muelmat_matrix blocks(mu_out.size(), mu_in.size(), rtepack::muelmat{0.0});
            bdrf(alpha, blocks, mu_out, mu_in);
            Matrix out(vdisort::stokes_dimension * mu_out.size(), vdisort::stokes_dimension * mu_in.size());
            for (Index i = 0; i < blocks.nrows(); ++i)
              for (Index j = 0; j < blocks.ncols(); ++j)
                for (Index so = 0; so < vdisort::stokes_dimension; ++so)
                  for (Index si = 0; si < vdisort::stokes_dimension; ++si)
                    out[4 * i + so, 4 * j + si] = blocks[i, j][so, si];
            return out;
          },
          "alpha"_a,
          "mu_out"_a,
          "mu_in"_a);
  generic_interface(vdisbdrf);
  py::implicitly_convertible<bdrf_func, vdisort::BDRF>();

  auto vvecs =
      py::bind_vector<std::vector<vdisort::BDRF>, py::rv_policy::reference_internal>(vdisort_nm, "ArrayOfBDRF");
  vvecs.doc() = "An array of polarized BDRF Fourier modes";
  generic_interface(vvecs);

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
  vdisort_nm.def("cox_munk_fourier_modes",
                 &vdisort::brdf::cox_munk_fourier_modes,
                 "wind_speed"_a,
                 "refractive_index"_a,
                 "shadowing"_a,
                 "number_of_modes"_a,
                 "azimuth_quadrature_points"_a = 100,
                 "Construct VDISORT-ready combined Fourier modes for a polarized Cox-Munk ocean");
  vdisort_nm.def(
      "fresnel_reflection",
      [](const Numeric incident_mu, const Complex refractive_index) {
        return vdisort::brdf::Fresnel{refractive_index}(incident_mu);
      },
      "incident_mu"_a,
      "refractive_index"_a = Complex{1.5, 0.0},
      "Evaluate the polarized Fresnel reflection matrix of a flat dielectric interface");
  vdisort_nm.def("fresnel_fourier_modes",
                 &vdisort::brdf::fresnel_fourier_modes,
                 "refractive_index"_a,
                 "number_of_modes"_a,
                 "Construct quadrature-normalized VDISORT Fourier modes for an ideal Fresnel surface");
  vdisort_nm.def("lambertian_fourier_modes",
                 &vdisort::brdf::lambertian_fourier_modes,
                 "albedo"_a,
                 "number_of_modes"_a,
                 "Construct exact fully depolarizing Lambertian VDISORT Fourier modes");
  vdisort_nm.def("hapke_fourier_modes",
                 &vdisort::brdf::hapke_fourier_modes,
                 "opposition_amplitude"_a,
                 "opposition_width"_a,
                 "single_scattering_albedo"_a,
                 "number_of_modes"_a,
                 "azimuth_quadrature_points"_a = 100,
                 "Construct fully depolarizing Hapke VDISORT Fourier modes");
  vdisort_nm.def("rpv_fourier_modes",
                 &vdisort::brdf::rpv_fourier_modes,
                 "rho0"_a,
                 "kappa"_a,
                 "asymmetry"_a,
                 "hotspot"_a,
                 "number_of_modes"_a,
                 "azimuth_quadrature_points"_a = 100,
                 "Construct fully depolarizing RPV VDISORT Fourier modes");
  vdisort_nm.def("ross_li_fourier_modes",
                 &vdisort::brdf::ross_li_fourier_modes,
                 "isotropic"_a,
                 "volumetric"_a,
                 "geometric"_a,
                 "hotspot_angle"_a,
                 "number_of_modes"_a,
                 "azimuth_quadrature_points"_a = 100,
                 "Construct fully depolarizing Ross-Li VDISORT Fourier modes");
  vdisort_nm.def("combine_fourier_modes",
                 &vdisort::brdf::combine_fourier_modes,
                 "first"_a,
                 "first_weight"_a,
                 "second"_a,
                 "second_weight"_a,
                 "Form a weighted sum of two polarized Fourier-mode surface models");
  vdisort_nm.def("fresnel_lambertian_fourier_modes",
                 &vdisort::brdf::fresnel_lambertian_fourier_modes,
                 "fresnel_fraction"_a,
                 "lambertian_albedo"_a,
                 "refractive_index"_a,
                 "number_of_modes"_a,
                 "Construct an ideal Fresnel/depolarizing-Lambertian VDISORT mixture");
  vdisort_nm.def("cox_munk_lambertian_fourier_modes",
                 &vdisort::brdf::cox_munk_lambertian_fourier_modes,
                 "cox_munk_fraction"_a,
                 "lambertian_albedo"_a,
                 "wind_speed"_a,
                 "refractive_index"_a,
                 "shadowing"_a,
                 "number_of_modes"_a,
                 "azimuth_quadrature_points"_a = 100,
                 "Construct a rough Fresnel Cox-Munk/depolarizing-Lambertian VDISORT mixture");
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
         const std::vector<DisortBDRF>& bdrf,
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
                                  bdrf,
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
      "BDRF_Fourier_modes"_a   = std::vector<DisortBDRF>{},
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
          [](disort::main_data& dis, const Vector& mu, const AscendingGrid& tau, const Vector& phi) {
            Tensor3             out(mu.size(), tau.size(), phi.size());
            disort::user_u_data data;
            for (Size t = 0; t < tau.size(); ++t)
              for (Size p = 0; p < phi.size(); ++p) {
                dis.u_user(data, tau[t], phi[p], mu);
                out[joker, t, p] = data.intensities;
              }
            return out;
          },
          "mu"_a,
          "tau"_a,
          "phi"_a,
          "Compute intensity at user polar-angle cosines using DISORT source reconstruction and formal ray integration")
      .def(
          "u_user_corr",
          [](disort::main_data& dis, const Vector& mu, const AscendingGrid& tau, const Vector& phi) {
            Tensor3             out(mu.size(), tau.size(), phi.size());
            disort::user_u_data data;
            disort::tms_data    tms;
            Vector              ims;
            for (Size t = 0; t < tau.size(); ++t)
              for (Size p = 0; p < phi.size(); ++p) {
                dis.u_user_corr(data, ims, tms, tau[t], phi[p], mu);
                out[joker, t, p] = data.intensities;
              }
            return out;
          },
          "mu"_a,
          "tau"_a,
          "phi"_a,
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
``s_poly_coeffs`` has shape ``[NLayers, Ncoeffs, 4]``.  Its coefficients are
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
         const std::vector<vdisort::BDRF>&                     bdrf,
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
                                   bdrf,
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
      "BDRF_Fourier_modes"_a              = std::vector<vdisort::BDRF>{},
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
  disort_settings.def_rw("bidirectional_reflectance_distribution_functions",
                         &DisortSettings::bidirectional_reflectance_distribution_functions,
                         ".. :class:`~pyarts3.arts.MatrixOfDisortBDRF`");
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

  py::class_<DisortFlux> df(m, "DisortFlux");
  generic_interface(df);
  df.def_rw(
      "freq_grid", &DisortFlux::freq_grid, "Frequency grid of the fluxes\n\n.. :class:`~pyarts3.arts.AscendingGrid`");
  df.def_rw("alt_grid",
            &DisortFlux::alt_grid,
            "Altitude grid of the fluxes (level values)\n\n.. :class:`~pyarts3.arts.DescendingGrid`");
  df.def_rw("up", &DisortFlux::up, "Upwelling flux (layer values)\n\n.. :class:`~pyarts3.arts.Matrix`");
  df.def_rw("down_diffuse",
            &DisortFlux::down_diffuse,
            "Downward diffuse flux (layer values)\n\n.. :class:`~pyarts3.arts.Matrix`");
  df.def_rw("down_direct",
            &DisortFlux::down_direct,
            "Downward direct flux (layer values)\n\n.. :class:`~pyarts3.arts.Matrix`");
  df.def_rw(
      "dfdt",
      &DisortFlux::dfdt,
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
