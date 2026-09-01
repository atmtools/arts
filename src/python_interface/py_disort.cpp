#include <disort-brdf.h>
#include <disort.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/vector.h>
#include <pydocs.h>
#include <python_interface.h>
#include <vdisort-brdf.h>
#include <vdisort.h>

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
  delta_m_transport.def_rw("tau", &vdisort::delta_m_transport_data::tau)
      .def_rw("omega", &vdisort::delta_m_transport_data::omega)
      .def_rw("phase_matrix", &vdisort::delta_m_transport_data::phase_matrix)
      .def_rw("beam_phase_matrix", &vdisort::delta_m_transport_data::beam_phase_matrix)
      .def_rw("source_coordinate_scale", &vdisort::delta_m_transport_data::source_coordinate_scale)
      .def_rw("source_coordinate_offset", &vdisort::delta_m_transport_data::source_coordinate_offset);
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
      .def_rw(
          "iterations", &disort::coupling_result::iterations, "Number of fixed-point iterations\n\n.. :class:`Index`")
      .def_rw("max_relative_change",
              &disort::coupling_result::max_relative_change,
              "Maximum relative interface update in the last iteration\n\n.. :class:`Numeric`")
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

The calling style mirrors :class:`cppdisort`, but scalar phase coefficients,
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
          "u_user",
          [](vdisort::main_data&                                   dis,
             const AscendingGrid&                                  tau,
             const Vector&                                         phi,
             const Vector&                                         mu,
             const vdisort::phase_matrix_data&                     phase_matrix,
             const std::optional<vdisort::beam_phase_matrix_data>& beam_phase_matrix) {
            rtepack::stokvec_tensor3 out(tau.size(), phi.size(), mu.size());
            dis.ungridded_u_user(
                out, tau, phi, mu, phase_matrix, beam_phase_matrix.value_or(vdisort::beam_phase_matrix_data{}));
            return out;
          },
          "tau"_a,
          "phi"_a,
          "mu"_a,
          "phase_matrix"_a,
          "beam_phase_matrix"_a.none() = py::none(),
          R"(Compute Stokes radiance at arbitrary nonzero polar-angle cosines.

The phase arrays contain the combined cosine/sine Mueller coefficients sampled
at the requested outgoing directions.  Their numerical shapes are
``[2, NFourier, NLayers, NUser, NQuad, 4, 4]`` and
``[2, NFourier, NLayers, NUser, 4, 4]`` for the direct beam.)")
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
  disort_settings.def_rw("quadrature_dimension", &DisortSettings::quadrature_dimension, ".. :class:`Index`");
  disort_settings.def_rw(
      "legendre_polynomial_dimension", &DisortSettings::legendre_polynomial_dimension, ".. :class:`Index`");
  disort_settings.def_rw("fourier_mode_dimension", &DisortSettings::fourier_mode_dimension, ".. :class:`Index`");
  disort_settings.def_rw("freq_grid", &DisortSettings::freq_grid, ".. :class:`AscendingGrid`");
  disort_settings.def_rw("alt_grid", &DisortSettings::alt_grid, ".. :class:`DescendingGrid`");
  disort_settings.def_rw("solar_azimuth_angle", &DisortSettings::solar_azimuth_angle, ".. :class:`Vector`");
  disort_settings.def_rw("solar_zenith_angle", &DisortSettings::solar_zenith_angle, ".. :class:`Vector`");
  disort_settings.def_rw("solar_source", &DisortSettings::solar_source, ".. :class:`Vector`");
  disort_settings.def_rw("bidirectional_reflectance_distribution_functions",
                         &DisortSettings::bidirectional_reflectance_distribution_functions,
                         ".. :class:`MatrixOfDisortBDRF`");
  disort_settings.def_rw("optical_thicknesses", &DisortSettings::optical_thicknesses, ".. :class:`Matrix`");
  disort_settings.def_rw("single_scattering_albedo", &DisortSettings::single_scattering_albedo, ".. :class:`Matrix`");
  disort_settings.def_rw("fractional_scattering", &DisortSettings::fractional_scattering, ".. :class:`Matrix`");
  disort_settings.def_rw("delta_m_peak_moments", &DisortSettings::delta_m_peak_moments, ".. :class:`Tensor3`");
  disort_settings.def_rw("source_polynomial", &DisortSettings::source_polynomial, ".. :class:`Tensor3`");
  disort_settings.def_rw("legendre_coefficients", &DisortSettings::legendre_coefficients, ".. :class:`Tensor3`");
  disort_settings.def_rw(
      "upward_boundary_condition", &DisortSettings::upward_boundary_condition, ".. :class:`Tensor3`");
  disort_settings.def_rw(
      "downward_boundary_condition", &DisortSettings::downward_boundary_condition, ".. :class:`Tensor3`");

  py::class_<DisortFlux> df(m, "DisortFlux");
  generic_interface(df);
  df.def_rw("freq_grid", &DisortFlux::freq_grid, "Frequency grid of the fluxes\n\n.. :class:`AscendingGrid`");
  df.def_rw(
      "alt_grid", &DisortFlux::alt_grid, "Altitude grid of the fluxes (level values)\n\n.. :class:`DescendingGrid`");
  df.def_rw("up", &DisortFlux::up, "Upwelling flux (layer values)\n\n.. :class:`Matrix`");
  df.def_rw("down_diffuse", &DisortFlux::down_diffuse, "Downward diffuse flux (layer values)\n\n.. :class:`Matrix`");
  df.def_rw("down_direct", &DisortFlux::down_direct, "Downward direct flux (layer values)\n\n.. :class:`Matrix`");
  df.def_rw("dfdt", &DisortFlux::dfdt, "Flux divergence (layer values)\n\n.. :class:`Matrix`");

  py::class_<DisortRadiance> dr(m, "DisortRadiance");
  generic_interface(dr);
  dr.def_rw("freq_grid", &DisortRadiance::freq_grid, "Frequency grid of the fluxes\n\n.. :class:`AscendingGrid`");
  dr.def_rw("alt_grid",
            &DisortRadiance::alt_grid,
            "Altitude grid of the fluxes (level values)\n\n.. :class:`DescendingGrid`");
  dr.def_rw("zen_grid", &DisortRadiance::zen_grid, "Zenith grid\n\n.. :class:`ZenGrid`");
  dr.def_rw("azi_grid", &DisortRadiance::azi_grid, "Azimuth grid\n\n.. :class:`AziGrid`");
  dr.def_rw("data", &DisortRadiance::data, "Radiance field (layer values)\n\n.. :class:`Tensor4`");
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize disort\n{}", e.what()));
}
}  // namespace Python
