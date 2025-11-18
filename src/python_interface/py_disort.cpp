#include <disort.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/vector.h>
#include <pydocs.h>
#include <python_interface.h>

#include <concepts>
#include <optional>

#include "configtypes.h"
#include "debug.h"
#include "hpy_arts.h"
#include "operators.h"
#include "sorting.h"

namespace Python {
using DisortBDRFOperator = CustomOperator<Matrix, const Vector&, const Vector&>;
using bdrf_func          = DisortBDRFOperator::func_t;

void py_disort(py::module_& m) try {
  auto disort_nm  = m.def_submodule("disort");
  disort_nm.doc() = "DISORT solver internal types";

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
      .def(
          "__call__",
          [](DisortBDRFOperator& f, const Vector& x, const Vector& y) {
            return f.f(x, y);
          },
          "x"_a,
          "y"_a);
  generic_interface(bdrfop);  // FIXME OLE
  py::implicitly_convertible<DisortBDRFOperator::func_t, DisortBDRFOperator>();

  py::class_<DisortBDRF> disbdrf(m, "DisortBDRF");
  disbdrf
      .def(
          "__init__",
          [](DisortBDRF* b, const DisortBDRFOperator& f) {
            new (b)
                DisortBDRF(DisortBDRF::func_t{[f](MatrixView mat,
                                                  const ConstVectorView& a,
                                                  const ConstVectorView& b) {
                  const Matrix out = f(Vector{a}, Vector{b});
                  if (out.shape() != mat.shape()) {
                    throw std::runtime_error(std::format(
                        "BDRF function returned wrong shape\n{:B,} vs {:B,}",
                        out.shape(),
                        mat.shape()));
                  }
                  mat = out;
                }});
          },
          py::keep_alive<0, 1>())
      .def("__call__",
           [](const DisortBDRF& bdrf, const Vector& a, const Vector& b) {
             Matrix out(a.size(), b.size());
             bdrf(out, a, b);
             return out;
           });
  generic_interface(disbdrf);
  py::implicitly_convertible<bdrf_func, DisortBDRF>();

  py::class_<MatrixOfDisortBDRF> mat_disbdrf(m, "MatrixOfDisortBDRF");
  generic_interface(mat_disbdrf);

  auto vecs  = py::bind_vector<std::vector<DisortBDRF>,
                               py::rv_policy::reference_internal>(disort_nm,
                                                                 "ArrayOfBDRF");
  vecs.doc() = "An array of BDRF functions";
  generic_interface(vecs);

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
      [](disort::main_data* n,
         const AscendingGrid& tau_arr,
         const Vector& omega_arr,
         const Index NQuad,
         const Matrix& Leg_coeffs_all,
         Numeric mu0,
         Numeric I0,
         Numeric phi0,
         const std::optional<Index> NLeg_,
         const std::optional<Index> NFourier_,
         const std::optional<Matrix>& b_pos,
         const std::optional<Matrix>& b_neg,
         const std::optional<Vector>& f_arr,
         const std::vector<DisortBDRF>& bdrf,
         const std::optional<Matrix>& s_poly_coeffs) {
        const Index NFourier = NFourier_.value_or(NQuad);
        const Index NLeg     = NLeg_.value_or(NQuad);
        const Index NLayers  = tau_arr.size();

        new (n)
            disort::main_data(NQuad,
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
                              phi0);
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
      "s_poly_coeffs"_a.none() = py::none());
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
          "flux",
          [](disort::main_data& dis, const AscendingGrid& tau) {
            Matrix out(3, tau.size());
            dis.ungridded_flux(out[0], out[1], out[2], tau);
            return out;
          },
          "tau"_a,
          "Compute the flux")
      .def(
          "pydisort_u",
          [](disort::main_data& dis, Vector tau_, const Vector& phi) {
            std::vector<Index> sorting(tau_.size());
            std::iota(sorting.begin(), sorting.end(), 0);
            bubble_sort_by(
                [&tau_](auto il, auto jl) { return tau_[il] > tau_[jl]; },
                sorting,
                tau_);

            AscendingGrid tau{std::move(tau_)};
            Tensor3 res(tau.size(), phi.size(), dis.mu().size());
            dis.ungridded_u(res, tau, phi);

            Tensor3 out(dis.mu().size(), tau.size(), phi.size());
            for (Size i = 0; i < tau.size(); i++) {
              out[joker, i, joker] = transpose(res[sorting[i]]);
            }
            return out;
          },
          "tau"_a,
          "phi"_a,
          "Compute the intensity")
      .def(
          "pydisort_flux_up",
          [](disort::main_data& dis, Vector tau_) {
            std::vector<Index> sorting(tau_.size());
            std::iota(sorting.begin(), sorting.end(), 0);
            bubble_sort_by(
                [&tau_](auto il, auto jl) { return tau_[il] > tau_[jl]; },
                sorting,
                tau_);

            AscendingGrid tau{std::move(tau_)};
            Matrix res(3, tau.size());
            dis.ungridded_flux(res[0], res[1], res[2], tau);

            Vector out(tau.size());
            for (Size i = 0; i < tau.size(); i++) {
              out[i] = res[0, sorting[i]];
            }
            return out;
          },
          "tau"_a,
          "Compute the upward flux")
      .def(
          "pydisort_flux_down",
          [](disort::main_data& dis, Vector tau_) {
            std::vector<Index> sorting(tau_.size());
            std::iota(sorting.begin(), sorting.end(), 0);
            bubble_sort_by(
                [&tau_](auto il, auto jl) { return tau_[il] > tau_[jl]; },
                sorting,
                tau_);

            AscendingGrid tau{std::move(tau_)};
            Matrix res(3, tau.size());
            dis.ungridded_flux(res[0], res[1], res[2], tau);

            ArrayOfVector out(2, Vector(tau.size()));
            for (Size i = 0; i < tau.size(); i++) {
              out[0][i] = res[1, sorting[i]];
              out[1][i] = res[2, sorting[i]];
            }
            return out;
          },
          "tau"_a,
          "Compute the downward flux");
  generic_interface(x);

  py::class_<DisortSettings> disort_settings(m, "DisortSettings");
  generic_interface(disort_settings);
  disort_settings.def_rw("quadrature_dimension",
                         &DisortSettings::quadrature_dimension,
                         ".. :class:`Index`");
  disort_settings.def_rw("legendre_polynomial_dimension",
                         &DisortSettings::legendre_polynomial_dimension,
                         ".. :class:`Index`");
  disort_settings.def_rw("fourier_mode_dimension",
                         &DisortSettings::fourier_mode_dimension,
                         ".. :class:`Index`");
  disort_settings.def_rw(
      "freq_grid", &DisortSettings::freq_grid, ".. :class:`AscendingGrid`");
  disort_settings.def_rw(
      "alt_grid", &DisortSettings::alt_grid, ".. :class:`DescendingGrid`");
  disort_settings.def_rw("solar_azimuth_angle",
                         &DisortSettings::solar_azimuth_angle,
                         ".. :class:`Vector`");
  disort_settings.def_rw("solar_zenith_angle",
                         &DisortSettings::solar_zenith_angle,
                         ".. :class:`Vector`");
  disort_settings.def_rw(
      "solar_source", &DisortSettings::solar_source, ".. :class:`Vector`");
  disort_settings.def_rw(
      "bidirectional_reflectance_distribution_functions",
      &DisortSettings::bidirectional_reflectance_distribution_functions,
      ".. :class:`MatrixOfDisortBDRF`");
  disort_settings.def_rw("optical_thicknesses",
                         &DisortSettings::optical_thicknesses,
                         ".. :class:`Matrix`");
  disort_settings.def_rw("single_scattering_albedo",
                         &DisortSettings::single_scattering_albedo,
                         ".. :class:`Matrix`");
  disort_settings.def_rw("fractional_scattering",
                         &DisortSettings::fractional_scattering,
                         ".. :class:`Matrix`");
  disort_settings.def_rw("source_polynomial",
                         &DisortSettings::source_polynomial,
                         ".. :class:`Tensor3`");
  disort_settings.def_rw("legendre_coefficients",
                         &DisortSettings::legendre_coefficients,
                         ".. :class:`Tensor3`");
  disort_settings.def_rw("positive_boundary_condition",
                         &DisortSettings::positive_boundary_condition,
                         ".. :class:`Tensor3`");
  disort_settings.def_rw("negative_boundary_condition",
                         &DisortSettings::negative_boundary_condition,
                         ".. :class:`Tensor3`");

  py::class_<DisortFlux> df(m, "DisortFlux");
  generic_interface(df);
  df.def_rw("freq_grid",
            &DisortFlux::freq_grid,
            "Frequency grid of the fluxes\n\n.. :class:`AscendingGrid`");
  df.def_rw(
      "alt_grid",
      &DisortFlux::alt_grid,
      "Altitude grid of the fluxes (level values)\n\n.. :class:`DescendingGrid`");
  df.def_rw("up",
            &DisortFlux::up,
            "Upwelling flux (layer values)\n\n.. :class:`Matrix`");
  df.def_rw("down_diffuse",
            &DisortFlux::down_diffuse,
            "Downward diffuse flux (layer values)\n\n.. :class:`Matrix`");
  df.def_rw("down_direct",
            &DisortFlux::down_direct,
            "Downward direct flux (layer values)\n\n.. :class:`Matrix`");

  py::class_<DisortRadiance> dr(m, "DisortRadiance");
  generic_interface(dr);
  dr.def_rw("freq_grid",
            &DisortRadiance::freq_grid,
            "Frequency grid of the fluxes\n\n.. :class:`AscendingGrid`");
  dr.def_rw(
      "alt_grid",
      &DisortRadiance::alt_grid,
      "Altitude grid of the fluxes (level values)\n\n.. :class:`DescendingGrid`");
  dr.def_rw("zen_grid",
            &DisortRadiance::zen_grid,
            "Zenith grid\n\n.. :class:`ZenGrid`");
  dr.def_rw("azi_grid",
            &DisortRadiance::azi_grid,
            "Azimuth grid\n\n.. :class:`AziGrid`");
  dr.def_rw("data",
            &DisortRadiance::data,
            "Radiance field (layer values)\n\n.. :class:`Matrix`");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize disort\n{}", e.what()));
}
}  // namespace Python
