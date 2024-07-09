#include <disort.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>

#include <memory>
#include <optional>

#include "configtypes.h"
#include "debug.h"
#include "matpack_iter.h"
#include "sorted_grid.h"
#include "sorting.h"

NB_MAKE_OPAQUE(std::vector<disort::BDRF>);

namespace Python {
using bdrf_func = std::function<Matrix(const Vector&, const Vector&)>;

void py_disort(py::module_& m) try {
  auto disort_nm = m.def_submodule("disort");

  py::class_<disort::BDRF>(disort_nm, "bdrf")
      .def(py::init<>())
      .def(
          "__init__",
          [](disort::BDRF* b, const bdrf_func& f) {
            new (b) disort::BDRF([f](ExhaustiveMatrixView mat,
                                     const ExhaustiveConstVectorView& a,
                                     const ExhaustiveConstVectorView& b) {
              const Matrix out = f(Vector{a}, Vector{b});
              if (out.shape() != mat.shape()) {
                throw std::runtime_error(
                    var_string("BDRF function returned wrong shape\n",
                               out.shape(),
                               " vs ",
                               mat.shape()));
              }
              mat = out;
            });
          },
          py::keep_alive<0, 1>())
      .def("__call__",
           [](const disort::BDRF& bdrf, const Vector& a, const Vector& b) {
             Matrix out(a.size(), b.size());
             bdrf(out, a, b);
             return out;
           });
  py::implicitly_convertible<bdrf_func, disort::BDRF>();

  py::bind_vector<std::vector<disort::BDRF>, py::rv_policy::reference_internal>(
      disort_nm, "ArrayOfBDRF");

  py::class_<disort::main_data> x(m, "cppdisort");
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
         const std::vector<disort::BDRF>& bdrf,
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
      "BDRF_Fourier_modes"_a   = std::vector<disort::BDRF>{},
      "s_poly_coeffs"_a.none() = py::none());
  x.def(
       "u",
       [](disort::main_data& dis, const AscendingGrid& tau, const Vector& phi) {
         Tensor3 out(tau.size(), phi.size(), dis.mu().size());
         dis.ungridded_u(out, tau, phi);
         return out;
       },
       "tau"_a,
       "phi"_a)
      .def(
          "flux",
          [](disort::main_data& dis, const AscendingGrid& tau) {
            Matrix out(3, tau.size());
            dis.ungridded_flux(out[0], out[1], out[2], tau);
            return out;
          },
          "tau"_a)
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
            for (Index i = 0; i < tau.size(); i++) {
              out(joker, i, joker) = transpose(res[sorting[i]]);
            }
            return out;
          },
          "tau"_a,
          "phi"_a)
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
            for (Index i = 0; i < tau.size(); i++) {
              out[i] = res(0, sorting[i]);
            }
            return out;
          },
          "tau"_a)
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
            for (Index i = 0; i < tau.size(); i++) {
              out[0][i] = res(1, sorting[i]);
              out[1][i] = res(2, sorting[i]);
            }
            return out;
          },
          "tau"_a);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize disort\n", e.what()));
}
}  // namespace Python
