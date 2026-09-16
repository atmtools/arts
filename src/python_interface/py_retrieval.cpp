#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <oem_settings.h>

#include "hpy_arts.h"
#include "python_interface.h"

namespace Python {
namespace {
LevenbergMarquardtSettings lm_settings_from_array(const std::array<Numeric, 6>& values) {
  LevenbergMarquardtSettings settings{.initial_damping           = values[0],
                                      .decrease_factor           = values[1],
                                      .increase_factor           = values[2],
                                      .maximum_damping           = values[3],
                                      .damping_threshold         = values[4],
                                      .convergence_damping_limit = values[5]};
  settings.validate();
  return settings;
}

template <typename T> void lm_setting_property(py::class_<LevenbergMarquardtSettings>& binding,
                                               const char*                             name,
                                               T LevenbergMarquardtSettings::* member,
                                               const char*                     description) {
  binding.def_prop_rw(
      name,
      [member](const LevenbergMarquardtSettings& settings) { return settings.*member; },
      [member](LevenbergMarquardtSettings& settings, T value) {
        // Validate before assignment so a failed edit preserves a usable object
        // and reports the named error before nanobind's implicit conversion.
        auto candidate    = settings;
        candidate.*member = value;
        candidate.validate();
        settings = candidate;
      },
      description);
}
}  // namespace

void py_retrieval(py::module_& m) try {
  py::class_<OptimalEstimationDiagnostics> diagnostics(m, "OptimalEstimationDiagnostics");
  generic_interface(diagnostics);
  diagnostics
      .def_rw("status",
              &OptimalEstimationDiagnostics::status,
              "Named OEM outcome; NotRun before inversion.\n\n.. :class:`~pyarts3.arts.OptimalEstimationStatus`")
      .def_rw("initial_cost",
              &OptimalEstimationDiagnostics::initial_cost,
              "Starting total cost per measurement.\n\n.. :class:`float`")
      .def_rw("final_cost",
              &OptimalEstimationDiagnostics::final_cost,
              "Final total cost per measurement.\n\n.. :class:`float`")
      .def_rw("measurement_cost",
              &OptimalEstimationDiagnostics::measurement_cost,
              "Final measurement cost per measurement.\n\n.. :class:`float`")
      .def_rw("iterations",
              &OptimalEstimationDiagnostics::iterations,
              "Number of completed outer iterations, zero when not run.\n\n.. :class:`int`")
      .def_rw("lm_ga_history",
              &OptimalEstimationDiagnostics::lm_ga_history,
              "Initial and updated LM damping values.\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("errors",
              &OptimalEstimationDiagnostics::errors,
              "Errors and warnings recorded by OEM.\n\n.. :class:`list[str]`");

  py::class_<OptimalEstimationData> data(m, "OptimalEstimationData");
  generic_interface(data);
  data.def_prop_ro("checked",
                   &OptimalEstimationData::checked,
                   "Whether input validation succeeded. Attribute replacement requires uncheck().\n\n.. :class:`bool`");
  data.def("check",
           &OptimalEstimationData::check,
           "jac_targets"_a.none() = py::none(),
           py::call_guard<py::gil_scoped_release>(),
           "Validate numerical inputs and optional finalized targets; report all problems or mark checked. "
           "In-place edits remain possible and require calling check() again when necessary.");
  data.def("uncheck", &OptimalEstimationData::uncheck, "Allow manual attribute replacement without changing values.");
  const auto member = [&]<typename T>(const char* name, T OptimalEstimationData::* field, const char* doc) {
    data.def_prop_rw(
        name,
        [field](OptimalEstimationData& value) -> T& { return value.*field; },
        [field, name](OptimalEstimationData& value, const T& replacement) {
          value.require_unchecked(std::format("replace member '{}'", name));
          value.*field = replacement;
        },
        py::for_getter(py::rv_policy::reference_internal),
        doc);
  };
  member(
      "covmat_diagonal_blocks",
      &OptimalEstimationData::covmat_diagonal_blocks,
      "Pending per-target covariance blocks used during target-based setup.\n\n.. :class:`~pyarts3.arts.JacobianTargetsDiagonalCovarianceMatrixMap`");
  member("measurement_vec",
         &OptimalEstimationData::measurement_vec,
         "Measurement vector.\n\n.. :class:`~pyarts3.arts.Vector`");
  member("model_state_vec_apriori",
         &OptimalEstimationData::model_state_vec_apriori,
         "A priori model state vector.\n\n.. :class:`~pyarts3.arts.Vector`");
  member("model_state_covmat",
         &OptimalEstimationData::model_state_covmat,
         "Covariance matrix of the model state.\n\n.. :class:`~pyarts3.arts.CovarianceMatrix`");
  member("measurement_vec_error_covmat",
         &OptimalEstimationData::measurement_vec_error_covmat,
         "Covariance matrix of the measurement vector error.\n\n.. :class:`~pyarts3.arts.CovarianceMatrix`");
  member("model_state_vec",
         &OptimalEstimationData::model_state_vec,
         "Current model state vector.\n\n.. :class:`~pyarts3.arts.Vector`");
  member("measurement_vec_fit",
         &OptimalEstimationData::measurement_vec_fit,
         "Fitted measurement vector.\n\n.. :class:`~pyarts3.arts.Vector`");
  member("measurement_jac",
         &OptimalEstimationData::measurement_jac,
         "Jacobian of the measurement operator.\n\n.. :class:`~pyarts3.arts.Matrix`");
  member("model_state_basis_mat",
         &OptimalEstimationData::model_state_basis_mat,
         "Basis matrix of the model state.\n\n.. :class:`~pyarts3.arts.BlockMatrix`");
  member("measurement_basis_mat",
         &OptimalEstimationData::measurement_basis_mat,
         "Basis matrix of the measurement.\n\n.. :class:`~pyarts3.arts.BlockMatrix`");
  member("model_state_covmat_normalization",
         &OptimalEstimationData::model_state_covmat_normalization,
         "Normalization vector for the model state covariance matrix.\n\n.. :class:`~pyarts3.arts.Vector`");
  member("measurement_vec_normalization",
         &OptimalEstimationData::measurement_vec_normalization,
         "Normalization vector for the measurement vector.\n\n.. :class:`~pyarts3.arts.Vector`");
  member("measurement_gain_mat",
         &OptimalEstimationData::measurement_gain_mat,
         "Gain matrix of the measurement.\n\n.. :class:`~pyarts3.arts.Matrix`");
  member("measurement_averaging_kernel",
         &OptimalEstimationData::measurement_averaging_kernel,
         "Averaging kernel of the measurement.\n\n.. :class:`~pyarts3.arts.Matrix`");
  member("observation_error_covmat",
         &OptimalEstimationData::observation_error_covmat,
         "Covariance matrix of the observation error.\n\n.. :class:`~pyarts3.arts.Matrix`");
  member("smoothing_error_covmat",
         &OptimalEstimationData::smoothing_error_covmat,
         "Covariance matrix of the smoothing error.\n\n.. :class:`~pyarts3.arts.Matrix`");
  member("diagnostics",
         &OptimalEstimationData::diagnostics,
         "Diagnostics information.\n\n.. :class:`~pyarts3.arts.OptimalEstimationDiagnostics`");
  member("basis_singular_values",
         &OptimalEstimationData::basis_singular_values,
         "Singular values of the basis matrix.\n\n.. :class:`~pyarts3.arts.Vector`");
  member("basis_lost_dofs",
         &OptimalEstimationData::basis_lost_dofs,
         "Degrees of freedom lost in the basis matrix.\n\n.. :class:`float`");
  member("basis_lost_information_bits",
         &OptimalEstimationData::basis_lost_information_bits,
         "Information bits lost in the basis matrix.\n\n.. :class:`float`");
  data.def("clear_auxiliary",
           &OptimalEstimationData::clear_auxiliary,
           "Release recomputable products and caches, preserving inputs, current state, bases and diagnostics.");
  data.def(
      "clear",
      [](OptimalEstimationData& value) {
        value.require_unchecked("clear data");
        value.clear();
      },
      "Reset all inputs and results.");

  const LevenbergMarquardtSettings       defaults;
  py::class_<LevenbergMarquardtSettings> lm(m, "LevenbergMarquardtSettings");
  lm.def(py::init<const LevenbergMarquardtSettings&>());
  lm.def(
      "__init__",
      [](LevenbergMarquardtSettings* settings, const std::array<Numeric, 6>& values) {
        new (settings) LevenbergMarquardtSettings(lm_settings_from_array(values));
      },
      "values"_a);
  py::implicitly_convertible<std::array<Numeric, 6>, LevenbergMarquardtSettings>();

  xml_interface(lm);
  lm.doc() = R"(Named Levenberg--Marquardt damping controls for OEM.

Pass this object as ``ws.oemCalc(settings=OptimalEstimationSettings(method="lm", lm=settings))``.
The same settings apply to ``lm_cg`` and the ``ml``/``ml_cg`` aliases.
Defaults provide an explicit starting configuration, with ordinary
convergence enabled only once damping reaches zero. They do not guarantee
convergence for every forward model. Use :meth:`describe` to inspect their
meaning and the current values.

Construction and every field assignment validate all controls.
Invalid edits raise an error naming the setting and preserve the previous
configuration. When increasing ``initial_damping`` beyond the current
maximum, raise ``maximum_damping`` first. To change several controls at
once, construct a replacement object with their named arguments.
:meth:`validate` checks the settings again. OEM accepts this type directly
and validates it before a Levenberg-Marquardt retrieval.
Python also accepts a six-value sequence directly or through the constructor.
Its order is initial_damping, decrease_factor, increase_factor,
maximum_damping, damping_threshold, convergence_damping_limit.
This input shorthand does not add a vector conversion to the C++ type.
See :ref:`sec-user-oem` for tuning guidance.
)";
  lm.def(
      "__init__",
      [](LevenbergMarquardtSettings* settings,
         Numeric                     initial_damping,
         Numeric                     decrease_factor,
         Numeric                     increase_factor,
         Numeric                     maximum_damping,
         Numeric                     damping_threshold,
         Numeric                     convergence_damping_limit,
         Index                       maximum_trials) {
        LevenbergMarquardtSettings value{
            .initial_damping           = initial_damping,
            .decrease_factor           = decrease_factor,
            .increase_factor           = increase_factor,
            .maximum_damping           = maximum_damping,
            .damping_threshold         = damping_threshold,
            .convergence_damping_limit = convergence_damping_limit,
            .maximum_trials            = maximum_trials,
        };
        value.validate();
        new (settings) LevenbergMarquardtSettings(value);
      },
      py::kw_only(),
      "initial_damping"_a           = defaults.initial_damping,
      "decrease_factor"_a           = defaults.decrease_factor,
      "increase_factor"_a           = defaults.increase_factor,
      "maximum_damping"_a           = defaults.maximum_damping,
      "damping_threshold"_a         = defaults.damping_threshold,
      "convergence_damping_limit"_a = defaults.convergence_damping_limit,
      "maximum_trials"_a            = defaults.maximum_trials,
      "Construct validated damping controls using named arguments.");
  lm_setting_property(lm,
                      "initial_damping",
                      &LevenbergMarquardtSettings::initial_damping,
                      R"(Initial damping gamma. Larger values restrain initial steps more strongly.
Must be finite, nonnegative, and no greater than ``maximum_damping``.
Zero starts with a Gauss--Newton step.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "decrease_factor",
                      &LevenbergMarquardtSettings::decrease_factor,
                      R"(Divisor used when the local model warrants decreasing damping.
Must be finite and greater than one. A larger value releases damping faster.
Not every accepted step causes a decrease.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "increase_factor",
                      &LevenbergMarquardtSettings::increase_factor,
                      R"(Multiplier used when an unsuccessful trial requires more damping.
Must be finite and greater than one. A larger value increases damping faster;
an increase from zero first restarts at ``damping_threshold``.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "maximum_damping",
                      &LevenbergMarquardtSettings::maximum_damping,
                      R"(Upper damping limit. Must be finite and positive, and no smaller than
``initial_damping`` or ``damping_threshold``.
Failure to obtain an acceptable step at this value stops the retrieval.
Check the forward model, Jacobian, and covariance scales before raising it.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "damping_threshold",
                      &LevenbergMarquardtSettings::damping_threshold,
                      R"(Positive restart damping and threshold for returning to Gauss--Newton.
A proposed decrease below this value sets damping to zero. A rejected trial
below this value restarts here. Must be finite and no greater than ``maximum_damping``.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "convergence_damping_limit",
                      &LevenbergMarquardtSettings::convergence_damping_limit,
                      R"(Largest updated damping at which the ordinary ``stop_dx`` test is enabled.
Must be finite and nonnegative. Zero waits until damping reaches zero.
A positive value permits convergence while damping still restrains steps,
which can hide a remaining distance to the minimum.

.. :class:`float`
)");
  lm_setting_property(
      lm,
      "maximum_trials",
      &LevenbergMarquardtSettings::maximum_trials,
      "Maximum linear-solve trials per outer iteration, including stationarity checks.\n\n.. :class:`int`");
  lm.def("validate",
         &LevenbergMarquardtSettings::validate,
         "Check the current named values and their coupled constraints.")
      .def("describe",
           &LevenbergMarquardtSettings::describe,
           "Explain the current values, their effects, and tuning tradeoffs.")
      .def("__repr__", &LevenbergMarquardtSettings::repr)
      .def("__str__", &LevenbergMarquardtSettings::repr)
      .def("__copy__", [](const LevenbergMarquardtSettings& value) { return value; })
      .def("__deepcopy__", [](const LevenbergMarquardtSettings& value, py::dict&) { return value; })
      .def("__getstate__",
           [](const LevenbergMarquardtSettings& value) {
             value.validate();
             return py::make_tuple(value.initial_damping,
                                   value.decrease_factor,
                                   value.increase_factor,
                                   value.maximum_damping,
                                   value.damping_threshold,
                                   value.convergence_damping_limit,
                                   value.maximum_trials);
           })
      .def("__setstate__",
           [](LevenbergMarquardtSettings*                                                    settings,
              const std::tuple<Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Index>& state) {
             const auto& [initial, decrease, increase, maximum, threshold, convergence, trials] = state;
             LevenbergMarquardtSettings value{.initial_damping           = initial,
                                              .decrease_factor           = decrease,
                                              .increase_factor           = increase,
                                              .maximum_damping           = maximum,
                                              .damping_threshold         = threshold,
                                              .convergence_damping_limit = convergence,
                                              .maximum_trials            = trials};
             value.validate();
             new (settings) LevenbergMarquardtSettings(value);
           });

  const OptimalEstimationSettings       calculation_defaults;
  py::class_<OptimalEstimationSettings> settings(m, "OptimalEstimationSettings");
  generic_interface(settings);
  settings.def(
      "__init__",
      [](OptimalEstimationSettings* self, const std::string& method) {
        OptimalEstimationSettings value{.method = to<OptimalEstimationMethod>(method)};
        value.validate();
        new (self) OptimalEstimationSettings(std::move(value));
      },
      "method"_a);
  py::implicitly_convertible<std::string, OptimalEstimationSettings>();

  settings.def(
      "__init__",
      [](OptimalEstimationSettings* self,
         OptimalEstimationMethod    method,
         Index                      max_iter,
         Numeric                    stop_dx,
         Numeric                    max_start_cost,
         Numeric                    cg_tolerance,
         Index                      cg_max_iter,
         LevenbergMarquardtSettings lm,
         Index                      display_progress,
         bool                       clear_matrices) {
        OptimalEstimationSettings value{.method           = method,
                                        .max_iter         = max_iter,
                                        .stop_dx          = stop_dx,
                                        .max_start_cost   = max_start_cost,
                                        .cg_tolerance     = cg_tolerance,
                                        .cg_max_iter      = cg_max_iter,
                                        .lm               = lm,
                                        .display_progress = display_progress,
                                        .clear_matrices   = clear_matrices};
        value.validate();
        new (self) OptimalEstimationSettings(std::move(value));
      },
      py::kw_only(),
      "method"_a           = calculation_defaults.method,
      "max_iter"_a         = calculation_defaults.max_iter,
      "stop_dx"_a          = calculation_defaults.stop_dx,
      "max_start_cost"_a   = calculation_defaults.max_start_cost,
      "cg_tolerance"_a     = calculation_defaults.cg_tolerance,
      "cg_max_iter"_a      = calculation_defaults.cg_max_iter,
      "lm"_a               = calculation_defaults.lm,
      "display_progress"_a = calculation_defaults.display_progress,
      "clear_matrices"_a   = calculation_defaults.clear_matrices);
  settings.def("validate",
               &OptimalEstimationSettings::validate,
               "Check all controls; oemCalc and oemCalcReduced also validate before calculation.");
  settings.def_rw("method",
                  &OptimalEstimationSettings::method,
                  "Algorithm and linear solver.\n\n.. :class:`~pyarts3.arts.OptimalEstimationMethod`");
  settings.def_rw("max_iter", &OptimalEstimationSettings::max_iter, "Maximum outer iterations.\n\n.. :class:`int`");
  settings.def_rw(
      "stop_dx", &OptimalEstimationSettings::stop_dx, "Positive convergence tolerance.\n\n.. :class:`float`");
  settings.def_rw("max_start_cost",
                  &OptimalEstimationSettings::max_start_cost,
                  "Maximum initial cost; infinity disables the cutoff.\n\n.. :class:`float`");
  settings.def_rw("cg_tolerance",
                  &OptimalEstimationSettings::cg_tolerance,
                  "Positive relative CG residual tolerance.\n\n.. :class:`float`");
  settings.def_rw("cg_max_iter",
                  &OptimalEstimationSettings::cg_max_iter,
                  "CG iteration limit; zero selects max(1000, 2 * dimension).\n\n.. :class:`int`");
  settings.def_rw(
      "lm", &OptimalEstimationSettings::lm, "LM damping and trial controls.\n\n.. :class:`~pyarts3.arts.LevenbergMarquardtSettings`");
  settings.def_rw(
      "display_progress", &OptimalEstimationSettings::display_progress, "Print progress when 1.\n\n.. :class:`int`");
  settings.def_rw("clear_matrices",
                  &OptimalEstimationSettings::clear_matrices,
                  "Release Jacobian and gain when true.\n\n.. :class:`bool`");
  settings
      .def("__getstate__",
           [](const OptimalEstimationSettings& value) {
             return py::make_tuple(std::string{toString(value.method)},
                                   value.max_iter,
                                   value.stop_dx,
                                   value.max_start_cost,
                                   value.cg_tolerance,
                                   value.cg_max_iter,
                                   value.lm,
                                   value.display_progress,
                                   value.clear_matrices);
           })
      .def("__setstate__",
           [](OptimalEstimationSettings* self,
              const std::
                  tuple<std::string, Index, Numeric, Numeric, Numeric, Index, LevenbergMarquardtSettings, Index, bool>&
                      state) {
             auto value = std::apply(
                 [](const std::string& method, auto... values) {
                   return OptimalEstimationSettings{to<OptimalEstimationMethod>(method), values...};
                 },
                 state);
             value.validate();
             new (self) OptimalEstimationSettings(std::move(value));
           });

  auto jtdcmm = py::bind_map<JacobianTargetsDiagonalCovarianceMatrixMap, py::rv_policy::reference_internal>(
      m, "JacobianTargetsDiagonalCovarianceMatrixMap");
  generic_interface(jtdcmm);

  py::class_<PairOfBlockMatrix> pobm(m, "PairOfBlockMatrix");
  generic_interface(pobm);
  pobm.def_rw("first", &PairOfBlockMatrix::first, "Matrix\n\n.. :class:`~pyarts3.arts.BlockMatrix`");
  pobm.def_rw("second", &PairOfBlockMatrix::second, "Inverse of Matrix\n\n.. :class:`~pyarts3.arts.BlockMatrix`");

  py::class_<JacobianTargetType> jtt(m, "JacobianTargetType");
  jtt.def_rw("value", &JacobianTargetType::target, "Target\n\n.. :class:`object`");
  generic_interface(jtt);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize retrieval\n{}", e.what()));
}
}  // namespace Python
