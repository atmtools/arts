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

void lm_setting_property(py::class_<LevenbergMarquardtSettings>& binding,
                         const char*                             name,
                         Numeric LevenbergMarquardtSettings::* member,
                         const char*                           description) {
  binding.def_prop_rw(
      name,
      [member](const LevenbergMarquardtSettings& settings) { return settings.*member; },
      [member](LevenbergMarquardtSettings& settings, Numeric value) {
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
  diagnostics.def_rw("status", &OptimalEstimationDiagnostics::status, "Named OEM outcome; NotRun before inversion.")
      .def_rw("initial_cost", &OptimalEstimationDiagnostics::initial_cost, "Starting total cost per measurement.")
      .def_rw("final_cost", &OptimalEstimationDiagnostics::final_cost, "Final total cost per measurement.")
      .def_rw("measurement_cost",
              &OptimalEstimationDiagnostics::measurement_cost,
              "Final measurement cost per measurement.")
      .def_rw("iterations",
              &OptimalEstimationDiagnostics::iterations,
              "Number of completed outer iterations, zero when not run.")
      .def_rw("lm_ga_history", &OptimalEstimationDiagnostics::lm_ga_history, "Initial and updated LM damping values.")
      .def_rw("errors", &OptimalEstimationDiagnostics::errors, "Errors and warnings recorded by OEM.");

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

Pass this object as ``ws.OEM(method="lm", lm_ga_settings=settings)``.
The same settings apply to ``lm_cg`` and the ``ml``/``ml_cg`` aliases.
Defaults provide an explicit starting configuration, with ordinary
convergence enabled only once damping reaches zero. They do not guarantee
convergence for every forward model. Use :meth:`describe` to inspect their
meaning and the current values.

Construction and every field assignment validate all six controls.
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
         Numeric                     convergence_damping_limit) {
        LevenbergMarquardtSettings value{
            .initial_damping           = initial_damping,
            .decrease_factor           = decrease_factor,
            .increase_factor           = increase_factor,
            .maximum_damping           = maximum_damping,
            .damping_threshold         = damping_threshold,
            .convergence_damping_limit = convergence_damping_limit,
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
                                   value.convergence_damping_limit);
           })
      .def("__setstate__",
           [](LevenbergMarquardtSettings*                                             settings,
              const std::tuple<Numeric, Numeric, Numeric, Numeric, Numeric, Numeric>& state) {
             const auto& [initial, decrease, increase, maximum, threshold, convergence] = state;
             LevenbergMarquardtSettings value{initial, decrease, increase, maximum, threshold, convergence};
             value.validate();
             new (settings) LevenbergMarquardtSettings(value);
           });

  auto jtdcmm = py::bind_map<JacobianTargetsDiagonalCovarianceMatrixMap, py::rv_policy::reference_internal>(
      m, "JacobianTargetsDiagonalCovarianceMatrixMap");
  generic_interface(jtdcmm);

  py::class_<PairOfBlockMatrix> pobm(m, "PairOfBlockMatrix");
  generic_interface(pobm);
  pobm.def_rw("first", &PairOfBlockMatrix::first, "Matrix\n\n.. :class:`BlockMatrix`");
  pobm.def_rw("second", &PairOfBlockMatrix::second, "Inverse of Matrix\n\n.. :class:`BlockMatrix`");

  py::class_<JacobianTargetType> jtt(m, "JacobianTargetType");
  jtt.def_rw("value", &JacobianTargetType::target, "Target\n\n.. :class:`object`");
  generic_interface(jtt);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize retrieval\n{}", e.what()));
}
}  // namespace Python
