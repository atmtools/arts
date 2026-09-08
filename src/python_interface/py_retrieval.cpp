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
void lm_setting_property(py::class_<OEMLMSettings>& binding,
                         const char*                name,
                         Numeric OEMLMSettings::* member,
                         const char*              description) {
  binding.def_prop_rw(
      name,
      [member](const OEMLMSettings& settings) { return settings.*member; },
      [member](OEMLMSettings& settings, Numeric value) {
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
  const OEMLMSettings       defaults;
  py::class_<OEMLMSettings> lm(m, "OEMLMSettings");
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
:meth:`validate` and conversion to :class:`Vector` check the settings again.
The object converts to the legacy six-element vector at the OEM boundary.
Use :meth:`as_vector` when storing the settings in a workspace variable.
See :ref:`sec-user-oem` for tuning guidance.
)";
  lm.def(
      "__init__",
      [](OEMLMSettings* settings,
         Numeric        initial_damping,
         Numeric        decrease_factor,
         Numeric        increase_factor,
         Numeric        maximum_damping,
         Numeric        damping_threshold,
         Numeric        convergence_damping_limit) {
        OEMLMSettings value{
            .initial_damping           = initial_damping,
            .decrease_factor           = decrease_factor,
            .increase_factor           = increase_factor,
            .maximum_damping           = maximum_damping,
            .damping_threshold         = damping_threshold,
            .convergence_damping_limit = convergence_damping_limit,
        };
        value.validate();
        new (settings) OEMLMSettings(value);
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
                      &OEMLMSettings::initial_damping,
                      R"(Initial damping gamma. Larger values restrain initial steps more strongly.
Must be finite, nonnegative, and no greater than ``maximum_damping``.
Zero starts with a Gauss--Newton step.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "decrease_factor",
                      &OEMLMSettings::decrease_factor,
                      R"(Divisor used when the local model warrants decreasing damping.
Must be finite and greater than one. A larger value releases damping faster.
Not every accepted step causes a decrease.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "increase_factor",
                      &OEMLMSettings::increase_factor,
                      R"(Multiplier used when an unsuccessful trial requires more damping.
Must be finite and greater than one. A larger value increases damping faster;
an increase from zero first restarts at ``damping_threshold``.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "maximum_damping",
                      &OEMLMSettings::maximum_damping,
                      R"(Upper damping limit. Must be finite and positive, and no smaller than
``initial_damping`` or ``damping_threshold``.
Failure to obtain an acceptable step at this value stops the retrieval.
Check the forward model, Jacobian, and covariance scales before raising it.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "damping_threshold",
                      &OEMLMSettings::damping_threshold,
                      R"(Positive restart damping and threshold for returning to Gauss--Newton.
A proposed decrease below this value sets damping to zero. A rejected trial
below this value restarts here. Must be finite and no greater than ``maximum_damping``.

.. :class:`float`
)");
  lm_setting_property(lm,
                      "convergence_damping_limit",
                      &OEMLMSettings::convergence_damping_limit,
                      R"(Largest updated damping at which the ordinary ``stop_dx`` test is enabled.
Must be finite and nonnegative. Zero waits until damping reaches zero.
A positive value permits convergence while damping still restrains steps,
which can hide a remaining distance to the minimum.

.. :class:`float`
)");
  lm.def("validate", &OEMLMSettings::validate, "Check the current named values and their coupled constraints.")
      .def("as_vector",
           &OEMLMSettings::as_vector,
           "Validate and return a copy in the legacy six-element lm_ga_settings order.")
      .def_static("from_vector",
                  &OEMLMSettings::from_vector,
                  "values"_a,
                  "Validate and import a legacy six-element lm_ga_settings vector.")
      .def("describe", &OEMLMSettings::describe, "Explain the current values, their effects, and tuning tradeoffs.")
      .def("__repr__", &OEMLMSettings::repr)
      .def("__str__", &OEMLMSettings::repr)
      .def("__copy__", [](const OEMLMSettings& value) { return value; })
      .def("__deepcopy__", [](const OEMLMSettings& value, py::dict&) { return value; })
      .def("__getstate__",
           [](const OEMLMSettings& value) {
             value.validate();
             return py::make_tuple(value.initial_damping,
                                   value.decrease_factor,
                                   value.increase_factor,
                                   value.maximum_damping,
                                   value.damping_threshold,
                                   value.convergence_damping_limit);
           })
      .def("__setstate__",
           [](OEMLMSettings* settings, const std::tuple<Numeric, Numeric, Numeric, Numeric, Numeric, Numeric>& state) {
             const auto& [initial, decrease, increase, maximum, threshold, convergence] = state;
             auto value =
                 OEMLMSettings::from_vector(Vector{initial, decrease, increase, maximum, threshold, convergence});
             new (settings) OEMLMSettings(value);
           });

  // The agenda parser explicitly constructs Vector from captured arguments;
  // direct workspace calls also need the registered implicit conversion.
  auto vector = py::borrow<py::class_<Vector>>(m.attr("Vector"));
  vector.def(
      "__init__",
      [](Vector* value, const OEMLMSettings& settings) {
        auto converted = settings.as_vector();
        new (value) Vector(std::move(converted));
      },
      "settings"_a,
      "Validate and convert named OEM damping controls to the legacy vector.");
  py::implicitly_convertible<OEMLMSettings, Vector>();

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
