#include <python_interface.h>

#include "atm.h"
#include "hpy_arts.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "nanobind/stl/bind_vector.h"
#include "nanobind/stl/function.h"
#include "nanobind/stl/variant.h"
#include "surf.h"

NB_MAKE_OPAQUE(std::vector<Jacobian::AtmTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::SurfaceTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::LineTarget>)

namespace Python {
template <typename T>
void target(py::class_<T>& cls) {
  cls.def_ro("type", &T::type);
  cls.def_ro("d", &T::d);
  cls.def_ro("target_pos", &T::target_pos);
  cls.def_ro("x_start", &T::x_start);
  cls.def_ro("x_size", &T::x_size);

  str_interface(cls);
}

void py_jac(py::module_& m) try {
  py::class_<Jacobian::AtmTarget> atm(m, "AtmTarget");
  target(atm);
  atm.def(
      "model_state",
      [](const Jacobian::AtmTarget& t, const AtmField& f) {
        Vector x(t.x_size);
        t.set_state(x, f, t.type);
        return x;
      },
      R"--(Invokes the internal model-state-setter

Warning:

    May segment-fault if there's a discrepancy between the data container and
    the target.  Such discrepancies should just happen if you have a manual
    setter for the target, or if you have changed the data container since you
    created the target.

Parameters
----------
t : AtmField
    The atmospheric data container

Returns
-------
x : Vector
    The local state vector
)--");

  py::class_<Jacobian::SurfaceTarget> surf(m, "SurfaceTarget");
  target(surf);
  surf.def(
      "model_state",
      [](const Jacobian::SurfaceTarget& t, const SurfaceField& f) {
        Vector x(t.x_size);
        t.set_state(x, f, t.type);
        return x;
      },
      R"--(Invokes the internal model-state-setter

Warning:

    May segment-fault if there's a discrepancy between the data container and
    the target.  Such discrepancies should just happen if you have a manual
    setter for the target, or if you have changed the data container since you
    created the target.

Parameters
----------
t : SurfaceField
    The surface data container

Returns
-------
x : Vector
    The local state vector
)--");

  py::class_<Jacobian::LineTarget> line(m, "LineTarget");
  target(line);
  line.def(
      "model_state",
      [](const Jacobian::LineTarget& t, const ArrayOfAbsorptionBand& f) {
        Vector x(t.x_size);
        t.set_state(x, f, t.type);
        return x;
      },
      R"--(Invokes the internal model-state-setter

Warning:

    May segment-fault if there's a discrepancy between the data container and
    the target.  Such discrepancies should just happen if you have a manual
    setter for the target, or if you have changed the data container since you
    created the target.

Parameters
----------
t : ArrayOfAbsorptionBand
    The absorption data container

Returns
-------
x : Vector
    The local state vector
)--");

  py::bind_vector<std::vector<Jacobian::AtmTarget>>(m, "ArrayOfAtmTargets");
  py::bind_vector<std::vector<Jacobian::SurfaceTarget>>(m,
                                                        "ArrayOfSurfaceTarget");
  py::bind_vector<std::vector<Jacobian::LineTarget>>(m, "ArrayOfLineTarget");

  py::class_<JacobianTargets> jacs(m, "JacobianTargets");
  workspace_group_interface(jacs);
  jacs.def_prop_rw(
      "atm",
      [](JacobianTargets& j) { return j.atm(); },
      [](JacobianTargets& j, std::vector<Jacobian::AtmTarget> t) {
        j.atm() = std::move(t);
      },
      "List of atmospheric targets");
  jacs.def_prop_rw(
      "surf",
      [](JacobianTargets& j) { return j.surf(); },
      [](JacobianTargets& j, std::vector<Jacobian::SurfaceTarget> t) {
        j.surf() = std::move(t);
      },
      "List of surface targets");
  jacs.def_prop_rw(
      "line",
      [](JacobianTargets& j) { return j.line(); },
      [](JacobianTargets& j, std::vector<Jacobian::LineTarget> t) {
        j.line() = std::move(t);
      },
      "List of line targets");
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize jac\n", e.what()));
}
}  // namespace Python