#include <predefined/predef.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

void internalCKDMT350(py::module_& predef) {
  auto m = predef.def_submodule("CKDMT350");
  m.doc() = "Contains implementations for using MT CKD version 3.50";

  m.def(
      "get_self_h2o",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::CKDMT350::compute_self_h2o(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes self broadening using MT CKD version 3.50

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_foreign_h2o",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::CKDMT350::compute_foreign_h2o(
            pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes foreign broadening using MT CKD version 3.50

Parameters:
    f_grid : Vector
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));
}

void py_predefined(py::module_& m) {
  auto predef = m.def_submodule("predef");
  predef.doc() = "Contains predefined absorption models";

  internalCKDMT350(predef);
}
}  // namespace Python