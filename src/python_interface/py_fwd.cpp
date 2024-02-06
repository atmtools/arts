#include <fwd.h>
#include <python_interface.h>

#include "debug.h"
#include "fwd_spectral_radiance.h"
#include "py_macros.h"
#include "rtepack.h"

namespace Python {
void py_fwd(py::module_& m) try {
  auto fwd = m.def_submodule("fwd");

  artsclass<SpectralRadianceOperator::PosDistance>(
      fwd, "SpectralRadiancePosDistance");

  artsclass<SpectralRadianceOperator::PathVector>(fwd,
                                                  "SpectralRadiancePathVector")
      .def("__getitem__",
           [](const SpectralRadianceOperator::PathVector& v, const Size i) {
             return v[i];
           })
      .def("__setitem__",
           [](SpectralRadianceOperator::PathVector& v,
              const Size i,
              const SpectralRadianceOperator::PosDistance& p) { v[i] = p; });

  artsclass<SpectralRadianceOperator>(m, "SpectralRadianceOperator")
      .def(py::init<>())
      .PythonInterfaceCopyValue(SpectralRadianceOperator)
      .PythonInterfaceWorkspaceVariableConversion(SpectralRadianceOperator)
      .PythonInterfaceBasicRepresentation(SpectralRadianceOperator)
      .PythonInterfaceFileIO(SpectralRadianceOperator)
      .def("geometric_atm",
           [](const SpectralRadianceOperator& fwd,
              const Vector3 pos,
              const Vector2 los) {
             const auto path = fwd.geometric_planar(pos[0], los[0]);
             return fwd.get_atm(path);
           })
      .def("geometric_planar",
           [](const SpectralRadianceOperator& fwd,
              const Numeric frequency,
              const Vector3 pos,
              const Vector2 los) {
             const auto path = fwd.geometric_planar(pos[0], los[0]);
             return fwd(frequency, los, path);
           })
      .def("geometric_planar",
           [](const SpectralRadianceOperator& fwd,
              const Vector& frequency,
              const Vector3 pos,
              const Vector2 los) {
             const auto path = fwd.geometric_planar(pos[0], los[0]);
             StokvecVector out(frequency.size());
             std::transform(
                 frequency.begin(),
                 frequency.end(),
                 out.begin(),
                 [&](const Numeric& f) { return fwd(f, los, path); });
             return out;
           })
      .def_property_readonly("altitude", &SpectralRadianceOperator::altitude)
      .PythonInterfaceWorkspaceDocumentation(SpectralRadianceOperator);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize fwd\n", e.what()));
}
}  // namespace Python