#include <fwd.h>
#include <python_interface.h>

#include "arts_omp.h"
#include "debug.h"
#include "py_macros.h"
#include "rtepack.h"

namespace Python {
void py_fwd(py::module_& m) try {
  artsclass<SpectralRadianceOperator>(m, "SpectralRadianceOperator")
      .def(py::init<>())
      .PythonInterfaceCopyValue(SpectralRadianceOperator)
      .PythonInterfaceWorkspaceVariableConversion(SpectralRadianceOperator)
      .PythonInterfaceBasicRepresentation(SpectralRadianceOperator)
      .PythonInterfaceFileIO(SpectralRadianceOperator)
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

            if (arts_omp_in_parallel()) {
             std::transform(
                 frequency.begin(),
                 frequency.end(),
                 out.begin(),
                 [&](const Numeric& f) { return fwd(f, los, path); });
                 } else {
                  String error;
                  #pragma omp parallel for
                  for (Index i=0; i<frequency.size(); i++) {
                    try {
                      out[i] = fwd(frequency[i], los, path);
                    } catch (std::exception& e) {
                      #pragma omp critical
                      error += e.what() + String("\n");
                    }
                  }
                  ARTS_USER_ERROR_IF(not error.empty(), error);
                 }
             return out;
           })
      .def_property_readonly("altitude", &SpectralRadianceOperator::altitude)
      .PythonInterfaceWorkspaceDocumentation(SpectralRadianceOperator);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize fwd\n", e.what()));
}
}  // namespace Python