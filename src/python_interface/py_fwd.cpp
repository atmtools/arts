#include <fwd.h>

#include <python_interface.h>

#include "debug.h"
#include "py_macros.h"

namespace Python {
void py_fwd(py::module_& m) try {
  artsclass<SpectralRadianceProfileOperator>(m, "SpectralRadianceProfileOperator")
      .def(py::init<>(), "From nothing")
      .PythonInterfaceCopyValue(SpectralRadianceProfileOperator)
      .PythonInterfaceWorkspaceVariableConversion(SpectralRadianceProfileOperator)
      .PythonInterfaceFileIO(SpectralRadianceProfileOperator)
      .PythonInterfaceBasicRepresentation(SpectralRadianceProfileOperator)
      .def(
          "planar",
          [](const SpectralRadianceProfileOperator& fwd_prof, Numeric f, Numeric za) {
            return fwd_prof.planar(f, za);
          },
          py::arg("f"),
          py::arg("za"),
          py::doc(R"--(Plane parallel computer of spectral radiance profile

Parameters
----------
f : float
    Frequency [Hz]
za : float
    Zenith angle [degrees]

Returns
-------
prof : ~pyarts.arts.Vector
    Profile of radiation
)--"))
      .def(
          "planar_par",
          [](const SpectralRadianceProfileOperator& fwd_prof, const Vector& f, Numeric za) {
            return fwd_prof.planar_par(f, za);
          },
          py::arg("f"),
          py::arg("za"),
          py::doc(R"--(Plane parallel computer of spectral radiance profile

Runs on parallel in frequency

Parameters
----------
f :  ~pyarts.arts.Vector
    Frequency [Hz]
za : float
    Zenith angle [degrees]

Returns
-------
prof : ~pyarts.arts.Matrix
    Profile of radiation
)--"))
    .def_readonly("altitude", &SpectralRadianceProfileOperator::altitude, py::doc(R"--(Altitude grid [m])--"))
      .PythonInterfaceWorkspaceDocumentation(SpectralRadianceProfileOperator);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize fwd\n", e.what()));
}
}  // namespace Python