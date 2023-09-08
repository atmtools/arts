#include <algorithm>
#include <lineshape.h>
#include <lineshapemodel.h>
#include <python_interface.h>
#include <zeemandata.h>

#include <memory>
#include <string>
#include <utility>

#include "absorptionlines.h"
#include "cia.h"
#include "debug.h"
#include "isotopologues.h"
#include "physics_funcs.h"
#include "py_macros.h"
#include "quantum_numbers.h"
#include "species_tags.h"

namespace Python {
void py_cia(py::module_& m) try {
  artsclass<CIARecord>(m, "CIARecord")
      .def(py::init([]() { return std::make_shared<CIARecord>(); }),
           "Empty record")
      .def(py::init([](const ArrayOfGriddedField2& data,
                       Species::Species spec1,
                       Species::Species spec2) {
             return std::make_shared<CIARecord>(data, spec1, spec2);
           }),
           "From values")
      .PythonInterfaceCopyValue(CIARecord)
      .PythonInterfaceWorkspaceVariableConversion(CIARecord)
      .PythonInterfaceBasicRepresentation(CIARecord)
      .PythonInterfaceFileIO(CIARecord)
      .def_property_readonly(
          "specs",
          [](const CIARecord& c) { return c.TwoSpecies(); },
          ":class:`list` The two species")
      .def_property_readonly(
          "data",
          [](const CIARecord& c) { return c.Data(); },
          ":class:`~pyarts.arts.ArrayOfGriddedField2` Data by bands")
      .def(
          "compute_abs",
          [](CIARecord& cia,
             Numeric T,
             Numeric P,
             Numeric X0,
             Numeric X1,
             const Vector& f,
             Numeric T_extrapolfac,
             Index robust) {
            Vector out(f.nelem(), 0);

            for (auto& cia_data: cia.Data()) {
              Vector result(f.nelem(), 0);
              cia_interpolation(result, f, T, cia_data, T_extrapolfac, robust);
              out += result;
            }

            out *= Math::pow2(number_density(P, T)) * X0 * X1;
            return out;
          },
          py::arg("T"),
          py::arg("P"),
          py::arg("X0"),
          py::arg("X1"),
          py::arg("f"),
          py::arg("T_extrapolfac") = 0.0,
          py::arg("robust") = 1,
          py::doc(
              R"--(Computes the collision-induced absorption in 1/m

Parameters
----------
T : Numeric
    Temperature [K]
P : Numeric
    Pressure [Pa]
X0 : Numeric
    VMR of species 1 [-]
X1 : Numeric
    VMR of species 2 [-]
f : Vector
    Frequency grid [Hz]
T_extrapolfac : Numeric, optional
    Extrapolation in temperature.  The default is 0
robust : Index, optional
    Returns NaN instead of throwing if it evaluates true.  The default is 1

Returns
-------
  abs : Vector
    Absorption profile [1/m]

)--"))
      .def(py::pickle(
          [](const CIARecord& t) {
            return py::make_tuple(t.Data(), t.TwoSpecies());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")
            auto out = std::make_shared<CIARecord>();
            out->Data() = t[0].cast<ArrayOfGriddedField2>();
            out->TwoSpecies() = t[1].cast<std::array<Species::Species, 2>>();
            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(CIARecord);

  PythonInterfaceWorkspaceArray(CIARecord);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize cia\n", e.what()));
}
}  // namespace Python
