#include <lineshapemodel.h>
#include <nanobind/stl/bind_vector.h>
#include <python_interface.h>
#include <zeemandata.h>

#include <algorithm>
#include <memory>
#include <string>
#include <utility>

#include "absorptionlines.h"
#include "cia.h"
#include "debug.h"
#include "hpy_arts.h"
#include "isotopologues.h"
#include "physics_funcs.h"
#include "py_macros.h"
#include "species_tags.h"

namespace Python {
void py_cia(py::module_& m) try {
  py::class_<CIARecord> cia(m, "CIARecord");
  workspace_group_interface(cia);
  cia.def(py::init<ArrayOfGriddedField2, SpeciesEnum, SpeciesEnum>())
      .def_prop_ro(
          "specs",
          [](const CIARecord& c) { return c.TwoSpecies(); },
          ":class:`list` The two species")
      .def_prop_ro(
          "data",
          [](const CIARecord& c) { return c.Data(); },
          ":class:`~pyarts.arts.ArrayOfGriddedField2` Data by bands")
      .def(
          "compute_abs",
          [](CIARecord& cia_,
             Numeric T,
             Numeric P,
             Numeric X0,
             Numeric X1,
             const Vector& f,
             Numeric T_extrapolfac,
             Index robust) {
            Vector out(f.nelem(), 0);

            for (auto& cia_data : cia_.Data()) {
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
          py::arg("robust")        = 1,
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

)--")
      .def("__getstate__",
           [](const CIARecord& t) {
             return std::make_tuple(t.Data(), t.TwoSpecies());
           })
      .def(
          "__setstate__",
          [](CIARecord* c,
             const std::tuple<ArrayOfGriddedField2, std::array<SpeciesEnum, 2>>&
                 state) {
            new (c) CIARecord();
            c->Data()       = std::get<0>(state);
            c->TwoSpecies() = std::get<1>(state);
          });

  auto acr = py::bind_vector<ArrayOfCIARecord, py::rv_policy::reference_internal>(m, "ArrayOfCIARecord");
  workspace_group_interface(acr);
  vector_interface(acr);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize cia\n", e.what()));
}
}  // namespace Python
