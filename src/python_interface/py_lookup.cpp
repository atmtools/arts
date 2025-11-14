#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/shared_ptr.h>
#include <python_interface.h>

#include "hpy_arts.h"

namespace Python {
void py_lookup(py::module_& m) try {
  py::class_<AbsorptionLookupTable> alt(m, "AbsorptionLookupTable");
  generic_interface(alt);
  alt.def_rw("f_grid",
             &AbsorptionLookupTable::f_grid,
             "The frequency grid in Hz\n\n.. :class:`AscendingGrid`");
  alt.def_rw(
      "log_p_grid",
      &AbsorptionLookupTable::log_p_grid,
      "The pressure grid in log Pa [same dimension as atm]\n\n.. :class:`DescendingGrid`");
  alt.def_rw(
      "t_pert",
      &AbsorptionLookupTable::t_pert,
      "The temperature perturbation grid in K [any number of elements or empty for nothing]\n\n.. :class:`AscendingGrid`");
  alt.def_rw(
      "w_pert",
      &AbsorptionLookupTable::w_pert,
      "The humidity perturbation grid in fractional units [any number of elements or empty for nothing]\n\n.. :class:`AscendingGrid`");
  alt.def_rw(
      "water_atmref",
      &AbsorptionLookupTable::water_atmref,
      "Local grids so that pressure interpolation may work\n\n.. :class:`Vector`");
  alt.def_rw(
      "t_atmref",
      &AbsorptionLookupTable::t_atmref,
      "Local grids so that pressure interpolation may work\n\n.. :class:`Vector`");
  alt.def_rw("xsec",
             &AbsorptionLookupTable::xsec,
             "The absorption cross section table\n\n.. :class:`Tensor4`");

  auto alts = py::bind_map<AbsorptionLookupTables>(m, "AbsorptionLookupTables");
  generic_interface(alts);

  alts.def(
      "spectral_propmat",
      [](const AbsorptionLookupTables& self,
         const AscendingGrid& f,
         const AtmPoint& atm,
         const SpeciesEnum& spec,
         const Index& no_negative_absorption,
         const Index& p_interp_order,
         const Index& t_interp_order,
         const Index& water_interp_order,
         const Index& f_interp_order,
         const Numeric& extpolfac,
         const py::kwargs&) {
        PropmatVector spectral_propmat(f.size());
        PropmatMatrix spectral_propmat_jac(0, f.size());
        JacobianTargets jac_targets{};

        spectral_propmatAddLookup(spectral_propmat,
                                  spectral_propmat_jac,
                                  f,
                                  jac_targets,
                                  spec,
                                  self,
                                  atm,
                                  no_negative_absorption,
                                  p_interp_order,
                                  t_interp_order,
                                  water_interp_order,
                                  f_interp_order,
                                  extpolfac);

        return spectral_propmat;
      },
      "f"_a,
      "atm"_a,
      "spec"_a                   = SpeciesEnum::Bath,
      "no_negative_absorption"_a = Index{1},
      "p_interp_order"_a         = Index{7},
      "t_interp_order"_a         = Index{7},
      "water_interp_order"_a     = Index{7},
      "f_interp_order"_a         = Index{7},
      "extpolfac"_a              = Numeric{0.5},
      "kwargs"_a                 = py::kwargs{},
      R"--(Computes the line-by-line model absorption in 1/m

Parameters
----------
f : AscendingGrid
    Frequency grid [Hz]
atm : AtmPoint
    Atmospheric point
spec : SpeciesEnum, optional
    Species to use.  Defaults to all.
no_negative_absorption : Index, optional
    If 1, the absorption is set to zero if it is negative. The default is 1.
p_interp_order : Index, optional
    Order of interpolation in pressure.  The default is 7.
t_interp_order : Index, optional
    Order of interpolation in temperature.  The default is 7.
water_interp_order : Index, optional
    Order of interpolation in water VMR.  The default is 7.
f_interp_order : Index, optional
    Order of interpolation in frequency.  The default is 7.
extpolfac : Numeric, optional
    How far the grids are allowed to be extended in relative distance of the outer two most grid points.  Default is 0.5.

Returns
-------
spectral_propmat : PropmatVector
    Propagation matrix by frequency [1/m]

)--");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize lookup\n{}", e.what()));
}
}  // namespace Python
