#include <python_interface.h>

#include "mc_antenna.h"
#include "py_macros.h"

namespace Python {
void py_mcantenna(py::module_& m) try {
  py::enum_<AntennaType>(m, "AntennaType")
      .value("ANTENNA_TYPE_PENCIL_BEAM",
             AntennaType::ANTENNA_TYPE_PENCIL_BEAM,
             "As pencil beam")
      .value("ANTENNA_TYPE_GAUSSIAN",
             AntennaType::ANTENNA_TYPE_GAUSSIAN,
             "As gaussian beam")
      .value("ANTENNA_TYPE_LOOKUP",
             AntennaType::ANTENNA_TYPE_LOOKUP,
             "As from a lookup")
      .PythonInterfaceCopyValue(AntennaType)
      .def("__getstate__",
           [](const AntennaType& self) {
             return std::make_tuple(static_cast<Index>(self));
           })
      .def("__setstate__",
           [](AntennaType* a, const std::tuple<Index>& state) {
             new (a) AntennaType{static_cast<AntennaType>(std::get<0>(state))};
           })
      .doc() = "An antenna type";

  py::class_<MCAntenna>(m, "MCAntenna")
      .def_rw("atype",
              &MCAntenna::atype,
              ":class:`~pyarts.arts.AntennaType` The antenna type")
      .def_rw("sigma_aa",
              &MCAntenna::sigma_aa,
              ":class:`float` The azimuthal half-width")
      .def_rw("sigma_za",
              &MCAntenna::sigma_za,
              ":class:`float` The zenith half-width")
      .def_rw("aa_grid",
              &MCAntenna::aa_grid,
              ":class:`~pyarts.arts.Vector` The azimuth grid")
      .def_rw("za_grid",
              &MCAntenna::za_grid,
              ":class:`~pyarts.arts.Vector` The zenith grid")
      .def_rw("G_lookup",
              &MCAntenna::G_lookup,
              ":class:`~pyarts.arts.Matrix` The lookup")
      .def("__getstate__",
           [](const MCAntenna& self) {
             return std::make_tuple(self.atype,
                                    self.sigma_aa,
                                    self.sigma_za,
                                    self.aa_grid,
                                    self.za_grid,
                                    self.G_lookup);
           })
      .def("__setstate__",
           [](MCAntenna* m,
              const std::
                  tuple<AntennaType, Numeric, Numeric, Vector, Vector, Matrix>&
                      state) {
             new (m) MCAntenna{};
             m->atype    = std::get<0>(state);
             m->sigma_aa = std::get<1>(state);
             m->sigma_za = std::get<2>(state);
             m->aa_grid  = std::get<3>(state);
             m->za_grid  = std::get<4>(state);
             m->G_lookup = std::get<5>(state);
           });
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize mcantenna\n", e.what()));
}
}  // namespace Python
