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
          }).def("__setstate__",
          [](AntennaType*a, const std::tuple<Index>& state) {
            new (a) AntennaType{static_cast<AntennaType>(std::get<0>(state))};
          }).doc() = "An antenna type";

  py::class_<MCAntenna>(m,"MCAntenna")
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
      .def(py::pickle(
          [](const MCAntenna& self) {
            return py::make_tuple(self.atype,
                                  self.sigma_aa,
                                  self.sigma_za,
                                  self.aa_grid,
                                  self.za_grid,
                                  self.G_lookup);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 6, "Invalid state!")

            auto out = std::make_shared<MCAntenna>();
            out->atype = t[0].cast<AntennaType>();
            out->sigma_aa = t[1].cast<Numeric>();
            out->sigma_za = t[2].cast<Numeric>();
            out->aa_grid = t[3].cast<Vector>();
            out->za_grid = t[4].cast<Vector>();
            out->G_lookup = t[5].cast<Matrix>();
            return out;
          }));
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize mcantenna\n", e.what()));
}
}  // namespace Python
