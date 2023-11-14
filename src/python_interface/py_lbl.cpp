#include <python_interface.h>

#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "lbl_temperature_model.h"
#include "py_macros.h"

#include <lbl.h>

namespace Python {
void py_lbl(py::module_& m) try {
  artsclass<lbl::temperature::data>(m, "TemperatureModel");

  artsclass<lbl::line_shape::species_model>(m, "LineShapeSpeciesModel")
  .def_readwrite("species", &lbl::line_shape::species_model::species, "The species")
  .def_readwrite("data", &lbl::line_shape::species_model::data, "The data");

  artsclass<lbl::line_shape::model>(m, "LineShapeModelFIXMENAMEODR")
  .def_readwrite("one_by_one", &lbl::line_shape::model::one_by_one, "If true, the lines are treated one by one")
  .def_readwrite("T0", &lbl::line_shape::model::T0, "The reference temperature")
  .def_readwrite("single_models", &lbl::line_shape::model::single_models, "The single models");

  artsclass<lbl::zeeman::model>(m, "ZeemanLineModel")
  .def_property("gl", &lbl::zeeman::model::gl, &lbl::zeeman::model::gl, ":class:`~pyarts.arts.Numeric` The lower level statistical weight")
  .def_property("gu", &lbl::zeeman::model::gu, &lbl::zeeman::model::gu, ":class:`~pyarts.arts.Numeric` The upper level statistical weight");

  artsclass<lbl::line>(m, "AbsorptionLine")
    .def_readwrite("a", &lbl::line::a, ":class:`~pyarts.arts.Numeric` The Einstein coefficient")
    .def_readwrite("f0", &lbl::line::f0, ":class:`~pyarts.arts.Numeric` The line center frequency [Hz]")
    .def_readwrite("e0", &lbl::line::e0, ":class:`~pyarts.arts.Numeric` The lower level energy [J]")
    .def_readwrite("gu", &lbl::line::gu, ":class:`~pyarts.arts.Numeric` The upper level statistical weight")
    .def_readwrite("gl", &lbl::line::gl, ":class:`~pyarts.arts.Numeric` The lower level statistical weight")
    .def_readwrite("z", &lbl::line::z, "The Zeeman model")
    .def_readwrite("ls", &lbl::line::ls, "The line shape model")
    .def_readwrite("qn", &lbl::line::qn, "The quantum numbers of this line");

  artsclass<lbl::band_data>(m, "AbsorptionBandData")
  .def_readwrite("lines", &lbl::band_data::lines, "The lines in the band")
  .def_readwrite("lineshape", &lbl::band_data::lineshape, "The lineshape type")
  .def_readwrite("linestrength", &lbl::band_data::linestrength, "The linestrength type")
  .def_readwrite("cutoff", &lbl::band_data::cutoff, "The cutoff type")
  .def_readwrite("cutoff_value", &lbl::band_data::cutoff_value, "The cutoff value [Hz]");

  artsclass<AbsorptionBand>(m, "AbsorptionBand")
      .def(py::init([]() { return std::make_shared<AbsorptionBand>(); }), "Default target")
      .def_readwrite("data", &AbsorptionBand::data, ":class:`~pyarts.arts.AbsorptionBandData`")
      .def_readwrite("key", &AbsorptionBand::data, ":class:`~pyarts.arts.QuantumIdentifier`")
      .PythonInterfaceBasicRepresentation(AbsorptionBand);

  artsarray<AbsorptionBands>(m, "AbsorptionBands");
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize lbl\n", e.what()));
}
}  // namespace Python