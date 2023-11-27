#include <lbl.h>
#include <python_interface.h>

#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "lbl_temperature_model.h"
#include "py_macros.h"

namespace Python {
void py_lbl(py::module_& m) try {
  artsclass<lbl::temperature::data>(m, "TemperatureModel")
      .def(py::init<lbl::temperature::model_type, Vector>())
      .def_property_readonly(
          "type",
          &lbl::temperature::data::Type,
          ":class:`~pyarts.arts.options.TemperatureModelType` The type of the model")
      .def_property_readonly("data",
                             &lbl::temperature::data::X,
                             ":class:`~pyarts.arts.Vector` The coefficients")
      .PythonInterfaceBasicRepresentation(lbl::temperature::data);

  using pair_vector_type =
      std::vector<std::pair<lbl::line_shape::variable, lbl::temperature::data>>;
  artsarray<pair_vector_type>(m, "LineShapeVariableTemperatureModelList")
      .def("get",
           [](const pair_vector_type& self,
              lbl::line_shape::variable x) -> lbl::temperature::data {
             for (auto& [var, data] : self) {
               if (var == x) {
                 return data;
               }
             }
             return {};
           })
      .def("set",
           [](pair_vector_type& self,
              lbl::line_shape::variable x,
              lbl::temperature::data y) {
             for (auto& [var, data] : self) {
               if (var == x) {
                 data = y;
                 return;
               }
             }
             self.emplace_back(x, y);
           })
      .PythonInterfaceBasicRepresentation(pair_vector_type);

  artsclass<lbl::line_shape::species_model>(m, "LineShapeSpeciesModel")
      .def_readwrite(
          "species", &lbl::line_shape::species_model::species, "The species")
      .def_readwrite("data", &lbl::line_shape::species_model::data, "The data")
      .PythonInterfaceBasicRepresentation(lbl::line_shape::species_model);

  artsarray<std::vector<lbl::line_shape::species_model>>(m,
                                                         "LineShapeModelList")
      .PythonInterfaceBasicRepresentation(
          std::vector<lbl::line_shape::species_model>);

  artsclass<lbl::line_shape::model>(m, "LineShapeModelFIXMENAMEODR")
      .def_readwrite("one_by_one",
                     &lbl::line_shape::model::one_by_one,
                     "If true, the lines are treated one by one")
      .def_readwrite(
          "T0", &lbl::line_shape::model::T0, "The reference temperature")
      .def_readwrite("single_models",
                     &lbl::line_shape::model::single_models,
                     "The single models")
      .PythonInterfaceBasicRepresentation(lbl::line_shape::model);

  artsclass<lbl::zeeman::model>(m, "ZeemanLineModel")
      .def_property(
          "on",
          [](const lbl::zeeman::model& z) { return z.active(); },
          [](lbl::zeeman::model& z, bool x) { z.active(x); },
          ":class:`~pyarts.arts.Bool` If True, the Zeeman effect is included")
      .def_property(
          "gl",
          &lbl::zeeman::model::gl,
          &lbl::zeeman::model::gl,
          ":class:`~pyarts.arts.Numeric` The lower level statistical weight")
      .def_property(
          "gu",
          &lbl::zeeman::model::gu,
          &lbl::zeeman::model::gu,
          ":class:`~pyarts.arts.Numeric` The upper level statistical weight")
      .PythonInterfaceBasicRepresentation(lbl::zeeman::model);

  artsclass<lbl::line>(m, "AbsorptionLine")
      .def_readwrite("a",
                     &lbl::line::a,
                     ":class:`~pyarts.arts.Numeric` The Einstein coefficient")
      .def_readwrite(
          "f0",
          &lbl::line::f0,
          ":class:`~pyarts.arts.Numeric` The line center frequency [Hz]")
      .def_readwrite("e0",
                     &lbl::line::e0,
                     ":class:`~pyarts.arts.Numeric` The lower level energy [J]")
      .def_readwrite(
          "gu",
          &lbl::line::gu,
          ":class:`~pyarts.arts.Numeric` The upper level statistical weight")
      .def_readwrite(
          "gl",
          &lbl::line::gl,
          ":class:`~pyarts.arts.Numeric` The lower level statistical weight")
      .def_readwrite("z", &lbl::line::z, "The Zeeman model")
      .def_readwrite("ls", &lbl::line::ls, "The line shape model")
      .def_readwrite("qn", &lbl::line::qn, "The quantum numbers of this line")
      .PythonInterfaceBasicRepresentation(lbl::line);

  artsarray<std::vector<lbl::line>>(m, "LineList")
      .PythonInterfaceBasicRepresentation(std::vector<lbl::line>);

  artsclass<lbl::band_data>(m, "AbsorptionBandData")
      .def_readwrite("lines", &lbl::band_data::lines, "The lines in the band")
      .def_readwrite(
          "lineshape", &lbl::band_data::lineshape, "The lineshape type")
      .def_readwrite("linestrength",
                     &lbl::band_data::linestrength,
                     "The linestrength type")
      .def_readwrite("cutoff", &lbl::band_data::cutoff, "The cutoff type")
      .def_readwrite("cutoff_value",
                     &lbl::band_data::cutoff_value,
                     "The cutoff value [Hz]")
      .PythonInterfaceBasicRepresentation(lbl::band_data);

  artsclass<AbsorptionBand>(m, "AbsorptionBand")
      .def(py::init([]() { return std::make_shared<AbsorptionBand>(); }),
           "Default target")
      .def_readwrite("data",
                     &AbsorptionBand::data,
                     ":class:`~pyarts.arts.AbsorptionBandData`")
      .def_readwrite("key",
                     &AbsorptionBand::data,
                     ":class:`~pyarts.arts.QuantumIdentifier`")
      .PythonInterfaceBasicRepresentation(AbsorptionBand);

  artsarray<AbsorptionBands>(m, "AbsorptionBands");
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize lbl\n", e.what()));
}
}  // namespace Python