#include <lbl.h>
#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <python_interface.h>

#include <memory>

#include "enums.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "lbl_lineshape_voigt_ecs.h"
#include "partfun.h"
#include "py_macros.h"
#include "quantum_numbers.h"

namespace Python {
void py_lbl(py::module_& m) try {
  auto lbl = m.def_submodule("lbl", "Line-by-line helper functions");

  artsclass<lbl::temperature::data>(m, "TemperatureModel")
      .def(py::init<LineShapeModelType, Vector>())
      .def_property_readonly(
          "type",
          &lbl::temperature::data::Type,
          ":class:`~pyarts.arts.options.TemperatureModelType` The type of the model")
      .def_property_readonly("data",
                             &lbl::temperature::data::X,
                             ":class:`~pyarts.arts.Vector` The coefficients")
      .PythonInterfaceBasicRepresentation(lbl::temperature::data);

  using pair_vector_type =
      std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>;
  artsarray<pair_vector_type>(m, "LineShapeVariableTemperatureModelList")
      .def("get",
           [](const pair_vector_type& self,
              LineShapeModelVariable x) -> lbl::temperature::data {
             for (auto& [var, data] : self) {
               if (var == x) {
                 return data;
               }
             }
             return {};
           })
      .def("set",
           [](pair_vector_type& self,
              LineShapeModelVariable x,
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
      .def("G0",
           py::vectorize(&lbl::line_shape::species_model::G0),
           "The G0 coefficient",
           py::arg("T0"),
           py::arg("T"),
           py::arg("P"))
      .def("Y",
           py::vectorize(&lbl::line_shape::species_model::Y),
           "The Y coefficient",
           py::arg("T0"),
           py::arg("T"),
           py::arg("P"))
      .def("D0",
           py::vectorize(&lbl::line_shape::species_model::D0),
           "The D0 coefficient",
           py::arg("T0"),
           py::arg("T"),
           py::arg("P"))
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
      .def("G0",
           &lbl::line_shape::model::G0,
           "The G0 coefficient",
           py::arg("atm"))
      .def("Y", &lbl::line_shape::model::Y, "The Y coefficient", py::arg("atm"))
      .def("D0",
           &lbl::line_shape::model::D0,
           "The D0 coefficient",
           py::arg("atm"))
      .def("remove",
           [](lbl::line_shape::model& self, LineShapeModelVariable x) {
             for (auto& specmod : self.single_models) {
               auto ptr = std::find_if(specmod.data.begin(),
                                       specmod.data.end(),
                                       [x](auto& y) { return y.first == x; });
               if (ptr != specmod.data.end()) specmod.data.erase(ptr);
             }
           })
      .def("remove_zeros",
           [](lbl::line_shape::model& self) { self.clear_zeroes(); })
      .PythonInterfaceBasicRepresentation(lbl::line_shape::model);

  artsclass<lbl::zeeman::model>(m, "ZeemanLineModel")
      .def_readwrite(
          "on",
          &lbl::zeeman::model::on,
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
      .def("s", py::vectorize(&lbl::line::s), "The line strength")
      .PythonInterfaceBasicRepresentation(lbl::line);

  artsarray<std::vector<lbl::line>>(m, "LineList")
      .PythonInterfaceBasicRepresentation(std::vector<lbl::line>);

  artsclass<lbl::band_data>(m, "AbsorptionBandData")
      .def_readwrite("lines", &lbl::band_data::lines, "The lines in the band")
      .def_readwrite(
          "lineshape", &lbl::band_data::lineshape, "The lineshape type")
      .def_readwrite("cutoff", &lbl::band_data::cutoff, "The cutoff type")
      .def_readwrite("cutoff_value",
                     &lbl::band_data::cutoff_value,
                     "The cutoff value [Hz]")
      .PythonInterfaceBasicRepresentation(lbl::band_data);

  artsclass<AbsorptionBand>(m, "AbsorptionBand")
      .def(py::init([]() { return std::make_shared<AbsorptionBand>(); }),
           "Default target")
      .PythonInterfaceCopyValue(AbsorptionBand)
      .PythonInterfaceWorkspaceVariableConversion(AbsorptionBand)
      .PythonInterfaceBasicRepresentation(AbsorptionBand)
      .PythonInterfaceFileIO(AbsorptionBand)
      .def_readwrite("data",
                     &AbsorptionBand::data,
                     ":class:`~pyarts.arts.AbsorptionBandData`")
      .def_readwrite("key",
                     &AbsorptionBand::key,
                     ":class:`~pyarts.arts.QuantumIdentifier`")
      .PythonInterfaceWorkspaceDocumentation(AbsorptionBand);

  artsarray<ArrayOfAbsorptionBand>(m, "ArrayOfAbsorptionBand")
      .def(py::init([]() { return std::make_shared<ArrayOfAbsorptionBand>(); }),
           "Default target")
      .PythonInterfaceCopyValue(ArrayOfAbsorptionBand)
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfAbsorptionBand)
      .PythonInterfaceBasicRepresentation(ArrayOfAbsorptionBand)
      .PythonInterfaceFileIO(ArrayOfAbsorptionBand)
      .def(
          "__getitem__",
          [](ArrayOfAbsorptionBand& x,
             const QuantumIdentifier& key) -> std::shared_ptr<lbl::band_data> {
            for (auto& v : x) {
              if (v.key == key) {
                return {&v.data, [](void*) {}};
              }
            }
            return std::make_shared<lbl::band_data>();
          },
          py::return_value_policy::reference_internal,
          py::keep_alive<0, 1>(),
          ":class:`~pyarts.arts.AbsorptionBandData`")
      .def(
          "__setitem__",
          [](ArrayOfAbsorptionBand& x,
             const QuantumIdentifier& key,
             const lbl::band_data& y) {
            for (auto& v : x) {
              if (v.key == key) {
                v.data = y;
                return;
              }
            }
            x.emplace_back(key, y);
          },
          py::return_value_policy::reference_internal,
          py::keep_alive<0, 1>(),
          ":class:`~pyarts.arts.AbsorptionBandData`")
      .PythonInterfaceWorkspaceDocumentation(ArrayOfAbsorptionBand);

  artsclass<LinemixingEcsData>(m, "LinemixingEcsData")
      .def(py::init([]() { return std::make_shared<LinemixingEcsData>(); }),
           "Default target")
      .PythonInterfaceCopyValue(LinemixingEcsData)
      .PythonInterfaceWorkspaceVariableConversion(LinemixingEcsData)
      .PythonInterfaceBasicRepresentation(LinemixingEcsData)
      .PythonInterfaceFileIO(LinemixingEcsData)
      .PythonInterfaceWorkspaceDocumentation(LinemixingEcsData);

  lbl.def(
      "equivalent_lines",
      [](const AbsorptionBand& band,
         const LinemixingEcsData& ecs_data,
         const AtmPoint& atm,
         const Vector& T) {
        lbl::voigt::ecs::ComputeData com_data({}, atm);

        const auto K = band.data.front().ls.one_by_one
                           ? band.data.front().ls.single_models.size()
                           : 1;
        const auto N = band.data.size();
        const auto M = T.size();

        auto eqv_str = ComplexTensor3(M, K, N);
        auto eqv_val = ComplexTensor3(M, K, N);

        equivalent_values(eqv_str,
                          eqv_val,
                          com_data,
                          band.key,
                          band.data,
                          ecs_data.data.at(band.key.Isotopologue()),
                          atm,
                          T);

        return std::pair{eqv_str, eqv_val};
      },
      "Compute equivalent lines for a given band",
      py::arg("band"),
      py::arg("ecs_data"),
      py::arg("atm"),
      py::arg("T"));

  lbl.def("Q",
          &PartitionFunctions::Q,
          "Partition function",
          py::arg("T"),
          py::arg("isotopologue"));
  lbl.def(
      "Q",
      [](const Vector& T, const SpeciesIsotopeRecord& isot) {
        Vector Q(T.size());
        std::transform(
            T.begin(), T.end(), Q.begin(), [&isot](Numeric t) -> Numeric {
              return PartitionFunctions::Q(t, isot);
            });
        return Q;
      },
      "Partition function",
      py::arg("T"),
      py::arg("isotopologue"));
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize lbl\n", e.what()));
}
}  // namespace Python