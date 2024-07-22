#include <lbl.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/pair.h>
#include <python_interface.h>

#include <iomanip>
#include <memory>
#include <stdexcept>

#include "enums.h"
#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "partfun.h"
#include "py_macros.h"
#include "quantum_numbers.h"

namespace Python {
void py_lbl(py::module_& m) try {
  auto lbl = m.def_submodule("lbl", "Line-by-line helper functions");

  py::class_<lbl::line_key> line_key(lbl, "line_key");
  line_key.def_rw("band", &lbl::line_key::band);
  line_key.def_rw("line", &lbl::line_key::line);
  line_key.def_rw("spec", &lbl::line_key::spec);
  line_key.def_rw("ls_var", &lbl::line_key::ls_var);
  line_key.def_rw("ls_coeff", &lbl::line_key::ls_coeff);
  line_key.def_rw("var", &lbl::line_key::var);
  str_interface(line_key);

  py::class_<lbl::temperature::data>(m, "TemperatureModel")
      .def(py::init<LineShapeModelType, Vector>(),
           "type"_a,
           "data"_a = Vector{0.0})
      .def_prop_rw(
          "type",
          &lbl::temperature::data::Type,
          [](lbl::temperature::data& self, LineShapeModelType x) {
            self = lbl::temperature::data{x, self.X()};
          },
          ":class:`~pyarts.arts.options.TemperatureModelType` The type of the model")
      .def_prop_rw(
          "data",
          [](lbl::temperature::data& self) { return self.X(); },
          [](lbl::temperature::data& self, const Vector& x) {
            self = lbl::temperature::data{self.Type(), x};
          },
          ":class:`~pyarts.arts.Vector` The coefficients")
      .PythonInterfaceBasicRepresentation(lbl::temperature::data)
      .def("__getstate__",
           [](const lbl::temperature::data& v) {
             return std::tuple<LineShapeModelType, Vector>{v.Type(), v.X()};
           })
      .def("__setstate__",
           [](lbl::temperature::data& v,
              const std::tuple<LineShapeModelType, Vector>& x) {
             new (&v) lbl::temperature::data{std::get<0>(x), std::get<1>(x)};
           });

  using pair_vector_type =
      std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>;
  auto lsvtml =
      py::bind_vector<pair_vector_type, py::rv_policy::reference_internal>(
          m, "LineShapeVariableTemperatureModelList")
          .def("__getitem__",
               [](pair_vector_type& self,
                  LineShapeModelVariable x) -> lbl::temperature::data& {
                 for (auto& [var, data] : self) {
                   if (var == x) {
                     return data;
                   }
                 }
                 throw std::out_of_range(var_string('"', x, '"'));
               })
          .def("__setitem__",
               [](pair_vector_type& self,
                  LineShapeModelVariable x,
                  const lbl::temperature::data& y) {
                 for (auto& [var, data] : self) {
                   if (var == x) {
                     data = y;
                     return;
                   }
                 }
                 self.emplace_back(x, y);
               })
          .PythonInterfaceBasicRepresentation(pair_vector_type);
  vector_interface(lsvtml);

  py::class_<lbl::line_shape::species_model>(m, "LineShapeSpeciesModel")
      .def_rw(
          "species", &lbl::line_shape::species_model::species, "The species")
      .def_rw("data", &lbl::line_shape::species_model::data, "The data")
      .def("__getitem__",
           [](py::object& x, const py::object& key) {
             return x.attr("data").attr("__getitem__")(key);
           })
      .def("__setitem__",
           [](py::object& x, const py::object& key, const py::object& val) {
             x.attr("data").attr("__setitem__")(key, val);
           })
      .def(
          "G0",
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.G0(t0, t, p);
                },
                T0,
                T,
                P);
          },
          "The G0 coefficient",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "Y",
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.Y(t0, t, p);
                },
                T0,
                T,
                P);
          },
          "The Y coefficient",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "D0",
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.D0(t0, t, p);
                },
                T0,
                T,
                P);
          },
          "The D0 coefficient",
          "T0"_a,
          "T"_a,
          "P"_a)
      .PythonInterfaceBasicRepresentation(lbl::line_shape::species_model);

  using line_shape_model_list = std::vector<lbl::line_shape::species_model>;
  auto lsml =
      py::bind_vector<line_shape_model_list, py::rv_policy::reference_internal>(
          m, "LineShapeModelList")
          .PythonInterfaceBasicRepresentation(line_shape_model_list);
  vector_interface(lsml);

  py::class_<lbl::line_shape::model>(m, "LineShapeModelFIXMENAMEODR")
      .def_rw("one_by_one",
              &lbl::line_shape::model::one_by_one,
              "If true, the lines are treated one by one")
      .def_rw("T0", &lbl::line_shape::model::T0, "The reference temperature")
      .def_rw("single_models",
              &lbl::line_shape::model::single_models,
              "The single models")
      .def("G0", &lbl::line_shape::model::G0, "The G0 coefficient", "atm"_a)
      .def("Y", &lbl::line_shape::model::Y, "The Y coefficient", "atm"_a)
      .def("D0", &lbl::line_shape::model::D0, "The D0 coefficient", "atm"_a)
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

  py::class_<lbl::zeeman::model>(m, "ZeemanLineModel")
      .def_rw(
          "on",
          &lbl::zeeman::model::on,
          ":class:`~pyarts.arts.Bool` If True, the Zeeman effect is included")
      .def_prop_rw(
          "gl",
          [](lbl::zeeman::model& z) { return z.gl(); },
          [](lbl::zeeman::model& z, Numeric g) { z.gl(g); },
          ":class:`~pyarts.arts.Numeric` The lower level statistical weight")
      .def_prop_rw(
          "gu",
          [](lbl::zeeman::model& z) { return z.gu(); },
          [](lbl::zeeman::model& z, Numeric g) { z.gu(g); },
          ":class:`~pyarts.arts.Numeric` The upper level statistical weight")
      .PythonInterfaceBasicRepresentation(lbl::zeeman::model);

  py::class_<lbl::line>(m, "AbsorptionLine")
      .def_rw("a",
              &lbl::line::a,
              ":class:`~pyarts.arts.Numeric` The Einstein coefficient")
      .def_rw("f0",
              &lbl::line::f0,
              ":class:`~pyarts.arts.Numeric` The line center frequency [Hz]")
      .def_rw("e0",
              &lbl::line::e0,
              ":class:`~pyarts.arts.Numeric` The lower level energy [J]")
      .def_rw(
          "gu",
          &lbl::line::gu,
          ":class:`~pyarts.arts.Numeric` The upper level statistical weight")
      .def_rw(
          "gl",
          &lbl::line::gl,
          ":class:`~pyarts.arts.Numeric` The lower level statistical weight")
      .def_rw("z", &lbl::line::z, "The Zeeman model")
      .def_rw("ls", &lbl::line::ls, "The line shape model")
      .def_rw("qn", &lbl::line::qn, "The quantum numbers of this line")
      .def(
          "s",
          [](const lbl::line& self, py::object& T, py::object& Q) {
            return vectorize(
                [&self](Numeric t, Numeric q) { return self.s(t, q); }, T, Q);
          },
          "T"_a,
          "Q"_a,
          "The line strength")
      .PythonInterfaceBasicRepresentation(lbl::line);

  auto ll = py::bind_vector<std::vector<lbl::line>,
                            py::rv_policy::reference_internal>(m, "LineList")
                .PythonInterfaceBasicRepresentation(std::vector<lbl::line>);
  vector_interface(ll);

  py::class_<lbl::band_data>(m, "AbsorptionBandData")
      .def_rw("lines", &lbl::band_data::lines, "The lines in the band")
      .def_rw("lineshape", &lbl::band_data::lineshape, "The lineshape type")
      .def_rw("cutoff", &lbl::band_data::cutoff, "The cutoff type")
      .def_rw("cutoff_value",
              &lbl::band_data::cutoff_value,
              "The cutoff value [Hz]")
      .PythonInterfaceBasicRepresentation(lbl::band_data);

  py::class_<AbsorptionBand> absd(m, "AbsorptionBand");
  workspace_group_interface(absd);
  absd.def_rw("data",
              &AbsorptionBand::data,
              ":class:`~pyarts.arts.AbsorptionBandData`")
      .def_rw("key",
              &AbsorptionBand::key,
              ":class:`~pyarts.arts.QuantumIdentifier`");

  auto aoab =
      py::bind_vector<ArrayOfAbsorptionBand, py::rv_policy::reference_internal>(
          m, "ArrayOfAbsorptionBand");
  workspace_group_interface(aoab);
  vector_interface(aoab);
  aoab.def(
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
          py::rv_policy::reference_internal,
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
          py::rv_policy::reference_internal,
          py::keep_alive<0, 1>(),
          ":class:`~pyarts.arts.AbsorptionBandData`");
  aoab.def("__getitem__",
           [](const ArrayOfAbsorptionBand& x,
              const lbl::line_key& key) -> Numeric { return key.get_value(x); })
      .def("__setitem__",
           [](ArrayOfAbsorptionBand& x, const lbl::line_key& key, Numeric v) {
             key.get_value(x) = v;
           });

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
      "band"_a,
      "ecs_data"_a,
      "atm"_a,
      "T"_a);

  lbl.def(
      "Q",
      [](py::object t, SpeciesIsotope isot) {
        return vectorize(
            [isot](Numeric T) { return PartitionFunctions::Q(T, isot); }, t);
      },
      "Partition function",
      "T"_a,
      "isotopologue"_a);

  py::class_<LinemixingEcsData> led(m, "LinemixingEcsData");
  workspace_group_interface(led);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize lbl\n", e.what()));
}
}  // namespace Python
