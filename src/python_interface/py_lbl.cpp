#include <lbl.h>
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>

#include <iomanip>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "partfun.h"
#include "quantum.h"

namespace Python {
void py_lbl(py::module_& m) try {
  auto lbl = m.def_submodule("lbl", "Line-by-line helper functions");

  py::class_<lbl::line_key> line_key(lbl, "line_key");
  line_key.def_rw("band", &lbl::line_key::band, "The band");
  line_key.def_rw("line", &lbl::line_key::line, "The line");
  line_key.def_rw("spec", &lbl::line_key::spec, "The species");
  line_key.def_rw("ls_var", &lbl::line_key::ls_var, "The line shape variable");
  line_key.def_rw(
      "ls_coeff", &lbl::line_key::ls_coeff, "The line shape coefficient");
  line_key.def_rw("var", &lbl::line_key::var, "The variable");
  line_key.doc() = "A key for a line";
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
          ":class:`~pyarts3.arts.TemperatureModelType` The type of the model")
      .def_prop_rw(
          "data",
          [](lbl::temperature::data& self) { return self.X(); },
          [](lbl::temperature::data& self, const Vector& x) {
            self = lbl::temperature::data{self.Type(), x};
          },
          ":class:`~pyarts3.arts.Vector` The coefficients")
      .def("__getstate__",
           [](const lbl::temperature::data& v) {
             return std::tuple<LineShapeModelType, Vector>{v.Type(), v.X()};
           })
      .def("__setstate__",
           [](lbl::temperature::data& v,
              const std::tuple<LineShapeModelType, Vector>& x) {
             new (&v) lbl::temperature::data{std::get<0>(x), std::get<1>(x)};
           })
      .doc() = "Temperature model";

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
                 throw std::out_of_range(std::format("\"{}\"", x));
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
               });
  lsvtml.doc() = "A list of line shape models with temperature coefficients";
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
      .doc() = "Line shape model for a species";

  using line_shape_model_list = std::vector<lbl::line_shape::species_model>;
  auto lsml =
      py::bind_vector<line_shape_model_list, py::rv_policy::reference_internal>(
          m, "LineShapeModelList");
  lsml.doc() = "A list of line shape models";
  vector_interface(lsml);

  py::class_<lbl::line_shape::model>(m, "LineShapeModel")
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
      .def(
          "remove",
          [](lbl::line_shape::model& self, LineShapeModelVariable x) {
            for (auto& specmod : self.single_models) {
              auto ptr = std::find_if(specmod.data.begin(),
                                      specmod.data.end(),
                                      [x](auto& y) { return y.first == x; });
              if (ptr != specmod.data.end()) specmod.data.erase(ptr);
            }
          },
          "x"_a,
          "Remove a variable")
      .def(
          "remove_zeros",
          [](lbl::line_shape::model& self) { self.clear_zeroes(); },
          "Remove zero coefficients")
      .doc() = "Line shape model";

  py::class_<lbl::zeeman::model>(m, "ZeemanLineModel")
      .def_rw(
          "on",
          &lbl::zeeman::model::on,
          ":class:`~pyarts3.arts.Bool` If True, the Zeeman effect is included")
      .def_prop_rw(
          "gl",
          [](lbl::zeeman::model& z) { return z.gl(); },
          [](lbl::zeeman::model& z, Numeric g) { z.gl(g); },
          ":class:`~pyarts3.arts.Numeric` The lower level statistical weight")
      .def_prop_rw(
          "gu",
          [](lbl::zeeman::model& z) { return z.gu(); },
          [](lbl::zeeman::model& z, Numeric g) { z.gu(g); },
          ":class:`~pyarts3.arts.Numeric` The upper level statistical weight")
      .def(
          "strengths",
          [](const lbl::zeeman::model& mod, const QuantumState& qn) {
            std::map<std::string, std::vector<double>> out;

            const Index Npi = mod.size(qn, lbl::zeeman::pol::pi);
            for (Index i = 0; i < Npi; i++) {
              out["pi"].push_back(mod.Strength(qn, lbl::zeeman::pol::pi, i));
            }

            const Index Nsp = mod.size(qn, lbl::zeeman::pol::sp);
            for (Index i = 0; i < Nsp; i++) {
              out["sp"].push_back(mod.Strength(qn, lbl::zeeman::pol::sp, i));
            }

            const Index Nsm = mod.size(qn, lbl::zeeman::pol::sm);
            for (Index i = 0; i < Nsm; i++) {
              out["sm"].push_back(mod.Strength(qn, lbl::zeeman::pol::sm, i));
            }

            return out;
          },
          "The number of Zeeman lines")
      .doc() = "Zeeman model";

  py::class_<lbl::line>(m, "AbsorptionLine")
      .def_rw("a",
              &lbl::line::a,
              ":class:`~pyarts3.arts.Numeric` The Einstein coefficient [1 / s]")
      .def_rw("f0",
              &lbl::line::f0,
              ":class:`~pyarts3.arts.Numeric` The line center frequency [Hz]")
      .def_rw("e0",
              &lbl::line::e0,
              ":class:`~pyarts3.arts.Numeric` The lower level energy [J]")
      .def_rw(
          "gu",
          &lbl::line::gu,
          ":class:`~pyarts3.arts.Numeric` The upper level statistical weight [-]")
      .def_rw(
          "gl",
          &lbl::line::gl,
          ":class:`~pyarts3.arts.Numeric` The lower level statistical weight [-]")
      .def_rw("z",
              &lbl::line::z,
              ":class:`~pyarts3.arts.ZeemanLineModel` The Zeeman model")
      .def_rw("ls",
              &lbl::line::ls,
              ":class:`~pyarts3.arts.LineShapeModel` The line shape model")
      .def_rw(
          "qn",
          &lbl::line::qn,
          ":class:`~pyarts3.arts.QuantumNumberLocalState` The local quantum numbers of this line")
      .def(
          "s",
          [](const lbl::line& self, py::object& T, py::object& Q) {
            return vectorize(
                [&self](Numeric t, Numeric q) { return self.s(t, q); }, T, Q);
          },
          "T"_a,
          "Q"_a,
          "The line strength")
      .def(
          "hitran_s",
          [](const lbl::line& self, const SpeciesIsotope& isot, Numeric T0) {
            return self.hitran_s(isot, T0);
          },
          "isot"_a,
          "T0"_a = 296.0,
          "The HITRAN-like line strength")
      .doc() = "A single absorption line";

  auto ll  = py::bind_vector<std::vector<lbl::line>,
                             py::rv_policy::reference_internal>(m, "LineList");
  ll.doc() = "A list of absorption lines";
  vector_interface(ll);

  py::class_<AbsorptionBand> ab(m, "AbsorptionBand");
  ab.def_rw("lines", &AbsorptionBand::lines, "The lines in the band")
      .def_rw("lineshape", &AbsorptionBand::lineshape, "The lineshape type")
      .def_rw("cutoff", &AbsorptionBand::cutoff, "The cutoff type")
      .def_rw("cutoff_value",
              &AbsorptionBand::cutoff_value,
              "The cutoff value [Hz]")
      .def(
          "keep_frequencies",
          [](AbsorptionBand& band, Vector2 freqs) {
            band.sort();
            auto l = band.active_lines(freqs[0], freqs[1]).second;
            std::vector<lbl::line> new_lines(l.begin(), l.end());
            band.lines = std::move(new_lines);
          },
          "freqs"_a,
          "Keep only the lines within the given frequency range")
      .def(
          "keep_hitran_s",
          [](AbsorptionBand& band,
             Numeric min_s,
             const SpeciesIsotope& isot,
             Numeric T0) {
            std::erase_if(band.lines, [&isot, &T0, &min_s](auto& line) {
              return line.hitran_s(isot, T0) < min_s;
            });
          },
          "min_s"_a,
          "isot"_a,
          "T0"_a = 296.0,
          "Keep only the lines with a stronger HITRAN-like line strength");
  generic_interface(ab);

  auto aoab = py::bind_map<AbsorptionBands, py::rv_policy::reference_internal>(
      m, "AbsorptionBands");
  generic_interface(aoab);
  aoab.def("__getitem__",
           [](const AbsorptionBands& x, const lbl::line_key& key) -> Numeric {
             return key.get_value(x);
           })
      .def("__setitem__",
           [](AbsorptionBands& x, const lbl::line_key& key, Numeric v) {
             key.get_value(x) = v;
           });
  aoab.def(
      "merge",
      [](AbsorptionBands& self,
         const AbsorptionBands& other) -> std::pair<Size, Size> {
        Size added = 0, updated = 0;

        for (const auto& [key, band] : other) {
          auto [it, newband] = self.try_emplace(key, band);

          if (newband) {
            added += band.size();
          } else {
            for (auto& line : band.lines) {
              if (it->second.merge(line)) {
                added++;
              } else {
                updated++;
              }
            }
          }
        }
        return {added, updated};
      },
      R"(Merge the other absorption bands into this

If the key in the other absorption bands already exists in this, the lines with the same local
quantum numbers are overwritten by those of the other absorption bands.

Parameters
----------
other : AbsorptionBands
    The other absorption bands to merge into this
)",
      "other"_a);
  aoab.def(
      "clear_linemixing",
      [](AbsorptionBands& self) {
        Size sum = 0;
        for (auto& [_, band] : self) {
          for (auto& line : band.lines) {
            for (auto& lsm : line.ls.single_models) {
              sum += lsm.remove_variables<LineShapeModelVariable::Y,
                                          LineShapeModelVariable::G,
                                          LineShapeModelVariable::DV>();
            }
          }
        }
        return sum;
      },
      R"(Clear the linemixing data from all bands by removing it from the inner line shape models

Returns
-------
int : The number of removed variables
)");

  aoab.def(
      "count_lines",
      [](const AbsorptionBands& self, SpeciesEnum spec) {
        Size n = 0;

        if (spec == SpeciesEnum::Bath) {
          for (auto& band : self | std::views::values) {
            n += band.size();
          }
          return n;
        }
        for (auto& [key, band] : self) {
          if (spec == key.isot.spec) n += band.size();
        }
        return n;
      },
      "spec"_a = SpeciesEnum::Bath,
      "Return the total number of lines");

  aoab.def(
      "remove_hitran_s",
      &lbl::keep_hitran_s,
      R"(Removes all lines with a weaker HITRAN-like line strength than those provided by the remove map.

Parameters
----------
remove : dict
    The species to keep with their respective minimum HITRAN-like line strengths
T0 : float
    The reference temperature. Defaults to 296.0.
)",
      "remove"_a,
      "T0"_a = 296.0);

  aoab.def(
      "percentile_hitran_s",
      [](const AbsorptionBands& self,
         const std::variant<Numeric, std::unordered_map<SpeciesEnum, Numeric>>&
             percentile,
         const Numeric T0) {
        return std::visit(
            [&](auto& i) { return lbl::percentile_hitran_s(self, i, T0); },
            percentile);
      },
      R"(Map of HITRAN linestrengths at a given percentile

.. note::

  The percentile is approximated by floating point arithmetic on the sorted HITRAN line strenght values.

Parameters
----------
percentile : float or dict
    The percentile to keep. If a float, the same percentile is used for all species. If a dict, the species are mapped to their respective percentiles.  Values must be [0, 100].
T0 : float
    The reference temperature. Defaults to 296.0.
)",
      "approximate_percentile"_a,
      "T0"_a = 296.0);

  aoab.def(
      "keep_hitran_s",
      [](AbsorptionBands& self,
         const std::variant<Numeric, std::unordered_map<SpeciesEnum, Numeric>>&
             percentile,
         const Numeric T0) {
        lbl::keep_hitran_s(
            self,
            std::visit(
                [&](auto& i) { return lbl::percentile_hitran_s(self, i, T0); },
                percentile),
            T0);
      },
      R"(Wraps calling percentile_hitran_s followed by remove_hitran_s.

Parameters
----------
percentile : float or dict
    See percentile_hitran_s.
T0 : float
    The reference temperature. Defaults to 296.0.
)",
      "approximate_percentile"_a,
      "T0"_a = 296.0);

  py::class_<LinemixingSingleEcsData> ed(m, "LinemixingSingleEcsData");
  ed.def_rw("scaling",
            &LinemixingSingleEcsData::scaling,
            ":class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("beta",
            &LinemixingSingleEcsData::beta,
            ":class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("lambda_",  // Fix name, not python
            &LinemixingSingleEcsData::lambda,
            ":class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("collisional_distance",
            &LinemixingSingleEcsData::collisional_distance,
            ":class:`~pyarts3.arts.TemperatureModel`");
  ed.def("Q",
         &LinemixingSingleEcsData::Q,
         "J"_a,
         "T"_a,
         "T0"_a,
         "energy"_a,
         R"(The Q coefficient for the ECS model)");
  ed.def("Omega",
         &LinemixingSingleEcsData::Omega,
         "T"_a,
         "T0"_a,
         "mass"_a,
         "other_mass"_a,
         "energy_x"_a,
         "energy_xm2"_a,
         R"(The Omega coefficient for the ECS model)");
  generic_interface(ed);

  auto lsed =
      py::bind_map<LinemixingSpeciesEcsData>(m, "LinemixingSpeciesEcsData");
  generic_interface(lsed);

  auto led = py::bind_map<LinemixingEcsData>(m, "LinemixingEcsData");
  generic_interface(led);

  lbl.def(
      "equivalent_lines",
      [](const AbsorptionBand& band,
         const QuantumIdentifier& qid,
         const LinemixingEcsData& ecs_data,
         const AtmPoint& atm,
         const Vector& T) {
        lbl::voigt::ecs::ComputeData com_data({}, atm);

        const auto K = band.front().ls.one_by_one
                           ? band.front().ls.single_models.size()
                           : 1;
        const auto N = band.size();
        const auto M = T.size();

        auto eqv_str = ComplexTensor3(M, K, N);
        auto eqv_val = ComplexTensor3(M, K, N);

        equivalent_values(eqv_str,
                          eqv_val,
                          com_data,
                          qid,
                          band,
                          ecs_data.at(qid.isot),
                          atm,
                          T);

        return std::pair{eqv_str, eqv_val};
      },
      "Compute equivalent lines for a given band",
      "band"_a,
      "qid"_a,
      "ecs_data"_a,
      "atm"_a,
      "T"_a);

  aoab.def(
      "propagation_matrix",
      [](const AbsorptionBands& self,
         const AscendingGrid& f,
         const AtmPoint& atm,
         const SpeciesEnum& spec,
         const PropagationPathPoint& path_point,
         const LinemixingEcsData& ecs_data,
         const Index& no_negative_absorption,
         const py::kwargs&) {
        PropmatVector propagation_matrix(f.size());
        StokvecVector nlte_vector(f.size());
        PropmatMatrix propagation_matrix_jacobian(0, f.size());
        StokvecMatrix nlte_matrix(0, f.size());
        JacobianTargets jacobian_targets{};

        propagation_matrixAddLines(propagation_matrix,
                                   nlte_vector,
                                   propagation_matrix_jacobian,
                                   nlte_matrix,
                                   f,
                                   jacobian_targets,
                                   spec,
                                   self,
                                   ecs_data,
                                   atm,
                                   path_point,
                                   no_negative_absorption);

        return propagation_matrix;
      },
      "f"_a,
      "atm"_a,
      "spec"_a                   = SpeciesEnum::Bath,
      "path_point"_a             = PropagationPathPoint{},
      "ecs_data"_a               = LinemixingEcsData{},
      "no_negative_absorption"_a = Index{1},
      "kwargs"_a                 = py::kwargs{},
      R"--(Computes the line-by-line model absorption in 1/m

The method accepts any number of kwargs to be compatible
with similar methods for computing the propagation matrix.

Parameters
----------
f : AscendingGrid
    Frequency grid [Hz]
atm : AtmPoint
    Atmospheric point
spec : SpeciesEnum, optional
    Species to use.  Defaults to all species.
path_point : PropagationPathPoint, optional
    The path point.  Default is POS [0, 0, 0], LOS [0, 0].
ecs_data : LinemixingEcsData, optional
    The ECS data.  Default is empty.
no_negative_absorption : Index, optional
    If 1, the absorption is set to zero if it is negative. The default is 1.

Returns
-------
propagation_matrix : PropmatVector
    Propagation matrix by frequency [1/m]

)--");

  py::class_<PartitionFunctionsData> partfun(m, "PartitionFunctionsData");
  generic_interface(partfun);
  partfun.doc() =
      "Data for partition functions, used in the line-by-line model";
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize lbl\n{}", e.what()));
}
}  // namespace Python
