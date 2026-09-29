#include <auto_wsm.h>
#include <enumsLineShapeModelVariable.h>
#include <enumsSpeciesEnum.h>
#include <hpy_arts.h>
#include <hpy_numpy.h>
#include <hpy_vector.h>
#include <isotopologues.h>
#include <lbl.h>
#include <lbl_data.h>
#include <lbl_lineshape_model.h>
#include <lbl_lineshape_voigt_ecs.h>
#include <lbl_lineshape_voigt_ecs_hadded.h>
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <partfun.h>
#include <partition_function_data.h>
#include <python_interface.h>
#include <quantum.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace Python {
void py_lbl(py::module_& m) try {
  auto lbl = m.def_submodule("lbl", "Line-by-line helper functions");

  auto lssmm = py::bind_map<lbl::line_shape::species_model::map_t, py::rv_policy::reference_internal>(
      lbl, "LineShapeSpeciesModelMap");
  lssmm.doc() = "A map from model variable to line shape models";
  generic_interface(lssmm);

  auto lsmm  = py::bind_map<lbl::line_shape::model::map_t, py::rv_policy::reference_internal>(lbl, "LineShapeModelMap");
  lsmm.doc() = "A map from species to species line shape models";
  generic_interface(lsmm);

  py::class_<lbl::line_key> line_key(lbl, "line_key");
  generic_interface(line_key);
  line_key.def_rw("band", &lbl::line_key::band, "The band\n\n.. :class:`~pyarts3.arts.QuantumIdentifier`");
  line_key.def_rw("line", &lbl::line_key::line, "The line\n\n.. :class:`int`");
  line_key.def_rw("spec", &lbl::line_key::spec, "The species\n\n.. :class:`int`");
  line_key.def_rw(
      "ls_var", &lbl::line_key::ls_var, "The line shape variable\n\n.. :class:`~pyarts3.arts.LineShapeModelVariable`");
  line_key.def_rw("ls_coeff",
                  &lbl::line_key::ls_coeff,
                  "The line shape coefficient\n\n.. :class:`~pyarts3.arts.LineShapeModelCoefficient`");
  line_key.def_rw("var", &lbl::line_key::var, "The variable\n\n.. :class:`~pyarts3.arts.LineByLineVariable`");
  line_key.doc() = "A key for a line";

  py::class_<lbl::temperature::data> tm(m, "TemperatureModel");
  generic_interface(tm);
  tm.def(py::init<LineShapeModelType, Vector>(), "type"_a, "data"_a = Vector{0.0})
      .def_prop_rw(
          "type",
          &lbl::temperature::data::Type,
          [](lbl::temperature::data& self, LineShapeModelType x) { self = lbl::temperature::data{x, self.X()}; },
          "The type of the model\n\n.. :class:`~pyarts3.arts.LineShapeModelType`")
      .def_prop_rw(
          "data",
          [](lbl::temperature::data& self) { return self.X(); },
          [](lbl::temperature::data& self, const Vector& x) { self = lbl::temperature::data{self.Type(), x}; },
          "The coefficients\n\n.. :class:`~pyarts3.arts.Vector`")

      .doc() = "Temperature model";

  py::class_<lbl::line_shape::species_model> lssm(m, "LineShapeSpeciesModel");
  generic_interface(lssm);
  lssm.def_rw("data",
              &lbl::line_shape::species_model::data,
              "The data\n\n.. :class:`dict[tuple[LineShapeModelVariable, TemperatureModel]]`")
      .def("__getitem__", [](py::object& x, const py::object& key) { return x.attr("data").attr("__getitem__")(key); })
      .def("__setitem__",
           [](py::object& x, const py::object& key, const py::object& val) {
             x.attr("data").attr("__setitem__")(key, val);
           })
      .def(
          "G0",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.G0(t0, t, p); }, T0, T, P);
          },
          R"(Computes the G0 coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The G0 coefficient(s) [Hz]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "G2",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.G2(t0, t, p); }, T0, T, P);
          },
          R"(Computes the G2 coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The G2 coefficient(s) [Hz]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "D0",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.D0(t0, t, p); }, T0, T, P);
          },
          R"(Computes the D0 coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The D0 coefficient(s) [Hz]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "D2",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.D2(t0, t, p); }, T0, T, P);
          },
          R"(Computes the D2 coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The D2 coefficient(s) [Hz]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "ETA",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.ETA(t0, t, p); }, T0, T, P);
          },
          R"(Computes the ETA coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The ETA coefficient(s) [dimensionless]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "G",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.G(t0, t, p); }, T0, T, P);
          },
          R"(Computes the G coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The G coefficient(s) [dimensionless]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "Y",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.Y(t0, t, p); }, T0, T, P);
          },
          R"(Computes the Y coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The Y coefficient(s) [dimensionless]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "DV",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.DV(t0, t, p); }, T0, T, P);
          },
          R"(Computes the DV coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The DV coefficient(s) [Hz]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .def(
          "FVC",
          [](const lbl::line_shape::species_model& self, py::object& T0, py::object& T, py::object& P) {
            return vectorize([&self](Numeric t0, Numeric t, Numeric p) { return self.FVC(t0, t, p); }, T0, T, P);
          },
          R"(Computes the FVC coefficient for the given conditions.

Parameters
----------
T0 : Numeric or array-like
    The reference temperature(s) [K]
T : Numeric or array-like
    The temperature(s) [K]
P : Numeric or array-like
    The pressure(s) [Pa]

Returns
-------
Numeric or array-like
    The FVC coefficient(s) [Hz]
)",
          "T0"_a,
          "T"_a,
          "P"_a)
      .doc() = "Line shape model for a species";

  using line_shape_model_list = std::vector<lbl::line_shape::species_model>;
  auto lsml  = py::bind_vector<line_shape_model_list, py::rv_policy::reference_internal>(m, "LineShapeModelList");
  lsml.doc() = "A list of line shape models";
  vector_interface(lsml);
  generic_interface(lsml);

  py::class_<lbl::line_shape::model> lsm(m, "LineShapeModel");
  generic_interface(lsm);
  lsm.def_rw("T0", &lbl::line_shape::model::T0, "The reference temperature [K]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("single_models",
              &lbl::line_shape::model::single_models,
              "The single models\n\n.. :class:`dict[SpeciesEnum, LineShapeSpeciesModel]`")
      .def("G0",
           &lbl::line_shape::model::G0,
           R"(Computes the G0 coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The G0 coefficient(s) [Hz]
)",
           "atm"_a)
      .def("D0",
           &lbl::line_shape::model::D0,
           R"(Computes the D0 coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The D0 coefficient(s) [Hz]
)",
           "atm"_a)
      .def("DV",
           &lbl::line_shape::model::DV,
           R"(Computes the DV coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The DV coefficient(s) [Hz]
)",
           "atm"_a)
      .def("D2",
           &lbl::line_shape::model::D2,
           R"(Computes the D2 coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The D2 coefficient(s) [Hz]
)",
           "atm"_a)
      .def("G2",
           &lbl::line_shape::model::G2,
           R"(Computes the G2 coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The G2 coefficient(s) [Hz]
)",
           "atm"_a)
      .def("FVC",
           &lbl::line_shape::model::FVC,
           R"(Computes the FVC coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The FVC coefficient(s) [Hz]
)",
           "atm"_a)
      .def("G",
           &lbl::line_shape::model::G,
           R"(Computes the G coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The G coefficient(s) [dimensionless]
)",
           "atm"_a)
      .def("Y",
           &lbl::line_shape::model::Y,
           R"(Computes the Y coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The Y coefficient(s) [dimensionless]
)",
           "atm"_a)
      .def("ETA",
           &lbl::line_shape::model::ETA,
           R"(Computes the ETA coefficient

Parameters
----------
atm : ~pyarts3.arts.AtmPoint
    The atmospheric point - must contain pressure and temperature and VMR of all species in the model

Returns
-------
Numeric or array-like
    The ETA coefficient(s) [dimensionless]
)",
           "atm"_a)
      .def(
          "remove",
          [](lbl::line_shape::model& self, LineShapeModelVariable x) {
            for (auto& mod : self.single_models | stdv::values) {
              auto ptr = mod.data.find(x);
              if (ptr != mod.data.end()) mod.data.erase(ptr);
            }
          },
          "x"_a,
          R"(Remove a type of variable from the line shape model.

Parameters
----------
x : LineShapeModelVariable
    The variable to remove
)")
      .def(
          "remove_zeros",
          [](lbl::line_shape::model& self) { self.clear_zeroes(); },
          "Remove zero coefficients from the line shape model")
      .doc() = "Line shape model";

  py::class_<lbl::zeeman::model> zlm(m, "ZeemanLineModel");
  generic_interface(zlm);
  zlm.def_rw("on", &lbl::zeeman::model::on, "If True, the Zeeman effect is included\n\n.. :class:`bool`")
      .def_prop_rw(
          "gl",
          [](lbl::zeeman::model& z) { return z.gl(); },
          [](lbl::zeeman::model& z, Numeric g) { z.gl(g); },
          "The lower level statistical weight\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_prop_rw(
          "gu",
          [](lbl::zeeman::model& z) { return z.gu(); },
          [](lbl::zeeman::model& z, Numeric g) { z.gu(g); },
          "The upper level statistical weight\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def(
          "strengths",
          [](const lbl::zeeman::model& mod, const QuantumState& qn) {
            std::map<std::string, std::vector<double>> out;

            const Index Npi = mod.size(qn, ZeemanPolarization::pi);
            for (Index i = 0; i < Npi; i++) { out["pi"].push_back(mod.Strength(qn, ZeemanPolarization::pi, i)); }

            const Index Nsp = mod.size(qn, ZeemanPolarization::sp);
            for (Index i = 0; i < Nsp; i++) { out["sp"].push_back(mod.Strength(qn, ZeemanPolarization::sp, i)); }

            const Index Nsm = mod.size(qn, ZeemanPolarization::sm);
            for (Index i = 0; i < Nsm; i++) { out["sm"].push_back(mod.Strength(qn, ZeemanPolarization::sm, i)); }

            return out;
          },
          "qn"_a,
          R"(The relative strengths of the Zeeman components for the given quantum numbers state.
Parameters
----------
qn : QuantumState
    The quantum numbers of the line. Must be an instance of the QuantumState class as provided by this module.
Returns
-------
dict[str, list[float]]
)")
      .doc() = "Zeeman model";

  py::class_<lbl::line> al(m, "AbsorptionLine");
  generic_interface(al);
  al.def_rw("a", &lbl::line::a, "The Einstein coefficient [1 / s]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("f0", &lbl::line::f0, "The line center frequency [Hz]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("e0", &lbl::line::e0, "The lower level energy [J]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("gu", &lbl::line::gu, "The upper level statistical weight [-]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("gl", &lbl::line::gl, "The lower level statistical weight [-]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("z", &lbl::line::z, "The Zeeman model\n\n.. :class:`~pyarts3.arts.ZeemanLineModel`")
      .def_rw("ls", &lbl::line::ls, "The line shape model\n\n.. :class:`~pyarts3.arts.LineShapeModel`")
      .def_rw("qn", &lbl::line::qn, "The local quantum numbers of this line\n\n.. :class:`~pyarts3.arts.QuantumState`")
      .def(
          "s",
          [](const lbl::line& self, py::object& T, py::object& Q) {
            return vectorize([&self](Numeric t, Numeric q) { return self.s(t, q); }, T, Q);
          },
          "T"_a,
          "Q"_a,
          R"(The line strength
Parameters
----------
T : Numeric or array-like
    The temperature(s) [K]
Q : Numeric or array-like
    The partition function(s) [dimensionless]

Returns
-------
Numeric or array-like
    The line strength(s)
)")
      .def(
          "hitran_s",
          [](const lbl::line& self, const SpeciesIsotope& isot, Numeric T0) { return self.hitran_s(isot, T0); },
          "isot"_a,
          "T0"_a = 296.0,
          R"(The HITRAN-like line strength
Parameters
----------
isot : SpeciesIsotope
    The species and isotope of the line.
T0 : Numeric
    The reference temperature [K]. Defaults to 296.0 K.

Returns
-------
Numeric
    The HITRAN-like line strength
)")
      .doc() = "A single absorption line";

  auto ll  = py::bind_vector<std::vector<lbl::line>, py::rv_policy::reference_internal>(m, "ArrayOfAbsorptionLine");
  ll.doc() = "A list of :class:`~pyarts3.arts.AbsorptionLine`";
  vector_interface(ll);
  generic_interface(ll);

  py::class_<AbsorptionBand> ab(m, "AbsorptionBand");
  generic_interface(ab);
  ab.def(
      "__getitem__",
      [](const py::object& x, py::object& i) { return x.attr("lines").attr("__getitem__")(i); },
      py::rv_policy::reference_internal);
  ab.def("__setitem__",
         [](py::object& x, const py::object& i, const py::object& v) { x.attr("lines").attr("__setitem__")(i, v); });
  ab.def("__len__", [](const AbsorptionBand& x) { return x.lines.size(); }, "Return the number of lines in the band");
  ab.def_rw("lines", &AbsorptionBand::lines, "The lines in the band\n\n.. :class:`~pyarts3.arts.ArrayOfAbsorptionLine`")
      .def_rw("lineshape",
              &AbsorptionBand::lineshape,
              "The lineshape type\n\n.. :class:`~pyarts3.arts.LineByLineLineshape`")
      .def_prop_rw(
          "cutoff",
          [](const AbsorptionBand& band) { return band.cutoff.type; },
          [](AbsorptionBand& band, LineByLineCutoffType x) { band.cutoff.type = x; },
          "The cutoff type\n\n.. :class:`~pyarts3.arts.LineByLineCutoffType`")
      .def_prop_rw(
          "cutoff_value",
          [](const AbsorptionBand& band) { return band.cutoff.value; },
          [](AbsorptionBand& band, Numeric x) { band.cutoff.value = x; },
          "The cutoff value [Hz]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def(
          "keep_frequencies",
          [](AbsorptionBand& band, Vector2 freqs) {
            band.sort();
            auto                   l = band.active_lines(freqs[0], freqs[1]).second;
            std::vector<lbl::line> new_lines(l.begin(), l.end());
            band.lines = std::move(new_lines);
          },
          "freqs"_a,
          "Keep only the lines within the given frequency range")
      .def(
          "keep_hitran_s",
          [](AbsorptionBand& band, Numeric min_s, const SpeciesIsotope& isot, Numeric T0) {
            std::erase_if(band.lines, [&isot, &T0, &min_s](auto& line) { return line.hitran_s(isot, T0) < min_s; });
          },
          "min_s"_a,
          "isot"_a,
          "T0"_a = 296.0,
          "Keep only the lines with a stronger HITRAN-like line strength");

  auto aoab = py::bind_map<AbsorptionBands, py::rv_policy::reference_internal>(m, "AbsorptionBands");
  generic_interface(aoab);
  aoab.def("__getitem__",
           [](const AbsorptionBands& x, const lbl::line_key& key) -> Numeric { return key.get_value(x); })
      .def("__setitem__", [](AbsorptionBands& x, const lbl::line_key& key, Numeric v) { key.get_value(x) = v; });
  aoab.def(
      "extract_species",
      [](const AbsorptionBands& x, const std::variant<SpeciesEnum, SpeciesIsotope>& vkey) -> AbsorptionBands {
        AbsorptionBands out;
        std::visit(
            [&]<typename T>(const T& key) {
              if constexpr (std::same_as<T, SpeciesEnum>) {
                for (auto& [k, v] : x) {
                  if (k.isot.spec == key) out.try_emplace(k, v);
                }
              } else {
                for (auto& [k, v] : x) {
                  if (k.isot == key) out.try_emplace(k, v);
                }
              }
            },
            vkey);

        return out;
      },
      "spec"_a,
      R"(Extract absorption bands for a given species or isotope.

Parameters
----------
vkey : SpeciesEnum or SpeciesIsotope
    The species or isotope to extract

Returns
-------
AbsorptionBands
    The extracted absorption bands
)");
  aoab.def(
      "merge",
      [](AbsorptionBands& self, const AbsorptionBands& other) -> std::pair<Size, Size> {
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
        using enum LineShapeModelVariable;

        Size sum = 0;
        for (auto& band : self | stdv::values) {
          for (auto& line : band.lines) {
            for (auto& lsm : line.ls.single_models | stdv::values) { sum += lsm.remove_variables<Y, G, DV>(); }
          }
        }

        return sum;
      },
      R"(Clear the linemixing data from all bands by removing it from the inner line shape models

Returns
-------
count : int
    The number of removed variables.
)");

  aoab.def(
      "count_lines",
      [](const AbsorptionBands& self, SpeciesEnum spec) {
        Size n = 0;

        if (spec == SpeciesEnum::Bath) {
          for (auto& band : self | stdv::values) n += band.size();
          return n;
        }
        for (auto& [key, band] : self) {
          if (spec == key.isot.spec) n += band.size();
        }
        return n;
      },
      "spec"_a = SpeciesEnum::Bath,
      "Return the total number of lines");

  aoab.def("remove_hitran_s",
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
      [](const AbsorptionBands&                                                 self,
         const std::variant<Numeric, std::unordered_map<SpeciesEnum, Numeric>>& percentile,
         const Numeric                                                          T0) {
        return std::visit([&](auto& i) { return lbl::percentile_hitran_s(self, i, T0); }, percentile);
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
      [](AbsorptionBands&                                                       self,
         const std::variant<Numeric, std::unordered_map<SpeciesEnum, Numeric>>& percentile,
         const Numeric                                                          T0) {
        lbl::keep_hitran_s(
            self, std::visit([&](auto& i) { return lbl::percentile_hitran_s(self, i, T0); }, percentile), T0);
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
  aoab.def(
      "keep_frequencies",
      [](AbsorptionBands& bands, const Numeric& fmin, const Numeric& fmax) {
        abs_bandsSelectFrequencyByLine(bands, fmin, fmax);
      },
      "fmin"_a = -std::numeric_limits<Numeric>::infinity(),
      "fmax"_a = std::numeric_limits<Numeric>::infinity(),
      R"(Keep the frequencies within the specified range

Parameters
----------
fmin : ~pyarts3.arts.Numeric
    Minimum frequency
fmax : ~pyarts3.arts.Numeric
    Maximum frequency
)");

  py::class_<LinemixingSingleEcsData> ed(m, "LinemixingSingleEcsData");
  generic_interface(ed);
  ed.def_rw("scaling", &LinemixingSingleEcsData::scaling, ".. :class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("beta", &LinemixingSingleEcsData::beta, ".. :class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("lambda_",  // Fix name, not python
            &LinemixingSingleEcsData::lambda,
            ".. :class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("collisional_distance",
            &LinemixingSingleEcsData::collisional_distance,
            ".. :class:`~pyarts3.arts.TemperatureModel`");
  ed.def("Q", &LinemixingSingleEcsData::Q, "J"_a, "T"_a, "T0"_a, "energy"_a, R"(The Q coefficient for the ECS model)");
  ed.def("Omega",
         &LinemixingSingleEcsData::Omega,
         "T"_a,
         "T0"_a,
         "mass"_a,
         "other_mass"_a,
         "energy_x"_a,
         "energy_xm2"_a,
         R"(The Omega coefficient for the ECS model)");

  auto lsed = py::bind_map<LinemixingSpeciesEcsData, py::rv_policy::reference_internal>(m, "LinemixingSpeciesEcsData");
  generic_interface(lsed);

  auto led = py::bind_map<LinemixingEcsData, py::rv_policy::reference_internal>(m, "LinemixingEcsData");
  generic_interface(led);

  namespace hadded = lbl::voigt::ecs::hadded;
  py::class_<hadded::rotational_line>(
      lbl, "hadded_rotational_line", "Prepared NH3 parallel-band rotational transition.")
      .def(
          "__init__",
          [](hadded::rotational_line* self, Index Ju, Index Jl, Index K, bool lower_antisymmetric) {
            using enum hadded::inversion;
            new (self) hadded::rotational_line{
                .upper = {.J = Ju, .K = K, .symmetry = lower_antisymmetric ? symmetric : antisymmetric},
                .lower = {.J = Jl, .K = K, .symmetry = lower_antisymmetric ? antisymmetric : symmetric}};
          },
          "Ju"_a,
          "Jl"_a,
          "K"_a,
          "lower_antisymmetric"_a = false,
          "A parallel-band line with equal upper/lower K and opposite inversion symmetries.")
      .def_prop_ro(
          "Ju",
          [](const hadded::rotational_line& line) { return line.upper.J; },
          "Upper-state angular momentum.\n\n.. :class:`~pyarts3.arts.Index`")
      .def_prop_ro(
          "Jl",
          [](const hadded::rotational_line& line) { return line.lower.J; },
          "Lower-state angular momentum.\n\n.. :class:`~pyarts3.arts.Index`")
      .def_prop_ro(
          "K",
          [](const hadded::rotational_line& line) { return line.lower.K; },
          "Body-fixed projection shared by both states.\n\n.. :class:`~pyarts3.arts.Index`")
      .def_prop_ro(
          "lower_antisymmetric",
          [](const hadded::rotational_line& line) { return line.lower.symmetry == hadded::inversion::antisymmetric; },
          "Whether the lower state has antisymmetric inversion symmetry.\n\n.. :class:`bool`");
  py::class_<hadded::collision_channel>(
      lbl, "hadded_collision_channel", "Signed angular channel of the NH3 collision basis.")
      .def(
          "__init__",
          [](hadded::collision_channel* self, Index L, Index Mi, Index Mf) {
            new (self) hadded::collision_channel{.L = L, .Mi = Mi, .Mf = Mf};
          },
          "L"_a,
          "Mi"_a,
          "Mf"_a,
          "A collision channel with signed body-fixed projections Mi and Mf (multiples of three).")
      .def_ro("L", &hadded::collision_channel::L, "Collision angular rank.\n\n .. :class:`~pyarts3.arts.Index`")
      .def_ro("Mi",
              &hadded::collision_channel::Mi,
              "Signed lower-state projection transfer.\n\n .. :class:`~pyarts3.arts.Index`")
      .def_ro("Mf",
              &hadded::collision_channel::Mf,
              "Signed upper-state projection transfer.\n\n .. :class:`~pyarts3.arts.Index`");
  py::class_<hadded::basis_data>(lbl, "hadded_basis_data", "Prepared NH3 collision rates and adiabatic factors.")
      .def(
          "__init__",
          [](hadded::basis_data* self, std::vector<hadded::collision_channel> channels, Vector Q, Vector Omega) {
            new (self) hadded::basis_data{std::move(channels), std::move(Q), std::move(Omega)};
          },
          "channels"_a,
          "Q"_a,
          "Omega"_a,
          R"(Prepared dynamical factors in channel order.

Q has the desired relaxation-matrix units (Hz for spectra); convert cross sections
in m^2 with number_density * mean_relative_speed / (2*pi). Omega is the paper's
factor >= 1. Supply every signed channel explicitly; omitted channels are zero.
No collision calibration or temperature dependence is supplied by this class.)")
      .def_rw("channels",
              &hadded::basis_data::channels,
              "Signed collision channels in basis order.\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("Q",
              &hadded::basis_data::Q,
              "Dynamical rates in channel order; Hz for spectra.\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("Omega",
              &hadded::basis_data::Omega,
              "Paper-I adiabatic factors in channel order, each at least one.\n\n.. :class:`~pyarts3.arts.Vector`");
  lbl.def("hadded_rotational_energy",
          &hadded::rotational_energy,
          "J"_a,
          "K"_a,
          "B"_a,
          "C"_a,
          "Rigid symmetric-top energy [J]; rotational constants B and C are in J.");
  lbl.def("hadded_reduced_dipole",
          &hadded::reduced_dipole,
          "line"_a,
          "Signed reduced dipole for populations containing the LOWER-state rotational degeneracy.");
  lbl.def(
      "hadded_adiabatic_factors",
      [](const Vector& gap, Numeric duration) {
        Vector Omega(gap.size());
        hadded::adiabatic_factors(Omega, gap, duration);
        return Omega;
      },
      "gap"_a,
      "duration"_a,
      "Paper-I Eq. 19 from prepared energy gaps [J] and collision duration [s]; returns Omega >= 1.");
  lbl.def(
      "hadded_relaxation_matrix_offdiagonal",
      [](const std::vector<hadded::rotational_line>& lines,
         const hadded::basis_data&                   basis,
         const Vector&                               e0,
         const Vector&                               Omega_line,
         Numeric                                     T,
         const Vector&                               widths) {
        const Size n = lines.size();
        ARTS_USER_ERROR_IF(widths.size() != n, "NH3 widths must have one entry per line")
        Matrix W(n, n, 0.0);
        for (Size i = 0; i < n; ++i) {
          ARTS_USER_ERROR_IF(not std::isfinite(widths[i]) or widths[i] < 0, "NH3 widths must be finite and nonnegative")
          W[i, i] = widths[i];
        }
        hadded::relaxation_matrix_offdiagonal(W, lines, basis, e0, Omega_line, T);
        return W;
      },
      "lines"_a,
      "basis"_a,
      "e0"_a,
      "Omega_line"_a,
      "T"_a,
      "widths"_a,
      R"(Return the prepared NH3 relaxation matrix, preserving supplied diagonal widths.

The four-term IOS angular kernel includes detailed balance and ECS energy
corrections. Setting all Omega factors to one gives IOS with detailed balance.
e0 [J] and Omega_line describe original lower states in line order; use one
consistent energy model for these energies and all adiabatic gaps. T is in K.
W uses row-rate storage W[from,to]; its units match Q and widths (Hz for spectra).
The caller supplies collision data and diagonal widths; no calibration or
truncated-band sum-rule adjustment is applied. Initialize Wigner tables first.)");

  lbl.def(
      "relaxation_matrix_profile",
      [](const Vector& frequency,
         const Vector& f0,
         const Matrix& W,
         const Vector& population,
         const Vector& dipole,
         Numeric       gd_fac) {
        const Size n = f0.size();
        ARTS_USER_ERROR_IF(W.nrows() != static_cast<Index>(n) or W.ncols() != static_cast<Index>(n) or
                               population.size() != n or dipole.size() != n,
                           "Inconsistent prepared relaxation-matrix profile dimensions")
        ARTS_USER_ERROR_IF(not std::isfinite(gd_fac) or gd_fac <= 0, "Doppler width factor must be positive and finite")
        for (Numeric f : frequency)
          ARTS_USER_ERROR_IF(not std::isfinite(f) or f <= 0, "Frequencies must be positive and finite")
        AtmPoint atm;
        atm.temperature = 296;
        atm.pressure    = 0;
        lbl::voigt::ecs::ComputeData data({}, atm);
        data.pop    = population;
        data.dip    = dipole;
        data.vmrs   = Vector{1.0};
        data.gd_fac = gd_fac;
        data.Ws.resize(1, n, n);
        for (Size i = 0; i < n; ++i) {
          ARTS_USER_ERROR_IF(not std::isfinite(f0[i]) or f0[i] <= 0 or not std::isfinite(population[i]) or
                                 population[i] < 0 or not std::isfinite(dipole[i]),
                             "Invalid prepared frequency, population, or dipole")
          ARTS_USER_ERROR_IF((W[i, i] < 0), "Diagonal relaxation widths must be nonnegative")
          for (Size j = 0; j < n; ++j) {
            ARTS_USER_ERROR_IF(not std::isfinite(W[i, j]), "Non-finite relaxation matrix")
            data.Ws[0, i, j] = Complex(i == j ? f0[i] : 0.0, W[i, j]);
          }
        }
        data.core_calc(frequency);
        return std::move(data.shape);
      },
      "frequency"_a,
      "f0"_a,
      "W"_a,
      "population"_a,
      "dipole"_a,
      "gd_fac"_a,
      R"(Evaluate the existing ECS equivalent-line Voigt profile from prepared inputs.

frequency and f0 are in Hz. W is a real relaxation matrix in Hz, in row-rate
storage W[from,to]; pressure shifts may be included in f0. population and signed
dipole must use the same state convention as W. gd_fac is the Gaussian 1/e
half-width divided by frequency. This function reuses the core equivalent-line
solver and Faddeeva profile and returns its raw complex shape.

For physical populations and dipoles, absorption [1/m] is real(shape)/sqrt(pi)
times n_abs * frequency * (1-exp(-h*frequency/(k*T))), where n_abs includes the
absorber abundance and isotope ratio. In the Hadded lower-state convention,
use population = gl*exp(-e0/(k*T))/Qpart and signed dipole magnitude
c*sqrt(A*gu/(8*pi*f0^3*gl)); these preserve catalogue isolated-line strengths.
No density, stimulated-emission, or abundance factor is applied here.)");

  lbl.def(
      "equivalent_lines",
      [](const AbsorptionBand&    band,
         const QuantumIdentifier& qid,
         const LinemixingEcsData& abs_ecs_data,
         const AtmPoint&          atm,
         const Vector&            T) {
        lbl::voigt::ecs::ComputeData com_data({}, atm);

        const auto K = band.front().ls.single_models.size();
        const auto N = band.size();
        const auto M = T.size();

        auto eqv_str = ComplexTensor3(M, K, N);
        auto eqv_val = ComplexTensor3(M, K, N);

        equivalent_values(eqv_str, eqv_val, com_data, qid, band, abs_ecs_data.at(qid.isot), atm, T);

        return std::pair{eqv_str, eqv_val};
      },
      "Compute equivalent lines for a given band",
      "band"_a,
      "qid"_a,
      "abs_ecs_data"_a,
      "atm"_a,
      "T"_a);

  aoab.def(
      "spectral_propmat",
      [](const AbsorptionBands&      self,
         const AscendingGrid&        f,
         const AtmPoint&             atm,
         const SpeciesEnum&          spec,
         const PropagationPathPoint& path_point,
         const LinemixingEcsData&    abs_ecs_data,
         const Index&                no_negative_absorption,
         const py::kwargs&) {
        PropmatVector   spectral_propmat(f.size());
        StokvecVector   nlte_vector(f.size());
        PropmatMatrix   spectral_propmat_jac(0, f.size());
        StokvecMatrix   nlte_matrix(0, f.size());
        JacobianTargets jac_targets{};

        spectral_propmatAddLines(spectral_propmat,
                                 nlte_vector,
                                 spectral_propmat_jac,
                                 nlte_matrix,
                                 f,
                                 jac_targets,
                                 spec,
                                 self,
                                 abs_ecs_data,
                                 atm,
                                 path_point,
                                 no_negative_absorption);

        return spectral_propmat;
      },
      "f"_a,
      "atm"_a,
      "spec"_a                   = SpeciesEnum::Bath,
      "path_point"_a             = PropagationPathPoint{},
      "abs_ecs_data"_a           = LinemixingEcsData{},
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
abs_ecs_data : LinemixingEcsData, optional
    The ECS data.  Default is empty.
no_negative_absorption : Index, optional
    If 1, the absorption is set to zero if it is negative. The default is 1.

Returns
-------
spectral_propmat : PropmatVector
    Propagation matrix by frequency [1/m]

)--");

  py::class_<PartitionFunctionsData> partfun(m, "PartitionFunctionsData");
  partfun.def_rw(
      "data", &PartitionFunctionsData::data, "The partition function data\n\n.. :class:`~pyarts3.arts.Matrix`");
  partfun.def_rw("type",
                 &PartitionFunctionsData::type,
                 "The type of partition function data\n\n.. :class:`~pyarts3.arts.PartitionFunctionsType`");
  generic_interface(partfun);
  partfun.doc() = "Data for partition functions, used in the line-by-line model";
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize lbl\n{}", e.what()));
}
}  // namespace Python
