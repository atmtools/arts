#include <enumsLineShapeModelVariable.h>
#include <enumsSpeciesEnum.h>
#include <hpy_arts.h>
#include <hpy_numpy.h>
#include <hpy_vector.h>
#include <isotopologues.h>
#include <lbl.h>
#include <lbl_data.h>
#include <lbl_lineshape_model.h>
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
#include <python_interface.h>
#include <quantum.h>

#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "matpack_mdspan_common_complex.h"
#include "matpack_mdspan_helpers_grid_t.h"

NB_MAKE_OPAQUE(lbl::line_shape::species_model::map_t)
NB_MAKE_OPAQUE(lbl::line_shape::model::map_t)

namespace Python {
void py_lbl(py::module_& m) try {
  auto lbl = m.def_submodule("lbl", "Line-by-line helper functions");

  auto lssmm = py::bind_map<lbl::line_shape::species_model::map_t>(
      lbl, "LineShapeSpeciesModelMap");
  lssmm.doc() = "A map from model variable to line shape models";
  generic_interface(lssmm);

  auto lsmm =
      py::bind_map<lbl::line_shape::model::map_t>(lbl, "LineShapeModelMap");
  lsmm.doc() = "A map from species to species line shape models";
  generic_interface(lsmm);

  py::class_<lbl::line_key> line_key(lbl, "line_key");
  generic_interface(line_key);
  line_key.def_rw("band",
                  &lbl::line_key::band,
                  "The band\n\n.. :class:`QuantumIdentifier`");
  line_key.def_rw("line", &lbl::line_key::line, "The line\n\n.. :class:`int`");
  line_key.def_rw(
      "spec", &lbl::line_key::spec, "The species\n\n.. :class:`int`");
  line_key.def_rw(
      "ls_var",
      &lbl::line_key::ls_var,
      "The line shape variable\n\n.. :class:`LineShapeModelVariable`");
  line_key.def_rw(
      "ls_coeff",
      &lbl::line_key::ls_coeff,
      "The line shape coefficient\n\n.. :class:`LineShapeModelCoefficient`");
  line_key.def_rw("var",
                  &lbl::line_key::var,
                  "The variable\n\n.. :class:`LineByLineVariable`");
  line_key.doc() = "A key for a line";

  py::class_<lbl::temperature::data> tm(m, "TemperatureModel");
  generic_interface(tm);
  tm.def(py::init<LineShapeModelType, Vector>(),
         "type"_a,
         "data"_a = Vector{0.0})
      .def_prop_rw(
          "type",
          &lbl::temperature::data::Type,
          [](lbl::temperature::data& self, LineShapeModelType x) {
            self = lbl::temperature::data{x, self.X()};
          },
          "The type of the model\n\n.. :class:`~pyarts3.arts.TemperatureModelType`")
      .def_prop_rw(
          "data",
          [](lbl::temperature::data& self) { return self.X(); },
          [](lbl::temperature::data& self, const Vector& x) {
            self = lbl::temperature::data{self.Type(), x};
          },
          "The coefficients\n\n.. :class:`~pyarts3.arts.Vector`")
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

  py::class_<lbl::line_shape::species_model> lssm(m, "LineShapeSpeciesModel");
  generic_interface(lssm);
  lssm.def_rw(
          "data",
          &lbl::line_shape::species_model::data,
          "The data\n\n.. :class:`dict[tuple[LineShapeModelVariable, TemperatureModel]]`")
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
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.G2(t0, t, p);
                },
                T0,
                T,
                P);
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
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.D2(t0, t, p);
                },
                T0,
                T,
                P);
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
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.ETA(t0, t, p);
                },
                T0,
                T,
                P);
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
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.G(t0, t, p);
                },
                T0,
                T,
                P);
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
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.DV(t0, t, p);
                },
                T0,
                T,
                P);
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
          [](const lbl::line_shape::species_model& self,
             py::object& T0,
             py::object& T,
             py::object& P) {
            return vectorize(
                [&self](Numeric t0, Numeric t, Numeric p) {
                  return self.FVC(t0, t, p);
                },
                T0,
                T,
                P);
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
  auto lsml =
      py::bind_vector<line_shape_model_list, py::rv_policy::reference_internal>(
          m, "LineShapeModelList");
  lsml.doc() = "A list of line shape models";
  vector_interface(lsml);
  generic_interface(lsml);

  py::class_<lbl::line_shape::model> lsm(m, "LineShapeModel");
  generic_interface(lsm);
  lsm.def_rw("T0",
              &lbl::line_shape::model::T0,
              "The reference temperature [K]\n\n.. :class:`Numeric`")
      .def_rw(
          "single_models",
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
  zlm.def_rw("on",
             &lbl::zeeman::model::on,
             "If True, the Zeeman effect is included\n\n.. :class:`bool`")
      .def_prop_rw(
          "gl",
          [](lbl::zeeman::model& z) { return z.gl(); },
          [](lbl::zeeman::model& z, Numeric g) { z.gl(g); },
          "The lower level statistical weight\n\n.. :class:`Numeric`")
      .def_prop_rw(
          "gu",
          [](lbl::zeeman::model& z) { return z.gu(); },
          [](lbl::zeeman::model& z, Numeric g) { z.gu(g); },
          "The upper level statistical weight\n\n.. :class:`Numeric`")
      .def(
          "strengths",
          [](const lbl::zeeman::model& mod, const QuantumState& qn) {
            std::map<std::string, std::vector<double>> out;

            const Index Npi = mod.size(qn, ZeemanPolarization::pi);
            for (Index i = 0; i < Npi; i++) {
              out["pi"].push_back(mod.Strength(qn, ZeemanPolarization::pi, i));
            }

            const Index Nsp = mod.size(qn, ZeemanPolarization::sp);
            for (Index i = 0; i < Nsp; i++) {
              out["sp"].push_back(mod.Strength(qn, ZeemanPolarization::sp, i));
            }

            const Index Nsm = mod.size(qn, ZeemanPolarization::sm);
            for (Index i = 0; i < Nsm; i++) {
              out["sm"].push_back(mod.Strength(qn, ZeemanPolarization::sm, i));
            }

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
  al.def_rw("a",
            &lbl::line::a,
            "The Einstein coefficient [1 / s]\n\n.. :class:`Numeric`")
      .def_rw("f0",
              &lbl::line::f0,
              "The line center frequency [Hz]\n\n.. :class:`Numeric`")
      .def_rw("e0",
              &lbl::line::e0,
              "The lower level energy [J]\n\n.. :class:`Numeric`")
      .def_rw("gu",
              &lbl::line::gu,
              "The upper level statistical weight [-]\n\n.. :class:`Numeric`")
      .def_rw("gl",
              &lbl::line::gl,
              "The lower level statistical weight [-]\n\n.. :class:`Numeric`")
      .def_rw("z",
              &lbl::line::z,
              "The Zeeman model\n\n.. :class:`~pyarts3.arts.ZeemanLineModel`")
      .def_rw(
          "ls",
          &lbl::line::ls,
          "The line shape model\n\n.. :class:`~pyarts3.arts.LineShapeModel`")
      .def_rw(
          "qn",
          &lbl::line::qn,
          "The local quantum numbers of this line\n\n.. :class:`~pyarts3.arts.QuantumState`")
      .def(
          "s",
          [](const lbl::line& self, py::object& T, py::object& Q) {
            return vectorize(
                [&self](Numeric t, Numeric q) { return self.s(t, q); }, T, Q);
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
          [](const lbl::line& self, const SpeciesIsotope& isot, Numeric T0) {
            return self.hitran_s(isot, T0);
          },
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

  auto ll = py::bind_vector<std::vector<lbl::line>,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfAbsorptionLine");
  ll.doc() = "A list of :class:`AbsorptionLine`";
  vector_interface(ll);
  generic_interface(ll);

  py::class_<AbsorptionBand> ab(m, "AbsorptionBand");
  generic_interface(ab);
  ab.def(
      "__getitem__",
      [](const py::object& x, py::object& i) {
        return x.attr("lines").attr("__getitem__")(i);
      },
      py::rv_policy::reference_internal);
  ab.def("__setitem__",
         [](py::object& x, const py::object& i, const py::object& v) {
           x.attr("lines").attr("__setitem__")(i, v);
         });
  ab.def(
      "__len__",
      [](const AbsorptionBand& x) { return x.lines.size(); },
      "Return the number of lines in the band");
  ab.def_rw("lines",
            &AbsorptionBand::lines,
            "The lines in the band\n\n.. :class:`ArrayOfAbsorptionLine`")
      .def_rw("lineshape",
              &AbsorptionBand::lineshape,
              "The lineshape type\n\n.. :class:`LineByLineLineshape`")
      .def_rw("cutoff",
              &AbsorptionBand::cutoff,
              "The cutoff type\n\n.. :class:`LineByLineCutoffType`")
      .def_rw("cutoff_value",
              &AbsorptionBand::cutoff_value,
              "The cutoff value [Hz]\n\n.. :class:`Numeric`")
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
      "extract_species",
      [](const AbsorptionBands& x,
         const std::variant<SpeciesEnum, SpeciesIsotope>& vkey)
          -> AbsorptionBands {
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
        using enum LineShapeModelVariable;

        Size sum = 0;
        for (auto& band : self | stdv::values) {
          for (auto& line : band.lines) {
            for (auto& lsm : line.ls.single_models | stdv::values) {
              sum += lsm.remove_variables<Y, G, DV>();
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
  ed.def_rw("scaling",
            &LinemixingSingleEcsData::scaling,
            ".. :class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("beta",
            &LinemixingSingleEcsData::beta,
            ".. :class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("lambda_",  // Fix name, not python
            &LinemixingSingleEcsData::lambda,
            ".. :class:`~pyarts3.arts.TemperatureModel`");
  ed.def_rw("collisional_distance",
            &LinemixingSingleEcsData::collisional_distance,
            ".. :class:`~pyarts3.arts.TemperatureModel`");
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

  auto lsed =
      py::bind_map<LinemixingSpeciesEcsData>(m, "LinemixingSpeciesEcsData");
  generic_interface(lsed);

  auto led = py::bind_map<LinemixingEcsData>(m, "LinemixingEcsData");
  generic_interface(led);

  lbl.def(
      "equivalent_lines",
      [](const AbsorptionBand& band,
         const QuantumIdentifier& qid,
         const LinemixingEcsData& abs_ecs_data,
         const AtmPoint& atm,
         const Vector& T) {
        lbl::voigt::ecs::ComputeData com_data({}, atm);

        const auto K = band.front().ls.single_models.size();
        const auto N = band.size();
        const auto M = T.size();

        auto eqv_str = ComplexTensor3(M, K, N);
        auto eqv_val = ComplexTensor3(M, K, N);

        equivalent_values(eqv_str,
                          eqv_val,
                          com_data,
                          qid,
                          band,
                          abs_ecs_data.at(qid.isot),
                          atm,
                          T);

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
      [](const AbsorptionBands& self,
         const AscendingGrid& f,
         const AtmPoint& atm,
         const SpeciesEnum& spec,
         const PropagationPathPoint& path_point,
         const LinemixingEcsData& abs_ecs_data,
         const Index& no_negative_absorption,
         const py::kwargs&) {
        PropmatVector spectral_propmat(f.size());
        StokvecVector nlte_vector(f.size());
        PropmatMatrix spectral_propmat_jac(0, f.size());
        StokvecMatrix nlte_matrix(0, f.size());
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
  generic_interface(partfun);
  partfun.doc() =
      "Data for partition functions, used in the line-by-line model";
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize lbl\n{}", e.what()));
}
}  // namespace Python
