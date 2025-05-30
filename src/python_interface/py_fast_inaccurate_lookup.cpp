#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>

#include <unordered_map>

#include "hpy_arts.h"

struct FastLinearGasOptics {
  Vector offset{};
  Vector slopeP{};
  Vector slopeT{};
};

struct FastNonLinearGasOptics {
  FastLinearGasOptics linear{};
  FastLinearGasOptics nonlinear{};
};

using FastGasOpticsUnit =
    std::variant<FastLinearGasOptics, FastNonLinearGasOptics>;

using FastGasOpticsUnitMap = std::unordered_map<SpeciesEnum, FastGasOpticsUnit>;

struct FastGasOptics {
  Numeric p0{1};           // Reference pressure, defaults to 1 Pa
  Numeric T0{250};         // Reference temperature, defaults to 250 K
  AscendingGrid f_grid{};  // Frequency grid
  FastGasOpticsUnitMap data{};
};

// FIX BELOW FOR FORMATTING

template <>
struct std::formatter<FastLinearGasOptics> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const FastLinearGasOptics& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(ctx,
                       "offset: ",
                       v.offset,
                       sep,
                       "slopeP: ",
                       v.slopeP,
                       sep,
                       "slopeT: ",
                       v.slopeT);
  }
};

template <>
struct std::formatter<FastNonLinearGasOptics> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const FastNonLinearGasOptics& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(
        ctx, "linear: ", v.linear, sep, "nonlinear: ", v.nonlinear);
  }
};

template <>
struct std::formatter<FastGasOptics> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const FastGasOptics& v, FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(ctx,
                       "p0: ",
                       v.p0,
                       sep,
                       "T0: ",
                       v.T0,
                       sep,
                       "f_grid: ",
                       v.f_grid,
                       sep,
                       "data: ",
                       v.data);
  }
};

// MOVE THE ABOVE TO PROPER FILES WHEN READY

namespace Python {
void py_fast_inaccurate_loopup(py::module_& m) {
  py::class_<FastLinearGasOptics> flgo(m, "FastLinearGasOptics");
  flgo.def(py::init<Vector, Vector, Vector>(),
           "offset"_a = Vector{},
           "slopeP"_a = Vector{},
           "slopeT"_a = Vector{})
      .def_rw("offset", &FastLinearGasOptics::offset)
      .def_rw("slopeP", &FastLinearGasOptics::slopeP)
      .def_rw("slopeT", &FastLinearGasOptics::slopeT);
  str_interface(flgo);

  py::class_<FastNonLinearGasOptics> fnlgo(m, "FastNonLinearGasOptics");
  fnlgo
      .def(py::init<FastLinearGasOptics, FastLinearGasOptics>(),
           "linear"_a    = FastLinearGasOptics{},
           "nonlinear"_a = FastLinearGasOptics{})
      .def_rw("linear", &FastNonLinearGasOptics::linear)
      .def_rw("nonlinear", &FastNonLinearGasOptics::nonlinear);
  str_interface(fnlgo);

  py::class_<FastGasOptics> fgo(m, "FastGasOptics");
  fgo.def(py::init<Numeric,
                   Numeric,
                   AscendingGrid,
                   std::unordered_map<SpeciesEnum, FastGasOpticsUnit>>(),
          "p0"_a     = Numeric{1},
          "T0"_a     = Numeric{250},
          "f_grid"_a = AscendingGrid{},
          "data"_a   = std::unordered_map<SpeciesEnum, FastGasOpticsUnit>{})
      .def_rw("p0", &FastGasOptics::p0)
      .def_rw("T0", &FastGasOptics::T0)
      .def_rw("f_grid", &FastGasOptics::f_grid)
      .def(
          "__getitem__",
          [](FastGasOptics& self, SpeciesEnum key) -> FastGasOpticsUnit& {
            return self.data.at(key);
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](FastGasOptics& self, SpeciesEnum key, FastGasOpticsUnit val) {
             self.data[key] = std::move(val);
           });
  str_interface(fgo);
}
}  // namespace Python
