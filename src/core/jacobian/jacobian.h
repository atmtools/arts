#pragma once

#include <atm.h>
#include <lbl.h>
#include <matpack.h>
#include <obsel.h>
#include <surf.h>

#include <concepts>
#include <functional>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <variant>
#include <vector>

#include "operators.h"

using AtmTargetSetState     = CustomOperator<void,
                                             ExhaustiveVectorView,
                                             const AtmField&,
                                             const AtmKeyVal&>;
using AtmTargetSetModel     = CustomOperator<void,
                                             AtmField&,
                                             const AtmKeyVal&,
                                             const ExhaustiveConstVectorView>;
using SurfaceTargetSetState = CustomOperator<void,
                                             ExhaustiveVectorView,
                                             const SurfaceField&,
                                             const SurfaceKeyVal&>;
using SurfaceTargetSetModel = CustomOperator<void,
                                             SurfaceField&,
                                             const SurfaceKeyVal&,
                                             const ExhaustiveConstVectorView>;
using LineTargetSetState    = CustomOperator<void,
                                             ExhaustiveVectorView,
                                             const AbsorptionBands&,
                                             const LblLineKey&>;
using LineTargetSetModel    = CustomOperator<void,
                                             AbsorptionBands&,
                                             const LblLineKey&,
                                             const ExhaustiveConstVectorView>;
using SensorTargetSetState  = CustomOperator<void,
                                             ExhaustiveVectorView,
                                             const ArrayOfSensorObsel&,
                                             const SensorKey&>;
using SensorTargetSetModel  = CustomOperator<void,
                                             ArrayOfSensorObsel&,
                                             const SensorKey&,
                                             const ExhaustiveConstVectorView>;

namespace Jacobian {
void default_atm_x_set(ExhaustiveVectorView, const AtmField&, const AtmKeyVal&);
void default_x_atm_set(AtmField&,
                       const AtmKeyVal&,
                       const ExhaustiveConstVectorView);

struct AtmTarget {
  AtmKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  AtmTargetSetState set_state{default_atm_x_set};

  AtmTargetSetModel set_model{default_x_atm_set};

  void update(AtmField& atm, const Vector& x) const;

  void update(Vector& x, const AtmField& atm) const;

  [[nodiscard]] bool is_wind() const;
};

void default_surf_x_set(ExhaustiveVectorView,
                        const SurfaceField&,
                        const SurfaceKeyVal&);
void default_x_surf_set(SurfaceField&,
                        const SurfaceKeyVal&,
                        const ExhaustiveConstVectorView);

struct SurfaceTarget {
  SurfaceKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  SurfaceTargetSetState set_state{default_surf_x_set};

  SurfaceTargetSetModel set_model{default_x_surf_set};

  void update(SurfaceField& surf, const Vector& x) const;

  void update(Vector& x, const SurfaceField& surf) const;
};

void default_line_x_set(ExhaustiveVectorView,
                        const AbsorptionBands&,
                        const LblLineKey&);
void default_x_line_set(AbsorptionBands&,
                        const LblLineKey&,
                        const ExhaustiveConstVectorView);

struct LineTarget {
  LblLineKey type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  LineTargetSetState set_state{default_line_x_set};

  LineTargetSetModel set_model{default_x_line_set};

  void update(AbsorptionBands&, const Vector&) const;

  void update(Vector&, const AbsorptionBands&) const;
};

void default_sensor_x_set(ExhaustiveVectorView,
                          const ArrayOfSensorObsel&,
                          const SensorKey&);
void default_x_sensor_set(ArrayOfSensorObsel&,
                          const SensorKey&,
                          const ExhaustiveConstVectorView);

struct SensorTarget {
  SensorKey type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  SensorTargetSetState set_state{default_sensor_x_set};

  SensorTargetSetModel set_model{default_x_sensor_set};

  void update(ArrayOfSensorObsel&, const Vector&) const;

  void update(Vector&, const ArrayOfSensorObsel&) const;
};

template <typename T>
concept target_type = requires(T a) {
  { a.d } -> std::same_as<Numeric&>;
  { a.target_pos } -> std::same_as<Size&>;
  { a.x_start } -> std::same_as<Size&>;
  { a.x_size } -> std::same_as<Size&>;
  a.type;
};

template <typename U, typename T>
concept target_comparable = target_type<T> and requires(T a, U b) {
  a.type == b;
  b == a.type;
};

template <typename U, typename... T>
concept valid_target = (std::same_as<U, T> or ...);

template <target_type... Targets>
struct targets_t {
  static constexpr Size N = sizeof...(Targets);

  std::tuple<std::vector<Targets>...> targets{};

  bool finalized{false};

  template <valid_target<Targets...> T>
  constexpr auto& target() {
    return std::get<std::vector<T>>(targets);
  }

  template <valid_target<Targets...> T>
  [[nodiscard]] constexpr const auto& target() const {
    return std::get<std::vector<T>>(targets);
  }

  template <valid_target<Targets...> T>
  constexpr auto find(const target_comparable<T> auto& t) const {
    auto out =
        std::pair{false, std::ranges::find_if(target<T>(), [&t](auto& v) {
                    return v.type == t;
                  })};
    out.first = out.second not_eq target<T>().end();
    return out;
  }

  template <valid_target<Targets...> T>
  constexpr auto find_all(const target_comparable<T> auto&... t) const {
    return std::array{find<T>(t)...};
  }

  template <valid_target<Targets...> T>
  constexpr bool contains(const target_comparable<T> auto& t) const {
    return find<T>(t).first;
  }

  template <valid_target<Targets...> T>
  constexpr Index target_position(const target_comparable<T> auto& t) const {
    auto pair = find<T>(t);
    return pair.first ? pair.second->target_pos : -1;
  }

  [[nodiscard]] constexpr Size x_size() const {
    ARTS_USER_ERROR_IF(target_count() != 0 and not finalized, "Not finalized.")

    const auto sz = [](const auto& x) { return x.x_size; };
    return (std::transform_reduce(target<Targets>().begin(),
                                  target<Targets>().end(),
                                  Size{0},
                                  std::plus<>{},
                                  sz) +
            ...);
  }

  [[nodiscard]] constexpr Size target_count() const {
    return (target<Targets>().size() + ...);
  }

  [[nodiscard]] constexpr bool any() const { return target_count() != 0; }

  void throwing_check(Size xsize) const {
    const auto t_size = target_count();

    ARTS_USER_ERROR_IF(
        xsize not_eq x_size(),
        "The size of the x-vector does not match the size of the targets.")

    ((std::ranges::for_each(
         target<Targets>(),
         [xsize, t_size](auto& a) {
           ARTS_USER_ERROR_IF((a.x_start + a.x_size) > xsize,
                              "The target {}"
                              " is out of bounds of the x-vector.  (xsize: {})",
                              a,
                              xsize)
           ARTS_USER_ERROR_IF(
               t_size <= a.target_pos,
               "The target {}"
               " is out of bounds of the target vector.  (t_size: {})",
               a,
               t_size)
         })),
     ...);
  }

  void clear() { ((target<Targets>().clear()), ...); }

  void zero_out_x() {
    ((std::ranges::for_each(target<Targets>(),
                            [](auto& a) {
                              a.x_start = 0;
                              a.x_size  = 0;
                            })),
     ...);
  }
};

struct Targets final : targets_t<AtmTarget, SurfaceTarget, LineTarget, SensorTarget> {
  [[nodiscard]] const std::vector<AtmTarget>& atm() const;
  [[nodiscard]] const std::vector<SurfaceTarget>& surf() const;
  [[nodiscard]] const std::vector<LineTarget>& line() const;
  [[nodiscard]] const std::vector<SensorTarget>& sensor() const;

  [[nodiscard]] std::vector<AtmTarget>& atm();
  [[nodiscard]] std::vector<SurfaceTarget>& surf();
  [[nodiscard]] std::vector<LineTarget>& line();
  [[nodiscard]] std::vector<SensorTarget>& sensor();

  //! Sets the sizes and x-positions of the targets.
  void finalize(const AtmField& atmospheric_field,
                const SurfaceField& surface_field,
                const AbsorptionBands& absorption_bands,
                const ArrayOfSensorObsel& measurement_sensor);
};

struct TargetType {
  using variant_t = std::variant<AtmKeyVal, SurfaceKeyVal, LblLineKey, SensorKey>;
  variant_t target;

  template <class AtmKeyValFunc, class SurfaceKeyValFunc, class LblLineKeyFunc, class SensorKeyFunc>
  [[nodiscard]] constexpr auto apply(const AtmKeyValFunc& ifatm,
                                     const SurfaceKeyValFunc& ifsurf,
                                     const LblLineKeyFunc& ifline,
                                     const SensorKeyFunc& ifsensor) const {
    return std::visit(
        [&](const auto& t) {
          using T = std::decay_t<decltype(t)>;
          if constexpr (std::same_as<T, AtmKeyVal>) {
            return ifatm(t);
          } else if constexpr (std::same_as<T, SurfaceKeyVal>) {
            return ifsurf(t);
          } else if constexpr (std::same_as<T, LblLineKey>) {
            return ifline(t);
          } else if constexpr (std::same_as<T, SensorKey>) {
            return ifsensor(t);
          }
        },
        target);
  }

  constexpr bool operator==(const TargetType&) const = default;

  [[nodiscard]] std::string type() const {
    return apply([](auto&) { return "AtmKeyVal"s; },
                 [](auto&) { return "SurfaceKeyVal"s; },
                 [](auto&) { return "LblLineKey"s; },
                 [](auto&) { return "SensorKey"s; });
  }
};
}  // namespace Jacobian

Numeric field_perturbation(const auto& f) {
  return std::transform_reduce(
      f.begin(),
      f.end(),
      0.0,
      [](auto a, auto b) { return std::max(a, b); },
      [](auto& m) { return m.first ? m.second->d : 0.0; });
}

using JacobianTargets    = Jacobian::Targets;
using JacobianTargetType = Jacobian::TargetType;

template <>
struct std::hash<JacobianTargetType> {
  std::size_t operator()(const JacobianTargetType& g) const {
    return std::hash<JacobianTargetType::variant_t>{}(g.target);
  }
};

template <>
struct std::formatter<JacobianTargetType> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const JacobianTargetType& jtt,
                              FmtContext& ctx) const {
    jtt.apply(
        [this, &ctx](const AtmKeyVal& key) {
          tags.format(ctx, "AtmKeyVal::"sv, key);
        },
        [this, &ctx](const SurfaceKeyVal& key) {
          tags.format(ctx, "SurfaceKeyVal::"sv, key);
        },
        [this, &ctx](const LblLineKey& key) {
          tags.format(ctx, "LblLineKey::"sv, key);
        },
        [this, &ctx](const SensorKey& key) {
          tags.format(ctx, "SensorKey::"sv, key);
        });

    return ctx.out();
  }
};

template <>
struct std::formatter<Jacobian::AtmTarget> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Jacobian::AtmTarget& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(ctx,
                       v.type,
                       ": "sv,
                       v.d,
                       sep,
                       "target_pos: "sv,
                       v.target_pos,
                       sep,
                       "x_start: "sv,
                       v.x_start,
                       sep,
                       "x_size: "sv,
                       v.x_size);
  }
};

template <>
struct std::formatter<Jacobian::SurfaceTarget> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Jacobian::SurfaceTarget& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(ctx,
                       v.type,
                       ": "sv,
                       v.d,
                       sep,
                       "target_pos: "sv,
                       v.target_pos,
                       sep,
                       "x_start: "sv,
                       v.x_start,
                       sep,
                       "x_size: "sv,
                       v.x_size);
  }
};

template <>
struct std::formatter<Jacobian::LineTarget> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Jacobian::LineTarget& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(ctx,
                       v.type,
                       ": "sv,
                       v.d,
                       sep,
                       "target_pos: "sv,
                       v.target_pos,
                       sep,
                       "x_start: "sv,
                       v.x_start,
                       sep,
                       "x_size: "sv,
                       v.x_size);
  }
};

template <>
struct std::formatter<Jacobian::SensorTarget> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Jacobian::SensorTarget& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(ctx,
                       v.type,
                       ": "sv,
                       v.d,
                       sep,
                       "target_pos: "sv,
                       v.target_pos,
                       sep,
                       "x_start: "sv,
                       v.x_start,
                       sep,
                       "x_size: "sv,
                       v.x_size);
  }
};

template <>
struct std::formatter<JacobianTargets> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const JacobianTargets& v, FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '{');

    const std::string_view sep = tags.sep(true);
    tags.format(ctx,
                R"("atm": )"sv,
                v.atm(),
                sep,
                R"("surf": )"sv,
                v.surf(),
                sep,
                R"("line": )"sv,
                v.line());

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};
