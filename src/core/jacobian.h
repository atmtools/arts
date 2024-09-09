#pragma once

#include <atm.h>
#include <lbl.h>
#include <matpack.h>
#include <surf.h>

#include <concepts>
#include <functional>
#include <limits>
#include <numeric>
#include <iosfwd>
#include <unordered_map>
#include <variant>
#include <vector>

namespace Jacobian {
struct AtmTarget {
  AtmKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  std::function<void(ExhaustiveVectorView, const AtmField&, const AtmKeyVal&)>
      set_state{[](ExhaustiveVectorView x,
                   const AtmField& atm,
                   const AtmKeyVal& key) {
        ARTS_USER_ERROR_IF(not atm.contains(key),
                           "Atmosphere does not contain key value {}",
                           key)

        auto xn = atm[key].flat_view();

        ARTS_USER_ERROR_IF(
            atm[key].flat_view().size() not_eq x.size(),
            "Problem with sizes.  \n"
            "Did you change your atmosphere since you set the jacobian targets?\n"
            "Did you forget to finalize the JacobianTargets?")

        x = xn;
      }};

  std::function<void(
      AtmField&, const AtmKeyVal&, const ExhaustiveConstVectorView)>
      set_model{[](AtmField& atm,
                   const AtmKeyVal& key,
                   const ExhaustiveConstVectorView x) {
        ARTS_USER_ERROR_IF(not atm.contains(key),
                           "Atmosphere does not contain key value {}",
                           key)

        auto xn = atm[key].flat_view();

        ARTS_USER_ERROR_IF(
            atm[key].flat_view().size() not_eq x.size(),
            "Problem with sizes.  \n"
            "Did you change your atmosphere since you set the jacobian targets?\n"
            "Did you forget to finalize the JacobianTargets?")

        xn = x;
      }};

  friend std::ostream& operator<<(std::ostream& os, const AtmTarget& target);

  void update(AtmField& atm, const Vector& x) const;

  void update(Vector& x, const AtmField& atm) const;
};

struct SurfaceTarget {
  SurfaceKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  std::function<void(
      ExhaustiveVectorView, const SurfaceField&, const SurfaceKeyVal&)>
      set_state{[](ExhaustiveVectorView x,
                   const SurfaceField& surf,
                   const SurfaceKeyVal& key) {
        ARTS_USER_ERROR_IF(
            not surf.contains(key), "Surface does not contain key value {}", key)

        auto xn = surf[key].flat_view();

        ARTS_USER_ERROR_IF(
            xn.size() not_eq x.size(),
            "Problem with sizes.\n"
            "Did you change your surface since you set the jacobian targets?\n"
            "Did you forget to finalize the JacobianTargets?")

        x = xn;
      }};

  std::function<void(
      SurfaceField&, const SurfaceKeyVal&, const ExhaustiveConstVectorView)>
      set_model{[](SurfaceField& surf,
                   const SurfaceKeyVal& key,
                   const ExhaustiveConstVectorView x) {
        ARTS_USER_ERROR_IF(
            not surf.contains(key), "Surface does not contain key value {}", key)

        auto xn = surf[key].flat_view();

        ARTS_USER_ERROR_IF(
            xn.size() not_eq x.size(),
            "Problem with sizes.\n"
            "Did you change your surface since you set the jacobian targets?\n"
            "Did you forget to finalize the JacobianTargets?")

        xn = x;
      }};

  friend std::ostream& operator<<(std::ostream& os,
                                  const SurfaceTarget& target);

  void update(SurfaceField& surf, const Vector& x) const;

  void update(Vector& x, const SurfaceField& surf) const;
};

struct LineTarget {
  LblLineKey type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  std::function<void(
      ExhaustiveVectorView, const ArrayOfAbsorptionBand&, const LblLineKey&)>
      set_state{[](ExhaustiveVectorView x,
                   const ArrayOfAbsorptionBand& bands,
                   const LblLineKey& key) { x = key.get_value(bands); }};

  std::function<void(ArrayOfAbsorptionBand&,
                     const LblLineKey&,
                     const ExhaustiveConstVectorView)>
      set_model{[](ArrayOfAbsorptionBand& bands,
                   const LblLineKey& key,
                   const ExhaustiveConstVectorView x) {
        ExhaustiveVectorView{key.get_value(bands)} = x;
      }};

  friend std::ostream& operator<<(std::ostream& os, const LineTarget&);

  void update(ArrayOfAbsorptionBand&, const Vector&) const;

  void update(Vector&, const ArrayOfAbsorptionBand&) const;
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
         [&](auto& a) {
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

struct Targets final : targets_t<AtmTarget, SurfaceTarget, LineTarget> {
  [[nodiscard]] const std::vector<AtmTarget>& atm() const;
  [[nodiscard]] const std::vector<SurfaceTarget>& surf() const;
  [[nodiscard]] const std::vector<LineTarget>& line() const;
  [[nodiscard]] std::vector<AtmTarget>& atm();
  [[nodiscard]] std::vector<SurfaceTarget>& surf();
  [[nodiscard]] std::vector<LineTarget>& line();

  //! Sets the sizes and x-positions of the targets.
  void finalize(const AtmField& atmospheric_field,
                const SurfaceField& surface_field,
                const ArrayOfAbsorptionBand& absorption_bands);

  friend std::ostream& operator<<(std::ostream& os, const Targets& targets);
};

struct TargetType {
  using variant_t = std::variant<AtmKeyVal, SurfaceKeyVal, LblLineKey>;
  variant_t target;

  template <class AtmKeyValFunc, class SurfaceKeyValFunc, class LblLineKeyFunc>
  [[nodiscard]] constexpr auto apply(const AtmKeyValFunc& ifatm,
                                     const SurfaceKeyValFunc& ifsurf,
                                     const LblLineKeyFunc& ifline) const {
    return std::visit(
        [&](const auto& t) {
          using T = std::decay_t<decltype(t)>;
          if constexpr (std::same_as<T, AtmKeyVal>) {
            return ifatm(t);
          } else if constexpr (std::same_as<T, SurfaceKeyVal>) {
            return ifsurf(t);
          } else if constexpr (std::same_as<T, LblLineKey>) {
            return ifline(t);
          }
        },
        target);
  }

  constexpr bool operator==(const TargetType&) const = default;

  [[nodiscard]] std::string type() const {
    return apply([](auto&) { return "AtmKeyVal"s; },
                 [](auto&) { return "SurfaceKeyVal"s; },
                 [](auto&) { return "LblLineKey"s; });
  }

  friend std::ostream& operator<<(std::ostream& os, const TargetType& tt) {
    tt.apply(
        [&os](const AtmKeyVal& key) { os << "AtmKeyVal::" << key; },
        [&os](const SurfaceKeyVal& key) { os << "SurfaceKeyVal::" << key; },
        [&os](const LblLineKey& key) { os << "LblLineKey::" << key; });
    return os;
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
