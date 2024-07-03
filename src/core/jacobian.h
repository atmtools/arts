#pragma once

#include <atm.h>
#include <lbl.h>
#include <matpack.h>
#include <surf.h>

#include <functional>
#include <limits>
#include <numeric>
#include <vector>

#include "array.h"
#include "debug.h"

namespace Jacobian {
struct AtmTarget {
  AtmKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  std::function<void(
      ExhaustiveVectorView, const AtmField&, const ExhaustiveConstVectorView&)>
      set{[](ExhaustiveVectorView xnew,
             const AtmField&,
             const ConstVectorView& xold) { xnew = xold; }};

  std::function<void(
      ExhaustiveVectorView, const AtmField&, const ExhaustiveConstVectorView&)>
      unset{[](VectorView xnew,
               const AtmField&,
               const ExhaustiveConstVectorView& xold) { xnew = xold; }};

  friend std::ostream& operator<<(std::ostream& os, const AtmTarget& target) {
    return os << "Atmosphere key value: " << target.type << ", starting at "
              << target.x_start << " of size " << target.x_size;
  }

  void update(AtmField& atm, const Vector& x) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (x_start + x_size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not atm.contains(type),
                       "Atmosphere does not contain key value ",
                       type,
                       '.')

    auto xnew   = atm[type].flat_view();
    auto xold_d = x.slice(x_start, x_size);
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.  \n"
        "Did you change your atmosphere since you set the jacobian targets?")

    set(xnew, atm, xold_d);
  }

  void update(Vector& x, const AtmField& atm) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (x_start + x_size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not atm.contains(type),
                       "Atmosphere does not contain key value ",
                       type,
                       '.')

    auto xnew   = x.slice(x_start, x_size);
    auto xold_d = atm[type].flat_view();
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.  \n"
        "Did you change your atmosphere since you set the jacobian targets?\n"
        "Did you forget to finalize the JacobianTargets?")

    unset(xnew, atm, xold_d);
  }
};

struct SurfaceTarget {
  SurfaceKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  std::function<void(ExhaustiveVectorView,
                     const SurfaceField&,
                     const ExhaustiveConstVectorView&)>
      set{[](ExhaustiveVectorView xnew,
             const SurfaceField&,
             const ConstVectorView& xold) { xnew = xold; }};

  std::function<void(ExhaustiveVectorView,
                     const SurfaceField&,
                     const ExhaustiveConstVectorView&)>
      unset{[](VectorView xnew,
               const SurfaceField&,
               const ExhaustiveConstVectorView& xold) { xnew = xold; }};

  friend std::ostream& operator<<(std::ostream& os,
                                  const SurfaceTarget& target) {
    return os << "Surface key value: " << target.type << ", starting at "
              << target.x_start << " of size " << target.x_size;
  }

  void update(SurfaceField& surf, const Vector& x) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (x_start + x_size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not surf.contains(type),
                       "Surface does not contain key value ",
                       type,
                       '.')

    auto xnew   = surf[type].flat_view();
    auto xold_d = x.slice(x_start, x_size);
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.\n"
        "Did you change your surface since you set the jacobian targets?\n"
        "Did you forget to finalize the JacobianTargets?")

    set(xnew, surf, xold_d);
  }

  void update(Vector& x, const SurfaceField& surf) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (x_start + x_size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not surf.contains(type),
                       "Surface does not contain key value ",
                       type,
                       '.')

    auto xnew   = x.slice(x_start, x_size);
    auto xold_d = surf[type].flat_view();
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.\n"
        "Did you change your surface since you set the jacobian targets?")

    unset(xnew, surf, xold_d);
  }
};

struct LineTarget {
  LblLineKey type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};

  Size x_start{std::numeric_limits<Size>::max()};

  Size x_size{std::numeric_limits<Size>::max()};

  std::function<void(ExhaustiveVectorView,
                     const ArrayOfAbsorptionBand&,
                     const ExhaustiveConstVectorView&)>
      set{[](ExhaustiveVectorView xnew,
             const ArrayOfAbsorptionBand&,
             const ConstVectorView& xold) { xnew = xold; }};

  std::function<void(ExhaustiveVectorView,
                     const ArrayOfAbsorptionBand&,
                     const ExhaustiveConstVectorView&)>
      unset{[](VectorView xnew,
               const ArrayOfAbsorptionBand&,
               const ExhaustiveConstVectorView& xold) { xnew = xold; }};

  friend std::ostream& operator<<(std::ostream& os, const LineTarget&) {
    return os << "Line key value: ";
  }

  void update(ArrayOfAbsorptionBand&, const Vector&) const {
    ARTS_ASSERT(false)
  }

  void update(Vector&, const ArrayOfAbsorptionBand&) const {
    ARTS_ASSERT(false)
  }
};

template <typename U, typename T>
concept target_comparable = requires(T a, U b) {
  a.type == b;
  b == a.type;
};

template <typename U, typename... T>
concept valid_target = (std::same_as<U, T> or ...);

template <typename... Targets>
struct targets_t {
  static constexpr Size N = sizeof...(Targets);

  std::tuple<std::vector<Targets>...> targets{};

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
           ARTS_USER_ERROR_IF(a.x_start + a.x_size >= xsize,
                              "The target ",
                              a,
                              " is out of bounds of the x-vector.")
           ARTS_USER_ERROR_IF(t_size >= a.target_pos,
                              "The target ",
                              a,
                              " is out of bounds of the target vector.")
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
  [[nodiscard]] const auto& atm() const { return target<AtmTarget>(); }
  [[nodiscard]] const auto& surf() const { return target<SurfaceTarget>(); }
  [[nodiscard]] const auto& line() const { return target<LineTarget>(); }

  //! Sets the sizes and x-positions of the targets.
  void finalize(const AtmField& atm_field) {
    zero_out_x();

    const Size Nt = target_count();
    for (Size i = 0; i < Nt; i++) {
      if ([&atm_field, this, pos = i]() -> bool {
            for (auto& t : target<AtmTarget>()) {
              if (t.target_pos == pos) {
                ARTS_USER_ERROR_IF(not atm_field.contains(t.type),
                                   "The target ",
                                   t,
                                   " is not in the atmosphere,"
                                   " but is required by jacobian target ",
                                   this)
                t.x_start = x_size();
                t.x_size  = atm_field[t.type].flat_view().size();
                return true;
              }
            }
            return false;
          }())
        continue;

      ARTS_ASSERT(false, "Missing impl for surface and line parameters")
      ARTS_USER_ERROR("Cannot reach here, there is no target position ", i)
    }

    throwing_check(x_size());
  }

  friend std::ostream& operator<<(std::ostream& os, const Targets& targets) {
    os << "Jacobian targets:\n";
    for (const auto& t : targets.atm()) {
      os << "  " << t << '\n';
    }
    for (const auto& t : targets.surf()) {
      os << "  " << t << '\n';
    }
    for (const auto& t : targets.line()) {
      os << "  " << t << '\n';
    }
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

using JacobianTargets = Jacobian::Targets;

template <>
struct std::formatter<Jacobian::AtmTarget> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Jacobian::AtmTarget& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();

    std::formatter<AtmKeyVal> type{};
    make_compat(type);
    type.format(v.type, ctx);
    std::format_to(ctx.out(),
                   ": [{}{}{}{}{}{}{}]",
                   v.d,
                   sep,
                   v.target_pos,
                   sep,
                   v.x_start,
                   sep,
                   v.x_size);

    return ctx.out();
  }
};

template <>
struct std::formatter<Jacobian::SurfaceTarget> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Jacobian::SurfaceTarget& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();

    std::formatter<SurfaceKeyVal> type{};
    make_compat(type);
    type.format(v.type, ctx);
    std::format_to(ctx.out(),
                   ": [{}{}{}{}{}{}{}]",
                   v.d,
                   sep,
                   v.target_pos,
                   sep,
                   v.x_start,
                   sep,
                   v.x_size);

    return ctx.out();
  }
};

template <>
struct std::formatter<Jacobian::LineTarget> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Jacobian::LineTarget& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();

    std::formatter<LblLineKey> type{};
    make_compat(type);
    type.format(v.type, ctx);
    std::format_to(ctx.out(),
                   ": [{}{}{}{}{}{}{}]",
                   v.d,
                   sep,
                   v.target_pos,
                   sep,
                   v.x_start,
                   sep,
                   v.x_size);

    return ctx.out();
  }
};

template <>
struct std::formatter<JacobianTargets> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const JacobianTargets& v, FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '{');

    const std::string_view sep = tags.sep(true);

    std::formatter<std::vector<Jacobian::AtmTarget>> atm{};
    std::formatter<std::vector<Jacobian::SurfaceTarget>> surf{};
    std::formatter<std::vector<Jacobian::LineTarget>> line{};

    std::format_to(ctx.out(), R"("atm": )");
    atm.format(v.atm(), ctx);
    std::format_to(ctx.out(), R"({}"surf": )", sep);
    surf.format(v.surf(), ctx);
    std::format_to(ctx.out(), R"({}"line": )", sep);
    line.format(v.line(), ctx);

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};
