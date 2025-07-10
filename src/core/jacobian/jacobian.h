#pragma once

#include <atm.h>
#include <lbl.h>
#include <matpack.h>
#include <obsel.h>
#include <operators.h>
#include <subsurface.h>
#include <surf.h>
#include <xml.h>

#include <concepts>
#include <limits>
#include <numeric>
#include <variant>
#include <vector>

#include "lbl_data.h"
#include "matpack_mdspan_view_t.h"

struct ErrorKey {
  Size y_start;

  Size y_size;

  constexpr bool operator==(const ErrorKey& other) const {
    return other.y_start == y_start and other.y_size == y_size;
  }
};

template <>
struct std::hash<ErrorKey> {
  std::size_t operator()(const ErrorKey& g) const {
    return std::hash<Size>{}(g.y_start) ^ (std::hash<Size>{}(g.y_size) << 32);
  }
};

namespace Jacobian {
using atm_vec     = std::function<Vector(Vector, const AtmField&)>;
using atm_mat     = std::function<Matrix(Matrix, const AtmField&)>;
using surf_vec    = std::function<Vector(Vector, const SurfaceField&)>;
using surf_mat    = std::function<Matrix(Matrix, const SurfaceField&)>;
using subsurf_vec = std::function<Vector(Vector, const SubsurfaceField&)>;
using subsurf_mat = std::function<Matrix(Matrix, const SubsurfaceField&)>;
using line_vec    = std::function<Vector(Vector, const AbsorptionBands&)>;
using line_mat    = std::function<Matrix(Matrix, const AbsorptionBands&)>;
using sensor_vec  = std::function<Vector(Vector, const ArrayOfSensorObsel&)>;
using sensor_mat  = std::function<Matrix(Matrix, const ArrayOfSensorObsel&)>;
using error_vec   = std::function<Vector(Vector, const ConstVectorView)>;
using error_mat   = std::function<Matrix(Matrix, const ConstVectorView)>;

/** The class that deals with Jacobian targets that are part of an AtmField object
 * 
 * The intent here is that the type should be able to match into any Key that can
 * be held by an AtmField object.  This includes things like atmospheric temperature,
 * pressure, VMR, Scattering species, etc.
 * 
 * As for all targets, we need to know where in a local list the target is
 * (target_pos), where in the global list the target starts (x_start), and how
 * many global elements the target takes up (x_size).
 * 
 * The set_state and set_model functions are used to update the target in
 * iteration and model space, respectively.  The default functions are provided
 * but can be overridden with advantage as they perform a lot of computations
 * that are not required when the type is known.
 */
struct AtmTarget {
  AtmKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};

  atm_vec transform_state{};
  atm_mat transform_jacobian{};
  atm_vec inverse_state{};

  void update_model(AtmField& atm, const ConstVectorView x) const;
  void update_state(VectorView x, const AtmField& atm) const;
  void update_jac(MatrixView x, const AtmField& atm) const;

  [[nodiscard]] bool is_wind() const;
};

/** The class that deals with Jacobian targets that are part of a SurfaceField object
 * 
 * The intent here is that the type should be able to match into any Key that can
 * be held by a SurfaceField object.  This includes things like surface temperature,
 * surface properties, etc.
 * 
 * As for all targets, we need to know where in a local vector the target is
 * (target_pos), where in the global vector the target starts (x_start), and how
 * many global elements the target takes up (x_size).
 * 
 * The set_state and set_model functions are used to update the target in
 * iteration and model space, respectively.  The default functions are provided
 * but can be overridden with advantage as they perform a lot of computations
 * that are not required when the type is known.
 */
struct SurfaceTarget {
  SurfaceKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};

  surf_vec transform_state{};
  surf_mat transform_jacobian{};
  surf_vec inverse_state{};

  void update_model(SurfaceField& surf, const ConstVectorView x) const;
  void update_state(VectorView x, const SurfaceField& surf) const;
  void update_jac(MatrixView x, const SurfaceField& surf) const;
};

/** The class that deals with Jacobian targets that are part of a SubsurfaceField object
 * 
 * The intent here is that the type should be able to match into any Key that can
 * be held by an SubsurfaceField object.  This includes things like subsurface temperature,
 * density, etc.
 * 
 * As for all targets, we need to know where in a local list the target is
 * (target_pos), where in the global list the target starts (x_start), and how
 * many global elements the target takes up (x_size).
 * 
 * The set_state and set_model functions are used to update the target in
 * iteration and model space, respectively.  The default functions are provided
 * but can be overridden with advantage as they perform a lot of computations
 * that are not required when the type is known.
 */
struct SubsurfaceTarget {
  SubsurfaceKeyVal type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};

  subsurf_vec transform_state{};
  subsurf_mat transform_jacobian{};
  subsurf_vec inverse_state{};

  void update_model(SubsurfaceField& surf, const ConstVectorView x) const;
  void update_state(VectorView x, const SubsurfaceField& surf) const;
  void update_jac(MatrixView x, const SubsurfaceField& surf) const;
};

/** The class that deals with Jacobian targets that are part of an AbsorptionBands object
 * 
 * The intent here is that the type should be able to match into any Key that can
 * be held by a AbsorptionBands object.  This includes things like einstein coefficient
 * pressure parameters, etc.
 * 
 * As for all targets, we need to know where in a local vector the target is
 * (target_pos), where in the global vector the target starts (x_start), and how
 * many global elements the target takes up (x_size).
 * 
 * The set_state and set_model functions are used to update the target in
 * iteration and model space, respectively.  The default functions are provided
 * but can be overridden with advantage as they perform a lot of computations
 * that are not required when the type is known.
 */
struct LineTarget {
  LblLineKey type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};

  line_vec transform_state{};
  line_mat transform_jacobian{};
  line_vec inverse_state{};

  void update_model(AbsorptionBands& surf, const ConstVectorView x) const;
  void update_state(VectorView x, const AbsorptionBands& surf) const;
  void update_jac(MatrixView x, const AbsorptionBands& surf) const;
};

/** The class that deals with Jacobian targets that are part of a ArrayOfSensorObsel object
 * 
 * The intent here is that the type should be able to match into any Key that can
 * be held by a SensorObsel object.  This includes things like sensor frequency grid,
 * sensor position, and sensor line of sight.
 * 
 * As for all targets, we need to know where in a local vector the target is
 * (target_pos), where in the global vector the target starts (x_start), and how
 * many global elements the target takes up (x_size).
 * 
 * The set_state and set_model functions are used to update the target in
 * iteration and model space, respectively.  The default functions are provided
 * but can be overridden with advantage as they perform a lot of computations
 * that are not required when the type is known.
 * 
 * NOTE:  Several SensorObsel objects share grids.  By default, shared grids
 * are updated together.
 */
struct SensorTarget {
  SensorKey type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};

  sensor_vec transform_state{};
  sensor_mat transform_jacobian{};
  sensor_vec inverse_state{};

  void update_model(ArrayOfSensorObsel& surf, const ConstVectorView x) const;
  void update_state(VectorView x, const ArrayOfSensorObsel& surf) const;
  void update_jac(MatrixView x, const ArrayOfSensorObsel& surf) const;
};

/** The class that deals with Jacobian targets that are not part of any input object
 * 
 * The intent here is that some properties of a retrieval process are best described
 * as an error on the measurement vector rather than as a property of the measurement
 * 
 * As for all targets, we need to know where in a local vector the target is
 * (target_pos), where in the global vector the target starts (x_start), and how
 * many global elements the target takes up (x_size).
 * 
 * The set_state and set_model functions are used to update the target in
 * iteration and model space, respectively.  Note that the default functions
 * will just throw an error if called.  It is necessary to override these functions
 * manually for each case.
 */
struct ErrorTarget {
  ErrorKey type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};

  error_vec transform_state{};
  error_mat transform_jacobian{};
  error_vec inverse_state{};

  void update_model(VectorView m, const ConstVectorView x) const;
  void update_state(VectorView x, const ConstVectorView m) const;
  void update_jac(MatrixView x, const ConstVectorView m) const;
};

template <typename T>
concept target_type = requires(T a) {
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

  using tup_t = std::tuple<std::vector<Targets>...>;
  tup_t targets{};

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
    if (target_count() != 0 and not finalized)
      throw std::runtime_error("Not finalized.");

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

    if (xsize != x_size())
      throw std::runtime_error(
          "The size of the x-vector does not match the size of the targets.");

    ((std::ranges::for_each(
         target<Targets>(),
         [xsize, t_size](auto& a) {
           if ((a.x_start + a.x_size) > xsize)
             throw std::runtime_error("x-vector out-of-bounds");
           if (t_size <= a.target_pos)
             throw std::runtime_error("target-vector out-of-bounds.");
         })),
     ...);
  }

  void clear() { ((target<Targets>().clear()), ...); }
};

struct Targets final : targets_t<AtmTarget,
                                 SurfaceTarget,
                                 SubsurfaceTarget,
                                 LineTarget,
                                 SensorTarget,
                                 ErrorTarget> {
  [[nodiscard]] const std::vector<AtmTarget>& atm() const;
  [[nodiscard]] const std::vector<SurfaceTarget>& surf() const;
  [[nodiscard]] const std::vector<SubsurfaceTarget>& subsurf() const;
  [[nodiscard]] const std::vector<LineTarget>& line() const;
  [[nodiscard]] const std::vector<SensorTarget>& sensor() const;
  [[nodiscard]] const std::vector<ErrorTarget>& error() const;

  [[nodiscard]] std::vector<AtmTarget>& atm();
  [[nodiscard]] std::vector<SurfaceTarget>& surf();
  [[nodiscard]] std::vector<SubsurfaceTarget>& subsurf();
  [[nodiscard]] std::vector<LineTarget>& line();
  [[nodiscard]] std::vector<SensorTarget>& sensor();
  [[nodiscard]] std::vector<ErrorTarget>& error();

  //! Sets the sizes and x-positions of the targets.
  void finalize(const AtmField& atmospheric_field,
                const SurfaceField& surface_field,
                const SubsurfaceField& subsurface_field,
                const AbsorptionBands& absorption_bands,
                const ArrayOfSensorObsel& measurement_sensor);

  AtmTarget& emplace_back(AtmKeyVal&& t, Numeric d);
  AtmTarget& emplace_back(const AtmKeyVal& t, Numeric d);
  LineTarget& emplace_back(LblLineKey&& t, Numeric d);
  LineTarget& emplace_back(const LblLineKey& t, Numeric d);
  SensorTarget& emplace_back(SensorKey&& t, Numeric d);
  SensorTarget& emplace_back(const SensorKey& t, Numeric d);
  SurfaceTarget& emplace_back(SurfaceKeyVal&& t, Numeric d);
  SurfaceTarget& emplace_back(const SurfaceKeyVal& t, Numeric d);
  SubsurfaceTarget& emplace_back(SubsurfaceKeyVal&& t, Numeric d);
  SubsurfaceTarget& emplace_back(const SubsurfaceKeyVal& t, Numeric d);
  ErrorTarget& emplace_back(ErrorKey&& t, Numeric d);
  ErrorTarget& emplace_back(const ErrorKey& t, Numeric d);
};

struct TargetType {
  using variant_t = std::variant<AtmKeyVal,
                                 SurfaceKeyVal,
                                 SubsurfaceKeyVal,
                                 LblLineKey,
                                 SensorKey,
                                 ErrorKey>;
  variant_t target;

  template <class AtmKeyValFunc,
            class SurfaceKeyValFunc,
            class SubsurfaceKeyValFunc,
            class LblLineKeyFunc,
            class SensorKeyFunc,
            class ErrorKeyFunc>
  [[nodiscard]] constexpr auto apply(const AtmKeyValFunc& ifatm,
                                     const SurfaceKeyValFunc& ifsurf,
                                     const SubsurfaceKeyValFunc& ifsubsurf,
                                     const LblLineKeyFunc& ifline,
                                     const SensorKeyFunc& ifsensor,
                                     const ErrorKeyFunc& iferror) const {
    return std::visit(
        [&](const auto& t) {
          using T = std::decay_t<decltype(t)>;
          if constexpr (std::same_as<T, AtmKeyVal>) {
            return ifatm(t);
          } else if constexpr (std::same_as<T, SurfaceKeyVal>) {
            return ifsurf(t);
          } else if constexpr (std::same_as<T, SubsurfaceKeyVal>) {
            return ifsubsurf(t);
          } else if constexpr (std::same_as<T, LblLineKey>) {
            return ifline(t);
          } else if constexpr (std::same_as<T, SensorKey>) {
            return ifsensor(t);
          } else if constexpr (std::same_as<T, ErrorKey>) {
            return iferror(t);
          }
        },
        target);
  }

  constexpr bool operator==(const TargetType&) const = default;

  [[nodiscard]] std::string type() const {
    return apply([](auto&) { return "AtmKeyVal"s; },
                 [](auto&) { return "SurfaceKeyVal"s; },
                 [](auto&) { return "SubsurfaceKeyVal"s; },
                 [](auto&) { return "LblLineKey"s; },
                 [](auto&) { return "SensorKey"s; },
                 [](auto&) { return "ErrorKey"s; });
  }
};

void polyfit(VectorView param,
             const ConstVectorView x,
             const ConstVectorView y);
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
struct std::formatter<ErrorKey> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ErrorKey& v, FmtContext& ctx) const {
    return tags.format(ctx, '[', v.y_start, ", "sv, v.y_start + v.y_size, ')');
  }
};

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
        [this, &ctx](const SubsurfaceKeyVal& key) {
          tags.format(ctx, "SubsurfaceKeyVal::"sv, key);
        },
        [this, &ctx](const LblLineKey& key) {
          tags.format(ctx, "LblLineKey::"sv, key);
        },
        [this, &ctx](const SensorKey& key) {
          tags.format(ctx, "SensorKey::"sv, key);
        },
        [this, &ctx](const ErrorKey& key) {
          tags.format(ctx, "ErrorKey::"sv, key);
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
struct std::formatter<Jacobian::ErrorTarget> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Jacobian::ErrorTarget& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(ctx,
                       v.type,
                       ": target_pos: "sv,
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
                v.line(),
                sep,
                R"("sensor": )"sv,
                v.sensor(),
                sep,
                R"("error": )"sv,
                v.error());

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};

template <>
struct xml_io_stream<JacobianTargetType> {
  static constexpr std::string_view type_name = "JacobianTargetType"sv;

  static void write(std::ostream& os,
                    const JacobianTargetType& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   JacobianTargetType& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<Jacobian::AtmTarget> {
  static constexpr std::string_view type_name = "AtmTarget"sv;

  static void write(std::ostream& os,
                    const Jacobian::AtmTarget& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   Jacobian::AtmTarget& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<Jacobian::LineTarget> {
  static constexpr std::string_view type_name = "LineTarget"sv;

  static void write(std::ostream& os,
                    const Jacobian::LineTarget& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   Jacobian::LineTarget& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<Jacobian::SurfaceTarget> {
  static constexpr std::string_view type_name = "SurfaceTarget"sv;

  static void write(std::ostream& os,
                    const Jacobian::SurfaceTarget& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   Jacobian::SurfaceTarget& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<Jacobian::SubsurfaceTarget> {
  static constexpr std::string_view type_name = "SubsurfaceTarget"sv;

  static void write(std::ostream& os,
                    const Jacobian::SubsurfaceTarget& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   Jacobian::SubsurfaceTarget& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<Jacobian::SensorTarget> {
  static constexpr std::string_view type_name = "SensorTarget"sv;

  static void write(std::ostream& os,
                    const Jacobian::SensorTarget& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   Jacobian::SensorTarget& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<ErrorKey> {
  static constexpr std::string_view type_name = "ErrorKey"sv;

  static void write(std::ostream& os,
                    const ErrorKey& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, ErrorKey& x, bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<Jacobian::ErrorTarget> {
  static constexpr std::string_view type_name = "ErrorTarget"sv;

  static void write(std::ostream& os,
                    const Jacobian::ErrorTarget& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   Jacobian::ErrorTarget& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<JacobianTargets> {
  static constexpr std::string_view type_name = "JacobianTargets"sv;

  static void write(std::ostream& os,
                    const JacobianTargets& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   JacobianTargets& x,
                   bifstream* pbifs = nullptr);
};
