#pragma once

#include <atm.h>
#include <lbl.h>
#include <matpack.h>
#include <sensor.h>
#include <subsurf.h>
#include <surf.h>
#include <xml.h>

#include <algorithm>
#include <boost/container_hash/hash.hpp>
#include <concepts>
#include <limits>
#include <variant>
#include <vector>

#include "atm_field.h"
#include "jac_log.h"
#include "jac_logrel.h"
#include "jac_polyfit.h"
#include "jac_rel.h"
#include "jac_rh.h"
#include "jacobian_names.h"
#include "xml_io_stream_functional.h"

struct ErrorKey {
  Size y_start;

  Size y_size;

  constexpr bool operator==(const ErrorKey& other) const { return other.y_start == y_start and other.y_size == y_size; }
};

template <> struct std::hash<ErrorKey> {
  static std::size_t operator()(const ErrorKey& g) {
    std::size_t seed = 0;
    boost::hash_combine(seed, g.y_start);
    boost::hash_combine(seed, g.y_size);
    return seed;
  }
};

namespace Jacobian {
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
  atm_mat inverse_jacobian{};
  atm_vec inverse_state{};

  bool      overlap{false};
  AtmKeyVal overlap_key{};

  void update_model(AtmField& y, cv x) const;
  void update_state(VectorView x, const AtmField& y) const;
  void update_jac(MatrixView dy, cv x, const AtmField& y) const;
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
  surf_mat inverse_jacobian{};
  surf_vec inverse_state{};

  bool          overlap{false};
  SurfaceKeyVal overlap_key{};

  void update_model(SurfaceField& y, cv x) const;
  void update_state(VectorView x, const SurfaceField& y) const;
  void update_jac(MatrixView dy, cv x, const SurfaceField& y) const;
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
  subsurf_mat inverse_jacobian{};
  subsurf_vec inverse_state{};

  bool             overlap{false};
  SubsurfaceKeyVal overlap_key{};

  void update_model(SubsurfaceField& y, cv x) const;
  void update_state(VectorView x, const SubsurfaceField& y) const;
  void update_jac(MatrixView dy, cv x, const SubsurfaceField& y) const;
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
  line_mat inverse_jacobian{};
  line_vec inverse_state{};

  bool       overlap{false};
  LblLineKey overlap_key{};

  void update_model(AbsorptionBands& y, cv x) const;
  void update_state(VectorView x, const AbsorptionBands& y) const;
  void update_jac(MatrixView dy, cv x, const AbsorptionBands& y) const;
};

/** The class that deals with Jacobian targets that are part of a ArrayOfSensor object
 * 
 * The intent here is that the type should be able to match into any Key that can
 * be held by a Sensor object.  This includes things like sensor frequency grid,
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
 * NOTE:  A target changes the frequency grid or the geometries of one sensor,
 * which all of its elements share.
 */
struct SensorTarget {
  SensorKey type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};

  sensor_vec transform_state{};
  sensor_mat inverse_jacobian{};
  sensor_vec inverse_state{};

  bool      overlap{false};
  SensorKey overlap_key{};

  void update_model(ArrayOfSensor& y, cv x) const;
  void update_state(VectorView x, const ArrayOfSensor& y) const;
  void update_jac(MatrixView dy, cv x, const ArrayOfSensor& y) const;
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
  error_mat inverse_jacobian{};
  error_vec inverse_state{};

  bool     overlap{false};
  ErrorKey overlap_key{};

  void update_model(VectorView y, cv x) const;
  void update_state(VectorView x, ConstVectorView y) const;
  void update_jac(MatrixView dy, cv x, cv y) const;
};

/** The class that deals with Jacobian targets that are part of a spectral field object
 *
 * The target is the data of one property of the field (AtmSpectralField,
 * SurfaceSpectralField or SubsurfaceSpectralField), and its model state is
 * the flat view of that data, so frequency is its slowest dimension.  The
 * flat weights of the data map the property at a frequency and position to
 * the model state.
 *
 * As for all targets, we need to know where in a local list the target is
 * (target_pos), where in the global list the target starts (x_start), and how
 * many global elements the target takes up (x_size).
 *
 * The targets have no transformations, so the field and the model state are
 * the same numbers.
 */
template <typename Tag, typename Field> struct SpectralTarget {
  Tag type;

  Numeric d{};

  Size target_pos{std::numeric_limits<Size>::max()};
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};

  void update_model(Field& y, cv x) const;
  void update_state(VectorView x, const Field& y) const;
};

using AtmSpectralTarget        = SpectralTarget<AtmSpectralPropertyTag, AtmSpectralField>;
using SurfaceSpectralTarget    = SpectralTarget<SurfaceSpectralPropertyTag, SurfaceSpectralField>;
using SubsurfaceSpectralTarget = SpectralTarget<SubsurfaceSpectralPropertyTag, SubsurfaceSpectralField>;

template <typename T>
concept target_type = requires(T a) {
  { a.target_pos } -> std::same_as<Size&>;
  { a.x_start } -> std::same_as<Size&>;
  { a.x_size } -> std::same_as<Size&>;
  a.type;
};

//! True if the target takes its model state from another target; targets that cannot overlap never do
template <target_type T> [[nodiscard]] bool is_overlap(const T& t) {
  if constexpr (requires { t.overlap; }) return t.overlap;
  return false;
}

//! True if the target shares its model state with other targets, either as an overlap or as the original
template <target_type T> [[nodiscard]] bool overlaps(const std::vector<T>& targets, const T& t) {
  if constexpr (requires { t.overlap_key; }) {
    return t.overlap or
           std::ranges::any_of(targets, [&t](const T& o) { return o.overlap and o.overlap_key == t.type; });
  }
  return false;
}

//! True if the target transforms between the field and the model state; targets without transformations never do
template <target_type T> [[nodiscard]] bool has_transformation(const T& t) {
  if constexpr (requires { t.transform_state; }) return t.transform_state or t.inverse_state or t.inverse_jacobian;
  return false;
}

struct Targets final {
  std::vector<AtmTarget>                atm;
  std::vector<AtmSpectralTarget>        atm_spectral;
  std::vector<SurfaceTarget>            surf;
  std::vector<SurfaceSpectralTarget>    surf_spectral;
  std::vector<SubsurfaceTarget>         subsurf;
  std::vector<SubsurfaceSpectralTarget> subsurf_spectral;
  std::vector<LineTarget>               line;
  std::vector<SensorTarget>             sensor;
  std::vector<ErrorTarget>              error;

  bool finalized{false};

  //! The number of targets when finalize was last called
  Size finalized_count{0};

  //! Calls f on each list of targets, in model state order
  void for_each_list(this auto&& self, auto&& f) {
    f(self.atm);
    f(self.atm_spectral);
    f(self.surf);
    f(self.surf_spectral);
    f(self.subsurf);
    f(self.subsurf_spectral);
    f(self.line);
    f(self.sensor);
    f(self.error);
  }

  //! True if finalize has been called and no target has been added since
  [[nodiscard]] bool is_finalized() const;

  [[nodiscard]] Size x_size() const;

  [[nodiscard]] Size target_count() const;

  void throwing_check(Size xsize) const;

  //! Sets the sizes and x-positions of the targets.
  void finalize(const AtmField&                atm_field,
                const AtmSpectralField&        atm_spectral_field,
                const SurfaceField&            surf_field,
                const SurfaceSpectralField&    surf_spectral_field,
                const SubsurfaceField&         subsurf_field,
                const SubsurfaceSpectralField& subsurf_spectral_field,
                const AbsorptionBands&         abs_bands,
                const ArrayOfSensor&           measurement_sensors);

  AtmTarget&                emplace_back(AtmKeyVal&& t, Numeric d = 0.0);
  AtmTarget&                emplace_back(const AtmKeyVal& t, Numeric d = 0.0);
  LineTarget&               emplace_back(LblLineKey&& t, Numeric d = 0.0);
  LineTarget&               emplace_back(const LblLineKey& t, Numeric d = 0.0);
  SensorTarget&             emplace_back(SensorKey&& t, Numeric d = 0.0);
  SensorTarget&             emplace_back(const SensorKey& t, Numeric d = 0.0);
  SurfaceTarget&            emplace_back(SurfaceKeyVal&& t, Numeric d = 0.0);
  SurfaceTarget&            emplace_back(const SurfaceKeyVal& t, Numeric d = 0.0);
  SubsurfaceTarget&         emplace_back(SubsurfaceKeyVal&& t, Numeric d = 0.0);
  SubsurfaceTarget&         emplace_back(const SubsurfaceKeyVal& t, Numeric d = 0.0);
  AtmSpectralTarget&        emplace_back(const AtmSpectralPropertyTag& t, Numeric d = 0.0);
  SurfaceSpectralTarget&    emplace_back(const SurfaceSpectralPropertyTag& t, Numeric d = 0.0);
  SubsurfaceSpectralTarget& emplace_back(const SubsurfaceSpectralPropertyTag& t, Numeric d = 0.0);
  ErrorTarget&              emplace_back(ErrorKey&& t, Numeric d = 0.0);
  ErrorTarget&              emplace_back(const ErrorKey& t, Numeric d = 0.0);

  [[nodiscard]] std::vector<AtmTarget>::const_iterator                find(const AtmKeyVal& t) const;
  [[nodiscard]] std::vector<SurfaceTarget>::const_iterator            find(const SurfaceKeyVal& t) const;
  [[nodiscard]] std::vector<SubsurfaceTarget>::const_iterator         find(const SubsurfaceKeyVal& t) const;
  [[nodiscard]] std::vector<AtmSpectralTarget>::const_iterator        find(const AtmSpectralPropertyTag& t) const;
  [[nodiscard]] std::vector<SurfaceSpectralTarget>::const_iterator    find(const SurfaceSpectralPropertyTag& t) const;
  [[nodiscard]] std::vector<SubsurfaceSpectralTarget>::const_iterator find(
      const SubsurfaceSpectralPropertyTag& t) const;
  [[nodiscard]] std::vector<LineTarget>::const_iterator   find(const LblLineKey& t) const;
  [[nodiscard]] std::vector<SensorTarget>::const_iterator find(const SensorKey& t) const;
  [[nodiscard]] std::vector<ErrorTarget>::const_iterator  find(const ErrorKey& t) const;

  [[nodiscard]] Index target_position(const AtmKeyVal& t) const;
  [[nodiscard]] Index target_position(const SurfaceKeyVal& t) const;
  [[nodiscard]] Index target_position(const SubsurfaceKeyVal& t) const;
  [[nodiscard]] Index target_position(const AtmSpectralPropertyTag& t) const;
  [[nodiscard]] Index target_position(const SurfaceSpectralPropertyTag& t) const;
  [[nodiscard]] Index target_position(const SubsurfaceSpectralPropertyTag& t) const;
  [[nodiscard]] Index target_position(const LblLineKey& t) const;
  [[nodiscard]] Index target_position(const SensorKey& t) const;
  [[nodiscard]] Index target_position(const ErrorKey& t) const;

  void               clear();
  [[nodiscard]] bool empty() const;
};

struct TargetType {
  using variant_t = std::variant<AtmKeyVal,
                                 AtmSpectralPropertyTag,
                                 SurfaceKeyVal,
                                 SurfaceSpectralPropertyTag,
                                 SubsurfaceKeyVal,
                                 SubsurfaceSpectralPropertyTag,
                                 LblLineKey,
                                 SensorKey,
                                 ErrorKey>;
  variant_t target;

  template <class AtmKeyValFunc,
            class AtmSpectralFunc,
            class SurfaceKeyValFunc,
            class SurfaceSpectralFunc,
            class SubsurfaceKeyValFunc,
            class SubsurfaceSpectralFunc,
            class LblLineKeyFunc,
            class SensorKeyFunc,
            class ErrorKeyFunc>
  [[nodiscard]] constexpr auto apply(const AtmKeyValFunc&          ifatm,
                                     const AtmSpectralFunc&        ifatm_spectral,
                                     const SurfaceKeyValFunc&      ifsurf,
                                     const SurfaceSpectralFunc&    ifsurf_spectral,
                                     const SubsurfaceKeyValFunc&   ifsubsurf,
                                     const SubsurfaceSpectralFunc& ifsubsurf_spectral,
                                     const LblLineKeyFunc&         ifline,
                                     const SensorKeyFunc&          ifsensor,
                                     const ErrorKeyFunc&           iferror) const {
    return std::visit(
        [&](const auto& t) {
          using T = std::decay_t<decltype(t)>;
          if constexpr (std::same_as<T, AtmKeyVal>) {
            return ifatm(t);
          } else if constexpr (std::same_as<T, AtmSpectralPropertyTag>) {
            return ifatm_spectral(t);
          } else if constexpr (std::same_as<T, SurfaceKeyVal>) {
            return ifsurf(t);
          } else if constexpr (std::same_as<T, SurfaceSpectralPropertyTag>) {
            return ifsurf_spectral(t);
          } else if constexpr (std::same_as<T, SubsurfaceKeyVal>) {
            return ifsubsurf(t);
          } else if constexpr (std::same_as<T, SubsurfaceSpectralPropertyTag>) {
            return ifsubsurf_spectral(t);
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

  [[nodiscard]] bool        operator==(const TargetType&) const;
  [[nodiscard]] std::string type() const;
};
}  // namespace Jacobian

using JacobianTargets    = Jacobian::Targets;
using JacobianTargetType = Jacobian::TargetType;

template <> struct std::formatter<ErrorKey> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const ErrorKey& v, FmtContext& ctx) const {
    return tags.format(ctx, "["sv, v.y_start, ", "sv, v.y_start + v.y_size, ")"sv);
  }
};

template <> struct std::hash<JacobianTargetType> {
  static std::size_t operator()(const JacobianTargetType& g) {
    return std::hash<JacobianTargetType::variant_t>{}(g.target);
  }
};

template <> struct std::formatter<JacobianTargetType> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const JacobianTargetType& jtt, FmtContext& ctx) const {
    jtt.apply([this, &ctx](const AtmKeyVal& key) { tags.format(ctx, "AtmKeyVal::"sv, key); },
              [this, &ctx](const AtmSpectralPropertyTag& key) { tags.format(ctx, "AtmSpectralPropertyTag::"sv, key); },
              [this, &ctx](const SurfaceKeyVal& key) { tags.format(ctx, "SurfaceKeyVal::"sv, key); },
              [this, &ctx](const SurfaceSpectralPropertyTag& key) {
                tags.format(ctx, "SurfaceSpectralPropertyTag::"sv, key);
              },
              [this, &ctx](const SubsurfaceKeyVal& key) { tags.format(ctx, "SubsurfaceKeyVal::"sv, key); },
              [this, &ctx](const SubsurfaceSpectralPropertyTag& key) {
                tags.format(ctx, "SubsurfaceSpectralPropertyTag::"sv, key);
              },
              [this, &ctx](const LblLineKey& key) { tags.format(ctx, "LblLineKey::"sv, key); },
              [this, &ctx](const SensorKey& key) { tags.format(ctx, "SensorKey::"sv, key); },
              [this, &ctx](const ErrorKey& key) { tags.format(ctx, "ErrorKey::"sv, key); });

    return ctx.out();
  }
};

template <Jacobian::target_type T> struct std::formatter<T> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const T& v, FmtContext& ctx) const {
    if (tags.short_str) return tags.format(ctx, v.type);

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

template <> struct std::formatter<JacobianTargets> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  template <class FmtContext> FmtContext::iterator format(const JacobianTargets& v, FmtContext& ctx) const {
    tags.add_if_bracket(ctx, "{"sv);

    const std::string_view sep = tags.sep();

    if (not v.atm.empty()) { tags.format(ctx, R"("atm": )"sv, v.atm, sep); }

    if (not v.atm_spectral.empty()) { tags.format(ctx, R"("atm_spectral": )"sv, v.atm_spectral, sep); }

    if (not v.surf.empty()) { tags.format(ctx, R"("surf": )"sv, v.surf, sep); }

    if (not v.surf_spectral.empty()) { tags.format(ctx, R"("surf_spectral": )"sv, v.surf_spectral, sep); }

    if (not v.subsurf.empty()) { tags.format(ctx, R"("subsurf": )"sv, v.subsurf, sep); }

    if (not v.subsurf_spectral.empty()) { tags.format(ctx, R"("subsurf_spectral": )"sv, v.subsurf_spectral, sep); }

    if (not v.line.empty()) tags.format(ctx, R"("line": )"sv, v.line, sep);

    if (not v.sensor.empty()) { tags.format(ctx, R"("sensor": )"sv, v.sensor, sep); }

    if (not v.error.empty()) { tags.format(ctx, R"("error": )"sv, v.error); }

    tags.add_if_bracket(ctx, "}"sv);
    return ctx.out();
  }
};

template <> struct xml_io_stream<JacobianTargetType> {
  static constexpr std::string_view type_name = "JacobianTargetType"sv;

  static void write(std::ostream&             os,
                    const JacobianTargetType& x,
                    bofstream*                pbofs = nullptr,
                    std::string_view          name  = ""sv);

  static void read(std::istream& is, JacobianTargetType& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream_name<Jacobian::AtmTarget> {
  static constexpr std::string_view name = "AtmTarget"sv;
};

template <> struct xml_io_stream_name<Jacobian::LineTarget> {
  static constexpr std::string_view name = "LineTarget"sv;
};

template <> struct xml_io_stream_name<Jacobian::SurfaceTarget> {
  static constexpr std::string_view name = "SurfaceTarget"sv;
};

template <> struct xml_io_stream_name<Jacobian::SubsurfaceTarget> {
  static constexpr std::string_view name = "SubsurfaceTarget"sv;
};

template <> struct xml_io_stream_name<Jacobian::AtmSpectralTarget> {
  static constexpr std::string_view name = "AtmSpectralTarget"sv;
};

template <> struct xml_io_stream_name<Jacobian::SurfaceSpectralTarget> {
  static constexpr std::string_view name = "SurfaceSpectralTarget"sv;
};

template <> struct xml_io_stream_name<Jacobian::SubsurfaceSpectralTarget> {
  static constexpr std::string_view name = "SubsurfaceSpectralTarget"sv;
};

template <> struct xml_io_stream_name<Jacobian::SensorTarget> {
  static constexpr std::string_view name = "SensorTarget"sv;
};

template <> struct xml_io_stream_name<ErrorKey> {
  static constexpr std::string_view name = "ErrorKey"sv;
};

template <> struct xml_io_stream_name<Jacobian::ErrorTarget> {
  static constexpr std::string_view name = "ErrorTarget"sv;
};

template <> struct xml_io_stream_aggregate<Jacobian::AtmTarget> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_aggregate<Jacobian::LineTarget> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_aggregate<Jacobian::SurfaceTarget> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_aggregate<Jacobian::SubsurfaceTarget> {
  static constexpr bool value = true;
};

template <typename Tag, typename Field> struct xml_io_stream_aggregate<Jacobian::SpectralTarget<Tag, Field>> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_aggregate<Jacobian::SensorTarget> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_aggregate<ErrorKey> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_aggregate<Jacobian::ErrorTarget> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream<JacobianTargets> {
  static constexpr std::string_view type_name = "JacobianTargets"sv;

  static void write(std::ostream&          os,
                    const JacobianTargets& x,
                    bofstream*             pbofs = nullptr,
                    std::string_view       name  = ""sv);

  static void read(std::istream& is, JacobianTargets& x, bifstream* pbifs = nullptr);
};
