#include "lbl_jpl.h"

#include <arts_conversions.h>
#include <fast_float/fast_float.h>
#include <jpl_species.h>
#include <partfun.h>
#include <quantum.h>

#include <charconv>
#include <istream>

#include "enumsLineShapeModelType.h"

namespace lbl {
namespace {
struct reader {
  std::string::const_iterator it;
  std::string::const_iterator end;

  reader(const std::string& s) : it(s.begin()), end(s.end()) {}

  template <typename T>
  constexpr T read_next(Size n, std::string_view error_context) {
    std::string_view orig(it, it + n);
    std::string_view sv = orig;
    skip(n);

    if constexpr (std::same_as<T, char>) {
      assert(n == 1);
      return orig[0];
    } else {
      while (sv.size() and sv.front() == ' ') sv.remove_prefix(1);
      while (sv.size() and sv.back() == ' ') sv.remove_suffix(1);

      T x{};
      if constexpr (std::same_as<T, double> or std::same_as<T, float>) {
        auto res = fast_float::from_chars(sv.data(), sv.data() + sv.size(), x);
        ARTS_USER_ERROR_IF(res.ec != std::errc{}, "Failed to parse {} from string \"{}\"", error_context, orig)
        ARTS_USER_ERROR_IF(
            res.ptr != sv.data() + sv.size(), "Failed to fully parse {} string \"{}\"", error_context, orig)
      } else {
        auto res = std::from_chars(sv.data(), sv.data() + sv.size(), x);
        ARTS_USER_ERROR_IF(res.ec != std::errc{}, "Failed to parse {} from string \"{}\"", error_context, orig)
        ARTS_USER_ERROR_IF(
            res.ptr != sv.data() + sv.size(), "Failed to fully parse {} string \"{}\"", error_context, orig)
      }

      return x;
    }
  }

  constexpr void skip(Size n) {
    ARTS_USER_ERROR_IF(it + n > end, "Unexpected end of string");
    it += n;
  }

  [[nodiscard]] constexpr bool             end_of_string() const { return it == end; }
  [[nodiscard]] constexpr std::string_view remaining_string() const { return {it, end}; }
};

bool read_jpl_entry(jpl_record& record, reader& data) {
  using namespace Conversion;

  record.f0 = data.read_next<Numeric>(13, "Central Frequency"sv) * 1e6;   // MHz to Hz
  record.df = data.read_next<Numeric>(8, "Frequency Deviation"sv) * 1e6;  // MHz to Hz (?)
  record.s =
      std::pow(10., data.read_next<Numeric>(8, "Line Intensity"sv)) / 1e12;  // log_10(nm2 MHz at T0) to ARTS units
  record.dr     = data.read_next<Index>(2, "dr"sv);
  record.E      = kaycm2joule(data.read_next<Numeric>(10, "Energy"sv));  // cm^-1 to J
  record.g_upp  = data.read_next<Index>(3, "Upper State Degeneracy"sv);
  record.jpl_id = Jpl::data_lookup(std::abs(data.read_next<Index>(7, "JPL ID"sv)));
  record.qnfmt  = data.read_next<Index>(4, "Quantum Number Format"sv);

  return true;
}

bool read_jpl_line(jpl_record& record, const std::string& linedata) try {
  reader data(linedata);

  return read_jpl_entry(record, data);
} catch (std::exception& e) {
  ARTS_USER_ERROR("Internal error:\n\n{}\n\nFailed to read JPL line record:\n\n{}", e.what(), linedata);
}
}  // namespace

jpl_data read_jpl_lines(std::istream& file) {
  jpl_data out;

  std::string linedata;
  bool        last_ok = true;

  while (std::getline(file, linedata)) { last_ok = read_jpl_line(last_ok ? out.emplace_back() : out.back(), linedata); }

  if (not last_ok) out.pop_back();

  return out;
}

jpl_data read_jpl_lines(std::istream&& file) { return read_jpl_lines(file); }

line jpl_record::from() const {
  using enum LineShapeModelVariable;
  using enum LineShapeModelType;

  line out;

  out.f0 = f0;
  out.e0 = E;
  out.gu = g_upp > 0 ? static_cast<Numeric>(g_upp) : -1;
  out.gl = -1.0;
  // Rescale so at T0, the line intensity matches the JPL value.
  // This is needed because ARTS defines the line intensities based
  // on built-in partition functions, which may differ from the JPL ones.
  out.a = einstein_a(
      s, out.gu, out.e0, out.f0, jpl_id.T0, Math::pow2(jpl_id.QT0) / PartitionFunctions::Q(jpl_id.T0, jpl_id.qid.isot));
  out.z.on                                  = false;
  out.z.mdata                               = {};
  out.qn                                    = {};
  out.ls                                    = line_shape::model{};
  out.ls.T0                                 = jpl_id.T0;
  out.ls.single_models["AIR"_spec].data[G0] = {T1, Vector{25e3, 0.75}};

  return out;
}
}  // namespace lbl