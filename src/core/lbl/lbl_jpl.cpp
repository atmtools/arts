#include "lbl_jpl.h"

#include <arts_conversions.h>
#include <fast_float/fast_float.h>
#include <jpl_species.h>
#include <mystring.h>
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

  template <typename T> constexpr T read_next(Size n, std::string_view error_context) {
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
        if (res.ec != std::errc{}) arts::user_error("Failed to parse {} from string \"{}\"", error_context, orig);
        if (res.ptr != sv.data() + sv.size())
          arts::user_error("Failed to fully parse {} string \"{}\"", error_context, orig);
      } else {
        auto res = std::from_chars(sv.data(), sv.data() + sv.size(), x);
        if (res.ec != std::errc{}) arts::user_error("Failed to parse {} from string \"{}\"", error_context, orig);
        if (res.ptr != sv.data() + sv.size())
          arts::user_error("Failed to fully parse {} string \"{}\"", error_context, orig);
      }

      return x;
    }
  }

  constexpr void skip(Size n) {
    if (it + n > end) arts::user_error("Unexpected end of string");
    it += n;
  }

  /*! The upper state degeneracy, a 3-character integer field

    Values above 999 are written with a letter for the hundreds, A = 10, so that
    A32 is 1032.  They are only read if parse_letters is true.  A field with an
    asterisk holds a degeneracy that does not even fit that way, and is read as
    an unknown degeneracy, -1.
  */
  Index read_degeneracy(bool parse_letters) {
    if (it + 3 > end) arts::user_error("Unexpected end of string");
    const std::string_view orig(it, it + 3);

    if (orig.find('*') != std::string_view::npos) {
      skip(3);
      return -1;
    }

    if (orig.front() >= 'A' and orig.front() <= 'Z') {
      const Index hundreds = orig.front() - 'A' + 10;
      skip(1);
      const Index rest = read_next<Index>(2, "Upper State Degeneracy"sv);
      if (not parse_letters)
        arts::user_error(
            "The upper state degeneracy \"{}\" has a letter for the hundreds.  Set parse_letters to read it as {}",
            orig,
            100 * hundreds + rest);
      return 100 * hundreds + rest;
    }

    return read_next<Index>(3, "Upper State Degeneracy"sv);
  }

  [[nodiscard]] constexpr bool             end_of_string() const { return it == end; }
  [[nodiscard]] constexpr std::string_view remaining_string() const { return {it, end}; }
};

bool read_jpl_entry(jpl_record& record, reader& data, std::span<const JplSpeciesInfo> species, bool parse_letters) {
  using namespace Conversion;

  record.f0 = data.read_next<Numeric>(13, "Central Frequency"sv) * 1e6;   // MHz to Hz
  record.df = data.read_next<Numeric>(8, "Frequency Deviation"sv) * 1e6;  // MHz to Hz (?)
  record.s =
      std::pow(10., data.read_next<Numeric>(8, "Line Intensity"sv)) / 1e12;  // log_10(nm2 MHz at T0) to ARTS units
  record.dr     = data.read_next<Index>(2, "dr"sv);
  record.E      = kaycm2joule(data.read_next<Numeric>(10, "Energy"sv));  // cm^-1 to J
  record.g_upp  = data.read_degeneracy(parse_letters);
  record.jpl_id = Jpl::data_lookup(std::abs(data.read_next<Index>(7, "JPL ID"sv)), species);
  record.qnfmt  = data.read_next<Index>(4, "Quantum Number Format"sv);

  return true;
}

bool read_jpl_line(jpl_record&                     record,
                   const std::string&              linedata,
                   std::span<const JplSpeciesInfo> species,
                   bool                            parse_letters) try {
  reader data(linedata);

  return read_jpl_entry(record, data, species, parse_letters);
} catch (std::exception& e) {
  arts::user_error("Failed to read JPL line record:\n\n{}\n\n{}", linedata, indent(e.what()));
}
}  // namespace

jpl_data read_jpl_lines(std::istream& file, std::span<const JplSpeciesInfo> species, bool parse_letters) {
  jpl_data out;

  std::string linedata;
  bool        last_ok = true;

  while (std::getline(file, linedata)) {
    last_ok = read_jpl_line(last_ok ? out.emplace_back() : out.back(), linedata, species, parse_letters);
  }

  if (not last_ok) out.pop_back();

  return out;
}

jpl_data read_jpl_lines(std::istream&& file, std::span<const JplSpeciesInfo> species, bool parse_letters) {
  return read_jpl_lines(file, species, parse_letters);
}

line jpl_record::from() const {
  using enum LineShapeModelVariable;
  using enum LineShapeModelType;

  line out;

  out.f0 = f0;
  out.e0 = E;
  out.gu = g_upp > 0 ? static_cast<Numeric>(g_upp) : -1;
  out.gl = -1.0;

  // The Einstein A coefficient follows from the JPL intensity with JPL's own
  // upper state degeneracy and partition function at T0.
  out.a = einstein_a(s, out.gu, out.e0, out.f0, jpl_id.T0, jpl_id.QT0);

  // JPL's degeneracies and partition functions often leave out factors that are
  // common to all levels, such as nuclear spin statistics, which the partition
  // functions of ARTS include.  Scaling the upper state degeneracy by the ratio
  // of the partition functions at T0 converts it to the convention of ARTS, so
  // that the line strength at T0 equals the JPL intensity.
  out.gu *= PartitionFunctions::Q(jpl_id.T0, jpl_id.qid.isot) / jpl_id.QT0;

  out.z.on                                  = false;
  out.z.mdata                               = {};
  out.qn                                    = {};
  out.ls                                    = line_shape::model{};
  out.ls.T0                                 = jpl_id.T0;
  out.ls.single_models["AIR"_spec].data[G0] = {T1, Vector{25e3, 0.75}};

  return out;
}
}  // namespace lbl