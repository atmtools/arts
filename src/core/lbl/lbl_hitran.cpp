#include "lbl_hitran.h"

#include <arts_conversions.h>
#include <fast_float/fast_float.h>
#include <hitran_species.h>
#include <partfun.h>
#include <quantum.h>

#include <charconv>

namespace lbl {
namespace {
struct reader {
  std::string::const_iterator it;
  std::string::const_iterator end;

  reader(const std::string& s) : it(s.begin()), end(s.end()) {}

  template <typename T>
  constexpr T read_next(Size n) {
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
        ARTS_USER_ERROR_IF(res.ec != std::errc{},
                           "Failed to parse value from string \"{}\"",
                           orig)
        ARTS_USER_ERROR_IF(res.ptr != sv.data() + sv.size(),
                           "Failed to fully parse string \"{}\"",
                           orig)
      } else {
        auto res = std::from_chars(sv.data(), sv.data() + sv.size(), x);
        ARTS_USER_ERROR_IF(res.ec != std::errc{},
                           "Failed to parse value from string \"{}\"",
                           orig)
        ARTS_USER_ERROR_IF(res.ptr != sv.data() + sv.size(),
                           "Failed to fully parse string \"{}\"",
                           orig)
      }

      return x;
    }
  }

  constexpr void skip(Size n) {
    ARTS_USER_ERROR_IF(it + n > end, "Unexpected end of string");
    it += n;
  }

  [[nodiscard]] constexpr bool end_of_string() const { return it == end; }
};

bool read_hitran_par_record(hitran_record& record,
                            const std::string& linedata,
                            const Numeric fmin) try {
  using namespace Conversion;

  reader data(linedata);

  const auto M = data.read_next<Index>(2);
  const auto I = data.read_next<char>(1);

  record.f0 = kaycm2freq(data.read_next<Numeric>(12));
  if (record.f0 < fmin) return false;

  // Set this after the frequency check to avoid unnecessary work
  record.qid = Hitran::id_from_lookup(M, I);

  record.S = kaycm_per_cmsquared2hz_per_msquared(data.read_next<Numeric>(10));
  record.A = data.read_next<Numeric>(10);
  record.gamma_air  = kaycm_per_atm2hz_per_pa(data.read_next<Numeric>(5));
  record.gamma_self = kaycm_per_atm2hz_per_pa(data.read_next<Numeric>(5));
  record.E          = kaycm2joule(data.read_next<Numeric>(10));
  record.n          = data.read_next<Numeric>(4);
  record.delta      = kaycm_per_atm2hz_per_pa(data.read_next<Numeric>(8));
  data.skip(79);
  record.g_upp = data.read_next<Numeric>(7);
  record.g_low = data.read_next<Numeric>(7);

  if (not data.end_of_string()) {
    const std::string_view remainder{data.it + 1, data.end};
    const std::string_view::const_iterator space =
        stdr::find_if(remainder, nonstd::isspace);
    ARTS_USER_ERROR_IF(space == remainder.end(),
                       "Failed to parse HITRAN Quantum numbers:\n\n{}",
                       remainder);

    record.qid.state =
        Quantum::from_hitran(std::string_view{remainder.begin(), space},
                             std::string_view{space + 1, remainder.end()});
  }

  return true;
} catch (std::exception& e) {
  ARTS_USER_ERROR(
      "Internal error:\n\n{}\n\nFailed to read HITRAN line record:\n\n{}",
      e.what(),
      linedata);
}
}  // namespace

hitran_data read_hitran_par(std::istream& file,
                            const Vector2& frequency_range) {
  hitran_data out;

  std::string linedata;
  bool last_ok = true;

  while (std::getline(file, linedata)) {
    last_ok = read_hitran_par_record(last_ok ? out.emplace_back() : out.back(),
                                     linedata,
                                     frequency_range[0]);

    if (last_ok and out.back().f0 > frequency_range[1]) {
      out.pop_back();
      break;
    }
  }

  if (not last_ok) out.pop_back();

  return out;
}

hitran_data read_hitran_par(std::istream&& file,
                            const Vector2& frequency_range) {
  return read_hitran_par(file, frequency_range);
}

line hitran_record::from(HitranLineStrengthOption ls,
                         QuantumState&& local,
                         bool do_zeeman) const {
  using enum LineShapeModelVariable;
  using enum LineShapeModelType;

  line l;
  l.a  = A;
  l.f0 = f0;
  l.e0 = E;
  l.gu = g_upp;
  l.gl = g_low;

  switch (ls) {
    case HitranLineStrengthOption::S:
      if (g_upp == 0.0) {
        l.gu = -1.0;
        l.gl = -1.0;
      }

      l.a = l.hitran_a(S, qid.isot);

      break;
    case HitranLineStrengthOption::A: break;
  }

  ARTS_USER_ERROR_IF(
      not std::isnormal(l.a) or not std::isnormal(l.gu),
      "Invalid Einstein coefficient {} or gu {} for full HITRAN RECORD: {}",
      l.a,
      l.gu,
      *this)

  if (do_zeeman) {
    l.z = lbl::zeeman::GetAdvancedModel(qid);
  } else {
    l.z = {};
  }
  l.z.on = false;

  // Set the line shape
  l.ls    = line_shape::model{};
  l.ls.T0 = 296.0;
  l.ls.single_models.reserve(2);

  auto& self    = l.ls.single_models[qid.isot.spec];
  auto& air     = l.ls.single_models[SpeciesEnum::Bath];
  self.data[G0] = {T1, Vector{gamma_self, n}};
  air.data[G0]  = {T1, Vector{gamma_air, n}};

  if (delta != 0) {
    self.data[D0] = {T0, Vector{delta}};
    air.data[D0]  = {T0, Vector{delta}};
  }

  l.qn = std::move(local);

  return l;
}
}  // namespace lbl