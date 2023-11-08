#pragma once

#include <arts_constants.h>
#include <atm.h>
#include <empty.h>
#include <isotopologues.h>
#include <matpack.h>
#include <physics_funcs.h>
#include <quantum_numbers.h>
#include <rtepack.h>
#include <sorting.h>

#include <Faddeeva/Faddeeva.hh>
#include <algorithm>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>

#include "configtypes.h"
#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "lbl_zeeman.h"
#include "matpack_view.h"

namespace lbl::voigt::lte {
struct single_shape {
  //! Line center after all adjustments, must be as_zeeman if zeeman effect is intended
  Numeric f0{};

  //! Inverse of the Doppler broadening factor (missing an ln2)
  Numeric inv_gd{};

  //! The imaginary part of the complex argument that goes into the Faddeeva function
  Numeric z_imag{};

  //! Linestrength but lacks N * f * (1 - exp(-hf/kt)) factor, must be as_zeeman if zeeman effect is intended, also has Constant::inv_sqrt_pi * inv_gd factor
  Complex s{};

  constexpr single_shape() = default;

  single_shape(const SpeciesIsotopeRecord& spec, const line&, const AtmPoint&);

  single_shape(const SpeciesIsotopeRecord& spec,
               const line&,
               const AtmPoint&,
               const Size);

  void as_zeeman(const line& line,
                 const Numeric H,
                 const zeeman::pol type,
                 const Index iz);

  [[nodiscard]] Complex z(Numeric f) const {
    return Complex{inv_gd * (f - f0), z_imag};
  }

  [[nodiscard]] Complex F(const Complex z_) const noexcept {
    return Faddeeva::w(z_);  // FIXME: Should factor be part of s?
  }

  [[nodiscard]] Complex F(const Numeric f) const noexcept { return F(z(f)); }

  [[nodiscard]] Complex operator()(const Numeric f) const noexcept {
    return s * F(f);
  }

  [[nodiscard]] Complex dF(const Complex z_, const Complex F_) const noexcept {
    return 2 * inv_gd * (1i * Constant::inv_pi * inv_gd - z_ * F_);
  }

  [[nodiscard]] Complex dF(const Numeric f) const noexcept {
    Complex z_ = z(f);
    return dF(z_, F(z_));
  }

 private:
  struct zFdF {
    Complex z, F, dF;
  };

  [[nodiscard]] zFdF all(const Numeric f) const noexcept {
    zFdF out;
    out.z = z(f);
    out.F = F(out.z);
    out.dF = dF(out.z, out.F);
    return out;
  }

 public:
  [[nodiscard]] Complex df(const Numeric f) const noexcept { return s * dF(f); }

  [[nodiscard]] Complex df0(const Complex ds_df0,
                            const Numeric f) const noexcept {
    const auto [z_, F_, dF_] = all(f);
    return ds_df0 * F_ - s * dF_;
  }

  [[nodiscard]] Complex dDV(const Complex dz_dDV,
                            const Numeric f) const noexcept {
    return s * dz_dDV * dF(f);
  }

  [[nodiscard]] Complex dD0(const Complex dz_dD0,
                            const Numeric f) const noexcept {
    return s * dz_dD0 * dF(f);
  }

  [[nodiscard]] Complex dG0(const Complex dz_dG0,
                            const Numeric f) const noexcept {
    return s * dz_dG0 * dF(f);
  }

  [[nodiscard]] Complex dH(const Complex dz_dH,
                           const Numeric f) const noexcept {
    return s * dz_dH * dF(f);
  }

  [[nodiscard]] Complex dVMR(const Complex ds_dVMR,
                             const Complex dz_dVMR,
                             const Numeric f) const noexcept {
    const auto [z_, F_, dF_] = all(f);
    return ds_dVMR * F_ + dz_dVMR * dF_;
  }

  [[nodiscard]] Complex dT(const Complex ds_dT,
                           const Complex dz_dT,
                           const Numeric f) const noexcept {
    const auto [z_, F_, dF_] = all(f);
    return ds_dT * F_ + dz_dT * dF_;
    //FIXME: invGD factor of F_ is missing
  }

  [[nodiscard]] Complex da(const Complex ds_da,
                           const Numeric f) const noexcept {
    return ds_da * F(f);
  }

  [[nodiscard]] Complex de0(const Complex ds_de0,
                            const Numeric f) const noexcept {
    return ds_de0 * F(f);
  }

  [[nodiscard]] Complex dG(const Complex ds_dG,
                           const Numeric f) const noexcept {
    return ds_dG * F(f);
  }

  [[nodiscard]] Complex dY(const Complex ds_dY,
                           const Numeric f) const noexcept {
    return ds_dY * F(f);
  }
};

struct line_pos {
  Size line;
  Size spec{std::numeric_limits<Size>::max()};
  Size iz{std::numeric_limits<Size>::max()};
};

inline Size count_lines(const band& bnd, const zeeman::pol type) {
  return std::transform_reduce(
      bnd.begin(), bnd.end(), Index{}, std::plus<>{}, [type](auto& line) {
        const Index factor =
            line.ls.one_by_one ? line.ls.single_models.size() : 1;
        return factor * line.z.size(line.qn.val, type);
      });
}

inline void zeeman_set_back(std::vector<single_shape>& lines,
                     std::vector<line_pos>& pos,
                     const single_shape& s,
                     const line& line,
                     const Numeric H,
                     const Size spec,
                     const zeeman::pol pol,
                     Size& last_single_shape_pos) {
  const auto line_nr = static_cast<Size>(pos.size() ? pos.back().line + 1 : 0);

  if (pol == zeeman::pol::no) {
    lines[last_single_shape_pos] = s;
    pos[last_single_shape_pos] = {line_nr, spec};
    ++last_single_shape_pos;
  } else {
    const auto nz = line.z.size(line.qn.val, pol);
    for (Index iz = 0; iz < nz; iz++) {
      lines[last_single_shape_pos] = s;
      lines[last_single_shape_pos].as_zeeman(line, H, pol, iz);
      pos[last_single_shape_pos] = {line_nr, spec, static_cast<Size>(iz)};
      ++last_single_shape_pos;
    }
  }
}

inline void lines_set(std::vector<single_shape>& lines,
               std::vector<line_pos>& pos,
               const SpeciesIsotopeRecord& spec,
               const line& line,
               const AtmPoint& atm,
               const zeeman::pol pol,
               Size& last_single_shape_pos) {
  const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  if (line.ls.one_by_one) {
    for (Size i = 0; i < line.ls.single_models.size(); ++i) {
      zeeman_set_back(lines,
                      pos,
                      {spec, line, atm, i},
                      line,
                      H,
                      i,
                      pol,
                      last_single_shape_pos);
    }
  } else {
    zeeman_set_back(lines,
                    pos,
                    {spec, line, atm},
                    line,
                    H,
                    std::numeric_limits<Size>::max(),
                    pol,
                    last_single_shape_pos);
  }
}

//! Helper for initializing the band_shape
inline void band_shape_helper(std::vector<single_shape>& lines,
                              std::vector<line_pos>& pos,
                              const SpeciesIsotopeRecord& spec,
                              const band& bnd,
                              const AtmPoint& atm,
                              const Numeric fmin,
                              const Numeric fmax,
                              const zeeman::pol pol) {
  lines.resize(count_lines(bnd, pol));
  pos.resize(lines.size());

  Size i{0};

  using enum CutoffType;
  switch (bnd.cutoff) {
    case None:
      for (auto& line : bnd) {
        lines_set(lines, pos, spec, line, atm, pol, i);
      }
      break;
    case Freq: {
      const auto by_range =
          std::views::filter([fmin, fmax, c = bnd.cutoff_value](auto& line) {
            return (line.f0 - c) > fmin and (line.f0 + c) < fmax;
          });
      for (auto& line : bnd | by_range) {
        lines_set(lines, pos, spec, line, atm, pol, i);
      }
    } break;
    case FINAL:
      ARTS_USER_ERROR("Bad state")
  }

  bubble_sort_by(
      [&](const Size l1, const Size l2) { return lines[l1].f0 < lines[l2].f0; },
      lines,
      pos);
}

inline std::pair<Index, Index> find_offset_and_count_of_frequency_range(
    const std::span<const single_shape> lines, Numeric f, Numeric cutoff) {
  if (cutoff < std::numeric_limits<Numeric>::infinity()) {
    auto first = lines.begin();

    struct limit {
      Numeric cutoff;
      constexpr bool operator()(const single_shape& line,
                                Numeric f) const noexcept {
        return cutoff < std::abs(f - line.f0);
      }
      constexpr bool operator()(Numeric f,
                                const single_shape& line) const noexcept {
        return cutoff < std::abs(f - line.f0);
      }
    };

    const auto [start_range, end_range] =
        std::equal_range(first, lines.end(), f, limit{cutoff});

    return {std::distance(first, start_range),
            std::distance(start_range, end_range)};
  }

  return {0, lines.size()};
}

namespace detail {
auto frequency_span(const auto& list, const Size start, const Size count) {
  std::span out{list};
  return out.subspan(start, count);
}
}  // namespace detail

template <typename... Ts>
auto frequency_spans(const Numeric cutoff,
                     const Numeric f,
                     const std::span<const single_shape>& lines,
                     const Ts&... lists) {
  ARTS_ASSERT(lines.size() == (static_cast<Size>(lists.size()) and ...))

  const auto [start, count] =
      find_offset_and_count_of_frequency_range(lines, f, cutoff);
  return std::tuple{detail::frequency_span(lines, start, count),
                    detail::frequency_span(lists, start, count)...};
}

inline std::tuple<std::span<const single_shape>, std::span<const line_pos>>
select_lines(const std::vector<single_shape>& lines,
             const std::vector<line_pos>& pos,
             const line_key& key) {
  ARTS_ASSERT(pos.size() == lines.size())

  const auto pred = [line = key.line](auto& a) { return a.line == line; };

  auto start_pos = std::find_if(pos.begin(), pos.end(), pred);

  if (start_pos == pos.end()) {
    return {{}, {}};
  }

  auto last_pos = std::find_if(pos.rbegin(), pos.rend(), pred);

  auto offset = std::distance(pos.begin(), start_pos);
  auto last_index = std::distance(last_pos, pos.rend());
  auto count = last_index - offset + 1;
  return {std::span<const single_shape>{lines}.subspan(offset, count),
          std::span<const line_pos>{pos}.subspan(offset, count)};
}

inline std::tuple<std::span<const single_shape>,
                  std::span<const line_pos>,
                  std::span<const Complex>>
select_lines(const std::vector<single_shape>& lines,
             const std::vector<line_pos>& pos,
             const line_key& key,
             const ExhaustiveConstComplexVectorView& cut) {
  ARTS_ASSERT(pos.size() == lines.size())

  const auto pred = [line = key.line](auto& a) { return a.line == line; };

  auto start_pos = std::find_if(pos.begin(), pos.end(), pred);

  if (start_pos == pos.end()) {
    return {{}, {}, {}};
  }

  auto last_pos = std::find_if(pos.rbegin(), pos.rend(), pred);

  auto offset = std::distance(pos.begin(), start_pos);
  auto last_index = std::distance(last_pos, pos.rend());
  auto count = last_index - offset + 1;
  return {std::span<const single_shape>{lines}.subspan(offset, count),
          std::span<const line_pos>{pos}.subspan(offset, count),
          std::span<const Complex>{cut}.subspan(offset, count)};
}

inline std::tuple<std::span<const single_shape>,
                  std::span<const line_pos>,
                  std::span<Complex>>
select_lines(const std::vector<single_shape>& lines,
             const std::vector<line_pos>& pos,
             const line_key& key,
             ExhaustiveComplexVectorView cut) {
  ARTS_ASSERT(pos.size() == lines.size())

  const auto pred = [line = key.line](auto& a) { return a.line == line; };

  auto start_pos = std::find_if(pos.begin(), pos.end(), pred);

  if (start_pos == pos.end()) {
    return {{}, {}, {}};
  }

  auto last_pos = std::find_if(pos.rbegin(), pos.rend(), pred);

  auto offset = std::distance(pos.begin(), start_pos);
  auto last_index = std::distance(last_pos, pos.rend());
  auto count = last_index - offset + 1;
  return {std::span<const single_shape>{lines}.subspan(offset, count),
          std::span<const line_pos>{pos}.subspan(offset, count),
          std::span<Complex>{cut}.subspan(offset, count)};
}

//! A band shape is a collection of single shapes.  The shapes are sorted by frequency.
struct band_shape {
  //! Line absorption shapes (lacking the f * (1 - exp(-hf/kt)) factor)
  std::vector<single_shape> lines{};
  Numeric cutoff{-1};

  [[nodiscard]] Size size() const noexcept { return lines.size(); }

  band_shape() = default;

  band_shape(std::vector<single_shape>&& ls, const Numeric cut) noexcept
      : lines(std::move(ls)), cutoff(cut) {}

  [[nodiscard]] Complex operator()(const Numeric f) const noexcept {
    return std::transform_reduce(
        lines.begin(), lines.end(), Complex{}, std::plus<>{}, [f](auto& ls) {
          return ls(f);
        });
  }

  [[nodiscard]] Complex df(const Numeric f) const {
    return std::transform_reduce(
        lines.begin(), lines.end(), Complex{}, std::plus<>{}, [f](auto& ls) {
          return ls.df(f);
        });
  }

  [[nodiscard]] Complex dH(const ExhaustiveConstComplexVectorView& dz_dH,
                           const Numeric f) const {
    ARTS_ASSERT(static_cast<Size>(dz_dH.size()) == lines.size())

    return std::transform_reduce(
        lines.begin(),
        lines.end(),
        dz_dH.begin(),
        Complex{},
        std::plus<>{},
        [f](auto& ls, auto& d) { return ls.dH(d, f); });
  }

  [[nodiscard]] Complex dT(const ExhaustiveConstComplexVectorView& ds_dT,
                           const ExhaustiveConstComplexVectorView& dz_dT,
                           const Numeric f) const {
    ARTS_ASSERT(ds_dT.size() == dz_dT.size())
    ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < lines.size(); ++i) {
      out += lines[i].dT(ds_dT[i], dz_dT[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex dVMR(const ExhaustiveConstComplexVectorView& ds_dVMR,
                             const ExhaustiveConstComplexVectorView& dz_dVMR,
                             const Numeric f) const {
    ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
    ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < lines.size(); ++i) {
      out += lines[i].dVMR(ds_dVMR[i], dz_dVMR[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex df0(const Complex ds_df0,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::f0) {
      return {};
    }

    const auto [s, ps] = select_lines(lines, pos, line_key);

    return std::transform_reduce(
        s.begin(),
        s.end(),
        ps.begin(),
        Complex{},
        std::plus<>{},
        [f, ds_df0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line) ? ls.df0(ds_df0, f) : Complex{};
        });
  }

  [[nodiscard]] Complex da(const Complex ds_da,
                           const Numeric f,
                           const std::vector<line_pos>& pos,
                           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::a) {
      return {};
    }

    const auto [s, ps] = select_lines(lines, pos, line_key);

    return std::transform_reduce(
        s.begin(),
        s.end(),
        ps.begin(),
        Complex{},
        std::plus<>{},
        [f, ds_da, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line) ? ls.da(ds_da, f) : Complex{};
        });
  }

  [[nodiscard]] Complex de0(const Complex ds_de0,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::a) {
      return {};
    }

    const auto [s, ps] = select_lines(lines, pos, line_key);

    return std::transform_reduce(
        s.begin(),
        s.end(),
        ps.begin(),
        Complex{},
        std::plus<>{},
        [f, ds_de0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line) ? ls.de0(ds_de0, f) : Complex{};
        });
  }

  [[nodiscard]] Complex dDV(const Complex dz_dDV,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::DV) {
      return {};
    }

    const auto [s, ps] = select_lines(lines, pos, line_key);

    return std::transform_reduce(
        s.begin(),
        s.end(),
        ps.begin(),
        Complex{},
        std::plus<>{},
        [f, dz_dDV, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dDV(dz_dDV, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dD0(const Complex dz_dD0,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::D0) {
      return {};
    }

    const auto [s, ps] = select_lines(lines, pos, line_key);

    return std::transform_reduce(
        s.begin(),
        s.end(),
        ps.begin(),
        Complex{},
        std::plus<>{},
        [f, dz_dD0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dD0(dz_dD0, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dG0(const Complex dz_dG0,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::G0) {
      return {};
    }

    const auto [s, ps] = select_lines(lines, pos, line_key);

    return std::transform_reduce(
        s.begin(),
        s.end(),
        ps.begin(),
        Complex{},
        std::plus<>{},
        [f, dz_dG0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dG0(dz_dG0, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dY(const Complex ds_dY,
                           const Numeric f,
                           const std::vector<line_pos>& pos,
                           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::Y) {
      return {};
    }

    const auto [s, ps] = select_lines(lines, pos, line_key);

    return std::transform_reduce(
        s.begin(),
        s.end(),
        ps.begin(),
        Complex{},
        std::plus<>{},
        [f, ds_dY, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dY(ds_dY, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dG(const Complex ds_dG,
                           const Numeric f,
                           const std::vector<line_pos>& pos,
                           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::G) {
      return {};
    }

    const auto [s, ps] = select_lines(lines, pos, line_key);

    return std::transform_reduce(
        s.begin(),
        s.end(),
        ps.begin(),
        Complex{},
        std::plus<>{},
        [f, ds_dG, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dG(ds_dG, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex operator()(const ExhaustiveConstComplexVectorView& cut,
                                   const Numeric f) const {
    const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
    return std::transform_reduce(s.begin(),
                                 s.end(),
                                 cs.begin(),
                                 Complex{},
                                 std::plus<>{},
                                 [f](auto& ls, auto& c) { return ls(f) - c; });
  }

  void operator()(ExhaustiveComplexVectorView cut) const {
    std::transform(
        lines.begin(),
        lines.end(),
        cut.begin(),
        [cutoff_freq = cutoff](auto& ls) { return ls(ls.f0 + cutoff_freq); });
  }

  [[nodiscard]] Complex df(const ExhaustiveConstComplexVectorView& cut,
                           const Numeric f) const {
    const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
    return std::transform_reduce(
        s.begin(),
        s.end(),
        cs.begin(),
        Complex{},
        std::plus<>{},
        [f](auto& ls, auto& c) { return ls.df(f) - c; });
  }

  void df(ExhaustiveComplexVectorView cut) const {
    std::transform(lines.begin(),
                   lines.end(),
                   cut.begin(),
                   [cutoff_freq = cutoff](auto& ls) {
                     return ls.df(ls.f0 + cutoff_freq);
                   });
  }

  [[nodiscard]] Complex dH(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView& dz_dH,
                           const Numeric f) const {
    ARTS_ASSERT(static_cast<Size>(dz_dH.size()) == lines.size())

    const auto [s, cs, dH] = frequency_spans(cutoff, f, lines, cut, dz_dH);

    Complex out{};  //! Fixme, use zip in C++ 23...
    for (Size i = 0; i < s.size(); ++i) {
      out += s[i].dH(dH[i], f) - cs[i];
    }

    return out;
  }

  void dH(ExhaustiveComplexVectorView cut,
          const ExhaustiveConstComplexVectorView& df0_dH) const {
    ARTS_ASSERT(static_cast<Size>(df0_dH.size()) == lines.size())

    std::transform(lines.begin(),
                   lines.end(),
                   df0_dH.begin(),
                   cut.begin(),
                   [cutoff_freq = cutoff](auto& ls, auto& d) {
                     return ls.dH(d, ls.f0 + cutoff_freq);
                   });
  }

  [[nodiscard]] Complex dT(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView& ds_dT,
                           const ExhaustiveConstComplexVectorView& dz_dT,
                           const Numeric f) const {
    ARTS_ASSERT(ds_dT.size() == dz_dT.size())
    ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

    Complex out{};  //! Fixme, use zip in C++ 23...

    const auto [s, cs, ds, dz] =
        frequency_spans(cutoff, f, lines, cut, ds_dT, dz_dT);

    for (Size i = 0; i < s.size(); ++i) {
      out += s[i].dT(ds[i], dz[i], f) - cs[i];
    }

    return out;
  }

  void dT(ExhaustiveComplexVectorView cut,
          const ExhaustiveConstComplexVectorView& ds_dT,
          const ExhaustiveConstComplexVectorView& dz_dT) const {
    ARTS_ASSERT(ds_dT.size() == dz_dT.size())
    ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

    for (Size i = 0; i < lines.size(); ++i) {
      cut[i] += lines[i].dT(ds_dT[i], dz_dT[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex dVMR(const ExhaustiveConstComplexVectorView& cut,
                             const ExhaustiveConstComplexVectorView& ds_dVMR,
                             const ExhaustiveConstComplexVectorView& dz_dVMR,
                             const Numeric f) const {
    ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
    ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

    Complex out{};  //! Fixme, use zip in C++ 23...

    const auto [s, cs, ds, dz] =
        frequency_spans(cutoff, f, lines, cut, ds_dVMR, dz_dVMR);

    for (Size i = 0; i < s.size(); ++i) {
      out += s[i].dVMR(ds[i], dz[i], f) - cs[i];
    }

    return out;
  }

  void dVMR(ExhaustiveComplexVectorView cut,
            const ExhaustiveConstComplexVectorView& ds_dVMR,
            const ExhaustiveConstComplexVectorView& dz_dVMR) const {
    ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
    ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

    for (Size i = 0; i < lines.size(); ++i) {
      cut[i] += lines[i].dVMR(ds_dVMR[i], dz_dVMR[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex df0(const ExhaustiveConstComplexVectorView& cut,
                            const Complex ds_df0,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::f0) {
      return {};
    }

    const auto [sp, psp, csp] = select_lines(lines, pos, line_key, cut);

    const auto [s, ps, cs] = frequency_spans(cutoff, f, sp, psp, csp);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < s.size(); i++) {
      if (ps[i].line == line_key.line) out += s[i].df0(ds_df0, f) - cs[i];
    }

    return out;
  }

  void df0(ExhaustiveComplexVectorView cut,
           const Complex ds_df0,
           const std::vector<line_pos>& pos,
           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::f0) {
      return;
    }

    const auto [s, ps, cs] = select_lines(lines, pos, line_key, cut);

    std::transform(
        s.begin(),
        s.end(),
        ps.begin(),
        cs.begin(),
        [ds_df0, &line_key, cutoff_freq = cutoff](auto& ls, auto& p) {
          return (p.line == line_key.line) ? ls.df0(ds_df0, ls.f0 + cutoff_freq)
                                           : Complex{};
        });
  }

  [[nodiscard]] Complex da(const ExhaustiveConstComplexVectorView& cut,
                           const Complex ds_da,
                           const Numeric f,
                           const std::vector<line_pos>& pos,
                           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::a) {
      return {};
    }

    const auto [sp, psp, csp] = select_lines(lines, pos, line_key, cut);

    const auto [s, ps, cs] = frequency_spans(cutoff, f, sp, psp, csp);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < s.size(); i++) {
      if (ps[i].line == line_key.line) out += s[i].da(ds_da, f) - cs[i];
    }

    return out;
  }

  void da(ExhaustiveComplexVectorView cut,
          const Complex ds_da,
          const std::vector<line_pos>& pos,
          const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::a) {
      return;
    }

    const auto [s, ps, cs] = select_lines(lines, pos, line_key, cut);

    std::transform(s.begin(),
                   s.end(),
                   ps.begin(),
                   cs.begin(),
                   [cutoff_freq = cutoff, ds_da, &line_key](auto& ls, auto& p) {
                     return (p.line == line_key.line)
                                ? ls.da(ds_da, ls.f0 + cutoff_freq)
                                : Complex{};
                   });
  }

  [[nodiscard]] Complex de0(const ExhaustiveConstComplexVectorView& cut,
                            const Complex ds_de0,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::a) {
      return {};
    }

    const auto [sp, psp, csp] = select_lines(lines, pos, line_key, cut);

    const auto [s, ps, cs] = frequency_spans(cutoff, f, sp, psp, csp);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < s.size(); i++) {
      if (ps[i].line == line_key.line) out += s[i].de0(ds_de0, f) - cs[i];
    }

    return out;
  }

  void de0(ExhaustiveComplexVectorView cut,
           const Complex ds_de0,
           const std::vector<line_pos>& pos,
           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.var != variable::a) {
      return;
    }

    const auto [s, ps, cs] = select_lines(lines, pos, line_key, cut);

    std::transform(
        s.begin(),
        s.end(),
        ps.begin(),
        cs.begin(),
        [cutoff_freq = cutoff, ds_de0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line) ? ls.de0(ds_de0, ls.f0 + cutoff_freq)
                                           : Complex{};
        });
  }

  [[nodiscard]] Complex dDV(const ExhaustiveConstComplexVectorView& cut,
                            const Complex dz_dDV,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::DV) {
      return {};
    }

    const auto [sp, psp, csp] = select_lines(lines, pos, line_key, cut);

    const auto [s, ps, cs] = frequency_spans(cutoff, f, sp, psp, csp);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < s.size(); i++) {
      if (ps[i].line == line_key.line and ps[i].spec == line_key.spec)
        out += s[i].dDV(dz_dDV, f) - cs[i];
    }

    return out;
  }

  void dDV(ExhaustiveComplexVectorView cut,
           const Complex dz_dDV,
           const std::vector<line_pos>& pos,
           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::DV) {
      return;
    }

    const auto [s, ps, cs] = select_lines(lines, pos, line_key, cut);

    std::transform(
        s.begin(),
        s.end(),
        ps.begin(),
        cs.begin(),
        [cutoff_freq = cutoff, dz_dDV, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dDV(dz_dDV, ls.f0 + cutoff_freq)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dD0(const ExhaustiveConstComplexVectorView& cut,
                            const Complex dz_dD0,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::D0) {
      return {};
    }

    const auto [sp, psp, csp] = select_lines(lines, pos, line_key, cut);

    const auto [s, ps, cs] = frequency_spans(cutoff, f, sp, psp, csp);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < s.size(); i++) {
      if (ps[i].line == line_key.line and ps[i].spec == line_key.spec)
        out += s[i].dD0(dz_dD0, f) - cs[i];
    }

    return out;
  }

  void dD0(ExhaustiveComplexVectorView cut,
           const Complex dz_dD0,
           const std::vector<line_pos>& pos,
           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::D0) {
      return;
    }

    const auto [s, ps, cs] = select_lines(lines, pos, line_key, cut);

    std::transform(
        s.begin(),
        s.end(),
        ps.begin(),
        cs.begin(),
        [cutoff_freq = cutoff, dz_dD0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dD0(dz_dD0, ls.f0 + cutoff_freq)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dG0(const ExhaustiveConstComplexVectorView& cut,
                            const Complex dz_dG0,
                            const Numeric f,
                            const std::vector<line_pos>& pos,
                            const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::G0) {
      return {};
    }

    const auto [sp, psp, csp] = select_lines(lines, pos, line_key, cut);

    const auto [s, ps, cs] = frequency_spans(cutoff, f, sp, psp, csp);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < s.size(); i++) {
      if (ps[i].line == line_key.line and ps[i].spec == line_key.spec)
        out += s[i].dG0(dz_dG0, f) - cs[i];
    }

    return out;
  }

  void dG0(ExhaustiveComplexVectorView cut,
           const Complex dz_dG0,
           const std::vector<line_pos>& pos,
           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::G0) {
      return;
    }

    const auto [s, ps, cs] = select_lines(lines, pos, line_key, cut);

    std::transform(
        s.begin(),
        s.end(),
        ps.begin(),
        cs.begin(),
        [cutoff_freq = cutoff, dz_dG0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dG0(dz_dG0, ls.f0 + cutoff_freq)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dY(const ExhaustiveConstComplexVectorView& cut,
                           const Complex ds_dY,
                           const Numeric f,
                           const std::vector<line_pos>& pos,
                           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::Y) {
      return {};
    }

    const auto [sp, psp, csp] = select_lines(lines, pos, line_key, cut);

    const auto [s, ps, cs] = frequency_spans(cutoff, f, sp, psp, csp);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < s.size(); i++) {
      if (ps[i].line == line_key.line and ps[i].spec == line_key.spec)
        out += s[i].dY(ds_dY, f) - cs[i];
    }

    return out;
  }

  void dY(ExhaustiveComplexVectorView cut,
          const Complex ds_dY,
          const std::vector<line_pos>& pos,
          const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::Y) {
      return;
    }

    const auto [s, ps, cs] = select_lines(lines, pos, line_key, cut);

    std::transform(
        s.begin(),
        s.end(),
        ps.begin(),
        cs.begin(),
        [cutoff_freq = cutoff, ds_dY, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dY(ds_dY, ls.f0 + cutoff_freq)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dG(const ExhaustiveConstComplexVectorView& cut,
                           const Complex ds_dG,
                           const Numeric f,
                           const std::vector<line_pos>& pos,
                           const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::G) {
      return {};
    }

    const auto [sp, psp, csp] = select_lines(lines, pos, line_key, cut);

    const auto [s, ps, cs] = frequency_spans(cutoff, f, sp, psp, csp);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < s.size(); i++) {
      if (ps[i].line == line_key.line and ps[i].spec == line_key.spec)
        out += s[i].dG(ds_dG, f) - cs[i];
    }

    return out;
  }

  void dG(ExhaustiveComplexVectorView cut,
          const Complex ds_dG,
          const std::vector<line_pos>& pos,
          const line_key& line_key) const {
    ARTS_ASSERT(pos.size() == lines.size())

    if (line_key.ls_var != line_shape::variable::G) {
      return;
    }

    const auto [s, ps, cs] = select_lines(lines, pos, line_key, cut);

    std::transform(
        s.begin(),
        s.end(),
        ps.begin(),
        cs.begin(),
        [cutoff_freq = cutoff, ds_dG, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dG(ds_dG, ls.f0 + cutoff_freq)
                     : Complex{};
        });
  }
};

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}

void calculate(PropmatVectorView pm,
               PropmatMatrixView dpm,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const band_key& bnd_qid,
               const band& bnd,
               const AtmPoint& atm_point,
               const Vector2 los,
               const zeeman::pol pol = zeeman::pol::no);
}  // namespace lbl::voigt::lte
