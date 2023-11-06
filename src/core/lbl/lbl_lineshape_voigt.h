#pragma once

#include <arts_constants.h>
#include <atm.h>
#include <matpack.h>

#include <Faddeeva/Faddeeva.hh>
#include <algorithm>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>

#include "configtypes.h"
#include "empty.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "lbl_zeeman.h"
#include "matpack_data.h"
#include "matpack_view.h"
#include "physics_funcs.h"
#include "quantum_numbers.h"
#include "sorting.h"

namespace lbl::voigt::lte {
struct single_shape {
  //! Linestrenght but lacks f * (1 - exp(-hf/kt)) factor, must be as_zeeman if zeeman effect is intended
  Complex s{};

  //! Line center after all adjustments, must be as_zeeman if zeeman effect is intended
  Numeric f0{};

  //! Inverse of the Doppler broadening factor (missing an ln2)
  Numeric inv_gd{};

  //! The imaginary part of the complex argument that goes into the Faddeeva function
  Numeric z_imag{};

  constexpr single_shape() = default;

  single_shape(const SpeciesIsotopeRecord& spec, const line&, const AtmPoint&);

  single_shape(const SpeciesIsotopeRecord& spec,
               const line&,
               const AtmPoint&,
               const Size);

  void as_zeeman(const line& line, const zeeman::pol type, const Index iz);

  [[nodiscard]] Complex z(Numeric f) const {
    return Complex{inv_gd * (f - f0), z_imag};
  }

  [[nodiscard]] Complex F(const Complex z_) const noexcept {
    return Constant::inv_sqrt_pi * inv_gd *
           Faddeeva::w(z_);  // FIXME: Should factor be part of s?
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
  [[nodiscard]] Complex df(Numeric f) const noexcept { return s * dF(f); }

  [[nodiscard]] Complex df0(const Complex ds_df0,
                            const Numeric f) const noexcept {
    const auto [z_, F_, dF_] = all(f);
    return ds_df0 * F_ - s * dF_;
  }

  [[nodiscard]] Complex dDV(const Numeric dX_dDV,
                            const Numeric f) const noexcept {
    return -s * dX_dDV * dF(f);
  }

  [[nodiscard]] Complex dD0(const Numeric dX_dD0,
                            const Numeric f) const noexcept {
    return -s * dX_dD0 * dF(f);
  }

  [[nodiscard]] Complex dG0(const Numeric dX_dG0,
                            const Numeric f) const noexcept {
    return s * Complex(0, dX_dG0) * dF(f);
  }

  [[nodiscard]] Complex dH(const Numeric df0_dH,
                           const Numeric f) const noexcept {
    return - s * df0_dH * dF(f);
  }

  [[nodiscard]] Complex dVMR(const Complex ds_dVMR,
                             const Numeric dG0_dVMR,
                             const Numeric dD0_dVMR,
                             const Numeric dDV_dVMR,
                             const Numeric f) const noexcept {
    const auto [z_, F_, dF_] = all(f);

    return ds_dVMR * F_ + Complex(-dD0_dVMR - dDV_dVMR, dG0_dVMR) * dF_;
  }

  [[nodiscard]] Complex dT(const Complex ds_dT,
                           const Numeric dG0_dT,
                           const Numeric dD0_dT,
                           const Numeric dDV_dT,
                           const Numeric T,
                           const Numeric f) const noexcept {
    const auto [z_, F_, dF_] = all(f);

    return ds_dT * F_ -
           (F_ * inv_gd + dF_ * z_) * (2 * T * (dD0_dT + dDV_dT) + f0) /
               (2 * T * inv_gd * f0) +
           Complex(-dD0_dT - dDV_dT, dG0_dT) * dF_;
  }

  [[nodiscard]] Complex da(const Complex ds_da,
                           const Numeric f) const noexcept {
    return ds_da * F(f);
  }

  [[nodiscard]] Complex de0(const Complex de0_da,
                            const Numeric f) const noexcept {
    return de0_da * F(f);
  }

  [[nodiscard]] Complex dG(const Numeric dX_dG, const Numeric f) const noexcept {
    return dX_dG * F(f);
  }

  [[nodiscard]] Complex dY(const Numeric dX_dY, const Numeric f) const noexcept {
    return -1i * dX_dY * F(f);
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

void zeeman_push_back(std::vector<single_shape>& lines,
                      auto& pos,
                      const single_shape& s,
                      const line& line,
                      const Size spec,
                      const zeeman::pol pol) {
  const auto line_nr = static_cast<Size>(pos.size() ? pos.back().line + 1 : 0);

  if (pol == zeeman::pol::no) {
    lines.push_back(s);
    pos.emplace_back(line_nr, spec);
  } else {
    const auto nz = line.z.size(line.qn.val, pol);
    for (Index iz = 0; iz < nz; iz++) {
      lines.push_back(s);
      lines.back().as_zeeman(line, pol, iz);
      pos.emplace_back(line_nr, spec, static_cast<Size>(iz));
    }
  }
}

void lines_push_back(std::vector<single_shape>& lines,
                     auto& pos,
                     const SpeciesIsotopeRecord& spec,
                     const line& line,
                     const AtmPoint& atm,
                     const zeeman::pol pol) {
  if (line.ls.one_by_one) {
    for (Size i = 0; i < line.ls.single_models.size(); ++i) {
      zeeman_push_back(lines, pos, {spec, line, atm, i}, line, i, pol);
    }
  } else {
    zeeman_push_back(lines,
                     pos,
                     {spec, line, atm},
                     line,
                     std::numeric_limits<Size>::max(),
                     pol);
  }
}

//! Helper for initializing the band_shape
template <bool keep_position = false>
std::conditional_t<keep_position,
                   std::pair<std::vector<single_shape>, std::vector<line_pos>>,
                   std::pair<std::vector<single_shape>, Empty>>
band_shape_helper(const SpeciesIsotopeRecord& spec,
                  const band& bnd,
                  const AtmPoint& atm,
                  const Numeric fmin,
                  const Numeric fmax,
                  const zeeman::pol pol) {
  std::pair<std::vector<single_shape>, std::vector<line_pos>> out;
  out.first.reserve(count_lines(bnd, pol));
  out.second.reserve(out.first.size());

  using enum CutoffType;
  switch (bnd.cutoff) {
    case None:
      for (auto& line : bnd) {
        lines_push_back(out.first, out.second, spec, line, atm, pol);
      }
      break;
    case Freq: {
      const auto by_range =
          std::views::filter([fmin, fmax, c = bnd.cutoff_value](auto& line) {
            return (line.f0 - c) > fmin and (line.f0 + c) < fmax;
          });
      for (auto& line : bnd | by_range) {
        lines_push_back(out.first, out.second, spec, line, atm, pol);
      }
    } break;
    case FINAL:
      ARTS_USER_ERROR("Bad state")
  }

  if constexpr (keep_position) {
    bubble_sort_by(
        [&](auto& l1, auto& l2) { return out.first[l1].f0 < out.first[l2].f0; },
        out.first,
        out.second);
  } else {
    std::ranges::sort(out.first,
                      [](auto& l1, auto& l2) { return l1.f0 < l2.f0; });
  }

  return out;
}

inline std::pair<Index, Index> find_offset_and_count_of_frequency_range(
    const std::span<const single_shape> lines, Numeric f, Numeric cutoff) {
  if (cutoff < std::numeric_limits<Numeric>::infinity()) {
    auto first = lines.begin();
    auto last = lines.end();

    auto start_range =
        std::lower_bound(first, last, f, [cutoff](auto& line, auto f) {
          return std::abs(cutoff - line.f0) < f;
        });

    auto end_range = start_range;
    while (end_range < last and std::abs(cutoff - end_range->f0) <= f) {
      ++end_range;
    }

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

inline std::pair<std::span<const single_shape>, std::span<const line_pos>>
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

//! A band shape is a collection of single shapes.  The shapes are sorted by frequency.
struct band_shape {
  std::vector<single_shape> lines{};
  Numeric cutoff{-1};

  band_shape() = default;

 public:
  band_shape(std::vector<single_shape>&& ls, const Numeric cut) noexcept
      : lines(ls), cutoff(cut) {}

  band_shape(const SpeciesIsotopeRecord& spec,
             const band& bnd,
             const AtmPoint& atm,
             const Numeric fmin = std::numeric_limits<Numeric>::lowest(),
             const Numeric fmax = std::numeric_limits<Numeric>::max(),
             const zeeman::pol pol = zeeman::pol::no)
      : lines(band_shape_helper<false>(spec, bnd, atm, fmin, fmax, pol).first),
        cutoff(bnd.get_cutoff_frequency()) {}

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

  [[nodiscard]] Complex dH(const ExhaustiveConstVectorView& df0_dH,
                           const Numeric f) const {
    ARTS_ASSERT(static_cast<Size>(df0_dH.size()) == lines.size())

    return std::transform_reduce(
        lines.begin(),
        lines.end(),
        df0_dH.begin(),
        Complex{},
        std::plus<>{},
        [f](auto& ls, auto& d) { return ls.dH(d, f); });
  }

  [[nodiscard]] Complex dT(const ExhaustiveConstVectorView& ds_dT,
                           const ExhaustiveConstVectorView& dG0_dT,
                           const ExhaustiveConstVectorView& dD0_dT,
                           const ExhaustiveConstVectorView& dDV_dT,
                           const Numeric T,
                           const Numeric f) const {
    ARTS_ASSERT(ds_dT.size() == dG0_dT.size())
    ARTS_ASSERT(ds_dT.size() == dD0_dT.size())
    ARTS_ASSERT(ds_dT.size() == dDV_dT.size())
    ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < lines.size(); ++i) {
      out += lines[i].dT(ds_dT[i], dG0_dT[i], dD0_dT[i], dDV_dT[i], T, f);
    }

    return out;
  }

  [[nodiscard]] Complex dT(const ExhaustiveConstVectorView& ds_dVMR,
                           const ExhaustiveConstVectorView& dG0_dVMR,
                           const ExhaustiveConstVectorView& dD0_dVMR,
                           const ExhaustiveConstVectorView& dDV_dVMR,
                           const Numeric f) const {
    ARTS_ASSERT(ds_dVMR.size() == dG0_dVMR.size())
    ARTS_ASSERT(ds_dVMR.size() == dD0_dVMR.size())
    ARTS_ASSERT(ds_dVMR.size() == dDV_dVMR.size())
    ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i = 0; i < lines.size(); ++i) {
      out += lines[i].dVMR(ds_dVMR[i], dG0_dVMR[i], dD0_dVMR[i], dDV_dVMR[i], f);
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

  [[nodiscard]] Complex dDV(const Numeric dX_dDV,
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
        [f, dX_dDV, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dDV(dX_dDV, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dD0(const Numeric dX_dD0,
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
        [f, dX_dD0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dD0(dX_dD0, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dG0(const Numeric dX_dG0,
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
        [f, dX_dG0, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dG0(dX_dG0, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dY(const Numeric dX_dY,
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
        [f, dX_dY, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dY(dX_dY, f)
                     : Complex{};
        });
  }

  [[nodiscard]] Complex dG(const Numeric dX_dG,
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
        [f, dX_dG, &line_key](auto& ls, auto& p) {
          return (p.line == line_key.line and p.spec == line_key.spec)
                     ? ls.dG(dX_dG, f)
                     : Complex{};
        });
  }
};
}  // namespace lbl::voigt::lte
