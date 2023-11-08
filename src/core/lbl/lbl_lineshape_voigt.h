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

  [[nodiscard]] constexpr Complex z(Numeric f) const noexcept {
    return Complex{inv_gd * (f - f0), z_imag};
  }

  [[nodiscard]] Complex F(const Complex z_) const noexcept ;

  [[nodiscard]] Complex F(const Numeric f) const noexcept;

  [[nodiscard]] Complex operator()(const Numeric f) const noexcept;

  [[nodiscard]] constexpr Complex dF(const Complex z_, const Complex F_) const noexcept {
    return 2 * inv_gd * (1i * Constant::inv_pi * inv_gd - z_ * F_);
  }

  [[nodiscard]] Complex dF(const Numeric f) const noexcept;

 private:
  struct zFdF {
    Complex z, F, dF;
  };

  [[nodiscard]] zFdF all(const Numeric f) const noexcept;

 public:
  [[nodiscard]] Complex df(const Numeric f) const noexcept;

  [[nodiscard]] Complex df0(const Complex ds_df0,
                            const Complex dz_df0,
                            const Numeric f) const noexcept ;

  [[nodiscard]] Complex dDV(const Complex dz_dDV,
                            const Numeric f) const noexcept;

  [[nodiscard]] Complex dD0(const Complex dz_dD0,
                            const Numeric f) const noexcept ;

  [[nodiscard]] Complex dG0(const Complex dz_dG0,
                            const Numeric f) const noexcept;

  [[nodiscard]] Complex dH(const Complex dz_dH,
                           const Numeric f) const noexcept;

  [[nodiscard]] Complex dVMR(const Complex ds_dVMR,
                             const Complex dz_dVMR,
                             const Numeric f) const noexcept;

  [[nodiscard]] Complex dT(const Complex ds_dT,
                           const Complex dz_dT,
                           const Numeric f) const noexcept;

  [[nodiscard]] Complex da(const Complex ds_da,
                           const Numeric f) const noexcept;

  [[nodiscard]] Complex de0(const Complex ds_de0,
                            const Numeric f) const noexcept;

  [[nodiscard]] Complex dG(const Complex ds_dG,
                           const Numeric f) const noexcept;

  [[nodiscard]] Complex dY(const Complex ds_dY,
                           const Numeric f) const noexcept;
};

struct line_pos {
  Size line;
  Size spec{std::numeric_limits<Size>::max()};
  Size iz{std::numeric_limits<Size>::max()};
};

Size count_lines(const band& bnd, const zeeman::pol type);

void zeeman_set_back(std::vector<single_shape>& lines,
                            std::vector<line_pos>& pos,
                            const single_shape& s,
                            const line& line,
                            const Numeric H,
                            const Size spec,
                            const zeeman::pol pol,
                            Size& last_single_shape_pos);

 void lines_set(std::vector<single_shape>& lines,
                      std::vector<line_pos>& pos,
                      const SpeciesIsotopeRecord& spec,
                      const line& line,
                      const AtmPoint& atm,
                      const zeeman::pol pol,
                      Size& last_single_shape_pos);

//! Helper for initializing the band_shape
inline void band_shape_helper(std::vector<single_shape>& lines,
                              std::vector<line_pos>& pos,
                              const SpeciesIsotopeRecord& spec,
                              const band& bnd,
                              const AtmPoint& atm,
                              const Numeric fmin,
                              const Numeric fmax,
                              const zeeman::pol pol);

constexpr std::pair<Index, Index> find_offset_and_count_of_frequency_range(
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
constexpr auto frequency_span(const auto& list, const Size start, const Size count) {
  std::span out{list};
  return out.subspan(start, count);
}
}  // namespace detail

template <typename... Ts>
constexpr auto frequency_spans(const Numeric cutoff,
                     const Numeric f,
                     const std::span<const single_shape>& lines,
                     const Ts&... lists) {
  ARTS_ASSERT(lines.size() == (static_cast<Size>(lists.size()) and ...))

  const auto [start, count] =
      find_offset_and_count_of_frequency_range(lines, f, cutoff);
  return std::tuple{detail::frequency_span(lines, start, count),
                    detail::frequency_span(lists, start, count)...};
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

  [[nodiscard]] Complex df0(const ExhaustiveConstComplexVectorView ds_df0,
                            const ExhaustiveConstComplexVectorView dz_df0,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += lines[i].df0(ds_df0[i], dz_df0[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex da(const ExhaustiveConstComplexVectorView ds_da,
                           const Numeric f,
                           const std::vector<Size>& filter) const {
    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += lines[i].da(ds_da[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex de0(const ExhaustiveConstComplexVectorView ds_de0,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += lines[i].de0(ds_de0[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex dDV(const ExhaustiveConstComplexVectorView dz_dDV,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += lines[i].dDV(dz_dDV[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex dD0(const ExhaustiveConstComplexVectorView dz_dD0,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += lines[i].dD0(dz_dD0[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex dG0(const ExhaustiveConstComplexVectorView dz_dG0,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += lines[i].dG0(dz_dG0[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex dY(const ExhaustiveConstComplexVectorView ds_dY,
                           const Numeric f,
                           const std::vector<Size>& filter) const {
    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += lines[i].dY(ds_dY[i], f);
    }

    return out;
  }

  [[nodiscard]] Complex dG(const ExhaustiveConstComplexVectorView ds_dG,
                           const Numeric f,
                           const std::vector<Size>& filter) const {
    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += lines[i].dG(ds_dG[i], f);
    }

    return out;
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
                            const ExhaustiveConstComplexVectorView ds_df0,
                            const ExhaustiveConstComplexVectorView dz_df0,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    const auto [s, cs, ds, dz] =
        frequency_spans(cutoff, f, lines, cut, ds_df0, dz_df0);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += s[i].df0(ds[i], dz[i], f) - cs[i];
    }

    return out;
  }

  void df0(ExhaustiveComplexVectorView cut,
           const ExhaustiveConstComplexVectorView ds_df0,
           const ExhaustiveConstComplexVectorView dz_df0,
           const std::vector<Size>& filter) const {
    for (Size i : filter) {
      cut[i] = lines[i].df0(ds_df0[i], dz_df0[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex da(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView ds_da,
                           const Numeric f,
                           const std::vector<Size>& filter) const {
    const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_da);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += s[i].da(ds[i], f) - cs[i];
    }

    return out;
  }

  void da(ExhaustiveComplexVectorView cut,
          const ExhaustiveComplexVectorView ds_da,
          const std::vector<Size>& filter) const {
    for (Size i : filter) {
      cut[i] = lines[i].da(ds_da[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex de0(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView ds_de0,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_de0);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += s[i].de0(ds[i], f) - cs[i];
    }

    return out;
  }

  void de0(ExhaustiveComplexVectorView cut,
           const ExhaustiveComplexVectorView ds_de0,
           const std::vector<Size>& filter) const {
    for (Size i : filter) {
      cut[i] = lines[i].de0(ds_de0[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex dDV(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView dz_dDV,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    const auto [s, cs, dz] = frequency_spans(cutoff, f, lines, cut, dz_dDV);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += s[i].dDV(dz[i], f) - cs[i];
    }

    return out;
  }

  void dDV(ExhaustiveComplexVectorView cut,
           const ExhaustiveComplexVectorView dz_dDV,
           const std::vector<Size>& filter) const {
    for (Size i : filter) {
      cut[i] = lines[i].dDV(dz_dDV[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex dD0(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView dz_dD0,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    const auto [s, cs, dz] = frequency_spans(cutoff, f, lines, cut, dz_dD0);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += s[i].dD0(dz[i], f) - cs[i];
    }

    return out;
  }

  void dD0(ExhaustiveComplexVectorView cut,
           const ExhaustiveComplexVectorView dz_dD0,
           const std::vector<Size>& filter) const {
    for (Size i : filter) {
      cut[i] = lines[i].dD0(dz_dD0[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex dG0(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView dz_dG0,
                            const Numeric f,
                            const std::vector<Size>& filter) const {
    const auto [s, cs, dz] = frequency_spans(cutoff, f, lines, cut, dz_dG0);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += s[i].dG0(dz[i], f) - cs[i];
    }

    return out;
  }

  void dG0(ExhaustiveComplexVectorView cut,
           const ExhaustiveComplexVectorView dz_dG0,
           const std::vector<Size>& filter) const {
    for (Size i : filter) {
      cut[i] = lines[i].dG0(dz_dG0[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex dY(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView ds_dY,
                           const Numeric f,
                           const std::vector<Size>& filter) const {
    const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_dY);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += s[i].dY(ds[i], f) - cs[i];
    }

    return out;
  }

  void dY(ExhaustiveComplexVectorView cut,
          const ExhaustiveComplexVectorView ds_dY,
          const std::vector<Size>& filter) const {
    for (Size i : filter) {
      cut[i] = lines[i].dY(ds_dY[i], lines[i].f0 + cutoff);
    }
  }

  [[nodiscard]] Complex dG(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView ds_dG,
                           const Numeric f,
                           const std::vector<Size>& filter) const {
    const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_dG);

    Complex out{};  //! Fixme, use zip in C++ 23...

    for (Size i : filter) {
      out += s[i].dG(ds[i], f) - cs[i];
    }

    return out;
  }

  void dG(ExhaustiveComplexVectorView cut,
          const ExhaustiveComplexVectorView ds_dG,
          const std::vector<Size>& filter) const {
    for (Size i : filter) {
      cut[i] = lines[i].dG(ds_dG[i], lines[i].f0 + cutoff);
    }
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
