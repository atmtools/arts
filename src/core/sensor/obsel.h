#pragma once

#include <matpack.h>
#include <rtepack.h>

#include <algorithm>
#include <tuple>

#include "format_tags.h"
#include "matpack_constexpr.h"
#include "sorted_grid.h"

namespace sensor {
struct PosLos {
  Vector3 pos;
  Vector2 los;

  constexpr auto operator==(const PosLos& other) const {
    return std::tie(pos, los) == std::tie(other.pos, other.los);
  };

  constexpr auto operator!=(const PosLos& other) const {
    return std::tie(pos, los) != std::tie(other.pos, other.los);
  };

  constexpr auto operator<=(const PosLos& other) const {
    return std::tie(pos, los) <= std::tie(other.pos, other.los);
  };

  constexpr auto operator>=(const PosLos& other) const {
    return std::tie(pos, los) >= std::tie(other.pos, other.los);
  };

  constexpr auto operator<(const PosLos& other) const {
    return std::tie(pos, los) < std::tie(other.pos, other.los);
  };

  constexpr auto operator>(const PosLos& other) const {
    return std::tie(pos, los) > std::tie(other.pos, other.los);
  };

  friend std::ostream& operator<<(std::ostream& os, const PosLos& poslos);
};

using PosLosVector = matpack::matpack_data<PosLos, 1>;

struct Obsel {
  //! Frequency grid weight, should likely sum to 1.0 - must have same size as f_grid
  Vector f_grid_w{};

  //! Frequency grid, must be ascending - must have same size as f_grid_w
  AscendingGrid f_grid{};

  //! Position and line of sight grid weight - must have same size as poslos_grid
  MuelmatVector poslos_grid_w{};

  //! Position and line of sight grid of the sensor - must have same size as poslos_grid_w
  PosLosVector poslos_grid{};

  //! Sampled polariazation state at antenna
  Stokvec polarization{1.0, 0.0, 0.0, 0.0};

  //! Constant indicating that the frequency or poslos is not found in the grid
  constexpr static Index dont_have = -1;

  /** Returns the state vector that should be sampled for a given pos-los and frequency.
   *
   * The dot-product of this output and the corresponding radiance Stokvec should
   * give the observed radiance for this weight.
   *
   * Returns [0,0,0,0] if either the frequency or the pos-los is not found.
   * 
   * @param f Frequency [Hz]
   * @param pl Pos-los [Altitude m, Latitude deg, Longitude deg], [Zenith deg, Azimuth deg]
   * @return Stokvec 
   */
  [[nodiscard]] Stokvec at(const Numeric& f, const PosLos& pl) const;

  /** As other at() function but with pure indices.
   *
   * Returns [0,0,0,0] if either the frequency or the pos-los indices are the value of dont_have. 
   * 
   * @param freq_ind Freqeuncy index [-1, f_grid.size()-1]
   * @param poslos_ind Frequency index [-1, poslos.size()-1]
   * @return Stokvec 
   */
  [[nodiscard]] Stokvec at(Index freq_ind, Index poslos_ind) const;

  /** Returns the index of the frequency in the frequency grid.
   *
   * Returns dont_have if the frequency is not found.
   * 
   * @param f Frequency [Hz]
   * @return Index 
   */
  [[nodiscard]] Index ind(const Numeric& f) const;

  /** Returns the index of the pos-los in the pos-los grid.
   *
   * Returns dont_have if the pos-los is not found.
   * 
   * @param poslos Pos-los [Altitude m, Latitude deg, Longitude deg], [Zenith deg, Azimuth deg]
   * @return Index 
   */
  [[nodiscard]] Index ind(const PosLos& poslos) const;

  friend std::ostream& operator<<(std::ostream& os, const Obsel& obsel);

  [[nodiscard]] bool ok() const;

  /** Set the frequency response to Dirac (single frequency, 1 weight)
   * 
   * @param f0 The frequency [Hz]
   */
  void set_frequency_dirac(const Numeric& f0);

  /** Set the frequency response to Boxcar (single frequency, 1 weight)
   *
   * A boxcar is a constant response over a range of frequencies.
   * The frequenct range is f0 +/- width/2.
   *
   * If the number of points N is 1, a dirac expression is used and width is ignored.
   *
   * If more N is larger than 1, the upper and lower edges of the boxcar are included.
   *
   * If N is less then 1, an error is thrown.
   * 
   * @param f0 The frequency [Hz]
   * @param width The width of the boxcar [Hz]
   * @param N The number of points in the boxcar.
   */
  void set_frequency_boxcar(const Numeric& f0,
                            const Numeric& width,
                            const Index& N);

  /** Set the frequency response to Boxcar (single frequency, 1 weight)
   *
   * A boxcar is a constant response over a range of frequencies.
   * The frequenct range is f0 +/- width/2.
   * 
   * @param f0 The frequency [Hz]
   * @param width The width of the boxcar [Hz]
   * @param f_grid An external frequency grid on which the boxcar element is placed
   * @param error_if_empty If true, an error is thrown if the grid is empty
   */
  void set_frequency_boxcar(const Numeric& f0,
                            const Numeric& width,
                            const AscendingGrid& f_grid,
                            const bool error_if_empty = true);

  /**  Set the frequency response to Gaussian
   *
   * Set the frequency response to a Gaussian function.
   * This means std::exp(- 4 * Constant::ln_2 * Math::pow2((f - f0) / fwhm))
   * for every valid frequency f.  The weights are normalized so that they sum to 1.
   *
   * There will be 2 * Nfwhm * (Nhwhm - 1) + 1 points in the frequency response.
   * The last and first are at f0 +- Nfwhm * fwhm / 2. In between, there is
   * a division so that there are Nhwhm points per inclusive fwhm / 2 range.
   *
   * As an example, f0 = 100 GHz, fwhm = 10 GHz, Nfwhm = 3, Nhwhm = 5 will
   * yield a frequency grid of: 
   * [85.0, 86.25, 87.5, 88.75, 90.0, 91.25, 92.5, 93.75, 95.0, 96.25, 97.5, 98.75, 100.0, 101.25, 102.5, 103.75, 105.0, 106.25, 107.5, 108.75, 110.0, 111.25, 112.5, 113.75, 115.0] GHz.
   * and a relative weight of the inverse of:
   * [512.0, 189.03, 76.11, 33.42, 16.0, 8.35, 4.76, 2.95, 2.0, 1.48, 1.19, 1.04, 1.0, 1.04, 1.19, 1.48, 2.0, 2.95, 4.76, 8.35, 16.0, 33.42, 76.11, 189.03, 512.0]
   * (these relative weights are the maximum of the distribution divided by the actual weights, used here just to demonstrate which of the 2^N points in the FWHM range we hit)
   * 
   * @param f0 The central frequency [Hz]
   * @param fwhm The full width at half maximum [Hz]
   * @param Nfwhm The number of points in the FWHM
   * @param Nhwhm The number of points per half width
   */
  void set_frequency_gaussian(const Numeric& f0,
                              const Numeric& fwhm,
                              const Index& Nfwhm = 3,
                              const Index& Nhwhm = 5);

  /**  Set the frequency response to Gaussian
   *
   * See other overloads for details on Gaussian distribution method.
   *
   * The provided frequency grid is used entirely.
   * 
   * @param f0 The central frequency [Hz]
   * @param fwhm The full width at half maximum [Hz]
   * @param f_grid The frequency grid on which the Gaussian is placed
   */
  void set_frequency_gaussian(const Numeric& f0,
                              const Numeric& fwhm,
                              const AscendingGrid& f_grid);

  /** Set the frequency distribution to that of an LO chain
   *
   * This is a local oscillator style channel selection frequency grid.
   * The f0s represents the LO chain.  The channels are placed in the
   * center of the resulting chain.  The length of f0s is the number of
   * applied LOs.  For example, if f0s = [100, 10, 1, 0.1] GHz there will be
   * channels at the end at all 8 combinations of 100 ± 10 ± 1 ± 0.1 GHz,
   * that is at [88.9, 89.1, 90.9, 91.1, 108.9, 109.1, 110.9, 111.1] GHz.
   *
   * The width is the width of each channel.  The number of points per
   * channel is N.  If N is 1, only the central frequency is used to
   * produce the channel.  If N is larger than 1, the channel is
   * represented by a boxcar function with the given width around the
   * central frequency of the LO chain.
   *
   * Finally, filter is given as a way to apply a bandpass filter to the
   * LO chain.  The filter is string of U, L, or anything else.  U means
   * that the upper sideband is filtered in.  L means that the lower
   * sideband is filtered in.  Anything else means that the full bandpass
   * is used.  So, for example, if f0s = [100, 10, 1, 0.1] GHz and filter
   * is "UUU" the resulting channels will be at [111.1] GHz.  If filter
   * instead was "UXX" the resulting channels would be at
   * [108.9, 109.1, 110.9, 111.1] GHz.  Leave it empty to use the full
   * bandpass.
   * 
   * @param f0s The frequency of the LOs
   * @param width The width of the channel [Hz]
   * @param N The number of points per channel
   * @param filter The sideband filter - must be one shorter than f0s or empty
   */
  void set_frequency_lochain(const DescendingGrid& f0s,
                             const Numeric& width,
                             const Index& N,
                             const String& filter       = {},
                             const Numeric& lower_width = -0.5,
                             const Numeric& upper_width = 0.5);

  /** Set the frequency distribution to that of an LO chain
   *
   * See other overload for details on the LO frequency grid distribution.
   *
   * All frequencies within the f_grid are mathced to fitting channels and added
   * to this element's frequency grid.  If any single channel has no frequencies,
   * then an error is thrown if error_if_empty is true.
   *
   * @param f0s The frequency of the LOs
   * @param width The width of the channel [Hz]
   * @param f_grid An external frequency grid on which the boxcar element is placed
   * @param filter The sideband filter - must be one shorter than f0s or empty
   * @param error_if_empty If true, an error is thrown if the grid is empty
   */
  void set_frequency_lochain(const DescendingGrid& f0s,
                             const Numeric& width,
                             const AscendingGrid& f_grid,
                             const String& filter       = {},
                             const Numeric& lower_width = -0.5,
                             const Numeric& upper_width = 0.5,
                             const bool error_if_empty  = true);

  //! Normalize the frequencies weights so that they sum to 1 - requires existing weights
  void normalize_frequency_weights();

  //! Remove weights above and below the cutoff - requires existing weights, returns true if anything was removed, cannot remove last element unless relative_cutoff is 1.0

  /** Remove weights above and below the cutoff - requires existing weights
   *
   * This will not work if any of the weights are negative or the weights are not normalized.
   *
   * In relative mode cutoff=0.1 keeps the top 90% of the weights.
   *
   * In absolute mode cutoff=0.1 removes weights with a contribution of less than 10%.
   * 
   * @param cutoff The cutoff value
   * @param relative Whether the cutoff is a percentage of the total weight or an absolute value
   * @return true if anything was removed
   * @return false if the cutoff is the same
   */
  bool cutoff_frequency_weights(const Numeric& cutoff, bool relative = true);
};

AscendingGrid& collect_f_grid(AscendingGrid& f_grid,
                              const Array<Obsel>& obsels,
                              const PosLos& poslos);
AscendingGrid collect_f_grid(const Array<Obsel>& obsels, const PosLos& poslos);

AscendingGrid& collect_f_grid(AscendingGrid& f_grid,
                              const Array<Obsel>& obsels);
AscendingGrid collect_f_grid(const Array<Obsel>& obsels);

PosLosVector& collect_poslos(PosLosVector& poslos, const Array<Obsel>& obsels);
PosLosVector collect_poslos(const Array<Obsel>& obsels);

Index max_frequency_size(const Array<Obsel>& obsels);

void sumup(Vector& out,
           const StokvecVector& in,
           const AscendingGrid& freqs,
           const Array<Obsel>& obsels,
           const PosLos& poslos);
void sumup(Vector& out,
           Matrix& out_jac,
           const StokvecVector& in,
           const StokvecMatrix& in_jac,
           const AscendingGrid& freqs,
           const Array<Obsel>& obsels,
           const PosLos& poslos);

void exhaustive_sumup(Vector& out,
                      const StokvecVector& in,
                      const Array<Obsel>& obsels,
                      const PosLos& poslos);
void exhaustive_sumup(Vector& out,
                      Matrix& out_jac,
                      const StokvecVector& in,
                      const StokvecMatrix& in_jac,
                      const Array<Obsel>& obsels,
                      const PosLos& poslos);

bool all_ok(const Array<Obsel>& obsels);

bool is_exhaustive_like(const Array<Obsel>& obsels);

std::ostream& operator<<(std::ostream& os, const Array<Obsel>& obsel);
}  // namespace sensor

using SensorPosLos       = sensor::PosLos;
using SensorPosLosVector = sensor::PosLosVector;
using SensorObsel        = sensor::Obsel;
using ArrayOfSensorObsel = Array<SensorObsel>;

template <>
struct std::formatter<SensorPosLos> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  constexpr void make_compat(std::formatter<Ts>&... xs) const {
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
  FmtContext::iterator format(const SensorPosLos& v, FmtContext& ctx) const {
    std::formatter<Vector3> fmt3{};
    std::formatter<Vector2> fmt2{};
    make_compat(fmt3, fmt2);

    fmt3.format(v.pos, ctx);

    if (tags.comma) std::ranges::copy(","sv, ctx.out());
    std::ranges::copy(" "sv, ctx.out());

    return fmt2.format(v.los, ctx);
  }
};
template <>

struct std::formatter<SensorObsel> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  constexpr void make_compat(std::formatter<Ts>&... xs) const {
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
  FmtContext::iterator format(const SensorObsel& v, FmtContext& ctx) const {
    std::formatter<Vector> f_grid_w{};
    std::formatter<AscendingGrid> f_grid{};
    std::formatter<MuelmatVector> poslos_grid_w{};
    std::formatter<SensorPosLosVector> poslos_grid{};
    std::formatter<Stokvec> polarization{};
    make_compat(f_grid_w, f_grid, poslos_grid_w, poslos_grid, polarization);

    std::ranges::copy("Obsel:"sv, ctx.out());

    std::ranges::copy("\n  frequency grid:                 "sv, ctx.out());
    f_grid.format(v.f_grid, ctx);

    std::ranges::copy("\n  pos-los grid:                   "sv, ctx.out());
    poslos_grid.format(v.poslos_grid, ctx);

    std::ranges::copy("\n  polarization:                   "sv, ctx.out());
    polarization.format(v.polarization, ctx);

    std::ranges::copy("\n  frequency grid weights:         "sv, ctx.out());
    f_grid_w.format(v.f_grid_w, ctx);

    std::ranges::copy("\n  pos-los grid polarized weigths: "sv, ctx.out());
    poslos_grid_w.format(v.poslos_grid_w, ctx);

    return ctx.out();
  }
};
