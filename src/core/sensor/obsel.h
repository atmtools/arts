#pragma once

#include <matpack.h>
#include <rtepack.h>

#include <tuple>

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
  //! Relevant frequencies
  Vector f_grid_w{};
  AscendingGrid f_grid{};

  //! Relevant propagation from antenna to sensor
  MuelmatVector poslos_grid_w{};
  PosLosVector poslos_grid{};

  //! Sampled polariazation state at antenna
  Stokvec polarization{1.0, 0.0, 0.0, 0.0};

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

  [[nodiscard]] Index ind(const PosLos& poslos) const;

  [[nodiscard]] Index ind(const Numeric& f) const;

  friend std::ostream& operator<<(std::ostream& os, const Obsel& obsel);

  [[nodiscard]] bool ok() const;
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

using SensorPosLos = sensor::PosLos;
using SensorPosLosVector = sensor::PosLosVector;
using SensorObsel = sensor::Obsel;
using ArrayOfSensorObsel = Array<SensorObsel>;
