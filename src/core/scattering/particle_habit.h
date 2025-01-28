#pragma once

#include <optional>
#include <variant>
#include <iostream>

#include "scattering/single_scattering_data.h"

namespace scattering {

using ParticleData = std::variant<
  SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>,
  SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>,
  SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>,
  SingleScatteringData<Numeric, Format::ARO, Representation::Spectral>
  >;

using ScatteringData = std::variant<
    std::vector<
        SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>>,
    std::vector<SingleScatteringData<Numeric,
                                     Format::TRO,
                                     Representation::Spectral>>,
    std::vector<
        SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>>,
    std::vector<SingleScatteringData<Numeric,
                                     Format::ARO,
                                     Representation::Spectral>>>;

/** A habit of scattering particles.
 *
 * Represents a set of particles that can be used to infer bulk scattering
 * properties at a given point in the atmosphere.
 */
class ParticleHabit {

  public:
  static ParticleHabit from_legacy_tro(std::vector<::SingleScatteringData> ssd_,
                                       std::vector<::ScatteringMetaData> meta_) {

    std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> ssd;
    ssd.reserve(ssd_.size());
    for (auto ind = 0; ind < Index(ssd_.size()); ++ind) {
      std::cout << "IND :: " << ind << std::endl;
      ssd.push_back(SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>::from_legacy_tro(ssd_[ind], meta_[ind]));
    }
    return ParticleHabit(ssd);
  }

  ParticleHabit() {}

  /** Create particle habit from TRO scattering data.
   *
   * @param t_grid The temperature grid over which the scattering data
   * is defined.
   * @param f_grid The frequency grid over which the scattering data is
   * defined
   * @param za_scat_grid The scattering zenith angle grid over which the
   * scattering data is defined.
   * @param particle_properties: The particle properties (size, mass, etc.) of
   * each particle.
   * @param scattering_data: Vector containing the single scattering data of
   * each particle.
   */
   template<Format format, Representation representation>
     ParticleHabit(std::vector<SingleScatteringData<Numeric, format, representation>> scattering_data_)
      : scattering_data(scattering_data_) {}


  template <Format format, Representation repr>
    void append_particle(const SingleScatteringData<Numeric, format, repr> &ssd) {
    auto append = [&ssd](auto vec) {vec.push_back(ssd);};
    std::visit(append, scattering_data);
  }

  Index size() {
    return std::visit([](auto vec){return vec.size();}, scattering_data);
  }

  ParticleData operator[](Index ind) {
    return std::visit([&ind](auto vec){return ParticleData(vec[ind]);}, scattering_data);
  }

 private:
  ScatteringData scattering_data;
};
}  // namespace scattering
