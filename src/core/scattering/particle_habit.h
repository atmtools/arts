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

using ScatteringData = std::vector<ParticleData>;


Format get_format(const ParticleData& pd) {
  return std::visit([](const auto& ssd) {return ssd.get_format();}, pd);
}

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
     {
       scattering_data.resize(scattering_data.size());
       std::transform(scattering_data_.begin(),
                      scattering_data_.end(),
                      scattering_data.begin(),
                      [](const auto& ssd) {return ssd;});
     }

   template<Format format, Representation representation>
     ParticleHabit(std::vector<SingleScatteringData<Numeric, format, representation>> scattering_data_,
                   ScatteringDataGrids grids_)
      : grids(grids_) {
       scattering_data.resize(scattering_data.size());
       std::transform(scattering_data_.begin(),
                      scattering_data_.end(),
                      scattering_data.begin(),
                      [](const auto& ssd) {return ssd;});
   }

  template <Format format, Representation repr>
    void append_particle(const SingleScatteringData<Numeric, format, repr> &ssd) {
    if (grids.has_value()) {
      scattering_data.push_back(ssd.regrid(grids));
    } else {
      scattering_data.push_back(ssd);
    }
  }




  ParticleHabit to_tro_spectral(const Vector& t_grid, const Vector& f_grid, Index l) {
    auto new_grids = ScatteringDataGrids(std::make_shared<Vector>(t_grid),
                                         std::make_shared<Vector>(f_grid));

    std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>> new_scat_data;
    auto regrid_and_transform = [&l, &new_grids](const auto& ssd)
      -> SingleScatteringData<Numeric, Format::TRO, Representation::Spectral> {
      if constexpr (ssd.get_format() == Format::ARO) {
        ARTS_USER_ERROR("Cannot convert scattering data from ARO format to TRO format.");
      } else {
        return ssd.to_spectral(l, 0).regrid(new_grids);
      }
    };
    std::transform(scattering_data.begin(),
                   scattering_data.end(),
                   new_scat_data.begin(),
                   [&regrid_and_transform](const ParticleData& pd) {return std::visit(regrid_and_transform, pd);});
    return ParticleHabit(new_scat_data, new_grids);
  }

  ParticleHabit to_tro_gridded(const Vector& t_grid, const Vector& f_grid, const Vector &za_scat_grid) {
    auto new_grids = ScatteringDataGrids(std::make_shared<const Vector>(t_grid),
                                         std::make_shared<const Vector>(f_grid),
                                         std::make_shared<const ZenithAngleGrid>(IrregularZenithAngleGrid(za_scat_grid)));

    auto regrid_and_transform = [&new_grids](const auto& ssd)
    -> SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>{
      if constexpr (ssd.get_format() == Format::ARO) {
        ARTS_USER_ERROR("Cannot convert scattering data from ARO format to TRO format.");
      } else if constexpr (ssd.get_representation() == Representation::Gridded) {
        return ssd.regrid(new_grids);
      } else {
        auto gridded = ssd.to_gridded();
        return gridded.regrid(new_grids);
      }
    };

    std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> new_scattering_data;
    new_scattering_data.reserve(scattering_data.size());
    for (const ParticleData& pd : scattering_data) {
      new_scattering_data.push_back(std::visit(regrid_and_transform, pd));
    }
    return ParticleHabit(new_scattering_data, new_grids);
  }

  ParticleHabit to_aro_spectral(const Vector& t_grid, const Vector& f_grid, Index l, Index m) {

    return {};
  }

  ParticleHabit to_aro_gridded(const Vector& t_grid, const Vector& f_grid, const Vector &za_scat_grid) {
    return {};
  }

  Index size() {
    return scattering_data.size();
  }

  ParticleData operator[](Index ind) {
    return scattering_data[ind];
  }

 private:

  ScatteringData scattering_data;
  std::optional<ScatteringDataGrids> grids;

};

}  // namespace scattering
