#pragma once

#include <optional>
#include <variant>
#include <iostream>

#include "scattering/single_scattering_data.h"
#include "scattering/psd.h"

namespace scattering {

using ParticleData = std::variant<
  SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>,
  SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>,
  SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>,
  SingleScatteringData<Numeric, Format::ARO, Representation::Spectral>
  >;

using ScatteringData = std::vector<ParticleData>;


template <Format format, Representation repr>
auto ssd_to_tro_spectral(const ScatteringDataGrids& new_grids,
                         Index l,
                         const SingleScatteringData<Numeric, format, repr>& ssd)
  -> SingleScatteringData<Numeric, Format::TRO, Representation::Spectral> {
      if constexpr (format == Format::ARO) {
        ARTS_USER_ERROR("Cannot convert scattering data from ARO format to TRO format.");
      } else {
        return ssd.to_spectral(l, 0).regrid(new_grids);
      }
};

template <Format format, Representation repr>
auto ssd_to_tro_gridded(const ScatteringDataGrids& new_grids,
                        const SingleScatteringData<Numeric, format, repr>& ssd)
  -> SingleScatteringData<Numeric, Format::TRO, Representation::Gridded> {
  if constexpr (format == Format::ARO) {
    ARTS_USER_ERROR("Cannot convert scattering data from ARO format to TRO format.");
  } else {
    if constexpr (repr == Representation::Gridded) {
      return ssd.regrid(new_grids);
    } else {
      return ssd.to_gridded().regrid(new_grids);
    }
  }
};


class ScatteringHabit;


/** A habit of scattering particles.
 *
 * Represents a set of particles that can be used to infer bulk scattering
 * properties at a given point in the atmosphere.
 */
class ParticleHabit {

  public:

  static ParticleHabit liquid_sphere(const StridedVectorView &t_grid,
                                     const StridedVectorView &f_grid,
                                     const StridedVectorView &diameters,
                                     const ZenithAngleGrid &za_scat_grid) {
    std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> ssd;
    ssd.reserve(diameters.size());
    for (auto ind = 0; ind < Index(diameters.size()); ++ind) {
      ssd.push_back(SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>::liquid_sphere(t_grid,
                                                                                                       f_grid,
                                                                                                       diameters[ind],
                                                                                                       za_scat_grid));
    }
    return ParticleHabit(ssd);
  }


  static ParticleHabit from_legacy_tro(std::vector<::SingleScatteringData> ssd_,
                                       std::vector<::ScatteringMetaData> meta_) {

    std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> ssd;
    ssd.reserve(ssd_.size());
    for (auto ind = 0; ind < Index(ssd_.size()); ++ind) {
      ssd.push_back(SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>::from_legacy_tro(ssd_[ind], meta_[ind]));
    }
    return ParticleHabit(ssd);
  }

  ParticleHabit() {}

  inline Vector get_sizes(SizeParameter param) const {
    Index n_particles = scattering_data.size();
    Vector sizes(n_particles);
    for (Index ind = 0; ind < n_particles; ++ind) {
      auto size = std::visit([&param, &ind](const auto& ssd) {return ssd.get_size(param);}, scattering_data[ind]);
      if (size.has_value()) {
        sizes[ind] = size.value();
      } else {
        ARTS_USER_ERROR("Encountered particle without size information.");
      }
    }
    return sizes;
  }

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
     ParticleHabit(std::vector<SingleScatteringData<Numeric, format, representation>> scattering_data_) {
       if (scattering_data_.size() > 0) {
         for (Index ind = 0; ind < scattering_data_.size(); ++ind) {
           scattering_data.push_back(scattering_data_[ind]);
         }
       }
     }

   template<Format format, Representation representation>
     ParticleHabit(std::vector<SingleScatteringData<Numeric, format, representation>> scattering_data_,
                   ScatteringDataGrids grids_)
       : grids(grids_) {
     if (scattering_data_.size() > 0) {
       scattering_data.resize(scattering_data_.size(), scattering_data_[0]);
       std::transform(scattering_data_.begin(),
                      scattering_data_.end(),
                      scattering_data.begin(),
                      [](const auto& ssd) {return ssd;});
     }
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
    std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>> new_scat_data;
    auto new_grids = ScatteringDataGrids(std::make_shared<Vector>(t_grid),
                                         std::make_shared<Vector>(f_grid));
    auto transform = [&new_grids, &l](const auto& ssd) {return ssd_to_tro_spectral(new_grids, l, ssd);};
    std::transform(scattering_data.begin(),
                   scattering_data.end(),
                   new_scat_data.begin(),
                   [&transform](const ParticleData& pd) {return std::visit(transform, pd);});
    return ParticleHabit(new_scat_data, new_grids);
  }

  ParticleHabit to_tro_gridded(const Vector& t_grid, const Vector& f_grid, const ZenithAngleGrid &za_scat_grid) {
    auto new_grids = ScatteringDataGrids(std::make_shared<const Vector>(t_grid),
                                         std::make_shared<const Vector>(f_grid),
                                         std::make_shared<const ZenithAngleGrid>(za_scat_grid));

    auto transform = [&new_grids](const auto& ssd) {return ssd_to_tro_gridded(new_grids, ssd);};
    std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> new_scattering_data;
    new_scattering_data.reserve(scattering_data.size());
    for (const ParticleData& pd : scattering_data) {
      new_scattering_data.push_back(std::visit(transform, pd));
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

  const ParticleData& operator[](Index ind) const {
    return scattering_data[ind];
  }

  ParticleData& operator[](Index ind) {
    return scattering_data[ind];
  }

  friend ScatteringHabit;

 private:

  ScatteringData scattering_data;
  std::optional<ScatteringDataGrids> grids;

};

}  // namespace scattering
