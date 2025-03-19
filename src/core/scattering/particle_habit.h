#pragma once

#include <optional>
#include <variant>
#include <iostream>

#include "scattering/single_scattering_data.h"
#include "scattering/sht.h"
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

template <Format format, Representation repr>
auto ssd_to_aro_gridded(const ScatteringDataGrids& new_grids,
                        const SingleScatteringData<Numeric, format, repr>& ssd)
  -> SingleScatteringData<Numeric, Format::ARO, Representation::Gridded> {
  if constexpr (format == Format::ARO) {
    if constexpr (repr == Representation::Gridded) {
      return ssd.regrid(new_grids);
    } else {
      return ssd.to_gridded().regrid(new_grids);
    }
  } else {
    if constexpr (repr == Representation::Gridded) {
      return ssd.to_lab_frame(new_grids);
    } else {
      return ssd.to_gridded().to_lab_frame(new_grids);
    }
  }
};


template <Format format, Representation repr>
auto ssd_to_aro_spectral(const ScatteringDataGrids& new_grids,
                         Index l,
                         Index m,
                         const SingleScatteringData<Numeric, format, repr>& ssd)
  -> SingleScatteringData<Numeric, Format::ARO, Representation::Spectral> {
  if constexpr (format == Format::ARO) {
    return ssd.to_spectral(l, m).regrid(new_grids);
  } else {
    if constexpr (repr == Representation::Gridded) {
      return ssd.to_lab_frame(new_grids).to_spectral(l, m);
    } else {
      return ssd.to_gridded().to_lab_frame(new_grids).to_spectral(l, m);
      }
  }
};


/*! Derives a and b for relationship mass = a * x^b

    The parameters a and b are derived by a fit including all data inside the
    size range [fit_start,fit_end].

    The vector x must have been checked to have at least 2 elements.

    An error is thrown if less than two data points are found inside
    [x_fit_start,x_fit_end].

    \param  masses      Size grid
    \param  mass        Particle masses
    \param  x_fit_start Start point of x-range to use for fitting
    \param  x_fit_end   Endpoint of x-range to use for fitting
    \return A pair containing the fitted parameters a and b

  \author Jana Mendrok, Patrick Eriksson, modified by Simon Pfreundschuh
  \date 2017-10-18, 2025-02-12

*/
std::pair<Numeric, Numeric> derive_scat_species_a_and_b(const Vector& sizes,
                                                        const Vector& masses,
                                                        const Numeric& fit_start,
                                                        const Numeric& fit_end);


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
      auto size = std::visit([&param](const auto& ssd) {return ssd.get_size(param);}, scattering_data[ind]);
      if (size.has_value()) {
        sizes[ind] = size.value();
      } else {
        ARTS_USER_ERROR("Encountered particle without size information.");
      }
    }
    return sizes;
  }

/* Workspace method: Doxygen documentation will be auto-generated */
  inline std::tuple<Vector, Numeric, Numeric> get_size_mass_info(SizeParameter size_parameter,
                                                                 const Numeric& fit_start = 0.0,
                                                                 const Numeric& fit_end = 1e9) {

    Index n_particles = scattering_data.size();

    Vector sizes(n_particles);
    Vector masses(n_particles);

    for (Index ind = 0; ind < n_particles; ++ind) {
      auto mass = std::visit([](const auto& ssd) {return ssd.get_mass();}, scattering_data[ind]);
      auto size = std::visit([&size_parameter](const auto& ssd) {return ssd.get_size(size_parameter);}, scattering_data[ind]);
      if (mass.has_value()) {
        masses[ind] = mass.value();
        sizes[ind] = size.value();
      } else {
        ARTS_USER_ERROR("Encountered particle without size information.");
      }
    }

    auto [a, b] = derive_scat_species_a_and_b(sizes, masses, fit_start, fit_end);
    return std::make_tuple(sizes, a, b);

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

  ParticleHabit(const ParticleHabit&) = default;

  template <Format format, Representation repr>
    void append_particle(const SingleScatteringData<Numeric, format, repr> &ssd) {
    if (grids.has_value()) {
      scattering_data.push_back(ssd.regrid(grids));
    } else {
      scattering_data.push_back(ssd);
    }
  }




  ParticleHabit to_tro_spectral(const Vector& t_grid, const Vector& f_grid, Index l) {

    std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>> new_scat_data{};
    auto new_grids = ScatteringDataGrids(std::make_shared<Vector>(t_grid),
                                         std::make_shared<Vector>(f_grid));
    auto transform = [&new_grids, &l](const auto& ssd) {return ssd_to_tro_spectral(new_grids, l, ssd);};
    for (size_t p_ind = 0; p_ind < scattering_data.size(); ++p_ind) {
      new_scat_data.push_back(std::visit(transform, scattering_data[p_ind]));
    }
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

  ParticleHabit to_aro_spectral(const Vector& t_grid,
                                const Vector& f_grid,
                                const Vector& za_inc_grid,
                                Index l,
                                Index m)
  {
    auto sht_ptr = sht::provider.get_instance_lm(l, m);
    auto aa_scat_grid = sht_ptr->get_azimuth_angle_grid();
    auto za_scat_grid = sht_ptr->get_zenith_angle_grid();
    auto new_grids = ScatteringDataGrids(std::make_shared<const Vector>(t_grid),
                                         std::make_shared<const Vector>(f_grid),
                                         std::make_shared<const Vector>(za_inc_grid),
                                         std::make_shared<const Vector>(aa_scat_grid),
                                         std::make_shared<const ZenithAngleGrid>(za_scat_grid));
    auto transform = [&new_grids, &l, &m](const auto& ssd) {return ssd_to_aro_spectral(new_grids, l, m, ssd);};
    std::vector<SingleScatteringData<Numeric, Format::ARO, Representation::Spectral>> new_scattering_data;
    new_scattering_data.reserve(scattering_data.size());
    for (const ParticleData& pd : scattering_data) {
      new_scattering_data.push_back(std::visit(transform, pd));
    }
    return ParticleHabit(new_scattering_data, new_grids);
  }

  ParticleHabit to_aro_gridded(const Vector& t_grid,
                               const Vector& f_grid,
                               const Vector& za_inc_grid,
                               const Vector& aa_scat_grid,
                               const Vector& za_scat_grid) {
    auto new_grids = ScatteringDataGrids(std::make_shared<const Vector>(t_grid),
                                         std::make_shared<const Vector>(f_grid),
                                         std::make_shared<const Vector>(aa_scat_grid),
                                         std::make_shared<const Vector>(za_scat_grid),
                                         std::make_shared<const ZenithAngleGrid>(za_scat_grid));
    auto transform = [&new_grids](const auto& ssd) {return ssd_to_aro_gridded(new_grids, ssd);};
    std::vector<SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>> new_scattering_data;
    new_scattering_data.reserve(scattering_data.size());
    for (const ParticleData& pd : scattering_data) {
      new_scattering_data.push_back(std::visit(transform, pd));
    }
    return ParticleHabit(new_scattering_data, new_grids);
  }

  Index size() const {
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
