/** \file particle_habit.h
 *
 * Contains the ParticleHabit class which represents a habit of particles that
 * scattering radiation. The habit is represented by  a collection of scattering
 * data for different particle sizes.
 *
 * @author Simon Pfreundschuh, 2020
 */
#pragma once

#include <random>

#include <scattering/single_scattering_data.h>
#include <scattering/particle.h>

namespace scattering {

// pxx :: export
/** Particle habit
 *
 * A particle habit represents a collection of scattering particles.
 */
class ParticleHabit {

public:

  /// Load particle from ARTS SSDB file.
  static ParticleHabit from_ssdb(std::string path);

  /// Create an empty ParticleHabit.
  ParticleHabit() {};

  /// Create a ParticleHabit from given particles
  ParticleHabit(const std::vector<scattering::Particle> &particles)
      : particles_(particles) {}

  static ParticleHabit liquid_spheres(
      math::Vector<double> f_grid,
      math::Vector<double> t_grid,
      math::Vector<double> lat_scat,
      math::Vector<double> radii
      ) {
      std::vector<Particle> particles;
      particles.reserve(radii.size());
      for (double radius : radii) {
          particles.emplace_back(
              Particle::liquid_sphere(f_grid, t_grid, lat_scat, radius)
              );
      }
      return ParticleHabit(particles);
  }

  static ParticleHabit deserialize(std::istream &input) {
      size_t n_particles;
      input.read(reinterpret_cast<char*>(&n_particles), sizeof(size_t));
      std::vector<Particle> particles{};
      particles.reserve(n_particles);
      for (size_t i = 0; i < n_particles; ++i) {
          particles.emplace_back(Particle::deserialize(input));
      }
      return ParticleHabit(particles);
  }

  std::ostream& serialize(std::ostream &output) const {
      size_t n_particles = size();
      output.write(reinterpret_cast<const char*>(&n_particles), sizeof(size_t));
      for (const Particle &part : particles_) {
          part.serialize(output);
      }
  }

  /// The number of particle that make up the particle habit.
  size_t size() const {
      return particles_.size();
  }

  /// Return vector contatining volume equivalent diameter of particles in the
  /// habit.
  math::Vector<double> get_d_eq() const {
    math::Vector<double> result(particles_.size());
    for (size_t i = 0; i < particles_.size(); ++i) {
      result[i] = particles_[i].get_d_eq();
    }
    return result;
  }

  /// Return vector contatining the maximum diameter of particles in the habit.
  math::Vector<double> get_d_max() const {
    math::Vector<double> result(particles_.size());
    for (size_t i = 0; i < particles_.size(); ++i) {
      result[i] = particles_[i].get_d_max();
    }
    return result;
  }

  /// Return vector contatining the mass of the particles in the habit.
  math::Vector<double> get_mass() const {
      math::Vector<double> result(particles_.size());
      for (size_t i = 0; i < particles_.size(); ++i) {
          result[i] = particles_[i].get_mass();
      }
      return result;
  }

  const std::vector<scattering::Particle> & get_particles() const {
      return particles_;
  }

  /// The minimum stokes dimension in the data.
  Index get_stokes_dim() const {
      auto stokes_dim = [](const scattering::Particle &p) -> Index {
          return p.get_stokes_dim();
      };
      std::vector<Index> stokes_dims;
      stokes_dims.reserve(particles_.size());
      std::transform(particles_.begin(),
                     particles_.end(),
                     stokes_dims.begin(),
                     stokes_dim);
      return *std::min_element(stokes_dims.begin(), stokes_dims.end());
  }

  /** Scattering data describing habit.
   * @return A vector containing the single scattering data describing the particles in
   * in the habit.
   */
  const SingleScatteringData &get_single_scattering_data(size_t index) const {
    return particles_[index].get_data();
  }

  // pxx :: hide
  /** Interpolate single scattering data along frequencies.
   *
   * @param f_grid The frequency grid to which to interpolate the data.
   * @return A new particle habit with the data interpolated to the given
   * frequency grid.
   */
  ParticleHabit interpolate_frequency(math::ConstVectorPtr<double> f_grid) {
      std::vector<scattering::Particle> new_data{};
      new_data.reserve(particles_.size());
      for (size_t i = 0; i < particles_.size(); ++i) {
          new_data.push_back(particles_[i].interpolate_frequency(f_grid));
      }
      return ParticleHabit(new_data);
  }

  ParticleHabit interpolate_frequency(math::Vector<double> f_grid) {
      auto f_grid_ptr = std::make_shared<math::Vector<double>>(f_grid);
      return interpolate_frequency(f_grid_ptr);
  }

  /** Regrid the data to SHTns-compatible grids.
   *
   * In order to perform SHT transforms, the data must be provided on a regular
   * Fejer quadrature grid. This method will interpolate the data to the grids
   * required for correct SHT transformation. The grids are computed by linear
   * interpolation of the existing grids to the corresponding Fejer grids with
   * the same number of grid nodes.
   *
   * @return A new ParticleHabit object with the data interpolated to the angular
   * grids required for SHT transforms.
   */
    ParticleHabit regrid() const {
        std::vector<scattering::Particle> new_particles{};
        new_particles.reserve(particles_.size());
        for (size_t i = 0; i < particles_.size(); ++i) {
            new_particles.push_back(particles_[i].regrid());
        }
        return ParticleHabit(new_particles);
    }

    /** Transform single scattering data in habit to spectral representation.
     *
     *
     */
    ParticleHabit to_spectral(Index l_max, Index m_max) {
        std::vector<scattering::Particle> new_particles{};
        new_particles.reserve(particles_.size());
        for (size_t i = 0; i < particles_.size(); ++i) {
            new_particles.push_back(particles_[i].to_spectral(l_max, m_max));
        }
        return ParticleHabit(new_particles);
    }

    ParticleHabit to_spectral() {
        std::vector<scattering::Particle> new_particles{};
        new_particles.reserve(particles_.size());
        for (size_t i = 0; i < particles_.size(); ++i) {
            new_particles.push_back(particles_[i].to_spectral());
        }
        return ParticleHabit(new_particles);
    }

    ParticleHabit set_stokes_dim(Index n) const {
        std::vector<scattering::Particle> new_particles{};
        new_particles.reserve(particles_.size());
        for (size_t i = 0; i < particles_.size(); ++i) {
            new_particles.push_back(particles_[i].set_stokes_dim(n));
        }
        return ParticleHabit(new_particles);
    }


    // pxx :: hide
    ParticleHabit to_gridded(math::ConstVectorPtr<double> lon_inc,
                             math::ConstVectorPtr<double> lat_inc,
                             math::ConstVectorPtr<double> lon_scat,
                             ConstLatitudeGridPtr<double> lat_scat) {
      std::vector<scattering::Particle> new_particles{};
      new_particles.reserve(particles_.size());
      for (size_t i = 0; i < particles_.size(); ++i) {
          new_particles.push_back(particles_[i].to_gridded(lon_inc,
                                                           lat_inc,
                                                           lon_scat,
                                                           lat_scat));
      }
      return ParticleHabit(new_particles);
    }

    ParticleHabit to_gridded(math::Vector<double> lon_inc,
                             math::Vector<double> lat_inc,
                             math::Vector<double> lon_scat,
                             math::Vector<double> lat_scat) {
      auto lon_inc_ptr = std::make_shared<math::Vector<double>>(lon_inc);
      auto lat_inc_ptr = std::make_shared<math::Vector<double>>(lat_inc);
      auto lon_scat_ptr = std::make_shared<math::Vector<double>>(lon_scat);
      auto lat_scat_ptr = std::make_shared<IrregularLatitudeGrid<double>>(lat_scat);
      return to_gridded(lon_inc_ptr, lat_inc_ptr, lon_scat_ptr, lat_scat_ptr);
    }

    // pxx :: hide
    ParticleHabit to_lab_frame(math::ConstVectorPtr<double> lat_inc,
                               math::ConstVectorPtr<double> lon_scat,
                               ConstLatitudeGridPtr<double> lat_scat,
                               Index stokes_dim) {
        std::vector<scattering::Particle> new_particles{};
        new_particles.reserve(particles_.size());
        for (size_t i = 0; i < particles_.size(); ++i) {
            new_particles.push_back(particles_[i].to_lab_frame(lat_inc, lon_scat, lat_scat, stokes_dim));
        }
        return ParticleHabit(new_particles);
    }

    ParticleHabit to_lab_frame(math::Vector<double> lat_inc,
                               math::Vector<double> lon_scat,
                               math::Vector<double> lat_scat,
                               Index stokes_dim) {
        auto lat_inc_ptr = std::make_shared<math::Vector<double>>(lat_inc);
        auto lon_scat_ptr = std::make_shared<math::Vector<double>>(lon_scat);
        auto lat_scat_ptr = std::make_shared<IrregularLatitudeGrid<double>>(lat_scat);
        return to_lab_frame(lat_inc_ptr, lon_scat_ptr, lat_scat_ptr, stokes_dim);
    }

    ParticleHabit to_lab_frame(Index n_lat_inc, Index n_lon_scat, Index stokes_dim) const {
        std::vector<scattering::Particle> new_particles{};
        new_particles.reserve(particles_.size());
      for (size_t i = 0; i < particles_.size(); ++i) {
          new_particles.push_back(particles_[i].to_lab_frame(n_lat_inc, n_lon_scat, stokes_dim));
      }
      return ParticleHabit(new_particles);
    }

    // pxx :: hide
    ParticleHabit downsample_scattering_angles(math::ConstVectorPtr<double> lon_scat,
                                               ConstLatitudeGridPtr<double> lat_scat) const {
        std::vector<scattering::Particle> new_particles{};
        new_particles.reserve(particles_.size());
        for (size_t i = 0; i < particles_.size(); ++i) {
            new_particles.push_back(particles_[i].downsample_scattering_angles(lon_scat,
                                                                               lat_scat));
        }
        return ParticleHabit(new_particles);
    }

    ParticleHabit downsample_scattering_angles(const math::Vector<double> &lon_scat,
                                               const math::Vector<double> &lat_scat) const {
      auto lon_scat_ptr = std::make_shared<math::Vector<double>>(lon_scat);
      auto lat_scat_ptr = std::make_shared<IrregularLatitudeGrid<double>>(lat_scat);
      return downsample_scattering_angles(lon_scat_ptr, lat_scat_ptr);
    }

    /** Calculate bulk scattering properties
     *
     * Calculates bulk scattering properties for a given temperature and particle number
     * distribution.
     *
     * @param The atmospheric temperature in K
     * @param pnd Vector containing the number of particles of each of the species in the habit.
     * @return The scattering properties corresponding to the sum of the particle in the habit
     * multiplied by the number given in pnd.
     */
    SingleScatteringData calculate_bulk_properties(
        double temperature,
        math::ConstVectorRef<double> pnd) {
        assert(static_cast<size_t>(pnd.size()) == particles_.size());

      auto temperature_vector = std::make_shared<math::Vector<double>>(1);
      (*temperature_vector)[0] = temperature;
      auto result = particles_[0].interpolate_temperature(temperature);
      result *= pnd[0];

      for (Index i = 1; i < pnd.size(); ++i) {
        auto data = particles_[i].interpolate_temperature(temperature);
        result += data * pnd[i];
      }
      return result;
    }

    math::Matrix<double> get_phase_matrix(
        double frequency,
        double temperature,
        double lon_inc,
        double lat_inc,
        double lon_scat,
        double lat_scat,
        math::ConstVectorRef<double> pnd,
        Index stokes_dim
        ) {
        math::Matrix<double> result;
        result = particles_[0].get_phase_matrix(
                frequency, temperature,
                lon_inc, lat_inc, lon_scat, lat_scat,
                stokes_dim
            ) * pnd[0];
        for (Index i = 1; i < pnd.size(); ++i) {
            auto data = particles_[i].get_phase_matrix(
                frequency, temperature,
                lon_inc, lat_inc, lon_scat, lat_scat,
                stokes_dim
                );
            result += data * pnd[i];
        }
        return result;
    }

    math::Matrix<double> get_extinction_matrix(
        double frequency,
        double temperature,
        double lon_inc,
        double lat_inc,
        math::ConstVectorRef<double> pnd,
        Index stokes_dim
        ) {
        math::Matrix<double> result = math::Matrix<double>::Constant(stokes_dim, stokes_dim, 0.0);
        for (Index i = 0; i < pnd.size(); ++i) {
            auto data = particles_[i].get_extinction_matrix(
                frequency, temperature,
                lon_inc, lat_inc,
                stokes_dim
                );
            result += data * pnd[i];
        }
        return result;
    }

    math::Vector<double> get_absorption_vector(
        double frequency,
        double temperature,
        double lon_inc,
        double lat_inc,
        math::ConstVectorRef<double> pnd,
        Index stokes_dim
        ) {
        math::Vector<double> result = math::Vector<double>::Constant(stokes_dim, 0.0);
        for (Index i = 0; i < pnd.size(); ++i) {
            auto data = particles_[i].get_absorption_vector(
                frequency, temperature,
                lon_inc, lat_inc,
                stokes_dim
                );
            result += data * pnd[i];
        }
        return result;
    }

    double get_phase_function_maximum_inc(
        double frequency,
        double temperature,
        double lon_scat,
        double lat_scat,
        math::ConstVectorRef<double> pnd) {
        double max = 0.0;
        for (Index i = 0; i < pnd.size(); ++i) {
            max += pnd[i] * particles_[i].get_phase_function_max_inc(
                frequency,
                temperature,
                lon_scat,
                lat_scat
                );
        }
        return max;
    }

    double get_phase_function_maximum_scat(
        double frequency,
        double temperature,
        double lon_scat,
        double lat_scat,
        math::ConstVectorRef<double> pnd) {
        double max = 0.0;
        for (Index i = 0; i < pnd.size(); ++i) {
            max += pnd[i] * particles_[i].get_phase_function_max_scat(
                frequency,
                temperature,
                lon_scat,
                lat_scat
                );
        }
        return max;
    }

    std::pair<double, double> sample_phase_matrix_inc(
        double frequency,
        double temperature,
        double lon_scat,
        double lat_scat,
        math::ConstVectorRef<double> pnd,
        Index stokes_dim
        ) {
        double norm = 0.0;
        for (Index i = 0; i < pnd.size(); ++i) {
            norm += pnd[i] * particles_[i].get_phase_function_max_inc(
                frequency,
                temperature,
                lon_scat,
                lat_scat
                );
        }

        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator{static_cast<std::default_random_engine::result_type>(seed)};
        std::uniform_real_distribution<double> lon_dist(0, 2 * M_PI);
        std::uniform_real_distribution<double> lat_dist(-1.0, 1.0);
        std::uniform_real_distribution<double> r_dist(0, 1.0);
        bool found = false;

        double lon = 0.0;
        double lat = 0.0;

        while (!found) {
            lon = lon_dist(generator);
            lat = acos(lat_dist(generator));

            math::Matrix<double> pha_mat = math::Matrix<double>::Zero(
                stokes_dim,
                stokes_dim
                );

            for (Index i = 0; i < pnd.size(); ++i) {
                pha_mat += pnd[i] * particles_[i].get_phase_matrix(
                    frequency,
                    temperature,
                    lon,
                    lat,
                    lon_scat,
                    lat_scat,
                    stokes_dim
                    );
            }
            double r = r_dist(generator);
            if (r * norm < pha_mat(0, 0)) found = true;
        }
        return std::make_pair(lon, lat);
    }

    std::pair<double, double> sample_phase_matrix_scat(
        double frequency,
        double temperature,
        double lon_inc,
        double lat_inc,
        math::ConstVectorRef<double> pnd,
        Index stokes_dim
        ) {
        double norm = 0.0;
        for (Index i = 0; i < pnd.size(); ++i) {
            norm += pnd[i] * particles_[i].get_phase_function_max_scat(
                frequency,
                temperature,
                lon_inc,
                lat_inc
                );
        }

        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator{static_cast<std::default_random_engine::result_type>(seed)};
        std::uniform_real_distribution<double> lon_dist(0, 2.0 * M_PI);
        std::uniform_real_distribution<double> lat_dist(-1.0, 1.0);
        std::uniform_real_distribution<double> r_dist(0, 1.0);
        bool found = false;

        double lon = 0.0;
        double lat = 0.0;

        while (!found) {
            lon = lon_dist(generator);
            lat = acos(lat_dist(generator));

            math::Matrix<double> pha_mat = math::Matrix<double>::Zero(
                stokes_dim,
                stokes_dim
                );
            for (Index i = 0; i < pnd.size(); ++i) {
                pha_mat += pnd[i] * particles_[i].get_phase_matrix(
                    frequency,
                    temperature,
                    lon_inc,
                    lat_inc,
                    lon,
                    lat,
                    stokes_dim
                    );
            }
            double r = r_dist(generator);
            if (r * norm < pha_mat(0, 0)) found = true;
        }
        return std::make_pair(lon, lat);
    }

    friend std::ostream& operator<<(std::ostream& out, const ParticleHabit&);

private:

    std::vector<scattering::Particle> particles_;
};


}
