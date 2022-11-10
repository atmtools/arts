/** \file arts_ssdb.h
 *
 * Provides an interface to read single scattering data from the ARTS single
 * scattering database.
 *
 * @author Simon Pfreundschuh, 2020
 */
#pragma once

#include <regex>
#include <set>
#include <utility>
#include <filesystem>

#include "netcdf.hpp"

#include <scattering/math.h>
#include <scattering/array.h>
#include <scattering/particle.h>
#include <scattering/particle_habit.h>

namespace scattering {

namespace arts_ssdb {

namespace detail {

////////////////////////////////////////////////////////////////////////////////
// Helper functions.
////////////////////////////////////////////////////////////////////////////////

/** Extract temperature and frequency from group name.
 * @param group_name The name of one of the NetCDF groups in an ASSDB particle
 * file.
 * @return Pair(temp, freq) containing the extracted temperatures and
 * frequency.
 */
std::pair<double, double> match_temp_and_freq(std::string group_name);

/** Extract particle metadata from filename.
 * @param std::filesystem::path Path obejct pointing to the NetCDF4 file to
 * read.
 * @return tuple (match, d_eq, d_max, m) containing
 *    - match: Flag indicating whether the filename matches the ARTS SSDB
 *      pattern.
 *    - d_eq: The volume-equivalent diameter
 *    - d_max: The maximum diameter.
 *    - m: The mass of the particle.
 */
std::tuple<bool, double, double, double> match_particle_properties(
    std::filesystem::path path);

/** Try to extract habit name from file path..
 *
 * Searches for the closest parent folder whose name matches the
 * pattern <habit_name>_Id<id> and returns the matching habit name.
 *
 * @param std::filesystem::path Path obejct pointing to the NetCDF4 file to
 * read.
 * @return The extracted habit
 */
std::string match_habit_name(
    std::filesystem::path path);

 /** Indirect sort w.r.t. equivalent diameter.
 *
 * @param d_eq Vector containing the water equivalent diameter of the particles.
 * @param d_max Vector containing the maximum diameter of the particles.
 * @param m Vector containing the masses of the particle.
 */
void sort_by_d_eq(std::vector<double> &d_eq,
                  std::vector<double> &d_max,
                  std::vector<double> &m);

}

////////////////////////////////////////////////////////////////////////////////
// Scattering data for given temperature and frequency.
////////////////////////////////////////////////////////////////////////////////
// pxx :: export
/** ARTS SSDB scattering data for given particle temperature and frequency.
 *
 * The scattering data for a given particle, temperature and frequency is
 * contained in a single NetCDF group. This class provides an interface to these
 * groups and provides access to the scattering data.
 */
class ScatteringData {
  // pxx :: hide
  // Determines format of scattering data. Called by constructor.
  void determine_format();
  // pxx :: hide
  /** Extract vector from NetCDF data.
   * @param name Name of the variable containing the vector.
   * @return math::Vector containing the data.
   */
  template <typename Float>
  math::Vector<Float> get_vector(std::string name);

 public:
  /** Create ScatteringData object from NetCDF group.
   *
   * Extracts temperature and frequency of scattering data and
   * determines the data format.
   */
  ScatteringData(netcdf4::Group group) : group_(group) {
    temperature_ = group.get_variable("temperature").read<double>();
    frequency_ = group.get_variable("frequency").read<double>();
    determine_format();
  }

  //
  // Particle properties
  //

  /// The particle type described by the data.
  ParticleType get_particle_type();
  /// The format the scattering data is in.
  DataFormat get_format() const { return format_; }
  /// The frequency in GHz for which this data is valid.
  double get_frequency() { return frequency_; }
  /// The temperature in K for which the data is valid.
  double get_temperature() { return temperature_; }
  /// The l_max value used in the SHT transformation of spectral data.
  Index get_l_max();
  /// Size of the lon. grid of incoming angles.
  Index get_n_lon_inc();
  /// Size of the lat. grid of incoming angles.
  Index get_n_lat_inc();
  /// Size of the lon. grid of scattering angles.
  Index get_n_lon_scat();
  /// Size of the lat. grid of scattering angles.
  Index get_n_lat_scat();

  /** Compatible SHT object.
   * @return A spherical harmonics transform object that can be used to
   * transform data in spectral format.
   */
  sht::SHT get_sht();
  /// Length-1 vector containing the frequency corresponding to the scattering data.
  math::Vector<double> get_f_grid() {return math::Vector<double>::Constant(1, frequency_);}
  /// Length-1 vector containing the temperature corresponding to the scattering data.
  math::Vector<double> get_t_grid() {return math::Vector<double>::Constant(1, temperature_);}
  /// The lon. component of the incoming angle grid of the scattering data.
  math::Vector<double> get_lon_inc();
  /// The lat. component of the incoming angle grid of the scattering data.
  math::Vector<double> get_lat_inc();
  /// The lon. component of the scattering angle grid of the scattering data.
  math::Vector<double> get_lon_scat() { return get_vector<double>("aa_scat"); }
  /// The lat. component of the scattering angle grid of the scattering data.
  math::Vector<double> get_lat_scat() { return get_vector<double>("za_scat"); }

  //
  // Extracting data
  //

  /** The phase matrix data in gridded format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * @return Rank-7 tensor containing the phase matrix in gridded format. Axes
   * correspond to frequency (1), temperature (2), incoming angle lon. (3) and
   * lat. (4),  scattering angle lon. (5) and lat. (6) and the stokes components
   * (7).
   */
  math::Tensor<double, 7> get_phase_matrix_data_gridded();

  /** The phase matrix data in spectral format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * @return Rank-6 complex tensor containing the phase matrix in gridded format.
   * Axes correspond to frequency (1), temperature (2), incoming angle lon. (3)
   * and lat. (4),  SHT components (5) nd the stokes components (6).
   */
  math::Tensor<std::complex<double>, 6> get_phase_matrix_data_spectral();

  /** The extinction matrix data in gridded format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * Note that axes corresponding to scattering angles are degenerate due to
   * the extinction matrix not exhibiting any scattering-angle dependency.
   *
   * @return Rank-7 tensor containing the phase matrix in gridded format. Axes
   * correspond to frequency (1), temperature (2), incoming angle lon. (3) and
   * lat. (4),  scattering angle lon. (5) and lat. (6) and the stokes components
   * (7).
   */
  math::Tensor<double, 7> get_extinction_matrix_data_gridded();

  /** The extinction matrix data in spectral format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * Note that axes corresponding to scattering angles are degenerate due to
   * the extinction matrix not exhibiting any scattering-angle dependency.
   *
   * @return Rank-6 complex tensor containing the phase matrix in gridded format.
   * Axes correspond to frequency (1), temperature (2), incoming angle lon. (3)
   * and lat. (4),  SHT components (5) nd the stokes components (6).
   */
  math::Tensor<std::complex<double>, 6> get_extinction_matrix_data_spectral();

  /** The absorption vector in gridded format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * Note that axes corresponding to scattering angles are degenerate due to
   * the absorption vector not exhibiting any scattering-angle dependency.
   *
   * @return Rank-7 tensor containing the phase matrix in gridded format. Axes
   * correspond to frequency (1), temperature (2), incoming angle lon. (3) and
   * lat. (4),  scattering angle lon. (5) and lat. (6) and the stokes components
   * (7).
   */
  math::Tensor<double, 7> get_absorption_vector_data_gridded();

  /** The absorption vector data in spectral format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * Note that axes corresponding to scattering angles are degenerate due to
   * the absorption vector not exhibiting any scattering-angle dependency.
   *
   * @return Rank-6 complex tensor containing the phase matrix in gridded format.
   * Axes correspond to frequency (1), temperature (2), incoming angle lon. (3)
   * and lat. (4),  SHT components (5) nd the stokes components (6).
   */
  math::Tensor<std::complex<double>, 6> get_absorption_vector_data_spectral();

  /** The backward scattering coefficient in gridded format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * Note that axes corresponding to scattering angles are degenerate due to
   * the backward scattering coefficient  not exhibiting any scattering-angle dependency.
   *
   * @return Rank-7 tensor containing the phase matrix in gridded format. Axes
   * correspond to frequency (1), temperature (2), incoming angle lon. (3) and
   * lat. (4),  scattering angle lon. (5) and lat. (6) and the stokes components
   * (7).
   */
  math::Tensor<double, 7> get_backward_scattering_coeff_data_gridded();

  /** The backward scattering coefficient data in spectral format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * Note that axes corresponding to scattering angles are degenerate due to
   * the backward scattering coefficient not exhibiting any scattering-angle dependency.
   *
   * @return Rank-6 complex tensor containing the phase matrix in gridded format.
   * Axes correspond to frequency (1), temperature (2), incoming angle lon. (3)
   * and lat. (4),  SHT components (5) nd the stokes components (6).
   */
  math::Tensor<std::complex<double>, 6> get_backward_scattering_coeff_data_spectral();

  /** The forward scattering coefficient in gridded format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * Note that axes corresponding to scattering angles are degenerate due to
   * the forward scattering coefficient  not exhibiting any scattering-angle dependency.
   *
   * @return Rank-7 tensor containing the phase matrix in gridded format. Axes
   * correspond to frequency (1), temperature (2), incoming angle lon. (3) and
   * lat. (4),  scattering angle lon. (5) and lat. (6) and the stokes components
   * (7).
   */
  math::Tensor<double, 7> get_forward_scattering_coeff_data_gridded();

  /** The forward scattering coefficient data in spectral format.
   *
   * Note that the axes corresponding to temperature and frequency are degenerate
   * because the data is specific to one given frequency and temperature.
   *
   * Note that axes corresponding to scattering angles are degenerate due to
   * the forward scattering coefficient not exhibiting any scattering-angle dependency.
   *
   * @return Rank-6 complex tensor containing the phase matrix in gridded format.
   * Axes correspond to frequency (1), temperature (2), incoming angle lon. (3)
   * and lat. (4),  SHT components (5) nd the stokes components (6).
   */
  math::Tensor<std::complex<double>, 6> get_forward_scattering_coeff_data_spectral();

  /// Conversion to SingleScatteringData in gridded format.
  operator SingleScatteringDataGridded<double>();
  /// Conversion to SingleScatteringData in spectral format.
  operator SingleScatteringDataSpectral<double>();
  /// Conversion to SingleScatteringData in native data format.
  operator SingleScatteringData();

 private:
  DataFormat format_;
  double temperature_, frequency_;
  netcdf4::Group group_;
  };

////////////////////////////////////////////////////////////////////////////////
// ParticleFile
////////////////////////////////////////////////////////////////////////////////
// pxx :: export
/** ARTS SSDB data for a specific particle.
 *
 * The scattering data for a particle a specific particle of given size is
 * contained in a single NetCDF file. This class provides an interface to these
 * files and access to the groups, which contain the scattering data for given
 * temperatures and frequencies.
 */
class ParticleFile {


  // pxx :: hide
  /// Parses temperatures and frequencies of the data in the NetCDF file.
  void parse_temps_and_freqs();
  // pxx :: hide
  /// Determine angular grids with the highest resolution.
  std::array<math::Vector<double>, 4> get_angular_grids_gridded();
  // pxx :: hide
  /// Determine angular grids with the highest resolution.
  std::tuple<math::Vector<double>, math::Vector<double>, Index> get_angular_grids_spectral();

 public:

  /// Iterator class providing access to scattering data for specific
  /// frequencies and temperatures.
  class DataIterator;

  /** Create particle file object from ARTS SSBD NetCDF file.
   * @param filename String containing the path to the file to read.
   */
  ParticleFile(std::string filename);

  /// The type of the particle, i.e. random or azimuthally-random orientation.
  ParticleType get_particle_type();

  /** The habit name.
   *
   * In the ARTS SSDB particle habits have names and IDs. If the name
   * has been successfully extracted from the file path this function
   * will return the matched name. Otherwise the empty string is returned.
   */
  std::string get_habit_name() { return habit_name_; }
  /// The volume equivalent diameter of the particle in microns.
  double get_d_eq() {return d_eq_;}
  /// The maximum diameter of the particle in microns.
  double get_d_max() {return d_max_;}
  /// The mass of the particle in kilo grams.
  double get_mass() {return mass_;}
  /// The frequencies at which data is available.
  const std::vector<double>& get_frequencies() {return freqs_;}
  /// The frequencies at which data is available as Eigen vector.
  math::Vector<double> get_f_grid() {
      return math::VectorMap<double>(freqs_.data(), freqs_.size());
  }
  /// The temperatures at which data is available.
  std::vector<double> get_temperatures() {return temps_;}
  /// The temperatures at which data is available as Eigen vector.
  math::Vector<double> get_t_grid() {
      return math::VectorMap<double>(temps_.data(), temps_.size());
  }

  /// Iterator pointing to first frequency-temperature pair for which
  /// for which data is available.
  DataIterator begin();
  /// Iterator pointing to the end of the data in the file.
  DataIterator end();

  /** Return scattering data for given frequency and temperature indices.
   * @param f_index The index of the frequency for which to return data.
   * @param t_index The index of the temperature for which to return data.
   */
  ScatteringData get_scattering_data(size_t f_index, size_t t_index);

  /** Convert data to SingleScatteringData.
   *
   * Note: This will interpolate all data in the file to the angular grids of
   * the data of the particle at lowest frequency and temperature.
   */
  operator SingleScatteringData();
  /** Convert data to SingleScatteringData.
   *
   * Note: This will interpolate all data in the file to the angular grids of
   * the data of the particle at lowest frequency and temperature.
   */
  SingleScatteringData to_single_scattering_data() { return *this; }
  /** Convert data to Particle.
   *
   * This function extracts the scattering data as SingleScatteringData object
   * and sets the source of the particle data to ARTS SSDB.
   */
  Particle to_particle();

 private:
  std::string habit_name_;
  double d_eq_, d_max_, mass_;
  std::vector<double> freqs_;
  std::vector<double> temps_;
  std::map<std::pair<double, double>, netcdf4::Group> group_map_;
  netcdf4::File file_handle_;
};

/// Iterator class for scattering data in ARTS SSDB file.
class ParticleFile::DataIterator {
 public:
  DataIterator(const ParticleFile *file,
               size_t f_index = 0,
               size_t t_index = 0);
  bool operator==(const DataIterator &other) const
  { return (t_index_ == other.t_index_) && (f_index_ == other.f_index_); }
  bool operator!=(const DataIterator &other) const { return !(*this == other); }
  DataIterator& operator++();
  ScatteringData operator*();

  /// The frequency corresponding to the data pointed to by the iterator.
  double get_frequency() { return file_->freqs_[f_index_]; }
  /// The temperature corresponding to the data pointed to by the iterator.
  double get_temperature() { return file_->temps_[t_index_]; }

  // iterator traits
  using difference_type = size_t;
  using value_type = ScatteringData;
  using pointer = const ScatteringData *;
  using reference = const ScatteringData &;
  using iterator_category = std::forward_iterator_tag;

 private:
  const ParticleFile *file_ = nullptr;
  size_t f_index_, t_index_;
};

////////////////////////////////////////////////////////////////////////////////
// Habit folder
////////////////////////////////////////////////////////////////////////////////
// pxx :: export
/** A folder describing a particle habit.
 *
 * Particle habits in the ARTS SSDB are represented by a folder containing NetCDF4
 * files for each available particle size. This class parses such a folder and
 * provides access to each particle of the habit.
 */
class HabitFolder {

  // pxx :: hide
  /// Parse files in folder.
  void parse_files();

public:

  /// Iterator over particles in folder.
  class DataIterator;

  /** Create particle HabitFolder object providing access to scattering data.
   * @param path String containing the path to the data to load.
   */
  HabitFolder(std::string path) : base_path_(path) { parse_files(); }

  /** The habit name.
   *
   * If the name * has been successfully extracted from the file path
   * this function will return the matched name. Otherwise the empty string is returned.
   */
  std::string get_habit_name() { return habit_name_; }
  /// Number of particle in habit.
  size_t get_n_particles() {return d_eq_.size();}
  /// Return vector of equivalent diameters of the particles in the habit.
  math::Vector<double> get_d_eq() {return d_eq_;}
  /// Return vector of maximum diameters of the particles in the habit.
  math::Vector<double> get_d_max() {return d_max_;}
  /// Return vector of masses of the particles in the habit.
  math::Vector<double> get_mass() {return mass_;}

  /// Iterator pointing to the particle with the smallest mass.
  DataIterator begin();
  /// Iterator pointing to the end of the data.
  DataIterator end();

  /// Convert to ParticleHabit object containing the scattering data.
  operator ParticleHabit();
  /// Convert to ParticleHabit object containing the scattering data.
  ParticleHabit to_particle_habit() { return *this; }

 private:
  std::string habit_name_;
  std::filesystem::path base_path_;
  math::Vector<double> d_eq_;
  math::Vector<double> d_max_;
  math::Vector<double> mass_;
  std::map<double, std::string> files_;
};

class HabitFolder::DataIterator {
public:
DataIterator(const HabitFolder *folder, size_t index = 0)
    : folder_(folder),
      index_(index) {}

    DataIterator& operator++() {index_++; return *this;}
    bool operator==(const DataIterator &other) const {return (index_ == other.index_);}
    bool operator!=(const DataIterator &other) const {return !(*this == other);}
    ParticleFile operator*() {
        double d_eq = folder_->d_eq_[index_];
        return ParticleFile(folder_->files_.find(d_eq)->second);
    }
    /// Return volume-equivalent diameter of particle the iterator points to.
    double get_d_eq() {return folder_->d_eq_[index_];}
    /// Return maximum diamater of the particle the iterator points to.
    double get_d_max() {return folder_->d_max_[index_];}
    /// Return mass of the particle the iterator points to.
    double get_mass() {return folder_->mass_[index_];}

    // iterator traits
    using difference_type = size_t;
    using value_type = ScatteringData;
    using pointer = const ScatteringData*;
    using reference = const ScatteringData&;
    using iterator_category = std::forward_iterator_tag;

private:
    const HabitFolder *folder_ = nullptr;
    size_t index_;
};

}  // namespace arts_ssdb
}
