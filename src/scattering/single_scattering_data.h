/** \file scattering_data.h
 *
 * Defines the SingleScatteringData class which provides a format-angostic
 * interface for scattering data.
 *
 * The SingleScatteringData class is essentially a container for the abstract
 * class SingleScatteringDataImpl, which describes the generic interface for
 * scattering data classes. The specific implementations for gridded and
 * spectral scattering data classes are also provided in this file.
 *
 * @author Simon Pfreundschuh, 2020
 */
#ifndef __ARTS_SCATTERING_SINGLE_SCATTERING_DATA__
#define __ARTS_SCATTERING_SINGLE_SCATTERING_DATA__

#include <scattering/single_scattering_data_impl.h>
#include <cassert>
#include <memory>

namespace scattering {

////////////////////////////////////////////////////////////////////////////////
// Format-agnostic single scattering data class
////////////////////////////////////////////////////////////////////////////////

// pxx :: export
/** Format-agnostic interface and container for single-scattering data.
 *
 * This class provides a format-agnostic interface and container for single
 * scattering data. This means that it can store scattering data in any format
 * and allows manipulating and combining this data with any other single
 * scattering data in any other format.
 *
 * Note that SingleScatteringData objects are lightweight object, i.e. when
 * copy-constructed the resulting object will share the scattering data with
 * the original object.
 */
class SingleScatteringData {
 public:

  //
  // Constructors
  //

  // pxx :: hide
  /** Create from existing pointer to implementation object.
   * @param data Pointer to existing format-specific scattering data object.
   */
  SingleScatteringData(SingleScatteringDataImpl *data) : data_(data) {}

  SingleScatteringData() {}

  /// The copy constructor creates a shallow copy of the given single scattering data
  /// object.
  SingleScatteringData(const SingleScatteringData &other) = default;

  // pxx :: hide
  /** Create scattering data obejct from gridded scattering data.
   *
   * This constructor creates a single scattering data object, whose data is
   * stored in gridded format.
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The incoming-angle longitude grid.
   * @param lat_inc The incoming-angle latitude grid.
   * @param lon_scat The scattering-angle longitude grid.
   * @param lat_scat The scattering-angle latitude grid.
   * @param phase_matrix Tensor containing the phase matrix data.
   * @param extinction_matrix Tensor containing the extinction matrix data.
   * @param absorption_vector Tensor containing the absorption vector data.
   * @param backward_scattering_coeff Tensor containing the backward_scattering
   * coefficients
   * @param forward_scattering_coeff Tensor containing the forward_scattering
   * coefficients.
   */
  SingleScatteringData(
      scattering::math::VectorPtr<double> f_grid,
      scattering::math::VectorPtr<double> t_grid,
      scattering::math::VectorPtr<double> lon_inc,
      scattering::math::VectorPtr<double> lat_inc,
      scattering::math::VectorPtr<double> lon_scat,
      std::shared_ptr<LatitudeGrid<double>> lat_scat,
      scattering::math::TensorPtr<double, 7> phase_matrix,
      scattering::math::TensorPtr<double, 7> extinction_matrix,
      scattering::math::TensorPtr<double, 7> absorption_vector,
      scattering::math::TensorPtr<double, 7> backward_scattering_coeff,
      scattering::math::TensorPtr<double, 7> forward_scattering_coeff);

  /** Create SingleScatteringData object from gridded scattering data.
   *
   * This constructor creates a single scattering data object, whose data is
   * stored in gridded format.
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The incoming-angle longitude grid.
   * @param lat_inc The incoming-angle latitude grid.
   * @param lon_scat The scattering-angle longitude grid.
   * @param lat_scat The scattering-angle latitude grid.
   * @param phase_matrix Tensor containing the phase matrix data.
   * @param extinction_matrix Tensor containing the extinction matrix data.
   * @param absorption_vector Tensor containing the absorption vector data.
   * @param backward_scattering_coeff Tensor containing the backward_scattering
   * coefficients
   * @param forward_scattering_coeff Tensor containing the forward_scattering
   * coefficients.
   */
  SingleScatteringData(
      scattering::math::Vector<double> f_grid,
      scattering::math::Vector<double> t_grid,
      scattering::math::Vector<double> lon_inc,
      scattering::math::Vector<double> lat_inc,
      scattering::math::Vector<double> lon_scat,
      scattering::math::Vector<double> lat_scat,
      scattering::math::Tensor<double, 7> phase_matrix,
      scattering::math::Tensor<double, 7> extinction_matrix,
      scattering::math::Tensor<double, 7> absorption_vector,
      scattering::math::Tensor<double, 7> backward_scattering_coeff,
      scattering::math::Tensor<double, 7> forward_scattering_coeff)
      : SingleScatteringData(
            std::make_shared<math::Vector<double>>(f_grid),
            std::make_shared<math::Vector<double>>(t_grid),
            std::make_shared<math::Vector<double>>(lon_inc),
            std::make_shared<math::Vector<double>>(lat_inc),
            std::make_shared<math::Vector<double>>(lon_scat),
            std::make_shared<IrregularLatitudeGrid<double>>(lat_scat),
            std::make_shared<math::Tensor<double, 7>>(phase_matrix),
            std::make_shared<math::Tensor<double, 7>>(extinction_matrix),
            std::make_shared<math::Tensor<double, 7>>(absorption_vector),
            std::make_shared<math::Tensor<double, 7>>(
                backward_scattering_coeff),
            std::make_shared<math::Tensor<double, 7>>(
                forward_scattering_coeff)) {}

  /** Create empty SingleScatteringData object.
   *
   * Creates an empty container for single scattering data to be stored in
   * in gridded format.
   *
   * @param f_grid The frequency grid of the scattering data.
   * @param t_grid The temperature grid of the scattering data.
   * @param lon_inc The incoming-angle longitude grid.
   * @param lat_inc The incoming-angle latitude grid.
   * @param lon_scat The scattering-angle longitude grid.
   * @param lat_scat The scattering-angle latitude grid.
   * @param type The type of particle the data corresponds to (random,
   *  azimuthally random, or general).
   */
  SingleScatteringData(scattering::math::Vector<double> f_grid,
                       scattering::math::Vector<double> t_grid,
                       scattering::math::Vector<double> lon_inc,
                       scattering::math::Vector<double> lat_inc,
                       scattering::math::Vector<double> lon_scat,
                       scattering::math::Vector<double> lat_scat,
                       ParticleType type);

  // pxx :: hide
  /** Create from spectral scattering data.
   *
   * This constructor creates a SingleScatteringData object from existing
   * spectral scattering data.
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The incoming-angle longitude grid.
   * @param lat_inc The incoming-angle latitude grid.
   * @param lon_scat The scattering-angle longitude grid.
   * @param lat_scat The scattering-angle latitude grid.
   * @param phase_matrix Tensor containing the phase matrix data.
   * @param extinction_matrix Tensor containing the extinction matrix data.
   * @param absorption_vector Tensor containing the absorption vector data.
   * @param backward_scattering_coeff Tensor containing the backward_scattering
   * coefficients
   * @param forward_scattering_coeff Tensor containing the forward_scattering
   * coefficients.
   */
  SingleScatteringData(
      scattering::math::VectorPtr<double> f_grid,
      scattering::math::VectorPtr<double> t_grid,
      scattering::math::VectorPtr<double> lon_inc,
      scattering::math::VectorPtr<double> lat_inc,
      std::shared_ptr<sht::SHT> sht_scat,
      scattering::math::TensorPtr<std::complex<double>, 6> phase_matrix,
      scattering::math::TensorPtr<std::complex<double>, 6> extinction_matrix,
      scattering::math::TensorPtr<std::complex<double>, 6> absorption_vector,
      scattering::math::TensorPtr<std::complex<double>, 6>
          backward_scattering_coeff,
      scattering::math::TensorPtr<std::complex<double>, 6>
          forward_scattering_coeff);

  /** Create from spectral scattering data.
   *
   * This constructor creates a SingleScatteringData object from existing
   * spectral scattering data. Parameters
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The incoming-angle longitude grid.
   * @param lat_inc The incoming-angle latitude grid.
   * @param phase_matrix Tensor containing the phase matrix data.
   * @param extinction_matrix Tensor containing the extinction matrix data.
   * @param absorption_vector Tensor containing the absorption vector data.
   * @param backward_scattering_coeff Tensor containing the backward_scattering
   * coefficients
   * @param forward_scattering_coeff Tensor containing the forward_scattering
   * coefficients.
   */
  SingleScatteringData(
      scattering::math::Vector<double> f_grid,
      scattering::math::Vector<double> t_grid,
      scattering::math::Vector<double> lon_inc,
      scattering::math::Vector<double> lat_inc,
      sht::SHT sht_scat,
      scattering::math::Tensor<std::complex<double>, 6> phase_matrix,
      scattering::math::Tensor<std::complex<double>, 6> extinction_matrix,
      scattering::math::Tensor<std::complex<double>, 6> absorption_vector,
      scattering::math::Tensor<std::complex<double>, 6> backward_scattering_coeff,
      scattering::math::Tensor<std::complex<double>, 6> forward_scattering_coeff);

  /** Create from spectral scattering data.
   *
   * This constructor creates a SingleScatteringData object from existing
   * spectral scattering data. Parameters of the SHT transform are inferred
   * assuming that l_max = m_max.
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The incoming-angle longitude grid.
   * @param lat_inc The incoming-angle latitude grid.
   * @param phase_matrix Tensor containing the phase matrix data.
   * @param extinction_matrix Tensor containing the extinction matrix data.
   * @param absorption_vector Tensor containing the absorption vector data.
   * @param backward_scattering_coeff Tensor containing the backward_scattering
   * coefficients
   * @param forward_scattering_coeff Tensor containing the forward_scattering
   * coefficients.
   */
  SingleScatteringData(
      scattering::math::Vector<double> f_grid,
      scattering::math::Vector<double> t_grid,
      scattering::math::Vector<double> lon_inc,
      scattering::math::Vector<double> lat_inc,
      scattering::math::Tensor<std::complex<double>, 6> phase_matrix,
      scattering::math::Tensor<std::complex<double>, 6> extinction_matrix,
      scattering::math::Tensor<std::complex<double>, 6> absorption_vector,
      scattering::math::Tensor<std::complex<double>, 6> backward_scattering_coeff,
      scattering::math::Tensor<std::complex<double>, 6> forward_scattering_coeff);

  /** Create empty SingleScatteringData object in spectral format.
   *
   * This creates an empty container for single scattering data in spectral format
   * expecting data from an SHT transform with given l_max parameter. The m_max parameter
   * of the transform is assumed to be the same as the m_max parameter.
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The incoming-angle longitude grid.
   * @param lat_inc The incoming-angle latitude grid.
   * @param l_max The l_max parameter of the SHT transform defining how many p
   * @param lat_scat The scattering-angle latitude grid.
   * @param phase_matrix Tensor containing the phase matrix data.
   * @param extinction_matrix Tensor containing the extinction matrix data.
   * @param absorption_vector Tensor containing the absorption vector data.
   * @param backward_scattering_coeff Tensor containing the backward_scattering
   * coefficients
   * @param forward_scattering_coeff Tensor containing the forward_scattering
   * coefficients.
   */
  SingleScatteringData(scattering::math::Vector<double> f_grid,
                       scattering::math::Vector<double> t_grid,
                       scattering::math::Vector<double> lon_inc,
                       scattering::math::Vector<double> lat_inc,
                       Index l_max,
                       ParticleType type);

  /// Perform a deep copy of the scattering data object.
  SingleScatteringData copy() const {return SingleScatteringData(data_->copy());}

  //
  // Getters and setters.
  //

  /// The frequency grid.
  const math::Vector<double>& get_f_grid() const { return data_->get_f_grid(); }
  /// The temperature grid.
  const math::Vector<double>& get_t_grid() const { return data_->get_t_grid(); }
  /// The incoming-angle longitude grid.
  math::Vector<double> get_lon_inc() const { return data_->get_lon_inc(); }
  /// The incoming-angle latitude grid.
  math::Vector<double> get_lat_inc() const { return data_->get_lat_inc(); }
  /// The scattering-angle longitude grid.
  math::Vector<double> get_lon_scat() const { return data_->get_lon_scat(); }
  /// The scattering-angle latitude grid.
  math::Vector<double> get_lat_scat() const { return data_->get_lat_scat(); }

  /// Size of frequency grid.
  math::Index get_n_freqs() const { return data_->get_n_freqs(); }
  /// Size of temperature grid.
  math::Index get_n_temps() const { return data_->get_n_temps(); }
  /// Size of incoming-angle longitude grid.
  math::Index get_n_lon_inc() const { return data_->get_n_lon_inc(); }
  /// Size of incoming-angle latitude grid.
  math::Index get_n_lat_inc() const { return data_->get_n_lat_inc(); }
  /// Size of scattering-angle longitude grid.
  math::Index get_n_lon_scat() const { return data_->get_n_lon_scat(); }
  /// Size of scattering-angle latitude grid.
  math::Index get_n_lat_scat() const { return data_->get_n_lat_scat(); }
  /// L-max parameter of SHT transform used for scattering angles.
  math::Index get_l_max_scat() const { return data_->get_l_max_scat(); }
  /// M-max parameter of SHT transform used for scattering angles.
  math::Index get_m_max_scat() const { return data_->get_m_max_scat(); }
  /// Stokes dimensions of scattering data.
  math::Index get_stokes_dim() const { return data_->get_stokes_dim(); }

  /// The type of particle: random, azimuthally-random or general.
  ParticleType get_particle_type() const { return data_->get_particle_type(); }
  /// The data format: Gridded, spectral or fully spectral.
  DataFormat get_data_format() const { return data_->get_data_format(); }

  /** Copies data from existing single scattering data object into this object.
   *
   * @param f_index Index along the frequency grid for which to set the data.
   * @param t_index Index along the temperature grid for which to set the data.
   * @param source The object from which to copy the data. The data that is
   * copied corresponds to the data at frequency and temperature indices (0, 0).
   */
  void set_data(Index f_index,
                Index t_index,
                const SingleScatteringData &other) {
      data_->set_data(f_index, t_index, *other.data_);
  }

  //
  // Phase matrix data.
  //

  /** Extract phase function.
   *
   * Converts scattering data to gridded format (if necessary) and returns
   * rank-6 tensor containing only the first component of the scattering data.
   *
   * @return Rank-6 tensor containing the first coefficient of gridded
   * phase-matrix data.
   */
  math::Tensor<double, 6> get_phase_function() const {
      return data_->get_phase_function();
  }

  /** Calculate maxima of phase function along incoming angles.
   *
   * @return Rank-4 tensor containing the maximum values of the phase
   * function calculated across the incoming angles.
   */
  math::Tensor<double, 4> get_phase_function_max_inc() const {
      return data_->get_phase_function_max_inc();
  }

  /** Calculate maxima of phase function along scattering angles.
   *
   * @return Rank-4 tensor containing the maximum values of the phase
   * function calculated across the scattering angles.
   */
  math::Tensor<double, 4> get_phase_function_max_scat() const {
      return data_->get_phase_function_max_scat();
  }

  /** Extract phase function.
   *
   * Converts scattering data to spectral format and return rank-5 tensor
   * containing only the first component of the scattering data.
   *
   * @return Rank-5 tensor containing the first coefficient of spectral
   * phase-matrix data.
   */
  math::Tensor<std::complex<double>, 5> get_phase_function_spectral() const {
      return data_->get_phase_function_spectral();
  }
  /** Extract phase matrix data.
   *
   * Converts scattering data to gridded format (if necessary) and returns
   * rank-7 tensor containing the full scattering data.
   *
   * @return Rank-7 tensor containing the gridded phase-matrix data.
   */
  math::Tensor<double, 7> get_phase_matrix_data() const {
      return data_->get_phase_matrix_data();
  }
  /** Extract phase matrix data.
   *
   * Converts scattering data to spectral format (if necessary) and returns
   * rank-6 tensor containing the full scattering data.
   *
   * @return Rank-6 tensor containing the spectral phase-matrix data.
   */
  math::Tensor<std::complex<double>, 6> get_phase_matrix_data_spectral() const {
      return data_->get_phase_matrix_data_spectral();
  }

  /** Calculate phase matrix.
   *
   * Converts phase matrix data to gridded format (if necessary), assembles
   * the phase matrix data and returns a rank-8 tensor  containing the
   * phase matrix for all frequencies, temperatures and angles. The last two
   * dimensions correspond to the rows and columns of the scattering
   * matrix, respectively.
   *
   * @param stokes_dim The number of stokes dimensions of the scattering matrix
   * to return.
   * @return Rank-8 tensor containing the phase matrix.
   */
  math::Tensor<double, 8> get_phase_matrix(Index stokes_dim) const {
      return data_->get_phase_matrix(stokes_dim);
  }

  /** Interpolate and extract phase matrix data.
   *
   * Calculates the Stokes phase matrix for given frequency, temperature,
   * and incoming and outgoing directions.
   *
   * @param frequency The frequency for which to extract the extinction matrix.
   * @param temperature The temperature for which to extract the extinction
   * matrix.
   * @param lon_inc The longitude comp. of the incoming angle in radians.
   * @param lat_inc The latitude comp. of the incoming angle in radians.
   * @param lon_scat The longitude comp. of the scat. angle in radians.
   * @param lat_scat The latitude comp. of the scat. angle in radians.
   * @stokes_dim The stokes dimension for which to extract
   */
  math::Matrix<double> get_phase_matrix(double frequency,
                                         double temperature,
                                         double lon_inc,
                                         double lat_inc,
                                         double lon_scat,
                                         double lat_scat,
                                         Index stokes_dim) {
      return data_->get_phase_matrix(
          frequency, temperature, lon_inc, lat_inc, lon_scat, lat_scat, stokes_dim
      );
  }

  //
  // Extinction matrix data.
  //

  /** Extract extinction coefficient.
   *
   * Converts extinction matrix data to gridded format (if necessary) and
   * returns rank-6 tensor containing the first coefficient of the
   * extinction matrix data.
   *
   * @return Rank-6 tensor containing the first component of the extinction matrix
   * data.
   */
  math::Tensor<double, 6> get_extinction_coeff() const {
      return data_->get_extinction_coeff();
  }

  /** Extract extinction matrix data.
   *
   * Converts extinction matrix data to gridded format (if necessary) and
   * returns rank-7 tensor containing the full extinction matrix data.
   *
   * @return Rank-7 tensor containing the extinction matrix
   * data.
   */
  math::Tensor<double, 7> get_extinction_matrix_data() const {
      return data_->get_extinction_matrix_data();
  }

  /** Extract extinction matrix.
   *
   * Converts extinction matrix data to gridded format (if necessary) and
   * returns rank-8 tensor containing extinction matrix data expanded to
   * full extinction matrices.
   *
   * The last two dimensions correspond to the rows and columns of the
   * extinction matrix, respectively.
   *
   *
   * @return Rank-7 tensor containing the extinction matrix data in compact
   * format.
   */
  math::Tensor<double, 8> get_extinction_matrix(Index stokes_dim) const {
    return data_->get_extinction_matrix(stokes_dim);
  }

  /** Calculate extinction matrix.
   *
   * Calculate extinction matrix for given frequency temperature,
   * incoming and scattering angles.
   *
   * @param frequency The frequency for which to extract the extinction
   * matrix
   * @param temperature The temperature for which to extract the extinction
   * matrix.
   * @param lon_inc The longitude comp. of the incoming angle in radians.
   * @param lat_inc The latitude comp. of the incoming angle in radians.
   * @stokes_dim The stokes dimension of the output matrix.
   * @return The extinction matrix for the given frequency, temperature
   * and angles.
   */
  math::Matrix<double> get_extinction_matrix(
      double frequency,
      double temperature,
      double lon_inc,
      double lat_inc,
      Index stokes_dim) const {
      return data_->get_extinction_matrix(
          frequency, temperature, lon_inc, lat_inc, stokes_dim
          );
  }

  //
  // Absorption vector.
  //

  /** Extract absorption coefficient.
   *
   * Converts absorption vector data to gridded format (if necessary) and
   * returns rank-6 tensor containing the first coefficient of the
   * absorption vector data.
   *
   * @return Rank-6 tensor containing the first component of the absorption
   * vector data.
   */
  math::Tensor<double, 6> get_absorption_coeff() const {
      return data_->get_absorption_coeff();
  }

  /** Extract absorption vector data.
   *
   * Converts absorption vector data to gridded format (if necessary) and
   * returns rank-7 tensor containing the absorption vector in compact
   * format.
   *
   * @return Rank-7 tensor containing the absorption vector data.
   */
  math::Tensor<double, 7> get_absorption_vector_data() const {
      return data_->get_absorption_vector_data();
  }

  /** Extract and expand absorption vector.
   *
   * Converts absorption vector data to gridded format (if necessary) and
   * returns rank-7 tensor containing the absorption vector data expanded
   * to full stokes vector form.
   *
   * @param stokes_dim The stokes dimensions to which to expand the data.
   * @return Rank-7 tensor containing the absorption vector data.
   */
  math::Tensor<double, 7> get_absorption_vector(Index stokes_dim) const {
    return data_->get_absorption_vector(stokes_dim);
  }

  /** Calculate absorption vector.
   *
   * Interpolate and expands absorption vector data for given frequency,
   * temperature, incoming and scattering angles.
   *
   * @param temperature The temperature for which to calculate the absorption vector.
   * @param frequency The frequency for which to calculate the absorption vector.
   * @param lon_inc The incoming-angle longitude.
   * @param lat_inc The incoming-angle latitude.
   * @param stokes_dim The stokes dimensions to which to expand the data.
   */
  math::Vector<double> get_absorption_vector(
      double temperature,
      double frequency,
      double lon_inc,
      double lat_inc,
      Index stokes_dim) const {
      return data_->get_absorption_vector(
          temperature,
          frequency,
          lon_inc,
          lat_inc,
          stokes_dim);
  }

  /** Forward scattering coefficient.
   *
   * Converts forward scattering coefficient to gridded format (if necessary) and
   * return the forward scattering coefficient data.
   *
   * @return Rank-7 tensor containing the forward scattering coefficient data.
   */
  math::Tensor<double, 7> get_forward_scattering_coeff() const {
    return data_->get_forward_scattering_coeff();
  }

  /** Backward scattering coefficient.
   *
   * Converts backward scattering coefficient to gridded format (if necessary) and
   * returns the forward scattering coefficient data.
   *
   * @return Rank-7 tensor containing the backward scattering coefficient data.
   */
  math::Tensor<double, 7> get_backward_scattering_coeff() const {
    return data_->get_backward_scattering_coeff();
  }

  /** Set the number of scattering coefficients.
   *
   * Reduces or expands the dimension of the scattering data corresponding
   * to the different Stokes coefficients.
   *
   * @param n The new dimensions of the Stokes-coefficient dimension.
   */
  void set_number_of_scattering_coeffs(Index n) {
      data_->set_number_of_scattering_coeffs(n);
  }

  /** Set the stokes dimension of the scattering data.
   *
   * Reduces scattering data to the minimum amount required to represent data
   * data for the given amount of stokes dimensions.
   *
   * @param n The number Stokes dimension to be included in the data.
   */
  void set_stokes_dim(Index stokes_dim) { data_->set_stokes_dim(stokes_dim); }

  //
  // Interpolation functions
  //

  /** Linear interpolation along frequency.
   * @param frequencies The frequency grid to which to interpolate the data.
   * @return New scattering data object interpolated to the given frequencies.
   */
  SingleScatteringData interpolate_frequency(
      math::Vector<double> frequencies) const {
    auto result = data_->interpolate_frequency(
        std::make_shared<math::Vector<double>>(frequencies));
    return SingleScatteringData(std::move(result));
  }
  // pxx :: hide
  SingleScatteringData interpolate_frequency(
      std::shared_ptr<math::Vector<double>> frequencies) const {
      auto result = data_->interpolate_frequency(frequencies);
      return SingleScatteringData(std::move(result));
  }

  /** Linear interpolation along temperature.
   * @param temperatures The temperature grid to which to interpolate the data.
   * @param extrapolate Whether or not to extrapolate the data at the boundaries.
   * @return New scattering data object interpolated to the given temperatures.
   */
  SingleScatteringData interpolate_temperature(
      math::Vector<double> temperatures,
      bool extrapolate=false) const {
    auto result = data_->interpolate_temperature(
        std::make_shared<math::Vector<double>>(temperatures),
        extrapolate);
    return SingleScatteringData(std::move(result));
  }

  // pxx :: hide
  SingleScatteringData interpolate_temperature(
      std::shared_ptr<math::Vector<double>> temperatures,
      bool extrapolate=false) const {
      auto result = data_->interpolate_temperature(temperatures, extrapolate);
      return SingleScatteringData(std::move(result));
  }

  /** Linear interpolation along angles.
   * @param lon_inc The incoming-angle longitudes to which to interpolate the data.
   * @param lat_inc The incoming-angle latitudes to which to interpolate the data.
   * @param lon_scat The scattering-angle longitudes to which to interpolate the data.
   * @param lat_scat The scattering-angle latitudes to which to interpolate the data.
   * @return New temperatures object interpolated to the given temperatures.
   */
  SingleScatteringData interpolate_angles(math::Vector<double> lon_inc,
                                          math::Vector<double> lat_inc,
                                          math::Vector<double> lon_scat,
                                          math::Vector<double> lat_scat) const {
    auto result = data_->interpolate_angles(
        std::make_shared<math::Vector<double>>(lon_inc),
        std::make_shared<math::Vector<double>>(lat_inc),
        std::make_shared<math::Vector<double>>(lon_scat),
        std::make_shared<IrregularLatitudeGrid<double>>(lat_scat));
    return SingleScatteringData(std::move(result));
  }

  // pxx :: hide
  SingleScatteringData interpolate_angles(std::shared_ptr<math::Vector<double>> lon_inc,
                                          std::shared_ptr<math::Vector<double>> lat_inc,
                                          std::shared_ptr<math::Vector<double>> lon_scat,
                                          std::shared_ptr<LatitudeGrid<double>> lat_scat) const {
      auto result = data_->interpolate_angles(lon_inc,
                                              lat_inc,
                                              lon_scat,
                                              lat_scat);
      return SingleScatteringData(std::move(result));
  }

  /** Downsample scattering angles by averaging.
   *
   * This function downsamples the angular resolution of scattering data but keeps
   * the angular integral constant by averaging over the original data.
   *
   * @param lon_scat The scattering-angle longitude grid to downsample the data to.
   * @param lat_scat The scattering-angle latitude grid to downsample the data to.
   * @return New SingleScatteringData object with the phase matrix data downsampled
   * along scattering angles.
   */
  SingleScatteringData downsample_scattering_angles(math::Vector<double> lon_scat,
                                                    math::Vector<double> lat_scat) const {
    auto result = data_->downsample_scattering_angles(
        std::make_shared<math::Vector<double>>(lon_scat),
        std::make_shared<IrregularLatitudeGrid<double>>(lat_scat)
        );
    return SingleScatteringData(std::move(result));
  }

  // pxx :: hide
  SingleScatteringData downsample_scattering_angles(
      std::shared_ptr<math::Vector<double>> lon_scat,
      std::shared_ptr<LatitudeGrid<double>> lat_scat) const {
    auto result =
        data_->downsample_scattering_angles(lon_scat, lat_scat);
    return SingleScatteringData(std::move(result));
  }

  // pxx :: hide
  SingleScatteringData downsample_lon_scat(
      std::shared_ptr<math::Vector<double>> lon_scat) const {
      auto result =
          data_->downsample_lon_scat(lon_scat);
      return SingleScatteringData(std::move(result));
  }

  /** Downsample scattering longitude angles by averaging.
   *
   * This function downsamples the angular resolution of longitude components
   * of the scattering angle but keeps the angular integral constant by averaging
   * over the original data.
   *
   * @param lon_scat The scattering-angle longitude grid to downsample the data to.
   * @return New SingleScatteringData object with the phase matrix data downsampled
   * along the longitude component of the scattering angles.
   */
  SingleScatteringData downsample_scattering_angles(math::Vector<double> lon_scat) const {
      auto result = data_->downsample_lon_scat(
          std::make_shared<math::Vector<double>>(lon_scat));
      return SingleScatteringData(std::move(result));
  }

  /// Regrid scattering data to shnts conform grids.
  SingleScatteringData regrid() const { return data_->regrid(); }

  //
  // Arithmetic operations.
  //

  /** Accumulate scattering data into this object.
   *
   * Converts (if necessary) scattering data in other into format of this data
   * and accumulates the data of other into the data of this object. Note
   * that this may require regridding for the gridded format.
   *
   * @param other The scattering data to accumulate into this object.
   */
  SingleScatteringData &operator+=(const SingleScatteringData &other) {
    data_->operator+=(other.data_.get());
    return *this;
  }
  /** Sum scattering data.
   *
   * Creates a deep copy of this object and accumulates the data of other
   * into it.
   *
   * @param other The right-hand operand of the sum.
   * @return A new and independent scattering data object containing the sum of
   * the data in this and other.
   */
  SingleScatteringData operator+(const SingleScatteringData &other) {
    return SingleScatteringData(data_->operator+(other.data_.get()));
  }

  // Scaling
  /** Scale scattering data.
   * @param c The scaling factor.
   */
  SingleScatteringData &operator*=(double c) {
    data_->operator*=(c);
    return *this;
  }
  /** Copy and scale scattering data.
   * Performs a deep copy of this object and scales the result by the given scaling
   * factor.
   * @param c The scaling factor.
   */
  SingleScatteringData operator*(double c) const {
    return SingleScatteringData(data_->operator*(c));
  }

  /** Normalize phase matrix data.
   * Normalizes the phase matrix data to integrate to the given value.
   * @param norm The value to which the phase matrix data should integrate.
   */
  void normalize(double norm) { data_->normalize(norm); }

  // Conversion
  inline SingleScatteringData to_gridded() const;
  inline SingleScatteringData to_spectral() const;
  inline SingleScatteringData to_spectral(Index l_max, Index m_max) const;
  inline SingleScatteringData to_spectral(Index l_max,
                                          Index m_max,
                                          Index n_lon,
                                          Index n_lat) const;

  // pxx :: hide
  SingleScatteringData to_lab_frame(std::shared_ptr<math::Vector<double>> lat_inc,
                                    std::shared_ptr<math::Vector<double>> lon_scat,
                                    std::shared_ptr<LatitudeGrid<double>> lat_scat,
                                    Index stokes_dim) const {
      return data_->to_lab_frame(lat_inc, lon_scat, lat_scat, stokes_dim);
  }

  SingleScatteringData to_lab_frame(Index n_lat_inc,
                                    Index n_lon_scat,
                                    Index stokes_dim) const {
    return data_->to_lab_frame(n_lat_inc, n_lon_scat, stokes_dim);
  }

 private:
  std::shared_ptr<SingleScatteringDataImpl> data_;
};

////////////////////////////////////////////////////////////////////////////////
// SingleScatteringData
////////////////////////////////////////////////////////////////////////////////

SingleScatteringData SingleScatteringData::to_gridded() const {
  return SingleScatteringData(
      new SingleScatteringDataGridded<double>(*data_->to_gridded()));
}

SingleScatteringData SingleScatteringData::to_spectral() const {
  return SingleScatteringData(
      new SingleScatteringDataSpectral<double>(*data_->to_spectral()));
}

SingleScatteringData SingleScatteringData::to_spectral(Index l_max,
                                                       Index m_max) const {
  return SingleScatteringData(new SingleScatteringDataSpectral<double>(
      *data_->to_spectral(l_max, m_max)));
}

SingleScatteringData SingleScatteringData::to_spectral(Index l_max,
                                                       Index m_max,
                                                       Index n_lon,
                                                       Index n_lat) const {
    return SingleScatteringData(new SingleScatteringDataSpectral<double>(
                                    *data_->to_spectral(l_max, m_max, n_lon, n_lat)));
}

}  // namespace scattering

#endif
