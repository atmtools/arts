/** \file scattering_data_impl.h
 *
 * Contains the definition of the SingleScatteringDataImple class which
 * defines a format-agnostic interface for scattering data.
 *
 * Also provides the SingleScatteringDataGridded and
 * SingleScatteringDataSpectral class, which provide implementations of
 * this interface for gridded and spectral scattering data.
 *
 * @author Simon Pfreundschuh, 2020
 */
#pragma once

#include <scattering/maths.h>
#include <scattering/interpolation.h>
#include <scattering/scattering_data_field.h>
#include <scattering/sht.h>
#include <scattering/stokes.h>

#include <cassert>
#include <memory>

namespace scattering {

////////////////////////////////////////////////////////////////////////////////
// Forward declarations
////////////////////////////////////////////////////////////////////////////////

template <typename Scalar>
class SingleScatteringDataGridded;
template <typename Scalar>
class SingleScatteringDataSpectral;
template <typename Scalar>
class SingleScatteringDataFullySpectral;

////////////////////////////////////////////////////////////////////////////////
// Interface for format-specific implementations.
////////////////////////////////////////////////////////////////////////////////

namespace detail {

/// The number of elements required to represent the phase matrix.
inline Index get_n_phase_matrix_elements(ParticleType type) {
  switch (type) {
    case ParticleType::Random:
      return 6;
    case ParticleType::AzimuthallyRandom:
      return 16;
    case ParticleType::General:
      return 16;
  }
  return 1;
}

/// The number of elements required to represent the extinction matrix.
inline Index get_n_extinction_matrix_elements(ParticleType type) {
  switch (type) {
    case ParticleType::Random:
      return 1;
    case ParticleType::AzimuthallyRandom:
      return 3;
    case ParticleType::General:
      return 4;
  }
  return 1;
}

/// The number of elements required to represent the absorption vector.
inline Index get_n_absorption_vector_elements(ParticleType type) {
  switch (type) {
    case ParticleType::Random:
      return 1;
    case ParticleType::AzimuthallyRandom:
      return 2;
    case ParticleType::General:
      return 4;
  }
  return 1;
}

/** Conditional deleter for pointers returned from a conversion.
 *
 * This custom deleter is used to avoid unnecessary copy when trying
 * to convert an object to a data representation, that it is already
 * in.
 */
template <typename T>
struct ConditionalDeleter {
  void operator()(T *ptr) {
    if (ptr && do_delete) {
      delete ptr;
    }
  }

  bool do_delete;
};

template <typename T>
using ConversionPtr = std::unique_ptr<T, ConditionalDeleter<T>>;

template <typename T>
ConversionPtr<T> make_conversion_ptr(T *t, bool do_delete) {
  return ConversionPtr<T>(t, ConditionalDeleter<T>{do_delete});
}

}  // namespace detail

/** Generic interface for scattering data classes.
 *
 * This abstract class defines the generic interface for class representing
 * different scattering data formats.
 *
 * Scattering data object are lightweight object, meaning that their default
 * behavior is that data the underlying data is not copied when an object
 * is copied.
 */
class SingleScatteringDataImpl {
 public:
  virtual ~SingleScatteringDataImpl() = default;

  /// Size of the frequency grid.
  virtual math::Index get_n_freqs() const = 0;
  /// Size of the temperature grid.
  virtual math::Index get_n_temps() const = 0;
  /// Size of the incoming-angle longitude grid.
  virtual math::Index get_n_lon_inc() const = 0;
  /// Size of the incoming-angle latitude grid.
  virtual math::Index get_n_lat_inc() const = 0;
  /// Size of the scattering-angle longitude grid.
  virtual math::Index get_n_lon_scat() const = 0;
  /// Size of the scattering-angle latitude grid.
  virtual math::Index get_n_lat_scat() const = 0;
  /// L-max parameter of SHT transform for scattering_angles.
  virtual math::Index get_l_max_scat() const = 0;
  /// M-max parameter of SHT transform for scattering-angles.
  virtual math::Index get_m_max_scat() const = 0;
  /// Stokes dimension of scattering data.
  virtual math::Index get_stokes_dim() const = 0;

  /// The frequency grid.
  virtual const math::Vector<double>& get_f_grid() const = 0;
  /// The temperature grid.
  virtual const math::Vector<double>& get_t_grid() const = 0;
  /// The incoming-angle longitude grid.
  virtual math::Vector<double> get_lon_inc() const = 0;
  /// The incoming-angle latitude grid.
  virtual math::Vector<double> get_lat_inc() const = 0;
  /// The scattering-angle longitude grid.
  virtual math::Vector<double> get_lon_scat() const = 0;
  /// The scattering-angle latitude grid.
  virtual math::Vector<double> get_lat_scat() const = 0;

  /// The type of particle: random, azimuthally-random or general.
  virtual ParticleType get_particle_type() const = 0;
  /// The data format: Gridded, spectral or fully spectral.
  virtual DataFormat get_data_format() const = 0;

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
  virtual math::Tensor<double, 6> get_phase_function() const = 0;

  /** Calculate maxima of phase function along incoming angles.
   *
   * @return Rank-4 tensor containing the maximum values of the phase
   * function calculated across the incoming angles.
   */
  virtual math::Tensor<double, 4> get_phase_function_max_inc() const = 0;

  /** Calculate maxima of phase function along scattering angles.
   *
   * @return Rank-4 tensor containing the maximum values of the phase
   * function calculated across the scattering angles.
   */
  virtual math::Tensor<double, 4> get_phase_function_max_scat() const = 0;

  /** Extract phase function.
   *
   * Converts scattering data to spectral format (if necessary) and returns
   * rank-5 tensor containing only the first component of the scattering data.
   *
   * @return Rank-5 tensor containing the first coefficient of spectral
   * phase-matrix data.
   */
  virtual math::Tensor<std::complex<double>, 5> get_phase_function_spectral() const = 0;

  /** Extract phase matrix data.
   *
   * Converts scattering data to gridded format (if necessary) and returns
   * rank-7 tensor containing the full scattering data.
   *
   * @return Rank-7 tensor containing the gridded phase-matrix data.
   */
  virtual math::Tensor<double, 7> get_phase_matrix_data() const = 0;

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
  virtual math::Matrix<double> get_phase_matrix(
      double frequency,
      double temperature,
      double lon_inc,
      double lat_inc,
      double lon_scat,
      double lat_scat,
      Index stokes_dim
      ) const = 0;

  /** Extract phase matrix data.
   *
   * Converts scattering data to spectral format (if necessary) and returns
   * rank-6 tensor containing the full scattering data.
   *
   * @return Rank-6 tensor containing the spectral phase-matrix data.
   */
  virtual math::Tensor<std::complex<double>, 6> get_phase_matrix_data_spectral() const = 0;

  /** Extract phase matrix data.
   *
   * Converts phase matrix data to gridded format (if necessary), assembles
   * phase matrix data and returns rank-8 tensor containing the phase matrix.
   *
   * The last two dimensions correspond to the rows and columns of the
   * phase matrix, respectively.
   *
   * @param stokes_dim The number of stokes dimensions of the scattering matrix
   * to return.
   * @return Rank-8 tensor containing the Stokes phase matrix.
   */
  virtual math::Tensor<double, 8> get_phase_matrix(Index stokes_dim) const = 0;

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
  virtual math::Tensor<double, 6> get_extinction_coeff() const = 0;

  /** Extract extinction matrix data.
   *
   * Converts extinction matrix data to gridded format (if necessary) and
   * returns rank-7 tensor containing the full extinction matrix data.
   *
   * @return Rank-7 tensor containing the extinction matrix
   * data.
   */
  virtual math::Tensor<double, 7> get_extinction_matrix_data() const = 0;

  /** Extract extinction matrix.
   *
   * Converts extinction matrix data to gridded format (if necessary) and
   * returns rank-8 tensor containing extinction matrix data expanded to
   * full extinction matrices.
   *
   * The last two dimensions correspond to the rows and columns of the
   * scattering matrix, respectively.
   *
   *
   * @return Rank-7 tensor containing the extinction matrix data in compact
   * format.
   */
  virtual math::Tensor<double, 8> get_extinction_matrix(Index stokes_dim) const = 0;

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
   * @param lon_scat The longitude comp. of the scat. angle in radians.
   * @param lat_scat The latitude comp. of the scat. angle in radians.
   * @stokes_dim The stokes dimension of the output matrix.
   * @return The extinction matrix for the given frequency, temperature
   * and angles.
   */
  virtual math::Matrix<double> get_extinction_matrix(
      double frequency,
      double temperature,
      double lon_inc,
      double lat_inc,
      Index stokes_dim) const = 0;

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
  virtual math::Tensor<double, 6> get_absorption_coeff() const = 0;

  /** Extract absorption vector data.
   *
   * Converts absorption vector data to gridded format (if necessary) and
   * returns rank-7 tensor containing the absorption vector in compact
   * format.
   *
   * @return Rank-7 tensor containing the absorption vector data.
   */
  virtual math::Tensor<double, 7> get_absorption_vector_data() const = 0;

  /** Extract and expand absorption vector.
   *
   * Converts absorption vector data to gridded format (if necessary) and
   * returns rank-7 tensor containing the absorption vector data expanded
   * to full stokes vector form.
   *
   * @param stokes_dim The stokes dimensions to which to expand the data.
   * @return Rank-7 tensor containing the absorption vector data.
   */
  virtual math::Tensor<double, 7> get_absorption_vector(Index stokes_dim) const = 0;

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
  virtual math::Vector<double> get_absorption_vector(
      double temperature,
      double frequency,
      double lon_inc,
      double lat_inc,
      Index stokes_dim) const = 0;

  //
  // Forward and backward scattering coefficient.
  //

  /** Forward scattering coefficient.
   *
   * Converts forward scattering coefficient to gridded format (if necessary) and
   * return the forward scattering coefficient data.
   *
   * @return Rank-7 tensor containing the forward scattering coefficient data.
   */
  virtual math::Tensor<double, 7> get_forward_scattering_coeff() const = 0;

  /** Backward scattering coefficient.
   *
   * Converts backward scattering coefficient to gridded format (if necessary) and
   * returns the forward scattering coefficient data.
   *
   * @return Rank-7 tensor containing the backward scattering coefficient data.
   */
  virtual math::Tensor<double, 7> get_backward_scattering_coeff() const = 0;

  /** Deep copy of scattering data.
   *
   * Return new scattering data object with the data copied from this
   * object.
   */
  virtual SingleScatteringDataImpl *copy() const = 0;

  //
  // Setting of scattering data.
  //

  /** Copies data from existing single scattering data object into this object.
   *
   * @param f_index Index along the frequency grid for which to set the data.
   * @param t_index Index along the temperature grid for which to set the data.
   * @param source The object from which to copy the data. The data that is
   * copied corresponds to the data at frequency and temperature indices (0, 0).
   */
  virtual void set_data(Index f_index,
                        Index t_index,
                        const SingleScatteringDataImpl &source) = 0;

  /** Set the number of scattering coefficients.
   *
   * Reduces or expands the dimension of the scattering data corresponding
   * to the different Stokes coefficients.
   *
   * @param n The new dimensions of the Stokes-coefficient dimension.
   */
  virtual void set_number_of_scattering_coeffs(Index n) = 0;

  /** Set the stokes dimension of the scattering data.
   *
   * Reduces scattering data to the minimum amount required to represent data
   * data for the given amount of stokes dimensions.
   *
   * @param n The number Stokes dimension to be included in the data.
   */
  virtual void set_stokes_dim(Index stokes_dim) = 0;


  //
  // Interpolation.
  //

  /** Linear interpolation along frequency.
   * @param frequencies The frequency grid to which to interpolate the data.
   * @return New scattering data object interpolated to the given frequencies.
   */
  virtual SingleScatteringDataImpl *interpolate_frequency(
      math::VectorPtr<double> frequencies) const = 0;

  /** Linear interpolation along temperature.
   * @param temperatures The temperature grid to which to interpolate the data.
   * @param extrapolate Whether or not to extrapolate the data at the boundaries.
   * @return New scattering data object interpolated to the given temperatures.
   */
  virtual SingleScatteringDataImpl *interpolate_temperature(
      math::VectorPtr<double> temperatures,
      bool extrapolate=false) const = 0;
  /** Linear interpolation along angles.
   * @param lon_inc The incoming-angle longitudes to which to interpolate the data.
   * @param lat_inc The incoming-angle latitudes to which to interpolate the data.
   * @param lon_scat The scattering-angle longitudes to which to interpolate the data.
   * @param lat_scat The scattering-angle latitudes to which to interpolate the data.
   * @return New temperatures object interpolated to the given temperatures.
   */
  virtual SingleScatteringDataImpl *interpolate_angles(
      math::VectorPtr<double> lon_inc,
      math::VectorPtr<double> lat_inc,
      math::VectorPtr<double> lon_scat,
      std::shared_ptr<LatitudeGrid<double>>  lat_scat) const = 0;

  /** Downsample scattering angles by averaging.
   *
   * This function downsamples the angular resolution of scattering data but keeps
   * the angular integral constant by averaging over the original data.
   *
   * @param lon_scat The scattering-angle longitude grid to downsample the data to.
   * @param lat_scat The scattering-angle latitude grid to downsample the data to.
   */
  virtual SingleScatteringDataImpl *downsample_scattering_angles(
      math::VectorPtr<double> lon_scat,
      std::shared_ptr<LatitudeGrid<double>> lat_scat) const = 0;

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
  virtual SingleScatteringDataImpl *downsample_lon_scat(
      math::VectorPtr<double> lon_scat) const = 0;

  /// Regrid scattering data to shnts conform grids.
  virtual SingleScatteringDataImpl *regrid() const = 0;

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
  virtual void operator+=(const SingleScatteringDataImpl *other) = 0;
  /** Sum scattering data.
   *
   * Creates a deep copy of this object and accumulates the data of other
   * into it.
   *
   * @param other The right-hand operand of the sum.
   * @return A new and independent scattering data object containing the sum of
   * the data in this and other.
   */
  virtual SingleScatteringDataImpl *operator+(
      const SingleScatteringDataImpl *other) = 0;

  // Scaling
  /** Scale scattering data.
   * @param c The scaling factor.
   */
  virtual void operator*=(double c) = 0;
  /** Copy and scale scattering data.
   * Performs a deep copy of this object and scales the result by the given scaling
   * factor.
   * @param c The scaling factor.
   */
  virtual SingleScatteringDataImpl *operator*(double c) const = 0;
  /** Normalize phase matrix data.
   * Normalizes the phase matrix data to integrate to the given value.
   * @param norm The value to which the phase matrix data should integrate.
   */
  virtual void normalize(double norm) = 0;

  //
  // Format conversion
  //

  /// Convert to gridded format if necessary.
  virtual detail::ConversionPtr<const SingleScatteringDataGridded<double>>
  to_gridded() const = 0;

  /** Convert to spectral format.
   *
   * Converts scattering data to spectral format with the given parameters
   * of the SHT transform. If data is already in spectral format, the data
   * is only truncated/expanded to match the given parameters.
   *
   * @param l_max The l_max parameter of the SHT transformation.
   * @param m_max The m_max parameter of the SHT transformation.
   * @param n_lon The size of the longitude grid for the SHT transformation.
   * @param n_lat The size of the latitude grid for the SHT transformation.
   */
  virtual detail::ConversionPtr<const SingleScatteringDataSpectral<double>>
  to_spectral(Index l_max, Index m_max, Index n_lon, Index n_lat) const = 0;

  /** Convert to spectral format.
   *
   * Converts scattering data to spectral format with the given parameters
   * of the SHT transform. If data is already in spectral format, the data
   * is only truncated/expanded to match the given parameters.
   *
   * @param l_max The l_max parameter of the SHT transformation.
   * @param m_max The m_max parameter of the SHT transformation.
   * @param n_lon The size of the longitude grid for the SHT transformation.
   * @param n_lat The size of the latitude grid for the SHT transformation.
   */
  virtual detail::ConversionPtr<const SingleScatteringDataSpectral<double>>
      to_spectral(Index l_max, Index m_max) const = 0;

  /** Convert to spectral format.
   *
   * Convert to spectral format an infer transformation parameters from data.
   */
  virtual detail::ConversionPtr<const SingleScatteringDataSpectral<double>>
  to_spectral() const = 0;

  //
  // Conversion of reference frame.
  //

  /** Convert data from scattering frame to lab-frame.
   *
   * @param lat_inc The incoming-angle latitude grid.
   * @param lon_scat The scattering-angle longitude grid.
   * @param lat_scat The scattering-angle latitude grid.
   */
  virtual SingleScatteringDataImpl *to_lab_frame(std::shared_ptr<math::Vector<double>> lat_inc,
                                                 std::shared_ptr<math::Vector<double>> lon_scat,
                                                 std::shared_ptr<LatitudeGrid<double>> lat_scat,
                                                 Index stokes_dimensions) const = 0;

  /** Convert data from scattering frame to lab-frame.
   *
   * @param n_lat_inc The size of incoming-angle latitude grid.
   * @param n_lon_scat The size of the scattering-angle longitude grid.
   * @param stokes_dimensions The number of stokes to include in the transformed data.
   */
  virtual SingleScatteringDataImpl *to_lab_frame(Index n_lat_inc,
                                                 Index n_lon_scat,
                                                 Index stokes_dimensions) const = 0;

};


////////////////////////////////////////////////////////////////////////////////
// ScatteringDataBase
////////////////////////////////////////////////////////////////////////////////

template <typename Scalar>
class SingleScatteringDataBase {
 public:
  using VectorPtr = std::shared_ptr<math::Vector<Scalar>>;
  using ScatteringCoeffPtr = std::shared_ptr<math::Tensor<Scalar, 4>>;

  SingleScatteringDataBase(math::VectorPtr<Scalar> f_grid,
                           math::VectorPtr<Scalar> t_grid)
      : n_freqs_(f_grid->size()),
        n_temps_(t_grid->size()),
        f_grid_(f_grid),
        t_grid_(t_grid) {}

 protected:
  size_t n_freqs_, n_temps_;
  VectorPtr f_grid_;
  VectorPtr t_grid_;
};

////////////////////////////////////////////////////////////////////////////////
// ScatteringDataGridded
////////////////////////////////////////////////////////////////////////////////

/** Single scattering data given in gridded format.
 *
 * This class implements the single scattering data interface defined by
 * SingleScatteringDataImpl for single scattering data in gridded format.
 */
template <typename Scalar>
class SingleScatteringDataGridded : public SingleScatteringDataBase<Scalar>,
                                    public SingleScatteringDataImpl {

  using SingleScatteringDataBase<Scalar>::f_grid_;
  using SingleScatteringDataBase<Scalar>::n_freqs_;
  using SingleScatteringDataBase<Scalar>::n_temps_;
  using SingleScatteringDataBase<Scalar>::t_grid_;

  /** Check that dimensions of phase-matrix tensor are consistent with grids.
   * @throws runtime_error when phase matrix data is inconsistent with grids.
   */
  void check_consistency_phase_matrix() {
      auto data_tensor = phase_matrix_.get_data();
      if (data_tensor.dimension(0) != get_n_freqs()) {
          throw std::runtime_error(
              "Frequency grid of phase matrix is inconsistent with provided "
              "data tensor."
              );
      }
      if (data_tensor.dimension(1) != get_n_temps()) {
          throw std::runtime_error(
              "Temperature grid of phase matrix is inconsistent with provided "
              "data tensor."
              );
      }
      if (data_tensor.dimension(2) != get_n_lon_inc()) {
          throw std::runtime_error(
              "Incoming-angle longitude grid of phase matrix is inconsistent "
              "with provided data tensor."
              );
      }
      if (data_tensor.dimension(3) != get_n_lat_inc()) {
          throw std::runtime_error(
              "Incoming-angle latitude grid of phase matrix is inconsistent "
              "with provided data tensor."
              );
      }
      if (data_tensor.dimension(4) != get_n_lon_scat()) {
          throw std::runtime_error(
              "Scattering-angle longitude grid of phase matrix is inconsistent "
              "with provided data tensor."
              );
      }
      if (data_tensor.dimension(5) != get_n_lat_scat()) {
          throw std::runtime_error(
              "Scattering-angle latitude grid of phase matrix is inconsistent "
              "with provided data tensor."
              );
      }
  }

  /** Check that dimensions of extinction-matrix tensor are consistent with grids.
   * @throws runtime_error when extinction-matrix data is inconsistent with grids.
   */
  void check_consistency_extinction_matrix() {
      auto data_tensor = extinction_matrix_.get_data();
      if (data_tensor.dimension(0) != get_n_freqs()) {
          throw std::runtime_error(
              "Frequency grid of extinction matrix is inconsistent with provided "
              "data tensor."
              );
      }
      if (data_tensor.dimension(1) != get_n_temps()) {
          throw std::runtime_error(
              "Temperature grid of extinction matrix is inconsistent with "
              "provided data tensor."
              );
      }
      if ((data_tensor.dimension(2) != 1) && (data_tensor.dimension(2) != get_n_lon_inc())) {
          throw std::runtime_error(
              "Incoming-angle longitude grid of extinction matrix is "
              "inconsistent with provided data tensor."
              );
      }
      if ((data_tensor.dimension(2) != 1) && (data_tensor.dimension(3) != get_n_lat_inc())) {
          throw std::runtime_error(
              "Incoming-angle latitude grid of extinction matrix is "
              "inconsistent with provided data tensor."
              );
      }
      if (data_tensor.dimension(4) != 1) {
          throw std::runtime_error(
              "Expected scattering-angle longitude grid of absorption vector"
              "to be 1."
              );
      }
      if (data_tensor.dimension(5) != 1) {
          throw std::runtime_error(
              "Expected scattering-angle latitude grid of absorption vector"
              "to be 1."
              );
      }
  }

  /** Check that dimensions of absorption-vector tensor are consistent with grids.
   * @throws runtime_error when extinction-matrix data is inconsistent with grids.
   */
  void check_consistency_absorption_vector() {
      auto data_tensor = absorption_vector_.get_data();
      if (data_tensor.dimension(0) != get_n_freqs()) {
          throw std::runtime_error(
              "Frequency grid of absorption vector is inconsistent with provided "
              "data tensor."
              );
      }
      if (data_tensor.dimension(1) != get_n_temps()) {
          throw std::runtime_error(
              "Temperature grid of absorption vector is inconsistent with "
              "provided data tensor."
              );
      }
      if ((data_tensor.dimension(2) != 1) && (data_tensor.dimension(2) != get_n_lon_inc())) {
          throw std::runtime_error(
              "Incoming-angle longitude grid of absorption vector is "
              "inconsistent with provided data tensor."
              );
      }
      if ((data_tensor.dimension(2) != 1) && (data_tensor.dimension(3) != get_n_lat_inc())) {
          throw std::runtime_error(
              "Incoming-angle latitude grid of absorption vector is "
              "inconsistent with provided data tensor."
              );
      }
      if (data_tensor.dimension(4) != 1) {
          throw std::runtime_error(
              "Expected scattering-angle longitude grid of absorption vector"
              "to be 1."
              );
      }
      if (data_tensor.dimension(5) != 1) {
          throw std::runtime_error(
              "Expected scattering-angle latitude grid of absorption vector"
              "to be 1."
              );
      }
  }

 public:
  using Vector = math::Vector<Scalar>;
  using VectorPtr = std::shared_ptr<Vector>;
  using LatitudeGridPtr = std::shared_ptr<const LatitudeGrid<Scalar>>;
  using DataTensor = math::Tensor<Scalar, 7>;
  using DataPtr = std::shared_ptr<DataTensor>;
  using ScatteringCoeff = math::Tensor<Scalar, 4>;
  using ScatteringCoeffPtr = std::shared_ptr<ScatteringCoeff>;
  using OtherScalar =
      std::conditional<std::is_same<Scalar, double>::value, float, double>;

  //
  // Constructors
  //

  /** Create gridded scattering data from existing data.
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param phase_matrix The phase matrix data.
   * @param extinction_matrix The extinction matrix data.
   * @param absorption_vector The absorption vector data.
   * @param backward_scattering_coeff The backward scattering coefficient data.
   * @param forward_scattering_coeff The forward scattering coefficient data.
   */
  SingleScatteringDataGridded(
      VectorPtr f_grid,
      VectorPtr t_grid,
      ScatteringDataFieldGridded<Scalar> phase_matrix,
      ScatteringDataFieldGridded<Scalar> extinction_matrix,
      ScatteringDataFieldGridded<Scalar> absorption_vector,
      ScatteringDataFieldGridded<Scalar> backward_scattering_coeff,
      ScatteringDataFieldGridded<Scalar> forward_scattering_coeff)

      : SingleScatteringDataBase<Scalar>(f_grid, t_grid),
        phase_matrix_(phase_matrix),
        extinction_matrix_(extinction_matrix),
        absorption_vector_(absorption_vector),
        backward_scattering_coeff_(backward_scattering_coeff),
        forward_scattering_coeff_(forward_scattering_coeff) {
    check_consistency_phase_matrix();
    check_consistency_extinction_matrix();
    check_consistency_absorption_vector();
  }

  /** Create gridded scattering data from existing data.
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The incoming-angle latitude grid.
   * @param lat_inc The incoming-angle longitude grid.
   * @param lon_scat The scattering_angle longitude grid.
   * @param lat_scat The scattering-angle latitude grid.
   * @param phase_matrix Tensor containing the phase-matrix data.
   * @param extinction_matrix Tensor containing the extinction-matrix data.
   * @param absorption_vector Tensor containing the absorption-vector data.
   * @param backward_scattering_coeff Tensor containing the backscattering
   * coefficient data.
   * @param forward_scattering_coeff Tensor containing the forward scattering
   * coeff. data.
   */
  SingleScatteringDataGridded(VectorPtr f_grid,
                              VectorPtr t_grid,
                              VectorPtr lon_inc,
                              VectorPtr lat_inc,
                              VectorPtr lon_scat,
                              LatitudeGridPtr lat_scat,
                              DataPtr phase_matrix,
                              DataPtr extinction_matrix,
                              DataPtr absorption_vector,
                              DataPtr backward_scattering_coeff,
                              DataPtr forward_scattering_coeff)
      : SingleScatteringDataBase<Scalar>(f_grid, t_grid),
        phase_matrix_{f_grid,
                      t_grid,
                      lon_inc,
                      lat_inc,
                      lon_scat,
                      lat_scat,
                      phase_matrix},
        extinction_matrix_{f_grid,
                           t_grid,
                           lon_inc,
                           lat_inc,
                           dummy_grid_,
                           dummy_lat_grid_,
                           extinction_matrix},
        absorption_vector_{f_grid,
                           t_grid,
                           lon_inc,
                           lat_inc,
                           dummy_grid_,
                           dummy_lat_grid_,
                           absorption_vector},
        backward_scattering_coeff_{f_grid,
                                   t_grid,
                                   lon_inc,
                                   lat_inc,
                                   dummy_grid_,
                                   dummy_lat_grid_,
                                   backward_scattering_coeff},
        forward_scattering_coeff_(f_grid,
                                  t_grid,
                                  lon_inc,
                                  lat_inc,
                                  dummy_grid_,
                                  dummy_lat_grid_,
                                  forward_scattering_coeff) {
    check_consistency_phase_matrix();
    check_consistency_extinction_matrix();
    check_consistency_absorption_vector();
  }

  //
  // Getters and setters.
  //

  ParticleType get_particle_type() const { return phase_matrix_.get_particle_type(); }
  DataFormat get_data_format() const { return phase_matrix_.get_data_format(); }

  const math::Vector<double> &get_f_grid() const { return *f_grid_; }
  const math::Vector<double> &get_t_grid() const { return *t_grid_; }
  math::Vector<double> get_lon_inc() const { return phase_matrix_.get_lon_inc(); }
  math::Vector<double> get_lat_inc() const { return phase_matrix_.get_lat_inc(); }
  math::Vector<double> get_lon_scat() const { return phase_matrix_.get_lon_scat(); }
  math::Vector<double> get_lat_scat() const { return phase_matrix_.get_lat_scat(); }

  math::Index get_n_freqs() const { return phase_matrix_.get_n_freqs(); }
  math::Index get_n_temps() const { return phase_matrix_.get_n_temps(); }
  math::Index get_n_lon_inc() const { return phase_matrix_.get_n_lon_inc(); }
  math::Index get_n_lat_inc() const { return phase_matrix_.get_n_lat_inc(); }
  math::Index get_n_lon_scat() const { return phase_matrix_.get_n_lon_scat(); }
  math::Index get_l_max_scat() const { return phase_matrix_.get_sht_scat_params()[0]; }
  math::Index get_m_max_scat() const { return phase_matrix_.get_sht_scat_params()[1]; }
  math::Index get_n_lat_scat() const { return phase_matrix_.get_n_lat_scat(); }
  math::Index get_stokes_dim() const { return phase_matrix_.get_stokes_dim(); }

  void set_data(Index f_index,
                Index t_index,
                const SingleScatteringDataImpl &other) {
    auto converted = other.to_gridded();
    phase_matrix_.set_data(f_index, t_index, converted->phase_matrix_);
    extinction_matrix_.set_data(f_index,
                                t_index,
                                converted->extinction_matrix_);
    absorption_vector_.set_data(f_index,
                                t_index,
                                converted->absorption_vector_);
    forward_scattering_coeff_.set_data(f_index,
                                       t_index,
                                       converted->forward_scattering_coeff_);
    backward_scattering_coeff_.set_data(f_index,
                                        t_index,
                                        converted->backward_scattering_coeff_);
  }

  math::Tensor<Scalar, 6> get_phase_function() const {
    return phase_matrix_.get_phase_function();
  }

  math::Tensor<double, 4> get_phase_function_max_inc() const {
      return phase_matrix_.get_phase_function_max_inc();
  }

  math::Tensor<double, 4> get_phase_function_max_scat() const {
      return phase_matrix_.get_phase_function_max_scat();
  }

  math::Tensor<std::complex<Scalar>, 5> get_phase_function_spectral() const {
    stokes::PhaseMatrix<ScatteringDataFieldSpectral<Scalar>> spectral = phase_matrix_.to_spectral();
    return spectral.get_phase_function();
  }

  math::Tensor<Scalar, 8> get_phase_matrix(Index stokes_dim) const {
      return phase_matrix_.get_phase_matrix(stokes_dim);
  }

  math::Tensor<Scalar, 7> get_phase_matrix_data() const {
      return phase_matrix_.get_data();
  }

  math::Matrix<double> get_phase_matrix(
      double frequency,
      double temperature,
      double lon_inc,
      double lat_inc,
      double lon_scat,
      double lat_scat,
      Index stokes_dim
      ) const {
      return phase_matrix_.get_phase_matrix(
          frequency, temperature, lon_inc, lat_inc, lon_scat, lat_scat, stokes_dim
          );
  }

  math::Tensor<std::complex<Scalar>, 6> get_phase_matrix_data_spectral() const {
      stokes::PhaseMatrix<ScatteringDataFieldSpectral<Scalar>> spectral = phase_matrix_.to_spectral();
      return spectral.get_data();
  }

  math::Tensor<Scalar, 7> get_extinction_matrix_data() const {
      return extinction_matrix_.get_data();
  }

  math::Tensor<Scalar, 6> get_extinction_coeff() const {
      return extinction_matrix_.get_extinction_coeff();
  }

  math::Tensor<Scalar, 8> get_extinction_matrix(Index stokes_dim) const {
    return extinction_matrix_.get_extinction_matrix(stokes_dim);
  }

  math::Matrix<double> get_extinction_matrix(
      double frequency,
      double temperature,
      double lon_inc,
      double lat_inc,
      Index stokes_dim) const {
      return extinction_matrix_.get_extinction_matrix(
          frequency, temperature, lon_inc, lat_inc, stokes_dim
          );
  }

  math::Tensor<Scalar, 7> get_absorption_vector_data() const {
      return absorption_vector_.get_data();
  }

  math::Tensor<Scalar, 6> get_absorption_coeff() const {
      return absorption_vector_.get_absorption_coeff();
  }

  math::Tensor<Scalar, 7> get_absorption_vector(Index stokes_dim) const {
    return absorption_vector_.get_absorption_vector(stokes_dim);
  }

  math::Vector<double> get_absorption_vector(
      double temperature,
      double frequency,
      double lon_inc,
      double lat_inc,
      Index stokes_dim) const {
      return absorption_vector_.get_absorption_vector(
          temperature,
          frequency,
          lon_inc,
          lat_inc,
          stokes_dim);
  }

  math::Tensor<Scalar, 7> get_backward_scattering_coeff() const {
    return backward_scattering_coeff_.get_data();
  }

  math::Tensor<Scalar, 7> get_forward_scattering_coeff() const {
      return forward_scattering_coeff_.get_data();
  }

  SingleScatteringDataGridded *copy() const {
    return new SingleScatteringDataGridded(f_grid_,
                                           t_grid_,
                                           phase_matrix_.copy(),
                                           extinction_matrix_.copy(),
                                           absorption_vector_.copy(),
                                           backward_scattering_coeff_.copy(),
                                           forward_scattering_coeff_.copy());
  }

  SingleScatteringDataImpl *interpolate_frequency(
      math::VectorPtr<double> frequencies) const {
    auto phase_matrix = ScatteringDataFieldGridded<Scalar>(
        phase_matrix_.interpolate_frequency(frequencies));
    auto extinction_matrix = ScatteringDataFieldGridded<Scalar>(
        extinction_matrix_.interpolate_frequency(frequencies));
    auto absorption_vector = ScatteringDataFieldGridded<Scalar>(
        absorption_vector_.interpolate_frequency(frequencies));
    auto backward_scattering_coeff = ScatteringDataFieldGridded<Scalar>(
        backward_scattering_coeff_.interpolate_frequency(frequencies));
    auto forward_scattering_coeff = ScatteringDataFieldGridded<Scalar>(
        forward_scattering_coeff_.interpolate_frequency(frequencies));

    return new SingleScatteringDataGridded(frequencies,
                                           t_grid_,
                                           phase_matrix,
                                           extinction_matrix,
                                           absorption_vector,
                                           backward_scattering_coeff,
                                           forward_scattering_coeff);
  }

  SingleScatteringDataImpl *interpolate_temperature(
      math::VectorPtr<Scalar> temperatures,
      bool extrapolate=false) const {
    auto phase_matrix = ScatteringDataFieldGridded<Scalar>(
        phase_matrix_.interpolate_temperature(temperatures, extrapolate));
    auto extinction_matrix = ScatteringDataFieldGridded<Scalar>(
        extinction_matrix_.interpolate_temperature(temperatures, extrapolate));
    auto absorption_vector = ScatteringDataFieldGridded<Scalar>(
        absorption_vector_.interpolate_temperature(temperatures, extrapolate));
    auto backward_scattering_coeff = ScatteringDataFieldGridded<Scalar>(
        backward_scattering_coeff_.interpolate_temperature(temperatures, extrapolate));
    auto forward_scattering_coeff = ScatteringDataFieldGridded<Scalar>(
        forward_scattering_coeff_.interpolate_temperature(temperatures, extrapolate));

    return new SingleScatteringDataGridded(f_grid_,
                                           temperatures,
                                           phase_matrix,
                                           extinction_matrix,
                                           absorption_vector,
                                           backward_scattering_coeff,
                                           forward_scattering_coeff);
  }

  SingleScatteringDataImpl *interpolate_angles(
      math::VectorPtr<Scalar> lon_inc,
      math::VectorPtr<Scalar> lat_inc,
      math::VectorPtr<Scalar> lon_scat,
      std::shared_ptr<LatitudeGrid<Scalar>> lat_scat) const {
    auto phase_matrix = ScatteringDataFieldGridded<Scalar>(
        phase_matrix_.interpolate_angles(lon_inc, lat_inc, lon_scat, lat_scat));
    auto extinction_matrix = ScatteringDataFieldGridded<Scalar>(
        extinction_matrix_.interpolate_angles(lon_inc,
                                              lat_inc,
                                              dummy_grid_,
                                              dummy_lat_grid_));
    auto absorption_vector = ScatteringDataFieldGridded<Scalar>(
        absorption_vector_.interpolate_angles(lon_inc,
                                              lat_inc,
                                              dummy_grid_,
                                              dummy_lat_grid_));
    auto backward_scattering_coeff = ScatteringDataFieldGridded<Scalar>(
        backward_scattering_coeff_.interpolate_angles(lon_inc,
                                                      lat_inc,
                                                      dummy_grid_,
                                                      dummy_lat_grid_));
    auto forward_scattering_coeff = ScatteringDataFieldGridded<Scalar>(
        forward_scattering_coeff_.interpolate_angles(lon_inc,
                                                     lat_inc,
                                                     dummy_grid_,
                                                     dummy_lat_grid_));

    return new SingleScatteringDataGridded(f_grid_,
                                           t_grid_,
                                           phase_matrix,
                                           extinction_matrix,
                                           absorption_vector,
                                           backward_scattering_coeff,
                                           forward_scattering_coeff);
  }

  SingleScatteringDataImpl *downsample_scattering_angles(
      math::VectorPtr<Scalar> lon_scat,
      std::shared_ptr<LatitudeGrid<Scalar>> lat_scat) const {
    auto phase_matrix = ScatteringDataFieldGridded<Scalar>(
        phase_matrix_.downsample_scattering_angles(lon_scat, lat_scat));

    return new SingleScatteringDataGridded(f_grid_,
                                           t_grid_,
                                           phase_matrix,
                                           extinction_matrix_,
                                           absorption_vector_,
                                           backward_scattering_coeff_,
                                           forward_scattering_coeff_);
  }

  SingleScatteringDataImpl *downsample_lon_scat(
      math::VectorPtr<Scalar> lon_scat) const {
    auto phase_matrix = ScatteringDataFieldGridded<Scalar>(
        phase_matrix_.downsample_lon_scat(lon_scat));

    return new SingleScatteringDataGridded(f_grid_,
                                           t_grid_,
                                           phase_matrix,
                                           extinction_matrix_,
                                           absorption_vector_,
                                           backward_scattering_coeff_,
                                           forward_scattering_coeff_);
  }

  void operator+=(const SingleScatteringDataImpl *other) {
    auto converted = other->to_gridded();
    phase_matrix_ += converted->phase_matrix_;
    extinction_matrix_ += converted->extinction_matrix_;
    absorption_vector_ += converted->absorption_vector_;
    backward_scattering_coeff_ += converted->backward_scattering_coeff_;
    forward_scattering_coeff_ += converted->forward_scattering_coeff_;
  }

  SingleScatteringDataImpl *operator+(const SingleScatteringDataImpl *other) {
    auto result = this->copy();
    result->operator+=(other);
    return result;
  }

  void operator*=(Scalar c) {
    phase_matrix_ *= c;
    extinction_matrix_ *= c;
    absorption_vector_ *= c;
    backward_scattering_coeff_ *= c;
    forward_scattering_coeff_ *= c;
  }

  SingleScatteringDataImpl *operator*(Scalar c) const {
    auto result = this->copy();
    result->operator*=(c);
    return result;
  }

  void normalize(Scalar norm) {
      phase_matrix_.normalize(norm);
  }

  void set_number_of_scattering_coeffs(Index n) {
    phase_matrix_.set_number_of_scattering_coeffs(n);
    extinction_matrix_.set_number_of_scattering_coeffs(n);
    absorption_vector_.set_number_of_scattering_coeffs(n);
    backward_scattering_coeff_.set_number_of_scattering_coeffs(n);
    forward_scattering_coeff_.set_number_of_scattering_coeffs(n);
  }

  detail::ConversionPtr<const SingleScatteringDataGridded> to_gridded() const {
    return detail::make_conversion_ptr(this, false);
  }

  detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>> to_spectral(
      std::shared_ptr<sht::SHT>) const;
  detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>>
  to_spectral() const;
  detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>> to_spectral(
      Index l_max,
      Index m_max) const;
  detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>> to_spectral(
      Index l_max,
      Index m_max,
      Index n_lon,
      Index n_lat) const;

  SingleScatteringDataImpl *to_lab_frame(std::shared_ptr<math::Vector<double>> lat_inc,
                                         std::shared_ptr<math::Vector<double>> lon_scat,
                                         std::shared_ptr<LatitudeGrid<double>> lat_scat,
                                         Index stokes_dim) const {
      auto phase_matrix =
          phase_matrix_.to_lab_frame(lat_inc, lon_scat, lat_scat, stokes_dim);
      return new SingleScatteringDataGridded(f_grid_,
                                             t_grid_,
                                             phase_matrix,
                                             extinction_matrix_,
                                             absorption_vector_,
                                             backward_scattering_coeff_,
                                             forward_scattering_coeff_);
  }

  SingleScatteringDataImpl *to_lab_frame(Index n_lat_inc,
                                         Index n_lon_scat,
                                         Index stokes_dimension) const {
    auto phase_matrix =
        phase_matrix_.to_lab_frame(n_lat_inc, n_lon_scat, stokes_dimension);
    return new SingleScatteringDataGridded(f_grid_,
                                           t_grid_,
                                           phase_matrix,
                                           extinction_matrix_,
                                           absorption_vector_,
                                           backward_scattering_coeff_,
                                           forward_scattering_coeff_);
  }

  void set_stokes_dim(Index stokes_dim) {
    phase_matrix_.set_stokes_dim(stokes_dim);
    extinction_matrix_.set_stokes_dim(stokes_dim);
    absorption_vector_.set_stokes_dim(stokes_dim);
  }

  // explicit operator SingeScatteringDataSpectral();
  // explicit operator SingeScatteringDataFullySpectral();

  SingleScatteringDataImpl* regrid() const {
      auto n_lon_inc = phase_matrix_.get_n_lon_inc();
      auto n_lat_inc = phase_matrix_.get_n_lat_inc();
      auto n_lon_scat = phase_matrix_.get_n_lon_scat();
      auto n_lat_scat = phase_matrix_.get_n_lat_scat();
      n_lon_inc = std::max<Index>(n_lon_inc - n_lon_inc % 2, 1);
      n_lat_inc = std::max<Index>(n_lat_inc - n_lat_inc % 2, 1);
      n_lon_scat = std::max<Index>(n_lon_scat - n_lon_scat % 2, 1);
      n_lat_scat = std::max<Index>(n_lat_scat - n_lat_scat % 2, 1);

      auto lon_inc = std::make_shared<Vector>(sht::SHT::get_longitude_grid(n_lon_inc));
      auto lat_inc = std::make_shared<Vector>(sht::SHT::get_latitude_grid(n_lat_inc));
      auto lon_scat = std::make_shared<Vector>(sht::SHT::get_longitude_grid(n_lon_scat));
      auto lat_scat = std::make_shared<sht::SHT::LatGrid>(sht::SHT::get_latitude_grid(n_lat_scat));

      auto phase_matrix = phase_matrix_.interpolate_angles(lon_inc, lat_inc, lon_scat, lat_scat);
      auto extinction_matrix = extinction_matrix_.interpolate_angles(lon_inc, lat_inc, dummy_grid_, dummy_lat_grid_);
      auto absorption_vector = absorption_vector_.interpolate_angles(lon_inc, lat_inc, dummy_grid_, dummy_lat_grid_);
      auto backward_scattering_coeff = backward_scattering_coeff_.interpolate_angles(lon_inc, lat_inc, dummy_grid_, dummy_lat_grid_);
      auto forward_scattering_coeff = forward_scattering_coeff_.interpolate_angles(lon_inc, lat_inc, dummy_grid_, dummy_lat_grid_);

      return new SingleScatteringDataGridded(f_grid_,
                                             t_grid_,
                                             phase_matrix,
                                             extinction_matrix,
                                             absorption_vector,
                                             backward_scattering_coeff,
                                             forward_scattering_coeff);
  }

 private:
  VectorPtr dummy_grid_ = std::make_shared<Vector>(Vector::Constant(1, 1));
  LatitudeGridPtr dummy_lat_grid_ = std::make_shared<IrregularLatitudeGrid<double>>(Vector::Constant(1, 0.5 * M_PI));
  stokes::PhaseMatrix<ScatteringDataFieldGridded<Scalar>> phase_matrix_;
  stokes::ExtinctionMatrix<ScatteringDataFieldGridded<Scalar>> extinction_matrix_;
  stokes::AbsorptionVector<ScatteringDataFieldGridded<Scalar>> absorption_vector_;
  ScatteringDataFieldGridded<Scalar> backward_scattering_coeff_;
  ScatteringDataFieldGridded<Scalar> forward_scattering_coeff_;
};

////////////////////////////////////////////////////////////////////////////////
// Spectral single scattering data.
////////////////////////////////////////////////////////////////////////////////

template <typename Scalar>
class SingleScatteringDataSpectral : public SingleScatteringDataBase<Scalar>,
                                     public SingleScatteringDataImpl {
  using SingleScatteringDataBase<Scalar>::f_grid_;
  using SingleScatteringDataBase<Scalar>::n_freqs_;
  using SingleScatteringDataBase<Scalar>::n_temps_;
  using SingleScatteringDataBase<Scalar>::t_grid_;

  /** Check that dimensions of phase-matrix tensor are consistent with grids.
   * @throws runtime_error when phase matrix data is inconsistent with grids.
   */
  void check_consistency_phase_matrix() {
      auto data_tensor = phase_matrix_.get_data();
      if (data_tensor.dimension(0) != get_n_freqs()) {
          throw std::runtime_error(
              "Frequency grid of phase matrix is inconsistent with provided "
              "data tensor."
              );
      }
      if (data_tensor.dimension(1) != get_n_temps()) {
          throw std::runtime_error(
              "Temperature grid of phase matrix is inconsistent with provided "
              "data tensor."
              );
      }
      if (data_tensor.dimension(2) != get_n_lon_inc()) {
          throw std::runtime_error(
              "Incoming-angle longitude grid of phase matrix is inconsistent "
              "with provided data tensor."
              );
      }
      if (data_tensor.dimension(3) != get_n_lat_inc()) {
          throw std::runtime_error(
              "Incoming-angle latitude grid of phase matrix is inconsistent "
              "with provided data tensor."
              );
      }
      if (data_tensor.dimension(4) != sht_scat_->get_n_spectral_coeffs()) {
          throw std::runtime_error(
              "Scattering-angle longitude grid of phase matrix is inconsistent "
              "with provided data tensor."
              );
      }
  }

  /** Check that dimensions of extinction-matrix tensor are consistent with grids.
   * @throws runtime_error when extinction-matrix data is inconsistent with grids.
   */
  void check_consistency_extinction_matrix() {
      auto data_tensor = extinction_matrix_.get_data();
      if (data_tensor.dimension(0) != get_n_freqs()) {
          throw std::runtime_error(
              "Frequency grid of extinction matrix is inconsistent with provided "
              "data tensor."
              );
      }
      if (data_tensor.dimension(1) != get_n_temps()) {
          throw std::runtime_error(
              "Temperature grid of extinction matrix is inconsistent with "
              "provided data tensor."
              );
      }
      if (data_tensor.dimension(2) != get_n_lon_inc()) {
          throw std::runtime_error(
              "Incoming-angle longitude grid of extinction matrix is "
              "inconsistent with provided data tensor."
              );
      }
      if (data_tensor.dimension(3) != get_n_lat_inc()) {
          throw std::runtime_error(
              "Incoming-angle latitude grid of extinction matrix is "
              "inconsistent with provided data tensor."
              );
      }
      if (data_tensor.dimension(4) != 1) {
          throw std::runtime_error(
              "Expected scattering-angle longitude grid of absorption vector"
              "to be 1."
              );
      }
  }

  /** Check that dimensions of absorption-vector tensor are consistent with grids.
   * @throws runtime_error when extinction-matrix data is inconsistent with grids.
   */
  void check_consistency_absorption_vector() {
      auto data_tensor = absorption_vector_.get_data();
      if (data_tensor.dimension(0) != get_n_freqs()) {
          throw std::runtime_error(
              "Frequency grid of absorption vector is inconsistent with provided "
              "data tensor."
              );
      }
      if (data_tensor.dimension(1) != get_n_temps()) {
          throw std::runtime_error(
              "Temperature grid of absorption vector is inconsistent with "
              "provided data tensor."
              );
      }
      if (data_tensor.dimension(2) != get_n_lon_inc()) {
          throw std::runtime_error(
              "Incoming-angle longitude grid of absorption vector is "
              "inconsistent with provided data tensor."
              );
      }
      if (data_tensor.dimension(3) != get_n_lat_inc()) {
          throw std::runtime_error(
              "Incoming-angle latitude grid of absorption vector is "
              "inconsistent with provided data tensor."
              );
      }
      if (data_tensor.dimension(4) != 1) {
          throw std::runtime_error(
              "Expected scattering-angle longitude grid of absorption vector"
              "to be 1."
              );
      }
  }

 public:
  using Vector = math::Vector<Scalar>;
  using VectorPtr = std::shared_ptr<Vector>;
  using ShtPtr = std::shared_ptr<sht::SHT>;
  using DataTensor = math::Tensor<std::complex<Scalar>, 6>;
  using DataPtr = std::shared_ptr<DataTensor>;
  using ScatteringCoeff = math::Tensor<Scalar, 4>;
  using ScatteringCoeffPtr = std::shared_ptr<ScatteringCoeff>;
  using OtherScalar =
      std::conditional<std::is_same<Scalar, double>::value, float, double>;

  SingleScatteringDataSpectral(
      VectorPtr f_grid,
      VectorPtr t_grid,
      ScatteringDataFieldSpectral<Scalar> phase_matrix,
      ScatteringDataFieldSpectral<Scalar> extinction_matrix,
      ScatteringDataFieldSpectral<Scalar> absorption_vector,
      ScatteringDataFieldSpectral<Scalar> backward_scattering_coeff,
      ScatteringDataFieldSpectral<Scalar> forward_scattering_coeff)

      : SingleScatteringDataBase<Scalar>(f_grid, t_grid),
        sht_scat_(std::make_shared<sht::SHT>(phase_matrix.get_sht_scat())),
        phase_matrix_(phase_matrix),
        extinction_matrix_(extinction_matrix),
        absorption_vector_(absorption_vector),
        backward_scattering_coeff_(backward_scattering_coeff),
        forward_scattering_coeff_(forward_scattering_coeff) {
            check_consistency_phase_matrix();
            check_consistency_extinction_matrix();
            check_consistency_absorption_vector();
        }

  SingleScatteringDataSpectral(VectorPtr f_grid,
                               VectorPtr t_grid,
                               VectorPtr lon_inc,
                               VectorPtr lat_inc,
                               ShtPtr sht_scat,
                               DataPtr phase_matrix,
                               DataPtr extinction_matrix,
                               DataPtr absorption_vector,
                               DataPtr backward_scattering_coeff,
                               DataPtr forward_scattering_coeff)
      : SingleScatteringDataBase<Scalar>(f_grid, t_grid),
        sht_scat_(sht_scat),
        phase_matrix_{f_grid, t_grid, lon_inc, lat_inc, sht_scat, phase_matrix},
        extinction_matrix_{f_grid,
                           t_grid,
                           lon_inc,
                           lat_inc,
                           sht_dummy_,
                           extinction_matrix},
        absorption_vector_{f_grid,
                           t_grid,
                           lon_inc,
                           lat_inc,
                           sht_dummy_,
                           absorption_vector},
        backward_scattering_coeff_{f_grid,
                                   t_grid,
                                   lon_inc,
                                   lat_inc,
                                   sht_dummy_,
                                   backward_scattering_coeff},
        forward_scattering_coeff_{f_grid,
                                  t_grid,
                                  lon_inc,
                                  lat_inc,
                                  sht_dummy_,
                                  forward_scattering_coeff} {

        }

  ParticleType get_particle_type() const { return phase_matrix_.get_particle_type(); }
  DataFormat get_data_format() const { return phase_matrix_.get_data_format(); }

  const math::Vector<double>& get_f_grid() const { return *f_grid_; }
  const math::Vector<double>& get_t_grid() const { return *t_grid_; }
  math::Vector<double> get_lon_inc() const { return phase_matrix_.get_lon_inc(); }
  math::Vector<double> get_lat_inc() const { return phase_matrix_.get_lat_inc(); }
  math::Vector<double> get_lon_scat() const { return phase_matrix_.get_lon_scat(); }
  math::Vector<double> get_lat_scat() const { return phase_matrix_.get_lat_scat(); }

  math::Index get_n_freqs() const { return phase_matrix_.get_n_freqs(); }
  math::Index get_n_temps() const { return phase_matrix_.get_n_temps(); }
  math::Index get_n_lon_inc() const { return phase_matrix_.get_n_lon_inc(); }
  math::Index get_n_lat_inc() const { return phase_matrix_.get_n_lat_inc(); }
  math::Index get_n_lon_scat() const { return phase_matrix_.get_n_lon_scat(); }
  math::Index get_l_max_scat() const { return phase_matrix_.get_sht_scat_params()[0]; }
  math::Index get_m_max_scat() const { return phase_matrix_.get_sht_scat_params()[1]; }
  math::Index get_n_lat_scat() const { return phase_matrix_.get_n_lat_scat(); }
  math::Index get_stokes_dim() const { return phase_matrix_.get_stokes_dim(); }

  void set_data(Index f_index,
                Index t_index,
                const SingleScatteringDataImpl &other) {
    auto converted = other.to_spectral();
    phase_matrix_.set_data(f_index, t_index, converted->phase_matrix_);
    extinction_matrix_.set_data(f_index,
                                t_index,
                                converted->extinction_matrix_);
    absorption_vector_.set_data(f_index,
                                t_index,
                                converted->absorption_vector_);
    backward_scattering_coeff_.set_data(f_index,
                                        t_index,
                                        converted->backward_scattering_coeff_);
    forward_scattering_coeff_.set_data(f_index,
                                       t_index,
                                       converted->forward_scattering_coeff_);
  }

  math::Tensor<Scalar, 6> get_phase_function() const {
    auto phase_matrix_gridded = phase_matrix_.to_gridded();
    stokes::PhaseMatrix<ScatteringDataFieldGridded<Scalar>> gridded = phase_matrix_.to_gridded();
    return gridded.get_phase_function();
  }

  math::Tensor<double, 4> get_phase_function_max_inc() const {
      return phase_matrix_.get_phase_function_max_inc();
  }

  math::Tensor<double, 4> get_phase_function_max_scat() const {
      return phase_matrix_.get_phase_function_max_scat();
  }

  math::Tensor<std::complex<Scalar>, 5> get_phase_function_spectral() const {
    return phase_matrix_.get_phase_function();
  }

  math::Tensor<Scalar, 8> get_phase_matrix(Index stokes_dim) const {
    stokes::PhaseMatrix<ScatteringDataFieldGridded<Scalar>> gridded =
        phase_matrix_.to_gridded();
    return gridded.get_phase_matrix(stokes_dim);
  }

  math::Matrix<double> get_phase_matrix(
      double frequency,
      double temperature,
      double lon_inc,
      double lat_inc,
      double lon_scat,
      double lat_scat,
      Index stokes_dim
      ) const {
      return phase_matrix_.get_phase_matrix(
          frequency, temperature, lon_inc, lat_inc, lon_scat, lat_scat, stokes_dim
          );
  }

  math::Tensor<Scalar, 7> get_phase_matrix_data() const {
      stokes::PhaseMatrix<ScatteringDataFieldGridded<Scalar>> gridded = phase_matrix_.to_gridded();
      return gridded.get_data();
  }

  math::Tensor<std::complex<Scalar>, 6> get_phase_matrix_data_spectral() const {
      return phase_matrix_.get_data();
  }

  math::Tensor<Scalar, 7> get_extinction_matrix_data() const {
      auto gridded = extinction_matrix_.to_gridded();
      return gridded.get_data();
}

math::Tensor<Scalar, 6> get_extinction_coeff() const {
    stokes::ExtinctionMatrix<ScatteringDataFieldGridded<Scalar>> extinction_matrix_gridded = extinction_matrix_.to_gridded();
    return extinction_matrix_gridded.get_extinction_coeff();
}

math::Tensor<Scalar, 8> get_extinction_matrix(Index stokes_dim) const {
    stokes::ExtinctionMatrix<ScatteringDataFieldGridded<Scalar>> extinction_matrix_gridded = extinction_matrix_.to_gridded();
    return extinction_matrix_gridded.get_extinction_matrix(stokes_dim);
  }

math::Matrix<double> get_extinction_matrix(
    double frequency,
    double temperature,
    double lon_inc,
    double lat_inc,
    Index stokes_dim) const {
    return extinction_matrix_.get_extinction_matrix(
        frequency, temperature, lon_inc, lat_inc, stokes_dim
        );
}

  math::Tensor<Scalar, 7> get_absorption_vector_data() const {
      auto absorption_vector_gridded = absorption_vector_.to_gridded();
      return absorption_vector_gridded.get_data();
  }

  math::Tensor<Scalar, 6> get_absorption_coeff() const {
      stokes::AbsorptionVector<ScatteringDataFieldGridded<Scalar>> absorption_vector_gridded = absorption_vector_.to_gridded();
      return absorption_vector_gridded.get_absorption_coeff();
  }

  math::Tensor<Scalar, 7> get_absorption_vector(Index stokes_dim) const {
    stokes::AbsorptionVector<ScatteringDataFieldGridded<Scalar>> absorption_vector_gridded = absorption_vector_.to_gridded();
    return absorption_vector_gridded.get_absorption_vector(stokes_dim);
  }

  math::Vector<double> get_absorption_vector(
      double temperature,
      double frequency,
      double lon_inc,
      double lat_inc,
      Index stokes_dim) const {
      return absorption_vector_.get_absorption_vector(
          temperature,
          frequency,
          lon_inc,
          lat_inc,
          stokes_dim);
  }

  math::Tensor<Scalar, 7> get_backward_scattering_coeff() const {
    auto backward_scattering_coeff_gridded =
        backward_scattering_coeff_.to_gridded();
    return backward_scattering_coeff_gridded.get_data();
  }

  math::Tensor<Scalar, 7> get_forward_scattering_coeff() const {
    auto forward_scattering_coeff_gridded =
        forward_scattering_coeff_.to_gridded();
    return forward_scattering_coeff_gridded.get_data();
  }

  SingleScatteringDataSpectral *copy() const {
    return new SingleScatteringDataSpectral(f_grid_,
                                            t_grid_,
                                            phase_matrix_.copy(),
                                            extinction_matrix_.copy(),
                                            absorption_vector_.copy(),
                                            backward_scattering_coeff_.copy(),
                                            forward_scattering_coeff_.copy());
  }

  SingleScatteringDataImpl *interpolate_frequency(
      math::VectorPtr<double> frequencies) const {
    auto phase_matrix = ScatteringDataFieldSpectral<Scalar>(
        phase_matrix_.interpolate_frequency(frequencies));
    auto extinction_matrix = ScatteringDataFieldSpectral<Scalar>(
        extinction_matrix_.interpolate_frequency(frequencies));
    auto absorption_vector = ScatteringDataFieldSpectral<Scalar>(
        absorption_vector_.interpolate_frequency(frequencies));
    auto backward_scattering_coeff = ScatteringDataFieldSpectral<Scalar>(
        backward_scattering_coeff_.interpolate_frequency(frequencies));
    auto forward_scattering_coeff = ScatteringDataFieldSpectral<Scalar>(
        forward_scattering_coeff_.interpolate_frequency(frequencies));
    return new SingleScatteringDataSpectral(frequencies,
                                            t_grid_,
                                            phase_matrix,
                                            extinction_matrix,
                                            absorption_vector,
                                            backward_scattering_coeff,
                                            forward_scattering_coeff);
  }

  SingleScatteringDataImpl *interpolate_temperature(
      math::VectorPtr<Scalar> temperatures,
      bool extrapolate=false) const {
    auto phase_matrix = ScatteringDataFieldSpectral<Scalar>(
        phase_matrix_.interpolate_temperature(temperatures, extrapolate));
    auto extinction_matrix = ScatteringDataFieldSpectral<Scalar>(
        extinction_matrix_.interpolate_temperature(temperatures, extrapolate));
    auto absorption_vector = ScatteringDataFieldSpectral<Scalar>(
        absorption_vector_.interpolate_temperature(temperatures, extrapolate));
    auto backward_scattering_coeff = ScatteringDataFieldSpectral<Scalar>(
        backward_scattering_coeff_.interpolate_temperature(temperatures, extrapolate));
    auto forward_scattering_coeff = ScatteringDataFieldSpectral<Scalar>(
        forward_scattering_coeff_.interpolate_temperature(temperatures, extrapolate));

    return new SingleScatteringDataSpectral(f_grid_,
                                            temperatures,
                                            phase_matrix,
                                            extinction_matrix,
                                            absorption_vector,
                                            backward_scattering_coeff,
                                            forward_scattering_coeff);
  }

  SingleScatteringDataImpl *interpolate_angles(
      math::VectorPtr<Scalar> lon_inc,
      math::VectorPtr<Scalar> lat_inc,
      math::VectorPtr<Scalar> /*lon_scat*/,
      std::shared_ptr<LatitudeGrid<Scalar>> /*lat_scat*/) const {
    auto phase_matrix = ScatteringDataFieldSpectral<Scalar>(
        phase_matrix_.interpolate_angles(lon_inc, lat_inc));
    auto extinction_matrix = ScatteringDataFieldSpectral<Scalar>(
        extinction_matrix_.interpolate_angles(lon_inc, lat_inc));
    auto absorption_vector = ScatteringDataFieldSpectral<Scalar>(
        absorption_vector_.interpolate_angles(lon_inc, lat_inc));
    auto backward_scattering_coeff = ScatteringDataFieldSpectral<Scalar>(
        backward_scattering_coeff_.interpolate_angles(lon_inc, lat_inc));
    auto forward_scattering_coeff = ScatteringDataFieldSpectral<Scalar>(
        forward_scattering_coeff_.interpolate_angles(lon_inc, lat_inc));

    return new SingleScatteringDataSpectral(f_grid_,
                                            t_grid_,
                                            phase_matrix,
                                            extinction_matrix,
                                            absorption_vector,
                                            backward_scattering_coeff,
                                            forward_scattering_coeff);
  }

  SingleScatteringDataImpl *downsample_scattering_angles(
      math::VectorPtr<Scalar> /*lon_scat*/,
      std::shared_ptr<LatitudeGrid<Scalar>> /*lat_scat*/) const {
      return copy();
  }

  SingleScatteringDataImpl *downsample_lon_scat(
      math::VectorPtr<Scalar> /*lon_scat*/) const {
      return copy();
  }

  void operator+=(const SingleScatteringDataImpl *other) {
    auto converted = other->to_spectral();
    phase_matrix_ += converted->phase_matrix_;
    extinction_matrix_ += converted->extinction_matrix_;
    absorption_vector_ += converted->absorption_vector_;
    backward_scattering_coeff_ += converted->backward_scattering_coeff_;
    forward_scattering_coeff_ += converted->forward_scattering_coeff_;
  }

  SingleScatteringDataImpl *operator+(const SingleScatteringDataImpl *other) {
    auto result = this->copy();
    result->operator+=(other);
    return result;
  }

  void operator*=(Scalar c) {
    phase_matrix_ *= c;
    extinction_matrix_ *= c;
    absorption_vector_ *= c;
    backward_scattering_coeff_ *= c;
    forward_scattering_coeff_ *= c;
  }

  SingleScatteringDataImpl *operator*(Scalar c) const {
    auto result = this->copy();
    result->operator*=(c);
    return result;
  }

  void normalize(Scalar norm) { phase_matrix_.normalize(norm); }

  void set_number_of_scattering_coeffs(Index n) {
      phase_matrix_.set_number_of_scattering_coeffs(n);
      extinction_matrix_.set_number_of_scattering_coeffs(n);
      absorption_vector_.set_number_of_scattering_coeffs(n);
      backward_scattering_coeff_.set_number_of_scattering_coeffs(n);
      forward_scattering_coeff_.set_number_of_scattering_coeffs(n);
  }

  SingleScatteringDataImpl* regrid() const {
      auto n_lon_inc = phase_matrix_.get_n_lon_inc();
      auto n_lat_inc = phase_matrix_.get_n_lat_inc();
      n_lon_inc = std::max<Index>(n_lon_inc - n_lon_inc % 2, 1);
      n_lat_inc = std::max<Index>(n_lat_inc - n_lat_inc % 2, 1);

      auto lon_inc = std::make_shared<Vector>(sht::SHT::get_longitude_grid(n_lon_inc));
      auto lat_inc = std::make_shared<Vector>(sht::SHT::get_latitude_grid(n_lat_inc));

      auto phase_matrix = phase_matrix_.interpolate_angles(lon_inc, lat_inc);
      auto extinction_matrix = extinction_matrix_.interpolate_angles(lon_inc, lat_inc);
      auto absorption_vector = absorption_vector_.interpolate_angles(lon_inc, lat_inc);
      auto backward_scattering_coeff = backward_scattering_coeff_.interpolate_angles(lon_inc, lat_inc);
      auto forward_scattering_coeff = forward_scattering_coeff_.interpolate_angles(lon_inc, lat_inc);

      return new SingleScatteringDataSpectral(f_grid_,
                                              t_grid_,
                                              phase_matrix,
                                              extinction_matrix,
                                              absorption_vector,
                                              backward_scattering_coeff,
                                              forward_scattering_coeff);
  }

  detail::ConversionPtr<const SingleScatteringDataGridded<Scalar>> to_gridded()
      const;

  detail::ConversionPtr<const SingleScatteringDataSpectral> to_spectral()
      const {
    return detail::make_conversion_ptr(this, false);
  }

  detail::ConversionPtr<const SingleScatteringDataSpectral> to_spectral(
      Index l_max,
      Index m_max) const;

  detail::ConversionPtr<const SingleScatteringDataSpectral> to_spectral(
      Index l_max,
      Index m_max,
      Index n_lat,
      Index n_lon) const;

  SingleScatteringDataImpl *to_lab_frame(std::shared_ptr<math::Vector<double>> lat_inc,
                                         std::shared_ptr<math::Vector<double>> lon_scat,
                                         std::shared_ptr<LatitudeGrid<double>> lat_scat,
                                         Index stokes_dim) const {
      stokes::PhaseMatrix<ScatteringDataFieldGridded<Scalar>> phase_matrix = phase_matrix_.to_gridded();
      auto phase_matrix_lab = phase_matrix.to_lab_frame(lat_inc, lon_scat, lat_scat, stokes_dim);
      auto n_lon_scat = lon_scat->size();
      Index l_max = (n_lon_scat - 2 + (n_lon_scat % 2)) / 2;
      auto sht_scat = std::make_shared<sht::SHT>(l_max, l_max, n_lon_scat, n_lon_scat);
      return new SingleScatteringDataSpectral(f_grid_,
                                              t_grid_,
                                              phase_matrix_lab.to_spectral(sht_scat),
                                              extinction_matrix_,
                                              absorption_vector_,
                                              backward_scattering_coeff_,
                                              forward_scattering_coeff_);
  }

  SingleScatteringDataImpl *to_lab_frame(Index n_lat_inc,
                                         Index n_lon_scat,
                                         Index stokes_dimension) const {
      stokes::PhaseMatrix<ScatteringDataFieldGridded<Scalar>> phase_matrix = phase_matrix_.to_gridded();
      auto phase_matrix_lab = phase_matrix.to_lab_frame(n_lat_inc, n_lon_scat, stokes_dimension);
      Index l_max = (n_lon_scat - 2 + (n_lon_scat % 2)) / 2;
      auto sht_scat = std::make_shared<sht::SHT>(l_max, l_max, n_lon_scat, n_lon_scat);
      return new SingleScatteringDataSpectral(f_grid_,
                                              t_grid_,
                                              phase_matrix_lab.to_spectral(sht_scat),
                                              extinction_matrix_,
                                              absorption_vector_,
                                              backward_scattering_coeff_,
                                              forward_scattering_coeff_);
  }

  void set_stokes_dim(Index stokes_dim) {
      phase_matrix_.set_stokes_dim(stokes_dim);
      extinction_matrix_.set_stokes_dim(stokes_dim);
      absorption_vector_.set_stokes_dim(stokes_dim);
  }
  // explicit operator SingeScatteringDataSpectral();
  // explicit operator SingeScatteringDataFullySpectral();

 private:
  ShtPtr sht_scat_;
  ShtPtr sht_dummy_ = std::make_shared<sht::SHT>(0, 0, 1, 1);
  stokes::PhaseMatrix<ScatteringDataFieldSpectral<Scalar>> phase_matrix_;
  stokes::ExtinctionMatrix<ScatteringDataFieldSpectral<Scalar>> extinction_matrix_;
  stokes::AbsorptionVector<ScatteringDataFieldSpectral<Scalar>> absorption_vector_;
  ScatteringDataFieldSpectral<Scalar> backward_scattering_coeff_;
  ScatteringDataFieldSpectral<Scalar> forward_scattering_coeff_;
};

// pxx :: hide
template <typename Scalar>
class SingleScatteringDataFullySpectral {};

////////////////////////////////////////////////////////////////////////////////
// Member function definitions.
////////////////////////////////////////////////////////////////////////////////

template <typename Scalar>
detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>>
SingleScatteringDataGridded<Scalar>::to_spectral(
    std::shared_ptr<sht::SHT> sht) const {
  using ReturnType = const SingleScatteringDataSpectral<Scalar>;
  auto phase_matrix =
      ScatteringDataFieldSpectral<Scalar>(phase_matrix_.to_spectral(sht));
  auto extinction_matrix =
      ScatteringDataFieldSpectral<Scalar>(extinction_matrix_.to_spectral());
  auto absorption_vector =
      ScatteringDataFieldSpectral<Scalar>(absorption_vector_.to_spectral());
  auto backward_scattering_coeff = ScatteringDataFieldSpectral<Scalar>(
      backward_scattering_coeff_.to_spectral());
  auto forward_scattering_coeff = ScatteringDataFieldSpectral<Scalar>(
      forward_scattering_coeff_.to_spectral());
  return detail::make_conversion_ptr<ReturnType>(
      new SingleScatteringDataSpectral<Scalar>(f_grid_,
                                               t_grid_,
                                               phase_matrix,
                                               extinction_matrix,
                                               absorption_vector,
                                               backward_scattering_coeff,
                                               forward_scattering_coeff),
      true);
}

template <typename Scalar>
detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>>
SingleScatteringDataGridded<Scalar>::to_spectral() const {
  auto sht_params = phase_matrix_.get_sht_scat_params();
  auto sht_scat = std::make_shared<sht::SHT>(sht_params[0],
                                             sht_params[1],
                                             phase_matrix_.get_n_lon_scat(),
                                             phase_matrix_.get_n_lat_scat());
  return to_spectral(sht_scat);
}

template <typename Scalar>
detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>>
SingleScatteringDataGridded<Scalar>::to_spectral(Index l_max,
                                                 Index m_max) const {
  auto sht_scat = std::make_shared<sht::SHT>(l_max,
                                             m_max,
                                             phase_matrix_.get_n_lon_scat(),
                                             phase_matrix_.get_n_lat_scat());
  return to_spectral(sht_scat);
}

template <typename Scalar>
    detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>>
    SingleScatteringDataGridded<Scalar>::to_spectral(Index l_max,
                                                     Index m_max,
                                                     Index n_lon,
                                                     Index n_lat) const {
    auto sht_scat = std::make_shared<sht::SHT>(l_max,
                                               m_max,
                                               n_lon,
                                               n_lat);
    return to_spectral(sht_scat);
}

template <typename Scalar>
detail::ConversionPtr<const SingleScatteringDataGridded<Scalar>>
SingleScatteringDataSpectral<Scalar>::to_gridded() const {
  using ReturnType = const SingleScatteringDataGridded<Scalar>;
  auto phase_matrix =
      ScatteringDataFieldGridded<Scalar>(phase_matrix_.to_gridded());
  auto extinction_matrix =
      ScatteringDataFieldGridded<Scalar>(extinction_matrix_.to_gridded());
  auto absorption_vector =
      ScatteringDataFieldGridded<Scalar>(absorption_vector_.to_gridded());
  auto backward_scattering_coeff = ScatteringDataFieldGridded<Scalar>(
      backward_scattering_coeff_.to_gridded());
  auto forward_scattering_coeff = ScatteringDataFieldGridded<Scalar>(
      forward_scattering_coeff_.to_gridded());
  auto lon_scat =
      std::make_shared<math::Vector<Scalar>>(phase_matrix_.get_lon_scat());
  auto lat_scat =
      std::make_shared<math::Vector<Scalar>>(phase_matrix_.get_lat_scat());
  return detail::make_conversion_ptr<ReturnType>(
      new SingleScatteringDataGridded<Scalar>(f_grid_,
                                              t_grid_,
                                              phase_matrix,
                                              extinction_matrix,
                                              absorption_vector,
                                              backward_scattering_coeff,
                                              forward_scattering_coeff),
      true);
}

template <typename Scalar>
detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>>
SingleScatteringDataSpectral<Scalar>::to_spectral(Index l_max,
                                                  Index m_max,
                                                  Index n_lon,
                                                  Index n_lat) const {
  using ReturnType = const SingleScatteringDataSpectral<Scalar>;
  auto phase_matrix = ScatteringDataFieldSpectral<Scalar>(
      phase_matrix_.to_spectral(l_max, m_max, n_lon, n_lat));
  auto extinction_matrix =
      ScatteringDataFieldSpectral<Scalar>(extinction_matrix_);
  auto absorption_vector =
      ScatteringDataFieldSpectral<Scalar>(absorption_vector_);
  auto backward_scattering_coeff =
      ScatteringDataFieldSpectral<Scalar>(backward_scattering_coeff_);
  auto forward_scattering_coeff =
      ScatteringDataFieldSpectral<Scalar>(forward_scattering_coeff_);
  return detail::make_conversion_ptr<ReturnType>(
      new SingleScatteringDataSpectral<Scalar>(f_grid_,
                                               t_grid_,
                                               phase_matrix,
                                               extinction_matrix,
                                               absorption_vector,
                                               backward_scattering_coeff,
                                               forward_scattering_coeff),
      true);
}

template <typename Scalar>
detail::ConversionPtr<const SingleScatteringDataSpectral<Scalar>>
SingleScatteringDataSpectral<Scalar>::to_spectral(Index l_max,
                                                  Index m_max) const {
  auto sht = phase_matrix_.get_sht_scat();
  auto n_lat = sht.get_n_latitudes();
  auto n_lon = sht.get_n_longitudes();
  return to_spectral(l_max, m_max, n_lon, n_lat);
}

}  // namespace scattering
