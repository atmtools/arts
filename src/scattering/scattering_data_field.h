/** \file scattering_data_field.h
 *
 * This files contains three classes used to represent scattering data fields in
 * three different formats. These classes are used as the basic container to
 * hold scattering data. The scattering data itself is typically multi-dimensional
 * (depending on the stokes dimension) and in general depends on frequency,
 * temperature, and incoming and scattering angles.
 *
 * The scattering data field classes differ in how the dependency w.r.t. to the
 * incoming and scattering angles is represented:
 *
 * The ScatteringDataGridded class represents the angle-dependencies using
 * two-dimensional grids over the  longitudinal and latitudinal components
 * of incoming and scattering angles.
 *
 * The ScatteringDataFieldSpectral represents the scattering-angle dpendency
 * using spherical harmonics, whereas the ScatteringDataFieldFullySpectral
 * also uses a spherical-harmonics transform to encode the dependency on
 * incoming angle.
 *
 * Finally, ScatteringDataFieldFullySpectral represents both the incoming
 * and scattering angle dependencies using spherical harmonics.
 *
 * @author Simon Pfreundschuh, 2020 - 2022
 */
#pragma once

#include <scattering/maths.h>
#include <scattering/integration.h>
#include <scattering/interpolation.h>
#include <scattering/sht.h>
#include <scattering/array.h>

#include <cassert>
#include <memory>

namespace scattering {

using math::Index;

// pxx :: export
enum class DataFormat {Gridded = 0, Spectral = 1, FullySpectral = 2 };
// pxx :: export
enum class ParticleType { Random = 0, AzimuthallyRandom = 1, General = 2 };

// pxx :: hide
template <typename Scalar>
class ScatteringDataFieldGridded;
template <typename Scalar>
class ScatteringDataFieldSpectral;
template <typename Scalar>
class ScatteringDataFieldFullySpectral;

/** ScatteringDataField base class.
 *
 * Holds information on the size of the angular grids and the
 * the type of scattering data.
 *
 */
class ScatteringDataFieldBase {
 public:
  /// Size of the frequency grid.
  Index get_n_freqs() const { return n_freqs_; }
  /// Size of the temperature grid.
  Index get_n_temps() const { return n_temps_; }
  /// Size of incoming longitude  angle grid.
  Index get_n_lon_inc() const { return n_lon_inc_; }
  /// Size of incoming latitude  angle grid.
  Index get_n_lat_inc() const { return n_lat_inc_; }
  /// Size of scattering longitude  angle grid.
  Index get_n_lon_scat() const { return n_lon_scat_; }
  /// Size of scattering latitude  angle grid.
  Index get_n_lat_scat() const { return n_lat_scat_; }

 protected:
  ScatteringDataFieldBase(Index n_freqs,
                          Index n_temps,
                          Index n_lon_inc,
                          Index n_lat_inc,
                          Index n_lon_scat,
                          Index n_lat_scat)
      : n_freqs_(n_freqs),
        n_temps_(n_temps),
        n_lon_inc_(n_lon_inc),
        n_lat_inc_(n_lat_inc),
        n_lon_scat_(n_lon_scat),
        n_lat_scat_(n_lat_scat) {}

 protected:
  Index n_freqs_;
  Index n_temps_;
  Index n_lon_inc_;
  Index n_lat_inc_;
  Index n_lon_scat_;
  Index n_lat_scat_;
  ParticleType type_;
};

////////////////////////////////////////////////////////////////////////////////
// Gridded format
////////////////////////////////////////////////////////////////////////////////
// pxx :: export
// pxx :: instance(["double"])
/** Gridded scattering data field.
 *
 * Holds scattering data in gridded format. The data is in this case given
 * in the form of a rank-7 tensor with the dimensions corresponding to the
 * following grids:
 *     1: Frequency
 *     2: Temperature
 *     3: Incoming azimuth angle
 *     4: Incoming zenith angle
 *     5: Scattering azimuth angle
 *     6: Scattering zenith angle
 *     7: Stokes-dimension of scattering data.
 */
template <typename Scalar_>
class ScatteringDataFieldGridded : public ScatteringDataFieldBase {
 public:

  using ScatteringDataFieldBase::n_freqs_;
  using ScatteringDataFieldBase::n_lat_inc_;
  using ScatteringDataFieldBase::n_lat_scat_;
  using ScatteringDataFieldBase::n_lon_inc_;
  using ScatteringDataFieldBase::n_lon_scat_;
  using ScatteringDataFieldBase::n_temps_;
  using ScatteringDataFieldBase::get_n_lon_inc;
  using ScatteringDataFieldBase::get_n_lat_inc;
  using ScatteringDataFieldBase::get_n_lon_scat;
  using ScatteringDataFieldBase::get_n_lat_scat;

  using Scalar = Scalar_;
  using Coefficient = Scalar;
  using Vector = math::Vector<Scalar>;
  using VectorMap = math::VectorMap<Scalar>;
  using VectorPtr = std::shared_ptr<const math::Vector<Scalar>>;
  using LatitudeGridPtr = std::shared_ptr<const LatitudeGrid<Scalar>>;
  using ConstVectorMap = math::ConstVectorMap<Scalar>;
  using Matrix = math::Matrix<Scalar>;
  using MatrixMap = math::MatrixMap<Scalar>;
  using ConstMatrixMap = math::ConstMatrixMap<Scalar>;
  using OneAngle = math::MatrixFixedRows<Scalar, 1>;
  using ThreeAngles = math::MatrixFixedRows<Scalar, 3>;
  using FourAngles = math::MatrixFixedRows<Scalar, 4>;

  template <math::Index rank>
  using Tensor = math::Tensor<Scalar, rank>;
  template <math::Index rank>
  using TensorMap = math::TensorMap<Scalar, rank>;
  template <math::Index rank>
  using ConstTensorMap = math::ConstTensorMap<Scalar, rank>;
  using DataTensor = math::Tensor<Scalar, 7>;
  using DataTensorPtr = std::shared_ptr<DataTensor>;

  static constexpr Index coeff_dim = 6;
  static constexpr Index rank = 7;

  // pxx :: hide
  /** Create gridded scattering data field.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @lon_inc The incoming azimuth angle
   * @lat_inc The incoming zenith angle
   * @lon_scat The scattering zenith angle
   * @lat_scat The scattering azimuth angle
   * @data The tensor containing the scattering data.
   */
  ScatteringDataFieldGridded(VectorPtr f_grid,
                             VectorPtr t_grid,
                             VectorPtr lon_inc,
                             VectorPtr lat_inc,
                             VectorPtr lon_scat,
                             LatitudeGridPtr lat_scat,
                             DataTensorPtr data)
      : ScatteringDataFieldBase(f_grid->size(),
                                t_grid->size(),
                                lon_inc->size(),
                                lat_inc->size(),
                                lon_scat->size(),
                                lat_scat->size()),
        f_grid_(f_grid),
        t_grid_(t_grid),
        lon_inc_(lon_inc),
        lat_inc_(lat_inc),
        lon_scat_(lon_scat),
        lat_scat_(lat_scat),
        f_grid_map_(f_grid->data(), n_freqs_),
        t_grid_map_(t_grid->data(), n_temps_),
        lon_inc_map_(lon_inc->data(), n_lon_inc_),
        lat_inc_map_(lat_inc->data(), n_lat_inc_),
        lon_scat_map_(lon_scat->data(), n_lon_scat_),
        lat_scat_map_(lat_scat->data(), n_lat_scat_),
        data_(data) {}

  /** Create gridded scattering data field.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @lon_inc The incoming azimuth angle
   * @lat_inc The incoming zenith angle
   * @lon_scat The scattering zenith angle
   * @lat_scat The scattering azimuth angle
   * @data The tensor containing the scattering data.
   */
  ScatteringDataFieldGridded(Vector f_grid,
                             Vector t_grid,
                             Vector &lon_inc,
                             Vector &lat_inc,
                             Vector &lon_scat,
                             Vector &lat_scat,
                             Tensor<7> &data)
      : ScatteringDataFieldBase(f_grid.size(),
                                t_grid.size(),
                                lon_inc.size(),
                                lat_inc.size(),
                                lon_scat.size(),
                                lat_scat.size()),
        f_grid_(std::make_shared<Vector>(f_grid)),
        t_grid_(std::make_shared<Vector>(t_grid)),
        lon_inc_(std::make_shared<Vector>(lon_inc)),
        lat_inc_(std::make_shared<Vector>(lat_inc)),
        lon_scat_(std::make_shared<Vector>(lon_scat)),
        lat_scat_(std::make_shared<IrregularLatitudeGrid<Scalar>>(lat_scat)),
        f_grid_map_(f_grid_->data(), n_freqs_),
        t_grid_map_(t_grid_->data(), n_temps_),
        lon_inc_map_(lon_inc_->data(), n_lon_inc_),
        lat_inc_map_(lat_inc_->data(), n_lat_inc_),
        lon_scat_map_(lon_scat_->data(), n_lon_scat_),
        lat_scat_map_(lat_scat_->data(), n_lat_scat_),
        data_(std::make_shared<DataTensor>(data)) {}

  /** Create empty gridded scattering data field.
   *
   * This constructor is useful to pre-allocate data for sequentially
   * loading scattering data from multiple files or that is defined
   * on different grids.
   *
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @lon_inc The incoming azimuth angle
   * @lat_inc The incoming zenith angle
   * @lon_scat The scattering zenith angle
   * @lat_scat The scattering azimuth angle
   * @n_element The number of scattering-data elements, e.g.
   * phase matrix components that need to be stored for this
   * data field.
   */
  ScatteringDataFieldGridded(Vector f_grid,
                             Vector t_grid,
                             Vector &lon_inc,
                             Vector &lat_inc,
                             Vector &lon_scat,
                             Vector &lat_scat,
                             Index n_elements)
      : ScatteringDataFieldBase(f_grid.size(),
                                t_grid.size(),
                                lon_inc.size(),
                                lat_inc.size(),
                                lon_scat.size(),
                                lat_scat.size()),
        f_grid_(std::make_shared<Vector>(f_grid)),
        t_grid_(std::make_shared<Vector>(t_grid)),
        lon_inc_(std::make_shared<Vector>(lon_inc)),
        lat_inc_(std::make_shared<Vector>(lat_inc)),
        lon_scat_(std::make_shared<Vector>(lon_scat)),
        lat_scat_(std::make_shared<IrregularLatitudeGrid<Scalar>>(lat_scat)),
        f_grid_map_(f_grid_->data(), n_freqs_),
        t_grid_map_(t_grid_->data(), n_temps_),
        lon_inc_map_(lon_inc_->data(), n_lon_inc_),
        lat_inc_map_(lat_inc_->data(), n_lat_inc_),
        lon_scat_map_(lon_scat_->data(), n_lon_scat_),
        lat_scat_map_(lat_scat_->data(), n_lat_scat_),
        data_(std::make_shared<DataTensor>(std::array<Index, 7>{n_freqs_,
                                                                n_temps_,
                                                                n_lon_inc_,
                                                                n_lat_inc_,
                                                                n_lon_scat_,
                                                                n_lat_scat_,
                                                                n_elements})) {}
  /// Shallow copy of the ScatteringDataField.
  ScatteringDataFieldGridded(const ScatteringDataFieldGridded &) = default;

  /// Enum representing the data format.
  DataFormat get_data_format() const { return DataFormat::Gridded; }

  /// The number of scattering-data coefficients.
  Index get_n_coeffs() const { return data_->dimension(6); }
  /// Largest SHT parameters satisfying shtns aliasing requirements for
  /// scattering angle.
  std::array<Index, 4> get_sht_scat_params() const {
      return sht::SHT::get_params(n_lon_scat_, n_lat_scat_);
  }
  /// Largest SHT parameters satisfying shtns aliasing requirements for
  /// incoming angle.
  std::array<Index, 4> get_sht_inc_params() const {
      return sht::SHT::get_params(n_lon_scat_, n_lat_scat_);
  }
  /// The frequency grid.
  const math::Vector<double> &get_f_grid() const { return *f_grid_; }
  /// The temperature grid.
  const math::Vector<double>& get_t_grid() const { return *t_grid_; }
  /// The incoming-angle longitude grid.
  math::Vector<double> get_lon_inc() const { return *lon_inc_; }
  /// The incoming-angle latitude grid.
  math::Vector<double> get_lat_inc() const { return *lat_inc_; }
  //// The scattering-angle longitude grid.
  math::Vector<double> get_lon_scat() const { return *lon_scat_; }
  //// The scattering-angle latitude grid.
  math::Vector<double> get_lat_scat() const { return *lat_scat_; }

  /// Deep copy of the scattering data.
  ScatteringDataFieldGridded copy() const {
    auto data_new = std::make_shared<DataTensor>(*data_);
    return ScatteringDataFieldGridded(f_grid_,
                                      t_grid_,
                                      lon_inc_,
                                      lat_inc_,
                                      lon_scat_,
                                      lat_scat_,
                                      data_new);
  }

  /** Set scattering data for given frequency and temperature index.
   *
   * This function copies the data from the given scattering data field
   * into the sub-tensor of this objects' data tensor identified by
   * the given frequency and temperature indices. The data is automatically
   * regridded to the scattering data grids of this object.
   *
   * This function is useful to combine scattering data at different
   * temperatures and frequencies that have different scattering grids.
   *
   * @frequency_index The index along the frequency dimension
   * @temperature_index The index along the temperature dimension.
   */
  void set_data(math::Index frequency_index,
                math::Index temperature_index,
                const ScatteringDataFieldGridded &other) {
    using Regridder = RegularRegridder<Scalar, 2, 3, 4, 5>;

    auto f_grid_other = other.f_grid_;
    auto lon_inc_other = other.lon_inc_;
    auto lat_inc_other = other.lat_inc_;
    auto lon_scat_other = other.lon_scat_;
    auto lat_scat_other = other.lat_scat_;

    auto regridder = Regridder(
        {*lon_inc_other, *lat_inc_other, *lon_scat_other, *lat_scat_other},
        {*lon_inc_, *lat_inc_, *lon_scat_, *lat_scat_});
    auto regridded = regridder.regrid(*other.data_);

    std::array<math::Index, 2> data_index = {frequency_index,
                                              temperature_index};
    std::array<math::Index, 2> input_index = {0, 0};
    math::tensor_index(*data_, data_index) =
        math::tensor_index(regridded, input_index);
  }

  Vector interpolate(Scalar frequency,
                     Scalar temperature,
                     Scalar lon_inc,
                     Scalar lat_inc,
                     Scalar lon_scat,
                     Scalar lat_scat) const {
    using Regridder = RegularRegridder<Scalar, 0, 1, 2, 3, 4, 5>;
    auto f_grid = std::make_shared<Vector>(Vector::Constant(1, frequency));
    auto t_grid = std::make_shared<Vector>(Vector::Constant(1, temperature));
    auto lon_inc_grid = std::make_shared<Vector>(Vector::Constant(1, lon_inc));
    auto lat_inc_grid = std::make_shared<Vector>(Vector::Constant(1, lat_inc));
    auto lon_scat_grid = std::make_shared<Vector>(Vector::Constant(1, lon_scat));
    auto lat_scat_grid = std::make_shared<IrregularLatitudeGrid<Scalar>>(
        Vector::Constant(1, lat_scat));
    Regridder regridder(
        {*f_grid_, *t_grid_, *lon_inc_, *lat_inc_, *lon_scat_, *lat_scat_},
        {*f_grid,
         *t_grid,
         *lon_inc_grid,
         *lat_inc_grid,
         *lon_scat_grid,
         *lat_scat_grid});
    auto data_interp = regridder.regrid(*data_);
    auto n_stokes = data_interp.dimension(6);
    math::Vector<Scalar> results{n_stokes};
    for (decltype(n_stokes) i = 0; i < n_stokes; ++i) {
        results(i) = data_interp.coeff({0, 0, 0, 0, 0, 0, i});
    }
    return results;
  }

  // pxx :: hide
  /** Linear interpolation along frequency dimension.
   * @param frequencies The frequency grid to which to interpolate the data
   * @return New scattering data field with the given frequencies as
   * frequency grid.
   */
  ScatteringDataFieldGridded interpolate_frequency(
      VectorPtr frequencies) const {
    using Regridder = RegularRegridder<Scalar, 0>;
    Regridder regridder({*f_grid_}, {*frequencies});
    auto dimensions_new = data_->dimensions();
    auto data_interp = regridder.regrid(*data_);
    dimensions_new[0] = frequencies->size();
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldGridded(frequencies,
                                      t_grid_,
                                      lon_inc_,
                                      lat_inc_,
                                      lon_scat_,
                                      lat_scat_,
                                      data_new);
  }

  /** Linear interpolation along frequency dimension.
   * @param frequencies The frequency grid to which to interpolate the data
   * @return New scattering data field with the given frequencies as
   * frequency grid.
   */
  ScatteringDataFieldGridded interpolate_frequency(
      const Vector &frequencies) const {
    return interpolate_frequency(std::make_shared<Vector>(frequencies));
  }

  // pxx :: hide
  /** Linear interpolation along temperature dimension.
   * @param temperature The temperature grid to which to interpolate the data
   * @return New scattering data field with the given temperatures as
   * temperature grid.
   */
  ScatteringDataFieldGridded interpolate_temperature(
      VectorPtr temperatures,
      bool extrapolate=false) const {
    using Regridder = RegularRegridder<Scalar, 1>;
    Regridder regridder({*t_grid_}, {*temperatures}, extrapolate);
    auto dimensions_new = data_->dimensions();
    auto data_interp = regridder.regrid(*data_);
    dimensions_new[1] = temperatures->size();
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldGridded(f_grid_,
                                      temperatures,
                                      lon_inc_,
                                      lat_inc_,
                                      lon_scat_,
                                      lat_scat_,
                                      data_new);
  }

  /** Linear interpolation along temperature dimension.
   * @param temperature The temperature grid to which to interpolate the data
   * @return New scattering data field with the given temperatures as
   * temperature grid.
   */
  ScatteringDataFieldGridded interpolate_temperature(
      const Vector &temperatures,
      bool extrapolate=false) const {
    return interpolate_temperature(std::make_shared<Vector>(temperatures), extrapolate);
  }

  // pxx :: hide
  /** Linear interpolation along angles.
   * @param lon_inc_new The incoming-angle longitude grid to which to interpolate
   * the data.
   * @param lat_inc_new The incoming-angle latitude grid to which to interpolate
   * the data.
   * @param lon_scat_new The scattering-angle longitude grid to which to interpolate
   * the data.
   * @param lat_scat_new The scattering-angle longitude grid to which to interpolate
   * the data.
   * @return New scattering data field with the data interpolated
   */
  ScatteringDataFieldGridded interpolate_angles(VectorPtr lon_inc_new,
                                                VectorPtr lat_inc_new,
                                                VectorPtr lon_scat_new,
                                                LatitudeGridPtr lat_scat_new) const {
    using Regridder = RegularRegridder<Scalar, 2, 3, 4, 5>;
    Regridder regridder(
        {*lon_inc_, *lat_inc_, *lon_scat_, *lat_scat_},
        {*lon_inc_new, *lat_inc_new, *lon_scat_new, *lat_scat_new});
    auto data_interp = regridder.regrid(*data_);
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldGridded(f_grid_,
                                      t_grid_,
                                      lon_inc_new,
                                      lat_inc_new,
                                      lon_scat_new,
                                      lat_scat_new,
                                      data_new);
  }

  /** Interpolate angular grids.
   * @param lon_inc_new The incoming-angle longitude grid to which to interpolate
   * the data
   * @param lat_inc_new The incoming-angle latitude angle grid to which to
   * interpolate the data
   * @param lon_scat_new The scattering-angle longitude grid to which to
   * interpolate the data.
   * @param lat_scat_new The scattering-angle latitude grid to which to
   * interpolate the data
   * @return New scattering data field with the given angles as angular grids.
   */
  ScatteringDataFieldGridded interpolate_angles(Vector lon_inc_new,
                                                Vector lat_inc_new,
                                                Vector lon_scat_new,
                                                Vector lat_scat_new) const {
    return interpolate_angles(std::make_shared<const Vector>(lon_inc_new),
                              std::make_shared<const Vector>(lat_inc_new),
                              std::make_shared<const Vector>(lon_scat_new),
                              std::make_shared<const IrregularLatitudeGrid<Scalar>>(lat_scat_new));
  }

  // pxx :: hide
  /** Downsample scattering angles by averaging.
   *
   * @param lon_scat_new The scattering-angle longitude grid to which to
   * downsample the data.
   * @param lat_scat_new The scattering-angle latitude grid to which to
   * interpolate the data
   * @param Whether or not to apply downsampling also to latitude angles
   * or only to longitude angles.
   * @return New scattering data field with the given scattering angles.
   */
  ScatteringDataFieldGridded downsample_scattering_angles(VectorPtr lon_scat_new,
                                                          LatitudeGridPtr lat_scat_new,
                                                          bool interpolate_latitudes=true) const {
      auto data_downsampled = downsample_dimension<4>(*data_, *lon_scat_, *lon_scat_new, 0.0, 2.0 * M_PI);
      Vector colatitudes = -lat_scat_->array().cos();
      Vector colatitudes_new = -lat_scat_new->array().cos();

      if (interpolate_latitudes) {
          using Regridder = RegularRegridder<Scalar, 5>;
          Regridder regridder({*lat_scat_}, {*lat_scat_new}, false);
          data_downsampled = regridder.regrid(data_downsampled);
      } else {
          data_downsampled = downsample_dimension<5>(data_downsampled, colatitudes, colatitudes_new, -1.0, 1.0);
      }

      return ScatteringDataFieldGridded(f_grid_,
                                        t_grid_,
                                        lon_inc_,
                                        lat_inc_,
                                        lon_scat_new,
                                        lat_scat_new,
                                        std::make_shared<DataTensor>(data_downsampled));
  }

  /** Downsample scattering angles by averaging.
   *
   * @param lon_scat_new The scattering-angle longitude grid to which to
   * downsample the data.
   * @param lat_scat_new The scattering-angle latitude grid to which to
   * interpolate the data
   * @return New scattering data field with the given scattering angles.
   */
  ScatteringDataFieldGridded downsample_scattering_angles(Vector lon_scat_new,
                                                          Vector lat_scat_new) const {
    return downsample_scattering_angles(std::make_shared<Vector>(lon_scat_new),
                                        std::make_shared<IrregularLatitudeGrid<Scalar>>(lat_scat_new));
  }

  // pxx :: hide
  /** Downsample scattering longitudes angles by averaging.
   *
   * @param lon_scat_new The scattering-angle longitude grid to which to
   * downsample the data.
   * @return New scattering data field with the data along the longitude-component
   * of the scattering angle downsampled to the given grid.
   */
  ScatteringDataFieldGridded downsample_lon_scat(VectorPtr lon_scat_new) const {
      auto data_downsampled = downsample_dimension<4>(*data_, *lon_scat_, *lon_scat_new, 0.0, 2.0 * M_PI);
      return ScatteringDataFieldGridded(f_grid_,
                                        t_grid_,
                                        lon_inc_,
                                        lat_inc_,
                                        lon_scat_new,
                                        lat_scat_,
                                        std::make_shared<DataTensor>(data_downsampled));
  }

  // pxx :: hide
  /** Regrid data to new grids.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @lon_inc The incoming-angle longitude grid.
   * @lat_inc The incoming-angle latitude grid.
   * @lon_scat The scattering-angle longitude grid.
   * @lat_scat The scattering-angle latitude grid.
   * @return A new ScatteringDataFieldGridded with the given grids.
   */
  ScatteringDataFieldGridded regrid(VectorPtr f_grid,
                                    VectorPtr t_grid,
                                    VectorPtr lon_inc,
                                    VectorPtr lat_inc,
                                    VectorPtr lon_scat,
                                    LatitudeGridPtr lat_scat) const {
    using Regridder = RegularRegridder<Scalar, 0, 1, 2, 3, 4, 5>;
    Regridder regridder(
        {*f_grid_, *t_grid_, *lon_inc_, *lat_inc_, *lon_scat_, *lat_scat_},
        {*f_grid, *t_grid, *lon_inc, *lat_inc, *lon_scat, *lat_scat});
    auto data_interp = regridder.regrid(*data_);
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldGridded(f_grid,
                                      t_grid,
                                      lon_inc,
                                      lat_inc,
                                      lon_scat,
                                      lat_scat,
                                      data_new);
  }

  /** Calculate scattering-angle integral.
   *
   * @return Rank-5 tensor containing the scattering-angle integrals of
   * the data tensor.
   */
  math::Tensor<Scalar, 5> integrate_scattering_angles() {
    math::IndexArray<5> dimensions = {n_freqs_,
                                       n_temps_,
                                       n_lon_inc_,
                                       n_lat_inc_,
                                       data_->dimension(6)};
    auto result = math::Tensor<Scalar, 5>(dimensions);
    for (auto i = math::DimensionCounter<5>{dimensions}; i; ++i) {
      auto matrix = math::get_submatrix<4, 5>(*data_, i.coordinates);
      result.coeffRef(i.coordinates) =
          integrate_angles<Scalar>(matrix, *lon_scat_, *lat_scat_);
    }
    return result;
  }

  /** Normalize w.r.t. scattering-angle integral.
   *
   * Normalization is performed in place, i.e. the object
   * is changed.
   *
   * @param value The value to normalize the integrals to.
   */
  void normalize(Scalar value) {
    auto integrals = integrate_scattering_angles();
    math::IndexArray<4> dimensions = {n_freqs_,
                                       n_temps_,
                                       n_lon_inc_,
                                       n_lat_inc_};
    for (auto i = math::DimensionCounter<4>{dimensions}; i; ++i) {
      for (Index j = 0; j < data_->dimension(6); ++j) {
        auto matrix_coords = concat<Index, 4, 1>(i.coordinates, {j});
        auto matrix = math::get_submatrix<4, 5>(*data_, matrix_coords);
        auto integral = integrals(concat<Index, 4, 1>(i.coordinates, {0}));
        if (integral != 0.0) {
            matrix *= value / integral;
        }
      }
    }
  }

  /** Accumulate scattering data into this object.
   *
   * Regrids the given scattering data field and accumulates its interpolated
   * data tensor into this object's data tensor.
   *
   * @param other The ScatteringDataField to accumulate into this.
   * @return Reference to this object.
   */
  ScatteringDataFieldGridded operator+=(
      const ScatteringDataFieldGridded &other) {
    auto regridded = other.regrid(f_grid_,
                                  t_grid_,
                                  lon_inc_,
                                  lat_inc_,
                                  lon_scat_,
                                  lat_scat_);
    *data_ += regridded.get_data();
    return *this;
  }

  /** Addition of scattering data.
   *
   * Creates a new scattering data object with the same grids as the
   * left-hand operand containing the sum of the scattering data in both
   * operands.
   *
   * @param other The right-hand summand.
   * @return The sum of the two scatttering data objects.
   */
  ScatteringDataFieldGridded operator+(
      const ScatteringDataFieldGridded &other) const {
    auto result = copy();
    result += other;
    return result;
  }

  /** In-place scaling of scattering data.
   *
   * @param c The scaling factor.
   * @return Reference to this object.
   */
  ScatteringDataFieldGridded operator*=(Scalar c) {
    (*data_) = c * (*data_);
    return *this;
  }

  /** Scale scattering data.
   *
   * @param c The scaling factor.
   * @return A new object containing the scaled scattering data.
   */
  ScatteringDataFieldGridded operator*(Scalar c) const {
    auto result = copy();
    result *= c;
    return result;
  }

  /** Set the number of scattering coefficients.
   *
   * De- or increases the number of scattering coefficients that are stored.
   * If the number of scattering coefficients is increased, the new elements are
   * set to 0.
   *
   * @param n The number of scattering coefficients to change the data to have.
   */
  void set_number_of_scattering_coeffs(Index n) {
    Index current_stokes_dim = data_->dimension(6);
    if (current_stokes_dim == n) {
      return;
    }
    auto new_dimensions = data_->dimensions();
    new_dimensions[6] = n;
    DataTensorPtr data_new = std::make_shared<DataTensor>(new_dimensions);
    math::copy(*data_new, *data_);
    data_ = data_new;
  }

  // pxx :: hide
  /** Convert gridded data to spectral format.
   * @param sht SHT instance to use for the transformation.
   * @return The scattering data field transformed to spectral format.
   */
  ScatteringDataFieldSpectral<Scalar> to_spectral(
      std::shared_ptr<sht::SHT> sht) const;

  /** Convert gridded data to spectral format.
   * @param l_max The maximum degree l to use in the SH expansion.
   * @param m_max The maximum order m to use in the SH expansion.
   * @return The scattering data field transformed to spectral format.
   */
  ScatteringDataFieldSpectral<Scalar> to_spectral(Index l_max,
                                                  Index m_max) const {
    std::shared_ptr<sht::SHT> sht =
        std::make_shared<sht::SHT>(l_max, m_max, n_lon_scat_, n_lat_scat_);
    return to_spectral(sht);
  }

  /** Convert gridded data to spectral format.
   * @param l_max The maximum degree l to use in the SH expansion.
   * @param m_max The maximum order m to use in the SH expansion.
   * @return The scattering data field transformed to spectral format.
   */
  ScatteringDataFieldSpectral<Scalar> to_spectral(Index l_max,
                                                  Index m_max,
                                                  Index n_lon,
                                                  Index n_lat) const {
      std::shared_ptr<sht::SHT> sht =
          std::make_shared<sht::SHT>(l_max,
                                     m_max,
                                     n_lon,
                                     n_lat);
      return to_spectral(sht);
  }

  /** Convert gridded data to spectral format.
   *
   * This version uses the highest possible values for the maximum order
   * and degree that fulfill the anti-aliasing conditions.
   * @return The scattering data field transformed to spectral format.
   */
  ScatteringDataFieldSpectral<Scalar> to_spectral() const {
    auto sht_params = get_sht_scat_params();
    return to_spectral(std::get<0>(sht_params), std::get<1>(sht_params));
  }

  /// The data tensor containing the scattering data.
  const DataTensor &get_data() const { return *data_; }

  /** Phase function maximum.
   *
   * Maximum of the phase function calculated across the dimensions
   * of incoming angles.
   */
      math::Tensor<Scalar, 4> get_phase_function_max_inc() const {
      Eigen::array<int, 2> dims({2, 3});
      return data_-> template chip<6>(0).maximum(dims);
  }

  /** Phase function maximum.
   *
   * Maximum of the phase function calculated across the dimensions
   * of scattering angles.
   */
  math::Tensor<Scalar, 4> get_phase_function_max_scat() const {
      std::array<int, 2> dims({4, 5});
      return data_-> template chip<6>(0).maximum(dims);
  }

 protected:

  VectorPtr f_grid_;
  VectorPtr t_grid_;
  VectorPtr lon_inc_;
  VectorPtr lat_inc_;
  VectorPtr lon_scat_;
  LatitudeGridPtr lat_scat_;

  ConstVectorMap f_grid_map_;
  ConstVectorMap t_grid_map_;
  ConstVectorMap lon_inc_map_;
  ConstVectorMap lat_inc_map_;
  ConstVectorMap lon_scat_map_;
  ConstVectorMap lat_scat_map_;

  DataTensorPtr data_;
};

////////////////////////////////////////////////////////////////////////////////
// Spectral format
////////////////////////////////////////////////////////////////////////////////

// pxx :: export
// pxx :: instance(["double"])
/** Scattering data in spectral format.
 *
 * Uses a spherical-harmonics expansion to represent the scattering-angle
 * dependency. The data tensor is therefore complex with the axes corresponding
 * to the following quantities:
 * 1. Frequency
 * 2. Temperature
 * 3. Incoming-angle longitude
 * 4. Incoming-angle latitude
 * 5. Scattering-angle SH coefficients.
 * 6. Stokes dimension of scattering data.
 */
template <typename Scalar_>
class ScatteringDataFieldSpectral : public ScatteringDataFieldBase {
 public:
  using ScatteringDataFieldBase::n_freqs_;
  using ScatteringDataFieldBase::n_lat_inc_;
  using ScatteringDataFieldBase::n_lat_scat_;
  using ScatteringDataFieldBase::n_lon_inc_;
  using ScatteringDataFieldBase::n_lon_scat_;
  using ScatteringDataFieldBase::n_temps_;
  using ScatteringDataFieldBase::type_;
  using ScatteringDataFieldBase::get_n_lon_inc;
  using ScatteringDataFieldBase::get_n_lat_inc;
  using ScatteringDataFieldBase::get_n_lon_scat;
  using ScatteringDataFieldBase::get_n_lat_scat;

  using Scalar = Scalar_;
  using Coefficient = std::complex<Scalar>;
  using Vector = math::Vector<Scalar>;
  using VectorMap = math::VectorMap<Scalar>;
  using VectorPtr = std::shared_ptr<const math::Vector<Scalar>>;
  using ConstVectorMap = math::ConstVectorMap<Scalar>;
  using Matrix = math::Matrix<Scalar>;
  using MatrixMap = math::MatrixMap<Scalar>;
  using ConstMatrixMap = math::ConstMatrixMap<Scalar>;
  using ShtPtr = std::shared_ptr<sht::SHT>;

  template <math::Index rank>
  using CmplxTensor = math::Tensor<std::complex<Scalar>, rank>;
  template <math::Index rank>
  using CmplxTensorMap = math::TensorMap<std::complex<Scalar>, rank>;
  template <math::Index rank>
  using ConstCmplxTensorMap = math::ConstTensorMap<std::complex<Scalar>, rank>;
  using DataTensor = math::Tensor<std::complex<Scalar>, 6>;
  using DataTensorPtr = std::shared_ptr<DataTensor>;

  static constexpr Index coeff_dim = 5;
  static constexpr Index rank = 6;

  // pxx :: hide
  /** Create spectral scattering data field.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The longitude grid for the incoming angles.
   * @param lat_inc The latitude grid for the incoming angles.
   * @param sht_scat The SH transform used to expand the scattering-angle
   * dependency.
   * @data The scattering data.
   */
  ScatteringDataFieldSpectral(VectorPtr f_grid,
                              VectorPtr t_grid,
                              VectorPtr lon_inc,
                              VectorPtr lat_inc,
                              ShtPtr sht_scat,
                              DataTensorPtr data)
      : ScatteringDataFieldBase(f_grid->size(),
                                t_grid->size(),
                                lon_inc->size(),
                                lat_inc->size(),
                                sht_scat->get_n_longitudes(),
                                sht_scat->get_n_latitudes()),
        f_grid_(f_grid),
        t_grid_(t_grid),
        lon_inc_(lon_inc),
        lat_inc_(lat_inc),
        sht_scat_(sht_scat),
        f_grid_map_(f_grid->data(), n_freqs_),
        t_grid_map_(t_grid->data(), n_temps_),
        lon_inc_map_(lon_inc->data(), n_lon_inc_),
        lat_inc_map_(lat_inc->data(), n_lat_inc_),
        data_(data) {}

  /** Create spectral scattering data field.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The longitude grid for the incoming angles.
   * @param lat_inc The latitude grid for the incoming angles.
   * @param sht_scat The SH transform used to expand the scattering-angle
   * dependency.
   * @data The scattering data.
   */
  ScatteringDataFieldSpectral(const Vector &f_grid,
                              const Vector &t_grid,
                              const Vector &lon_inc,
                              const Vector &lat_inc,
                              const sht::SHT &sht_scat,
                              const DataTensor &data)
      : ScatteringDataFieldBase(f_grid.size(),
                                t_grid.size(),
                                lon_inc.size(),
                                lat_inc.size(),
                                sht_scat.get_n_longitudes(),
                                sht_scat.get_n_latitudes()),
        f_grid_(std::make_shared<Vector>(f_grid)),
        t_grid_(std::make_shared<Vector>(t_grid)),
        lon_inc_(std::make_shared<Vector>(lon_inc)),
        lat_inc_(std::make_shared<Vector>(lat_inc)),
        sht_scat_(std::make_shared<sht::SHT>(sht_scat)),
        f_grid_map_(f_grid_->data(), n_freqs_),
        t_grid_map_(t_grid_->data(), n_temps_),
        lon_inc_map_(lon_inc_->data(), n_freqs_),
        lat_inc_map_(lat_inc_->data(), n_temps_),
        data_(std::make_shared<DataTensor>(data)) {}

  /** Create empty scattering data field.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param lon_inc The longitude grid for the incoming angles.
   * @param lat_inc The latitude grid for the incoming angles.
   * @param sht_scat The SH transform used to expand the scattering-angle
   * dependency.
   * @data The scattering data.
   */
  ScatteringDataFieldSpectral(const Vector &f_grid,
                              const Vector &t_grid,
                              const Vector &lon_inc,
                              const Vector &lat_inc,
                              const sht::SHT &sht_scat,
                              Index n_elements)
      : ScatteringDataFieldBase(f_grid.size(),
                                t_grid.size(),
                                lon_inc.size(),
                                lat_inc.size(),
                                sht_scat.get_n_longitudes(),
                                sht_scat.get_n_latitudes()),
        f_grid_(std::make_shared<Vector>(f_grid)),
        t_grid_(std::make_shared<Vector>(t_grid)),
        lon_inc_(std::make_shared<Vector>(lon_inc)),
        lat_inc_(std::make_shared<Vector>(lat_inc)),
        sht_scat_(std::make_shared<sht::SHT>(sht_scat)),
        f_grid_map_(f_grid_->data(), n_freqs_),
        t_grid_map_(t_grid_->data(), n_temps_),
        lon_inc_map_(lon_inc_->data(), n_freqs_),
        lat_inc_map_(lat_inc_->data(), n_temps_),
        data_(std::make_shared<DataTensor>(
            std::array<Index, 6>{n_freqs_,
                                 n_temps_,
                                 n_lon_inc_,
                                 n_lat_inc_,
                                 sht_scat.get_n_spectral_coeffs(),
                                 n_elements})) {}

  /// Shallow copy of the ScatteringDataField.
  ScatteringDataFieldSpectral(const ScatteringDataFieldSpectral &) = default;

  /// Deep copy of the scattering data.
  ScatteringDataFieldSpectral copy() const {
    auto data_new = std::make_shared<DataTensor>(*data_);
    return ScatteringDataFieldSpectral(f_grid_,
                                       t_grid_,
                                       lon_inc_,
                                       lat_inc_,
                                       sht_scat_,
                                       data_new);
  }

  /// Enum representing the data format.
  constexpr DataFormat get_data_format() const { return DataFormat::Spectral; }

  /// The number of scattering-data coefficients.
  Index get_n_coeffs() const { return data_->dimension(5); }
  /// The frequency grid.
  /// Parameters of SHT transformation used to transform
  /// scattering angle.
  std::array<Index, 4> get_sht_scat_params() const {
      return {sht_scat_->get_l_max(),
              sht_scat_->get_m_max(),
              sht_scat_->get_n_longitudes(),
              sht_scat_->get_n_latitudes()};
  }
  /// Largest SHT parameters satisfying shtns aliasing requirements for
  /// incoming angle.
  std::array<Index, 4> get_sht_inc_params() const {
      return sht::SHT::get_params(n_lon_inc_, n_lat_inc_);
  }
  const math::Vector<double>& get_f_grid() const { return *f_grid_; }
  /// The temperature grid.
  const math::Vector<double>& get_t_grid() const { return *t_grid_; }
  /// The incoming-angle longitude grid.
  math::Vector<double> get_lon_inc() const { return *lon_inc_; }
  /// The incoming-angle latitude grid.
  math::Vector<double> get_lat_inc() const { return *lat_inc_; }
  /// The scattering-angle longitude grid.
  math::Vector<double> get_lon_scat() const {
    return sht_scat_->get_longitude_grid();
  }
  /// The scattering-angle latitude grid.
  math::Vector<double> get_lat_scat() const {
    return sht_scat_->get_latitude_grid();
  }

  /// The SHT object used to transform the data along the scattering
  /// angle.
  sht::SHT &get_sht_scat() const { return *sht_scat_; }

  /// The raw scattering data.
  const DataTensor &get_data() const { return *data_; }

  /** Set scattering data for given frequency and temperature index.
   *
   * This function copies the data from the given scattering data field
   * into the sub-tensor of this objects' data tensor identified by
   * the given frequency and temperature indices. The data is automatically
   * regridded to the scattering data grids of this object.
   *
   * This function is useful to combine scattering data at different
   * temperatures and frequencies that have different scattering grids.
   *
   * @frequency_index The index along the frequency dimension
   * @temperature_index The index along the temperature dimension.
   */
  void set_data(math::Index frequency_index,
                math::Index temperature_index,
                const ScatteringDataFieldSpectral &other) {
    using Regridder = RegularRegridder<Scalar, 2, 3>;

    auto lon_inc_other = other.lon_inc_;
    auto lat_inc_other = other.lat_inc_;
    auto regridder =
        Regridder({*lon_inc_other, *lat_inc_other}, {*lon_inc_, *lat_inc_});
    auto regridded = regridder.regrid(*other.data_);

    std::array<math::Index, 2> data_index = {frequency_index,
                                          temperature_index};
    std::array<math::Index, 2> input_index = {0, 0};

    math::IndexArray<3> dimensions_loop = {n_lon_inc_,
                                            n_lat_inc_,
                                            data_->dimension(5)};
    auto data_map = math::tensor_index(*data_, data_index);
    auto other_data_map = math::tensor_index(regridded, input_index);
    for (math::DimensionCounter<3> i{dimensions_loop}; i; ++i) {
      auto result = math::get_subvector<2>(data_map, i.coordinates);
      auto in_l = math::get_subvector<2>(data_map, i.coordinates);
      auto in_r = math::get_subvector<2>(other_data_map, i.coordinates);
      result = sht::SHT::add_coeffs(*sht_scat_, in_l, *other.sht_scat_, in_r);
    }
  }

  Vector interpolate(Scalar frequency,
                     Scalar temperature,
                     Scalar lon_inc,
                     Scalar lat_inc,
                     Scalar lon_scat,
                     Scalar lat_scat) const {
    using Regridder = RegularRegridder<Scalar, 0, 1, 2, 3>;
    auto f_grid = std::make_shared<Vector>(Vector::Constant(1, frequency));
    auto t_grid = std::make_shared<Vector>(Vector::Constant(1, temperature));
    auto lon_inc_grid = std::make_shared<Vector>(Vector::Constant(1, lon_inc));
    auto lat_inc_grid = std::make_shared<Vector>(Vector::Constant(1, lat_inc));
    auto lon_scat_grid = std::make_shared<Vector>(Vector::Constant(1, lon_scat));
    auto lat_scat_grid = std::make_shared<IrregularLatitudeGrid<Scalar>>(Vector::Constant(1, lat_scat)) ;
    Regridder regridder({*f_grid_,
                         *t_grid_,
                         *lon_inc_,
                         *lat_inc_},
                        {*f_grid,
                         *t_grid,
                         *lon_inc_grid,
                         *lat_inc_grid});
    auto data_interp = regridder.regrid(*data_);

    auto n_stokes = data_->dimension(5);
    math::IndexArray<5> dimensions_loop = {1, 1, 1, 1, n_stokes};

    using DataTensor = math::Tensor<Scalar, 7>;
    Vector results{n_stokes};
    for (math::DimensionCounter<5> i{dimensions_loop}; i; ++i) {
        auto synthesized = sht_scat_->evaluate(math::get_subvector<4>(data_interp, i.coordinates),
                                               lon_scat,
                                               lat_scat);
        results[i.coordinates[4]] = synthesized;
    }
    return results;
  }

  /** Interpolate data along frequency dimension.
   * @param frequencies The frequencies to which to interpolate the data.
   * @return New ScatteringDataFieldSpectral with the data interpolated
   * to the given frequencies.
   */
  // pxx :: hide
  ScatteringDataFieldSpectral interpolate_frequency(
      VectorPtr frequencies) const {
    using Regridder = RegularRegridder<Scalar, 0>;
    Regridder regridder({*f_grid_}, {*frequencies});
    auto dimensions_new = data_->dimensions();
    auto data_interp = regridder.regrid(*data_);
    dimensions_new[0] = frequencies->size();
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldSpectral(frequencies,
                                       t_grid_,
                                       lon_inc_,
                                       lat_inc_,
                                       sht_scat_,
                                       data_new);
  }

  /** Interpolate data along frequency dimension.
   * @param frequencies The frequencies to which to interpolate the data.
   * @return New ScatteringDataFieldSpectral with the data interpolated
   * to the given frequencies.
   */
  ScatteringDataFieldSpectral interpolate_frequency(
      const Vector &frequencies) const {
    return interpolate_frequency(std::make_shared<Vector>(frequencies));
  }

  // pxx :: hide
  /** Interpolate data along temperature dimension.
   * @param temperatures The temperatures to which to interpolate the data.
   * @return New ScatteringDataFieldSpectral with the data interpolated
   * to the given temperatures.
   */
  ScatteringDataFieldSpectral interpolate_temperature(
      VectorPtr temperatures,
      bool extrapolate=false) const {
    using Regridder = RegularRegridder<Scalar, 1>;
    Regridder regridder({*t_grid_}, {*temperatures}, extrapolate);
    auto dimensions_new = data_->dimensions();
    auto data_interp = regridder.regrid(*data_);
    dimensions_new[1] = temperatures->size();
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldSpectral(f_grid_,
                                       temperatures,
                                       lon_inc_,
                                       lat_inc_,
                                       sht_scat_,
                                       data_new);
  }

  /** Interpolate data along temperature dimension.
   * @param temperatures The temperatures to which to interpolate the data.
   * @return New ScatteringDataFieldSpectral with the data interpolated
   * to the given temperatures.
   */
  ScatteringDataFieldSpectral interpolate_temperature(
      const Vector &temperatures,
      bool extrapolate=false) const {
      return interpolate_temperature(std::make_shared<Vector>(temperatures),
                                     extrapolate);
  }

  // pxx :: hide
  /** Interpolate data along incoming angles.
   * @param lon_inc_new The incoming-angle longitudes to which to interpolate the data.
   * @param lat_inc_new The incoming-angle latitudes to which to interpolate the data.
   * @return New ScatteringDataFieldSpectral with the data interpolated
   * to the given incoming angles.
   */
  ScatteringDataFieldSpectral interpolate_angles(VectorPtr lon_inc_new,
                                                 VectorPtr lat_inc_new) const {
    using Regridder = RegularRegridder<Scalar, 2, 3>;
    Regridder regridder({*lon_inc_, *lat_inc_}, {*lon_inc_new, *lat_inc_new});
    auto dimensions_new = data_->dimensions();
    dimensions_new[2] = lon_inc_new->size();
    dimensions_new[3] = lat_inc_new->size();
    auto data_new = std::make_shared<DataTensor>(DataTensor(dimensions_new));
    regridder.regrid(*data_new, *data_);
    return ScatteringDataFieldSpectral(f_grid_,
                                       t_grid_,
                                       lon_inc_new,
                                       lat_inc_new,
                                       sht_scat_,
                                       data_new);
  }

  /** Interpolate data along incoming angles.
   * @param lon_inc_new The incoming-angle longitudes to which to interpolate the data.
   * @param lat_inc_new The incoming-angle latitudes to which to interpolate the data.
   * @return New ScatteringDataFieldSpectral with the data interpolated
   * to the given incoming angles.
   */
  ScatteringDataFieldSpectral interpolate_angles(Vector lon_inc_new,
                                                 Vector lat_inc_new) const {
    return interpolate_angles(std::make_shared<Vector>(lon_inc_new),
                              std::make_shared<Vector>(lat_inc_new));
  }

  /** Regrid data to new grids.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @lon_inc The incoming azimuth angle
   * @lat_inc The incoming zenith angle
   * @return A new ScatteringDataFieldSpectral with the given grids.
   */
  // pxx :: hide
  ScatteringDataFieldSpectral regrid(VectorPtr f_grid,
                                     VectorPtr t_grid,
                                     VectorPtr lon_inc,
                                     VectorPtr lat_inc) const {
    using Regridder = RegularRegridder<Scalar, 0, 1, 2, 3>;
    Regridder regridder({*f_grid_, *t_grid_, *lon_inc_, *lat_inc_},
                        {*f_grid, *t_grid, *lon_inc, *lat_inc});
    auto data_interp = regridder.regrid(*data_);
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldSpectral(f_grid,
                                       t_grid,
                                       lon_inc,
                                       lat_inc,
                                       sht_scat_,
                                       data_new);
  }

  /** Calculate scattering-angle integral.
   *
   * @return Rank-5 tensor containing the scattering-angle integrals of
   * the data tensor.
   */
  math::Tensor<Scalar, 5> integrate_scattering_angles() const {
      math::Tensor<std::complex<Scalar>, 5> result = data_->template chip<4>(0);
      return result.real() * sqrt(4.0 * M_PI);
  }

  /** Normalize w.r.t. scattering-angle integral.
   *
   * Normalization is performed in place, i.e. the object
   * is changed.
   *
   * @param value The value to normalize the integrals to.
   */
  void normalize(Scalar value) {
    auto integrals = integrate_scattering_angles();
    math::IndexArray<4> dimensions = {n_freqs_,
                                       n_temps_,
                                       n_lon_inc_,
                                       n_lat_inc_};
    for (auto i = math::DimensionCounter<4>{dimensions}; i; ++i) {
      auto matrix = math::get_submatrix<4, 5>(*data_, i.coordinates);
      auto integral = integrals(concat<Index, 4, 1>(i.coordinates, {0}));
      if (integral != 0.0) {
          matrix *= value / integral;
      }
    }
  }

  /** Accumulate scattering data into this object.
   *
   * Regrids the given scattering data field and accumulates its interpolated
   * data tensor into this object's data tensor.
   *
   * @param other The ScatteringDataField to accumulate into this.
   * @return Reference to this object.
   */
  ScatteringDataFieldSpectral &operator+=(
      const ScatteringDataFieldSpectral &other) {
    auto regridded = other.regrid(f_grid_, t_grid_, lon_inc_, lat_inc_);
    math::IndexArray<5> dimensions_loop = {n_freqs_,
                                            n_temps_,
                                            n_lon_inc_,
                                            n_lat_inc_,
                                            data_->dimension(5)};
    for (math::DimensionCounter<5> i{dimensions_loop}; i; ++i) {
      auto result = math::get_subvector<4>(*data_, i.coordinates);
      auto in_l = math::get_subvector<4>(*data_, i.coordinates);
      auto in_r = math::get_subvector<4>(*regridded.data_, i.coordinates);
      result =
          sht::SHT::add_coeffs(*sht_scat_, in_l, *regridded.sht_scat_, in_r);
    }
    return *this;
  }

  /** Add scattering data fields.
   *
   * Regrids the given scattering data field to the grids of this object and
   * computes the sum of the two scattering data fields.
   *
   * @param other The ScatteringDataField to accumulate into this.
   * @return Reference to this object.
   */
  ScatteringDataFieldSpectral operator+(
      const ScatteringDataFieldSpectral &other) {
    auto result = copy();
    result += other;
    return result;
  }

  /** In-place scaling scattering data.
   *
   * @param c The scaling factor.
   * @return Reference to this object.
   */
  ScatteringDataFieldSpectral &operator*=(Scalar c) {
    (*data_) = c * (*data_);
    return *this;
  }

  /** Scale scattering data.
   *
   * @param c The scaling factor.
   * @return A new object containing the scaled scattering data.
   */
  ScatteringDataFieldSpectral operator*(Scalar c) const {
    auto result = copy();
    result *= c;
    return result;
  }

  /** Set the number of scattering coefficients.
   *
   * De- or increases the number of scattering coefficients that are stored.
   * If the number of scattering coefficients is increased, the new elements are
   * set to 0.
   *
   * @param n The number of scattering coefficients to change the data to have.
   */
  void set_number_of_scattering_coeffs(Index n) {
      Index current_stokes_dim = data_->dimension(5);
      if (current_stokes_dim == n) {
          return;
      }
      auto new_dimensions = data_->dimensions();
      new_dimensions[5] = n;
      DataTensorPtr data_new = std::make_shared<DataTensor>(new_dimensions);
      math::copy(*data_new, *data_);
      data_ = data_new;
  }

  /** Convert data to SHT representation with other parameters.
   *
   * @param sht_other SHT object representing SHT-representation to convert to.
   * @returns A new ScatteringDataFieldSpectral object containg the scattering
   * data of this object in the requested representation.
   */
  ScatteringDataFieldSpectral to_spectral(ShtPtr sht_other) const {
    auto new_dimensions = data_->dimensions();
    new_dimensions[4] = sht_other->get_n_spectral_coeffs();
    auto data_new_ =
        std::make_shared<DataTensor>(DataTensor(new_dimensions).setZero());
    auto result = ScatteringDataFieldSpectral(f_grid_,
                                              t_grid_,
                                              lon_inc_,
                                              lat_inc_,
                                              sht_other,
                                              data_new_);
    result += *this;
    return result;
  }

  /** Convert data to SHT representation with other parameters.
   *
   * @param l_max The l_max parameter of the SHT-representation to convert to.
   * @param m_max The m_max parameter of the SHT-representation to convert to.
   * @returns A new ScatteringDataFieldSpectral object containg the scattering
   * data of this object in the requested representation.
   */
  ScatteringDataFieldSpectral to_spectral(Index l_max, Index m_max) const {
    auto n_lat = sht_scat_->get_n_latitudes();
    auto n_lon = sht_scat_->get_n_longitudes();
    return to_spectral(std::make_shared<sht::SHT>(l_max, m_max, n_lon, n_lat));
  }

  /** Convert data to SHT representation with other parameters.
   *
   * @param l_max The l_max parameter of the SHT-representation to convert to.
   * @param m_max The m_max parameter of the SHT-representation to convert to.
   * @param n_lon Size of the scattering-angle longitude grid of the new
   * representation.
   * @param n_lat Size of the scattering-angle latitude of the new
   * representation.
   * @returns A new ScatteringDataFieldSpectral object containg the scattering
   * data of this object in the requested representation.
   */
  ScatteringDataFieldSpectral to_spectral(Index l_max,
                                          Index m_max,
                                          Index n_lon,
                                          Index n_lat) const {
      return to_spectral(std::make_shared<sht::SHT>(l_max, m_max, n_lon, n_lat));
  }

  /** Convert data to SHT representation with other parameters.
   *
   * @param l_max The l_max parameter of the SHT-representation to convert to.
   * @returns A new ScatteringDataFieldSpectral object containg the scattering
   * data of this object in the requested representation.
   */
  ScatteringDataFieldSpectral to_spectral(Index l_max) const {
    return to_spectral(l_max, l_max);
  }

  /** Convert data gridded representation.
   *
   * @returns A new ScatteringDataFieldSpectral object containg the scattering
   * data converted to gridded representation using the native scattering-angle
   * latitude * and longitude grids of the SHT transformation.
   */
  ScatteringDataFieldGridded<Scalar> to_gridded() const;

  ScatteringDataFieldGridded<Scalar> to_gridded(Index n_lon,
                                                Index n_lat) const {
      auto sht = std::make_shared<sht::SHT>(sht_scat_->get_l_max(),
                                            sht_scat_->get_m_max(),
                                            n_lon,
                                            n_lat);
      return to_spectral(sht).to_gridded();
  }

  /** Convert data to fully-spectral representation.
   *
   * @param sht SHT object to use to perform the incoming-angle transform.
   * @returns Representation of this scattering data converted to fully-spectral
   * format.
   */
  ScatteringDataFieldFullySpectral<Scalar> to_fully_spectral(ShtPtr sht) const;

  /** Convert data to fully-spectral representation.
   *
   * @param l_max The l_max parameter to use for the SHT transform of the incoming-angle
   * grid.
   * @param m_max The m_max parameter to use for the SHT transform of the incoming-angle
   * grid.
   * @returns Representation of this scattering data converted to fully-spectral
   * format.
   */
  ScatteringDataFieldFullySpectral<Scalar> to_fully_spectral(
      Index l_max,
      Index m_max) const {
    std::shared_ptr<sht::SHT> sht =
        std::make_shared<sht::SHT>(l_max, m_max, n_lon_inc_, n_lat_inc_);
    return to_fully_spectral(sht);
  }

  /** Convert data to fully-spectral representation.
   *
   * Automatically determines the largest l_max and m_max values that satisfy
   * the shtns aliasing requirements and uses the to transform the data
   * to fully-spectral format.
   *
   * @returns Representation of this scattering data converted to fully-spectral
   * format.
   */
  ScatteringDataFieldFullySpectral<Scalar> to_fully_spectral() {
    auto sht_params = get_sht_inc_params();
    return to_fully_spectral(std::get<0>(sht_params), std::get<1>(sht_params));
  }

  /** Phase function maximum.
   *
   * Maximum of the phase function calculated across the dimensions
   * of incoming angles.
   */
  math::Tensor<Scalar, 4> get_phase_function_max_inc() const {
      return to_gridded().get_phase_function_max_inc();
  }

  /** Phase function maximum.
   *
   * Maximum of the phase function calculated across the dimensions
   * of scattering angles.
   */
  math::Tensor<Scalar, 4> get_phase_function_max_scat() const {
      Eigen::array<int, 2> dims({4, 5});
      return to_gridded().get_phase_function_max_scat();
  }


 protected:

  VectorPtr f_grid_;
  VectorPtr t_grid_;
  VectorPtr lon_inc_;
  VectorPtr lat_inc_;
  VectorPtr lon_scat_;
  VectorPtr lat_scat_;
  ShtPtr sht_scat_;

  ConstVectorMap f_grid_map_;
  ConstVectorMap t_grid_map_;
  ConstVectorMap lon_inc_map_;
  ConstVectorMap lat_inc_map_;

  DataTensorPtr data_;
};

// pxx :: export
// pxx :: instance(["double"])
////////////////////////////////////////////////////////////////////////////////
// Fully-spectral format
////////////////////////////////////////////////////////////////////////////////
/** Fully-spectral scattering data.
 *
 * Represents scattering data fields whose dependency to both the incoming and
 * the scattering angles is represented using SHs.
 *
 */
/** Scattering data in fully-spectral format.
 *
 * Uses spherical-harmonics expansion to represent the scattering- as well
 * as the incoming-angle dependency. The data tensor is therefore complex with
 * the axes corresponding to the following quantities:
 * 1. Frequency
 * 2. Temperature
 * 3. Incoming-angle SH coefficients.
 * 4. Scattering-angle SH coefficients.
 * 5. Stokes dimension of scattering data.
 */
template <typename Scalar>
class ScatteringDataFieldFullySpectral : public ScatteringDataFieldBase {
 public:
  using ScatteringDataFieldBase::n_freqs_;
  using ScatteringDataFieldBase::n_lat_inc_;
  using ScatteringDataFieldBase::n_lat_scat_;
  using ScatteringDataFieldBase::n_lon_inc_;
  using ScatteringDataFieldBase::n_lon_scat_;
  using ScatteringDataFieldBase::n_temps_;
  using ScatteringDataFieldBase::type_;
  using ScatteringDataFieldBase::get_n_lon_inc;
  using ScatteringDataFieldBase::get_n_lat_inc;
  using ScatteringDataFieldBase::get_n_lon_scat;
  using ScatteringDataFieldBase::get_n_lat_scat;

  using Vector = math::Vector<Scalar>;
  using VectorMap = math::VectorMap<Scalar>;
  using VectorPtr = std::shared_ptr<const math::Vector<Scalar>>;
  using ConstVectorMap = math::ConstVectorMap<Scalar>;
  using Matrix = math::Matrix<Scalar>;
  using MatrixMap = math::MatrixMap<Scalar>;
  using ConstMatrixMap = math::ConstMatrixMap<Scalar>;
  using ShtPtr = std::shared_ptr<sht::SHT>;

  template <math::Index rank>
  using CmplxTensor = math::Tensor<std::complex<Scalar>, rank>;
  template <math::Index rank>
  using CmplxTensorMap = math::TensorMap<std::complex<Scalar>, rank>;
  template <math::Index rank>
  using ConstCmplxTensorMap = math::ConstTensorMap<std::complex<Scalar>, rank>;
  using DataTensor = math::Tensor<std::complex<Scalar>, 5>;
  using DataTensorPtr = std::shared_ptr<DataTensor>;

  // pxx :: hide
  /** Create scattering data field.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param sht_inc The SH transform used to expand the incoming-angle
   * dependency.
   * @param sht_scat The SH transform used to expand the scattering-angle
   * dependency.
   * @data The scattering data.
   */
  ScatteringDataFieldFullySpectral(VectorPtr f_grid,
                                   VectorPtr t_grid,
                                   ShtPtr sht_inc,
                                   ShtPtr sht_scat,
                                   DataTensorPtr data)
      : ScatteringDataFieldBase(f_grid->size(),
                                t_grid->size(),
                                sht_inc->get_n_longitudes(),
                                sht_inc->get_n_latitudes(),
                                sht_scat->get_n_longitudes(),
                                sht_scat->get_n_latitudes()),
        f_grid_(f_grid),
        t_grid_(t_grid),
        sht_inc_(sht_inc),
        sht_scat_(sht_scat),
        f_grid_map_(f_grid->data(), n_freqs_),
        t_grid_map_(t_grid->data(), n_temps_),
        data_(data) {}

  /** Create scattering data field.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param sht_inc The SH transform used to expand the incoming-angle
   * dependency.
   * @param sht_scat The SH transform used to expand the scattering-angle
   * dependency.
   * @data The scattering data.
   */
  ScatteringDataFieldFullySpectral(const Vector &f_grid,
                                   const Vector &t_grid,
                                   const sht::SHT &sht_inc,
                                   const sht::SHT &sht_scat,
                                   const DataTensor &data)
      : ScatteringDataFieldBase(f_grid.size(),
                                t_grid.size(),
                                sht_inc.get_n_longitudes(),
                                sht_inc.get_n_latitudes(),
                                sht_scat.get_n_longitudes(),
                                sht_scat.get_n_latitudes()),
        f_grid_(std::make_shared<Vector>(f_grid)),
        t_grid_(std::make_shared<Vector>(t_grid)),
        sht_inc_(std::make_shared<sht::SHT>(sht_inc)),
        sht_scat_(std::make_shared<sht::SHT>(sht_scat)),
        f_grid_map_(f_grid_->data(), n_freqs_),
        t_grid_map_(t_grid_->data(), n_temps_),
        data_(std::make_shared<DataTensor>(data)) {}

  /** Create empty scattering data field.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @param sht_inc The SH transform used to expand the incoming-angle
   * dependency.
   * @param sht_scat The SH transform used to expand the scattering-angle
   * dependency.
   * @data The scattering data.
   */
  ScatteringDataFieldFullySpectral(const Vector &f_grid,
                                   const Vector &t_grid,
                                   const sht::SHT &sht_inc,
                                   const sht::SHT &sht_scat,
                                   Index n_elements)
      : ScatteringDataFieldBase(f_grid.size(),
                                t_grid.size(),
                                sht_inc.get_n_longitudes(),
                                sht_inc.get_n_latitudes(),
                                sht_scat.get_n_longitudes(),
                                sht_scat.get_n_latitudes()),
        f_grid_(std::make_shared<Vector>(f_grid)),
        t_grid_(std::make_shared<Vector>(t_grid)),
        sht_inc_(std::make_shared<sht::SHT>(sht_inc)),
        sht_scat_(std::make_shared<sht::SHT>(sht_scat)),
        f_grid_map_(f_grid_->data(), n_freqs_),
        t_grid_map_(t_grid_->data(), n_temps_),
        data_(std::make_shared<DataTensor>(
            std::array<Index, 5>{f_grid.size(),
                                 t_grid.size(),
                                 sht_inc.get_n_spectral_coeffs_cmplx(),
                                 sht_scat.get_n_spectral_coeffs(),
                                 n_elements})) {}

  /// Shallow copy of the ScatteringDataField.
  ScatteringDataFieldFullySpectral(const ScatteringDataFieldFullySpectral &) =
      default;

  /// Deep copy of the scattering data.
  ScatteringDataFieldFullySpectral copy() const {
    auto data_new = std::make_shared<DataTensor>(*data_);
    return ScatteringDataFieldFullySpectral(f_grid_,
                                            t_grid_,
                                            sht_inc_,
                                            sht_scat_,
                                            data_new);
  }

  /// Return enum identifying the data format.
  constexpr DataFormat get_data_format() const { return DataFormat::FullySpectral; }

  /// The number of scattering-data coefficients.
  Index get_n_coeffs() const { return data_->dimension(4); }
  /// Parameters of SHT transformation used to transform
  /// scattering angle.
  std::array<Index, 4> get_sht_inc_params() const {
      return {sht_inc_->get_l_max(),
              sht_inc_->get_m_max(),
              sht_inc_->get_n_longitudes(),
              sht_inc_->get_n_latitudes()};
  }
  /// The frequency grid.
  /// Parameters of SHT transformation used to transform
  /// scattering angle.
  std::array<Index, 4> get_sht_scat_params() const {
      return {sht_scat_->get_l_max(),
              sht_scat_->get_m_max(),
              sht_scat_->get_n_longitudes(),
              sht_scat_->get_n_latitudes()};
  }
  /// The frequency grid.
  const math::Vector<double>& get_f_grid() { return *f_grid_; }
  /// The temperature grid.
  const math::Vector<double>& get_t_grid() { return *t_grid_; }
  /// The incoming-angle longitude grid.
  math::Vector<double> get_lon_inc() { return sht_inc_->get_longitude_grid(); }
  /// The incoming-angle latitude grid.
  math::Vector<double> get_lat_inc() { return sht_inc_->get_latitude_grid(); }
  /// The scattering-angle longitude grid.
  math::Vector<double> get_lon_scat() {
    return sht_scat_->get_longitude_grid();
  }
  /// The scattering-angle latitude grid.
  math::Vector<double> get_lat_scat() {
    return sht_scat_->get_latitude_grid();
  }

  /// The SHT object representing the SHT transform used to transform
  /// the data along the incoming angles.
  sht::SHT &get_sht_inc() const { return *sht_inc_; }
  /// The SHT object representing the SHT transform used to transform
  /// the data along the scattering angles.
  sht::SHT &get_sht_scat() const { return *sht_scat_; }

  /** Set scattering data for given frequency and temperature index.
   *
   * This function copies the data from the given scattering data field
   * into the sub-tensor of this objects' data tensor identified by
   * the given frequency and temperature indices. The data is automatically
   * regridded to the scattering data grids of this object.
   *
   * This function is useful to combine scattering data at different
   * temperatures and frequencies that have different scattering grids.
   *
   * @frequency_index The index along the frequency dimension
   * @temperature_index The index along the temperature dimension.
   */
  void set_data(math::Index frequency_index,
                math::Index temperature_index,
                const ScatteringDataFieldFullySpectral &other) {
    std::array<math::Index, 2> data_index = {frequency_index,
                                              temperature_index};
    std::array<math::Index, 2> input_index = {0, 0};
    auto data_map = math::tensor_index(*data_, data_index);
    auto other_data_map = math::tensor_index(*other.data_, input_index);

    math::IndexArray<1> dimensions_loop = {data_->dimension(5)};
    for (math::DimensionCounter<1> i{dimensions_loop}; i; ++i) {
      auto result = math::get_submatrix<0, 1>(data_map, i.coordinates);
      auto in_l = math::get_submatrix<0, 1>(data_map, i.coordinates);
      auto in_r = math::get_submatrix<0, 1>(other_data_map, i.coordinates);
      result = sht::SHT::add_coeffs(*sht_inc_,
                                    *sht_scat_,
                                    in_l,
                                    *other.sht_inc_,
                                    *other.sht_scat_,
                                    in_r);
    }
  }

  Vector interpolate(Scalar frequency,
                     Scalar temperature,
                     Scalar lon_inc,
                     Scalar lat_inc,
                     Scalar lon_scat,
                     Scalar lat_scat)  const {
      return to_spectral().interpolate(
          frequency, temperature, lon_inc, lat_inc,
          lon_scat, lat_scat);
  }

  // pxx :: hide
  /** Linear interpolation along frequency dimension.
   * @param frequencies The frequencies to which to interpolate the data.
   * @return New ScatteringDataFieldFullySpectral with the data interpolated
   * to the given frequencies.
   */
  ScatteringDataFieldFullySpectral interpolate_frequency(
      VectorPtr frequencies) const {
    using Regridder = RegularRegridder<Scalar, 0>;
    Regridder regridder({*f_grid_}, {*frequencies});
    auto dimensions_new = data_->dimensions();
    auto data_interp = regridder.regrid(*data_);
    dimensions_new[0] = frequencies->size();
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldFullySpectral(frequencies,
                                            t_grid_,
                                            sht_inc_,
                                            sht_scat_,
                                            data_new);
  }

  /** Linear interpolation along frequency dimension.
   * @param frequencies The frequencies to which to interpolate the data.
   * @return New ScatteringDataFieldFullySpectral with the data interpolated
   * to the given frequencies.
   */
  ScatteringDataFieldFullySpectral interpolate_frequency(
      const Vector &frequencies) const {
    return interpolate_frequency(std::make_shared<Vector>(frequencies));
  }

  // pxx :: hide
  /** Linear interpolation along temperature dimension.
   * @param temperatures The temperatures to which to interpolate the data.
   * @return New ScatteringDataFieldFullySpectral with the data interpolated
   * to the given temperatures.
   */
  ScatteringDataFieldFullySpectral interpolate_temperature(
      VectorPtr temperatures,
      bool extrapolate=false) const {
    using Regridder = RegularRegridder<Scalar, 1>;
    Regridder regridder({*t_grid_}, {*temperatures}, extrapolate);
    auto dimensions_new = data_->dimensions();
    auto data_interp = regridder.regrid(*data_);
    dimensions_new[1] = temperatures->size();
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    ;
    return ScatteringDataFieldFullySpectral(f_grid_,
                                            temperatures,
                                            sht_inc_,
                                            sht_scat_,
                                            data_new);
  }

  /** Linear interpolation along temperature dimension.
   * @param temperatures The temperatures to which to interpolate the data.
   * @return New ScatteringDataFieldFullySpectral with the data interpolated
   * to the given temperatures.
   */
  ScatteringDataFieldFullySpectral interpolate_temperature(
      const Vector &temperatures,
      bool extrapolate=false) const {
      return interpolate_temperature(std::make_shared<Vector>(temperatures),
                                     extrapolate);
  }

  /** Regrid data to new frequency and temperature grids.
   * @param f_grid The frequency grid.
   * @param t_grid The temperature grid.
   * @return A new ScatteringDataFieldFullySpectral with the given grids.
   */
  // pxx :: hide
  ScatteringDataFieldFullySpectral regrid(VectorPtr f_grid,
                                          VectorPtr t_grid) const {
    using Regridder = RegularRegridder<Scalar, 0, 1>;
    Regridder regridder({*f_grid_, *t_grid_}, {*f_grid, *t_grid});
    auto data_interp = regridder.regrid(*data_);
    auto data_new = std::make_shared<DataTensor>(std::move(data_interp));
    return ScatteringDataFieldFullySpectral(f_grid,
                                            t_grid,
                                            sht_inc_,
                                            sht_scat_,
                                            data_new);
  }

  /** Accumulate scattering data into this object.
   *
   * Regrids the given scattering data field and accumulates its interpolated
   * data tensor into this object's data tensor.
   *
   * @param other The ScatteringDataField to accumulate into this.
   * @return Reference to this object.
   */
  ScatteringDataFieldFullySpectral &operator+=(
      const ScatteringDataFieldFullySpectral &other) {
    auto regridded = other.regrid(f_grid_, t_grid_);
    math::IndexArray<3> dimensions_loop = {n_freqs_,
                                            n_temps_,
                                            data_->dimension(4)};
    for (math::DimensionCounter<3> i{dimensions_loop}; i; ++i) {
      auto result = math::get_submatrix<2, 3>(*data_, i.coordinates);
      auto in_l = math::get_submatrix<2, 3>(*data_, i.coordinates);
      auto in_r = math::get_submatrix<2, 3>(*regridded.data_, i.coordinates);
      result = sht::SHT::add_coeffs(*sht_inc_,
                                    *sht_scat_,
                                    in_l,
                                    *regridded.sht_inc_,
                                    *regridded.sht_scat_,
                                    in_r);
    }
    return *this;
  }

  /** Add scattering data fields.
   *
   * Regrids the given scattering data field to the grids of this object and
   * computes the sum of the two scattering data fields.
   *
   * @param other The other scattering data field to add to this.
   * @return A new ScatteringDataFieldFullySpectral object representing the sum
   * of the two other fields.
   */
  ScatteringDataFieldFullySpectral operator+(
      const ScatteringDataFieldFullySpectral &other) const {
    auto result = copy();
    result += other;
    return result;
  }

  /** In-place scaling scattering data.
   *
   * @param c The scaling factor.
   * @return Reference to this object.
   */
  ScatteringDataFieldFullySpectral &operator*=(Scalar c) {
    *data_ = c * (*data_);
    return *this;
  }

  /** Scale scattering data.
   *
   * @param c The scaling factor.
   * @return A new object containing the scaled scattering data.
   */
  ScatteringDataFieldFullySpectral operator*(Scalar c) const {
    auto result = copy();
    result *= c;
    return result;
  }

  /** Set the number of scattering coefficients.
   *
   * De- or increases the number of scattering coefficients that are stored.
   * If the number of scattering coefficients is increased, the new elements are
   * set to 0.
   *
   * @param n The number of scattering coefficients to change the data to have.
   */
  void set_number_of_scattering_coeffs(Index n) {
      Index current_stokes_dim = data_->dimension(4);
      if (current_stokes_dim == n) {
          return;
      }
      auto new_dimensions = data_->dimensions();
      new_dimensions[4] = n;
      DataTensorPtr data_new = std::make_shared<DataTensor>(new_dimensions);
      math::copy(*data_new, *data_);
      data_ = data_new;
  }

  /// Convert to spectral scattering data format using native SHT parameters.
  ScatteringDataFieldSpectral<Scalar> to_spectral() const;

  /** Change SHT representation of incoming angles.
   *
   * @param sht_other SHT-object representing the transformation to which
   * to convert the data.
   * @return Scattering data field in fully-spectral format but with different
   * transformation parameters.
   */
  ScatteringDataFieldSpectral<Scalar> to_spectral(ShtPtr sht_other) const {
    auto new_dimensions = data_->dimensions();
    new_dimensions[3] = sht_other->get_n_spectral_coeffs();
    auto data_new_ =
        std::make_shared<DataTensor>(DataTensor(new_dimensions).setZero());
    auto result = ScatteringDataFieldFullySpectral(f_grid_,
                                                   t_grid_,
                                                   sht_inc_,
                                                   sht_other,
                                                   data_new_);
    result += *this;
    return result.to_spectral();
  }

  /** Change SHT representation of incoming angles.
   *
   * @param l_max The l_max parameter of the SHT transform to use to transform
   * the incoming angles.
   * @param m_max The m_max parameter of the SHT transform to use to transform
   * the incoming angles.
   * @return Scattering data field in fully-spectral format but with different
   * transformation parameters.
   */
  ScatteringDataFieldSpectral<Scalar> to_spectral(Index l_max,
                                                  Index m_max) const {
    auto n_lat = sht_scat_->get_n_latitudes();
    auto n_lon = sht_scat_->get_n_longitudes();
    auto sht_other = std::make_shared<sht::SHT>(l_max, m_max, n_lon, n_lat);
    return to_spectral(sht_other);
  }

  /** Change SHT representation of incoming angles.
   *
   * @param l_max The l_max parameter of the SHT transform to use to transform
   * the incoming angles.
   * @param m_max The m_max parameter of the SHT transform to use to transform
   * the incoming angles.
   * @param n_lon The size of the longitude-angle grid to use for the new
   * SHT transform.
   * @param n_lat The size of the latitude-angle grid to use for the new
   * SHT transform.
   * @return Scattering data field in fully-spectral format but with different
   * transformation parameters.
   */
  ScatteringDataFieldSpectral<Scalar> to_spectral(Index l_max,
                                                  Index m_max,
                                                  Index n_lon,
                                                  Index n_lat) const {
      auto sht_other = std::make_shared<sht::SHT>(l_max, m_max, n_lon, n_lat);
      return to_spectral(sht_other);
  }

  const DataTensor &get_data() const { return *data_; }


  /** Phase function maximum.
   *
   * Maximum of the phase function calculated across the dimensions
   * of incoming angles.
   */
  math::Tensor<Scalar, 4> get_phase_function_max_inc() const {
      return to_spectral().get_phase_function_max_inc();
  }

  /** Phase function maximum.
   *
   * Maximum of the phase function calculated across the dimensions
   * of scattering angles.
   */
  math::Tensor<Scalar, 4> get_phase_function_max_scat() const {
      std::array<int, 2> dims({4, 5});
      return to_spectral().get_phase_function_max_scat();
  }

 protected:
  VectorPtr f_grid_;
  VectorPtr t_grid_;
  VectorPtr lon_inc_;
  VectorPtr lat_inc_;
  ShtPtr sht_inc_;
  ShtPtr sht_scat_;

  ConstVectorMap f_grid_map_;
  ConstVectorMap t_grid_map_;

  DataTensorPtr data_;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation of conversion methods.
////////////////////////////////////////////////////////////////////////////////

template <typename Scalar>
ScatteringDataFieldSpectral<Scalar>
ScatteringDataFieldGridded<Scalar>::to_spectral(
    std::shared_ptr<sht::SHT> sht) const {
  math::IndexArray<5> dimensions_loop = {n_freqs_,
                                          n_temps_,
                                          n_lon_inc_,
                                          n_lat_inc_,
                                          data_->dimension(6)};
  math::IndexArray<6> dimensions_new = {n_freqs_,
                                         n_temps_,
                                         n_lon_inc_,
                                         n_lat_inc_,
                                         sht->get_n_spectral_coeffs(),
                                         data_->dimension(6)};
  using CmplxDataTensor = math::Tensor<std::complex<Scalar>, 6>;
  auto data_new = std::make_shared<CmplxDataTensor>(dimensions_new);
  for (math::DimensionCounter<5> i{dimensions_loop}; i; ++i) {
    math::get_subvector<4>(*data_new, i.coordinates) =
        sht->transform(math::get_submatrix<4, 5>(*data_, i.coordinates));
  }
  return ScatteringDataFieldSpectral<Scalar>(f_grid_,
                                             t_grid_,
                                             lon_inc_,
                                             lat_inc_,
                                             sht,
                                             data_new);
}

template <typename Scalar>
ScatteringDataFieldGridded<Scalar>
ScatteringDataFieldSpectral<Scalar>::to_gridded() const {
  math::IndexArray<5> dimensions_loop = {n_freqs_,
                                          n_temps_,
                                          n_lon_inc_,
                                          n_lat_inc_,
                                          data_->dimension(5)};
  math::IndexArray<7> dimensions_new = {n_freqs_,
                                         n_temps_,
                                         n_lon_inc_,
                                         n_lat_inc_,
                                         sht_scat_->get_n_longitudes(),
                                         sht_scat_->get_n_latitudes(),
                                         data_->dimension(5)};
  using Vector = math::Vector<Scalar>;
  using DataTensor = math::Tensor<Scalar, 7>;
  auto data_new = std::make_shared<DataTensor>(dimensions_new);
  for (math::DimensionCounter<5> i{dimensions_loop}; i; ++i) {
      auto synthesized = sht_scat_->synthesize(math::get_subvector<4>(*data_, i.coordinates));
      math::get_submatrix<4, 5>(*data_new, i.coordinates) = synthesized;
  }
  auto lon_scat_ = std::make_shared<Vector>(sht_scat_->get_longitude_grid());
  auto lat_scat_ = std::make_shared<sht::SHT::LatGrid>(sht_scat_->get_latitude_grid());
  return ScatteringDataFieldGridded<Scalar>(f_grid_,
                                            t_grid_,
                                            lon_inc_,
                                            lat_inc_,
                                            lon_scat_,
                                            lat_scat_,
                                            data_new);
}

template <typename Scalar>
ScatteringDataFieldFullySpectral<Scalar>
ScatteringDataFieldSpectral<Scalar>::to_fully_spectral(
    std::shared_ptr<sht::SHT> sht) const {
  math::IndexArray<4> dimensions_loop = {n_freqs_,
                                          n_temps_,
                                          data_->dimension(4),
                                          data_->dimension(5)};
  math::IndexArray<5> dimensions_new = {n_freqs_,
                                         n_temps_,
                                         sht->get_n_spectral_coeffs_cmplx(),
                                         data_->dimension(4),
                                         data_->dimension(5)};
  using CmplxDataTensor = math::Tensor<std::complex<Scalar>, 5>;
  auto data_new = std::make_shared<CmplxDataTensor>(dimensions_new);
  for (math::DimensionCounter<4> i{dimensions_loop}; i; ++i) {
    math::get_subvector<2>(*data_new, i.coordinates) =
        sht->transform_cmplx(math::get_submatrix<2, 3>(*data_, i.coordinates));
  }
  return ScatteringDataFieldFullySpectral<Scalar>(f_grid_,
                                                  t_grid_,
                                                  sht,
                                                  sht_scat_,
                                                  data_new);
}

template <typename Scalar>
ScatteringDataFieldSpectral<Scalar>
ScatteringDataFieldFullySpectral<Scalar>::to_spectral() const {
  math::IndexArray<4> dimensions_loop = {n_freqs_,
                                          n_temps_,
                                          data_->dimension(3),
                                          data_->dimension(4)};
  math::IndexArray<6> dimensions_new = {n_freqs_,
                                         n_temps_,
                                         sht_inc_->get_n_longitudes(),
                                         sht_inc_->get_n_latitudes(),
                                         data_->dimension(3),
                                         data_->dimension(4)};
  using CmplxDataTensor = math::Tensor<std::complex<Scalar>, 6>;
  auto data_new = std::make_shared<CmplxDataTensor>(dimensions_new);
  for (math::DimensionCounter<4> i{dimensions_loop}; i; ++i) {
    math::get_submatrix<2, 3>(*data_new, i.coordinates) =
        sht_inc_->synthesize_cmplx(
            math::get_subvector<2>(*data_, i.coordinates));
  }

  auto lon_inc_ = std::make_shared<Vector>(sht_inc_->get_longitude_grid());
  auto lat_inc_ = std::make_shared<sht::SHT::LatGrid>(sht_inc_->get_latitude_grid());

  return ScatteringDataFieldSpectral<Scalar>(f_grid_,
                                             t_grid_,

                                             lon_inc_,
                                             lat_inc_,
                                             sht_scat_,
                                             data_new);
}

}  // namespace scattering
