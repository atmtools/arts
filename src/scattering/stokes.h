/** \file stokes.h
 *
 * Defines class to represent the different scattering quantities: The absorption
 * vector, the extinction matrix and the phase matrix. All of these use slightly
 * different compact storage formats and need to be assembled in different ways.
 * The functions defined in this file take care of this.
 *
 * @author Simon Pfreundschuh, 2020
 */
#pragma once

#include "scattering/maths.h"
#include "scattering/scattering_data_field.h"

namespace scattering {
namespace stokes {

/// [1, 1] element of phase matrix stored in compact format.
template <typename VectorType>
auto f11(const VectorType &v) {
  return v[0];
}

/// [1, 2] element of phase matrix stored in compact format.
template <typename VectorType>
auto f12(const VectorType &v) {
  return v[1];
}

/// [2, 2] element of phase matrix stored in compact format.
template <typename VectorType>
auto f22(const VectorType &v) {
  return v[2];
}

/// [3, 3] element of phase matrix stored in compact format.
template <typename VectorType>
auto f33(const VectorType &v) {
  return v[3];
}

/// [3, 4] element of phase matrix stored in compact format.
template <typename VectorType>
auto f34(const VectorType &v) {
  return v[4];
}

/// [4, 4] element of phase matrix stored in compact format.
template <typename VectorType>
auto f44(const VectorType &v) {
  return v[5];
}

////////////////////////////////////////////////////////////////////////////////
// Calculation of rotation coefficients.
////////////////////////////////////////////////////////////////////////////////

/** Calculate angle between incoming and outgoing directions in the scattering
 *  plane.
 *
 * @param lon_inc The incoming-angle longitude component in radians.
 * @param lat_inc The incoming-angle latitude component in radians.
 * @param lon_scat The outgoing (scattering) angle longitude component in
 * radians.
 * @param lat_scat The outgoing (scattering) angle longitude component in
 * radians.
 * @return The angle between the incoming and outgoing directions in
 * radians.
 */
template <typename Scalar>
Scalar scattering_angle(Scalar lon_inc,
                        Scalar lat_inc,
                        Scalar lon_scat,
                        Scalar lat_scat) {
    Scalar cos_theta = cos(lat_inc) * cos(lat_scat) +
        sin(lat_inc) * sin(lat_scat) * cos(lon_scat - lon_inc);
    return math::save_acos(cos_theta);
}

/** Calculate rotation coefficients for scattering matrix.
 *
 * This method calculates the rotation coefficients that are required to
 * transform the phase matrix of a randomly-oriented particle to the Stokes
 * phase  matrix, which describes its scattering behavior w.r.t. to the
 * laboratory frame. This equation calculates the angle Theta and the
 * coefficients C_1, * C_2, S_1, S_2 as defined in equation (4.16) of
 * "Scattering, Absorption, and Emission of Light by Small Particles."
 *
 * @param lon_inc The longitude component of the incoming angle in radians.
 * @param lat_inc The latitude component of the incoming angle in radians.
 * @param lon_scat The longitude component of the scattering angle in radians.
 * @param lat_scat The latitude component of the scattering angle in radians.
 */
template <typename Scalar>
std::array<Scalar, 5> rotation_coefficients(Scalar lon_inc,
                                            Scalar lat_inc,
                                            Scalar lon_scat,
                                            Scalar lat_scat) {
  Scalar cos_theta = cos(lat_inc) * cos(lat_scat) +
                     sin(lat_inc) * sin(lat_scat) * cos(lon_scat - lon_inc);
  Scalar theta = math::save_acos(cos_theta);
  if ((math::small(abs(lon_scat - lon_inc))) ||
      (math::equal(abs(lon_scat - lon_inc), 2.0 * M_PI))) {
    theta = abs(lat_inc - lat_scat);
  } else if ((math::equal(lon_scat - lon_inc, 2.0 * M_PI))) {
    theta = lat_scat + lat_inc;
    if (theta > M_PI) {
      theta = 2.0 * M_PI - theta;
    }
  }

  Scalar sigma_1, sigma_2;

  if (math::small(lat_inc)) {
      sigma_1 = lon_scat - lon_inc;
      sigma_2 = 0.0;
  } else if (math::equal(lat_inc, M_PI)) {
      sigma_1 = lon_scat - lon_inc;
      sigma_2 = M_PI;
  } else if (math::small(lat_scat)) {
      sigma_1 = 0.0;
      sigma_2 = M_PI + lon_scat - lon_inc;
  } else if (math::equal(lat_scat, M_PI)) {
      sigma_1 = M_PI;
      sigma_2 = lon_scat - lon_inc;
  } else {
      sigma_1 = math::save_acos((cos(lat_scat) - cos(lat_inc) * cos_theta) /
                                (sin(lat_inc) * sin(theta)));
      sigma_2 = math::save_acos((cos(lat_inc) - cos(lat_scat) * cos_theta) /
                                (sin(lat_scat) * sin(theta)));
  }

  Scalar c_1 = cos(2.0 * sigma_1);
  Scalar c_2 = cos(2.0 * sigma_2);
  Scalar s_1 = sin(2.0 * sigma_1);
  Scalar s_2 = sin(2.0 * sigma_2);

  return {theta, c_1, c_2, s_1, s_2};
}

////////////////////////////////////////////////////////////////////////////////
// Manipulation of compact phase matrix format.
////////////////////////////////////////////////////////////////////////////////

template <typename Derived>
/** Converter class for the compact phase-function.
 *
 * This helper struct implements the functionality to convert the phase
 * function of spherical particles given in compact format to the full
 * Stokes phase matrix.
 */
struct CompactFormatBase {
  template <typename TensorType, typename VectorType, typename MatrixType>

  /** Expand and transform compact phase function coefficients.
   *
   * @param output Tensor to which to write the results.
   * @scat_angs The scattering angles corresponding to the angular grid
   * of the desired output.
   * @lon_scat The longitude component of the scattering angle in radians.
   * @rotation_coeffs The precalculated rotation coefficients.
   */
  static void expand_and_transform(TensorType &output,
                                   const TensorType &input,
                                   const VectorType &scat_angs,
                                   const VectorType &lon_scat,
                                   const MatrixType &rotation_coeffs) {
    using Scalar = typename TensorType::Scalar;
    using CoefficientVectorMap = Eigen::Map<
        const Eigen::Matrix<Scalar, 1, Derived::n_coeffs, Eigen::RowMajor>>;
    using PhaseMatrixMap = Eigen::Map<Eigen::Matrix<Scalar,
                                                    Derived::stokes_dim,
                                                    Derived::stokes_dim,
                                                    Eigen::RowMajor>>;

    math::IndexArray<6> dimensions_loop = {input.dimension(0),
                                            input.dimension(1),
                                            input.dimension(2),
                                            input.dimension(3),
                                            input.dimension(4),
                                            input.dimension(5)};
    math::Index scat_ang_index = 0;
    math::Index input_index = 0;
    math::Index output_index = 0;
    for (math::DimensionCounter<6> i{dimensions_loop}; i; ++i) {
      PhaseMatrixMap output_matrix(output.data() + output_index);
      CoefficientVectorMap input_vector(input.data() + input_index);

      scat_ang_index %= rotation_coeffs.rows();
      auto scat_ang = scat_angs[scat_ang_index];

      if (math::small(scat_ang) || math::equal(scat_ang, M_PI)) {
        Derived::expand(output_matrix, input_vector);
      } else {
        bool lon_scat_gt_pi = lon_scat[i.coordinates[4]] > M_PI;
        Derived::expand_and_transform(output_matrix,
                                      input_vector,
                                      rotation_coeffs(scat_ang_index, 0),
                                      rotation_coeffs(scat_ang_index, 1),
                                      rotation_coeffs(scat_ang_index, 2),
                                      rotation_coeffs(scat_ang_index, 3),
                                      lon_scat_gt_pi);
      }
      scat_ang_index++;
      output_index += output.dimension(6);
      input_index += input.dimension(6);
    }
  }


  /** Expand and transform single angle.
    *
    * @param output Matrix to write the expanded stokes phase matrix to.
    * @input The compressed phase function coefficients.
    * @lon_inc Longitude-component of the incoming angle in radians.
    * @lat_inc Latitude-component of the incoming angle in radians.
    * @lon_scat Longitude-component of the scattering angle in radians.
    * @lat_scat Latitude-component of the scattering angle.
    */
  template <typename MatrixType, typename VectorType, typename Scalar>
  static void expand_and_transform(MatrixType &output,
                                   const VectorType &input,
                                   Scalar lon_scat,
                                   std::array<Scalar, 5> rotation_coeffs) {
    bool lon_scat_gt_pi = lon_scat > M_PI;
    auto scat_ang = rotation_coeffs[0];
    if (math::small(scat_ang) || math::equal(scat_ang, M_PI)) {
      Derived::expand(output, input);
    } else {
      Derived::expand_and_transform(output,
                                    input,
                                    rotation_coeffs[1],
                                    rotation_coeffs[2],
                                    rotation_coeffs[3],
                                    rotation_coeffs[4],
                                    lon_scat_gt_pi);
    }
  }
};

template <math::Index stokes_dim_> struct CompactFormat;

/// Specialized compact format class for stokes dimension 4.
template <>
struct CompactFormat<4> : public CompactFormatBase<CompactFormat<4>> {

  using CompactFormatBase<CompactFormat<4>>::expand_and_transform;
  static constexpr math::Index n_coeffs = 6;
  static constexpr math::Index stokes_dim = 4;

  /// Sort coefficients from compact vector into Stokes phase matrix.
  template <typename MatrixType, typename VectorType>
  static void expand(MatrixType &output, const VectorType &input) {
    output(0, 0) = input[0];
    output(0, 1) = f12(input);
    output(1, 0) = f12(input);
    output(1, 1) = f22(input);
    output(0, 2) = 0.0;
    output(1, 2) = 0.0;
    output(2, 0) = 0.0;
    output(2, 1) = 0.0;
    output(2, 2) = f33(input);
    output(0, 3) = 0.0;
    output(1, 3) = 0.0;
    output(2, 3) = f34(input);
    output(3, 0) = 0.0;
    output(3, 1) = 0.0;
    output(3, 2) = -f34(input);
    output(3, 3) = f44(input);
  }

  /// Sort and transform coefficients from compact vector into Stokes phase
  /// matrix.
  template <typename MatrixType, typename VectorType, typename Scalar>
  static void expand_and_transform(MatrixType &output,
                                   const VectorType &input,
                                   Scalar c_1,
                                   Scalar c_2,
                                   Scalar s_1,
                                   Scalar s_2,
                                   bool lon_scat_gt_pi) {
    output(0, 0) = input[0];

    output(0, 1) = c_1 * f12(input);
    output(1, 0) = c_2 * f12(input);
    output(1, 1) = c_1 * c_2 * f22(input) - s_1 * s_2 * f33(input);

    output(0, 2) = s_1 * f12(input);
    output(1, 2) = s_1 * c_2 * f22(input) + c_1 * s_2 * f33(input);
    output(2, 0) = -s_2 * f12(input);
    output(2, 1) = -c_1 * s_2 * f22(input) - s_1 * c_2 * f33(input);
    output(2, 2) = -s_1 * s_2 * f22(input) + c_1 * c_2 * f33(input);

    if (lon_scat_gt_pi) {
      output(0, 2) *= -1.0;
      output(1, 2) *= -1.0;
      output(2, 0) *= -1.0;
      output(2, 1) *= -1.0;
    }

    output(0, 3) = 0.0;
    output(1, 3) = s_2 * f34(input);
    output(2, 3) = c_2 * f34(input);
    output(3, 0) = 0.0;
    output(3, 1) = s_1 * f34(input);
    output(3, 2) = -c_1 * f34(input);
    output(3, 3) = f44(input);

    if (lon_scat_gt_pi) {
      output(1, 3) *= -1.0;
      output(3, 1) *= -1.0;
    }
  }

};

/// Specialized compact format class for stokes dimension 3.
template <>
struct CompactFormat<3> : public CompactFormatBase<CompactFormat<3>> {

  using CompactFormatBase<CompactFormat<3>>::expand_and_transform;
  static constexpr math::Index n_coeffs = 6;
  static constexpr math::Index stokes_dim = 3;

  /// Sort coefficients from compact vector into Stokes phase matrix.
  template <typename MatrixType, typename VectorType>
  static void expand(MatrixType &output, const VectorType &input) {
    output(0, 0) = input[0];
    output(0, 1) = f12(input);
    output(1, 0) = f12(input);
    output(1, 1) = f22(input);
    output(0, 2) = 0.0;
    output(1, 2) = 0.0;
    output(2, 0) = 0.0;
    output(2, 1) = 0.0;
    output(2, 2) = f33(input);
  }

  /// Sort and transform coefficients from compact vector into Stokes phase
  /// matrix.
  template <typename MatrixType, typename VectorType, typename Scalar>
  static void expand_and_transform(MatrixType &output,
                                   const VectorType &input,
                                   Scalar c_1,
                                   Scalar c_2,
                                   Scalar s_1,
                                   Scalar s_2,
                                   bool lon_scat_gt_pi) {

    output(0, 0) = input[0];

    output(0, 1) = c_1 * f12(input);
    output(1, 0) = c_2 * f12(input);
    output(1, 1) = c_1 * c_2 * f22(input) - s_1 * s_2 * f33(input);

    output(0, 2) = s_1 * f12(input);
    output(1, 2) = s_1 * c_2 * f22(input) + c_1 * s_2 * f33(input);
    output(2, 0) = -s_2 * f12(input);
    output(2, 1) = -c_1 * s_2 * f22(input) - s_1 * c_2 * f33(input);
    output(2, 2) = -s_1 * s_2 * f22(input) + c_1 * c_2 * f33(input);

    if (lon_scat_gt_pi) {
        output(0, 2) *= -1.0;
        output(1, 2) *= -1.0;
        output(2, 0) *= -1.0;
        output(2, 1) *= -1.0;
    }
  }
};

/// Specialized compact format class for stokes dimension 2.
template <>
struct CompactFormat<2> : public CompactFormatBase<CompactFormat<2>> {

    using CompactFormatBase<CompactFormat<2>>::expand_and_transform;
    static constexpr math::Index n_coeffs = 4;
    static constexpr math::Index stokes_dim = 2;

  /// Sort coefficients from compact vector into Stokes phase matrix.
  template <typename MatrixType, typename VectorType>
  static void expand(MatrixType &output, const VectorType &input) {
    output(0, 0) = input[0];
    output(0, 1) = f12(input);
    output(1, 0) = f12(input);
    output(1, 1) = f22(input);
  }

  /// Sort and transform coefficients from compact vector into Stokes phase
  /// matrix.
  template <typename MatrixType, typename VectorType, typename Scalar>
  static void expand_and_transform(MatrixType &output,
                                   const VectorType &input,
                                   Scalar c_1,
                                   Scalar c_2,
                                   Scalar s_1,
                                   Scalar s_2,
                                   bool /*lon_scat_gt_pi*/) {
    output(0, 0) = input[0];

    output(0, 1) = c_1 * f12(input);
    output(1, 0) = c_2 * f12(input);
    output(1, 1) = c_1 * c_2 * f22(input) - s_1 * s_2 * f33(input);
  }
};

/// Specialized compact format class for stokes dimension 1.
template <>
struct CompactFormat<1> : public CompactFormatBase<CompactFormat<1>> {

  using CompactFormatBase<CompactFormat<1>>::expand_and_transform;
  static constexpr math::Index n_coeffs = 1;
  static constexpr math::Index stokes_dim = 1;

  /// Sort coefficients from compact vector into Stokes phase matrix.
  template <typename MatrixType, typename VectorType>
  static void expand(MatrixType &output, const VectorType &input) {
    output(0, 0) = input[0];
  }

  /// Sort and transform coefficients from compact vector into Stokes phase
  /// matrix.
  template <typename MatrixType, typename VectorType, typename Scalar>
  static void expand_and_transform(MatrixType &output,
                                   const VectorType &input,
                                   Scalar /*c_1*/,
                                   Scalar /*c_2*/,
                                   Scalar /*s_1*/,
                                   Scalar /*s_2*/,
                                   bool /*lon_scat_gt_pi*/) {
    output(0, 0) = input[0];
  }
};

////////////////////////////////////////////////////////////////////////////////
// Phase matrix
////////////////////////////////////////////////////////////////////////////////

// pxx :: hide
// pxx :: instance(["scattering::ScatteringDataFieldGridded<double>"])
/** The Phase matrix class.
 *
 * This class provides an interface to a scattering data field containing
 * coefficients of the phase matrix. It takes care of correctly handling the
 * the different compression schemes and reference frames that are
 * used for the scattering data.
 */
template <typename Base>
    class PhaseMatrix : public Base {


  using typename Base::Coefficient;
  using typename Base::DataTensor;
  using typename Base::DataTensorPtr;
  using typename Base::Matrix;
  using typename Base::Vector;
  using typename Base::VectorPtr;

  using Base::coeff_dim;
  using Base::data_;
  using Base::f_grid_;
  using Base::t_grid_;
  using Base::lon_inc_;
  using Base::lat_inc_;
  using Base::lon_scat_;
  using Base::lat_scat_;
  using Base::rank;

  using Base::n_temps_;
  using Base::n_freqs_;
  using Base::n_lon_inc_;
  using Base::n_lat_inc_;
  using Base::n_lon_scat_;
  using Base::n_lat_scat_;

 public:

  using PhaseFunction = math::Tensor<Coefficient, rank - 1>;
  using ScatteringMatrix = math::Tensor<Coefficient, rank + 1>;
  using Scalar = typename Base::Scalar;
  using Base::copy;

  /// Perfect forwarding constructor.
  template <typename... Args>
  PhaseMatrix(Args... args) : Base(std::forward<Args>(args)...) {}

  ParticleType get_particle_type() const {
      if (n_lon_inc_ > 1) {
          return ParticleType::General;
      }
      if ((n_lon_inc_ > 1) || (n_lat_inc_ > 1) || (n_lon_scat_ > 1)) {
          return ParticleType::AzimuthallyRandom;
      }
      return ParticleType::Random;
  }

  /// Determine stokes dimension of data.
  math::Index get_stokes_dim() const {
    auto n_coeffs = Base::get_n_coeffs();
    auto type = get_particle_type();
    if (type == ParticleType::Random) {
      if (n_coeffs == 6) {
        return 4;
      } else if (n_coeffs == 4) {
        return 3;
      }
      return 1;
    } else {
      return sqrt(n_coeffs);
    }
  }

  /** In-place reduction of scattering data to given stokes dimension.
   *
   * This function discards any stored data that exceeds a given stokes
   * dimension. This can help speed up scattering data preparation
   * in simulation that don't require all stokes components.
   *
   * @param n The maximum required Stokes dimension.
   */
  void set_stokes_dim(math::Index n) {
      auto stokes_dim_out = n;
      auto stokes_dim = get_stokes_dim();
      auto dimensions_new = data_->dimensions();

    if (get_particle_type() == ParticleType::Random) {
      if (stokes_dim_out == 1) {
        dimensions_new[coeff_dim] = 1;
      } else if (stokes_dim_out == 2) {
        dimensions_new[coeff_dim] = 4;
      } else {
        dimensions_new[coeff_dim] = 6;
      }
      auto data_new = std::make_shared<DataTensor>(dimensions_new);
      math::copy(*data_new, *data_);
      data_ = data_new;
    } else {
        dimensions_new[coeff_dim] = stokes_dim_out * stokes_dim_out;
        auto data_new = std::make_shared<DataTensor>(dimensions_new);

        auto stokes_dim_min = std::min(stokes_dim_out, stokes_dim);
        for (Index i = 0; i < stokes_dim_min; ++i) {
            for (Index j = 0; j < stokes_dim_min; ++j) {
                data_new->template chip<coeff_dim>(i * stokes_dim_out + j)
                    = data_->template chip<coeff_dim>(i * stokes_dim  + j);
            }
        }
        data_ = data_new;
    }
  }

  /** Transform phase matrix data to lab frame.
   *
   * The reference frame for scattering data of spherical particles is the
   * scattering plane. This function converts in this format to the full
   * but flattened Stokes phase matrix format defined w.r.t. the lab
   * frame.
   *
   * @param lat_inc_new Pointer to the new latitude component
   * grid of the incoming angle.
   * @param lon_scat_new Pointer to the new latitude component grid
   * of the scattering angle.
   * @param lat_scat_new Pointer to the new latitude angle grid.
   * @param stokes_dim The number of stokes components to retain.
   *
   * @return A new PhaseMatrix object containing the transformed data.
   */
  PhaseMatrix to_lab_frame(VectorPtr lat_inc_new,
                           VectorPtr lon_scat_new,
                           std::shared_ptr<LatitudeGrid<Scalar>> lat_scat_new,
                           Index stokes_dim) const {
    if ((n_lat_inc_ > 1) || (n_lon_scat_ > 1)) {
        return copy();
    }
    auto n_lat_inc = lat_inc_new->size();
    auto n_lon_scat = lon_scat_new->size();
    auto n_lat_scat = lat_scat_new->size();

    auto scat_ang_new = Vector(n_lat_inc * n_lon_scat * n_lat_scat);
    math::MatrixFixedRows<Scalar, 4> rotation_coeffs{
        n_lat_inc * n_lon_scat * n_lat_scat,
        4};

    Index index = 0;
    for (Index i = 0; i < n_lat_inc; ++i) {
      for (Index j = 0; j < n_lon_scat; ++j) {
        for (Index k = 0; k < n_lat_scat; ++k) {
          auto lat_inc = lat_inc_new->operator[](i);
          auto lon_scat = lon_scat_new->operator[](j);
          auto lat_scat = lat_scat_new->operator[](k);
          auto coeffs = rotation_coefficients(0.0, lat_inc, lon_scat, lat_scat);
          scat_ang_new[index] = coeffs[0];
          rotation_coeffs(index, 0) = coeffs[1];
          rotation_coeffs(index, 1) = coeffs[2];
          rotation_coeffs(index, 2) = coeffs[3];
          rotation_coeffs(index, 3) = coeffs[4];
          index += 1;
        }
      }
    }
    Index stokes_dim_in = get_stokes_dim();

    // Interpolate data to scattering angles.
    using Regridder = RegularRegridder<Scalar, 5>;
    Regridder regridder({*lat_scat_}, {scat_ang_new});
    math::IndexArray<7> dimensions_interp = {n_freqs_,
                                              n_temps_,
                                              1,
                                              n_lat_inc,
                                              n_lon_scat,
                                              n_lat_scat,
                                              Base::get_n_coeffs()};
    auto data_interp = regridder.regrid(*data_);
    data_interp.resize(dimensions_interp);

    stokes_dim = std::min(stokes_dim_in, stokes_dim);

    // Tensor to hold results.
    math::IndexArray<7> dimensions_new = {n_freqs_,
                                           n_temps_,
                                           1,
                                           n_lat_inc,
                                           n_lon_scat,
                                           n_lat_scat,
                                           stokes_dim * stokes_dim};
    auto data_new = std::make_shared<DataTensor>(dimensions_new);

    if (stokes_dim == 1) {
      CompactFormat<1>::expand_and_transform(*data_new,
                                             data_interp,
                                             scat_ang_new,
                                             *lon_scat_new,
                                             rotation_coeffs);
    } else if (stokes_dim == 2) {
      CompactFormat<2>::expand_and_transform(*data_new,
                                             data_interp,
                                             scat_ang_new,
                                             *lon_scat_new,
                                             rotation_coeffs);
    } else if (stokes_dim == 3) {
      CompactFormat<3>::expand_and_transform(*data_new,
                                             data_interp,
                                             scat_ang_new,
                                             *lon_scat_new,
                                             rotation_coeffs);
    } else {
      CompactFormat<4>::expand_and_transform(*data_new,
                                             data_interp,
                                             scat_ang_new,
                                             *lon_scat_new,
                                             rotation_coeffs);
    }
    return PhaseMatrix(f_grid_,
                       t_grid_,
                       lon_inc_,
                       lat_inc_new,
                       lon_scat_new,
                       lat_scat_new,
                       data_new);
  }

  /** Transform phase matrix data to lab frame.
   *
   * The reference frame for scattering data of spherical particles is the
   * scattering plane. This function converts in this format to the full
   * but flattened Stokes phase matrix format defined w.r.t. the lab
   * frame.
   *
   * @param lat_inc The new latitude component grid of the incoming angle.
   * @param lon_scat The new latitude component grid of the scattering angle.
   * @param lat_scat The new latitude component grid of the scattering angle.
   * @param stokes_dim The number of stokes components to retain.
   *
   * @return A new PhaseMatrix object containing the transformed data.
   */
  PhaseMatrix to_lab_frame(Vector lat_inc,
                           Vector lon_scat,
                           Vector lat_scat,
                           Index stokes_dim) const {
    return to_lab_frame(std::make_shared<Vector>(lat_inc),
                        std::make_shared<Vector>(lon_scat),
                        std::make_shared<IrregularLatitudeGrid<Scalar>>(lat_scat),
                        stokes_dim);
  }

  /** Transform phase matrix data to lab frame.
   *
   * The reference frame for scattering data of spherical particles is the
   * scattering plane. This function converts in this format to the full
   * but flattened Stokes phase matrix format defined w.r.t. the lab
   * frame.
   *
   * This method uses Fejer-quadrature grids that are compatible with
   * those expected by RT4.
   *
   * @param n_lat_inc The number of grid points of the latitude component  grid
   * of the incoming angle.
   * @param n_lon_scat The number of grid points of the longitude component
   * grid of the scattering angles.
   * @param stokes_dim The number of stokes components to retain.
   *
   * @return A new PhaseMatrix object containing the transformed data.
   */
  PhaseMatrix to_lab_frame(Index n_lat_inc,
                           Index n_lon_scat,
                           Index stokes_dim) const {
      if ((n_lat_inc_ > 1) || (n_lon_scat_ > 1)) {
          return copy();
      }
      auto lon_scat_new =
          std::make_shared<Vector>(sht::SHT::get_longitude_grid(n_lon_scat));
      auto lat_scat_new =
          std::make_shared<sht::SHT::LatGrid>(sht::SHT::get_latitude_grid(n_lon_scat));
      auto lat_inc_new =
          std::make_shared<Vector>(sht::SHT::get_latitude_grid(n_lat_inc));
      return to_lab_frame(lat_inc_new, lon_scat_new, lat_scat_new, stokes_dim);
  }

  /// Return phase function data for the full scattering data field.
  PhaseFunction get_phase_function() const { return data_->template chip<rank - 1>(0); }

  /** Return phase data for the full scattering data grid.
   *
   * This function provides access to the full phase matrix data on the
   * frequency, temperature and angular grids of the underlying scattering
   * data field.
   *
   * @param stokes_dim: The required stokes dimension.
   * @return A tensor with the first dimensions corresponding to those of
   * the scattering data field and the last two to those of the stokes matrix.
   */
  ScatteringMatrix get_phase_matrix(Index stokes_dim) const {
    stokes_dim = std::min(stokes_dim, get_stokes_dim());
    auto output_dimensions = math::get_dimensions<rank + 1>(*data_);
    output_dimensions[rank - 1] = stokes_dim;
    output_dimensions[rank] = stokes_dim;
    return data_->reshape(output_dimensions);
  }

  /** Extract a single phase matrix.
   *
   * This function extracts a single phase matrix for a given frequency,
   * temperature, incoming and scattering angles.
   *
   * The phase matrix is automatically transformed to lab frame and and
   * assembled into full Stokes phase-matrix format.
   *
   * @param frequency The frequency for which to extract the phase matrix.
   * @param temperature The temperature for which to extract the phase matrix.
   * @param lon_inc: The longitude component of the incoming angle in radians.
   * @param lat_inc: The latitude component of the incoming angle in radians
   * @param lon_scat: The longitude component of the scattering angle in
   * radians.
   * @param lat_scat: The latitude component of the scattering angle.
   * @param stokes_dim: The required Stokes dimension.
   */
  Matrix get_phase_matrix(
      Scalar frequency,
      Scalar temperature,
      Scalar lon_inc,
      Scalar lat_inc,
      Scalar lon_scat,
      Scalar lat_scat,
      Index stokes_dim) const {

      Index stokes_dim_in = get_stokes_dim();
      Matrix result(stokes_dim, stokes_dim);
      stokes_dim = std::min(stokes_dim_in, stokes_dim);
      auto d_azimuth = lon_scat - lon_inc;

      auto p_type = get_particle_type();
      if (p_type == ParticleType::Random) {
          auto rot_coeffs = rotation_coefficients(
              lon_inc, lat_inc, lon_scat, lat_scat
              );
          auto scat_ang = rot_coeffs[0];
          auto data = this->interpolate(frequency, temperature, 0.0, 0.0, 0.0, scat_ang);
          if (stokes_dim == 1) {
              CompactFormat<1>::expand_and_transform(
                  result, data,
                  lon_scat, rot_coeffs);
          } else if (stokes_dim == 2) {
              CompactFormat<2>::expand_and_transform(
                  result, data,
                  lon_scat, rot_coeffs);
          } else if (stokes_dim == 3) {
              CompactFormat<3>::expand_and_transform(
                  result, data,
                  lon_scat, rot_coeffs);
          } else if (stokes_dim == 4) {
              CompactFormat<4>::expand_and_transform(
                  result, data,
                  lon_scat, rot_coeffs);
          }
      } else if (p_type == ParticleType::AzimuthallyRandom) {
          auto data = this->interpolate(frequency, temperature, 0, lat_inc, abs(d_azimuth), lat_scat);
          for (Index i = 0; i < stokes_dim * stokes_dim; ++i) {
              result(i / stokes_dim, i % stokes_dim) = data[i];
          }
          if (d_azimuth < 0) {
              if (stokes_dim > 2) {
                  result(0, 2) *= -1.0;
                  result(1, 2) *= -1.0;
                  result(2, 0) *= -1.0;
                  result(2, 1) *= -1.0;
              } if (stokes_dim > 3) {
                  result(0, 2) *= -1.0;
                  result(1, 2) *= -1.0;
                  result(2, 0) *= -1.0;
                  result(2, 1) *= -1.0;
              }
          }
      } else {
          auto data = this->interpolate(frequency, temperature, 0, lat_inc, d_azimuth, lat_scat);
          for (Index i = 0; i < stokes_dim * stokes_dim; ++i) {
              result(i / stokes_dim, i % stokes_dim) = data[i];
          }
      }
      return result;
  }

};

////////////////////////////////////////////////////////////////////////////////
// Extinction matrix
////////////////////////////////////////////////////////////////////////////////

// pxx :: hide
// pxx :: instance(["scattering::ScatteringDataFieldGridded<double>"])
/** Extinction matrix
 *
 * Interface class for extracting scttering matrix data from scattering data field.
 * Similar to the PhaseMatrix class, this class wraps around the ScatteringDataField
 * class that is used to store the raw scatterign coefficients and takes care of
 * assembling the compact representation into the full extinction matrix.
 */
template <typename Base>
class ExtinctionMatrix : public Base {


  using typename Base::Coefficient;

  using typename Base::Vector;
  using typename Base::Matrix;
  using typename Base::DataTensor;
  using typename Base::DataTensorPtr;

  using Base::coeff_dim;
  using Base::data_;
  using Base::f_grid_;
  using Base::t_grid_;
  using Base::lon_inc_;
  using Base::lat_inc_;
  using Base::lon_scat_;
  using Base::lat_scat_;
  using Base::rank;

  using Base::n_temps_;
  using Base::n_freqs_;
  using Base::n_lon_inc_;
  using Base::n_lat_inc_;
  using Base::n_lon_scat_;
  using Base::n_lat_scat_;


  /** Expand extinction matrix coefficients from compact format.
   *
   * @param coefficient Vector-view on the coefficient in compact format.
   * @param type The particle type.
   * @param stokes_dim The required stokes dimension.
   * @param stokes_dim_min The minimum of the required dimension and the stokes
   * dimensions available in the scattering data.
   * @return The expanded extinction matrix.
   */
  template<typename VectorType>
  static Matrix expand_coefficients(const VectorType& coefficients,
                                    ParticleType type,
                                    Index stokes_dim,
                                    Index stokes_dim_min) {
    Matrix matrix = Matrix::Constant(stokes_dim, stokes_dim, 0.0);

    matrix(0, 0) = coefficients(0);
    if (stokes_dim >= 2) {
      matrix(1, 1) = coefficients(0);
    }
    if (stokes_dim >= 3) {
      matrix(2, 2) = coefficients(0);
    }
    if (stokes_dim >= 4) {
      matrix(3, 3) = coefficients(0);
    }

    if (type == ParticleType::AzimuthallyRandom) {
      if (stokes_dim > 1) {
        if (stokes_dim_min > 1) {
          matrix(0, 1) = coefficients(1);
          matrix(1, 0) = coefficients(1);
        }
      }
      if (stokes_dim > 3) {
        if (stokes_dim_min > 2) {
          matrix(2, 3) = coefficients(2);
          matrix(3, 2) = -coefficients(2);
        }
      }
    }
    return matrix;
  }

 public:

  using ExtinctionCoefficient = math::Tensor<Coefficient, rank - 1>;
  using ExtinctionMatrixData = math::Tensor<Coefficient, rank + 1>;
  using Scalar = typename Base::Scalar;

  using Base::copy;

  /// Perfect forwarding constructor.
  template <typename... Args>
  ExtinctionMatrix(Args... args) : Base(std::forward<Args>(args)...) {}

  /// Determine stokes dimension of data.
  math::Index get_stokes_dim() const {
    auto n_coeffs = Base::get_n_coeffs();
    auto type = get_particle_type();
    if (type == ParticleType::Random) {
      return 4;
    } else {
      if (n_coeffs >= 3) {
        return 4;
      }
      if (n_coeffs >= 2) {
        return 2;
      }
      return 1;
    }
  }

  /// Determine stokes dimension of data.
  ParticleType get_particle_type() const {
    if (data_->dimension(coeff_dim) == 1) {
      return ParticleType::Random;
    } else {
      return ParticleType::AzimuthallyRandom;
    }
  }

  /** In-place reduction of scattering data to given stokes dimension.
   *
   * This function discards any stored data that exceeds a given stokes
   * dimension. This can help speed up scattering data preparation
   * in simulation that don't require all stokes components.
   *
   * @param n The maximum required Stokes dimension.
   */
  void set_stokes_dim(math::Index n) {
    auto stokes_dim = std::min(n, get_stokes_dim());
    auto dimensions_new = data_->dimensions();

    if (get_particle_type() == ParticleType::Random) {
      return;
    } else {
      if (stokes_dim == 1) {
        dimensions_new[coeff_dim] = 1;
      } else if (stokes_dim == 2) {
        dimensions_new[coeff_dim] = 2;
      } else if (stokes_dim == 3) {
        dimensions_new[coeff_dim] = 3;
      }
      auto data_new = std::make_shared<DataTensor>(dimensions_new);
      math::copy(*data_new, *data_);
      data_ = data_new;
    }
  }

  /** Get extinction coefficient data in native format.
   *
   * This function returns the extinction coefficient in the format
   * of the underlying scattering data field.
   *
   * @return Tensor containing the extinction coefficient for all frequencies,
   * temperature and angles of the underlying scattering data.
   */
  ExtinctionCoefficient get_extinction_coeff() const {
      return data_->template chip<rank - 1>(0);
  }

  /** Get extinction matrix data.
   *
   * This assembles the extinction matrix data from its compressed
   * form into matrix form.
   *
   * Note: This function does not automatically transform the format of the
   * underlying scattering data field to gridded format.
   *
   * @return The tensor containing the full extinction matrix data with
   * the last dimension expanded into the dimensions of the extinction
   * matrix.
   */
  ExtinctionMatrixData get_extinction_matrix(Index stokes_dim) const {
      std::array<math::Index, rank> dimensions = data_->dimensions();
      dimensions[rank - 1] = stokes_dim * stokes_dim;

      auto stokes_dim_min = std::min(stokes_dim, get_stokes_dim());


      auto new_data = math::Tensor<Coefficient, rank>(dimensions);
      new_data.setZero();

      auto type = get_particle_type();
      if (type == ParticleType::Random) {
          new_data.template chip<rank - 1>(0) = data_->template chip<rank-1>(0);
          if (stokes_dim >= 2) {
              new_data.template chip<rank - 1>(stokes_dim + 1) = data_->template chip<rank-1>(0);
          }
          if (stokes_dim >= 3) {
              new_data.template chip<rank - 1>(2 * stokes_dim + 2) = data_->template chip<rank-1>(0);
          }
          if (stokes_dim >= 4) {
              new_data.template chip<rank - 1>(3 * stokes_dim + 3) = data_->template chip<rank-1>(0);
          }
      } else {
          new_data.template chip<rank - 1>(0) = data_->template chip<rank-1>(0);
          if (stokes_dim >= 2) {
              new_data.template chip<rank - 1>(stokes_dim + 1) = data_->template chip<rank-1>(0);
              if (stokes_dim_min > 1) {
                new_data.template chip<rank - 1>(2) = data_->template chip<rank-1>(1);
                new_data.template chip<rank - 1>(stokes_dim) = data_->template chip<rank-1>(1);
              }
          }
          if (stokes_dim >= 3) {
              new_data.template chip<rank - 1>(2 * stokes_dim + 2) = data_->template chip<rank-1>(0);
          }
          if (stokes_dim >= 4) {
              new_data.template chip<rank - 1>(3 * stokes_dim + 3) = data_->template chip<rank-1>(0);
              if (stokes_dim_min > 2) {
                new_data.template chip<rank - 1>(2 * stokes_dim + 3) =
                    data_->template chip<rank - 1>(2);
                new_data.template chip<rank - 1>(3 * stokes_dim + 2) =
                    -data_->template chip<rank - 1>(2);
              }
          }
      }
      auto output_dimensions = math::get_dimensions<rank + 1>(*data_);
      output_dimensions[rank - 1] = stokes_dim;
      output_dimensions[rank] = stokes_dim;
      return new_data.reshape(output_dimensions);
  }

  /** Calculate extinction matrix.
   *
   * Interpolate and assemble extinction matrix for given frequency
   * temperature and incoming and scattering angles.
   *
   * @param frequency The frequency for which to extract the extinction matrix.
   * @param temperature The temperature for which to extract the extinction
   * matrix.
   * @param lon_inc The longitude comp. of the incoming angle in radians.
   * @param lat_inc The latitude comp. of the incoming angle in radians.
   * @stokes_dim The stokes dimension for which to extract
   *
   * @return Matrix containing the extinction matrix.
   */
  Matrix get_extinction_matrix(Scalar frequency,
                               Scalar temperature,
                               Scalar lon_inc,
                               Scalar lat_inc,
                               Index stokes_dim) const {
    auto stokes_dim_min = std::min(stokes_dim, get_stokes_dim());
    auto data = this->interpolate(temperature,
                                  frequency,
                                  lon_inc,
                                  lat_inc,
                                  0.0,
                                  0.0);
    return expand_coefficients(data,
                               get_particle_type(),
                               stokes_dim,
                               stokes_dim_min);
  }
};

////////////////////////////////////////////////////////////////////////////////
// The absorption vector
////////////////////////////////////////////////////////////////////////////////

/** Absorption vector
 *
 * Interface class to for scattering data fields containing absorption
 * vector data.
 */
template <typename Base>
class AbsorptionVector : public Base {


  using typename Base::Coefficient;

  using typename Base::Vector;
  using typename Base::DataTensor;
  using typename Base::DataTensorPtr;

  using Base::coeff_dim;
  using Base::data_;
  using Base::f_grid_;
  using Base::t_grid_;
  using Base::lon_inc_;
  using Base::lat_inc_;
  using Base::lon_scat_;
  using Base::lat_scat_;
  using Base::rank;

  using Base::n_temps_;
  using Base::n_freqs_;
  using Base::n_lon_inc_;
  using Base::n_lat_inc_;
  using Base::n_lon_scat_;
  using Base::n_lat_scat_;

 public:

  using AbsorptionCoefficient = math::Tensor<Coefficient, rank - 1>;
  using AbsorptionVectorData = math::Tensor<Coefficient, rank>;
  using Scalar = typename Base::Scalar;

  using Base::copy;

  /// Perfect forwarding constructor.
  template <typename... Args>
  AbsorptionVector(Args... args) : Base(std::forward<Args>(args)...) {}

  ParticleType get_particle_type() const {
    if (data_->dimension(coeff_dim) == 1) {
      return ParticleType::Random;
    } else {
      return ParticleType::AzimuthallyRandom;
    }
  }

  /// Determine stokes dimension of data.
  math::Index get_stokes_dim() const {
    auto n_coeffs = Base::get_n_coeffs();
    auto type = get_particle_type();
    if (type == ParticleType::Random) {
      return 4;
    } else {
      if (n_coeffs >= 3) {
        return 4;
      }
      if (n_coeffs >= 2) {
        return 2;
      }
      return 1;
    }
  }

  /** In-place reduction of scattering data to given stokes dimension.
   *
   * This function discards any stored data that exceeds a given stokes
   * dimension. This can help speed up scattering data preparation
   * in simulation that don't require all stokes components.
   *
   * @param n The maximum required Stokes dimension.
   */
  void set_stokes_dim(math::Index n) {
    auto stokes_dim = std::min(n, get_stokes_dim());
    auto dimensions_new = data_->dimensions();

    if (get_particle_type() == ParticleType::Random) {
      return;
    } else {
      if (stokes_dim == 1) {
        dimensions_new[coeff_dim] = 1;
      } else if (stokes_dim >= 2) {
        dimensions_new[coeff_dim] = 2;
      }
      auto data_new = std::make_shared<DataTensor>(dimensions_new);
      math::copy(*data_new, *data_);
      data_ = data_new;
    }
  }

  /** Get absorption coefficient data in native format.
   *
   * This function returns the absorption coefficient in the format
   * of the underlying scattering data field.
   *
   * @return Tensor containing the absorption coefficient for all frequencies,
   * temperature and angles of the underlying scattering data.
   */
  AbsorptionCoefficient get_absorption_coeff() const { return data_->template chip<rank - 1>(0); }

  /** Get absorption vector data.
   *
   * This assembles the absorption vector from its compressed
   * form into vector form.
   *
   * Note: This function does not automatically transform the format of the
   * underlying scattering data field to gridded format.
   *
   * @return A tensor containing the full absorption vector data with
   * the first dimensions corresponding to frequency, temperature and
   * angular grids and the coefficient of the absorption vector along
   * the last dimension.
   */
  AbsorptionVectorData get_absorption_vector(Index stokes_dim) const {
      std::array<math::Index, rank> dimensions = data_->dimensions();
      dimensions[rank - 1] = stokes_dim;
      stokes_dim = std::min(get_stokes_dim(), stokes_dim);

      auto new_data = math::Tensor<Coefficient, rank>(dimensions);
      new_data.setZero();

      new_data.template chip<rank - 1>(0) = data_->template chip<rank-1>(0);
      auto type = get_particle_type();
      if ((type == ParticleType::AzimuthallyRandom) && (stokes_dim > 1)) {
          new_data.template chip<rank - 1>(1) = data_->template chip<rank-1>(1);
      }
      return new_data;
  }

  /** Calculate absorption vector.
   *
   * Interpolate and assemble extinction matrix for given frequency,
   * temperature and incoming and scattering angles.
   *
   * @param frequency The frequency for which to extract the extinction matrix.
   * @param temperature The temperature for which to extract the extinction
   * matrix.
   * @param lon_inc The longitude comp. of the incoming angle in radians.
   * @param lat_inc The latitude comp. of the incoming angle in radians.
   * @param lon_scat The longitude comp. of the scat. angle in radians.
   * @param lat_scat The latitude comp. of the scat. angle in radians.
   * @stokes_dim The stokes dimension for which to extract
   *
   * @return Matrix containing the extinction matrix.
   */
  Vector get_absorption_vector(Scalar frequency,
                               Scalar temperature,
                               Scalar lon_inc,
                               Scalar lat_inc,
                               Index stokes_dim) const {
      auto stokes_dim_min = std::min(stokes_dim, get_stokes_dim());
      auto data = this->interpolate(temperature,
                                    frequency,
                                    lon_inc,
                                    lat_inc,
                                    0.0,
                                    0.0);
      Vector result(stokes_dim);
      result.setZero();
      result[0] = data[0];
      auto type = get_particle_type();
      if ((type == ParticleType::AzimuthallyRandom) && (stokes_dim_min > 1)) {
          result[1] = data[1];
      }
      return result;
  }
};


}  // namespace stokes
}  // namespace scattering
