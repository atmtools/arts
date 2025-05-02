#pragma once

#include <matpack.h>

#include <memory>

#include "arts_constants.h"
#include "interpolation.h"
#include "sht.h"

#include "scattering/utils.h"

namespace scattering {

using sht::SHT;

enum class Format { TRO, ARO, General };
std::ostream &operator<<(std::ostream &out, Format format);

enum class Representation { Gridded, Spectral, DoublySpectral };
std::ostream &operator<<(std::ostream &out, Representation repr);

namespace detail {

template <typename Scalar>
constexpr bool equal(Scalar a, Scalar b, Scalar epsilon = 1e-6) {
  return std::abs(a - b) <=
         ((std::abs(a) > std::abs(b) ? std::abs(a) : std::abs(b)) * epsilon);
}

template <typename Scalar>
constexpr bool small(Scalar a, Scalar epsilon = 1e-6) {
  return std::abs(a) < epsilon;
}

template <typename Scalar>
constexpr Scalar save_acos(Scalar a, Scalar epsilon = 1e-6) {
  if (equal(a, 1.0, epsilon)) {
    return 0.0;
  }
  if (equal(a, -1.0, epsilon)) {
    return Constant::pi;
  }
  return acos(a);
}

/** Calculate angle between incoming and outgoing directions in the scattering
 *  plane.
 *
 * @param aa_inc The incoming-angle azimuth-angle component in degree.
 * @param za_inc The incoming-angle zenith-angle component in degree.
 * @param aa_scat The outgoing (scattering) angle azimuth-angle component in
 * degree.
 * @param za_scat The outgoing (scattering) angle azimuth-angle component in
 * degree.
 * @return The angle between the incoming and outgoing directions in
 * degree.
 */
template <typename Scalar>
Scalar scattering_angle(Scalar aa_inc,
                        Scalar za_inc,
                        Scalar aa_scat,
                        Scalar za_scat) {
  Scalar cos_theta = cos(za_inc) * cos(za_scat) +
                     sin(za_inc) * sin(za_scat) * cos(aa_scat - aa_inc);
  return save_acos(cos_theta);
}

/** Calculate rotation coefficients for scattering matrix.
 *
 * This method calculates the rotation coefficients that are required to
 * transform the scattering matrix of a randomly-oriented particle to the
 * phase matrix, which describes its scattering properties w.r.t. to the
 * laboratory frame. This equation calculates the angle Theta and the
 * coefficients C_1, * C_2, S_1, S_2 as defined in equation (4.16) in
 * "Scattering, Absorption, and Emission of Light by Small Particles."
 *
 * @param aa_inc The azimuth-angle component of the incoming angle in degree.
 * @param za_inc The zenith-angle component of the incoming angle in degree.
 * @param aa_scat The azimuth-angle component of the scattering angle in degree.
 * @param za_scat The zenith-angle component of the scattering angle in degree.
 */
template <typename Scalar>
std::array<Scalar, 5> rotation_coefficients(Scalar aa_inc_d,
                                            Scalar za_inc_d,
                                            Scalar aa_scat_d,
                                            Scalar za_scat_d) {
  Scalar aa_inc  = Conversion::deg2rad(aa_inc_d);
  Scalar za_inc  = Conversion::deg2rad(za_inc_d);
  Scalar aa_scat = Conversion::deg2rad(aa_scat_d);
  Scalar za_scat = Conversion::deg2rad(za_scat_d);

  Scalar cos_theta = cos(za_inc) * cos(za_scat) +
                     sin(za_inc) * sin(za_scat) * cos(aa_scat - aa_inc);
  Scalar theta = save_acos(cos_theta);
  if ((small(abs(aa_scat - aa_inc))) ||
      (equal(abs(aa_scat - aa_inc), 2.0 * pi_v<Scalar>))) {
    theta = abs(za_inc - za_scat);
  } else if ((equal(aa_scat - aa_inc, 2.0 * pi_v<Scalar>))) {
    theta = za_scat + za_inc;
    if (theta > pi_v<Scalar>) {
      theta = 2.0 * pi_v<Scalar> - theta;
    }
  }
  if (small(theta)) {
    return {theta, 1.0, -1.0, 0.0, 0.0};
  } else if (equal(theta, pi_v<Scalar>)) {
    return {theta, 1.0, 1.0, 0.0, 0.0};
  }

  Scalar sigma_1, sigma_2;

  if (small(za_inc)) {
    sigma_1 = aa_scat - aa_inc;
    sigma_2 = 0.0;
  } else if (equal(za_inc, pi_v<Scalar>)) {
    sigma_1 = aa_scat - aa_inc;
    sigma_2 = pi_v<Scalar>;
  } else if (small(za_scat)) {
    sigma_1 = 0.0;
    sigma_2 = pi_v<Scalar> + aa_scat - aa_inc;
  } else if (equal(za_scat, pi_v<Scalar>)) {
    sigma_1 = pi_v<Scalar>;
    sigma_2 = aa_scat - aa_inc;
  } else {
    sigma_1 = save_acos((cos(za_scat) - cos(za_inc) * cos_theta) /
                        (sin(za_inc) * sin(theta)));
    sigma_2 = save_acos((cos(za_inc) - cos(za_scat) * cos_theta) /
                        (sin(za_scat) * sin(theta)));
  }

  Scalar c_1 = cos(2.0 * sigma_1);
  Scalar c_2 = cos(2.0 * sigma_2);
  Scalar s_1 = sin(2.0 * sigma_1);
  Scalar s_2 = sin(2.0 * sigma_2);

  return {theta, c_1, c_2, s_1, s_2};
}

/// [1, 1] element of phase matrix stored in compact format.
template <typename VectorType>
constexpr auto f11(const VectorType &v) {
  return v[0];
}

/// [1, 2] element of phase matrix stored in compact format.
template <typename VectorType>
constexpr auto f12(const VectorType &v) {
  return v[1];
}

/// [2, 2] element of phase matrix stored in compact format.
template <typename VectorType>
constexpr auto f22(const VectorType &v) {
  return v[2];
}

/// [3, 3] element of phase matrix stored in compact format.
template <typename VectorType>
constexpr auto f33(const VectorType &v) {
  return v[3];
}

/// [3, 4] element of phase matrix stored in compact format.
template <typename VectorType>
constexpr auto f34(const VectorType &v) {
  return v[4];
}

/// [4, 4] element of phase matrix stored in compact format.
template <typename VectorType>
constexpr auto f44(const VectorType &v) {
  return v[5];
}


template <typename Scalar, bool strided>
void expand_and_transform(StridedVectorView output,
                          const StridedConstVectorView &input,
                          const std::array<Scalar, 5> rotation_coefficients,
                          bool delta_aa_gt_180) {
  Scalar c_1 = std::get<1>(rotation_coefficients);
  Scalar c_2 = std::get<2>(rotation_coefficients);
  Scalar s_1 = std::get<3>(rotation_coefficients);
  Scalar s_2 = std::get<4>(rotation_coefficients);

  // Stokes dim 1
  output[0] = f11(input);

  // Stokes dim 2 and higher.
  output[1] = c_1 * f12(input);
  output[4] = c_2 * f12(input);
  output[5] = c_1 * c_2 * f22(input) - s_1 * s_2 * f33(input);

  // Stokes dim 3 and higher.
  output[2]  = s_1 * f12(input);
  output[6]  = s_1 * c_2 * f22(input) + c_1 * s_2 * f33(input);
  output[8]  = -s_2 * f12(input);
  output[9]  = -c_1 * s_2 * f22(input) - s_1 * c_2 * f33(input);
  output[10] = -s_1 * s_2 * f22(input) + c_1 * c_2 * f33(input);

  if (delta_aa_gt_180) {
    output[2] *= -1.0;
    output[6] *= -1.0;
    output[8] *= -1.0;
    output[9] *= -1.0;
  }

  // Stokes dim 4 and higher.
  output[3]  = 0.0;
  output[7]  = s_2 * f34(input);
  output[11] = c_2 * f34(input);
  output[12] = 0.0;
  output[13] = s_1 * f34(input);
  output[14] = -c_1 * f34(input);
  output[15] = f44(input);

  if (delta_aa_gt_180) {
    output[7]  *= -1.0;
    output[13] *= -1.0;
  }
}

/** Number of stored phase matrix elements.
 *
 * Returns the number of phase matrix elements that are required for
 * a phase matrix in a given format and stokes dimension.
 *
 * @param format The phase matrix data format.
 * @return The number of elements that are required to be stored.
 */
constexpr Index get_n_mat_elems(Format format) {
  // Compact format used for phase matrix data in TRO format.
  if (format == Format::TRO) {
    return 6;
  }
  // All matrix elements stored for data in other formats.
  return 16;
}

}  // namespace detail

/** Expand phase matrix from compressed coefficient form. */
inline Matrix expand_phase_matrix(const StridedConstVectorView &compact) {
  Matrix mat{4, 4};
  mat[0, 0] = detail::f11(compact);
  mat[0, 1] = detail::f12(compact);
  mat[1, 0] = detail::f12(compact);
  mat[1, 1] = detail::f22(compact);
  mat[2, 2] = detail::f33(compact);
  mat[2, 3] = detail::f34(compact);
  mat[3, 2] = detail::f34(compact);
  mat[3, 3] = detail::f33(compact);
  return mat;
}

inline ComplexMatrix expand_phase_matrix(const StridedConstComplexVectorView &compact) {
  ComplexMatrix mat{4, 4};
  mat[0, 0] = detail::f11(compact);
  mat[0, 1] = detail::f12(compact);
  mat[1, 0] = detail::f12(compact);
  mat[1, 1] = detail::f22(compact);
  mat[2, 2] = detail::f33(compact);
  mat[2, 3] = detail::f34(compact);
  mat[3, 2] = detail::f34(compact);
  mat[3, 3] = detail::f33(compact);
  return mat;
}

/// The grid over which the scattering data is defined.
struct ScatteringDataGrids {
  ScatteringDataGrids(std::shared_ptr<const Vector> t_grid_,
                      std::shared_ptr<const Vector> f_grid_)
      : t_grid(t_grid_),
        f_grid(f_grid_),
        aa_inc_grid(nullptr),
        za_inc_grid(nullptr),
        aa_scat_grid(nullptr),
        za_scat_grid(nullptr) {}

  ScatteringDataGrids(std::shared_ptr<const Vector> t_grid_,
                      std::shared_ptr<const Vector> f_grid_,
                      std::shared_ptr<const ZenithAngleGrid> za_scat_grid_)
      : t_grid(t_grid_),
        f_grid(f_grid_),
        aa_inc_grid(nullptr),
        za_inc_grid(nullptr),
        aa_scat_grid(nullptr),
        za_scat_grid(za_scat_grid_) {}

  ScatteringDataGrids(std::shared_ptr<const Vector> t_grid_,
                      std::shared_ptr<const Vector> f_grid_,
                      std::shared_ptr<const Vector> za_inc_grid_,
                      std::shared_ptr<const Vector> delta_aa_grid_,
                      std::shared_ptr<const ZenithAngleGrid> za_scat_grid_)
      : t_grid(t_grid_),
        f_grid(f_grid_),
        aa_inc_grid(nullptr),
        za_inc_grid(za_inc_grid_),
        aa_scat_grid(delta_aa_grid_),
        za_scat_grid(za_scat_grid_) {}

  std::shared_ptr<const Vector> t_grid;
  std::shared_ptr<const Vector> f_grid;
  std::shared_ptr<const Vector> aa_inc_grid;
  std::shared_ptr<const Vector> za_inc_grid;
  std::shared_ptr<const Vector> aa_scat_grid;
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid;
};

struct RegridWeights {
  ArrayOfGridPos t_grid_weights;
  ArrayOfGridPos f_grid_weights;
  ArrayOfGridPos aa_inc_grid_weights;
  ArrayOfGridPos za_inc_grid_weights;
  ArrayOfGridPos aa_scat_grid_weights;
  ArrayOfGridPos za_scat_grid_weights;
};

inline RegridWeights calc_regrid_weights(
    std::shared_ptr<const Vector> t_grid,
    std::shared_ptr<const Vector> f_grid,
    std::shared_ptr<const Vector> aa_inc_grid,
    std::shared_ptr<const Vector> za_inc_grid,
    std::shared_ptr<const Vector> aa_scat_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
    ScatteringDataGrids new_grids) {
  RegridWeights res{};

  if (!t_grid) {
    ARTS_USER_ERROR(
        "The old t_grid must be provided for calculating regridding weights.");
  }
  if (!new_grids.t_grid) {
    ARTS_USER_ERROR(
        "The new t_grid must be provided for calculating regridding weights.");
  }
  if (!f_grid) {
    ARTS_USER_ERROR(
        "The old f_grid must be provided for calculating regridding weights.");
  }
  if (!new_grids.f_grid) {
    ARTS_USER_ERROR(
        "The new f_grid must be provided for calculating regridding weights.");
  }

  res.t_grid_weights = ArrayOfGridPos(new_grids.t_grid->size());
  gridpos(res.t_grid_weights, *t_grid, *new_grids.t_grid, 1e99);
  res.f_grid_weights = ArrayOfGridPos(new_grids.f_grid->size());
  gridpos(res.f_grid_weights, *f_grid, *new_grids.f_grid, 1e99);

  if ((aa_inc_grid) && (new_grids.aa_inc_grid)) {
    res.aa_inc_grid_weights = ArrayOfGridPos(new_grids.aa_inc_grid->size());
    gridpos(
        res.aa_inc_grid_weights, *aa_inc_grid, *new_grids.aa_inc_grid, 1e99);
  }
  if ((za_inc_grid) && (new_grids.za_inc_grid)) {
    res.za_inc_grid_weights = ArrayOfGridPos(new_grids.za_inc_grid->size());
    gridpos(
        res.za_inc_grid_weights, *za_inc_grid, *new_grids.za_inc_grid, 1e99);
  }
  if ((aa_scat_grid) && (new_grids.aa_scat_grid)) {
    res.aa_scat_grid_weights = ArrayOfGridPos(new_grids.aa_scat_grid->size());
    gridpos(
        res.aa_scat_grid_weights, *aa_scat_grid, *new_grids.aa_scat_grid, 1e99);
  }
  if ((za_scat_grid) && (new_grids.za_scat_grid)) {
    res.za_scat_grid_weights = ArrayOfGridPos(std::visit(
        [](const auto &grd) { return grd.angles.size(); }, *new_grids.za_scat_grid));
    gridpos(res.za_scat_grid_weights,
            std::visit([](const auto &grd) { return static_cast<Vector>(grd.angles); },
                       *za_scat_grid),
            std::visit([](const auto &grd) { return static_cast<Vector>(grd.angles); },
                       *new_grids.za_scat_grid),
            1e99);
  }
  return res;
}

template <std::floating_point Scalar, Format format>
class BackscatterMatrixData : public matpack::data_t<Scalar, 3> {
 public:
  /// The number of stokes coefficients.
  static constexpr Index n_stokes_coeffs = detail::get_n_mat_elems(format);
  using CoeffVector = Eigen::Matrix<Scalar, 1, n_stokes_coeffs>;

  BackscatterMatrixData(std::shared_ptr<const Vector> t_grid,
                        std::shared_ptr<const Vector> f_grid)
      : matpack::data_t<Scalar, 3>(
            t_grid->size(), f_grid->size(), n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid) {
    matpack::data_t<Scalar, 3>::operator=(0.0);
  }

  BackscatterMatrixData &operator=(const matpack::data_t<Scalar, 3> &data) {
    ARTS_USER_ERROR_IF(
        data.shape()[0] != n_temps_,
        "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[1] != n_freqs_,
        "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[2] != n_stokes_coeffs,
        "Provided backscatter coefficient data do not match expected number of stokes coefficients.");
    this->template data_t<Scalar, 3>::operator=(data);
    return *this;
  }

  std::shared_ptr<const Vector> get_t_grid_ptr() const {
    return t_grid_;
  }

  std::shared_ptr<const Vector> get_f_grid_ptr() const {
    return f_grid_;
  }

  constexpr matpack::view_t<CoeffVector, 2> get_coeff_vector_view() {
    return matpack::view_t<CoeffVector, 2>{matpack::mdview_t<CoeffVector, 2>(
        reinterpret_cast<CoeffVector *>(this->data_handle()),
        std::array<Index, 2>{this->extent(0), this->extent(1)})};
  }

  constexpr matpack::view_t<const CoeffVector, 2> get_const_coeff_vector_view() const {
    return matpack::view_t<const CoeffVector, 2>{matpack::mdview_t<const CoeffVector, 2>(
        reinterpret_cast<const CoeffVector *>(this->data_handle()),
        std::array<Index, 2>{this->extent(0), this->extent(1)})};
  }

  BackscatterMatrixData<Scalar, Format::TRO> extract_stokes_coeffs() const {
    constexpr Index n_stokes_coeffs_new = detail::get_n_mat_elems(format);
    BackscatterMatrixData<Scalar, format> result(t_grid_, f_grid_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_s = 0; i_s < n_stokes_coeffs_new; ++i_s) {
          result[i_t, i_f, i_s] = this->operator[](i_t, i_f, i_s);
        }
      }
    }
    return result;
  }

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids,
                               const RegridWeights &weights) const {
    BackscatterMatrixData result(grids.t_grid, grids.f_grid);
    auto coeffs_this = get_const_coeff_vector_view();
    auto coeffs_res  = result.get_coeff_vector_view();
    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos gp_t  = weights.t_grid_weights[i_t];
      Numeric w_t_l = gp_t.fd[1];
      Numeric w_t_r = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos gp_f  = weights.f_grid_weights[i_f];
        Numeric w_f_l = gp_f.fd[1];
        Numeric w_f_r = gp_f.fd[0];
        coeffs_res[i_t, i_f] =
            (w_t_l * w_f_l * coeffs_this[gp_t.idx, gp_f.idx] +
             w_t_l * w_f_r * coeffs_this[gp_t.idx, gp_f.idx + 1] +
             w_t_r * w_f_l * coeffs_this[gp_t.idx + 1, gp_f.idx] +
             w_t_r * w_f_r * coeffs_this[gp_t.idx + 1, gp_f.idx + 1]);
      }
    }
    return result;
  }

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

  BackscatterMatrixData to_gridded() const {
    return *this;
  }

 protected:
  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;
  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;
};  // namespace scattering

template <std::floating_point Scalar>
class BackscatterMatrixData<Scalar, Format::ARO>
    : public matpack::data_t<Scalar, 4> {
 public:
  /// The number of stokes coefficients.
  static constexpr Index n_stokes_coeffs = detail::get_n_mat_elems(Format::ARO);
  using CoeffVector = Eigen::Matrix<Scalar, 1, n_stokes_coeffs>;

  BackscatterMatrixData(std::shared_ptr<const Vector> t_grid,
                        std::shared_ptr<const Vector> f_grid,
                        std::shared_ptr<const Vector> za_inc_grid)
      : matpack::data_t<Scalar, 4>(t_grid->size(),
                                   f_grid->size(),
                                   za_inc_grid->size(),
                                   n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid) {
    matpack::data_t<Scalar, 4>::operator=(0.0);
  }

  BackscatterMatrixData(const BackscatterMatrixData<Scalar, Format::TRO> &bsmat,
                        std::shared_ptr<const Vector> za_inc_grid)
      : matpack::data_t<Scalar, 4>(bsmat.get_t_grid_ptr()->size(),
                                   bsmat.get_f_grid_ptr()->size(),
                                   za_inc_grid->size(),
                                   n_stokes_coeffs),
        n_temps_(bsmat.get_t_grid_ptr()->size()),
        t_grid_(bsmat.get_t_grid_ptr()),
        n_freqs_(bsmat.get_f_grid_ptr()->size()),
        f_grid_(bsmat.get_f_grid_ptr()),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid) {
      for (Index ind = 0; ind < n_za_inc_; ++ind) {
        matpack::data_t<Scalar, 4>::operator[](ind) = bsmat;
      }
  }

  constexpr matpack::view_t<CoeffVector, 3> get_coeff_vector_view() {
    return matpack::view_t<CoeffVector, 3>{matpack::mdview_t<CoeffVector, 3>(
        reinterpret_cast<CoeffVector *>(this->data_handle()),
        std::array<Index, 3>{
            this->extent(0), this->extent(1), this->extent(2)})};
  }

  constexpr matpack::view_t<const CoeffVector, 3> get_const_coeff_vector_view() const {
    return matpack::view_t<const CoeffVector, 3>{matpack::mdview_t<const CoeffVector, 3>(
        reinterpret_cast<const CoeffVector *>(this->data_handle()),
        std::array<Index, 3>{
            this->extent(0), this->extent(1), this->extent(2)})};
  }

  BackscatterMatrixData &operator=(const matpack::data_t<Scalar, 4> &data) {
    ARTS_USER_ERROR_IF(
        data.shape()[0] != n_temps_,
        "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[1] != n_freqs_,
        "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[2] != n_za_inc_,
        "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
    ARTS_USER_ERROR_IF(
        data.shape()[3] != n_stokes_coeffs,
        "Provided backscatter coefficient data do not match expected number of stokes coefficients.");
    this->template data_t<Scalar, 4>::operator=(data);
    return this;
  }

  BackscatterMatrixData<Scalar, Format::ARO> extract_stokes_coeffs() const {
    constexpr Index n_stokes_coeffs_new = detail::get_n_mat_elems(Format::ARO);
    BackscatterMatrixData<Scalar, Format::ARO> result(
        t_grid_, f_grid_, za_inc_grid_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          for (Index i_s = 0; i_s < n_stokes_coeffs_new; ++i_s) {
            result[i_t, i_f, i_za_inc, i_s] =
                this->operator[](i_t, i_f, i_za_inc, i_s);
          }
        }
      }
    }
    return result;
  }

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids,
                               const RegridWeights &weights) const {
    BackscatterMatrixData result(grids.t_grid, grids.f_grid, grids.za_inc_grid);
    auto coeffs_this = get_const_coeff_vector_view();
    auto coeffs_res  = result.get_coeff_vector_view();
    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos gp_t  = weights.t_grid_weights[i_t];
      Numeric w_t_l = gp_t.fd[1];
      Numeric w_t_r = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos gp_f  = weights.f_grid_weights[i_f];
        Numeric w_f_l = gp_f.fd[1];
        Numeric w_f_r = gp_f.fd[0];
        for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size();
             ++i_za_inc) {
          GridPos gp_za_inc  = weights.za_inc_grid_weights[i_za_inc];
          Numeric w_za_inc_l = gp_za_inc.fd[1];
          Numeric w_za_inc_r = gp_za_inc.fd[0];
          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
            coeffs_res[i_t, i_f, i_za_inc] =
                (w_t_l * w_f_l * w_za_inc_l *
                     coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx] +
                 w_t_l * w_f_l * w_za_inc_r *
                     coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx + 1] +
                 w_t_l * w_f_r * w_za_inc_l *
                     coeffs_this[gp_t.idx, gp_f.idx + 1, gp_za_inc.idx] +
                 w_t_l * w_f_r * w_za_inc_r *
                     coeffs_this[gp_t.idx, gp_f.idx + 1, gp_za_inc.idx + 1] +

                 w_t_r * w_f_l * w_za_inc_l *
                     coeffs_this[gp_t.idx + 1, gp_f.idx, gp_za_inc.idx] +
                 w_t_r * w_f_l * w_za_inc_r *
                     coeffs_this[gp_t.idx + 1, gp_f.idx, gp_za_inc.idx + 1] +
                 w_t_r * w_f_r * w_za_inc_l *
                     coeffs_this[gp_t.idx + 1, gp_f.idx + 1, gp_za_inc.idx] +
                 w_t_r * w_f_r * w_za_inc_r *
                     coeffs_this[gp_t.idx + 1,
                                 gp_f.idx + 1,
                                 gp_za_inc.idx + 1]);
          }
        }
      }
    }
    return result;
  }

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, za_inc_grid_, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

 protected:
  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;
  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;
  /// The size of the incoming zenith angle grid.
  Index n_za_inc_;
  /// The incoming zenith angle grid.
  std::shared_ptr<const Vector> za_inc_grid_;
};

template <std::floating_point Scalar, Format format>
using ForwardscatterMatrixData = BackscatterMatrixData<Scalar, format>;

///////////////////////////////////////////////////////////////////////////////
// Phase matrix data
///////////////////////////////////////////////////////////////////////////////

template <std::floating_point Scalar,
          Format format,
          Representation representation>
class PhaseMatrixData;

template <std::floating_point Scalar>
class PhaseMatrixData<Scalar, Format::TRO, Representation::Gridded>
    : public matpack::data_t<Scalar, 4> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::data_t<Scalar, 4>::resize;
  using matpack::data_t<Scalar, 4>::reshape;

 public:
  /// Spectral transform of this phase matrix.
  using PhaseMatrixDataSpectral =
      PhaseMatrixData<Scalar, Format::TRO, Representation::Spectral>;
  using PhaseMatrixDataLabFrame =
      PhaseMatrixData<Scalar, Format::ARO, Representation::Gridded>;

  /// The tensor type used to store the phase matrix data.
  using TensorType = matpack::data_t<Scalar, 4>;

  /// The number of stokes coefficients.
  static constexpr Index n_stokes_coeffs = detail::get_n_mat_elems(Format::TRO);
  using CoeffVector = Eigen::Matrix<Scalar, 1, n_stokes_coeffs>;

  using matpack::data_t<Scalar, 4>::operator[];

  PhaseMatrixData() {}
  /** Create a new PhaseMatrixData container.
   *
   * Creates a container to hold phase matrix data for the
   * provided grids. The phase matrix data in the container is
   * initialized to 0.
   *
   * @param t_grid: A pointer to the temperature grid over which the
   * data is defined.
   * @param f_grid: A pointer to the frequency grid over which the
   * data is defined.
   * @param za_scat_grid: A pointer to the scattering zenith-angle grid
   * over which the data is defined.
   *
   */
  PhaseMatrixData(std::shared_ptr<const Vector> t_grid,
                  std::shared_ptr<const Vector> f_grid,
                  std::shared_ptr<const ZenithAngleGrid> za_scat_grid)
      : matpack::data_t<Scalar, 4>(t_grid->size(),
                                   f_grid->size(),
                                   grid_size(*za_scat_grid),
                                   n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_scat_(grid_size(*za_scat_grid)),
        za_scat_grid_(za_scat_grid) {
    TensorType::operator=(0.0);
  }

  template <typename OtherScalar, Format format, Representation repr>
  PhaseMatrixData(const PhaseMatrixData<OtherScalar, format, repr> &other) {
    // Other must be TRO particle
    if (format != Format::TRO) {
      ARTS_USER_ERROR(
          "Phase matrix data in TRO format can only be constructed "
          " from data that is also in TRO format.");
    }

    // Extract required stokes parameters.
    auto other_stokes = other.extract_stokes_coeffs();

    if constexpr (repr == Representation::Gridded) {
      t_grid_       = other_stokes.get_t_grid();
      f_grid_       = other_stokes.get_f_grid();
      za_scat_grid_ = other_stokes.get_za_scat_grid();
      TensorType::operator=(other_stokes);
    } else {
      auto other_gridded = other.to_gridded();
      t_grid_            = other_gridded.get_t_grid();
      f_grid_            = other_gridded.get_f_grid();
      za_scat_grid_      = other_gridded.get_za_scat_grid();
      TensorType::operator=(other_gridded);
    }
    n_temps_   = t_grid_->size();
    n_freqs_   = f_grid_->size();
    n_za_scat_ = grid_size(*za_scat_grid_);
  }

  PhaseMatrixData &operator=(const matpack::data_t<Scalar, 4> &data) {
    ARTS_USER_ERROR_IF(
        data.shape()[0] != n_temps_,
        "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[1] != n_freqs_,
        "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[2] != n_za_scat_,
        "Provided backscatter coefficient data do not match expected number of scattering zenith angles.");
    ARTS_USER_ERROR_IF(
        data.shape()[3] != n_stokes_coeffs,
        "Provided backscatter coefficient data do not match expected number of stokes coefficients.");
    this->template data_t<Scalar, 4>::operator=(data);
    return *this;
  }

  std::shared_ptr<const Vector> get_t_grid() const { return t_grid_; }
  std::shared_ptr<const Vector> get_f_grid() const { return f_grid_; }
  std::shared_ptr<const ZenithAngleGrid> get_za_scat_grid() const {
    return za_scat_grid_;
  }

  /** Transform phase matrix to spectral format.
   *
   * @param Pointer to the SHT to use for the transformation.
   */
  PhaseMatrixDataSpectral to_spectral(std::shared_ptr<SHT> sht) const {
    ARTS_ASSERT(sht->get_n_azimuth_angles() == 1);

    // Regrid phase matrix along zenith-angles to ensure that it is on the grid
    // expected by SHT.
    if (!std::holds_alternative<FejerGrid>(*za_scat_grid_)) {
      auto new_grids = ScatteringDataGrids(t_grid_,
                                           f_grid_,
                                           std::make_shared<ZenithAngleGrid>(sht->get_zenith_angle_grid()));
      return regrid(new_grids).to_spectral(sht);
    }

    ARTS_ASSERT(sht->get_n_zenith_angles() == n_za_scat_);
    PhaseMatrixDataSpectral result(t_grid_, f_grid_, sht);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
          result[i_t, i_f, joker, i_s] =
              sht->transform(this->operator[](i_t, Range(i_f, 1), joker, i_s));
        }
      }
    }
    return result;
  }

  PhaseMatrixDataSpectral to_spectral(Index degree, Index order) const {
    auto sht_ptr = sht::provider.get_instance_lm(degree, order);
    return to_spectral(sht_ptr);
  }

  PhaseMatrixDataSpectral to_spectral() const {
    return to_spectral(sht::provider.get_instance(1, n_za_scat_));
  }

  PhaseMatrixDataLabFrame to_lab_frame(
      std::shared_ptr<const Vector> za_inc_grid,
      std::shared_ptr<const Vector> delta_aa_grid,
      std::shared_ptr<const ZenithAngleGrid> za_scat_grid_new) const {
    PhaseMatrixDataLabFrame result(
        t_grid_, f_grid_, za_inc_grid, delta_aa_grid, za_scat_grid_new);

    for (Size i_za_inc = 0; i_za_inc < za_inc_grid->size(); ++i_za_inc) {
      GridPos angle_interp;
      for (Size i_delta_aa = 0; i_delta_aa < delta_aa_grid->size();
           ++i_delta_aa) {
        for (Index i_za_scat = 0; i_za_scat < grid_size(*za_scat_grid_new);
             ++i_za_scat) {
          std::array<Scalar, 5> coeffs = detail::rotation_coefficients(
              0.0,
              (*za_inc_grid)[i_za_inc],
              (*delta_aa_grid)[i_delta_aa],
              grid_vector(*za_scat_grid_new)[i_za_scat]);

          // On the fly interpolation of stokes components and expansion.
          Scalar scat_angle = Conversion::rad2deg(std::get<0>(coeffs));
          gridpos(angle_interp, grid_vector(*za_scat_grid_), scat_angle, 1e99);

          Tensor3 scat_mat_interpd(n_temps_, n_freqs_, 4);
          Vector pm_comps(n_stokes_coeffs);
          for (Index i_t = 0; i_t < n_temps_; ++i_t) {
            for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
              for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
                pm_comps[i_s] =
                    (angle_interp.fd[1] *
                         this->operator[](i_t, i_f, angle_interp.idx, i_s) +
                     angle_interp.fd[0] *
                         this->operator[](i_t, i_f, angle_interp.idx + 1, i_s));
              }
              detail::expand_and_transform<Scalar, false>(
                  result[i_t, i_f, i_za_inc, i_delta_aa, i_za_scat, joker],
                  pm_comps,
                  coeffs,
                  (*delta_aa_grid)[i_delta_aa] > 180.0);
            }
          }
        }
      }
    }
    return result;
  }

  BackscatterMatrixData<Scalar, Format::TRO> extract_backscatter_matrix() {
    BackscatterMatrixData<Scalar, Format::TRO> result(t_grid_, f_grid_);
    GridPos interp = find_interp_weights(grid_vector(*za_scat_grid_), 180.0);
    //gridpos(interp, grid_vector(*za_scat_grid_), 180.0);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
          result[i_t, i_f, i_s] =
              (interp.fd[1] * this->operator[](i_t, i_f, interp.idx, i_s) +
               interp.fd[0] * this->operator[](i_t, i_f, interp.idx + 1, i_s));
        }
      }
    }
    return result;
  }

  ForwardscatterMatrixData<Scalar, Format::TRO>
  extract_forwardscatter_matrix() {
    ForwardscatterMatrixData<Scalar, Format::TRO> result(t_grid_, f_grid_);
    GridPos interp;
    gridpos(interp, grid_vector(*za_scat_grid_), 0.0, 1e99);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
          result[i_t, i_f, i_s] =
              (interp.fd[1] * this->operator[](i_t, i_f, interp.idx, i_s) +
               interp.fd[0] * this->operator[](i_t, i_f, interp.idx + 1, i_s));
        }
      }
    }
    return result;
  }

  constexpr matpack::view_t<CoeffVector, 3> get_coeff_vector_view() {
    return matpack::view_t<CoeffVector, 3>{matpack::mdview_t<CoeffVector, 3>(
        reinterpret_cast<CoeffVector *>(this->data_handle()),
        std::array<Index, 3>{
            this->extent(0), this->extent(1), this->extent(2)})};
  }

  constexpr matpack::view_t<const CoeffVector, 3> get_const_coeff_vector_view() const {
    return matpack::view_t<const CoeffVector, 3>{matpack::mdview_t<const CoeffVector, 3>(
        reinterpret_cast<const CoeffVector *>(this->data_handle()),
        std::array<Index, 3>{
            this->extent(0), this->extent(1), this->extent(2)})};
  }

  /** Calculate scattering-angle integral.
   *
   * Integrates the phase matrix over the scattering angles.
   * @return A Tensor3 containing the integral of the phase matrix
   * data with temperatures along the first axis, frequencies along
   * the second and stokes elements along the third.
   */
  Tensor3 integrate_phase_matrix() {
    Tensor3 results(this->extent(0), this->extent(1), n_stokes_coeffs);
    auto result_vec =
        matpack::view_t<CoeffVector, 2>(matpack::mdview_t<CoeffVector, 2>(
            reinterpret_cast<CoeffVector *>(results.data_handle()),
            std::array<Index, 2>{this->extent(0), this->extent(1)}));
    auto this_vec = get_coeff_vector_view();
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        result_vec[i_t, i_f] =
            2.0 * pi_v<Scalar> *
            integrate_zenith_angle(this_vec[i_t, i_f, joker], *za_scat_grid_);
      }
    }
    return results;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A new phase matrix data object containing only data required
   * for the requested stokes dimensions.
   */
  PhaseMatrixData<Scalar, Format::TRO, Representation::Gridded>
  extract_stokes_coeffs() const {
    PhaseMatrixData<Scalar, Format::TRO, Representation::Gridded> result(
        t_grid_, f_grid_, za_scat_grid_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_scat = 0; i_za_scat < n_za_scat_; ++i_za_scat) {
          for (Index i_s = 0; i_s < result.n_stokes_coeffs; ++i_s) {
            result[i_t, i_f, i_za_scat, i_s] =
                this->operator[](i_t, i_f, i_za_scat, i_s);
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids,
                         const RegridWeights &weights) const {
    PhaseMatrixData result(grids.t_grid, grids.f_grid, grids.za_scat_grid);
    auto coeffs_this = get_const_coeff_vector_view();
    auto coeffs_res  = result.get_coeff_vector_view();
    for (Index i_t = 0; i_t < static_cast<Index>(weights.t_grid_weights.size());
         ++i_t) {
      GridPos gp_t  = weights.t_grid_weights[i_t];
      Numeric w_t_l = gp_t.fd[1];
      Numeric w_t_r = gp_t.fd[0];
      for (Index i_f = 0;
           i_f < static_cast<Index>(weights.f_grid_weights.size());
           ++i_f) {
        GridPos gp_f  = weights.f_grid_weights[i_f];
        Numeric w_f_l = gp_f.fd[1];
        Numeric w_f_r = gp_f.fd[0];
        for (Index i_za_scat = 0;
             i_za_scat <
             static_cast<Index>(weights.za_scat_grid_weights.size());
             ++i_za_scat) {
          GridPos gp_za_scat  = weights.za_scat_grid_weights[i_za_scat];
          Numeric w_za_scat_l = gp_za_scat.fd[1];
          Numeric w_za_scat_r = gp_za_scat.fd[0];

          coeffs_res[i_t, i_f, i_za_scat].setZero();

          if (w_t_l > 0.0) {
            if (w_f_l > 0.0) {
              coeffs_res[i_t, i_f, i_za_scat] +=
                  w_t_l * w_f_l * w_za_scat_l *
                  coeffs_this[gp_t.idx, gp_f.idx, gp_za_scat.idx];
              coeffs_res[i_t, i_f, i_za_scat] +=
                  w_t_l * w_f_l * w_za_scat_r *
                  coeffs_this[gp_t.idx, gp_f.idx, gp_za_scat.idx + 1];
            }
            if (w_f_r > 0.0) {
              coeffs_res[i_t, i_f, i_za_scat] +=
                  w_t_l * w_f_r * w_za_scat_l *
                  coeffs_this[gp_t.idx, gp_f.idx + 1, gp_za_scat.idx];
              coeffs_res[i_t, i_f, i_za_scat] +=
                  w_t_l * w_f_r * w_za_scat_r *
                  coeffs_this[gp_t.idx, gp_f.idx + 1, gp_za_scat.idx + 1];
            }
          }
          if (w_t_r > 0.0) {
            if (w_f_l > 0.0) {
              coeffs_res[i_t, i_f, i_za_scat] +=
                  w_t_r * w_f_l * w_za_scat_l *
                  coeffs_this[gp_t.idx + 1, gp_f.idx, gp_za_scat.idx];
              coeffs_res[i_t, i_f, i_za_scat] +=
                  w_t_r * w_f_l * w_za_scat_r *
                  coeffs_this[gp_t.idx + 1, gp_f.idx, gp_za_scat.idx + 1];
            }
            if (w_f_r > 0.0) {
              coeffs_res[i_t, i_f, i_za_scat] +=
                  w_t_r * w_f_r * w_za_scat_l *
                  coeffs_this[gp_t.idx + 1, gp_f.idx + 1, gp_za_scat.idx];
              coeffs_res[i_t, i_f, i_za_scat] +=
                  w_t_r * w_f_r * w_za_scat_r *
                  coeffs_this[gp_t.idx + 1, gp_f.idx + 1, gp_za_scat.idx + 1];
            }
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, za_scat_grid_, grids);
    return regrid(grids, weights);
  }

 protected:
  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;

  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;

  /// The number of scattering zenith angles.
  Index n_za_scat_;
  /// The zenith angle grid.
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid_;
};

template <std::floating_point Scalar, Representation repr>
class PhaseMatrixData<Scalar, Format::TRO, repr>
    : public matpack::data_t<std::complex<Scalar>, 4> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::data_t<std::complex<Scalar>, 4>::resize;
  using matpack::data_t<std::complex<Scalar>, 4>::reshape;

 public:
  /// Gridded transform of this phase matrix.
  using PhaseMatrixDataGridded =
      PhaseMatrixData<Scalar, Format::TRO, Representation::Gridded>;
  using PhaseMatrixDataSpectral =
      PhaseMatrixData<Scalar, Format::TRO, Representation::Spectral>;
  using PhaseMatrixDataLabFrame =
      PhaseMatrixData<Scalar, Format::ARO, Representation::Gridded>;
  using PhaseMatrixDataDoublySpectral =
      PhaseMatrixData<Scalar, Format::TRO, Representation::DoublySpectral>;
  using TensorType = matpack::data_t<std::complex<Scalar>, 4>;

  /// The number of stokes coefficients.
  static constexpr Index n_stokes_coeffs = detail::get_n_mat_elems(Format::TRO);
  using CoeffVector = Eigen::Matrix<std::complex<Scalar>, 1, n_stokes_coeffs>;

  PhaseMatrixData() {}
  /** Create a new PhaseMatrixData container.
   *
   * Creates a container to hold phase matrix data for the
   * provided grids and SHT transform. The phase matrix data
   * in the container is initialized to 0.
   *
   * @param t_grid: A pointer to the temperature grid over which the
   * data is defined.
   * @param f_grid: A pointer to the frequency grid over which the
   * data is defined.
   * @param za_scat_grid: A pointer to the scattering zenith-angle grid
   * over which the data is defined.
   *
   */
  PhaseMatrixData(std::shared_ptr<const Vector> t_grid,
                  std::shared_ptr<const Vector> f_grid,
                  std::shared_ptr<SHT> sht)
      : TensorType(t_grid->size(),
                   f_grid->size(),
                   sht->get_n_spectral_coeffs(),
                   n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_spectral_coeffs_(sht->get_n_spectral_coeffs()),
        sht_(sht) {
    TensorType::operator=(std::complex<Scalar>(0.0, 0.0));
  }

  template <typename OtherScalar, Format format, Representation other_repr>
  PhaseMatrixData(
      const PhaseMatrixData<OtherScalar, format, other_repr> &other) {
    // Other must be TRO particle
    if (format != Format::TRO) {
      ARTS_USER_ERROR(
          "Phase matrix data in TRO format can only be constructed "
          " from data that is also in TRO format.");
    }

    t_grid_  = other.get_t_grid();
    f_grid_  = other.get_f_grid();
    t_grid_  = other.get_t_grid();
    f_grid_  = other.get_f_grid();
    n_temps_ = t_grid_->size();
    n_freqs_ = f_grid_->size();

    // Extract required stokes parameters.
    auto other_stokes = other.extract_stokes_coeffs();

    if constexpr (other_repr == Representation::Gridded) {
      auto other_spectral = other_stokes.to_spectral();
      TensorType::operator=(other_spectral);
      sht_ = other_spectral.get_sht();
      TensorType::operator=(other_spectral);
    } else {
      sht_ = other_stokes.get_sht();
      TensorType::operator=(other_stokes);
    }
    n_spectral_coeffs_ = sht_->get_n_spectral_coeffs();
  }

  PhaseMatrixData &operator=(
      const matpack::data_t<std::complex<Scalar>, 4> &data) {
    ARTS_USER_ERROR_IF(
        data.shape()[0] != n_temps_,
        "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[1] != n_freqs_,
        "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[2] != n_spectral_coeffs_,
        "Provided backscatter coefficient data do not match expected number of SHT coefficients.");
    ARTS_USER_ERROR_IF(
        data.shape()[3] != n_stokes_coeffs,
        "Provided backscatter coefficient data do not match expected number of stokes coefficients.");
    this->template data_t<std::complex<Scalar>, 4>::operator=(data);
    return *this;
  }

  constexpr matpack::view_t<CoeffVector, 3> get_coeff_vector_view() {
    return matpack::view_t<CoeffVector, 3>{matpack::mdview_t<CoeffVector, 3>(
        reinterpret_cast<CoeffVector *>(this->data_handle()),
        std::array<Index, 3>{
            this->extent(0), this->extent(1), this->extent(2)})};
  }

  constexpr matpack::view_t<const CoeffVector, 3> get_const_coeff_vector_view() const {
    return matpack::view_t<const CoeffVector, 3>{matpack::mdview_t<const CoeffVector, 3>(
        reinterpret_cast<const CoeffVector *>(this->data_handle()),
        std::array<Index, 3>{
            this->extent(0), this->extent(1), this->extent(2)})};
  }

  std::shared_ptr<const Vector> get_t_grid() const { return t_grid_; }
  std::shared_ptr<const Vector> get_f_grid() const { return f_grid_; }
  std::shared_ptr<SHT> get_sht() const { return sht_; }

  /** Transform phase matrix to gridded format.
   *
   * @param Pointer to the SHT to use for the transformation.
   */
  PhaseMatrixDataGridded to_gridded() const {
    ARTS_ASSERT(sht_->get_n_spectral_coeffs() == n_spectral_coeffs_);

    auto za_grid =
        std::make_shared<ZenithAngleGrid>(sht_->get_zenith_angle_grid());
    PhaseMatrixDataGridded result(t_grid_, f_grid_, za_grid);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
          result[i_t, i_f, joker, i_s] =
              sht_->synthesize(this->operator[](i_t, i_f, joker, i_s))[0];
        }
      }
    }
    return result;
  }

  PhaseMatrixDataLabFrame to_lab_frame(
      std::shared_ptr<const Vector> za_inc_grid,
      std::shared_ptr<const Vector> delta_aa_grid,
      std::shared_ptr<const ZenithAngleGrid> za_scat_grid_new) const {
    return to_gridded().to_lab_frame(
        za_inc_grid, delta_aa_grid, za_scat_grid_new);
  }

  BackscatterMatrixData<Scalar, Format::TRO> extract_backscatter_matrix()
      const {
    BackscatterMatrixData<Scalar, Format::TRO> result(t_grid_, f_grid_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
          result[i_t, i_f, i_s] =
              sht_->evaluate(this->operator[](i_t, i_f, joker, i_s),
                             Conversion::deg2rad(0.0),
                             Conversion::deg2rad(180.0));
        }
      }
    }
    return result;
  }

  ForwardscatterMatrixData<Scalar, Format::TRO> extract_forwardscatter_matrix()
      const {
    ForwardscatterMatrixData<Scalar, Format::TRO> result(t_grid_, f_grid_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
          result[i_t, i_f, i_s] =
              sht_->evaluate(this->operator[](i_t, i_f, joker, i_s), 0.0, 0.0);
        }
      }
    }
    return result;
  }

  /** Calculate scattering-angle integral.
   *
   * Integrates the phase matrix over the scattering angles.
   * @return A Tensor3 containing the integral of the phase matrix
   * data with temperatures along the first axis, frequencies along
   * the second and stokes elements along the third.
   */
  Tensor3 integrate_phase_matrix() {
    Tensor3 results(this->extent(0), this->extent(1), n_stokes_coeffs);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
          results[i_t, i_f, i_s] = this->operator[](i_t, i_f, 0, i_s).real() *
                                   sqrt(4.0 * pi_v<Scalar>);
        }
      }
    }
    return results;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A new phase matrix data object containing only data required
   * for the requested stokes dimensions.
   */
  PhaseMatrixData<Scalar, Format::TRO, Representation::Spectral>
  extract_stokes_coeffs() const {
    PhaseMatrixData<Scalar, Format::TRO, Representation::Spectral> result(
        t_grid_, f_grid_, sht_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
          for (Index i_s = 0; i_s < result.n_stokes_coeffs; ++i_s) {
            result[i_t, i_f, i_sht, i_s] =
                this->operator[](i_t, i_f, i_sht, i_s);
          }
        }
      }
    }
    return result;
  }

  /// Conversion from doubly-spectral format to spectral is just a
  /// copy for TRO format.
  PhaseMatrixDataSpectral to_spectral(std::shared_ptr<SHT>) const {
    return *this;
  }

  /// Conversion from spectral to doubly-spectral format is just a copy for TRO format.
  PhaseMatrixDataDoublySpectral to_doubly_spectral(std::shared_ptr<SHT>) {
    return *this;
  }

  PhaseMatrixData to_spectral(Index l_new, Index m_new) const {
    auto sht_new = sht::provider.get_instance_lm(l_new, m_new);
    PhaseMatrixData pm_new(t_grid_, f_grid_, sht_new);
    Index f_grid_size = f_grid_->size();
    Index t_grid_size = t_grid_->size();
    for (Index f_ind = 0; f_ind < f_grid_size; ++f_ind) {
      for (Index t_ind = 0; t_ind < t_grid_size; ++t_ind) {
        for (Index coeff_ind = 0;
             coeff_ind < std::min(this->extent(3), pm_new.extent(3)); ++coeff_ind) {
          pm_new[t_ind, f_ind, coeff_ind] =
            this->operator[](t_ind, f_ind, coeff_ind);
        }
      }
    }
    return pm_new;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids,
                         const RegridWeights weights) const {
    PhaseMatrixData result(grids.t_grid, grids.f_grid, sht_);
    auto coeffs_this = get_const_coeff_vector_view();
    auto coeffs_res  = result.get_coeff_vector_view();
    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos gp_t  = weights.t_grid_weights[i_t];
      Numeric w_t_l = gp_t.fd[1];
      Numeric w_t_r = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos gp_f  = weights.f_grid_weights[i_f];
        Numeric w_f_l = gp_f.fd[1];
        Numeric w_f_r = gp_f.fd[0];

        for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
          coeffs_res[i_t, i_f, i_sht].setZero();
        }

        if (w_t_l > 0.0) {
          if (w_f_l > 0.0) {
            for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
              coeffs_res[i_t, i_f, i_sht] +=
                  w_t_l * w_f_l * coeffs_this[gp_t.idx, gp_f.idx, i_sht];
            }
          }
          if (w_f_r > 0.0) {
            for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
              coeffs_res[i_t, i_f, i_sht] +=
                  w_t_l * w_f_r * coeffs_this[gp_t.idx, gp_f.idx + 1, i_sht];
            }
          }
        }
        if (w_t_r > 0.0) {
          if (w_f_l > 0.0) {
            for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
              coeffs_res[i_t, i_f, i_sht] +=
                  w_t_r * w_f_l * coeffs_this[gp_t.idx + 1, gp_f.idx, i_sht];
            }
          }
          if (w_f_r > 0.0) {
            for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
              coeffs_res[i_t, i_f, i_sht] +=
                  w_t_r * w_f_r *
                  coeffs_this[gp_t.idx + 1, gp_f.idx + 1, i_sht];
            }
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

 protected:
  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;

  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;

  /// Number of SHT coefficients.
  Index n_spectral_coeffs_;

  std::shared_ptr<SHT> sht_;
};

///////////////////////////////////////////////////////////////////////////////
// ARO format
///////////////////////////////////////////////////////////////////////////////

template <std::floating_point Scalar>
class PhaseMatrixData<Scalar, Format::ARO, Representation::Gridded>
    : public matpack::data_t<Scalar, 6> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::data_t<Scalar, 6>::resize;
  using matpack::data_t<Scalar, 6>::reshape;

 public:
  /// Spectral transform of this phase matrix.
  using PhaseMatrixDataSpectral =
      PhaseMatrixData<Scalar, Format::ARO, Representation::Spectral>;

  /// The number of stokes coefficients.
  static constexpr Index n_stokes_coeffs = detail::get_n_mat_elems(Format::ARO);
  using CoeffVector = Eigen::Matrix<Scalar, 1, n_stokes_coeffs>;

  PhaseMatrixData() {}
  /** Create a new PhaseMatrixData container.
   *
   * Creates a container to hold phase matrix data for the
   * provided grids. The phase matrix data in the container is
   * initialized to 0.
   *
   * @param t_grid: A pointer to the temperature grid over which the
   * data is defined.
   * @param f_grid: A pointer to the frequency grid over which the
   * data is defined.
   * @param za_inc_grid: A pointer to the incoming zenith-angle grid
   * over which the data is defined.
   * @param delta_aa_grid: A pointer to the azimuth angle difference
   * grid over which the data is defined.
   * @param za_scat_grid: A pointer to the scattering zenith-angle grid
   * over which the data is defined.
   *
   */
  PhaseMatrixData(std::shared_ptr<const Vector> t_grid,
                  std::shared_ptr<const Vector> f_grid,
                  std::shared_ptr<const Vector> za_inc_grid,
                  std::shared_ptr<const Vector> delta_aa_grid,
                  std::shared_ptr<const ZenithAngleGrid> za_scat_grid)
      : matpack::data_t<Scalar, 6>(t_grid->size(),
                                   f_grid->size(),
                                   za_inc_grid->size(),
                                   delta_aa_grid->size(),
                                   grid_size(*za_scat_grid),
                                   n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid),
        n_delta_aa_(delta_aa_grid->size()),
        delta_aa_grid_(delta_aa_grid),
        n_za_scat_(grid_size(*za_scat_grid)),
        za_scat_grid_(za_scat_grid) {
    matpack::data_t<Scalar, 6>::operator=(0.0);
  }

  constexpr matpack::view_t<CoeffVector, 5> get_coeff_vector_view() {
    return matpack::view_t<CoeffVector, 5>{matpack::mdview_t<CoeffVector, 5>(
        reinterpret_cast<CoeffVector *>(this->data_handle()),
        std::array<Index, 5>{this->extent(0),
                             this->extent(1),
                             this->extent(2),
                             this->extent(3),
                             this->extent(4)})};
  }

  constexpr matpack::view_t<const CoeffVector, 5> get_const_coeff_vector_view() const {
    return matpack::view_t<const CoeffVector, 5>{matpack::mdview_t<const CoeffVector, 5>(
        reinterpret_cast<const CoeffVector *>(this->data_handle()),
        std::array<Index, 5>{this->extent(0),
                             this->extent(1),
                             this->extent(2),
                             this->extent(3),
                             this->extent(4)})};
  }

  PhaseMatrixData &operator=(const matpack::data_t<Scalar, 6> &data) {
    ARTS_USER_ERROR_IF(
        data.shape()[0] != n_temps_,
        "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[1] != n_freqs_,
        "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[2] != n_za_inc_,
        "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
    ARTS_USER_ERROR_IF(
        data.shape()[3] != n_delta_aa_,
        "Provided backscatter coefficient data do not match expected number of scattering azimuth angles.");
    ARTS_USER_ERROR_IF(
        data.shape()[4] != n_za_scat_,
        "Provided backscatter coefficient data do not match expected number of scattering zenith angles.");
    ARTS_USER_ERROR_IF(
        data.shape()[5] != n_stokes_coeffs,
        "Provided backscatter coefficient data do not match expected number of stokes coefficients.");
    this->template data_t<Scalar, 6>::operator=(data);
    return *this;
  }

  /** Transform phase matrix to spectral format.
   *
   * @param Pointer to the SHT to use for the transformation.
   */
  PhaseMatrixDataSpectral to_spectral(std::shared_ptr<SHT> sht) const {
    ARTS_ASSERT(sht->get_n_azimuth_angles() == n_delta_aa_);
    ARTS_ASSERT(sht->get_n_zenith_angles() == n_za_scat_);

    PhaseMatrixDataSpectral result(t_grid_, f_grid_, za_inc_grid_, sht);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
            result[i_t, i_f, i_za_inc, joker, i_s] = sht->transform(
                this->operator[](i_t, i_f, i_za_inc, joker, joker, i_s));
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixDataSpectral to_spectral(Index degree, Index order) const {
    auto sht_ptr = sht::provider.get_instance_lm(degree, order);
    return to_spectral(sht_ptr);
  }

  PhaseMatrixDataSpectral to_spectral() const {
    return to_spectral(sht::provider.get_instance(n_delta_aa_, n_za_scat_));
  }

  BackscatterMatrixData<Scalar, Format::ARO> extract_backscatter_matrix() {
    BackscatterMatrixData<Scalar, Format::ARO> result(
        t_grid_, f_grid_, za_inc_grid_);
    GridPos za_scat_interp, delta_aa_interp;
    gridpos(delta_aa_interp, *delta_aa_grid_, 180.0, 1e99);
    Vector weights(4);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          auto za_inc = (*za_inc_grid_)[i_za_inc];
          gridpos(za_scat_interp,
                  grid_vector(*za_scat_grid_),
                  180.0 - za_inc,
                  1e99);
          interpweights(weights, delta_aa_interp, za_scat_interp);

          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
            auto mat = this->operator[](i_t, i_f, i_za_inc, joker, joker, i_s);
            result[i_t, i_f, i_za_inc, i_s] =
                interp(weights, mat, delta_aa_interp, za_scat_interp);
          }
        }
      }
    }
    return result;
  }

  ForwardscatterMatrixData<Scalar, Format::ARO>
  extract_forwardscatter_matrix() {
    BackscatterMatrixData<Scalar, Format::ARO> result(
        t_grid_, f_grid_, za_inc_grid_);
    GridPos za_scat_interp, delta_aa_interp;
    gridpos(delta_aa_interp, *delta_aa_grid_, 0.0, 1e99);
    Vector weights(4);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          auto za_inc = (*za_inc_grid_)[i_za_inc];
          gridpos(za_scat_interp, grid_vector(*za_scat_grid_), za_inc, 1e99);
          interpweights(weights, delta_aa_interp, za_scat_interp);

          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
            auto mat = this->operator[](i_t, i_f, i_za_inc, joker, joker, i_s);
            result[i_t, i_f, i_za_inc, i_s] =
                interp(weights, mat, delta_aa_interp, za_scat_interp);
          }
        }
      }
    }
    return result;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A new phase matrix data object containing only data required
   * for the requested stokes dimensions.
   */
  PhaseMatrixData<Scalar, Format::ARO, Representation::Gridded>
  extract_stokes_coeffs() const {
    PhaseMatrixData<Scalar, Format::ARO, Representation::Gridded> result(
        t_grid_, f_grid_, za_inc_grid_, delta_aa_grid_, za_scat_grid_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          for (Index i_aa_scat = 0; i_aa_scat < n_delta_aa_; ++i_aa_scat) {
            for (Index i_za_scat = 0; i_za_scat < n_za_scat_; ++i_za_scat) {
              for (Index i_s = 0; i_s < result.n_stokes_coeffs; ++i_s) {
                result[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat, i_s] =
                    this->operator[](
                        i_t, i_f, i_za_inc, i_aa_scat, i_za_scat, i_s);
              }
            }
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids,
                         const RegridWeights &weights) const {
    PhaseMatrixData result(grids.t_grid,
                           grids.f_grid,
                           grids.za_inc_grid,
                           grids.aa_scat_grid,
                           grids.za_scat_grid);
    auto coeffs_this = get_const_coeff_vector_view();
    auto coeffs_res  = result.get_coeff_vector_view();

    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos gp_t  = weights.t_grid_weights[i_t];
      Numeric w_t_l = gp_t.fd[1];
      Numeric w_t_r = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos gp_f  = weights.f_grid_weights[i_f];
        Numeric w_f_l = gp_f.fd[1];
        Numeric w_f_r = gp_f.fd[0];
        for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size();
             ++i_za_inc) {
          GridPos gp_za_inc  = weights.za_inc_grid_weights[i_za_inc];
          Numeric w_za_inc_l = gp_za_inc.fd[1];
          Numeric w_za_inc_r = gp_za_inc.fd[0];
          for (Size i_aa_scat = 0;
               i_aa_scat < weights.aa_scat_grid_weights.size();
               ++i_aa_scat) {
            GridPos gp_aa_scat  = weights.aa_scat_grid_weights[i_aa_scat];
            Numeric w_aa_scat_l = gp_aa_scat.fd[1];
            Numeric w_aa_scat_r = gp_aa_scat.fd[0];
            for (Size i_za_scat = 0;
                 i_za_scat < weights.za_scat_grid_weights.size();
                 ++i_za_scat) {
              GridPos gp_za_scat  = weights.za_scat_grid_weights[i_za_scat];
              Numeric w_za_scat_l = gp_za_scat.fd[1];
              Numeric w_za_scat_r = gp_za_scat.fd[0];

              coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat].setZero();

              if (w_t_l > 0.0) {
                if (w_f_l > 0.0) {
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx + 1];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx + 1];

                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx + 1];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx + 1];
                }
                if (w_f_r > 0.0) {
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx + 1];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx + 1];

                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx + 1];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx + 1];
                }
              }
              if (w_t_r > 0.0) {
                if (w_f_l > 0.0) {
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx + 1];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx + 1];

                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx + 1];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx + 1];
                }
                if (w_f_r > 0.0) {
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx + 1];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx + 1];

                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx,
                                  gp_za_scat.idx + 1];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx + 1,
                                  gp_f.idx + 1,
                                  gp_za_inc.idx + 1,
                                  gp_aa_scat.idx + 1,
                                  gp_za_scat.idx + 1];
                }
              }
            }
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const {
    auto weights = calc_regrid_weights(t_grid_,
                                       f_grid_,
                                       nullptr,
                                       za_inc_grid_,
                                       delta_aa_grid_,
                                       za_scat_grid_,
                                       grids);
    return regrid(grids, weights);
  }

 protected:
  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;

  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;

  /// The number of incoming zenith angles.
  Index n_za_inc_;
  /// The incoming angle grid.
  std::shared_ptr<const Vector> za_inc_grid_;

  /// The number of angles in the azimuth difference grid.
  Index n_delta_aa_;
  /// The azimuth difference grid.
  std::shared_ptr<const Vector> delta_aa_grid_;

  /// The number of scattering zenith angles.
  Index n_za_scat_;
  /// The zenith angle grid.
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid_;
};

template <std::floating_point Scalar>
class PhaseMatrixData<Scalar, Format::ARO, Representation::Spectral>
    : public matpack::data_t<std::complex<Scalar>, 5> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::data_t<std::complex<Scalar>, 5>::resize;
  using matpack::data_t<std::complex<Scalar>, 5>::reshape;

 public:
  /// Spectral transform of this phase matrix.
  using PhaseMatrixDataGridded =
      PhaseMatrixData<Scalar, Format::ARO, Representation::Gridded>;

  /// The number of stokes coefficients.
  static constexpr Index n_stokes_coeffs = detail::get_n_mat_elems(Format::ARO);
  using CoeffVector = Eigen::Matrix<std::complex<Scalar>, 1, n_stokes_coeffs>;

  PhaseMatrixData() {}
  /** Create a new PhaseMatrixData container.
   *
   * Creates a container to hold phase matrix data for the
   * provided grids. The phase matrix data in the container is
   * initialized to 0.
   *
   * @param t_grid: A pointer to the temperature grid over which the
   * data is defined.
   * @param f_grid: A pointer to the frequency grid over which the
   * data is defined.
   * @param za_inc_grid: A pointer to the incoming zenith-angle grid
   * over which the data is defined.
   * @param sht: A shared pointer to the SHT object used to transform
   * the phase matrix data.
   */
  PhaseMatrixData(std::shared_ptr<const Vector> t_grid,
                  std::shared_ptr<const Vector> f_grid,
                  std::shared_ptr<const Vector> za_inc_grid,
                  std::shared_ptr<SHT> sht)
      : matpack::data_t<std::complex<Scalar>, 5>(t_grid->size(),
                                                 f_grid->size(),
                                                 za_inc_grid->size(),
                                                 sht->get_n_spectral_coeffs(),
                                                 n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid),
        n_spectral_coeffs_(sht->get_n_spectral_coeffs()),
        sht_(sht) {
    matpack::data_t<std::complex<Scalar>, 5>::operator=(0.0);
  }

  PhaseMatrixData &operator=(
      const matpack::data_t<std::complex<Scalar>, 5> &data) {
    ARTS_USER_ERROR_IF(
        data.shape()[0] != n_temps_,
        "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[1] != n_freqs_,
        "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[2] != n_za_inc_,
        "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
    ARTS_USER_ERROR_IF(
        data.shape()[3] != n_spectral_coeffs_,
        "Provided backscatter coefficient data do not match expected number of SHT coefficients.");
    ARTS_USER_ERROR_IF(
        data.shape()[4] != n_stokes_coeffs,
        "Provided backscatter coefficient data do not match expected number of stokes coefficients.");
    this->template data_t<std::complex<Scalar>, 4>::operator=(data);
    return *this;
  }

  constexpr matpack::view_t<CoeffVector, 4> get_coeff_vector_view() {
    return matpack::view_t<CoeffVector, 4>{matpack::mdview_t<CoeffVector, 4>(
        reinterpret_cast<CoeffVector *>(this->data_handle()),
        std::array<Index, 4>{this->extent(0),
                             this->extent(1),
                             this->extent(2),
                             this->extent(3)})};
  }

  constexpr matpack::view_t<const CoeffVector, 4> get_const_coeff_vector_view() const {
    return matpack::view_t<const CoeffVector, 4>{matpack::mdview_t<const CoeffVector, 4>(
        reinterpret_cast<const CoeffVector *>(this->data_handle()),
        std::array<Index, 4>{this->extent(0),
                             this->extent(1),
                             this->extent(2),
                             this->extent(3)})};
  }

  /** Transform phase matrix to gridded format.
   *
   * @param Pointer to the SHT to use for the transformation.
   */
  PhaseMatrixDataGridded to_gridded() const {
    PhaseMatrixDataGridded result(t_grid_,
                                  f_grid_,
                                  za_inc_grid_,
                                  sht_->get_aa_grid_ptr(),
                                  sht_->get_za_grid_ptr());

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
            result[i_t, i_f, i_za_inc, joker, joker, i_s] = sht_->synthesize(
                this->operator[](i_t, i_f, i_za_inc, joker, i_s));
          }
        }
      }
    }
    return result;
  }

  /** Transform phase matixr to spectral format.
   *
   * @param Pointer to the SHT to use for the transformation.
   */
  PhaseMatrixData to_spectral(Index l_new, Index m_new) const {
    auto sht_new = sht::provider.get_instance_lm(l_new, m_new);
    PhaseMatrixData pm_new(t_grid_, f_grid_, za_inc_grid_, sht_new);
    for (Index f_ind = 0; f_ind < f_grid_->size(); ++f_ind) {
      for (Index t_ind = 0; t_ind < t_grid_->size(); ++t_ind) {
        for (Index za_inc_ind = 0; za_inc_ind < za_inc_grid_->size();
             ++za_inc_ind) {
          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {

          for (Index coeff_ind = 0; coeff_ind < std::min(this->extent(4), pm_new.extent(4)); ++coeff_ind) {
            pm_new[t_ind, f_ind, za_inc_ind, joker, i_s] = sht::add_coeffs(*sht_new,
                                                                           pm_new[t_ind, f_ind, za_inc_ind, joker, i_s],
                                                                           *sht_,
                                                                           (*this)[t_ind, f_ind, za_inc_ind, joker, i_s]);
            }
          }
        }
      }
    }
    return pm_new;
  }

  BackscatterMatrixData<Scalar, Format::ARO> extract_backscatter_matrix() const {
    BackscatterMatrixData<Scalar, Format::ARO> result(
        t_grid_, f_grid_, za_inc_grid_);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          auto za_inc    = (*za_inc_grid_)[i_za_inc];
          Scalar za_scat = 180.0 - za_inc;
          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
            auto coeffs = this->operator[](i_t, i_f, i_za_inc, joker, i_s);
            result[i_t, i_f, i_za_inc, i_s] =
                sht_->evaluate(coeffs,
                               Conversion::deg2rad(180.0),
                               Conversion::deg2rad(za_scat));
          }
        }
      }
    }
    return result;
  }

  ForwardscatterMatrixData<Scalar, Format::ARO>
  extract_forwardscatter_matrix() {
    BackscatterMatrixData<Scalar, Format::ARO> result(
        t_grid_, f_grid_, za_inc_grid_);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          auto za_inc = (*za_inc_grid_)[i_za_inc];
          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
            auto coeffs = this->operator[](i_t, i_f, i_za_inc, joker, i_s);
            result[i_t, i_f, i_za_inc, i_s] = sht_->evaluate(
                coeffs, Conversion::deg2rad(0.0), Conversion::deg2rad(za_inc));
          }
        }
      }
    }
    return result;
  }

  Tensor4 integrate_phase_matrix() {
    Tensor4 results(n_temps_, n_freqs_, n_za_inc_, n_stokes_coeffs);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          for (Index i_s = 0; i_s < n_stokes_coeffs; ++i_s) {
            results[i_t, i_f, i_za_inc, i_s] =
                this->operator[](i_t, i_f, i_za_inc, 0, i_s).real() *
                sqrt(4.0 * pi_v<Scalar>);
          }
        }
      }
    }
    return results;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A new phase matrix data object containing only data required
   * for the requested stokes dimensions.
   */
  PhaseMatrixData<Scalar, Format::ARO, Representation::Gridded>
  extract_stokes_coeffs() const {
    PhaseMatrixData<Scalar, Format::ARO, Representation::Spectral> result(
        t_grid_, f_grid_, za_inc_grid_, sht_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
            for (Index i_s = 0; i_s < result.n_stokes_coeffs; ++i_s) {
              result[i_t, i_f, i_za_inc, i_sht, i_s] =
                  this->operator[](i_t, i_f, i_za_inc, i_sht, i_s);
            }
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids,
                         const RegridWeights weights) const {
    PhaseMatrixData result(grids.t_grid, grids.f_grid, grids.za_inc_grid, sht_);
    auto coeffs_this = get_const_coeff_vector_view();
    auto coeffs_res  = result.get_coeff_vector_view();

    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos gp_t  = weights.t_grid_weights[i_t];
      Numeric w_t_l = gp_t.fd[1];
      Numeric w_t_r = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos gp_f  = weights.f_grid_weights[i_f];
        Numeric w_f_l = gp_f.fd[1];
        Numeric w_f_r = gp_f.fd[0];
        for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size();
             ++i_za_inc) {
          GridPos gp_za_inc  = weights.za_inc_grid_weights[i_za_inc];
          Numeric w_za_inc_l = gp_za_inc.fd[1];
          Numeric w_za_inc_r = gp_za_inc.fd[0];

          for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
            coeffs_res[i_t, i_f, i_za_inc, i_sht].setZero();
          }

          if (w_t_l > 0.0) {
            if (w_f_l > 0.0) {
              for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_l * w_f_l * w_za_inc_l *
                    coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, i_sht];
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_l * w_f_l * w_za_inc_r *
                    coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx + 1, i_sht];
              }
            }
            if (w_f_r > 0.0) {
              for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_l * w_f_r * w_za_inc_l *
                    coeffs_this[gp_t.idx, gp_f.idx + 1, gp_za_inc.idx, i_sht];
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_l * w_f_r * w_za_inc_r *
                    coeffs_this
                        [gp_t.idx, gp_f.idx + 1, gp_za_inc.idx + 1, i_sht];
              }
            }
          }
          if (w_t_r > 0.0) {
            if (w_f_l > 0.0) {
              for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_r * w_f_l * w_za_inc_l *
                    coeffs_this[gp_t.idx + 1, gp_f.idx, gp_za_inc.idx, i_sht];
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_r * w_f_l * w_za_inc_r *
                    coeffs_this
                        [gp_t.idx + 1, gp_f.idx, gp_za_inc.idx + 1, i_sht];
              }
            }
            if (w_f_r > 0.0) {
              for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_r * w_f_r * w_za_inc_l *
                    coeffs_this
                        [gp_t.idx + 1, gp_f.idx + 1, gp_za_inc.idx, i_sht];
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_r * w_f_r * w_za_inc_r *
                    coeffs_this
                        [gp_t.idx + 1, gp_f.idx + 1, gp_za_inc.idx + 1, i_sht];
              }
            }
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const {
    auto weights = calc_regrid_weights(t_grid_,
                                       f_grid_,
                                       nullptr,
                                       za_inc_grid_,
                                       nullptr,
                                       nullptr,
                                       grids);
    return regrid(grids, weights);
  }

 protected:
  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;
  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;
  /// The number of incoming zenith angles.
  Index n_za_inc_;
  /// The incoming zenith angle grid.
  std::shared_ptr<const Vector> za_inc_grid_;
  /// The number of SHT coefficients.
  Index n_spectral_coeffs_;
  /// The incoming zenith angle grid.
  std::shared_ptr<SHT> sht_;
};

}  // namespace scattering
