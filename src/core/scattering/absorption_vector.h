#ifndef SCATTERING_ABSORPTION_VECTOR_H_
#define SCATTERING_ABSORPTION_VECTOR_H_

#include <matpack.h>

#include <memory>

#include "phase_matrix.h"
#include "sht.h"

namespace scattering {

namespace absorption {

/** Number of stored extinction matrix elements.
 *
 * Returns the number of extinction matrix elements that stored for
 * an extinction matrix for a given format and stokes dimension.
 *
 * @param format The extinction matrix data format.
 * @return The number of elements that are required to be stored.
 */
constexpr Index get_n_mat_elems(Format format) {
  // Compact format used for extinction matrix data in TRO format.
  if (format == Format::TRO) {
    return 1;
  }
  return 2;
}

}  // namespace absorption

template <std::floating_point Scalar,
          Format format,
          Representation representation>
class AbsorptionVectorData;

template <std::floating_point Scalar, Representation repr>
class AbsorptionVectorData<Scalar, Format::TRO, repr>
    : public matpack::data_t<Scalar, 3> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::data_t<Scalar, 3>::resize;
  using matpack::data_t<Scalar, 3>::reshape;

 public:
  using AbsorptionVectorDataLabFrame =
      AbsorptionVectorData<Scalar, Format::ARO, Representation::Gridded>;

  constexpr static Index n_stokes_coeffs =
      absorption::get_n_mat_elems(Format::TRO);
  using CoeffVector = Eigen::Matrix<Scalar, 1, n_stokes_coeffs>;

  using matpack::data_t<Scalar, 3>::operator[];

  AbsorptionVectorData() {};

  /** Create a new AbsorptionVectorData container.
   *
   * Creates a container to hold extinction matrix data for the
   * provided grids. The extinction matrix data in the container is
   * initialized to 0.
   *
   * @param t_grid: A pointer to the temperature grid over which the
   * data is defined.
   * @param f_grid: A pointer to the frequency grid over which the
   * data is defined.
   */
  AbsorptionVectorData(std::shared_ptr<const Vector> t_grid,
                       std::shared_ptr<const Vector> f_grid)
      : matpack::data_t<Scalar, 3>(
            t_grid->size(), f_grid->size(), n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid) {
    matpack::data_t<Scalar, 3>::operator=(0.0);
  }

  AbsorptionVectorData& operator=(const matpack::data_t<Scalar, 3>& data) {
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

  constexpr matpack::view_t<CoeffVector, 2> get_coeff_vector_view() {
    return matpack::view_t<CoeffVector, 2>{matpack::mdview_t<CoeffVector, 2>(
        reinterpret_cast<CoeffVector*>(this->data_handle()),
        std::array<Index, 2>{this->extent(0), this->extent(1)})};
  }

  constexpr matpack::view_t<const CoeffVector, 2> get_const_coeff_vector_view() const {
    return matpack::view_t<const CoeffVector, 2>{matpack::mdview_t<const CoeffVector, 2>(
        reinterpret_cast<const CoeffVector*>(this->data_handle()),
        std::array<Index, 2>{this->extent(0), this->extent(1)})};
  }

  AbsorptionVectorData<Scalar, Format::ARO, repr> to_lab_frame(
      std::shared_ptr<const Vector> za_inc_grid) const {
    AbsorptionVectorData<Scalar, Format::ARO, repr> av_new{
        t_grid_, f_grid_, za_inc_grid};
    for (Size t_ind = 0; t_ind < t_grid_->size(); ++t_ind) {
      for (Size f_ind = 0; f_ind < f_grid_->size(); ++f_ind) {
        for (Size za_inc_ind = 0; za_inc_ind < za_inc_grid->size();
             ++za_inc_ind) {
          av_new[t_ind, f_ind, za_inc_ind, 0] =
              this->operator[](t_ind, f_ind, 0);
        }
      }
    }
    return av_new;
  }

  AbsorptionVectorData<Scalar, Format::TRO, Representation::Spectral>
  to_spectral() const {
    AbsorptionVectorData<Scalar, Format::TRO, Representation::Spectral> avd_new(
        t_grid_, f_grid_);
    reinterpret_cast<matpack::data_t<Scalar, 3>&>(avd_new) = *this;
    return avd_new;
  }

  AbsorptionVectorData regrid(const ScatteringDataGrids grids,
                              const RegridWeights weights) const  {
    AbsorptionVectorData result(grids.t_grid, grids.f_grid);
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
        coeffs_res[i_t, i_f] =
            (w_t_l * w_f_l * coeffs_this[gp_t.idx, gp_f.idx] +
             w_t_l * w_f_r * coeffs_this[gp_t.idx, gp_f.idx + 1] +
             w_t_r * w_f_l * coeffs_this[gp_t.idx + 1, gp_f.idx] +
             w_t_r * w_f_r * coeffs_this[gp_t.idx + 1, gp_f.idx + 1]);
      }
    }
    return result;
  }

  AbsorptionVectorData regrid(const ScatteringDataGrids& grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

  AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded> to_gridded() const {
    auto abs_new = AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded>(t_grid_, f_grid_);
    abs_new = *this;
    return abs_new;
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
};

template <std::floating_point Scalar, Representation repr>
class AbsorptionVectorData<Scalar, Format::ARO, repr>
    : public matpack::data_t<Scalar, 4> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::data_t<Scalar, 4>::resize;
  using matpack::data_t<Scalar, 4>::reshape;

 public:
  constexpr static Index n_stokes_coeffs =
      absorption::get_n_mat_elems(Format::ARO);
  using CoeffVector = Eigen::Matrix<Scalar, 1, n_stokes_coeffs>;

  AbsorptionVectorData() {};
  /** Create a new AbsorptionVectorData container.
   *
   * Creates a container to hold extinction matrix data for the
   * provided grids. The extinction matrix data in the container is
   * initialized to 0.
   *
   * @param t_grid: A pointer to the temperature grid over which the
   * data is defined.
   * @param f_grid: A pointer to the frequency grid over which the
   * data is defined.
   */
  AbsorptionVectorData(std::shared_ptr<const Vector> t_grid,
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

  AbsorptionVectorData& operator=(const matpack::data_t<Scalar, 4>& data) {
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
    return *this;
  }

  constexpr matpack::view_t<CoeffVector, 3> get_coeff_vector_view() {
    return matpack::view_t<CoeffVector, 3>{matpack::mdview_t<CoeffVector, 3>(
        reinterpret_cast<CoeffVector*>(this->data_handle()),
        std::array<Index, 3>{
            this->extent(0), this->extent(1), this->extent(2)})};
  }

  constexpr matpack::view_t<const CoeffVector, 3> get_const_coeff_vector_view() const {
    return matpack::view_t<const CoeffVector, 3>{matpack::mdview_t<const CoeffVector, 3>(
        reinterpret_cast<const CoeffVector*>(this->data_handle()),
        std::array<Index, 3>{
            this->extent(0), this->extent(1), this->extent(2)})};
  }

  AbsorptionVectorData<Scalar, Format::ARO, Representation::Spectral>
  to_spectral() const {
    AbsorptionVectorData<Scalar, Format::ARO, Representation::Spectral> avd_new(
        t_grid_, f_grid_, za_inc_grid_);
    reinterpret_cast<matpack::data_t<Scalar, 4>&>(avd_new) = *this;
    return avd_new;
  }

  AbsorptionVectorData regrid(const ScatteringDataGrids& grids,
                              const RegridWeights& weights) const {
    AbsorptionVectorData result(
        grids.t_grid,
        grids.f_grid,
        std::make_shared<Vector>(grid_vector(*grids.za_inc_grid)));
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
        for (Index i_za_inc = 0;
             i_za_inc < static_cast<Index>(weights.za_inc_grid_weights.size());
             ++i_za_inc) {
          GridPos gp_za_inc  = weights.za_inc_grid_weights[i_za_inc];
          Numeric w_za_inc_l = gp_za_inc.fd[1];
          Numeric w_za_inc_r = gp_za_inc.fd[0];
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
                   coeffs_this[gp_t.idx + 1, gp_f.idx + 1, gp_za_inc.idx + 1]);
        }
      }
    }
    return result;
  }

  AbsorptionVectorData regrid(const ScatteringDataGrids& grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, za_inc_grid_, nullptr, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

  AbsorptionVectorData<Numeric, Format::ARO, Representation::Gridded> to_gridded() const {
    auto abs_new = AbsorptionVectorData<Numeric, Format::ARO, Representation::Gridded>(t_grid_, f_grid_, za_inc_grid_);
    abs_new = *this;
    return abs_new;
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

}  // namespace scattering

#endif  // ABSORPTION_VECTOR_H_
