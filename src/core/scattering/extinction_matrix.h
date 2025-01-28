#ifndef SCATTERING_EXTINCTION_MATRIX_H_
#define SCATTERING_EXTINCTION_MATRIX_H_
#include <matpack.h>

#include <memory>

#include "phase_matrix.h"
#include "sht.h"

namespace scattering {

namespace extinction {

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
  // All matrix elements stored for data in other formats.
  return 3;
}

}  // namespace extinction

template <std::floating_point Scalar,
          Format format,
          Representation representation>
class ExtinctionMatrixData;

template <std::floating_point Scalar, Representation repr>
class ExtinctionMatrixData<Scalar, Format::TRO, repr>
    : public matpack::data_t<Scalar, 3> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::data_t<Scalar, 3>::resize;
  using matpack::data_t<Scalar, 3>::reshape;

 public:
  /// Spectral transform of this extinction matrix.
  using ExtinctionMatrixDataSpectral =
      ExtinctionMatrixData<Scalar, Format::TRO, Representation::Spectral>;
  using ExtinctionMatrixDataLabFrame =
      ExtinctionMatrixData<Scalar, Format::ARO, Representation::Gridded>;

  constexpr static Index n_stokes_coeffs =
      extinction::get_n_mat_elems(Format::TRO);
  using CoeffVector = Eigen::Matrix<Scalar, 1, n_stokes_coeffs>;

  using matpack::data_t<Scalar, 3>::operator[];

  ExtinctionMatrixData() = default;
  /** Create a new ExtinctionMatrixData container.
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
  ExtinctionMatrixData(std::shared_ptr<const Vector> t_grid,
                       std::shared_ptr<const Vector> f_grid)
      : matpack::data_t<Scalar, 3>(
            t_grid->size(), f_grid->size(), n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid) {
    matpack::data_t<Scalar, 3>::operator=(0.0);
  }

  ExtinctionMatrixData& operator=(const matpack::data_t<Scalar, 3>& data) {
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

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A new extinction matrix data object containing only data required
   * for the requested stokes dimensions.
   */
  ExtinctionMatrixData<Scalar, Format::TRO, repr> extract_stokes_coeffs()
      const {
    ExtinctionMatrixData<Scalar, Format::TRO, repr> result(t_grid_, f_grid_);
    result = this->operator()(joker, joker, result.n_stokes_coeffs);
    return result;
  }

  ExtinctionMatrixData<Scalar, Format::TRO, Representation::Spectral>
  to_spectral() {
    ExtinctionMatrixData<Scalar, Format::TRO, Representation::Spectral> emd_new{
        t_grid_, f_grid_};
    reinterpret_cast<matpack::data_t<Scalar, 3>&>(emd_new) = *this;
    return emd_new;
  }

  ExtinctionMatrixData<Scalar, Format::ARO, repr> to_lab_frame(
      std::shared_ptr<const Vector> za_inc_grid) {
    ExtinctionMatrixData<Scalar, Format::ARO, repr> em_new{
        t_grid_, f_grid_, za_inc_grid};
    for (Size t_ind = 0; t_ind < t_grid_->size(); ++t_ind) {
      for (Size f_ind = 0; f_ind < f_grid_->size(); ++f_ind) {
        for (Size za_inc_ind = 0; za_inc_ind < za_inc_grid->size();
             ++za_inc_ind) {
          em_new[t_ind, f_ind, za_inc_ind, 0] =
              this->operator[](t_ind, f_ind, 0);
        }
      }
    }
    return em_new;
  }

  ExtinctionMatrixData regrid(const ScatteringDataGrids& grids,
                              const RegridWeights& weights) {
    ExtinctionMatrixData result(grids.t_grid, grids.f_grid);
    auto coeffs_this = get_coeff_vector_view();
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
class ExtinctionMatrixData<Scalar, Format::ARO, repr>
    : public matpack::data_t<Scalar, 4> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::data_t<Scalar, 4>::resize;
  using matpack::data_t<Scalar, 4>::reshape;

 public:
  /// Spectral transform of this extinction matrix.
  using ExtinctionMatrixDataSpectral =
      ExtinctionMatrixData<Scalar, Format::ARO, Representation::Spectral>;
  using ExtinctionMatrixDataLabFrame =
      ExtinctionMatrixData<Scalar, Format::ARO, Representation::Gridded>;

  constexpr static Index n_stokes_coeffs =
      extinction::get_n_mat_elems(Format::ARO);
  using CoeffVector = Eigen::Matrix<Scalar, 1, n_stokes_coeffs>;

  ExtinctionMatrixData() {};
  /** Create a new ExtinctionMatrixData container.
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
  ExtinctionMatrixData(std::shared_ptr<const Vector> t_grid,
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

  ExtinctionMatrixData& operator=(const matpack::data_t<Scalar, 4>& data) {
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

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A new extinction matrix data object containing only data required
   * for the requested stokes dimensions.
   */
  PhaseMatrixData<Scalar, Format::TRO, Representation::Spectral>
  extract_stokes_coeffs() const {
    ExtinctionMatrixData<Scalar, Format::ARO, repr> result(
        t_grid_, f_grid_, za_inc_grid_);
    result = this->operator[](joker, joker, joker, result.n_stokes_coeffs);
    return result;
  }

  ExtinctionMatrixData<Scalar, Format::ARO, Representation::Spectral>
  to_spectral();

  ExtinctionMatrixData regrid(const ScatteringDataGrids& grids,
                              const RegridWeights& weights) {
    ExtinctionMatrixData result(grids.t_grid, grids.f_grid, grids.za_inc_grid);
    auto coeffs_this = get_coeff_vector_view();
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

template <std::floating_point Scalar, Representation repr>
ExtinctionMatrixData<Scalar, Format::ARO, Representation::Spectral>
ExtinctionMatrixData<Scalar, Format::ARO, repr>::to_spectral() {
  ExtinctionMatrixData<Scalar, Format::ARO, Representation::Spectral> emd_new{
      t_grid_, f_grid_, za_inc_grid_};
  emd_new = reinterpret_cast<matpack::data_t<Scalar, 4>&>(*this);
  return emd_new;
}

}  // namespace scattering

#endif  // SCATTERING_EXTINCTION_MATRIX_H_
