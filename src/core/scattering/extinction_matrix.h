#ifndef SCATTERING_EXTINCTION_MATRIX_H_
#define SCATTERING_EXTINCTION_MATRIX_H_
#include <memory>

#include "matpack/matpack_data.h"
#include "matpack/matpack_eigen.h"
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
 * @param stokes_dim The number of stokes elements that are stored.
 * @return The number of elements that are required to be stored.
 */
constexpr Index get_n_mat_elems(Format format, Index stokes_dim) {
  // Compact format used for extinction matrix data in TRO format.
  if (format == Format::TRO) {
    return 1;
  }
  // All matrix elements stored for data in other formats.
  switch (stokes_dim) {
    case 1:
      return 1;
    case 2:
      return 2;
    case 3:
      return 3;
    case 4:
      return 3;
  }
  return 0;
}

}  // namespace extinction

template <std::floating_point Scalar,
          Format format,
          Representation representation,
          Index stokes_dim>
class ExtinctionMatrixData;

template <std::floating_point Scalar, Representation repr, Index stokes_dim>
class ExtinctionMatrixData<Scalar, Format::TRO, repr, stokes_dim>
    : public matpack::matpack_data<Scalar, 3> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::matpack_data<Scalar, 3>::resize;
  using matpack::matpack_data<Scalar, 3>::reshape;

 public:
  /// Spectral transform of this extinction matrix.
  using ExtinctionMatrixDataSpectral =
      ExtinctionMatrixData<Scalar,
                           Format::TRO,
                           Representation::Spectral,
                           stokes_dim>;
  using ExtinctionMatrixDataLabFrame =
      ExtinctionMatrixData<Scalar,
                           Format::ARO,
                           Representation::Gridded,
                           stokes_dim>;

  constexpr static Index n_stokes_coeffs =
      extinction::get_n_mat_elems(Format::TRO, stokes_dim);
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
                       std::shared_ptr<const Vector> f_grid)
      : matpack::matpack_data<Scalar, 3>(
            t_grid->size(), f_grid->size(), n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid) {
    matpack::matpack_data<Scalar, 3>::operator=(0.0);
  }

  ExtinctionMatrixData& operator=(
      const matpack::matpack_data<Scalar, 3>& data) {
    ARTS_USER_ERROR_IF(
        data.shape()[0] != n_temps_,
        "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[1] != n_freqs_,
        "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[2] != n_stokes_coeffs,
        "Provided backscatter coefficient data do not match expected number of stokes coefficients.");
    this->template matpack_data<Scalar, 3>::operator=(data);
    return *this;
  }

  constexpr matpack::matpack_view<CoeffVector, 2, false, false>
  get_coeff_vector_view() {
    return matpack::matpack_view<CoeffVector, 2, false, false>(
        reinterpret_cast<CoeffVector*>(this->data_handle()),
        {this->extent(0), this->extent(1)});
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @tparam new_stokes_dim The stokes dimensions for which to extract the
   * relevant data.
   * @return A new extinction matrix data object containing only data required
   * for the requested stokes dimensions.
   */
  template <Index new_stokes_dim>
  PhaseMatrixData<Scalar, Format::TRO, repr, new_stokes_dim>
  extract_stokes_coeffs() const {
    ARTS_ASSERT(new_stokes_dim <= stokes_dim);
    ExtinctionMatrixData<Scalar, Format::TRO, repr, new_stokes_dim> result(
        t_grid_, f_grid_);
    result = this->operator()(joker, joker, result.n_stokes_coeffs);
    return result;
  }

  ExtinctionMatrixData regrid(const ScatteringDataGrids& grids,
                              const RegridWeights& weights) {
    ExtinctionMatrixData result(grids.t_grid, grids.f_grid);
    auto coeffs_this = get_coeff_vector_view();
    auto coeffs_res = result.get_coeff_vector_view();
    for (Index i_t = 0; i_t < static_cast<Index>(weights.t_grid_weights.size()); ++i_t) {
      GridPos gp_t = weights.t_grid_weights[i_t];
      Numeric w_t_l = gp_t.fd[1];
      Numeric w_t_r = gp_t.fd[0];
      for (Index i_f = 0; i_f < static_cast<Index>(weights.f_grid_weights.size()); ++i_f) {
        GridPos gp_f = weights.f_grid_weights[i_f];
        Numeric w_f_l = gp_f.fd[1];
        Numeric w_f_r = gp_f.fd[0];
        coeffs_res(i_t, i_f) =
            (w_t_l * w_f_l * coeffs_this(gp_t.idx, gp_f.idx) +
             w_t_l * w_f_r * coeffs_this(gp_t.idx, gp_f.idx + 1) +
             w_t_r * w_f_l * coeffs_this(gp_t.idx + 1, gp_f.idx) +
             w_t_r * w_f_r * coeffs_this(gp_t.idx + 1, gp_f.idx + 1));
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

template <std::floating_point Scalar, Representation repr, Index stokes_dim>
class ExtinctionMatrixData<Scalar, Format::ARO, repr, stokes_dim>
    : public matpack::matpack_data<Scalar, 4> {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using matpack::matpack_data<Scalar, 4>::resize;
  using matpack::matpack_data<Scalar, 4>::reshape;

 public:
  /// Spectral transform of this extinction matrix.
  using ExtinctionMatrixDataSpectral =
      ExtinctionMatrixData<Scalar,
                           Format::ARO,
                           Representation::Spectral,
                           stokes_dim>;
  using ExtinctionMatrixDataLabFrame =
      ExtinctionMatrixData<Scalar,
                           Format::ARO,
                           Representation::Gridded,
                           stokes_dim>;

  constexpr static Index n_stokes_coeffs =
      extinction::get_n_mat_elems(Format::ARO, stokes_dim);
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
      : matpack::matpack_data<Scalar, 4>(t_grid->size(),
                                         f_grid->size(),
                                         za_inc_grid->size(),
                                         n_stokes_coeffs),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid) {
    matpack::matpack_data<Scalar, 4>::operator=(0.0);
  }

  ExtinctionMatrixData& operator=(
      const matpack::matpack_data<Scalar, 4>& data) {
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
    this->template matpack_data<Scalar, 4>::operator=(data);
    return *this;
  }

  constexpr matpack::matpack_view<CoeffVector, 3, false, false>
  get_coeff_vector_view() {
    return matpack::matpack_view<CoeffVector, 3, false, false>(
        reinterpret_cast<CoeffVector*>(this->data_handle()),
        {this->extent(0), this->extent(1), this->extent(2)});
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @tparam new_stokes_dim The stokes dimensions for which to extract the
   * relevant data.
   * @return A new extinction matrix data object containing only data required
   * for the requested stokes dimensions.
   */
  template <Index new_stokes_dim>
  PhaseMatrixData<Scalar, Format::TRO, Representation::Spectral, new_stokes_dim>
  extract_stokes_coeffs() const {
    ARTS_ASSERT(new_stokes_dim <= stokes_dim);
    ExtinctionMatrixData<Scalar, Format::ARO, repr, new_stokes_dim> result(
        t_grid_, f_grid_, za_inc_grid_);
    result = this->operator()(joker, joker, joker, result.n_stokes_coeffs);
    return result;
  }

  ExtinctionMatrixData regrid(const ScatteringDataGrids& grids,
                              const RegridWeights& weights) {
    ExtinctionMatrixData result(grids.t_grid, grids.f_grid, grids.za_inc_grid);
    auto coeffs_this = get_coeff_vector_view();
    auto coeffs_res = result.get_coeff_vector_view();
    for (Index i_t = 0; i_t < static_cast<Index>(weights.t_grid_weights.size()); ++i_t) {
      GridPos gp_t = weights.t_grid_weights[i_t];
      Numeric w_t_l = gp_t.fd[1];
      Numeric w_t_r = gp_t.fd[0];
      for (Index i_f = 0; i_f < static_cast<Index>(weights.f_grid_weights.size()); ++i_f) {
        GridPos gp_f = weights.f_grid_weights[i_f];
        Numeric w_f_l = gp_f.fd[1];
        Numeric w_f_r = gp_f.fd[0];
        for (Index i_za_inc = 0; i_za_inc < static_cast<Index>(weights.za_inc_grid_weights.size());
             ++i_za_inc) {
          GridPos gp_za_inc = weights.za_inc_grid_weights[i_za_inc];
          Numeric w_za_inc_l = gp_za_inc.fd[1];
          Numeric w_za_inc_r = gp_za_inc.fd[0];
          coeffs_res(i_t, i_f, i_za_inc) =
              (w_t_l * w_f_l * w_za_inc_l *
                   coeffs_this(gp_t.idx, gp_f.idx, gp_za_inc.idx) +
               w_t_l * w_f_l * w_za_inc_r *
                   coeffs_this(gp_t.idx, gp_f.idx, gp_za_inc.idx + 1) +
               w_t_l * w_f_r * w_za_inc_l *
                   coeffs_this(gp_t.idx, gp_f.idx + 1, gp_za_inc.idx) +
               w_t_l * w_f_r * w_za_inc_r *
                   coeffs_this(gp_t.idx, gp_f.idx + 1, gp_za_inc.idx + 1) +

               w_t_r * w_f_l * w_za_inc_l *
                   coeffs_this(gp_t.idx + 1, gp_f.idx, gp_za_inc.idx) +
               w_t_r * w_f_l * w_za_inc_r *
                   coeffs_this(gp_t.idx + 1, gp_f.idx, gp_za_inc.idx + 1) +
               w_t_r * w_f_r * w_za_inc_l *
                   coeffs_this(gp_t.idx + 1, gp_f.idx + 1, gp_za_inc.idx) +
               w_t_r * w_f_r * w_za_inc_r *
                   coeffs_this(gp_t.idx + 1, gp_f.idx + 1, gp_za_inc.idx + 1));
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

}  // namespace scattering

#endif  // SCATTERING_EXTINCTION_MATRIX_H_
