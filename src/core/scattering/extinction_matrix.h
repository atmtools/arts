#ifndef SCATTERING_EXTINCTION_MATRIX_H_
#define SCATTERING_EXTINCTION_MATRIX_H_
#include <matpack.h>
#include <rtepack.h>

#include <memory>

#include "phase_matrix.h"
#include "sht.h"

namespace scattering {

template <Format format, Representation representation> class ExtinctionMatrixData;

/** The extinction matrices [t, f] of totally randomly oriented particles
 *
 * Their extinction matrix is K11 times the identity: A of each
 * compact_planar_propmat is K11, and B = W = 0.
 */
template <Representation repr> class ExtinctionMatrixData<Format::TRO, repr>
    : public rtepack::compact_planar_propmat_matrix {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::compact_planar_propmat_matrix::reshape;
  using rtepack::compact_planar_propmat_matrix::resize;

 public:
  /// Spectral transform of this extinction matrix.
  using ExtinctionMatrixDataSpectral = ExtinctionMatrixData<Format::TRO, Representation::Spectral>;
  using ExtinctionMatrixDataLabFrame = ExtinctionMatrixData<Format::ARO, Representation::Gridded>;
  /// The tensor type used to store the extinction matrix data.
  using TensorType = rtepack::compact_planar_propmat_matrix;

  using TensorType::operator[];

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
  ExtinctionMatrixData(std::shared_ptr<const Vector> t_grid, std::shared_ptr<const Vector> f_grid)
      : TensorType(t_grid->size(), f_grid->size(), rtepack::compact_planar_propmat{}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid) {}

  ExtinctionMatrixData& operator=(const TensorType& data) {
    ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                       "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_,
                       "Provided backscatter coefficient data do not match frequency grid.");
    TensorType::operator=(data);
    return *this;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: TRO extinction holds K11 only.
   */
  [[nodiscard]] ExtinctionMatrixData<Format::TRO, repr> extract_stokes_coeffs() const { return *this; }

  [[nodiscard]] ExtinctionMatrixData<Format::TRO, Representation::Spectral> to_spectral() const {
    ExtinctionMatrixData<Format::TRO, Representation::Spectral> emd_new{t_grid_, f_grid_};
    emd_new = *this;
    return emd_new;
  }

  [[nodiscard]] ExtinctionMatrixData<Format::ARO, repr> to_lab_frame(std::shared_ptr<const Vector> za_inc_grid) const {
    ExtinctionMatrixData<Format::ARO, repr> em_new{t_grid_, f_grid_, za_inc_grid};
    for (Size t_ind = 0; t_ind < t_grid_->size(); ++t_ind) {
      for (Size f_ind = 0; f_ind < f_grid_->size(); ++f_ind) {
        for (Size za_inc_ind = 0; za_inc_ind < za_inc_grid->size(); ++za_inc_ind) {
          em_new[t_ind, f_ind, za_inc_ind] = this->operator[](t_ind, f_ind);
        }
      }
    }
    return em_new;
  }

  [[nodiscard]] ExtinctionMatrixData regrid(const ScatteringDataGrids& grids, const RegridWeights& weights) const {
    ExtinctionMatrixData result(grids.t_grid, grids.f_grid);
    regrid_with(
        static_cast<TensorType&>(result), static_cast<const TensorType&>(*this), weights.t_lags, weights.f_lags);
    return result;
  }

  [[nodiscard]] ExtinctionMatrixData regrid(const ScatteringDataGrids& grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

  [[nodiscard]] ExtinctionMatrixData<Format::TRO, Representation::Gridded> to_gridded() const {
    auto new_ext = ExtinctionMatrixData<Format::TRO, Representation::Gridded>(t_grid_, f_grid_);
    new_ext      = *this;
    return new_ext;
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

/** The extinction matrices [t, f, za_inc] of azimuthally randomly oriented particles
 *
 * In the laboratory (meridional) basis A, B and W of each
 * compact_planar_propmat are K11, K12 and K34.
 */
template <Representation repr> class ExtinctionMatrixData<Format::ARO, repr>
    : public rtepack::compact_planar_propmat_tensor3 {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::compact_planar_propmat_tensor3::reshape;
  using rtepack::compact_planar_propmat_tensor3::resize;

 public:
  /// Spectral transform of this extinction matrix.
  using ExtinctionMatrixDataSpectral = ExtinctionMatrixData<Format::ARO, Representation::Spectral>;
  using ExtinctionMatrixDataLabFrame = ExtinctionMatrixData<Format::ARO, Representation::Gridded>;
  /// The tensor type used to store the extinction matrix data.
  using TensorType = rtepack::compact_planar_propmat_tensor3;

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
                       std::shared_ptr<const Vector> f_grid,
                       std::shared_ptr<const Vector> za_inc_grid)
      : TensorType(t_grid->size(), f_grid->size(), za_inc_grid->size(), rtepack::compact_planar_propmat{}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid) {}

  ExtinctionMatrixData& operator=(const TensorType& data) {
    ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                       "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_,
                       "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(data.shape()[2] != n_za_inc_,
                       "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
    TensorType::operator=(data);
    return *this;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: ARO extinction holds K11, K12 and K34 only.
   */
  [[nodiscard]] ExtinctionMatrixData extract_stokes_coeffs() const { return *this; }

  [[nodiscard]] ExtinctionMatrixData<Format::ARO, Representation::Spectral> to_spectral() const;

  //! The same data, as the extinction matrix of azimuthal Fourier-mode optics (it has no azimuth dependence)
  [[nodiscard]] ExtinctionMatrixData<Format::ARO, Representation::Fourier> to_fourier() const {
    ExtinctionMatrixData<Format::ARO, Representation::Fourier> emd_new{t_grid_, f_grid_, za_inc_grid_};
    emd_new = *this;
    return emd_new;
  }

  [[nodiscard]] ExtinctionMatrixData regrid(const ScatteringDataGrids& grids, const RegridWeights& weights) const {
    ExtinctionMatrixData result(grids.t_grid, grids.f_grid, grids.za_inc_grid);
    regrid_with(static_cast<TensorType&>(result),
                static_cast<const TensorType&>(*this),
                weights.t_lags,
                weights.f_lags,
                weights.za_inc_lags);
    return result;
  }

  [[nodiscard]] ExtinctionMatrixData regrid(const ScatteringDataGrids& grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, za_inc_grid_, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

  [[nodiscard]] ExtinctionMatrixData<Format::ARO, Representation::Gridded> to_gridded() const {
    auto new_ext = ExtinctionMatrixData<Format::ARO, Representation::Gridded>(t_grid_, f_grid_, za_inc_grid_);
    new_ext      = *this;
    return new_ext;
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

template <Representation repr> [[nodiscard]] ExtinctionMatrixData<Format::ARO, Representation::Spectral>
ExtinctionMatrixData<Format::ARO, repr>::to_spectral() const {
  ExtinctionMatrixData<Format::ARO, Representation::Spectral> emd_new{t_grid_, f_grid_, za_inc_grid_};
  emd_new = *this;
  return emd_new;
}

}  // namespace scattering

#endif  // SCATTERING_EXTINCTION_MATRIX_H_
