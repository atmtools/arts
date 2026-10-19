#pragma once

#include <arts_omp.h>
#include <legendre.h>
#include <matpack.h>
#include <rtepack.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>

#include "arts_constants.h"
#include "integration.h"
#include "sht.h"
#include "tro_legendre.h"
#include "utils.h"

namespace scattering {

using sht::SHT;

enum class Format : char { TRO, ARO, General };
std::ostream &operator<<(std::ostream &out, Format format);

/** The representation of scattering data
 *
 * Gridded: on angular grids.  Spectral: spherical-harmonics (SHT) coefficients
 * over the scattering directions (Legendre coefficients for TRO).  Fourier:
 * azimuthal Fourier modes of the laboratory-frame phase matrix at given zenith
 * angles, the form plane-parallel solvers consume (ARO only).
 */
enum class Representation : char { Gridded, Spectral, DoublySpectral, Fourier };
std::ostream &operator<<(std::ostream &out, Representation repr);

namespace detail {
/** Calculate rotation coefficients for scattering matrix.
 *
 * This method calculates the rotation coefficients that are required to
 * transform the scattering matrix of a randomly-oriented particle to the
 * phase matrix, which describes its scattering properties w.r.t. to the
 * laboratory frame. This equation calculates the angle Theta and the
 * coefficients C_1, * C_2, S_1, S_2 as defined in equation (4.16) in
 * "Scattering, Absorption, and Emission of Light by Small Particles."
 *
 * The coefficients are those of rtepack::scattering_plane_rotation, so the
 * phase matrix is rtepack::rotated(F, C_1, S_1, C_2, S_2).  With x = East,
 * y = North and z = Up, the direction (za, aa) is
 * k = (sin za sin aa, sin za cos aa, cos za), in the meridional basis
 * rtepack::meridional_basis(k, h) of ARTS's propagation matrices; at a pole
 * h = (cos aa, -sin aa, 0), the limit along the ray's meridian.  Only the
 * side |aa_scat - aa_inc| <= 180 deg is computed: the phase matrix on the
 * other side is the mirror image (rtepack::mirror) of that at 360 deg
 * minus the azimuth difference, which the callers apply.
 *
 * @param aa_inc The azimuth-angle component of the incoming angle in degree.
 * @param za_inc The zenith-angle component of the incoming angle in degree.
 * @param aa_scat The azimuth-angle component of the scattering angle in degree.
 * @param za_scat The zenith-angle component of the scattering angle in degree.
 */
std::array<Numeric, 5> rotation_coefficients(Numeric aa_inc_d, Numeric za_inc_d, Numeric aa_scat_d, Numeric za_scat_d);
}  // namespace detail

/// The grid over which the scattering data is defined.
struct ScatteringDataGrids {
  ScatteringDataGrids(std::shared_ptr<const Vector> t_grid_, std::shared_ptr<const Vector> f_grid_);

  ScatteringDataGrids(std::shared_ptr<const Vector>          t_grid_,
                      std::shared_ptr<const Vector>          f_grid_,
                      std::shared_ptr<const ZenithAngleGrid> za_scat_grid_);

  ScatteringDataGrids(std::shared_ptr<const Vector>          t_grid_,
                      std::shared_ptr<const Vector>          f_grid_,
                      std::shared_ptr<const Vector>          za_inc_grid_,
                      std::shared_ptr<const Vector>          delta_aa_grid_,
                      std::shared_ptr<const ZenithAngleGrid> za_scat_grid_);

  std::shared_ptr<const Vector>          t_grid;
  std::shared_ptr<const Vector>          f_grid;
  std::shared_ptr<const Vector>          aa_inc_grid;
  std::shared_ptr<const Vector>          za_inc_grid;
  std::shared_ptr<const Vector>          aa_scat_grid;
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

RegridWeights calc_regrid_weights(std::shared_ptr<const Vector>          t_grid,
                                  std::shared_ptr<const Vector>          f_grid,
                                  std::shared_ptr<const Vector>          aa_inc_grid,
                                  std::shared_ptr<const Vector>          za_inc_grid,
                                  std::shared_ptr<const Vector>          aa_scat_grid,
                                  std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                  ScatteringDataGrids                    new_grids);

template <Format format> class BackscatterMatrixData;

/** The scattering matrix at one scattering angle of a totally randomly oriented scatterer [t, f]
 *
 * The forward- and backscattering matrices, in the scattering-plane basis.
 */
template <> class BackscatterMatrixData<Format::TRO> : public rtepack::compact_planar_muelmat_matrix {
 public:
  using TensorType = rtepack::compact_planar_muelmat_matrix;

  BackscatterMatrixData(std::shared_ptr<const Vector> t_grid, std::shared_ptr<const Vector> f_grid)
      : TensorType(t_grid->size(), f_grid->size(), rtepack::compact_planar_muelmat{}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid) {}

  BackscatterMatrixData &operator=(const TensorType &data);

  [[nodiscard]] std::shared_ptr<const Vector> get_t_grid_ptr() const { return t_grid_; }

  [[nodiscard]] std::shared_ptr<const Vector> get_f_grid_ptr() const { return f_grid_; }

  [[nodiscard]] BackscatterMatrixData extract_stokes_coeffs() const { return *this; }

  [[nodiscard]] BackscatterMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const;

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids) const;

  BackscatterMatrixData to_gridded() const { return *this; }

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

/** The laboratory-frame phase matrix at one scattering direction of an azimuthally randomly oriented scatterer
 *
 * The forward- and backscattering Mueller matrices [t, f, za_inc].
 */
template <> class BackscatterMatrixData<Format::ARO> : public rtepack::muelmat_tensor3 {
 public:
  using TensorType = rtepack::muelmat_tensor3;

  BackscatterMatrixData(std::shared_ptr<const Vector> t_grid,
                        std::shared_ptr<const Vector> f_grid,
                        std::shared_ptr<const Vector> za_inc_grid)
      : TensorType(t_grid->size(), f_grid->size(), za_inc_grid->size(), rtepack::muelmat{0.0}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid) {}

  BackscatterMatrixData(const BackscatterMatrixData<Format::TRO> &bsmat, std::shared_ptr<const Vector> za_inc_grid);

  BackscatterMatrixData &operator=(const TensorType &data);

  BackscatterMatrixData extract_stokes_coeffs() const { return *this; }

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const;

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids) const;

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

template <Format format> using ForwardscatterMatrixData = BackscatterMatrixData<format>;

///////////////////////////////////////////////////////////////////////////////
// Phase matrix data
///////////////////////////////////////////////////////////////////////////////

template <Format format, Representation representation> class PhaseMatrixData;

/** The scattering matrix of a totally randomly oriented scatterer at a scattering angle
 *
 * scattering_matrix(theta, F) sets F [t, f] to the scattering-plane Mueller matrix at scattering angle theta
 * [rad]; see tro_lab_frame.  It may keep scratch between calls.
 */
using ScatteringMatrixFunction = std::function<void(Numeric, rtepack::compact_planar_muelmat_matrix &)>;

/** Phase matrix in the laboratory frame of a totally randomly oriented scatterer.
 *
 * Rotates the scattering matrix into the laboratory frame for every incidence
 * zenith angle, azimuth difference and scattering zenith angle (Mishchenko et
 * al., 2002, Eq. 4.16).  This is the conversion between the two frames for
 * every representation of such a scatterer; the representation only supplies
 * its scattering matrix at the scattering angle of each direction pair:
 *
 * scattering_matrix(theta, F) sets F[i_t, i_f, joker] to the six independent
 * elements [F11, F12, F22, F33, F34, F44] at scattering angle theta [rad] for
 * every temperature of t_grid and frequency of f_grid.
 *
 * The result is as accurate as scattering_matrix: exact for a closed form,
 * interpolated for gridded data.
 *
 * @param t_grid The temperature grid
 * @param f_grid The frequency grid
 * @param za_inc_grid The incidence zenith angles [deg]
 * @param delta_aa_grid The azimuth differences [deg]
 * @param za_scat_grid The scattering zenith angles [deg]
 * @param scattering_matrix The scattering matrix as above
 */
PhaseMatrixData<Format::ARO, Representation::Gridded> tro_lab_frame(std::shared_ptr<const Vector> t_grid,
                                                                    std::shared_ptr<const Vector> f_grid,
                                                                    std::shared_ptr<const Vector> za_inc_grid,
                                                                    std::shared_ptr<const Vector> delta_aa_grid,
                                                                    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                                                    const ScatteringMatrixFunction &scattering_matrix);

/** The azimuthal Fourier modes of the laboratory-frame phase matrix of a totally randomly oriented scatterer.
 *
 * At fixed incidence and scattering zenith angles, the laboratory-frame
 * phase matrix of tro_lab_frame is an analytic periodic function of the
 * azimuth difference, except where the scattering plane is undefined: at
 * exact forward (delta = 0, za_inc = za_scat) and backward (delta = 180 deg)
 * scattering it has a kink unless the scattering matrix is regular there, as
 * a physical one is (F12 = F34 = 0, F22 = F33 forward and F22 = -F33
 * backward), and a Legendre series of data is only to its accuracy.  Z is
 * mirror symmetric, Z(360 deg - delta) = D Z(delta) D with D = diag(1, 1,
 * -1, -1), so the elements with exactly one index in {U, V} have sine modes
 * only and the others cosine modes only, from [0, 180] deg alone.  The modes
 * are integrated with a Gauss-Legendre rule on [0, 180] deg, whose ends those
 * directions are, doubling the nodes until the modes no longer change to
 * rounding (1e-13 of the largest phase-matrix element of the direction pair);
 * a pair that has not converged at 4096 azimuths is an error.
 *
 * scattering_matrix is as for tro_lab_frame, and must be exact: the modes
 * are as accurate as it is.  phase_integral [t, f] is its 2 pi int F11
 * dcos(Theta), the phase integral at every incidence angle; the producer
 * knows it exactly (sqrt(4 pi) a_0 of a Legendre series, the scattering
 * coefficient of a normalised closed form).
 *
 * @param t_grid The temperature grid
 * @param f_grid The frequency grid
 * @param za_inc_grid The incidence zenith angles [deg]
 * @param za_scat_grid The scattering zenith angles [deg]
 * @param max_mode The highest mode M >= 0
 * @param phase_integral [t, f], 2 pi int F11 dcos(Theta) of scattering_matrix
 * @param scattering_matrix As for tro_lab_frame
 */
PhaseMatrixData<Format::ARO, Representation::Fourier> tro_lab_frame_fourier_modes(
    std::shared_ptr<const Vector>          t_grid,
    std::shared_ptr<const Vector>          f_grid,
    std::shared_ptr<const Vector>          za_inc_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
    Index                                  max_mode,
    const ConstMatrixView                 &phase_integral,
    const ScatteringMatrixFunction        &scattering_matrix);

/** The phase matrix of a totally randomly oriented scatterer on a scattering-angle grid
 *
 * The scattering-plane Mueller matrices [t, f, za_scat], in compact form.
 */
template <> class PhaseMatrixData<Format::TRO, Representation::Gridded>
    : public rtepack::compact_planar_muelmat_tensor3 {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::compact_planar_muelmat_tensor3::reshape;
  using rtepack::compact_planar_muelmat_tensor3::resize;

 public:
  /// Spectral transform of this phase matrix.
  using PhaseMatrixDataSpectral = PhaseMatrixData<Format::TRO, Representation::Spectral>;
  using PhaseMatrixDataLabFrame = PhaseMatrixData<Format::ARO, Representation::Gridded>;

  /// The tensor type used to store the phase matrix data.
  using TensorType = rtepack::compact_planar_muelmat_tensor3;

  using TensorType::operator[];

  PhaseMatrixData()                                   = default;
  PhaseMatrixData(const PhaseMatrixData &)            = default;
  PhaseMatrixData(PhaseMatrixData &&)                 = default;
  PhaseMatrixData &operator=(const PhaseMatrixData &) = default;
  PhaseMatrixData &operator=(PhaseMatrixData &&)      = default;

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
  PhaseMatrixData(std::shared_ptr<const Vector>          t_grid,
                  std::shared_ptr<const Vector>          f_grid,
                  std::shared_ptr<const ZenithAngleGrid> za_scat_grid)
      : TensorType(t_grid->size(), f_grid->size(), grid_size(*za_scat_grid), rtepack::compact_planar_muelmat{}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_scat_(grid_size(*za_scat_grid)),
        za_scat_grid_(za_scat_grid) {}

  template <Format format, Representation repr> PhaseMatrixData(const PhaseMatrixData<format, repr> &other) {
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

  PhaseMatrixData &operator=(const TensorType &data);

  std::shared_ptr<const Vector>          get_t_grid() const { return t_grid_; }
  std::shared_ptr<const Vector>          get_f_grid() const { return f_grid_; }
  std::shared_ptr<const ZenithAngleGrid> get_za_scat_grid() const { return za_scat_grid_; }

  /** The Legendre series of this phase matrix to the given degree.
   *
   * The coefficients are those of the function the gridded data define (see
   * tro_legendre.h): linear in the scattering angle between the nodes and
   * constant beyond the first and the last.  The projection is exact for that
   * function; whether the degree resolves it is what legendre_report tells.
   *
   * @param degree The highest Legendre degree, >= 0
   * @param order Must be 0: a TRO phase matrix depends on the scattering angle only
   */
  PhaseMatrixDataSpectral to_spectral(Index degree, Index order = 0) const;

  /** How well a Legendre series represents this phase matrix, for every temperature and frequency.
   *
   * normalisation_error is left NaN; SingleScatteringData fills it.
   *
   * @param spectral A series on the same temperature and frequency grids, e.g. to_spectral(degree)
   */
  LegendreReport legendre_report(const PhaseMatrixDataSpectral &spectral) const;

  /** The phase matrix in the laboratory frame.
   *
   * The scattering matrix is interpolated linearly in scattering angle, and
   * held constant outside the range of the scattering-angle grid.
   */
  PhaseMatrixDataLabFrame to_lab_frame(std::shared_ptr<const Vector>          za_inc_grid,
                                       std::shared_ptr<const Vector>          delta_aa_grid,
                                       std::shared_ptr<const ZenithAngleGrid> za_scat_grid_new) const;

  BackscatterMatrixData<Format::TRO> extract_backscatter_matrix();

  ForwardscatterMatrixData<Format::TRO> extract_forwardscatter_matrix();

  /** Calculate scattering-angle integral.
   *
   * Integrates the phase matrix over the scattering angles.
   * @return The integral of the phase matrix [t, f]
   */
  rtepack::compact_planar_muelmat_matrix integrate_phase_matrix();

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: TRO data hold all six elements.
   */
  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const;

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const;

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

/** The Legendre series of the phase matrix of a totally randomly oriented scatterer.
 *
 * Holds, for every temperature and frequency, the coefficients a_l,
 * l = 0..degree, on the orthonormal Y_l0 of the scattering-plane phase matrix
 * (see tro_legendre.h): a compact_planar_muelmat per degree, [t, f, l].  A TRO
 * phase matrix depends on the scattering angle only, so its coefficients are
 * real.
 */
template <> class PhaseMatrixData<Format::TRO, Representation::Spectral>
    : public rtepack::compact_planar_muelmat_tensor3 {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::compact_planar_muelmat_tensor3::reshape;
  using rtepack::compact_planar_muelmat_tensor3::resize;

 public:
  /// Gridded transform of this phase matrix.
  using PhaseMatrixDataGridded  = PhaseMatrixData<Format::TRO, Representation::Gridded>;
  using PhaseMatrixDataSpectral = PhaseMatrixData<Format::TRO, Representation::Spectral>;
  using PhaseMatrixDataLabFrame = PhaseMatrixData<Format::ARO, Representation::Gridded>;
  using TensorType              = rtepack::compact_planar_muelmat_tensor3;

  PhaseMatrixData()                                   = default;
  PhaseMatrixData(const PhaseMatrixData &)            = default;
  PhaseMatrixData(PhaseMatrixData &&)                 = default;
  PhaseMatrixData &operator=(const PhaseMatrixData &) = default;
  PhaseMatrixData &operator=(PhaseMatrixData &&)      = default;

  /** A zero series of the given degree.
   *
   * @param t_grid The temperature grid
   * @param f_grid The frequency grid
   * @param degree The highest Legendre degree, >= 0
   */
  PhaseMatrixData(std::shared_ptr<const Vector> t_grid, std::shared_ptr<const Vector> f_grid, Index degree)
      : TensorType(t_grid->size(), f_grid->size(), checked_degree(degree) + 1, rtepack::compact_planar_muelmat{}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        degree_(degree) {}

  PhaseMatrixData &operator=(const TensorType &data);

  std::shared_ptr<const Vector> get_t_grid() const { return t_grid_; }
  std::shared_ptr<const Vector> get_f_grid() const { return f_grid_; }
  Index                         get_degree() const { return degree_; }

  /** The coefficients [degree + 1] at one temperature and frequency */
  CompactPlanarMuelmatVector coefficients(Index i_t, Index i_f) const {
    return CompactPlanarMuelmatVector{(*this)[i_t, i_f, joker]};
  }

  /** The series at the scattering angles of a grid, exactly
   *
   * @param za_scat_grid The scattering angles [deg]
   */
  PhaseMatrixDataGridded to_gridded(std::shared_ptr<const ZenithAngleGrid> za_scat_grid) const;

  /** The series at the 2 degree + 2 Gauss-Legendre nodes, exactly */
  PhaseMatrixDataGridded to_gridded() const {
    return to_gridded(std::make_shared<const ZenithAngleGrid>(GaussLegendreGrid(2 * degree_ + 2)));
  }

  /** The laboratory-frame phase matrix, the series evaluated at the exact scattering angle of every direction pair */
  PhaseMatrixDataLabFrame to_lab_frame(std::shared_ptr<const Vector>          za_inc_grid,
                                       std::shared_ptr<const Vector>          delta_aa_grid,
                                       std::shared_ptr<const ZenithAngleGrid> za_scat_grid_new) const;

  /** The azimuthal Fourier modes m = 0..max_mode of the laboratory-frame phase matrix of the series
   *
   * See tro_lab_frame_fourier_modes; the series is evaluated exactly at every scattering angle.
   */
  PhaseMatrixData<Format::ARO, Representation::Fourier> to_lab_frame_fourier_modes(
      std::shared_ptr<const Vector>          za_inc_grid,
      std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
      Index                                  max_mode) const;

  /** The scattering matrix as tro_lab_frame and tro_lab_frame_fourier_modes take it, from the series */
  ScatteringMatrixFunction scattering_matrix_function() const;

  BackscatterMatrixData<Format::TRO> extract_backscatter_matrix() const { return extract_at(180.0); }

  ForwardscatterMatrixData<Format::TRO> extract_forwardscatter_matrix() const { return extract_at(0.0); }

  /** Calculate scattering-angle integral.
   *
   * Integrates the phase matrix over the scattering angles: 2 pi int F dcos(Theta) = sqrt(4 pi) a_0.
   * @return The integral of the phase matrix [t, f]
   */
  rtepack::compact_planar_muelmat_matrix integrate_phase_matrix() const;

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: TRO data hold all six elements.
   */
  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  /** The series truncated to a lower degree.
   *
   * @param degree The new degree, at most this series' degree: a truncated
   * series does not tell its higher coefficients
   * @param order Must be 0
   */
  PhaseMatrixData to_spectral(Index degree, Index order = 0) const;

  /** Interpolate the coefficients linearly in temperature and frequency.
   *
   * The angular dependence is the series and needs no interpolation; the
   * zenith-angle grids of grids are not used.
   */
  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights weights) const;

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const;

 protected:
  static Index checked_degree(Index degree);

  /** The series at one scattering angle [deg] */
  BackscatterMatrixData<Format::TRO> extract_at(Numeric angle) const;

  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;

  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;

  /// The highest Legendre degree.
  Index degree_;
};

///////////////////////////////////////////////////////////////////////////////
// ARO format
///////////////////////////////////////////////////////////////////////////////

/** The laboratory-frame phase matrix of an azimuthally randomly oriented scatterer on angular grids
 *
 * The Mueller matrices [t, f, za_inc, delta_aa, za_scat].
 */
template <> class PhaseMatrixData<Format::ARO, Representation::Gridded> : public rtepack::muelmat_tensor5 {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::muelmat_tensor5::reshape;
  using rtepack::muelmat_tensor5::resize;

 public:
  /// Spectral transform of this phase matrix.
  using PhaseMatrixDataSpectral = PhaseMatrixData<Format::ARO, Representation::Spectral>;
  /// Azimuthal Fourier modes of this phase matrix.
  using PhaseMatrixDataFourier = PhaseMatrixData<Format::ARO, Representation::Fourier>;
  /// The tensor type used to store the phase matrix data.
  using TensorType = rtepack::muelmat_tensor5;

  PhaseMatrixData()                                   = default;
  PhaseMatrixData(const PhaseMatrixData &)            = default;
  PhaseMatrixData(PhaseMatrixData &&)                 = default;
  PhaseMatrixData &operator=(const PhaseMatrixData &) = default;
  PhaseMatrixData &operator=(PhaseMatrixData &&)      = default;

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
  PhaseMatrixData(std::shared_ptr<const Vector>          t_grid,
                  std::shared_ptr<const Vector>          f_grid,
                  std::shared_ptr<const Vector>          za_inc_grid,
                  std::shared_ptr<const Vector>          delta_aa_grid,
                  std::shared_ptr<const ZenithAngleGrid> za_scat_grid)
      : TensorType(t_grid->size(),
                   f_grid->size(),
                   za_inc_grid->size(),
                   delta_aa_grid->size(),
                   grid_size(*za_scat_grid),
                   rtepack::muelmat{0.0}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid),
        n_delta_aa_(delta_aa_grid->size()),
        delta_aa_grid_(delta_aa_grid),
        n_za_scat_(grid_size(*za_scat_grid)),
        za_scat_grid_(za_scat_grid) {}

  std::shared_ptr<const Vector>          get_t_grid() const { return t_grid_; }
  std::shared_ptr<const Vector>          get_f_grid() const { return f_grid_; }
  std::shared_ptr<const Vector>          get_za_inc_grid() const { return za_inc_grid_; }
  std::shared_ptr<const Vector>          get_delta_aa_grid() const { return delta_aa_grid_; }
  std::shared_ptr<const ZenithAngleGrid> get_za_scat_grid() const { return za_scat_grid_; }

  PhaseMatrixData &operator=(const TensorType &data);

  /** Transform phase matrix to spectral format.
   *
   * @param Pointer to the SHT to use for the transformation.
   */
  PhaseMatrixDataSpectral to_spectral(std::shared_ptr<SHT> sht) const;

  PhaseMatrixDataSpectral to_spectral(Index degree, Index order) const;

  PhaseMatrixDataSpectral to_spectral() const;

  /** The azimuthal Fourier modes of this phase matrix, m = 0..max_mode.
   *
   * The modes are those of the function the gridded data define in the
   * azimuth difference: linear between the nodes of the azimuth-difference
   * grid, which must span one period (its last angle 360 deg above its
   * first).  The integrals are exact for that function.  The zenith-angle
   * grids are kept.  The phase integral is that of the data linear in the
   * scattering zenith angle between the nodes, exactly; NaN unless the nodes
   * span [0, 180] deg.
   *
   * @param max_mode The highest mode, >= 0
   */
  PhaseMatrixDataFourier to_fourier(Index max_mode) const;

  BackscatterMatrixData<Format::ARO> extract_backscatter_matrix();

  ForwardscatterMatrixData<Format::ARO> extract_forwardscatter_matrix();

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: ARO data hold all sixteen elements.
   */
  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const;

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const;

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

/** The spherical-harmonic transform of the laboratory-frame phase matrix of an azimuthally randomly oriented scatterer
 *
 * The coefficients [t, f, za_inc, nlm] of the SHT over the scattering
 * directions, a specmat per coefficient: the complex coefficients of the 16
 * elements of the Mueller matrix.
 */
template <> class PhaseMatrixData<Format::ARO, Representation::Spectral> : public rtepack::specmat_tensor4 {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::specmat_tensor4::reshape;
  using rtepack::specmat_tensor4::resize;

 public:
  /// Spectral transform of this phase matrix.
  using PhaseMatrixDataGridded = PhaseMatrixData<Format::ARO, Representation::Gridded>;
  /// The tensor type used to store the phase matrix data.
  using TensorType = rtepack::specmat_tensor4;

  PhaseMatrixData()                                   = default;
  PhaseMatrixData(const PhaseMatrixData &)            = default;
  PhaseMatrixData(PhaseMatrixData &&)                 = default;
  PhaseMatrixData &operator=(const PhaseMatrixData &) = default;
  PhaseMatrixData &operator=(PhaseMatrixData &&)      = default;

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
                  std::shared_ptr<SHT>          sht)
      : TensorType(
            t_grid->size(), f_grid->size(), za_inc_grid->size(), sht->get_n_spectral_coeffs(), rtepack::specmat{0.0}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid),
        n_spectral_coeffs_(sht->get_n_spectral_coeffs()),
        sht_(sht) {}

  PhaseMatrixData &operator=(const TensorType &data);

  std::shared_ptr<const Vector> get_t_grid() const { return t_grid_; }
  std::shared_ptr<const Vector> get_f_grid() const { return f_grid_; }
  std::shared_ptr<const Vector> get_za_inc_grid() const { return za_inc_grid_; }
  std::shared_ptr<SHT>          get_sht() const { return sht_; }

  /** The azimuthal Fourier modes m = 0..max_mode at the scattering zenith angles of a grid, exactly
   *
   * The series has orders up to the SHT's m_max in the azimuth difference, so
   * equally spaced azimuths, more than m_max + max_mode of them, integrate
   * the modes exactly, and the modes above m_max are zero.  The series is
   * evaluated exactly at each scattering zenith angle.  The phase integral is
   * that of the series (integrate_phase_matrix).
   *
   * @param za_scat_grid The scattering zenith angles [deg]
   * @param max_mode The highest mode, >= 0
   */
  PhaseMatrixData<Format::ARO, Representation::Fourier> to_fourier(std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                                                   Index max_mode) const;

  /** Transform phase matrix to gridded format.
   *
   * @param Pointer to the SHT to use for the transformation.
   */
  PhaseMatrixDataGridded to_gridded() const;

  /** Transform phase matixr to spectral format.
   *
   * @param Pointer to the SHT to use for the transformation.
   */
  PhaseMatrixData to_spectral(Index l_new, Index m_new) const;

  BackscatterMatrixData<Format::ARO> extract_backscatter_matrix() const;

  ForwardscatterMatrixData<Format::ARO> extract_forwardscatter_matrix();

  /** The integral over all scattering directions of the phase matrix [t, f, za_inc]
   *
   * sqrt(4 pi) times the l = 0 coefficient on the orthonormal Y_00; a
   * degree-0 SHT holds the value itself as its coefficient.
   */
  rtepack::muelmat_tensor3 integrate_phase_matrix() const;

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: ARO data hold all sixteen elements.
   */
  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights weights) const;

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const;

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

/** The azimuthal Fourier modes of the laboratory-frame phase matrix of an azimuthally randomly oriented scatterer.
 *
 * For every temperature, frequency, incidence zenith angle and scattering
 * zenith angle, the phase matrix as a function of the azimuth difference
 * Delta = aa_scat - aa_inc is
 *
 *   Z(Delta) = sum_{m = 0}^{M} [C_m cos(m Delta) + S_m sin(m Delta)],
 *
 * so C_0 is the azimuthal mean, and S_0 = 0.  This is the form the
 * plane-parallel solvers (VDISORT, RT4) consume at their streams.  The data
 * are the Mueller matrices [t, f, za_inc, za_scat, m, 2 (C, S)].
 *
 * The zenith directions are those of the grids, not interpolated: the modes
 * are known there and nowhere else.
 *
 * The data also carry the phase integral sigma(za_inc) = int Z11 dOmega over
 * all scattering directions, [t, f, za_inc], which the modes at a few
 * scattering zenith angles cannot give: the scattering coefficient implied by
 * the phase matrix, to be compared with K11 - a1.  Its producers give it
 * exactly, or NaN where their data do not cover all scattering directions.
 */
template <> class PhaseMatrixData<Format::ARO, Representation::Fourier> : public rtepack::muelmat_tensor6 {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::muelmat_tensor6::reshape;
  using rtepack::muelmat_tensor6::resize;

 public:
  using PhaseMatrixDataGridded = PhaseMatrixData<Format::ARO, Representation::Gridded>;
  using TensorType             = rtepack::muelmat_tensor6;

  PhaseMatrixData()                                   = default;
  PhaseMatrixData(const PhaseMatrixData &)            = default;
  PhaseMatrixData(PhaseMatrixData &&)                 = default;
  PhaseMatrixData &operator=(const PhaseMatrixData &) = default;
  PhaseMatrixData &operator=(PhaseMatrixData &&)      = default;

  /** Zero modes m = 0..max_mode on the given grids */
  PhaseMatrixData(std::shared_ptr<const Vector>          t_grid,
                  std::shared_ptr<const Vector>          f_grid,
                  std::shared_ptr<const Vector>          za_inc_grid,
                  std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                  Index                                  max_mode)
      : TensorType(t_grid->size(),
                   f_grid->size(),
                   za_inc_grid->size(),
                   grid_size(*za_scat_grid),
                   checked_mode(max_mode) + 1,
                   2,
                   rtepack::muelmat{0.0}),
        t_grid_(t_grid),
        f_grid_(f_grid),
        za_inc_grid_(za_inc_grid),
        za_scat_grid_(za_scat_grid),
        max_mode_(max_mode),
        phase_integral_(t_grid->size(), f_grid->size(), za_inc_grid->size(), 0.0) {}

  PhaseMatrixData &operator=(const rtepack::muelmat_tensor6 &data);

  std::shared_ptr<const Vector>          get_t_grid() const { return t_grid_; }
  std::shared_ptr<const Vector>          get_f_grid() const { return f_grid_; }
  std::shared_ptr<const Vector>          get_za_inc_grid() const { return za_inc_grid_; }
  std::shared_ptr<const ZenithAngleGrid> get_za_scat_grid() const { return za_scat_grid_; }
  Index                                  get_max_mode() const { return max_mode_; }

  //! The phase integral int Z11 dOmega [t, f, za_inc] (see the class)
  const Tensor3 &get_phase_integral() const { return phase_integral_; }

  void set_phase_integral(Tensor3 phase_integral);

  PhaseMatrixData &operator+=(const PhaseMatrixData &other);

  PhaseMatrixData &operator*=(Numeric x);

  /** The phase matrix at the azimuth differences [deg] of a grid, the series evaluated exactly */
  PhaseMatrixDataGridded to_gridded(std::shared_ptr<const Vector> delta_aa_grid) const;

  /** The phase matrix at 2 max_mode + 3 azimuth differences over [-180, 180] deg, which resolve the modes */
  PhaseMatrixDataGridded to_gridded() const {
    return to_gridded(std::make_shared<const Vector>(nlinspace(-180.0, 180.0, 2 * max_mode_ + 3)));
  }

  /** The modes truncated to a lower highest mode.
   *
   * @param max_mode At most this data's highest mode: truncated modes do not tell the higher ones
   */
  PhaseMatrixData to_fourier(Index max_mode) const;

  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  /** Interpolate the modes linearly in temperature and frequency, and select zenith angles.
   *
   * The zenith directions are not interpolated: every incidence and
   * scattering zenith angle of grids (when given) must be one of this data's,
   * whose modes and phase integral are then taken as they are.
   */
  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const;

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const;

 protected:
  static Index checked_mode(Index max_mode);

  //! The index in nodes of every angle [deg] of wanted, which must be among them
  static std::vector<Index> node_indices(const StridedConstVectorView &wanted,
                                         const StridedConstVectorView &nodes,
                                         const char                   *kind);

  std::shared_ptr<const Vector>          t_grid_;
  std::shared_ptr<const Vector>          f_grid_;
  std::shared_ptr<const Vector>          za_inc_grid_;
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid_;
  Index                                  max_mode_{0};
  //! int Z11 dOmega over all scattering directions [t, f, za_inc]
  Tensor3 phase_integral_;
};

}  // namespace scattering
