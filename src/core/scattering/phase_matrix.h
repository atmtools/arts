#pragma once

#include <arts_omp.h>
#include <legendre.h>
#include <matpack.h>
#include <rtepack.h>

#include <algorithm>
#include <limits>
#include <map>
#include <memory>

#include "arts_constants.h"
#include "integration.h"
#include "sht.h"
#include "tro_legendre.h"
#include "utils.h"

namespace scattering {

using sht::SHT;

enum class Format { TRO, ARO, General };
std::ostream &operator<<(std::ostream &out, Format format);

/** The representation of scattering data
 *
 * Gridded: on angular grids.  Spectral: spherical-harmonics (SHT) coefficients
 * over the scattering directions (Legendre coefficients for TRO).  Fourier:
 * azimuthal Fourier modes of the laboratory-frame phase matrix at given zenith
 * angles, the form plane-parallel solvers consume (ARO only).
 */
enum class Representation { Gridded, Spectral, DoublySpectral, Fourier };
std::ostream &operator<<(std::ostream &out, Representation repr);

namespace detail {

template <typename Scalar> constexpr bool equal(Scalar a, Scalar b, Scalar epsilon = 1e-6) {
  return std::abs(a - b) <= ((std::abs(a) > std::abs(b) ? std::abs(a) : std::abs(b)) * epsilon);
}

template <typename Scalar> constexpr bool small(Scalar a, Scalar epsilon = 1e-6) { return std::abs(a) < epsilon; }

/** acos of a cosine that rounding may have put just outside [-1, 1]
 *
 * Only clamps: snapping cosines within a relative 1e-6 of +-1 to 0 or pi, as
 * this did, moves angles by up to 1.4e-3 rad and makes the laboratory-frame
 * phase matrix jump there.
 */
template <typename Scalar> constexpr Scalar save_acos(Scalar a) { return std::acos(std::clamp<Scalar>(a, -1.0, 1.0)); }

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
template <typename Scalar> Scalar scattering_angle(Scalar aa_inc, Scalar za_inc, Scalar aa_scat, Scalar za_scat) {
  Scalar cos_theta = cos(za_inc) * cos(za_scat) + sin(za_inc) * sin(za_scat) * cos(aa_scat - aa_inc);
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
inline std::array<Numeric, 5> rotation_coefficients(Numeric aa_inc_d,
                                                    Numeric za_inc_d,
                                                    Numeric aa_scat_d,
                                                    Numeric za_scat_d) {
  Numeric aa_inc  = Conversion::deg2rad(aa_inc_d);
  Numeric za_inc  = Conversion::deg2rad(za_inc_d);
  Numeric aa_scat = Conversion::deg2rad(aa_scat_d);
  Numeric za_scat = Conversion::deg2rad(za_scat_d);

  Numeric cos_theta = cos(za_inc) * cos(za_scat) + sin(za_inc) * sin(za_scat) * cos(aa_scat - aa_inc);
  Numeric theta     = save_acos(cos_theta);
  if ((small(std::abs(aa_scat - aa_inc))) || (equal(std::abs(aa_scat - aa_inc), 2.0 * pi_v<Numeric>))) {
    theta = std::abs(za_inc - za_scat);
  } else if ((equal(aa_scat - aa_inc, 2.0 * pi_v<Numeric>))) {
    theta = za_scat + za_inc;
    if (theta > pi_v<Numeric>) { theta = 2.0 * pi_v<Numeric> - theta; }
  }
  const Numeric delta_aa = std::abs(std::remainder(static_cast<Numeric>(aa_scat - aa_inc), Constant::two_pi));
  const auto    ray      = [](Numeric za, Numeric aa) {
    const Numeric s = std::sin(za);
    return rtepack::meridional_basis(Vector3{s * std::sin(aa), s * std::cos(aa), std::cos(za)},
                                     Vector3{std::cos(aa), -std::sin(aa), 0.0});
  };
  const auto [c_1, s_1, c_2, s_2] = rtepack::scattering_plane_rotation(ray(static_cast<Numeric>(za_inc), 0.0),
                                                                       ray(static_cast<Numeric>(za_scat), delta_aa));

  return {theta,
          static_cast<Numeric>(c_1),
          static_cast<Numeric>(c_2),
          static_cast<Numeric>(s_1),
          static_cast<Numeric>(s_2)};
}

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

/** The scattering matrix at one scattering angle of a totally randomly oriented scatterer [t, f]
 *
 * The forward- and backscattering matrices, in the scattering-plane basis.
 */
template <Format format> class BackscatterMatrixData : public rtepack::compact_planar_muelmat_matrix {
 public:
  using TensorType = rtepack::compact_planar_muelmat_matrix;

  BackscatterMatrixData(std::shared_ptr<const Vector> t_grid, std::shared_ptr<const Vector> f_grid)
      : TensorType(t_grid->size(), f_grid->size(), rtepack::compact_planar_muelmat{}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid) {}

  BackscatterMatrixData &operator=(const TensorType &data) {
    ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                       "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_,
                       "Provided backscatter coefficient data do not match frequency grid.");
    TensorType::operator=(data);
    return *this;
  }

  std::shared_ptr<const Vector> get_t_grid_ptr() const { return t_grid_; }

  std::shared_ptr<const Vector> get_f_grid_ptr() const { return f_grid_; }

  BackscatterMatrixData extract_stokes_coeffs() const { return *this; }

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const {
    BackscatterMatrixData result(grids.t_grid, grids.f_grid);
    const TensorType     &coeffs_this = *this;
    TensorType           &coeffs_res  = result;
    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos    gp_t  = weights.t_grid_weights[i_t];
      Numeric    w_t_l = gp_t.fd[1];
      Numeric    w_t_r = gp_t.fd[0];
      const auto i_t_l = std::clamp<Index>(gp_t.idx, 0, n_temps_ - 1);
      const auto i_t_r = std::min(i_t_l + 1, n_temps_ - 1);
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos    gp_f      = weights.f_grid_weights[i_f];
        Numeric    w_f_l     = gp_f.fd[1];
        Numeric    w_f_r     = gp_f.fd[0];
        const auto i_f_l     = std::clamp<Index>(gp_f.idx, 0, n_freqs_ - 1);
        const auto i_f_r     = std::min(i_f_l + 1, n_freqs_ - 1);
        coeffs_res[i_t, i_f] = (w_t_l * w_f_l * coeffs_this[i_t_l, i_f_l] + w_t_l * w_f_r * coeffs_this[i_t_l, i_f_r] +
                                w_t_r * w_f_l * coeffs_this[i_t_r, i_f_l] + w_t_r * w_f_r * coeffs_this[i_t_r, i_f_r]);
      }
    }
    return result;
  }

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

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
};  // namespace scattering

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

  BackscatterMatrixData(const BackscatterMatrixData<Format::TRO> &bsmat, std::shared_ptr<const Vector> za_inc_grid)
      : BackscatterMatrixData(bsmat.get_t_grid_ptr(), bsmat.get_f_grid_ptr(), za_inc_grid) {
    for (Index i_t = 0; i_t < n_temps_; ++i_t)
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) this->operator[](i_t, i_f, joker) = bsmat[i_t, i_f].expand();
  }

  BackscatterMatrixData &operator=(const TensorType &data) {
    ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                       "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_,
                       "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(data.shape()[2] != n_za_inc_,
                       "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
    TensorType::operator=(data);
    return *this;
  }

  BackscatterMatrixData extract_stokes_coeffs() const { return *this; }

  BackscatterMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const {
    BackscatterMatrixData result(grids.t_grid, grids.f_grid, grids.za_inc_grid);
    const TensorType     &coeffs_this = *this;
    TensorType           &coeffs_res  = result;
    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos     gp_t    = weights.t_grid_weights[i_t];
      const Index t_upper = std::min<Index>(gp_t.idx + 1, coeffs_this.extent(0) - 1);
      Numeric     w_t_l   = gp_t.fd[1];
      Numeric     w_t_r   = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos     gp_f    = weights.f_grid_weights[i_f];
        const Index f_upper = std::min<Index>(gp_f.idx + 1, coeffs_this.extent(1) - 1);
        Numeric     w_f_l   = gp_f.fd[1];
        Numeric     w_f_r   = gp_f.fd[0];
        for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size(); ++i_za_inc) {
          GridPos     gp_za_inc    = weights.za_inc_grid_weights[i_za_inc];
          const Index za_inc_upper = std::min<Index>(gp_za_inc.idx + 1, coeffs_this.extent(2) - 1);
          Numeric     w_za_inc_l   = gp_za_inc.fd[1];
          Numeric     w_za_inc_r   = gp_za_inc.fd[0];
          coeffs_res[i_t, i_f, i_za_inc] =
              (w_t_l * w_f_l * w_za_inc_l * coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx] +
               w_t_l * w_f_l * w_za_inc_r * coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper] +
               w_t_l * w_f_r * w_za_inc_l * coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx] +
               w_t_l * w_f_r * w_za_inc_r * coeffs_this[gp_t.idx, f_upper, za_inc_upper] +

               w_t_r * w_f_l * w_za_inc_l * coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx] +
               w_t_r * w_f_l * w_za_inc_r * coeffs_this[t_upper, gp_f.idx, za_inc_upper] +
               w_t_r * w_f_r * w_za_inc_l * coeffs_this[t_upper, f_upper, gp_za_inc.idx] +
               w_t_r * w_f_r * w_za_inc_r * coeffs_this[t_upper, f_upper, za_inc_upper]);
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

template <Format format> using ForwardscatterMatrixData = BackscatterMatrixData<format>;

///////////////////////////////////////////////////////////////////////////////
// Phase matrix data
///////////////////////////////////////////////////////////////////////////////

template <Format format, Representation representation> class PhaseMatrixData;

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
template <typename ScatteringMatrix>
PhaseMatrixData<Format::ARO, Representation::Gridded> tro_lab_frame(std::shared_ptr<const Vector> t_grid,
                                                                    std::shared_ptr<const Vector> f_grid,
                                                                    std::shared_ptr<const Vector> za_inc_grid,
                                                                    std::shared_ptr<const Vector> delta_aa_grid,
                                                                    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                                                    ScatteringMatrix &&scattering_matrix);

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
template <typename ScatteringMatrix> PhaseMatrixData<Format::ARO, Representation::Fourier> tro_lab_frame_fourier_modes(
    std::shared_ptr<const Vector>          t_grid,
    std::shared_ptr<const Vector>          f_grid,
    std::shared_ptr<const Vector>          za_inc_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
    Index                                  max_mode,
    const ConstMatrixView                 &phase_integral,
    ScatteringMatrix                     &&scattering_matrix);

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

  PhaseMatrixData &operator=(const TensorType &data) {
    ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                       "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_,
                       "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(
        data.shape()[2] != n_za_scat_,
        "Provided backscatter coefficient data do not match expected number of scattering zenith angles.");
    TensorType::operator=(data);
    return *this;
  }

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

  BackscatterMatrixData<Format::TRO> extract_backscatter_matrix() {
    BackscatterMatrixData<Format::TRO> result(t_grid_, f_grid_);
    GridPos                            interp = find_interp_weights(grid_vector(*za_scat_grid_), 180.0);
    //gridpos(interp, grid_vector(*za_scat_grid_), 180.0);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        result[i_t, i_f] =
            interp.fd[1] * (*this)[i_t, i_f, interp.idx] + interp.fd[0] * (*this)[i_t, i_f, interp.idx + 1];
      }
    }
    return result;
  }

  ForwardscatterMatrixData<Format::TRO> extract_forwardscatter_matrix() {
    ForwardscatterMatrixData<Format::TRO> result(t_grid_, f_grid_);
    GridPos                               interp;
    gridpos(interp, grid_vector(*za_scat_grid_), 0.0, 1e99);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        result[i_t, i_f] =
            interp.fd[1] * (*this)[i_t, i_f, interp.idx] + interp.fd[0] * (*this)[i_t, i_f, interp.idx + 1];
      }
    }
    return result;
  }

  /** Calculate scattering-angle integral.
   *
   * Integrates the phase matrix over the scattering angles.
   * @return The integral of the phase matrix [t, f]
   */
  rtepack::compact_planar_muelmat_matrix integrate_phase_matrix() {
    rtepack::compact_planar_muelmat_matrix results(n_temps_, n_freqs_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        results[i_t, i_f] =
            2.0 * pi_v<Numeric> * integrate_zenith_angle(this->operator[](i_t, i_f, joker), *za_scat_grid_);
      }
    }
    return results;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: TRO data hold all six elements.
   */
  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const {
    PhaseMatrixData   result(grids.t_grid, grids.f_grid, grids.za_scat_grid);
    const TensorType &coeffs_this = *this;
    TensorType       &coeffs_res  = result;
    for (Index i_t = 0; i_t < static_cast<Index>(weights.t_grid_weights.size()); ++i_t) {
      GridPos     gp_t  = weights.t_grid_weights[i_t];
      Numeric     w_t_l = gp_t.fd[1];
      Numeric     w_t_r = gp_t.fd[0];
      const Index it0   = std::clamp<Index>(gp_t.idx, 0, n_temps_ - 1);
      const Index it1   = std::min(it0 + 1, n_temps_ - 1);
      for (Index i_f = 0; i_f < static_cast<Index>(weights.f_grid_weights.size()); ++i_f) {
        GridPos     gp_f  = weights.f_grid_weights[i_f];
        Numeric     w_f_l = gp_f.fd[1];
        Numeric     w_f_r = gp_f.fd[0];
        const Index if0   = std::clamp<Index>(gp_f.idx, 0, n_freqs_ - 1);
        const Index if1   = std::min(if0 + 1, n_freqs_ - 1);
        for (Index i_za_scat = 0; i_za_scat < static_cast<Index>(weights.za_scat_grid_weights.size()); ++i_za_scat) {
          GridPos     gp_za_scat  = weights.za_scat_grid_weights[i_za_scat];
          Numeric     w_za_scat_l = gp_za_scat.fd[1];
          Numeric     w_za_scat_r = gp_za_scat.fd[0];
          const Index iza0        = std::clamp<Index>(gp_za_scat.idx, 0, n_za_scat_ - 1);
          const Index iza1        = std::min(iza0 + 1, n_za_scat_ - 1);

          coeffs_res[i_t, i_f, i_za_scat] = rtepack::compact_planar_muelmat{};

          if (w_t_l > 0.0) {
            if (w_f_l > 0.0) {
              coeffs_res[i_t, i_f, i_za_scat] += w_t_l * w_f_l * w_za_scat_l * coeffs_this[it0, if0, iza0];
              coeffs_res[i_t, i_f, i_za_scat] += w_t_l * w_f_l * w_za_scat_r * coeffs_this[it0, if0, iza1];
            }
            if (w_f_r > 0.0) {
              coeffs_res[i_t, i_f, i_za_scat] += w_t_l * w_f_r * w_za_scat_l * coeffs_this[it0, if1, iza0];
              coeffs_res[i_t, i_f, i_za_scat] += w_t_l * w_f_r * w_za_scat_r * coeffs_this[it0, if1, iza1];
            }
          }
          if (w_t_r > 0.0) {
            if (w_f_l > 0.0) {
              coeffs_res[i_t, i_f, i_za_scat] += w_t_r * w_f_l * w_za_scat_l * coeffs_this[it1, if0, iza0];
              coeffs_res[i_t, i_f, i_za_scat] += w_t_r * w_f_l * w_za_scat_r * coeffs_this[it1, if0, iza1];
            }
            if (w_f_r > 0.0) {
              coeffs_res[i_t, i_f, i_za_scat] += w_t_r * w_f_r * w_za_scat_l * coeffs_this[it1, if1, iza0];
              coeffs_res[i_t, i_f, i_za_scat] += w_t_r * w_f_r * w_za_scat_r * coeffs_this[it1, if1, iza1];
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

  PhaseMatrixData &operator=(const TensorType &data) {
    ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                       "Provided phase matrix coefficients do not match temperature grid.");
    ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_, "Provided phase matrix coefficients do not match frequency grid.");
    ARTS_USER_ERROR_IF(data.shape()[2] != degree_ + 1,
                       "Provided phase matrix coefficients do not match the Legendre degree {}.",
                       degree_);
    TensorType::operator=(data);
    return *this;
  }

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
  PhaseMatrixDataGridded to_gridded(std::shared_ptr<const ZenithAngleGrid> za_scat_grid) const {
    PhaseMatrixDataGridded result(t_grid_, f_grid_, za_scat_grid);
    const Vector           angles{grid_vector(*za_scat_grid)};
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        result[i_t, i_f, joker] = tro_legendre::evaluate((*this)[i_t, i_f, joker], angles);
      }
    }
    return result;
  }

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
  auto scattering_matrix_function() const {
    Vector norm(degree_ + 1);
    for (Index l = 0; l <= degree_; ++l) norm[l] = std::sqrt(static_cast<Numeric>(2 * l + 1) / (4.0 * pi_v<Numeric>));
    return [this, p = Vector(degree_ + 1), norm = std::move(norm)](
               Numeric theta, rtepack::compact_planar_muelmat_matrix &scattering_matrix) mutable {
      Legendre::legendre_polynomials(p, std::clamp<Numeric>(std::cos(theta), -1.0, 1.0));
      for (Index l = 0; l <= degree_; ++l) p[l] *= norm[l];
      for (Index i_t = 0; i_t < n_temps_; ++i_t) {
        for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
          rtepack::compact_planar_muelmat value{};
          for (Index l = 0; l <= degree_; ++l) value += (*this)[i_t, i_f, l] * p[l];
          scattering_matrix[i_t, i_f] = value;
        }
      }
    };
  }

  BackscatterMatrixData<Format::TRO> extract_backscatter_matrix() const { return extract_at(180.0); }

  ForwardscatterMatrixData<Format::TRO> extract_forwardscatter_matrix() const { return extract_at(0.0); }

  /** Calculate scattering-angle integral.
   *
   * Integrates the phase matrix over the scattering angles: 2 pi int F dcos(Theta) = sqrt(4 pi) a_0.
   * @return The integral of the phase matrix [t, f]
   */
  rtepack::compact_planar_muelmat_matrix integrate_phase_matrix() const {
    rtepack::compact_planar_muelmat_matrix results(n_temps_, n_freqs_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        results[i_t, i_f] = this->operator[](i_t, i_f, 0) * sqrt(4.0 * pi_v<Numeric>);
      }
    }
    return results;
  }

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
  PhaseMatrixData to_spectral(Index degree, Index order = 0) const {
    ARTS_USER_ERROR_IF(order != 0,
                       "A TRO phase matrix depends on the scattering angle only, so its spectral form has order 0; "
                       "got order {}",
                       order)
    ARTS_USER_ERROR_IF(degree > degree_,
                       "The Legendre series has degree {}, so it cannot give the coefficients to degree {}; convert "
                       "the gridded data to that degree instead",
                       degree_,
                       degree)
    PhaseMatrixData out(t_grid_, f_grid_, degree);
    for (Index i_t = 0; i_t < n_temps_; ++i_t)
      for (Index i_f = 0; i_f < n_freqs_; ++i_f)
        for (Index l = 0; l <= degree; ++l) out[i_t, i_f, l] = this->operator[](i_t, i_f, l);
    return out;
  }

  /** Interpolate the coefficients linearly in temperature and frequency.
   *
   * The angular dependence is the series and needs no interpolation; the
   * zenith-angle grids of grids are not used.
   */
  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights weights) const {
    PhaseMatrixData   result(grids.t_grid, grids.f_grid, degree_);
    const TensorType &coeffs_this = *this;
    TensorType       &coeffs_res  = result;
    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos     gp_t    = weights.t_grid_weights[i_t];
      const Index t_upper = std::min<Index>(gp_t.idx + 1, coeffs_this.extent(0) - 1);
      Numeric     w_t_l   = gp_t.fd[1];
      Numeric     w_t_r   = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos     gp_f    = weights.f_grid_weights[i_f];
        const Index f_upper = std::min<Index>(gp_f.idx + 1, coeffs_this.extent(1) - 1);
        Numeric     w_f_l   = gp_f.fd[1];
        Numeric     w_f_r   = gp_f.fd[0];

        for (Index l = 0; l <= degree_; ++l) {
          coeffs_res[i_t, i_f, l] = rtepack::compact_planar_muelmat{};
          if (w_t_l > 0.0) {
            if (w_f_l > 0.0) coeffs_res[i_t, i_f, l] += w_t_l * w_f_l * coeffs_this[gp_t.idx, gp_f.idx, l];
            if (w_f_r > 0.0) coeffs_res[i_t, i_f, l] += w_t_l * w_f_r * coeffs_this[gp_t.idx, f_upper, l];
          }
          if (w_t_r > 0.0) {
            if (w_f_l > 0.0) coeffs_res[i_t, i_f, l] += w_t_r * w_f_l * coeffs_this[t_upper, gp_f.idx, l];
            if (w_f_r > 0.0) coeffs_res[i_t, i_f, l] += w_t_r * w_f_r * coeffs_this[t_upper, f_upper, l];
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
  static Index checked_degree(Index degree) {
    ARTS_USER_ERROR_IF(degree < 0, "The Legendre degree must be >= 0, got {}", degree)
    return degree;
  }

  /** The series at one scattering angle [deg] */
  BackscatterMatrixData<Format::TRO> extract_at(Numeric angle) const {
    BackscatterMatrixData<Format::TRO> result(t_grid_, f_grid_);
    const Vector                       angles{angle};
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        result[i_t, i_f] = tro_legendre::evaluate((*this)[i_t, i_f, joker], angles)[0];
      }
    }
    return result;
  }

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

  PhaseMatrixData &operator=(const TensorType &data) {
    ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                       "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_,
                       "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(data.shape()[2] != n_za_inc_,
                       "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
    ARTS_USER_ERROR_IF(
        data.shape()[3] != n_delta_aa_,
        "Provided backscatter coefficient data do not match expected number of scattering azimuth angles.");
    ARTS_USER_ERROR_IF(
        data.shape()[4] != n_za_scat_,
        "Provided backscatter coefficient data do not match expected number of scattering zenith angles.");
    TensorType::operator=(data);
    return *this;
  }

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

  BackscatterMatrixData<Format::ARO> extract_backscatter_matrix() {
    BackscatterMatrixData<Format::ARO> result(t_grid_, f_grid_, za_inc_grid_);
    GridPos                            za_scat_interp, delta_aa_interp;
    gridpos(delta_aa_interp, *delta_aa_grid_, 180.0, 1e99);
    Vector     weights(4);
    const auto elements = rtepack::elements(*this);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          auto za_inc = (*za_inc_grid_)[i_za_inc];
          gridpos(za_scat_interp, grid_vector(*za_scat_grid_), 180.0 - za_inc, 1e99);
          interpweights(weights, delta_aa_interp, za_scat_interp);

          for (Index i_s = 0; i_s < 16; ++i_s) {
            auto mat                             = elements[i_t, i_f, i_za_inc, joker, joker, i_s];
            result[i_t, i_f, i_za_inc].data[i_s] = interp(weights, mat, delta_aa_interp, za_scat_interp);
          }
        }
      }
    }
    return result;
  }

  ForwardscatterMatrixData<Format::ARO> extract_forwardscatter_matrix() {
    BackscatterMatrixData<Format::ARO> result(t_grid_, f_grid_, za_inc_grid_);
    GridPos                            za_scat_interp, delta_aa_interp;
    gridpos(delta_aa_interp, *delta_aa_grid_, 0.0, 1e99);
    Vector     weights(4);
    const auto elements = rtepack::elements(*this);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          auto za_inc = (*za_inc_grid_)[i_za_inc];
          gridpos(za_scat_interp, grid_vector(*za_scat_grid_), za_inc, 1e99);
          interpweights(weights, delta_aa_interp, za_scat_interp);

          for (Index i_s = 0; i_s < 16; ++i_s) {
            auto mat                             = elements[i_t, i_f, i_za_inc, joker, joker, i_s];
            result[i_t, i_f, i_za_inc].data[i_s] = interp(weights, mat, delta_aa_interp, za_scat_interp);
          }
        }
      }
    }
    return result;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: ARO data hold all sixteen elements.
   */
  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const {
    PhaseMatrixData   result(grids.t_grid, grids.f_grid, grids.za_inc_grid, grids.aa_scat_grid, grids.za_scat_grid);
    const TensorType &coeffs_this = *this;
    TensorType       &coeffs_res  = result;

    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos     gp_t    = weights.t_grid_weights[i_t];
      const Index t_upper = std::min<Index>(gp_t.idx + 1, coeffs_this.extent(0) - 1);
      Numeric     w_t_l   = gp_t.fd[1];
      Numeric     w_t_r   = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos     gp_f    = weights.f_grid_weights[i_f];
        const Index f_upper = std::min<Index>(gp_f.idx + 1, coeffs_this.extent(1) - 1);
        Numeric     w_f_l   = gp_f.fd[1];
        Numeric     w_f_r   = gp_f.fd[0];
        for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size(); ++i_za_inc) {
          GridPos     gp_za_inc    = weights.za_inc_grid_weights[i_za_inc];
          const Index za_inc_upper = std::min<Index>(gp_za_inc.idx + 1, coeffs_this.extent(2) - 1);
          Numeric     w_za_inc_l   = gp_za_inc.fd[1];
          Numeric     w_za_inc_r   = gp_za_inc.fd[0];
          for (Size i_aa_scat = 0; i_aa_scat < weights.aa_scat_grid_weights.size(); ++i_aa_scat) {
            GridPos     gp_aa_scat    = weights.aa_scat_grid_weights[i_aa_scat];
            const Index aa_scat_upper = std::min<Index>(gp_aa_scat.idx + 1, coeffs_this.extent(3) - 1);
            Numeric     w_aa_scat_l   = gp_aa_scat.fd[1];
            Numeric     w_aa_scat_r   = gp_aa_scat.fd[0];
            for (Size i_za_scat = 0; i_za_scat < weights.za_scat_grid_weights.size(); ++i_za_scat) {
              GridPos     gp_za_scat    = weights.za_scat_grid_weights[i_za_scat];
              const Index za_scat_upper = std::min<Index>(gp_za_scat.idx + 1, coeffs_this.extent(4) - 1);
              Numeric     w_za_scat_l   = gp_za_scat.fd[1];
              Numeric     w_za_scat_r   = gp_za_scat.fd[0];

              coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] = rtepack::muelmat{0.0};

              if (w_t_l > 0.0) {
                if (w_f_l > 0.0) {
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, gp_aa_scat.idx, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, gp_aa_scat.idx, za_scat_upper];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, aa_scat_upper, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, aa_scat_upper, za_scat_upper];

                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, gp_aa_scat.idx, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, gp_aa_scat.idx, za_scat_upper];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, aa_scat_upper, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, aa_scat_upper, za_scat_upper];
                }
                if (w_f_r > 0.0) {
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, gp_aa_scat.idx, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, gp_aa_scat.idx, za_scat_upper];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, aa_scat_upper, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, aa_scat_upper, za_scat_upper];

                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[gp_t.idx, f_upper, za_inc_upper, gp_aa_scat.idx, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[gp_t.idx, f_upper, za_inc_upper, gp_aa_scat.idx, za_scat_upper];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[gp_t.idx, f_upper, za_inc_upper, aa_scat_upper, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_l * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[gp_t.idx, f_upper, za_inc_upper, aa_scat_upper, za_scat_upper];
                }
              }
              if (w_t_r > 0.0) {
                if (w_f_l > 0.0) {
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, gp_aa_scat.idx, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, gp_aa_scat.idx, za_scat_upper];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, aa_scat_upper, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, aa_scat_upper, za_scat_upper];

                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[t_upper, gp_f.idx, za_inc_upper, gp_aa_scat.idx, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[t_upper, gp_f.idx, za_inc_upper, gp_aa_scat.idx, za_scat_upper];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[t_upper, gp_f.idx, za_inc_upper, aa_scat_upper, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[t_upper, gp_f.idx, za_inc_upper, aa_scat_upper, za_scat_upper];
                }
                if (w_f_r > 0.0) {
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[t_upper, f_upper, gp_za_inc.idx, gp_aa_scat.idx, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[t_upper, f_upper, gp_za_inc.idx, gp_aa_scat.idx, za_scat_upper];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[t_upper, f_upper, gp_za_inc.idx, aa_scat_upper, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[t_upper, f_upper, gp_za_inc.idx, aa_scat_upper, za_scat_upper];

                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                      coeffs_this[t_upper, f_upper, za_inc_upper, gp_aa_scat.idx, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                      coeffs_this[t_upper, f_upper, za_inc_upper, gp_aa_scat.idx, za_scat_upper];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                      coeffs_this[t_upper, f_upper, za_inc_upper, aa_scat_upper, gp_za_scat.idx];
                  coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                      w_t_r * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                      coeffs_this[t_upper, f_upper, za_inc_upper, aa_scat_upper, za_scat_upper];
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
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, za_inc_grid_, delta_aa_grid_, za_scat_grid_, grids);
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

template <typename ScatteringMatrix>
PhaseMatrixData<Format::ARO, Representation::Gridded> tro_lab_frame(std::shared_ptr<const Vector> t_grid,
                                                                    std::shared_ptr<const Vector> f_grid,
                                                                    std::shared_ptr<const Vector> za_inc_grid,
                                                                    std::shared_ptr<const Vector> delta_aa_grid,
                                                                    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                                                    ScatteringMatrix &&scattering_matrix) {
  PhaseMatrixData<Format::ARO, Representation::Gridded> result(
      t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);

  const auto                             za_scat = grid_vector(*za_scat_grid);
  rtepack::compact_planar_muelmat_matrix tro(t_grid->size(), f_grid->size());
  for (Size i_za_inc = 0; i_za_inc < za_inc_grid->size(); ++i_za_inc) {
    for (Size i_delta_aa = 0; i_delta_aa < delta_aa_grid->size(); ++i_delta_aa) {
      // In [0, 360) deg, so that delta_aa > 180 tells the side of the principal plane for any grid, e.g. [-180, 180]
      Numeric delta_aa = std::fmod((*delta_aa_grid)[i_delta_aa], Numeric{360});
      if (delta_aa < 0) delta_aa += 360;
      for (Size i_za_scat = 0; i_za_scat < za_scat.size(); ++i_za_scat) {
        const auto [theta, c_1, c_2, s_1, s_2] =
            detail::rotation_coefficients(0.0, (*za_inc_grid)[i_za_inc], delta_aa, za_scat[i_za_scat]);
        scattering_matrix(theta, tro);
        for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
          for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
            // Beyond 180 deg, the mirror image of the other side of the principal plane
            const rtepack::muelmat Z                          = rtepack::rotated(tro[i_t, i_f], c_1, s_1, c_2, s_2);
            result[i_t, i_f, i_za_inc, i_delta_aa, i_za_scat] = delta_aa > 180.0 ? rtepack::mirror(Z) : Z;
          }
        }
      }
    }
  }
  return result;
}

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

  PhaseMatrixData &operator=(const TensorType &data) {
    ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                       "Provided backscatter coefficient data do not match temperature grid.");
    ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_,
                       "Provided backscatter coefficient data do not match frequency grid.");
    ARTS_USER_ERROR_IF(data.shape()[2] != n_za_inc_,
                       "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
    ARTS_USER_ERROR_IF(data.shape()[3] != n_spectral_coeffs_,
                       "Provided backscatter coefficient data do not match expected number of SHT coefficients.");
    TensorType::operator=(data);
    return *this;
  }

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
  PhaseMatrixDataGridded to_gridded() const {
    PhaseMatrixDataGridded result(t_grid_, f_grid_, za_inc_grid_, sht_->get_aa_grid_ptr(), sht_->get_za_grid_ptr());
    auto                   elements = rtepack::elements(result);
    const auto             coeffs   = rtepack::elements(*this);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          for (Index i_s = 0; i_s < 16; ++i_s) {
            elements[i_t, i_f, i_za_inc, joker, joker, i_s] = sht_->synthesize(coeffs[i_t, i_f, i_za_inc, joker, i_s]);
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
    auto            sht_new = sht::provider.get_instance_lm(l_new, m_new);
    PhaseMatrixData pm_new(t_grid_, f_grid_, za_inc_grid_, sht_new);
    auto            coeffs_new = rtepack::elements(pm_new);
    const auto      coeffs     = rtepack::elements(*this);
    for (Size f_ind = 0; f_ind < f_grid_->size(); ++f_ind) {
      for (Size t_ind = 0; t_ind < t_grid_->size(); ++t_ind) {
        for (Size za_inc_ind = 0; za_inc_ind < za_inc_grid_->size(); ++za_inc_ind) {
          for (Index i_s = 0; i_s < 16; ++i_s) {
            coeffs_new[t_ind, f_ind, za_inc_ind, joker, i_s] =
                sht::add_coeffs(*sht_new,
                                coeffs_new[t_ind, f_ind, za_inc_ind, joker, i_s],
                                *sht_,
                                coeffs[t_ind, f_ind, za_inc_ind, joker, i_s]);
          }
        }
      }
    }
    return pm_new;
  }

  BackscatterMatrixData<Format::ARO> extract_backscatter_matrix() const {
    BackscatterMatrixData<Format::ARO> result(t_grid_, f_grid_, za_inc_grid_);
    const auto                         elements = rtepack::elements(*this);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          auto    za_inc  = (*za_inc_grid_)[i_za_inc];
          Numeric za_scat = 180.0 - za_inc;
          for (Index i_s = 0; i_s < 16; ++i_s) {
            auto coeffs = elements[i_t, i_f, i_za_inc, joker, i_s];
            result[i_t, i_f, i_za_inc].data[i_s] =
                sht_->evaluate(coeffs, Conversion::deg2rad(180.0), Conversion::deg2rad(za_scat));
          }
        }
      }
    }
    return result;
  }

  ForwardscatterMatrixData<Format::ARO> extract_forwardscatter_matrix() {
    BackscatterMatrixData<Format::ARO> result(t_grid_, f_grid_, za_inc_grid_);
    const auto                         elements = rtepack::elements(*this);

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          auto za_inc = (*za_inc_grid_)[i_za_inc];
          for (Index i_s = 0; i_s < 16; ++i_s) {
            auto coeffs = elements[i_t, i_f, i_za_inc, joker, i_s];
            result[i_t, i_f, i_za_inc].data[i_s] =
                sht_->evaluate(coeffs, Conversion::deg2rad(0.0), Conversion::deg2rad(za_inc));
          }
        }
      }
    }
    return result;
  }

  /** The integral over all scattering directions of the phase matrix [t, f, za_inc]
   *
   * sqrt(4 pi) times the l = 0 coefficient on the orthonormal Y_00; a
   * degree-0 SHT holds the value itself as its coefficient.
   */
  rtepack::muelmat_tensor3 integrate_phase_matrix() const {
    const Numeric            factor = sht_->get_l_max() == 0 ? 4.0 * pi_v<Numeric> : std::sqrt(4.0 * pi_v<Numeric>);
    rtepack::muelmat_tensor3 results(n_temps_, n_freqs_, n_za_inc_);
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
          const rtepack::specmat &a_0 = this->operator[](i_t, i_f, i_za_inc, 0);
          for (Index i_s = 0; i_s < 16; ++i_s) {
            results[i_t, i_f, i_za_inc].data[i_s] = a_0.data[i_s].real() * factor;
          }
        }
      }
    }
    return results;
  }

  /** Extract single scattering data for given stokes dimension.
   *
   * @return A copy: ARO data hold all sixteen elements.
   */
  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights weights) const {
    PhaseMatrixData   result(grids.t_grid, grids.f_grid, grids.za_inc_grid, sht_);
    const TensorType &coeffs_this = *this;
    TensorType       &coeffs_res  = result;

    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      GridPos     gp_t    = weights.t_grid_weights[i_t];
      const Index t_upper = std::min<Index>(gp_t.idx + 1, coeffs_this.extent(0) - 1);
      Numeric     w_t_l   = gp_t.fd[1];
      Numeric     w_t_r   = gp_t.fd[0];
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        GridPos     gp_f    = weights.f_grid_weights[i_f];
        const Index f_upper = std::min<Index>(gp_f.idx + 1, coeffs_this.extent(1) - 1);
        Numeric     w_f_l   = gp_f.fd[1];
        Numeric     w_f_r   = gp_f.fd[0];
        for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size(); ++i_za_inc) {
          GridPos     gp_za_inc    = weights.za_inc_grid_weights[i_za_inc];
          const Index za_inc_upper = std::min<Index>(gp_za_inc.idx + 1, coeffs_this.extent(2) - 1);
          Numeric     w_za_inc_l   = gp_za_inc.fd[1];
          Numeric     w_za_inc_r   = gp_za_inc.fd[0];

          for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
            coeffs_res[i_t, i_f, i_za_inc, i_sht] = rtepack::specmat{0.0};
          }

          if (w_t_l > 0.0) {
            if (w_f_l > 0.0) {
              for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_l * w_f_l * w_za_inc_l * coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, i_sht];
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_l * w_f_l * w_za_inc_r * coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, i_sht];
              }
            }
            if (w_f_r > 0.0) {
              for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_l * w_f_r * w_za_inc_l * coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, i_sht];
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_l * w_f_r * w_za_inc_r * coeffs_this[gp_t.idx, f_upper, za_inc_upper, i_sht];
              }
            }
          }
          if (w_t_r > 0.0) {
            if (w_f_l > 0.0) {
              for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_r * w_f_l * w_za_inc_l * coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, i_sht];
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_r * w_f_l * w_za_inc_r * coeffs_this[t_upper, gp_f.idx, za_inc_upper, i_sht];
              }
            }
            if (w_f_r > 0.0) {
              for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_r * w_f_r * w_za_inc_l * coeffs_this[t_upper, f_upper, gp_za_inc.idx, i_sht];
                coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                    w_t_r * w_f_r * w_za_inc_r * coeffs_this[t_upper, f_upper, za_inc_upper, i_sht];
              }
            }
          }
        }
      }
    }
    return result;
  }

  PhaseMatrixData regrid(const ScatteringDataGrids &grids) const {
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

  PhaseMatrixData &operator=(const rtepack::muelmat_tensor6 &data) {
    ARTS_USER_ERROR_IF(data.shape() != this->shape(),
                       "The Fourier modes have shape {:B,}, but these grids need {:B,}",
                       data.shape(),
                       this->shape())
    TensorType::operator=(data);
    return *this;
  }

  std::shared_ptr<const Vector>          get_t_grid() const { return t_grid_; }
  std::shared_ptr<const Vector>          get_f_grid() const { return f_grid_; }
  std::shared_ptr<const Vector>          get_za_inc_grid() const { return za_inc_grid_; }
  std::shared_ptr<const ZenithAngleGrid> get_za_scat_grid() const { return za_scat_grid_; }
  Index                                  get_max_mode() const { return max_mode_; }

  //! The phase integral int Z11 dOmega [t, f, za_inc] (see the class)
  const Tensor3 &get_phase_integral() const { return phase_integral_; }

  void set_phase_integral(Tensor3 phase_integral) {
    ARTS_USER_ERROR_IF(phase_integral.shape() != phase_integral_.shape(),
                       "The phase integral has shape {:B,}, but these grids need {:B,}",
                       phase_integral.shape(),
                       phase_integral_.shape())
    phase_integral_ = std::move(phase_integral);
  }

  PhaseMatrixData &operator+=(const PhaseMatrixData &other) {
    ARTS_USER_ERROR_IF(other.shape() != this->shape(),
                       "Cannot add Fourier modes of shape {:B,} to modes of shape {:B,}",
                       other.shape(),
                       this->shape())
    static_cast<TensorType &>(*this) += static_cast<const TensorType &>(other);
    phase_integral_                  += other.phase_integral_;
    return *this;
  }

  PhaseMatrixData &operator*=(Numeric x) {
    static_cast<TensorType &>(*this) *= x;
    phase_integral_                  *= x;
    return *this;
  }

  /** The phase matrix at the azimuth differences [deg] of a grid, the series evaluated exactly */
  PhaseMatrixDataGridded to_gridded(std::shared_ptr<const Vector> delta_aa_grid) const {
    PhaseMatrixDataGridded result(t_grid_, f_grid_, za_inc_grid_, delta_aa_grid, za_scat_grid_);
    const auto [nt, nf, ni, ns, nm, ncs] = this->shape();
    for (Size k = 0; k < delta_aa_grid->size(); ++k) {
      const Numeric delta = Conversion::deg2rad((*delta_aa_grid)[k]);
      for (Index m = 0; m < nm; ++m) {
        const Numeric c = std::cos(static_cast<Numeric>(m) * delta), s = std::sin(static_cast<Numeric>(m) * delta);
        for (Index i_t = 0; i_t < nt; ++i_t)
          for (Index i_f = 0; i_f < nf; ++i_f)
            for (Index i_i = 0; i_i < ni; ++i_i)
              for (Index i_s = 0; i_s < ns; ++i_s) {
                const rtepack::muelmat &C = this->operator[](i_t, i_f, i_i, i_s, m, 0);
                const rtepack::muelmat &S = this->operator[](i_t, i_f, i_i, i_s, m, 1);
                for (Index e = 0; e < 16; ++e) result[i_t, i_f, i_i, k, i_s].data[e] += c * C.data[e] + s * S.data[e];
              }
      }
    }
    return result;
  }

  /** The phase matrix at 2 max_mode + 3 azimuth differences over [-180, 180] deg, which resolve the modes */
  PhaseMatrixDataGridded to_gridded() const {
    return to_gridded(std::make_shared<const Vector>(nlinspace(-180.0, 180.0, 2 * max_mode_ + 3)));
  }

  /** The modes truncated to a lower highest mode.
   *
   * @param max_mode At most this data's highest mode: truncated modes do not tell the higher ones
   */
  PhaseMatrixData to_fourier(Index max_mode) const {
    ARTS_USER_ERROR_IF(max_mode > max_mode_,
                       "The Fourier modes go to m = {}, so they cannot give the modes to m = {}; convert the gridded "
                       "data to that mode instead",
                       max_mode_,
                       max_mode)
    PhaseMatrixData out(t_grid_, f_grid_, za_inc_grid_, za_scat_grid_, max_mode);
    const auto [nt, nf, ni, ns, nm, ncs] = out.shape();
    for (Index i_t = 0; i_t < nt; ++i_t)
      for (Index i_f = 0; i_f < nf; ++i_f)
        for (Index i_i = 0; i_i < ni; ++i_i)
          for (Index i_s = 0; i_s < ns; ++i_s)
            for (Index m = 0; m < nm; ++m) out[i_t, i_f, i_i, i_s, m] = this->operator[](i_t, i_f, i_i, i_s, m);
    out.phase_integral_ = phase_integral_;
    return out;
  }

  PhaseMatrixData extract_stokes_coeffs() const { return *this; }

  /** Interpolate the modes linearly in temperature and frequency, and select zenith angles.
   *
   * The zenith directions are not interpolated: every incidence and
   * scattering zenith angle of grids (when given) must be one of this data's,
   * whose modes and phase integral are then taken as they are.
   */
  PhaseMatrixData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const {
    auto            new_inc  = grids.za_inc_grid ? grids.za_inc_grid : za_inc_grid_;
    auto            new_scat = grids.za_scat_grid ? grids.za_scat_grid : za_scat_grid_;
    const auto      inc      = node_indices(*new_inc, *za_inc_grid_, "incidence");
    const auto      scat     = node_indices(grid_vector(*new_scat), grid_vector(*za_scat_grid_), "scattering");
    PhaseMatrixData result(grids.t_grid, grids.f_grid, new_inc, new_scat, max_mode_);
    const Index     nt = this->extent(0), nf = this->extent(1);
    for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
      const GridPos gp_t = weights.t_grid_weights[i_t];
      const Index   t0 = std::clamp<Index>(gp_t.idx, 0, nt - 1), t1 = std::min<Index>(t0 + 1, nt - 1);
      for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
        const GridPos gp_f = weights.f_grid_weights[i_f];
        const Index   f0 = std::clamp<Index>(gp_f.idx, 0, nf - 1), f1 = std::min<Index>(f0 + 1, nf - 1);
        const std::array<std::tuple<Numeric, Index, Index>, 4> corners{{{gp_t.fd[1] * gp_f.fd[1], t0, f0},
                                                                        {gp_t.fd[1] * gp_f.fd[0], t0, f1},
                                                                        {gp_t.fd[0] * gp_f.fd[1], t1, f0},
                                                                        {gp_t.fd[0] * gp_f.fd[0], t1, f1}}};
        for (const auto &[w, it, jf] : corners) {
          if (w == 0.0) continue;
          for (Size i_i = 0; i_i < inc.size(); ++i_i) {
            result.phase_integral_[i_t, i_f, i_i] += w * phase_integral_[it, jf, inc[i_i]];
            for (Size i_s = 0; i_s < scat.size(); ++i_s) {
              auto       out = result[i_t, i_f, i_i, i_s];
              const auto in  = this->operator[](it, jf, inc[i_i], scat[i_s]);
              for (auto [o, x] : std::views::zip(out | by_elem, in | by_elem))
                for (auto [oe, xe] : std::views::zip(o.data, x.data)) oe += w * xe;
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
  static Index checked_mode(Index max_mode) {
    ARTS_USER_ERROR_IF(max_mode < 0, "The highest Fourier mode must be >= 0, got {}", max_mode)
    return max_mode;
  }

  //! The index in nodes of every angle [deg] of wanted, which must be among them
  static std::vector<Index> node_indices(const StridedConstVectorView &wanted,
                                         const StridedConstVectorView &nodes,
                                         const char                   *kind) {
    std::vector<Index> out(wanted.size());
    for (Size i = 0; i < wanted.size(); ++i) {
      const auto it = stdr::find_if(nodes, [x = wanted[i]](Numeric y) { return std::abs(x - y) <= 1e-9; });
      ARTS_USER_ERROR_IF(it == nodes.end(),
                         "ARO Fourier modes are known on their own {} zenith angles only, {:B,} deg, which do not "
                         "include {} deg; give the data on the wanted angles",
                         kind,
                         nodes,
                         wanted[i])
      out[i] = static_cast<Index>(it - nodes.begin());
    }
    return out;
  }

  std::shared_ptr<const Vector>          t_grid_;
  std::shared_ptr<const Vector>          f_grid_;
  std::shared_ptr<const Vector>          za_inc_grid_;
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid_;
  Index                                  max_mode_{0};
  //! int Z11 dOmega over all scattering directions [t, f, za_inc]
  Tensor3 phase_integral_;
};

template <typename ScatteringMatrix> PhaseMatrixData<Format::ARO, Representation::Fourier> tro_lab_frame_fourier_modes(
    std::shared_ptr<const Vector>          t_grid,
    std::shared_ptr<const Vector>          f_grid,
    std::shared_ptr<const Vector>          za_inc_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
    Index                                  max_mode,
    const ConstMatrixView                 &phase_integral,
    ScatteringMatrix                     &&scattering_matrix) {
  PhaseMatrixData<Format::ARO, Representation::Fourier> out(t_grid, f_grid, za_inc_grid, za_scat_grid, max_mode);
  ARTS_USER_ERROR_IF(phase_integral.nrows() != static_cast<Index>(t_grid->size()) or
                         phase_integral.ncols() != static_cast<Index>(f_grid->size()),
                     "The phase integral has shape [{}, {}], but there are {} temperatures and {} frequencies",
                     phase_integral.nrows(),
                     phase_integral.ncols(),
                     t_grid->size(),
                     f_grid->size()) {
    Tensor3 integral(t_grid->size(), f_grid->size(), za_inc_grid->size());
    for (Index i_t = 0; i_t < integral.extent(0); ++i_t)
      for (Index i_f = 0; i_f < integral.extent(1); ++i_f) integral[i_t, i_f, joker] = phase_integral[i_t, i_f];
    out.set_phase_integral(std::move(integral));
  }

  const Vector za_scat{grid_vector(*za_scat_grid)};
  for (const Numeric za : *za_inc_grid)
    ARTS_USER_ERROR_IF(not(za >= 0.0 and za <= 180.0), "Incidence zenith angles must be in [0, 180] deg, got {}", za)
  for (const Numeric za : za_scat)
    ARTS_USER_ERROR_IF(not(za >= 0.0 and za <= 180.0), "Scattering zenith angles must be in [0, 180] deg, got {}", za)

  constexpr Index   max_nodes = 4096;
  constexpr Numeric tolerance = 1e-13;
  const Index       nt = t_grid->size(), nf = f_grid->size(), M = max_mode, nset = nt * nf;

  // Z(2 pi - delta) = D Z(delta) D with D = diag(1, 1, -1, -1), the side flip rtepack::mirror: the
  // elements with exactly one index in {U, V} are odd in delta and have sine modes only, the others cosine modes
  // only.  So [0, pi] gives the modes: C_m = (2 - delta_m0) / 2 sum_k w_k Z(phi_k) cos(m phi_k) and S_m = sum_k
  // w_k Z(phi_k) sin(m phi_k) for the n-point Gauss-Legendre rule (w_k, phi_k) on [0, pi], whose ends are the
  // only directions where Z can have a kink (exact forward and backward scattering).
  constexpr std::array<bool, 16> odd{
      false, false, true, true, false, false, true, true, true, true, false, false, true, true, false, false};

  //! The rule of n nodes on [0, pi], and its weights of every mode, [n, M + 1]
  struct rule {
    Vector phi;
    Matrix cos_weight, sin_weight;
  };

  Index start = 16;
  while (start < M + 8) start *= 2;

  // [set, m, (cos, sin), element]
  const Index nmode = (M + 1) * 2 * 16;

  // The direction pairs are independent; each thread has its own copy of the scattering matrix (which may keep
  // scratch), its own rules and its own work arrays, and every pair writes its own part of out
  std::string error;
#pragma omp parallel if (not arts_omp_in_parallel())
  {
    auto                                   local_matrix = scattering_matrix;
    std::map<Index, rule>                  rules;
    rtepack::compact_planar_muelmat_matrix F(nt, nf);
    std::vector<Numeric>                   coarse(nset * nmode), fine(nset * nmode);

    const auto get_rule = [&](Index n) -> const rule & {
      auto [it, inserted] = rules.try_emplace(n);
      if (inserted) {
        Vector x(n), w(n);
        Legendre::GaussLegendre(x, w);
        rule &r = it->second;
        r.phi   = Vector(n);
        r.cos_weight.resize(n, M + 1);
        r.sin_weight.resize(n, M + 1);
        for (Index k = 0; k < n; k++) {
          r.phi[k] = pi_v<Numeric> * 0.5 * (x[k] + 1.0);
          for (Index m = 0; m <= M; m++) {
            const Numeric dm   = static_cast<Numeric>(m);
            r.cos_weight[k, m] = (m == 0 ? 0.5 : 1.0) * w[k] * std::cos(dm * r.phi[k]);
            r.sin_weight[k, m] = w[k] * std::sin(dm * r.phi[k]);
          }
        }
      }
      return it->second;
    };

#pragma omp for collapse(2) schedule(dynamic)
    for (Size ii = 0; ii < za_inc_grid->size(); ii++) {
      for (Size is = 0; is < za_scat.size(); is++) {
        try {
          Numeric    scale = 0.0;
          const auto modes = [&](Index n, std::vector<Numeric> &c) {
            const rule &r = get_rule(n);
            stdr::fill(c, 0.0);
            for (Index k = 0; k < n; k++) {
              const Numeric delta = Conversion::rad2deg(r.phi[k]);
              const auto [theta, c_1, c_2, s_1, s_2] =
                  detail::rotation_coefficients(0.0, (*za_inc_grid)[ii], delta, za_scat[is]);
              local_matrix(theta, F);
              for (Index set = 0; set < nset; set++) {
                // Its 16 elements, row-major
                const auto Z = rtepack::rotated(F[set / nf, set % nf], c_1, s_1, c_2, s_2).data;
                for (Index e = 0; e < 16; e++) scale = std::max<Numeric>(scale, std::abs(Z[e]));
                Numeric *cs = c.data() + set * nmode;
                for (Index m = 0; m <= M; m++, cs += 32) {
                  const Numeric wc = r.cos_weight[k, m], ws = r.sin_weight[k, m];
                  for (Index e = 0; e < 16; e++) {
                    cs[e]      += wc * Z[e];
                    cs[16 + e] += ws * Z[e];
                  }
                }
              }
            }
            // The parts that vanish by the symmetry
            for (Index j = 0; j < nset * (M + 1); j++)
              for (Index e = 0; e < 16; e++) c[j * 32 + (odd[e] ? 0 : 16) + e] = 0.0;
          };

          Index n = start;
          modes(n, coarse);
          while (true) {
            n *= 2;
            modes(n, fine);
            Numeric change = 0.0;
            for (Size j = 0; j < fine.size(); j++) change = std::max(change, std::abs(fine[j] - coarse[j]));
            if (change <= tolerance * scale) break;
            ARTS_USER_ERROR_IF(2 * n > max_nodes,
                               "The azimuthal Fourier modes of the laboratory-frame phase matrix at za_inc = {} deg "
                               "and za_scat = {} deg did not converge with {} azimuths (they still change by {:.3e} "
                               "of the largest phase-matrix element); the scattering matrix is too sharp in the "
                               "scattering angle to resolve",
                               (*za_inc_grid)[ii],
                               za_scat[is],
                               2 * n,
                               change / scale)
            std::swap(coarse, fine);
          }

          for (Index i_t = 0; i_t < nt; i_t++)
            for (Index i_f = 0; i_f < nf; i_f++)
              for (Index m = 0; m <= M; m++)
                for (Index cs = 0; cs < 2; cs++)
                  for (Index e = 0; e < 16; e++)
                    out[i_t, i_f, ii, is, m, cs].data[e] = fine[(i_t * nf + i_f) * nmode + (m * 2 + cs) * 16 + e];
        } catch (const std::exception &e) {
#pragma omp critical(tro_lab_frame_fourier_modes)
          if (error.empty()) error = e.what();
        }
      }
    }
  }
  ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
  return out;
}

}  // namespace scattering
