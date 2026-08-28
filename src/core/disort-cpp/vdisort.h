#pragma once

#include <matpack.h>
#include <operators.h>
#include <rtepack_mueller_matrix.h>
#include <rtepack_spectral_matrix.h>
#include <rtepack_stokes_vector.h>

#include <array>
#include <functional>
#include <span>
#include <vector>

/**
 * Vector DISORT core.
 *
 * The class deliberately follows disort::main_data's naming and update model.
 * Every semantic difference from the scalar class is enclosed in a
 * `VDISORT CHANGE BEGIN/END` marker.  Ordinary bookkeeping which is identical
 * to DISORT is left unmarked.
 *
 * Conventions follow Lin et al. (2022), doi:10.3389/frsen.2022.880768:
 *   - Stokes order is [I, Q, U, V].
 *   - positive mu is upward and the positive quadrature points precede the
 *     corresponding negative points;
 *   - phase_matrix[alpha,m,l,i,j] is an rtepack Mueller block mapping incident
 *     stream j to outgoing stream i and already contains the combined matrices
 *     of Eq. (81);
 *   - alpha=0 is the combined cosine equation and alpha=1 the combined sine
 *     equation;
 *   - boundary arrays use the same combined-component convention.
 */
namespace vdisort {
inline constexpr Index stokes_dimension = 4;
inline constexpr Index cosine_mode      = 0;
inline constexpr Index sine_mode        = 1;

using phase_matrix_data      = rtepack::muelmat_tensor5;
using beam_phase_matrix_data = rtepack::muelmat_tensor4;

struct phase_matrix_fourier_coefficients {
  rtepack::muelmat_matrix cosine;
  rtepack::muelmat_matrix sine;
};

/** Split a complex scattering phase matrix into real Fourier coefficients.
 *
 * For C = A + iB, this uses A as the cosine coefficient and -B as the sine
 * coefficient.  No factor of two is applied.  Both outputs preserve the input
 * [frequency, spectral-coefficient] shape.
 */
[[nodiscard]] phase_matrix_fourier_coefficients phase_matrix_fourier_split(
    const rtepack::specmat_matrix_const_view& phase_matrix);

// VDISORT CHANGE BEGIN: scalar radiances become one Stokes vector per stream.
struct u_data {
  rtepack::stokvec_vector intensities;  // [NQuad]
};

/** Scratch and result storage for the formal solution at user directions. */
struct user_u_data {
  rtepack::stokvec_vector intensities;  // [NUser]
};

struct u0_data {
  rtepack::stokvec_vector u0;  // [NQuad]
};
// VDISORT CHANGE END

struct flux_data {
  rtepack::stokvec_vector u0;  // scratch [NQuad]
};

struct flux_values {
  Numeric up{};
  Numeric down_diffuse{};
  Numeric down_direct{};
  Numeric dfdt{};
};

/** Laboratory-frame Mueller phase function used by the delta-M correction.
 * Arguments are layer, outgoing (mu, phi), and incident (mu, phi).
 */
using lab_phase_function = std::function<rtepack::muelmat(Index, Numeric, Numeric, Numeric, Numeric)>;
using lab_phase_pair_convolution_function =
    std::function<rtepack::muelmat(Index, Index, Numeric, Numeric, Numeric, Numeric)>;

/** Cached, fully Mueller-valued Nakajima-Tanaka delta-M correction.
 *
 * `removed_phase` is the normalized angular shape R of the removed peak.  The
 * IMS cache evaluates 2R - (1/4pi) integral R(out,mid) R(mid,beam) dOmega_mid,
 * with all matrices expressed in the laboratory Stokes frame.  Consequently
 * reference-plane rotations belong in the supplied callbacks and are retained
 * through the Mueller product.
 */
class delta_m_correction_cache {
  AscendingGrid    physical_tau_;
  AscendingGrid    scaled_tau_;
  Vector           omega_;
  Vector           fraction_;
  Vector           scale_;
  Vector           user_mu_;
  Vector           phi_;
  Numeric          mu0_{};
  Numeric          phi0_{};
  rtepack::stokvec beam_{};

  rtepack::muelmat_tensor3 tms_operator_;  // [layer, phi, user]
  rtepack::muelmat_tensor3 ims_operator_;  // [boundary, phi, user]
  Vector                   ims_scalar_;    // [boundary]
  Vector                   ims_mu0_;       // [boundary]

 public:
  delta_m_correction_cache() = default;
  delta_m_correction_cache(AscendingGrid                       physical_tau,
                           Vector                              omega,
                           Vector                              fraction,
                           Numeric                             mu0,
                           Numeric                             phi0,
                           rtepack::stokvec                    beam,
                           Vector                              user_mu,
                           Vector                              phi,
                           lab_phase_function                  original_phase,
                           lab_phase_function                  transport_phase,
                           lab_phase_function                  removed_phase,
                           Index                               intermediate_mu          = 32,
                           Index                               intermediate_phi         = 64,
                           lab_phase_pair_convolution_function removed_pair_convolution = {});

  /** Return the cached TMS+IMS Stokes correction at one cached azimuth. */
  [[nodiscard]] rtepack::stokvec_vector evaluate(Numeric tau, Index phi_index) const;

  [[nodiscard]] const AscendingGrid& physical_tau() const { return physical_tau_; }
  [[nodiscard]] const AscendingGrid& scaled_tau() const { return scaled_tau_; }
  [[nodiscard]] const Vector&        user_mu() const { return user_mu_; }
  [[nodiscard]] const Vector&        phi() const { return phi_; }
};

// VDISORT CHANGE BEGIN: a BRDF Fourier mode is a 4x4 block operator.
/** Polarized BRDF Fourier mode in the combined representation.
 *
 * Each callback fills a Mueller-block matrix with shape [n_out, n_in].  Block
 * (i,j) maps the incident Stokes vector at mu_in[j] to the reflected Stokes
 * vector at mu_out[i].  Values are the reflection matrices R^m in Eq. (123),
 * before the pi*w*mu quadrature factor in Eq. (124).
 */
struct BDRF {
  using func_t = CustomOperator<void, rtepack::muelmat_matrix_view, const ConstVectorView&, const ConstVectorView&>;

  func_t cosine;
  func_t sine;

  void operator()(Index                        alpha,
                  rtepack::muelmat_matrix_view out,
                  const ConstVectorView&       mu_out,
                  const ConstVectorView&       mu_in) const;
};
// VDISORT CHANGE END

/** Convert numeric ordinary phase-matrix cosine/sine Fourier coefficients to the
 * combined VDISORT matrices of Lin et al. (2022), Eqs. (81)-(82).
 *
 * Inputs have shape [NFourier, NLayers, NQuad, NQuad, 4, 4].  The result has
 * shape [2, NFourier, NLayers, NQuad, NQuad, 4, 4].  This transform assumes
 * randomly oriented particles and therefore the block structure in Eq. (75).
 */
[[nodiscard]] Tensor7 combine_phase_matrices(const Tensor6& cosine, const Tensor6& sine);

/** Native rtepack overload of combine_phase_matrices.
 *
 * Inputs have shape [NFourier, NLayers, NQuad, NQuad] of Mueller blocks and
 * the result has shape [2, NFourier, NLayers, NQuad, NQuad].
 */
[[nodiscard]] phase_matrix_data combine_phase_matrices(const rtepack::muelmat_tensor4& cosine,
                                                       const rtepack::muelmat_tensor4& sine);

/** Beam-angle counterpart of combine_phase_matrices.
 *
 * Inputs have shape [NFourier, NLayers, NQuad, 4, 4], with the incident
 * direction fixed at -mu0.  The result has shape
 * [2, NFourier, NLayers, NQuad, 4, 4].
 */
[[nodiscard]] Tensor6 combine_beam_phase_matrices(const Tensor5& cosine, const Tensor5& sine);

/** Native rtepack overload of combine_beam_phase_matrices. */
[[nodiscard]] beam_phase_matrix_data combine_beam_phase_matrices(const rtepack::muelmat_tensor3& cosine,
                                                                 const rtepack::muelmat_tensor3& sine);

/** The main data structure for the polarized VDISORT algorithm.
 *
 * This is the direct full-eigenproblem formulation (paper Eq. 87), retaining
 * complex conjugate eigenpairs.  It intentionally favors a transparent port
 * of disort::main_data over the optional reduced eigenproblem of Eq. (96).
 */
class main_data {
  Index NLayers{0};
  Index NQuad{0};
  Index NFourier{0};
  Index N{0};
  Index Nscoeffs{0};
  Index NBDRF{0};
  Index NState{0};
  Index NHalfState{0};
  bool  has_source_poly{false};
  bool  has_beam_source{false};

  //! User inputs
  AscendingGrid tau_arr{};    // [NLayers]
  Vector        omega_arr{};  // [NLayers]

  // VDISORT CHANGE BEGIN: all radiation-bearing inputs carry Stokes and mode dimensions.
  rtepack::stokvec_matrix  source_poly_coeffs{};        // [NLayers, Nscoeffs]
  Vector                   source_coordinate_scale{};   // [NLayers], x = offset + scale*tau
  Vector                   source_coordinate_offset{};  // [NLayers]
  phase_matrix_data        phase_matrix{};              // [2, NFourier, NLayers, NQuad, NQuad] of 4x4 blocks
  rtepack::stokvec_tensor3 boundary_up{};               // [2, NFourier, N]
  rtepack::stokvec_tensor3 boundary_down{};             // [2, NFourier, N]
  std::vector<BDRF>        brdf_fourier_modes{};        // [NBDRF]
  Numeric                  mu0{};
  rtepack::stokvec         beam_stokes{};  // irradiance Stokes vector S_b
  Numeric                  phi0{};
  beam_phase_matrix_data   beam_phase_matrix{};  // [2, NFourier, NLayers, NQuad] of 4x4 blocks
  // VDISORT CHANGE END

  //! Derived values
  Vector mu_arr{};      // [NQuad]
  Vector inv_mu_arr{};  // [NQuad]
  Vector W{};           // [N]

  // VDISORT CHANGE BEGIN: the homogeneous solution is complex (paper Sec. 3.3.1).
  ComplexTensor5             G_collect{};       // [2, NFourier, NLayers, NState, NState]
  ComplexTensor4             K_collect{};       // [2, NFourier, NLayers, NState]
  ComplexTensor4             GC_collect{};      // [2, NFourier, NLayers, NState]
  rtepack::stokvec_tensor4   B_collect{};       // [2, NFourier, NLayers, NQuad]
  rtepack::stokvec_tensor5   source_collect{};  // [2, NFourier, NLayers, NQuad, Nscoeffs]
  rtepack::stokvec_tensor4   um{};              // [NLayers, 2, NFourier, NQuad], layer-bottom fields
  std::vector<unsigned char> top_anchored{};    // one flag per eigenmode collection
  // VDISORT CHANGE END

  [[nodiscard]] Index   state_index(Index stream, Index stokes) const;
  [[nodiscard]] Index   anchor_index(Index alpha, Index m, Index layer, Index eigen) const;
  [[nodiscard]] Numeric layer_top(Index layer) const;
  [[nodiscard]] Complex homogeneous(Index alpha, Index m, Index layer, Index state, Index eigen, Numeric tau) const;
  [[nodiscard]] Numeric particular(Index alpha, Index m, Index layer, Index state, Numeric tau) const;
  void                  combined_field(MatrixView out, Numeric tau) const;

 public:
  main_data()                            = default;
  main_data(const main_data&)            = default;
  main_data(main_data&&)                 = default;
  main_data& operator=(const main_data&) = default;
  main_data& operator=(main_data&&)      = default;

  main_data(Index NLayers, Index NQuad, Index NFourier, Index Nscoeffs, Index NBDRF);

  // VDISORT CHANGE BEGIN: compatibility constructor for numeric 4x4 trailing dimensions.
  main_data(Index             NQuad,
            Index             NFourier,
            AscendingGrid     tau_arr,
            Vector            omega_arr,
            Tensor7           phase_matrix,
            Tensor4           boundary_up,
            Tensor4           boundary_down,
            Tensor3           source_poly_coeffs,
            std::vector<BDRF> brdf_fourier_modes,
            Numeric           mu0,
            Vector            beam_stokes,
            Numeric           phi0,
            Tensor6           beam_phase_matrix        = {},
            Vector            source_coordinate_scale  = {},
            Vector            source_coordinate_offset = {});

  // Native constructor: phase operators are stored as rtepack Mueller blocks.
  main_data(Index                    NQuad,
            Index                    NFourier,
            AscendingGrid            tau_arr,
            Vector                   omega_arr,
            phase_matrix_data        phase_matrix,
            rtepack::stokvec_tensor3 boundary_up,
            rtepack::stokvec_tensor3 boundary_down,
            rtepack::stokvec_matrix  source_poly_coeffs,
            std::vector<BDRF>        brdf_fourier_modes,
            Numeric                  mu0,
            rtepack::stokvec         beam_stokes,
            Numeric                  phi0,
            beam_phase_matrix_data   beam_phase_matrix        = {},
            Vector                   source_coordinate_scale  = {},
            Vector                   source_coordinate_offset = {});
  // VDISORT CHANGE END

  [[nodiscard]] Index tau_index(Numeric tau) const;

  void u(u_data& data, Numeric tau, Numeric phi) const;

  /** Radiance at arbitrary nonzero polar-angle cosines.
   *
   * The directional phase inputs contain the combined VDISORT matrices at
   * the requested outgoing directions.  Their shapes are
   * [2, NFourier, NLayers, NUser, NQuad] and, for the direct beam,
   * [2, NFourier, NLayers, NUser].  Supplying these values explicitly avoids
   * attempting to reconstruct unsampled polarized phase matrices from the
   * quadrature grid.
   *
   * The discrete-ordinate field is formally integrated along each requested
   * ray.  This call is on demand and does not add work to u(), gridded_u(), or
   * the flux methods.
   */
  void u_user(user_u_data&                  data,
              Numeric                       tau,
              Numeric                       phi,
              const ConstVectorView&        user_mu,
              const phase_matrix_data&      user_phase_matrix,
              const beam_phase_matrix_data& user_beam_phase_matrix = {}) const;
  void u0(u0_data& data, Numeric tau) const;

  [[nodiscard]] Numeric                     flux_up(flux_data&, Numeric tau) const;
  [[nodiscard]] std::pair<Numeric, Numeric> flux_down(flux_data&, Numeric tau) const;
  [[nodiscard]] flux_values                 flux(flux_data&, Numeric tau) const;

  void gridded_u(Tensor4View out, const Vector& phi) const;
  void gridded_flux(VectorView up, VectorView down, VectorView down_direct) const;
  void ungridded_u(Tensor4View out, const AscendingGrid& tau, const Vector& phi) const;
  void ungridded_flux(VectorView up, VectorView down, VectorView down_direct, const AscendingGrid& tau) const;

  [[nodiscard]] rtepack::stokvec_tensor3_const_view layer_um(Size layer) const;

  void diagonalize();
  void transmission();
  void source_function();
  void solve_for_coefs();
  void rad_field();

  void check_input_size() const;
  void check_input_value() const;
  void update_all();

  [[nodiscard]] bool has_complex_eigensolutions(Numeric tolerance = 1e-12) const;

  [[nodiscard]] const Vector&        mu() const { return mu_arr; }
  [[nodiscard]] const Vector&        weights() const { return W; }
  [[nodiscard]] const AscendingGrid& tau() const { return tau_arr; }
  [[nodiscard]] const Vector&        omega() const { return omega_arr; }

  // VDISORT CHANGE BEGIN: polarized accessors replace scalar coefficient accessors.
  [[nodiscard]] const phase_matrix_data&        all_phase_matrices() const { return phase_matrix; }
  [[nodiscard]] rtepack::muelmat_tensor5_view   all_phase_matrices() { return phase_matrix; }
  [[nodiscard]] const rtepack::stokvec_tensor3& upward_boundary() const { return boundary_up; }
  [[nodiscard]] rtepack::stokvec_tensor3_view   upward_boundary() { return boundary_up; }
  [[nodiscard]] const rtepack::stokvec_tensor3& downward_boundary() const { return boundary_down; }
  [[nodiscard]] rtepack::stokvec_tensor3_view   downward_boundary() { return boundary_down; }
  [[nodiscard]] const rtepack::stokvec_matrix&  source_poly() const { return source_poly_coeffs; }
  [[nodiscard]] rtepack::stokvec_matrix_view    source_poly() { return source_poly_coeffs; }
  [[nodiscard]] std::span<const BDRF>           brdf_modes() const { return brdf_fourier_modes; }
  [[nodiscard]] std::span<BDRF>                 brdf_modes() { return brdf_fourier_modes; }
  [[nodiscard]] const rtepack::stokvec&         beam_source() const { return beam_stokes; }
  [[nodiscard]] rtepack::stokvec&               beam_source() { return beam_stokes; }
  void                                          set_beam_source(rtepack::stokvec beam);
  // VDISORT CHANGE END

  [[nodiscard]] Numeric  solar_zenith() const { return mu0; }
  [[nodiscard]] Numeric& solar_zenith() { return mu0; }
  [[nodiscard]] Numeric  beam_azimuth() const { return phi0; }
  [[nodiscard]] Numeric& beam_azimuth() { return phi0; }
};
}  // namespace vdisort
