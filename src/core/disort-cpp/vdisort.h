#pragma once

#include <debug.h>
#include <matpack.h>
#include <operators.h>
#include <rtepack_mueller_matrix.h>
#include <rtepack_spectral_matrix.h>
#include <rtepack_stokes_vector.h>

#include <array>
#include <functional>
#include <span>
#include <vector>

#include "common.h"

/**
 * Vector DISORT core.
 *
 * The class follows disort::main_data's naming and update model so equivalent
 * scalar and polarized solver stages can be compared directly.
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

struct u_data {
  rtepack::stokvec_vector intensities;  // [NQuad]
};

/** Result storage for the formal solution at user directions. */
struct user_u_data {
  rtepack::stokvec_vector intensities;  // [NUser]
};

struct u0_data {
  rtepack::stokvec_vector u0;  // [NQuad]
};

struct flux_data {
  rtepack::stokvec_vector u0;  // scratch [NQuad]
};

using flux_values    = disort_common::flux_values;
using ims_convention = disort_common::ims_convention;

/** Laboratory-frame Mueller phase function used by the delta-M correction.
 * Arguments are layer, outgoing (mu, phi), and incident (mu, phi).
 */
using lab_phase_function = std::function<rtepack::muelmat(Index, Numeric, Numeric, Numeric, Numeric)>;
using lab_phase_pair_convolution_function =
    std::function<rtepack::muelmat(Index, Index, Numeric, Numeric, Numeric, Numeric)>;

/** Cached, fully Mueller-valued Nakajima-Tanaka delta-M correction.
 *
 * TMS replaces the single scattering of the beam that the transport solver
 * carries by that of the original phase matrix.  The cache holds the
 * original part, from `original_phase`.  The transport part is what the
 * solver used: the combined beam phase matrices into the correction's
 * directions, synthesised from their NFourier modes at each azimuth, which
 * TMS(), gridded_TMS() and evaluate() take as `transport_beam_phase`
 * (main_data passes its own, or those given to u_user()).  For F11 this is
 * the whole transport matrix, since its series has azimuthal modes
 * m < NFourier only; the polarized elements also contain the rotation
 * angles, which are not band-limited, so subtracting the full transport
 * matrix would add their azimuthal truncation back as an error.
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

  rtepack::muelmat_tensor3 tms_operator_;  // [layer, phi, user], the original part of TMS
  rtepack::muelmat_tensor3 ims_operator_;  // [boundary, phi, user]
  Vector                   ims_scalar_;    // [boundary]
  Vector                   ims_mu0_;       // [boundary]

  //! The TMS source of a layer into a user direction at a cached azimuth: original minus transport part
  [[nodiscard]] rtepack::stokvec tms_source(const beam_phase_matrix_data& transport_beam_phase,
                                            Index                         layer,
                                            Index                         phi_index,
                                            Index                         user) const;

  //! Check that the transport beam phase matrices are [2, NFourier, NLayers, NUser] or empty
  void check_transport_beam_phase(const beam_phase_matrix_data& transport_beam_phase) const;

 public:
  /** Construct an empty correction cache. */
  delta_m_correction_cache() = default;

  /** Precompute polarized TMS and IMS operators on the requested user-angle grid. */
  delta_m_correction_cache(AscendingGrid                       physical_tau,
                           Vector                              omega,
                           Vector                              fraction,
                           Numeric                             mu0,
                           Numeric                             phi0,
                           rtepack::stokvec                    beam,
                           Vector                              user_mu,
                           Vector                              phi,
                           lab_phase_function                  original_phase,
                           lab_phase_function                  removed_phase,
                           Index                               intermediate_mu          = 32,
                           Index                               intermediate_phi         = 64,
                           lab_phase_pair_convolution_function removed_pair_convolution = {});

  /** Evaluate the cached Mueller-valued TMS correction at one physical depth and cached azimuth.
   *
   * transport_beam_phase is [2, NFourier, NLayers, NUser]: the solver's combined beam phase matrices into
   * the correction's directions.  Empty means that the solver carries no beam single scattering.
   */
  void TMS(rtepack::stokvec_vector&      tms,
           Numeric                       tau,
           Index                         phi_index,
           const beam_phase_matrix_data& transport_beam_phase) const;

  /** Evaluate the cached Mueller-valued IMS correction at one physical depth and cached azimuth. */
  void IMS(rtepack::stokvec_vector& ims,
           Numeric                  tau,
           Index                    phi_index,
           ims_convention           convention = ims_convention::disort) const;

  /** Evaluate TMS at every physical layer bottom, cached azimuth, and user direction (see TMS()). */
  void gridded_TMS(rtepack::stokvec_tensor3_view tms, const beam_phase_matrix_data& transport_beam_phase) const;

  /** Evaluate IMS at every physical layer bottom, cached azimuth, and user direction. */
  void gridded_IMS(rtepack::stokvec_tensor3_view ims, ims_convention convention = ims_convention::disort) const;

  /** Return the cached TMS+IMS Stokes correction at one cached azimuth (see TMS()). */
  [[nodiscard]] rtepack::stokvec_vector evaluate(Numeric                       tau,
                                                 Index                         phi_index,
                                                 const beam_phase_matrix_data& transport_beam_phase,
                                                 ims_convention convention = ims_convention::disort) const;

  /** Map one physical optical depth to the transport solver's scaled coordinate. */
  [[nodiscard]] Numeric scaled_tau(Numeric physical_tau) const;

  /** Return the physical, unscaled layer-bottom optical depths. */
  [[nodiscard]] const AscendingGrid& physical_tau() const { return physical_tau_; }

  /** Return the delta-M-scaled layer-bottom optical depths. */
  [[nodiscard]] const AscendingGrid& scaled_tau() const { return scaled_tau_; }

  /** Return the signed user-direction cosines cached by the correction. */
  [[nodiscard]] const Vector& user_mu() const { return user_mu_; }

  /** Return the user azimuths cached by the correction. */
  [[nodiscard]] const Vector& phi() const { return phi_; }
};

/** A surface's reflection into upward directions, in VDISORT's combined representation.
 *
 * brdf [2, NBRDF, n_out, N]: the combined cosine (alpha = 0) and sine
 *   (alpha = 1) Fourier modes R^{alpha,m}(mu_o, mu_j) of the polarized BRDF,
 *   m < NBRDF, from the downward streams of cosine mu_j into the upward
 *   directions of cosine mu_o; Mueller block (o, j) maps the incident Stokes
 *   vector to the reflected one.  These are the reflection matrices R^m of
 *   Lin et al. (2022), Eq. (123), before the pi w_j mu_j quadrature factor of
 *   Eq. (124), which the solver applies.  Higher modes reflect nothing.
 * beam_brdf [2, NBRDF, n_out]: the same from the direct beam, at the beam
 *   cosine beam_mu; beam_mu is NaN when it was computed without a beam, and a
 *   solver with a beam refuses a surface whose beam_mu is not its own.
 * specular [n_out], or empty: the mirror part R(mu_o).  Radiance propagating
 *   down at the cosine mu is reflected into the upward direction of the same
 *   cosine and azimuth, I_up(mu) = R(mu) I_down(mu), in every Fourier mode and
 *   both combined systems.  A delta distribution in direction has no kernel
 *   on the streams, so the solver applies it itself: R(mu_i) between equal
 *   streams, and at a user direction to the downward radiance in the partner
 *   direction.  It cannot reflect the direct beam (its mirror image is an
 *   upward beam), so a solver with a beam refuses it, and it is never for a
 *   beam (beam_mu is NaN).
 *
 * The shapes, and whether there is a specular part, are those it is
 * constructed with.  Its data is set in place, through the views or set(), so
 * nothing is reallocated; beam_mu is set only with the beam reflection it
 * belongs to.  A surface without modes or specular part reflects nothing,
 * whatever its directions.
 *
 * On the solver's upward streams (n_out = N) it is main_data's surface, which
 * the settings set for each frequency.  At user directions it is the
 * reflection u_user() needs there: like() makes one for the modes, streams,
 * beam and specular part of the solver's surface, or brdf:: evaluates a
 * surface model there.
 */
class surface_reflection {
  rtepack::muelmat_tensor4 kernel{};       // [2, NBRDF, n_out, N]
  rtepack::muelmat_tensor3 beam_kernel{};  // [2, NBRDF, n_out]
  rtepack::muelmat_vector  mirror{};       // [n_out] or [0]
  Numeric                  beam_cosine{std::numeric_limits<Numeric>::quiet_NaN()};

 public:
  /** No reflection */
  surface_reflection() = default;

  /** Zero reflection of NBRDF modes from N streams into n_out directions, for the beam cosine beam_mu (NaN: none),
   * with a (zero) specular part if specular */
  surface_reflection(Index   NBRDF,
                     Index   n_out,
                     Index   N,
                     Numeric beam_mu  = std::numeric_limits<Numeric>::quiet_NaN(),
                     bool    specular = false);

  /** The reflection brdf [2, NBRDF, n_out, N] and beam_brdf [2, NBRDF, n_out], for the beam cosine beam_mu (NaN:
   * none), and the specular part [n_out] or none */
  surface_reflection(rtepack::muelmat_tensor4 brdf,
                     rtepack::muelmat_tensor3 beam_brdf,
                     Numeric                  beam_mu  = std::numeric_limits<Numeric>::quiet_NaN(),
                     rtepack::muelmat_vector  specular = {});

  /** Zero reflection of this one's modes from its streams into n_out directions, for its beam cosine, with a (zero)
   * specular part if it has one */
  [[nodiscard]] surface_reflection like(Index n_out) const;

  /** The number of Fourier modes of the BRDF, NBRDF */
  [[nodiscard]] Index modes() const { return kernel.extent(1); }

  /** The number of upward directions reflected into, n_out */
  [[nodiscard]] Index directions() const { return kernel.extent(2); }

  /** The number of downward streams reflected from, N */
  [[nodiscard]] Index streams() const { return kernel.extent(3); }

  /** Whether there is a specular part */
  [[nodiscard]] bool has_specular() const { return not mirror.empty(); }

  /** The beam cosine of the beam reflection, NaN if computed without a beam */
  [[nodiscard]] Numeric beam_mu() const { return beam_cosine; }

  //! [2, NBRDF, n_out, N]
  [[nodiscard]] rtepack::muelmat_tensor4_const_view brdf() const { return kernel; }
  [[nodiscard]] rtepack::muelmat_tensor4_view       brdf() { return kernel; }

  //! [2, NBRDF, n_out]
  [[nodiscard]] rtepack::muelmat_tensor3_const_view beam_brdf() const { return beam_kernel; }
  [[nodiscard]] rtepack::muelmat_tensor3_view       beam_brdf() { return beam_kernel; }

  //! [n_out] or [0]
  [[nodiscard]] rtepack::muelmat_vector_const_view specular() const { return mirror; }
  [[nodiscard]] rtepack::muelmat_vector_view       specular() { return mirror; }

  //! The specular part into direction i
  [[nodiscard]] const rtepack::muelmat& specular(Index i) const { return mirror[i]; }

  /** Set the reflection in place, for the beam cosine beam_mu (NaN: none); the shapes must be this one's */
  void set(const rtepack::muelmat_tensor4_const_view& brdf,
           const rtepack::muelmat_tensor3_const_view& beam_brdf,
           Numeric                                    beam_mu,
           const rtepack::muelmat_vector_const_view&  specular);

  /** Set the reflection in place from one of the same shapes */
  void set(const surface_reflection& other);
};

/** Convert ordinary phase-matrix cosine/sine Fourier coefficients to the
 * combined VDISORT matrices of Lin et al. (2022), Eqs. (81)-(82).
 *
 * Inputs have shape [NFourier, NLayers, NQuad, NQuad] of Mueller blocks and
 * the result has shape [2, NFourier, NLayers, NQuad, NQuad].
 */
[[nodiscard]] phase_matrix_data combine_phase_matrices(const rtepack::muelmat_tensor4& cosine,
                                                       const rtepack::muelmat_tensor4& sine);

/** Beam-angle counterpart of combine_phase_matrices.
 *
 * Inputs are the ordinary cosine/sine coefficients C^m(mu_i, -mu0) and
 * S^m(mu_i, -mu0) with shape [NFourier, NLayers, NQuad] of Mueller blocks.
 * The beam has only cosine azimuth terms, so this is not Eq. (81): the
 * cosine system takes rows I, Q of C^m and rows U, V of S^m, and the sine
 * system rows I, Q of S^m and rows U, V of C^m.
 */
[[nodiscard]] beam_phase_matrix_data combine_beam_phase_matrices(const rtepack::muelmat_tensor3& cosine,
                                                                 const rtepack::muelmat_tensor3& sine);

/** Solver-ready result of an explicitly specified polarized delta-M transform. */
struct delta_m_transport_data {
  AscendingGrid          tau;
  Vector                 omega;
  phase_matrix_data      phase_matrix;
  beam_phase_matrix_data beam_phase_matrix;
  Vector                 source_coordinate_scale;
  Vector                 source_coordinate_offset;
};

/** Apply a caller-defined polarized delta-M split to VDISORT inputs.
 *
 * `removed_phase_matrix` and `removed_beam_phase_matrix` contain the
 * normalized Mueller forward peak R in the same combined representation and
 * on the same angular grids as their original counterparts.  This function
 * applies
 *
 *   dt' = (1 - omega*f) dt,
 *   omega' = omega*(1-f)/(1-omega*f),
 *   P' = (P-f*R)/(1-f),
 *
 * and returns the affine coordinate map needed to evaluate physical source
 * polynomials on the scaled depth grid.  It validates sizes, ranges, and
 * finite values, but deliberately does not infer or certify the physical
 * normalization of the caller's removed Mueller peak.
 */
[[nodiscard]] delta_m_transport_data delta_m_preprocess(const AscendingGrid&          physical_tau,
                                                        const Vector&                 physical_omega,
                                                        const Vector&                 fraction,
                                                        const phase_matrix_data&      original_phase_matrix,
                                                        const phase_matrix_data&      removed_phase_matrix,
                                                        const beam_phase_matrix_data& original_beam_phase_matrix = {},
                                                        const beam_phase_matrix_data& removed_beam_phase_matrix  = {});

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

  // Radiation-bearing inputs carry Stokes and combined-mode dimensions.
  rtepack::stokvec_matrix  source_poly_coeffs{};        // [NLayers, Nscoeffs], Stokes source function B
  Vector                   source_coordinate_scale{};   // [NLayers], x = offset + scale*tau
  Vector                   source_coordinate_offset{};  // [NLayers]
  phase_matrix_data        phase_matrix{};              // [2, NFourier, NLayers, NQuad, NQuad] of 4x4 blocks
  rtepack::stokvec_tensor3 boundary_up{};               // [2, NFourier, N], beneath the surface, transmitted with 1 - R
  rtepack::stokvec_tensor3 boundary_down{};             // [2, NFourier, N]
  surface_reflection       surface_refl{};              // [2, NBDRF, N, N], [2, NBDRF, N] and [N] or [0]
  Numeric                  mu0{};
  rtepack::stokvec         beam_stokes{};  // irradiance Stokes vector S_b
  Numeric                  phi0{};
  beam_phase_matrix_data   beam_phase_matrix{};  // [2, NFourier, NLayers, NQuad] of 4x4 blocks

  //! Derived values
  Vector                  mu_arr{};                          // [NQuad]
  Vector                  inv_mu_arr{};                      // [NQuad]
  Vector                  W{};                               // [N]
  Vector                  half_range_barycentric_weights{};  // [N]
  rtepack::stokvec_matrix scaled_source_poly_coeffs{};       // [NLayers, Nscoeffs], transport emission

  /* Scratch of update_all(), sized by the constructor and reused by every
     update.  What is accumulated or written only conditionally is reset
     where it is used; everything else is overwritten. */
  Matrix               transport{};              // [NState, NState], A of one Fourier mode and layer
  Matrix               transport_work{};         // [NState, NState], a copy for the destructive solves
  Matrix               real_eigenvectors{};      // [NState, NState], a conjugate pair as Re and Im columns
  Vector               eigenvalues_real{};       // [NState]
  Vector               eigenvalues_imag{};       // [NState]
  ComplexMatrix        eigenvectors{};           // [NState, NState]
  ComplexVector        eigenvalues{};            // [NState]
  std::vector<Index>   eigen_order{};            // [NState]
  std::vector<Index>   eigen_partner{};          // [NState], the conjugate partner before sorting, or -1
  std::vector<Index>   eigen_position{};         // [NState], the sorted position of each eigenpair
  Matrix               half_sum{};               // [NHalfState, NHalfState], alpha + beta Delta
  Matrix               half_difference{};        // [NHalfState, NHalfState], alpha - beta Delta
  Matrix               half_product{};           // [NHalfState, NHalfState], their product
  Matrix               half_eigenvectors{};      // [NHalfState, NHalfState], a conjugate pair as Re and Im columns
  Matrix               half_images{};            // [NHalfState, NHalfState], half_sum times the eigenvectors
  Vector               half_eigenvalues_real{};  // [NHalfState]
  Vector               half_eigenvalues_imag{};  // [NHalfState]
  diagonalize_workdata half_eigen_work{};        // [NHalfState]
  diagonalize_workdata eigen_work{};             // [NState]
  solve_workdata       solve_work{};             // [NState]
  Vector               state_rhs{};              // [NState], right-hand sides of the solves
  Vector               state_next{};             // [NState], the source polynomial's next coefficient
  matpack::band_matrix boundary_lhs{};           // [NLayers NState, NLayers NState], banded and real
  Vector               boundary_rhs{};           // [NLayers NState]
  ComplexTensor3       mode_exponential{};       // [NLayers, 2, NState], at each layer's top and bottom
  Matrix               surface_operator{};       // [NHalfState, NHalfState], R of one mode with its quadrature weights
  Vector               direct_reflection{};      // [NHalfState]
  Matrix               layer_field{};            // [2 NFourier, NState]

  // The homogeneous solution is complex (paper Sec. 3.3.1).
  ComplexTensor5                    G_collect{};                // [2, NFourier, NLayers, NState, NState]
  ComplexTensor4                    K_collect{};                // [2, NFourier, NLayers, NState]
  ComplexTensor4                    GC_collect{};               // [2, NFourier, NLayers, NState]
  rtepack::stokvec_tensor4          B_collect{};                // [2, NFourier, NLayers, NQuad]
  rtepack::stokvec_tensor5          source_collect{};           // [2, NFourier, NLayers, NQuad, Nscoeffs]
  rtepack::stokvec_tensor4          um{};                       // [NLayers, 2, NFourier, NQuad], layer-bottom fields
  std::vector<unsigned char>        top_anchored{};             // one flag per eigenmode collection
  std::vector<Index>                conjugate_partner{};        // per eigenmode collection: its conjugate, or -1
  std::vector<std::array<Index, 2>> conservative_pair_index{};  // [NLayers], stable m=0 cosine pair, or {-1,-1}
  Vector                            conservative_pair_kappa{};  // [NLayers]

  /** Flatten a quadrature-stream and Stokes-component pair into a state index. */
  [[nodiscard]] Index state_index(Index stream, Index stokes) const;
  /** Flatten the mode, layer, and eigenvector indices used by the anchoring flags. */
  [[nodiscard]] Index anchor_index(Index alpha, Index m, Index layer, Index eigen) const;
  /** Return the physical optical depth at the top of a layer. */
  [[nodiscard]] Numeric layer_top(Index layer) const;
  /** Evaluate one anchored complex homogeneous eigenmode component. */
  [[nodiscard]] Complex homogeneous(Index alpha, Index m, Index layer, Index state, Index eigen, Numeric tau) const;
  /** Evaluate the beam and polynomial particular solution for one state. */
  [[nodiscard]] Numeric particular(Index alpha, Index m, Index layer, Index state, Numeric tau) const;
  /** Assemble all combined Fourier and Stokes components on quadrature streams. */
  void combined_field(MatrixView out, Numeric tau) const;
  /** Formally integrate the cached modal solution along arbitrary user directions. */
  void user_fourier_modes(ComplexTensor4&                 out,
                          const AscendingGrid&            tau,
                          const ConstVectorView&          user_mu,
                          const phase_matrix_data&        user_phase_matrix,
                          const beam_phase_matrix_data&   user_beam_phase_matrix,
                          const rtepack::stokvec_tensor3& user_boundary,
                          const surface_reflection&       user_surface) const;
  /** Whether a layer's combined phase matrix of one Fourier mode is mirror symmetric about the horizontal
   * plane, P(-mu, -mu') = D P(mu, mu') D and P(-mu, mu') = D P(mu, -mu') D with D = diag(1, 1, -1, -1),
   * to rounding. */
  [[nodiscard]] bool mirror_symmetric(Index alpha, Index m, Index layer) const;
  /** Diagonalize the mirror-symmetric transport operator through its half-size product problem.
   * Returns false, leaving the eigensolution to the full problem, if a squared eigenvalue is zero. */
  [[nodiscard]] bool half_size_eigensolution();
  /** Validate that a correction cache uses this solver's scaled depth grid and optionally its streams. */
  void check_correction_compatibility(const delta_m_correction_cache& correction, bool quadrature_angles) const;

 public:
  /** Construct an empty, uninitialized solver. */
  main_data() = default;
  /** Copy a fully initialized solver and all of its cached solution data. */
  main_data(const main_data&) = default;
  /** Move a solver and its cached solution data. */
  main_data(main_data&&) = default;
  /** Copy-assign a solver and all of its cached solution data. */
  main_data& operator=(const main_data&) = default;
  /** Move-assign a solver and its cached solution data. */
  main_data& operator=(main_data&&) = default;

  /** Allocate solver input and work arrays without computing a solution; the surface has NBDRF modes and, if
   * specular, a specular part. */
  main_data(Index NLayers, Index NQuad, Index NFourier, Index Nscoeffs, Index NBDRF, bool specular = false);

  /** Construct and solve a polarized discrete-ordinate problem.
   *
   * Phase, boundary, beam, and source arrays use Mueller or Stokes blocks in
   * the combined cosine/sine representation documented for this class.
   */
  main_data(Index                     NQuad,
            Index                     NFourier,
            AscendingGrid             tau_arr,
            Vector                    omega_arr,
            phase_matrix_data         phase_matrix,
            rtepack::stokvec_tensor3  boundary_up,
            rtepack::stokvec_tensor3  boundary_down,
            rtepack::stokvec_matrix   source_poly_coeffs,
            const surface_reflection& surface,
            Numeric                   mu0,
            rtepack::stokvec          beam_stokes,
            Numeric                   phi0,
            beam_phase_matrix_data    beam_phase_matrix        = {},
            Vector                    source_coordinate_scale  = {},
            Vector                    source_coordinate_offset = {});

  /** Return the layer containing a physical optical depth. */
  [[nodiscard]] Index tau_index(Numeric tau) const;

  /** Evaluate polarized TMS from a compatible physical-phase correction cache on the quadrature streams. */
  void TMS(rtepack::stokvec_vector&        tms,
           Numeric                         physical_tau,
           Index                           phi_index,
           const delta_m_correction_cache& correction) const;

  /** Evaluate polarized IMS from a compatible physical-phase correction cache. */
  void IMS(rtepack::stokvec_vector&        ims,
           Numeric                         physical_tau,
           Index                           phi_index,
           const delta_m_correction_cache& correction,
           ims_convention                  convention = ims_convention::disort) const;

  /** Evaluate polarized TMS at every layer bottom and cached azimuth on the quadrature streams. */
  void gridded_TMS(rtepack::stokvec_tensor3_view tms, const delta_m_correction_cache& correction) const;

  /** Evaluate polarized IMS at every layer bottom and cached correction angle. */
  void gridded_IMS(rtepack::stokvec_tensor3_view   ims,
                   const delta_m_correction_cache& correction,
                   ims_convention                  convention = ims_convention::disort) const;

  /** Evaluate the Stokes radiance on all quadrature streams at one point. */
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
   *
   * Over a reflecting surface the upward user directions need the surface's
   * reflection into them, user_surface [2, NBDRF, NUser, N] (with the beam
   * reflection for the solver's beam, and the specular part [NUser] if the
   * surface has one); only the upward rows are read, so it may be evaluated
   * at |mu| for all.  Otherwise this throws.  Over a surface with a specular
   * part (e.g. brdf::fresnel_fourier_modes), the upward radiance at a user
   * cosine mu reflects the downward radiance at -mu at the surface, so every
   * upward user cosine must also come with -mu (and its phase matrices).
   *
   * user_boundary [2, NFourier, NUser], optional: the boundary radiance at
   * the start of each user ray, in the combined cosine/sine Fourier modes of
   * boundary_up and boundary_down, i.e. the radiance beneath the surface for
   * an upward direction (which the surface transmits with 1 - R, as
   * boundary_up) and the radiance incident at the top for a downward one.
   * Without it these are interpolated from the streams' boundary values
   * (barycentric in mu), which is exact only for boundaries that are
   * polynomials of low degree in mu.
   */
  void u_user(user_u_data&                    data,
              Numeric                         tau,
              Numeric                         phi,
              const ConstVectorView&          user_mu,
              const phase_matrix_data&        user_phase_matrix,
              const beam_phase_matrix_data&   user_beam_phase_matrix = {},
              const rtepack::stokvec_tensor3& user_boundary          = {},
              const surface_reflection&       user_surface           = {}) const;

  /** Evaluate and apply cached polarized IMS/TMS corrections at arbitrary user directions. */
  void u_user_corr(user_u_data&                    data,
                   rtepack::stokvec_vector&        ims,
                   rtepack::stokvec_vector&        tms,
                   Numeric                         physical_tau,
                   Index                           phi_index,
                   const delta_m_correction_cache& correction,
                   const phase_matrix_data&        user_phase_matrix,
                   const beam_phase_matrix_data&   user_beam_phase_matrix = {},
                   const surface_reflection&       user_surface           = {},
                   ims_convention                  convention             = ims_convention::disort) const;

  /** Radiance at arbitrary polar angles, optical depths, and azimuths.
   *
   * The output contains Stokes blocks with shape [NTau, NPhi, NUser].
   * User-direction phase projections are formed once per call and the
   * layerwise modal exponentials are integrated analytically.  This is the
   * bulk, test-facing counterpart of DISORT's TERPEV/TERPSO/USRINT path.
   * Over a reflecting surface, upward user cosines need the surface's
   * reflection into them, and over a specular one their downward partners;
   * user_boundary is optional; all as for u_user().
   */
  void ungridded_u_user(rtepack::stokvec_tensor3_view   out,
                        const AscendingGrid&            tau,
                        const Vector&                   phi,
                        const ConstVectorView&          user_mu,
                        const phase_matrix_data&        user_phase_matrix,
                        const beam_phase_matrix_data&   user_beam_phase_matrix = {},
                        const rtepack::stokvec_tensor3& user_boundary          = {},
                        const surface_reflection&       user_surface           = {}) const;

  /** Evaluate the azimuth-independent quadrature-stream field at one depth. */
  void u0(u0_data& data, Numeric tau) const;

  /** Evaluate and apply cached polarized IMS/TMS corrections on the quadrature streams. */
  void u_corr(u_data&                         data,
              rtepack::stokvec_vector&        ims,
              rtepack::stokvec_vector&        tms,
              Numeric                         physical_tau,
              Index                           phi_index,
              const delta_m_correction_cache& correction,
              ims_convention                  convention = ims_convention::disort) const;

  /** Return all flux components and the flux divergence at one optical depth. */
  [[nodiscard]] flux_values flux(flux_data&, Numeric tau) const;

  /** Evaluate radiances at every layer bottom for all quadrature streams and azimuths. */
  void gridded_u(Tensor4View out, const Vector& phi) const;
  /** Evaluate corrected radiances at all layer bottoms, cached azimuths, and quadrature streams. */
  void gridded_u_corr(Tensor4View                     out,
                      rtepack::stokvec_tensor3_view   tms,
                      rtepack::stokvec_tensor3_view   ims,
                      const delta_m_correction_cache& correction,
                      ims_convention                  convention = ims_convention::disort) const;
  /** Evaluate every flux quantity at every layer bottom. */
  void gridded_flux(VectorView up, VectorView down_diffuse, VectorView down_direct, VectorView dfdt) const;
  /** Evaluate quadrature-stream radiances on an arbitrary ascending optical-depth grid. */
  void ungridded_u(Tensor4View out, const AscendingGrid& tau, const Vector& phi) const;
  /** Evaluate every flux quantity on an arbitrary ascending optical-depth grid. */
  void ungridded_flux(
      VectorView up, VectorView down_diffuse, VectorView down_direct, VectorView dfdt, const AscendingGrid& tau) const;

  /** Return the cached combined Fourier field at the bottom of one layer. */
  [[nodiscard]] rtepack::stokvec_tensor3_const_view layer_um(Size layer) const;

  /** Assemble and diagonalize each layer and Fourier-mode transport operator.
   *
   * The operator is diagonalized in real arithmetic, and each conjugate pair
   * is expanded to two complex columns.  Without scattering the eigenvectors
   * are the states.  A mirror-symmetric operator is diagonalized through its
   * half-size product problem, except in the cosine system of order zero,
   * which holds the energy-conservation pair.  An exactly conservative
   * operator is represented by an explicitly constructed null vector and
   * generalized Jordan vector.
   */
  void diagonalize();
  /** Derive source-coordinate scales and internally scaled source coefficients. */
  void set_scales();
  /** Prepare transmission data; complex VDISORT modes are evaluated from anchored exponentials. */
  void transmission();
  /** Compute polynomial particular solutions for internal thermal sources. */
  void source_function();
  /** Solve the banded multilayer boundary and continuity system, in real arithmetic: a conjugate pair of
   * modes contributes the real and imaginary parts of one of them. */
  void solve_for_coefs();
  /** Cache the combined Fourier fields at every layer bottom. */
  void rad_field();

  /** Validate all input-array dimensions against the configured solver sizes. */
  void check_input_size() const;
  /** Validate physical ranges and angular conventions of solver inputs. */
  void check_input_value() const;
  /** Recompute every derived quantity and solution cache after input changes. */
  void update_all();

  /** Report whether any retained eigenvalue has a significant imaginary component. */
  [[nodiscard]] bool has_complex_eigensolutions(Numeric tolerance = 1e-12) const;

  /** Return signed quadrature cosines, with upward streams first. */
  [[nodiscard]] const Vector& mu() const { return mu_arr; }
  /** Return positive-hemisphere double-Gauss weights. */
  [[nodiscard]] const Vector& weights() const { return W; }
  /** Return physical layer-bottom optical depths. */
  [[nodiscard]] const AscendingGrid& tau() const { return tau_arr; }
  /** Replace the physical layer-bottom optical depths, keeping their number. */
  template <matpack::exact_md<Numeric, 1> T> void tau(T&& x) {
    if (static_cast<Index>(x.size()) != NLayers) arts::user_error("Invalid size for tau: {}", x.size());
    tau_arr = std::forward<T>(x);
  }
  /** Return physical single-scattering albedos. */
  [[nodiscard]] const Vector& omega() const { return omega_arr; }
  /** Return writable physical single-scattering albedos. */
  [[nodiscard]] VectorView omega() { return omega_arr; }

  /** Return a read-only view of all combined phase-matrix Fourier coefficients. */
  [[nodiscard]] const phase_matrix_data& all_phase_matrices() const { return phase_matrix; }
  /** Return a writable view of all combined phase-matrix Fourier coefficients. */
  [[nodiscard]] rtepack::muelmat_tensor5_view all_phase_matrices() { return phase_matrix; }
  /** Return a read-only view of all combined beam phase-matrix Fourier coefficients. */
  [[nodiscard]] const beam_phase_matrix_data& all_beam_phase_matrices() const { return beam_phase_matrix; }
  /** Return a writable view of all combined beam phase-matrix Fourier coefficients. */
  [[nodiscard]] rtepack::muelmat_tensor4_view all_beam_phase_matrices() { return beam_phase_matrix; }
  /** Return the radiance beneath the surface on the upward streams, combined Fourier coefficients.
   *
   * The surface transmits it with 1 - R, R the surface's reflection on the
   * streams (BRDF modes and specular part).  For a Planck radiance B(T) in
   * the I component of the cosine system of mode 0 this is Kirchhoff's
   * emission B(T) (e_I - R e_I); without a surface it is the radiance itself.
   */
  [[nodiscard]] const rtepack::stokvec_tensor3& upward_boundary() const { return boundary_up; }
  /** Return the writable radiance beneath the surface, see the const overload. */
  [[nodiscard]] rtepack::stokvec_tensor3_view upward_boundary() { return boundary_up; }
  /** Return the prescribed downward boundary Fourier coefficients. */
  [[nodiscard]] const rtepack::stokvec_tensor3& downward_boundary() const { return boundary_down; }
  /** Return writable prescribed downward boundary Fourier coefficients. */
  [[nodiscard]] rtepack::stokvec_tensor3_view downward_boundary() { return boundary_down; }
  /** Return the input Stokes source-function polynomial coefficients. */
  [[nodiscard]] const rtepack::stokvec_matrix& source_poly() const { return source_poly_coeffs; }
  /** Return writable Stokes source-function polynomial coefficients. */
  [[nodiscard]] rtepack::stokvec_matrix_view source_poly() { return source_poly_coeffs; }
  /** Return the surface's reflection on the streams. */
  [[nodiscard]] const surface_reflection& surface() const { return surface_refl; }
  /** Set the surface's reflection on the streams in place - [2, NBDRF, N, N], [2, NBDRF, N] and [N] or [0] as the
   * solver's, for the beam cosine beam_mu. */
  void surface(const rtepack::muelmat_tensor4_const_view& brdf,
               const rtepack::muelmat_tensor3_const_view& beam_brdf,
               Numeric                                    beam_mu,
               const rtepack::muelmat_vector_const_view&  specular) {
    surface_refl.set(brdf, beam_brdf, beam_mu, specular);
  }
  /** Set the surface's reflection on the streams in place from one of the solver's shapes. */
  void surface(const surface_reflection& x) { surface_refl.set(x); }
  /** Return the incident beam irradiance Stokes vector. */
  [[nodiscard]] const rtepack::stokvec& beam_source() const { return beam_stokes; }
  /** Return a writable incident beam irradiance Stokes vector. */
  [[nodiscard]] rtepack::stokvec& beam_source() { return beam_stokes; }
  /** Replace the beam Stokes vector and update the beam-presence flag. */
  void set_beam_source(rtepack::stokvec beam);

  /** Return the cosine of the incident beam zenith angle. */
  [[nodiscard]] Numeric solar_zenith() const { return mu0; }
  /** Return a writable cosine of the incident beam zenith angle. */
  [[nodiscard]] Numeric& solar_zenith() { return mu0; }
  /** Return the incident beam azimuth in radians. */
  [[nodiscard]] Numeric beam_azimuth() const { return phi0; }
  /** Return a writable incident beam azimuth in radians. */
  [[nodiscard]] Numeric& beam_azimuth() { return phi0; }
};
}  // namespace vdisort
