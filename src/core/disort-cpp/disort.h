#pragma once

#include <matpack.h>
#include <operators.h>
#include <xml.h>

#include <format>
#include <iosfwd>
#include <string_view>

#include "common.h"
#include "disort-eigen.h"

namespace disort {
struct radiances {
  AscendingGrid  freq_grid;  // nf
  DescendingGrid alt_grid;   // level; nl
  AziGrid        azi_grid;   // naa
  ZenGrid        zen_grid;   // nza
  Tensor4        data;       // nf, nl - 1, naa, nza

  /** Replace all grids and resize the radiance tensor to their Cartesian product. */
  void resize(AscendingGrid freq_grid, DescendingGrid alt_grid, AziGrid azi_grid, ZenGrid zen_grid);

  /** Merge radiances defined on adjacent or overlapping compatible grids. */
  [[nodiscard]] radiances combine(const radiances& other) const;

  /** Reorder the zenith dimension to match the supplied solver-stream ordering. */
  void sort(const Vector& solver_mu);
};

struct fluxes {
  AscendingGrid  freq_grid;     // nf
  DescendingGrid alt_grid;      // level; nl
  Matrix         up;            // nf, nl - 1
  Matrix         down_diffuse;  // nf, nl - 1
  Matrix         down_direct;   // nf, nl - 1
  Matrix         dfdt;          // nf, nl - 1

  /** Replace both grids and resize all four flux matrices. */
  void resize(AscendingGrid freq_grid, DescendingGrid alt_grid);

  /** Merge fluxes defined on adjacent or overlapping compatible grids. */
  [[nodiscard]] fluxes combine(const fluxes& other) const;
};

/** A surface's reflection into upward directions, in DISORT's Fourier convention.
 *
 * brdf [NBRDF, n_out, N]: the Fourier modes R_m(mu_o, mu_j) of the BRDF (see
 *   brdf::fourier_modes), m < NBRDF, from the downward streams of cosine mu_j
 *   into the upward directions of cosine mu_o, without the solver's quadrature
 *   factor (1 + delta_m0) W_j mu_j.  Higher modes reflect nothing.
 * beam_brdf [NBRDF, n_out]: R_m(mu_o, mu0) from the direct beam, at the beam
 *   cosine beam_mu; beam_mu is NaN when it was computed without a beam, and a
 *   solver with a beam refuses a surface whose beam_mu is not its own.
 *
 * The shapes are those it is constructed with.  Its data is set in place,
 * through the views or set(), so nothing is reallocated; beam_mu is set only
 * with the beam reflection it belongs to.  A surface without modes reflects
 * nothing, whatever its directions.
 *
 * On the solver's upward streams (n_out = N) it is main_data's surface, which
 * the settings set for each frequency.  At user directions it is the
 * reflection u_user() needs there: like() makes one for the modes, streams and
 * beam of the solver's surface, or brdf:: evaluates a surface model there.
 */
class surface_reflection {
  Tensor3 kernel{};       // [NBRDF, n_out, N]
  Matrix  beam_kernel{};  // [NBRDF, n_out]
  Numeric beam_cosine{std::numeric_limits<Numeric>::quiet_NaN()};

 public:
  /** No reflection */
  surface_reflection() = default;

  /** Zero reflection of NBRDF modes from N streams into n_out directions, for the beam cosine beam_mu (NaN: none) */
  surface_reflection(Index NBRDF, Index n_out, Index N, Numeric beam_mu = std::numeric_limits<Numeric>::quiet_NaN());

  /** The reflection brdf [NBRDF, n_out, N] and beam_brdf [NBRDF, n_out], for the beam cosine beam_mu (NaN: none) */
  surface_reflection(Tensor3 brdf, Matrix beam_brdf, Numeric beam_mu = std::numeric_limits<Numeric>::quiet_NaN());

  /** Zero reflection of this one's modes from its streams into n_out directions, for its beam cosine */
  [[nodiscard]] surface_reflection like(Index n_out) const;

  /** The number of Fourier modes, NBRDF */
  [[nodiscard]] Index modes() const { return kernel.npages(); }

  /** The number of upward directions reflected into, n_out */
  [[nodiscard]] Index directions() const { return kernel.nrows(); }

  /** The number of downward streams reflected from, N */
  [[nodiscard]] Index streams() const { return kernel.ncols(); }

  /** The beam cosine of the beam reflection, NaN if computed without a beam */
  [[nodiscard]] Numeric beam_mu() const { return beam_cosine; }

  //! [NBRDF, n_out, N]
  [[nodiscard]] ConstTensor3View brdf() const { return kernel; }
  [[nodiscard]] Tensor3View      brdf() { return kernel; }

  //! [NBRDF, n_out]
  [[nodiscard]] ConstMatrixView beam_brdf() const { return beam_kernel; }
  [[nodiscard]] MatrixView      beam_brdf() { return beam_kernel; }

  /** Set the reflection in place, for the beam cosine beam_mu (NaN: none); the shapes must be this one's */
  void set(const ConstTensor3View& brdf, const ConstMatrixView& beam_brdf, Numeric beam_mu);

  /** Set the reflection in place from one of the same shapes */
  void set(const surface_reflection& other);
};

/** A delta-M removed forward peak represented by its total fraction and
 * normalized Legendre moments.  Classical delta-M uses moments equal to one. */
struct delta_m_scaling {
  Vector fraction;
  Matrix moments;
};

/** Construct the DISORT 4 delta-M-plus Gaussian forward-peak model.
 * The input must contain phase moments through degree `nleg + 1`. */
delta_m_scaling delta_m_plus(ConstMatrixView phase_moments, Index nleg);

struct mathscr_v_data {
  Matrix         src;
  Matrix         G;
  solve_workdata solve_work;
  Vector         k1, k2;
  Vector         cvec;

  /** Allocate ordinary-source scratch arrays for Nk streams and Nc coefficients. */
  mathscr_v_data(const Index Nk = 0, const Index Nc = 0)
      : src(Nk, Nc), G(Nk, Nk), solve_work(Nk), k1(Nk), k2(Nk), cvec(Nc) {}
  /** Copy ordinary-source scratch storage. */
  mathscr_v_data(const mathscr_v_data&) = default;
  /** Move ordinary-source scratch storage. */
  mathscr_v_data(mathscr_v_data&&) = default;
  /** Copy-assign ordinary-source scratch storage. */
  mathscr_v_data& operator=(const mathscr_v_data&) = default;
  /** Move-assign ordinary-source scratch storage. */
  mathscr_v_data& operator=(mathscr_v_data&&) = default;

  /** Resize every ordinary-source work array consistently. */
  void resize(const Index Nk, const Index Nc) {
    src.resize(Nk, Nc);
    G.resize(Nk, Nk);
    solve_work.resize(Nk);
    k1.resize(Nk);
    k2.resize(Nk);
    cvec.resize(Nc);
  }
};

struct u_data {
  Matrix         exponent;
  Matrix         um;
  mathscr_v_data src;
  Vector         intensities;
};

/** Reusable scratch storage for DISORT user-angle radiances. */
struct user_u_data {
  mathscr_v_data source;
  Vector         barycentric_weights;
  Vector         particular;
  Vector         intensities;
};

struct u0_data {
  mathscr_v_data src;
  Vector         exponent;
  Vector         u0;
};

struct tms_data {
  Vector nu;
  Vector TMS;
  Matrix mathscr_B;
};

struct flux_data {
  u0_data u0;
};

using flux_values = disort_common::flux_values;

struct coupling_result {
  Index   iterations{0};
  Numeric max_relative_change{0.0};
  bool    converged{false};
};

using ims_convention = disort_common::ims_convention;

/** The main data structure for the DISORT algorithm
 * 
 * This pre-allocates all the memory needed for the computations
 *
 * Note to cite the original author of the algorithm when this is presented:
 * https://github.com/LDEO-CREW/Pythonic-DISORT
 *
 * Also note to present the original author of the new Legendre Gauss algorithm:
 * https://doi.org/10.1137/140954969
 */
class main_data {
  Index NLayers{0};
  Index NQuad{0};
  Index NLeg{0};
  Index NFourier{0};
  Index N{0};
  Index Nscoeffs{0};
  Index NLeg_all{0};
  Index NBDRF{0};
  bool  has_source_poly{false};
  bool  has_beam_source{false};

  //! User inputs
  AscendingGrid      tau_arr{};             // [NLayers]
  Vector             omega_arr{};           // [NLayers]
  Vector             f_arr{};               // [NLayers]
  Matrix             delta_m_peak{};        // [NLayers, NLeg]
  Matrix             source_poly_coeffs{};  // [NLayers, Nscoeffs]
  Matrix             Leg_coeffs_all{};      // [NLayers, NLeg_all]
  Matrix             boundary_up{};         // [NFourier, N], beneath the surface, transmitted with 1 - R
  Matrix             boundary_down{};       // [NFourier, N]
  surface_reflection surface_refl{};        // [NBDRF, N, N] and [NBDRF, N] on the streams
  Numeric            mu0{};
  Numeric            I0{};
  Numeric            phi0{};

  //! Derived values
  Vector             scale_tau{};                   // [NLayers]
  Vector             scaled_omega_arr{};            // [NLayers]
  Vector             scaled_tau_arr_with_0{};       // [NLayers + 1]
  Matrix             scaled_source_poly_coeffs{};   // [NLayers, Nscoeffs], emission polynomial in scaled tau
  Vector             mu_arr{};                      // [NQuad]
  Vector             inv_mu_arr{};                  // [NQuad]
  Vector             W{};                           // [N]
  Matrix             Leg_coeffs_residue_avg{};      // [NLayers + 1, NLeg_all]
  Vector             IMS_scalar{};                  // [NLayers + 1]
  Matrix             weighted_scaled_Leg_coeffs{};  // [NLayers, NLeg]
  Matrix             weighted_Leg_coeffs_all{};     // [NLayers, NLeg_all]
  Tensor4            GC_collect{};               // [NFourier, NLayers, NQuad, NQuad], cached basis-column contributions
  Tensor3            C_collect{};                // [NFourier, NLayers, NQuad], homogeneous-basis coefficients
  Tensor4            G_collect{};                // [NFourier, NLayers, NQuad, NQuad]
  Tensor3            K_collect{};                // [NFourier, NLayers, NQuad]
  Tensor3            expK_collect{};             // [NFourier, NLayers, NQuad]
  Tensor3            exponent{};                 // [NLayers, NFourier, NQuad]
  Tensor3            um{};                       // [NLayers, NFourier, NQuad]
  Tensor3            B_collect{};                // [NFourier, NLayers, NQuad]
  Tensor3            source_collect{};           // [NLayers, NQuad, Nscoeffs], polynomial particular solution
  Tensor3            transport_matrix{};         // [NLayers, NQuad, NQuad], zeroth-mode first-order operator
  std::vector<Index> conservative_pair_index{};  // [NLayers], first column of a stable near-zero pair, or -1
  Vector             conservative_pair_kappa{};  // [NLayers]
  Numeric            I0_orig{};
  Vector             scaled_mu0{};  // [NLayers + 1]

  //! Internal compute data
  Index  n{};                           // NQuad * NLayers;
  Vector RHS{};                         // [n]
  Vector jvec{};                        // [NQuad]
  Vector fac{};                         // [NLeg]
  Vector weighted_asso_Leg_coeffs_l{};  // [NLeg]
  Vector asso_leg_term_mu0{};           // [NLeg]
  Vector X_temp{};                      // [NLeg]
  Vector mathscr_X_pos{};               // [N]
  Vector E_Lm1L{};                      // [N]
  Vector E_lm1l{};                      // [N]
  Vector E_llp1{};                      // [N]
  Vector BDRF_RHS_contribution{};       // [N]
  Vector SRCB{};                        // [NQuad]
  Matrix SRC0{};                        // [NLayers, NQuad]
  Matrix SRC1{};                        // [NLayers, NQuad]
  Matrix Gml{};                         // [NQuad, NQuad]
  Matrix BDRF_LHS{};                    // [N, NQuad]
  Matrix R{};                           // [N, N]
  Matrix D_pos{};                       // [N, N]
  Matrix D_neg{};                       // [N, N]
  Matrix apb{};                         // [N, N]
  Matrix amb{};                         // [N, N]
  Matrix sqr{};                         // [N, N]
  Matrix asso_leg_term_pos{};           // [N, NLeg]
  Matrix asso_leg_term_neg{};           // [N, NLeg]
  Matrix D_temp{};                      // [N, NLeg]

  //! [NQuad / 2]
  solve_workdata solve_work{};

  //! [4 * N] + [N, N]
  real_diagonalize_workdata diag_work{};

  //! [n, 9 * N - 2] + [n / 2]
  matpack::band_matrix LHSB{};

  //! [NQuad, Nscoeffs] + [NQuad, NQuad] + 3 * [Nquad] + [Nscoeffs]
  mathscr_v_data comp_data{};

  /** Cache angle-dependent terms required to evaluate TMS at multiple depths. */
  void prepare_TMS(tms_data& data, Numeric phi, const ConstVectorView& mu) const;
  /** Require classical delta-M peak moments before evaluating scalar IMS/TMS. */
  void check_classical_delta_m_correction() const;
  /** Evaluate prepared TMS data at one optical depth. */
  void evaluate_TMS(tms_data& data, Numeric tau, const ConstVectorView& mu) const;

  /** Evaluate one scalar homogeneous eigenmode component. */
  [[nodiscard]] Numeric homogeneous(Index m, Index layer, Index state, Index eigen, Numeric tau) const;
  /** Assemble the complete homogeneous field for one Fourier mode and layer. */
  void homogeneous_field(VectorView out, Index m, Index layer, Numeric tau) const;
  /** Evaluate the polynomial-source particular solution for one stream. */
  [[nodiscard]] Numeric source_particular(Index layer, Index state, Numeric tau) const;
  /** Evaluate the beam and polynomial particular solution for one stream. */
  [[nodiscard]] Numeric particular(Index m, Index layer, Index state, Numeric tau) const;

 public:
  friend struct std::formatter<main_data>;

  /** Construct an empty, uninitialized solver. */
  main_data() = default;
  /** Copy a solver and all cached solution data. */
  main_data(const main_data&) = default;
  /** Move a solver and all cached solution data. */
  main_data(main_data&&) = default;
  /** Copy-assign a solver and all cached solution data. */
  main_data& operator=(const main_data&) = default;
  /** Move-assign a solver and all cached solution data. */
  main_data& operator=(main_data&&) = default;

  /** Allocate solver input and work arrays without computing a solution. */
  main_data(const Index NLayers,
            const Index NQuad,
            const Index NLeg,
            const Index NFourier,
            const Index Nscoeffs,
            const Index NLeg_all,
            const Index NBDRF);

  /** Construct and solve a scalar discrete-ordinate problem from complete inputs. */
  main_data(const Index               NQuad,
            const Index               NLeg,
            const Index               NFourier,
            AscendingGrid             tau_arr,
            Vector                    omega_arr,
            Matrix                    Leg_coeffs_all,
            Matrix                    boundary_up,
            Matrix                    boundary_down,
            Vector                    f_arr,
            Matrix                    source_poly_coeffs,
            const surface_reflection& surface,
            Numeric                   mu0,
            Numeric                   I0,
            Numeric                   phi0,
            Matrix                    delta_m_peak = {});

  /** Return the layer containing a physical optical depth.
    *
    * Throws if tau is out-of-bounds
    *
    * Safe for parallel use
    *
    * @param tau The point-wise optical thickness
    * @return Layer index, assigning an interface to the layer below it
    */
  [[nodiscard]] Index tau_index(const Numeric tau) const;

  /** Get the TMS correction factor
    *
    * Called by u_corr, exists for testing purposes
    *
    * Safe for parallel use if data is unique for each call
    *
    * @param data Compute data, the main output is in data.TMS
    * @param tau The point-wise optical thickness
    * @param phi The azimuthal angle of observation [0, 2 * pi)
    */
  void TMS(tms_data& data, const Numeric tau, const Numeric phi) const;
  /** TMS correction at arbitrary signed polar-angle cosines. */
  void TMS(tms_data& data, Numeric tau, Numeric phi, const ConstVectorView& mu) const;
  /** Evaluate TMS at every layer bottom, azimuth, and quadrature direction. */
  void gridded_TMS(Tensor3View tms, const Vector& phi) const;

  /** Get the IMS correction factor
    *
    * Called by u_corr, exists for testing purposes
    *
    * Safe for parallel use if ims is unique for each call
    *
    * @param ims Compute data
    * @param tau The point-wise optical thickness
    * @param phi The azimuthal angle of observation [0, 2 * pi)
    */
  void IMS(Vector& ims, Numeric tau, Numeric phi, ims_convention convention = ims_convention::disort) const;
  /** IMS correction at arbitrary signed polar-angle cosines.
   *
   * The result has the same size as mu and is zero in upward directions.
   */
  void IMS(Vector&                ims,
           Numeric                tau,
           Numeric                phi,
           const ConstVectorView& mu,
           ims_convention         convention = ims_convention::disort) const;
  /** Evaluate IMS at every layer bottom, azimuth, and quadrature direction. */
  void gridded_IMS(Tensor3View ims, const Vector& phi, ims_convention convention = ims_convention::disort) const;

  /** Spectral radiance at a given tau and phi
    *
    * The computations are done at the quadrature points
    *
    * Safe for parallel use if u_data is unique per thread
    * 
    * @param data Compute data, the spectral radiance is in data.intensities
    * @param tau The point-wise optical thickness 
    * @param phi The azimuthal angle of observation [0, 2 * pi)
    */
  void u(u_data& data, const Numeric tau, const Numeric phi) const;

  /** Radiance at arbitrary nonzero polar-angle cosines.
   *
   * Reconstructs the angular source from the discrete-ordinate solution and
   * formally integrates it along each ray, following DISORT's TERPEV,
   * TERPSO, and USRINT path.
   *
   * user_surface [NBDRF, NUser, N]: the surface's reflection into the user
   * directions (only the upward ones are read, so it may be evaluated at
   * |mu| for all).  It is required, with the beam's cosine, for an upward
   * user direction over a reflecting surface; otherwise this throws.
   */
  void u_user(user_u_data&              data,
              Numeric                   tau,
              Numeric                   phi,
              const ConstVectorView&    user_mu,
              const surface_reflection& user_surface = {}) const;

  /** Corrected radiance at arbitrary nonzero polar-angle cosines; user_surface as for u_user(). */
  void u_user_corr(user_u_data&              data,
                   Vector&                   ims,
                   tms_data&                 tms,
                   Numeric                   tau,
                   Numeric                   phi,
                   const ConstVectorView&    user_mu,
                   const surface_reflection& user_surface = {}) const;

  /** Spectral radiance at a given tau, only for the 0th Fourier mode
    *
    * The computations are done at the quadrature points
    *
    * Safe for parallel use if u0_data is unique per thread
    * 
    * @param data Compute data, the spectral radiance is in data.u0
    * @param tau The point-wise optical thickness 
    */
  void u0(u0_data& data, const Numeric tau) const;

  /** Spectral radiance at a given tau and phi
    *
    * The computations are done at the quadrature points
    *
    * Both the IMS and TMS corrections are applied
    *
    * Safe for parallel use if u_data, ims, and tms_data are unique per thread
    * 
    * @param u_data Compute data, the corrected spectral radiance is in data.intensities
    * @param ims The IMS correction compute data
    * @param tms_data The TMS correction compute data
    * @param tau The point-wise optical thickness 
    * @param phi The azimuthal angle of observation [0, 2 * pi)
    */
  void u_corr(u_data&        u_data,
              Vector&        ims,
              tms_data&      tms_data,
              Numeric        tau,
              Numeric        phi,
              ims_convention convention = ims_convention::disort) const;
  /** Evaluate IMS/TMS-corrected radiances at all layer bottoms and quadrature directions. */
  void gridded_u_corr(Tensor3View    u_data,
                      Tensor3View    tms,
                      Tensor3View    ims,
                      const Vector&  phi,
                      ims_convention convention = ims_convention::disort) const;

  /** Compute all flux quantities at a given tau from one zeroth-mode evaluation. */
  [[nodiscard]] flux_values flux(flux_data&, Numeric tau) const;

  /** Return the cached Fourier radiances at the bottom of one layer. */
  [[nodiscard]] ConstMatrixView layer_um(Size l) const;
  /** Evaluate radiances at every layer bottom for all quadrature streams and azimuths. */
  void gridded_u(Tensor3View, const Vector& phi) const;
  /** Evaluate every flux quantity at every layer bottom. */
  void gridded_flux(VectorView up, VectorView down_diffuse, VectorView down_direct, VectorView dfdt) const;

  /** Evaluate quadrature-stream radiances on an arbitrary ascending optical-depth grid. */
  void ungridded_u(Tensor3View out, const AscendingGrid& tau, const Vector& phi) const;
  /** Evaluate every flux quantity on an arbitrary ascending optical-depth grid. */
  void ungridded_flux(
      VectorView up, VectorView down_diffuse, VectorView down_direct, VectorView dfdt, const AscendingGrid& tau) const;

  /** Computes the IMS correction factors
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - omega_arr
    * - tau_arr
    * - f_arr
    * - mu0
    * - Leg_coeffs_all
    *
    * Modifies:
    * - scaled_mu0
    * - Leg_coeffs_residue_avg
    * - IMS_scalar
    */
  void set_ims_factors();

  /** Set the scaled objects
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    *
    * Depends on:
    * - omega_arr
    * - tau_arr
    * - f_arr
    * - Leg_coeffs_all
    *
    * Modifies:
    * - scale_tau
    * - scaled_tau_arr_with_0
    * - weighted_scaled_Leg_coeffs
    * - scaled_omega_arr
    */
  void set_scales();

  /** Diagonalize the system of equations
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - weighted_scaled_Leg_coeffs 
    * - scaled_omega_arr 
    * - mu0 
    * - I0
    *
    * Modifies:
    * - G_collect 
    * - K_collect 
    * - B_collect 
    */
  void diagonalize();

  /** Compute the transmission coefficients
    *
    * This function depends on diagonalize.
    * 
    * Depends on:
    * - K_collect
    *
    * Modifies:
    * - exp_K
    */
  void transmission();

  /** Cache the Fourier radiance field at every layer bottom. */
  void rad_field();

  /** Compute the source function
    *
    * This function depends on diagonalize.
    * 
    * Depends on:
    * - G_collect[0]
    * - K_collect[0]
    * - source_poly_coeffs 
    * - omega_arr
    * - tau_arr
    * - inv_mu_arr
    *
    * Modifies:
    * - SRC0
    * - SRC1
    */
  void source_function();

  /** Solves the system of equations
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - G_collect 
    * - K_collect 
    * - B_collect 
    * - source_poly_coeffs 
    * - boundary_up 
    * - boundary_down 
    * - tau_arr 
    * - scaled_tau_arr_with_0 
    * - surface_refl
    * - mu0 
    * - I0 
    *
    * Modifies:
    * - GC_collect 
    */
  void solve_for_coefs();

  /** Set the weighted Leg coeffs all object
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - Leg_coeffs_all
    *
    * Modifies:
    * - weighted_Leg_coeffs_all
    */
  void set_weighted_Leg_coeffs_all();

  /** Set the beam origin
    *
    * If you have manually changed any of the "Depends on" parameters,
    * you should call this function to update the internal data,
    * otherwise, the results will be incorrect.
    *
    * The "Modifies" parameters are updated by this function.
    *
    * Note that it is generally better to call "update_all" than
    * this method, as it will update all values correctly.
    *
    * Not safe for parallel use.
    * 
    * Depends on:
    * - boundary_up
    *
    * Modifies:
    * - I0_orig
    * - I0
    * 
    * @param I0 The new beam intensity
    */
  void set_beam_source(const Numeric I0);

  //! Checks the input sizes
  void check_input_size() const;

  //! Checks the input values
  void check_input_value() const;

  /** Updates all the internal data
    *
    * This function should be called after changing any of the input values
    * to ensure that the internal data is consistent.
    *
    * Additionally, this calls check_input_value to ensure that the input
    * values are valid.
    *
    * Not safe for parallel use.
    * 
    * @param I0 The new beam intensity if it should be changed, otherwise -1
    */
  void update_all(const Numeric I0 = -1);

  //! The angles of quadrature - NQuad
  [[nodiscard]] auto&& mu() const { return mu_arr; }

  //! The weights of quadrature - N
  [[nodiscard]] auto&& weights() const { return W; }

  //! The optical thicknesses grid - NLayers
  [[nodiscard]] auto&& tau() const { return tau_arr; }

  //! The cumulative delta-scaled optical thicknesses, including zero - NLayers + 1
  [[nodiscard]] auto&& scaled_tau() const { return scaled_tau_arr_with_0; }

  //! The single scattering albedo - NLayers
  [[nodiscard]] auto&& omega() const { return omega_arr; }

  //! The fractional scattering into the peak - NLayers or 0
  [[nodiscard]] auto&& f() const { return f_arr; }

  //! Legendre coefficients of the normalized removed forward peak - NLayers x NLeg
  [[nodiscard]] auto&& delta_m_peak_moments() const { return delta_m_peak; }

  //! Polynomial coefficients of isotropic internal sources - NLayers x Nscoeffs or 0 x 0
  [[nodiscard]] auto&& source_poly() const { return source_poly_coeffs; }

  //! Legendre coefficients of the scattering phase function (unweighted) - NLayers x NLeg_all
  [[nodiscard]] auto&& all_legendre_coeffs() const { return Leg_coeffs_all; }

  //! The radiance beneath the surface on the upward streams, Fourier modes - NFourier x N.  The surface transmits
  //! it with 1 - R, which for a Planck radiance B(T) is Kirchhoff's emission (1 - albedo(mu)) B(T)
  [[nodiscard]] auto&& upward_boundary() const { return boundary_up; }

  //! The radiance incident at the top on the downward streams, Fourier modes - NFourier x N
  [[nodiscard]] auto&& downward_boundary() const { return boundary_down; }

  //! The surface's reflection on the streams
  [[nodiscard]] const surface_reflection& surface() const { return surface_refl; }

  //! The cosine of the beam source zenith angle
  [[nodiscard]] auto&& solar_zenith() const { return mu0; }

  //! The intensity of the beam source
  [[nodiscard]] auto&& beam_source() const { return I0; }

  //! The azimuthal angle of the beam source
  [[nodiscard]] auto&& beam_azimuth() const { return phi0; }

  //! Set the optical thicknesses grid - NLayers
  template <matpack::exact_md<Numeric, 1> T> void tau(T&& x) {
    ARTS_USER_ERROR_IF(x.size() != tau_arr.size(), "Invalid size for tau: {}", x.size());
    tau_arr = std::forward<T>(x);
  }

  //! The single scattering albedo - NLayers
  [[nodiscard]] auto omega() { return VectorView{omega_arr}; }

  //! The fractional scattering into the peak - NLayers or 0
  [[nodiscard]] auto f() { return VectorView{f_arr}; }

  //! Legendre coefficients of the normalized removed forward peak - NLayers x NLeg
  [[nodiscard]] auto delta_m_peak_moments() { return MatrixView{delta_m_peak}; }

  //! Polynomial coefficients of isotropic internal sources - NLayers x Nscoeffs or 0 x 0
  [[nodiscard]] auto source_poly() { return MatrixView{source_poly_coeffs}; }

  //! Legendre coefficients of the scattering phase function (unweighted) - NLayers x NLeg_all
  [[nodiscard]] auto all_legendre_coeffs() { return MatrixView{Leg_coeffs_all}; }

  //! The radiance beneath the surface on the upward streams, Fourier modes - NFourier x N.  The surface transmits
  //! it with 1 - R, which for a Planck radiance B(T) is Kirchhoff's emission (1 - albedo(mu)) B(T)
  [[nodiscard]] auto upward_boundary() { return MatrixView{boundary_up}; }

  //! The radiance incident at the top on the downward streams, Fourier modes - NFourier x N
  [[nodiscard]] auto downward_boundary() { return MatrixView{boundary_down}; }

  //! Set the surface's reflection on the streams in place - [NBDRF, N, N] and [NBDRF, N], for the beam cosine beam_mu
  void surface(const ConstTensor3View& brdf, const ConstMatrixView& beam_brdf, Numeric beam_mu) {
    surface_refl.set(brdf, beam_brdf, beam_mu);
  }

  //! Set the surface's reflection on the streams in place from one of the solver's shapes
  void surface(const surface_reflection& x) { surface_refl.set(x); }

  //! The cosine of the beam source zenith angle
  [[nodiscard]] Numeric& solar_zenith() { return mu0; }

  //! The azimuthal angle of the beam source
  [[nodiscard]] Numeric& beam_azimuth() { return phi0; }

  //! Get weights on a grid
  [[nodiscard]] ZenGriddedField1 gridded_weights() const;
};

/** Iteratively exchange interface boundary conditions between two DISORT models
  *
  * The atmosphere lower boundary is updated from the subsurface top interface
  * field, and the subsurface upper boundary is updated from the atmosphere
  * bottom interface field. The exchanged quantities are the Fourier-mode
  * boundary arrays expected by main_data::positive_boundary and
  * main_data::negative_boundary.  The atmosphere's surface reflection R is
  * the interface's: both directions cross it with 1 - R, the atmosphere's
  * radiance here and the subsurface's in the atmosphere's solver.
  *
  * Both models must use the same quadrature and Fourier dimensions.
  *
  * @param atmosphere Atmospheric DISORT model. Its positive boundary is updated.
  * @param subsurface Subsurface DISORT model. Its negative boundary is updated.
  * @param tolerance Relative convergence tolerance for exchanged boundaries.
  * @param max_iterations Maximum number of fixed-point iterations.
  * @param relaxation Under-relaxation factor in (0, 1].
  */
[[nodiscard]] coupling_result couple(main_data& atmosphere,
                                     main_data& subsurface,
                                     Numeric    tolerance      = 1e-6,
                                     Index      max_iterations = 16,
                                     Numeric    relaxation     = 1.0);
}  // namespace disort

struct DisortSettings {
  Index quadrature_dimension{0};
  Index legendre_polynomial_dimension{0};
  Index fourier_mode_dimension{0};

  // Grids
  AscendingGrid  freq_grid{};
  DescendingGrid alt_grid{};  // levels not layers

  // freq_grid.size()
  Vector solar_azimuth_angle{};

  // freq_grid.size()
  Vector solar_zenith_angle{};

  // freq_grid.size()
  Vector solar_source{};

  // freq_grid.size() x nbrdf x quadrature_dimension / 2 x quadrature_dimension / 2, the Fourier modes of the
  // surface's BRDF between the streams (disort::surface_reflection)
  Tensor4 surface_brdf{};

  // freq_grid.size() x nbrdf x quadrature_dimension / 2, from the direct beam at the cosine surface_beam_mu
  Tensor3 surface_beam_brdf{};

  // freq_grid.size(), the beam cosine of surface_beam_brdf, NaN where it was computed without a beam
  Vector surface_beam_mu{};

  // freq_grid.size() x [alt_grid.size() - 1]
  Matrix optical_thicknesses{};

  // freq_grid.size() x [alt_grid.size() - 1]
  Matrix single_scattering_albedo{};

  // freq_grid.size() x [alt_grid.size() - 1]
  Matrix fractional_scattering{};

  // freq_grid.size() x [alt_grid.size() - 1] x legendre_polynomial_dimension
  Tensor3 delta_m_peak_moments{};

  // freq_grid.size() x [alt_grid.size() - 1] x nsrc
  Tensor3 source_polynomial{};

  // freq_grid.size() x [alt_grid.size() - 1] x legendre_polynomial_dimension_full <- last is unknown at construction, must be larger or equal to legendre_polynomial_dimension
  Tensor3 legendre_coefficients{};

  // freq_grid.size() x fourier_mode_dimension x quadrature_dimension / 2, the radiance beneath the surface, which
  // DISORT transmits with 1 - R (Kirchhoff's emission for a Planck radiance)
  Tensor3 upward_boundary_condition{};

  // freq_grid.size() x fourier_mode_dimension x quadrature_dimension / 2.
  Tensor3 downward_boundary_condition{};

  /** Construct empty DISORT settings. */
  DisortSettings() = default;

  /** Replace grids and dimensions, then resize every frequency- and layer-dependent input. */
  void resize(Index          quadrature_dimension,
              Index          legendre_polynomial_dimension,
              Index          fourier_mode_dimension,
              AscendingGrid  freq_grid,
              DescendingGrid alt_grid_);

  /** Return the number of configured frequencies. */
  [[nodiscard]] Index frequency_count() const { return freq_grid.size(); }
  /** Return the number of atmospheric layers. */
  [[nodiscard]] Index layer_count() const { return alt_grid.size() - 1; }

  /** Allocate a scalar solver with the dimensions represented by these settings. */
  [[nodiscard]] disort::main_data init() const;
  /** Populate and solve an existing scalar solver for frequency index iv. */
  disort::main_data& set(disort::main_data&, Index iv) const;
#ifdef ENABLE_CDISORT
  /** Populate an existing scalar solver using cdisort-compatible truncation conventions. */
  disort::main_data& set_cdisort(disort::main_data&, Index iv) const;
#endif
  /** Return whether every member agrees with the grids and dimensions. */
  [[nodiscard]] bool ok() const;
  /** Validate settings dimensions, grids, and physical input ranges. */
  void check() const;
};

template <> struct std::formatter<disort::main_data> {
  format_tags tags;

  /** Return the mutable nested formatter used by generic formatting helpers. */
  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  /** Return the read-only nested formatter used by generic formatting helpers. */
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  /** Parse ARTS format tags for a scalar DISORT solver. */
  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  /** Format solver dimensions or the complete scalar DISORT state. */
  template <class FmtContext> FmtContext::iterator format(const disort::main_data& v, FmtContext& ctx) const {
    const auto sep = tags.sep();
    return tags.format(ctx,
                       "NLayers: "sv,
                       v.NLayers,
                       ", NQuad: "sv,
                       v.NQuad,
                       ", NLeg: "sv,
                       v.NLeg,
                       ", NFourier: "sv,
                       v.NFourier,
                       ", N: "sv,
                       v.N,
                       ", Nscoeffs: "sv,
                       v.Nscoeffs,
                       ", NLeg_all: "sv,
                       v.NLeg_all,
                       ", NBDRF: "sv,
                       v.NBDRF,
                       sep,
                       "Has source: "sv,
                       v.has_source_poly,
                       ", Has beam source: "sv,
                       v.has_beam_source,
                       sep,
                       "tau_arr: "sv,
                       v.tau_arr,
                       sep,
                       "omega_arr: "sv,
                       v.omega_arr,
                       sep,
                       "f_arr: "sv,
                       v.f_arr,
                       sep,
                       "delta_m_peak: "sv,
                       v.delta_m_peak,
                       sep,
                       "source_poly_coeffs: "sv,
                       v.source_poly_coeffs,
                       sep,
                       "Leg_coeffs_all: "sv,
                       v.Leg_coeffs_all,
                       sep,
                       "boundary_up: "sv,
                       v.boundary_up,
                       sep,
                       "boundary_down: "sv,
                       v.boundary_down,
                       sep,
                       "surface_refl.brdf: "sv,
                       v.surface_refl.brdf(),
                       sep,
                       "mu0: "sv,
                       v.mu0,
                       sep,
                       "I0: "sv,
                       v.I0,
                       sep,
                       "phi0: "sv,
                       v.phi0,
                       sep,
                       "scale_tau: "sv,
                       v.scale_tau,
                       sep,
                       "scaled_omega_arr: "sv,
                       v.scaled_omega_arr,
                       sep,
                       "scaled_tau_arr_with_0: "sv,
                       v.scaled_tau_arr_with_0,
                       sep,
                       "mu_arr: "sv,
                       v.mu_arr,
                       sep,
                       "inv_mu_arr: "sv,
                       v.inv_mu_arr,
                       sep,
                       "W: "sv,
                       v.W,
                       sep,
                       "Leg_coeffs_residue_avg: "sv,
                       v.Leg_coeffs_residue_avg,
                       sep,
                       "weighted_scaled_Leg_coeffs: "sv,
                       v.weighted_scaled_Leg_coeffs,
                       sep,
                       "weighted_Leg_coeffs_all: "sv,
                       v.weighted_Leg_coeffs_all,
                       sep,
                       "GC_collect: "sv,
                       v.GC_collect,
                       sep,
                       "G_collect: "sv,
                       v.G_collect,
                       sep,
                       "K_collect: "sv,
                       v.K_collect,
                       sep,
                       "B_collect: "sv,
                       v.B_collect,
                       sep,
                       "I0_orig: "sv,
                       v.I0_orig,
                       sep,
                       "IMS_scalar: "sv,
                       v.IMS_scalar,
                       sep,
                       "scaled_mu0: "sv,
                       v.scaled_mu0,
                       sep,
                       "n: "sv,
                       v.n,
                       sep,
                       "RHS: "sv,
                       v.RHS,
                       sep,
                       "jvec: "sv,
                       v.jvec,
                       sep,
                       "fac: "sv,
                       v.fac,
                       sep,
                       "weighted_asso_Leg_coeffs_l: "sv,
                       v.weighted_asso_Leg_coeffs_l,
                       sep,
                       "asso_leg_term_mu0: "sv,
                       v.asso_leg_term_mu0,
                       sep,
                       "X_temp: "sv,
                       v.X_temp,
                       sep,
                       "mathscr_X_pos: "sv,
                       v.mathscr_X_pos,
                       sep,
                       "E_Lm1L: "sv,
                       v.E_Lm1L,
                       sep,
                       "E_lm1l: "sv,
                       v.E_lm1l,
                       sep,
                       "E_llp1: "sv,
                       v.E_llp1,
                       sep,
                       "BDRF_RHS_contribution: "sv,
                       v.BDRF_RHS_contribution,
                       sep,
                       "Gml: "sv,
                       v.Gml,
                       sep,
                       "BDRF_LHS: "sv,
                       v.BDRF_LHS,
                       sep,
                       "R: "sv,
                       v.R,
                       sep,
                       "D_pos: "sv,
                       v.D_pos,
                       sep,
                       "D_neg: "sv,
                       v.D_neg,
                       sep,
                       "apb: "sv,
                       v.apb,
                       sep,
                       "amb: "sv,
                       v.amb,
                       sep,
                       "sqr: "sv,
                       v.sqr,
                       sep,
                       "asso_leg_term_pos: "sv,
                       v.asso_leg_term_pos,
                       sep,
                       "asso_leg_term_neg: "sv,
                       v.asso_leg_term_neg,
                       sep,
                       "D_temp: "sv,
                       v.D_temp);
  }
};

template <> struct std::formatter<DisortSettings> {
  format_tags tags;

  /** Return the mutable nested formatter used by generic formatting helpers. */
  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  /** Return the read-only nested formatter used by generic formatting helpers. */
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  /** Parse ARTS format tags for DISORT settings. */
  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  /** Format settings dimensions, validation diagnostics, or complete values. */
  template <class FmtContext> FmtContext::iterator format(const DisortSettings& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return tags.format(
          ctx,
          R"-x-(quadrature_dimension:          )-x-"sv,
          v.quadrature_dimension,
          R"-x-(
legendre_polynomial_dimension: )-x-"sv,
          v.legendre_polynomial_dimension,
          R"-x-(
fourier_mode_dimension:        )-x-"sv,
          v.fourier_mode_dimension,
          R"-x-(
freq_grid.size():              )-x-"sv,
          v.freq_grid.size(),
          R"-x-(
alt_grid.size():               )-x-"sv,
          v.alt_grid.size(),
          R"-x-(
nsrc:                          )-x-"sv,
          v.source_polynomial.ncols(),
          R"-x-(
nbrdf:                         )-x-"sv,
          v.surface_brdf.npages(),
          R"-x-(

solar_source.shape():                                     )-x-"sv,
          v.solar_source.shape(),
          R"-x-( - should be freq_grid.size().
solar_zenith_angle.shape():                               )-x-"sv,
          v.solar_zenith_angle.shape(),
          R"-x-( - should be freq_grid.size().
solar_azimuth_angle.shape():                              )-x-"sv,
          v.solar_azimuth_angle.shape(),
          R"-x-( - should be freq_grid.size().
surface_brdf.shape():                                     )-x-"sv,
          v.surface_brdf.shape(),
          R"-x-( - should be freq_grid.size() x nbrdf x (quadrature_dimension / 2) x (quadrature_dimension / 2).
surface_beam_brdf.shape():                                )-x-"sv,
          v.surface_beam_brdf.shape(),
          R"-x-( - should be freq_grid.size() x nbrdf x (quadrature_dimension / 2).
surface_beam_mu.shape():                                  )-x-"sv,
          v.surface_beam_mu.shape(),
          R"-x-( - should be freq_grid.size().
optical_thicknesses.shape():                              )-x-"sv,
          v.optical_thicknesses.shape(),
          R"-x-( - should be freq_grid.size() x [alt_grid.size() - 1].
single_scattering_albedo.shape():                         )-x-"sv,
          v.single_scattering_albedo.shape(),
          R"-x-( - should be freq_grid.size() x [alt_grid.size() - 1].
fractional_scattering.shape():                            )-x-"sv,
          v.fractional_scattering.shape(),
          R"-x-( - should be freq_grid.size() x [alt_grid.size() - 1].
delta_m_peak_moments.shape():                             )-x-"sv,
          v.delta_m_peak_moments.shape(),
          R"-x-( - should be freq_grid.size() x [alt_grid.size() - 1] x nleg.
source_polynomial.shape():                                )-x-"sv,
          v.source_polynomial.shape(),
          R"-x-( - should be freq_grid.size() x [alt_grid.size() - 1] x nsrc.
legendre_coefficients.shape():                            )-x-"sv,
          v.legendre_coefficients.shape(),
          R"-x-( - should be freq_grid.size() x [alt_grid.size() - 1] x nleg.
upward_boundary_condition.shape():                        )-x-"sv,
          v.upward_boundary_condition.shape(),
          R"-x-( - should be freq_grid.size() x [alt_grid.size() - 1] x (quadrature_dimension / 2).
downward_boundary_condition.shape():                      )-x-"sv,
          v.downward_boundary_condition.shape(),
          R"-x-( - should be freq_grid.size() x [alt_grid.size() - 1] x (quadrature_dimension / 2).)-x-"sv);
    }
    return tags.format(ctx,
                       R"-x-(quadrature_dimension:          )-x-"sv,
                       v.quadrature_dimension,
                       R"-x-(
legendre_polynomial_dimension: )-x-"sv,
                       v.legendre_polynomial_dimension,
                       R"-x-(
fourier_mode_dimension:        )-x-"sv,
                       v.fourier_mode_dimension,
                       R"-x-(
freq_grid:                     )-x-"sv,
                       v.freq_grid,
                       R"-x-(
alt_grid:                      )-x-"sv,
                       v.alt_grid,
                       R"-x-(
nsrc:                          )-x-"sv,
                       v.source_polynomial.ncols(),
                       R"-x-(
nbrdf:                         )-x-"sv,
                       v.surface_brdf.npages(),
                       R"-x-(

solar_source:
)-x-"sv,
                       v.solar_source,
                       R"-x-(
solar_zenith_angle:
)-x-"sv,
                       v.solar_zenith_angle,
                       R"-x-(
solar_azimuth_angle:
)-x-"sv,
                       v.solar_azimuth_angle,
                       R"-x-(
surface_brdf:
)-x-"sv,
                       v.surface_brdf,
                       R"-x-(
surface_beam_brdf:
)-x-"sv,
                       v.surface_beam_brdf,
                       R"-x-(
surface_beam_mu:
)-x-"sv,
                       v.surface_beam_mu,
                       R"-x-(
optical_thicknesses:
)-x-"sv,
                       v.optical_thicknesses,
                       R"-x-(
single_scattering_albedo:
)-x-"sv,
                       v.single_scattering_albedo,
                       R"-x-(
fractional_scattering:
)-x-"sv,
                       v.fractional_scattering,
                       R"-x-(
delta_m_peak_moments:
)-x-"sv,
                       v.delta_m_peak_moments,
                       R"-x-(
source_polynomial:
)-x-"sv,
                       v.source_polynomial,
                       R"-x-(
legendre_coefficients:
)-x-"sv,
                       v.legendre_coefficients,
                       R"-x-(
upward_boundary_condition:
)-x-"sv,
                       v.upward_boundary_condition,
                       R"-x-(
downward_boundary_condition:
)-x-"sv,
                       v.downward_boundary_condition);
  }
};

template <> struct xml_io_stream<DisortSettings> {
  static constexpr std::string_view type_name = "DisortSettings"sv;

  /** Write DISORT settings to XML or the associated binary stream. */
  static void write(std::ostream&         os,
                    const DisortSettings& x,
                    bofstream*            pbofs = nullptr,
                    std::string_view      name  = ""sv);

  /** Read DISORT settings from XML or the associated binary stream. */
  static void read(std::istream& is, DisortSettings& x, bifstream* pbifs = nullptr);
};

using FluxComponentsProfile = disort::fluxes;

template <> struct xml_io_stream_name<FluxComponentsProfile> {
  constexpr static std::string_view name = "FluxComponentsProfile"sv;
};

template <> struct xml_io_stream_aggregate<FluxComponentsProfile> {
  constexpr static bool value = true;
};

template <> struct std::formatter<FluxComponentsProfile> {
  format_tags tags;

  /** Return the mutable nested formatter used by generic formatting helpers. */
  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  /** Return the read-only nested formatter used by generic formatting helpers. */
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  /** Parse ARTS format tags for flux output. */
  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  /** Format flux dimensions or complete gridded flux values. */
  template <class FmtContext> FmtContext::iterator format(const FluxComponentsProfile& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return tags.format(ctx,
                         "freq_grid.shape:    "sv,
                         v.freq_grid.shape(),
                         "\n"sv,
                         "alt_grid.shape:     "sv,
                         v.alt_grid.shape(),
                         "\n"sv,
                         "down_diffuse.shape: "sv,
                         v.down_diffuse.shape(),
                         "\n"sv,
                         "down_direct.shape:  "sv,
                         v.down_direct.shape(),
                         "\n"sv,
                         "up.shape:           "sv,
                         v.up.shape(),
                         "\n"sv,
                         "dfdt.shape:         "sv,
                         v.dfdt.shape());
    }

    auto sep = tags.sep();
    return tags.format(
        ctx, v.freq_grid, sep, v.alt_grid, sep, v.up, sep, v.down_diffuse, sep, v.down_direct, sep, v.dfdt);
  }
};

using DisortRadiance = disort::radiances;

template <> struct xml_io_stream_name<DisortRadiance> {
  constexpr static std::string_view name = "DisortRadiance"sv;
};

template <> struct xml_io_stream_aggregate<DisortRadiance> {
  constexpr static bool value = true;
};

template <> struct std::formatter<DisortRadiance> {
  format_tags tags;

  /** Return the mutable nested formatter used by generic formatting helpers. */
  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  /** Return the read-only nested formatter used by generic formatting helpers. */
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  /** Parse ARTS format tags for radiance output. */
  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  /** Format radiance dimensions or complete gridded radiance values. */
  template <class FmtContext> FmtContext::iterator format(const DisortRadiance& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return tags.format(ctx,
                         "freq_grid.shape: "sv,
                         v.freq_grid.shape(),
                         "\n"sv,
                         "alt_grid.shape:  "sv,
                         v.alt_grid.shape(),
                         "\n"sv,
                         "azi_grid.shape:  "sv,
                         v.azi_grid.shape(),
                         "\n"sv,
                         "zen_grid.shape:  "sv,
                         v.zen_grid.shape(),
                         "\n"sv,
                         "data.shape:      "sv,
                         v.data.shape());
    }

    auto sep = tags.sep();
    return tags.format(ctx, v.freq_grid, sep, v.alt_grid, sep, v.azi_grid, sep, v.zen_grid, sep, v.data);
  }
};
