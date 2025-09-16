#pragma once

#include <matpack.h>
#include <operators.h>
#include <xml.h>

#include <format>
#include <iosfwd>
#include <string_view>

#include "disort-eigen.h"

namespace disort {
struct radiances {
  AscendingGrid frequency_grid;  // nf
  DescendingGrid altitude_grid;  // level; nl
  AzimuthGrid azimuth_grid;      // naa
  ZenithGrid zenith_grid;        // nza
  Tensor4 data;                  // nf, nl - 1, naa, nza

  void resize(AscendingGrid frequency_grid,
              DescendingGrid altitude_grid,
              AzimuthGrid azimuth_grid,
              ZenithGrid zenith_grid);

  void sort(const Vector& solver_mu);
};

struct fluxes {
  AscendingGrid frequency_grid;  // nf
  DescendingGrid altitude_grid;  // level; nl
  Matrix up;                     // nf, nl - 1
  Matrix down_diffuse;           // nf, nl - 1
  Matrix down_direct;            // nf, nl - 1

  void resize(AscendingGrid frequency_grid, DescendingGrid altitude_grid);
};

struct BDRF {
  using func_t = CustomOperator<void,
                                MatrixView,
                                const ConstVectorView&,
                                const ConstVectorView&>;
  func_t f;

  void operator()(MatrixView x,
                  const ConstVectorView& a,
                  const ConstVectorView& b) const;

  Matrix operator()(const Vector& a, const Vector& b) const;
};

struct mathscr_v_data {
  Matrix src;
  Matrix G;
  solve_workdata solve_work;
  Vector k1, k2;
  Vector cvec;

  mathscr_v_data(const Index Nk = 0, const Index Nc = 0)
      : src(Nk, Nc), G(Nk, Nk), solve_work(Nk), k1(Nk), k2(Nk), cvec(Nc) {}
  mathscr_v_data(const mathscr_v_data&)            = default;
  mathscr_v_data(mathscr_v_data&&)                 = default;
  mathscr_v_data& operator=(const mathscr_v_data&) = default;
  mathscr_v_data& operator=(mathscr_v_data&&)      = default;

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
  Matrix exponent;
  Matrix um;
  mathscr_v_data src;
  Vector intensities;
};

struct u0_data {
  mathscr_v_data src;
  Vector exponent;
  Vector u0;
};

struct tms_data {
  Vector nu;
  Vector TMS;
  Matrix mathscr_B;
  Matrix contribution_from_other_layers_pos;
  Matrix contribution_from_other_layers_neg;
};

struct flux_data {
  Vector exponent;
  mathscr_v_data src;
  Vector u0_pos;
  Vector u0_neg;
};

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
  bool has_source_poly{false};
  bool has_beam_source{false};

  //! User inputs
  AscendingGrid tau_arr{};                 // [NLayers]
  Vector omega_arr{};                      // [NLayers]
  Vector f_arr{};                          // [NLayers]
  Matrix source_poly_coeffs{};             // [NLayers, Nscoeffs]
  Matrix Leg_coeffs_all{};                 // [NLayers, NLeg_all]
  Matrix b_pos{};                          // [NFourier, N]
  Matrix b_neg{};                          // [NFourier, N]
  std::vector<BDRF> brdf_fourier_modes{};  // [NBDRF]
  Numeric mu0{};
  Numeric I0{};
  Numeric phi0{};

  //! Derived values
  Vector scale_tau{};                   // [NLayers]
  Vector scaled_omega_arr{};            // [NLayers]
  Vector scaled_tau_arr_with_0{};       // [NLayers + 1]
  Vector mu_arr{};                      // [NQuad]
  Vector inv_mu_arr{};                  // [NQuad]
  Vector W{};                           // [N]
  Vector Leg_coeffs_residue_avg{};      // [NLeg_all]
  Matrix weighted_scaled_Leg_coeffs{};  // [NLayers, NLeg]
  Matrix weighted_Leg_coeffs_all{};     // [NLayers, NLeg_all]
  Tensor4 GC_collect{};                 // [NFourier, NLayers, NQuad, NQuad]
  Tensor4 G_collect{};                  // [NFourier, NLayers, NQuad, NQuad]
  Tensor3 K_collect{};                  // [NFourier, NLayers, NQuad]
  Tensor3 expK_collect{};               // [NFourier, NLayers, NQuad]
  Tensor3 B_collect{};                  // [NFourier, NLayers, NQuad]
  Numeric I0_orig{};
  Numeric f_avg{};
  Numeric omega_avg{};
  Numeric scaled_mu0{};

  //! Internal compute data
  Index n{};                            // NQuad * NLayers;
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
  Matrix mathscr_D_neg{};               // [N, N]
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

 public:
  friend struct std::formatter<main_data>;

  main_data()                            = default;
  main_data(const main_data&)            = default;
  main_data(main_data&&)                 = default;
  main_data& operator=(const main_data&) = default;
  main_data& operator=(main_data&&)      = default;

  main_data(const Index NLayers,
            const Index NQuad,
            const Index NLeg,
            const Index NFourier,
            const Index Nscoeffs,
            const Index NLeg_all,
            const Index NBDRF);

  main_data(const Index NQuad,
            const Index NLeg,
            const Index NFourier,
            AscendingGrid tau_arr,
            Vector omega_arr,
            Matrix Leg_coeffs_all,
            Matrix b_pos,
            Matrix b_neg,
            Vector f_arr,
            Matrix source_poly_coeffs,
            std::vector<BDRF> brdf_fourier_modes,
            Numeric mu0,
            Numeric I0,
            Numeric phi0);

  /** Get the index of the tau value closest to the given tau
    *
    * Throws if tau is out-of-bounds
    *
    * Safe for parallel use
    *
    * @param tau The point-wise optical thickness
    * @return Index of the tau value closest to the given tau
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
  void IMS(Vector& ims, const Numeric tau, const Numeric phi) const;
  void gridded_IMS(Tensor3View ims, const Vector& phi) const;

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
  void u_corr(u_data& u_data,
              Vector& ims,
              tms_data& tms_data,
              const Numeric tau,
              const Numeric phi) const;
  void gridded_u_corr(Tensor3View u_data,
                      Tensor3View tms,
                      Tensor3View ims,
                      const Vector& phi) const;

  /** Compute the upward flux at a given tau
    *
    * Safe for parallel use if flux_data is unique per thread
    * 
    * @param tau The point-wise optical thickness 
    * @return Numeric value of the upward flux
    */
  [[nodiscard]] Numeric flux_up(flux_data&, const Numeric tau) const;

  void gridded_u(Tensor3View, const Vector& phi) const;
  void gridded_flux(VectorView up,
                    VectorView down,
                    VectorView down_direct) const;

  void ungridded_u(Tensor3View out,
                   const AscendingGrid& tau,
                   const Vector& phi) const;
  void ungridded_flux(VectorView flux_up,
                      VectorView flux_do,
                      VectorView flux_dd,
                      const AscendingGrid& tau) const;

  /** Compute the downward flux at a given tau
    *
    * Safe for parallel use if flux_data is unique per thread
    * 
    * @param tau The point-wise optical thickness 
    * @return std::pair<Numeric, Numeric> Diffuse and direct downward flux, respectively
    */
  [[nodiscard]] std::pair<Numeric, Numeric> flux_down(flux_data&,
                                                      const Numeric tau) const;

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
    * - omega_avg
    * - f_avg
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
    * - b_pos 
    * - b_neg 
    * - tau_arr 
    * - scaled_tau_arr_with_0 
    * - brdf_fourier_modes 
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
    * - b_pos
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

  //! The single scattering albedo - NLayers
  [[nodiscard]] auto&& omega() const { return omega_arr; }

  //! The fractional scattering into the peak - NLayers or 0
  [[nodiscard]] auto&& f() const { return f_arr; }

  //! Polynomial coefficients of isotropic internal sources - NLayers x Nscoeffs or 0 x 0
  [[nodiscard]] auto&& source_poly() const { return source_poly_coeffs; }

  //! Legendre coefficients of the scattering phase function (unweighted) - NLayers x NLeg_all
  [[nodiscard]] auto&& all_legendre_coeffs() const { return Leg_coeffs_all; }

  //! Positive Fourier coefficients of the beam source - NFourier x N
  [[nodiscard]] auto&& positive_boundary() const { return b_pos; }

  //! Negative Fourier coefficients of the beam source - NFourier x N
  [[nodiscard]] auto&& negative_boundary() const { return b_neg; }

  //! Fourier modes of the BRDF - NBDRF
  [[nodiscard]] auto&& brdf_modes() const { return brdf_fourier_modes; }

  //! The cosine of the beam source zenith angle
  [[nodiscard]] auto&& solar_zenith() const { return mu0; }

  //! The intensity of the beam source
  [[nodiscard]] auto&& beam_source() const { return I0; }

  //! The azimuthal angle of the beam source
  [[nodiscard]] auto&& beam_azimuth() const { return phi0; }

  //! Set the optical thicknesses grid - NLayers
  template <matpack::exact_md<Numeric, 1> T>
  void tau(T&& x) {
    ARTS_USER_ERROR_IF(
        x.size() != tau_arr.size(), "Invalid size for tau: {}", x.size());
    tau_arr = std::forward<T>(x);
  }

  //! The single scattering albedo - NLayers
  [[nodiscard]] auto omega() { return VectorView{omega_arr}; }

  //! The fractional scattering into the peak - NLayers or 0
  [[nodiscard]] auto f() { return VectorView{f_arr}; }

  //! Polynomial coefficients of isotropic internal sources - NLayers x Nscoeffs or 0 x 0
  [[nodiscard]] auto source_poly() { return MatrixView{source_poly_coeffs}; }

  //! Legendre coefficients of the scattering phase function (unweighted) - NLayers x NLeg_all
  [[nodiscard]] auto all_legendre_coeffs() {
    return MatrixView{Leg_coeffs_all};
  }

  //! Positive Fourier coefficients of the beam source - NFourier x N
  [[nodiscard]] auto positive_boundary() { return MatrixView{b_pos}; }

  //! Negative Fourier coefficients of the beam source - NFourier x N
  [[nodiscard]] auto negative_boundary() { return MatrixView{b_neg}; }

  //! Fourier modes of the BRDF - NBDRF
  [[nodiscard]] auto brdf_modes() { return std::span{brdf_fourier_modes}; }

  //! The cosine of the beam source zenith angle
  [[nodiscard]] Numeric& solar_zenith() { return mu0; }

  /** The intensity of the beam source
    * 
    * This function is disabled, use set_beam_source instead, or if
    * it is part of a larger update, use update_all.
    *
    * If you really want the beam source intensity, you can use the
    * const version of this function by first const-casting your object.
    */
  [[nodiscard]] Numeric& beam_source() = delete;

  //! The azimuthal angle of the beam source
  [[nodiscard]] Numeric& beam_azimuth() { return phi0; }

  //! Get weights on a grid
  [[nodiscard]] ZenithGriddedField1 gridded_weights() const;
};
}  // namespace disort

using DisortBDRF         = disort::BDRF;
using MatrixOfDisortBDRF = matpack::data_t<DisortBDRF, 2>;

struct DisortSettings {
  Index quadrature_dimension{0};
  Index legendre_polynomial_dimension{0};
  Index fourier_mode_dimension{0};

  // Grids
  AscendingGrid frequency_grid{};
  DescendingGrid altitude_grid{};  // levels not layers

  // frequency_grid.size()
  Vector solar_azimuth_angle{};

  // frequency_grid.size()
  Vector solar_zenith_angle{};

  // frequency_grid.size()
  Vector solar_source{};

  // frequency_grid.size() x nbrdf
  MatrixOfDisortBDRF bidirectional_reflectance_distribution_functions{};

  // frequency_grid.size() x [altitude_grid.size() - 1]
  Matrix optical_thicknesses{};

  // frequency_grid.size() x [altitude_grid.size() - 1]
  Matrix single_scattering_albedo{};

  // frequency_grid.size() x [altitude_grid.size() - 1]
  Matrix fractional_scattering{};

  // frequency_grid.size() x [altitude_grid.size() - 1] x nsrc
  Tensor3 source_polynomial{};

  // frequency_grid.size() x [altitude_grid.size() - 1] x legendre_polynomial_dimension_full <- last is unknown at construction, must be larger or equal to legendre_polynomial_dimension
  Tensor3 legendre_coefficients{};

  // frequency_grid.size() x fourier_mode_dimension x quadrature_dimension / 2
  Tensor3 positive_boundary_condition{};

  // frequency_grid.size() x fourier_mode_dimension x quadrature_dimension / 2.
  Tensor3 negative_boundary_condition{};

  DisortSettings() = default;

  void resize(Index quadrature_dimension,
              Index legendre_polynomial_dimension,
              Index fourier_mode_dimension,
              AscendingGrid frequency_grid,
              DescendingGrid altitude_grid);

  [[nodiscard]] Index frequency_count() const { return frequency_grid.size(); }
  [[nodiscard]] Index layer_count() const { return altitude_grid.size() - 1; }

  [[nodiscard]] disort::main_data init() const;
  disort::main_data& set(disort::main_data&, Index iv) const;
#ifdef ENABLE_CDISORT
  disort::main_data& set_cdisort(disort::main_data&, Index iv) const;
#endif
  void check() const;
};

template <>
struct std::formatter<DisortBDRF> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }
  template <class FmtContext>
  FmtContext::iterator format(const DisortBDRF&, FmtContext& ctx) const {
    return tags.format(ctx, "BDRF"sv);
  }
};

template <>
struct std::formatter<disort::main_data> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  template <class FmtContext>
  FmtContext::iterator format(const disort::main_data& v,
                              FmtContext& ctx) const {
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
                       "source_poly_coeffs: "sv,
                       v.source_poly_coeffs,
                       sep,
                       "Leg_coeffs_all: "sv,
                       v.Leg_coeffs_all,
                       sep,
                       "b_pos: "sv,
                       v.b_pos,
                       sep,
                       "b_neg: "sv,
                       v.b_neg,
                       sep,
                       "brdf_fourier_modes: "sv,
                       // v.brdf_fourier_modes,
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
                       "f_avg: "sv,
                       v.f_avg,
                       sep,
                       "omega_avg: "sv,
                       v.omega_avg,
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
                       "mathscr_D_neg: "sv,
                       v.mathscr_D_neg,
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

template <>
struct std::formatter<DisortSettings> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const DisortSettings& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return tags.format(
          ctx,
          R"-x-(quadrature_dimension:          ")-x-"sv,
          v.quadrature_dimension,
          R"-x-(
legendre_polynomial_dimension: ")-x-"sv,
          v.legendre_polynomial_dimension,
          R"-x-(
fourier_mode_dimension:        ")-x-"sv,
          v.fourier_mode_dimension,
          R"-x-(
frequency_grid.size():         ")-x-"sv,
          v.frequency_grid.size(),
          R"-x-(
altitude_grid.size():          ")-x-"sv,
          v.altitude_grid.size(),
          R"-x-(
nsrc:                          ")-x-"sv,
          v.source_polynomial.ncols(),
          R"-x-(
nbrdf:                         ")-x-"sv,
          v.bidirectional_reflectance_distribution_functions.ncols(),
          R"-x-(

solar_source.shape():                                     ")-x-"sv,
          v.solar_source.shape(),
          R"-x-( - should be frequency_grid.size().
solar_zenith_angle.shape():                               ")-x-"sv,
          v.solar_zenith_angle.shape(),
          R"-x-( - should be frequency_grid.size().
solar_azimuth_angle.shape():                              ")-x-"sv,
          v.solar_azimuth_angle.shape(),
          R"-x-( - should be frequency_grid.size().
bidirectional_reflectance_distribution_functions.shape(): ")-x-"sv,
          v.bidirectional_reflectance_distribution_functions.shape(),
          R"-x-( - should be frequency_grid.size() x nbrdf.
optical_thicknesses.shape():                              ")-x-"sv,
          v.optical_thicknesses.shape(),
          R"-x-( - should be frequency_grid.size() x [altitude_grid.size() - 1].
single_scattering_albedo.shape():                         ")-x-"sv,
          v.single_scattering_albedo.shape(),
          R"-x-( - should be frequency_grid.size() x [altitude_grid.size() - 1].
fractional_scattering.shape():                            ")-x-"sv,
          v.fractional_scattering.shape(),
          R"-x-( - should be frequency_grid.size() x [altitude_grid.size() - 1].
source_polynomial.shape():                                ")-x-"sv,
          v.source_polynomial.shape(),
          R"-x-( - should be frequency_grid.size() x [altitude_grid.size() - 1] x nsrc.
legendre_coefficients.shape():                            ")-x-"sv,
          v.legendre_coefficients.shape(),
          R"-x-( - should be frequency_grid.size() x [altitude_grid.size() - 1] x nleg.
positive_boundary_condition.shape():                      ")-x-"sv,
          v.positive_boundary_condition.shape(),
          R"-x-( - should be frequency_grid.size() x [altitude_grid.size() - 1] x (quadrature_dimension / 2).
negative_boundary_condition.shape():                      ")-x-"sv,
          v.negative_boundary_condition.shape(),
          R"-x-( - should be frequency_grid.size() x [altitude_grid.size() - 1] x (quadrature_dimension / 2).)-x-"sv);
    }
    return tags.format(
        ctx,
        R"-x-(quadrature_dimension:          ")-x-"sv,
        v.quadrature_dimension,
        R"-x-(
legendre_polynomial_dimension: ")-x-"sv,
        v.legendre_polynomial_dimension,
        R"-x-(
fourier_mode_dimension:        ")-x-"sv,
        v.fourier_mode_dimension,
        R"-x-(
frequency_grid:                ")-x-"sv,
        v.frequency_grid,
        R"-x-(
altitude_grid:                 ")-x-"sv,
        v.altitude_grid,
        R"-x-(
nsrc:                          ")-x-"sv,
        v.source_polynomial.ncols(),
        R"-x-(
nbrdf:                         ")-x-"sv,
        v.bidirectional_reflectance_distribution_functions.ncols(),
        R"-x-(

solar_source:
")-x-"sv,
        v.solar_source,
        R"-x-(
solar_zenith_angle:
")-x-"sv,
        v.solar_zenith_angle,
        R"-x-(
solar_azimuth_angle:
")-x-"sv,
        v.solar_azimuth_angle,
        R"-x-(
bidirectional_reflectance_distribution_functions:
")-x-"sv,
        v.bidirectional_reflectance_distribution_functions,
        R"-x-(
optical_thicknesses:
")-x-"sv,
        v.optical_thicknesses,
        R"-x-(
single_scattering_albedo:
")-x-"sv,
        v.single_scattering_albedo,
        R"-x-(
fractional_scattering:
")-x-"sv,
        v.fractional_scattering,
        R"-x-(
source_polynomial:
")-x-"sv,
        v.source_polynomial,
        R"-x-(
legendre_coefficients:
")-x-"sv,
        v.legendre_coefficients,
        R"-x-(
positive_boundary_condition:
")-x-"sv,
        v.positive_boundary_condition,
        R"-x-(
negative_boundary_condition:
")-x-"sv,
        v.negative_boundary_condition);
  }
};

template <>
struct xml_io_stream<DisortBDRF> {
  static constexpr std::string_view type_name = "DisortBDRF"sv;

  static void write(std::ostream& os,
                    const DisortBDRF& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, DisortBDRF& x, bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<DisortSettings> {
  static constexpr std::string_view type_name = "DisortSettings"sv;

  static void write(std::ostream& os,
                    const DisortSettings& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   DisortSettings& x,
                   bifstream* pbifs = nullptr);
};

using DisortFlux = disort::fluxes;

template <>
struct xml_io_stream_name<DisortFlux> {
  constexpr static std::string_view name = "DisortFlux"sv;
};

template <>
struct xml_io_stream_aggregate<DisortFlux> {
  constexpr static bool value = true;
};

template <>
struct std::formatter<DisortFlux> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  template <class FmtContext>
  FmtContext::iterator format(const DisortFlux& v, FmtContext& ctx) const {
    auto sep = tags.sep();
    return tags.format(ctx,
                       v.frequency_grid,
                       sep,
                       v.altitude_grid,
                       sep,
                       v.up,
                       sep,
                       v.down_diffuse,
                       sep,
                       v.down_direct);
  }
};

using DisortRadiance = disort::radiances;

template <>
struct xml_io_stream_name<DisortRadiance> {
  constexpr static std::string_view name = "DisortRadiance"sv;
};

template <>
struct xml_io_stream_aggregate<DisortRadiance> {
  constexpr static bool value = true;
};

template <>
struct std::formatter<DisortRadiance> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  template <class FmtContext>
  FmtContext::iterator format(const DisortRadiance& v, FmtContext& ctx) const {
    auto sep = tags.sep();
    return tags.format(ctx,
                       v.frequency_grid,
                       sep,
                       v.altitude_grid,
                       sep,
                       v.azimuth_grid,
                       sep,
                       v.zenith_grid,
                       sep,
                       v.data);
  }
};
