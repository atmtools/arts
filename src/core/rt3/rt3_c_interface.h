#pragma once

#include <cstdint>

/* The ISO_C_BINDING entry points of 3rdparty/polradtran/rt3_c_interface.f90,
   available when ENABLE_RT3=ON (ARTS_HAS_RT3).

   Each wraps the Fortran routine of the same name (without the rt3_ and
   RT3_ prefixes).  Integers are int64_t, the one-character flags, the
   LOGICAL and the scalars are passed by value, a COMPLEX*16 as its real and
   imaginary parts.  Arrays are the Fortran column-major arrays: a Fortran
   A(d1, ..., dk) is the row-major [dk, ..., d1].  Arrays that the Fortran
   routine only reads are const here.

   None of the Fortran code is reentrant (COMMON blocks, SAVEd FFT tables
   and static local arrays): the caller must serialise all calls. */
extern "C" {
// radtran3.f
void rt3_radtran(std::int64_t        nstokes,
                 std::int64_t        nummu,
                 std::int64_t        aziorder,
                 double              max_delta_tau,
                 std::int64_t        src_code,
                 char                quad_type,
                 char                deltam,
                 double              direct_flux,
                 double              direct_mu,
                 double              ground_temp,
                 char                ground_type,
                 double              ground_albedo,
                 double              ground_index_re,
                 double              ground_index_im,
                 double              sky_temp,
                 double              wavelength,
                 std::int64_t        num_layers,
                 double*             height,
                 double*             temperatures,
                 double*             gas_extinct,
                 std::int64_t        nsl,
                 double*             scat_extinct,
                 double*             scat_scatter,
                 const std::int64_t* scat_nlegen,
                 std::int64_t        ldcoef,
                 double*             scat_coef,
                 const std::int64_t* scatlayers,
                 std::int64_t        noutlevels,
                 const std::int64_t* outlevels,
                 double*             mu_values,
                 double*             up_flux,
                 double*             down_flux,
                 double*             up_rad,
                 double*             down_rad);

// radscat3.f
void rt3_get_scat_set(char          deltam,
                      std::int64_t  nummu,
                      std::int64_t  nlegin,
                      const double* coefin,
                      double        extin,
                      double        scatin,
                      std::int64_t* nlegen,
                      double*       coef,
                      double*       extinction,
                      double*       scatter);
void rt3_scattering(std::int64_t  nummu,
                    std::int64_t  aziorder,
                    std::int64_t  nstokes,
                    const double* mu_values,
                    const double* quad_weights,
                    std::int64_t  numlegendre,
                    const double* legendre_coef,
                    std::int64_t  scat_num,
                    double*       scatbuf);
void rt3_direct_scattering(std::int64_t  nummu,
                           std::int64_t  aziorder,
                           std::int64_t  nstokes,
                           const double* mu_values,
                           std::int64_t  numlegendre,
                           const double* legendre_coef,
                           double        direct_mu,
                           std::int64_t  scat_num,
                           double*       directbuf);
void rt3_get_scattering(std::int64_t  nstokes,
                        std::int64_t  nummu,
                        std::int64_t  mode,
                        std::int64_t  aziorder,
                        std::int64_t  scat_num,
                        const double* scatbuf,
                        double*       scatter_matrix);
//! scat is (nstokes, nummu, nstokes, nummu, 4): makes parts 3 and 4 from 2
//! and 1
void rt3_scatter_symmetry(std::int64_t nstokes, std::int64_t nummu, double* scat);
//! STOPs if the phase function is not normalised
void rt3_check_norm(std::int64_t nstokes, std::int64_t nummu, const double* quad_weights, const double* scatter_matrix);
void rt3_get_direct(std::int64_t  nstokes,
                    std::int64_t  nummu,
                    std::int64_t  mode,
                    std::int64_t  aziorder,
                    std::int64_t  scat_num,
                    const double* directbuf,
                    double*       direct_vector);
//! dosum[6] is output
void rt3_number_sums(std::int64_t nstokes, std::int64_t nlegen, const double* coef, std::int64_t* dosum);
void rt3_sum_legendre(
    std::int64_t nlegen, const double* coef, double x, const std::int64_t* dosum, double* phase_matrix);
void rt3_rotate_phase_matrix(const double* phase_matrix1,
                             double        mu1,
                             double        mu2,
                             double        delphi,
                             double        cos_scat,
                             double*       phase_matrix2,
                             std::int64_t  nstokes);
//! matrix1 and matrix2 may be the same array
void rt3_matrix_symmetry(std::int64_t nstokes, const double* matrix1, double* matrix2);
void rt3_fourier_matrix(
    std::int64_t aziorder, std::int64_t numpts, std::int64_t nstokes, const double* real_matrix, double* basis_matrix);
//! Both vectors are input or output by direction; real_vector is
//! overwritten when order > 0 (in-place FFT)
void rt3_fourier_basis(std::int64_t numbasis,
                       std::int64_t order,
                       std::int64_t numpts,
                       std::int64_t direction,
                       double*      basis_vector,
                       double*      real_vector);
//! In place; isign +1 real to complex conjugate, -1 back; STOPs above
//! n = 512
void rt3_fft1dr(double* data, std::int64_t n, std::int64_t isign);
//! Complex FFT of n points in place (data holds 2 n values)
void rt3_fftc(double* data, std::int64_t n, const double* phase);
//! nyquist[2] is output for isign > 0, nyquist[0] input otherwise
void rt3_fixreal(double* data, double* nyquist, std::int64_t n, std::int64_t isign, const double* phase);
//! phase holds 4 nmax values
void rt3_makephase(double* phase, std::int64_t nmax);
void rt3_combine_phase_modes(std::int64_t  nstokes,
                             std::int64_t  aziorder,
                             std::int64_t  m,
                             double        tmp,
                             const double* basis_matrix,
                             double*       out_matrix);

// radintg3.f
void rt3_initialize(std::int64_t  nstokes,
                    std::int64_t  nummu,
                    std::int64_t  n,
                    double        delta_z,
                    const double* mu_values,
                    double        extinction,
                    double        albedo,
                    const double* phase_function,
                    double*       reflect,
                    double*       trans);
void rt3_initial_source(std::int64_t  nstokes,
                        std::int64_t  nummu,
                        std::int64_t  n,
                        double        delta_z,
                        const double* mu_values,
                        double        extinction,
                        const double* source_vector,
                        double*       source);
void rt3_nonscatter_layer(std::int64_t  nstokes,
                          std::int64_t  nummu,
                          std::int64_t  mode,
                          double        deltatau,
                          const double* mu_values,
                          double        planck0,
                          double        planck1,
                          double*       reflect,
                          double*       trans,
                          double*       source);
//! Also overwrites reflect, trans, exp_source and lin_source (doubling scratch)
void rt3_doubling_integration(std::int64_t n,
                              std::int64_t num_doubles,
                              std::int64_t src_code,
                              bool         symmetric,
                              double*      reflect,
                              double*      trans,
                              double*      exp_source,
                              double       expfactor,
                              double*      lin_source,
                              double       linfactor,
                              double*      t_reflect,
                              double*      t_trans,
                              double*      t_source);
void rt3_combine_layers(std::int64_t  n,
                        const double* reflect1,
                        const double* trans1,
                        const double* source1,
                        const double* reflect2,
                        const double* trans2,
                        const double* source2,
                        double*       out_reflect,
                        double*       out_trans,
                        double*       out_source);
void rt3_internal_radiance(std::int64_t  n,
                           const double* upreflect,
                           const double* uptrans,
                           const double* upsource,
                           const double* downreflect,
                           const double* downtrans,
                           const double* downsource,
                           const double* intoprad,
                           const double* inbottomrad,
                           double*       uprad,
                           double*       downrad);

// radutil3.f
void rt3_lambert_surface(std::int64_t  nstokes,
                         std::int64_t  nummu,
                         std::int64_t  mode,
                         const double* mu_values,
                         const double* quad_weights,
                         double        ground_albedo,
                         double*       reflect,
                         double*       trans,
                         double*       source);
void rt3_lambert_radiance(std::int64_t nstokes,
                          std::int64_t nummu,
                          std::int64_t mode,
                          std::int64_t src_code,
                          double       ground_albedo,
                          double       ground_temp,
                          double       wavelength,
                          double       direct_sfc_flux,
                          double*      radiance);
void rt3_fresnel_surface(std::int64_t  nstokes,
                         std::int64_t  nummu,
                         const double* mu_values,
                         double        index_re,
                         double        index_im,
                         double*       reflect,
                         double*       trans,
                         double*       source);
void rt3_fresnel_radiance(std::int64_t  nstokes,
                          std::int64_t  nummu,
                          std::int64_t  mode,
                          const double* mu_values,
                          double        index_re,
                          double        index_im,
                          double        ground_temp,
                          double        wavelength,
                          double*       radiance);
void rt3_thermal_radiance(std::int64_t nstokes,
                          std::int64_t nummu,
                          std::int64_t mode,
                          double       temperature,
                          double       albedo,
                          double       wavelength,
                          double*      radiance);
void rt3_double_gauss_quadrature(std::int64_t num, double* abscissas, double* weights);
void rt3_gauss_legendre_quadrature(std::int64_t num, double* abscissas, double* weights);
void rt3_lobatto_quadrature(std::int64_t num, double* abscissas, double* weights);
}
