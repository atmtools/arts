#pragma once

#include <cstdint>

/* The ISO_C_BINDING entry points of 3rdparty/polradtran/rt4_c_interface.f90,
   available when ENABLE_RT4=ON (ARTS_HAS_RT4).  ARTS's RT4 is the C++
   port; these are the Fortran routines it ports, for the tests that
   compare the two (cpp.fast.rt4-radtrano-test).

   Each wraps the Fortran routine of the same name (without the rt4_
   prefix).  Integers are int64_t, the one-character flags and the scalars
   are passed by value, a COMPLEX*16 as its real and imaginary parts.
   Arrays are the Fortran column-major arrays: a Fortran A(d1, ..., dk) is
   the row-major [dk, ..., d1].  Arrays that the Fortran routine only reads
   are const here.

   None of the Fortran code is reentrant (COMMON blocks and static local
   arrays): the caller must serialise all calls. */
extern "C" {
// radtran4.f
void rt4_radtrano(std::int64_t nstokes,
                  std::int64_t nummu,
                  std::int64_t nuummu,
                  double       max_delta_tau,
                  char         quad_type,
                  double       ground_temp,
                  char         ground_type,
                  double       ground_albedo,
                  double       ground_index_re,
                  double       ground_index_im,
                  double*      ground_reflec,
                  double*      surf_reflect,
                  double*      gnd_radiance,
                  double       sky_temp,
                  double       wavelength,
                  std::int64_t num_layers,
                  double*      height,
                  double*      temperatures,
                  double*      gas_extinct,
                  std::int64_t nsl,
                  double*      scatlayers,
                  double*      extinct_matrix,
                  double*      emis_vector,
                  double*      scatter_matrix,
                  double*      mu_values,
                  double*      up_rad,
                  double*      down_rad);

// radintg4.f
void rt4_initialize(std::int64_t  nstokes,
                    std::int64_t  nummu,
                    double        delta_z,
                    const double* mu_values,
                    const double* quad_weights,
                    double        gas_extinct,
                    const double* extinct_matrix,
                    const double* scatter_matrix,
                    double*       reflect,
                    double*       trans);
void rt4_initial_source(std::int64_t  nstokes,
                        std::int64_t  nummu,
                        double        delta_z,
                        const double* mu_values,
                        double        planck,
                        const double* emis_vector,
                        double        gas_extinct,
                        double*       source);
void rt4_nonscatter_layer(std::int64_t  nstokes,
                          std::int64_t  nummu,
                          double        deltatau,
                          const double* mu_values,
                          double        planck0,
                          double        planck1,
                          double*       reflect,
                          double*       trans,
                          double*       source);
void rt4_internal_radiance(std::int64_t  n,
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
//! Also overwrites reflect, trans and lin_source (doubling scratch)
void rt4_doubling_integration(std::int64_t n,
                              std::int64_t num_doubles,
                              bool         symmetric,
                              double*      reflect,
                              double*      trans,
                              double*      lin_source,
                              double       linfactor,
                              double*      t_reflect,
                              double*      t_trans,
                              double*      t_source);
void rt4_combine_layers(std::int64_t  n,
                        const double* reflect1,
                        const double* trans1,
                        const double* source1,
                        const double* reflect2,
                        const double* trans2,
                        const double* source2,
                        double*       out_reflect,
                        double*       out_trans,
                        double*       out_source);

// radutil4.f
void rt4_lambert_surface(std::int64_t  nstokes,
                         std::int64_t  nummu,
                         std::int64_t  mode,
                         const double* mu_values,
                         const double* quad_weights,
                         double        ground_albedo,
                         double*       reflect,
                         double*       trans,
                         double*       source);
void rt4_lambert_radiance(std::int64_t nstokes,
                          std::int64_t nummu,
                          double       ground_albedo,
                          double       ground_temp,
                          double       wavelength,
                          double*      radiance);
void rt4_fresnel_surface(std::int64_t  nstokes,
                         std::int64_t  nummu,
                         const double* mu_values,
                         double        index_re,
                         double        index_im,
                         double*       reflect,
                         double*       trans,
                         double*       source);
void rt4_fresnel_radiance(std::int64_t  nstokes,
                          std::int64_t  nummu,
                          const double* mu_values,
                          double        index_re,
                          double        index_im,
                          double        ground_temp,
                          double        wavelength,
                          double*       radiance);
void rt4_specular_surface(std::int64_t  nstokes,
                          std::int64_t  nummu,
                          const double* ground_reflec,
                          double*       reflect,
                          double*       trans,
                          double*       source);
void rt4_specular_radiance(std::int64_t  nstokes,
                           std::int64_t  nummu,
                           const double* ground_reflec,
                           double        ground_temp,
                           double        wavelength,
                           double*       radiance);
void rt4_external_surface(std::int64_t  nstokes,
                          std::int64_t  nummu,
                          const double* surf_refl,
                          const double* radiance,
                          double*       reflect,
                          double*       trans,
                          double*       source);
void rt4_thermal_radiance(
    std::int64_t nstokes, std::int64_t nummu, double temperature, double albedo, double wavelength, double* radiance);
//! PLANCK_FUNCTION with radiance units ('R')
double rt4_planck_function(double temp, double wavelength);
void   rt4_double_gauss_quadrature(std::int64_t num, double* abscissas, double* weights);
void   rt4_gauss_legendre_quadrature(std::int64_t num, double* abscissas, double* weights);
void   rt4_lobatto_quadrature(std::int64_t num, double* abscissas, double* weights);

// radmat.f
void rt4_mcopy(std::int64_t n, std::int64_t m, const double* matrix1, double* matrix2);
void rt4_mzero(std::int64_t n, std::int64_t m, double* matrix1);
void rt4_midentity(std::int64_t n, double* matrix);
}
