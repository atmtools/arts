/*
 * Copyright (c) 2010-2021 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

/** \file shtns.h
 \brief shtns.h is the definition file for SHTns : include this file in your source code to use SHTns.
**/

#ifdef __cplusplus
	#include <complex>
	typedef std::complex<double> cplx;		///< double precision complex number data type
	typedef std::complex<float> cplx_f;		///< single precision (float) complex number data type
#else
	#include <complex.h>
	typedef complex double cplx;			///< double precision complex number data type
	typedef complex float cplx_f;			///< single precision (float) complex number data type
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/// SHTns interface version (loosely follow versions) allowing simple version checks: (major << 16) | (minor << 8) | patchlevel
/// should be increased at least each time this file changes.
#define SHTNS_INTERFACE 0x30500

/// pointer to data structure describing an SHT, returned by shtns_init() or shtns_create().
typedef struct shtns_info* shtns_cfg;

/// pointer to data structure describing a rotation, returned by shtns_rotation_create().
typedef struct shtns_rot_* shtns_rot;

/// different Spherical Harmonic normalizations.
/// see also section \ref norm for details.
enum shtns_norm {
	sht_orthonormal,	///< orthonormalized spherical harmonics (default).
	sht_fourpi,			///< Geodesy and spectral analysis : 4.pi normalization.
	sht_schmidt,		///< Schmidt semi-normalized : 4.pi/(2l+1)
	sht_for_rotations	///< \internal normalization to generate wigner-d matrices internally, with otherwise limited functionality.
};
#define SHT_NO_CS_PHASE (256*4)		///< don't include Condon-Shortley phase (add to last argument of \ref shtns_create)
#define SHT_REAL_NORM (256*8)		///< use a "real" normalization. (add to last argument of \ref shtns_create)

/// different SHT types and algorithms
enum shtns_type {
	sht_gauss,		///< use <b>gaussian grid</b> and quadrature. Allows on-the-fly or matrix-based algorithms.
	sht_auto,		///< use fastest algorithm and grid. This currently implies using the Gauss grid and is thus the same as sht_gauss.
	sht_reg_fast,	///< quick initialization of a <b>regular grid</b>. The grid is the same as with sht_reg_dct
	sht_reg_dct,	///< slow initialization of a regular grid (self-tuning). The grid is equispaced and avoids the poles (Féjer quadrature).
	sht_quick_init, ///< gauss grid, with minimum initialization time (useful for pre/post-processing)
	sht_reg_poles,	///< quick initialization of a <b>regular grid including poles</b> (Clenshaw-Curtis quadrature). Useful for vizualisation.
	sht_gauss_fly	///< legendre polynomials are recomputed on-the-fly for each transform (may be faster on some machines, saves memory and bandwidth).
};
#define SHT_NATIVE_LAYOUT 0			///< Tells shtns_init to use \ref native
#define SHT_THETA_CONTIGUOUS 256	///< use \ref theta_fast
#define SHT_PHI_CONTIGUOUS (256*2)	///< use \ref phi_fast
#define SHT_SOUTH_POLE_FIRST (256*32)	///< latitudinal data are stored starting from south pole.

#define SHT_SCALAR_ONLY (256*16)	///< don't compute vector matrices. (add to flags in shtns_set_grid)
#define SHT_LOAD_SAVE_CFG (256*64)	///< try to load and save the config. (add to flags in shtns_set_grid)
#define SHT_ALLOW_GPU (256*128)		///< allows to use a GPU. This needs special care because the same plan cannot be used simultaneously by different threads anymore.
#define SHT_ALLOW_PADDING (256*256)	///< allows SHTns to add extra space between lines of the spatial array to avoid cache bank conflicts.


#ifndef SHTNS_PRIVATE
struct shtns_info {		// allow read-only access to some data (useful for optimization and helper macros)
	const unsigned int nlm;			///< total number of (l,m) spherical harmonics components.
	const unsigned short lmax;		///< maximum degree (lmax) of spherical harmonics.
	const unsigned short mmax;		///< maximum order (mmax*mres) of spherical harmonics.
	const unsigned short mres;		///< the periodicity along the phi axis.
	const unsigned short nlat_2;	///< half of spatial points in Theta direction (using (shtns.nlat+1)/2 allows odd shtns.nlat.)
	const unsigned int nlat;		///< number of spatial points in Theta direction (latitude) ...
	const unsigned int nphi;		///< number of spatial points in Phi direction (longitude)
	const unsigned int nspat;		///< number of real numbers that must be allocated in a spatial field.
	const unsigned short *const li;	///< degree l for given mode index (size nlm) : li[lm]
	const unsigned short *const mi;	///< order m for given mode index (size nlm) : mi[lm]
	const double *const ct;			///< cos(theta) array (size nlat)
	const double *const st;			///< sin(theta) array (size nlat)
	const unsigned int nlat_padded;	///< number of spatial points in Theta direction, including padding.
	const unsigned int nlm_cplx;	///< number of complex coefficients to represent a complex-valued spatial field.
};
#endif

// MACROS //

/*! \name Access to spherical harmonic components
 * The following macros give access to single spherical harmonic coefficient or perform loops spanning all of them.
**/
///@{
///LiM(shtns, l,im) : macro returning array index for given l and im, corresponding to config shtns.
#define LiM(shtns, l,im) ( (((im)*(2*shtns->lmax + 2 - ((im)+1)*shtns->mres))>>1) + (l) )
/// LM(shtns, l,m) : macro returning array index for given l and m, corresponding to config shtns.
#define LM(shtns, l,m) ( (((((unsigned short)(m))/shtns->mres)*(2*shtns->lmax + 2 - ((m)+shtns->mres)))>>1) + (l) )
/// LM_LOOP( shtns, action ) : macro that performs "action" for every (l,m), with lm set, but neither l, m nor im.
/// \c lm must be a declared int and is the loop counter and the SH array index. more info : \ref spec_data
#define LM_LOOP( shtns, action ) { int lm=0; do { action } while(++lm < shtns->nlm); }
/// LM_L_LOOP : loop over all (l,im) and perform "action"  : l and lm are defined (but NOT m and im).
/// \c lm and \c m must be declared int's. \c lm is the loop counter and SH array index, while \c l is the SH degree. more info : \ref spec_data
#define LM_L_LOOP( shtns, action ) { int lm=0; do { int l=shtns->li[lm]; action } while(++lm < shtns->nlm); }
/// LM_cplx(shtns, l,m) : macro returning array index for SH representation of complex-valued fields for given l and m. (mres=1 only!)
#define LM_cplx(shtns, l, m) ( ((l) <= shtns->mmax) ? (l)*((l)+1)+(m) : shtns->mmax*(2*(l) - shtns->mmax) + (l)+(m) )
///@}


/// total number of 'doubles' required for a spatial field (includes FFTW reserved space).
/// only the first shtns.nlat*shtns.nphi are real spatial data, the remaining is used by the Fourier Transform. more info : \ref spat
#define NSPAT_ALLOC(shtns) (shtns->nspat)

// HELPER MACROS //

/// phi angle value in degrees for given index ip.
#define PHI_DEG(shtns, ip) (360./(shtns->nphi*shtns->mres))*(ip)
/// phi angle value in radians for given index ip.
#define PHI_RAD(shtns, ip) (2.*M_PI/(shtns->nphi*shtns->mres))*(ip)


// FUNCTIONS //

/// compute number of spherical harmonics modes (l,m) to represent a real-valued spatial field for given size parameters. Does not require any previous setup.
long nlm_calc(long lmax, long mmax, long mres);

/// compute number of spherical harmonics modes (l,m) to represent a complex-valued spatial field, for given size parameters. Does not require any previous setup.
long nlm_cplx_calc(long lmax, long mmax, long mres);

void shtns_verbose(int);			///< controls output during initialization: 0=no output (default), 1=some output, 2=degug (if compiled in)
void shtns_print_version(void);		///< print version information to stdout.
const char* shtns_get_build_info(void);	///< get version and build information as a string (the same as the one printed by shtns_print_version())

#ifndef SWIG

void shtns_print_cfg(shtns_cfg);	///< print information about given config to stdout.

/// \name initialization
///@{
/// Simple initialization of the spherical harmonic transforms of given size. Calls \ref shtns_create and \ref shtns_set_grid_auto.
shtns_cfg shtns_init(enum shtns_type flags, int lmax, int mmax, int mres, int nlat, int nphi);
/// Defines the sizes of the spectral description. Use for advanced initialization.
shtns_cfg shtns_create(int lmax, int mmax, int mres, enum shtns_norm norm);
/// Precompute everything for a given spatial grid. Use for advanced initialization, after \ref shtns_create.
int shtns_set_grid(shtns_cfg, enum shtns_type flags, double eps, int nlat, int nphi);
/// Precompute everything and choose the optimal nlat and nphi for a given non-linear order.
int shtns_set_grid_auto(shtns_cfg, enum shtns_type flags, double eps, int nl_order, int *nlat, int *nphi);
/// Copy a given config but allow a different (smaller) mmax and the possibility to enable/disable fft.
shtns_cfg shtns_create_with_grid(shtns_cfg, int mmax, int nofft);
/// Enables multi-thread transform using OpenMP with num_threads (if available). Returns number of threads that will be used.
int shtns_use_threads(int num_threads);
/// Selects the gpu device (device_id % Num_devices). Must be called BEFORE any initialization. Internally calls cudaSetDevice(). Returns the actual device or -1 when no device found.
int shtns_use_gpu(int device_id);

void shtns_reset(void);				///< destroy all configs, free memory, and go back to initial state.
void shtns_destroy(shtns_cfg);		///< free memory of given config, which cannot be used afterwards.
void shtns_unset_grid(shtns_cfg);	///< unset the grid.

/// set the use of Robert form. If robert != 0, the vector synthesis returns a field multiplied by sin(theta), while the analysis divides by sin(theta) before the transform.
void shtns_robert_form(shtns_cfg, int robert);

void* shtns_malloc(size_t bytes);	///< alloc appropriate memory (pinned for gpu, aligned for avx, ...). Use \ref shtns_free to free it.
void shtns_free(void* p);			///< free memory allocated with \ref shtns_malloc


///@}

/// \name special values
///@{
double sh00_1(shtns_cfg);	///< return the spherical harmonic representation of 1 (l=0,m=0)
double sh10_ct(shtns_cfg);	///< return the spherical harmonic representation of cos(theta) (l=1,m=0)
double sh11_st(shtns_cfg);	///< return the spherical harmonic representation of sin(theta)*cos(phi) (l=1,m=1)
double shlm_e1(shtns_cfg, int l, int m);		///< return the l,m SH coefficient corresponding to unit energy.
/// fill the given array with Gauss weights. returns the number of weights written (0 if not a Gauss grid).
int shtns_gauss_wts(shtns_cfg, double *wts);

///@}

/// \name Rotation functions
///@{
/// Rotate a SH representation Qlm around the z-axis by angle alpha (in radians),
/// which is the same as rotating the reference frame by angle -alpha.
/// Result is stored in Rlm (which can be the same array as Qlm).
void SH_Zrotate(shtns_cfg, cplx *Qlm, double alpha, cplx *Rlm);
/// \deprecated Rotate SH representation around Y axis by alpha (in radians).
void SH_Yrotate(shtns_cfg, cplx *Qlm, double alpha, cplx *Rlm);
/// \deprecated Rotate SH representation around Y axis by 90 degrees.
void SH_Yrotate90(shtns_cfg, cplx *Qlm, cplx *Rlm);
/// \deprecated Rotate SH representation around X axis by 90 degrees.
void SH_Xrotate90(shtns_cfg, cplx *Qlm, cplx *Rlm);


/// Creation of a spherical harmonic rotation object for degrees and orders up to lmax and mmax respectively.
shtns_rot shtns_rotation_create(const int lmax, const int mmax, int norm);
/// Release memory allocated for the rotation object.
void shtns_rotation_destroy(shtns_rot r);
/// Defines a rotation parameterized by the 3 intrinsic Euler angles with ZYZ convention.
void shtns_rotation_set_angles_ZYZ(shtns_rot r, double alpha, double beta, double gamma);
/// Defines a rotation parameterized by the 3 intrinsic Euler angles with ZXZ convention.
void shtns_rotation_set_angles_ZXZ(shtns_rot r, double alpha, double beta, double gamma);
/// Defines a rotation of angle theta around axis of cartesian coordinates (Vx,Vy,Vz).
void shtns_rotation_set_angle_axis(shtns_rot r, double theta, double Vx, double Vy, double Vz);
/// Generate spherical-harmonic rotation matrix for given degree l and orthonormal convention (Wigner-d matrix)
/// \param[out] mx is an (2*l+1)*(2*l+1) array that will be filled with the Wigner-d matrix elements (rotation matrix along Y-axis in orthonormal spherical harmonic space).
/// \return 0 if error, or 2*l+1 (size of the square matrix) otherwise.
int shtns_rotation_wigner_d_matrix(shtns_rot r, const int l, double* mx);
/// apply rotation to the **orthonormal** spherical harmonic expansion Zlm of a complex-valued field, store the result into Rlm (can be the same as Zlm)
void shtns_rotation_apply_cplx(shtns_rot r, cplx* Zlm, cplx* Rlm);
/// apply rotation to the **orthonormal** spherical harmonic expansion Qlm of a real field, store the result into Rlm (can be the same as Qlm)
void shtns_rotation_apply_real(shtns_rot r, cplx* Qlm, cplx* Rlm);

///@}

/// \name Generation of Legendre associated functions
///@{
/// Compute values of legendre polynomials noramalized for spherical harmonics,
/// for a range of l=m..lmax, at given m and x, using stable recurrence.
/// Requires a previous call to \ref shtns_create().
/// Output compatible with the GSL function gsl_sf_legendre_sphPlm_array(lmax, m, x, yl)
/// \param[in] lmax maximum degree computed
/// \param[in] im = m/MRES with m the SH order
/// \param[in] x argument, x=cos(theta).
/// \param[out] yl is a double array of size (lmax-m+1) filled with the values (of increasing degree l).
int legendre_sphPlm_array(shtns_cfg shtns, const int lmax, const int im, const double x, double *yl);
int legendre_sphPlm_deriv_array(shtns_cfg shtns, const int lmax, const int im, const double x, const double sint, double *yl, double *dyl);
///@}

/// \name Special operator functions
///@{
/// compute the matrix (stored in mx, a double array of size 2*NLM) required
/// to multiply an SH representation by cos(theta) using \ref SH_mul_mx.
void mul_ct_matrix(shtns_cfg, double* mx);
/// compute the matrix (stored in mx, a double array of size 2*NLM) required
/// to apply sin(theta)*d/dtheta to an SH representation using \ref SH_mul_mx.
void st_dt_matrix(shtns_cfg, double* mx);
/// Apply a matrix involving l+1 and l-1 to an SH representation Qlm. Result stored in Rlm (must be different from Qlm).
void SH_mul_mx(shtns_cfg, double* mx, cplx *Qlm, cplx *Rlm);
///@}

/** \addtogroup sht Spherical Harmonic transform functions.
 * All these function perform a global spherical harmonic transform.
 * Their first argument is a shtns_cfg variable (which is a pointer to a \ref shtns_info struct)
 * obtained by a previous call to \ref shtns_create or \ref shtns_init.
 * \see \ref spat \see \ref spec \see \ref vsh
 */
///@{

/// \name Scalar transforms
///@{
/// transform the scalar field Vr into its spherical harmonic representation Qlm.
/// \param[in] shtns = a configuration created by \ref shtns_create with a grid set by \ref shtns_set_grid or \ref shtns_set_grid_auto
/// \param[in] Vr = spatial scalar field : double array of size shtns->nspat; NOT GUARANTEED TO BE PRESERVED.
/// \param[out] Qlm = spherical harmonics coefficients : cplx array of size shtns->nlm.
void spat_to_SH(shtns_cfg shtns, double *Vr, cplx *Qlm);
/// transform the spherical harmonic coefficients Qlm into its spatial representation Vr.
/// \param[in] shtns = a configuration created by \ref shtns_create with a grid set by \ref shtns_set_grid or \ref shtns_set_grid_auto
/// \param[in] Qlm = spherical harmonics coefficients : cplx array of size shtns->nlm.
/// \param[out] Vr = spatial scalar field : double array of size shtns->nspat.
void SH_to_spat(shtns_cfg shtns, cplx *Qlm, double *Vr);
/// complex scalar synthesis.
/// \param[in] shtns = a configuration created by \ref shtns_create with a grid set by \ref shtns_set_grid or \ref shtns_set_grid_auto
/// \param[in] alm[l*(l+1)+m] is the SH coefficient of order l and degree m (with -l <= m <= l) [total of (LMAX+1)^2 coefficients]
/// \param[out] z = complex spatial field
void SH_to_spat_cplx(shtns_cfg shtns, cplx *alm, cplx *z);
/// complex scalar analysis.
/// \param[in] shtns = a configuration created by \ref shtns_create with a grid set by \ref shtns_set_grid or \ref shtns_set_grid_auto
/// \param[in] z = complex spatial field
/// \param[out] alm[l*(l+1)+m] is the SH coefficient of order l and degree m (with -l <= m <= l) [total of (LMAX+1)^2 coefficients]
void spat_cplx_to_SH(shtns_cfg shtns, cplx *z, cplx *alm);
///@}

/// \name 2D vector transforms
///@{
/// transform the theta and phi components (Vt,Vp) of a vector into its spheroidal-toroidal spherical harmonic representation (Slm,Tlm). \see \ref vsh
void spat_to_SHsphtor(shtns_cfg, double *Vt, double *Vp, cplx *Slm, cplx *Tlm);
/// transform spheroidal-toroidal spherical harmonic coefficients (Slm,Tlm) to the spatial theta and phi components (Vt,Vp). \see \ref vsh
void SHsphtor_to_spat(shtns_cfg, cplx *Slm, cplx *Tlm, double *Vt, double *Vp);
/// transform spheroidal spherical harmonic coefficients Slm to the spatial theta and phi components (Vt,Vp), effectively computing the gradient of S. \see \ref vsh
void SHsph_to_spat(shtns_cfg, cplx *Slm, double *Vt, double *Vp);
/// transform toroidal spherical harmonic coefficients Tlm to the spatial theta and phi components (Vt,Vp). \see \ref vsh
void SHtor_to_spat(shtns_cfg, cplx *Tlm, double *Vt, double *Vp);

/// transform the theta and phi components (Vt,Vp) of a complex valued vector into its spheroidal-toroidal spherical harmonic representation (Slm,Tlm). \see \ref vsh
void spat_cplx_to_SHsphtor(shtns_cfg, cplx *Vt, cplx *Vp, cplx *Slm, cplx *Tlm);
/// transform spheroidal-toroidal spherical harmonic coefficients (Slm,Tlm) to the spatial theta and phi complex-valued components (Vt,Vp). \see \ref vsh
void SHsphtor_to_spat_cplx(shtns_cfg, cplx *Slm, cplx *Tlm, cplx *Vt, cplx *Vp);
///@}
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat
#define SH_to_grad_spat(shtns, S,Gt,Gp) SHsph_to_spat(shtns, S, Gt, Gp)

/// \name 3D transforms (combine scalar and vector)
///@{
/// 3D vector transform from spherical coordinates to radial-spheroidal-toroidal spectral components (see \ref vsh_def).
/// They should be prefered over separate calls to scalar and 2D vector transforms as they can be significantly faster.
void spat_to_SHqst(shtns_cfg, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm);
/// 3D vector transform from radial-spheroidal-toroidal spectral components (see \ref vsh_def) to spherical coordinates.
/// They should be prefered over separate calls to scalar and 2D vector transforms as they can be significantly faster.
void SHqst_to_spat(shtns_cfg, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp);

void spat_cplx_to_SHqst(shtns_cfg, cplx *Vr, cplx *Vt, cplx *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm);
void SHqst_to_spat_cplx(shtns_cfg, cplx *Qlm, cplx *Slm, cplx *Tlm, cplx *Vr, cplx *Vt, cplx *Vp);
///@}

/// \name Truncated transforms at given degree l
/// wiht l <= lmax used for setup.
///@{
void spat_to_SH_l(shtns_cfg, double *Vr, cplx *Qlm, int ltr);
void SH_to_spat_l(shtns_cfg, cplx *Qlm, double *Vr, int ltr);

void SHsphtor_to_spat_l(shtns_cfg, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, int ltr);
void SHsph_to_spat_l(shtns_cfg, cplx *Slm, double *Vt, double *Vp, int ltr);
void SHtor_to_spat_l(shtns_cfg, cplx *Tlm, double *Vt, double *Vp, int ltr);
void spat_to_SHsphtor_l(shtns_cfg, double *Vt, double *Vp, cplx *Slm, cplx *Tlm, int ltr);

void spat_to_SHqst_l(shtns_cfg, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr);
void SHqst_to_spat_l(shtns_cfg, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int ltr);
///@}
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat_l
#define SH_to_grad_spat_l(shtns, S,Gt,Gp,ltr) SHsph_to_spat_l(shtns, S, Gt, Gp, ltr)

/// \name Legendre transform at given m (no fft) and truncated at given degree l <= lmax
/// The input and output arrays contain only the specified m=im*mres, that is spatial size is nlat
/// and spectral size is lmax+1-m, containing only the coefficients of the given m.
///@{
void spat_to_SH_ml(shtns_cfg, int im, cplx *Vr, cplx *Ql, int ltr);
void SH_to_spat_ml(shtns_cfg, int im, cplx *Ql, cplx *Vr, int ltr);

void spat_to_SHsphtor_ml(shtns_cfg, int im, cplx *Vt, cplx *Vp, cplx *Sl, cplx *Tl, int ltr);
void SHsphtor_to_spat_ml(shtns_cfg, int im, cplx *Sl, cplx *Tl, cplx *Vt, cplx *Vp, int ltr);
void SHsph_to_spat_ml(shtns_cfg, int im, cplx *Sl, cplx *Vt, cplx *Vp, int ltr);
void SHtor_to_spat_ml(shtns_cfg, int im, cplx *Tl, cplx *Vt, cplx *Vp, int ltr);

void spat_to_SHqst_ml(shtns_cfg, int im, cplx *Vr, cplx *Vt, cplx *Vp, cplx *Ql, cplx *Sl, cplx *Tl, int ltr);
void SHqst_to_spat_ml(shtns_cfg, int im, cplx *Ql, cplx *Sl, cplx *Tl, cplx *Vr, cplx *Vt, cplx *Vp, int ltr);
///@}
/// Compute the spatial representation of the gradient of a scalar SH field. Alias for \ref SHsph_to_spat_l
#define SH_to_grad_spat_ml(shtns, im, S,Gt,Gp,ltr) SHsph_to_spat_ml(shtns, im, S, Gt, Gp, ltr)

///@}

/// \name Local and partial evalutions of a SH representation :
/// Does not require a call to \ref shtns_set_grid_auto
///@{
double SH_to_point(shtns_cfg, cplx *Qlm, double cost, double phi);
cplx SH_to_point_cplx(shtns_cfg, cplx *alm, double cost, double phi);
void SH_to_grad_point(shtns_cfg, cplx *DrSlm, cplx *Slm,
					double cost, double phi, double *vr, double *vt, double *vp);
void SHqst_to_point(shtns_cfg, cplx *Qlm, cplx *Slm, cplx *Tlm,
					double cost, double phi, double *vr, double *vt, double *vp);

void SH_to_lat(shtns_cfg shtns, cplx *Qlm, double cost,
					double *vr, int nphi, int ltr, int mtr);
void SHqst_to_lat(shtns_cfg, cplx *Qlm, cplx *Slm, cplx *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr);
///@}


#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */
