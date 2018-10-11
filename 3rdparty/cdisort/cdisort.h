/************************************************************************
 * $Id: cdisort.h 2887 2013-03-11 09:19:28Z robert.buras $
 ************************************************************************/

/*
 *   Copyright (c) 2011 by Timothy E. Dowling
 *   
 *   This file is part of cdisort.
 *
 *   cdisort is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   cdisort is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with cdisort.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __cdisort_h
#define __cdisort_h

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

#if HAVE_BRDF
#include "ocean.h"
#include "ambralsfor.h"
#endif

/*------------*
 * Structures *
 *------------*/
/*
 * See DISORT.txt for details on inputs and outputs.
 */
typedef struct {
  /* DISORT CONTROL FLAGS */
  int
    usrtau,    /* TRUE=> radiant quantities returned at user-specified optical depths  */
    usrang,    /* TRUE=> radiant quantities returned at user-specified polar angles    */
    ibcnd,     /* SPECIAL_BC => return only albedo and transmis., see Ref S2           */
    lamber,    /* TRUE=> isotropically reflecting bottom boundary, FALSE=>bi-dir.      */
    planck,    /* TRUE=>incl. thermal emission                                         */
    spher,     /* TRUE=>pseudo-spherical geometry, otherwise plane-parallel            */
    onlyfl,    /* FALSE=>return intensities in addition to other radiant quantities    */
    prnt[5],   /* Print flags: 0=input variables (except pmom), 1=fluxes,              */
               /*   2=intensities, 3=transmis. and albedo, 4=pmom                      */
    brdf_type, /* type of BRDF, can be Hapke, RPV, Cox&Munk, Ambrals                   */
    quiet,     /* quiet output                                                         */
    intensity_correction,      /* apply intensity correction                           */
    old_intensity_correction,  /* use original intensity correction routine            */
    general_source,      /* Include solution for a general user specified source term. */         
    output_uum; /* TRUE=> uum is returned as a seperate output                         */
} disort_flag;

typedef struct {
  /* DISORT OUTPUT RADIANT QUANTITIES */
  double
    rfldir,    /* Direct-beam flux (w/o delta-M scaling),                              */
    rfldn,     /* Diffuse down-flux (tot.-direct-beam; w/o delta-M scaling),           */
    flup,      /* Diffuse up-flux                                                      */
    dfdt,      /* Flux divergence, d(net flux)/d(optical depth)                        */
    uavg,      /* Mean intensity, incl. direct beam (not corr. for delta-M scaling)    */
    uavgdn,    /* Mean diffuse downward intensity, not incl. direct beam               */
	       /* (not corr. for delta-M scaling)                                      */
    uavgup,    /* Mean diffuse downward intensity, not incl. direct beam               */
	       /* (not corr. for delta-M scaling)                                      */
    uavgso;    /* Mean diffuse direct solar, that is the direct beam                   */
	       /* (not corr. for delta-M scaling)                                      */
} disort_radiant;

typedef struct {
  /* DISORT TOP AND BOTTOM BOUNDARY CONDITIONS */
  double
    fbeam,     /* Intensity of incident parallel beam at top boundary                  */
    umu0,      /* Polar angle cosine of incident beam (positive)                       */
    phi0,      /* Azimuth angle of incident beam (0 to 360 deg)                        */
    fisot,     /* Intensity of top-boundary isotropic illumination                     */
    fluor,     /* Intensity of bottom-boundary isotropic illumination                  */
    ttemp,     /* Temperature [K] of top boundary                                      */
    btemp,     /* Temperature [K] of bottom boundary                                   */
    temis,     /* Emissivity of top boundary. Needed if planck = TRUE                  */
    albedo;    /* Albedo of bottom boundary, needed if lamber = TRUE                   */
} disort_bc;

typedef struct {
  /* RPV BRDF specifications */
  double
    rho0,     /* BRDF rpv: rho0       */
    k,        /* BRDF rpv: k          */
    theta,    /* BRDF rpv: theta      */
    sigma,    /* BRDF rpv snow: sigma */
    t1,       /* BRDF rpv snow: t1    */
    t2,       /* BRDF rpv snow: t2    */
    scale;    /* BRDF rpv: scale      */
} rpv_brdf_spec;

typedef struct {
  /* Ambrals BRDF specifications */
  double
    iso,      /* BRDF ambrals: iso    */
    vol,      /* BRDF ambrals: vol    */
    geo;      /* BRDF ambrals: geo    */
} ambrals_brdf_spec;

typedef struct {
  /* Cox and Munk BRDF specifications */
  double
    u10,      /* BRDF C&M: u10        */
    pcl,      /* BRDF C&M: pcl        */
    xsal;     /* BRDF C&M: sal        */
} cam_brdf_spec;

typedef struct {
  /* brdf types */
  rpv_brdf_spec
    *rpv;      /* specification for rpv BRDF      */
#if HAVE_BRDF
  ambrals_brdf_spec
    *ambrals;  /* specification for ambrals BRDF  */
  cam_brdf_spec
    *cam;       /* specification for Cox&Munk BRDF */
#endif
} disort_brdf;


typedef struct {
  char
    header[128];
  disort_flag
    flag;
  disort_bc
    bc;
  disort_brdf
    brdf;
  int
    nlyr,      /* Number of computational layers                                       */
    nmom,      /* Number of phase function moments (not including the zeroth)          */
    nstr,      /* Number of streams (computational polar angles). Even and >= 2        */
    nmom_nstr, /* IMAX(ds.nmom,ds.nstr), used to set size of 1st dimension of PMOM     */
    ntau,      /* Number of computational optical depths                               */
    numu,      /* Number of computational polar angles                                 */
    nphi,      /* Number of azimuthal angles at which to return intensities            */
    nphase;    /* number of angles (grid points)                                       */
  double
    wvnmlo,    /* Wavenumber [cm^-1] lower range, used for Planck function             */
    wvnmhi,    /* Wavenumber [cm^-1] upper range, used for Planck function             */
    accur,     /* Convergence criteria for azimuthal (Fourier cosine) series           */
    radius;    /* Radius of the body of interest. Only used if spher=TRUE. Must be     */
               /* the same units as zd.                                                */
  double
    *dtauc,     /* Optical depths of computational layers, DTAUC(lc)                   */
    *ssalb,     /* Single-scatter albedos of computational layers, SSALB(lc)           */
    *pmom,      /* Coefficients (moments) of Legendre polynomials, PMOM(k,lc)          */
    *temper,    /* Temperatures [K] of levels, TEMPER(lev)                             */
    *utau,      /* Optical depths, in increasing order, on layer boundaries, UTAU(lu)  */
    *umu,       /* Cosines of polar angles, increasing order, UMU(iu)                  */
    *phi,       /* Azimuthal output angles [deg], used when onlyfl = FALSE             */
    *zd,        /* The altitude of levels, used when flag.spher = TRUE                 */
    *mu_phase,  /* values of scattering angles for which phase function given          */
    *phase,     /* phase function as a function of scattering angles                   */
    *gensrc,    /* User specified general source at computational angles               */
    *gensrcu;   /* User specified general source at user angles                        */
} disort_state;

typedef struct {
  disort_radiant
   *rad;       /* See typedef disort_radiant                                           */
  double
    *albmed,   /* Albedo of medium, ALBMED(iu) (ds.flag.ibcnd = SPECIAL_BC case only)  */
    *trnmed,   /* Transmissivity of medium, TRNMED(iu) (ds.flag.ibcnd = SPECIAL_BC)    */
    *uu,       /* Intensity, UU(iu,lu,j) (if ds.flag.onlyfl = FALSE; zero otherwise)   */
    *u0u,      /* Az.av.Int. U0U(iu,lu,j) (if ds.flag.onlyfl = FALSE; zero otherwise)  */
    *uum;      /* Intensity, UUM(iu,lu,j) (if ds.flag.output_uum = TRUE; not 
		  allocated space or used otherwise)                                   */
} disort_output;

typedef struct {
  double
    zero,
    one;
} disort_pair;

typedef struct {
  double
    zero,
    one,
    alpha;
} disort_triplet;

/*
 * Definitions specific to twostr()
 */

#define TWOSTR_NERR 22

typedef struct {
  double
    super,
    on,
    sub;
} twostr_diag;

typedef struct {
  double
    xb_0d,  /* x-sub-zero-sup-minus in expansion of pseudo-spherical beam source, Eq. KST(22) */
    xb_0u,  /* x-sub-zero-sup-plus  "                                                         */
    xb_1d,  /* x-sub-one-sup-minus  "                                                         */ 
    xb_1u,  /* x-sub-one-sup-plus   "                                                         */
    xp_0,   /* x-sub-zero in expansion of thermal source func.; Eq. KST(22), has no (mu) dep. */
    xp_1,   /* x-sub-one  "                                                                   */
    yb_0d,  /* y-sub-zero-sup-minus in Eq. KST(23), solution for pseudo-spherical beam source */
    yb_0u,  /* y-sub-zero-sup-plus  "                                                         */
    yb_1d,  /* y-sub-one-sup-minus  "                                                         */
    yb_1u,  /* y-sub-one-sup-plus   "                                                         */
    yp_0d,  /* y-sub-zero-sup-minus in Eq. KST(23), solution for thermal source               */
    yp_0u,  /* y-sub-zero-sup-plus  "                                                         */
    yp_1d,  /* y-sub-one-sup-minus  "                                                         */
    yp_1u,  /* y-sub-one-sup-plus   "                                                         */
    zb_a,   /* Alfa coefficient in Eq. KST(22) for pseudo-spherical beam source               */
    zp_a;   /* Alfa coefficient in Eq. KST(22) for thermal source                             */
} twostr_xyz;

/*
 * Array shift macros
 * Using unit-offset shift macros to match Fortran version
 *
 * NOTE: ARRAY(iq,jq) is defined locally instead of here, because its size is different
 *       in different subroutines.
 */
#define A(i,j)           a[i-1+(j-1)*lda]
#define AA(j,k)          aa[j-1+(k-1)*ia]
#define ABD(i,j)         abd[i-1+(j-1)*lda]
#define ALBMED(iu)       out->albmed[iu-1]
#define AMB(iq,jq)       ab[iq-1+(jq-1)*(ds->nstr/2)].zero
#define APB(iq,jq)       ab[iq-1+(jq-1)*(ds->nstr/2)].one

#define B(iq)            b[iq-1]
#define BDR(iq,jq)       bdr[iq-1+(jq)*(ds->nstr/2)]
#define BEM(iq)          bem[iq-1]

#define CBAND(irow,ncol) cband[irow-1+(ncol-1)*(9*(ds->nstr/2)-2)]
#define CC(iq,jq)        cc[iq-1+(jq-1)*ds->nstr]
#define CH(lc)           ch[lc-1]
#define CHTAU(ls)        chtau[ls]
#define CMU(iq)          cmu[iq-1]
#define CWT(iq)          cwt[iq-1]

#define DFDT(lu)         out->rad[lu-1].dfdt
#define DIAG(i)          diag[i-1].on
#define DTAUC(lc)        ds->dtauc[lc-1]
#define DTAU_C(lc)       dtau_c[lc-1]
#define DTAUCPR(lc)      dtaucpr[lc-1]

#define EMU(iu)          emu[iu-1]
#define EVAL(j)          eval[j-1]
#define EVEC(j,k)        evec[j-1+(k-1)*ievec]
#define EVECC(iq,jq)     evecc[iq-1+(jq-1)*ds->nstr]
#define EXPBEA(lc)       expbea[lc]

#define FLDIR(lu)        fl[lu-1].zero
#define FLDN(lu)         fl[lu-1].one
#define FLUP(lu)         out->rad[lu-1].flup
#define FLYR(lc)         flyr[lc-1]

#define GC(iq,jq,lc)     gc[iq-1+(jq-1+(lc-1)*ds->nstr)*ds->nstr]
#define GENSRC(maz,lc,iq)  ds->gensrc[iq-1+(lc-1+maz*ds->nlyr)*ds->nstr]
#define GENSRCU(maz,lc,iu) ds->gensrcu[iu-1+(lc-1+maz*ds->nlyr)*ds->numu]
#define GG(lc)           gg[lc-1]
#define GGPRIM(lc)       ggprim[lc-1]
#define GL(k,lc)         gl[k+(lc-1)*(ds->nstr+1)]
#define GMU(k)           gmu[k-1]
#define GU(iu,iq,lc)     gu[iu-1+(iq-1+(lc-1)*ds->nstr)*ds->numu]
#define GWT(k)           gwt[k-1]

#define IERROR(i)        ierror[i-1]
#define IPVT(k)          ipvt[k-1]

#define KK(iq,lc)        kk[iq-1+(lc-1)*ds->nstr]

#define LAYRU(lu)        layru[lu-1]
#define LL(iq,lc)        ll[iq-1+(lc-1)*ds->nstr]

#define MU(i)            mu[i-1]

#define OMEGA(lyr)       omega[lyr-1]
#define OPRIM(lc)        oprim[lc-1]

#define PKAG(lc)         pkag[lc]
#define PKAGC(lc)        pkagc[lc-1]
#define PHASA(lc)        phasa[lc-1]
#define PHASE(lc)        phase[lc-1]
#define PHASM(lc)        phasm[lc-1]
#define PHAST(lc)        phast[lc-1]
#define PHI(j)           ds->phi[j-1]
#define PHIRAD(jp)       phirad[jp-1]
#define PMOM(k,lc)       ds->pmom[k+(lc-1)*(ds->nmom_nstr+1)]
#define PRNTU0(i)        prntu0[i-1]
#define PSI0(iq)         psi[iq-1].zero
#define PSI1(iq)         psi[iq-1].one

#define RFLDIR(lu)       out->rad[lu-1].rfldir
#define RFLDN(lu)        out->rad[lu-1].rfldn
#define RMU(iu,iq)       rmu[iu-1+(iq)*ds->numu]
#define RR(lc)           rr[lc-1]

#define SSALB(lc)        ds->ssalb[lc-1]
#define SUBD(i)          diag[i-1].sub
#define SUPERD(i)        diag[i-1].super
#define SX(i)            sx[i-1]
#define SY(i)            sy[i-1]

#define TAU(lc)          tau[lc]
#define TAUC(lc)         tauc[lc]
#define TAUCPR(lc)       taucpr[lc]
#define TEMPER(lc)       ds->temper[lc]
#define TRNMED(iu)       out->trnmed[iu-1]

#define U0C(iq,lu)       u0c[iq-1+(lu-1)*ds->nstr]
#define U0U(iu,lu)       out->u0u[iu-1+(lu-1)*ds->numu]
#define UAVG(lu)         out->rad[lu-1].uavg
#define UAVGDN(lu)       out->rad[lu-1].uavgdn
#define UAVGUP(lu)       out->rad[lu-1].uavgup
#define UAVGSO(lu)       out->rad[lu-1].uavgso
#define UMU(iu)          ds->umu[iu-1]
#define UTAU(lu)         ds->utau[lu-1]
#define UTAUPR(lu)       utaupr[lu-1]
#define UUM(iu,lu)       uum[iu-1+(lu-1)*ds->numu]
#define UU(iu,lu,j)      out->uu[iu-1+(lu-1+(j-1)*ds->ntau)*ds->numu]
#define OUT_UUM(iu,lu,j) out->uum[iu-1+(lu-1+(j)*ds->ntau)*ds->numu] /* No -i behind j as mazim starts at 0, aky */

#define WK(iq)           wk[iq-1]

#define XBA(lc)          xba[lc]
#define XB0(iq,lc)       xb[iq-1+(lc-1)*ds->nstr].zero
#define XB1(iq,lc)       xb[iq-1+(lc-1)*ds->nstr].one
#define XB_0D(lc)        ts[lc-1].xb_0d
#define XB_0U(lc)        ts[lc-1].xb_0u
#define XB_1D(lc)        ts[lc-1].xb_1d
#define XB_1U(lc)        ts[lc-1].xb_1u
#define XP_0(lc)         ts[lc-1].xp_0
#define XP_1(lc)         ts[lc-1].xp_1
#define XR0(lc)          xr[lc-1].zero
#define XR1(lc)          xr[lc-1].one

#define YB_0D(lc)        ts[lc-1].yb_0d
#define YB_0U(lc)        ts[lc-1].yb_0u
#define YB_1D(lc)        ts[lc-1].yb_1d
#define YB_1U(lc)        ts[lc-1].yb_1u
#define YLM(l,i)         ylm[l+(i-1)*(maxmu+1)]
#define YLM0(iq)         ylm0[iq]
#define YLMC(l,iq)       ylmc[l+(iq-1)*(ds->nstr+1)]
#define YLMU(l,iu)       ylmu[l+(iu-1)*(ds->nstr+1)]
#define YP_0D(lc)        ts[lc-1].yp_0d
#define YP_0U(lc)        ts[lc-1].yp_0u
#define YP_1D(lc)        ts[lc-1].yp_1d
#define YP_1U(lc)        ts[lc-1].yp_1u

#define Z(j)             z[j-1]
#define Z0(iu)           zee[iu-1].zero
#define Z1(iq)           zee[iq-1].one
#define Z0U(iu,lc)       zu[iu-1+(lc-1)*ds->numu].zero
#define Z1U(iu,lc)       zu[iu-1+(lc-1)*ds->numu].one
#define ZB0U(iu,lc)      zbu[iu-1+(lc-1)*ds->numu].zero
#define ZB1U(iu,lc)      zbu[iu-1+(lc-1)*ds->numu].one
#define ZBAU(iu,lc)      zbu[iu-1+(lc-1)*ds->numu].alpha
#define ZB_A(lc)         ts[lc-1].zb_a
#define ZBEAM(iu,lc)     zbeam[iu-1+(lc-1)*ds->numu]
#define ZBEAMA(lc)       zbeama[lc-1]
#define ZBEAM0(iq,lc)    zbeamsp[iq-1+(lc-1)*ds->nstr].zero
#define ZBEAM1(iq,lc)    zbeamsp[iq-1+(lc-1)*ds->nstr].one
#define ZBS0(iq)         zbs[iq-1].zero
#define ZBS1(iq)         zbs[iq-1].one
#define ZD(j)            zd[j]
#define ZJ(j)            zj[j-1]
#define ZJG(j)           zjg[j-1]
#define ZJU(j)           zju[j-1]
#define ZGU(iu,lc)       zgu[iu-1+(lc-1)*ds->numu]
#define ZP_A(lc)         ts[lc-1].zp_a
#define ZPLK0(iq,lc)     plk[iq-1+(lc-1)*ds->nstr].zero
#define ZPLK1(iq,lc)     plk[iq-1+(lc-1)*ds->nstr].one
#define ZZ(iq,lc)        zz[iq-1+(lc-1)*ds->nstr]
#define ZZG(iq,lc)       zzg[iq-1+(lc-1)*ds->nstr]

/* BDE stuff */
#define MUP(it)          mu_phase[it-1]
#define PHASR(lc)        phasr[lc-1]
#define PHAS2(it,lc)     phas2[it-1+(lc-1)*nphase]
#define DSPHASE(it,lc)   ds->phase[it-1+(lc-1)*ds->nphase]
#define F_PHAS2_ABS(it)  f_phas2_abs[it-1]
#define MU_EQ(i,lu)      mu_eq[i-1+(lu-1)*nf]
#define NEG_PHAS(i,lu)   neg_phas[i-1+(lu-1)*nf]
#define NORM_PHAS(lu)    norm_phas[lu-1]

/* setout.f, inter.f stuff */
#define SDTAUC(i)        sdtauc[i-1]
#define SUTAU(i)         sutau[i-1]
#define ZOUT(i)          zout[i-1]
#define TAUINT(i)        tauint[i-1]
#define XARR(i)          xarr[i-1]
#define YARR(i)          yarr[i-1]

/*
 * Logical
 */
#define TRUE  1
#define FALSE 0

#define FIRST_IPHAS          1
#define ISOTROPIC            1
#define RAYLEIGH             2
#define HENYEY_GREENSTEIN    3
#define HAZE_GARCIA_SIEWERT  4
#define CLOUD_GARCIA_SIEWERT 5
#define LAST_IPHAS           5

#define GENERAL_BC 0
#define SPECIAL_BC 1

#define TOP_ILLUM 1
#define BOT_ILLUM 2

#define DS_WARNING 0
#define DS_ERROR   1

#define VERBOSE 0
#define QUIET   1

#define BRDF_NONE   0
#define BRDF_RPV    1  /* don't change these numbers as they are */
#define BRDF_CAM    2  /* used by Fortran code which of course   */
#define BRDF_AMB    3  /* has no access to this header file      */
#define BRDF_HAPKE  4

/*defined for new option names brdf_cam for cox_and_munk_sal,pcl,u10,uphi*/
#define BRDF_CAM_NN   4
#define BRDF_CAM_SAL  0 
#define BRDF_CAM_PCL  1
#define BRDF_CAM_U10  2
#define BRDF_CAM_UPHI 3

/*
 * NMUG : Number of angle cosine quadrature points on (-1,1) for integrating bidirectional reflectivity
 *        to get directional emissivity (it is necessary to use a quadrature set distinct from the 
 *        computational angles, because the computational angles may not be dense enough---ds->nstr
 *        may be too small---to give an accurate approximation for the integration).
 */
#define NMUG 50

/*
 * Mathematical
 */
#if !defined(M_E)
#  define M_E         2.7182818284590452354
#  define M_LOG2E     1.4426950408889634074
#  define M_LOG10E    0.43429448190325182765
#  define M_LN2       0.69314718055994530942
#  define M_LN10      2.30258509299404568402
#  define M_PI        3.14159265358979323846
#  define M_PI_2      1.57079632679489661923
#  define M_PI_4      0.78539816339744830962
#  define M_1_PI      0.31830988618379067154
#  define M_2_PI      0.63661977236758134308
#  define M_2_SQRTPI  1.12837916709551257390
#  define M_SQRT2     1.41421356237309504880
#  define M_SQRT1_2   0.70710678118654752440
#endif

#define DEG (M_PI/180.)

#define SQR(x) ({ \
          const double _x = (double)(x); \
          _x*_x; })

#define MIN(x,y) ({ \
         const double _x = (double)(x); \
         const double _y = (double)(y); \
         _x < _y ? _x : _y; })

#define MAX(x,y) ({ \
         const double _x = (double)(x); \
         const double _y = (double)(y); \
         _x > _y ? _x : _y; })

#define LIMIT_RANGE(min,x,max) ({ \
         const double _min = (double)(min); \
         const double _x   = (double)(x);   \
         const double _max = (double)(max); \
         _x < _min ? _min : ( _x > _max ? _max : _x ); })

#define IMIN(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i < _j ? _i : _j; })

#define IMAX(i,j) ({ \
         const int _i = (int)(i); \
         const int _j = (int)(j); \
         _i > _j ? _i : _j; })

#define F77_SIGN(a,b) ((b) >= 0. ? fabs(a) : -fabs(a))

/*---------------------*
 * Function prototypes *
 *---------------------*/

void c_disort(disort_state  *ds,
              disort_output *out);

double c_bidir_reflectivity ( double       wvnmlo,
			      double       wvnmhi,
			      double       mu,
			      double       mup,
			      double       dphi,
			      int          brdf_type,
			      disort_brdf *brdf,
			      int          callnum );

double c_bidir_reflectivity_hapke ( double wvnmlo,
				    double wvnmhi,
				    double mu,
				    double mup,
				    double dphi );

double c_bidir_reflectivity_rpv ( rpv_brdf_spec *brdf,
				  double         mu1,
				  double         mu2,
				  double         phi,
				  double         badmu );

double c_dref(double       wvnmlo,
              double       wvnmhi,
              double       mu,
	      int          brdf_type,
	      disort_brdf *brdf,
	      int          callnum );

void c_getmom(int    iphas,
             double  gg,
             int     nmom,
             double *pmom);

void c_asymmetric_matrix(double *aa,
                         double *evec,
                         double *eval,
                         int     m,
                         int     ia,
                         int     ievec,
                         int    *ier,
                         double *wk);

void c_intensity_components(disort_state *ds,
                            double       *gc,
                            double       *kk,
                            int          *layru,
                            double       *ll,
                            int           lyrcut,
                            int           mazim,
                            int           ncut,
                            int           nn,
                            double       *taucpr,
                            double       *utaupr,
                            double       *zz,
                            disort_pair  *plk,
                            double       *uum);

void c_fluxes(disort_state  *ds,
              disort_output *out,
              double        *ch,
              double        *cmu,
              double        *cwt,
              double        *gc,
              double        *kk,
              int           *layru,
              double        *ll,
              int            lyrcut,
              int            ncut,
              int            nn,
              int            prntu0,
              double        *taucpr,
              double        *utaupr,
              disort_pair   *xr,
	      disort_pair   *zbeamsp,
              double        *zbeama,
              double        *zz,
              double        *zzg,
              disort_pair   *plk,
              disort_pair   *fl,
              double        *u0c);

void c_intensity_correction(disort_state  *ds,
                            disort_output *out,
                            double         dither,
                            double        *flyr,
                            int           *layru,
                            int            lyrcut,
                            int            ncut,
                            double        *oprim,
                            double        *phasa,
                            double        *phast,
                            double        *phasm,
                            double        *phirad,
                            double        *tauc,
                            double        *taucpr,
                            double        *utaupr);

void c_new_intensity_correction(disort_state  *ds,
				disort_output *out,
				double         dither,
				double        *flyr,
				int           *layru,
				int            lyrcut,
				int            ncut,
				double        *oprim,
				double        *phasa,
				double        *phast,
				double        *phasm,
				double        *phirad,
				double        *tauc,
				double        *taucpr,
				double        *utaupr);

void prep_double_scat_integr (int           nphase,
			      int           ntau,
			      int           nf,
			      double       *mu_phase,
			      double       *phas2,
			      double       *mu_eq,
			      int          *neg_phas,
			      double       *norm_phas);

double c_secondary_scat(disort_state *ds,
                        int           iu,
                        int           lu,
                        double        ctheta,
                        double       *flyr,
                        int           layru,
                        double       *tauc);

double c_new_secondary_scat(disort_state *ds,
			    int           iu,
			    int           lu,
			    int           it,
			    double        ctheta,
			    double       *flyr,
			    int           layru,
			    double       *tauc,
			    int           nf,
			    double       *phas2,
			    double       *mu_eq,
			    int          *neg_phas,
			    double        norm_phas);

double calc_phase_squared (int           nphase,
			   int           lu,
			   double        ctheta,
			   int           nf,
			   double       *mu_phase,
			   double       *phas2,
			   double       *mu_eq,
			   int          *neg_phas,
			   double        norm_phas);

void c_disort_set(disort_state *ds,
                  double       *ch,
                  double       *chtau,
                  double       *cmu,
                  double       *cwt,
                  int           deltam,
                  double       *dtaucpr,
                  double       *expbea,
                  double       *flyr,
                  double       *gl,
                  int          *layru,
                  int          *lyrcut,
                  int          *ncut,
                  int          *nn,
                  int          *corint,
                  double       *oprim,
                  double       *tauc,
                  double       *taucpr,
                  double       *utaupr);

void c_set_matrix(disort_state *ds,
                  double       *bdr,
                  double       *cband,
                  double       *cmu,
                  double       *cwt,
                  double        delm0,
                  double       *dtaucpr,
                  double       *gc,
                  double       *kk,
                  int           lyrcut,
                  int          *ncol,
                  int           ncut,
                  double       *taucpr,
                  double       *wk);

double c_single_scat(double   dither,
                     int      layru,
                     int      nlyr,
                     double  *phase,
                     double  *omega,
                     double  *tau,
                     double   umu,
                     double   umu0,
                     double   utau,
                     double   fbeam);

void c_solve_eigen(disort_state *ds,
                   int           lc,
                   disort_pair  *ab,
                   double       *array,
                   double       *cmu,
                   double       *cwt,
                   double       *gl,
                   int           mazim,
                   int           nn,
                   double       *ylmc,
                   double       *cc,
                   double       *evecc,
                   double       *eval,
                   double       *kk,
                   double       *gc,
                   double       *wk);

void c_solve0(disort_state *ds,
              double       *b,
              double       *bdr,
              double       *bem,
              double        bplanck,
              double       *cband,
              double       *cmu,
              double       *cwt,
              double       *expbea,
              int          *ipvt,
              double       *ll,
              int           lyrcut,
              int           mazim,
              int           ncol,
              int           ncut,
              int           nn,
              double        tplanck,
              double       *taucpr,
              double       *z,
	      disort_pair  *zbeamsp,
	      double       *zbeama,
              double       *zz,
              double       *zzg,
              disort_pair  *plk);

void c_surface_bidir(disort_state *ds,
                     double        delm0,
                     double       *cmu,
                     int           mazim,
                     int           nn,
                     double       *bdr,
                     double       *emu,
                     double       *bem,
                     double       *rmu,
		     int           callnum);

void c_interp_eigenvec(disort_state *ds,
                       int           lc,
                       double       *cwt,
                       double       *evecc,
                       double       *gl,
                       double       *gu,
                       int           mazim,
                       int           nn,
                       double       *wk,
                       double       *ylmc,
                       double       *ylmu);

void c_interp_source(disort_state   *ds,
                     int             lc,
                     double         *cwt,
                     double          delm0,
                     double         *gl,
                     int             mazim,
                     double         *oprim,
                     double         *ylm0,
                     double         *ylmc,
                     double         *ylmu,
                     disort_pair    *psi,
                     disort_pair    *xr,
                     disort_pair    *zee,
                     double         *zj,
                     double         *zjg,
                     double         *zbeam,
		     disort_triplet *zbu,
		     disort_pair    *zbs,
		     double          zbsa,
		     double         *zgu,
                     disort_pair    *zu);

void c_set_coefficients_beam_source(disort_state *ds,
				    double       *ch,
				    double       *chtau,
				    double       *cmu, 
				    double        delm0,
				    double        fbeam,
				    double       *gl,
				    int           lc,
				    int           mazim,
				    int           nstr,
				    double       *taucpr,
				    double       *xba,
				    disort_pair  *xb,
				    double       *ylm0,
				    double       *ylmc,
				    double       *zj);

void c_interp_coefficients_beam_source(disort_state   *ds,
				       double         *chtau,
				       double          delm0,
				       double          fbeam,
				       double         *gl,
				       int             lc,
				       int             mazim,
				       int             nstr,
				       int             numu,
				       double         *taucpr,
				       disort_triplet *zbu,
				       double         *xba,
				       double         *zj,
				       double         *ylm0,
				       double         *ylmu);

void c_upbeam(disort_state *ds,
              int           lc,
              double       *array,
              double       *cc,
              double       *cmu,
              double        delm0,
              double       *gl,
              int          *ipvt,
              int           mazim,
              int           nn,
              double       *wk,
              double       *ylm0,
              double       *ylmc,
              double       *zj,
              double       *zz);

void c_upbeam_pseudo_spherical(disort_state *ds,
			       int           lc,
			       double       *array, 
			       double       *cc,
			       double       *cmu, 
			       int          *ipvt, 
			       int           nn,
			       double       *wk,
			       disort_pair  *xb,
			       double       *xba, 
			       disort_pair  *zbs,
			       double       *zbsa,
			       disort_pair  *zbeamsp,
			       double       *zbeama);

void c_upbeam_general_source(disort_state *ds,
			     int           lc,
 			     int           maz,
			     double       *array,
			     double       *cc,
			     int          *ipvt,
			     int           nn,
			     double       *wk,
			     double       *zjg,
			     double       *zzg);

void c_upisot(disort_state *ds,
              int           lc,
              double       *array,
              double       *cc,
              double       *cmu,
              int          *ipvt,
              int           nn,
              double       *oprim,
              double       *wk,
              disort_pair  *xr,
              disort_pair  *zee,
              disort_pair  *plk);

void c_user_intensities(disort_state   *ds,
                        double          bplanck,
                        double         *cmu,
                        double         *cwt,
                        double          delm0,
                        double         *dtaucpr,
                        double         *emu,
                        double         *expbea,
                        double         *gc,
                        double         *gu,
                        double         *kk,
                        int            *layru,
                        double         *ll,
                        int             lyrcut,
                        int             mazim,
                        int             ncut,
                        int             nn,
                        double         *rmu,
                        double         *taucpr,
                        double          tplanck,
                        double         *utaupr,
                        double         *wk,
			disort_triplet *zbu,
                        double         *zbeam,
			disort_pair    *zbeamsp,
                        double         *zbeama,
                        double         *zgu,
                        disort_pair    *zu,
                        double         *zz,
                        double         *zzg,
                        disort_pair    *plk,
                        double         *uum);

double c_xi_func(double umu1,
                 double umu2,
                 double tau);

void c_check_inputs(disort_state *ds,
                    int           scat_yes,
                    int           deltam,
                    int           corint,
                    double       *tauc,
		    int           callnum);

void c_legendre_poly(int     nmu,
                     int     m,
                     int     maxmu,
                     int     twonm1,
                     double *mu,
                     double *ylm);

void c_print_avg_intensities(disort_state *ds,
			     disort_output *out);

void c_print_inputs(disort_state *ds,
                    double       *dtaucpr,
                    int           scat_yes,
                    int           deltam,
                    int           corint,
                    double       *flyr,
                    int           lyrcut,
                    double       *oprim,
                    double       *tauc,
                    double       *taucpr);

void c_print_intensities(disort_state  *ds,
                         disort_output *out);

void c_gaussian_quadrature(int    m,
                           double *gmu,
                           double *gwt);

double c_ratio(double a,
               double b);

int c_fcmp(double x1,
           double x2);

void c_self_test(int            compare,
                 int           *prntu0,
                 disort_state  *ds,
                 disort_output *out);

void c_albtrans(disort_state  *ds,
                disort_output *out,
                disort_pair   *ab,
                double        *array,
                double        *b,
                double        *bdr,
                double        *cband,
                double        *cc,
                double        *cmu,
                double        *cwt,
                double        *dtaucpr,
                double        *eval,
                double        *evecc,
                double        *gl,
                double        *gc,
                double        *gu,
                int           *ipvt,
                double        *kk,
                double        *ll,
                int            nn,
                double        *taucpr,
                double        *ylmc,
                double        *ylmu,
                double        *z,
                double        *wk);

void c_albtrans_intensity(disort_state *ds,
			  disort_output *out,
                          double       *gu,
                          double       *kk,
                          double       *ll,
                          int           nn,
                          double       *taucpr,
                          double       *wk);

void c_print_albtrans(disort_state  *ds,
                      disort_output *out);

void c_solve1(disort_state *ds,
              double       *cband,
              int           ihom,
              int          *ipvt,
              int           ncol,
              int           ncut,
              int           nn,
              double       *b,
              double       *ll);

void c_albtrans_spherical(disort_state *ds,
                          double       *cmu,
                          double       *cwt,
                          double       *gc,
                          double       *kk,
                          double       *ll,
                          int           nn,
                          double       *taucpr,
                          double       *sflup,
                          double       *sfldn);

void c_errmsg(char *messag,
              int   type);

int c_write_bad_var(int   quiet,
                    char *varnam);

int c_write_too_small_dim(int   quiet,
                          char *dimnam,
                          int   minval);

void c_sgbco(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             double *rcond,
             double *z);

void c_sgbfa(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             int    *info);

void c_sgbsl(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             double *b,
             int     job);

void c_sgeco(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             double *rcond,
             double *z);

void c_sgefa(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             int    *info);

void c_sgesl(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             double *b,
             int     job);

double c_sasum(int     n,
               double *sx);

void c_saxpy(int     n,
             double  sa,
             double *sx,
             double *sy);

double c_sdot(int     n,
              double *sx,
              double *sy);

void c_sscal(int    n,
            double  sa,
            double *sx);

int c_isamax(int     n,
             double *sx);

void c_twostr(disort_state  *ds,
              disort_output *out,
              int            deltam,
              double        *gg,
              int           *ierror,
              double         radius);

double c_chapman(int     lc,
                 double  taup,
                 double *tauc,
                 int     nlyr,
                 double *zd,
                 double *dtau_c,
                 double  zenang,
                 double  r);

double c_chapman_simpler(int     lc,
                 double  taup,
                 int     nlyr,
                 double *zd,
                 double *dtau_c,
                 double  zenang,
                 double  r);

void c_twostr_check_inputs(disort_state *ds,
                           double       *gg,
                           int          *ierror,
                           double       *tauc);

void c_twostr_fluxes(disort_state  *ds,
                   twostr_xyz      *ts,
                   double          *ch,
                   double           cmu,
                   double          *kk,
                   int             *layru,
                   double          *ll,
                   int              lyrcut,
                   int              ncut,
                   double          *oprim,
                   double          *rr,
                   double          *taucpr,
                   double          *utaupr,
                   disort_output   *out,
                   double          *u0c,
                   disort_pair     *fl);

void c_twostr_solns(disort_state *ds,
                    double       *ch,
                    double       *chtau,
                    double        cmu,
                    int           ncut,
                    double       *oprim,
                    double       *pkag,
                    double       *pkagc,
                    double       *taucpr,
                    double       *ggprim,
                    double       *kk,
                    double       *rr,
                    twostr_xyz   *ts);

void c_twostr_print_inputs(disort_state *ds,
                           int           deltam,
                           double       *flyr,
                           double       *gg,
                           int           lyrcut,
                           double       *oprim,
                           double       *tauc,
                           double       *taucpr);

void c_twostr_set(disort_state *ds,
                  double       *bplanck,
                  double       *ch,
                  double       *chtau,
                  double       *cmu,
                  int           deltam,
                  double       *dtaucpr,
                  double       *expbea,
                  double       *flyr,
                  double       *gg,
                  double       *ggprim,
                  int          *layru,
                  int          *lyrcut,
                  int          *ncut,
                  int          *nn,
                  double       *oprim,
                  double       *pkag,
                  double       *pkagc,
                  double        radius,
                  double       *tauc,
                  double       *taucpr,
                  double       *tplanck,
                  double       *utaupr);

void c_twostr_solve_bc(disort_state *ds,
                       twostr_xyz   *ts,
                       double        bplanck,
                       double       *cband,
                       double        cmu,
                       double       *expbea,
                       int           lyrcut,
                       int           nn,
                       int           ncut,
                       double        tplanck,
                       double       *taucpr,
                       double       *kk,
                       double       *rr,
                       int          *ipvt,
                       double       *b,
                       double       *ll,
                       twostr_diag  *diag);

double c_planck_func1(double wnumlo,
                      double wnumhi,
                      double t);

double c_planck_func2(double wnumlo,
                      double wnumhi,
                      double t);

void c_disort_state_alloc(disort_state *ds);

void c_disort_state_free(disort_state *ds);

void c_disort_out_alloc(disort_state  *ds,
                        disort_output *out);

void c_disort_out_free(disort_state  *ds,
                       disort_output *out);

void c_twostr_state_alloc(disort_state *ds);

void c_twostr_state_free(disort_state *ds);

void c_twostr_out_alloc(disort_state  *ds,
                        disort_output *out);

void c_twostr_out_free(disort_state  *ds,
                       disort_output *out);

double *c_dbl_vector(int  nl, 
                     int  nh,
		     char *name);
int *c_int_vector(int  nl, 
		  int  nh,
		  char *name);

void print_test(disort_state  *ds_calc,
                disort_output *calc,
                disort_state  *ds_good,
                disort_output *good);

void c_free_dbl_vector(double *m, 
                       int     nl, 
                       int     nh);

int c_setout( float *sdtauc,
	      int    nlyr,
	      int    ntau,
	      float *sutau,
	      float *z,
	      float *zout );

double c_inter( int     npoints,
		int     itype,
		double  arg,
		float  *xarr,
		double *yarr,
		double *hh );

int c_gaussian_quadrature_test(int nstr, float *sza, double umu0);

void disort_test01(void);
void disort_test02(void);
void disort_test03(void);
void disort_test04(void);
void disort_test05(void);
void disort_test06(void);
void disort_test07(void);
void disort_test08(void);
void disort_test09(void);
void disort_test10(void);
void disort_test11(void);
void disort_test12(void);
void disort_test13(void);
void disort_test14(void);

/* * * * * * * * * * * * * * * * * * * * * * * end of cdisort.h * * * * * * * * * * * * * * * * * * */

#endif
