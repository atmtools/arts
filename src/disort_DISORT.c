/* DISORT.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__101 = 101;
static integer c__100 = 100;
static integer c__200 = 200;
static integer c__48 = 48;
static integer c__3 = 3;
static integer c__1000 = 1000;
static integer c__49 = 49;
static integer c__2304 = 2304;
static integer c__4900 = 4900;
static integer c__2352 = 2352;
static integer c__4800 = 4800;
static integer c_b51 = 230400;
static integer c__24 = 24;
static integer c__576 = 576;
static integer c__214 = 214;
static integer c__9600 = 9600;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__0 = 0;
static integer c__10 = 10;

/* ~~~~~~~~~~~~ */
/* VERSION 1.2 */
/* ~~~~~~~~~~~~ */
/* Subroutine */ int disort_(integer *nlyr, doublereal *dtauc, doublereal *
	ssalb, doublereal *pmom, doublereal *temper, doublereal *wvnmlo, 
	doublereal *wvnmhi, logical *usrtau, integer *ntau, doublereal *utau, 
	integer *nstr, logical *usrang, integer *numu, doublereal *umu, 
	integer *nphi, doublereal *phi, integer *ibcnd, doublereal *fbeam, 
	doublereal *umu0, doublereal *phi0, doublereal *fisot, logical *
	lamber, doublereal *albedo, doublereal *hl, doublereal *btemp, 
	doublereal *ttemp, doublereal *temis, logical *deltam, logical *plank,
	 logical *onlyfl, doublereal *accur, logical *prnt, char *header, 
	integer *maxcly, integer *maxulv, integer *maxumu, integer *maxcmu, 
	integer *maxphi, doublereal *rfldir, doublereal *rfldn, doublereal *
	flup, doublereal *dfdt, doublereal *uavg, doublereal *uu, doublereal *
	u0u, doublereal *albmed, doublereal *trnmed, ftnlen header_len)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    address a__1[2];
    integer dtauc_dim1, flup_dim1, phi_dim1, pmom_dim1, pmom_dim2, 
	    pmom_offset, rfldir_dim1, rfldn_dim1, ssalb_dim1, temper_dim1, 
	    u0u_dim1, u0u_dim2, u0u_offset, umu_dim1, utau_dim1, uu_dim1, 
	    uu_dim2, uu_dim3, uu_offset, i__1, i__2, i__3, i__4, i__5, i__6, 
	    i__7, i__8[2], i__9, i__10, i__11, i__12;
    doublereal d__1, d__2, d__3, d__4, d__5;
    char ch__1[136];

    /* Builtin functions */
    double asin(doublereal);
    integer s_rnge(char *, integer, char *, integer);
    double sqrt(doublereal);
    integer i_len(char *, ftnlen), s_wsfe(cilist *);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer do_fio(integer *, char *, ftnlen), e_wsfe();
    double cos(doublereal);

    /* Local variables */
    static doublereal pkag[101], fldn[200], eval[24];
    static integer ncol;
    static doublereal tauc[101];
    static integer ncos;
    static doublereal ylmc[2352]	/* was [49][48] */, hlpr[49];
    static integer ncut;
    static doublereal flyr[100];
    static integer ipvt[4800];
    static doublereal ylmu[2352]	/* was [49][48] */, delm0, zplk0[4800]
	    	/* was [48][100] */, zplk1[4800]	/* was [48][100] */, 
	    b[4800], cband[1027200]	/* was [214][4800] */;
    static integer j, l;
    static doublereal evecc[2304]	/* was [48][48] */, z__[4800], evald[
	    24], zbeam[4800]	/* was [48][100] */, fldir[200];
    static integer mazim;
    static doublereal array[2304]	/* was [48][48] */;
    extern doublereal ratio_(doublereal *, doublereal *);
    static integer kconv;
    static doublereal azerr, oprim[100];
    static integer layru[200];
    static doublereal z0[48], z1[48];
    extern doublereal r1mach_(integer *);
    static doublereal cc[2304]	/* was [48][48] */;
    extern /* Subroutine */ int solve0_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     logical *, doublereal *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal gc[230400]	/* was [48][48][100] */;
    static integer lc;
    static doublereal gl[4900]	/* was [49][100] */, pi;
    static integer iq;
    static doublereal kk[4800]	/* was [48][100] */;
    static integer nn;
    static doublereal gu[230400]	/* was [48][48][100] */;
    static integer iu;
    static doublereal ll[4800]	/* was [48][100] */, eveccd[576]	/* 
	    was [24][24] */;
    static integer lu, ns;
    static doublereal expbea[101], wk[48], bplank, phirad[3], zj[48], angcos;
    extern /* Subroutine */ int chekin_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , integer *, doublereal *, integer *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , logical *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), upbeam_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal dither, dtaucp[100];
    static logical compar;
    extern /* Subroutine */ int albtrn_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal zz[4800]	/* was [48][100] */, cosphi;
    extern doublereal plkavg_(doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int soleig_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal tplank, taucpr[101];
    extern /* Subroutine */ int cmpint_(doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    logical *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), pravin_(
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *);
    static doublereal azterm;
    extern /* Subroutine */ int lepoly_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *), fluxes_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, logical *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), setdis_(doublereal *, 
	    doublereal *, logical *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, logical *, integer *, logical *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, logical *, integer *, logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     logical *);
    static doublereal u0c[9600]	/* was [48][200] */;
    extern /* Subroutine */ int surfac_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, logical *, 
	    doublereal *, logical *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), prtinp_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, logical *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, logical *
	    , logical *, logical *, doublereal *, doublereal *, logical *, 
	    doublereal *, doublereal *, doublereal *, integer *, logical *), 
	    slftst_(doublereal *, doublereal *, doublereal *, logical *, 
	    doublereal *, doublereal *, doublereal *, integer *, logical *, 
	    integer *, logical *, integer *, integer *, integer *, integer *, 
	    logical *, doublereal *, doublereal *, doublereal *, logical *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, logical *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal utaupr[200];
    extern /* Subroutine */ int prtint_(doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *);
    static logical lyrcut;
    extern /* Subroutine */ int setmtx_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *), terpev_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), terpso_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, logical *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), upisot_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), usrint_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     integer *, doublereal *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    logical *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal xr0[100], xr1[100];
    extern /* Subroutine */ int zeroal_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , integer *, doublereal *, integer *, doublereal *), zeroit_(
	    doublereal *, integer *);
    static doublereal z0u[4800]	/* was [48][100] */, z1u[4800]	/* was [48][
	    100] */, aad[576]	/* was [24][24] */, amb[576]	/* was [24][
	    24] */, apb[576]	/* was [24][24] */, bem[24], bdr[600]	/* 
	    was [24][25] */, cmu[48], dum;
    static integer lev;
    static doublereal rpd;
    static integer naz;
    static doublereal sgn, emu[48], psi[48], wkd[48], cwt[48], rmu[1200]	
	    /* was [48][25] */, uum[9600]	/* was [48][200] */, sqt[1000]
	    , ylm0[49];

    /* Fortran I/O blocks */
    static cilist io___9 = { 0, 6, 0, "(//,1X,100('*'),/,A,/,1X,100('*'))", 0 
	    };


/* ******************************************************************* */
/*       Plane-parallel discrete ordinates radiative transfer program */
/*             ( see DISORT.doc for complete documentation ) */
/* ******************************************************************* */

/* +------------------------------------------------------------------+ */
/*  Calling Tree (omitting calls to ERRMSG): */
/*  (routines in parentheses are not in this file) */

/*  DISORT-+-(R1MACH) */
/*         +-SLFTST-+-(TSTBAD) */
/*         +-ZEROIT */
/*         +-CHEKIN-+-(WRTBAD) */
/*         |        +-(WRTDIM) */
/*         |        +-DREF */
/*         +-ZEROAL */
/*         +-SETDIS-+-QGAUSN (1)-+-(D1MACH) */
/*         +-PRTINP */
/*         +-ALBTRN-+-LEPOLY (2) */
/*         |        +-ZEROIT */
/*         |        +-SOLEIG (3)-+-ASYMTX-+-(D1MACH) */
/*         |        +-TERPEV */
/*         |        +-SETMTX (4)--ZEROIT */
/*         |        +-(SGBCO) */
/*         |        +-SOLVE1-+-ZEROIT */
/*         |        |        +-(SGBSL) */
/*         |        +-ALTRIN */
/*         |        +-SPALTR */
/*         |        +-PRALTR */
/*         +-PLKAVG-+-(R1MACH) */
/*         +-LEPOLY see 2 */
/*         +-SURFAC-+-QGAUSN see 1 */
/*         |        +-LEPOLY see 2 */
/*         |        +-ZEROIT */
/*         +-SOLEIG see 3 */
/*         +-UPBEAM-+-(SGECO) */
/*         |        +-(SGESL) */
/*         +-UPISOT-+-(SGECO) */
/*         |        +-(SGESL) */
/*         +-TERPEV */
/*         +-TERPSO */
/*         +-SETMTX see 4 */
/*         +-SOLVE0-+-ZEROIT */
/*         |        +-(SGBCO) */
/*         |        +-(SGBSL) */
/*         +-FLUXES--ZEROIT */
/*         +-USRINT */
/*         +-CMPINT */
/*         +-PRAVIN */
/*         +-RATIO--(R1MACH) */
/*         +-PRTINT */

/* *** Intrinsic Functions used in DISORT package which take */
/*     non-negligible amount of time: */

/*    EXP :  Called by- ALBTRN, ALTRIN, CMPINT, FLUXES, SETDIS, */
/*                      SETMTX, SPALTR, USRINT, PLKAVG */

/*    SQRT : Called by- ASYMTX, SOLEIG */

/* +-------------------------------------------------------------------+ */

/*  Index conventions (for all DO-loops and all variable descriptions): */

/*     IU     :  for user polar angles */

/*  IQ,JQ,KQ  :  for computational polar angles ('quadrature angles') */

/*   IQ/2     :  for half the computational polar angles (just the ones */
/*               in either 0-90 degrees, or 90-180 degrees) */

/*     J      :  for user azimuthal angles */

/*     K,L    :  for Legendre expansion coefficients or, alternatively, */
/*               subscripts of associated Legendre polynomials */

/*     LU     :  for user levels */

/*     LC     :  for computational layers (each having a different */
/*               single-scatter albedo and/or phase function) */

/*    LEV     :  for computational levels */

/*    MAZIM   :  for azimuthal components in Fourier cosine expansion */
/*               of intensity and phase function */

/* +------------------------------------------------------------------+ */

/*               I N T E R N A L    V A R I A B L E S */

/*   AMB(IQ/2,IQ/2)    First matrix factor in reduced eigenvalue problem */
/*                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG) */

/*   APB(IQ/2,IQ/2)    Second matrix factor in reduced eigenvalue problem */
/*                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG) */

/*   ARRAY(IQ,IQ)      Scratch matrix for SOLEIG, UPBEAM and UPISOT */
/*                     (see each subroutine for definition) */

/*   B()               Right-hand side vector of Eq. SC(5) going into */
/*                     SOLVE0,1;  returns as solution vector */
/*                     vector  L, the constants of integration */

/*   BDR(IQ/2,0:IQ/2)  Bottom-boundary bidirectional reflectivity for a */
/*                     given azimuthal component.  First index always */
/*                     refers to a computational angle.  Second index: */
/*                     if zero, refers to incident beam angle UMU0; */
/*                     if non-zero, refers to a computational angle. */

/*   BEM(IQ/2)         Bottom-boundary directional emissivity at compu- */
/*                     tational angles. */

/*   BPLANK            Intensity emitted from bottom boundary */

/*   CBAND()           Matrix of left-hand side of the linear system */
/*                     Eq. SC(5), scaled by Eq. SC(12);  in banded */
/*                     form required by LINPACK solution routines */

/*   CC(IQ,IQ)         C-sub-IJ in Eq. SS(5) */

/*   CMU(IQ)           Computational polar angles (Gaussian) */

/*   CWT(IQ)           Quadrature weights corresponding to CMU */

/*   DELM0             Kronecker delta, delta-sub-M0, where M = MAZIM */
/*                     is the number of the Fourier component in the */
/*                     azimuth cosine expansion */

/*   DITHER            Small quantity subtracted from single-scattering */
/*                     albedos of unity, in order to avoid using special */
/*                     case formulas;  prevents an eigenvalue of exactly */
/*                     zero from occurring, which would cause an */
/*                     immediate overflow */

/*   DTAUCP(LC)        Computational-layer optical depths (delta-M-scaled */
/*                     if DELTAM = TRUE, otherwise equal to DTAUC) */

/*   EMU(IU)           Bottom-boundary directional emissivity at user */
/*                     angles. */

/*   EVAL(IQ)          Temporary storage for eigenvalues of Eq. SS(12) */

/*   EVECC(IQ,IQ)      Complete eigenvectors of SS(7) on return from */
/*                     SOLEIG; stored permanently in  GC */

/*   EXPBEA(LC)        Transmission of direct beam in delta-M optical */
/*                     depth coordinates */

/*   FLYR(LC)          Truncated fraction in delta-M method */

/*   GL(K,LC)          Phase function Legendre polynomial expansion */
/*                     coefficients, calculated from PMOM by */
/*                     including single-scattering albedo, factor */
/*                     2K+1, and (if DELTAM=TRUE) the delta-M */
/*                     scaling */

/*   GC(IQ,IQ,LC)      Eigenvectors at polar quadrature angles, */
/*                     g  in Eq. SC(1) */

/*   GU(IU,IQ,LC)      Eigenvectors interpolated to user polar angles */
/*                     ( g  in Eqs. SC(3) and S1(8-9), i.e. */
/*                       G without the L factor ) */

/*   HLPR()            Legendre coefficients of bottom bidirectional */
/*                     reflectivity (after inclusion of 2K+1 factor) */

/*   IPVT(LC*IQ)       Integer vector of pivot indices for LINPACK */
/*                     routines */

/*   KK(IQ,LC)         Eigenvalues of coeff. matrix in Eq. SS(7) */

/*   KCONV             Counter in azimuth convergence test */

/*   LAYRU(LU)         Computational layer in which user output level */
/*                     UTAU(LU) is located */

/*   LL(IQ,LC)         Constants of integration L in Eq. SC(1), */
/*                     obtained by solving scaled version of Eq. SC(5) */

/*   LYRCUT            TRUE, radiation is assumed zero below layer */
/*                     NCUT because of almost complete absorption */

/*   NAZ               Number of azimuthal components considered */

/*   NCUT              Computational layer number in which absorption */
/*                     optical depth first exceeds ABSCUT */

/*   OPRIM(LC)         Single scattering albedo after delta-M scaling */

/*   PASS1             TRUE on first entry, FALSE thereafter */

/*   PKAG(0:LC)        Integrated Planck function for internal emission */

/*   PSI(IQ)           Sum just after square bracket in  Eq. SD(9) */

/*   RMU(IU,0:IQ)      Bottom-boundary bidirectional reflectivity for a */
/*                     given azimuthal component.  First index always */
/*                     refers to a user angle.  Second index: */
/*                     if zero, refers to incident beam angle UMU0; */
/*                     if non-zero, refers to a computational angle. */

/*   SQT(k)            Square root of k (used only in LEPOLY for */
/*                     computing associated Legendre polynomials) */

/*   TAUC(0:LC)        Cumulative optical depth (un-delta-M-scaled) */

/*   TAUCPR(0:LC)      Cumulative optical depth (delta-M-scaled if */
/*                     DELTAM = TRUE, otherwise equal to TAUC) */

/*   TPLANK            Intensity emitted from top boundary */

/*   UUM(IU,LU)        Expansion coefficients when the intensity */
/*                     (u-super-M) is expanded in Fourier cosine series */
/*                     in azimuth angle */

/*   U0C(IQ,LU)        Azimuthally-averaged intensity */

/*   UTAUPR(LU)        Optical depths of user output levels in delta-M */
/*                     coordinates;  equal to  UTAU(LU) if no delta-M */

/*   WK()              scratch array */

/*   XR0(LC)           X-sub-zero in expansion of thermal source func- */
/*                     tion preceding Eq. SS(14) (has no mu-dependence) */

/*   XR1(LC)           X-sub-one in expansion of thermal source func- */
/*                     tion;  see  Eqs. SS(14-16) */

/*   YLM0(L)           Normalized associated Legendre polynomial */
/*                     of subscript L at the beam angle (not saved */
/*                     as function of superscipt M) */

/*   YLMC(L,IQ)        Normalized associated Legendre polynomial */
/*                     of subscript L at the computational angles */
/*                     (not saved as function of superscipt M) */

/*   YLMU(L,IU)        Normalized associated Legendre polynomial */
/*                     of subscript L at the user angles */
/*                     (not saved as function of superscipt M) */

/*   Z()               scratch array used in SOLVE0, ALBTRN to solve */
/*                     a linear system for the constants of integration */

/*   Z0(IQ)            Solution vectors Z-sub-zero of Eq. SS(16) */

/*   Z0U(IU,LC)        Z-sub-zero in Eq. SS(16) interpolated to user */
/*                     angles from an equation derived from SS(16) */

/*   Z1(IQ)            Solution vectors Z-sub-one  of Eq. SS(16) */

/*   Z1U(IU,LC)        Z-sub-one in Eq. SS(16) interpolated to user */
/*                     angles from an equation derived from SS(16) */

/*   ZBEAM(IU,LC)      Particular solution for beam source */

/*   ZJ(IQ)            Right-hand side vector  X-sub-zero in */
/*                     Eq. SS(19), also the solution vector */
/*                     Z-sub-zero after solving that system */

/*   ZZ(IQ,LC)         Permanent storage for the beam source vectors ZJ */

/*   ZPLK0(IQ,LC)      Permanent storage for the thermal source */
/*                     vectors  Z0  obtained by solving  Eq. SS(16) */

/*   ZPLK1(IQ,LC)      Permanent storage for the thermal source */
/*                     vectors  Z1  obtained by solving  Eq. SS(16) */

/* +-------------------------------------------------------------------+ */

/*  LOCAL SYMBOLIC DIMENSIONS (have big effect on storage requirements): */

/*       MXCLY  = Max no. of computational layers */
/*       MXULV  = Max no. of output levels */
/*       MXCMU  = Max no. of computation polar angles */
/*       MXUMU  = Max no. of output polar angles */
/*       MXPHI  = Max no. of output azimuthal angles */
/*       MXSQT  = Max no. of square roots of integers (for LEPOLY) */
/* +-------------------------------------------------------------------+ */
/*     .. Parameters .. */
/*      PARAMETER ( MXCLY = 6, MXULV = 5, MXCMU = 48, MXUMU = 10, */
/*    &          MXPHI = 3, MI = MXCMU / 2, MI9M2 = 9*MI - 2, */
/*    &          NNLYRI = MXCMU*MXCLY, MXSQT = 1000 ) */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    temper_dim1 = *maxcly - 0 + 1;
    ssalb_dim1 = *maxcly;
    dtauc_dim1 = *maxcly;
    flup_dim1 = *maxulv;
    rfldn_dim1 = *maxulv;
    rfldir_dim1 = *maxulv;
    utau_dim1 = *maxulv;
    u0u_dim1 = *maxumu;
    u0u_dim2 = *maxulv;
    u0u_offset = 1 + u0u_dim1 * 1;
    umu_dim1 = *maxumu;
    pmom_dim1 = *maxcmu - 0 + 1;
    pmom_dim2 = *maxcly;
    pmom_offset = 0 + pmom_dim1 * 1;
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_dim3 = *maxphi;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2 * 1);
    phi_dim1 = *maxphi;

    /* Function Body */
    if (pass1) {
	pi = asin(1.) * 2.;
	dither = r1mach_(&c__4) * 10.;
/*                            ** Must dither more on Cray (14-digit prec) */
	if (dither < 1e-10) {
	    dither *= 10.;
	}
	rpd = pi / 180.;
	for (ns = 1; ns <= 1000; ++ns) {
	    sqt[(i__1 = ns - 1) < 1000 && 0 <= i__1 ? i__1 : s_rnge("sqt", 
		    i__1, "disort_", (ftnlen)391)] = sqrt((doublereal) ns);
/* L5: */
	}
/*                            ** Set input values for self-test. */
/*                            ** Be sure SLFTST sets all print flags off. */
	compar = FALSE_;
	slftst_(accur, albedo, btemp, deltam, &dtauc[(i__1 = 0) < 1 * 
		dtauc_dim1 ? i__1 : s_rnge("dtauc", i__1, "disort_", (ftnlen)
		397)], fbeam, fisot, ibcnd, lamber, nlyr, plank, nphi, numu, 
		nstr, ntau, onlyfl, &phi[(i__2 = 0) < 1 * phi_dim1 ? i__2 : 
		s_rnge("phi", i__2, "disort_", (ftnlen)397)], phi0, &pmom[(
		i__3 = pmom_dim1 - pmom_offset) < 1 * pmom_dim1 * pmom_dim2 &&
		 0 <= i__3 ? i__3 : s_rnge("pmom", i__3, "disort_", (ftnlen)
		397)], prnt, &ssalb[(i__4 = 0) < 1 * ssalb_dim1 ? i__4 : 
		s_rnge("ssalb", i__4, "disort_", (ftnlen)397)], temis, &
		temper[(i__5 = 0) < 1 * temper_dim1 ? i__5 : s_rnge("temper", 
		i__5, "disort_", (ftnlen)397)], ttemp, &umu[(i__6 = 0) < 1 * 
		umu_dim1 ? i__6 : s_rnge("umu", i__6, "disort_", (ftnlen)397)]
		, usrang, usrtau, &utau[(i__7 = 0) < 1 * utau_dim1 ? i__7 : 
		s_rnge("utau", i__7, "disort_", (ftnlen)397)], umu0, wvnmhi, 
		wvnmlo, &compar, &dum, &dum, &dum, &dum);
    }
L10:
    if (! pass1 && i_len(header, (ftnlen)127) != 0) {
	s_wsfe(&io___9);
/* Writing concatenation */
	i__8[0] = 9, a__1[0] = " DISORT: ";
	i__8[1] = 127, a__1[1] = header;
	s_cat(ch__1, a__1, i__8, &c__2, (ftnlen)136);
	do_fio(&c__1, ch__1, (ftnlen)136);
	e_wsfe();
    }
/*                                  ** Calculate cumulative optical depth */
/*                                  ** and dither single-scatter albedo */
/*                                  ** to improve numerical behavior of */
/*                                  ** eigenvalue/vector computation */
    zeroit_(tauc, &c__101);
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	if (ssalb[(i__2 = lc - 1) < 1 * ssalb_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("ssalb", i__2, "disort_", (ftnlen)420)] == 1.) {
	    ssalb[(i__3 = lc - 1) < 1 * ssalb_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("ssalb", i__3, "disort_", (ftnlen)420)] = 1. - 
		    dither;
	}
	tauc[(i__2 = lc) < 101 && 0 <= i__2 ? i__2 : s_rnge("tauc", i__2, 
		"disort_", (ftnlen)421)] = tauc[(i__3 = lc - 1) < 101 && 0 <= 
		i__3 ? i__3 : s_rnge("tauc", i__3, "disort_", (ftnlen)421)] + 
		dtauc[(i__4 = lc - 1) < 1 * dtauc_dim1 && 0 <= i__4 ? i__4 : 
		s_rnge("dtauc", i__4, "disort_", (ftnlen)421)];
/* L20: */
    }
/*                                ** Check input dimensions and variables */
    chekin_(nlyr, dtauc, ssalb, pmom, temper, wvnmlo, wvnmhi, usrtau, ntau, 
	    utau, nstr, usrang, numu, umu, nphi, phi, ibcnd, fbeam, umu0, 
	    phi0, fisot, lamber, albedo, hl, btemp, ttemp, temis, plank, 
	    onlyfl, accur, tauc, maxcly, maxulv, maxumu, maxcmu, maxphi, &
	    c__100, &c__200, &c__48, &c__48, &c__3, &c__1000);
/*                                 ** Zero internal and output arrays */
    i__1 = *maxumu * *maxulv;
    i__2 = *maxumu * *maxulv * *maxphi;
    zeroal_(&c__100, &expbea[1], flyr, oprim, &taucpr[1], xr0, xr1, &c__48, 
	    cmu, cwt, psi, wk, z0, z1, zj, &c__49, hlpr, ylm0, &c__2304, 
	    array, cc, evecc, &c__4900, gl, &c__2352, ylmc, &c__2352, ylmu, &
	    c__4800, kk, ll, zz, zplk0, zplk1, &c_b51, gc, &c__200, layru, 
	    utaupr, &c_b51, gu, &c__4800, z0u, z1u, zbeam, &c__24, eval, &
	    c__576, amb, apb, &c__4800, ipvt, z__, maxulv, rfldir, rfldn, 
	    flup, uavg, dfdt, maxumu, albmed, trnmed, &i__1, u0u, &i__2, uu);
/*                                 ** Perform various setup operations */
    setdis_(cmu, cwt, deltam, dtauc, dtaucp, expbea, fbeam, flyr, gl, hl, 
	    hlpr, ibcnd, lamber, layru, &lyrcut, maxumu, maxcmu, &c__48, &
	    ncut, nlyr, ntau, &nn, nstr, plank, numu, onlyfl, oprim, pmom, 
	    ssalb, tauc, taucpr, utau, utaupr, umu, umu0, usrtau, usrang);
/*                                 ** Print input information */
    if (prnt[0]) {
	prtinp_(nlyr, dtauc, dtaucp, ssalb, pmom, temper, wvnmlo, wvnmhi, 
		ntau, utau, nstr, numu, umu, nphi, phi, ibcnd, fbeam, umu0, 
		phi0, fisot, lamber, albedo, hl, btemp, ttemp, temis, deltam, 
		plank, onlyfl, accur, flyr, &lyrcut, oprim, tauc, taucpr, 
		maxcmu, &prnt[6]);
    }
/*                              ** Handle special case for getting albedo */
/*                              ** and transmissivity of medium for many */
/*                              ** beam angles at once */
    if (*ibcnd == 1) {
	albtrn_(albedo, amb, apb, array, b, bdr, cband, cc, cmu, cwt, dtaucp, 
		eval, evecc, gl, gc, gu, ipvt, kk, ll, nlyr, &nn, nstr, numu, 
		prnt, taucpr, umu, u0u, wk, ylmc, ylmu, z__, aad, evald, 
		eveccd, wkd, &c__24, &c__214, maxulv, maxumu, &c__48, &c__48, 
		&c__4800, sqt, albmed, trnmed);
	return 0;
    }
/*                                   ** Calculate Planck functions */
    if (! (*plank)) {
	bplank = 0.;
	tplank = 0.;
	zeroit_(pkag, &c__101);
    } else {
	tplank = *temis * plkavg_(wvnmlo, wvnmhi, ttemp);
	bplank = plkavg_(wvnmlo, wvnmhi, btemp);
	i__1 = *nlyr;
	for (lev = 0; lev <= i__1; ++lev) {
	    pkag[(i__2 = lev) < 101 && 0 <= i__2 ? i__2 : s_rnge("pkag", i__2,
		     "disort_", (ftnlen)499)] = plkavg_(wvnmlo, wvnmhi, &
		    temper[(i__3 = lev) < 1 * temper_dim1 && 0 <= i__3 ? i__3 
		    : s_rnge("temper", i__3, "disort_", (ftnlen)499)]);
/* L30: */
	}
    }
/* ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  ======= */
/*           (EQ STWJ 5) */
    kconv = 0;
    naz = *nstr - 1;
/*                                    ** Azimuth-independent case */
    if (*fbeam == 0. || (d__1 = 1. - *umu0, abs(d__1)) < 1e-5 || *onlyfl || *
	    numu == 1 && (d__2 = 1. - umu[(i__1 = 0) < 1 * umu_dim1 ? i__1 : 
	    s_rnge("umu", i__1, "disort_", (ftnlen)512)], abs(d__2)) < 1e-5 ||
	     *numu == 1 && (d__3 = umu[(i__2 = 0) < 1 * umu_dim1 ? i__2 : 
	    s_rnge("umu", i__2, "disort_", (ftnlen)512)] + 1., abs(d__3)) < 
	    1e-5 || *numu == 2 && (d__4 = umu[(i__3 = 0) < 1 * umu_dim1 ? 
	    i__3 : s_rnge("umu", i__3, "disort_", (ftnlen)512)] + 1., abs(
	    d__4)) < 1e-5 && (d__5 = 1. - umu[(i__4 = 1) < 1 * umu_dim1 ? 
	    i__4 : s_rnge("umu", i__4, "disort_", (ftnlen)512)], abs(d__5)) < 
	    1e-5) {
	naz = 0;
    }
    i__1 = naz;
    for (mazim = 0; mazim <= i__1; ++mazim) {
	if (mazim == 0) {
	    delm0 = 1.;
	}
	if (mazim > 0) {
	    delm0 = 0.;
	}
/*                             ** Get normalized associated Legendre */
/*                             ** polynomials for */
/*                             ** (a) incident beam angle cosine */
/*                             ** (b) computational and user polar angle */
/*                             **     cosines */
	if (*fbeam > 0.) {
	    ncos = 1;
	    angcos = -(*umu0);
	    i__2 = *nstr - 1;
	    lepoly_(&ncos, &mazim, &c__48, &i__2, &angcos, sqt, ylm0);
	}
	if (! (*onlyfl) && *usrang) {
	    i__2 = *nstr - 1;
	    lepoly_(numu, &mazim, &c__48, &i__2, umu, sqt, ylmu);
	}
	i__2 = *nstr - 1;
	lepoly_(&nn, &mazim, &c__48, &i__2, cmu, sqt, ylmc);
/*                       ** Get normalized associated Legendre polys. */
/*                       ** with negative arguments from those with */
/*                       ** positive arguments; Dave/Armstrong Eq. (15) */
	sgn = -1.;
	i__2 = *nstr - 1;
	for (l = mazim; l <= i__2; ++l) {
	    sgn = -sgn;
	    i__3 = *nstr;
	    for (iq = nn + 1; iq <= i__3; ++iq) {
		ylmc[(i__4 = l + iq * 49 - 49) < 2352 && 0 <= i__4 ? i__4 : 
			s_rnge("ylmc", i__4, "disort_", (ftnlen)555)] = sgn * 
			ylmc[(i__5 = l + (iq - nn) * 49 - 49) < 2352 && 0 <= 
			i__5 ? i__5 : s_rnge("ylmc", i__5, "disort_", (ftnlen)
			555)];
/* L40: */
	    }
/* L50: */
	}
/*                                 ** Specify users bottom reflectivity */
/*                                 ** and emissivity properties */
	if (! lyrcut) {
	    surfac_(albedo, &delm0, fbeam, hlpr, lamber, &c__24, &mazim, &
		    c__48, &c__48, &nn, numu, nstr, onlyfl, umu, usrang, ylm0,
		     ylmc, ylmu, bdr, emu, bem, rmu, sqt);
	}
/* ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  ============= */
	i__2 = ncut;
	for (lc = 1; lc <= i__2; ++lc) {
/*                        ** Solve eigenfunction problem in Eq. STWJ(8B); */
/*                        ** return eigenvalues and eigenvectors */
	    soleig_(amb, apb, array, cmu, cwt, &gl[(i__3 = lc * 49 - 49) < 
		    4900 && 0 <= i__3 ? i__3 : s_rnge("gl", i__3, "disort_", (
		    ftnlen)575)], &c__24, &mazim, &c__48, &nn, nstr, ylmc, cc,
		     evecc, eval, &kk[(i__4 = lc * 48 - 48) < 4800 && 0 <= 
		    i__4 ? i__4 : s_rnge("kk", i__4, "disort_", (ftnlen)575)],
		     &gc[(i__5 = (lc * 48 + 1) * 48 - 2352) < 230400 && 0 <= 
		    i__5 ? i__5 : s_rnge("gc", i__5, "disort_", (ftnlen)575)],
		     aad, eveccd, evald, wkd);
/*                                  ** Calculate particular solutions of */
/*                                  ** Eq.SS(18) for incident beam source */
	    if (*fbeam > 0.) {
		upbeam_(array, cc, cmu, &delm0, fbeam, &gl[(i__3 = lc * 49 - 
			49) < 4900 && 0 <= i__3 ? i__3 : s_rnge("gl", i__3, 
			"disort_", (ftnlen)582)], ipvt, &mazim, &c__48, &nn, 
			nstr, &pi, umu0, wk, ylm0, ylmc, zj, &zz[(i__4 = lc * 
			48 - 48) < 4800 && 0 <= i__4 ? i__4 : s_rnge("zz", 
			i__4, "disort_", (ftnlen)582)]);
	    }
/*                              ** Calculate particular solutions of */
/*                              ** Eq. SS(15) for thermal emission source */
	    if (*plank && mazim == 0) {
		xr1[(i__3 = lc - 1) < 100 && 0 <= i__3 ? i__3 : s_rnge("xr1", 
			i__3, "disort_", (ftnlen)592)] = 0.;
		if (dtaucp[(i__3 = lc - 1) < 100 && 0 <= i__3 ? i__3 : s_rnge(
			"dtaucp", i__3, "disort_", (ftnlen)594)] > 0.) {
		    xr1[(i__4 = lc - 1) < 100 && 0 <= i__4 ? i__4 : s_rnge(
			    "xr1", i__4, "disort_", (ftnlen)594)] = (pkag[(
			    i__5 = lc) < 101 && 0 <= i__5 ? i__5 : s_rnge(
			    "pkag", i__5, "disort_", (ftnlen)594)] - pkag[(
			    i__6 = lc - 1) < 101 && 0 <= i__6 ? i__6 : s_rnge(
			    "pkag", i__6, "disort_", (ftnlen)594)]) / dtaucp[(
			    i__7 = lc - 1) < 100 && 0 <= i__7 ? i__7 : s_rnge(
			    "dtaucp", i__7, "disort_", (ftnlen)594)];
		}
		xr0[(i__3 = lc - 1) < 100 && 0 <= i__3 ? i__3 : s_rnge("xr0", 
			i__3, "disort_", (ftnlen)597)] = pkag[(i__4 = lc - 1) 
			< 101 && 0 <= i__4 ? i__4 : s_rnge("pkag", i__4, 
			"disort_", (ftnlen)597)] - xr1[(i__5 = lc - 1) < 100 
			&& 0 <= i__5 ? i__5 : s_rnge("xr1", i__5, "disort_", (
			ftnlen)597)] * taucpr[(i__6 = lc - 1) < 101 && 0 <= 
			i__6 ? i__6 : s_rnge("taucpr", i__6, "disort_", (
			ftnlen)597)];
		upisot_(array, cc, cmu, ipvt, &c__48, &nn, nstr, &oprim[(i__3 
			= lc - 1) < 100 && 0 <= i__3 ? i__3 : s_rnge("oprim", 
			i__3, "disort_", (ftnlen)599)], wk, &xr0[(i__4 = lc - 
			1) < 100 && 0 <= i__4 ? i__4 : s_rnge("xr0", i__4, 
			"disort_", (ftnlen)599)], &xr1[(i__5 = lc - 1) < 100 
			&& 0 <= i__5 ? i__5 : s_rnge("xr1", i__5, "disort_", (
			ftnlen)599)], z0, z1, &zplk0[(i__6 = lc * 48 - 48) < 
			4800 && 0 <= i__6 ? i__6 : s_rnge("zplk0", i__6, 
			"disort_", (ftnlen)599)], &zplk1[(i__7 = lc * 48 - 48)
			 < 4800 && 0 <= i__7 ? i__7 : s_rnge("zplk1", i__7, 
			"disort_", (ftnlen)599)]);
	    }
	    if (! (*onlyfl) && *usrang) {
/*                                            ** Interpolate eigenvectors */
/*                                            ** to user angles */
		terpev_(cwt, evecc, &gl[(i__3 = lc * 49 - 49) < 4900 && 0 <= 
			i__3 ? i__3 : s_rnge("gl", i__3, "disort_", (ftnlen)
			610)], &gu[(i__4 = (lc * 48 + 1) * 48 - 2352) < 
			230400 && 0 <= i__4 ? i__4 : s_rnge("gu", i__4, "dis\
ort_", (ftnlen)610)], &mazim, &c__48, &c__48, &nn, nstr, numu, wk, ylmc, ylmu)
			;
/*                                            ** Interpolate source terms */
/*                                            ** to user angles */
		terpso_(cwt, &delm0, fbeam, &gl[(i__3 = lc * 49 - 49) < 4900 
			&& 0 <= i__3 ? i__3 : s_rnge("gl", i__3, "disort_", (
			ftnlen)616)], &mazim, &c__48, plank, numu, nstr, &
			oprim[(i__4 = lc - 1) < 100 && 0 <= i__4 ? i__4 : 
			s_rnge("oprim", i__4, "disort_", (ftnlen)616)], &pi, 
			ylm0, ylmc, ylmu, psi, &xr0[(i__5 = lc - 1) < 100 && 
			0 <= i__5 ? i__5 : s_rnge("xr0", i__5, "disort_", (
			ftnlen)616)], &xr1[(i__6 = lc - 1) < 100 && 0 <= i__6 
			? i__6 : s_rnge("xr1", i__6, "disort_", (ftnlen)616)],
			 z0, zj, &zbeam[(i__7 = lc * 48 - 48) < 4800 && 0 <= 
			i__7 ? i__7 : s_rnge("zbeam", i__7, "disort_", (
			ftnlen)616)], &z0u[(i__9 = lc * 48 - 48) < 4800 && 0 
			<= i__9 ? i__9 : s_rnge("z0u", i__9, "disort_", (
			ftnlen)616)], &z1u[(i__10 = lc * 48 - 48) < 4800 && 0 
			<= i__10 ? i__10 : s_rnge("z1u", i__10, "disort_", (
			ftnlen)616)]);
	    }
/* L60: */
	}
/* ===================  END LOOP ON COMPUTATIONAL LAYERS  =============== */
/*                      ** Set coefficient matrix of equations combining */
/*                      ** boundary and layer interface conditions */
	setmtx_(bdr, cband, cmu, cwt, &delm0, dtaucp, gc, kk, lamber, &lyrcut,
		 &c__24, &c__214, &c__48, &ncol, &ncut, &c__4800, &nn, nstr, 
		taucpr, wk);
/*                      ** Solve for constants of integration in homo- */
/*                      ** geneous solution (general boundary conditions) */
	solve0_(b, bdr, bem, &bplank, cband, cmu, cwt, expbea, fbeam, fisot, 
		ipvt, lamber, ll, &lyrcut, &mazim, &c__24, &c__214, &c__48, &
		ncol, &ncut, &nn, nstr, &c__4800, &pi, &tplank, taucpr, umu0, 
		z__, zz, zplk0, zplk1);
/*                                  ** Compute upward and downward fluxes */
	if (mazim == 0) {
	    fluxes_(cmu, cwt, fbeam, gc, kk, layru, ll, &lyrcut, maxulv, &
		    c__48, &c__200, &ncut, &nn, nstr, ntau, &pi, prnt, ssalb, 
		    taucpr, umu0, utau, utaupr, xr0, xr1, zz, zplk0, zplk1, 
		    dfdt, flup, fldn, fldir, rfldir, rfldn, uavg, u0c);
	}
	if (*onlyfl) {
	    if (*maxumu >= *nstr) {
/*                                     ** Save azimuthal-avg intensities */
/*                                     ** at quadrature angles */
		i__2 = *ntau;
		for (lu = 1; lu <= i__2; ++lu) {
		    i__3 = *nstr;
		    for (iq = 1; iq <= i__3; ++iq) {
			u0u[(i__4 = iq + lu * u0u_dim1 - u0u_offset) < 1 * 
				u0u_dim1 * u0u_dim2 && 0 <= i__4 ? i__4 : 
				s_rnge("u0u", i__4, "disort_", (ftnlen)660)] =
				 u0c[(i__5 = iq + lu * 48 - 49) < 9600 && 0 <=
				 i__5 ? i__5 : s_rnge("u0c", i__5, "disort_", 
				(ftnlen)660)];
/* L70: */
		    }
/* L80: */
		}
	    }
	    goto L170;
	}
	zeroit_(uum, &c__9600);
	if (*usrang) {
/*                                     ** Compute azimuthal intensity */
/*                                     ** components at user angles */
	    usrint_(&bplank, cmu, cwt, &delm0, dtaucp, emu, expbea, fbeam, 
		    fisot, gc, gu, kk, lamber, layru, ll, &lyrcut, &mazim, &
		    c__48, &c__200, &c__48, &ncut, nlyr, &nn, nstr, plank, 
		    numu, ntau, &pi, rmu, taucpr, &tplank, umu, umu0, utaupr, 
		    wk, zbeam, z0u, z1u, zz, zplk0, zplk1, uum);
	} else {
/*                                     ** Compute azimuthal intensity */
/*                                     ** components at quadrature angles */
	    cmpint_(fbeam, gc, kk, layru, ll, &lyrcut, &mazim, &c__48, &
		    c__200, &c__48, &ncut, &nn, nstr, plank, ntau, taucpr, 
		    umu0, utaupr, zz, zplk0, zplk1, uum);
	}
	if (mazim == 0) {
/*                               ** Save azimuthally averaged intensities */
	    i__2 = *ntau;
	    for (lu = 1; lu <= i__2; ++lu) {
		i__3 = *numu;
		for (iu = 1; iu <= i__3; ++iu) {
		    u0u[(i__4 = iu + lu * u0u_dim1 - u0u_offset) < 1 * 
			    u0u_dim1 * u0u_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			    "u0u", i__4, "disort_", (ftnlen)701)] = uum[(i__5 
			    = iu + lu * 48 - 49) < 9600 && 0 <= i__5 ? i__5 : 
			    s_rnge("uum", i__5, "disort_", (ftnlen)701)];
		    i__4 = *nphi;
		    for (j = 1; j <= i__4; ++j) {
			uu[(i__5 = iu + (lu + j * uu_dim2) * uu_dim1 - 
				uu_offset) < 1 * uu_dim1 * uu_dim2 * uu_dim3 
				&& 0 <= i__5 ? i__5 : s_rnge("uu", i__5, 
				"disort_", (ftnlen)704)] = uum[(i__6 = iu + 
				lu * 48 - 49) < 9600 && 0 <= i__6 ? i__6 : 
				s_rnge("uum", i__6, "disort_", (ftnlen)704)];
/* L90: */
		    }
/* L100: */
		}
/* L110: */
	    }
/*                              ** Print azimuthally averaged intensities */
/*                              ** at user angles */
	    if (prnt[3]) {
		pravin_(umu, numu, maxumu, utau, ntau, u0u);
	    }
	    if (naz > 0) {
		zeroit_(phirad, &c__3);
		i__2 = *nphi;
		for (j = 1; j <= i__2; ++j) {
		    phirad[(i__3 = j - 1) < 3 && 0 <= i__3 ? i__3 : s_rnge(
			    "phirad", i__3, "disort_", (ftnlen)719)] = rpd * (
			    phi[(i__4 = j - 1) < 1 * phi_dim1 && 0 <= i__4 ? 
			    i__4 : s_rnge("phi", i__4, "disort_", (ftnlen)719)
			    ] - *phi0);
/* L120: */
		}
	    }
	} else {
/*                                ** Increment intensity by current */
/*                                ** azimuthal component (Fourier */
/*                                ** cosine series);  Eq SD(2) */
	    azerr = 0.;
	    i__2 = *nphi;
	    for (j = 1; j <= i__2; ++j) {
		cosphi = cos(mazim * phirad[(i__3 = j - 1) < 3 && 0 <= i__3 ? 
			i__3 : s_rnge("phirad", i__3, "disort_", (ftnlen)733)]
			);
		i__3 = *ntau;
		for (lu = 1; lu <= i__3; ++lu) {
		    i__4 = *numu;
		    for (iu = 1; iu <= i__4; ++iu) {
			azterm = uum[(i__5 = iu + lu * 48 - 49) < 9600 && 0 <=
				 i__5 ? i__5 : s_rnge("uum", i__5, "disort_", 
				(ftnlen)738)] * cosphi;
			uu[(i__5 = iu + (lu + j * uu_dim2) * uu_dim1 - 
				uu_offset) < 1 * uu_dim1 * uu_dim2 * uu_dim3 
				&& 0 <= i__5 ? i__5 : s_rnge("uu", i__5, 
				"disort_", (ftnlen)739)] = uu[(i__6 = iu + (
				lu + j * uu_dim2) * uu_dim1 - uu_offset) < 1 *
				 uu_dim1 * uu_dim2 * uu_dim3 && 0 <= i__6 ? 
				i__6 : s_rnge("uu", i__6, "disort_", (ftnlen)
				739)] + azterm;
/* Computing MAX */
			d__4 = abs(azterm);
			d__5 = (d__1 = uu[(i__5 = iu + (lu + j * uu_dim2) * 
				uu_dim1 - uu_offset) < 1 * uu_dim1 * uu_dim2 *
				 uu_dim3 && 0 <= i__5 ? i__5 : s_rnge("uu", 
				i__5, "disort_", (ftnlen)740)], abs(d__1));
			d__2 = azerr, d__3 = ratio_(&d__4, &d__5);
			azerr = max(d__2,d__3);
/* L130: */
		    }
/* L140: */
		}
/* L150: */
	    }
	    if (azerr <= *accur) {
		++kconv;
	    }
	    if (kconv >= 2) {
		goto L170;
	    }
	}
/* L160: */
    }
/* ===================  END LOOP ON AZIMUTHAL COMPONENTS  =============== */
/*                                          ** Print intensities */
L170:
    if (prnt[4] && ! (*onlyfl)) {
	prtint_(uu, utau, ntau, umu, numu, phi, nphi, maxulv, maxumu);
    }
    if (pass1) {
/*                                    ** Compare test case results with */
/*                                    ** correct answers and abort if bad */
	compar = TRUE_;
	slftst_(accur, albedo, btemp, deltam, &dtauc[(i__1 = 0) < 1 * 
		dtauc_dim1 ? i__1 : s_rnge("dtauc", i__1, "disort_", (ftnlen)
		770)], fbeam, fisot, ibcnd, lamber, nlyr, plank, nphi, numu, 
		nstr, ntau, onlyfl, &phi[(i__2 = 0) < 1 * phi_dim1 ? i__2 : 
		s_rnge("phi", i__2, "disort_", (ftnlen)770)], phi0, &pmom[(
		i__3 = pmom_dim1 - pmom_offset) < 1 * pmom_dim1 * pmom_dim2 &&
		 0 <= i__3 ? i__3 : s_rnge("pmom", i__3, "disort_", (ftnlen)
		770)], prnt, &ssalb[(i__4 = 0) < 1 * ssalb_dim1 ? i__4 : 
		s_rnge("ssalb", i__4, "disort_", (ftnlen)770)], temis, &
		temper[(i__5 = 0) < 1 * temper_dim1 ? i__5 : s_rnge("temper", 
		i__5, "disort_", (ftnlen)770)], ttemp, &umu[(i__6 = 0) < 1 * 
		umu_dim1 ? i__6 : s_rnge("umu", i__6, "disort_", (ftnlen)770)]
		, usrang, usrtau, &utau[(i__7 = 0) < 1 * utau_dim1 ? i__7 : 
		s_rnge("utau", i__7, "disort_", (ftnlen)770)], umu0, wvnmhi, 
		wvnmlo, &compar, &flup[(i__9 = 0) < 1 * flup_dim1 ? i__9 : 
		s_rnge("flup", i__9, "disort_", (ftnlen)770)], &rfldir[(i__10 
		= 0) < 1 * rfldir_dim1 ? i__10 : s_rnge("rfldir", i__10, 
		"disort_", (ftnlen)770)], &rfldn[(i__11 = 0) < 1 * rfldn_dim1 
		? i__11 : s_rnge("rfldn", i__11, "disort_", (ftnlen)770)], &
		uu[(i__12 = (uu_dim2 + 1) * uu_dim1 + 1 - uu_offset) < 1 * 
		uu_dim1 * uu_dim2 * uu_dim3 && 0 <= i__12 ? i__12 : s_rnge(
		"uu", i__12, "disort_", (ftnlen)770)]);
	pass1 = FALSE_;
	goto L10;
    }
    return 0;
} /* disort_ */

/* Subroutine */ int asymtx_(doublereal *aa, doublereal *evec, doublereal *
	eval, integer *m, integer *ia, integer *ievec, integer *ier, 
	doublereal *wkd, doublereal *aad, doublereal *evecd, doublereal *
	evald)
{
    /* Initialized data */

    static doublereal c1 = .4375;
    static doublereal c2 = .5;
    static doublereal c3 = .75;
    static doublereal c4 = .95;
    static doublereal c5 = 16.;
    static doublereal c6 = 256.;
    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* System generated locals */
    integer aa_dim1, aa_dim2, aa_offset, eval_dim1, evec_dim1, evec_dim2, 
	    evec_offset, aad_dim1, aad_dim2, aad_offset, evald_dim1, 
	    evecd_dim1, evecd_dim2, evecd_offset, i__1, i__2, i__3, i__4, 
	    i__5;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal repl, f, g, h__;
    static integer i__, j, k, l, n;
    static doublereal p, q, r__, s, t, scale, w, x, y, z__, rnorm;
    extern doublereal d1mach_(integer *);
    static integer n1, n2, ka, lb, ii, in;
    static doublereal uu, vv, discri;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static logical notlas, noconv;
    static doublereal col;
    static integer kkk, lll;
    static doublereal sgn, tol, row;

/*    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ====== */
/*       Solves eigenfunction problem for real asymmetric matrix */
/*       for which it is known a priori that the eigenvalues are real. */
/*       This is an adaptation of a subroutine EIGRF in the IMSL */
/*       library to use real instead of complex arithmetic, accounting */
/*       for the known fact that the eigenvalues and eigenvectors in */
/*       the discrete ordinate solution are real.  Other changes include */
/*       putting all the called subroutines in-line, deleting the */
/*       performance index calculation, updating many DO-loops */
/*       to Fortran77, and in calculating the machine precision */
/*       TOL instead of specifying it in a data statement. */
/*       EIGRF is based primarily on EISPACK routines.  The matrix is */
/*       first balanced using the Parlett-Reinsch algorithm.  Then */
/*       the Martin-Wilkinson algorithm is applied. */
/*       There is a statement 'J  = WKD( I )' that converts a double */
/*       precision variable to an integer variable, that seems dangerous */
/*       to us in principle, but seems to work fine in practice. */
/*       References: */
/*          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving */
/*             Matrix Eigenvalue Problems, in Cowell, ed., 1984: */
/*             Sources and Development of Mathematical Software, */
/*             Prentice-Hall, Englewood Cliffs, NJ */
/*         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation */
/*             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304 */
/*         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem, */
/*             Clarendon Press, Oxford */

/*   I N P U T    V A R I A B L E S: */

/*       AA    :  input asymmetric matrix, destroyed after solved */
/*        M    :  order of  AA */
/*       IA    :  first dimension of  AA */
/*    IEVEC    :  first dimension of  EVEC */

/*   O U T P U T    V A R I A B L E S: */

/*       EVEC  :  (unnormalized) eigenvectors of  AA */
/*                ( column J corresponds to EVAL(J) ) */

/*       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M ) */

/*       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge; */
/*                in that case eigenvalues IER+1,IER+2,...,M  are */
/*                correct but eigenvalues 1,...,IER are set to zero. */

/*   S C R A T C H   V A R I A B L E S: */

/*       WKD   :  work area ( dimension at least 2*M ) */
/*       AAD   :  double precision stand-in for AA */
/*       EVECD :  double precision stand-in for EVEC */
/*       EVALD :  double precision stand-in for EVAL */

/*   Called by- SOLEIG */
/*   Calls- D1MACH, ERRMSG */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    evald_dim1 = *m;
    eval_dim1 = *m;
    evecd_dim1 = *ia;
    evecd_dim2 = *m;
    evecd_offset = 1 + evecd_dim1 * 1;
    aad_dim1 = *ia;
    aad_dim2 = *m;
    aad_offset = 1 + aad_dim1 * 1;
    aa_dim1 = *ia;
    aa_dim2 = *m;
    aa_offset = 1 + aa_dim1 * 1;
    evec_dim1 = *ievec;
    evec_dim2 = *m;
    evec_offset = 1 + evec_dim1 * 1;

    /* Function Body */
    *ier = 0;
    tol = d1mach_(&c__4);
    if (*m < 1 || *ia < *m || *ievec < *m) {
	errmsg_("ASYMTX--bad input variable(s)", &c_true, (ftnlen)29);
    }
/*                           ** Handle 1x1 and 2x2 special cases */
    if (*m == 1) {
	eval[(i__1 = 0) < 1 * eval_dim1 ? i__1 : s_rnge("eval", i__1, "asymt\
x_", (ftnlen)896)] = aa[(i__2 = aa_dim1 + 1 - aa_offset) < 1 * aa_dim1 * 
		aa_dim2 && 0 <= i__2 ? i__2 : s_rnge("aa", i__2, "asymtx_", (
		ftnlen)896)];
	evec[(i__1 = evec_dim1 + 1 - evec_offset) < 1 * evec_dim1 * evec_dim2 
		&& 0 <= i__1 ? i__1 : s_rnge("evec", i__1, "asymtx_", (ftnlen)
		897)] = 1.;
	return 0;
    } else if (*m == 2) {
/* Computing 2nd power */
	d__1 = aa[(i__1 = aa_dim1 + 1 - aa_offset) < 1 * aa_dim1 * aa_dim2 && 
		0 <= i__1 ? i__1 : s_rnge("aa", i__1, "asymtx_", (ftnlen)902)]
		 - aa[(i__2 = (aa_dim1 << 1) + 2 - aa_offset) < 1 * aa_dim1 * 
		aa_dim2 && 0 <= i__2 ? i__2 : s_rnge("aa", i__2, "asymtx_", (
		ftnlen)902)];
	discri = d__1 * d__1 + aa[(i__3 = (aa_dim1 << 1) + 1 - aa_offset) < 1 
		* aa_dim1 * aa_dim2 && 0 <= i__3 ? i__3 : s_rnge("aa", i__3, 
		"asymtx_", (ftnlen)902)] * 4. * aa[(i__4 = aa_dim1 + 2 - 
		aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__4 ? i__4 : 
		s_rnge("aa", i__4, "asymtx_", (ftnlen)902)];
	if (discri < 0.) {
	    errmsg_("ASYMTX--complex evals in 2x2 case", &c_true, (ftnlen)33);
	}
	sgn = one;
	if (aa[(i__1 = aa_dim1 + 1 - aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 
		<= i__1 ? i__1 : s_rnge("aa", i__1, "asymtx_", (ftnlen)909)] <
		 aa[(i__2 = (aa_dim1 << 1) + 2 - aa_offset) < 1 * aa_dim1 * 
		aa_dim2 && 0 <= i__2 ? i__2 : s_rnge("aa", i__2, "asymtx_", (
		ftnlen)909)]) {
	    sgn = -one;
	}
	eval[(i__1 = 0) < 1 * eval_dim1 ? i__1 : s_rnge("eval", i__1, "asymt\
x_", (ftnlen)911)] = (aa[(i__2 = aa_dim1 + 1 - aa_offset) < 1 * aa_dim1 * 
		aa_dim2 && 0 <= i__2 ? i__2 : s_rnge("aa", i__2, "asymtx_", (
		ftnlen)911)] + aa[(i__3 = (aa_dim1 << 1) + 2 - aa_offset) < 1 
		* aa_dim1 * aa_dim2 && 0 <= i__3 ? i__3 : s_rnge("aa", i__3, 
		"asymtx_", (ftnlen)911)] + sgn * sqrt(discri)) * .5;
	eval[(i__1 = 1) < 1 * eval_dim1 ? i__1 : s_rnge("eval", i__1, "asymt\
x_", (ftnlen)912)] = (aa[(i__2 = aa_dim1 + 1 - aa_offset) < 1 * aa_dim1 * 
		aa_dim2 && 0 <= i__2 ? i__2 : s_rnge("aa", i__2, "asymtx_", (
		ftnlen)912)] + aa[(i__3 = (aa_dim1 << 1) + 2 - aa_offset) < 1 
		* aa_dim1 * aa_dim2 && 0 <= i__3 ? i__3 : s_rnge("aa", i__3, 
		"asymtx_", (ftnlen)912)] - sgn * sqrt(discri)) * .5;
	evec[(i__1 = evec_dim1 + 1 - evec_offset) < 1 * evec_dim1 * evec_dim2 
		&& 0 <= i__1 ? i__1 : s_rnge("evec", i__1, "asymtx_", (ftnlen)
		913)] = 1.;
	evec[(i__1 = (evec_dim1 << 1) + 2 - evec_offset) < 1 * evec_dim1 * 
		evec_dim2 && 0 <= i__1 ? i__1 : s_rnge("evec", i__1, "asymtx_"
		, (ftnlen)914)] = 1.;
	if (aa[(i__1 = aa_dim1 + 1 - aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 
		<= i__1 ? i__1 : s_rnge("aa", i__1, "asymtx_", (ftnlen)916)] 
		== aa[(i__2 = (aa_dim1 << 1) + 2 - aa_offset) < 1 * aa_dim1 * 
		aa_dim2 && 0 <= i__2 ? i__2 : s_rnge("aa", i__2, "asymtx_", (
		ftnlen)916)] && (aa[(i__3 = aa_dim1 + 2 - aa_offset) < 1 * 
		aa_dim1 * aa_dim2 && 0 <= i__3 ? i__3 : s_rnge("aa", i__3, 
		"asymtx_", (ftnlen)916)] == 0. || aa[(i__4 = (aa_dim1 << 1) + 
		1 - aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__4 ? i__4 : 
		s_rnge("aa", i__4, "asymtx_", (ftnlen)916)] == 0.)) {
	    rnorm = (d__1 = aa[(i__1 = aa_dim1 + 1 - aa_offset) < 1 * aa_dim1 
		    * aa_dim2 && 0 <= i__1 ? i__1 : s_rnge("aa", i__1, "asym\
tx_", (ftnlen)919)], abs(d__1)) + (d__2 = aa[(i__2 = (aa_dim1 << 1) + 1 - 
		    aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__2 ? i__2 : 
		    s_rnge("aa", i__2, "asymtx_", (ftnlen)919)], abs(d__2)) + 
		    (d__3 = aa[(i__3 = aa_dim1 + 2 - aa_offset) < 1 * aa_dim1 
		    * aa_dim2 && 0 <= i__3 ? i__3 : s_rnge("aa", i__3, "asym\
tx_", (ftnlen)919)], abs(d__3)) + (d__4 = aa[(i__4 = (aa_dim1 << 1) + 2 - 
		    aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__4 ? i__4 : 
		    s_rnge("aa", i__4, "asymtx_", (ftnlen)919)], abs(d__4));
	    w = tol * rnorm;
	    evec[(i__1 = evec_dim1 + 2 - evec_offset) < 1 * evec_dim1 * 
		    evec_dim2 && 0 <= i__1 ? i__1 : s_rnge("evec", i__1, 
		    "asymtx_", (ftnlen)922)] = aa[(i__2 = aa_dim1 + 2 - 
		    aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__2 ? i__2 : 
		    s_rnge("aa", i__2, "asymtx_", (ftnlen)922)] / w;
	    evec[(i__1 = (evec_dim1 << 1) + 1 - evec_offset) < 1 * evec_dim1 *
		     evec_dim2 && 0 <= i__1 ? i__1 : s_rnge("evec", i__1, 
		    "asymtx_", (ftnlen)923)] = -aa[(i__2 = (aa_dim1 << 1) + 1 
		    - aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__2 ? i__2 :
		     s_rnge("aa", i__2, "asymtx_", (ftnlen)923)] / w;
	} else {
	    evec[(i__1 = evec_dim1 + 2 - evec_offset) < 1 * evec_dim1 * 
		    evec_dim2 && 0 <= i__1 ? i__1 : s_rnge("evec", i__1, 
		    "asymtx_", (ftnlen)927)] = aa[(i__2 = aa_dim1 + 2 - 
		    aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__2 ? i__2 : 
		    s_rnge("aa", i__2, "asymtx_", (ftnlen)927)] / (eval[(i__3 
		    = 0) < 1 * eval_dim1 ? i__3 : s_rnge("eval", i__3, "asym\
tx_", (ftnlen)927)] - aa[(i__4 = (aa_dim1 << 1) + 2 - aa_offset) < 1 * 
		    aa_dim1 * aa_dim2 && 0 <= i__4 ? i__4 : s_rnge("aa", i__4,
		     "asymtx_", (ftnlen)927)]);
	    evec[(i__1 = (evec_dim1 << 1) + 1 - evec_offset) < 1 * evec_dim1 *
		     evec_dim2 && 0 <= i__1 ? i__1 : s_rnge("evec", i__1, 
		    "asymtx_", (ftnlen)928)] = aa[(i__2 = (aa_dim1 << 1) + 1 
		    - aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__2 ? i__2 :
		     s_rnge("aa", i__2, "asymtx_", (ftnlen)928)] / (eval[(
		    i__3 = 1) < 1 * eval_dim1 ? i__3 : s_rnge("eval", i__3, 
		    "asymtx_", (ftnlen)928)] - aa[(i__4 = aa_dim1 + 1 - 
		    aa_offset) < 1 * aa_dim1 * aa_dim2 && 0 <= i__4 ? i__4 : 
		    s_rnge("aa", i__4, "asymtx_", (ftnlen)928)]);
	}
	return 0;
    }
/*                               ** Convert single-prec. matrix to double */
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (k = 1; k <= i__2; ++k) {
	    aad[(i__3 = j + k * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, "asym\
tx_", (ftnlen)940)] = aa[(i__4 = j + k * aa_dim1 - aa_offset) < 1 * aa_dim1 * 
		    aa_dim2 && 0 <= i__4 ? i__4 : s_rnge("aa", i__4, "asymtx_"
		    , (ftnlen)940)];
/* L10: */
	}
/* L20: */
    }
/*                                ** Initialize output variables */
    *ier = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	evald[(i__2 = i__ - 1) < 1 * evald_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"evald", i__2, "asymtx_", (ftnlen)949)] = zero;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    evecd[(i__3 = i__ + j * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecd", i__3, "asymtx_", (ftnlen)952)] = zero;
/* L30: */
	}
	evecd[(i__2 = i__ + i__ * evecd_dim1 - evecd_offset) < 1 * evecd_dim1 
		* evecd_dim2 && 0 <= i__2 ? i__2 : s_rnge("evecd", i__2, 
		"asymtx_", (ftnlen)955)] = one;
/* L40: */
    }
/*                  ** Balance the input matrix and reduce its norm by */
/*                  ** diagonal similarity transformation stored in WK; */
/*                  ** then search for rows isolating an eigenvalue */
/*                  ** and push them down */
    rnorm = zero;
    l = 1;
    k = *m;
L50:
    kkk = k;
    for (j = kkk; j >= 1; --j) {
	row = zero;
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ != j) {
		row += (d__1 = aad[(i__2 = j + i__ * aad_dim1 - aad_offset) < 
			1 * aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			"aad", i__2, "asymtx_", (ftnlen)975)], abs(d__1));
	    }
/* L60: */
	}
	if (row == zero) {
	    wkd[k - 1] = (doublereal) j;
	    if (j != k) {
		i__1 = k;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    repl = aad[(i__2 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			    "aad", i__2, "asymtx_", (ftnlen)986)];
		    aad[(i__2 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			    "aad", i__2, "asymtx_", (ftnlen)987)] = aad[(i__3 
			    = i__ + k * aad_dim1 - aad_offset) < 1 * aad_dim1 
			    * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", 
			    i__3, "asymtx_", (ftnlen)987)];
		    aad[(i__2 = i__ + k * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			    "aad", i__2, "asymtx_", (ftnlen)988)] = repl;
/* L70: */
		}
		i__1 = *m;
		for (i__ = l; i__ <= i__1; ++i__) {
		    repl = aad[(i__2 = j + i__ * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			    "aad", i__2, "asymtx_", (ftnlen)992)];
		    aad[(i__2 = j + i__ * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			    "aad", i__2, "asymtx_", (ftnlen)993)] = aad[(i__3 
			    = k + i__ * aad_dim1 - aad_offset) < 1 * aad_dim1 
			    * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", 
			    i__3, "asymtx_", (ftnlen)993)];
		    aad[(i__2 = k + i__ * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			    "aad", i__2, "asymtx_", (ftnlen)994)] = repl;
/* L80: */
		}
	    }
	    --k;
	    goto L50;
	}
/* L90: */
    }
/*                                ** Search for columns isolating an */
/*                                ** eigenvalue and push them left */
L100:
    lll = l;
    i__1 = k;
    for (j = lll; j <= i__1; ++j) {
	col = zero;
	i__2 = k;
	for (i__ = l; i__ <= i__2; ++i__) {
	    if (i__ != j) {
		col += (d__1 = aad[(i__3 = i__ + j * aad_dim1 - aad_offset) < 
			1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"aad", i__3, "asymtx_", (ftnlen)1016)], abs(d__1));
	    }
/* L110: */
	}
	if (col == zero) {
	    wkd[l - 1] = (doublereal) j;
	    if (j != l) {
		i__2 = k;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    repl = aad[(i__3 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			    "aad", i__3, "asymtx_", (ftnlen)1027)];
		    aad[(i__3 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			    "aad", i__3, "asymtx_", (ftnlen)1028)] = aad[(
			    i__4 = i__ + l * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			    "aad", i__4, "asymtx_", (ftnlen)1028)];
		    aad[(i__3 = i__ + l * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			    "aad", i__3, "asymtx_", (ftnlen)1029)] = repl;
/* L120: */
		}
		i__2 = *m;
		for (i__ = l; i__ <= i__2; ++i__) {
		    repl = aad[(i__3 = j + i__ * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			    "aad", i__3, "asymtx_", (ftnlen)1033)];
		    aad[(i__3 = j + i__ * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			    "aad", i__3, "asymtx_", (ftnlen)1034)] = aad[(
			    i__4 = l + i__ * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			    "aad", i__4, "asymtx_", (ftnlen)1034)];
		    aad[(i__3 = l + i__ * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			    "aad", i__3, "asymtx_", (ftnlen)1035)] = repl;
/* L130: */
		}
	    }
	    ++l;
	    goto L100;
	}
/* L140: */
    }
/*                           ** Balance the submatrix in rows L through K */
    i__1 = k;
    for (i__ = l; i__ <= i__1; ++i__) {
	wkd[i__ - 1] = one;
/* L150: */
    }
L160:
    noconv = FALSE_;
    i__1 = k;
    for (i__ = l; i__ <= i__1; ++i__) {
	col = zero;
	row = zero;
	i__2 = k;
	for (j = l; j <= i__2; ++j) {
	    if (j != i__) {
		col += (d__1 = aad[(i__3 = j + i__ * aad_dim1 - aad_offset) < 
			1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"aad", i__3, "asymtx_", (ftnlen)1064)], abs(d__1));
		row += (d__1 = aad[(i__3 = i__ + j * aad_dim1 - aad_offset) < 
			1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"aad", i__3, "asymtx_", (ftnlen)1065)], abs(d__1));
	    }
/* L170: */
	}
	f = one;
	g = row / c5;
	h__ = col + row;
L180:
	if (col < g) {
	    f *= c5;
	    col *= c6;
	    goto L180;
	}
	g = row * c5;
L190:
	if (col >= g) {
	    f /= c5;
	    col /= c6;
	    goto L190;
	}
/*                                                ** Now balance */
	if ((col + row) / f < c4 * h__) {
	    wkd[i__ - 1] *= f;
	    noconv = TRUE_;
	    i__2 = *m;
	    for (j = l; j <= i__2; ++j) {
		aad[(i__3 = i__ + j * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
			aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, 
			"asymtx_", (ftnlen)1101)] = aad[(i__4 = i__ + j * 
			aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 
			<= i__4 ? i__4 : s_rnge("aad", i__4, "asymtx_", (
			ftnlen)1101)] / f;
/* L200: */
	    }
	    i__2 = k;
	    for (j = 1; j <= i__2; ++j) {
		aad[(i__3 = j + i__ * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
			aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, 
			"asymtx_", (ftnlen)1105)] = aad[(i__4 = j + i__ * 
			aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 
			<= i__4 ? i__4 : s_rnge("aad", i__4, "asymtx_", (
			ftnlen)1105)] * f;
/* L210: */
	    }
	}
/* L220: */
    }
    if (noconv) {
	goto L160;
    }
/*                                   ** Is A already in Hessenberg form? */
    if (k - 1 < l + 1) {
	goto L370;
    }
/*                                   ** Transfer A to a Hessenberg form */
    i__1 = k - 1;
    for (n = l + 1; n <= i__1; ++n) {
	h__ = zero;
	wkd[n + *m - 1] = zero;
	scale = zero;
/*                                                 ** Scale column */
	i__2 = k;
	for (i__ = n; i__ <= i__2; ++i__) {
	    scale += (d__1 = aad[(i__3 = i__ + (n - 1) * aad_dim1 - 
		    aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 
		    : s_rnge("aad", i__3, "asymtx_", (ftnlen)1125)], abs(d__1)
		    );
/* L230: */
	}
	if (scale != zero) {
	    i__2 = n;
	    for (i__ = k; i__ >= i__2; --i__) {
		wkd[i__ + *m - 1] = aad[(i__3 = i__ + (n - 1) * aad_dim1 - 
			aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? 
			i__3 : s_rnge("aad", i__3, "asymtx_", (ftnlen)1131)] /
			 scale;
/* Computing 2nd power */
		d__1 = wkd[i__ + *m - 1];
		h__ += d__1 * d__1;
/* L240: */
	    }
	    d__1 = sqrt(h__);
	    g = -d_sign(&d__1, &wkd[n + *m - 1]);
	    h__ -= wkd[n + *m - 1] * g;
	    wkd[n + *m - 1] -= g;
/*                                            ** Form (I-(U*UT)/H)*A */
	    i__2 = *m;
	    for (j = n; j <= i__2; ++j) {
		f = zero;
		i__3 = n;
		for (i__ = k; i__ >= i__3; --i__) {
		    f += wkd[i__ + *m - 1] * aad[(i__4 = i__ + j * aad_dim1 - 
			    aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= 
			    i__4 ? i__4 : s_rnge("aad", i__4, "asymtx_", (
			    ftnlen)1144)];
/* L250: */
		}
		i__3 = k;
		for (i__ = n; i__ <= i__3; ++i__) {
		    aad[(i__4 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			    "aad", i__4, "asymtx_", (ftnlen)1148)] = aad[(
			    i__5 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__5 ? i__5 : s_rnge(
			    "aad", i__5, "asymtx_", (ftnlen)1148)] - wkd[i__ 
			    + *m - 1] * f / h__;
/* L260: */
		}
/* L270: */
	    }
/*                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H) */
	    i__2 = k;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		f = zero;
		i__3 = n;
		for (j = k; j >= i__3; --j) {
		    f += wkd[j + *m - 1] * aad[(i__4 = i__ + j * aad_dim1 - 
			    aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= 
			    i__4 ? i__4 : s_rnge("aad", i__4, "asymtx_", (
			    ftnlen)1158)];
/* L280: */
		}
		i__3 = k;
		for (j = n; j <= i__3; ++j) {
		    aad[(i__4 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			    "aad", i__4, "asymtx_", (ftnlen)1162)] = aad[(
			    i__5 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__5 ? i__5 : s_rnge(
			    "aad", i__5, "asymtx_", (ftnlen)1162)] - wkd[j + *
			    m - 1] * f / h__;
/* L290: */
		}
/* L300: */
	    }
	    wkd[n + *m - 1] = scale * wkd[n + *m - 1];
	    aad[(i__2 = n + (n - 1) * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asym\
tx_", (ftnlen)1168)] = scale * g;
	}
/* L310: */
    }
    i__1 = l;
    for (n = k - 2; n >= i__1; --n) {
	n1 = n + 1;
	n2 = n + 2;
	f = aad[(i__2 = n + 1 + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asymtx_", 
		(ftnlen)1179)];
	if (f != zero) {
	    f *= wkd[n + 1 + *m - 1];
	    i__2 = k;
	    for (i__ = n + 2; i__ <= i__2; ++i__) {
		wkd[i__ + *m - 1] = aad[(i__3 = i__ + n * aad_dim1 - 
			aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? 
			i__3 : s_rnge("aad", i__3, "asymtx_", (ftnlen)1186)];
/* L320: */
	    }
	    if (n + 1 <= k) {
		i__2 = *m;
		for (j = 1; j <= i__2; ++j) {
		    g = zero;
		    i__3 = k;
		    for (i__ = n + 1; i__ <= i__3; ++i__) {
			g += wkd[i__ + *m - 1] * evecd[(i__4 = i__ + j * 
				evecd_dim1 - evecd_offset) < 1 * evecd_dim1 * 
				evecd_dim2 && 0 <= i__4 ? i__4 : s_rnge("eve\
cd", i__4, "asymtx_", (ftnlen)1196)];
/* L330: */
		    }
		    g /= f;
		    i__3 = k;
		    for (i__ = n + 1; i__ <= i__3; ++i__) {
			evecd[(i__4 = i__ + j * evecd_dim1 - evecd_offset) < 
				1 * evecd_dim1 * evecd_dim2 && 0 <= i__4 ? 
				i__4 : s_rnge("evecd", i__4, "asymtx_", (
				ftnlen)1202)] = evecd[(i__5 = i__ + j * 
				evecd_dim1 - evecd_offset) < 1 * evecd_dim1 * 
				evecd_dim2 && 0 <= i__5 ? i__5 : s_rnge("eve\
cd", i__5, "asymtx_", (ftnlen)1202)] + g * wkd[i__ + *m - 1];
/* L340: */
		    }
/* L350: */
		}
	    }
	}
/* L360: */
    }
L370:
    n = 1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = n; j <= i__2; ++j) {
	    rnorm += (d__1 = aad[(i__3 = i__ + j * aad_dim1 - aad_offset) < 1 
		    * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", 
		    i__3, "asymtx_", (ftnlen)1221)], abs(d__1));
/* L380: */
	}
	n = i__;
	if (i__ < l || i__ > k) {
	    evald[(i__2 = i__ - 1) < 1 * evald_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("evald", i__2, "asymtx_", (ftnlen)1226)] = aad[(
		    i__3 = i__ + i__ * aad_dim1 - aad_offset) < 1 * aad_dim1 *
		     aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, "asy\
mtx_", (ftnlen)1226)];
	}
/* L390: */
    }
    n = k;
    t = zero;
/*                                      ** Search for next eigenvalues */
L400:
    if (n < l) {
	goto L550;
    }
    in = 0;
    n1 = n - 1;
    n2 = n - 2;
/*                          ** Look for single small sub-diagonal element */
L410:
    i__1 = n;
    for (i__ = l; i__ <= i__1; ++i__) {
	lb = n + l - i__;
	if (lb == l) {
	    goto L430;
	}
	s = (d__1 = aad[(i__2 = lb - 1 + (lb - 1) * aad_dim1 - aad_offset) < 
		1 * aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", 
		i__2, "asymtx_", (ftnlen)1248)], abs(d__1)) + (d__2 = aad[(
		i__3 = lb + lb * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, "asymtx_", 
		(ftnlen)1248)], abs(d__2));
	if (s == zero) {
	    s = rnorm;
	}
	if ((d__1 = aad[(i__2 = lb + (lb - 1) * aad_dim1 - aad_offset) < 1 * 
		aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, 
		"asymtx_", (ftnlen)1252)], abs(d__1)) <= tol * s) {
	    goto L430;
	}
/* L420: */
    }
L430:
    x = aad[(i__1 = n + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 
	    && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", (ftnlen)1258)
	    ];
    if (lb == n) {
/*                                        ** One eigenvalue found */
	aad[(i__1 = n + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 
		&& 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", (ftnlen)
		1262)] = x + t;
	evald[(i__1 = n - 1) < 1 * evald_dim1 && 0 <= i__1 ? i__1 : s_rnge(
		"evald", i__1, "asymtx_", (ftnlen)1263)] = aad[(i__2 = n + n *
		 aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= 
		i__2 ? i__2 : s_rnge("aad", i__2, "asymtx_", (ftnlen)1263)];
	n = n1;
	goto L400;
    }
    y = aad[(i__1 = n1 + n1 * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
	    aad_dim2 && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", (
	    ftnlen)1269)];
    w = aad[(i__1 = n + n1 * aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 
	    && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", (ftnlen)1270)
	    ] * aad[(i__2 = n1 + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
	    aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asymtx_", (
	    ftnlen)1270)];
    if (lb == n1) {
/*                                        ** Two eigenvalues found */
	p = (y - x) * c2;
/* Computing 2nd power */
	d__1 = p;
	q = d__1 * d__1 + w;
	z__ = sqrt((abs(q)));
	aad[(i__1 = n + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 
		&& 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", (ftnlen)
		1277)] = x + t;
	x = aad[(i__1 = n + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		aad_dim2 && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", 
		(ftnlen)1278)];
	aad[(i__1 = n1 + n1 * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		aad_dim2 && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", 
		(ftnlen)1279)] = y + t;
/*                                        ** Real pair */
	z__ = p + d_sign(&z__, &p);
	evald[(i__1 = n1 - 1) < 1 * evald_dim1 && 0 <= i__1 ? i__1 : s_rnge(
		"evald", i__1, "asymtx_", (ftnlen)1282)] = x + z__;
	evald[(i__1 = n - 1) < 1 * evald_dim1 && 0 <= i__1 ? i__1 : s_rnge(
		"evald", i__1, "asymtx_", (ftnlen)1283)] = evald[(i__2 = n1 - 
		1) < 1 * evald_dim1 && 0 <= i__2 ? i__2 : s_rnge("evald", 
		i__2, "asymtx_", (ftnlen)1283)];
	if (z__ != zero) {
	    evald[(i__1 = n - 1) < 1 * evald_dim1 && 0 <= i__1 ? i__1 : 
		    s_rnge("evald", i__1, "asymtx_", (ftnlen)1285)] = x - w / 
		    z__;
	}
	x = aad[(i__1 = n + n1 * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		aad_dim2 && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", 
		(ftnlen)1287)];
/*                                  ** Employ scale factor in case */
/*                                  ** X and Z are very small */
	r__ = sqrt(x * x + z__ * z__);
	p = x / r__;
	q = z__ / r__;
/*                                             ** Row modification */
	i__1 = *m;
	for (j = n1; j <= i__1; ++j) {
	    z__ = aad[(i__2 = n1 + j * aad_dim1 - aad_offset) < 1 * aad_dim1 *
		     aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asy\
mtx_", (ftnlen)1295)];
	    aad[(i__2 = n1 + j * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asym\
tx_", (ftnlen)1296)] = q * z__ + p * aad[(i__3 = n + j * aad_dim1 - 
		    aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 
		    : s_rnge("aad", i__3, "asymtx_", (ftnlen)1296)];
	    aad[(i__2 = n + j * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asym\
tx_", (ftnlen)1297)] = q * aad[(i__3 = n + j * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", 
		    i__3, "asymtx_", (ftnlen)1297)] - p * z__;
/* L440: */
	}
/*                                             ** Column modification */
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__ = aad[(i__2 = i__ + n1 * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", 
		    i__2, "asymtx_", (ftnlen)1301)];
	    aad[(i__2 = i__ + n1 * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asym\
tx_", (ftnlen)1302)] = q * z__ + p * aad[(i__3 = i__ + n * aad_dim1 - 
		    aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 
		    : s_rnge("aad", i__3, "asymtx_", (ftnlen)1302)];
	    aad[(i__2 = i__ + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asym\
tx_", (ftnlen)1303)] = q * aad[(i__3 = i__ + n * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", 
		    i__3, "asymtx_", (ftnlen)1303)] - p * z__;
/* L450: */
	}
/*                                          ** Accumulate transformations */
	i__1 = k;
	for (i__ = l; i__ <= i__1; ++i__) {
	    z__ = evecd[(i__2 = i__ + n1 * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__2 ? i__2 : s_rnge(
		    "evecd", i__2, "asymtx_", (ftnlen)1307)];
	    evecd[(i__2 = i__ + n1 * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__2 ? i__2 : s_rnge(
		    "evecd", i__2, "asymtx_", (ftnlen)1308)] = q * z__ + p * 
		    evecd[(i__3 = i__ + n * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecd", i__3, "asymtx_", (ftnlen)1308)];
	    evecd[(i__2 = i__ + n * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__2 ? i__2 : s_rnge(
		    "evecd", i__2, "asymtx_", (ftnlen)1309)] = q * evecd[(
		    i__3 = i__ + n * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecd", i__3, "asymtx_", (ftnlen)1309)] - p * z__;
/* L460: */
	}
	n = n2;
	goto L400;
    }
    if (in == 30) {
/*                    ** No convergence after 30 iterations; set error */
/*                    ** indicator to the index of the current eigenvalue */
	*ier = n;
	goto L700;
    }
/*                                                  ** Form shift */
    if (in == 10 || in == 20) {
	t += x;
	i__1 = n;
	for (i__ = l; i__ <= i__1; ++i__) {
	    aad[(i__2 = i__ + i__ * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asym\
tx_", (ftnlen)1332)] = aad[(i__3 = i__ + i__ * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", 
		    i__3, "asymtx_", (ftnlen)1332)] - x;
/* L470: */
	}
	s = (d__1 = aad[(i__1 = n + n1 * aad_dim1 - aad_offset) < 1 * 
		aad_dim1 * aad_dim2 && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, 
		"asymtx_", (ftnlen)1335)], abs(d__1)) + (d__2 = aad[(i__2 = 
		n1 + n2 * aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 && 
		0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asymtx_", (ftnlen)
		1335)], abs(d__2));
	x = c3 * s;
	y = x;
/* Computing 2nd power */
	d__1 = s;
	w = -c1 * (d__1 * d__1);
    }
    ++in;
/*                ** Look for two consecutive small sub-diagonal elements */
    i__1 = n2;
    for (j = lb; j <= i__1; ++j) {
	i__ = n2 + lb - j;
	z__ = aad[(i__2 = i__ + i__ * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asymtx_", 
		(ftnlen)1349)];
	r__ = x - z__;
	s = y - z__;
	p = (r__ * s - w) / aad[(i__2 = i__ + 1 + i__ * aad_dim1 - aad_offset)
		 < 1 * aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad",
		 i__2, "asymtx_", (ftnlen)1352)] + aad[(i__3 = i__ + (i__ + 1)
		 * aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= 
		i__3 ? i__3 : s_rnge("aad", i__3, "asymtx_", (ftnlen)1352)];
	q = aad[(i__2 = i__ + 1 + (i__ + 1) * aad_dim1 - aad_offset) < 1 * 
		aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, 
		"asymtx_", (ftnlen)1353)] - z__ - r__ - s;
	r__ = aad[(i__2 = i__ + 2 + (i__ + 1) * aad_dim1 - aad_offset) < 1 * 
		aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, 
		"asymtx_", (ftnlen)1354)];
	s = abs(p) + abs(q) + abs(r__);
	p /= s;
	q /= s;
	r__ /= s;
	if (i__ == lb) {
	    goto L490;
	}
	uu = (d__1 = aad[(i__2 = i__ + (i__ - 1) * aad_dim1 - aad_offset) < 1 
		* aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", 
		i__2, "asymtx_", (ftnlen)1362)], abs(d__1)) * (abs(q) + abs(
		r__));
	vv = abs(p) * ((d__1 = aad[(i__2 = i__ - 1 + (i__ - 1) * aad_dim1 - 
		aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : 
		s_rnge("aad", i__2, "asymtx_", (ftnlen)1363)], abs(d__1)) + 
		abs(z__) + (d__2 = aad[(i__3 = i__ + 1 + (i__ + 1) * aad_dim1 
		- aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : 
		s_rnge("aad", i__3, "asymtx_", (ftnlen)1363)], abs(d__2)));
	if (uu <= tol * vv) {
	    goto L490;
	}
/* L480: */
    }
L490:
    aad[(i__1 = i__ + 2 + i__ * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
	    aad_dim2 && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asymtx_", (
	    ftnlen)1371)] = zero;
    i__1 = n;
    for (j = i__ + 3; j <= i__1; ++j) {
	aad[(i__2 = j + (j - 2) * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asymtx_", 
		(ftnlen)1374)] = zero;
	aad[(i__2 = j + (j - 3) * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, "asymtx_", 
		(ftnlen)1375)] = zero;
/* L500: */
    }
/*             ** Double QR step involving rows K to N and columns M to N */
    i__1 = n1;
    for (ka = i__; ka <= i__1; ++ka) {
	notlas = ka != n1;
	if (ka == i__) {
	    d__1 = sqrt(p * p + q * q + r__ * r__);
	    s = d_sign(&d__1, &p);
	    if (lb != i__) {
		aad[(i__2 = ka + (ka - 1) * aad_dim1 - aad_offset) < 1 * 
			aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			"aad", i__2, "asymtx_", (ftnlen)1388)] = -aad[(i__3 = 
			ka + (ka - 1) * aad_dim1 - aad_offset) < 1 * aad_dim1 
			* aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, 
			"asymtx_", (ftnlen)1388)];
	    }
	} else {
	    p = aad[(i__2 = ka + (ka - 1) * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", 
		    i__2, "asymtx_", (ftnlen)1392)];
	    q = aad[(i__2 = ka + 1 + (ka - 1) * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", 
		    i__2, "asymtx_", (ftnlen)1393)];
	    r__ = zero;
	    if (notlas) {
		r__ = aad[(i__2 = ka + 2 + (ka - 1) * aad_dim1 - aad_offset) <
			 1 * aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			"aad", i__2, "asymtx_", (ftnlen)1396)];
	    }
	    x = abs(p) + abs(q) + abs(r__);
	    if (x == zero) {
		goto L540;
	    }
	    p /= x;
	    q /= x;
	    r__ /= x;
	    d__1 = sqrt(p * p + q * q + r__ * r__);
	    s = d_sign(&d__1, &p);
	    aad[(i__2 = ka + (ka - 1) * aad_dim1 - aad_offset) < 1 * aad_dim1 
		    * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge("aad", i__2, 
		    "asymtx_", (ftnlen)1406)] = -s * x;
	}
	p += s;
	x = p / s;
	y = q / s;
	z__ = r__ / s;
	q /= p;
	r__ /= p;
/*                                              ** Row modification */
	i__2 = *m;
	for (j = ka; j <= i__2; ++j) {
	    p = aad[(i__3 = ka + j * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, "asym\
tx_", (ftnlen)1419)] + q * aad[(i__4 = ka + 1 + j * aad_dim1 - aad_offset) < 
		    1 * aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge(
		    "aad", i__4, "asymtx_", (ftnlen)1419)];
	    if (notlas) {
		p += r__ * aad[(i__3 = ka + 2 + j * aad_dim1 - aad_offset) < 
			1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"aad", i__3, "asymtx_", (ftnlen)1423)];
		aad[(i__3 = ka + 2 + j * aad_dim1 - aad_offset) < 1 * 
			aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"aad", i__3, "asymtx_", (ftnlen)1424)] = aad[(i__4 = 
			ka + 2 + j * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
			aad_dim2 && 0 <= i__4 ? i__4 : s_rnge("aad", i__4, 
			"asymtx_", (ftnlen)1424)] - p * z__;
	    }
	    aad[(i__3 = ka + 1 + j * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, "asym\
tx_", (ftnlen)1428)] = aad[(i__4 = ka + 1 + j * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge("aad", 
		    i__4, "asymtx_", (ftnlen)1428)] - p * y;
	    aad[(i__3 = ka + j * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, "asym\
tx_", (ftnlen)1429)] = aad[(i__4 = ka + j * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge("aad", 
		    i__4, "asymtx_", (ftnlen)1429)] - p * x;
/* L510: */
	}
/*                                                 ** Column modification */
/* Computing MIN */
	i__3 = n, i__4 = ka + 3;
	i__2 = min(i__3,i__4);
	for (ii = 1; ii <= i__2; ++ii) {
	    p = x * aad[(i__3 = ii + ka * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", 
		    i__3, "asymtx_", (ftnlen)1434)] + y * aad[(i__4 = ii + (
		    ka + 1) * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__4 ? i__4 : s_rnge("aad", i__4, "asym\
tx_", (ftnlen)1434)];
	    if (notlas) {
		p += z__ * aad[(i__3 = ii + (ka + 2) * aad_dim1 - aad_offset) 
			< 1 * aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : 
			s_rnge("aad", i__3, "asymtx_", (ftnlen)1438)];
		aad[(i__3 = ii + (ka + 2) * aad_dim1 - aad_offset) < 1 * 
			aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"aad", i__3, "asymtx_", (ftnlen)1439)] = aad[(i__4 = 
			ii + (ka + 2) * aad_dim1 - aad_offset) < 1 * aad_dim1 
			* aad_dim2 && 0 <= i__4 ? i__4 : s_rnge("aad", i__4, 
			"asymtx_", (ftnlen)1439)] - p * r__;
	    }
	    aad[(i__3 = ii + (ka + 1) * aad_dim1 - aad_offset) < 1 * aad_dim1 
		    * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, 
		    "asymtx_", (ftnlen)1443)] = aad[(i__4 = ii + (ka + 1) * 
		    aad_dim1 - aad_offset) < 1 * aad_dim1 * aad_dim2 && 0 <= 
		    i__4 ? i__4 : s_rnge("aad", i__4, "asymtx_", (ftnlen)1443)
		    ] - p * q;
	    aad[(i__3 = ii + ka * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__3 ? i__3 : s_rnge("aad", i__3, "asym\
tx_", (ftnlen)1444)] = aad[(i__4 = ii + ka * aad_dim1 - aad_offset) < 1 * 
		    aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge("aad", 
		    i__4, "asymtx_", (ftnlen)1444)] - p;
/* L520: */
	}
/*                                          ** Accumulate transformations */
	i__2 = k;
	for (ii = l; ii <= i__2; ++ii) {
	    p = x * evecd[(i__3 = ii + ka * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecd", i__3, "asymtx_", (ftnlen)1449)] + y * evecd[(
		    i__4 = ii + (ka + 1) * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__4 ? i__4 : s_rnge(
		    "evecd", i__4, "asymtx_", (ftnlen)1449)];
	    if (notlas) {
		p += z__ * evecd[(i__3 = ii + (ka + 2) * evecd_dim1 - 
			evecd_offset) < 1 * evecd_dim1 * evecd_dim2 && 0 <= 
			i__3 ? i__3 : s_rnge("evecd", i__3, "asymtx_", (
			ftnlen)1453)];
		evecd[(i__3 = ii + (ka + 2) * evecd_dim1 - evecd_offset) < 1 *
			 evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"evecd", i__3, "asymtx_", (ftnlen)1454)] = evecd[(
			i__4 = ii + (ka + 2) * evecd_dim1 - evecd_offset) < 1 
			* evecd_dim1 * evecd_dim2 && 0 <= i__4 ? i__4 : 
			s_rnge("evecd", i__4, "asymtx_", (ftnlen)1454)] - p * 
			r__;
	    }
	    evecd[(i__3 = ii + (ka + 1) * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecd", i__3, "asymtx_", (ftnlen)1458)] = evecd[(i__4 = 
		    ii + (ka + 1) * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__4 ? i__4 : s_rnge(
		    "evecd", i__4, "asymtx_", (ftnlen)1458)] - p * q;
	    evecd[(i__3 = ii + ka * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecd", i__3, "asymtx_", (ftnlen)1459)] = evecd[(i__4 = 
		    ii + ka * evecd_dim1 - evecd_offset) < 1 * evecd_dim1 * 
		    evecd_dim2 && 0 <= i__4 ? i__4 : s_rnge("evecd", i__4, 
		    "asymtx_", (ftnlen)1459)] - p;
/* L530: */
	}
L540:
	;
    }
    goto L410;
/*                     ** All evals found, now backsubstitute real vector */
L550:
    if (rnorm != zero) {
	for (n = *m; n >= 1; --n) {
	    n2 = n;
	    aad[(i__1 = n + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
		    aad_dim2 && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, "asym\
tx_", (ftnlen)1472)] = one;
	    for (i__ = n - 1; i__ >= 1; --i__) {
		w = aad[(i__1 = i__ + i__ * aad_dim1 - aad_offset) < 1 * 
			aad_dim1 * aad_dim2 && 0 <= i__1 ? i__1 : s_rnge(
			"aad", i__1, "asymtx_", (ftnlen)1475)] - evald[(i__2 =
			 n - 1) < 1 * evald_dim1 && 0 <= i__2 ? i__2 : s_rnge(
			"evald", i__2, "asymtx_", (ftnlen)1475)];
		if (w == zero) {
		    w = tol * rnorm;
		}
		r__ = aad[(i__1 = i__ + n * aad_dim1 - aad_offset) < 1 * 
			aad_dim1 * aad_dim2 && 0 <= i__1 ? i__1 : s_rnge(
			"aad", i__1, "asymtx_", (ftnlen)1479)];
		i__1 = n - 1;
		for (j = n2; j <= i__1; ++j) {
		    r__ += aad[(i__2 = i__ + j * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			    "aad", i__2, "asymtx_", (ftnlen)1482)] * aad[(
			    i__3 = j + n * aad_dim1 - aad_offset) < 1 * 
			    aad_dim1 * aad_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			    "aad", i__3, "asymtx_", (ftnlen)1482)];
/* L560: */
		}
		aad[(i__1 = i__ + n * aad_dim1 - aad_offset) < 1 * aad_dim1 * 
			aad_dim2 && 0 <= i__1 ? i__1 : s_rnge("aad", i__1, 
			"asymtx_", (ftnlen)1485)] = -r__ / w;
		n2 = i__;
/* L570: */
	    }
/* L580: */
	}
/*                      ** End backsubstitution vectors of isolated evals */
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ < l || i__ > k) {
		i__2 = *m;
		for (j = i__; j <= i__2; ++j) {
		    evecd[(i__3 = i__ + j * evecd_dim1 - evecd_offset) < 1 * 
			    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : 
			    s_rnge("evecd", i__3, "asymtx_", (ftnlen)1496)] = 
			    aad[(i__4 = i__ + j * aad_dim1 - aad_offset) < 1 *
			     aad_dim1 * aad_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			    "aad", i__4, "asymtx_", (ftnlen)1496)];
/* L590: */
		}
	    }
/* L600: */
	}
/*                                   ** Multiply by transformation matrix */
	if (k != 0) {
	    i__1 = l;
	    for (j = *m; j >= i__1; --j) {
		i__2 = k;
		for (i__ = l; i__ <= i__2; ++i__) {
		    z__ = zero;
		    i__3 = min(j,k);
		    for (n = l; n <= i__3; ++n) {
			z__ += evecd[(i__4 = i__ + n * evecd_dim1 - 
				evecd_offset) < 1 * evecd_dim1 * evecd_dim2 &&
				 0 <= i__4 ? i__4 : s_rnge("evecd", i__4, 
				"asymtx_", (ftnlen)1511)] * aad[(i__5 = n + j 
				* aad_dim1 - aad_offset) < 1 * aad_dim1 * 
				aad_dim2 && 0 <= i__5 ? i__5 : s_rnge("aad", 
				i__5, "asymtx_", (ftnlen)1511)];
/* L610: */
		    }
		    evecd[(i__3 = i__ + j * evecd_dim1 - evecd_offset) < 1 * 
			    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : 
			    s_rnge("evecd", i__3, "asymtx_", (ftnlen)1514)] = 
			    z__;
/* L620: */
		}
/* L630: */
	    }
	}
    }
    i__1 = k;
    for (i__ = l; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    evecd[(i__3 = i__ + j * evecd_dim1 - evecd_offset) < 1 * 
		    evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecd", i__3, "asymtx_", (ftnlen)1527)] = evecd[(i__4 = 
		    i__ + j * evecd_dim1 - evecd_offset) < 1 * evecd_dim1 * 
		    evecd_dim2 && 0 <= i__4 ? i__4 : s_rnge("evecd", i__4, 
		    "asymtx_", (ftnlen)1527)] * wkd[i__ - 1];
/* L640: */
	}
/* L650: */
    }
/*                           ** Interchange rows if permutations occurred */
    for (i__ = l - 1; i__ >= 1; --i__) {
	j = (integer) wkd[i__ - 1];
	if (i__ != j) {
	    i__1 = *m;
	    for (n = 1; n <= i__1; ++n) {
		repl = evecd[(i__2 = i__ + n * evecd_dim1 - evecd_offset) < 1 
			* evecd_dim1 * evecd_dim2 && 0 <= i__2 ? i__2 : 
			s_rnge("evecd", i__2, "asymtx_", (ftnlen)1539)];
		evecd[(i__2 = i__ + n * evecd_dim1 - evecd_offset) < 1 * 
			evecd_dim1 * evecd_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			"evecd", i__2, "asymtx_", (ftnlen)1540)] = evecd[(
			i__3 = j + n * evecd_dim1 - evecd_offset) < 1 * 
			evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"evecd", i__3, "asymtx_", (ftnlen)1540)];
		evecd[(i__2 = j + n * evecd_dim1 - evecd_offset) < 1 * 
			evecd_dim1 * evecd_dim2 && 0 <= i__2 ? i__2 : s_rnge(
			"evecd", i__2, "asymtx_", (ftnlen)1541)] = repl;
/* L660: */
	    }
	}
/* L670: */
    }
    i__1 = *m;
    for (i__ = k + 1; i__ <= i__1; ++i__) {
	j = (integer) wkd[i__ - 1];
	if (i__ != j) {
	    i__2 = *m;
	    for (n = 1; n <= i__2; ++n) {
		repl = evecd[(i__3 = i__ + n * evecd_dim1 - evecd_offset) < 1 
			* evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : 
			s_rnge("evecd", i__3, "asymtx_", (ftnlen)1556)];
		evecd[(i__3 = i__ + n * evecd_dim1 - evecd_offset) < 1 * 
			evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"evecd", i__3, "asymtx_", (ftnlen)1557)] = evecd[(
			i__4 = j + n * evecd_dim1 - evecd_offset) < 1 * 
			evecd_dim1 * evecd_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			"evecd", i__4, "asymtx_", (ftnlen)1557)];
		evecd[(i__3 = j + n * evecd_dim1 - evecd_offset) < 1 * 
			evecd_dim1 * evecd_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"evecd", i__3, "asymtx_", (ftnlen)1558)] = repl;
/* L680: */
	    }
	}
/* L690: */
    }
/*                         ** Put results into output arrays */
L700:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	eval[(i__2 = j - 1) < 1 * eval_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"eval", i__2, "asymtx_", (ftnlen)1570)] = evald[(i__3 = j - 1)
		 < 1 * evald_dim1 && 0 <= i__3 ? i__3 : s_rnge("evald", i__3, 
		"asymtx_", (ftnlen)1570)];
	i__2 = *m;
	for (k = 1; k <= i__2; ++k) {
	    evec[(i__3 = j + k * evec_dim1 - evec_offset) < 1 * evec_dim1 * 
		    evec_dim2 && 0 <= i__3 ? i__3 : s_rnge("evec", i__3, 
		    "asymtx_", (ftnlen)1573)] = evecd[(i__4 = j + k * 
		    evecd_dim1 - evecd_offset) < 1 * evecd_dim1 * evecd_dim2 
		    && 0 <= i__4 ? i__4 : s_rnge("evecd", i__4, "asymtx_", (
		    ftnlen)1573)];
/* L710: */
	}
/* L720: */
    }
    return 0;
} /* asymtx_ */

/* Subroutine */ int cmpint_(doublereal *fbeam, doublereal *gc, doublereal *
	kk, integer *layru, doublereal *ll, logical *lyrcut, integer *mazim, 
	integer *mxcmu, integer *mxulv, integer *mxumu, integer *ncut, 
	integer *nn, integer *nstr, logical *plank, integer *ntau, doublereal 
	*taucpr, doublereal *umu0, doublereal *utaupr, doublereal *zz, 
	doublereal *zplk0, doublereal *zplk1, doublereal *uum)
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, utaupr_dim1, uum_dim1, uum_dim2, uum_offset, 
	    zplk0_dim1, zplk0_offset, zplk1_dim1, zplk1_offset, zz_dim1, 
	    zz_offset, i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    double exp(doublereal);

    /* Local variables */
    static doublereal zint;
    static integer iq, jq, lu, lyu;

/*          Calculates the Fourier intensity components at the quadrature */
/*          angles for azimuthal expansion terms (MAZIM) in Eq. SD(2) */


/*    I N P U T    V A R I A B L E S: */

/*       KK      :  Eigenvalues of coeff. matrix in Eq. SS(7) */

/*       GC      :  Eigenvectors at polar quadrature angles in Eq. SC(1) */

/*       LL      :  Constants of integration in Eq. SC(1), obtained */
/*                  by solving scaled version of Eq. SC(5); */
/*                  exponential term of Eq. SC(12) not included */

/*       LYRCUT  :  Logical flag for truncation of computational layer */

/*       MAZIM   :  Order of azimuthal component */

/*       NCUT    :  Number of computational layer where absorption */
/*                  optical depth exceeds ABSCUT */

/*       NN      :  Order of double-Gauss quadrature (NSTR/2) */

/*       TAUCPR  :  Cumulative optical depth (delta-M-scaled) */

/*       UTAUPR  :  Optical depths of user output levels in delta-M */
/*                  coordinates;  equal to UTAU if no delta-M */

/*       ZZ      :  Beam source vectors in Eq. SS(19) */

/*       ZPLK0   :  Thermal source vectors Z0, by solving Eq. SS(16) */

/*       ZPLK1   :  Thermal source vectors Z1, by solving Eq. SS(16) */

/*       (Remainder are 'DISORT' input variables) */


/*    O U T P U T   V A R I A B L E S: */

/*       UUM     :  Fourier components of the intensity in Eq. SD(12) */
/*                    (at polar quadrature angles) */


/*    I N T E R N A L   V A R I A B L E S: */

/*       FACT    :  EXP( - UTAUPR / UMU0 ) */
/*       ZINT    :  intensity of M=0 case, in Eq. SC(1) */

/*   Called by- DISORT */
/* +-------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*                                       ** Loop over user levels */
    /* Parameter adjustments */
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1 * 1;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1 * 1;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1 * 1;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    utaupr_dim1 = *mxulv;
    uum_dim1 = *mxumu;
    uum_dim2 = *mxulv;
    uum_offset = 1 + uum_dim1 * 1;

    /* Function Body */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	lyu = layru[lu - 1];
	if (*lyrcut && lyu > *ncut) {
	    goto L40;
	}
	i__2 = *nstr;
	for (iq = 1; iq <= i__2; ++iq) {
	    zint = 0.;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1 - gc_offset] * 
			ll[jq + lyu * ll_dim1 - ll_offset] * exp(-kk[jq + lyu 
			* kk_dim1 - kk_offset] * (utaupr[(i__4 = lu - 1) < 1 *
			 utaupr_dim1 && 0 <= i__4 ? i__4 : s_rnge("utaupr", 
			i__4, "cmpint_", (ftnlen)1671)] - taucpr[lyu]));
/* L10: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1 - gc_offset] * 
			ll[jq + lyu * ll_dim1 - ll_offset] * exp(-kk[jq + lyu 
			* kk_dim1 - kk_offset] * (utaupr[(i__4 = lu - 1) < 1 *
			 utaupr_dim1 && 0 <= i__4 ? i__4 : s_rnge("utaupr", 
			i__4, "cmpint_", (ftnlen)1677)] - taucpr[lyu - 1]));
/* L20: */
	    }
	    uum[(i__3 = iq + lu * uum_dim1 - uum_offset) < 1 * uum_dim1 * 
		    uum_dim2 && 0 <= i__3 ? i__3 : s_rnge("uum", i__3, "cmpi\
nt_", (ftnlen)1682)] = zint;
	    if (*fbeam > 0.) {
		uum[(i__4 = iq + lu * uum_dim1 - uum_offset) < 1 * uum_dim1 * 
			uum_dim2 && 0 <= i__4 ? i__4 : s_rnge("uum", i__4, 
			"cmpint_", (ftnlen)1684)] = zint + zz[iq + lyu * 
			zz_dim1 - zz_offset] * exp(-utaupr[(i__3 = lu - 1) < 
			1 * utaupr_dim1 && 0 <= i__3 ? i__3 : s_rnge("utaupr",
			 i__3, "cmpint_", (ftnlen)1684)] / *umu0);
	    }
	    if (*plank && *mazim == 0) {
		uum[(i__3 = iq + lu * uum_dim1 - uum_offset) < 1 * uum_dim1 * 
			uum_dim2 && 0 <= i__3 ? i__3 : s_rnge("uum", i__3, 
			"cmpint_", (ftnlen)1687)] = uum[(i__4 = iq + lu * 
			uum_dim1 - uum_offset) < 1 * uum_dim1 * uum_dim2 && 0 
			<= i__4 ? i__4 : s_rnge("uum", i__4, "cmpint_", (
			ftnlen)1687)] + zplk0[iq + lyu * zplk0_dim1 - 
			zplk0_offset] + zplk1[iq + lyu * zplk1_dim1 - 
			zplk1_offset] * utaupr[(i__5 = lu - 1) < 1 * 
			utaupr_dim1 && 0 <= i__5 ? i__5 : s_rnge("utaupr", 
			i__5, "cmpint_", (ftnlen)1687)];
	    }
/* L30: */
	}
L40:
	;
    }
    return 0;
} /* cmpint_ */

/* Subroutine */ int fluxes_(doublereal *cmu, doublereal *cwt, doublereal *
	fbeam, doublereal *gc, doublereal *kk, integer *layru, doublereal *ll,
	 logical *lyrcut, integer *maxulv, integer *mxcmu, integer *mxulv, 
	integer *ncut, integer *nn, integer *nstr, integer *ntau, doublereal *
	pi, logical *prnt, doublereal *ssalb, doublereal *taucpr, doublereal *
	umu0, doublereal *utau, doublereal *utaupr, doublereal *xr0, 
	doublereal *xr1, doublereal *zz, doublereal *zplk0, doublereal *zplk1,
	 doublereal *dfdt, doublereal *flup, doublereal *fldn, doublereal *
	fldir, doublereal *rfldir, doublereal *rfldn, doublereal *uavg, 
	doublereal *u0c)
{
    /* System generated locals */
    integer layru_dim1, cmu_dim1, cwt_dim1, dfdt_dim1, fldir_dim1, fldn_dim1, 
	    flup_dim1, gc_dim1, gc_dim2, gc_offset, kk_dim1, kk_offset, 
	    ll_dim1, ll_offset, rfldir_dim1, rfldn_dim1, u0c_dim1, u0c_dim2, 
	    u0c_offset, uavg_dim1, utau_dim1, utaupr_dim1, zplk0_dim1, 
	    zplk0_offset, zplk1_dim1, zplk1_offset, zz_dim1, zz_offset, i__1, 
	    i__2, i__3, i__4, i__5, i__6, i__7;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_rnge(char *, integer, char *, integer);
    double exp(doublereal), acos(doublereal);

    /* Local variables */
    static doublereal fact, fnet, zint;
    static integer iq, jq, lu;
    static doublereal dirint, fdntot, plsorc;
    extern /* Subroutine */ int zeroit_(doublereal *, integer *);
    static integer lyu;
    static doublereal ang1, ang2;

    /* Fortran I/O blocks */
    static cilist io___139 = { 0, 6, 0, "(//,21X,A,/,2A,/,2A,/)", 0 };
    static cilist io___150 = { 0, 6, 0, "(F10.4,I7,1P,7E12.3,E14.3)", 0 };
    static cilist io___151 = { 0, 6, 0, "(//,2A)", 0 };
    static cilist io___152 = { 0, 6, 0, "(/,A,F10.4,//,2A)", 0 };
    static cilist io___155 = { 0, 6, 0, "(2(0P,F16.4,F13.5,1P,E14.3))", 0 };


/*       Calculates the radiative fluxes, mean intensity, and flux */
/*       derivative with respect to optical depth from the m=0 intensity */
/*       components (the azimuthally-averaged intensity) */


/*    I N P U T     V A R I A B L E S: */

/*       CMU      :  Abscissae for Gauss quadrature over angle cosine */

/*       CWT      :  Weights for Gauss quadrature over angle cosine */

/*       GC       :  Eigenvectors at polar quadrature angles, SC(1) */

/*       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7) */

/*       LAYRU    :  Layer number of user level UTAU */

/*       LL       :  Constants of integration in Eq. SC(1), obtained */
/*                   by solving scaled version of Eq. SC(5); */
/*                   exponential term of Eq. SC(12) not included */

/*       LYRCUT   :  Logical flag for truncation of comput. layer */

/*       NN       :  Order of double-Gauss quadrature (NSTR/2) */

/*       NCUT     :  Number of computational layer where absorption */
/*                   optical depth exceeds ABSCUT */

/*       TAUCPR   :  Cumulative optical depth (delta-M-scaled) */

/*       UTAUPR   :  Optical depths of user output levels in delta-M */
/*                   coordinates;  equal to UTAU if no delta-M */

/*       XR0      :  Expansion of thermal source function in Eq. SS(14) */

/*       XR1      :  Expansion of thermal source function Eqs. SS(16) */

/*       ZZ       :  Beam source vectors in Eq. SS(19) */

/*       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16) */

/*       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16) */

/*       (remainder are DISORT input variables) */


/*    O U T P U T     V A R I A B L E S: */

/*       U0C      :  Azimuthally averaged intensities */
/*                   ( at polar quadrature angles ) */

/*       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables) */


/*    I N T E R N A L       V A R I A B L E S: */

/*       DIRINT   :  Direct intensity attenuated */
/*       FDNTOT   :  Total downward flux (direct + diffuse) */
/*       FLDIR    :  Direct-beam flux (delta-M scaled) */
/*       FLDN     :  Diffuse down-flux (delta-M scaled) */
/*       FNET     :  Net flux (total-down - diffuse-up) */
/*       FACT     :  EXP( - UTAUPR / UMU0 ) */
/*       PLSORC   :  Planck source function (thermal) */
/*       ZINT     :  Intensity of m = 0 case, in Eq. SC(1) */

/*   Called by- DISORT */
/*   Calls- ZEROIT */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    uavg_dim1 = *maxulv;
    rfldn_dim1 = *maxulv;
    rfldir_dim1 = *maxulv;
    flup_dim1 = *maxulv;
    dfdt_dim1 = *maxulv;
    utau_dim1 = *maxulv;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1 * 1;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1 * 1;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1 * 1;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    cwt_dim1 = *mxcmu;
    cmu_dim1 = *mxcmu;
    u0c_dim1 = *mxcmu;
    u0c_dim2 = *mxulv;
    u0c_offset = 1 + u0c_dim1 * 1;
    fldir_dim1 = *mxulv;
    fldn_dim1 = *mxulv;
    utaupr_dim1 = *mxulv;
    layru_dim1 = *mxulv;

    /* Function Body */
    if (prnt[1]) {
	s_wsfe(&io___139);
	do_fio(&c__1, "<----------------------- FLUXES ---------------------\
-->", (ftnlen)56);
	do_fio(&c__1, "   Optical  Compu    Downward    Downward    Downward\
     ", (ftnlen)58);
	do_fio(&c__1, " Upward                    Mean      Planck   d(Net F\
lux)", (ftnlen)57);
	do_fio(&c__1, "     Depth  Layer      Direct     Diffuse       Total\
     ", (ftnlen)58);
	do_fio(&c__1, "Diffuse         Net   Intensity      Source   / d(Op \
Dep)", (ftnlen)57);
	e_wsfe();
    }
/*                                        ** Zero DISORT output arrays */
    i__1 = *mxulv * *mxcmu;
    zeroit_(u0c, &i__1);
    zeroit_(fldir, mxulv);
    zeroit_(fldn, mxulv);
/*                                        ** Loop over user levels */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	lyu = layru[(i__2 = lu - 1) < 1 * layru_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("layru", i__2, "fluxes_", (ftnlen)1821)];
	if (*lyrcut && lyu > *ncut) {
/*                                                ** No radiation reaches */
/*                                                ** this level */
	    fdntot = 0.;
	    fnet = 0.;
	    plsorc = 0.;
	    goto L70;
	}
	if (*fbeam > 0.) {
	    fact = exp(-utaupr[(i__2 = lu - 1) < 1 * utaupr_dim1 && 0 <= i__2 
		    ? i__2 : s_rnge("utaupr", i__2, "fluxes_", (ftnlen)1836)] 
		    / *umu0);
	    dirint = *fbeam * fact;
	    fldir[(i__2 = lu - 1) < 1 * fldir_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("fldir", i__2, "fluxes_", (ftnlen)1838)] = *umu0 * 
		    (*fbeam * fact);
	    rfldir[(i__3 = lu - 1) < 1 * rfldir_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("rfldir", i__3, "fluxes_", (ftnlen)1839)] = *umu0 *
		     *fbeam * exp(-utau[(i__2 = lu - 1) < 1 * utau_dim1 && 0 
		    <= i__2 ? i__2 : s_rnge("utau", i__2, "fluxes_", (ftnlen)
		    1839)] / *umu0);
	} else {
	    dirint = 0.;
	    fldir[(i__2 = lu - 1) < 1 * fldir_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("fldir", i__2, "fluxes_", (ftnlen)1844)] = 0.;
	    rfldir[(i__2 = lu - 1) < 1 * rfldir_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("rfldir", i__2, "fluxes_", (ftnlen)1845)] = 0.;
	}
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    zint = 0.;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1 - gc_offset] * 
			ll[jq + lyu * ll_dim1 - ll_offset] * exp(-kk[jq + lyu 
			* kk_dim1 - kk_offset] * (utaupr[(i__4 = lu - 1) < 1 *
			 utaupr_dim1 && 0 <= i__4 ? i__4 : s_rnge("utaupr", 
			i__4, "fluxes_", (ftnlen)1855)] - taucpr[lyu]));
/* L10: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1 - gc_offset] * 
			ll[jq + lyu * ll_dim1 - ll_offset] * exp(-kk[jq + lyu 
			* kk_dim1 - kk_offset] * (utaupr[(i__4 = lu - 1) < 1 *
			 utaupr_dim1 && 0 <= i__4 ? i__4 : s_rnge("utaupr", 
			i__4, "fluxes_", (ftnlen)1861)] - taucpr[lyu - 1]));
/* L20: */
	    }
	    u0c[(i__3 = iq + lu * u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * 
		    u0c_dim2 && 0 <= i__3 ? i__3 : s_rnge("u0c", i__3, "flux\
es_", (ftnlen)1866)] = zint;
	    if (*fbeam > 0.) {
		u0c[(i__3 = iq + lu * u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * 
			u0c_dim2 && 0 <= i__3 ? i__3 : s_rnge("u0c", i__3, 
			"fluxes_", (ftnlen)1868)] = zint + zz[iq + lyu * 
			zz_dim1 - zz_offset] * fact;
	    }
	    u0c[(i__3 = iq + lu * u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * 
		    u0c_dim2 && 0 <= i__3 ? i__3 : s_rnge("u0c", i__3, "flux\
es_", (ftnlen)1870)] = u0c[(i__4 = iq + lu * u0c_dim1 - u0c_offset) < 1 * 
		    u0c_dim1 * u0c_dim2 && 0 <= i__4 ? i__4 : s_rnge("u0c", 
		    i__4, "fluxes_", (ftnlen)1870)] + zplk0[iq + lyu * 
		    zplk0_dim1 - zplk0_offset] + zplk1[iq + lyu * zplk1_dim1 
		    - zplk1_offset] * utaupr[(i__5 = lu - 1) < 1 * 
		    utaupr_dim1 && 0 <= i__5 ? i__5 : s_rnge("utaupr", i__5, 
		    "fluxes_", (ftnlen)1870)];
	    uavg[(i__3 = lu - 1) < 1 * uavg_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "uavg", i__3, "fluxes_", (ftnlen)1872)] = uavg[(i__4 = lu 
		    - 1) < 1 * uavg_dim1 && 0 <= i__4 ? i__4 : s_rnge("uavg", 
		    i__4, "fluxes_", (ftnlen)1872)] + cwt[(i__5 = *nn + 1 - 
		    iq - 1) < 1 * cwt_dim1 && 0 <= i__5 ? i__5 : s_rnge("cwt",
		     i__5, "fluxes_", (ftnlen)1872)] * u0c[(i__6 = iq + lu * 
		    u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * u0c_dim2 && 0 <= 
		    i__6 ? i__6 : s_rnge("u0c", i__6, "fluxes_", (ftnlen)1872)
		    ];
	    fldn[(i__3 = lu - 1) < 1 * fldn_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "fldn", i__3, "fluxes_", (ftnlen)1873)] = fldn[(i__4 = lu 
		    - 1) < 1 * fldn_dim1 && 0 <= i__4 ? i__4 : s_rnge("fldn", 
		    i__4, "fluxes_", (ftnlen)1873)] + cwt[(i__5 = *nn + 1 - 
		    iq - 1) < 1 * cwt_dim1 && 0 <= i__5 ? i__5 : s_rnge("cwt",
		     i__5, "fluxes_", (ftnlen)1873)] * cmu[(i__6 = *nn + 1 - 
		    iq - 1) < 1 * cmu_dim1 && 0 <= i__6 ? i__6 : s_rnge("cmu",
		     i__6, "fluxes_", (ftnlen)1873)] * u0c[(i__7 = iq + lu * 
		    u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * u0c_dim2 && 0 <= 
		    i__7 ? i__7 : s_rnge("u0c", i__7, "fluxes_", (ftnlen)1873)
		    ];
/* L30: */
	}
	i__2 = *nstr;
	for (iq = *nn + 1; iq <= i__2; ++iq) {
	    zint = 0.;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1 - gc_offset] * 
			ll[jq + lyu * ll_dim1 - ll_offset] * exp(-kk[jq + lyu 
			* kk_dim1 - kk_offset] * (utaupr[(i__4 = lu - 1) < 1 *
			 utaupr_dim1 && 0 <= i__4 ? i__4 : s_rnge("utaupr", 
			i__4, "fluxes_", (ftnlen)1883)] - taucpr[lyu]));
/* L40: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1 - gc_offset] * 
			ll[jq + lyu * ll_dim1 - ll_offset] * exp(-kk[jq + lyu 
			* kk_dim1 - kk_offset] * (utaupr[(i__4 = lu - 1) < 1 *
			 utaupr_dim1 && 0 <= i__4 ? i__4 : s_rnge("utaupr", 
			i__4, "fluxes_", (ftnlen)1889)] - taucpr[lyu - 1]));
/* L50: */
	    }
	    u0c[(i__3 = iq + lu * u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * 
		    u0c_dim2 && 0 <= i__3 ? i__3 : s_rnge("u0c", i__3, "flux\
es_", (ftnlen)1894)] = zint;
	    if (*fbeam > 0.) {
		u0c[(i__3 = iq + lu * u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * 
			u0c_dim2 && 0 <= i__3 ? i__3 : s_rnge("u0c", i__3, 
			"fluxes_", (ftnlen)1896)] = zint + zz[iq + lyu * 
			zz_dim1 - zz_offset] * fact;
	    }
	    u0c[(i__3 = iq + lu * u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * 
		    u0c_dim2 && 0 <= i__3 ? i__3 : s_rnge("u0c", i__3, "flux\
es_", (ftnlen)1898)] = u0c[(i__4 = iq + lu * u0c_dim1 - u0c_offset) < 1 * 
		    u0c_dim1 * u0c_dim2 && 0 <= i__4 ? i__4 : s_rnge("u0c", 
		    i__4, "fluxes_", (ftnlen)1898)] + zplk0[iq + lyu * 
		    zplk0_dim1 - zplk0_offset] + zplk1[iq + lyu * zplk1_dim1 
		    - zplk1_offset] * utaupr[(i__5 = lu - 1) < 1 * 
		    utaupr_dim1 && 0 <= i__5 ? i__5 : s_rnge("utaupr", i__5, 
		    "fluxes_", (ftnlen)1898)];
	    uavg[(i__3 = lu - 1) < 1 * uavg_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "uavg", i__3, "fluxes_", (ftnlen)1900)] = uavg[(i__4 = lu 
		    - 1) < 1 * uavg_dim1 && 0 <= i__4 ? i__4 : s_rnge("uavg", 
		    i__4, "fluxes_", (ftnlen)1900)] + cwt[(i__5 = iq - *nn - 
		    1) < 1 * cwt_dim1 && 0 <= i__5 ? i__5 : s_rnge("cwt", 
		    i__5, "fluxes_", (ftnlen)1900)] * u0c[(i__6 = iq + lu * 
		    u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * u0c_dim2 && 0 <= 
		    i__6 ? i__6 : s_rnge("u0c", i__6, "fluxes_", (ftnlen)1900)
		    ];
	    flup[(i__3 = lu - 1) < 1 * flup_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "flup", i__3, "fluxes_", (ftnlen)1901)] = flup[(i__4 = lu 
		    - 1) < 1 * flup_dim1 && 0 <= i__4 ? i__4 : s_rnge("flup", 
		    i__4, "fluxes_", (ftnlen)1901)] + cwt[(i__5 = iq - *nn - 
		    1) < 1 * cwt_dim1 && 0 <= i__5 ? i__5 : s_rnge("cwt", 
		    i__5, "fluxes_", (ftnlen)1901)] * cmu[(i__6 = iq - *nn - 
		    1) < 1 * cmu_dim1 && 0 <= i__6 ? i__6 : s_rnge("cmu", 
		    i__6, "fluxes_", (ftnlen)1901)] * u0c[(i__7 = iq + lu * 
		    u0c_dim1 - u0c_offset) < 1 * u0c_dim1 * u0c_dim2 && 0 <= 
		    i__7 ? i__7 : s_rnge("u0c", i__7, "fluxes_", (ftnlen)1901)
		    ];
/* L60: */
	}
	flup[(i__2 = lu - 1) < 1 * flup_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"flup", i__2, "fluxes_", (ftnlen)1906)] = *pi * 2. * flup[(
		i__3 = lu - 1) < 1 * flup_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"flup", i__3, "fluxes_", (ftnlen)1906)];
	fldn[(i__2 = lu - 1) < 1 * fldn_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"fldn", i__2, "fluxes_", (ftnlen)1907)] = *pi * 2. * fldn[(
		i__3 = lu - 1) < 1 * fldn_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"fldn", i__3, "fluxes_", (ftnlen)1907)];
	fdntot = fldn[(i__2 = lu - 1) < 1 * fldn_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("fldn", i__2, "fluxes_", (ftnlen)1908)] + fldir[(i__3 =
		 lu - 1) < 1 * fldir_dim1 && 0 <= i__3 ? i__3 : s_rnge("fldir"
		, i__3, "fluxes_", (ftnlen)1908)];
	fnet = fdntot - flup[(i__2 = lu - 1) < 1 * flup_dim1 && 0 <= i__2 ? 
		i__2 : s_rnge("flup", i__2, "fluxes_", (ftnlen)1909)];
	rfldn[(i__2 = lu - 1) < 1 * rfldn_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"rfldn", i__2, "fluxes_", (ftnlen)1910)] = fdntot - rfldir[(
		i__3 = lu - 1) < 1 * rfldir_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"rfldir", i__3, "fluxes_", (ftnlen)1910)];
	uavg[(i__2 = lu - 1) < 1 * uavg_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"uavg", i__2, "fluxes_", (ftnlen)1911)] = (*pi * 2. * uavg[(
		i__3 = lu - 1) < 1 * uavg_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"uavg", i__3, "fluxes_", (ftnlen)1911)] + dirint) / (*pi * 4.)
		;
	plsorc = xr0[lyu - 1] + xr1[lyu - 1] * utaupr[(i__2 = lu - 1) < 1 * 
		utaupr_dim1 && 0 <= i__2 ? i__2 : s_rnge("utaupr", i__2, 
		"fluxes_", (ftnlen)1912)];
	dfdt[(i__2 = lu - 1) < 1 * dfdt_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"dfdt", i__2, "fluxes_", (ftnlen)1913)] = (1. - ssalb[lyu - 1]
		) * 4. * *pi * (uavg[(i__3 = lu - 1) < 1 * uavg_dim1 && 0 <= 
		i__3 ? i__3 : s_rnge("uavg", i__3, "fluxes_", (ftnlen)1913)] 
		- plsorc);
L70:
	if (prnt[1]) {
	    s_wsfe(&io___150);
	    do_fio(&c__1, (char *)&utau[(i__2 = lu - 1) < 1 * utau_dim1 && 0 
		    <= i__2 ? i__2 : s_rnge("utau", i__2, "fluxes_", (ftnlen)
		    1917)], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&lyu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rfldir[(i__3 = lu - 1) < 1 * rfldir_dim1 &&
		     0 <= i__3 ? i__3 : s_rnge("rfldir", i__3, "fluxes_", (
		    ftnlen)1917)], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&rfldn[(i__4 = lu - 1) < 1 * rfldn_dim1 && 
		    0 <= i__4 ? i__4 : s_rnge("rfldn", i__4, "fluxes_", (
		    ftnlen)1917)], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&fdntot, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&flup[(i__5 = lu - 1) < 1 * flup_dim1 && 0 
		    <= i__5 ? i__5 : s_rnge("flup", i__5, "fluxes_", (ftnlen)
		    1917)], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&fnet, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&uavg[(i__6 = lu - 1) < 1 * uavg_dim1 && 0 
		    <= i__6 ? i__6 : s_rnge("uavg", i__6, "fluxes_", (ftnlen)
		    1917)], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&plsorc, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&dfdt[(i__7 = lu - 1) < 1 * dfdt_dim1 && 0 
		    <= i__7 ? i__7 : s_rnge("dfdt", i__7, "fluxes_", (ftnlen)
		    1917)], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
/* L80: */
    }
    if (prnt[2]) {
	s_wsfe(&io___151);
	do_fio(&c__1, " ******** AZIMUTHALLY AVERAGED ", (ftnlen)31);
	do_fio(&c__1, "INTENSITIES ( at polar quadrature angles ) *******", (
		ftnlen)50);
	e_wsfe();
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    s_wsfe(&io___152);
	    do_fio(&c__1, " Optical depth =", (ftnlen)16);
	    do_fio(&c__1, (char *)&utau[(i__2 = lu - 1) < 1 * utau_dim1 && 0 
		    <= i__2 ? i__2 : s_rnge("utau", i__2, "fluxes_", (ftnlen)
		    1931)], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, "     Angle (deg)   cos(Angle)     Intensity", (
		    ftnlen)43);
	    do_fio(&c__1, "     Angle (deg)   cos(Angle)     Intensity", (
		    ftnlen)43);
	    e_wsfe();
	    i__2 = *nn;
	    for (iq = 1; iq <= i__2; ++iq) {
		ang1 = 180. / *pi * acos(cmu[(i__3 = (*nn << 1) - iq) < 1 * 
			cmu_dim1 && 0 <= i__3 ? i__3 : s_rnge("cmu", i__3, 
			"fluxes_", (ftnlen)1937)]);
		ang2 = 180. / *pi * acos(cmu[(i__3 = iq - 1) < 1 * cmu_dim1 &&
			 0 <= i__3 ? i__3 : s_rnge("cmu", i__3, "fluxes_", (
			ftnlen)1938)]);
		s_wsfe(&io___155);
		do_fio(&c__1, (char *)&ang1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&cmu[(i__3 = (*nn << 1) - iq) < 1 * 
			cmu_dim1 && 0 <= i__3 ? i__3 : s_rnge("cmu", i__3, 
			"fluxes_", (ftnlen)1939)], (ftnlen)sizeof(doublereal))
			;
		do_fio(&c__1, (char *)&u0c[(i__4 = iq + lu * u0c_dim1 - 
			u0c_offset) < 1 * u0c_dim1 * u0c_dim2 && 0 <= i__4 ? 
			i__4 : s_rnge("u0c", i__4, "fluxes_", (ftnlen)1939)], 
			(ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&ang2, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&cmu[(i__5 = iq - 1) < 1 * cmu_dim1 && 
			0 <= i__5 ? i__5 : s_rnge("cmu", i__5, "fluxes_", (
			ftnlen)1939)], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&u0c[(i__6 = iq + *nn + lu * u0c_dim1 - 
			u0c_offset) < 1 * u0c_dim1 * u0c_dim2 && 0 <= i__6 ? 
			i__6 : s_rnge("u0c", i__6, "fluxes_", (ftnlen)1939)], 
			(ftnlen)sizeof(doublereal));
		e_wsfe();
/* L90: */
	    }
/* L100: */
	}
    }
    return 0;
} /* fluxes_ */

/* Subroutine */ int setdis_(doublereal *cmu, doublereal *cwt, logical *
	deltam, doublereal *dtauc, doublereal *dtaucp, doublereal *expbea, 
	doublereal *fbeam, doublereal *flyr, doublereal *gl, doublereal *hl, 
	doublereal *hlpr, integer *ibcnd, logical *lamber, integer *layru, 
	logical *lyrcut, integer *maxumu, integer *maxcmu, integer *mxcmu, 
	integer *ncut, integer *nlyr, integer *ntau, integer *nn, integer *
	nstr, logical *plank, integer *numu, logical *onlyfl, doublereal *
	oprim, doublereal *pmom, doublereal *ssalb, doublereal *tauc, 
	doublereal *taucpr, doublereal *utau, doublereal *utaupr, doublereal *
	umu, doublereal *umu0, logical *usrtau, logical *usrang)
{
    /* Initialized data */

    static doublereal abscut = 10.;

    /* System generated locals */
    integer cmu_dim1, cwt_dim1, gl_dim1, gl_offset, hl_dim1, hlpr_dim1, 
	    pmom_dim1, pmom_offset, umu_dim1, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal f;
    static integer k, lc, iq, iu, lu;
    static doublereal abstau;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen), qgausn_(
	    integer *, doublereal *, doublereal *);

/*          Perform miscellaneous setting-up operations */

/*    INPUT :  all are DISORT input variables (see DOC file) */


/*    O U T P U T     V A R I A B L E S: */

/*       NTAU,UTAU   if USRTAU = FALSE (defined in DISORT.doc) */
/*       NUMU,UMU    if USRANG = FALSE (defined in DISORT.doc) */

/*       CMU,CWT     computational polar angles and */
/*                   corresponding quadrature weights */

/*       EXPBEA      transmission of direct beam */

/*       FLYR        truncated fraction in delta-M method */

/*       GL          phase function Legendre coefficients multiplied */
/*                   by (2L+1) and single-scatter albedo */

/*       HLPR        Legendre moments of surface bidirectional */
/*                   reflectivity, times 2K+1 */

/*       LAYRU       Computational layer in which UTAU falls */

/*       LYRCUT      flag as to whether radiation will be zeroed */
/*                   below layer NCUT */

/*       NCUT        computational layer where absorption */
/*                   optical depth first exceeds  ABSCUT */

/*       NN          NSTR / 2 */

/*       OPRIM       delta-M-scaled single-scatter albedo */

/*       TAUCPR      delta-M-scaled optical depth */

/*       UTAUPR      delta-M-scaled version of  UTAU */

/*   Called by- DISORT */
/*   Calls- QGAUSN, ERRMSG */
/* --------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    umu_dim1 = *maxumu;
    pmom_dim1 = *maxcmu - 0 + 1;
    pmom_offset = 0 + pmom_dim1 * 1;
    hl_dim1 = *maxcmu - 0 + 1;
    hlpr_dim1 = *mxcmu - 0 + 1;
    gl_dim1 = *mxcmu - 0 + 1;
    gl_offset = 0 + gl_dim1 * 1;
    cwt_dim1 = *mxcmu;
    cmu_dim1 = *mxcmu;

    /* Function Body */
    if (! (*usrtau)) {
/*                              ** Set output levels at computational */
/*                              ** layer boundaries */
	*ntau = *nlyr + 1;
	i__1 = *ntau - 1;
	for (lc = 0; lc <= i__1; ++lc) {
	    utau[lc] = tauc[lc];
/* L10: */
	}
    }
/*                        ** Apply delta-M scaling and move description */
/*                        ** of computational layers to local variables */
    expbea[0] = 1.;
    taucpr[0] = 0.;
    abstau = 0.;
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	pmom[lc * pmom_dim1 - pmom_offset] = 1.;
	if (abstau < abscut) {
	    *ncut = lc;
	}
	abstau += (1. - ssalb[lc - 1]) * dtauc[lc - 1];
	if (! (*deltam)) {
	    oprim[lc - 1] = ssalb[lc - 1];
	    dtaucp[lc - 1] = dtauc[lc - 1];
	    taucpr[lc] = tauc[lc];
	    i__2 = *nstr - 1;
	    for (k = 0; k <= i__2; ++k) {
		gl[k + lc * gl_dim1 - gl_offset] = ((k << 1) + 1) * oprim[lc 
			- 1] * pmom[k + lc * pmom_dim1 - pmom_offset];
/* L20: */
	    }
	    f = 0.;
	} else {
/*                                    ** Do delta-M transformation */
	    f = pmom[*nstr + lc * pmom_dim1 - pmom_offset];
	    oprim[lc - 1] = ssalb[lc - 1] * (1. - f) / (1. - f * ssalb[lc - 1]
		    );
	    dtaucp[lc - 1] = (1. - f * ssalb[lc - 1]) * dtauc[lc - 1];
	    taucpr[lc] = taucpr[lc - 1] + dtaucp[lc - 1];
	    i__2 = *nstr - 1;
	    for (k = 0; k <= i__2; ++k) {
		gl[k + lc * gl_dim1 - gl_offset] = ((k << 1) + 1) * oprim[lc 
			- 1] * (pmom[k + lc * pmom_dim1 - pmom_offset] - f) / 
			(1. - f);
/* L30: */
	    }
	}
	flyr[lc - 1] = f;
	expbea[lc] = 0.;
	if (*fbeam > 0.) {
	    expbea[lc] = exp(-taucpr[lc] / *umu0);
	}
/* L40: */
    }
/*                      ** If no thermal emission, cut off medium below */
/*                      ** absorption optical depth = ABSCUT ( note that */
/*                      ** delta-M transformation leaves absorption */
/*                      ** optical depth invariant ).  Not worth the */
/*                      ** trouble for one-layer problems, though. */
    *lyrcut = FALSE_;
    if (abstau >= abscut && ! (*plank) && *ibcnd != 1 && *nlyr > 1) {
	*lyrcut = TRUE_;
    }
    if (! (*lyrcut)) {
	*ncut = *nlyr;
    }
/*                             ** Set arrays defining location of user */
/*                             ** output levels within delta-M-scaled */
/*                             ** computational mesh */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	i__2 = *nlyr;
	for (lc = 1; lc <= i__2; ++lc) {
	    if (utau[lu - 1] >= tauc[lc - 1] && utau[lu - 1] <= tauc[lc]) {
		goto L60;
	    }
/* L50: */
	}
	lc = *nlyr;
L60:
	utaupr[lu - 1] = utau[lu - 1];
	if (*deltam) {
	    utaupr[lu - 1] = taucpr[lc - 1] + (1. - ssalb[lc - 1] * flyr[lc - 
		    1]) * (utau[lu - 1] - tauc[lc - 1]);
	}
	layru[lu - 1] = lc;
/* L70: */
    }
/*                      ** Calculate computational polar angle cosines */
/*                      ** and associated quadrature weights for Gaussian */
/*                      ** quadrature on the interval (0,1) (upward) */
    *nn = *nstr / 2;
    qgausn_(nn, cmu, cwt);
/*                                  ** Downward (neg) angles and weights */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	cmu[(i__2 = iq + *nn - 1) < 1 * cmu_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"cmu", i__2, "setdis_", (ftnlen)2132)] = -cmu[(i__3 = iq - 1) 
		< 1 * cmu_dim1 && 0 <= i__3 ? i__3 : s_rnge("cmu", i__3, 
		"setdis_", (ftnlen)2132)];
	cwt[(i__2 = iq + *nn - 1) < 1 * cwt_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"cwt", i__2, "setdis_", (ftnlen)2133)] = cwt[(i__3 = iq - 1) <
		 1 * cwt_dim1 && 0 <= i__3 ? i__3 : s_rnge("cwt", i__3, "set\
dis_", (ftnlen)2133)];
/* L80: */
    }
    if (*fbeam > 0.) {
/*                               ** Compare beam angle to comput. angles */
	i__1 = *nn;
	for (iq = 1; iq <= i__1; ++iq) {
	    if ((d__1 = *umu0 - cmu[(i__2 = iq - 1) < 1 * cmu_dim1 && 0 <= 
		    i__2 ? i__2 : s_rnge("cmu", i__2, "setdis_", (ftnlen)2141)
		    ], abs(d__1)) / *umu0 < 1e-4) {
		errmsg_("SETDIS--beam angle=computational angle; change NSTR",
			 &c_true, (ftnlen)51);
	    }
/* L90: */
	}
    }
    if (! (*usrang) || *onlyfl && *maxumu >= *nstr) {
/*                                   ** Set output polar angles to */
/*                                   ** computational polar angles */
	*numu = *nstr;
	i__1 = *nn;
	for (iu = 1; iu <= i__1; ++iu) {
	    umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "umu", i__2, "setdis_", (ftnlen)2157)] = -cmu[(i__3 = *nn 
		    + 1 - iu - 1) < 1 * cmu_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "cmu", i__3, "setdis_", (ftnlen)2157)];
/* L100: */
	}
	i__1 = *nstr;
	for (iu = *nn + 1; iu <= i__1; ++iu) {
	    umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "umu", i__2, "setdis_", (ftnlen)2161)] = cmu[(i__3 = iu - 
		    *nn - 1) < 1 * cmu_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "cmu", i__3, "setdis_", (ftnlen)2161)];
/* L110: */
	}
    }
    if (*usrang && *ibcnd == 1) {
/*                               ** Shift positive user angle cosines to */
/*                               ** upper locations and put negatives */
/*                               ** in lower locations */
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    umu[(i__2 = iu + *numu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("umu", i__2, "setdis_", (ftnlen)2173)] = umu[(i__3 
		    = iu - 1) < 1 * umu_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "umu", i__3, "setdis_", (ftnlen)2173)];
/* L120: */
	}
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "umu", i__2, "setdis_", (ftnlen)2177)] = -umu[(i__3 = (*
		    numu << 1) + 1 - iu - 1) < 1 * umu_dim1 && 0 <= i__3 ? 
		    i__3 : s_rnge("umu", i__3, "setdis_", (ftnlen)2177)];
/* L130: */
	}
	*numu <<= 1;
    }
    if (! (*lyrcut) && ! (*lamber)) {
	i__1 = *nstr;
	for (k = 0; k <= i__1; ++k) {
	    hlpr[(i__2 = k) < 1 * hlpr_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "hlpr", i__2, "setdis_", (ftnlen)2188)] = ((k << 1) + 1) *
		     hl[(i__3 = k) < 1 * hl_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "hl", i__3, "setdis_", (ftnlen)2188)];
/* L140: */
	}
    }
    return 0;
} /* setdis_ */

/* Subroutine */ int setmtx_(doublereal *bdr, doublereal *cband, doublereal *
	cmu, doublereal *cwt, doublereal *delm0, doublereal *dtaucp, 
	doublereal *gc, doublereal *kk, logical *lamber, logical *lyrcut, 
	integer *mi, integer *mi9m2, integer *mxcmu, integer *ncol, integer *
	ncut, integer *nnlyri, integer *nn, integer *nstr, doublereal *taucpr,
	 doublereal *wk)
{
    /* System generated locals */
    integer bdr_dim1, bdr_dim2, bdr_offset, cband_dim1, cband_dim2, 
	    cband_offset, cmu_dim1, cwt_dim1, gc_dim1, gc_dim2, gc_offset, 
	    kk_dim1, kk_offset, wk_dim1, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    double exp(doublereal);

    /* Local variables */
    static integer jcol;
    static doublereal expa;
    static integer irow, k, nncol, lc, iq, jq, nshift;
    extern /* Subroutine */ int zeroit_(doublereal *, integer *);
    static integer lda, ncd;
    static doublereal sum;

/*        Calculate coefficient matrix for the set of equations */
/*        obtained from the boundary conditions and the continuity- */
/*        of-intensity-at-layer-interface equations;  store in the */
/*        special banded-matrix format required by LINPACK routines */


/*    I N P U T      V A R I A B L E S: */

/*       BDR      :  surface bidirectional reflectivity */

/*       CMU,CWT     abscissae, weights for Gauss quadrature */
/*                   over angle cosine */

/*       DELM0    :  Kronecker delta, delta-sub-m0 */

/*       GC       :  Eigenvectors at polar quadrature angles, SC(1) */

/*       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7) */

/*       LYRCUT   :  Logical flag for truncation of computational layers */

/*       NN       :  Number of streams in a hemisphere (NSTR/2) */

/*       NCUT     :  Total number of computational layers considered */

/*       TAUCPR   :  Cumulative optical depth (delta-M-scaled) */

/*       (remainder are DISORT input variables) */


/*   O U T P U T     V A R I A B L E S: */

/*       CBAND    :  Left-hand side matrix of linear system Eq. SC(5), */
/*                   scaled by Eq. SC(12); in banded form required */
/*                   by LINPACK solution routines */

/*       NCOL     :  Number of columns in CBAND */


/*   I N T E R N A L    V A R I A B L E S: */

/*       IROW     :  Points to row in CBAND */
/*       JCOL     :  Points to position in layer block */
/*       LDA      :  Row dimension of CBAND */
/*       NCD      :  Number of diagonals below or above main diagonal */
/*       NSHIFT   :  For positioning number of rows in band storage */
/*       WK       :  Temporary storage for EXP evaluations */


/*   BAND STORAGE */

/*      LINPACK requires band matrices to be input in a special */
/*      form where the elements of each diagonal are moved up or */
/*      down (in their column) so that each diagonal becomes a row. */
/*      (The column locations of diagonal elements are unchanged.) */

/*      Example:  if the original matrix is */

/*          11 12 13  0  0  0 */
/*          21 22 23 24  0  0 */
/*           0 32 33 34 35  0 */
/*           0  0 43 44 45 46 */
/*           0  0  0 54 55 56 */
/*           0  0  0  0 65 66 */

/*      then its LINPACK input form would be: */

/*           *  *  *  +  +  +  , * = not used */
/*           *  * 13 24 35 46  , + = used for pivoting */
/*           * 12 23 34 45 56 */
/*          11 22 33 44 55 66 */
/*          21 32 43 54 65  * */

/*      If A is a band matrix, the following program segment */
/*      will convert it to the form (ABD) required by LINPACK */
/*      band-matrix routines: */

/*               N  = (column dimension of A, ABD) */
/*               ML = (band width below the diagonal) */
/*               MU = (band width above the diagonal) */
/*               M = ML + MU + 1 */
/*               DO J = 1, N */
/*                  I1 = MAX(1, J-MU) */
/*                  I2 = MIN(N, J+ML) */
/*                  DO I = I1, I2 */
/*                     K = I - J + M */
/*                     ABD(K,J) = A(I,J) */
/*                  END DO */
/*               END DO */

/*      This uses rows  ML+1  through  2*ML+MU+1  of ABD. */
/*      The total number of rows needed in ABD is  2*ML+MU+1 . */
/*      In the example above, N = 6, ML = 1, MU = 2, and the */
/*      row dimension of ABD must be >= 5. */


/*   Called by- DISORT, ALBTRN */
/*   Calls- ZEROIT */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    bdr_dim1 = *mi;
    bdr_dim2 = *mi - 0 + 1;
    bdr_offset = 1 + bdr_dim1 * 0;
    wk_dim1 = *mxcmu;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    cwt_dim1 = *mxcmu;
    cmu_dim1 = *mxcmu;
    cband_dim1 = *mi9m2;
    cband_dim2 = *nnlyri;
    cband_offset = 1 + cband_dim1 * 1;

    /* Function Body */
    i__1 = *mi9m2 * *nnlyri;
    zeroit_(cband, &i__1);
    ncd = *nn * 3 - 1;
    lda = ncd * 3 + 1;
    nshift = lda - (*nstr << 1) + 1;
    *ncol = 0;
/*                         ** Use continuity conditions of Eq. STWJ(17) */
/*                         ** to form coefficient matrix in STWJ(20); */
/*                         ** employ scaling transformation STWJ(22) */
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    wk[(i__3 = iq - 1) < 1 * wk_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "wk", i__3, "setmtx_", (ftnlen)2340)] = exp(kk[iq + lc * 
		    kk_dim1 - kk_offset] * dtaucp[lc - 1]);
/* L10: */
	}
	jcol = 0;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ++(*ncol);
	    irow = nshift - jcol;
	    i__3 = *nstr;
	    for (jq = 1; jq <= i__3; ++jq) {
		cband[(i__4 = irow + *nstr + *ncol * cband_dim1 - 
			cband_offset) < 1 * cband_dim1 * cband_dim2 && 0 <= 
			i__4 ? i__4 : s_rnge("cband", i__4, "setmtx_", (
			ftnlen)2351)] = gc[jq + (iq + lc * gc_dim2) * gc_dim1 
			- gc_offset];
		cband[(i__4 = irow + *ncol * cband_dim1 - cband_offset) < 1 * 
			cband_dim1 * cband_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			"cband", i__4, "setmtx_", (ftnlen)2352)] = -gc[jq + (
			iq + lc * gc_dim2) * gc_dim1 - gc_offset] * wk[(i__5 =
			 iq - 1) < 1 * wk_dim1 && 0 <= i__5 ? i__5 : s_rnge(
			"wk", i__5, "setmtx_", (ftnlen)2352)];
		++irow;
/* L20: */
	    }
	    ++jcol;
/* L30: */
	}
	i__2 = *nstr;
	for (iq = *nn + 1; iq <= i__2; ++iq) {
	    ++(*ncol);
	    irow = nshift - jcol;
	    i__3 = *nstr;
	    for (jq = 1; jq <= i__3; ++jq) {
		cband[(i__4 = irow + *nstr + *ncol * cband_dim1 - 
			cband_offset) < 1 * cband_dim1 * cband_dim2 && 0 <= 
			i__4 ? i__4 : s_rnge("cband", i__4, "setmtx_", (
			ftnlen)2367)] = gc[jq + (iq + lc * gc_dim2) * gc_dim1 
			- gc_offset] * wk[(i__5 = *nstr + 1 - iq - 1) < 1 * 
			wk_dim1 && 0 <= i__5 ? i__5 : s_rnge("wk", i__5, 
			"setmtx_", (ftnlen)2367)];
		cband[(i__4 = irow + *ncol * cband_dim1 - cband_offset) < 1 * 
			cband_dim1 * cband_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			"cband", i__4, "setmtx_", (ftnlen)2369)] = -gc[jq + (
			iq + lc * gc_dim2) * gc_dim1 - gc_offset];
		++irow;
/* L40: */
	    }
	    ++jcol;
/* L50: */
	}
/* L60: */
    }
/*                  ** Use top boundary condition of STWJ(20a) for */
/*                  ** first layer */
    jcol = 0;
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	expa = exp(kk[iq + kk_dim1 - kk_offset] * taucpr[1]);
	irow = nshift - jcol + *nn;
	for (jq = *nn; jq >= 1; --jq) {
	    cband[(i__2 = irow + (jcol + 1) * cband_dim1 - cband_offset) < 1 *
		     cband_dim1 * cband_dim2 && 0 <= i__2 ? i__2 : s_rnge(
		    "cband", i__2, "setmtx_", (ftnlen)2388)] = gc[jq + (iq + 
		    gc_dim2) * gc_dim1 - gc_offset] * expa;
	    ++irow;
/* L70: */
	}
	++jcol;
/* L80: */
    }
    i__1 = *nstr;
    for (iq = *nn + 1; iq <= i__1; ++iq) {
	irow = nshift - jcol + *nn;
	for (jq = *nn; jq >= 1; --jq) {
	    cband[(i__2 = irow + (jcol + 1) * cband_dim1 - cband_offset) < 1 *
		     cband_dim1 * cband_dim2 && 0 <= i__2 ? i__2 : s_rnge(
		    "cband", i__2, "setmtx_", (ftnlen)2402)] = gc[jq + (iq + 
		    gc_dim2) * gc_dim1 - gc_offset];
	    ++irow;
/* L90: */
	}
	++jcol;
/* L100: */
    }
/*                           ** Use bottom boundary condition of */
/*                           ** STWJ(20c) for last layer */
    nncol = *ncol - *nstr;
    jcol = 0;
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	++nncol;
	irow = nshift - jcol + *nstr;
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    if (*lyrcut || *lamber && *delm0 == 0.) {
/*                          ** No azimuthal-dependent intensity if Lam- */
/*                          ** bert surface; no intensity component if */
/*                          ** truncated bottom layer */
		cband[(i__3 = irow + nncol * cband_dim1 - cband_offset) < 1 * 
			cband_dim1 * cband_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"cband", i__3, "setmtx_", (ftnlen)2427)] = gc[jq + (
			iq + *ncut * gc_dim2) * gc_dim1 - gc_offset];
	    } else {
		sum = 0.;
		i__3 = *nn;
		for (k = 1; k <= i__3; ++k) {
		    sum += cwt[(i__4 = k - 1) < 1 * cwt_dim1 && 0 <= i__4 ? 
			    i__4 : s_rnge("cwt", i__4, "setmtx_", (ftnlen)
			    2434)] * cmu[(i__5 = k - 1) < 1 * cmu_dim1 && 0 <=
			     i__5 ? i__5 : s_rnge("cmu", i__5, "setmtx_", (
			    ftnlen)2434)] * bdr[(i__6 = jq - *nn + k * 
			    bdr_dim1 - bdr_offset) < 1 * bdr_dim1 * bdr_dim2 
			    && 0 <= i__6 ? i__6 : s_rnge("bdr", i__6, "setmt\
x_", (ftnlen)2434)] * gc[*nn + 1 - k + (iq + *ncut * gc_dim2) * gc_dim1 - 
			    gc_offset];
/* L110: */
		}
		cband[(i__3 = irow + nncol * cband_dim1 - cband_offset) < 1 * 
			cband_dim1 * cband_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"cband", i__3, "setmtx_", (ftnlen)2438)] = gc[jq + (
			iq + *ncut * gc_dim2) * gc_dim1 - gc_offset] - (*
			delm0 + 1.) * sum;
	    }
	    ++irow;
/* L120: */
	}
	++jcol;
/* L130: */
    }
    i__1 = *nstr;
    for (iq = *nn + 1; iq <= i__1; ++iq) {
	++nncol;
	irow = nshift - jcol + *nstr;
	expa = wk[(i__2 = *nstr + 1 - iq - 1) < 1 * wk_dim1 && 0 <= i__2 ? 
		i__2 : s_rnge("wk", i__2, "setmtx_", (ftnlen)2455)];
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    if (*lyrcut || *lamber && *delm0 == 0.) {
		cband[(i__3 = irow + nncol * cband_dim1 - cband_offset) < 1 * 
			cband_dim1 * cband_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"cband", i__3, "setmtx_", (ftnlen)2461)] = gc[jq + (
			iq + *ncut * gc_dim2) * gc_dim1 - gc_offset] * expa;
	    } else {
		sum = 0.;
		i__3 = *nn;
		for (k = 1; k <= i__3; ++k) {
		    sum += cwt[(i__4 = k - 1) < 1 * cwt_dim1 && 0 <= i__4 ? 
			    i__4 : s_rnge("cwt", i__4, "setmtx_", (ftnlen)
			    2468)] * cmu[(i__5 = k - 1) < 1 * cmu_dim1 && 0 <=
			     i__5 ? i__5 : s_rnge("cmu", i__5, "setmtx_", (
			    ftnlen)2468)] * bdr[(i__6 = jq - *nn + k * 
			    bdr_dim1 - bdr_offset) < 1 * bdr_dim1 * bdr_dim2 
			    && 0 <= i__6 ? i__6 : s_rnge("bdr", i__6, "setmt\
x_", (ftnlen)2468)] * gc[*nn + 1 - k + (iq + *ncut * gc_dim2) * gc_dim1 - 
			    gc_offset];
/* L140: */
		}
		cband[(i__3 = irow + nncol * cband_dim1 - cband_offset) < 1 * 
			cband_dim1 * cband_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			"cband", i__3, "setmtx_", (ftnlen)2472)] = (gc[jq + (
			iq + *ncut * gc_dim2) * gc_dim1 - gc_offset] - (*
			delm0 + 1.) * sum) * expa;
	    }
	    ++irow;
/* L150: */
	}
	++jcol;
/* L160: */
    }
    return 0;
} /* setmtx_ */

/* Subroutine */ int soleig_(doublereal *amb, doublereal *apb, doublereal *
	array, doublereal *cmu, doublereal *cwt, doublereal *gl, integer *mi, 
	integer *mazim, integer *mxcmu, integer *nn, integer *nstr, 
	doublereal *ylmc, doublereal *cc, doublereal *evecc, doublereal *eval,
	 doublereal *kk, doublereal *gc, doublereal *aad, doublereal *eveccd, 
	doublereal *evald, doublereal *wkd)
{
    /* System generated locals */
    integer amb_dim1, amb_dim2, amb_offset, apb_dim1, apb_dim2, apb_offset, 
	    array_dim1, array_offset, cc_dim1, cc_dim2, cc_offset, cmu_dim1, 
	    cwt_dim1, eval_dim1, evecc_dim1, evecc_dim2, evecc_offset, 
	    gc_dim1, gc_dim2, gc_offset, gl_dim1, kk_dim1, ylmc_dim1, 
	    ylmc_dim2, ylmc_offset, aad_dim1, aad_offset, eveccd_dim1, 
	    eveccd_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer), s_wsfe(cilist *), 
	    do_fio(integer *, char *, ftnlen), e_wsfe();
    double sqrt(doublereal);

    /* Local variables */
    static doublereal beta;
    static integer l;
    static doublereal alpha;
    static integer iq, jq, kq;
    static doublereal gpmigm, gpplgm;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen), asymtx_(
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer ier;
    static doublereal sum;

    /* Fortran I/O blocks */
    static cilist io___184 = { 0, 6, 0, "(//,A,I4,A)", 0 };


/*         Solves eigenvalue/vector problem necessary to construct */
/*         homogeneous part of discrete ordinate solution; STWJ(8b) */
/*         ** NOTE ** Eigenvalue problem is degenerate when single */
/*                    scattering albedo = 1;  present way of doing it */
/*                    seems numerically more stable than alternative */
/*                    methods that we tried */


/*   I N P U T     V A R I A B L E S: */

/*       GL     :  Delta-M scaled Legendre coefficients of phase function */
/*                 (including factors 2l+1 and single-scatter albedo) */

/*       CMU    :  Computational polar angle cosines */

/*       CWT    :  Weights for quadrature over polar angle cosine */

/*       MAZIM  :  Order of azimuthal component */

/*       NN     :  Half the total number of streams */

/*       YLMC   :  Normalized associated Legendre polynomial */
/*                 at the quadrature angles CMU */

/*       (remainder are DISORT input variables) */


/*   O U T P U T    V A R I A B L E S: */

/*       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18) */

/*       EVAL   :  NN eigenvalues of Eq. SS(12) on return from ASYMTX */
/*                 but then square roots taken */

/*       EVECC  :  NN eigenvectors  (G+) - (G-)  on return */
/*                 from ASYMTX ( column j corresponds to EVAL(j) ) */
/*                 but then  (G+) + (G-)  is calculated from SS(10), */
/*                 G+  and  G-  are separated, and  G+  is stacked on */
/*                 top of  G-  to form NSTR eigenvectors of SS(7) */

/*       GC     :  Permanent storage for all NSTR eigenvectors, but */
/*                 in an order corresponding to KK */

/*       KK     :  Permanent storage for all NSTR eigenvalues of SS(7), */
/*                 but re-ordered with negative values first ( square */
/*                 roots of EVAL taken and negatives added ) */


/*   I N T E R N A L   V A R I A B L E S: */

/*       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced */
/*                    eigenvalue problem */
/*       ARRAY   :  Complete coefficient matrix of reduced eigenvalue */
/*                    problem: (alfa+beta)*(alfa-beta) */
/*       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11)) */
/*       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11)) */
/*       WKD     :  Scratch array required by ASYMTX */

/*   Called by- DISORT, ALBTRN */
/*   Calls- ASYMTX, ERRMSG */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*                             ** Calculate quantities in Eqs. SS(5-6) */
    /* Parameter adjustments */
    eveccd_dim1 = *mi;
    eveccd_offset = 1 + eveccd_dim1 * 1;
    aad_dim1 = *mi;
    aad_offset = 1 + aad_dim1 * 1;
    eval_dim1 = *mi;
    array_dim1 = *mi;
    array_offset = 1 + array_dim1 * 1;
    apb_dim1 = *mi;
    apb_dim2 = *mi;
    apb_offset = 1 + apb_dim1 * 1;
    amb_dim1 = *mi;
    amb_dim2 = *mi;
    amb_offset = 1 + amb_dim1 * 1;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * 1;
    kk_dim1 = *mxcmu;
    evecc_dim1 = *mxcmu;
    evecc_dim2 = *mxcmu;
    evecc_offset = 1 + evecc_dim1 * 1;
    cc_dim1 = *mxcmu;
    cc_dim2 = *mxcmu;
    cc_offset = 1 + cc_dim1 * 1;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_dim2 = *mxcmu;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    gl_dim1 = *mxcmu - 0 + 1;
    cwt_dim1 = *mxcmu;
    cmu_dim1 = *mxcmu;

    /* Function Body */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    sum = 0.;
	    i__3 = *nstr - 1;
	    for (l = *mazim; l <= i__3; ++l) {
		sum += gl[(i__4 = l) < 1 * gl_dim1 && 0 <= i__4 ? i__4 : 
			s_rnge("gl", i__4, "soleig_", (ftnlen)2588)] * ylmc[(
			i__5 = l + iq * ylmc_dim1 - ylmc_offset) < 1 * 
			ylmc_dim1 * ylmc_dim2 && 0 <= i__5 ? i__5 : s_rnge(
			"ylmc", i__5, "soleig_", (ftnlen)2588)] * ylmc[(i__6 =
			 l + jq * ylmc_dim1 - ylmc_offset) < 1 * ylmc_dim1 * 
			ylmc_dim2 && 0 <= i__6 ? i__6 : s_rnge("ylmc", i__6, 
			"soleig_", (ftnlen)2588)];
/* L10: */
	    }
	    cc[(i__3 = iq + jq * cc_dim1 - cc_offset) < 1 * cc_dim1 * cc_dim2 
		    && 0 <= i__3 ? i__3 : s_rnge("cc", i__3, "soleig_", (
		    ftnlen)2591)] = sum * .5 * cwt[(i__4 = jq - 1) < 1 * 
		    cwt_dim1 && 0 <= i__4 ? i__4 : s_rnge("cwt", i__4, "sole\
ig_", (ftnlen)2591)];
/* L20: */
	}
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
/*                             ** Fill remainder of array using symmetry */
/*                             ** relations  C(-mui,muj) = C(mui,-muj) */
/*                             ** and        C(-mui,-muj) = C(mui,muj) */
	    cc[(i__3 = iq + *nn + jq * cc_dim1 - cc_offset) < 1 * cc_dim1 * 
		    cc_dim2 && 0 <= i__3 ? i__3 : s_rnge("cc", i__3, "soleig_"
		    , (ftnlen)2600)] = cc[(i__4 = iq + (jq + *nn) * cc_dim1 - 
		    cc_offset) < 1 * cc_dim1 * cc_dim2 && 0 <= i__4 ? i__4 : 
		    s_rnge("cc", i__4, "soleig_", (ftnlen)2600)];
	    cc[(i__3 = iq + *nn + (jq + *nn) * cc_dim1 - cc_offset) < 1 * 
		    cc_dim1 * cc_dim2 && 0 <= i__3 ? i__3 : s_rnge("cc", i__3,
		     "soleig_", (ftnlen)2601)] = cc[(i__4 = iq + jq * cc_dim1 
		    - cc_offset) < 1 * cc_dim1 * cc_dim2 && 0 <= i__4 ? i__4 :
		     s_rnge("cc", i__4, "soleig_", (ftnlen)2601)];
/*                                       ** Get factors of coeff. matrix */
/*                                       ** of reduced eigenvalue problem */
	    alpha = cc[(i__3 = iq + jq * cc_dim1 - cc_offset) < 1 * cc_dim1 * 
		    cc_dim2 && 0 <= i__3 ? i__3 : s_rnge("cc", i__3, "soleig_"
		    , (ftnlen)2606)] / cmu[(i__4 = iq - 1) < 1 * cmu_dim1 && 
		    0 <= i__4 ? i__4 : s_rnge("cmu", i__4, "soleig_", (ftnlen)
		    2606)];
	    beta = cc[(i__3 = iq + (jq + *nn) * cc_dim1 - cc_offset) < 1 * 
		    cc_dim1 * cc_dim2 && 0 <= i__3 ? i__3 : s_rnge("cc", i__3,
		     "soleig_", (ftnlen)2607)] / cmu[(i__4 = iq - 1) < 1 * 
		    cmu_dim1 && 0 <= i__4 ? i__4 : s_rnge("cmu", i__4, "sole\
ig_", (ftnlen)2607)];
	    amb[(i__3 = iq + jq * amb_dim1 - amb_offset) < 1 * amb_dim1 * 
		    amb_dim2 && 0 <= i__3 ? i__3 : s_rnge("amb", i__3, "sole\
ig_", (ftnlen)2608)] = alpha - beta;
	    apb[(i__3 = iq + jq * apb_dim1 - apb_offset) < 1 * apb_dim1 * 
		    apb_dim2 && 0 <= i__3 ? i__3 : s_rnge("apb", i__3, "sole\
ig_", (ftnlen)2609)] = alpha + beta;
/* L30: */
	}
	amb[(i__2 = iq + iq * amb_dim1 - amb_offset) < 1 * amb_dim1 * 
		amb_dim2 && 0 <= i__2 ? i__2 : s_rnge("amb", i__2, "soleig_", 
		(ftnlen)2613)] = amb[(i__3 = iq + iq * amb_dim1 - amb_offset) 
		< 1 * amb_dim1 * amb_dim2 && 0 <= i__3 ? i__3 : s_rnge("amb", 
		i__3, "soleig_", (ftnlen)2613)] - 1. / cmu[(i__4 = iq - 1) < 
		1 * cmu_dim1 && 0 <= i__4 ? i__4 : s_rnge("cmu", i__4, "sole\
ig_", (ftnlen)2613)];
	apb[(i__2 = iq + iq * apb_dim1 - apb_offset) < 1 * apb_dim1 * 
		apb_dim2 && 0 <= i__2 ? i__2 : s_rnge("apb", i__2, "soleig_", 
		(ftnlen)2614)] = apb[(i__3 = iq + iq * apb_dim1 - apb_offset) 
		< 1 * apb_dim1 * apb_dim2 && 0 <= i__3 ? i__3 : s_rnge("apb", 
		i__3, "soleig_", (ftnlen)2614)] - 1. / cmu[(i__4 = iq - 1) < 
		1 * cmu_dim1 && 0 <= i__4 ? i__4 : s_rnge("cmu", i__4, "sole\
ig_", (ftnlen)2614)];
/* L40: */
    }
/*                      ** Finish calculation of coefficient matrix of */
/*                      ** reduced eigenvalue problem:  get matrix */
/*                      ** product (alfa+beta)*(alfa-beta); SS(12) */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    sum = 0.;
	    i__3 = *nn;
	    for (kq = 1; kq <= i__3; ++kq) {
		sum += apb[(i__4 = iq + kq * apb_dim1 - apb_offset) < 1 * 
			apb_dim1 * apb_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			"apb", i__4, "soleig_", (ftnlen)2626)] * amb[(i__5 = 
			kq + jq * amb_dim1 - amb_offset) < 1 * amb_dim1 * 
			amb_dim2 && 0 <= i__5 ? i__5 : s_rnge("amb", i__5, 
			"soleig_", (ftnlen)2626)];
/* L50: */
	    }
	    array[iq + jq * array_dim1 - array_offset] = sum;
/* L60: */
	}
/* L70: */
    }
/*                      ** Find (real) eigenvalues and eigenvectors */
    asymtx_(array, evecc, eval, nn, mi, mxcmu, &ier, wkd, aad, eveccd, evald);
    if (ier > 0) {
	s_wsfe(&io___184);
	do_fio(&c__1, " ASYMTX--eigenvalue no. ", (ftnlen)24);
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	do_fio(&c__1, "  didnt converge.  Lower-numbered eigenvalues wrong.", 
		(ftnlen)52);
	e_wsfe();
	errmsg_("ASYMTX--convergence problems", &c_true, (ftnlen)28);
    }
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	eval[(i__2 = iq - 1) < 1 * eval_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"eval", i__2, "soleig_", (ftnlen)2649)] = sqrt((d__1 = eval[(
		i__3 = iq - 1) < 1 * eval_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"eval", i__3, "soleig_", (ftnlen)2649)], abs(d__1)));
	kk[(i__2 = iq + *nn - 1) < 1 * kk_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"kk", i__2, "soleig_", (ftnlen)2650)] = eval[(i__3 = iq - 1) <
		 1 * eval_dim1 && 0 <= i__3 ? i__3 : s_rnge("eval", i__3, 
		"soleig_", (ftnlen)2650)];
/*                                      ** Add negative eigenvalue */
	kk[(i__2 = *nn + 1 - iq - 1) < 1 * kk_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("kk", i__2, "soleig_", (ftnlen)2652)] = -eval[(i__3 = 
		iq - 1) < 1 * eval_dim1 && 0 <= i__3 ? i__3 : s_rnge("eval", 
		i__3, "soleig_", (ftnlen)2652)];
/* L80: */
    }
/*                          ** Find eigenvectors (G+) + (G-) from SS(10) */
/*                          ** and store temporarily in APB array */
    i__1 = *nn;
    for (jq = 1; jq <= i__1; ++jq) {
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    sum = 0.;
	    i__3 = *nn;
	    for (kq = 1; kq <= i__3; ++kq) {
		sum += amb[(i__4 = iq + kq * amb_dim1 - amb_offset) < 1 * 
			amb_dim1 * amb_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			"amb", i__4, "soleig_", (ftnlen)2663)] * evecc[(i__5 =
			 kq + jq * evecc_dim1 - evecc_offset) < 1 * 
			evecc_dim1 * evecc_dim2 && 0 <= i__5 ? i__5 : s_rnge(
			"evecc", i__5, "soleig_", (ftnlen)2663)];
/* L90: */
	    }
	    apb[(i__3 = iq + jq * apb_dim1 - apb_offset) < 1 * apb_dim1 * 
		    apb_dim2 && 0 <= i__3 ? i__3 : s_rnge("apb", i__3, "sole\
ig_", (ftnlen)2666)] = sum / eval[(i__4 = jq - 1) < 1 * eval_dim1 && 0 <= 
		    i__4 ? i__4 : s_rnge("eval", i__4, "soleig_", (ftnlen)
		    2666)];
/* L100: */
	}
/* L110: */
    }
    i__1 = *nn;
    for (jq = 1; jq <= i__1; ++jq) {
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    gpplgm = apb[(i__3 = iq + jq * apb_dim1 - apb_offset) < 1 * 
		    apb_dim1 * apb_dim2 && 0 <= i__3 ? i__3 : s_rnge("apb", 
		    i__3, "soleig_", (ftnlen)2677)];
	    gpmigm = evecc[(i__3 = iq + jq * evecc_dim1 - evecc_offset) < 1 * 
		    evecc_dim1 * evecc_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecc", i__3, "soleig_", (ftnlen)2678)];
/*                                ** Recover eigenvectors G+,G- from */
/*                                ** their sum and difference; stack them */
/*                                ** to get eigenvectors of full system */
/*                                ** SS(7) (JQ = eigenvector number) */
	    evecc[(i__3 = iq + jq * evecc_dim1 - evecc_offset) < 1 * 
		    evecc_dim1 * evecc_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecc", i__3, "soleig_", (ftnlen)2684)] = (gpplgm + 
		    gpmigm) * .5;
	    evecc[(i__3 = iq + *nn + jq * evecc_dim1 - evecc_offset) < 1 * 
		    evecc_dim1 * evecc_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecc", i__3, "soleig_", (ftnlen)2685)] = (gpplgm - 
		    gpmigm) * .5;
/*                                ** Eigenvectors corresponding to */
/*                                ** negative eigenvalues (corresp. to */
/*                                ** reversing sign of 'k' in SS(10) ) */
	    gpplgm = -gpplgm;
	    evecc[(i__3 = iq + (jq + *nn) * evecc_dim1 - evecc_offset) < 1 * 
		    evecc_dim1 * evecc_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecc", i__3, "soleig_", (ftnlen)2691)] = (gpplgm + 
		    gpmigm) * .5;
	    evecc[(i__3 = iq + *nn + (jq + *nn) * evecc_dim1 - evecc_offset) <
		     1 * evecc_dim1 * evecc_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "evecc", i__3, "soleig_", (ftnlen)2692)] = (gpplgm - 
		    gpmigm) * .5;
	    gc[(i__3 = iq + *nn + (jq + *nn) * gc_dim1 - gc_offset) < 1 * 
		    gc_dim1 * gc_dim2 && 0 <= i__3 ? i__3 : s_rnge("gc", i__3,
		     "soleig_", (ftnlen)2693)] = evecc[(i__4 = iq + jq * 
		    evecc_dim1 - evecc_offset) < 1 * evecc_dim1 * evecc_dim2 
		    && 0 <= i__4 ? i__4 : s_rnge("evecc", i__4, "soleig_", (
		    ftnlen)2693)];
	    gc[(i__3 = *nn + 1 - iq + (jq + *nn) * gc_dim1 - gc_offset) < 1 * 
		    gc_dim1 * gc_dim2 && 0 <= i__3 ? i__3 : s_rnge("gc", i__3,
		     "soleig_", (ftnlen)2694)] = evecc[(i__4 = iq + *nn + jq *
		     evecc_dim1 - evecc_offset) < 1 * evecc_dim1 * evecc_dim2 
		    && 0 <= i__4 ? i__4 : s_rnge("evecc", i__4, "soleig_", (
		    ftnlen)2694)];
	    gc[(i__3 = iq + *nn + (*nn + 1 - jq) * gc_dim1 - gc_offset) < 1 * 
		    gc_dim1 * gc_dim2 && 0 <= i__3 ? i__3 : s_rnge("gc", i__3,
		     "soleig_", (ftnlen)2695)] = evecc[(i__4 = iq + (jq + *nn)
		     * evecc_dim1 - evecc_offset) < 1 * evecc_dim1 * 
		    evecc_dim2 && 0 <= i__4 ? i__4 : s_rnge("evecc", i__4, 
		    "soleig_", (ftnlen)2695)];
	    gc[(i__3 = *nn + 1 - iq + (*nn + 1 - jq) * gc_dim1 - gc_offset) < 
		    1 * gc_dim1 * gc_dim2 && 0 <= i__3 ? i__3 : s_rnge("gc", 
		    i__3, "soleig_", (ftnlen)2696)] = evecc[(i__4 = iq + *nn 
		    + (jq + *nn) * evecc_dim1 - evecc_offset) < 1 * 
		    evecc_dim1 * evecc_dim2 && 0 <= i__4 ? i__4 : s_rnge(
		    "evecc", i__4, "soleig_", (ftnlen)2696)];
/* L120: */
	}
/* L130: */
    }
    return 0;
} /* soleig_ */

/* Subroutine */ int solve0_(doublereal *b, doublereal *bdr, doublereal *bem, 
	doublereal *bplank, doublereal *cband, doublereal *cmu, doublereal *
	cwt, doublereal *expbea, doublereal *fbeam, doublereal *fisot, 
	integer *ipvt, logical *lamber, doublereal *ll, logical *lyrcut, 
	integer *mazim, integer *mi, integer *mi9m2, integer *mxcmu, integer *
	ncol, integer *ncut, integer *nn, integer *nstr, integer *nnlyri, 
	doublereal *pi, doublereal *tplank, doublereal *taucpr, doublereal *
	umu0, doublereal *z__, doublereal *zz, doublereal *zplk0, doublereal *
	zplk1)
{
    /* System generated locals */
    integer b_dim1, bdr_dim1, bdr_dim2, bdr_offset, bem_dim1, cband_dim1, 
	    cband_offset, cmu_dim1, cwt_dim1, ll_dim1, ll_offset, zplk0_dim1, 
	    zplk0_offset, zplk1_dim1, zplk1_offset, zz_dim1, zz_offset, i__1, 
	    i__2, i__3, i__4, i__5;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer ipnt;
    extern /* Subroutine */ int sgbco_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static doublereal rcond;
    extern /* Subroutine */ int sgbsl_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *);
    static integer lc, iq, jq, it;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen), zeroit_(
	    doublereal *, integer *);
    static integer ncd;
    static doublereal sum;

/*        Construct right-hand side vector B for general boundary */
/*        conditions STWJ(17) and solve system of equations obtained */
/*        from the boundary conditions and the continuity-of- */
/*        intensity-at-layer-interface equations. */
/*        Thermal emission contributes only in azimuthal independence. */


/*    I N P U T      V A R I A B L E S: */

/*       BDR      :  Surface bidirectional reflectivity */

/*       BEM      :  Surface bidirectional emissivity */

/*       BPLANK   :  Bottom boundary thermal emission */

/*       CBAND    :  Left-hand side matrix of linear system Eq. SC(5), */
/*                   scaled by Eq. SC(12); in banded form required */
/*                   by LINPACK solution routines */

/*       CMU,CWT     Abscissae, weights for Gauss quadrature */
/*                   over angle cosine */

/*       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0) */

/*       LYRCUT   :  Logical flag for truncation of computational layers */

/*       MAZIM    :  Order of azimuthal component */

/*       NCOL     :  Number of columns in CBAND */

/*       NN       :  Order of double-Gauss quadrature (NSTR/2) */

/*       NCUT     :  Total number of computational layers considered */

/*       TPLANK   :  Top boundary thermal emission */

/*       TAUCPR   :  Cumulative optical depth (delta-M-scaled) */

/*       ZZ       :  Beam source vectors in Eq. SS(19) */

/*       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16) */

/*       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16) */

/*       (remainder are DISORT input variables) */


/*    O U T P U T     V A R I A B L E S: */

/*       B        :  Right-hand side vector of Eq. SC(5) going into */
/*                   SGBSL; returns as solution vector of Eq. SC(12), */
/*                   constants of integration without exponential term */

/*      LL        :  Permanent storage for B, but re-ordered */


/*   I N T E R N A L    V A R I A B L E S: */

/*       IPVT     :  Integer vector of pivot indices */
/*       IT       :  Pointer for position in  B */
/*       NCD      :  Number of diagonals below or above main diagonal */
/*       RCOND    :  Indicator of singularity for CBAND */
/*       Z        :  Scratch array required by SGBCO */

/*   Called by- DISORT */
/*   Calls- ZEROIT, SGBCO, ERRMSG, SGBSL */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    /* Parameter adjustments */
    bem_dim1 = *mi;
    bdr_dim1 = *mi;
    bdr_dim2 = *mi - 0 + 1;
    bdr_offset = 1 + bdr_dim1 * 0;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1 * 1;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1 * 1;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1 * 1;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    cwt_dim1 = *mxcmu;
    cmu_dim1 = *mxcmu;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1 * 1;
    b_dim1 = *nnlyri;

    /* Function Body */
    zeroit_(b, nnlyri);
/*                              ** Construct B,  STWJ(20a,c) for */
/*                              ** parallel beam + bottom reflection + */
/*                              ** thermal emission at top and/or bottom */
    if (*mazim > 0 && *fbeam > 0.) {
/*                                         ** Azimuth-dependent case */
/*                                         ** (never called if FBEAM = 0) */
	if (*lyrcut || *lamber) {
/*               ** No azimuthal-dependent intensity for Lambert surface; */
/*               ** no intensity component for truncated bottom layer */
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
/*                                                  ** Top boundary */
		b[(i__2 = iq - 1) < 1 * b_dim1 && 0 <= i__2 ? i__2 : s_rnge(
			"b", i__2, "solve0_", (ftnlen)2821)] = -zz[*nn + 1 - 
			iq + zz_dim1 - zz_offset];
/*                                                  ** Bottom boundary */
		b[(i__2 = *ncol - *nn + iq - 1) < 1 * b_dim1 && 0 <= i__2 ? 
			i__2 : s_rnge("b", i__2, "solve0_", (ftnlen)2824)] = 
			-zz[iq + *nn + *ncut * zz_dim1 - zz_offset] * expbea[*
			ncut];
/* L10: */
	    }
	} else {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
		b[(i__2 = iq - 1) < 1 * b_dim1 && 0 <= i__2 ? i__2 : s_rnge(
			"b", i__2, "solve0_", (ftnlen)2833)] = -zz[*nn + 1 - 
			iq + zz_dim1 - zz_offset];
		sum = 0.;
		i__2 = *nn;
		for (jq = 1; jq <= i__2; ++jq) {
		    sum += cwt[(i__3 = jq - 1) < 1 * cwt_dim1 && 0 <= i__3 ? 
			    i__3 : s_rnge("cwt", i__3, "solve0_", (ftnlen)
			    2837)] * cmu[(i__4 = jq - 1) < 1 * cmu_dim1 && 0 
			    <= i__4 ? i__4 : s_rnge("cmu", i__4, "solve0_", (
			    ftnlen)2837)] * bdr[(i__5 = iq + jq * bdr_dim1 - 
			    bdr_offset) < 1 * bdr_dim1 * bdr_dim2 && 0 <= 
			    i__5 ? i__5 : s_rnge("bdr", i__5, "solve0_", (
			    ftnlen)2837)] * zz[*nn + 1 - jq + *ncut * zz_dim1 
			    - zz_offset] * expbea[*ncut];
/* L20: */
		}
		b[(i__2 = *ncol - *nn + iq - 1) < 1 * b_dim1 && 0 <= i__2 ? 
			i__2 : s_rnge("b", i__2, "solve0_", (ftnlen)2841)] = 
			sum;
		if (*fbeam > 0.) {
		    b[(i__2 = *ncol - *nn + iq - 1) < 1 * b_dim1 && 0 <= i__2 
			    ? i__2 : s_rnge("b", i__2, "solve0_", (ftnlen)
			    2842)] = sum + (bdr[(i__3 = iq - bdr_offset) < 1 *
			     bdr_dim1 * bdr_dim2 && 0 <= i__3 ? i__3 : s_rnge(
			    "bdr", i__3, "solve0_", (ftnlen)2842)] * *umu0 * *
			    fbeam / *pi - zz[iq + *nn + *ncut * zz_dim1 - 
			    zz_offset]) * expbea[*ncut];
		}
/* L30: */
	    }
	}
/*                             ** Continuity condition for layer */
/*                             ** interfaces of Eq. STWJ(20b) */
	it = *nn;
	i__1 = *ncut - 1;
	for (lc = 1; lc <= i__1; ++lc) {
	    i__2 = *nstr;
	    for (iq = 1; iq <= i__2; ++iq) {
		++it;
		b[(i__3 = it - 1) < 1 * b_dim1 && 0 <= i__3 ? i__3 : s_rnge(
			"b", i__3, "solve0_", (ftnlen)2857)] = (zz[iq + (lc + 
			1) * zz_dim1 - zz_offset] - zz[iq + lc * zz_dim1 - 
			zz_offset]) * expbea[lc];
/* L40: */
	    }
/* L50: */
	}
    } else {
/*                                   ** Azimuth-independent case */
	if (*fbeam == 0.) {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
/*                                      ** Top boundary */
		b[(i__2 = iq - 1) < 1 * b_dim1 && 0 <= i__2 ? i__2 : s_rnge(
			"b", i__2, "solve0_", (ftnlen)2871)] = -zplk0[*nn + 1 
			- iq + zplk0_dim1 - zplk0_offset] + *fisot + *tplank;
/* L60: */
	    }
	    if (*lyrcut) {
/*                               ** No intensity component for truncated */
/*                               ** bottom layer */
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
/*                                      ** Bottom boundary */
		    b[(i__2 = *ncol - *nn + iq - 1) < 1 * b_dim1 && 0 <= i__2 
			    ? i__2 : s_rnge("b", i__2, "solve0_", (ftnlen)
			    2882)] = -zplk0[iq + *nn + *ncut * zplk0_dim1 - 
			    zplk0_offset] - zplk1[iq + *nn + *ncut * 
			    zplk1_dim1 - zplk1_offset] * taucpr[*ncut];
/* L70: */
		}
	    } else {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    sum = 0.;
		    i__2 = *nn;
		    for (jq = 1; jq <= i__2; ++jq) {
			sum += cwt[(i__3 = jq - 1) < 1 * cwt_dim1 && 0 <= 
				i__3 ? i__3 : s_rnge("cwt", i__3, "solve0_", (
				ftnlen)2894)] * cmu[(i__4 = jq - 1) < 1 * 
				cmu_dim1 && 0 <= i__4 ? i__4 : s_rnge("cmu", 
				i__4, "solve0_", (ftnlen)2894)] * bdr[(i__5 = 
				iq + jq * bdr_dim1 - bdr_offset) < 1 * 
				bdr_dim1 * bdr_dim2 && 0 <= i__5 ? i__5 : 
				s_rnge("bdr", i__5, "solve0_", (ftnlen)2894)] 
				* (zplk0[*nn + 1 - jq + *ncut * zplk0_dim1 - 
				zplk0_offset] + zplk1[*nn + 1 - jq + *ncut * 
				zplk1_dim1 - zplk1_offset] * taucpr[*ncut]);
/* L80: */
		    }
		    b[(i__2 = *ncol - *nn + iq - 1) < 1 * b_dim1 && 0 <= i__2 
			    ? i__2 : s_rnge("b", i__2, "solve0_", (ftnlen)
			    2899)] = sum * 2. + bem[(i__3 = iq - 1) < 1 * 
			    bem_dim1 && 0 <= i__3 ? i__3 : s_rnge("bem", i__3,
			     "solve0_", (ftnlen)2899)] * *bplank - zplk0[iq + 
			    *nn + *ncut * zplk0_dim1 - zplk0_offset] - zplk1[
			    iq + *nn + *ncut * zplk1_dim1 - zplk1_offset] * 
			    taucpr[*ncut];
/* L90: */
		}
	    }
/*                             ** Continuity condition for layer */
/*                             ** interfaces, STWJ(20b) */
	    it = *nn;
	    i__1 = *ncut - 1;
	    for (lc = 1; lc <= i__1; ++lc) {
		i__2 = *nstr;
		for (iq = 1; iq <= i__2; ++iq) {
		    ++it;
		    b[(i__3 = it - 1) < 1 * b_dim1 && 0 <= i__3 ? i__3 : 
			    s_rnge("b", i__3, "solve0_", (ftnlen)2913)] = 
			    zplk0[iq + (lc + 1) * zplk0_dim1 - zplk0_offset] 
			    - zplk0[iq + lc * zplk0_dim1 - zplk0_offset] + (
			    zplk1[iq + (lc + 1) * zplk1_dim1 - zplk1_offset] 
			    - zplk1[iq + lc * zplk1_dim1 - zplk1_offset]) * 
			    taucpr[lc];
/* L100: */
		}
/* L110: */
	    }
	} else {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
		b[(i__2 = iq - 1) < 1 * b_dim1 && 0 <= i__2 ? i__2 : s_rnge(
			"b", i__2, "solve0_", (ftnlen)2924)] = -zz[*nn + 1 - 
			iq + zz_dim1 - zz_offset] - zplk0[*nn + 1 - iq + 
			zplk0_dim1 - zplk0_offset] + *fisot + *tplank;
/* L120: */
	    }
	    if (*lyrcut) {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    b[(i__2 = *ncol - *nn + iq - 1) < 1 * b_dim1 && 0 <= i__2 
			    ? i__2 : s_rnge("b", i__2, "solve0_", (ftnlen)
			    2931)] = -zz[iq + *nn + *ncut * zz_dim1 - 
			    zz_offset] * expbea[*ncut] - zplk0[iq + *nn + *
			    ncut * zplk0_dim1 - zplk0_offset] - zplk1[iq + *
			    nn + *ncut * zplk1_dim1 - zplk1_offset] * taucpr[*
			    ncut];
/* L130: */
		}
	    } else {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    sum = 0.;
		    i__2 = *nn;
		    for (jq = 1; jq <= i__2; ++jq) {
			sum += cwt[(i__3 = jq - 1) < 1 * cwt_dim1 && 0 <= 
				i__3 ? i__3 : s_rnge("cwt", i__3, "solve0_", (
				ftnlen)2943)] * cmu[(i__4 = jq - 1) < 1 * 
				cmu_dim1 && 0 <= i__4 ? i__4 : s_rnge("cmu", 
				i__4, "solve0_", (ftnlen)2943)] * bdr[(i__5 = 
				iq + jq * bdr_dim1 - bdr_offset) < 1 * 
				bdr_dim1 * bdr_dim2 && 0 <= i__5 ? i__5 : 
				s_rnge("bdr", i__5, "solve0_", (ftnlen)2943)] 
				* (zz[*nn + 1 - jq + *ncut * zz_dim1 - 
				zz_offset] * expbea[*ncut] + zplk0[*nn + 1 - 
				jq + *ncut * zplk0_dim1 - zplk0_offset] + 
				zplk1[*nn + 1 - jq + *ncut * zplk1_dim1 - 
				zplk1_offset] * taucpr[*ncut]);
/* L140: */
		    }
		    b[(i__2 = *ncol - *nn + iq - 1) < 1 * b_dim1 && 0 <= i__2 
			    ? i__2 : s_rnge("b", i__2, "solve0_", (ftnlen)
			    2949)] = sum * 2. + (bdr[(i__3 = iq - bdr_offset) 
			    < 1 * bdr_dim1 * bdr_dim2 && 0 <= i__3 ? i__3 : 
			    s_rnge("bdr", i__3, "solve0_", (ftnlen)2949)] * *
			    umu0 * *fbeam / *pi - zz[iq + *nn + *ncut * 
			    zz_dim1 - zz_offset]) * expbea[*ncut] + bem[(i__4 
			    = iq - 1) < 1 * bem_dim1 && 0 <= i__4 ? i__4 : 
			    s_rnge("bem", i__4, "solve0_", (ftnlen)2949)] * *
			    bplank - zplk0[iq + *nn + *ncut * zplk0_dim1 - 
			    zplk0_offset] - zplk1[iq + *nn + *ncut * 
			    zplk1_dim1 - zplk1_offset] * taucpr[*ncut];
/* L150: */
		}
	    }
	    it = *nn;
	    i__1 = *ncut - 1;
	    for (lc = 1; lc <= i__1; ++lc) {
		i__2 = *nstr;
		for (iq = 1; iq <= i__2; ++iq) {
		    ++it;
		    b[(i__3 = it - 1) < 1 * b_dim1 && 0 <= i__3 ? i__3 : 
			    s_rnge("b", i__3, "solve0_", (ftnlen)2966)] = (zz[
			    iq + (lc + 1) * zz_dim1 - zz_offset] - zz[iq + lc 
			    * zz_dim1 - zz_offset]) * expbea[lc] + zplk0[iq + 
			    (lc + 1) * zplk0_dim1 - zplk0_offset] - zplk0[iq 
			    + lc * zplk0_dim1 - zplk0_offset] + (zplk1[iq + (
			    lc + 1) * zplk1_dim1 - zplk1_offset] - zplk1[iq + 
			    lc * zplk1_dim1 - zplk1_offset]) * taucpr[lc];
/* L160: */
		}
/* L170: */
	    }
	}
    }
/*                     ** Find L-U (lower/upper triangular) decomposition */
/*                     ** of band matrix CBAND and test if it is nearly */
/*                     ** singular (note: CBAND is destroyed) */
/*                     ** (CBAND is in LINPACK packed format) */
    rcond = 0.;
    ncd = *nn * 3 - 1;
    sgbco_(cband, mi9m2, ncol, &ncd, &ncd, ipvt, &rcond, z__);
    if (rcond + 1. == 1.) {
	errmsg_("SOLVE0--SGBCO says matrix near singular", &c_false, (ftnlen)
		39);
    }
/*                   ** Solve linear system with coeff matrix CBAND */
/*                   ** and R.H. side(s) B after CBAND has been L-U */
/*                   ** decomposed.  Solution is returned in B. */
    sgbsl_(cband, mi9m2, ncol, &ncd, &ncd, ipvt, b, &c__0);
/*                   ** Zero CBAND (it may contain 'foreign' */
/*                   ** elements upon returning from LINPACK); */
/*                   ** necessary to prevent errors */
    i__1 = *mi9m2 * *nnlyri;
    zeroit_(cband, &i__1);
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	ipnt = lc * *nstr - *nn;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ll[*nn + 1 - iq + lc * ll_dim1 - ll_offset] = b[(i__3 = ipnt + 1 
		    - iq - 1) < 1 * b_dim1 && 0 <= i__3 ? i__3 : s_rnge("b", 
		    i__3, "solve0_", (ftnlen)3005)];
	    ll[iq + *nn + lc * ll_dim1 - ll_offset] = b[(i__3 = iq + ipnt - 1)
		     < 1 * b_dim1 && 0 <= i__3 ? i__3 : s_rnge("b", i__3, 
		    "solve0_", (ftnlen)3006)];
/* L180: */
	}
/* L190: */
    }
    return 0;
} /* solve0_ */

/* Subroutine */ int surfac_(doublereal *albedo, doublereal *delm0, 
	doublereal *fbeam, doublereal *hlpr, logical *lamber, integer *mi, 
	integer *mazim, integer *mxcmu, integer *mxumu, integer *nn, integer *
	numu, integer *nstr, logical *onlyfl, doublereal *umu, logical *
	usrang, doublereal *ylm0, doublereal *ylmc, doublereal *ylmu, 
	doublereal *bdr, doublereal *emu, doublereal *bem, doublereal *rmu, 
	doublereal *sqt)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    integer bdr_dim1, bdr_dim2, bdr_offset, bem_dim1, emu_dim1, hlpr_dim1, 
	    rmu_dim1, rmu_dim2, rmu_offset, ylm0_dim1, ylmc_dim1, ylmc_dim2, 
	    ylmc_offset, ylmu_dim1, ylmu_dim2, ylmu_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal dref, ylmg[1010]	/* was [101][10] */;
    static integer k, jg, iq, jq, iu;
    extern /* Subroutine */ int qgausn_(integer *, doublereal *, doublereal *)
	    , errmsg_(char *, logical *, ftnlen), lepoly_(integer *, integer *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *),
	     zeroit_(doublereal *, integer *);
    static doublereal sgn, gmu[10], gwt[10], sum;

/*       Specifies user's surface bidirectional properties, STWJ(21) */


/*   I N P U T     V A R I A B L E S: */

/*       DELM0  :  Kronecker delta, delta-sub-m0 */

/*       HLPR   :  Legendre moments of surface bidirectional reflectivity */
/*                 (with 2K+1 factor included) */

/*       MAZIM  :  Order of azimuthal component */

/*       NN     :  Order of double-Gauss quadrature (NSTR/2) */

/*       YLM0   :  Normalized associated Legendre polynomial */
/*                 at the beam angle */

/*       YLMC   :  Normalized associated Legendre polynomials */
/*                 at the quadrature angles */

/*       YLMU   :  Normalized associated Legendre polynomials */
/*                 at the user angles */

/*       SQT(k) :  Square root of k */

/*       (remainder are DISORT input variables) */


/*    O U T P U T     V A R I A B L E S: */

/*       BDR :  Surface bidirectional reflectivity (computational angles) */

/*       RMU :  Surface bidirectional reflectivity (user angles) */

/*       BEM :  Surface directional emissivity (computational angles) */

/*       EMU :  Surface directional emissivity (user angles) */


/*    I N T E R N A L     V A R I A B L E S: */

/*       DREF      Directional reflectivity */

/*       NMUG   :  Number of angle cosine quadrature points on (0,1) for */
/*                 integrating bidirectional reflectivity to get */
/*                 directional emissivity (it is necessary to use a */
/*                 quadrature set distinct from the computational */
/*                 angles, because the computational angles may not be */
/*                 dense enough--NSTR may be too small--to give an */
/*                 accurate approximation for the integration). */

/*       GMU    :  The NMUG angle cosine quadrature points on (0,1) */
/*       GWT    :  The NMUG angle cosine quadrature weights on (0,1) */

/*       YLMG   :  Normalized associated Legendre polynomials */
/*                 at the NMUG quadrature angles */

/*   Called by- DISORT */
/*   Calls- QGAUSN, LEPOLY, ZEROIT, ERRMSG */
/* +-------------------------------------------------------------------+ */
/*     .. Parameters .. */
/*                             ** CAUTION:  Do not increase MAXSTR */
/*                             **           without checking if this */
/*                             **           would require a larger */
/*                             **           dimension for SQT */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    /* Parameter adjustments */
    bem_dim1 = *mi;
    bdr_dim1 = *mi;
    bdr_dim2 = *mi - 0 + 1;
    bdr_offset = 1 + bdr_dim1 * 0;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_dim2 = *mxcmu;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylm0_dim1 = *mxcmu - 0 + 1;
    hlpr_dim1 = *mxcmu - 0 + 1;
    rmu_dim1 = *mxumu;
    rmu_dim2 = *mi - 0 + 1;
    rmu_offset = 1 + rmu_dim1 * 0;
    emu_dim1 = *mxumu;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_dim2 = *mxumu;
    ylmu_offset = 0 + ylmu_dim1 * 1;

    /* Function Body */
    if (pass1) {
	pass1 = FALSE_;
	qgausn_(&c__10, gmu, gwt);
	lepoly_(&c__10, &c__0, &c__100, &c__100, gmu, sqt, ylmg);
/*                       ** Convert Legendre polys. to negative GMU */
	sgn = -1.;
	for (k = 0; k <= 100; ++k) {
	    sgn = -sgn;
	    for (jg = 1; jg <= 10; ++jg) {
		ylmg[(i__1 = k + jg * 101 - 101) < 1010 && 0 <= i__1 ? i__1 : 
			s_rnge("ylmg", i__1, "surfac_", (ftnlen)3136)] = sgn *
			 ylmg[(i__2 = k + jg * 101 - 101) < 1010 && 0 <= i__2 
			? i__2 : s_rnge("ylmg", i__2, "surfac_", (ftnlen)3136)
			];
/* L10: */
	    }
/* L20: */
	}
    }
    i__1 = *mi * (*mi + 1);
    zeroit_(bdr, &i__1);
    zeroit_(bem, mi);
    if (*lamber && *mazim == 0) {
	i__1 = *nn;
	for (iq = 1; iq <= i__1; ++iq) {
	    bem[(i__2 = iq - 1) < 1 * bem_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "bem", i__2, "surfac_", (ftnlen)3151)] = 1. - *albedo;
	    i__2 = *nn;
	    for (jq = 0; jq <= i__2; ++jq) {
		bdr[(i__3 = iq + jq * bdr_dim1 - bdr_offset) < 1 * bdr_dim1 * 
			bdr_dim2 && 0 <= i__3 ? i__3 : s_rnge("bdr", i__3, 
			"surfac_", (ftnlen)3154)] = *albedo;
/* L30: */
	    }
/* L40: */
	}
    } else if (! (*lamber)) {
/*                                  ** Compute surface bidirectional */
/*                                  ** properties at computational angles */
	i__1 = *nn;
	for (iq = 1; iq <= i__1; ++iq) {
	    i__2 = *nn;
	    for (jq = 1; jq <= i__2; ++jq) {
		sum = 0.;
		i__3 = *nstr - 1;
		for (k = *mazim; k <= i__3; ++k) {
		    sum += hlpr[(i__4 = k) < 1 * hlpr_dim1 && 0 <= i__4 ? 
			    i__4 : s_rnge("hlpr", i__4, "surfac_", (ftnlen)
			    3169)] * ylmc[(i__5 = k + iq * ylmc_dim1 - 
			    ylmc_offset) < 1 * ylmc_dim1 * ylmc_dim2 && 0 <= 
			    i__5 ? i__5 : s_rnge("ylmc", i__5, "surfac_", (
			    ftnlen)3169)] * ylmc[(i__6 = k + (jq + *nn) * 
			    ylmc_dim1 - ylmc_offset) < 1 * ylmc_dim1 * 
			    ylmc_dim2 && 0 <= i__6 ? i__6 : s_rnge("ylmc", 
			    i__6, "surfac_", (ftnlen)3169)];
/* L50: */
		}
		bdr[(i__3 = iq + jq * bdr_dim1 - bdr_offset) < 1 * bdr_dim1 * 
			bdr_dim2 && 0 <= i__3 ? i__3 : s_rnge("bdr", i__3, 
			"surfac_", (ftnlen)3173)] = (2. - *delm0) * sum;
/* L60: */
	    }
	    if (*fbeam > 0.) {
		sum = 0.;
		i__2 = *nstr - 1;
		for (k = *mazim; k <= i__2; ++k) {
		    sum += hlpr[(i__3 = k) < 1 * hlpr_dim1 && 0 <= i__3 ? 
			    i__3 : s_rnge("hlpr", i__3, "surfac_", (ftnlen)
			    3182)] * ylmc[(i__4 = k + iq * ylmc_dim1 - 
			    ylmc_offset) < 1 * ylmc_dim1 * ylmc_dim2 && 0 <= 
			    i__4 ? i__4 : s_rnge("ylmc", i__4, "surfac_", (
			    ftnlen)3182)] * ylm0[(i__5 = k) < 1 * ylm0_dim1 &&
			     0 <= i__5 ? i__5 : s_rnge("ylm0", i__5, "surfac_"
			    , (ftnlen)3182)];
/* L70: */
		}
		bdr[(i__2 = iq - bdr_offset) < 1 * bdr_dim1 * bdr_dim2 && 0 <=
			 i__2 ? i__2 : s_rnge("bdr", i__2, "surfac_", (ftnlen)
			3185)] = (2. - *delm0) * sum;
	    }
/* L80: */
	}
	if (*mazim == 0) {
	    if (*nstr > 100) {
		errmsg_("SURFAC--parameter MAXSTR too small", &c_true, (
			ftnlen)34);
	    }
/*                              ** Integrate bidirectional reflectivity */
/*                              ** at reflection polar angles CMU and */
/*                              ** incident angles GMU to get */
/*                              ** directional emissivity at */
/*                              ** computational angles CMU. */
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
		dref = 0.;
		for (jg = 1; jg <= 10; ++jg) {
		    sum = 0.;
		    i__2 = *nstr - 1;
		    for (k = 0; k <= i__2; ++k) {
			sum += hlpr[(i__3 = k) < 1 * hlpr_dim1 && 0 <= i__3 ? 
				i__3 : s_rnge("hlpr", i__3, "surfac_", (
				ftnlen)3210)] * ylmc[(i__4 = k + iq * 
				ylmc_dim1 - ylmc_offset) < 1 * ylmc_dim1 * 
				ylmc_dim2 && 0 <= i__4 ? i__4 : s_rnge("ylmc",
				 i__4, "surfac_", (ftnlen)3210)] * ylmg[(i__5 
				= k + jg * 101 - 101) < 1010 && 0 <= i__5 ? 
				i__5 : s_rnge("ylmg", i__5, "surfac_", (
				ftnlen)3210)];
/* L90: */
		    }
		    dref += gwt[(i__2 = jg - 1) < 10 && 0 <= i__2 ? i__2 : 
			    s_rnge("gwt", i__2, "surfac_", (ftnlen)3214)] * 
			    2. * gmu[(i__3 = jg - 1) < 10 && 0 <= i__3 ? i__3 
			    : s_rnge("gmu", i__3, "surfac_", (ftnlen)3214)] * 
			    sum;
/* L100: */
		}
		bem[(i__2 = iq - 1) < 1 * bem_dim1 && 0 <= i__2 ? i__2 : 
			s_rnge("bem", i__2, "surfac_", (ftnlen)3218)] = 1. - 
			dref;
/* L110: */
	    }
	}
    }
/*                                       ** Compute surface bidirectional */
/*                                       ** properties at user angles */
    if (! (*onlyfl) && *usrang) {
	zeroit_(emu, mxumu);
	i__1 = *mxumu * (*mi + 1);
	zeroit_(rmu, &i__1);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    if (umu[iu - 1] > 0.) {
		if (*lamber && *mazim == 0) {
		    i__2 = *nn;
		    for (iq = 0; iq <= i__2; ++iq) {
			rmu[(i__3 = iu + iq * rmu_dim1 - rmu_offset) < 1 * 
				rmu_dim1 * rmu_dim2 && 0 <= i__3 ? i__3 : 
				s_rnge("rmu", i__3, "surfac_", (ftnlen)3240)] 
				= *albedo;
/* L120: */
		    }
		    emu[(i__2 = iu - 1) < 1 * emu_dim1 && 0 <= i__2 ? i__2 : 
			    s_rnge("emu", i__2, "surfac_", (ftnlen)3243)] = 
			    1. - *albedo;
		} else if (! (*lamber)) {
		    i__2 = *nn;
		    for (iq = 1; iq <= i__2; ++iq) {
			sum = 0.;
			i__3 = *nstr - 1;
			for (k = *mazim; k <= i__3; ++k) {
			    sum += hlpr[(i__4 = k) < 1 * hlpr_dim1 && 0 <= 
				    i__4 ? i__4 : s_rnge("hlpr", i__4, "surf\
ac_", (ftnlen)3252)] * ylmu[(i__5 = k + iu * ylmu_dim1 - ylmu_offset) < 1 * 
				    ylmu_dim1 * ylmu_dim2 && 0 <= i__5 ? i__5 
				    : s_rnge("ylmu", i__5, "surfac_", (ftnlen)
				    3252)] * ylmc[(i__6 = k + (iq + *nn) * 
				    ylmc_dim1 - ylmc_offset) < 1 * ylmc_dim1 *
				     ylmc_dim2 && 0 <= i__6 ? i__6 : s_rnge(
				    "ylmc", i__6, "surfac_", (ftnlen)3252)];
/* L130: */
			}
			rmu[(i__3 = iu + iq * rmu_dim1 - rmu_offset) < 1 * 
				rmu_dim1 * rmu_dim2 && 0 <= i__3 ? i__3 : 
				s_rnge("rmu", i__3, "surfac_", (ftnlen)3256)] 
				= (2. - *delm0) * sum;
/* L140: */
		    }
		    if (*fbeam > 0.) {
			sum = 0.;
			i__2 = *nstr - 1;
			for (k = *mazim; k <= i__2; ++k) {
			    sum += hlpr[(i__3 = k) < 1 * hlpr_dim1 && 0 <= 
				    i__3 ? i__3 : s_rnge("hlpr", i__3, "surf\
ac_", (ftnlen)3265)] * ylmu[(i__4 = k + iu * ylmu_dim1 - ylmu_offset) < 1 * 
				    ylmu_dim1 * ylmu_dim2 && 0 <= i__4 ? i__4 
				    : s_rnge("ylmu", i__4, "surfac_", (ftnlen)
				    3265)] * ylm0[(i__5 = k) < 1 * ylm0_dim1 
				    && 0 <= i__5 ? i__5 : s_rnge("ylm0", i__5,
				     "surfac_", (ftnlen)3265)];
/* L150: */
			}
			rmu[(i__2 = iu - rmu_offset) < 1 * rmu_dim1 * 
				rmu_dim2 && 0 <= i__2 ? i__2 : s_rnge("rmu", 
				i__2, "surfac_", (ftnlen)3268)] = (2. - *
				delm0) * sum;
		    }
		    if (*mazim == 0) {
/*                               ** Integrate bidirectional reflectivity */
/*                               ** at reflection angles UMU and */
/*                               ** incident angles GMU to get */
/*                               ** directional emissivity at */
/*                               ** user angles UMU. */
			dref = 0.;
			for (jg = 1; jg <= 10; ++jg) {
			    sum = 0.;
			    i__2 = *nstr - 1;
			    for (k = 0; k <= i__2; ++k) {
				sum += hlpr[(i__3 = k) < 1 * hlpr_dim1 && 0 <=
					 i__3 ? i__3 : s_rnge("hlpr", i__3, 
					"surfac_", (ftnlen)3286)] * ylmu[(
					i__4 = k + iu * ylmu_dim1 - 
					ylmu_offset) < 1 * ylmu_dim1 * 
					ylmu_dim2 && 0 <= i__4 ? i__4 : 
					s_rnge("ylmu", i__4, "surfac_", (
					ftnlen)3286)] * ylmg[(i__5 = k + jg * 
					101 - 101) < 1010 && 0 <= i__5 ? i__5 
					: s_rnge("ylmg", i__5, "surfac_", (
					ftnlen)3286)];
/* L160: */
			    }
			    dref += gwt[(i__2 = jg - 1) < 10 && 0 <= i__2 ? 
				    i__2 : s_rnge("gwt", i__2, "surfac_", (
				    ftnlen)3290)] * 2. * gmu[(i__3 = jg - 1) <
				     10 && 0 <= i__3 ? i__3 : s_rnge("gmu", 
				    i__3, "surfac_", (ftnlen)3290)] * sum;
/* L170: */
			}
			emu[(i__2 = iu - 1) < 1 * emu_dim1 && 0 <= i__2 ? 
				i__2 : s_rnge("emu", i__2, "surfac_", (ftnlen)
				3294)] = 1. - dref;
		    }
		}
	    }
/* L180: */
	}
    }
    return 0;
} /* surfac_ */

/* Subroutine */ int terpev_(doublereal *cwt, doublereal *evecc, doublereal *
	gl, doublereal *gu, integer *mazim, integer *mxcmu, integer *mxumu, 
	integer *nn, integer *nstr, integer *numu, doublereal *wk, doublereal 
	*ylmc, doublereal *ylmu)
{
    /* System generated locals */
    integer cwt_dim1, evecc_dim1, evecc_dim2, evecc_offset, gl_dim1, gu_dim1, 
	    gu_dim2, gu_offset, wk_dim1, ylmc_dim1, ylmc_dim2, ylmc_offset, 
	    ylmu_dim1, ylmu_dim2, ylmu_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer l, iq, jq, iu;
    static doublereal sum;

/*         Interpolate eigenvectors to user angles; Eq SD(8) */
/*   Called by- DISORT, ALBTRN */
/* --------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
    /* Parameter adjustments */
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_dim2 = *mxcmu;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    wk_dim1 = *mxcmu;
    gl_dim1 = *mxcmu - 0 + 1;
    evecc_dim1 = *mxcmu;
    evecc_dim2 = *mxcmu;
    evecc_offset = 1 + evecc_dim1 * 1;
    cwt_dim1 = *mxcmu;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_dim2 = *mxumu;
    ylmu_offset = 0 + ylmu_dim1 * 1;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * 1;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr - 1;
	for (l = *mazim; l <= i__2; ++l) {
/*                                   ** Inner sum in SD(8) times all */
/*                                   ** factors in outer sum but PLM(mu) */
	    sum = 0.;
	    i__3 = *nstr;
	    for (jq = 1; jq <= i__3; ++jq) {
		sum += cwt[(i__4 = jq - 1) < 1 * cwt_dim1 && 0 <= i__4 ? i__4 
			: s_rnge("cwt", i__4, "terpev_", (ftnlen)3342)] * 
			ylmc[(i__5 = l + jq * ylmc_dim1 - ylmc_offset) < 1 * 
			ylmc_dim1 * ylmc_dim2 && 0 <= i__5 ? i__5 : s_rnge(
			"ylmc", i__5, "terpev_", (ftnlen)3342)] * evecc[(i__6 
			= jq + iq * evecc_dim1 - evecc_offset) < 1 * 
			evecc_dim1 * evecc_dim2 && 0 <= i__6 ? i__6 : s_rnge(
			"evecc", i__6, "terpev_", (ftnlen)3342)];
/* L10: */
	    }
	    wk[(i__3 = l) < 1 * wk_dim1 && 0 <= i__3 ? i__3 : s_rnge("wk", 
		    i__3, "terpev_", (ftnlen)3345)] = gl[(i__4 = l) < 1 * 
		    gl_dim1 && 0 <= i__4 ? i__4 : s_rnge("gl", i__4, "terpev_"
		    , (ftnlen)3345)] * .5 * sum;
/* L20: */
	}
/*                                    ** Finish outer sum in SD(8) */
/*                                    ** and store eigenvectors */
	i__2 = *numu;
	for (iu = 1; iu <= i__2; ++iu) {
	    sum = 0.;
	    i__3 = *nstr - 1;
	    for (l = *mazim; l <= i__3; ++l) {
		sum += wk[(i__4 = l) < 1 * wk_dim1 && 0 <= i__4 ? i__4 : 
			s_rnge("wk", i__4, "terpev_", (ftnlen)3354)] * ylmu[(
			i__5 = l + iu * ylmu_dim1 - ylmu_offset) < 1 * 
			ylmu_dim1 * ylmu_dim2 && 0 <= i__5 ? i__5 : s_rnge(
			"ylmu", i__5, "terpev_", (ftnlen)3354)];
/* L30: */
	    }
	    if (iq <= *nn) {
		gu[(i__3 = iu + (iq + *nn) * gu_dim1 - gu_offset) < 1 * 
			gu_dim1 * gu_dim2 && 0 <= i__3 ? i__3 : s_rnge("gu", 
			i__3, "terpev_", (ftnlen)3357)] = sum;
	    }
	    if (iq > *nn) {
		gu[(i__3 = iu + (*nstr + 1 - iq) * gu_dim1 - gu_offset) < 1 * 
			gu_dim1 * gu_dim2 && 0 <= i__3 ? i__3 : s_rnge("gu", 
			i__3, "terpev_", (ftnlen)3358)] = sum;
	    }
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* terpev_ */

/* Subroutine */ int terpso_(doublereal *cwt, doublereal *delm0, doublereal *
	fbeam, doublereal *gl, integer *mazim, integer *mxcmu, logical *plank,
	 integer *numu, integer *nstr, doublereal *oprim, doublereal *pi, 
	doublereal *ylm0, doublereal *ylmc, doublereal *ylmu, doublereal *psi,
	 doublereal *xr0, doublereal *xr1, doublereal *z0, doublereal *zj, 
	doublereal *zbeam, doublereal *z0u, doublereal *z1u)
{
    /* System generated locals */
    integer cwt_dim1, gl_dim1, psi_dim1, ylm0_dim1, ylmc_dim1, ylmc_dim2, 
	    ylmc_offset, ylmu_dim1, ylmu_offset, z0_dim1, zj_dim1, i__1, i__2,
	     i__3, i__4, i__5;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal fact, psum;
    static integer iq, jq, iu;
    static doublereal sum;

/*         Interpolates source functions to user angles */


/*    I N P U T      V A R I A B L E S: */

/*       CWT    :  Weights for Gauss quadrature over angle cosine */

/*       DELM0  :  Kronecker delta, delta-sub-m0 */

/*       GL     :  Delta-M scaled Legendre coefficients of phase function */
/*                 (including factors 2L+1 and single-scatter albedo) */

/*       MAZIM  :  Order of azimuthal component */

/*       OPRIM  :  Single scattering albedo */

/*       XR0    :  Expansion of thermal source function */

/*       XR1    :  Expansion of thermal source function Eqs.SS(14-16) */

/*       YLM0   :  Normalized associated Legendre polynomial */
/*                 at the beam angle */

/*       YLMC   :  Normalized associated Legendre polynomial */
/*                 at the quadrature angles */

/*       YLMU   :  Normalized associated Legendre polynomial */
/*                 at the user angles */

/*       Z0     :  Solution vectors Z-sub-zero of Eq. SS(16) */

/*       ZJ     :  Solution vector Z-sub-zero after solving Eq. SS(19) */

/*       (remainder are DISORT input variables) */


/*    O U T P U T     V A R I A B L E S: */

/*       ZBEAM  :  Incident-beam source function at user angles */

/*       Z0U,Z1U:  Components of a linear-in-optical-depth-dependent */
/*                 source (approximating the Planck emission source) */


/*   I N T E R N A L    V A R I A B L E S: */

/*       PSI    :  Sum just after square bracket in  Eq. SD(9) */

/*   Called by- DISORT */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
    /* Parameter adjustments */
    zj_dim1 = *mxcmu;
    z0_dim1 = *mxcmu;
    psi_dim1 = *mxcmu;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1 * 1;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_dim2 = *mxcmu;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylm0_dim1 = *mxcmu - 0 + 1;
    gl_dim1 = *mxcmu - 0 + 1;
    cwt_dim1 = *mxcmu;

    /* Function Body */
    if (*fbeam > 0.) {
/*                                  ** Beam source terms; Eq. SD(9) */
	i__1 = *nstr - 1;
	for (iq = *mazim; iq <= i__1; ++iq) {
	    psum = 0.;
	    i__2 = *nstr;
	    for (jq = 1; jq <= i__2; ++jq) {
		psum += cwt[(i__3 = jq - 1) < 1 * cwt_dim1 && 0 <= i__3 ? 
			i__3 : s_rnge("cwt", i__3, "terpso_", (ftnlen)3450)] *
			 ylmc[(i__4 = iq + jq * ylmc_dim1 - ylmc_offset) < 1 *
			 ylmc_dim1 * ylmc_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			"ylmc", i__4, "terpso_", (ftnlen)3450)] * zj[(i__5 = 
			jq - 1) < 1 * zj_dim1 && 0 <= i__5 ? i__5 : s_rnge(
			"zj", i__5, "terpso_", (ftnlen)3450)];
/* L10: */
	    }
	    psi[(i__2 = iq) < 1 * psi_dim1 && 0 <= i__2 ? i__2 : s_rnge("psi",
		     i__2, "terpso_", (ftnlen)3453)] = gl[(i__3 = iq) < 1 * 
		    gl_dim1 && 0 <= i__3 ? i__3 : s_rnge("gl", i__3, "terpso_"
		    , (ftnlen)3453)] * .5 * psum;
/* L20: */
	}
	fact = (2. - *delm0) * *fbeam / (*pi * 4.);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    sum = 0.;
	    i__2 = *nstr - 1;
	    for (iq = *mazim; iq <= i__2; ++iq) {
		sum += ylmu[iq + iu * ylmu_dim1 - ylmu_offset] * (psi[(i__3 = 
			iq) < 1 * psi_dim1 && 0 <= i__3 ? i__3 : s_rnge("psi",
			 i__3, "terpso_", (ftnlen)3463)] + fact * gl[(i__4 = 
			iq) < 1 * gl_dim1 && 0 <= i__4 ? i__4 : s_rnge("gl", 
			i__4, "terpso_", (ftnlen)3463)] * ylm0[(i__5 = iq) < 
			1 * ylm0_dim1 && 0 <= i__5 ? i__5 : s_rnge("ylm0", 
			i__5, "terpso_", (ftnlen)3463)]);
/* L30: */
	    }
	    zbeam[iu - 1] = sum;
/* L40: */
	}
    }
    if (*plank && *mazim == 0) {
/*                                   ** Thermal source terms, STWJ(27c) */
	i__1 = *nstr - 1;
	for (iq = *mazim; iq <= i__1; ++iq) {
	    psum = 0.;
	    i__2 = *nstr;
	    for (jq = 1; jq <= i__2; ++jq) {
		psum += cwt[(i__3 = jq - 1) < 1 * cwt_dim1 && 0 <= i__3 ? 
			i__3 : s_rnge("cwt", i__3, "terpso_", (ftnlen)3481)] *
			 ylmc[(i__4 = iq + jq * ylmc_dim1 - ylmc_offset) < 1 *
			 ylmc_dim1 * ylmc_dim2 && 0 <= i__4 ? i__4 : s_rnge(
			"ylmc", i__4, "terpso_", (ftnlen)3481)] * z0[(i__5 = 
			jq - 1) < 1 * z0_dim1 && 0 <= i__5 ? i__5 : s_rnge(
			"z0", i__5, "terpso_", (ftnlen)3481)];
/* L50: */
	    }
	    psi[(i__2 = iq) < 1 * psi_dim1 && 0 <= i__2 ? i__2 : s_rnge("psi",
		     i__2, "terpso_", (ftnlen)3484)] = gl[(i__3 = iq) < 1 * 
		    gl_dim1 && 0 <= i__3 ? i__3 : s_rnge("gl", i__3, "terpso_"
		    , (ftnlen)3484)] * .5 * psum;
/* L60: */
	}
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    sum = 0.;
	    i__2 = *nstr - 1;
	    for (iq = *mazim; iq <= i__2; ++iq) {
		sum += ylmu[iq + iu * ylmu_dim1 - ylmu_offset] * psi[(i__3 = 
			iq) < 1 * psi_dim1 && 0 <= i__3 ? i__3 : s_rnge("psi",
			 i__3, "terpso_", (ftnlen)3492)];
/* L70: */
	    }
	    z0u[iu - 1] = sum + (1. - *oprim) * *xr0;
	    z1u[iu - 1] = *xr1;
/* L80: */
	}
    }
    return 0;
} /* terpso_ */

/* Subroutine */ int upbeam_(doublereal *array, doublereal *cc, doublereal *
	cmu, doublereal *delm0, doublereal *fbeam, doublereal *gl, integer *
	ipvt, integer *mazim, integer *mxcmu, integer *nn, integer *nstr, 
	doublereal *pi, doublereal *umu0, doublereal *wk, doublereal *ylm0, 
	doublereal *ylmc, doublereal *zj, doublereal *zz)
{
    /* System generated locals */
    integer array_dim1, array_dim2, array_offset, cc_dim1, cc_dim2, cc_offset,
	     cmu_dim1, gl_dim1, ylm0_dim1, ylmc_dim1, ylmc_offset, zj_dim1, 
	    zz_dim1, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int sgeco_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static doublereal rcond;
    extern /* Subroutine */ int sgesl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *);
    static integer iq, jq;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static integer job;
    static doublereal sum;

/*         Finds the incident-beam particular solution of SS(18) */


/*   I N P U T    V A R I A B L E S: */

/*       CC     :  C-sub-ij in Eq. SS(5) */

/*       CMU    :  Abscissae for Gauss quadrature over angle cosine */

/*       DELM0  :  Kronecker delta, delta-sub-m0 */

/*       GL     :  Delta-M scaled Legendre coefficients of phase function */
/*                 (including factors 2L+1 and single-scatter albedo) */

/*       MAZIM  :  Order of azimuthal component */

/*       YLM0   :  Normalized associated Legendre polynomial */
/*                 at the beam angle */

/*       YLMC   :  Normalized associated Legendre polynomial */
/*                 at the quadrature angles */

/*       (remainder are DISORT input variables) */


/*   O U T P U T    V A R I A B L E S: */

/*       ZJ     :  Right-hand side vector X-sub-zero in SS(19); also the */
/*                 solution vector Z-sub-zero after solving that system */

/*       ZZ     :  Permanent storage for ZJ, but re-ordered */


/*   I N T E R N A L    V A R I A B L E S: */

/*       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19) */
/*       IPVT   :  Integer vector of pivot indices required by LINPACK */
/*       WK     :  Scratch array required by LINPACK */

/*   Called by- DISORT */
/*   Calls- SGECO, ERRMSG, SGESL */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    /* Parameter adjustments */
    zz_dim1 = *mxcmu;
    zj_dim1 = *mxcmu;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylm0_dim1 = *mxcmu - 0 + 1;
    gl_dim1 = *mxcmu - 0 + 1;
    cmu_dim1 = *mxcmu;
    cc_dim1 = *mxcmu;
    cc_dim2 = *mxcmu;
    cc_offset = 1 + cc_dim1 * 1;
    array_dim1 = *mxcmu;
    array_dim2 = *mxcmu;
    array_offset = 1 + array_dim1 * 1;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    array[(i__3 = iq + jq * array_dim1 - array_offset) < 1 * 
		    array_dim1 * array_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "array", i__3, "upbeam_", (ftnlen)3579)] = -cc[(i__4 = iq 
		    + jq * cc_dim1 - cc_offset) < 1 * cc_dim1 * cc_dim2 && 0 
		    <= i__4 ? i__4 : s_rnge("cc", i__4, "upbeam_", (ftnlen)
		    3579)];
/* L10: */
	}
	array[(i__2 = iq + iq * array_dim1 - array_offset) < 1 * array_dim1 * 
		array_dim2 && 0 <= i__2 ? i__2 : s_rnge("array", i__2, "upbe\
am_", (ftnlen)3582)] = cmu[(i__3 = iq - 1) < 1 * cmu_dim1 && 0 <= i__3 ? i__3 
		: s_rnge("cmu", i__3, "upbeam_", (ftnlen)3582)] / *umu0 + 1. 
		+ array[(i__4 = iq + iq * array_dim1 - array_offset) < 1 * 
		array_dim1 * array_dim2 && 0 <= i__4 ? i__4 : s_rnge("array", 
		i__4, "upbeam_", (ftnlen)3582)];
	sum = 0.;
	i__2 = *nstr - 1;
	for (k = *mazim; k <= i__2; ++k) {
	    sum += gl[(i__3 = k) < 1 * gl_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		    "gl", i__3, "upbeam_", (ftnlen)3586)] * ylmc[k + iq * 
		    ylmc_dim1 - ylmc_offset] * ylm0[(i__4 = k) < 1 * 
		    ylm0_dim1 && 0 <= i__4 ? i__4 : s_rnge("ylm0", i__4, 
		    "upbeam_", (ftnlen)3586)];
/* L20: */
	}
	zj[(i__2 = iq - 1) < 1 * zj_dim1 && 0 <= i__2 ? i__2 : s_rnge("zj", 
		i__2, "upbeam_", (ftnlen)3589)] = (2. - *delm0) * *fbeam * 
		sum / (*pi * 4.);
/* L30: */
    }
/*                  ** Find L-U (lower/upper triangular) decomposition */
/*                  ** of ARRAY and see if it is nearly singular */
/*                  ** (NOTE:  ARRAY is altered) */
    rcond = 0.;
    sgeco_(array, mxcmu, nstr, ipvt, &rcond, wk);
    if (rcond + 1. == 1.) {
	errmsg_("UPBEAM--SGECO says matrix near singular", &c_false, (ftnlen)
		39);
    }
/*                ** Solve linear system with coeff matrix ARRAY */
/*                ** (assumed already L-U decomposed) and R.H. side(s) */
/*                ** ZJ;  return solution(s) in ZJ */
    job = 0;
    sgesl_(array, mxcmu, nstr, ipvt, zj, &job);
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zz[(i__2 = iq + *nn - 1) < 1 * zz_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"zz", i__2, "upbeam_", (ftnlen)3610)] = zj[(i__3 = iq - 1) < 
		1 * zj_dim1 && 0 <= i__3 ? i__3 : s_rnge("zj", i__3, "upbeam_"
		, (ftnlen)3610)];
	zz[(i__2 = *nn + 1 - iq - 1) < 1 * zz_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("zz", i__2, "upbeam_", (ftnlen)3611)] = zj[(i__3 = iq 
		+ *nn - 1) < 1 * zj_dim1 && 0 <= i__3 ? i__3 : s_rnge("zj", 
		i__3, "upbeam_", (ftnlen)3611)];
/* L40: */
    }
    return 0;
} /* upbeam_ */

/* Subroutine */ int upisot_(doublereal *array, doublereal *cc, doublereal *
	cmu, integer *ipvt, integer *mxcmu, integer *nn, integer *nstr, 
	doublereal *oprim, doublereal *wk, doublereal *xr0, doublereal *xr1, 
	doublereal *z0, doublereal *z1, doublereal *zplk0, doublereal *zplk1)
{
    /* System generated locals */
    integer array_dim1, array_dim2, array_offset, cc_dim1, cc_dim2, cc_offset,
	     cmu_dim1, z0_dim1, z1_dim1, zplk0_dim1, zplk1_dim1, i__1, i__2, 
	    i__3, i__4;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    extern /* Subroutine */ int sgeco_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static doublereal rcond;
    extern /* Subroutine */ int sgesl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *);
    static integer iq, jq;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);

/*       Finds the particular solution of thermal radiation of SS(15) */


/*    I N P U T     V A R I A B L E S: */

/*       CC     :  C-sub-ij in Eq. SS(5) */

/*       CMU    :  Abscissae for Gauss quadrature over angle cosine */

/*       OPRIM  :  Delta-M scaled single scattering albedo */

/*       XR0    :  Expansion of thermal source function */

/*       XR1    :  Expansion of thermal source function Eqs. SS(14-16) */

/*       (remainder are DISORT input variables) */


/*    O U T P U T    V A R I A B L E S: */

/*       Z0     :  Solution vectors Z-sub-zero of Eq. SS(16) */

/*       Z1     :  Solution vectors Z-sub-one  of Eq. SS(16) */

/*       ZPLK0, :  Permanent storage for Z0,Z1, but re-ordered */
/*        ZPLK1 */


/*   I N T E R N A L    V A R I A B L E S: */

/*       ARRAY  :  Coefficient matrix in left-hand side of EQ. SS(16) */
/*       IPVT   :  Integer vector of pivot indices required by LINPACK */
/*       WK     :  Scratch array required by LINPACK */

/*   Called by- DISORT */
/*   Calls- SGECO, ERRMSG, SGESL */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    /* Parameter adjustments */
    zplk1_dim1 = *mxcmu;
    zplk0_dim1 = *mxcmu;
    z1_dim1 = *mxcmu;
    z0_dim1 = *mxcmu;
    cmu_dim1 = *mxcmu;
    cc_dim1 = *mxcmu;
    cc_dim2 = *mxcmu;
    cc_offset = 1 + cc_dim1 * 1;
    array_dim1 = *mxcmu;
    array_dim2 = *mxcmu;
    array_offset = 1 + array_dim1 * 1;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    array[(i__3 = iq + jq * array_dim1 - array_offset) < 1 * 
		    array_dim1 * array_dim2 && 0 <= i__3 ? i__3 : s_rnge(
		    "array", i__3, "upisot_", (ftnlen)3685)] = -cc[(i__4 = iq 
		    + jq * cc_dim1 - cc_offset) < 1 * cc_dim1 * cc_dim2 && 0 
		    <= i__4 ? i__4 : s_rnge("cc", i__4, "upisot_", (ftnlen)
		    3685)];
/* L10: */
	}
	array[(i__2 = iq + iq * array_dim1 - array_offset) < 1 * array_dim1 * 
		array_dim2 && 0 <= i__2 ? i__2 : s_rnge("array", i__2, "upis\
ot_", (ftnlen)3688)] = array[(i__3 = iq + iq * array_dim1 - array_offset) < 1 
		* array_dim1 * array_dim2 && 0 <= i__3 ? i__3 : s_rnge("array"
		, i__3, "upisot_", (ftnlen)3688)] + 1.;
	z1[(i__2 = iq - 1) < 1 * z1_dim1 && 0 <= i__2 ? i__2 : s_rnge("z1", 
		i__2, "upisot_", (ftnlen)3690)] = *xr1;
	z0[(i__2 = iq - 1) < 1 * z0_dim1 && 0 <= i__2 ? i__2 : s_rnge("z0", 
		i__2, "upisot_", (ftnlen)3691)] = (1. - *oprim) * *xr0 + cmu[(
		i__3 = iq - 1) < 1 * cmu_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"cmu", i__3, "upisot_", (ftnlen)3691)] * z1[(i__4 = iq - 1) < 
		1 * z1_dim1 && 0 <= i__4 ? i__4 : s_rnge("z1", i__4, "upisot_"
		, (ftnlen)3691)];
/* L20: */
    }
/*                       ** Solve linear equations: same as in UPBEAM, */
/*                       ** except ZJ replaced by Z0 */
    rcond = 0.;
    sgeco_(array, mxcmu, nstr, ipvt, &rcond, wk);
    if (rcond + 1. == 1.) {
	errmsg_("UPISOT--SGECO says matrix near singular", &c_false, (ftnlen)
		39);
    }
    sgesl_(array, mxcmu, nstr, ipvt, z0, &c__0);
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zplk0[(i__2 = iq + *nn - 1) < 1 * zplk0_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("zplk0", i__2, "upisot_", (ftnlen)3706)] = z0[(i__3 = 
		iq - 1) < 1 * z0_dim1 && 0 <= i__3 ? i__3 : s_rnge("z0", i__3,
		 "upisot_", (ftnlen)3706)];
	zplk1[(i__2 = iq + *nn - 1) < 1 * zplk1_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("zplk1", i__2, "upisot_", (ftnlen)3707)] = z1[(i__3 = 
		iq - 1) < 1 * z1_dim1 && 0 <= i__3 ? i__3 : s_rnge("z1", i__3,
		 "upisot_", (ftnlen)3707)];
	zplk0[(i__2 = *nn + 1 - iq - 1) < 1 * zplk0_dim1 && 0 <= i__2 ? i__2 :
		 s_rnge("zplk0", i__2, "upisot_", (ftnlen)3708)] = z0[(i__3 = 
		iq + *nn - 1) < 1 * z0_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"z0", i__3, "upisot_", (ftnlen)3708)];
	zplk1[(i__2 = *nn + 1 - iq - 1) < 1 * zplk1_dim1 && 0 <= i__2 ? i__2 :
		 s_rnge("zplk1", i__2, "upisot_", (ftnlen)3709)] = z1[(i__3 = 
		iq + *nn - 1) < 1 * z1_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"z1", i__3, "upisot_", (ftnlen)3709)];
/* L30: */
    }
    return 0;
} /* upisot_ */

/* Subroutine */ int usrint_(doublereal *bplank, doublereal *cmu, doublereal *
	cwt, doublereal *delm0, doublereal *dtaucp, doublereal *emu, 
	doublereal *expbea, doublereal *fbeam, doublereal *fisot, doublereal *
	gc, doublereal *gu, doublereal *kk, logical *lamber, integer *layru, 
	doublereal *ll, logical *lyrcut, integer *mazim, integer *mxcmu, 
	integer *mxulv, integer *mxumu, integer *ncut, integer *nlyr, integer 
	*nn, integer *nstr, logical *plank, integer *numu, integer *ntau, 
	doublereal *pi, doublereal *rmu, doublereal *taucpr, doublereal *
	tplank, doublereal *umu, doublereal *umu0, doublereal *utaupr, 
	doublereal *wk, doublereal *zbeam, doublereal *z0u, doublereal *z1u, 
	doublereal *zz, doublereal *zplk0, doublereal *zplk1, doublereal *uum)
{
    /* System generated locals */
    integer cmu_dim1, cwt_dim1, emu_dim1, gc_dim1, gc_dim2, gc_offset, 
	    gu_dim1, gu_dim2, gu_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, rmu_dim1, rmu_offset, utaupr_dim1, uum_dim1, uum_dim2, 
	    uum_offset, wk_dim1, z0u_dim1, z0u_offset, z1u_dim1, z1u_offset, 
	    zbeam_dim1, zbeam_offset, zplk0_dim1, zplk0_offset, zplk1_dim1, 
	    zplk1_offset, zz_dim1, zz_offset, i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    double exp(doublereal);

    /* Local variables */
    static doublereal fact, dtau, expn, dtau1, dtau2, denom;
    static integer lc, iq, jq, iu, lu;
    static doublereal bnddfu, bnddir, bndint, palint, dfuint;
    static integer lyrend;
    static logical negumu;
    static doublereal plkint;
    static integer lyrstr;
    static doublereal sgn;
    static integer lyu;
    static doublereal exp0, exp1, exp2;

/*       Computes intensity components at user output angles */
/*       for azimuthal expansion terms in Eq. SD(2) */


/*   I N P U T    V A R I A B L E S: */

/*       BPLANK :  Integrated Planck function for emission from */
/*                 bottom boundary */

/*       CMU    :  Abscissae for Gauss quadrature over angle cosine */

/*       CWT    :  Weights for Gauss quadrature over angle cosine */

/*       DELM0  :  Kronecker delta, delta-sub-M0 */

/*       EMU    :  Surface directional emissivity (user angles) */

/*       EXPBEA :  Transmission of incident beam, EXP(-TAUCPR/UMU0) */

/*       GC     :  Eigenvectors at polar quadrature angles, SC(1) */

/*       GU     :  Eigenvectors interpolated to user polar angles */
/*                    (i.e., G in Eq. SC(1) ) */

/*       KK     :  Eigenvalues of coeff. matrix in Eq. SS(7) */

/*       LAYRU  :  Layer number of user level UTAU */

/*       LL     :  Constants of integration in Eq. SC(1), obtained */
/*                 by solving scaled version of Eq. SC(5); */
/*                 exponential term of Eq. SC(12) not included */

/*       LYRCUT :  Logical flag for truncation of computational layer */

/*       MAZIM  :  Order of azimuthal component */

/*       NCUT   :  Total number of computational layers considered */

/*       NN     :  Order of double-Gauss quadrature (NSTR/2) */

/*       RMU    :  Surface bidirectional reflectivity (user angles) */

/*       TAUCPR :  Cumulative optical depth (delta-M-Scaled) */

/*       TPLANK :  Integrated Planck function for emission from */
/*                 top boundary */

/*       UTAUPR :  Optical depths of user output levels in delta-M */
/*                 coordinates;  equal to UTAU if no delta-M */

/*       Z0U    :  Z-sub-zero in Eq. SS(16) interpolated to user */
/*                 angles from an equation derived from SS(16) */

/*       Z1U    :  Z-sub-one in Eq. SS(16) interpolated to user */
/*                 angles from an equation derived from SS(16) */

/*       ZZ     :  Beam source vectors in Eq. SS(19) */

/*       ZPLK0  :  Thermal source vectors Z0, by solving Eq. SS(16) */

/*       ZPLK1  :  Thermal source vectors Z1, by solving Eq. SS(16) */

/*       ZBEAM  :  Incident-beam source vectors */

/*       (Remainder are DISORT input variables) */


/*    O U T P U T    V A R I A B L E S: */

/*       UUM    :  Azimuthal components of the intensity in EQ. STWJ(5) */


/*    I N T E R N A L    V A R I A B L E S: */

/*       BNDDIR :  Direct intensity down at the bottom boundary */
/*       BNDDFU :  Diffuse intensity down at the bottom boundary */
/*       BNDINT :  Intensity attenuated at both boundaries, STWJ(25-6) */
/*       DTAU   :  Optical depth of a computational layer */
/*       LYREND :  End layer of integration */
/*       LYRSTR :  Start layer of integration */
/*       PALINT :  Intensity component from parallel beam */
/*       PLKINT :  Intensity component from planck source */
/*       WK     :  Scratch vector for saving EXP evaluations */

/*       All the exponential factors ( EXP1, EXPN,... etc.) */
/*       come from the substitution of constants of integration in */
/*       Eq. SC(12) into Eqs. S1(8-9).  They all have negative */
/*       arguments so there should never be overflow problems. */

/*   Called by- DISORT */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*                          ** Incorporate constants of integration into */
/*                          ** interpolated eigenvectors */
    /* Parameter adjustments */
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1 * 1;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1 * 1;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1 * 1;
    wk_dim1 = *mxcmu;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    cwt_dim1 = *mxcmu;
    cmu_dim1 = *mxcmu;
    utaupr_dim1 = *mxulv;
    uum_dim1 = *mxumu;
    uum_dim2 = *mxulv;
    uum_offset = 1 + uum_dim1 * 1;
    z1u_dim1 = *mxumu;
    z1u_offset = 1 + z1u_dim1 * 1;
    z0u_dim1 = *mxumu;
    z0u_offset = 1 + z0u_dim1 * 1;
    zbeam_dim1 = *mxumu;
    zbeam_offset = 1 + zbeam_dim1 * 1;
    rmu_dim1 = *mxumu;
    rmu_offset = 1 + rmu_dim1 * 0;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2 * 1);
    emu_dim1 = *mxumu;

    /* Function Body */
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	i__2 = *nstr;
	for (iq = 1; iq <= i__2; ++iq) {
	    i__3 = *numu;
	    for (iu = 1; iu <= i__3; ++iu) {
		gu[iu + (iq + lc * gu_dim2) * gu_dim1 - gu_offset] *= ll[iq + 
			lc * ll_dim1 - ll_offset];
/* L10: */
	    }
/* L20: */
	}
/* L30: */
    }
/*                           ** Loop over levels at which intensities */
/*                           ** are desired ('user output levels') */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	if (*fbeam > 0.) {
	    exp0 = exp(-utaupr[(i__2 = lu - 1) < 1 * utaupr_dim1 && 0 <= i__2 
		    ? i__2 : s_rnge("utaupr", i__2, "usrint_", (ftnlen)3862)] 
		    / *umu0);
	}
	lyu = layru[lu - 1];
/*                              ** Loop over polar angles at which */
/*                              ** intensities are desired */
	i__2 = *numu;
	for (iu = 1; iu <= i__2; ++iu) {
	    if (*lyrcut && lyu > *ncut) {
		goto L150;
	    }
	    negumu = umu[iu - 1] < 0.;
	    if (negumu) {
		lyrstr = 1;
		lyrend = lyu - 1;
		sgn = -1.;
	    } else {
		lyrstr = lyu + 1;
		lyrend = *ncut;
		sgn = 1.;
	    }
/*                          ** For downward intensity, integrate from top */
/*                          ** to LYU-1 in Eq. S1(8); for upward, */
/*                          ** integrate from bottom to LYU+1 in S1(9) */
	    palint = 0.;
	    plkint = 0.;
	    i__3 = lyrend;
	    for (lc = lyrstr; lc <= i__3; ++lc) {
		dtau = dtaucp[lc - 1];
		exp1 = exp((utaupr[(i__4 = lu - 1) < 1 * utaupr_dim1 && 0 <= 
			i__4 ? i__4 : s_rnge("utaupr", i__4, "usrint_", (
			ftnlen)3894)] - taucpr[lc - 1]) / umu[iu - 1]);
		exp2 = exp((utaupr[(i__4 = lu - 1) < 1 * utaupr_dim1 && 0 <= 
			i__4 ? i__4 : s_rnge("utaupr", i__4, "usrint_", (
			ftnlen)3895)] - taucpr[lc]) / umu[iu - 1]);
		if (*plank && *mazim == 0) {
		    plkint += sgn * (z0u[iu + lc * z0u_dim1 - z0u_offset] * (
			    exp1 - exp2) + z1u[iu + lc * z1u_dim1 - 
			    z1u_offset] * ((taucpr[lc - 1] + umu[iu - 1]) * 
			    exp1 - (taucpr[lc] + umu[iu - 1]) * exp2));
		}
		if (*fbeam > 0.) {
		    denom = umu[iu - 1] / *umu0 + 1.;
		    if (abs(denom) < 1e-4) {
/*                                                   ** L'Hospital limit */
			expn = dtau / *umu0 * exp0;
		    } else {
			expn = (exp1 * expbea[lc - 1] - exp2 * expbea[lc]) * 
				sgn / denom;
		    }
		    palint += zbeam[iu + lc * zbeam_dim1 - zbeam_offset] * 
			    expn;
		}
/*                                                   ** KK is negative */
		i__4 = *nn;
		for (iq = 1; iq <= i__4; ++iq) {
		    wk[(i__5 = iq - 1) < 1 * wk_dim1 && 0 <= i__5 ? i__5 : 
			    s_rnge("wk", i__5, "usrint_", (ftnlen)3923)] = 
			    exp(kk[iq + lc * kk_dim1 - kk_offset] * dtau);
		    denom = umu[iu - 1] * kk[iq + lc * kk_dim1 - kk_offset] + 
			    1.;
		    if (abs(denom) < 1e-4) {
/*                                                   ** L'Hospital limit */
			expn = dtau / umu[iu - 1] * exp2;
		    } else {
			expn = sgn * (exp1 * wk[(i__5 = iq - 1) < 1 * wk_dim1 
				&& 0 <= i__5 ? i__5 : s_rnge("wk", i__5, 
				"usrint_", (ftnlen)3932)] - exp2) / denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1 - 
			    gu_offset] * expn;
/* L40: */
		}
/*                                                   ** KK is positive */
		i__4 = *nstr;
		for (iq = *nn + 1; iq <= i__4; ++iq) {
		    denom = umu[iu - 1] * kk[iq + lc * kk_dim1 - kk_offset] + 
			    1.;
		    if (abs(denom) < 1e-4) {
/*                                                   ** L'Hospital limit */
			expn = -dtau / umu[iu - 1] * exp1;
		    } else {
			expn = sgn * (exp1 - exp2 * wk[(i__5 = *nstr + 1 - iq 
				- 1) < 1 * wk_dim1 && 0 <= i__5 ? i__5 : 
				s_rnge("wk", i__5, "usrint_", (ftnlen)3951)]) 
				/ denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1 - 
			    gu_offset] * expn;
/* L50: */
		}
/* L60: */
	    }
/*                           ** Calculate contribution from user */
/*                           ** output level to next computational level */
	    dtau1 = utaupr[(i__3 = lu - 1) < 1 * utaupr_dim1 && 0 <= i__3 ? 
		    i__3 : s_rnge("utaupr", i__3, "usrint_", (ftnlen)3964)] - 
		    taucpr[lyu - 1];
	    dtau2 = utaupr[(i__3 = lu - 1) < 1 * utaupr_dim1 && 0 <= i__3 ? 
		    i__3 : s_rnge("utaupr", i__3, "usrint_", (ftnlen)3965)] - 
		    taucpr[lyu];
	    if (abs(dtau1) < 1e-6 && negumu) {
		goto L90;
	    }
	    if (abs(dtau2) < 1e-6 && ! negumu) {
		goto L90;
	    }
	    if (negumu) {
		exp1 = exp(dtau1 / umu[iu - 1]);
	    }
	    if (! negumu) {
		exp2 = exp(dtau2 / umu[iu - 1]);
	    }
	    if (*fbeam > 0.) {
		denom = umu[iu - 1] / *umu0 + 1.;
		if (abs(denom) < 1e-4) {
		    expn = dtau1 / *umu0 * exp0;
		} else if (negumu) {
		    expn = (exp0 - expbea[lyu - 1] * exp1) / denom;
		} else {
		    expn = (exp0 - expbea[lyu] * exp2) / denom;
		}
		palint += zbeam[iu + lyu * zbeam_dim1 - zbeam_offset] * expn;
	    }
/*                                                   ** KK is negative */
	    dtau = dtaucp[lyu - 1];
	    i__3 = *nn;
	    for (iq = 1; iq <= i__3; ++iq) {
		denom = umu[iu - 1] * kk[iq + lyu * kk_dim1 - kk_offset] + 1.;
		if (abs(denom) < 1e-4) {
		    expn = -dtau2 / umu[iu - 1] * exp2;
		} else if (negumu) {
		    expn = (exp(-kk[iq + lyu * kk_dim1 - kk_offset] * dtau2) 
			    - exp(kk[iq + lyu * kk_dim1 - kk_offset] * dtau) *
			     exp1) / denom;
		} else {
		    expn = (exp(-kk[iq + lyu * kk_dim1 - kk_offset] * dtau2) 
			    - exp2) / denom;
		}
		palint += gu[iu + (iq + lyu * gu_dim2) * gu_dim1 - gu_offset] 
			* expn;
/* L70: */
	    }
/*                                                   ** KK is positive */
	    i__3 = *nstr;
	    for (iq = *nn + 1; iq <= i__3; ++iq) {
		denom = umu[iu - 1] * kk[iq + lyu * kk_dim1 - kk_offset] + 1.;
		if (abs(denom) < 1e-4) {
		    expn = -dtau1 / umu[iu - 1] * exp1;
		} else if (negumu) {
		    expn = (exp(-kk[iq + lyu * kk_dim1 - kk_offset] * dtau1) 
			    - exp1) / denom;
		} else {
		    expn = (exp(-kk[iq + lyu * kk_dim1 - kk_offset] * dtau1) 
			    - exp(-kk[iq + lyu * kk_dim1 - kk_offset] * dtau) 
			    * exp2) / denom;
		}
		palint += gu[iu + (iq + lyu * gu_dim2) * gu_dim1 - gu_offset] 
			* expn;
/* L80: */
	    }
	    if (*plank && *mazim == 0) {
		if (negumu) {
		    expn = exp1;
		    fact = taucpr[lyu - 1] + umu[iu - 1];
		} else {
		    expn = exp2;
		    fact = taucpr[lyu] + umu[iu - 1];
		}
		plkint = plkint + z0u[iu + lyu * z0u_dim1 - z0u_offset] * (1. 
			- expn) + z1u[iu + lyu * z1u_dim1 - z1u_offset] * (
			utaupr[(i__3 = lu - 1) < 1 * utaupr_dim1 && 0 <= i__3 
			? i__3 : s_rnge("utaupr", i__3, "usrint_", (ftnlen)
			4059)] + umu[iu - 1] - fact * expn);
	    }
/*                            ** Calculate intensity components */
/*                            ** attenuated at both boundaries. */
/*                            ** NOTE: no azimuthal intensity */
/*                            ** component for isotropic surface */
L90:
	    bndint = 0.;
	    if (negumu && *mazim == 0) {
		bndint = (*fisot + *tplank) * exp(utaupr[(i__3 = lu - 1) < 1 *
			 utaupr_dim1 && 0 <= i__3 ? i__3 : s_rnge("utaupr", 
			i__3, "usrint_", (ftnlen)4073)] / umu[iu - 1]);
	    } else if (! negumu) {
		if (*lyrcut || *lamber && *mazim > 0) {
		    goto L140;
		}
		i__3 = *nstr;
		for (jq = *nn + 1; jq <= i__3; ++jq) {
		    wk[(i__4 = jq - 1) < 1 * wk_dim1 && 0 <= i__4 ? i__4 : 
			    s_rnge("wk", i__4, "usrint_", (ftnlen)4082)] = 
			    exp(-kk[jq + *nlyr * kk_dim1 - kk_offset] * 
			    dtaucp[*nlyr - 1]);
/* L100: */
		}
		bnddfu = 0.;
		for (iq = *nn; iq >= 1; --iq) {
		    dfuint = 0.;
		    i__3 = *nn;
		    for (jq = 1; jq <= i__3; ++jq) {
			dfuint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1 - 
				gc_offset] * ll[jq + *nlyr * ll_dim1 - 
				ll_offset];
/* L110: */
		    }
		    i__3 = *nstr;
		    for (jq = *nn + 1; jq <= i__3; ++jq) {
			dfuint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1 - 
				gc_offset] * ll[jq + *nlyr * ll_dim1 - 
				ll_offset] * wk[(i__4 = jq - 1) < 1 * wk_dim1 
				&& 0 <= i__4 ? i__4 : s_rnge("wk", i__4, 
				"usrint_", (ftnlen)4095)];
/* L120: */
		    }
		    if (*fbeam > 0.) {
			dfuint += zz[iq + *nlyr * zz_dim1 - zz_offset] * 
				expbea[*nlyr];
		    }
		    dfuint += *delm0 * (zplk0[iq + *nlyr * zplk0_dim1 - 
			    zplk0_offset] + zplk1[iq + *nlyr * zplk1_dim1 - 
			    zplk1_offset] * taucpr[*nlyr]);
		    bnddfu += (*delm0 + 1.) * rmu[iu + (*nn + 1 - iq) * 
			    rmu_dim1 - rmu_offset] * cmu[(i__3 = *nn + 1 - iq 
			    - 1) < 1 * cmu_dim1 && 0 <= i__3 ? i__3 : s_rnge(
			    "cmu", i__3, "usrint_", (ftnlen)4104)] * cwt[(
			    i__4 = *nn + 1 - iq - 1) < 1 * cwt_dim1 && 0 <= 
			    i__4 ? i__4 : s_rnge("cwt", i__4, "usrint_", (
			    ftnlen)4104)] * dfuint;
/* L130: */
		}
		bnddir = 0.;
		if (*fbeam > 0.) {
		    bnddir = *umu0 * *fbeam / *pi * rmu[iu - rmu_offset] * 
			    expbea[*nlyr];
		}
		bndint = (bnddfu + bnddir + *delm0 * emu[(i__4 = iu - 1) < 1 *
			 emu_dim1 && 0 <= i__4 ? i__4 : s_rnge("emu", i__4, 
			"usrint_", (ftnlen)4112)] * *bplank) * exp((utaupr[(
			i__3 = lu - 1) < 1 * utaupr_dim1 && 0 <= i__3 ? i__3 :
			 s_rnge("utaupr", i__3, "usrint_", (ftnlen)4112)] - 
			taucpr[*nlyr]) / umu[iu - 1]);
	    }
L140:
	    uum[(i__3 = iu + lu * uum_dim1 - uum_offset) < 1 * uum_dim1 * 
		    uum_dim2 && 0 <= i__3 ? i__3 : s_rnge("uum", i__3, "usri\
nt_", (ftnlen)4118)] = palint + plkint + bndint;
L150:
	    ;
	}
/* L160: */
    }
    return 0;
} /* usrint_ */

/* ****************************************************************** */
/* ********** DISORT service routines ******************************* */
/* ****************************************************************** */
/* Subroutine */ int chekin_(integer *nlyr, doublereal *dtauc, doublereal *
	ssalb, doublereal *pmom, doublereal *temper, doublereal *wvnmlo, 
	doublereal *wvnmhi, logical *usrtau, integer *ntau, doublereal *utau, 
	integer *nstr, logical *usrang, integer *numu, doublereal *umu, 
	integer *nphi, doublereal *phi, integer *ibcnd, doublereal *fbeam, 
	doublereal *umu0, doublereal *phi0, doublereal *fisot, logical *
	lamber, doublereal *albedo, doublereal *hl, doublereal *btemp, 
	doublereal *ttemp, doublereal *temis, logical *plank, logical *onlyfl,
	 doublereal *accur, doublereal *tauc, integer *maxcly, integer *
	maxulv, integer *maxumu, integer *maxcmu, integer *maxphi, integer *
	mxcly, integer *mxulv, integer *mxumu, integer *mxcmu, integer *mxphi,
	 integer *mxsqt)
{
    /* System generated locals */
    integer dtauc_dim1, phi_dim1, pmom_dim1, pmom_dim2, pmom_offset, 
	    ssalb_dim1, tauc_dim1, temper_dim1, umu_dim1, utau_dim1, i__1, 
	    i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    extern doublereal dref_(doublereal *, doublereal *, integer *);
    static integer irmu, j, k, lc, iu, lu;
    static doublereal flxalb;
    extern logical wrtbad_(char *, ftnlen);
    static logical inperr;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    extern logical wrtdim_(char *, integer *, ftnlen);
    static integer numsqt;
    static doublereal rmu;

/*           Checks the input dimensions and variables */
/*   Calls- WRTBAD, WRTDIM, DREF, ERRMSG */
/*   Called by- DISORT */
/* -------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    temper_dim1 = *maxcly - 0 + 1;
    ssalb_dim1 = *maxcly;
    dtauc_dim1 = *maxcly;
    utau_dim1 = *maxulv;
    umu_dim1 = *maxumu;
    pmom_dim1 = *maxcmu - 0 + 1;
    pmom_dim2 = *maxcly;
    pmom_offset = 0 + pmom_dim1 * 1;
    phi_dim1 = *maxphi;
    tauc_dim1 = *mxcly - 0 + 1;

    /* Function Body */
    inperr = FALSE_;
    if (*nlyr < 1) {
	inperr = wrtbad_("NLYR", (ftnlen)4);
    }
    if (*nlyr > *maxcly) {
	inperr = wrtbad_("MAXCLY", (ftnlen)6);
    }
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	if (dtauc[(i__2 = lc - 1) < 1 * dtauc_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("dtauc", i__2, "chekin_", (ftnlen)4193)] < 0.) {
	    inperr = wrtbad_("DTAUC", (ftnlen)5);
	}
	if (ssalb[(i__2 = lc - 1) < 1 * ssalb_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("ssalb", i__2, "chekin_", (ftnlen)4195)] < 0. || ssalb[
		(i__3 = lc - 1) < 1 * ssalb_dim1 && 0 <= i__3 ? i__3 : s_rnge(
		"ssalb", i__3, "chekin_", (ftnlen)4195)] > 1.) {
	    inperr = wrtbad_("SSALB", (ftnlen)5);
	}
	if (*plank && *ibcnd != 1) {
	    if (lc == 1 && temper[(i__2 = 0) < 1 * temper_dim1 ? i__2 : 
		    s_rnge("temper", i__2, "chekin_", (ftnlen)4200)] < 0.) {
		inperr = wrtbad_("TEMPER", (ftnlen)6);
	    }
	    if (temper[(i__2 = lc) < 1 * temper_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("temper", i__2, "chekin_", (ftnlen)4203)] < 0.) {
		inperr = wrtbad_("TEMPER", (ftnlen)6);
	    }
	}
	i__2 = *nstr;
	for (k = 0; k <= i__2; ++k) {
	    if (pmom[(i__3 = k + lc * pmom_dim1 - pmom_offset) < 1 * 
		    pmom_dim1 * pmom_dim2 && 0 <= i__3 ? i__3 : s_rnge("pmom",
		     i__3, "chekin_", (ftnlen)4209)] < -1. || pmom[(i__4 = k 
		    + lc * pmom_dim1 - pmom_offset) < 1 * pmom_dim1 * 
		    pmom_dim2 && 0 <= i__4 ? i__4 : s_rnge("pmom", i__4, 
		    "chekin_", (ftnlen)4209)] > 1.) {
		inperr = wrtbad_("PMOM", (ftnlen)4);
	    }
/* L10: */
	}
/* L20: */
    }
    if (*ibcnd == 1) {
	if (*maxulv < 2) {
	    inperr = wrtbad_("MAXULV", (ftnlen)6);
	}
    } else if (*usrtau) {
	if (*ntau < 1) {
	    inperr = wrtbad_("NTAU", (ftnlen)4);
	}
	if (*maxulv < *ntau) {
	    inperr = wrtbad_("MAXULV", (ftnlen)6);
	}
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    if ((d__1 = utau[(i__2 = lu - 1) < 1 * utau_dim1 && 0 <= i__2 ? 
		    i__2 : s_rnge("utau", i__2, "chekin_", (ftnlen)4229)] - 
		    tauc[(i__3 = *nlyr) < 1 * tauc_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("tauc", i__3, "chekin_", (ftnlen)4229)], abs(d__1))
		     <= 1e-4) {
		utau[(i__4 = lu - 1) < 1 * utau_dim1 && 0 <= i__4 ? i__4 : 
			s_rnge("utau", i__4, "chekin_", (ftnlen)4229)] = tauc[
			(i__5 = *nlyr) < 1 * tauc_dim1 && 0 <= i__5 ? i__5 : 
			s_rnge("tauc", i__5, "chekin_", (ftnlen)4229)];
	    }
	    if (utau[(i__2 = lu - 1) < 1 * utau_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("utau", i__2, "chekin_", (ftnlen)4232)] < 0. || 
		    utau[(i__3 = lu - 1) < 1 * utau_dim1 && 0 <= i__3 ? i__3 :
		     s_rnge("utau", i__3, "chekin_", (ftnlen)4232)] > tauc[(
		    i__4 = *nlyr) < 1 * tauc_dim1 && 0 <= i__4 ? i__4 : 
		    s_rnge("tauc", i__4, "chekin_", (ftnlen)4232)]) {
		inperr = wrtbad_("UTAU", (ftnlen)4);
	    }
/* L30: */
	}
    } else {
	if (*maxulv < *nlyr + 1) {
	    inperr = wrtbad_("MAXULV", (ftnlen)6);
	}
    }
    if (*nstr < 2 || *nstr % 2 != 0) {
	inperr = wrtbad_("NSTR", (ftnlen)4);
    }
    if (*nstr == 2) {
	errmsg_("CHEKIN--2 streams not recommended; use specialized 2-stream\
 code instead", &c_false, (ftnlen)72);
    }
    if (*nstr > *maxcmu) {
	inperr = wrtbad_("MAXCMU", (ftnlen)6);
    }
    if (*usrang) {
	if (*numu < 0) {
	    inperr = wrtbad_("NUMU", (ftnlen)4);
	}
	if (! (*onlyfl) && *numu == 0) {
	    inperr = wrtbad_("NUMU", (ftnlen)4);
	}
	if (*numu > *maxumu) {
	    inperr = wrtbad_("MAXUMU", (ftnlen)6);
	}
	if (*ibcnd == 1 && *numu << 1 > *maxumu) {
	    inperr = wrtbad_("MAXUMU", (ftnlen)6);
	}
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    if (umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("umu", i__2, "chekin_", (ftnlen)4265)] < -1. || 
		    umu[(i__3 = iu - 1) < 1 * umu_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("umu", i__3, "chekin_", (ftnlen)4265)] > 1. || umu[
		    (i__4 = iu - 1) < 1 * umu_dim1 && 0 <= i__4 ? i__4 : 
		    s_rnge("umu", i__4, "chekin_", (ftnlen)4265)] == 0.) {
		inperr = wrtbad_("UMU", (ftnlen)3);
	    }
	    if (*ibcnd == 1 && umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= 
		    i__2 ? i__2 : s_rnge("umu", i__2, "chekin_", (ftnlen)4268)
		    ] < 0.) {
		inperr = wrtbad_("UMU", (ftnlen)3);
	    }
	    if (iu > 1) {
		if (umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : 
			s_rnge("umu", i__2, "chekin_", (ftnlen)4273)] < umu[(
			i__3 = iu - 2) < 1 * umu_dim1 && 0 <= i__3 ? i__3 : 
			s_rnge("umu", i__3, "chekin_", (ftnlen)4273)]) {
		    inperr = wrtbad_("UMU", (ftnlen)3);
		}
	    }
/* L40: */
	}
    } else {
	if (*maxumu < *nstr) {
	    inperr = wrtbad_("MAXUMU", (ftnlen)6);
	}
    }
    if (! (*onlyfl) && *ibcnd != 1) {
	if (*nphi <= 0) {
	    inperr = wrtbad_("NPHI", (ftnlen)4);
	}
	if (*nphi > *maxphi) {
	    inperr = wrtbad_("MAXPHI", (ftnlen)6);
	}
	i__1 = *nphi;
	for (j = 1; j <= i__1; ++j) {
	    if (phi[(i__2 = j - 1) < 1 * phi_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("phi", i__2, "chekin_", (ftnlen)4294)] < 0. || phi[
		    (i__3 = j - 1) < 1 * phi_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("phi", i__3, "chekin_", (ftnlen)4294)] > 360.) {
		inperr = wrtbad_("PHI", (ftnlen)3);
	    }
/* L50: */
	}
    }
    if (*ibcnd < 0 || *ibcnd > 1) {
	inperr = wrtbad_("IBCND", (ftnlen)5);
    }
    if (*ibcnd == 0) {
	if (*fbeam < 0.) {
	    inperr = wrtbad_("FBEAM", (ftnlen)5);
	}
	if (*fbeam > 0. && (*umu0 <= 0. || *umu0 > 1.)) {
	    inperr = wrtbad_("UMU0", (ftnlen)4);
	}
	if (*fbeam > 0. && (*phi0 < 0. || *phi0 > 360.)) {
	    inperr = wrtbad_("PHI0", (ftnlen)4);
	}
	if (*fisot < 0.) {
	    inperr = wrtbad_("FISOT", (ftnlen)5);
	}
	if (*lamber) {
	    if (*albedo < 0. || *albedo > 1.) {
		inperr = wrtbad_("ALBEDO", (ftnlen)6);
	    }
	} else {
/*                    ** Make sure flux albedo at dense mesh of incident */
/*                    ** angles does not assume unphysical values */
	    for (irmu = 0; irmu <= 100; ++irmu) {
		rmu = irmu * .01;
		flxalb = dref_(&rmu, hl, nstr);
		if (flxalb < 0. || flxalb > 1.) {
		    inperr = wrtbad_("HL", (ftnlen)2);
		}
/* L60: */
	    }
	}
    } else if (*ibcnd == 1) {
	if (*albedo < 0. || *albedo > 1.) {
	    inperr = wrtbad_("ALBEDO", (ftnlen)6);
	}
    }
    if (*plank && *ibcnd != 1) {
	if (*wvnmlo < 0. || *wvnmhi < *wvnmlo) {
	    inperr = wrtbad_("WVNMLO,HI", (ftnlen)9);
	}
	if (*temis < 0. || *temis > 1.) {
	    inperr = wrtbad_("TEMIS", (ftnlen)5);
	}
	if (*btemp < 0.) {
	    inperr = wrtbad_("BTEMP", (ftnlen)5);
	}
	if (*ttemp < 0.) {
	    inperr = wrtbad_("TTEMP", (ftnlen)5);
	}
    }
    if (*accur < 0. || *accur > .01) {
	inperr = wrtbad_("ACCUR", (ftnlen)5);
    }
    if (*mxcly < *nlyr) {
	inperr = wrtdim_("MXCLY", nlyr, (ftnlen)5);
    }
    if (*ibcnd != 1) {
	if (*usrtau && *mxulv < *ntau) {
	    inperr = wrtdim_("MXULV", ntau, (ftnlen)5);
	}
	if (! (*usrtau) && *mxulv < *nlyr + 1) {
	    i__1 = *nlyr + 1;
	    inperr = wrtdim_("MXULV", &i__1, (ftnlen)5);
	}
    } else {
	if (*mxulv < 2) {
	    inperr = wrtdim_("MXULV", &c__2, (ftnlen)5);
	}
    }
    if (*mxcmu < *nstr) {
	inperr = wrtdim_("MXCMU", nstr, (ftnlen)5);
    }
    if (*usrang && *mxumu < *numu) {
	inperr = wrtdim_("MXUMU", numu, (ftnlen)5);
    }
    if (*usrang && *ibcnd == 1 && *mxumu < *numu << 1) {
	i__1 = *numu << 1;
	inperr = wrtdim_("MXUMU", &i__1, (ftnlen)5);
    }
    if (! (*usrang) && *mxumu < *nstr) {
	inperr = wrtdim_("MXUMU", nstr, (ftnlen)5);
    }
    if (! (*onlyfl) && *ibcnd != 1 && *mxphi < *nphi) {
	inperr = wrtdim_("MXPHI", nphi, (ftnlen)5);
    }
    numsqt = max(100,*nstr) << 1;
    if (*mxsqt < numsqt) {
	inperr = wrtdim_("MXSQT", &numsqt, (ftnlen)5);
    }
    if (inperr) {
	errmsg_("DISORT--input and/or dimension errors", &c_true, (ftnlen)37);
    }
    if (*plank) {
	i__1 = *nlyr;
	for (lc = 1; lc <= i__1; ++lc) {
	    if ((d__1 = temper[(i__2 = lc) < 1 * temper_dim1 && 0 <= i__2 ? 
		    i__2 : s_rnge("temper", i__2, "chekin_", (ftnlen)4400)] - 
		    temper[(i__3 = lc - 1) < 1 * temper_dim1 && 0 <= i__3 ? 
		    i__3 : s_rnge("temper", i__3, "chekin_", (ftnlen)4400)], 
		    abs(d__1)) > 20.) {
		errmsg_("CHEKIN--vertical temperature step may be too large \
for good accuracy", &c_false, (ftnlen)68);
	    }
/* L70: */
	}
    }
    return 0;
} /* chekin_ */

doublereal dref_(doublereal *mu, doublereal *hl, integer *nstr)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    integer hl_dim1, i__1, i__2, i__3;
    doublereal ret_val;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal c__[100];
    static integer l;
    static doublereal cl, pl;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static doublereal plm1, plm2;

/*        Exact flux albedo for given angle of incidence, given */
/*        a bidirectional reflectivity characterized by its */
/*        Legendre coefficients ( NOTE** these will only agree */
/*        with bottom-boundary albedos calculated by DISORT in */
/*        the limit as number of streams go to infinity, because */
/*        DISORT evaluates the integral 'CL' only approximately, */
/*        by quadrature, while this routine calculates it exactly.) */

/*  INPUT :   MU     Cosine of incidence angle */

/*            HL     Legendre coefficients of bidirectional reflectivity */

/*          NSTR     Number of elements of HL to consider */


/*  INTERNAL VARIABLES (P-sub-L is the L-th Legendre polynomial) : */

/*       CL      Integral from 0 to 1 of  MU * P-sub-L(MU) */
/*                   (vanishes for  L = 3, 5, 7, ... ) */
/*       PL      P-sub-L */
/*       PLM1    P-sub-(L-1) */
/*       PLM2    P-sub-(L-2) */

/*   Called by- CHEKIN */
/*   Calls- ERRMSG */
/* +-------------------------------------------------------------------+ */
/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    hl_dim1 = *nstr - 0 + 1;

    /* Function Body */
/*     .. */
    if (pass1) {
	pass1 = FALSE_;
	cl = .125;
	c__[1] = cl * 10.;
	for (l = 4; l <= 100; l += 2) {
	    cl = -cl * (l - 3) / (l + 2);
	    c__[(i__1 = l - 1) < 100 && 0 <= i__1 ? i__1 : s_rnge("c", i__1, 
		    "dref_", (ftnlen)4486)] = ((l << 1) + 1) * 2. * cl;
/* L10: */
	}
    }
    if (*nstr < 2 || abs(*mu) > 1.) {
	errmsg_("DREF--input argument error(s)", &c_true, (ftnlen)29);
    }
    if (*nstr > 100) {
	errmsg_("DREF--parameter MAXTRM too small", &c_true, (ftnlen)32);
    }
    ret_val = hl[(i__1 = 0) < 1 * hl_dim1 ? i__1 : s_rnge("hl", i__1, "dref_",
	     (ftnlen)4499)] - hl[(i__2 = 1) < 1 * hl_dim1 ? i__2 : s_rnge(
	    "hl", i__2, "dref_", (ftnlen)4499)] * 2. * *mu;
    plm2 = 1.;
    plm1 = -(*mu);
    i__1 = *nstr - 1;
    for (l = 2; l <= i__1; ++l) {
/*                                ** Legendre polynomial recurrence */
	pl = (((l << 1) - 1) * (-(*mu)) * plm1 - (l - 1) * plm2) / l;
	if (l % 2 == 0) {
	    ret_val += c__[(i__2 = l - 1) < 100 && 0 <= i__2 ? i__2 : s_rnge(
		    "c", i__2, "dref_", (ftnlen)4508)] * hl[(i__3 = l) < 1 * 
		    hl_dim1 && 0 <= i__3 ? i__3 : s_rnge("hl", i__3, "dref_", 
		    (ftnlen)4508)] * pl;
	}
	plm2 = plm1;
	plm1 = pl;
/* L20: */
    }
    if (ret_val < 0. || ret_val > 1.) {
	errmsg_("DREF--albedo value not in (0,1)", &c_false, (ftnlen)31);
    }
    return ret_val;
} /* dref_ */

/* Subroutine */ int lepoly_(integer *nmu, integer *m, integer *maxmu, 
	integer *twonm1, doublereal *mu, doublereal *sqt, doublereal *ylm)
{
    /* System generated locals */
    integer ylm_dim1, ylm_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, l;
    static doublereal tmp1, tmp2;

/*       Computes the normalized associated Legendre polynomial, */
/*       defined in terms of the associated Legendre polynomial */
/*       Plm = P-sub-l-super-m as */

/*             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU) */

/*       for fixed order m and all degrees from l = m to TWONM1. */
/*       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available */
/*       from a prior call to the routine. */

/*       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of */
/*                  High-Order Associated Legendre Polynomials, */
/*                  J. Quant. Spectrosc. Radiat. Transfer 10, */
/*                  557-562, 1970.  (hereafter D/A) */

/*       METHOD: Varying degree recurrence relationship. */

/*       NOTES: */
/*       (1) The D/A formulas are transformed by setting M=n-1; L=k-1. */
/*       (2) Assumes that routine is called first with  M = 0, then with */
/*           M = 1, etc. up to  M = TWONM1. */
/*       (3) Loops are written in such a way as to vectorize. */


/*  I N P U T     V A R I A B L E S: */

/*       NMU    :  Number of arguments of YLM */

/*       M      :  Order of YLM */

/*       MAXMU  :  First dimension of YLM */

/*       TWONM1 :  Max degree of YLM */

/*       MU(i)  :  Arguments of YLM (i = 1 to NMU) */

/*       SQT(k) :  Square root of k */

/*       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist */
/*       from a prior call. */


/*  O U T P U T     V A R I A B L E: */

/*       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre */
/*                   polynomials evaluated at argument MU(i) */

/*   Called by- DISORT, ALBTRN, SURFAC */
/*   Calls- ERRMSG */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
    /* Parameter adjustments */
    ylm_dim1 = *maxmu - 0 + 1;
    ylm_offset = 0 + ylm_dim1 * 1;

    /* Function Body */
    if (*m == 0) {
/*                             ** Upward recurrence for ordinary */
/*                             ** Legendre polynomials */
	i__1 = *nmu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ylm[i__ * ylm_dim1 - ylm_offset] = 1.;
	    ylm[i__ * ylm_dim1 + 1 - ylm_offset] = mu[i__ - 1];
/* L20: */
	}
	i__1 = *twonm1;
	for (l = 2; l <= i__1; ++l) {
	    i__2 = *nmu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ylm[l + i__ * ylm_dim1 - ylm_offset] = (((l << 1) - 1) * mu[
			i__ - 1] * ylm[l - 1 + i__ * ylm_dim1 - ylm_offset] - 
			(l - 1) * ylm[l - 2 + i__ * ylm_dim1 - ylm_offset]) / 
			l;
/* L30: */
	    }
/* L40: */
	}
    } else {
	i__1 = *nmu;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*                               ** Y-sub-m-super-m; derived from */
/*                               ** D/A Eqs. (11,12) */
/* Computing 2nd power */
	    d__1 = mu[i__ - 1];
	    ylm[*m + i__ * ylm_dim1 - ylm_offset] = -sqt[(*m << 1) - 2] / sqt[
		    (*m << 1) - 1] * sqrt(1. - d__1 * d__1) * ylm[*m - 1 + 
		    i__ * ylm_dim1 - ylm_offset];
/*                              ** Y-sub-(m+1)-super-m; derived from */
/*                              ** D/A Eqs.(13,14) using Eqs.(11,12) */
	    ylm[*m + 1 + i__ * ylm_dim1 - ylm_offset] = sqt[*m * 2] * mu[i__ 
		    - 1] * ylm[*m + i__ * ylm_dim1 - ylm_offset];
/* L50: */
	}
/*                                   ** Upward recurrence; D/A EQ.(10) */
	i__1 = *twonm1;
	for (l = *m + 2; l <= i__1; ++l) {
	    tmp1 = sqt[l - *m - 1] * sqt[l + *m - 1];
	    tmp2 = sqt[l - *m - 2] * sqt[l + *m - 2];
	    i__2 = *nmu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ylm[l + i__ * ylm_dim1 - ylm_offset] = (((l << 1) - 1) * mu[
			i__ - 1] * ylm[l - 1 + i__ * ylm_dim1 - ylm_offset] - 
			tmp2 * ylm[l - 2 + i__ * ylm_dim1 - ylm_offset]) / 
			tmp1;
/* L60: */
	    }
/* L70: */
	}
    }
    return 0;
} /* lepoly_ */

doublereal plkavg_(doublereal *wnumlo, doublereal *wnumhi, doublereal *t)
{
    /* Initialized data */

    static doublereal c2 = 1.438786;
    static doublereal sigma = 5.67032e-8;
    static doublereal vcut = 1.5;
    static doublereal vcp[7] = { 10.25,5.7,3.9,2.9,2.3,1.9,0. };
    static doublereal pi = 0.;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double asin(doublereal), log(doublereal), exp(doublereal);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal conc;
    static integer mmax;
    static doublereal vmax, plkconst, a, b, d__[2];
    static integer i__, k, m, n;
    static doublereal p[2], v[2], epsil, boltz;
    extern doublereal r1mach_(integer *);
    static doublereal hh, ex, mv, sigdpi, oldval;
    static integer smallv;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static doublereal speedlight, del, val, exm, vsq, val0;

/*        Computes Planck function integrated between two wavenumbers */

/*  INPUT :  WNUMLO : Lower wavenumber (inv cm) of spectral interval */

/*           WNUMHI : Upper wavenumber */

/*           T      : Temperature (K) */

/*  OUTPUT : PLKAVG : Integrated Planck function ( Watts/sq m ) */
/*                      = Integral (WNUMLO to WNUMHI) of */
/*                        2h c**2  nu**3 / ( EXP(hc nu/kT) - 1) */
/*                        (where h=Plancks constant, c=speed of */
/*                         light, nu=wavenumber, T=temperature, */
/*                         and k = Boltzmann constant) */

/*  Reference : Specifications of the Physical World: New Value */
/*                 of the Fundamental Constants, Dimensions/N.B.S., */
/*                 Jan. 1974 */

/*  Method :  For WNUMLO close to WNUMHI, a Simpson-rule quadrature */
/*            is done to avoid ill-conditioning; otherwise */

/*            (1)  For WNUMLO or WNUMHI small, */
/*                 integral(0 to WNUMLO/HI) is calculated by expanding */
/*                 the integrand in a power series and integrating */
/*                 term by term; */

/*            (2)  Otherwise, integral(WNUMLO/HI to INFINITY) is */
/*                 calculated by expanding the denominator of the */
/*                 integrand in powers of the exponential and */
/*                 integrating term by term. */

/*  Accuracy :  At least 6 significant digits, assuming the */
/*              physical constants are infinitely accurate */

/*  ERRORS WHICH ARE NOT TRAPPED: */

/*      * power or exponential series may underflow, giving no */
/*        significant digits.  This may or may not be of concern, */
/*        depending on the application. */

/*      * Simpson-rule special case is skipped when denominator of */
/*        integrand will cause overflow.  In that case the normal */
/*        procedure is used, which may be inaccurate if the */
/*        wavenumber limits (WNUMLO, WNUMHI) are close together. */

/*  LOCAL VARIABLES */

/*        A1,2,... :  Power series coefficients */
/*        C2       :  h * c / k, in units cm*K (h = Plancks constant, */
/*                      c = speed of light, k = Boltzmann constant) */
/*        D(I)     :  Exponential series expansion of integral of */
/*                       Planck function from WNUMLO (i=1) or WNUMHI */
/*                       (i=2) to infinity */
/*        EPSIL    :  SMALLEST NUMBER SUCH THAT 1+EPSIL .GT. 1 on */
/*                       computer */
/*        EX       :  EXP( - V(I) ) */
/*        EXM      :  EX**M */
/*        MMAX     :  No. of terms to take in exponential series */
/*        MV       :  Multiples of V(I) */
/*        P(I)     :  Power series expansion of integral of */
/*                       Planck function from zero to WNUMLO (I=1) or */
/*                       WNUMHI (I=2) */
/*        PI       :  3.14159... */
/*        SIGMA    :  Stefan-Boltzmann constant (W/m**2/K**4) */
/*        SIGDPI   :  SIGMA / PI */
/*        SMALLV   :  Number of times the power series is used (0,1,2) */
/*        V(I)     :  C2 * (WNUMLO(I=1) or WNUMHI(I=2)) / temperature */
/*        VCUT     :  Power-series cutoff point */
/*        VCP      :  Exponential series cutoff points */
/*        VMAX     :  Largest allowable argument of EXP function */

/*   Called by- DISORT */
/*   Calls- R1MACH, ERRMSG */
/* ---------------------------------------------------------------------- */
/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
    if (pi == 0.) {
	pi = asin(1.) * 2.;
	vmax = log(r1mach_(&c__2));
	epsil = r1mach_(&c__4);
	sigdpi = sigma / pi;
/* Computing 4th power */
	d__1 = pi, d__1 *= d__1;
	conc = 15. / (d__1 * d__1);
    }
    if (*t < 0. || *wnumhi < *wnumlo || *wnumlo < 0.) {
	errmsg_("PLKAVG--temperature or wavenums. wrong", &c_true, (ftnlen)38)
		;
    }
    if (*t < 1e-4) {
	ret_val = 0.;
	return ret_val;
    }
    v[0] = c2 * *wnumlo / *t;
    v[1] = c2 * *wnumhi / *t;
/*     (Added by CE) Included the monochromatic Planck function needed for */
/*     monochromatic arts calculations, this is especially included fo use in */
/*     arts, so give wave number in si units */
    if (*wnumhi == *wnumlo) {
	plkconst = 6.6262e-34;
	speedlight = 299792458.;
	boltz = 1.380662e-23;
	a = plkconst * 2 * speedlight;
	b = plkconst * speedlight / boltz;
	ret_val = a * *wnumlo * *wnumlo * *wnumlo / (exp(b * *wnumlo / *t) - 
		1);
	return ret_val;
    }
    if (v[0] > epsil && v[1] < vmax && (*wnumhi - *wnumlo) / *wnumhi < .01) {
/*                          ** Wavenumbers are very close.  Get integral */
/*                          ** by iterating Simpson rule to convergence. */
	hh = v[1] - v[0];
	oldval = 0.;
/* Computing 3rd power */
	d__1 = v[0];
/* Computing 3rd power */
	d__2 = v[1];
	val0 = d__1 * (d__1 * d__1) / (exp(v[0]) - 1) + d__2 * (d__2 * d__2) /
		 (exp(v[1]) - 1);
	for (n = 1; n <= 10; ++n) {
	    del = hh / (n << 1);
	    val = val0;
	    i__1 = (n << 1) - 1;
	    for (k = 1; k <= i__1; ++k) {
		d__1 = v[0] + k * del;
/* Computing 3rd power */
		d__2 = d__1;
		val += (k % 2 + 1 << 1) * (d__2 * (d__2 * d__2) / (exp(d__1) 
			- 1));
/* L10: */
	    }
	    val = del / 3. * val;
	    if ((d__1 = (val - oldval) / val, abs(d__1)) <= 1e-6) {
		goto L30;
	    }
	    oldval = val;
/* L20: */
	}
	errmsg_("PLKAVG--Simpson rule didnt converge", &c_false, (ftnlen)35);
L30:
/* Computing 4th power */
	d__1 = *t, d__1 *= d__1;
	ret_val = sigdpi * (d__1 * d__1) * conc * val;
	return ret_val;
    }
/*                          *** General case *** */
    smallv = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	if (v[(i__1 = i__ - 1) < 2 && 0 <= i__1 ? i__1 : s_rnge("v", i__1, 
		"plkavg_", (ftnlen)4853)] < vcut) {
/*                                   ** Use power series */
	    ++smallv;
/* Computing 2nd power */
	    d__1 = v[(i__1 = i__ - 1) < 2 && 0 <= i__1 ? i__1 : s_rnge("v", 
		    i__1, "plkavg_", (ftnlen)4856)];
	    vsq = d__1 * d__1;
	    p[(i__1 = i__ - 1) < 2 && 0 <= i__1 ? i__1 : s_rnge("p", i__1, 
		    "plkavg_", (ftnlen)4857)] = conc * vsq * v[(i__2 = i__ - 
		    1) < 2 && 0 <= i__2 ? i__2 : s_rnge("v", i__2, "plkavg_", 
		    (ftnlen)4857)] * (v[(i__3 = i__ - 1) < 2 && 0 <= i__3 ? 
		    i__3 : s_rnge("v", i__3, "plkavg_", (ftnlen)4857)] * (v[(
		    i__4 = i__ - 1) < 2 && 0 <= i__4 ? i__4 : s_rnge("v", 
		    i__4, "plkavg_", (ftnlen)4857)] * (vsq * (vsq * (vsq * 
		    -7.5156325156325161e-8 + 3.6743092298647855e-6) - 
		    1.9841269841269841e-4) + .016666666666666666) - .125) + 
		    .33333333333333331);
	} else {
/*                      ** Use exponential series */
	    mmax = 0;
/*                                ** Find upper limit of series */
L40:
	    ++mmax;
	    if (v[(i__1 = i__ - 1) < 2 && 0 <= i__1 ? i__1 : s_rnge("v", i__1,
		     "plkavg_", (ftnlen)4868)] < vcp[(i__2 = mmax - 1) < 7 && 
		    0 <= i__2 ? i__2 : s_rnge("vcp", i__2, "plkavg_", (ftnlen)
		    4868)]) {
		goto L40;
	    }
	    ex = exp(-v[(i__1 = i__ - 1) < 2 && 0 <= i__1 ? i__1 : s_rnge(
		    "v", i__1, "plkavg_", (ftnlen)4870)]);
	    exm = 1.;
	    d__[(i__1 = i__ - 1) < 2 && 0 <= i__1 ? i__1 : s_rnge("d", i__1, 
		    "plkavg_", (ftnlen)4872)] = 0.;
	    i__1 = mmax;
	    for (m = 1; m <= i__1; ++m) {
		mv = m * v[(i__2 = i__ - 1) < 2 && 0 <= i__2 ? i__2 : s_rnge(
			"v", i__2, "plkavg_", (ftnlen)4875)];
		exm = ex * exm;
/* Computing 4th power */
		i__4 = m, i__4 *= i__4;
		d__[(i__2 = i__ - 1) < 2 && 0 <= i__2 ? i__2 : s_rnge("d", 
			i__2, "plkavg_", (ftnlen)4877)] = d__[(i__3 = i__ - 1)
			 < 2 && 0 <= i__3 ? i__3 : s_rnge("d", i__3, "plkavg_"
			, (ftnlen)4877)] + exm * (mv * (mv * (mv + 3.) + 6.) 
			+ 6.) / (i__4 * i__4);
/* L50: */
	    }
	    d__[(i__1 = i__ - 1) < 2 && 0 <= i__1 ? i__1 : s_rnge("d", i__1, 
		    "plkavg_", (ftnlen)4881)] = conc * d__[(i__2 = i__ - 1) < 
		    2 && 0 <= i__2 ? i__2 : s_rnge("d", i__2, "plkavg_", (
		    ftnlen)4881)];
	}
/* L60: */
    }
/*                              ** Handle ill-conditioning */
    if (smallv == 2) {
/*                                    ** WNUMLO and WNUMHI both small */
	ret_val = p[1] - p[0];
    } else if (smallv == 1) {
/*                                    ** WNUMLO small, WNUMHI large */
	ret_val = 1. - p[0] - d__[1];
    } else {
/*                                    ** WNUMLO and WNUMHI both large */
	ret_val = d__[0] - d__[1];
    }
/* Computing 4th power */
    d__1 = *t, d__1 *= d__1;
    ret_val = sigdpi * (d__1 * d__1) * ret_val;
    if (ret_val == 0.) {
	errmsg_("PLKAVG--returns zero; possible underflow", &c_false, (ftnlen)
		40);
    }
    return ret_val;
} /* plkavg_ */

/* Subroutine */ int pravin_(doublereal *umu, integer *numu, integer *maxumu, 
	doublereal *utau, integer *ntau, doublereal *u0u)
{
    /* System generated locals */
    integer u0u_dim1, u0u_dim2, u0u_offset, umu_dim1, utau_dim1, i__1, i__2, 
	    i__3, i__4, i__5;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer iumin, iumax, npass, iu, np, lu, lenfmt;

    /* Fortran I/O blocks */
    static cilist io___305 = { 0, 6, 0, "(//,A)", 0 };
    static cilist io___308 = { 0, 6, 0, "(/,A,/,A)", 0 };
    static cilist io___312 = { 0, 6, 0, "(/,10X,8F14.5)", 0 };
    static cilist io___315 = { 0, 6, 0, "(0P,F10.4,1P,8E14.4)", 0 };


/*        Print azimuthally averaged intensities at user angles */
/*   Called by- DISORT */
/*     LENFMT   Max number of polar angle cosines UMU that can be */
/*              printed on one line, as set in FORMAT statement */
/* -------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    umu_dim1 = *numu;
    u0u_dim1 = *maxumu;
    u0u_dim2 = *ntau;
    u0u_offset = 1 + u0u_dim1 * 1;
    utau_dim1 = *ntau;

    /* Function Body */
    if (*numu < 1) {
	return 0;
    }
    s_wsfe(&io___305);
    do_fio(&c__1, " *******  AZIMUTHALLY AVERAGED INTENSITIES (at user polar\
 angles)  ********", (ftnlen)75);
    e_wsfe();
    lenfmt = 8;
    npass = (*numu - 1) / lenfmt + 1;
    s_wsfe(&io___308);
    do_fio(&c__1, "   Optical   Polar Angle Cosines", (ftnlen)32);
    do_fio(&c__1, "     Depth", (ftnlen)10);
    e_wsfe();
    i__1 = npass;
    for (np = 1; np <= i__1; ++np) {
	iumin = lenfmt * (np - 1) + 1;
/* Computing MIN */
	i__2 = lenfmt * np;
	iumax = min(i__2,*numu);
	s_wsfe(&io___312);
	i__2 = iumax;
	for (iu = iumin; iu <= i__2; ++iu) {
	    do_fio(&c__1, (char *)&umu[(i__3 = iu - 1) < 1 * umu_dim1 && 0 <= 
		    i__3 ? i__3 : s_rnge("umu", i__3, "pravin_", (ftnlen)4955)
		    ], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	i__3 = *ntau;
	for (lu = 1; lu <= i__3; ++lu) {
	    s_wsfe(&io___315);
	    do_fio(&c__1, (char *)&utau[(i__2 = lu - 1) < 1 * utau_dim1 && 0 
		    <= i__2 ? i__2 : s_rnge("utau", i__2, "pravin_", (ftnlen)
		    4958)], (ftnlen)sizeof(doublereal));
	    i__4 = iumax;
	    for (iu = iumin; iu <= i__4; ++iu) {
		do_fio(&c__1, (char *)&u0u[(i__5 = iu + lu * u0u_dim1 - 
			u0u_offset) < 1 * u0u_dim1 * u0u_dim2 && 0 <= i__5 ? 
			i__5 : s_rnge("u0u", i__5, "pravin_", (ftnlen)4958)], 
			(ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
/* L10: */
	}
/* L20: */
    }
    return 0;
} /* pravin_ */

/* Subroutine */ int prtinp_(integer *nlyr, doublereal *dtauc, doublereal *
	dtaucp, doublereal *ssalb, doublereal *pmom, doublereal *temper, 
	doublereal *wvnmlo, doublereal *wvnmhi, integer *ntau, doublereal *
	utau, integer *nstr, integer *numu, doublereal *umu, integer *nphi, 
	doublereal *phi, integer *ibcnd, doublereal *fbeam, doublereal *umu0, 
	doublereal *phi0, doublereal *fisot, logical *lamber, doublereal *
	albedo, doublereal *hl, doublereal *btemp, doublereal *ttemp, 
	doublereal *temis, logical *deltam, logical *plank, logical *onlyfl, 
	doublereal *accur, doublereal *flyr, logical *lyrcut, doublereal *
	oprim, doublereal *tauc, doublereal *taucpr, integer *maxcmu, logical 
	*prtmom)
{
    /* System generated locals */
    integer hl_dim1, pmom_dim1, pmom_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer j, k, lc, iu, lu;
    static doublereal yessct;

    /* Fortran I/O blocks */
    static cilist io___316 = { 0, 6, 0, "(/,A,I4,A,I4)", 0 };
    static cilist io___317 = { 0, 6, 0, "(I4,A,10F10.4,/,(26X,10F10.4))", 0 };
    static cilist io___319 = { 0, 6, 0, "(I4,A,10F9.5,/,(31X,10F9.5))", 0 };
    static cilist io___321 = { 0, 6, 0, "(I4,A,10F9.2,/,(28X,10F9.2))", 0 };
    static cilist io___323 = { 0, 6, 0, "(A)", 0 };
    static cilist io___324 = { 0, 6, 0, "(A,I2)", 0 };
    static cilist io___325 = { 0, 6, 0, "(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P\
,E11.3)", 0 };
    static cilist io___326 = { 0, 6, 0, "(A,0P,F8.4)", 0 };
    static cilist io___327 = { 0, 6, 0, "(A,/,(10X,10F9.5))", 0 };
    static cilist io___329 = { 0, 6, 0, "(A,2F14.4,/,A,F10.2,A,F10.2,A,F8.4)",
	     0 };
    static cilist io___330 = { 0, 6, 0, "(A)", 0 };
    static cilist io___331 = { 0, 6, 0, "(A,0P,F8.4)", 0 };
    static cilist io___332 = { 0, 6, 0, "(A)", 0 };
    static cilist io___333 = { 0, 6, 0, "(A)", 0 };
    static cilist io___334 = { 0, 6, 0, "(A)", 0 };
    static cilist io___335 = { 0, 6, 0, "(A)", 0 };
    static cilist io___336 = { 0, 6, 0, "(A)", 0 };
    static cilist io___337 = { 0, 6, 0, "(A,1P,E11.2)", 0 };
    static cilist io___338 = { 0, 6, 0, "(A)", 0 };
    static cilist io___339 = { 0, 6, 0, "(/,37X,A,3(/,2A))", 0 };
    static cilist io___340 = { 0, 6, 0, "(/,37X,A,3(/,2A))", 0 };
    static cilist io___343 = { 0, 6, 0, "(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5\
,F9.4,F14.3)", 0 };
    static cilist io___344 = { 0, 6, 0, "(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5\
,F9.4)", 0 };
    static cilist io___345 = { 0, 6, 0, "(85X,F14.3)", 0 };
    static cilist io___346 = { 0, 6, 0, "(/,A)", 0 };
    static cilist io___347 = { 0, 6, 0, "(I6,10F11.6,/,(6X,10F11.6))", 0 };


/*        Print values of input variables */
/*   Called by- DISORT */
/* -------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
    /* Parameter adjustments */
    hl_dim1 = *maxcmu - 0 + 1;
    pmom_dim1 = *maxcmu - 0 + 1;
    pmom_offset = 0 + pmom_dim1 * 1;

    /* Function Body */
    s_wsfe(&io___316);
    do_fio(&c__1, " No. streams =", (ftnlen)14);
    do_fio(&c__1, (char *)&(*nstr), (ftnlen)sizeof(integer));
    do_fio(&c__1, "     No. computational layers =", (ftnlen)31);
    do_fio(&c__1, (char *)&(*nlyr), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ibcnd != 1) {
	s_wsfe(&io___317);
	do_fio(&c__1, (char *)&(*ntau), (ftnlen)sizeof(integer));
	do_fio(&c__1, " User optical depths :", (ftnlen)22);
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    do_fio(&c__1, (char *)&utau[lu - 1], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (! (*onlyfl)) {
	s_wsfe(&io___319);
	do_fio(&c__1, (char *)&(*numu), (ftnlen)sizeof(integer));
	do_fio(&c__1, " User polar angle cosines :", (ftnlen)27);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    do_fio(&c__1, (char *)&umu[iu - 1], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (! (*onlyfl) && *ibcnd != 1) {
	s_wsfe(&io___321);
	do_fio(&c__1, (char *)&(*nphi), (ftnlen)sizeof(integer));
	do_fio(&c__1, " User azimuthal angles :", (ftnlen)24);
	i__1 = *nphi;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&phi[j - 1], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    if (! (*plank) || *ibcnd == 1) {
	s_wsfe(&io___323);
	do_fio(&c__1, " No thermal emission", (ftnlen)20);
	e_wsfe();
    }
    s_wsfe(&io___324);
    do_fio(&c__1, " Boundary condition flag: IBCND =", (ftnlen)33);
    do_fio(&c__1, (char *)&(*ibcnd), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ibcnd == 0) {
	s_wsfe(&io___325);
	do_fio(&c__1, "    Incident beam with intensity =", (ftnlen)34);
	do_fio(&c__1, (char *)&(*fbeam), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, " and polar angle cosine = ", (ftnlen)26);
	do_fio(&c__1, (char *)&(*umu0), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "  and azimuth angle =", (ftnlen)21);
	do_fio(&c__1, (char *)&(*phi0), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "    plus isotropic incident intensity =", (ftnlen)39);
	do_fio(&c__1, (char *)&(*fisot), (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (*lamber) {
	    s_wsfe(&io___326);
	    do_fio(&c__1, "    Bottom albedo (Lambertian) =", (ftnlen)32);
	    do_fio(&c__1, (char *)&(*albedo), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	if (! (*lamber)) {
	    s_wsfe(&io___327);
	    do_fio(&c__1, "    Legendre coeffs of bottom bidirectional refle\
ctivity :", (ftnlen)58);
	    i__1 = *nstr;
	    for (k = 0; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&hl[(i__2 = k) < 1 * hl_dim1 && 0 <= 
			i__2 ? i__2 : s_rnge("hl", i__2, "prtinp_", (ftnlen)
			5031)], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
	}
	if (*plank) {
	    s_wsfe(&io___329);
	    do_fio(&c__1, "    Thermal emission in wavenumber interval :", (
		    ftnlen)45);
	    do_fio(&c__1, (char *)&(*wvnmlo), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*wvnmhi), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, "    Bottom temperature =", (ftnlen)24);
	    do_fio(&c__1, (char *)&(*btemp), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, "    Top temperature =", (ftnlen)21);
	    do_fio(&c__1, (char *)&(*ttemp), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, "    Top emissivity =", (ftnlen)20);
	    do_fio(&c__1, (char *)&(*temis), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else if (*ibcnd == 1) {
	s_wsfe(&io___330);
	do_fio(&c__1, "    Isotropic illumination from top and bottom", (
		ftnlen)46);
	e_wsfe();
	s_wsfe(&io___331);
	do_fio(&c__1, "    Bottom albedo (Lambertian) =", (ftnlen)32);
	do_fio(&c__1, (char *)&(*albedo), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*deltam) {
	s_wsfe(&io___332);
	do_fio(&c__1, " Uses delta-M method", (ftnlen)20);
	e_wsfe();
    }
    if (! (*deltam)) {
	s_wsfe(&io___333);
	do_fio(&c__1, " Does not use delta-M method", (ftnlen)28);
	e_wsfe();
    }
    if (*ibcnd == 1) {
	s_wsfe(&io___334);
	do_fio(&c__1, " Calculate albedo and transmissivity of medium vs. in\
cident beam angle", (ftnlen)70);
	e_wsfe();
    } else if (*onlyfl) {
	s_wsfe(&io___335);
	do_fio(&c__1, " Calculate fluxes and azim-averaged intensities only", 
		(ftnlen)52);
	e_wsfe();
    } else {
	s_wsfe(&io___336);
	do_fio(&c__1, " Calculate fluxes and intensities", (ftnlen)33);
	e_wsfe();
    }
    s_wsfe(&io___337);
    do_fio(&c__1, " Relative convergence criterion for azimuth series =", (
	    ftnlen)52);
    do_fio(&c__1, (char *)&(*accur), (ftnlen)sizeof(doublereal));
    e_wsfe();
    if (*lyrcut) {
	s_wsfe(&io___338);
	do_fio(&c__1, " Sets radiation = 0 below absorption optical depth 10",
		 (ftnlen)53);
	e_wsfe();
    }
/*                                    ** Print layer variables */
/*                                    ** (to read, skip every other line) */
    if (*plank) {
	s_wsfe(&io___339);
	do_fio(&c__1, "<------------- Delta-M --------------->", (ftnlen)39);
	do_fio(&c__1, "                   Total    Single                   \
        ", (ftnlen)61);
	do_fio(&c__1, "Total    Single", (ftnlen)15);
	do_fio(&c__1, "       Optical   Optical   Scatter   Truncated   ", (
		ftnlen)49);
	do_fio(&c__1, "Optical   Optical   Scatter    Asymm", (ftnlen)36);
	do_fio(&c__1, "         Depth     Depth    Albedo    Fraction     ", (
		ftnlen)51);
	do_fio(&c__1, "Depth     Depth    Albedo   Factor   Temperature", (
		ftnlen)48);
	e_wsfe();
    }
    if (! (*plank)) {
	s_wsfe(&io___340);
	do_fio(&c__1, "<------------- Delta-M --------------->", (ftnlen)39);
	do_fio(&c__1, "                   Total    Single                   \
        ", (ftnlen)61);
	do_fio(&c__1, "Total    Single", (ftnlen)15);
	do_fio(&c__1, "       Optical   Optical   Scatter   Truncated   ", (
		ftnlen)49);
	do_fio(&c__1, "Optical   Optical   Scatter    Asymm", (ftnlen)36);
	do_fio(&c__1, "         Depth     Depth    Albedo    Fraction     ", (
		ftnlen)51);
	do_fio(&c__1, "Depth     Depth    Albedo   Factor", (ftnlen)34);
	e_wsfe();
    }
    yessct = 0.;
    i__2 = *nlyr;
    for (lc = 1; lc <= i__2; ++lc) {
	yessct += ssalb[lc - 1];
/*                                       ** f90 nonadvancing I/O would */
/*                                       ** simplify this a lot (also the */
/*                                       ** two WRITEs above) */
	if (*plank) {
	    s_wsfe(&io___343);
	    do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dtauc[lc - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&tauc[lc], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ssalb[lc - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&flyr[lc - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&dtaucp[lc - 1], (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&taucpr[lc], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&oprim[lc - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pmom[lc * pmom_dim1 + 1 - pmom_offset], (
		    ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&temper[lc - 1], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsfe();
	}
	if (! (*plank)) {
	    s_wsfe(&io___344);
	    do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dtauc[lc - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&tauc[lc], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ssalb[lc - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&flyr[lc - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&dtaucp[lc - 1], (ftnlen)sizeof(doublereal))
		    ;
	    do_fio(&c__1, (char *)&taucpr[lc], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&oprim[lc - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&pmom[lc * pmom_dim1 + 1 - pmom_offset], (
		    ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
/* L10: */
    }
    if (*plank) {
	s_wsfe(&io___345);
	do_fio(&c__1, (char *)&temper[*nlyr], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*prtmom && yessct > 0.) {
	s_wsfe(&io___346);
	do_fio(&c__1, " Layer   Phase Function Moments", (ftnlen)31);
	e_wsfe();
	i__2 = *nlyr;
	for (lc = 1; lc <= i__2; ++lc) {
	    if (ssalb[lc - 1] > 0.) {
		s_wsfe(&io___347);
		do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
		i__1 = *nstr;
		for (k = 0; k <= i__1; ++k) {
		    do_fio(&c__1, (char *)&pmom[k + lc * pmom_dim1 - 
			    pmom_offset], (ftnlen)sizeof(doublereal));
		}
		e_wsfe();
	    }
/* L20: */
	}
    }
    return 0;
} /* prtinp_ */

/* Subroutine */ int prtint_(doublereal *uu, doublereal *utau, integer *ntau, 
	doublereal *umu, integer *numu, doublereal *phi, integer *nphi, 
	integer *maxulv, integer *maxumu)
{
    /* System generated locals */
    integer uu_dim1, uu_dim2, uu_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer jmin, jmax, j, npass, iu, np, lu, lenfmt;

    /* Fortran I/O blocks */
    static cilist io___348 = { 0, 6, 0, "(//,A)", 0 };
    static cilist io___351 = { 0, 6, 0, "(/,A,/,A,/,A)", 0 };
    static cilist io___356 = { 0, 6, 0, "(/,18X,10F11.2)", 0 };
    static cilist io___358 = { 0, 6, 0, "(F10.4,F8.4,1P,10E11.3)", 0 };
    static cilist io___359 = { 0, 6, 0, "(10X,F8.4,1P,10E11.3)", 0 };
    static cilist io___361 = { 0, 6, 0, "(10X,F8.4,1P,10E11.3)", 0 };


/*         Prints the intensity at user polar and azimuthal angles */
/*     All arguments are DISORT input or output variables */
/*   Called by- DISORT */
/*     LENFMT   Max number of azimuth angles PHI that can be printed */
/*                on one line, as set in FORMAT statement */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2 * 1);

    /* Function Body */
    if (*nphi < 1) {
	return 0;
    }
    s_wsfe(&io___348);
    do_fio(&c__1, " *********  I N T E N S I T I E S  *********", (ftnlen)44);
    e_wsfe();
    lenfmt = 10;
    npass = (*nphi - 1) / lenfmt + 1;
    s_wsfe(&io___351);
    do_fio(&c__1, "             Polar   Azimuth angles (degrees)", (ftnlen)45)
	    ;
    do_fio(&c__1, "   Optical   Angle", (ftnlen)18);
    do_fio(&c__1, "    Depth   Cosine", (ftnlen)18);
    e_wsfe();
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	i__2 = npass;
	for (np = 1; np <= i__2; ++np) {
	    jmin = lenfmt * (np - 1) + 1;
/* Computing MIN */
	    i__3 = lenfmt * np;
	    jmax = min(i__3,*nphi);
	    s_wsfe(&io___356);
	    i__3 = jmax;
	    for (j = jmin; j <= i__3; ++j) {
		do_fio(&c__1, (char *)&phi[j - 1], (ftnlen)sizeof(doublereal))
			;
	    }
	    e_wsfe();
	    if (np == 1) {
		s_wsfe(&io___358);
		do_fio(&c__1, (char *)&utau[lu - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&umu[0], (ftnlen)sizeof(doublereal));
		i__3 = jmax;
		for (j = jmin; j <= i__3; ++j) {
		    do_fio(&c__1, (char *)&uu[(lu + j * uu_dim2) * uu_dim1 + 
			    1 - uu_offset], (ftnlen)sizeof(doublereal));
		}
		e_wsfe();
	    }
	    if (np > 1) {
		s_wsfe(&io___359);
		do_fio(&c__1, (char *)&umu[0], (ftnlen)sizeof(doublereal));
		i__3 = jmax;
		for (j = jmin; j <= i__3; ++j) {
		    do_fio(&c__1, (char *)&uu[(lu + j * uu_dim2) * uu_dim1 + 
			    1 - uu_offset], (ftnlen)sizeof(doublereal));
		}
		e_wsfe();
	    }
	    i__3 = *numu;
	    for (iu = 2; iu <= i__3; ++iu) {
		s_wsfe(&io___361);
		do_fio(&c__1, (char *)&umu[iu - 1], (ftnlen)sizeof(doublereal)
			);
		i__4 = jmax;
		for (j = jmin; j <= i__4; ++j) {
		    do_fio(&c__1, (char *)&uu[iu + (lu + j * uu_dim2) * 
			    uu_dim1 - uu_offset], (ftnlen)sizeof(doublereal));
		}
		e_wsfe();
/* L10: */
	    }
/* L20: */
	}
/* L30: */
    }
    return 0;
} /* prtint_ */

/* Subroutine */ int qgausn_(integer *m, doublereal *gmu, doublereal *gwt)
{
    /* Initialized data */

    static doublereal pi = 0.;
    static integer maxit = 1000;
    static doublereal one = 1.;
    static doublereal two = 2.;

    /* System generated locals */
    integer gmu_dim1, gwt_dim1, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double asin(doublereal);
    integer s_rnge(char *, integer, char *, integer);
    double tan(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal cona;
    static integer iter;
    static doublereal prod, p2pri;
    static integer k;
    static doublereal p, t, x;
    extern doublereal d1mach_(integer *);
    static doublereal en;
    static integer nn;
    static doublereal xi;
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static doublereal pm1;
    static integer np1;
    static doublereal pm2;
    static integer lim;
    static doublereal tol, tmp, ppr, nnp1;

/*       Compute weights and abscissae for ordinary Gaussian quadrature */
/*       on the interval (0,1);  that is, such that */
/*           sum(i=1 to M) ( GWT(i) f(GMU(i)) ) */
/*       is a good approximation to */
/*           integral(0 to 1) ( f(x) dx ) */
/*   INPUT :    M       order of quadrature rule */
/*   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M) */
/*             GWT(I)   array of weights (I = 1 TO M) */
/*   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical */
/*                   Integration, Academic Press, New York, pp. 87, 1975 */
/*   METHOD:  Compute the abscissae as roots of the Legendre */
/*            polynomial P-sub-M using a cubically convergent */
/*            refinement of Newton's method.  Compute the */
/*            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note */
/*            that Newton's method can very easily diverge; only a */
/*            very good initial guess can guarantee convergence. */
/*            The initial guess used here has never led to divergence */
/*            even for M up to 1000. */
/*   ACCURACY:  relative error no better than TOL or computer */
/*              precision (machine epsilon), whichever is larger */
/*   INTERNAL VARIABLES: */
/*    ITER      : number of Newton Method iterations */
/*    MAXIT     : maximum allowed iterations of Newton Method */
/*    PM2,PM1,P : 3 successive Legendre polynomials */
/*    PPR       : derivative of Legendre polynomial */
/*    P2PRI     : 2nd derivative of Legendre polynomial */
/*    TOL       : convergence criterion for Legendre poly root iteration */
/*    X,XI      : successive iterates in cubically-convergent version */
/*                of Newtons Method (seeking roots of Legendre poly.) */
/*   Called by- SETDIS, SURFAC */
/*   Calls- D1MACH, ERRMSG */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    gwt_dim1 = *m;
    gmu_dim1 = *m;

    /* Function Body */
    if (pi == 0.) {
	pi = asin(1.) * 2.;
	tol = d1mach_(&c__4) * 10.;
    }
    if (*m < 1) {
	errmsg_("QGAUSN--Bad value of M", &c_true, (ftnlen)22);
    }
    if (*m == 1) {
	gmu[(i__1 = 0) < 1 * gmu_dim1 ? i__1 : s_rnge("gmu", i__1, "qgausn_", 
		(ftnlen)5303)] = .5;
	gwt[(i__1 = 0) < 1 * gwt_dim1 ? i__1 : s_rnge("gwt", i__1, "qgausn_", 
		(ftnlen)5304)] = 1.;
	return 0;
    }
    en = (doublereal) (*m);
    np1 = *m + 1;
    nnp1 = (doublereal) (*m * np1);
/* Computing 3rd power */
    i__1 = *m;
    cona = (doublereal) (*m - 1) / (i__1 * (i__1 * i__1) << 3);
    lim = *m / 2;
    i__1 = lim;
    for (k = 1; k <= i__1; ++k) {
/*                                        ** Initial guess for k-th root */
/*                                        ** of Legendre polynomial, from */
/*                                        ** Davis/Rabinowitz (2.7.3.3a) */
	t = ((k << 2) - 1) * pi / ((*m << 2) + 2);
	x = cos(t + cona / tan(t));
	iter = 0;
/*                                        ** Upward recurrence for */
/*                                        ** Legendre polynomials */
L10:
	++iter;
	pm2 = one;
	pm1 = x;
	i__2 = *m;
	for (nn = 2; nn <= i__2; ++nn) {
	    p = (((nn << 1) - 1) * x * pm1 - (nn - 1) * pm2) / nn;
	    pm2 = pm1;
	    pm1 = p;
/* L20: */
	}
/*                                              ** Newton Method */
/* Computing 2nd power */
	d__1 = x;
	tmp = one / (one - d__1 * d__1);
	ppr = en * (pm2 - x * p) * tmp;
	p2pri = (two * x * ppr - nnp1 * p) * tmp;
	xi = x - p / ppr * (one + p / ppr * p2pri / (two * ppr));
/*                                              ** Check for convergence */
	if ((d__1 = xi - x, abs(d__1)) > tol) {
	    if (iter > maxit) {
		errmsg_("QGAUSN--max iteration count", &c_true, (ftnlen)27);
	    }
	    x = xi;
	    goto L10;
	}
/*                             ** Iteration finished--calculate weights, */
/*                             ** abscissae for (-1,1) */
	gmu[(i__2 = k - 1) < 1 * gmu_dim1 && 0 <= i__2 ? i__2 : s_rnge("gmu", 
		i__2, "qgausn_", (ftnlen)5354)] = -x;
/* Computing 2nd power */
	d__1 = en * pm2;
	gwt[(i__2 = k - 1) < 1 * gwt_dim1 && 0 <= i__2 ? i__2 : s_rnge("gwt", 
		i__2, "qgausn_", (ftnlen)5355)] = two / (tmp * (d__1 * d__1));
	gmu[(i__2 = np1 - k - 1) < 1 * gmu_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"gmu", i__2, "qgausn_", (ftnlen)5356)] = -gmu[(i__3 = k - 1) <
		 1 * gmu_dim1 && 0 <= i__3 ? i__3 : s_rnge("gmu", i__3, "qga\
usn_", (ftnlen)5356)];
	gwt[(i__2 = np1 - k - 1) < 1 * gwt_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"gwt", i__2, "qgausn_", (ftnlen)5357)] = gwt[(i__3 = k - 1) < 
		1 * gwt_dim1 && 0 <= i__3 ? i__3 : s_rnge("gwt", i__3, "qgau\
sn_", (ftnlen)5357)];
/* L30: */
    }
/*                                    ** Set middle abscissa and weight */
/*                                    ** for rules of odd order */
    if (*m % 2 != 0) {
	gmu[(i__1 = lim) < 1 * gmu_dim1 && 0 <= i__1 ? i__1 : s_rnge("gmu", 
		i__1, "qgausn_", (ftnlen)5363)] = 0.;
	prod = one;
	i__1 = *m;
	for (k = 3; k <= i__1; k += 2) {
	    prod = prod * k / (k - 1);
/* L40: */
	}
/* Computing 2nd power */
	d__1 = prod;
	gwt[(i__1 = lim) < 1 * gwt_dim1 && 0 <= i__1 ? i__1 : s_rnge("gwt", 
		i__1, "qgausn_", (ftnlen)5370)] = two / (d__1 * d__1);
    }
/*                                        ** Convert from (-1,1) to (0,1) */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	gmu[(i__2 = k - 1) < 1 * gmu_dim1 && 0 <= i__2 ? i__2 : s_rnge("gmu", 
		i__2, "qgausn_", (ftnlen)5375)] = gmu[(i__3 = k - 1) < 1 * 
		gmu_dim1 && 0 <= i__3 ? i__3 : s_rnge("gmu", i__3, "qgausn_", 
		(ftnlen)5375)] * .5 + .5;
	gwt[(i__2 = k - 1) < 1 * gwt_dim1 && 0 <= i__2 ? i__2 : s_rnge("gwt", 
		i__2, "qgausn_", (ftnlen)5376)] = gwt[(i__3 = k - 1) < 1 * 
		gwt_dim1 && 0 <= i__3 ? i__3 : s_rnge("gwt", i__3, "qgausn_", 
		(ftnlen)5376)] * .5;
/* L50: */
    }
    return 0;
} /* qgausn_ */

doublereal ratio_(doublereal *a, doublereal *b)
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double d_lg10(doublereal *), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal absa, absb, huge__, powa, powb, tiny;
    extern doublereal r1mach_(integer *);
    static doublereal powmin, powmax;

/*        Calculate ratio  A/B  with over- and under-flow protection */
/*        (thanks to Prof. Jeff Dozier for some suggestions here). */
/*        Since this routine takes two logs, it is no speed demon, */
/*        but it is invaluable for comparing results from two runs */
/*        of a program under development. */

/*        NOTE:  In Fortran90, built-in functions TINY and HUGE */
/*               can replace the R1MACH calls. */
/* --------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. */
    if (pass1) {
	tiny = r1mach_(&c__1);
	huge__ = r1mach_(&c__2);
	powmax = d_lg10(&huge__);
	powmin = d_lg10(&tiny);
	pass1 = FALSE_;
    }
    if (*a == 0.) {
	if (*b == 0.) {
	    ret_val = 1.;
	} else {
	    ret_val = 0.;
	}
    } else if (*b == 0.) {
	ret_val = d_sign(&huge__, a);
    } else {
	absa = abs(*a);
	absb = abs(*b);
	powa = d_lg10(&absa);
	powb = d_lg10(&absb);
	if (absa < tiny && absb < tiny) {
	    ret_val = 1.;
	} else if (powa - powb >= powmax) {
	    ret_val = huge__;
	} else if (powa - powb <= powmin) {
	    ret_val = tiny;
	} else {
	    ret_val = absa / absb;
	}
/*                      ** DONT use old trick of determining sign */
/*                      ** from A*B because A*B may (over/under)flow */
	if (*a > 0. && *b < 0. || *a < 0. && *b > 0.) {
	    ret_val = -ret_val;
	}
    }
    return ret_val;
} /* ratio_ */

/* Subroutine */ int slftst_(doublereal *accur, doublereal *albedo, 
	doublereal *btemp, logical *deltam, doublereal *dtauc, doublereal *
	fbeam, doublereal *fisot, integer *ibcnd, logical *lamber, integer *
	nlyr, logical *plank, integer *nphi, integer *numu, integer *nstr, 
	integer *ntau, logical *onlyfl, doublereal *phi, doublereal *phi0, 
	doublereal *pmom, logical *prnt, doublereal *ssalb, doublereal *temis,
	 doublereal *temper, doublereal *ttemp, doublereal *umu, logical *
	usrang, logical *usrtau, doublereal *utau, doublereal *umu0, 
	doublereal *wvnmhi, doublereal *wvnmlo, logical *compar, doublereal *
	flup, doublereal *rfldir, doublereal *rfldn, doublereal *uu)
{
    /* Initialized data */

    static doublereal acc = 1e-4;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal phis, umus, phi0s;
    static integer i__;
    static doublereal umu0s;
    static integer n, nphis, ntaus;
    static doublereal pmoms[5], utaus;
    static logical prnts[7];
    static integer nlyrs, numus, nstrs;
    static doublereal error1, error2, error3, error4;
    static logical ok;
    static doublereal albeds, fbeams;
    static logical lambes;
    static integer ibcnds;
    static logical deltas;
    static doublereal accurs, dtaucs;
    extern logical tstbad_(char *, doublereal *, ftnlen);
    static doublereal ssalbs;
    static logical planks;
    static doublereal btemps, tempes[2];
    extern /* Subroutine */ int errmsg_(char *, logical *, ftnlen);
    static doublereal temiss, fisots;
    static logical onlyfs, usrans;
    static doublereal ttemps;
    static logical usrtas;
    static doublereal wvnmhs, wvnmls;

/*       If  COMPAR = FALSE, save user input values that would otherwise */
/*       be destroyed and replace them with input values for self-test. */
/*       If  COMPAR = TRUE, compare self-test case results with correct */
/*       answers and restore user input values if test is passed. */

/*       (See file 'DISORT.doc' for variable definitions.) */


/*     I N T E R N A L    V A R I A B L E S: */

/*         ACC     Relative accuracy required for passing self-test */

/*         ERRORn  Relative errors in DISORT output variables */

/*         OK      Logical variable for determining failure of self-test */

/*         All variables ending in 'S' are temporary 'S'torage for input */

/*   Called by- DISORT */
/*   Calls- TSTBAD, ERRMSG */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    if (! (*compar)) {
/*                                     ** Save user input values */
	nlyrs = *nlyr;
	dtaucs = *dtauc;
	ssalbs = *ssalb;
	for (n = 0; n <= 4; ++n) {
	    pmoms[(i__1 = n) < 5 && 0 <= i__1 ? i__1 : s_rnge("pmoms", i__1, 
		    "slftst_", (ftnlen)5561)] = pmom[n];
/* L10: */
	}
	nstrs = *nstr;
	usrans = *usrang;
	numus = *numu;
	umus = *umu;
	usrtas = *usrtau;
	ntaus = *ntau;
	utaus = *utau;
	nphis = *nphi;
	phis = *phi;
	ibcnds = *ibcnd;
	fbeams = *fbeam;
	umu0s = *umu0;
	phi0s = *phi0;
	fisots = *fisot;
	lambes = *lamber;
	albeds = *albedo;
	deltas = *deltam;
	onlyfs = *onlyfl;
	accurs = *accur;
	planks = *plank;
	wvnmls = *wvnmlo;
	wvnmhs = *wvnmhi;
	btemps = *btemp;
	ttemps = *ttemp;
	temiss = *temis;
	tempes[0] = temper[0];
	tempes[1] = temper[1];
	for (i__ = 1; i__ <= 7; ++i__) {
	    prnts[(i__1 = i__ - 1) < 7 && 0 <= i__1 ? i__1 : s_rnge("prnts", 
		    i__1, "slftst_", (ftnlen)5593)] = prnt[i__ - 1];
/* L20: */
	}
/*                                     ** Set input values for self-test */
	*nstr = 4;
	*nlyr = 1;
	*dtauc = 1.;
	*ssalb = .9;
/*                          ** Haze L moments */
	pmom[0] = 1.;
	pmom[1] = .8042;
	pmom[2] = .646094;
	pmom[3] = .481851;
	pmom[4] = .359056;
	*usrang = TRUE_;
	*numu = 1;
	*umu = .5;
	*usrtau = TRUE_;
	*ntau = 1;
	*utau = .5;
	*nphi = 1;
	*phi = 90.;
	*ibcnd = 0;
	*fbeam = 3.14159265;
	*umu0 = .866;
	*phi0 = 0.;
	*fisot = 1.;
	*lamber = TRUE_;
	*albedo = .7;
	*deltam = TRUE_;
	*onlyfl = FALSE_;
	*accur = 1e-4;
	*plank = TRUE_;
	*wvnmlo = 0.;
	*wvnmhi = 5e4;
	*btemp = 300.;
	*ttemp = 100.;
	*temis = .8;
	temper[0] = 210.;
	temper[1] = 200.;
	for (i__ = 1; i__ <= 7; ++i__) {
	    prnt[i__ - 1] = FALSE_;
/* L30: */
	}
    } else {
/*                                    ** Compare test case results with */
/*                                    ** correct answers and abort if bad */
	ok = TRUE_;
	error1 = (*uu - 47.86005) / 47.86005;
	error2 = (*rfldir - 1.527286) / 1.527286;
	error3 = (*rfldn - 28.37223) / 28.37223;
	error4 = (*flup - 152.5853) / 152.5853;
	if (abs(error1) > acc) {
	    ok = tstbad_("UU", &error1, (ftnlen)2);
	}
	if (abs(error2) > acc) {
	    ok = tstbad_("RFLDIR", &error2, (ftnlen)6);
	}
	if (abs(error3) > acc) {
	    ok = tstbad_("RFLDN", &error3, (ftnlen)5);
	}
	if (abs(error4) > acc) {
	    ok = tstbad_("FLUP", &error4, (ftnlen)4);
	}
	if (! ok) {
	    errmsg_("DISORT--self-test failed", &c_true, (ftnlen)24);
	}
/*                                      ** Restore user input values */
	*nlyr = nlyrs;
	*dtauc = dtaucs;
	*ssalb = ssalbs;
	for (n = 0; n <= 4; ++n) {
	    pmom[n] = pmoms[(i__1 = n) < 5 && 0 <= i__1 ? i__1 : s_rnge("pmo\
ms", i__1, "slftst_", (ftnlen)5664)];
/* L40: */
	}
	*nstr = nstrs;
	*usrang = usrans;
	*numu = numus;
	*umu = umus;
	*usrtau = usrtas;
	*ntau = ntaus;
	*utau = utaus;
	*nphi = nphis;
	*phi = phis;
	*ibcnd = ibcnds;
	*fbeam = fbeams;
	*umu0 = umu0s;
	*phi0 = phi0s;
	*fisot = fisots;
	*lamber = lambes;
	*albedo = albeds;
	*deltam = deltas;
	*onlyfl = onlyfs;
	*accur = accurs;
	*plank = planks;
	*wvnmlo = wvnmls;
	*wvnmhi = wvnmhs;
	*btemp = btemps;
	*ttemp = ttemps;
	*temis = temiss;
	temper[0] = tempes[0];
	temper[1] = tempes[1];
	for (i__ = 1; i__ <= 7; ++i__) {
	    prnt[i__ - 1] = prnts[(i__1 = i__ - 1) < 7 && 0 <= i__1 ? i__1 : 
		    s_rnge("prnts", i__1, "slftst_", (ftnlen)5696)];
/* L50: */
	}
    }
    return 0;
} /* slftst_ */

/* Subroutine */ int zeroal_(integer *nd1, doublereal *expbea, doublereal *
	flyr, doublereal *oprim, doublereal *taucpr, doublereal *xr0, 
	doublereal *xr1, integer *nd2, doublereal *cmu, doublereal *cwt, 
	doublereal *psi, doublereal *wk, doublereal *z0, doublereal *z1, 
	doublereal *zj, integer *nd3, doublereal *hlpr, doublereal *ylm0, 
	integer *nd4, doublereal *array, doublereal *cc, doublereal *evecc, 
	integer *nd5, doublereal *gl, integer *nd6, doublereal *ylmc, integer 
	*nd7, doublereal *ylmu, integer *nd8, doublereal *kk, doublereal *ll, 
	doublereal *zz, doublereal *zplk0, doublereal *zplk1, integer *nd9, 
	doublereal *gc, integer *nd10, integer *layru, doublereal *utaupr, 
	integer *nd11, doublereal *gu, integer *nd12, doublereal *z0u, 
	doublereal *z1u, doublereal *zbeam, integer *nd13, doublereal *eval, 
	integer *nd14, doublereal *amb, doublereal *apb, integer *nd15, 
	integer *ipvt, doublereal *z__, integer *nd16, doublereal *rfldir, 
	doublereal *rfldn, doublereal *flup, doublereal *uavg, doublereal *
	dfdt, integer *nd17, doublereal *albmed, doublereal *trnmed, integer *
	nd18, doublereal *u0u, integer *nd19, doublereal *uu)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer n;

/*         ZERO ARRAYS; NDn is dimension of all arrays following */
/*         it in the argument list */
/*   Called by- DISORT */
/* -------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
    i__1 = *nd1;
    for (n = 1; n <= i__1; ++n) {
	expbea[n - 1] = 0.;
	flyr[n - 1] = 0.;
	oprim[n - 1] = 0.;
	taucpr[n - 1] = 0.;
	xr0[n - 1] = 0.;
	xr1[n - 1] = 0.;
/* L10: */
    }
    i__1 = *nd2;
    for (n = 1; n <= i__1; ++n) {
	cmu[n - 1] = 0.;
	cwt[n - 1] = 0.;
	psi[n - 1] = 0.;
	wk[n - 1] = 0.;
	z0[n - 1] = 0.;
	z1[n - 1] = 0.;
	zj[n - 1] = 0.;
/* L20: */
    }
    i__1 = *nd3;
    for (n = 1; n <= i__1; ++n) {
	hlpr[n - 1] = 0.;
	ylm0[n - 1] = 0.;
/* L30: */
    }
    i__1 = *nd4;
    for (n = 1; n <= i__1; ++n) {
	array[n - 1] = 0.;
	cc[n - 1] = 0.;
	evecc[n - 1] = 0.;
/* L40: */
    }
    i__1 = *nd5;
    for (n = 1; n <= i__1; ++n) {
	gl[n - 1] = 0.;
/* L50: */
    }
    i__1 = *nd6;
    for (n = 1; n <= i__1; ++n) {
	ylmc[n - 1] = 0.;
/* L60: */
    }
    i__1 = *nd7;
    for (n = 1; n <= i__1; ++n) {
	ylmu[n - 1] = 0.;
/* L70: */
    }
    i__1 = *nd8;
    for (n = 1; n <= i__1; ++n) {
	kk[n - 1] = 0.;
	ll[n - 1] = 0.;
	zz[n - 1] = 0.;
	zplk0[n - 1] = 0.;
	zplk1[n - 1] = 0.;
/* L80: */
    }
    i__1 = *nd9;
    for (n = 1; n <= i__1; ++n) {
	gc[n - 1] = 0.;
/* L90: */
    }
    i__1 = *nd10;
    for (n = 1; n <= i__1; ++n) {
	layru[n - 1] = 0;
	utaupr[n - 1] = 0.;
/* L100: */
    }
    i__1 = *nd11;
    for (n = 1; n <= i__1; ++n) {
	gu[n - 1] = 0.;
/* L110: */
    }
    i__1 = *nd12;
    for (n = 1; n <= i__1; ++n) {
	z0u[n - 1] = 0.;
	z1u[n - 1] = 0.;
	zbeam[n - 1] = 0.;
/* L120: */
    }
    i__1 = *nd13;
    for (n = 1; n <= i__1; ++n) {
	eval[n - 1] = 0.;
/* L130: */
    }
    i__1 = *nd14;
    for (n = 1; n <= i__1; ++n) {
	amb[n - 1] = 0.;
	apb[n - 1] = 0.;
/* L140: */
    }
    i__1 = *nd15;
    for (n = 1; n <= i__1; ++n) {
	ipvt[n - 1] = 0;
	z__[n - 1] = 0.;
/* L150: */
    }
    i__1 = *nd16;
    for (n = 1; n <= i__1; ++n) {
	rfldir[n - 1] = 0.;
	rfldn[n - 1] = 0.;
	flup[n - 1] = 0.;
	uavg[n - 1] = 0.;
	dfdt[n - 1] = 0.;
/* L160: */
    }
    i__1 = *nd17;
    for (n = 1; n <= i__1; ++n) {
	albmed[n - 1] = 0.;
	trnmed[n - 1] = 0.;
/* L170: */
    }
    i__1 = *nd18;
    for (n = 1; n <= i__1; ++n) {
	u0u[n - 1] = 0.;
/* L180: */
    }
    i__1 = *nd19;
    for (n = 1; n <= i__1; ++n) {
	uu[n - 1] = 0.;
/* L190: */
    }
    return 0;
} /* zeroal_ */

/* Subroutine */ int zeroit_(doublereal *a, integer *length)
{
    /* System generated locals */
    integer a_dim1, i__1, i__2;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer l;

/*         Zeros a real array A having LENGTH elements */
/* -------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
    /* Parameter adjustments */
    a_dim1 = *length;

    /* Function Body */
    i__1 = *length;
    for (l = 1; l <= i__1; ++l) {
	a[(i__2 = l - 1) < 1 * a_dim1 && 0 <= i__2 ? i__2 : s_rnge("a", i__2, 
		"zeroit_", (ftnlen)5882)] = 0.;
/* L10: */
    }
    return 0;
} /* zeroit_ */

/* ****************************************************************** */
/* ********** end of DISORT service routines ************************ */
/* ****************************************************************** */
/* ****************************************************************** */
/* ********** IBCND=1 special case routines ************************* */
/* ****************************************************************** */
/* Subroutine */ int albtrn_(doublereal *albedo, doublereal *amb, doublereal *
	apb, doublereal *array, doublereal *b, doublereal *bdr, doublereal *
	cband, doublereal *cc, doublereal *cmu, doublereal *cwt, doublereal *
	dtaucp, doublereal *eval, doublereal *evecc, doublereal *gl, 
	doublereal *gc, doublereal *gu, integer *ipvt, doublereal *kk, 
	doublereal *ll, integer *nlyr, integer *nn, integer *nstr, integer *
	numu, logical *prnt, doublereal *taucpr, doublereal *umu, doublereal *
	u0u, doublereal *wk, doublereal *ylmc, doublereal *ylmu, doublereal *
	z__, doublereal *aad, doublereal *evald, doublereal *eveccd, 
	doublereal *wkd, integer *mi, integer *mi9m2, integer *maxulv, 
	integer *maxumu, integer *mxcmu, integer *mxumu, integer *nnlyri, 
	doublereal *sqt, doublereal *albmed, doublereal *trnmed)
{
    /* System generated locals */
    integer albmed_dim1, amb_dim1, amb_offset, apb_dim1, apb_offset, 
	    array_dim1, array_offset, bdr_dim1, bdr_offset, cband_dim1, 
	    cband_offset, cc_dim1, cc_offset, evecc_dim1, evecc_offset, 
	    gc_dim1, gc_dim2, gc_offset, gl_dim1, gl_offset, gu_dim1, gu_dim2,
	     gu_offset, kk_dim1, kk_offset, ll_dim1, ll_offset, trnmed_dim1, 
	    u0u_dim1, u0u_dim2, u0u_offset, umu_dim1, ylmc_dim1, ylmc_dim2, 
	    ylmc_offset, ylmu_dim1, ylmu_offset, aad_dim1, aad_offset, 
	    eveccd_dim1, eveccd_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    double exp(doublereal);

    /* Local variables */
    static integer ncol, ncut;
    static doublereal delm0;
    static integer l;
    extern /* Subroutine */ int sgbco_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *);
    static doublereal rcond;
    static integer mazim;
    static doublereal fisot;
    extern /* Subroutine */ int solve1_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *);
    static integer lc, iq, iu;
    static logical lamber;
    static doublereal sphalb;
    extern /* Subroutine */ int soleig_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), altrin_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    , errmsg_(char *, logical *, ftnlen), lepoly_(integer *, integer *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *),
	     praltr_(doublereal *, integer *, doublereal *, doublereal *), 
	    spaltr_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), terpev_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal sphtrn;
    extern /* Subroutine */ int zeroit_(doublereal *, integer *);
    static logical lyrcut;
    extern /* Subroutine */ int setmtx_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static integer ncd;
    static doublereal sgn;

/*    DISORT special case to get only albedo and transmissivity */
/*    of entire medium as a function of incident beam angle */
/*    (many simplifications because boundary condition is just */
/*    isotropic illumination, there are no thermal sources, and */
/*    particular solutions do not need to be computed).  See */
/*    Ref. S2 and references therein for details. */
/*    The basic idea is as follows.  The reciprocity principle leads to */
/*    the following relationships for a plane-parallel, vertically */
/*    inhomogeneous medium lacking thermal (or other internal) sources: */

/*       albedo(theta) = u_0(theta) for unit-intensity isotropic */
/*                       illumination at *top* boundary */

/*       trans(theta) =  u_0(theta) for unit-intensity isotropic */
/*                       illumination at *bottom* boundary */

/*    where */

/*       albedo(theta) = albedo for beam incidence at angle theta */
/*       trans(theta) = transmissivity for beam incidence at angle theta */
/*       u_0(theta) = upward azim-avg intensity at top boundary */
/*                    at angle theta */
/*   O U T P U T    V A R I A B L E S: */

/*       ALBMED(IU)   Albedo of the medium as a function of incident */
/*                    beam angle cosine UMU(IU) */

/*       TRNMED(IU)   Transmissivity of the medium as a function of */
/*                    incident beam angle cosine UMU(IU) */
/*    I N T E R N A L   V A R I A B L E S: */
/*       NCD         number of diagonals below/above main diagonal */
/*       RCOND       estimate of the reciprocal condition of matrix */
/*                   CBAND; for system  CBAND*X = B, relative */
/*                   perturbations in CBAND and B of size epsilon may */
/*                   cause relative perturbations in X of size */
/*                   epsilon/RCOND.  If RCOND is so small that */
/*                          1.0 + RCOND .EQ. 1.0 */
/*                   is true, then CBAND may be singular to working */
/*                   precision. */
/*       CBAND       Left-hand side matrix of linear system Eq. SC(5), */
/*                   scaled by Eq. SC(12); in banded form required */
/*                   by LINPACK solution routines */
/*       NCOL        number of columns in CBAND matrix */
/*       IPVT        INTEGER vector of pivot indices */
/*       (most others documented in DISORT) */
/*   Called by- DISORT */
/*   Calls- LEPOLY, ZEROIT, SGBCO, SOLEIG, TERPEV, SETMTX, SOLVE1, */
/*          ALTRIN, SPALTR, PRALTR */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    eveccd_dim1 = *mi;
    eveccd_offset = 1 + eveccd_dim1 * 1;
    aad_dim1 = *mi;
    aad_offset = 1 + aad_dim1 * 1;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    apb_dim1 = *mi;
    apb_offset = 1 + apb_dim1 * 1;
    amb_dim1 = *mi;
    amb_offset = 1 + amb_dim1 * 1;
    trnmed_dim1 = *maxumu;
    albmed_dim1 = *maxumu;
    u0u_dim1 = *maxumu;
    u0u_dim2 = *maxulv;
    u0u_offset = 1 + u0u_dim1 * 1;
    umu_dim1 = *maxumu;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1 * 1;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_dim2 = *mxcmu;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    gl_dim1 = *mxcmu - 0 + 1;
    gl_offset = 0 + gl_dim1 * 1;
    evecc_dim1 = *mxcmu;
    evecc_offset = 1 + evecc_dim1 * 1;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1 * 1;
    array_dim1 = *mxcmu;
    array_offset = 1 + array_dim1 * 1;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2 * 1);
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1 * 1;

    /* Function Body */
    mazim = 0;
    delm0 = 1.;
/*                    ** Set DISORT variables that are ignored in this */
/*                    ** special case but are needed below in argument */
/*                    ** lists of subroutines shared with general case */
    ncut = *nlyr;
    lyrcut = FALSE_;
    fisot = 1.;
    lamber = TRUE_;
/*                          ** Get Legendre polynomials for computational */
/*                          ** and user polar angle cosines */
    i__1 = *nstr - 1;
    lepoly_(numu, &mazim, mxcmu, &i__1, umu, sqt, ylmu);
    i__1 = *nstr - 1;
    lepoly_(nn, &mazim, mxcmu, &i__1, cmu, sqt, ylmc);
/*                       ** Evaluate Legendre polynomials with negative */
/*                       ** arguments from those with positive arguments; */
/*                       ** Dave/Armstrong Eq. (15) */
    sgn = -1.;
    i__1 = *nstr - 1;
    for (l = mazim; l <= i__1; ++l) {
	sgn = -sgn;
	i__2 = *nstr;
	for (iq = *nn + 1; iq <= i__2; ++iq) {
	    ylmc[(i__3 = l + iq * ylmc_dim1 - ylmc_offset) < 1 * ylmc_dim1 * 
		    ylmc_dim2 && 0 <= i__3 ? i__3 : s_rnge("ylmc", i__3, 
		    "albtrn_", (ftnlen)6032)] = sgn * ylmc[(i__4 = l + (iq - *
		    nn) * ylmc_dim1 - ylmc_offset) < 1 * ylmc_dim1 * 
		    ylmc_dim2 && 0 <= i__4 ? i__4 : s_rnge("ylmc", i__4, 
		    "albtrn_", (ftnlen)6032)];
/* L10: */
	}
/* L20: */
    }
/*                                  ** Zero out bottom reflectivity */
/*                                  ** (ALBEDO is used only in analytic */
/*                                  ** formulae involving ALBEDO = 0 */
/*                                  ** solutions; Eqs 16-17 of Ref S2) */
    i__1 = *mi * (*mi + 1);
    zeroit_(bdr, &i__1);
/* ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  ============= */
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
/*                        ** Solve eigenfunction problem in Eq. STWJ(8b) */
	soleig_(amb, apb, array, cmu, cwt, &gl[lc * gl_dim1 - gl_offset], mi, 
		&mazim, mxcmu, nn, nstr, ylmc, cc, evecc, eval, &kk[lc * 
		kk_dim1 + 1 - kk_offset], &gc[(lc * gc_dim2 + 1) * gc_dim1 + 
		1 - gc_offset], aad, eveccd, evald, wkd);
/*                          ** Interpolate eigenvectors to user angles */
	terpev_(cwt, evecc, &gl[lc * gl_dim1 - gl_offset], &gu[(lc * gu_dim2 
		+ 1) * gu_dim1 + 1 - gu_offset], &mazim, mxcmu, mxumu, nn, 
		nstr, numu, wk, ylmc, ylmu);
/* L30: */
    }
/* ===================  END LOOP ON COMPUTATIONAL LAYERS  =============== */
/*                      ** Set coefficient matrix (CBAND) of equations */
/*                      ** combining boundary and layer interface */
/*                      ** conditions (in band-storage mode required by */
/*                      ** LINPACK routines) */
    setmtx_(bdr, cband, cmu, cwt, &delm0, dtaucp, gc, kk, &lamber, &lyrcut, 
	    mi, mi9m2, mxcmu, &ncol, &ncut, nnlyri, nn, nstr, taucpr, wk);
/*                      ** LU-decompose the coeff. matrix (LINPACK) */
    ncd = *nn * 3 - 1;
    sgbco_(cband, mi9m2, &ncol, &ncd, &ncd, ipvt, &rcond, z__);
    if (rcond + 1. == 1.) {
	errmsg_("ALBTRN--SGBCO says matrix near singular", &c_false, (ftnlen)
		39);
    }
/*                             ** First, illuminate from top; if only */
/*                             ** one layer, this will give us everything */
/*                             ** Solve for constants of integration in */
/*                             ** homogeneous solution */
    solve1_(b, cband, &fisot, &c__1, ipvt, ll, mi9m2, mxcmu, &ncol, nlyr, nn, 
	    nnlyri, nstr);
/*                             ** Compute azimuthally-averaged intensity */
/*                             ** at user angles; gives albedo if multi- */
/*                             ** layer (Eq. 9 of Ref S2); gives both */
/*                             ** albedo and transmissivity if single */
/*                             ** layer (Eqs. 3-4 of Ref S2) */
    altrin_(gu, kk, ll, mxcmu, mxumu, maxumu, nlyr, nn, nstr, numu, taucpr, 
	    umu, u0u, wk);
/*                               ** Get beam-incidence albedos from */
/*                               ** reciprocity principle */
    i__1 = *numu / 2;
    for (iu = 1; iu <= i__1; ++iu) {
	albmed[(i__2 = iu - 1) < 1 * albmed_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		"albmed", i__2, "albtrn_", (ftnlen)6102)] = u0u[(i__3 = iu + *
		numu / 2 + u0u_dim1 - u0u_offset) < 1 * u0u_dim1 * u0u_dim2 &&
		 0 <= i__3 ? i__3 : s_rnge("u0u", i__3, "albtrn_", (ftnlen)
		6102)];
/* L40: */
    }
    if (*nlyr == 1) {
	i__1 = *numu / 2;
	for (iu = 1; iu <= i__1; ++iu) {
/*                               ** Get beam-incidence transmissivities */
/*                               ** from reciprocity principle (1 layer); */
/*                               ** flip them end over end to correspond */
/*                               ** to positive UMU instead of negative */
	    trnmed[(i__3 = iu - 1) < 1 * trnmed_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("trnmed", i__3, "albtrn_", (ftnlen)6114)] = u0u[(
		    i__4 = *numu / 2 + 1 - iu + (u0u_dim1 << 1) - u0u_offset) 
		    < 1 * u0u_dim1 * u0u_dim2 && 0 <= i__4 ? i__4 : s_rnge(
		    "u0u", i__4, "albtrn_", (ftnlen)6114)] + exp(-taucpr[*
		    nlyr] / umu[(i__2 = iu + *numu / 2 - 1) < 1 * umu_dim1 && 
		    0 <= i__2 ? i__2 : s_rnge("umu", i__2, "albtrn_", (ftnlen)
		    6114)]);
/* L50: */
	}
    } else {
/*                             ** Second, illuminate from bottom */
/*                             ** (if multiple layers) */
	solve1_(b, cband, &fisot, &c__2, ipvt, ll, mi9m2, mxcmu, &ncol, nlyr, 
		nn, nnlyri, nstr);
	altrin_(gu, kk, ll, mxcmu, mxumu, maxumu, nlyr, nn, nstr, numu, 
		taucpr, umu, u0u, wk);
/*                               ** Get beam-incidence transmissivities */
/*                               ** from reciprocity principle */
	i__1 = *numu / 2;
	for (iu = 1; iu <= i__1; ++iu) {
	    trnmed[(i__3 = iu - 1) < 1 * trnmed_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("trnmed", i__3, "albtrn_", (ftnlen)6132)] = u0u[(
		    i__4 = iu + *numu / 2 + u0u_dim1 - u0u_offset) < 1 * 
		    u0u_dim1 * u0u_dim2 && 0 <= i__4 ? i__4 : s_rnge("u0u", 
		    i__4, "albtrn_", (ftnlen)6132)] + exp(-taucpr[*nlyr] / 
		    umu[(i__2 = iu + *numu / 2 - 1) < 1 * umu_dim1 && 0 <= 
		    i__2 ? i__2 : s_rnge("umu", i__2, "albtrn_", (ftnlen)6132)
		    ]);
/* L60: */
	}
    }
    if (*albedo > 0.) {
/*                             ** Get spherical albedo and transmissivity */
	if (*nlyr == 1) {
	    spaltr_(cmu, cwt, gc, kk, ll, mxcmu, nlyr, nn, nstr, taucpr, &
		    sphalb, &sphtrn);
	} else {
	    spaltr_(cmu, cwt, gc, kk, ll, mxcmu, nlyr, nn, nstr, taucpr, &
		    sphtrn, &sphalb);
	}
/*                                ** Ref. S2, Eqs. 16-17 (these eqs. have */
/*                                ** a simple physical interpretation */
/*                                ** like that of adding-doubling eqs.) */
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    albmed[(i__2 = iu - 1) < 1 * albmed_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("albmed", i__2, "albtrn_", (ftnlen)6157)] = albmed[
		    (i__3 = iu - 1) < 1 * albmed_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("albmed", i__3, "albtrn_", (ftnlen)6157)] + *
		    albedo / (1. - *albedo * sphalb) * sphtrn * trnmed[(i__4 =
		     iu - 1) < 1 * trnmed_dim1 && 0 <= i__4 ? i__4 : s_rnge(
		    "trnmed", i__4, "albtrn_", (ftnlen)6157)];
	    trnmed[(i__2 = iu - 1) < 1 * trnmed_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("trnmed", i__2, "albtrn_", (ftnlen)6160)] = trnmed[
		    (i__3 = iu - 1) < 1 * trnmed_dim1 && 0 <= i__3 ? i__3 : 
		    s_rnge("trnmed", i__3, "albtrn_", (ftnlen)6160)] + *
		    albedo / (1. - *albedo * sphalb) * sphalb * trnmed[(i__4 =
		     iu - 1) < 1 * trnmed_dim1 && 0 <= i__4 ? i__4 : s_rnge(
		    "trnmed", i__4, "albtrn_", (ftnlen)6160)];
/* L70: */
	}
    }
/*                          ** Return UMU to all positive values, to */
/*                          ** agree with ordering in ALBMED, TRNMED */
    *numu /= 2;
    i__1 = *numu;
    for (iu = 1; iu <= i__1; ++iu) {
	umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : s_rnge("umu",
		 i__2, "albtrn_", (ftnlen)6169)] = umu[(i__3 = iu + *numu - 1)
		 < 1 * umu_dim1 && 0 <= i__3 ? i__3 : s_rnge("umu", i__3, 
		"albtrn_", (ftnlen)6169)];
/* L80: */
    }
    if (prnt[5]) {
	praltr_(umu, numu, albmed, trnmed);
    }
    return 0;
} /* albtrn_ */

/* Subroutine */ int altrin_(doublereal *gu, doublereal *kk, doublereal *ll, 
	integer *mxcmu, integer *mxumu, integer *maxumu, integer *nlyr, 
	integer *nn, integer *nstr, integer *numu, doublereal *taucpr, 
	doublereal *umu, doublereal *u0u, doublereal *wk)
{
    /* System generated locals */
    integer gu_dim1, gu_dim2, gu_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, u0u_dim1, u0u_offset, umu_dim1, wk_dim1, i__1, i__2, 
	    i__3, i__4;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);
    double exp(doublereal);

    /* Local variables */
    static doublereal dtau, expn, denom;
    static integer iumin, iumax, lc, iq, iu, lu;
    static doublereal mu, palint, utaupr[2], sgn, exp1, exp2;

/*       Computes azimuthally-averaged intensity at top and bottom */
/*       of medium (related to albedo and transmission of medium by */
/*       reciprocity principles; see Ref S2).  User polar angles are */
/*       used as incident beam angles. (This is a very specialized */
/*       version of USRINT) */

/*       ** NOTE **  User input values of UMU (assumed positive) are */
/*                   temporarily in upper locations of  UMU  and */
/*                   corresponding negatives are in lower locations */
/*                   (this makes GU come out right).  I.e. the contents */
/*                   of the temporary UMU array are: */

/*                     -UMU(NUMU),..., -UMU(1), UMU(1),..., UMU(NUMU) */


/*   I N P U T    V A R I A B L E S: */

/*       GU     :  Eigenvectors interpolated to user polar angles */
/*                   (i.e., g in Eq. SC(1) ) */

/*       KK     :  Eigenvalues of coeff. matrix in Eq. SS(7) */

/*       LL     :  Constants of integration in Eq. SC(1), obtained */
/*                   by solving scaled version of Eq. SC(5); */
/*                   exponential term of Eq. SC(12) not included */

/*       NN     :  Order of double-Gauss quadrature (NSTR/2) */

/*       TAUCPR :  Cumulative optical depth (delta-M-scaled) */

/*       (remainder are DISORT input variables) */


/*   O U T P U T    V A R I A B L E: */

/*       U0U  :    Diffuse azimuthally-averaged intensity at top and */
/*                 bottom of medium (directly transmitted component, */
/*                 corresponding to BNDINT in USRINT, is omitted). */


/*   I N T E R N A L    V A R I A B L E S: */

/*       DTAU   :  Optical depth of a computational layer */
/*       PALINT :  Non-boundary-forced intensity component */
/*       UTAUPR :  Optical depths of user output levels (delta-M scaled) */
/*       WK     :  Scratch vector for saving 'EXP' evaluations */
/*       All the exponential factors (i.e., EXP1, EXPN,... etc.) */
/*       come from the substitution of constants of integration in */
/*       Eq. SC(12) into Eqs. S1(8-9).  All have negative arguments. */

/*   Called by- ALBTRN */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    wk_dim1 = *mxcmu;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2 * 1);
    u0u_dim1 = *maxumu;
    u0u_offset = 1 + u0u_dim1 * 1;
    umu_dim1 = *maxumu;

    /* Function Body */
    utaupr[0] = 0.;
    utaupr[1] = taucpr[*nlyr];
    for (lu = 1; lu <= 2; ++lu) {
	if (lu == 1) {
	    iumin = *numu / 2 + 1;
	    iumax = *numu;
	    sgn = 1.;
	} else {
	    iumin = 1;
	    iumax = *numu / 2;
	    sgn = -1.;
	}
/*                                   ** Loop over polar angles at which */
/*                                   ** albedos/transmissivities desired */
/*                                   ** ( upward angles at top boundary, */
/*                                   ** downward angles at bottom ) */
	i__1 = iumax;
	for (iu = iumin; iu <= i__1; ++iu) {
	    mu = umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : 
		    s_rnge("umu", i__2, "altrin_", (ftnlen)6282)];
/*                                     ** Integrate from top to bottom */
/*                                     ** computational layer */
	    palint = 0.;
	    i__2 = *nlyr;
	    for (lc = 1; lc <= i__2; ++lc) {
		dtau = taucpr[lc] - taucpr[lc - 1];
		exp1 = exp((utaupr[(i__3 = lu - 1) < 2 && 0 <= i__3 ? i__3 : 
			s_rnge("utaupr", i__3, "altrin_", (ftnlen)6290)] - 
			taucpr[lc - 1]) / mu);
		exp2 = exp((utaupr[(i__3 = lu - 1) < 2 && 0 <= i__3 ? i__3 : 
			s_rnge("utaupr", i__3, "altrin_", (ftnlen)6291)] - 
			taucpr[lc]) / mu);
/*                                      ** KK is negative */
		i__3 = *nn;
		for (iq = 1; iq <= i__3; ++iq) {
		    wk[(i__4 = iq - 1) < 1 * wk_dim1 && 0 <= i__4 ? i__4 : 
			    s_rnge("wk", i__4, "altrin_", (ftnlen)6296)] = 
			    exp(kk[iq + lc * kk_dim1 - kk_offset] * dtau);
		    denom = mu * kk[iq + lc * kk_dim1 - kk_offset] + 1.;
		    if (abs(denom) < 1e-4) {
/*                                                   ** L'Hospital limit */
			expn = dtau / mu * exp2;
		    } else {
			expn = (exp1 * wk[(i__4 = iq - 1) < 1 * wk_dim1 && 0 
				<= i__4 ? i__4 : s_rnge("wk", i__4, "altrin_",
				 (ftnlen)6305)] - exp2) * sgn / denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1 - 
			    gu_offset] * ll[iq + lc * ll_dim1 - ll_offset] * 
			    expn;
/* L10: */
		}
/*                                        ** KK is positive */
		i__3 = *nstr;
		for (iq = *nn + 1; iq <= i__3; ++iq) {
		    denom = mu * kk[iq + lc * kk_dim1 - kk_offset] + 1.;
		    if (abs(denom) < 1e-4) {
			expn = -dtau / mu * exp1;
		    } else {
			expn = (exp1 - exp2 * wk[(i__4 = *nstr + 1 - iq - 1) <
				 1 * wk_dim1 && 0 <= i__4 ? i__4 : s_rnge(
				"wk", i__4, "altrin_", (ftnlen)6324)]) * sgn /
				 denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1 - 
			    gu_offset] * ll[iq + lc * ll_dim1 - ll_offset] * 
			    expn;
/* L20: */
		}
/* L30: */
	    }
	    u0u[iu + lu * u0u_dim1 - u0u_offset] = palint;
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* altrin_ */

/* Subroutine */ int praltr_(doublereal *umu, integer *numu, doublereal *
	albmed, doublereal *trnmed)
{
    /* System generated locals */
    integer albmed_dim1, trnmed_dim1, umu_dim1, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(), 
	    s_rnge(char *, integer, char *, integer);
    double acos(doublereal);

    /* Local variables */
    static integer iu;

    /* Fortran I/O blocks */
    static cilist io___466 = { 0, 6, 0, "(///,A,//,A)", 0 };
    static cilist io___468 = { 0, 6, 0, "(0P,F13.4,F20.6,F12.5,1P,E17.4)", 0 }
	    ;


/*        Print planar albedo and transmissivity of medium */
/*        as a function of incident beam angle */
/*   Called by- ALBTRN */
/* -------------------------------------------------------------------- */
/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    trnmed_dim1 = *numu;
    albmed_dim1 = *numu;
    umu_dim1 = *numu;

    /* Function Body */
    s_wsfe(&io___466);
    do_fio(&c__1, " *******  Flux Albedo and/or Transmissivity of entire med\
ium  ********", (ftnlen)70);
    do_fio(&c__1, " Beam Zen Ang   cos(Beam Zen Ang)      Albedo   Transmiss\
ivity", (ftnlen)62);
    e_wsfe();
    i__1 = *numu;
    for (iu = 1; iu <= i__1; ++iu) {
	s_wsfe(&io___468);
	d__1 = acos(umu[(i__2 = iu - 1) < 1 * umu_dim1 && 0 <= i__2 ? i__2 : 
		s_rnge("umu", i__2, "praltr_", (ftnlen)6379)]) * 
		57.295779578552292;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&umu[(i__3 = iu - 1) < 1 * umu_dim1 && 0 <= 
		i__3 ? i__3 : s_rnge("umu", i__3, "praltr_", (ftnlen)6379)], (
		ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&albmed[(i__4 = iu - 1) < 1 * albmed_dim1 && 0 
		<= i__4 ? i__4 : s_rnge("albmed", i__4, "praltr_", (ftnlen)
		6379)], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&trnmed[(i__5 = iu - 1) < 1 * trnmed_dim1 && 0 
		<= i__5 ? i__5 : s_rnge("trnmed", i__5, "praltr_", (ftnlen)
		6379)], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L10: */
    }
    return 0;
} /* praltr_ */

/* Subroutine */ int solve1_(doublereal *b, doublereal *cband, doublereal *
	fisot, integer *ihom, integer *ipvt, doublereal *ll, integer *mi9m2, 
	integer *mxcmu, integer *ncol, integer *ncut, integer *nn, integer *
	nnlyri, integer *nstr)
{
    /* System generated locals */
    integer b_dim1, cband_dim1, cband_offset, ll_dim1, ll_offset, i__1, i__2, 
	    i__3;

    /* Builtin functions */
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static integer ipnt, i__;
    extern /* Subroutine */ int sgbsl_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *);
    static integer lc, iq;
    extern /* Subroutine */ int zeroit_(doublereal *, integer *);
    static integer ncd;

/*        Construct right-hand side vector B for isotropic incidence */
/*        (only) on either top or bottom boundary and solve system */
/*        of equations obtained from the boundary conditions and the */
/*        continuity-of-intensity-at-layer-interface equations */

/*     I N P U T      V A R I A B L E S: */

/*       CBAND    :  Left-hand side matrix of banded linear system */
/*                   Eq. SC(5), scaled by Eq. SC(12); assumed already */
/*                   in LU-decomposed form, ready for LINPACK solver */

/*       IHOM     :  Direction of illumination flag (1, top; 2, bottom) */

/*       NCOL     :  Number of columns in CBAND */

/*       NN       :  Order of double-Gauss quadrature (NSTR/2) */

/*       (remainder are DISORT input variables) */


/*    O U T P U T     V A R I A B L E S: */

/*       B        :  Right-hand side vector of Eq. SC(5) going into */
/*                   SGBSL; returns as solution vector of Eq. */
/*                   SC(12), constants of integration without */
/*                   exponential term */

/*       LL      :   permanent storage for B, but re-ordered */


/*    I N T E R N A L    V A R I A B L E S: */

/*       IPVT     :  INTEGER vector of pivot indices */
/*       NCD      :  Number of diagonals below or above main diagonal */

/*   Called by- ALBTRN */
/*   Calls- ZEROIT, ERRMSG, SGBSL */
/* +-------------------------------------------------------------------+ */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    /* Parameter adjustments */
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1 * 1;
    b_dim1 = *nnlyri;

    /* Function Body */
    zeroit_(b, nnlyri);
    if (*ihom == 1) {
/*                             ** Because there are no beam or emission */
/*                             ** sources, remainder of B array is zero */
	i__1 = *nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[(i__2 = i__ - 1) < 1 * b_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "b", i__2, "solve1_", (ftnlen)6454)] = *fisot;
	    b[(i__2 = *ncol - *nn + i__ - 1) < 1 * b_dim1 && 0 <= i__2 ? i__2 
		    : s_rnge("b", i__2, "solve1_", (ftnlen)6455)] = 0.;
/* L10: */
	}
    } else if (*ihom == 2) {
	i__1 = *nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[(i__2 = i__ - 1) < 1 * b_dim1 && 0 <= i__2 ? i__2 : s_rnge(
		    "b", i__2, "solve1_", (ftnlen)6461)] = 0.;
	    b[(i__2 = *ncol - *nn + i__ - 1) < 1 * b_dim1 && 0 <= i__2 ? i__2 
		    : s_rnge("b", i__2, "solve1_", (ftnlen)6462)] = *fisot;
/* L20: */
	}
    }
    ncd = *nn * 3 - 1;
    sgbsl_(cband, mi9m2, ncol, &ncd, &ncd, ipvt, b, &c__0);
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	ipnt = lc * *nstr - *nn;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ll[*nn + 1 - iq + lc * ll_dim1 - ll_offset] = b[(i__3 = ipnt + 1 
		    - iq - 1) < 1 * b_dim1 && 0 <= i__3 ? i__3 : s_rnge("b", 
		    i__3, "solve1_", (ftnlen)6476)];
	    ll[iq + *nn + lc * ll_dim1 - ll_offset] = b[(i__3 = iq + ipnt - 1)
		     < 1 * b_dim1 && 0 <= i__3 ? i__3 : s_rnge("b", i__3, 
		    "solve1_", (ftnlen)6477)];
/* L30: */
	}
/* L40: */
    }
    return 0;
} /* solve1_ */

/* Subroutine */ int spaltr_(doublereal *cmu, doublereal *cwt, doublereal *gc,
	 doublereal *kk, doublereal *ll, integer *mxcmu, integer *nlyr, 
	integer *nn, integer *nstr, doublereal *taucpr, doublereal *sflup, 
	doublereal *sfldn)
{
    /* System generated locals */
    integer cmu_dim1, cwt_dim1, gc_dim1, gc_dim2, gc_offset, kk_dim1, 
	    kk_offset, ll_dim1, ll_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double exp(doublereal);
    integer s_rnge(char *, integer, char *, integer);

    /* Local variables */
    static doublereal zint;
    static integer iq, jq;

/*       Calculates spherical albedo and transmissivity for the entire */
/*       medium from the m=0 intensity components */
/*       (this is a very specialized version of FLUXES) */


/*    I N P U T    V A R I A B L E S: */

/*       CMU,CWT    Abscissae, weights for Gauss quadrature */
/*                  over angle cosine */

/*       KK      :  Eigenvalues of coeff. matrix in eq. SS(7) */

/*       GC      :  Eigenvectors at polar quadrature angles, SC(1) */

/*       LL      :  Constants of integration in eq. SC(1), obtained */
/*                  by solving scaled version of Eq. SC(5); */
/*                  exponential term of Eq. SC(12) not included */

/*       NN      :  Order of double-Gauss quadrature (NSTR/2) */

/*       (remainder are DISORT input variables) */


/*    O U T P U T   V A R I A B L E S: */

/*       SFLUP   :  Up-flux at top (equivalent to spherical albedo due to */
/*                  reciprocity).  For illumination from below it gives */
/*                  spherical transmissivity */

/*       SFLDN   :  Down-flux at bottom (for single layer, equivalent to */
/*                  spherical transmissivity due to reciprocity) */


/*    I N T E R N A L   V A R I A B L E S: */

/*       ZINT    :  Intensity of m=0 case, in Eq. SC(1) */

/*   Called by- ALBTRN */
/* +-------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    cwt_dim1 = *mxcmu;
    cmu_dim1 = *mxcmu;

    /* Function Body */
    *sflup = 0.;
    i__1 = *nstr;
    for (iq = *nn + 1; iq <= i__1; ++iq) {
	zint = 0.;
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + gc_dim2) * gc_dim1 - gc_offset] * ll[jq + 
		    ll_dim1 - ll_offset] * exp(kk[jq + kk_dim1 - kk_offset] * 
		    taucpr[1]);
/* L10: */
	}
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + gc_dim2) * gc_dim1 - gc_offset] * ll[jq + 
		    ll_dim1 - ll_offset];
/* L20: */
	}
	*sflup += cwt[(i__2 = iq - *nn - 1) < 1 * cwt_dim1 && 0 <= i__2 ? 
		i__2 : s_rnge("cwt", i__2, "spaltr_", (ftnlen)6563)] * cmu[(
		i__3 = iq - *nn - 1) < 1 * cmu_dim1 && 0 <= i__3 ? i__3 : 
		s_rnge("cmu", i__3, "spaltr_", (ftnlen)6563)] * zint;
/* L30: */
    }
    *sfldn = 0.;
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zint = 0.;
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1 - gc_offset] * 
		    ll[jq + *nlyr * ll_dim1 - ll_offset];
/* L40: */
	}
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1 - gc_offset] * 
		    ll[jq + *nlyr * ll_dim1 - ll_offset] * exp(-kk[jq + *nlyr 
		    * kk_dim1 - kk_offset] * (taucpr[*nlyr] - taucpr[*nlyr - 
		    1]));
/* L50: */
	}
	*sfldn += cwt[(i__2 = *nn + 1 - iq - 1) < 1 * cwt_dim1 && 0 <= i__2 ? 
		i__2 : s_rnge("cwt", i__2, "spaltr_", (ftnlen)6583)] * cmu[(
		i__3 = *nn + 1 - iq - 1) < 1 * cmu_dim1 && 0 <= i__3 ? i__3 : 
		s_rnge("cmu", i__3, "spaltr_", (ftnlen)6583)] * zint;
/* L60: */
    }
    *sflup *= 2.;
    *sfldn *= 2.;
    return 0;
} /* spaltr_ */

#ifdef __cplusplus
	}
#endif
