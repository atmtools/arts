/* DISORT.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__7 = 7;
static integer c__6 = 6;
static integer c__5 = 5;
static integer c__10 = 10;
static integer c__48 = 48;
static integer c__3 = 3;
static integer c__1000 = 1000;
static integer c__49 = 49;
static integer c__2304 = 2304;
static integer c__294 = 294;
static integer c__2352 = 2352;
static integer c__490 = 490;
static integer c__288 = 288;
static integer c__13824 = 13824;
static integer c__2880 = 2880;
static integer c__60 = 60;
static integer c__24 = 24;
static integer c__576 = 576;
static integer c__214 = 214;
static integer c__50 = 50;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__0 = 0;
static integer c__100 = 100;

/* ~~~~~~~~~~~~ */
/* VERSION 1.2 */
/* ~~~~~~~~~~~~ */
/* Subroutine */ int disort_(nlyr, dtauc, ssalb, pmom, temper, wvnmlo, wvnmhi,
	 usrtau, ntau, utau, nstr, usrang, numu, umu, nphi, phi, ibcnd, fbeam,
	 umu0, phi0, fisot, lamber, albedo, hl, btemp, ttemp, temis, deltam, 
	plank, onlyfl, accur, prnt, header, maxcly, maxulv, maxumu, maxcmu, 
	maxphi, rfldir, rfldn, flup, dfdt, uavg, uu, u0u, albmed, trnmed, 
	header_len)
integer *nlyr;
real *dtauc, *ssalb, *pmom, *temper, *wvnmlo, *wvnmhi;
logical *usrtau;
integer *ntau;
real *utau;
integer *nstr;
logical *usrang;
integer *numu;
real *umu;
integer *nphi;
real *phi;
integer *ibcnd;
real *fbeam, *umu0, *phi0, *fisot;
logical *lamber;
real *albedo, *hl, *btemp, *ttemp, *temis;
logical *deltam, *plank, *onlyfl;
real *accur;
logical *prnt;
char *header;
integer *maxcly, *maxulv, *maxumu, *maxcmu, *maxphi;
real *rfldir, *rfldn, *flup, *dfdt, *uavg, *uu, *u0u, *albmed, *trnmed;
ftnlen header_len;
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    address a__1[2];
    integer pmom_dim1, pmom_offset, u0u_dim1, u0u_offset, uu_dim1, uu_dim2, 
	    uu_offset, i__1[2], i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5;
    char ch__1[136];

    /* Builtin functions */
    double asin(), sqrt();
    integer i_len(), s_wsfe();
    /* Subroutine */ int s_cat();
    integer do_fio(), e_wsfe();
    double cos();

    /* Local variables */
    static real pkag[7], fldn[5], eval[24];
    static integer ncol;
    static real tauc[7];
    static integer ncos;
    static real ylmc[2352]	/* was [49][48] */, hlpr[49];
    static integer ncut;
    static real flyr[6];
    static integer ipvt[288];
    static real ylmu[490]	/* was [49][10] */, delm0, zplk0[288]	/* 
	    was [48][6] */, zplk1[288]	/* was [48][6] */, b[288], cband[
	    61632]	/* was [214][288] */;
    static integer j, l;
    static real evecc[2304]	/* was [48][48] */, z__[288];
    static doublereal evald[24];
    static real zbeam[60]	/* was [10][6] */, fldir[5];
    static integer mazim;
    static real array[2304]	/* was [48][48] */;
    extern doublereal ratio_();
    static integer kconv;
    static real azerr, oprim[6];
    static integer layru[5];
    static real z0[48], z1[48];
    extern doublereal r1mach_();
    static real cc[2304]	/* was [48][48] */;
    extern /* Subroutine */ int solve0_();
    static real gc[13824]	/* was [48][48][6] */;
    static integer lc;
    static real gl[294]	/* was [49][6] */, pi;
    static integer iq;
    static real kk[288]	/* was [48][6] */;
    static integer nn;
    static real gu[2880]	/* was [10][48][6] */;
    static integer iu;
    static real ll[288]	/* was [48][6] */;
    static doublereal eveccd[576]	/* was [24][24] */;
    static integer lu, ns;
    static real expbea[7], wk[48], bplank, phirad[3], zj[48], angcos;
    extern /* Subroutine */ int chekin_(), upbeam_();
    static real dither, dtaucp[6];
    static logical compar;
    extern /* Subroutine */ int albtrn_();
    static real zz[288]	/* was [48][6] */, cosphi;
    extern doublereal plkavg_();
    extern /* Subroutine */ int soleig_();
    static real tplank, taucpr[7];
    extern /* Subroutine */ int cmpint_(), pravin_();
    static real azterm;
    extern /* Subroutine */ int lepoly_(), fluxes_(), setdis_();
    static real u0c[240]	/* was [48][5] */;
    extern /* Subroutine */ int surfac_(), prtinp_(), slftst_();
    static real utaupr[5];
    extern /* Subroutine */ int prtint_();
    static logical lyrcut;
    extern /* Subroutine */ int setmtx_(), terpev_(), terpso_(), upisot_(), 
	    usrint_();
    static real xr0[6], xr1[6];
    extern /* Subroutine */ int zeroal_(), zeroit_();
    static real z0u[60]	/* was [10][6] */, z1u[60]	/* was [10][6] */;
    static doublereal aad[576]	/* was [24][24] */;
    static real amb[576]	/* was [24][24] */, apb[576]	/* was [24][
	    24] */, bem[24], bdr[600]	/* was [24][25] */, cmu[48], dum;
    static integer lev;
    static real rpd;
    static integer naz;
    static real sgn, emu[10], psi[48];
    static doublereal wkd[48];
    static real cwt[48], rmu[250]	/* was [10][25] */, uum[50]	/* 
	    was [10][5] */, sqt[1000], ylm0[49];

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
    --prnt;
    --ssalb;
    --dtauc;
    --uavg;
    --dfdt;
    --flup;
    --rfldn;
    --rfldir;
    --utau;
    --trnmed;
    --albmed;
    u0u_dim1 = *maxumu;
    u0u_offset = 1 + u0u_dim1 * 1;
    u0u -= u0u_offset;
    --umu;
    pmom_dim1 = *maxcmu - 0 + 1;
    pmom_offset = 0 + pmom_dim1 * 1;
    pmom -= pmom_offset;
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2 * 1);
    uu -= uu_offset;
    --phi;

    /* Function Body */
    if (pass1) {
	pi = asin((float)1.) * (float)2.;
	dither = r1mach_(&c__4) * (float)10.;
/*                            ** Must dither more on Cray (14-digit prec) */
	if (dither < (float)1e-10) {
	    dither *= (float)10.;
	}
	rpd = pi / (float)180.;
	for (ns = 1; ns <= 1000; ++ns) {
	    sqt[ns - 1] = sqrt((real) ns);
/* L5: */
	}
/*                            ** Set input values for self-test. */
/*                            ** Be sure SLFTST sets all print flags off. */
	compar = FALSE_;
	slftst_(accur, albedo, btemp, deltam, &dtauc[1], fbeam, fisot, ibcnd, 
		lamber, nlyr, plank, nphi, numu, nstr, ntau, onlyfl, &phi[1], 
		phi0, &pmom[pmom_dim1], &prnt[1], &ssalb[1], temis, temper, 
		ttemp, &umu[1], usrang, usrtau, &utau[1], umu0, wvnmhi, 
		wvnmlo, &compar, &dum, &dum, &dum, &dum);
    }
L10:
    if (! pass1 && i_len(header, (ftnlen)127) != 0) {
	s_wsfe(&io___9);
/* Writing concatenation */
	i__1[0] = 9, a__1[0] = " DISORT: ";
	i__1[1] = 127, a__1[1] = header;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)136);
	do_fio(&c__1, ch__1, (ftnlen)136);
	e_wsfe();
    }
/*                                  ** Calculate cumulative optical depth */
/*                                  ** and dither single-scatter albedo */
/*                                  ** to improve numerical behavior of */
/*                                  ** eigenvalue/vector computation */
    zeroit_(tauc, &c__7);
    i__2 = *nlyr;
    for (lc = 1; lc <= i__2; ++lc) {
	if (ssalb[lc] == (float)1.) {
	    ssalb[lc] = (float)1. - dither;
	}
	tauc[lc] = tauc[lc - 1] + dtauc[lc];
/* L20: */
    }
/*                                ** Check input dimensions and variables */
    chekin_(nlyr, &dtauc[1], &ssalb[1], &pmom[pmom_offset], temper, wvnmlo, 
	    wvnmhi, usrtau, ntau, &utau[1], nstr, usrang, numu, &umu[1], nphi,
	     &phi[1], ibcnd, fbeam, umu0, phi0, fisot, lamber, albedo, hl, 
	    btemp, ttemp, temis, plank, onlyfl, accur, tauc, maxcly, maxulv, 
	    maxumu, maxcmu, maxphi, &c__6, &c__5, &c__10, &c__48, &c__3, &
	    c__1000);
/*                                 ** Zero internal and output arrays */
    i__2 = *maxumu * *maxulv;
    i__3 = *maxumu * *maxulv * *maxphi;
    zeroal_(&c__6, &expbea[1], flyr, oprim, &taucpr[1], xr0, xr1, &c__48, cmu,
	     cwt, psi, wk, z0, z1, zj, &c__49, hlpr, ylm0, &c__2304, array, 
	    cc, evecc, &c__294, gl, &c__2352, ylmc, &c__490, ylmu, &c__288, 
	    kk, ll, zz, zplk0, zplk1, &c__13824, gc, &c__5, layru, utaupr, &
	    c__2880, gu, &c__60, z0u, z1u, zbeam, &c__24, eval, &c__576, amb, 
	    apb, &c__288, ipvt, z__, maxulv, &rfldir[1], &rfldn[1], &flup[1], 
	    &uavg[1], &dfdt[1], maxumu, &albmed[1], &trnmed[1], &i__2, &u0u[
	    u0u_offset], &i__3, &uu[uu_offset]);
/*                                 ** Perform various setup operations */
    setdis_(cmu, cwt, deltam, &dtauc[1], dtaucp, expbea, fbeam, flyr, gl, hl, 
	    hlpr, ibcnd, lamber, layru, &lyrcut, maxumu, maxcmu, &c__48, &
	    ncut, nlyr, ntau, &nn, nstr, plank, numu, onlyfl, oprim, &pmom[
	    pmom_offset], &ssalb[1], tauc, taucpr, &utau[1], utaupr, &umu[1], 
	    umu0, usrtau, usrang);
/*                                 ** Print input information */
    if (prnt[1]) {
	prtinp_(nlyr, &dtauc[1], dtaucp, &ssalb[1], &pmom[pmom_offset], 
		temper, wvnmlo, wvnmhi, ntau, &utau[1], nstr, numu, &umu[1], 
		nphi, &phi[1], ibcnd, fbeam, umu0, phi0, fisot, lamber, 
		albedo, hl, btemp, ttemp, temis, deltam, plank, onlyfl, accur,
		 flyr, &lyrcut, oprim, tauc, taucpr, maxcmu, &prnt[7]);
    }
/*                              ** Handle special case for getting albedo */
/*                              ** and transmissivity of medium for many */
/*                              ** beam angles at once */
    if (*ibcnd == 1) {
	albtrn_(albedo, amb, apb, array, b, bdr, cband, cc, cmu, cwt, dtaucp, 
		eval, evecc, gl, gc, gu, ipvt, kk, ll, nlyr, &nn, nstr, numu, 
		&prnt[1], taucpr, &umu[1], &u0u[u0u_offset], wk, ylmc, ylmu, 
		z__, aad, evald, eveccd, wkd, &c__24, &c__214, maxulv, maxumu,
		 &c__48, &c__10, &c__288, sqt, &albmed[1], &trnmed[1]);
	return 0;
    }
/*                                   ** Calculate Planck functions */
    if (! (*plank)) {
	bplank = (float)0.;
	tplank = (float)0.;
	zeroit_(pkag, &c__7);
    } else {
	tplank = *temis * plkavg_(wvnmlo, wvnmhi, ttemp);
	bplank = plkavg_(wvnmlo, wvnmhi, btemp);
	i__2 = *nlyr;
	for (lev = 0; lev <= i__2; ++lev) {
	    pkag[lev] = plkavg_(wvnmlo, wvnmhi, &temper[lev]);
/* L30: */
	}
    }
/* ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  ======= */
/*           (EQ STWJ 5) */
    kconv = 0;
    naz = *nstr - 1;
/*                                    ** Azimuth-independent case */
    if (*fbeam == (float)0. || (r__1 = (float)1. - *umu0, dabs(r__1)) < (
	    float)1e-5 || *onlyfl || *numu == 1 && (r__2 = (float)1. - umu[1],
	     dabs(r__2)) < (float)1e-5 || *numu == 1 && (r__3 = umu[1] + (
	    float)1., dabs(r__3)) < (float)1e-5 || *numu == 2 && (r__4 = umu[
	    1] + (float)1., dabs(r__4)) < (float)1e-5 && (r__5 = (float)1. - 
	    umu[2], dabs(r__5)) < (float)1e-5) {
	naz = 0;
    }
    i__2 = naz;
    for (mazim = 0; mazim <= i__2; ++mazim) {
	if (mazim == 0) {
	    delm0 = (float)1.;
	}
	if (mazim > 0) {
	    delm0 = (float)0.;
	}
/*                             ** Get normalized associated Legendre */
/*                             ** polynomials for */
/*                             ** (a) incident beam angle cosine */
/*                             ** (b) computational and user polar angle */
/*                             **     cosines */
	if (*fbeam > (float)0.) {
	    ncos = 1;
	    angcos = -(*umu0);
	    i__3 = *nstr - 1;
	    lepoly_(&ncos, &mazim, &c__48, &i__3, &angcos, sqt, ylm0);
	}
	if (! (*onlyfl) && *usrang) {
	    i__3 = *nstr - 1;
	    lepoly_(numu, &mazim, &c__48, &i__3, &umu[1], sqt, ylmu);
	}
	i__3 = *nstr - 1;
	lepoly_(&nn, &mazim, &c__48, &i__3, cmu, sqt, ylmc);
/*                       ** Get normalized associated Legendre polys. */
/*                       ** with negative arguments from those with */
/*                       ** positive arguments; Dave/Armstrong Eq. (15) */
	sgn = (float)-1.;
	i__3 = *nstr - 1;
	for (l = mazim; l <= i__3; ++l) {
	    sgn = -sgn;
	    i__4 = *nstr;
	    for (iq = nn + 1; iq <= i__4; ++iq) {
		ylmc[l + iq * 49 - 49] = sgn * ylmc[l + (iq - nn) * 49 - 49];
/* L40: */
	    }
/* L50: */
	}
/*                                 ** Specify users bottom reflectivity */
/*                                 ** and emissivity properties */
	if (! lyrcut) {
	    surfac_(albedo, &delm0, fbeam, hlpr, lamber, &c__24, &mazim, &
		    c__48, &c__10, &nn, numu, nstr, onlyfl, &umu[1], usrang, 
		    ylm0, ylmc, ylmu, bdr, emu, bem, rmu, sqt);
	}
/* ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  ============= */
	i__3 = ncut;
	for (lc = 1; lc <= i__3; ++lc) {
/*                        ** Solve eigenfunction problem in Eq. STWJ(8B); */
/*                        ** return eigenvalues and eigenvectors */
	    soleig_(amb, apb, array, cmu, cwt, &gl[lc * 49 - 49], &c__24, &
		    mazim, &c__48, &nn, nstr, ylmc, cc, evecc, eval, &kk[lc * 
		    48 - 48], &gc[(lc * 48 + 1) * 48 - 2352], aad, eveccd, 
		    evald, wkd);
/*                                  ** Calculate particular solutions of */
/*                                  ** Eq.SS(18) for incident beam source */
	    if (*fbeam > (float)0.) {
		upbeam_(array, cc, cmu, &delm0, fbeam, &gl[lc * 49 - 49], 
			ipvt, &mazim, &c__48, &nn, nstr, &pi, umu0, wk, ylm0, 
			ylmc, zj, &zz[lc * 48 - 48]);
	    }
/*                              ** Calculate particular solutions of */
/*                              ** Eq. SS(15) for thermal emission source */
	    if (*plank && mazim == 0) {
		xr1[lc - 1] = (float)0.;
		if (dtaucp[lc - 1] > (float)0.) {
		    xr1[lc - 1] = (pkag[lc] - pkag[lc - 1]) / dtaucp[lc - 1];
		}
		xr0[lc - 1] = pkag[lc - 1] - xr1[lc - 1] * taucpr[lc - 1];
		upisot_(array, cc, cmu, ipvt, &c__48, &nn, nstr, &oprim[lc - 
			1], wk, &xr0[lc - 1], &xr1[lc - 1], z0, z1, &zplk0[lc 
			* 48 - 48], &zplk1[lc * 48 - 48]);
	    }
	    if (! (*onlyfl) && *usrang) {
/*                                            ** Interpolate eigenvectors */
/*                                            ** to user angles */
		terpev_(cwt, evecc, &gl[lc * 49 - 49], &gu[(lc * 48 + 1) * 10 
			- 490], &mazim, &c__48, &c__10, &nn, nstr, numu, wk, 
			ylmc, ylmu);
/*                                            ** Interpolate source terms */
/*                                            ** to user angles */
		terpso_(cwt, &delm0, fbeam, &gl[lc * 49 - 49], &mazim, &c__48,
			 plank, numu, nstr, &oprim[lc - 1], &pi, ylm0, ylmc, 
			ylmu, psi, &xr0[lc - 1], &xr1[lc - 1], z0, zj, &zbeam[
			lc * 10 - 10], &z0u[lc * 10 - 10], &z1u[lc * 10 - 10])
			;
	    }
/* L60: */
	}
/* ===================  END LOOP ON COMPUTATIONAL LAYERS  =============== */
/*                      ** Set coefficient matrix of equations combining */
/*                      ** boundary and layer interface conditions */
	setmtx_(bdr, cband, cmu, cwt, &delm0, dtaucp, gc, kk, lamber, &lyrcut,
		 &c__24, &c__214, &c__48, &ncol, &ncut, &c__288, &nn, nstr, 
		taucpr, wk);
/*                      ** Solve for constants of integration in homo- */
/*                      ** geneous solution (general boundary conditions) */
	solve0_(b, bdr, bem, &bplank, cband, cmu, cwt, expbea, fbeam, fisot, 
		ipvt, lamber, ll, &lyrcut, &mazim, &c__24, &c__214, &c__48, &
		ncol, &ncut, &nn, nstr, &c__288, &pi, &tplank, taucpr, umu0, 
		z__, zz, zplk0, zplk1);
/*                                  ** Compute upward and downward fluxes */
	if (mazim == 0) {
	    fluxes_(cmu, cwt, fbeam, gc, kk, layru, ll, &lyrcut, maxulv, &
		    c__48, &c__5, &ncut, &nn, nstr, ntau, &pi, &prnt[1], &
		    ssalb[1], taucpr, umu0, &utau[1], utaupr, xr0, xr1, zz, 
		    zplk0, zplk1, &dfdt[1], &flup[1], fldn, fldir, &rfldir[1],
		     &rfldn[1], &uavg[1], u0c);
	}
	if (*onlyfl) {
	    if (*maxumu >= *nstr) {
/*                                     ** Save azimuthal-avg intensities */
/*                                     ** at quadrature angles */
		i__3 = *ntau;
		for (lu = 1; lu <= i__3; ++lu) {
		    i__4 = *nstr;
		    for (iq = 1; iq <= i__4; ++iq) {
			u0u[iq + lu * u0u_dim1] = u0c[iq + lu * 48 - 49];
/* L70: */
		    }
/* L80: */
		}
	    }
	    goto L170;
	}
	zeroit_(uum, &c__50);
	if (*usrang) {
/*                                     ** Compute azimuthal intensity */
/*                                     ** components at user angles */
	    usrint_(&bplank, cmu, cwt, &delm0, dtaucp, emu, expbea, fbeam, 
		    fisot, gc, gu, kk, lamber, layru, ll, &lyrcut, &mazim, &
		    c__48, &c__5, &c__10, &ncut, nlyr, &nn, nstr, plank, numu,
		     ntau, &pi, rmu, taucpr, &tplank, &umu[1], umu0, utaupr, 
		    wk, zbeam, z0u, z1u, zz, zplk0, zplk1, uum);
	} else {
/*                                     ** Compute azimuthal intensity */
/*                                     ** components at quadrature angles */
	    cmpint_(fbeam, gc, kk, layru, ll, &lyrcut, &mazim, &c__48, &c__5, 
		    &c__10, &ncut, &nn, nstr, plank, ntau, taucpr, umu0, 
		    utaupr, zz, zplk0, zplk1, uum);
	}
	if (mazim == 0) {
/*                               ** Save azimuthally averaged intensities */
	    i__3 = *ntau;
	    for (lu = 1; lu <= i__3; ++lu) {
		i__4 = *numu;
		for (iu = 1; iu <= i__4; ++iu) {
		    u0u[iu + lu * u0u_dim1] = uum[iu + lu * 10 - 11];
		    i__5 = *nphi;
		    for (j = 1; j <= i__5; ++j) {
			uu[iu + (lu + j * uu_dim2) * uu_dim1] = uum[iu + lu * 
				10 - 11];
/* L90: */
		    }
/* L100: */
		}
/* L110: */
	    }
/*                              ** Print azimuthally averaged intensities */
/*                              ** at user angles */
	    if (prnt[4]) {
		pravin_(&umu[1], numu, maxumu, &utau[1], ntau, &u0u[
			u0u_offset]);
	    }
	    if (naz > 0) {
		zeroit_(phirad, &c__3);
		i__3 = *nphi;
		for (j = 1; j <= i__3; ++j) {
		    phirad[j - 1] = rpd * (phi[j] - *phi0);
/* L120: */
		}
	    }
	} else {
/*                                ** Increment intensity by current */
/*                                ** azimuthal component (Fourier */
/*                                ** cosine series);  Eq SD(2) */
	    azerr = (float)0.;
	    i__3 = *nphi;
	    for (j = 1; j <= i__3; ++j) {
		cosphi = cos(mazim * phirad[j - 1]);
		i__4 = *ntau;
		for (lu = 1; lu <= i__4; ++lu) {
		    i__5 = *numu;
		    for (iu = 1; iu <= i__5; ++iu) {
			azterm = uum[iu + lu * 10 - 11] * cosphi;
			uu[iu + (lu + j * uu_dim2) * uu_dim1] += azterm;
/* Computing MAX */
			r__4 = dabs(azterm);
			r__5 = (r__1 = uu[iu + (lu + j * uu_dim2) * uu_dim1], 
				dabs(r__1));
			r__2 = azerr, r__3 = ratio_(&r__4, &r__5);
			azerr = dmax(r__2,r__3);
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
    if (prnt[5] && ! (*onlyfl)) {
	prtint_(&uu[uu_offset], &utau[1], ntau, &umu[1], numu, &phi[1], nphi, 
		maxulv, maxumu);
    }
    if (pass1) {
/*                                    ** Compare test case results with */
/*                                    ** correct answers and abort if bad */
	compar = TRUE_;
	slftst_(accur, albedo, btemp, deltam, &dtauc[1], fbeam, fisot, ibcnd, 
		lamber, nlyr, plank, nphi, numu, nstr, ntau, onlyfl, &phi[1], 
		phi0, &pmom[pmom_dim1], &prnt[1], &ssalb[1], temis, temper, 
		ttemp, &umu[1], usrang, usrtau, &utau[1], umu0, wvnmhi, 
		wvnmlo, &compar, &flup[1], &rfldir[1], &rfldn[1], &uu[(
		uu_dim2 + 1) * uu_dim1 + 1]);
	pass1 = FALSE_;
	goto L10;
    }
    return 0;
} /* disort_ */

/* Subroutine */ int asymtx_(aa, evec, eval, m, ia, ievec, ier, wkd, aad, 
	evecd, evald)
real *aa, *evec, *eval;
integer *m, *ia, *ievec, *ier;
doublereal *wkd, *aad, *evecd, *evald;
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
    integer aa_dim1, aa_offset, evec_dim1, evec_offset, aad_dim1, aad_offset, 
	    evecd_dim1, evecd_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(), d_sign();

    /* Local variables */
    static doublereal repl, f, g, h__;
    static integer i__, j, k, l, n;
    static doublereal p, q, r__, s, t, scale, w, x, y, z__, rnorm;
    extern doublereal d1mach_();
    static integer n1, n2, ka, lb, ii, in;
    static doublereal uu, vv, discri;
    extern /* Subroutine */ int errmsg_();
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
    --evald;
    --eval;
    evecd_dim1 = *ia;
    evecd_offset = 1 + evecd_dim1 * 1;
    evecd -= evecd_offset;
    aad_dim1 = *ia;
    aad_offset = 1 + aad_dim1 * 1;
    aad -= aad_offset;
    aa_dim1 = *ia;
    aa_offset = 1 + aa_dim1 * 1;
    aa -= aa_offset;
    evec_dim1 = *ievec;
    evec_offset = 1 + evec_dim1 * 1;
    evec -= evec_offset;
    --wkd;

    /* Function Body */
    *ier = 0;
    tol = d1mach_(&c__4);
    if (*m < 1 || *ia < *m || *ievec < *m) {
	errmsg_("ASYMTX--bad input variable(s)", &c_true, (ftnlen)29);
    }
/*                           ** Handle 1x1 and 2x2 special cases */
    if (*m == 1) {
	eval[1] = aa[aa_dim1 + 1];
	evec[evec_dim1 + 1] = (float)1.;
	return 0;
    } else if (*m == 2) {
/* Computing 2nd power */
	r__1 = aa[aa_dim1 + 1] - aa[(aa_dim1 << 1) + 2];
	discri = r__1 * r__1 + aa[(aa_dim1 << 1) + 1] * (float)4. * aa[
		aa_dim1 + 2];
	if (discri < (float)0.) {
	    errmsg_("ASYMTX--complex evals in 2x2 case", &c_true, (ftnlen)33);
	}
	sgn = one;
	if (aa[aa_dim1 + 1] < aa[(aa_dim1 << 1) + 2]) {
	    sgn = -one;
	}
	eval[1] = (aa[aa_dim1 + 1] + aa[(aa_dim1 << 1) + 2] + sgn * sqrt(
		discri)) * (float).5;
	eval[2] = (aa[aa_dim1 + 1] + aa[(aa_dim1 << 1) + 2] - sgn * sqrt(
		discri)) * (float).5;
	evec[evec_dim1 + 1] = (float)1.;
	evec[(evec_dim1 << 1) + 2] = (float)1.;
	if (aa[aa_dim1 + 1] == aa[(aa_dim1 << 1) + 2] && (aa[aa_dim1 + 2] == (
		float)0. || aa[(aa_dim1 << 1) + 1] == (float)0.)) {
	    rnorm = (r__1 = aa[aa_dim1 + 1], dabs(r__1)) + (r__2 = aa[(
		    aa_dim1 << 1) + 1], dabs(r__2)) + (r__3 = aa[aa_dim1 + 2],
		     dabs(r__3)) + (r__4 = aa[(aa_dim1 << 1) + 2], dabs(r__4))
		    ;
	    w = tol * rnorm;
	    evec[evec_dim1 + 2] = aa[aa_dim1 + 2] / w;
	    evec[(evec_dim1 << 1) + 1] = -aa[(aa_dim1 << 1) + 1] / w;
	} else {
	    evec[evec_dim1 + 2] = aa[aa_dim1 + 2] / (eval[1] - aa[(aa_dim1 << 
		    1) + 2]);
	    evec[(evec_dim1 << 1) + 1] = aa[(aa_dim1 << 1) + 1] / (eval[2] - 
		    aa[aa_dim1 + 1]);
	}
	return 0;
    }
/*                               ** Convert single-prec. matrix to double */
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (k = 1; k <= i__2; ++k) {
	    aad[j + k * aad_dim1] = aa[j + k * aa_dim1];
/* L10: */
	}
/* L20: */
    }
/*                                ** Initialize output variables */
    *ier = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	evald[i__] = zero;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    evecd[i__ + j * evecd_dim1] = zero;
/* L30: */
	}
	evecd[i__ + i__ * evecd_dim1] = one;
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
		row += (d__1 = aad[j + i__ * aad_dim1], abs(d__1));
	    }
/* L60: */
	}
	if (row == zero) {
	    wkd[k] = (doublereal) j;
	    if (j != k) {
		i__1 = k;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    repl = aad[i__ + j * aad_dim1];
		    aad[i__ + j * aad_dim1] = aad[i__ + k * aad_dim1];
		    aad[i__ + k * aad_dim1] = repl;
/* L70: */
		}
		i__1 = *m;
		for (i__ = l; i__ <= i__1; ++i__) {
		    repl = aad[j + i__ * aad_dim1];
		    aad[j + i__ * aad_dim1] = aad[k + i__ * aad_dim1];
		    aad[k + i__ * aad_dim1] = repl;
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
		col += (d__1 = aad[i__ + j * aad_dim1], abs(d__1));
	    }
/* L110: */
	}
	if (col == zero) {
	    wkd[l] = (doublereal) j;
	    if (j != l) {
		i__2 = k;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    repl = aad[i__ + j * aad_dim1];
		    aad[i__ + j * aad_dim1] = aad[i__ + l * aad_dim1];
		    aad[i__ + l * aad_dim1] = repl;
/* L120: */
		}
		i__2 = *m;
		for (i__ = l; i__ <= i__2; ++i__) {
		    repl = aad[j + i__ * aad_dim1];
		    aad[j + i__ * aad_dim1] = aad[l + i__ * aad_dim1];
		    aad[l + i__ * aad_dim1] = repl;
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
	wkd[i__] = one;
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
		col += (d__1 = aad[j + i__ * aad_dim1], abs(d__1));
		row += (d__1 = aad[i__ + j * aad_dim1], abs(d__1));
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
	    wkd[i__] *= f;
	    noconv = TRUE_;
	    i__2 = *m;
	    for (j = l; j <= i__2; ++j) {
		aad[i__ + j * aad_dim1] /= f;
/* L200: */
	    }
	    i__2 = k;
	    for (j = 1; j <= i__2; ++j) {
		aad[j + i__ * aad_dim1] *= f;
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
	wkd[n + *m] = zero;
	scale = zero;
/*                                                 ** Scale column */
	i__2 = k;
	for (i__ = n; i__ <= i__2; ++i__) {
	    scale += (d__1 = aad[i__ + (n - 1) * aad_dim1], abs(d__1));
/* L230: */
	}
	if (scale != zero) {
	    i__2 = n;
	    for (i__ = k; i__ >= i__2; --i__) {
		wkd[i__ + *m] = aad[i__ + (n - 1) * aad_dim1] / scale;
/* Computing 2nd power */
		d__1 = wkd[i__ + *m];
		h__ += d__1 * d__1;
/* L240: */
	    }
	    d__1 = sqrt(h__);
	    g = -d_sign(&d__1, &wkd[n + *m]);
	    h__ -= wkd[n + *m] * g;
	    wkd[n + *m] -= g;
/*                                            ** Form (I-(U*UT)/H)*A */
	    i__2 = *m;
	    for (j = n; j <= i__2; ++j) {
		f = zero;
		i__3 = n;
		for (i__ = k; i__ >= i__3; --i__) {
		    f += wkd[i__ + *m] * aad[i__ + j * aad_dim1];
/* L250: */
		}
		i__3 = k;
		for (i__ = n; i__ <= i__3; ++i__) {
		    aad[i__ + j * aad_dim1] -= wkd[i__ + *m] * f / h__;
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
		    f += wkd[j + *m] * aad[i__ + j * aad_dim1];
/* L280: */
		}
		i__3 = k;
		for (j = n; j <= i__3; ++j) {
		    aad[i__ + j * aad_dim1] -= wkd[j + *m] * f / h__;
/* L290: */
		}
/* L300: */
	    }
	    wkd[n + *m] = scale * wkd[n + *m];
	    aad[n + (n - 1) * aad_dim1] = scale * g;
	}
/* L310: */
    }
    i__1 = l;
    for (n = k - 2; n >= i__1; --n) {
	n1 = n + 1;
	n2 = n + 2;
	f = aad[n + 1 + n * aad_dim1];
	if (f != zero) {
	    f *= wkd[n + 1 + *m];
	    i__2 = k;
	    for (i__ = n + 2; i__ <= i__2; ++i__) {
		wkd[i__ + *m] = aad[i__ + n * aad_dim1];
/* L320: */
	    }
	    if (n + 1 <= k) {
		i__2 = *m;
		for (j = 1; j <= i__2; ++j) {
		    g = zero;
		    i__3 = k;
		    for (i__ = n + 1; i__ <= i__3; ++i__) {
			g += wkd[i__ + *m] * evecd[i__ + j * evecd_dim1];
/* L330: */
		    }
		    g /= f;
		    i__3 = k;
		    for (i__ = n + 1; i__ <= i__3; ++i__) {
			evecd[i__ + j * evecd_dim1] += g * wkd[i__ + *m];
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
	    rnorm += (d__1 = aad[i__ + j * aad_dim1], abs(d__1));
/* L380: */
	}
	n = i__;
	if (i__ < l || i__ > k) {
	    evald[i__] = aad[i__ + i__ * aad_dim1];
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
	s = (d__1 = aad[lb - 1 + (lb - 1) * aad_dim1], abs(d__1)) + (d__2 = 
		aad[lb + lb * aad_dim1], abs(d__2));
	if (s == zero) {
	    s = rnorm;
	}
	if ((d__1 = aad[lb + (lb - 1) * aad_dim1], abs(d__1)) <= tol * s) {
	    goto L430;
	}
/* L420: */
    }
L430:
    x = aad[n + n * aad_dim1];
    if (lb == n) {
/*                                        ** One eigenvalue found */
	aad[n + n * aad_dim1] = x + t;
	evald[n] = aad[n + n * aad_dim1];
	n = n1;
	goto L400;
    }
    y = aad[n1 + n1 * aad_dim1];
    w = aad[n + n1 * aad_dim1] * aad[n1 + n * aad_dim1];
    if (lb == n1) {
/*                                        ** Two eigenvalues found */
	p = (y - x) * c2;
/* Computing 2nd power */
	d__1 = p;
	q = d__1 * d__1 + w;
	z__ = sqrt((abs(q)));
	aad[n + n * aad_dim1] = x + t;
	x = aad[n + n * aad_dim1];
	aad[n1 + n1 * aad_dim1] = y + t;
/*                                        ** Real pair */
	z__ = p + d_sign(&z__, &p);
	evald[n1] = x + z__;
	evald[n] = evald[n1];
	if (z__ != zero) {
	    evald[n] = x - w / z__;
	}
	x = aad[n + n1 * aad_dim1];
/*                                  ** Employ scale factor in case */
/*                                  ** X and Z are very small */
	r__ = sqrt(x * x + z__ * z__);
	p = x / r__;
	q = z__ / r__;
/*                                             ** Row modification */
	i__1 = *m;
	for (j = n1; j <= i__1; ++j) {
	    z__ = aad[n1 + j * aad_dim1];
	    aad[n1 + j * aad_dim1] = q * z__ + p * aad[n + j * aad_dim1];
	    aad[n + j * aad_dim1] = q * aad[n + j * aad_dim1] - p * z__;
/* L440: */
	}
/*                                             ** Column modification */
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__ = aad[i__ + n1 * aad_dim1];
	    aad[i__ + n1 * aad_dim1] = q * z__ + p * aad[i__ + n * aad_dim1];
	    aad[i__ + n * aad_dim1] = q * aad[i__ + n * aad_dim1] - p * z__;
/* L450: */
	}
/*                                          ** Accumulate transformations */
	i__1 = k;
	for (i__ = l; i__ <= i__1; ++i__) {
	    z__ = evecd[i__ + n1 * evecd_dim1];
	    evecd[i__ + n1 * evecd_dim1] = q * z__ + p * evecd[i__ + n * 
		    evecd_dim1];
	    evecd[i__ + n * evecd_dim1] = q * evecd[i__ + n * evecd_dim1] - p 
		    * z__;
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
	    aad[i__ + i__ * aad_dim1] -= x;
/* L470: */
	}
	s = (d__1 = aad[n + n1 * aad_dim1], abs(d__1)) + (d__2 = aad[n1 + n2 *
		 aad_dim1], abs(d__2));
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
	z__ = aad[i__ + i__ * aad_dim1];
	r__ = x - z__;
	s = y - z__;
	p = (r__ * s - w) / aad[i__ + 1 + i__ * aad_dim1] + aad[i__ + (i__ + 
		1) * aad_dim1];
	q = aad[i__ + 1 + (i__ + 1) * aad_dim1] - z__ - r__ - s;
	r__ = aad[i__ + 2 + (i__ + 1) * aad_dim1];
	s = abs(p) + abs(q) + abs(r__);
	p /= s;
	q /= s;
	r__ /= s;
	if (i__ == lb) {
	    goto L490;
	}
	uu = (d__1 = aad[i__ + (i__ - 1) * aad_dim1], abs(d__1)) * (abs(q) + 
		abs(r__));
	vv = abs(p) * ((d__1 = aad[i__ - 1 + (i__ - 1) * aad_dim1], abs(d__1))
		 + abs(z__) + (d__2 = aad[i__ + 1 + (i__ + 1) * aad_dim1], 
		abs(d__2)));
	if (uu <= tol * vv) {
	    goto L490;
	}
/* L480: */
    }
L490:
    aad[i__ + 2 + i__ * aad_dim1] = zero;
    i__1 = n;
    for (j = i__ + 3; j <= i__1; ++j) {
	aad[j + (j - 2) * aad_dim1] = zero;
	aad[j + (j - 3) * aad_dim1] = zero;
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
		aad[ka + (ka - 1) * aad_dim1] = -aad[ka + (ka - 1) * aad_dim1]
			;
	    }
	} else {
	    p = aad[ka + (ka - 1) * aad_dim1];
	    q = aad[ka + 1 + (ka - 1) * aad_dim1];
	    r__ = zero;
	    if (notlas) {
		r__ = aad[ka + 2 + (ka - 1) * aad_dim1];
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
	    aad[ka + (ka - 1) * aad_dim1] = -s * x;
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
	    p = aad[ka + j * aad_dim1] + q * aad[ka + 1 + j * aad_dim1];
	    if (notlas) {
		p += r__ * aad[ka + 2 + j * aad_dim1];
		aad[ka + 2 + j * aad_dim1] -= p * z__;
	    }
	    aad[ka + 1 + j * aad_dim1] -= p * y;
	    aad[ka + j * aad_dim1] -= p * x;
/* L510: */
	}
/*                                                 ** Column modification */
/* Computing MIN */
	i__3 = n, i__4 = ka + 3;
	i__2 = min(i__3,i__4);
	for (ii = 1; ii <= i__2; ++ii) {
	    p = x * aad[ii + ka * aad_dim1] + y * aad[ii + (ka + 1) * 
		    aad_dim1];
	    if (notlas) {
		p += z__ * aad[ii + (ka + 2) * aad_dim1];
		aad[ii + (ka + 2) * aad_dim1] -= p * r__;
	    }
	    aad[ii + (ka + 1) * aad_dim1] -= p * q;
	    aad[ii + ka * aad_dim1] -= p;
/* L520: */
	}
/*                                          ** Accumulate transformations */
	i__2 = k;
	for (ii = l; ii <= i__2; ++ii) {
	    p = x * evecd[ii + ka * evecd_dim1] + y * evecd[ii + (ka + 1) * 
		    evecd_dim1];
	    if (notlas) {
		p += z__ * evecd[ii + (ka + 2) * evecd_dim1];
		evecd[ii + (ka + 2) * evecd_dim1] -= p * r__;
	    }
	    evecd[ii + (ka + 1) * evecd_dim1] -= p * q;
	    evecd[ii + ka * evecd_dim1] -= p;
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
	    aad[n + n * aad_dim1] = one;
	    for (i__ = n - 1; i__ >= 1; --i__) {
		w = aad[i__ + i__ * aad_dim1] - evald[n];
		if (w == zero) {
		    w = tol * rnorm;
		}
		r__ = aad[i__ + n * aad_dim1];
		i__1 = n - 1;
		for (j = n2; j <= i__1; ++j) {
		    r__ += aad[i__ + j * aad_dim1] * aad[j + n * aad_dim1];
/* L560: */
		}
		aad[i__ + n * aad_dim1] = -r__ / w;
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
		    evecd[i__ + j * evecd_dim1] = aad[i__ + j * aad_dim1];
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
			z__ += evecd[i__ + n * evecd_dim1] * aad[n + j * 
				aad_dim1];
/* L610: */
		    }
		    evecd[i__ + j * evecd_dim1] = z__;
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
	    evecd[i__ + j * evecd_dim1] *= wkd[i__];
/* L640: */
	}
/* L650: */
    }
/*                           ** Interchange rows if permutations occurred */
    for (i__ = l - 1; i__ >= 1; --i__) {
	j = (integer) wkd[i__];
	if (i__ != j) {
	    i__1 = *m;
	    for (n = 1; n <= i__1; ++n) {
		repl = evecd[i__ + n * evecd_dim1];
		evecd[i__ + n * evecd_dim1] = evecd[j + n * evecd_dim1];
		evecd[j + n * evecd_dim1] = repl;
/* L660: */
	    }
	}
/* L670: */
    }
    i__1 = *m;
    for (i__ = k + 1; i__ <= i__1; ++i__) {
	j = (integer) wkd[i__];
	if (i__ != j) {
	    i__2 = *m;
	    for (n = 1; n <= i__2; ++n) {
		repl = evecd[i__ + n * evecd_dim1];
		evecd[i__ + n * evecd_dim1] = evecd[j + n * evecd_dim1];
		evecd[j + n * evecd_dim1] = repl;
/* L680: */
	    }
	}
/* L690: */
    }
/*                         ** Put results into output arrays */
L700:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	eval[j] = evald[j];
	i__2 = *m;
	for (k = 1; k <= i__2; ++k) {
	    evec[j + k * evec_dim1] = evecd[j + k * evecd_dim1];
/* L710: */
	}
/* L720: */
    }
    return 0;
} /* asymtx_ */

/* Subroutine */ int cmpint_(fbeam, gc, kk, layru, ll, lyrcut, mazim, mxcmu, 
	mxulv, mxumu, ncut, nn, nstr, plank, ntau, taucpr, umu0, utaupr, zz, 
	zplk0, zplk1, uum)
real *fbeam, *gc, *kk;
integer *layru;
real *ll;
logical *lyrcut;
integer *mazim, *mxcmu, *mxulv, *mxumu, *ncut, *nn, *nstr;
logical *plank;
integer *ntau;
real *taucpr, *umu0, *utaupr, *zz, *zplk0, *zplk1, *uum;
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, uum_dim1, uum_offset, zplk0_dim1, zplk0_offset, 
	    zplk1_dim1, zplk1_offset, zz_dim1, zz_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static real zint;
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
    --layru;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1 * 1;
    zplk1 -= zplk1_offset;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1 * 1;
    zplk0 -= zplk0_offset;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1 * 1;
    zz -= zz_offset;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    gc -= gc_offset;
    --utaupr;
    uum_dim1 = *mxumu;
    uum_offset = 1 + uum_dim1 * 1;
    uum -= uum_offset;

    /* Function Body */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	lyu = layru[lu];
	if (*lyrcut && lyu > *ncut) {
	    goto L40;
	}
	i__2 = *nstr;
	for (iq = 1; iq <= i__2; ++iq) {
	    zint = (float)0.;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu]));
/* L10: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu - 1]));
/* L20: */
	    }
	    uum[iq + lu * uum_dim1] = zint;
	    if (*fbeam > (float)0.) {
		uum[iq + lu * uum_dim1] = zint + zz[iq + lyu * zz_dim1] * exp(
			-utaupr[lu] / *umu0);
	    }
	    if (*plank && *mazim == 0) {
		uum[iq + lu * uum_dim1] = uum[iq + lu * uum_dim1] + zplk0[iq 
			+ lyu * zplk0_dim1] + zplk1[iq + lyu * zplk1_dim1] * 
			utaupr[lu];
	    }
/* L30: */
	}
L40:
	;
    }
    return 0;
} /* cmpint_ */

/* Subroutine */ int fluxes_(cmu, cwt, fbeam, gc, kk, layru, ll, lyrcut, 
	maxulv, mxcmu, mxulv, ncut, nn, nstr, ntau, pi, prnt, ssalb, taucpr, 
	umu0, utau, utaupr, xr0, xr1, zz, zplk0, zplk1, dfdt, flup, fldn, 
	fldir, rfldir, rfldn, uavg, u0c)
real *cmu, *cwt, *fbeam, *gc, *kk;
integer *layru;
real *ll;
logical *lyrcut;
integer *maxulv, *mxcmu, *mxulv, *ncut, *nn, *nstr, *ntau;
real *pi;
logical *prnt;
real *ssalb, *taucpr, *umu0, *utau, *utaupr, *xr0, *xr1, *zz, *zplk0, *zplk1, 
	*dfdt, *flup, *fldn, *fldir, *rfldir, *rfldn, *uavg, *u0c;
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, u0c_dim1, u0c_offset, zplk0_dim1, zplk0_offset, 
	    zplk1_dim1, zplk1_offset, zz_dim1, zz_offset, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();
    double exp(), acos();

    /* Local variables */
    static real fact, fnet, zint;
    static integer iq, jq, lu;
    static real dirint, fdntot, plsorc;
    extern /* Subroutine */ int zeroit_();
    static integer lyu;
    static real ang1, ang2;

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
    --uavg;
    --rfldn;
    --rfldir;
    --flup;
    --dfdt;
    --utau;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1 * 1;
    zplk1 -= zplk1_offset;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1 * 1;
    zplk0 -= zplk0_offset;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1 * 1;
    zz -= zz_offset;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    gc -= gc_offset;
    --cwt;
    --cmu;
    u0c_dim1 = *mxcmu;
    u0c_offset = 1 + u0c_dim1 * 1;
    u0c -= u0c_offset;
    --fldir;
    --fldn;
    --utaupr;
    --layru;
    --prnt;
    --ssalb;
    --xr0;
    --xr1;

    /* Function Body */
    if (prnt[2]) {
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
    zeroit_(&u0c[u0c_offset], &i__1);
    zeroit_(&fldir[1], mxulv);
    zeroit_(&fldn[1], mxulv);
/*                                        ** Loop over user levels */
    i__1 = *ntau;
    for (lu = 1; lu <= i__1; ++lu) {
	lyu = layru[lu];
	if (*lyrcut && lyu > *ncut) {
/*                                                ** No radiation reaches */
/*                                                ** this level */
	    fdntot = (float)0.;
	    fnet = (float)0.;
	    plsorc = (float)0.;
	    goto L70;
	}
	if (*fbeam > (float)0.) {
	    fact = exp(-utaupr[lu] / *umu0);
	    dirint = *fbeam * fact;
	    fldir[lu] = *umu0 * (*fbeam * fact);
	    rfldir[lu] = *umu0 * *fbeam * exp(-utau[lu] / *umu0);
	} else {
	    dirint = (float)0.;
	    fldir[lu] = (float)0.;
	    rfldir[lu] = (float)0.;
	}
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    zint = (float)0.;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu]));
/* L10: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu - 1]));
/* L20: */
	    }
	    u0c[iq + lu * u0c_dim1] = zint;
	    if (*fbeam > (float)0.) {
		u0c[iq + lu * u0c_dim1] = zint + zz[iq + lyu * zz_dim1] * 
			fact;
	    }
	    u0c[iq + lu * u0c_dim1] = u0c[iq + lu * u0c_dim1] + zplk0[iq + 
		    lyu * zplk0_dim1] + zplk1[iq + lyu * zplk1_dim1] * utaupr[
		    lu];
	    uavg[lu] += cwt[*nn + 1 - iq] * u0c[iq + lu * u0c_dim1];
	    fldn[lu] += cwt[*nn + 1 - iq] * cmu[*nn + 1 - iq] * u0c[iq + lu * 
		    u0c_dim1];
/* L30: */
	}
	i__2 = *nstr;
	for (iq = *nn + 1; iq <= i__2; ++iq) {
	    zint = (float)0.;
	    i__3 = *nn;
	    for (jq = 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu]));
/* L40: */
	    }
	    i__3 = *nstr;
	    for (jq = *nn + 1; jq <= i__3; ++jq) {
		zint += gc[iq + (jq + lyu * gc_dim2) * gc_dim1] * ll[jq + lyu 
			* ll_dim1] * exp(-kk[jq + lyu * kk_dim1] * (utaupr[lu]
			 - taucpr[lyu - 1]));
/* L50: */
	    }
	    u0c[iq + lu * u0c_dim1] = zint;
	    if (*fbeam > (float)0.) {
		u0c[iq + lu * u0c_dim1] = zint + zz[iq + lyu * zz_dim1] * 
			fact;
	    }
	    u0c[iq + lu * u0c_dim1] = u0c[iq + lu * u0c_dim1] + zplk0[iq + 
		    lyu * zplk0_dim1] + zplk1[iq + lyu * zplk1_dim1] * utaupr[
		    lu];
	    uavg[lu] += cwt[iq - *nn] * u0c[iq + lu * u0c_dim1];
	    flup[lu] += cwt[iq - *nn] * cmu[iq - *nn] * u0c[iq + lu * 
		    u0c_dim1];
/* L60: */
	}
	flup[lu] = *pi * (float)2. * flup[lu];
	fldn[lu] = *pi * (float)2. * fldn[lu];
	fdntot = fldn[lu] + fldir[lu];
	fnet = fdntot - flup[lu];
	rfldn[lu] = fdntot - rfldir[lu];
	uavg[lu] = (*pi * (float)2. * uavg[lu] + dirint) / (*pi * (float)4.);
	plsorc = xr0[lyu] + xr1[lyu] * utaupr[lu];
	dfdt[lu] = ((float)1. - ssalb[lyu]) * (float)4. * *pi * (uavg[lu] - 
		plsorc);
L70:
	if (prnt[2]) {
	    s_wsfe(&io___150);
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&lyu, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&rfldir[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&rfldn[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&fdntot, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&flup[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&fnet, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&uavg[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&plsorc, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dfdt[lu], (ftnlen)sizeof(real));
	    e_wsfe();
	}
/* L80: */
    }
    if (prnt[3]) {
	s_wsfe(&io___151);
	do_fio(&c__1, " ******** AZIMUTHALLY AVERAGED ", (ftnlen)31);
	do_fio(&c__1, "INTENSITIES ( at polar quadrature angles ) *******", (
		ftnlen)50);
	e_wsfe();
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    s_wsfe(&io___152);
	    do_fio(&c__1, " Optical depth =", (ftnlen)16);
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	    do_fio(&c__1, "     Angle (deg)   cos(Angle)     Intensity", (
		    ftnlen)43);
	    do_fio(&c__1, "     Angle (deg)   cos(Angle)     Intensity", (
		    ftnlen)43);
	    e_wsfe();
	    i__2 = *nn;
	    for (iq = 1; iq <= i__2; ++iq) {
		ang1 = (float)180. / *pi * acos(cmu[(*nn << 1) - iq + 1]);
		ang2 = (float)180. / *pi * acos(cmu[iq]);
		s_wsfe(&io___155);
		do_fio(&c__1, (char *)&ang1, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&cmu[(*nn << 1) - iq + 1], (ftnlen)
			sizeof(real));
		do_fio(&c__1, (char *)&u0c[iq + lu * u0c_dim1], (ftnlen)
			sizeof(real));
		do_fio(&c__1, (char *)&ang2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&cmu[iq], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&u0c[iq + *nn + lu * u0c_dim1], (ftnlen)
			sizeof(real));
		e_wsfe();
/* L90: */
	    }
/* L100: */
	}
    }
    return 0;
} /* fluxes_ */

/* Subroutine */ int setdis_(cmu, cwt, deltam, dtauc, dtaucp, expbea, fbeam, 
	flyr, gl, hl, hlpr, ibcnd, lamber, layru, lyrcut, maxumu, maxcmu, 
	mxcmu, ncut, nlyr, ntau, nn, nstr, plank, numu, onlyfl, oprim, pmom, 
	ssalb, tauc, taucpr, utau, utaupr, umu, umu0, usrtau, usrang)
real *cmu, *cwt;
logical *deltam;
real *dtauc, *dtaucp, *expbea, *fbeam, *flyr, *gl, *hl, *hlpr;
integer *ibcnd;
logical *lamber;
integer *layru;
logical *lyrcut;
integer *maxumu, *maxcmu, *mxcmu, *ncut, *nlyr, *ntau, *nn, *nstr;
logical *plank;
integer *numu;
logical *onlyfl;
real *oprim, *pmom, *ssalb, *tauc, *taucpr, *utau, *utaupr, *umu, *umu0;
logical *usrtau, *usrang;
{
    /* Initialized data */

    static real abscut = (float)10.;

    /* System generated locals */
    integer gl_dim1, gl_offset, pmom_dim1, pmom_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static real f;
    static integer k, lc, iq, iu, lu;
    static real abstau;
    extern /* Subroutine */ int errmsg_(), qgausn_();

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
    --dtauc;
    --dtaucp;
    --flyr;
    --layru;
    --umu;
    pmom_dim1 = *maxcmu - 0 + 1;
    pmom_offset = 0 + pmom_dim1 * 1;
    pmom -= pmom_offset;
    gl_dim1 = *mxcmu - 0 + 1;
    gl_offset = 0 + gl_dim1 * 1;
    gl -= gl_offset;
    --cwt;
    --cmu;
    --oprim;
    --ssalb;
    --utau;
    --utaupr;

    /* Function Body */
    if (! (*usrtau)) {
/*                              ** Set output levels at computational */
/*                              ** layer boundaries */
	*ntau = *nlyr + 1;
	i__1 = *ntau - 1;
	for (lc = 0; lc <= i__1; ++lc) {
	    utau[lc + 1] = tauc[lc];
/* L10: */
	}
    }
/*                        ** Apply delta-M scaling and move description */
/*                        ** of computational layers to local variables */
    expbea[0] = (float)1.;
    taucpr[0] = (float)0.;
    abstau = (float)0.;
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	pmom[lc * pmom_dim1] = (float)1.;
	if (abstau < abscut) {
	    *ncut = lc;
	}
	abstau += ((float)1. - ssalb[lc]) * dtauc[lc];
	if (! (*deltam)) {
	    oprim[lc] = ssalb[lc];
	    dtaucp[lc] = dtauc[lc];
	    taucpr[lc] = tauc[lc];
	    i__2 = *nstr - 1;
	    for (k = 0; k <= i__2; ++k) {
		gl[k + lc * gl_dim1] = ((k << 1) + 1) * oprim[lc] * pmom[k + 
			lc * pmom_dim1];
/* L20: */
	    }
	    f = (float)0.;
	} else {
/*                                    ** Do delta-M transformation */
	    f = pmom[*nstr + lc * pmom_dim1];
	    oprim[lc] = ssalb[lc] * ((float)1. - f) / ((float)1. - f * ssalb[
		    lc]);
	    dtaucp[lc] = ((float)1. - f * ssalb[lc]) * dtauc[lc];
	    taucpr[lc] = taucpr[lc - 1] + dtaucp[lc];
	    i__2 = *nstr - 1;
	    for (k = 0; k <= i__2; ++k) {
		gl[k + lc * gl_dim1] = ((k << 1) + 1) * oprim[lc] * (pmom[k + 
			lc * pmom_dim1] - f) / ((float)1. - f);
/* L30: */
	    }
	}
	flyr[lc] = f;
	expbea[lc] = (float)0.;
	if (*fbeam > (float)0.) {
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
	    if (utau[lu] >= tauc[lc - 1] && utau[lu] <= tauc[lc]) {
		goto L60;
	    }
/* L50: */
	}
	lc = *nlyr;
L60:
	utaupr[lu] = utau[lu];
	if (*deltam) {
	    utaupr[lu] = taucpr[lc - 1] + ((float)1. - ssalb[lc] * flyr[lc]) *
		     (utau[lu] - tauc[lc - 1]);
	}
	layru[lu] = lc;
/* L70: */
    }
/*                      ** Calculate computational polar angle cosines */
/*                      ** and associated quadrature weights for Gaussian */
/*                      ** quadrature on the interval (0,1) (upward) */
    *nn = *nstr / 2;
    qgausn_(nn, &cmu[1], &cwt[1]);
/*                                  ** Downward (neg) angles and weights */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	cmu[iq + *nn] = -cmu[iq];
	cwt[iq + *nn] = cwt[iq];
/* L80: */
    }
    if (*fbeam > (float)0.) {
/*                               ** Compare beam angle to comput. angles */
	i__1 = *nn;
	for (iq = 1; iq <= i__1; ++iq) {
	    if ((r__1 = *umu0 - cmu[iq], dabs(r__1)) / *umu0 < (float)1e-4) {
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
	    umu[iu] = -cmu[*nn + 1 - iu];
/* L100: */
	}
	i__1 = *nstr;
	for (iu = *nn + 1; iu <= i__1; ++iu) {
	    umu[iu] = cmu[iu - *nn];
/* L110: */
	}
    }
    if (*usrang && *ibcnd == 1) {
/*                               ** Shift positive user angle cosines to */
/*                               ** upper locations and put negatives */
/*                               ** in lower locations */
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    umu[iu + *numu] = umu[iu];
/* L120: */
	}
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    umu[iu] = -umu[(*numu << 1) + 1 - iu];
/* L130: */
	}
	*numu <<= 1;
    }
    if (! (*lyrcut) && ! (*lamber)) {
	i__1 = *nstr;
	for (k = 0; k <= i__1; ++k) {
	    hlpr[k] = ((k << 1) + 1) * hl[k];
/* L140: */
	}
    }
    return 0;
} /* setdis_ */

/* Subroutine */ int setmtx_(bdr, cband, cmu, cwt, delm0, dtaucp, gc, kk, 
	lamber, lyrcut, mi, mi9m2, mxcmu, ncol, ncut, nnlyri, nn, nstr, 
	taucpr, wk)
real *bdr, *cband, *cmu, *cwt, *delm0, *dtaucp, *gc, *kk;
logical *lamber, *lyrcut;
integer *mi, *mi9m2, *mxcmu, *ncol, *ncut, *nnlyri, *nn, *nstr;
real *taucpr, *wk;
{
    /* System generated locals */
    integer bdr_dim1, bdr_offset, cband_dim1, cband_offset, gc_dim1, gc_dim2, 
	    gc_offset, kk_dim1, kk_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static integer jcol;
    static real expa;
    static integer irow, k, nncol, lc, iq, jq, nshift;
    extern /* Subroutine */ int zeroit_();
    static integer lda, ncd;
    static real sum;

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
    --dtaucp;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    bdr -= bdr_offset;
    --wk;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    gc -= gc_offset;
    --cwt;
    --cmu;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1 * 1;
    cband -= cband_offset;

    /* Function Body */
    i__1 = *mi9m2 * *nnlyri;
    zeroit_(&cband[cband_offset], &i__1);
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
	    wk[iq] = exp(kk[iq + lc * kk_dim1] * dtaucp[lc]);
/* L10: */
	}
	jcol = 0;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ++(*ncol);
	    irow = nshift - jcol;
	    i__3 = *nstr;
	    for (jq = 1; jq <= i__3; ++jq) {
		cband[irow + *nstr + *ncol * cband_dim1] = gc[jq + (iq + lc * 
			gc_dim2) * gc_dim1];
		cband[irow + *ncol * cband_dim1] = -gc[jq + (iq + lc * 
			gc_dim2) * gc_dim1] * wk[iq];
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
		cband[irow + *nstr + *ncol * cband_dim1] = gc[jq + (iq + lc * 
			gc_dim2) * gc_dim1] * wk[*nstr + 1 - iq];
		cband[irow + *ncol * cband_dim1] = -gc[jq + (iq + lc * 
			gc_dim2) * gc_dim1];
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
	expa = exp(kk[iq + kk_dim1] * taucpr[1]);
	irow = nshift - jcol + *nn;
	for (jq = *nn; jq >= 1; --jq) {
	    cband[irow + (jcol + 1) * cband_dim1] = gc[jq + (iq + gc_dim2) * 
		    gc_dim1] * expa;
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
	    cband[irow + (jcol + 1) * cband_dim1] = gc[jq + (iq + gc_dim2) * 
		    gc_dim1];
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
	    if (*lyrcut || *lamber && *delm0 == (float)0.) {
/*                          ** No azimuthal-dependent intensity if Lam- */
/*                          ** bert surface; no intensity component if */
/*                          ** truncated bottom layer */
		cband[irow + nncol * cband_dim1] = gc[jq + (iq + *ncut * 
			gc_dim2) * gc_dim1];
	    } else {
		sum = (float)0.;
		i__3 = *nn;
		for (k = 1; k <= i__3; ++k) {
		    sum += cwt[k] * cmu[k] * bdr[jq - *nn + k * bdr_dim1] * 
			    gc[*nn + 1 - k + (iq + *ncut * gc_dim2) * gc_dim1]
			    ;
/* L110: */
		}
		cband[irow + nncol * cband_dim1] = gc[jq + (iq + *ncut * 
			gc_dim2) * gc_dim1] - (*delm0 + (float)1.) * sum;
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
	expa = wk[*nstr + 1 - iq];
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    if (*lyrcut || *lamber && *delm0 == (float)0.) {
		cband[irow + nncol * cband_dim1] = gc[jq + (iq + *ncut * 
			gc_dim2) * gc_dim1] * expa;
	    } else {
		sum = (float)0.;
		i__3 = *nn;
		for (k = 1; k <= i__3; ++k) {
		    sum += cwt[k] * cmu[k] * bdr[jq - *nn + k * bdr_dim1] * 
			    gc[*nn + 1 - k + (iq + *ncut * gc_dim2) * gc_dim1]
			    ;
/* L140: */
		}
		cband[irow + nncol * cband_dim1] = (gc[jq + (iq + *ncut * 
			gc_dim2) * gc_dim1] - (*delm0 + (float)1.) * sum) * 
			expa;
	    }
	    ++irow;
/* L150: */
	}
	++jcol;
/* L160: */
    }
    return 0;
} /* setmtx_ */

/* Subroutine */ int soleig_(amb, apb, array, cmu, cwt, gl, mi, mazim, mxcmu, 
	nn, nstr, ylmc, cc, evecc, eval, kk, gc, aad, eveccd, evald, wkd)
real *amb, *apb, *array, *cmu, *cwt, *gl;
integer *mi, *mazim, *mxcmu, *nn, *nstr;
real *ylmc, *cc, *evecc, *eval, *kk, *gc;
doublereal *aad, *eveccd, *evald, *wkd;
{
    /* System generated locals */
    integer amb_dim1, amb_offset, apb_dim1, apb_offset, array_dim1, 
	    array_offset, cc_dim1, cc_offset, evecc_dim1, evecc_offset, 
	    gc_dim1, gc_offset, ylmc_dim1, ylmc_offset, aad_dim1, aad_offset, 
	    eveccd_dim1, eveccd_offset, i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();
    double sqrt();

    /* Local variables */
    static real beta;
    static integer l;
    static real alpha;
    static integer iq, jq, kq;
    static real gpmigm, gpplgm;
    extern /* Subroutine */ int errmsg_(), asymtx_();
    static integer ier;
    static real sum;

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
    --evald;
    eveccd_dim1 = *mi;
    eveccd_offset = 1 + eveccd_dim1 * 1;
    eveccd -= eveccd_offset;
    aad_dim1 = *mi;
    aad_offset = 1 + aad_dim1 * 1;
    aad -= aad_offset;
    --eval;
    array_dim1 = *mi;
    array_offset = 1 + array_dim1 * 1;
    array -= array_offset;
    apb_dim1 = *mi;
    apb_offset = 1 + apb_dim1 * 1;
    apb -= apb_offset;
    amb_dim1 = *mi;
    amb_offset = 1 + amb_dim1 * 1;
    amb -= amb_offset;
    --wkd;
    gc_dim1 = *mxcmu;
    gc_offset = 1 + gc_dim1 * 1;
    gc -= gc_offset;
    --kk;
    evecc_dim1 = *mxcmu;
    evecc_offset = 1 + evecc_dim1 * 1;
    evecc -= evecc_offset;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1 * 1;
    cc -= cc_offset;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylmc -= ylmc_offset;
    --cwt;
    --cmu;

    /* Function Body */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    sum = (float)0.;
	    i__3 = *nstr - 1;
	    for (l = *mazim; l <= i__3; ++l) {
		sum += gl[l] * ylmc[l + iq * ylmc_dim1] * ylmc[l + jq * 
			ylmc_dim1];
/* L10: */
	    }
	    cc[iq + jq * cc_dim1] = sum * (float).5 * cwt[jq];
/* L20: */
	}
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
/*                             ** Fill remainder of array using symmetry */
/*                             ** relations  C(-mui,muj) = C(mui,-muj) */
/*                             ** and        C(-mui,-muj) = C(mui,muj) */
	    cc[iq + *nn + jq * cc_dim1] = cc[iq + (jq + *nn) * cc_dim1];
	    cc[iq + *nn + (jq + *nn) * cc_dim1] = cc[iq + jq * cc_dim1];
/*                                       ** Get factors of coeff. matrix */
/*                                       ** of reduced eigenvalue problem */
	    alpha = cc[iq + jq * cc_dim1] / cmu[iq];
	    beta = cc[iq + (jq + *nn) * cc_dim1] / cmu[iq];
	    amb[iq + jq * amb_dim1] = alpha - beta;
	    apb[iq + jq * apb_dim1] = alpha + beta;
/* L30: */
	}
	amb[iq + iq * amb_dim1] -= (float)1. / cmu[iq];
	apb[iq + iq * apb_dim1] -= (float)1. / cmu[iq];
/* L40: */
    }
/*                      ** Finish calculation of coefficient matrix of */
/*                      ** reduced eigenvalue problem:  get matrix */
/*                      ** product (alfa+beta)*(alfa-beta); SS(12) */
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    sum = (float)0.;
	    i__3 = *nn;
	    for (kq = 1; kq <= i__3; ++kq) {
		sum += apb[iq + kq * apb_dim1] * amb[kq + jq * amb_dim1];
/* L50: */
	    }
	    array[iq + jq * array_dim1] = sum;
/* L60: */
	}
/* L70: */
    }
/*                      ** Find (real) eigenvalues and eigenvectors */
    asymtx_(&array[array_offset], &evecc[evecc_offset], &eval[1], nn, mi, 
	    mxcmu, &ier, &wkd[1], &aad[aad_offset], &eveccd[eveccd_offset], &
	    evald[1]);
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
	eval[iq] = sqrt((r__1 = eval[iq], dabs(r__1)));
	kk[iq + *nn] = eval[iq];
/*                                      ** Add negative eigenvalue */
	kk[*nn + 1 - iq] = -eval[iq];
/* L80: */
    }
/*                          ** Find eigenvectors (G+) + (G-) from SS(10) */
/*                          ** and store temporarily in APB array */
    i__1 = *nn;
    for (jq = 1; jq <= i__1; ++jq) {
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    sum = (float)0.;
	    i__3 = *nn;
	    for (kq = 1; kq <= i__3; ++kq) {
		sum += amb[iq + kq * amb_dim1] * evecc[kq + jq * evecc_dim1];
/* L90: */
	    }
	    apb[iq + jq * apb_dim1] = sum / eval[jq];
/* L100: */
	}
/* L110: */
    }
    i__1 = *nn;
    for (jq = 1; jq <= i__1; ++jq) {
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    gpplgm = apb[iq + jq * apb_dim1];
	    gpmigm = evecc[iq + jq * evecc_dim1];
/*                                ** Recover eigenvectors G+,G- from */
/*                                ** their sum and difference; stack them */
/*                                ** to get eigenvectors of full system */
/*                                ** SS(7) (JQ = eigenvector number) */
	    evecc[iq + jq * evecc_dim1] = (gpplgm + gpmigm) * (float).5;
	    evecc[iq + *nn + jq * evecc_dim1] = (gpplgm - gpmigm) * (float).5;
/*                                ** Eigenvectors corresponding to */
/*                                ** negative eigenvalues (corresp. to */
/*                                ** reversing sign of 'k' in SS(10) ) */
	    gpplgm = -gpplgm;
	    evecc[iq + (jq + *nn) * evecc_dim1] = (gpplgm + gpmigm) * (float)
		    .5;
	    evecc[iq + *nn + (jq + *nn) * evecc_dim1] = (gpplgm - gpmigm) * (
		    float).5;
	    gc[iq + *nn + (jq + *nn) * gc_dim1] = evecc[iq + jq * evecc_dim1];
	    gc[*nn + 1 - iq + (jq + *nn) * gc_dim1] = evecc[iq + *nn + jq * 
		    evecc_dim1];
	    gc[iq + *nn + (*nn + 1 - jq) * gc_dim1] = evecc[iq + (jq + *nn) * 
		    evecc_dim1];
	    gc[*nn + 1 - iq + (*nn + 1 - jq) * gc_dim1] = evecc[iq + *nn + (
		    jq + *nn) * evecc_dim1];
/* L120: */
	}
/* L130: */
    }
    return 0;
} /* soleig_ */

/* Subroutine */ int solve0_(b, bdr, bem, bplank, cband, cmu, cwt, expbea, 
	fbeam, fisot, ipvt, lamber, ll, lyrcut, mazim, mi, mi9m2, mxcmu, ncol,
	 ncut, nn, nstr, nnlyri, pi, tplank, taucpr, umu0, z__, zz, zplk0, 
	zplk1)
real *b, *bdr, *bem, *bplank, *cband, *cmu, *cwt, *expbea, *fbeam, *fisot;
integer *ipvt;
logical *lamber;
real *ll;
logical *lyrcut;
integer *mazim, *mi, *mi9m2, *mxcmu, *ncol, *ncut, *nn, *nstr, *nnlyri;
real *pi, *tplank, *taucpr, *umu0, *z__, *zz, *zplk0, *zplk1;
{
    /* System generated locals */
    integer bdr_dim1, bdr_offset, cband_dim1, cband_offset, ll_dim1, 
	    ll_offset, zplk0_dim1, zplk0_offset, zplk1_dim1, zplk1_offset, 
	    zz_dim1, zz_offset, i__1, i__2;

    /* Local variables */
    static integer ipnt;
    extern /* Subroutine */ int sgbco_();
    static real rcond;
    extern /* Subroutine */ int sgbsl_();
    static integer lc, iq, jq, it;
    extern /* Subroutine */ int errmsg_(), zeroit_();
    static integer ncd;
    static real sum;

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
    --ipvt;
    --bem;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    bdr -= bdr_offset;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1 * 1;
    zplk1 -= zplk1_offset;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1 * 1;
    zplk0 -= zplk0_offset;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1 * 1;
    zz -= zz_offset;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    ll -= ll_offset;
    --cwt;
    --cmu;
    --z__;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1 * 1;
    cband -= cband_offset;
    --b;

    /* Function Body */
    zeroit_(&b[1], nnlyri);
/*                              ** Construct B,  STWJ(20a,c) for */
/*                              ** parallel beam + bottom reflection + */
/*                              ** thermal emission at top and/or bottom */
    if (*mazim > 0 && *fbeam > (float)0.) {
/*                                         ** Azimuth-dependent case */
/*                                         ** (never called if FBEAM = 0) */
	if (*lyrcut || *lamber) {
/*               ** No azimuthal-dependent intensity for Lambert surface; */
/*               ** no intensity component for truncated bottom layer */
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
/*                                                  ** Top boundary */
		b[iq] = -zz[*nn + 1 - iq + zz_dim1];
/*                                                  ** Bottom boundary */
		b[*ncol - *nn + iq] = -zz[iq + *nn + *ncut * zz_dim1] * 
			expbea[*ncut];
/* L10: */
	    }
	} else {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
		b[iq] = -zz[*nn + 1 - iq + zz_dim1];
		sum = (float)0.;
		i__2 = *nn;
		for (jq = 1; jq <= i__2; ++jq) {
		    sum += cwt[jq] * cmu[jq] * bdr[iq + jq * bdr_dim1] * zz[*
			    nn + 1 - jq + *ncut * zz_dim1] * expbea[*ncut];
/* L20: */
		}
		b[*ncol - *nn + iq] = sum;
		if (*fbeam > (float)0.) {
		    b[*ncol - *nn + iq] = sum + (bdr[iq] * *umu0 * *fbeam / *
			    pi - zz[iq + *nn + *ncut * zz_dim1]) * expbea[*
			    ncut];
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
		b[it] = (zz[iq + (lc + 1) * zz_dim1] - zz[iq + lc * zz_dim1]) 
			* expbea[lc];
/* L40: */
	    }
/* L50: */
	}
    } else {
/*                                   ** Azimuth-independent case */
	if (*fbeam == (float)0.) {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
/*                                      ** Top boundary */
		b[iq] = -zplk0[*nn + 1 - iq + zplk0_dim1] + *fisot + *tplank;
/* L60: */
	    }
	    if (*lyrcut) {
/*                               ** No intensity component for truncated */
/*                               ** bottom layer */
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
/*                                      ** Bottom boundary */
		    b[*ncol - *nn + iq] = -zplk0[iq + *nn + *ncut * 
			    zplk0_dim1] - zplk1[iq + *nn + *ncut * zplk1_dim1]
			     * taucpr[*ncut];
/* L70: */
		}
	    } else {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    sum = (float)0.;
		    i__2 = *nn;
		    for (jq = 1; jq <= i__2; ++jq) {
			sum += cwt[jq] * cmu[jq] * bdr[iq + jq * bdr_dim1] * (
				zplk0[*nn + 1 - jq + *ncut * zplk0_dim1] + 
				zplk1[*nn + 1 - jq + *ncut * zplk1_dim1] * 
				taucpr[*ncut]);
/* L80: */
		    }
		    b[*ncol - *nn + iq] = sum * (float)2. + bem[iq] * *bplank 
			    - zplk0[iq + *nn + *ncut * zplk0_dim1] - zplk1[iq 
			    + *nn + *ncut * zplk1_dim1] * taucpr[*ncut];
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
		    b[it] = zplk0[iq + (lc + 1) * zplk0_dim1] - zplk0[iq + lc 
			    * zplk0_dim1] + (zplk1[iq + (lc + 1) * zplk1_dim1]
			     - zplk1[iq + lc * zplk1_dim1]) * taucpr[lc];
/* L100: */
		}
/* L110: */
	    }
	} else {
	    i__1 = *nn;
	    for (iq = 1; iq <= i__1; ++iq) {
		b[iq] = -zz[*nn + 1 - iq + zz_dim1] - zplk0[*nn + 1 - iq + 
			zplk0_dim1] + *fisot + *tplank;
/* L120: */
	    }
	    if (*lyrcut) {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    b[*ncol - *nn + iq] = -zz[iq + *nn + *ncut * zz_dim1] * 
			    expbea[*ncut] - zplk0[iq + *nn + *ncut * 
			    zplk0_dim1] - zplk1[iq + *nn + *ncut * zplk1_dim1]
			     * taucpr[*ncut];
/* L130: */
		}
	    } else {
		i__1 = *nn;
		for (iq = 1; iq <= i__1; ++iq) {
		    sum = (float)0.;
		    i__2 = *nn;
		    for (jq = 1; jq <= i__2; ++jq) {
			sum += cwt[jq] * cmu[jq] * bdr[iq + jq * bdr_dim1] * (
				zz[*nn + 1 - jq + *ncut * zz_dim1] * expbea[*
				ncut] + zplk0[*nn + 1 - jq + *ncut * 
				zplk0_dim1] + zplk1[*nn + 1 - jq + *ncut * 
				zplk1_dim1] * taucpr[*ncut]);
/* L140: */
		    }
		    b[*ncol - *nn + iq] = sum * (float)2. + (bdr[iq] * *umu0 *
			     *fbeam / *pi - zz[iq + *nn + *ncut * zz_dim1]) * 
			    expbea[*ncut] + bem[iq] * *bplank - zplk0[iq + *
			    nn + *ncut * zplk0_dim1] - zplk1[iq + *nn + *ncut 
			    * zplk1_dim1] * taucpr[*ncut];
/* L150: */
		}
	    }
	    it = *nn;
	    i__1 = *ncut - 1;
	    for (lc = 1; lc <= i__1; ++lc) {
		i__2 = *nstr;
		for (iq = 1; iq <= i__2; ++iq) {
		    ++it;
		    b[it] = (zz[iq + (lc + 1) * zz_dim1] - zz[iq + lc * 
			    zz_dim1]) * expbea[lc] + zplk0[iq + (lc + 1) * 
			    zplk0_dim1] - zplk0[iq + lc * zplk0_dim1] + (
			    zplk1[iq + (lc + 1) * zplk1_dim1] - zplk1[iq + lc 
			    * zplk1_dim1]) * taucpr[lc];
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
    rcond = (float)0.;
    ncd = *nn * 3 - 1;
    sgbco_(&cband[cband_offset], mi9m2, ncol, &ncd, &ncd, &ipvt[1], &rcond, &
	    z__[1]);
    if (rcond + (float)1. == (float)1.) {
	errmsg_("SOLVE0--SGBCO says matrix near singular", &c_false, (ftnlen)
		39);
    }
/*                   ** Solve linear system with coeff matrix CBAND */
/*                   ** and R.H. side(s) B after CBAND has been L-U */
/*                   ** decomposed.  Solution is returned in B. */
    sgbsl_(&cband[cband_offset], mi9m2, ncol, &ncd, &ncd, &ipvt[1], &b[1], &
	    c__0);
/*                   ** Zero CBAND (it may contain 'foreign' */
/*                   ** elements upon returning from LINPACK); */
/*                   ** necessary to prevent errors */
    i__1 = *mi9m2 * *nnlyri;
    zeroit_(&cband[cband_offset], &i__1);
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	ipnt = lc * *nstr - *nn;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ll[*nn + 1 - iq + lc * ll_dim1] = b[ipnt + 1 - iq];
	    ll[iq + *nn + lc * ll_dim1] = b[iq + ipnt];
/* L180: */
	}
/* L190: */
    }
    return 0;
} /* solve0_ */

/* Subroutine */ int surfac_(albedo, delm0, fbeam, hlpr, lamber, mi, mazim, 
	mxcmu, mxumu, nn, numu, nstr, onlyfl, umu, usrang, ylm0, ylmc, ylmu, 
	bdr, emu, bem, rmu, sqt)
real *albedo, *delm0, *fbeam, *hlpr;
logical *lamber;
integer *mi, *mazim, *mxcmu, *mxumu, *nn, *numu, *nstr;
logical *onlyfl;
real *umu;
logical *usrang;
real *ylm0, *ylmc, *ylmu, *bdr, *emu, *bem, *rmu, *sqt;
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    integer bdr_dim1, bdr_offset, rmu_dim1, rmu_offset, ylmc_dim1, 
	    ylmc_offset, ylmu_dim1, ylmu_offset, i__1, i__2, i__3;

    /* Local variables */
    static real dref, ylmg[1010]	/* was [101][10] */;
    static integer k, jg, iq, jq, iu;
    extern /* Subroutine */ int qgausn_(), errmsg_(), lepoly_(), zeroit_();
    static real sgn, gmu[10], gwt[10], sum;

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
    --bem;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    bdr -= bdr_offset;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylmc -= ylmc_offset;
    rmu_dim1 = *mxumu;
    rmu_offset = 1 + rmu_dim1 * 0;
    rmu -= rmu_offset;
    --emu;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1 * 1;
    ylmu -= ylmu_offset;
    --umu;
    --sqt;

    /* Function Body */
    if (pass1) {
	pass1 = FALSE_;
	qgausn_(&c__10, gmu, gwt);
	lepoly_(&c__10, &c__0, &c__100, &c__100, gmu, &sqt[1], ylmg);
/*                       ** Convert Legendre polys. to negative GMU */
	sgn = (float)-1.;
	for (k = 0; k <= 100; ++k) {
	    sgn = -sgn;
	    for (jg = 1; jg <= 10; ++jg) {
		ylmg[k + jg * 101 - 101] = sgn * ylmg[k + jg * 101 - 101];
/* L10: */
	    }
/* L20: */
	}
    }
    i__1 = *mi * (*mi + 1);
    zeroit_(&bdr[bdr_offset], &i__1);
    zeroit_(&bem[1], mi);
    if (*lamber && *mazim == 0) {
	i__1 = *nn;
	for (iq = 1; iq <= i__1; ++iq) {
	    bem[iq] = (float)1. - *albedo;
	    i__2 = *nn;
	    for (jq = 0; jq <= i__2; ++jq) {
		bdr[iq + jq * bdr_dim1] = *albedo;
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
		sum = (float)0.;
		i__3 = *nstr - 1;
		for (k = *mazim; k <= i__3; ++k) {
		    sum += hlpr[k] * ylmc[k + iq * ylmc_dim1] * ylmc[k + (jq 
			    + *nn) * ylmc_dim1];
/* L50: */
		}
		bdr[iq + jq * bdr_dim1] = ((float)2. - *delm0) * sum;
/* L60: */
	    }
	    if (*fbeam > (float)0.) {
		sum = (float)0.;
		i__2 = *nstr - 1;
		for (k = *mazim; k <= i__2; ++k) {
		    sum += hlpr[k] * ylmc[k + iq * ylmc_dim1] * ylm0[k];
/* L70: */
		}
		bdr[iq] = ((float)2. - *delm0) * sum;
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
		dref = (float)0.;
		for (jg = 1; jg <= 10; ++jg) {
		    sum = (float)0.;
		    i__2 = *nstr - 1;
		    for (k = 0; k <= i__2; ++k) {
			sum += hlpr[k] * ylmc[k + iq * ylmc_dim1] * ylmg[k + 
				jg * 101 - 101];
/* L90: */
		    }
		    dref += gwt[jg - 1] * (float)2. * gmu[jg - 1] * sum;
/* L100: */
		}
		bem[iq] = (float)1. - dref;
/* L110: */
	    }
	}
    }
/*                                       ** Compute surface bidirectional */
/*                                       ** properties at user angles */
    if (! (*onlyfl) && *usrang) {
	zeroit_(&emu[1], mxumu);
	i__1 = *mxumu * (*mi + 1);
	zeroit_(&rmu[rmu_offset], &i__1);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    if (umu[iu] > (float)0.) {
		if (*lamber && *mazim == 0) {
		    i__2 = *nn;
		    for (iq = 0; iq <= i__2; ++iq) {
			rmu[iu + iq * rmu_dim1] = *albedo;
/* L120: */
		    }
		    emu[iu] = (float)1. - *albedo;
		} else if (! (*lamber)) {
		    i__2 = *nn;
		    for (iq = 1; iq <= i__2; ++iq) {
			sum = (float)0.;
			i__3 = *nstr - 1;
			for (k = *mazim; k <= i__3; ++k) {
			    sum += hlpr[k] * ylmu[k + iu * ylmu_dim1] * ylmc[
				    k + (iq + *nn) * ylmc_dim1];
/* L130: */
			}
			rmu[iu + iq * rmu_dim1] = ((float)2. - *delm0) * sum;
/* L140: */
		    }
		    if (*fbeam > (float)0.) {
			sum = (float)0.;
			i__2 = *nstr - 1;
			for (k = *mazim; k <= i__2; ++k) {
			    sum += hlpr[k] * ylmu[k + iu * ylmu_dim1] * ylm0[
				    k];
/* L150: */
			}
			rmu[iu] = ((float)2. - *delm0) * sum;
		    }
		    if (*mazim == 0) {
/*                               ** Integrate bidirectional reflectivity */
/*                               ** at reflection angles UMU and */
/*                               ** incident angles GMU to get */
/*                               ** directional emissivity at */
/*                               ** user angles UMU. */
			dref = (float)0.;
			for (jg = 1; jg <= 10; ++jg) {
			    sum = (float)0.;
			    i__2 = *nstr - 1;
			    for (k = 0; k <= i__2; ++k) {
				sum += hlpr[k] * ylmu[k + iu * ylmu_dim1] * 
					ylmg[k + jg * 101 - 101];
/* L160: */
			    }
			    dref += gwt[jg - 1] * (float)2. * gmu[jg - 1] * 
				    sum;
/* L170: */
			}
			emu[iu] = (float)1. - dref;
		    }
		}
	    }
/* L180: */
	}
    }
    return 0;
} /* surfac_ */

/* Subroutine */ int terpev_(cwt, evecc, gl, gu, mazim, mxcmu, mxumu, nn, 
	nstr, numu, wk, ylmc, ylmu)
real *cwt, *evecc, *gl, *gu;
integer *mazim, *mxcmu, *mxumu, *nn, *nstr, *numu;
real *wk, *ylmc, *ylmu;
{
    /* System generated locals */
    integer evecc_dim1, evecc_offset, gu_dim1, gu_offset, ylmc_dim1, 
	    ylmc_offset, ylmu_dim1, ylmu_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer l, iq, jq, iu;
    static real sum;

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
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylmc -= ylmc_offset;
    --wk;
    evecc_dim1 = *mxcmu;
    evecc_offset = 1 + evecc_dim1 * 1;
    evecc -= evecc_offset;
    --cwt;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1 * 1;
    ylmu -= ylmu_offset;
    gu_dim1 = *mxumu;
    gu_offset = 1 + gu_dim1 * 1;
    gu -= gu_offset;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr - 1;
	for (l = *mazim; l <= i__2; ++l) {
/*                                   ** Inner sum in SD(8) times all */
/*                                   ** factors in outer sum but PLM(mu) */
	    sum = (float)0.;
	    i__3 = *nstr;
	    for (jq = 1; jq <= i__3; ++jq) {
		sum += cwt[jq] * ylmc[l + jq * ylmc_dim1] * evecc[jq + iq * 
			evecc_dim1];
/* L10: */
	    }
	    wk[l + 1] = gl[l] * (float).5 * sum;
/* L20: */
	}
/*                                    ** Finish outer sum in SD(8) */
/*                                    ** and store eigenvectors */
	i__2 = *numu;
	for (iu = 1; iu <= i__2; ++iu) {
	    sum = (float)0.;
	    i__3 = *nstr - 1;
	    for (l = *mazim; l <= i__3; ++l) {
		sum += wk[l + 1] * ylmu[l + iu * ylmu_dim1];
/* L30: */
	    }
	    if (iq <= *nn) {
		gu[iu + (iq + *nn) * gu_dim1] = sum;
	    }
	    if (iq > *nn) {
		gu[iu + (*nstr + 1 - iq) * gu_dim1] = sum;
	    }
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* terpev_ */

/* Subroutine */ int terpso_(cwt, delm0, fbeam, gl, mazim, mxcmu, plank, numu,
	 nstr, oprim, pi, ylm0, ylmc, ylmu, psi, xr0, xr1, z0, zj, zbeam, z0u,
	 z1u)
real *cwt, *delm0, *fbeam, *gl;
integer *mazim, *mxcmu;
logical *plank;
integer *numu, *nstr;
real *oprim, *pi, *ylm0, *ylmc, *ylmu, *psi, *xr0, *xr1, *z0, *zj, *zbeam, *
	z0u, *z1u;
{
    /* System generated locals */
    integer ylmc_dim1, ylmc_offset, ylmu_dim1, ylmu_offset, i__1, i__2;

    /* Local variables */
    static real fact, psum;
    static integer iq, jq, iu;
    static real sum;

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
    --zj;
    --z0;
    --psi;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1 * 1;
    ylmu -= ylmu_offset;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylmc -= ylmc_offset;
    --cwt;
    --zbeam;
    --z0u;
    --z1u;

    /* Function Body */
    if (*fbeam > (float)0.) {
/*                                  ** Beam source terms; Eq. SD(9) */
	i__1 = *nstr - 1;
	for (iq = *mazim; iq <= i__1; ++iq) {
	    psum = (float)0.;
	    i__2 = *nstr;
	    for (jq = 1; jq <= i__2; ++jq) {
		psum += cwt[jq] * ylmc[iq + jq * ylmc_dim1] * zj[jq];
/* L10: */
	    }
	    psi[iq + 1] = gl[iq] * (float).5 * psum;
/* L20: */
	}
	fact = ((float)2. - *delm0) * *fbeam / (*pi * (float)4.);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    sum = (float)0.;
	    i__2 = *nstr - 1;
	    for (iq = *mazim; iq <= i__2; ++iq) {
		sum += ylmu[iq + iu * ylmu_dim1] * (psi[iq + 1] + fact * gl[
			iq] * ylm0[iq]);
/* L30: */
	    }
	    zbeam[iu] = sum;
/* L40: */
	}
    }
    if (*plank && *mazim == 0) {
/*                                   ** Thermal source terms, STWJ(27c) */
	i__1 = *nstr - 1;
	for (iq = *mazim; iq <= i__1; ++iq) {
	    psum = (float)0.;
	    i__2 = *nstr;
	    for (jq = 1; jq <= i__2; ++jq) {
		psum += cwt[jq] * ylmc[iq + jq * ylmc_dim1] * z0[jq];
/* L50: */
	    }
	    psi[iq + 1] = gl[iq] * (float).5 * psum;
/* L60: */
	}
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    sum = (float)0.;
	    i__2 = *nstr - 1;
	    for (iq = *mazim; iq <= i__2; ++iq) {
		sum += ylmu[iq + iu * ylmu_dim1] * psi[iq + 1];
/* L70: */
	    }
	    z0u[iu] = sum + ((float)1. - *oprim) * *xr0;
	    z1u[iu] = *xr1;
/* L80: */
	}
    }
    return 0;
} /* terpso_ */

/* Subroutine */ int upbeam_(array, cc, cmu, delm0, fbeam, gl, ipvt, mazim, 
	mxcmu, nn, nstr, pi, umu0, wk, ylm0, ylmc, zj, zz)
real *array, *cc, *cmu, *delm0, *fbeam, *gl;
integer *ipvt, *mazim, *mxcmu, *nn, *nstr;
real *pi, *umu0, *wk, *ylm0, *ylmc, *zj, *zz;
{
    /* System generated locals */
    integer array_dim1, array_offset, cc_dim1, cc_offset, ylmc_dim1, 
	    ylmc_offset, i__1, i__2;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int sgeco_();
    static real rcond;
    extern /* Subroutine */ int sgesl_();
    static integer iq, jq;
    extern /* Subroutine */ int errmsg_();
    static integer job;
    static real sum;

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
    --ipvt;
    --zz;
    --zj;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylmc -= ylmc_offset;
    --wk;
    --cmu;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1 * 1;
    cc -= cc_offset;
    array_dim1 = *mxcmu;
    array_offset = 1 + array_dim1 * 1;
    array -= array_offset;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    array[iq + jq * array_dim1] = -cc[iq + jq * cc_dim1];
/* L10: */
	}
	array[iq + iq * array_dim1] = cmu[iq] / *umu0 + (float)1. + array[iq 
		+ iq * array_dim1];
	sum = (float)0.;
	i__2 = *nstr - 1;
	for (k = *mazim; k <= i__2; ++k) {
	    sum += gl[k] * ylmc[k + iq * ylmc_dim1] * ylm0[k];
/* L20: */
	}
	zj[iq] = ((float)2. - *delm0) * *fbeam * sum / (*pi * (float)4.);
/* L30: */
    }
/*                  ** Find L-U (lower/upper triangular) decomposition */
/*                  ** of ARRAY and see if it is nearly singular */
/*                  ** (NOTE:  ARRAY is altered) */
    rcond = (float)0.;
    sgeco_(&array[array_offset], mxcmu, nstr, &ipvt[1], &rcond, &wk[1]);
    if (rcond + (float)1. == (float)1.) {
	errmsg_("UPBEAM--SGECO says matrix near singular", &c_false, (ftnlen)
		39);
    }
/*                ** Solve linear system with coeff matrix ARRAY */
/*                ** (assumed already L-U decomposed) and R.H. side(s) */
/*                ** ZJ;  return solution(s) in ZJ */
    job = 0;
    sgesl_(&array[array_offset], mxcmu, nstr, &ipvt[1], &zj[1], &job);
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zz[iq + *nn] = zj[iq];
	zz[*nn + 1 - iq] = zj[iq + *nn];
/* L40: */
    }
    return 0;
} /* upbeam_ */

/* Subroutine */ int upisot_(array, cc, cmu, ipvt, mxcmu, nn, nstr, oprim, wk,
	 xr0, xr1, z0, z1, zplk0, zplk1)
real *array, *cc, *cmu;
integer *ipvt, *mxcmu, *nn, *nstr;
real *oprim, *wk, *xr0, *xr1, *z0, *z1, *zplk0, *zplk1;
{
    /* System generated locals */
    integer array_dim1, array_offset, cc_dim1, cc_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int sgeco_();
    static real rcond;
    extern /* Subroutine */ int sgesl_();
    static integer iq, jq;
    extern /* Subroutine */ int errmsg_();

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
    --ipvt;
    --zplk1;
    --zplk0;
    --z1;
    --z0;
    --wk;
    --cmu;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1 * 1;
    cc -= cc_offset;
    array_dim1 = *mxcmu;
    array_offset = 1 + array_dim1 * 1;
    array -= array_offset;

    /* Function Body */
    i__1 = *nstr;
    for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *nstr;
	for (jq = 1; jq <= i__2; ++jq) {
	    array[iq + jq * array_dim1] = -cc[iq + jq * cc_dim1];
/* L10: */
	}
	array[iq + iq * array_dim1] += (float)1.;
	z1[iq] = *xr1;
	z0[iq] = ((float)1. - *oprim) * *xr0 + cmu[iq] * z1[iq];
/* L20: */
    }
/*                       ** Solve linear equations: same as in UPBEAM, */
/*                       ** except ZJ replaced by Z0 */
    rcond = (float)0.;
    sgeco_(&array[array_offset], mxcmu, nstr, &ipvt[1], &rcond, &wk[1]);
    if (rcond + (float)1. == (float)1.) {
	errmsg_("UPISOT--SGECO says matrix near singular", &c_false, (ftnlen)
		39);
    }
    sgesl_(&array[array_offset], mxcmu, nstr, &ipvt[1], &z0[1], &c__0);
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zplk0[iq + *nn] = z0[iq];
	zplk1[iq + *nn] = z1[iq];
	zplk0[*nn + 1 - iq] = z0[iq + *nn];
	zplk1[*nn + 1 - iq] = z1[iq + *nn];
/* L30: */
    }
    return 0;
} /* upisot_ */

/* Subroutine */ int usrint_(bplank, cmu, cwt, delm0, dtaucp, emu, expbea, 
	fbeam, fisot, gc, gu, kk, lamber, layru, ll, lyrcut, mazim, mxcmu, 
	mxulv, mxumu, ncut, nlyr, nn, nstr, plank, numu, ntau, pi, rmu, 
	taucpr, tplank, umu, umu0, utaupr, wk, zbeam, z0u, z1u, zz, zplk0, 
	zplk1, uum)
real *bplank, *cmu, *cwt, *delm0, *dtaucp, *emu, *expbea, *fbeam, *fisot, *gc,
	 *gu, *kk;
logical *lamber;
integer *layru;
real *ll;
logical *lyrcut;
integer *mazim, *mxcmu, *mxulv, *mxumu, *ncut, *nlyr, *nn, *nstr;
logical *plank;
integer *numu, *ntau;
real *pi, *rmu, *taucpr, *tplank, *umu, *umu0, *utaupr, *wk, *zbeam, *z0u, *
	z1u, *zz, *zplk0, *zplk1, *uum;
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, gu_dim1, gu_dim2, gu_offset, kk_dim1,
	     kk_offset, ll_dim1, ll_offset, rmu_dim1, rmu_offset, uum_dim1, 
	    uum_offset, z0u_dim1, z0u_offset, z1u_dim1, z1u_offset, 
	    zbeam_dim1, zbeam_offset, zplk0_dim1, zplk0_offset, zplk1_dim1, 
	    zplk1_offset, zz_dim1, zz_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static real fact, dtau, expn, dtau1, dtau2, denom;
    static integer lc, iq, jq, iu, lu;
    static real bnddfu, bnddir, bndint, palint, dfuint;
    static integer lyrend;
    static logical negumu;
    static real plkint;
    static integer lyrstr;
    static real sgn;
    static integer lyu;
    static real exp0, exp1, exp2;

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
    --dtaucp;
    --layru;
    zplk1_dim1 = *mxcmu;
    zplk1_offset = 1 + zplk1_dim1 * 1;
    zplk1 -= zplk1_offset;
    zplk0_dim1 = *mxcmu;
    zplk0_offset = 1 + zplk0_dim1 * 1;
    zplk0 -= zplk0_offset;
    zz_dim1 = *mxcmu;
    zz_offset = 1 + zz_dim1 * 1;
    zz -= zz_offset;
    --wk;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    gc -= gc_offset;
    --cwt;
    --cmu;
    --utaupr;
    uum_dim1 = *mxumu;
    uum_offset = 1 + uum_dim1 * 1;
    uum -= uum_offset;
    z1u_dim1 = *mxumu;
    z1u_offset = 1 + z1u_dim1 * 1;
    z1u -= z1u_offset;
    z0u_dim1 = *mxumu;
    z0u_offset = 1 + z0u_dim1 * 1;
    z0u -= z0u_offset;
    zbeam_dim1 = *mxumu;
    zbeam_offset = 1 + zbeam_dim1 * 1;
    zbeam -= zbeam_offset;
    rmu_dim1 = *mxumu;
    rmu_offset = 1 + rmu_dim1 * 0;
    rmu -= rmu_offset;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2 * 1);
    gu -= gu_offset;
    --emu;
    --umu;

    /* Function Body */
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	i__2 = *nstr;
	for (iq = 1; iq <= i__2; ++iq) {
	    i__3 = *numu;
	    for (iu = 1; iu <= i__3; ++iu) {
		gu[iu + (iq + lc * gu_dim2) * gu_dim1] *= ll[iq + lc * 
			ll_dim1];
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
	if (*fbeam > (float)0.) {
	    exp0 = exp(-utaupr[lu] / *umu0);
	}
	lyu = layru[lu];
/*                              ** Loop over polar angles at which */
/*                              ** intensities are desired */
	i__2 = *numu;
	for (iu = 1; iu <= i__2; ++iu) {
	    if (*lyrcut && lyu > *ncut) {
		goto L150;
	    }
	    negumu = umu[iu] < (float)0.;
	    if (negumu) {
		lyrstr = 1;
		lyrend = lyu - 1;
		sgn = (float)-1.;
	    } else {
		lyrstr = lyu + 1;
		lyrend = *ncut;
		sgn = (float)1.;
	    }
/*                          ** For downward intensity, integrate from top */
/*                          ** to LYU-1 in Eq. S1(8); for upward, */
/*                          ** integrate from bottom to LYU+1 in S1(9) */
	    palint = (float)0.;
	    plkint = (float)0.;
	    i__3 = lyrend;
	    for (lc = lyrstr; lc <= i__3; ++lc) {
		dtau = dtaucp[lc];
		exp1 = exp((utaupr[lu] - taucpr[lc - 1]) / umu[iu]);
		exp2 = exp((utaupr[lu] - taucpr[lc]) / umu[iu]);
		if (*plank && *mazim == 0) {
		    plkint += sgn * (z0u[iu + lc * z0u_dim1] * (exp1 - exp2) 
			    + z1u[iu + lc * z1u_dim1] * ((taucpr[lc - 1] + 
			    umu[iu]) * exp1 - (taucpr[lc] + umu[iu]) * exp2));
		}
		if (*fbeam > (float)0.) {
		    denom = umu[iu] / *umu0 + (float)1.;
		    if (dabs(denom) < (float)1e-4) {
/*                                                   ** L'Hospital limit */
			expn = dtau / *umu0 * exp0;
		    } else {
			expn = (exp1 * expbea[lc - 1] - exp2 * expbea[lc]) * 
				sgn / denom;
		    }
		    palint += zbeam[iu + lc * zbeam_dim1] * expn;
		}
/*                                                   ** KK is negative */
		i__4 = *nn;
		for (iq = 1; iq <= i__4; ++iq) {
		    wk[iq] = exp(kk[iq + lc * kk_dim1] * dtau);
		    denom = umu[iu] * kk[iq + lc * kk_dim1] + (float)1.;
		    if (dabs(denom) < (float)1e-4) {
/*                                                   ** L'Hospital limit */
			expn = dtau / umu[iu] * exp2;
		    } else {
			expn = sgn * (exp1 * wk[iq] - exp2) / denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1] * expn;
/* L40: */
		}
/*                                                   ** KK is positive */
		i__4 = *nstr;
		for (iq = *nn + 1; iq <= i__4; ++iq) {
		    denom = umu[iu] * kk[iq + lc * kk_dim1] + (float)1.;
		    if (dabs(denom) < (float)1e-4) {
/*                                                   ** L'Hospital limit */
			expn = -dtau / umu[iu] * exp1;
		    } else {
			expn = sgn * (exp1 - exp2 * wk[*nstr + 1 - iq]) / 
				denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1] * expn;
/* L50: */
		}
/* L60: */
	    }
/*                           ** Calculate contribution from user */
/*                           ** output level to next computational level */
	    dtau1 = utaupr[lu] - taucpr[lyu - 1];
	    dtau2 = utaupr[lu] - taucpr[lyu];
	    if (dabs(dtau1) < (float)1e-6 && negumu) {
		goto L90;
	    }
	    if (dabs(dtau2) < (float)1e-6 && ! negumu) {
		goto L90;
	    }
	    if (negumu) {
		exp1 = exp(dtau1 / umu[iu]);
	    }
	    if (! negumu) {
		exp2 = exp(dtau2 / umu[iu]);
	    }
	    if (*fbeam > (float)0.) {
		denom = umu[iu] / *umu0 + (float)1.;
		if (dabs(denom) < (float)1e-4) {
		    expn = dtau1 / *umu0 * exp0;
		} else if (negumu) {
		    expn = (exp0 - expbea[lyu - 1] * exp1) / denom;
		} else {
		    expn = (exp0 - expbea[lyu] * exp2) / denom;
		}
		palint += zbeam[iu + lyu * zbeam_dim1] * expn;
	    }
/*                                                   ** KK is negative */
	    dtau = dtaucp[lyu];
	    i__3 = *nn;
	    for (iq = 1; iq <= i__3; ++iq) {
		denom = umu[iu] * kk[iq + lyu * kk_dim1] + (float)1.;
		if (dabs(denom) < (float)1e-4) {
		    expn = -dtau2 / umu[iu] * exp2;
		} else if (negumu) {
		    expn = (exp(-kk[iq + lyu * kk_dim1] * dtau2) - exp(kk[iq 
			    + lyu * kk_dim1] * dtau) * exp1) / denom;
		} else {
		    expn = (exp(-kk[iq + lyu * kk_dim1] * dtau2) - exp2) / 
			    denom;
		}
		palint += gu[iu + (iq + lyu * gu_dim2) * gu_dim1] * expn;
/* L70: */
	    }
/*                                                   ** KK is positive */
	    i__3 = *nstr;
	    for (iq = *nn + 1; iq <= i__3; ++iq) {
		denom = umu[iu] * kk[iq + lyu * kk_dim1] + (float)1.;
		if (dabs(denom) < (float)1e-4) {
		    expn = -dtau1 / umu[iu] * exp1;
		} else if (negumu) {
		    expn = (exp(-kk[iq + lyu * kk_dim1] * dtau1) - exp1) / 
			    denom;
		} else {
		    expn = (exp(-kk[iq + lyu * kk_dim1] * dtau1) - exp(-kk[iq 
			    + lyu * kk_dim1] * dtau) * exp2) / denom;
		}
		palint += gu[iu + (iq + lyu * gu_dim2) * gu_dim1] * expn;
/* L80: */
	    }
	    if (*plank && *mazim == 0) {
		if (negumu) {
		    expn = exp1;
		    fact = taucpr[lyu - 1] + umu[iu];
		} else {
		    expn = exp2;
		    fact = taucpr[lyu] + umu[iu];
		}
		plkint = plkint + z0u[iu + lyu * z0u_dim1] * ((float)1. - 
			expn) + z1u[iu + lyu * z1u_dim1] * (utaupr[lu] + umu[
			iu] - fact * expn);
	    }
/*                            ** Calculate intensity components */
/*                            ** attenuated at both boundaries. */
/*                            ** NOTE: no azimuthal intensity */
/*                            ** component for isotropic surface */
L90:
	    bndint = (float)0.;
	    if (negumu && *mazim == 0) {
		bndint = (*fisot + *tplank) * exp(utaupr[lu] / umu[iu]);
	    } else if (! negumu) {
		if (*lyrcut || *lamber && *mazim > 0) {
		    goto L140;
		}
		i__3 = *nstr;
		for (jq = *nn + 1; jq <= i__3; ++jq) {
		    wk[jq] = exp(-kk[jq + *nlyr * kk_dim1] * dtaucp[*nlyr]);
/* L100: */
		}
		bnddfu = (float)0.;
		for (iq = *nn; iq >= 1; --iq) {
		    dfuint = (float)0.;
		    i__3 = *nn;
		    for (jq = 1; jq <= i__3; ++jq) {
			dfuint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1] * 
				ll[jq + *nlyr * ll_dim1];
/* L110: */
		    }
		    i__3 = *nstr;
		    for (jq = *nn + 1; jq <= i__3; ++jq) {
			dfuint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1] * 
				ll[jq + *nlyr * ll_dim1] * wk[jq];
/* L120: */
		    }
		    if (*fbeam > (float)0.) {
			dfuint += zz[iq + *nlyr * zz_dim1] * expbea[*nlyr];
		    }
		    dfuint += *delm0 * (zplk0[iq + *nlyr * zplk0_dim1] + 
			    zplk1[iq + *nlyr * zplk1_dim1] * taucpr[*nlyr]);
		    bnddfu += (*delm0 + (float)1.) * rmu[iu + (*nn + 1 - iq) *
			     rmu_dim1] * cmu[*nn + 1 - iq] * cwt[*nn + 1 - iq]
			     * dfuint;
/* L130: */
		}
		bnddir = (float)0.;
		if (*fbeam > (float)0.) {
		    bnddir = *umu0 * *fbeam / *pi * rmu[iu] * expbea[*nlyr];
		}
		bndint = (bnddfu + bnddir + *delm0 * emu[iu] * *bplank) * exp(
			(utaupr[lu] - taucpr[*nlyr]) / umu[iu]);
	    }
L140:
	    uum[iu + lu * uum_dim1] = palint + plkint + bndint;
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
/* Subroutine */ int chekin_(nlyr, dtauc, ssalb, pmom, temper, wvnmlo, wvnmhi,
	 usrtau, ntau, utau, nstr, usrang, numu, umu, nphi, phi, ibcnd, fbeam,
	 umu0, phi0, fisot, lamber, albedo, hl, btemp, ttemp, temis, plank, 
	onlyfl, accur, tauc, maxcly, maxulv, maxumu, maxcmu, maxphi, mxcly, 
	mxulv, mxumu, mxcmu, mxphi, mxsqt)
integer *nlyr;
real *dtauc, *ssalb, *pmom, *temper, *wvnmlo, *wvnmhi;
logical *usrtau;
integer *ntau;
real *utau;
integer *nstr;
logical *usrang;
integer *numu;
real *umu;
integer *nphi;
real *phi;
integer *ibcnd;
real *fbeam, *umu0, *phi0, *fisot;
logical *lamber;
real *albedo, *hl, *btemp, *ttemp, *temis;
logical *plank, *onlyfl;
real *accur, *tauc;
integer *maxcly, *maxulv, *maxumu, *maxcmu, *maxphi, *mxcly, *mxulv, *mxumu, *
	mxcmu, *mxphi, *mxsqt;
{
    /* System generated locals */
    integer pmom_dim1, pmom_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    extern doublereal dref_();
    static integer irmu, j, k, lc, iu, lu;
    static real flxalb;
    extern logical wrtbad_();
    static logical inperr;
    extern /* Subroutine */ int errmsg_();
    extern logical wrtdim_();
    static integer numsqt;
    static real rmu;

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
    --ssalb;
    --dtauc;
    --utau;
    --umu;
    pmom_dim1 = *maxcmu - 0 + 1;
    pmom_offset = 0 + pmom_dim1 * 1;
    pmom -= pmom_offset;
    --phi;

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
	if (dtauc[lc] < (float)0.) {
	    inperr = wrtbad_("DTAUC", (ftnlen)5);
	}
	if (ssalb[lc] < (float)0. || ssalb[lc] > (float)1.) {
	    inperr = wrtbad_("SSALB", (ftnlen)5);
	}
	if (*plank && *ibcnd != 1) {
	    if (lc == 1 && temper[0] < (float)0.) {
		inperr = wrtbad_("TEMPER", (ftnlen)6);
	    }
	    if (temper[lc] < (float)0.) {
		inperr = wrtbad_("TEMPER", (ftnlen)6);
	    }
	}
	i__2 = *nstr;
	for (k = 0; k <= i__2; ++k) {
	    if (pmom[k + lc * pmom_dim1] < (float)-1. || pmom[k + lc * 
		    pmom_dim1] > (float)1.) {
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
	    if ((r__1 = utau[lu] - tauc[*nlyr], dabs(r__1)) <= (float)1e-4) {
		utau[lu] = tauc[*nlyr];
	    }
	    if (utau[lu] < (float)0. || utau[lu] > tauc[*nlyr]) {
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
	    if (umu[iu] < (float)-1. || umu[iu] > (float)1. || umu[iu] == (
		    float)0.) {
		inperr = wrtbad_("UMU", (ftnlen)3);
	    }
	    if (*ibcnd == 1 && umu[iu] < (float)0.) {
		inperr = wrtbad_("UMU", (ftnlen)3);
	    }
	    if (iu > 1) {
		if (umu[iu] < umu[iu - 1]) {
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
	    if (phi[j] < (float)0. || phi[j] > (float)360.) {
		inperr = wrtbad_("PHI", (ftnlen)3);
	    }
/* L50: */
	}
    }
    if (*ibcnd < 0 || *ibcnd > 1) {
	inperr = wrtbad_("IBCND", (ftnlen)5);
    }
    if (*ibcnd == 0) {
	if (*fbeam < (float)0.) {
	    inperr = wrtbad_("FBEAM", (ftnlen)5);
	}
	if (*fbeam > (float)0. && (*umu0 <= (float)0. || *umu0 > (float)1.)) {
	    inperr = wrtbad_("UMU0", (ftnlen)4);
	}
	if (*fbeam > (float)0. && (*phi0 < (float)0. || *phi0 > (float)360.)) 
		{
	    inperr = wrtbad_("PHI0", (ftnlen)4);
	}
	if (*fisot < (float)0.) {
	    inperr = wrtbad_("FISOT", (ftnlen)5);
	}
	if (*lamber) {
	    if (*albedo < (float)0. || *albedo > (float)1.) {
		inperr = wrtbad_("ALBEDO", (ftnlen)6);
	    }
	} else {
/*                    ** Make sure flux albedo at dense mesh of incident */
/*                    ** angles does not assume unphysical values */
	    for (irmu = 0; irmu <= 100; ++irmu) {
		rmu = irmu * (float).01;
		flxalb = dref_(&rmu, hl, nstr);
		if (flxalb < (float)0. || flxalb > (float)1.) {
		    inperr = wrtbad_("HL", (ftnlen)2);
		}
/* L60: */
	    }
	}
    } else if (*ibcnd == 1) {
	if (*albedo < (float)0. || *albedo > (float)1.) {
	    inperr = wrtbad_("ALBEDO", (ftnlen)6);
	}
    }
    if (*plank && *ibcnd != 1) {
	if (*wvnmlo < (float)0. || *wvnmhi <= *wvnmlo) {
	    inperr = wrtbad_("WVNMLO,HI", (ftnlen)9);
	}
	if (*temis < (float)0. || *temis > (float)1.) {
	    inperr = wrtbad_("TEMIS", (ftnlen)5);
	}
	if (*btemp < (float)0.) {
	    inperr = wrtbad_("BTEMP", (ftnlen)5);
	}
	if (*ttemp < (float)0.) {
	    inperr = wrtbad_("TTEMP", (ftnlen)5);
	}
    }
    if (*accur < (float)0. || *accur > (float).01) {
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
	    if ((r__1 = temper[lc] - temper[lc - 1], dabs(r__1)) > (float)20.)
		     {
		errmsg_("CHEKIN--vertical temperature step may be too large \
for good accuracy", &c_false, (ftnlen)68);
	    }
/* L70: */
	}
    }
    return 0;
} /* chekin_ */

doublereal dref_(mu, hl, nstr)
real *mu, *hl;
integer *nstr;
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static real c__[100];
    static integer l;
    static real cl, pl;
    extern /* Subroutine */ int errmsg_();
    static real plm1, plm2;

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
/*     .. */
    if (pass1) {
	pass1 = FALSE_;
	cl = (float).125;
	c__[1] = cl * (float)10.;
	for (l = 4; l <= 100; l += 2) {
	    cl = -cl * (l - 3) / (l + 2);
	    c__[l - 1] = ((l << 1) + 1) * (float)2. * cl;
/* L10: */
	}
    }
    if (*nstr < 2 || dabs(*mu) > (float)1.) {
	errmsg_("DREF--input argument error(s)", &c_true, (ftnlen)29);
    }
    if (*nstr > 100) {
	errmsg_("DREF--parameter MAXTRM too small", &c_true, (ftnlen)32);
    }
    ret_val = hl[0] - hl[1] * (float)2. * *mu;
    plm2 = (float)1.;
    plm1 = -(*mu);
    i__1 = *nstr - 1;
    for (l = 2; l <= i__1; ++l) {
/*                                ** Legendre polynomial recurrence */
	pl = (((l << 1) - 1) * (-(*mu)) * plm1 - (l - 1) * plm2) / l;
	if (l % 2 == 0) {
	    ret_val += c__[l - 1] * hl[l] * pl;
	}
	plm2 = plm1;
	plm1 = pl;
/* L20: */
    }
    if (ret_val < (float)0. || ret_val > (float)1.) {
	errmsg_("DREF--albedo value not in (0,1)", &c_false, (ftnlen)31);
    }
    return ret_val;
} /* dref_ */

/* Subroutine */ int lepoly_(nmu, m, maxmu, twonm1, mu, sqt, ylm)
integer *nmu, *m, *maxmu, *twonm1;
real *mu, *sqt, *ylm;
{
    /* System generated locals */
    integer ylm_dim1, ylm_offset, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__, l;
    static real tmp1, tmp2;

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
    ylm -= ylm_offset;
    --mu;
    --sqt;

    /* Function Body */
    if (*m == 0) {
/*                             ** Upward recurrence for ordinary */
/*                             ** Legendre polynomials */
	i__1 = *nmu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ylm[i__ * ylm_dim1] = (float)1.;
	    ylm[i__ * ylm_dim1 + 1] = mu[i__];
/* L20: */
	}
	i__1 = *twonm1;
	for (l = 2; l <= i__1; ++l) {
	    i__2 = *nmu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ylm[l + i__ * ylm_dim1] = (((l << 1) - 1) * mu[i__] * ylm[l - 
			1 + i__ * ylm_dim1] - (l - 1) * ylm[l - 2 + i__ * 
			ylm_dim1]) / l;
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
	    r__1 = mu[i__];
	    ylm[*m + i__ * ylm_dim1] = -sqt[(*m << 1) - 1] / sqt[*m * 2] * 
		    sqrt((float)1. - r__1 * r__1) * ylm[*m - 1 + i__ * 
		    ylm_dim1];
/*                              ** Y-sub-(m+1)-super-m; derived from */
/*                              ** D/A Eqs.(13,14) using Eqs.(11,12) */
	    ylm[*m + 1 + i__ * ylm_dim1] = sqt[(*m << 1) + 1] * mu[i__] * ylm[
		    *m + i__ * ylm_dim1];
/* L50: */
	}
/*                                   ** Upward recurrence; D/A EQ.(10) */
	i__1 = *twonm1;
	for (l = *m + 2; l <= i__1; ++l) {
	    tmp1 = sqt[l - *m] * sqt[l + *m];
	    tmp2 = sqt[l - *m - 1] * sqt[l + *m - 1];
	    i__2 = *nmu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ylm[l + i__ * ylm_dim1] = (((l << 1) - 1) * mu[i__] * ylm[l - 
			1 + i__ * ylm_dim1] - tmp2 * ylm[l - 2 + i__ * 
			ylm_dim1]) / tmp1;
/* L60: */
	    }
/* L70: */
	}
    }
    return 0;
} /* lepoly_ */

doublereal plkavg_(wnumlo, wnumhi, t)
real *wnumlo, *wnumhi, *t;
{
    /* Initialized data */

    static real c2 = (float)1.438786;
    static real sigma = (float)5.67032e-8;
    static real vcut = (float)1.5;
    static real vcp[7] = { (float)10.25,(float)5.7,(float)3.9,(float)2.9,(
	    float)2.3,(float)1.9,(float)0. };
    static real pi = (float)0.;

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double asin(), log(), exp();

    /* Local variables */
    static real conc;
    static integer mmax;
    static real vmax, d__[2];
    static integer i__, k, m, n;
    static real p[2], v[2], epsil;
    extern doublereal r1mach_();
    static real hh, ex, mv, sigdpi, oldval;
    static integer smallv;
    extern /* Subroutine */ int errmsg_();
    static real del, val, exm, vsq, val0;

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
    if (pi == (float)0.) {
	pi = asin((float)1.) * (float)2.;
	vmax = log(r1mach_(&c__2));
	epsil = r1mach_(&c__4);
	sigdpi = sigma / pi;
/* Computing 4th power */
	r__1 = pi, r__1 *= r__1;
	conc = (float)15. / (r__1 * r__1);
    }
    if (*t < (float)0. || *wnumhi <= *wnumlo || *wnumlo < (float)0.) {
	errmsg_("PLKAVG--temperature or wavenums. wrong", &c_true, (ftnlen)38)
		;
    }
    if (*t < (float)1e-4) {
	ret_val = (float)0.;
	return ret_val;
    }
    v[0] = c2 * *wnumlo / *t;
    v[1] = c2 * *wnumhi / *t;
    if (v[0] > epsil && v[1] < vmax && (*wnumhi - *wnumlo) / *wnumhi < (float)
	    .01) {
/*                          ** Wavenumbers are very close.  Get integral */
/*                          ** by iterating Simpson rule to convergence. */
	hh = v[1] - v[0];
	oldval = (float)0.;
/* Computing 3rd power */
	r__1 = v[0];
/* Computing 3rd power */
	r__2 = v[1];
	val0 = r__1 * (r__1 * r__1) / (exp(v[0]) - 1) + r__2 * (r__2 * r__2) /
		 (exp(v[1]) - 1);
	for (n = 1; n <= 10; ++n) {
	    del = hh / (n << 1);
	    val = val0;
	    i__1 = (n << 1) - 1;
	    for (k = 1; k <= i__1; ++k) {
		r__1 = v[0] + k * del;
/* Computing 3rd power */
		r__2 = r__1;
		val += (k % 2 + 1 << 1) * (r__2 * (r__2 * r__2) / (exp(r__1) 
			- 1));
/* L10: */
	    }
	    val = del / (float)3. * val;
	    if ((r__1 = (val - oldval) / val, dabs(r__1)) <= (float)1e-6) {
		goto L30;
	    }
	    oldval = val;
/* L20: */
	}
	errmsg_("PLKAVG--Simpson rule didnt converge", &c_false, (ftnlen)35);
L30:
/* Computing 4th power */
	r__1 = *t, r__1 *= r__1;
	ret_val = sigdpi * (r__1 * r__1) * conc * val;
	return ret_val;
    }
/*                          *** General case *** */
    smallv = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	if (v[i__ - 1] < vcut) {
/*                                   ** Use power series */
	    ++smallv;
/* Computing 2nd power */
	    r__1 = v[i__ - 1];
	    vsq = r__1 * r__1;
	    p[i__ - 1] = conc * vsq * v[i__ - 1] * (v[i__ - 1] * (v[i__ - 1] *
		     (vsq * (vsq * (vsq * (float)-7.5156325156325161e-8 + (
		    float)3.6743092298647855e-6) - (float)
		    1.9841269841269841e-4) + (float).016666666666666666) - (
		    float).125) + (float).33333333333333331);
	} else {
/*                      ** Use exponential series */
	    mmax = 0;
/*                                ** Find upper limit of series */
L40:
	    ++mmax;
	    if (v[i__ - 1] < vcp[mmax - 1]) {
		goto L40;
	    }
	    ex = exp(-v[i__ - 1]);
	    exm = (float)1.;
	    d__[i__ - 1] = (float)0.;
	    i__1 = mmax;
	    for (m = 1; m <= i__1; ++m) {
		mv = m * v[i__ - 1];
		exm = ex * exm;
/* Computing 4th power */
		i__2 = m, i__2 *= i__2;
		d__[i__ - 1] += exm * (mv * (mv * (mv + (float)3.) + (float)
			6.) + (float)6.) / (i__2 * i__2);
/* L50: */
	    }
	    d__[i__ - 1] = conc * d__[i__ - 1];
	}
/* L60: */
    }
/*                              ** Handle ill-conditioning */
    if (smallv == 2) {
/*                                    ** WNUMLO and WNUMHI both small */
	ret_val = p[1] - p[0];
    } else if (smallv == 1) {
/*                                    ** WNUMLO small, WNUMHI large */
	ret_val = (float)1. - p[0] - d__[1];
    } else {
/*                                    ** WNUMLO and WNUMHI both large */
	ret_val = d__[0] - d__[1];
    }
/* Computing 4th power */
    r__1 = *t, r__1 *= r__1;
    ret_val = sigdpi * (r__1 * r__1) * ret_val;
    if (ret_val == (float)0.) {
	errmsg_("PLKAVG--returns zero; possible underflow", &c_false, (ftnlen)
		40);
    }
    return ret_val;
} /* plkavg_ */

/* Subroutine */ int pravin_(umu, numu, maxumu, utau, ntau, u0u)
real *umu;
integer *numu, *maxumu;
real *utau;
integer *ntau;
real *u0u;
{
    /* System generated locals */
    integer u0u_dim1, u0u_offset, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer iumin, iumax, npass, iu, np, lu, lenfmt;

    /* Fortran I/O blocks */
    static cilist io___300 = { 0, 6, 0, "(//,A)", 0 };
    static cilist io___303 = { 0, 6, 0, "(/,A,/,A)", 0 };
    static cilist io___307 = { 0, 6, 0, "(/,10X,8F14.5)", 0 };
    static cilist io___310 = { 0, 6, 0, "(0P,F10.4,1P,8E14.4)", 0 };


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
    --umu;
    u0u_dim1 = *maxumu;
    u0u_offset = 1 + u0u_dim1 * 1;
    u0u -= u0u_offset;
    --utau;

    /* Function Body */
    if (*numu < 1) {
	return 0;
    }
    s_wsfe(&io___300);
    do_fio(&c__1, " *******  AZIMUTHALLY AVERAGED INTENSITIES (at user polar\
 angles)  ********", (ftnlen)75);
    e_wsfe();
    lenfmt = 8;
    npass = (*numu - 1) / lenfmt + 1;
    s_wsfe(&io___303);
    do_fio(&c__1, "   Optical   Polar Angle Cosines", (ftnlen)32);
    do_fio(&c__1, "     Depth", (ftnlen)10);
    e_wsfe();
    i__1 = npass;
    for (np = 1; np <= i__1; ++np) {
	iumin = lenfmt * (np - 1) + 1;
/* Computing MIN */
	i__2 = lenfmt * np;
	iumax = min(i__2,*numu);
	s_wsfe(&io___307);
	i__2 = iumax;
	for (iu = iumin; iu <= i__2; ++iu) {
	    do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
	}
	e_wsfe();
	i__2 = *ntau;
	for (lu = 1; lu <= i__2; ++lu) {
	    s_wsfe(&io___310);
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	    i__3 = iumax;
	    for (iu = iumin; iu <= i__3; ++iu) {
		do_fio(&c__1, (char *)&u0u[iu + lu * u0u_dim1], (ftnlen)
			sizeof(real));
	    }
	    e_wsfe();
/* L10: */
	}
/* L20: */
    }
    return 0;
} /* pravin_ */

/* Subroutine */ int prtinp_(nlyr, dtauc, dtaucp, ssalb, pmom, temper, wvnmlo,
	 wvnmhi, ntau, utau, nstr, numu, umu, nphi, phi, ibcnd, fbeam, umu0, 
	phi0, fisot, lamber, albedo, hl, btemp, ttemp, temis, deltam, plank, 
	onlyfl, accur, flyr, lyrcut, oprim, tauc, taucpr, maxcmu, prtmom)
integer *nlyr;
real *dtauc, *dtaucp, *ssalb, *pmom, *temper, *wvnmlo, *wvnmhi;
integer *ntau;
real *utau;
integer *nstr, *numu;
real *umu;
integer *nphi;
real *phi;
integer *ibcnd;
real *fbeam, *umu0, *phi0, *fisot;
logical *lamber;
real *albedo, *hl, *btemp, *ttemp, *temis;
logical *deltam, *plank, *onlyfl;
real *accur, *flyr;
logical *lyrcut;
real *oprim, *tauc, *taucpr;
integer *maxcmu;
logical *prtmom;
{
    /* System generated locals */
    integer pmom_dim1, pmom_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer j, k, lc, iu, lu;
    static real yessct;

    /* Fortran I/O blocks */
    static cilist io___311 = { 0, 6, 0, "(/,A,I4,A,I4)", 0 };
    static cilist io___312 = { 0, 6, 0, "(I4,A,10F10.4,/,(26X,10F10.4))", 0 };
    static cilist io___314 = { 0, 6, 0, "(I4,A,10F9.5,/,(31X,10F9.5))", 0 };
    static cilist io___316 = { 0, 6, 0, "(I4,A,10F9.2,/,(28X,10F9.2))", 0 };
    static cilist io___318 = { 0, 6, 0, "(A)", 0 };
    static cilist io___319 = { 0, 6, 0, "(A,I2)", 0 };
    static cilist io___320 = { 0, 6, 0, "(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P\
,E11.3)", 0 };
    static cilist io___321 = { 0, 6, 0, "(A,0P,F8.4)", 0 };
    static cilist io___322 = { 0, 6, 0, "(A,/,(10X,10F9.5))", 0 };
    static cilist io___324 = { 0, 6, 0, "(A,2F14.4,/,A,F10.2,A,F10.2,A,F8.4)",
	     0 };
    static cilist io___325 = { 0, 6, 0, "(A)", 0 };
    static cilist io___326 = { 0, 6, 0, "(A,0P,F8.4)", 0 };
    static cilist io___327 = { 0, 6, 0, "(A)", 0 };
    static cilist io___328 = { 0, 6, 0, "(A)", 0 };
    static cilist io___329 = { 0, 6, 0, "(A)", 0 };
    static cilist io___330 = { 0, 6, 0, "(A)", 0 };
    static cilist io___331 = { 0, 6, 0, "(A)", 0 };
    static cilist io___332 = { 0, 6, 0, "(A,1P,E11.2)", 0 };
    static cilist io___333 = { 0, 6, 0, "(A)", 0 };
    static cilist io___334 = { 0, 6, 0, "(/,37X,A,3(/,2A))", 0 };
    static cilist io___335 = { 0, 6, 0, "(/,37X,A,3(/,2A))", 0 };
    static cilist io___338 = { 0, 6, 0, "(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5\
,F9.4,F14.3)", 0 };
    static cilist io___339 = { 0, 6, 0, "(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5\
,F9.4)", 0 };
    static cilist io___340 = { 0, 6, 0, "(85X,F14.3)", 0 };
    static cilist io___341 = { 0, 6, 0, "(/,A)", 0 };
    static cilist io___342 = { 0, 6, 0, "(I6,10F11.6,/,(6X,10F11.6))", 0 };


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
    --dtauc;
    --dtaucp;
    --ssalb;
    --utau;
    --umu;
    --phi;
    --flyr;
    --oprim;
    pmom_dim1 = *maxcmu - 0 + 1;
    pmom_offset = 0 + pmom_dim1 * 1;
    pmom -= pmom_offset;

    /* Function Body */
    s_wsfe(&io___311);
    do_fio(&c__1, " No. streams =", (ftnlen)14);
    do_fio(&c__1, (char *)&(*nstr), (ftnlen)sizeof(integer));
    do_fio(&c__1, "     No. computational layers =", (ftnlen)31);
    do_fio(&c__1, (char *)&(*nlyr), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ibcnd != 1) {
	s_wsfe(&io___312);
	do_fio(&c__1, (char *)&(*ntau), (ftnlen)sizeof(integer));
	do_fio(&c__1, " User optical depths :", (ftnlen)22);
	i__1 = *ntau;
	for (lu = 1; lu <= i__1; ++lu) {
	    do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    if (! (*onlyfl)) {
	s_wsfe(&io___314);
	do_fio(&c__1, (char *)&(*numu), (ftnlen)sizeof(integer));
	do_fio(&c__1, " User polar angle cosines :", (ftnlen)27);
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    if (! (*onlyfl) && *ibcnd != 1) {
	s_wsfe(&io___316);
	do_fio(&c__1, (char *)&(*nphi), (ftnlen)sizeof(integer));
	do_fio(&c__1, " User azimuthal angles :", (ftnlen)24);
	i__1 = *nphi;
	for (j = 1; j <= i__1; ++j) {
	    do_fio(&c__1, (char *)&phi[j], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    if (! (*plank) || *ibcnd == 1) {
	s_wsfe(&io___318);
	do_fio(&c__1, " No thermal emission", (ftnlen)20);
	e_wsfe();
    }
    s_wsfe(&io___319);
    do_fio(&c__1, " Boundary condition flag: IBCND =", (ftnlen)33);
    do_fio(&c__1, (char *)&(*ibcnd), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ibcnd == 0) {
	s_wsfe(&io___320);
	do_fio(&c__1, "    Incident beam with intensity =", (ftnlen)34);
	do_fio(&c__1, (char *)&(*fbeam), (ftnlen)sizeof(real));
	do_fio(&c__1, " and polar angle cosine = ", (ftnlen)26);
	do_fio(&c__1, (char *)&(*umu0), (ftnlen)sizeof(real));
	do_fio(&c__1, "  and azimuth angle =", (ftnlen)21);
	do_fio(&c__1, (char *)&(*phi0), (ftnlen)sizeof(real));
	do_fio(&c__1, "    plus isotropic incident intensity =", (ftnlen)39);
	do_fio(&c__1, (char *)&(*fisot), (ftnlen)sizeof(real));
	e_wsfe();
	if (*lamber) {
	    s_wsfe(&io___321);
	    do_fio(&c__1, "    Bottom albedo (Lambertian) =", (ftnlen)32);
	    do_fio(&c__1, (char *)&(*albedo), (ftnlen)sizeof(real));
	    e_wsfe();
	}
	if (! (*lamber)) {
	    s_wsfe(&io___322);
	    do_fio(&c__1, "    Legendre coeffs of bottom bidirectional refle\
ctivity :", (ftnlen)58);
	    i__1 = *nstr;
	    for (k = 0; k <= i__1; ++k) {
		do_fio(&c__1, (char *)&hl[k], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
	if (*plank) {
	    s_wsfe(&io___324);
	    do_fio(&c__1, "    Thermal emission in wavenumber interval :", (
		    ftnlen)45);
	    do_fio(&c__1, (char *)&(*wvnmlo), (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&(*wvnmhi), (ftnlen)sizeof(real));
	    do_fio(&c__1, "    Bottom temperature =", (ftnlen)24);
	    do_fio(&c__1, (char *)&(*btemp), (ftnlen)sizeof(real));
	    do_fio(&c__1, "    Top temperature =", (ftnlen)21);
	    do_fio(&c__1, (char *)&(*ttemp), (ftnlen)sizeof(real));
	    do_fio(&c__1, "    Top emissivity =", (ftnlen)20);
	    do_fio(&c__1, (char *)&(*temis), (ftnlen)sizeof(real));
	    e_wsfe();
	}
    } else if (*ibcnd == 1) {
	s_wsfe(&io___325);
	do_fio(&c__1, "    Isotropic illumination from top and bottom", (
		ftnlen)46);
	e_wsfe();
	s_wsfe(&io___326);
	do_fio(&c__1, "    Bottom albedo (Lambertian) =", (ftnlen)32);
	do_fio(&c__1, (char *)&(*albedo), (ftnlen)sizeof(real));
	e_wsfe();
    }
    if (*deltam) {
	s_wsfe(&io___327);
	do_fio(&c__1, " Uses delta-M method", (ftnlen)20);
	e_wsfe();
    }
    if (! (*deltam)) {
	s_wsfe(&io___328);
	do_fio(&c__1, " Does not use delta-M method", (ftnlen)28);
	e_wsfe();
    }
    if (*ibcnd == 1) {
	s_wsfe(&io___329);
	do_fio(&c__1, " Calculate albedo and transmissivity of medium vs. in\
cident beam angle", (ftnlen)70);
	e_wsfe();
    } else if (*onlyfl) {
	s_wsfe(&io___330);
	do_fio(&c__1, " Calculate fluxes and azim-averaged intensities only", 
		(ftnlen)52);
	e_wsfe();
    } else {
	s_wsfe(&io___331);
	do_fio(&c__1, " Calculate fluxes and intensities", (ftnlen)33);
	e_wsfe();
    }
    s_wsfe(&io___332);
    do_fio(&c__1, " Relative convergence criterion for azimuth series =", (
	    ftnlen)52);
    do_fio(&c__1, (char *)&(*accur), (ftnlen)sizeof(real));
    e_wsfe();
    if (*lyrcut) {
	s_wsfe(&io___333);
	do_fio(&c__1, " Sets radiation = 0 below absorption optical depth 10",
		 (ftnlen)53);
	e_wsfe();
    }
/*                                    ** Print layer variables */
/*                                    ** (to read, skip every other line) */
    if (*plank) {
	s_wsfe(&io___334);
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
	s_wsfe(&io___335);
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
    yessct = (float)0.;
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
	yessct += ssalb[lc];
/*                                       ** f90 nonadvancing I/O would */
/*                                       ** simplify this a lot (also the */
/*                                       ** two WRITEs above) */
	if (*plank) {
	    s_wsfe(&io___338);
	    do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dtauc[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&tauc[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ssalb[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&flyr[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dtaucp[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&taucpr[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&oprim[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&pmom[lc * pmom_dim1 + 1], (ftnlen)sizeof(
		    real));
	    do_fio(&c__1, (char *)&temper[lc - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	}
	if (! (*plank)) {
	    s_wsfe(&io___339);
	    do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dtauc[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&tauc[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ssalb[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&flyr[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dtaucp[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&taucpr[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&oprim[lc], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&pmom[lc * pmom_dim1 + 1], (ftnlen)sizeof(
		    real));
	    e_wsfe();
	}
/* L10: */
    }
    if (*plank) {
	s_wsfe(&io___340);
	do_fio(&c__1, (char *)&temper[*nlyr], (ftnlen)sizeof(real));
	e_wsfe();
    }
    if (*prtmom && yessct > (float)0.) {
	s_wsfe(&io___341);
	do_fio(&c__1, " Layer   Phase Function Moments", (ftnlen)31);
	e_wsfe();
	i__1 = *nlyr;
	for (lc = 1; lc <= i__1; ++lc) {
	    if (ssalb[lc] > (float)0.) {
		s_wsfe(&io___342);
		do_fio(&c__1, (char *)&lc, (ftnlen)sizeof(integer));
		i__2 = *nstr;
		for (k = 0; k <= i__2; ++k) {
		    do_fio(&c__1, (char *)&pmom[k + lc * pmom_dim1], (ftnlen)
			    sizeof(real));
		}
		e_wsfe();
	    }
/* L20: */
	}
    }
    return 0;
} /* prtinp_ */

/* Subroutine */ int prtint_(uu, utau, ntau, umu, numu, phi, nphi, maxulv, 
	maxumu)
real *uu, *utau;
integer *ntau;
real *umu;
integer *numu;
real *phi;
integer *nphi, *maxulv, *maxumu;
{
    /* System generated locals */
    integer uu_dim1, uu_dim2, uu_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer jmin, jmax, j, npass, iu, np, lu, lenfmt;

    /* Fortran I/O blocks */
    static cilist io___343 = { 0, 6, 0, "(//,A)", 0 };
    static cilist io___346 = { 0, 6, 0, "(/,A,/,A,/,A)", 0 };
    static cilist io___351 = { 0, 6, 0, "(/,18X,10F11.2)", 0 };
    static cilist io___353 = { 0, 6, 0, "(F10.4,F8.4,1P,10E11.3)", 0 };
    static cilist io___354 = { 0, 6, 0, "(10X,F8.4,1P,10E11.3)", 0 };
    static cilist io___356 = { 0, 6, 0, "(10X,F8.4,1P,10E11.3)", 0 };


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
    --utau;
    --umu;
    --phi;
    uu_dim1 = *maxumu;
    uu_dim2 = *maxulv;
    uu_offset = 1 + uu_dim1 * (1 + uu_dim2 * 1);
    uu -= uu_offset;

    /* Function Body */
    if (*nphi < 1) {
	return 0;
    }
    s_wsfe(&io___343);
    do_fio(&c__1, " *********  I N T E N S I T I E S  *********", (ftnlen)44);
    e_wsfe();
    lenfmt = 10;
    npass = (*nphi - 1) / lenfmt + 1;
    s_wsfe(&io___346);
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
	    s_wsfe(&io___351);
	    i__3 = jmax;
	    for (j = jmin; j <= i__3; ++j) {
		do_fio(&c__1, (char *)&phi[j], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    if (np == 1) {
		s_wsfe(&io___353);
		do_fio(&c__1, (char *)&utau[lu], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&umu[1], (ftnlen)sizeof(real));
		i__3 = jmax;
		for (j = jmin; j <= i__3; ++j) {
		    do_fio(&c__1, (char *)&uu[(lu + j * uu_dim2) * uu_dim1 + 
			    1], (ftnlen)sizeof(real));
		}
		e_wsfe();
	    }
	    if (np > 1) {
		s_wsfe(&io___354);
		do_fio(&c__1, (char *)&umu[1], (ftnlen)sizeof(real));
		i__3 = jmax;
		for (j = jmin; j <= i__3; ++j) {
		    do_fio(&c__1, (char *)&uu[(lu + j * uu_dim2) * uu_dim1 + 
			    1], (ftnlen)sizeof(real));
		}
		e_wsfe();
	    }
	    i__3 = *numu;
	    for (iu = 2; iu <= i__3; ++iu) {
		s_wsfe(&io___356);
		do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
		i__4 = jmax;
		for (j = jmin; j <= i__4; ++j) {
		    do_fio(&c__1, (char *)&uu[iu + (lu + j * uu_dim2) * 
			    uu_dim1], (ftnlen)sizeof(real));
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

/* Subroutine */ int qgausn_(m, gmu, gwt)
integer *m;
real *gmu, *gwt;
{
    /* Initialized data */

    static real pi = (float)0.;
    static integer maxit = 1000;
    static doublereal one = 1.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double asin(), tan(), cos();

    /* Local variables */
    static real cona;
    static integer iter;
    static doublereal prod, p2pri;
    static integer k;
    static doublereal p;
    static real t;
    static doublereal x;
    extern doublereal d1mach_();
    static doublereal en;
    static integer nn;
    static doublereal xi;
    extern /* Subroutine */ int errmsg_();
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
    --gwt;
    --gmu;

    /* Function Body */
    if (pi == (float)0.) {
	pi = asin((float)1.) * (float)2.;
	tol = d1mach_(&c__4) * (float)10.;
    }
    if (*m < 1) {
	errmsg_("QGAUSN--Bad value of M", &c_true, (ftnlen)22);
    }
    if (*m == 1) {
	gmu[1] = (float).5;
	gwt[1] = (float)1.;
	return 0;
    }
    en = (doublereal) (*m);
    np1 = *m + 1;
    nnp1 = (doublereal) (*m * np1);
/* Computing 3rd power */
    i__1 = *m;
    cona = (real) (*m - 1) / (i__1 * (i__1 * i__1) << 3);
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
	gmu[k] = -x;
/* Computing 2nd power */
	d__1 = en * pm2;
	gwt[k] = two / (tmp * (d__1 * d__1));
	gmu[np1 - k] = -gmu[k];
	gwt[np1 - k] = gwt[k];
/* L30: */
    }
/*                                    ** Set middle abscissa and weight */
/*                                    ** for rules of odd order */
    if (*m % 2 != 0) {
	gmu[lim + 1] = (float)0.;
	prod = one;
	i__1 = *m;
	for (k = 3; k <= i__1; k += 2) {
	    prod = prod * k / (k - 1);
/* L40: */
	}
/* Computing 2nd power */
	d__1 = prod;
	gwt[lim + 1] = two / (d__1 * d__1);
    }
/*                                        ** Convert from (-1,1) to (0,1) */
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	gmu[k] = gmu[k] * (float).5 + (float).5;
	gwt[k] *= (float).5;
/* L50: */
    }
    return 0;
} /* qgausn_ */

doublereal ratio_(a, b)
real *a, *b;
{
    /* Initialized data */

    static logical pass1 = TRUE_;

    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double r_lg10(), r_sign();

    /* Local variables */
    static real absa, absb, huge__, powa, powb, tiny;
    extern doublereal r1mach_();
    static real powmin, powmax;

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
	powmax = r_lg10(&huge__);
	powmin = r_lg10(&tiny);
	pass1 = FALSE_;
    }
    if (*a == (float)0.) {
	if (*b == (float)0.) {
	    ret_val = (float)1.;
	} else {
	    ret_val = (float)0.;
	}
    } else if (*b == (float)0.) {
	ret_val = r_sign(&huge__, a);
    } else {
	absa = dabs(*a);
	absb = dabs(*b);
	powa = r_lg10(&absa);
	powb = r_lg10(&absb);
	if (absa < tiny && absb < tiny) {
	    ret_val = (float)1.;
	} else if (powa - powb >= powmax) {
	    ret_val = huge__;
	} else if (powa - powb <= powmin) {
	    ret_val = tiny;
	} else {
	    ret_val = absa / absb;
	}
/*                      ** DONT use old trick of determining sign */
/*                      ** from A*B because A*B may (over/under)flow */
	if (*a > (float)0. && *b < (float)0. || *a < (float)0. && *b > (float)
		0.) {
	    ret_val = -ret_val;
	}
    }
    return ret_val;
} /* ratio_ */

/* Subroutine */ int slftst_(accur, albedo, btemp, deltam, dtauc, fbeam, 
	fisot, ibcnd, lamber, nlyr, plank, nphi, numu, nstr, ntau, onlyfl, 
	phi, phi0, pmom, prnt, ssalb, temis, temper, ttemp, umu, usrang, 
	usrtau, utau, umu0, wvnmhi, wvnmlo, compar, flup, rfldir, rfldn, uu)
real *accur, *albedo, *btemp;
logical *deltam;
real *dtauc, *fbeam, *fisot;
integer *ibcnd;
logical *lamber;
integer *nlyr;
logical *plank;
integer *nphi, *numu, *nstr, *ntau;
logical *onlyfl;
real *phi, *phi0, *pmom;
logical *prnt;
real *ssalb, *temis, *temper, *ttemp, *umu;
logical *usrang, *usrtau;
real *utau, *umu0, *wvnmhi, *wvnmlo;
logical *compar;
real *flup, *rfldir, *rfldn, *uu;
{
    /* Initialized data */

    static real acc = (float)1e-4;

    static real phis, umus, phi0s;
    static integer i__;
    static real umu0s;
    static integer n, nphis, ntaus;
    static real pmoms[5], utaus;
    static logical prnts[7];
    static integer nlyrs, numus, nstrs;
    static real error1, error2, error3, error4;
    static logical ok;
    static real albeds, fbeams;
    static logical lambes;
    static integer ibcnds;
    static logical deltas;
    static real accurs, dtaucs;
    extern logical tstbad_();
    static real ssalbs;
    static logical planks;
    static real btemps, tempes[2];
    extern /* Subroutine */ int errmsg_();
    static real temiss, fisots;
    static logical onlyfs, usrans;
    static real ttemps;
    static logical usrtas;
    static real wvnmhs, wvnmls;

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
    /* Parameter adjustments */
    --prnt;

    /* Function Body */
    if (! (*compar)) {
/*                                     ** Save user input values */
	nlyrs = *nlyr;
	dtaucs = *dtauc;
	ssalbs = *ssalb;
	for (n = 0; n <= 4; ++n) {
	    pmoms[n] = pmom[n];
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
	    prnts[i__ - 1] = prnt[i__];
/* L20: */
	}
/*                                     ** Set input values for self-test */
	*nstr = 4;
	*nlyr = 1;
	*dtauc = (float)1.;
	*ssalb = (float).9;
/*                          ** Haze L moments */
	pmom[0] = (float)1.;
	pmom[1] = (float).8042;
	pmom[2] = (float).646094;
	pmom[3] = (float).481851;
	pmom[4] = (float).359056;
	*usrang = TRUE_;
	*numu = 1;
	*umu = (float).5;
	*usrtau = TRUE_;
	*ntau = 1;
	*utau = (float).5;
	*nphi = 1;
	*phi = (float)90.;
	*ibcnd = 0;
	*fbeam = (float)3.14159265;
	*umu0 = (float).866;
	*phi0 = (float)0.;
	*fisot = (float)1.;
	*lamber = TRUE_;
	*albedo = (float).7;
	*deltam = TRUE_;
	*onlyfl = FALSE_;
	*accur = (float)1e-4;
	*plank = TRUE_;
	*wvnmlo = (float)0.;
	*wvnmhi = (float)5e4;
	*btemp = (float)300.;
	*ttemp = (float)100.;
	*temis = (float).8;
	temper[0] = (float)210.;
	temper[1] = (float)200.;
	for (i__ = 1; i__ <= 7; ++i__) {
	    prnt[i__] = FALSE_;
/* L30: */
	}
    } else {
/*                                    ** Compare test case results with */
/*                                    ** correct answers and abort if bad */
	ok = TRUE_;
	error1 = (*uu - (float)47.86005) / (float)47.86005;
	error2 = (*rfldir - (float)1.527286) / (float)1.527286;
	error3 = (*rfldn - (float)28.37223) / (float)28.37223;
	error4 = (*flup - (float)152.5853) / (float)152.5853;
	if (dabs(error1) > acc) {
	    ok = tstbad_("UU", &error1, (ftnlen)2);
	}
	if (dabs(error2) > acc) {
	    ok = tstbad_("RFLDIR", &error2, (ftnlen)6);
	}
	if (dabs(error3) > acc) {
	    ok = tstbad_("RFLDN", &error3, (ftnlen)5);
	}
	if (dabs(error4) > acc) {
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
	    pmom[n] = pmoms[n];
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
	    prnt[i__] = prnts[i__ - 1];
/* L50: */
	}
    }
    return 0;
} /* slftst_ */

/* Subroutine */ int zeroal_(nd1, expbea, flyr, oprim, taucpr, xr0, xr1, nd2, 
	cmu, cwt, psi, wk, z0, z1, zj, nd3, hlpr, ylm0, nd4, array, cc, evecc,
	 nd5, gl, nd6, ylmc, nd7, ylmu, nd8, kk, ll, zz, zplk0, zplk1, nd9, 
	gc, nd10, layru, utaupr, nd11, gu, nd12, z0u, z1u, zbeam, nd13, eval, 
	nd14, amb, apb, nd15, ipvt, z__, nd16, rfldir, rfldn, flup, uavg, 
	dfdt, nd17, albmed, trnmed, nd18, u0u, nd19, uu)
integer *nd1;
real *expbea, *flyr, *oprim, *taucpr, *xr0, *xr1;
integer *nd2;
real *cmu, *cwt, *psi, *wk, *z0, *z1, *zj;
integer *nd3;
real *hlpr, *ylm0;
integer *nd4;
real *array, *cc, *evecc;
integer *nd5;
real *gl;
integer *nd6;
real *ylmc;
integer *nd7;
real *ylmu;
integer *nd8;
real *kk, *ll, *zz, *zplk0, *zplk1;
integer *nd9;
real *gc;
integer *nd10, *layru;
real *utaupr;
integer *nd11;
real *gu;
integer *nd12;
real *z0u, *z1u, *zbeam;
integer *nd13;
real *eval;
integer *nd14;
real *amb, *apb;
integer *nd15, *ipvt;
real *z__;
integer *nd16;
real *rfldir, *rfldn, *flup, *uavg, *dfdt;
integer *nd17;
real *albmed, *trnmed;
integer *nd18;
real *u0u;
integer *nd19;
real *uu;
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
    /* Parameter adjustments */
    --uu;
    --u0u;
    --trnmed;
    --albmed;
    --dfdt;
    --uavg;
    --flup;
    --rfldn;
    --rfldir;
    --z__;
    --ipvt;
    --apb;
    --amb;
    --eval;
    --zbeam;
    --z1u;
    --z0u;
    --gu;
    --utaupr;
    --layru;
    --gc;
    --zplk1;
    --zplk0;
    --zz;
    --ll;
    --kk;
    --ylmu;
    --ylmc;
    --gl;
    --evecc;
    --cc;
    --array;
    --ylm0;
    --hlpr;
    --zj;
    --z1;
    --z0;
    --wk;
    --psi;
    --cwt;
    --cmu;
    --xr1;
    --xr0;
    --taucpr;
    --oprim;
    --flyr;
    --expbea;

    /* Function Body */
    i__1 = *nd1;
    for (n = 1; n <= i__1; ++n) {
	expbea[n] = (float)0.;
	flyr[n] = (float)0.;
	oprim[n] = (float)0.;
	taucpr[n] = (float)0.;
	xr0[n] = (float)0.;
	xr1[n] = (float)0.;
/* L10: */
    }
    i__1 = *nd2;
    for (n = 1; n <= i__1; ++n) {
	cmu[n] = (float)0.;
	cwt[n] = (float)0.;
	psi[n] = (float)0.;
	wk[n] = (float)0.;
	z0[n] = (float)0.;
	z1[n] = (float)0.;
	zj[n] = (float)0.;
/* L20: */
    }
    i__1 = *nd3;
    for (n = 1; n <= i__1; ++n) {
	hlpr[n] = (float)0.;
	ylm0[n] = (float)0.;
/* L30: */
    }
    i__1 = *nd4;
    for (n = 1; n <= i__1; ++n) {
	array[n] = (float)0.;
	cc[n] = (float)0.;
	evecc[n] = (float)0.;
/* L40: */
    }
    i__1 = *nd5;
    for (n = 1; n <= i__1; ++n) {
	gl[n] = (float)0.;
/* L50: */
    }
    i__1 = *nd6;
    for (n = 1; n <= i__1; ++n) {
	ylmc[n] = (float)0.;
/* L60: */
    }
    i__1 = *nd7;
    for (n = 1; n <= i__1; ++n) {
	ylmu[n] = (float)0.;
/* L70: */
    }
    i__1 = *nd8;
    for (n = 1; n <= i__1; ++n) {
	kk[n] = (float)0.;
	ll[n] = (float)0.;
	zz[n] = (float)0.;
	zplk0[n] = (float)0.;
	zplk1[n] = (float)0.;
/* L80: */
    }
    i__1 = *nd9;
    for (n = 1; n <= i__1; ++n) {
	gc[n] = (float)0.;
/* L90: */
    }
    i__1 = *nd10;
    for (n = 1; n <= i__1; ++n) {
	layru[n] = 0;
	utaupr[n] = (float)0.;
/* L100: */
    }
    i__1 = *nd11;
    for (n = 1; n <= i__1; ++n) {
	gu[n] = (float)0.;
/* L110: */
    }
    i__1 = *nd12;
    for (n = 1; n <= i__1; ++n) {
	z0u[n] = (float)0.;
	z1u[n] = (float)0.;
	zbeam[n] = (float)0.;
/* L120: */
    }
    i__1 = *nd13;
    for (n = 1; n <= i__1; ++n) {
	eval[n] = (float)0.;
/* L130: */
    }
    i__1 = *nd14;
    for (n = 1; n <= i__1; ++n) {
	amb[n] = (float)0.;
	apb[n] = (float)0.;
/* L140: */
    }
    i__1 = *nd15;
    for (n = 1; n <= i__1; ++n) {
	ipvt[n] = 0;
	z__[n] = (float)0.;
/* L150: */
    }
    i__1 = *nd16;
    for (n = 1; n <= i__1; ++n) {
	rfldir[n] = (float)0.;
	rfldn[n] = (float)0.;
	flup[n] = (float)0.;
	uavg[n] = (float)0.;
	dfdt[n] = (float)0.;
/* L160: */
    }
    i__1 = *nd17;
    for (n = 1; n <= i__1; ++n) {
	albmed[n] = (float)0.;
	trnmed[n] = (float)0.;
/* L170: */
    }
    i__1 = *nd18;
    for (n = 1; n <= i__1; ++n) {
	u0u[n] = (float)0.;
/* L180: */
    }
    i__1 = *nd19;
    for (n = 1; n <= i__1; ++n) {
	uu[n] = (float)0.;
/* L190: */
    }
    return 0;
} /* zeroal_ */

/* Subroutine */ int zeroit_(a, length)
real *a;
integer *length;
{
    /* System generated locals */
    integer i__1;

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
    --a;

    /* Function Body */
    i__1 = *length;
    for (l = 1; l <= i__1; ++l) {
	a[l] = (float)0.;
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
/* Subroutine */ int albtrn_(albedo, amb, apb, array, b, bdr, cband, cc, cmu, 
	cwt, dtaucp, eval, evecc, gl, gc, gu, ipvt, kk, ll, nlyr, nn, nstr, 
	numu, prnt, taucpr, umu, u0u, wk, ylmc, ylmu, z__, aad, evald, eveccd,
	 wkd, mi, mi9m2, maxulv, maxumu, mxcmu, mxumu, nnlyri, sqt, albmed, 
	trnmed)
real *albedo, *amb, *apb, *array, *b, *bdr, *cband, *cc, *cmu, *cwt, *dtaucp, 
	*eval, *evecc, *gl, *gc, *gu;
integer *ipvt;
real *kk, *ll;
integer *nlyr, *nn, *nstr, *numu;
logical *prnt;
real *taucpr, *umu, *u0u, *wk, *ylmc, *ylmu, *z__;
doublereal *aad, *evald, *eveccd, *wkd;
integer *mi, *mi9m2, *maxulv, *maxumu, *mxcmu, *mxumu, *nnlyri;
real *sqt, *albmed, *trnmed;
{
    /* System generated locals */
    integer amb_dim1, amb_offset, apb_dim1, apb_offset, array_dim1, 
	    array_offset, bdr_dim1, bdr_offset, cband_dim1, cband_offset, 
	    cc_dim1, cc_offset, evecc_dim1, evecc_offset, gc_dim1, gc_dim2, 
	    gc_offset, gl_dim1, gl_offset, gu_dim1, gu_dim2, gu_offset, 
	    kk_dim1, kk_offset, ll_dim1, ll_offset, u0u_dim1, u0u_offset, 
	    ylmc_dim1, ylmc_offset, ylmu_dim1, ylmu_offset, aad_dim1, 
	    aad_offset, eveccd_dim1, eveccd_offset, i__1, i__2;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static integer ncol, ncut;
    static real delm0;
    static integer l;
    extern /* Subroutine */ int sgbco_();
    static real rcond;
    static integer mazim;
    static real fisot;
    extern /* Subroutine */ int solve1_();
    static integer lc, iq, iu;
    static logical lamber;
    static real sphalb;
    extern /* Subroutine */ int soleig_(), altrin_(), errmsg_(), lepoly_(), 
	    praltr_(), spaltr_(), terpev_();
    static real sphtrn;
    extern /* Subroutine */ int zeroit_();
    static logical lyrcut;
    extern /* Subroutine */ int setmtx_();
    static integer ncd;
    static real sgn;

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
    --dtaucp;
    --ipvt;
    --prnt;
    eveccd_dim1 = *mi;
    eveccd_offset = 1 + eveccd_dim1 * 1;
    eveccd -= eveccd_offset;
    --evald;
    aad_dim1 = *mi;
    aad_offset = 1 + aad_dim1 * 1;
    aad -= aad_offset;
    --eval;
    bdr_dim1 = *mi;
    bdr_offset = 1 + bdr_dim1 * 0;
    bdr -= bdr_offset;
    apb_dim1 = *mi;
    apb_offset = 1 + apb_dim1 * 1;
    apb -= apb_offset;
    amb_dim1 = *mi;
    amb_offset = 1 + amb_dim1 * 1;
    amb -= amb_offset;
    --trnmed;
    --albmed;
    u0u_dim1 = *maxumu;
    u0u_offset = 1 + u0u_dim1 * 1;
    u0u -= u0u_offset;
    --umu;
    --wkd;
    ylmu_dim1 = *mxcmu - 0 + 1;
    ylmu_offset = 0 + ylmu_dim1 * 1;
    ylmu -= ylmu_offset;
    ylmc_dim1 = *mxcmu - 0 + 1;
    ylmc_offset = 0 + ylmc_dim1 * 1;
    ylmc -= ylmc_offset;
    --wk;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    gc -= gc_offset;
    gl_dim1 = *mxcmu - 0 + 1;
    gl_offset = 0 + gl_dim1 * 1;
    gl -= gl_offset;
    evecc_dim1 = *mxcmu;
    evecc_offset = 1 + evecc_dim1 * 1;
    evecc -= evecc_offset;
    --cwt;
    --cmu;
    cc_dim1 = *mxcmu;
    cc_offset = 1 + cc_dim1 * 1;
    cc -= cc_offset;
    array_dim1 = *mxcmu;
    array_offset = 1 + array_dim1 * 1;
    array -= array_offset;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2 * 1);
    gu -= gu_offset;
    --z__;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1 * 1;
    cband -= cband_offset;
    --b;
    --sqt;

    /* Function Body */
    mazim = 0;
    delm0 = (float)1.;
/*                    ** Set DISORT variables that are ignored in this */
/*                    ** special case but are needed below in argument */
/*                    ** lists of subroutines shared with general case */
    ncut = *nlyr;
    lyrcut = FALSE_;
    fisot = (float)1.;
    lamber = TRUE_;
/*                          ** Get Legendre polynomials for computational */
/*                          ** and user polar angle cosines */
    i__1 = *nstr - 1;
    lepoly_(numu, &mazim, mxcmu, &i__1, &umu[1], &sqt[1], &ylmu[ylmu_offset]);
    i__1 = *nstr - 1;
    lepoly_(nn, &mazim, mxcmu, &i__1, &cmu[1], &sqt[1], &ylmc[ylmc_offset]);
/*                       ** Evaluate Legendre polynomials with negative */
/*                       ** arguments from those with positive arguments; */
/*                       ** Dave/Armstrong Eq. (15) */
    sgn = (float)-1.;
    i__1 = *nstr - 1;
    for (l = mazim; l <= i__1; ++l) {
	sgn = -sgn;
	i__2 = *nstr;
	for (iq = *nn + 1; iq <= i__2; ++iq) {
	    ylmc[l + iq * ylmc_dim1] = sgn * ylmc[l + (iq - *nn) * ylmc_dim1];
/* L10: */
	}
/* L20: */
    }
/*                                  ** Zero out bottom reflectivity */
/*                                  ** (ALBEDO is used only in analytic */
/*                                  ** formulae involving ALBEDO = 0 */
/*                                  ** solutions; Eqs 16-17 of Ref S2) */
    i__1 = *mi * (*mi + 1);
    zeroit_(&bdr[bdr_offset], &i__1);
/* ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  ============= */
    i__1 = *nlyr;
    for (lc = 1; lc <= i__1; ++lc) {
/*                        ** Solve eigenfunction problem in Eq. STWJ(8b) */
	soleig_(&amb[amb_offset], &apb[apb_offset], &array[array_offset], &
		cmu[1], &cwt[1], &gl[lc * gl_dim1], mi, &mazim, mxcmu, nn, 
		nstr, &ylmc[ylmc_offset], &cc[cc_offset], &evecc[evecc_offset]
		, &eval[1], &kk[lc * kk_dim1 + 1], &gc[(lc * gc_dim2 + 1) * 
		gc_dim1 + 1], &aad[aad_offset], &eveccd[eveccd_offset], &
		evald[1], &wkd[1]);
/*                          ** Interpolate eigenvectors to user angles */
	terpev_(&cwt[1], &evecc[evecc_offset], &gl[lc * gl_dim1], &gu[(lc * 
		gu_dim2 + 1) * gu_dim1 + 1], &mazim, mxcmu, mxumu, nn, nstr, 
		numu, &wk[1], &ylmc[ylmc_offset], &ylmu[ylmu_offset]);
/* L30: */
    }
/* ===================  END LOOP ON COMPUTATIONAL LAYERS  =============== */
/*                      ** Set coefficient matrix (CBAND) of equations */
/*                      ** combining boundary and layer interface */
/*                      ** conditions (in band-storage mode required by */
/*                      ** LINPACK routines) */
    setmtx_(&bdr[bdr_offset], &cband[cband_offset], &cmu[1], &cwt[1], &delm0, 
	    &dtaucp[1], &gc[gc_offset], &kk[kk_offset], &lamber, &lyrcut, mi, 
	    mi9m2, mxcmu, &ncol, &ncut, nnlyri, nn, nstr, taucpr, &wk[1]);
/*                      ** LU-decompose the coeff. matrix (LINPACK) */
    ncd = *nn * 3 - 1;
    sgbco_(&cband[cband_offset], mi9m2, &ncol, &ncd, &ncd, &ipvt[1], &rcond, &
	    z__[1]);
    if (rcond + (float)1. == (float)1.) {
	errmsg_("ALBTRN--SGBCO says matrix near singular", &c_false, (ftnlen)
		39);
    }
/*                             ** First, illuminate from top; if only */
/*                             ** one layer, this will give us everything */
/*                             ** Solve for constants of integration in */
/*                             ** homogeneous solution */
    solve1_(&b[1], &cband[cband_offset], &fisot, &c__1, &ipvt[1], &ll[
	    ll_offset], mi9m2, mxcmu, &ncol, nlyr, nn, nnlyri, nstr);
/*                             ** Compute azimuthally-averaged intensity */
/*                             ** at user angles; gives albedo if multi- */
/*                             ** layer (Eq. 9 of Ref S2); gives both */
/*                             ** albedo and transmissivity if single */
/*                             ** layer (Eqs. 3-4 of Ref S2) */
    altrin_(&gu[gu_offset], &kk[kk_offset], &ll[ll_offset], mxcmu, mxumu, 
	    maxumu, nlyr, nn, nstr, numu, taucpr, &umu[1], &u0u[u0u_offset], &
	    wk[1]);
/*                               ** Get beam-incidence albedos from */
/*                               ** reciprocity principle */
    i__1 = *numu / 2;
    for (iu = 1; iu <= i__1; ++iu) {
	albmed[iu] = u0u[iu + *numu / 2 + u0u_dim1];
/* L40: */
    }
    if (*nlyr == 1) {
	i__1 = *numu / 2;
	for (iu = 1; iu <= i__1; ++iu) {
/*                               ** Get beam-incidence transmissivities */
/*                               ** from reciprocity principle (1 layer); */
/*                               ** flip them end over end to correspond */
/*                               ** to positive UMU instead of negative */
	    trnmed[iu] = u0u[*numu / 2 + 1 - iu + (u0u_dim1 << 1)] + exp(
		    -taucpr[*nlyr] / umu[iu + *numu / 2]);
/* L50: */
	}
    } else {
/*                             ** Second, illuminate from bottom */
/*                             ** (if multiple layers) */
	solve1_(&b[1], &cband[cband_offset], &fisot, &c__2, &ipvt[1], &ll[
		ll_offset], mi9m2, mxcmu, &ncol, nlyr, nn, nnlyri, nstr);
	altrin_(&gu[gu_offset], &kk[kk_offset], &ll[ll_offset], mxcmu, mxumu, 
		maxumu, nlyr, nn, nstr, numu, taucpr, &umu[1], &u0u[
		u0u_offset], &wk[1]);
/*                               ** Get beam-incidence transmissivities */
/*                               ** from reciprocity principle */
	i__1 = *numu / 2;
	for (iu = 1; iu <= i__1; ++iu) {
	    trnmed[iu] = u0u[iu + *numu / 2 + u0u_dim1] + exp(-taucpr[*nlyr] /
		     umu[iu + *numu / 2]);
/* L60: */
	}
    }
    if (*albedo > (float)0.) {
/*                             ** Get spherical albedo and transmissivity */
	if (*nlyr == 1) {
	    spaltr_(&cmu[1], &cwt[1], &gc[gc_offset], &kk[kk_offset], &ll[
		    ll_offset], mxcmu, nlyr, nn, nstr, taucpr, &sphalb, &
		    sphtrn);
	} else {
	    spaltr_(&cmu[1], &cwt[1], &gc[gc_offset], &kk[kk_offset], &ll[
		    ll_offset], mxcmu, nlyr, nn, nstr, taucpr, &sphtrn, &
		    sphalb);
	}
/*                                ** Ref. S2, Eqs. 16-17 (these eqs. have */
/*                                ** a simple physical interpretation */
/*                                ** like that of adding-doubling eqs.) */
	i__1 = *numu;
	for (iu = 1; iu <= i__1; ++iu) {
	    albmed[iu] += *albedo / ((float)1. - *albedo * sphalb) * sphtrn * 
		    trnmed[iu];
	    trnmed[iu] += *albedo / ((float)1. - *albedo * sphalb) * sphalb * 
		    trnmed[iu];
/* L70: */
	}
    }
/*                          ** Return UMU to all positive values, to */
/*                          ** agree with ordering in ALBMED, TRNMED */
    *numu /= 2;
    i__1 = *numu;
    for (iu = 1; iu <= i__1; ++iu) {
	umu[iu] = umu[iu + *numu];
/* L80: */
    }
    if (prnt[6]) {
	praltr_(&umu[1], numu, &albmed[1], &trnmed[1]);
    }
    return 0;
} /* albtrn_ */

/* Subroutine */ int altrin_(gu, kk, ll, mxcmu, mxumu, maxumu, nlyr, nn, nstr,
	 numu, taucpr, umu, u0u, wk)
real *gu, *kk, *ll;
integer *mxcmu, *mxumu, *maxumu, *nlyr, *nn, *nstr, *numu;
real *taucpr, *umu, *u0u, *wk;
{
    /* System generated locals */
    integer gu_dim1, gu_dim2, gu_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, u0u_dim1, u0u_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static real dtau, expn, denom;
    static integer iumin, iumax, lc, iq, iu, lu;
    static real mu, palint, utaupr[2], sgn, exp1, exp2;

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
    --wk;
    ll_dim1 = *mxcmu;
    ll_offset = 1 + ll_dim1 * 1;
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    kk -= kk_offset;
    gu_dim1 = *mxumu;
    gu_dim2 = *mxcmu;
    gu_offset = 1 + gu_dim1 * (1 + gu_dim2 * 1);
    gu -= gu_offset;
    u0u_dim1 = *maxumu;
    u0u_offset = 1 + u0u_dim1 * 1;
    u0u -= u0u_offset;
    --umu;

    /* Function Body */
    utaupr[0] = (float)0.;
    utaupr[1] = taucpr[*nlyr];
    for (lu = 1; lu <= 2; ++lu) {
	if (lu == 1) {
	    iumin = *numu / 2 + 1;
	    iumax = *numu;
	    sgn = (float)1.;
	} else {
	    iumin = 1;
	    iumax = *numu / 2;
	    sgn = (float)-1.;
	}
/*                                   ** Loop over polar angles at which */
/*                                   ** albedos/transmissivities desired */
/*                                   ** ( upward angles at top boundary, */
/*                                   ** downward angles at bottom ) */
	i__1 = iumax;
	for (iu = iumin; iu <= i__1; ++iu) {
	    mu = umu[iu];
/*                                     ** Integrate from top to bottom */
/*                                     ** computational layer */
	    palint = (float)0.;
	    i__2 = *nlyr;
	    for (lc = 1; lc <= i__2; ++lc) {
		dtau = taucpr[lc] - taucpr[lc - 1];
		exp1 = exp((utaupr[lu - 1] - taucpr[lc - 1]) / mu);
		exp2 = exp((utaupr[lu - 1] - taucpr[lc]) / mu);
/*                                      ** KK is negative */
		i__3 = *nn;
		for (iq = 1; iq <= i__3; ++iq) {
		    wk[iq] = exp(kk[iq + lc * kk_dim1] * dtau);
		    denom = mu * kk[iq + lc * kk_dim1] + (float)1.;
		    if (dabs(denom) < (float)1e-4) {
/*                                                   ** L'Hospital limit */
			expn = dtau / mu * exp2;
		    } else {
			expn = (exp1 * wk[iq] - exp2) * sgn / denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1] * ll[iq 
			    + lc * ll_dim1] * expn;
/* L10: */
		}
/*                                        ** KK is positive */
		i__3 = *nstr;
		for (iq = *nn + 1; iq <= i__3; ++iq) {
		    denom = mu * kk[iq + lc * kk_dim1] + (float)1.;
		    if (dabs(denom) < (float)1e-4) {
			expn = -dtau / mu * exp1;
		    } else {
			expn = (exp1 - exp2 * wk[*nstr + 1 - iq]) * sgn / 
				denom;
		    }
		    palint += gu[iu + (iq + lc * gu_dim2) * gu_dim1] * ll[iq 
			    + lc * ll_dim1] * expn;
/* L20: */
		}
/* L30: */
	    }
	    u0u[iu + lu * u0u_dim1] = palint;
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* altrin_ */

/* Subroutine */ int praltr_(umu, numu, albmed, trnmed)
real *umu;
integer *numu;
real *albmed, *trnmed;
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();
    double acos();

    /* Local variables */
    static integer iu;

    /* Fortran I/O blocks */
    static cilist io___461 = { 0, 6, 0, "(///,A,//,A)", 0 };
    static cilist io___463 = { 0, 6, 0, "(0P,F13.4,F20.6,F12.5,1P,E17.4)", 0 }
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
    --trnmed;
    --albmed;
    --umu;

    /* Function Body */
    s_wsfe(&io___461);
    do_fio(&c__1, " *******  Flux Albedo and/or Transmissivity of entire med\
ium  ********", (ftnlen)70);
    do_fio(&c__1, " Beam Zen Ang   cos(Beam Zen Ang)      Albedo   Transmiss\
ivity", (ftnlen)62);
    e_wsfe();
    i__1 = *numu;
    for (iu = 1; iu <= i__1; ++iu) {
	s_wsfe(&io___463);
	r__1 = acos(umu[iu]) * (float)57.295779578552292;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&umu[iu], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&albmed[iu], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&trnmed[iu], (ftnlen)sizeof(real));
	e_wsfe();
/* L10: */
    }
    return 0;
} /* praltr_ */

/* Subroutine */ int solve1_(b, cband, fisot, ihom, ipvt, ll, mi9m2, mxcmu, 
	ncol, ncut, nn, nnlyri, nstr)
real *b, *cband, *fisot;
integer *ihom, *ipvt;
real *ll;
integer *mi9m2, *mxcmu, *ncol, *ncut, *nn, *nnlyri, *nstr;
{
    /* System generated locals */
    integer cband_dim1, cband_offset, ll_dim1, ll_offset, i__1, i__2;

    /* Local variables */
    static integer ipnt, i__;
    extern /* Subroutine */ int sgbsl_();
    static integer lc, iq;
    extern /* Subroutine */ int zeroit_();
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
    ll -= ll_offset;
    --ipvt;
    cband_dim1 = *mi9m2;
    cband_offset = 1 + cband_dim1 * 1;
    cband -= cband_offset;
    --b;

    /* Function Body */
    zeroit_(&b[1], nnlyri);
    if (*ihom == 1) {
/*                             ** Because there are no beam or emission */
/*                             ** sources, remainder of B array is zero */
	i__1 = *nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__] = *fisot;
	    b[*ncol - *nn + i__] = (float)0.;
/* L10: */
	}
    } else if (*ihom == 2) {
	i__1 = *nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__] = (float)0.;
	    b[*ncol - *nn + i__] = *fisot;
/* L20: */
	}
    }
    ncd = *nn * 3 - 1;
    sgbsl_(&cband[cband_offset], mi9m2, ncol, &ncd, &ncd, &ipvt[1], &b[1], &
	    c__0);
    i__1 = *ncut;
    for (lc = 1; lc <= i__1; ++lc) {
	ipnt = lc * *nstr - *nn;
	i__2 = *nn;
	for (iq = 1; iq <= i__2; ++iq) {
	    ll[*nn + 1 - iq + lc * ll_dim1] = b[ipnt + 1 - iq];
	    ll[iq + *nn + lc * ll_dim1] = b[iq + ipnt];
/* L30: */
	}
/* L40: */
    }
    return 0;
} /* solve1_ */

/* Subroutine */ int spaltr_(cmu, cwt, gc, kk, ll, mxcmu, nlyr, nn, nstr, 
	taucpr, sflup, sfldn)
real *cmu, *cwt, *gc, *kk, *ll;
integer *mxcmu, *nlyr, *nn, *nstr;
real *taucpr, *sflup, *sfldn;
{
    /* System generated locals */
    integer gc_dim1, gc_dim2, gc_offset, kk_dim1, kk_offset, ll_dim1, 
	    ll_offset, i__1, i__2;

    /* Builtin functions */
    double exp();

    /* Local variables */
    static real zint;
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
    ll -= ll_offset;
    kk_dim1 = *mxcmu;
    kk_offset = 1 + kk_dim1 * 1;
    kk -= kk_offset;
    gc_dim1 = *mxcmu;
    gc_dim2 = *mxcmu;
    gc_offset = 1 + gc_dim1 * (1 + gc_dim2 * 1);
    gc -= gc_offset;
    --cwt;
    --cmu;

    /* Function Body */
    *sflup = (float)0.;
    i__1 = *nstr;
    for (iq = *nn + 1; iq <= i__1; ++iq) {
	zint = (float)0.;
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + gc_dim2) * gc_dim1] * ll[jq + ll_dim1] * 
		    exp(kk[jq + kk_dim1] * taucpr[1]);
/* L10: */
	}
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + gc_dim2) * gc_dim1] * ll[jq + ll_dim1];
/* L20: */
	}
	*sflup += cwt[iq - *nn] * cmu[iq - *nn] * zint;
/* L30: */
    }
    *sfldn = (float)0.;
    i__1 = *nn;
    for (iq = 1; iq <= i__1; ++iq) {
	zint = (float)0.;
	i__2 = *nn;
	for (jq = 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1] * ll[jq + *nlyr 
		    * ll_dim1];
/* L40: */
	}
	i__2 = *nstr;
	for (jq = *nn + 1; jq <= i__2; ++jq) {
	    zint += gc[iq + (jq + *nlyr * gc_dim2) * gc_dim1] * ll[jq + *nlyr 
		    * ll_dim1] * exp(-kk[jq + *nlyr * kk_dim1] * (taucpr[*
		    nlyr] - taucpr[*nlyr - 1]));
/* L50: */
	}
	*sfldn += cwt[*nn + 1 - iq] * cmu[*nn + 1 - iq] * zint;
/* L60: */
    }
    *sflup *= (float)2.;
    *sfldn *= (float)2.;
    return 0;
} /* spaltr_ */

