/************************************************************************
 * $Id: cdisort.c 2966 2013-07-24 08:58:48Z svn-kylling $
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
 
/*
 * DISORT: Discrete Ordinates Radiative Transfer Program
 *
 * Version 2.1.2 (rewritten in C from Fortran DISORT.f, v 2.1)
 *
 * The C-version has the following extensions NOT present in the original fortran version:
 *  - Correction of intensity fields by Buras-Emde algorithm, included by Robert Buras.
 *  - Pseudospherical geometry for direct beam source, included by Arve Kylling.
 *  - Solution for a general source term, included by Arve Kylling.
 *
 * This file, cdisort.c, contains the source code for disort(), twostr(), and supporting subroutines.
 * The corresponding header file is cdisort.h.
 *
 * See DISORT.txt for a full description and history of the Fortran version.
 *
 *-----------------------REFERENCES (cited in code using the acronyms shown)------------------------------------
 *
 *     BDE: Buras R, Dowling T, Emde C, 201X, ...
 *    DGIS: Devaux C, Grandjean P, Ishiguro Y, Siewert CE, 1979, 
 *              On Multi-Region Problems in Radiative Transfer, Astrophys. Space Sci. 62, 225-233
 *      GS: Garcia RDM, Siewert CE, 1985, Benchmark Results in Radiative Transfer,
 *              Transport Theory and Statistical Physics 14, 437-483
 *      KS: Kylling A, Stamnes K, 1992, Efficient yet accurate solution of the linear transport
 *              equation in the presence of internal sources: The exponential-linear-in-depth
 *              approximation, J. Comp. Phys., 102, 265-276
 *       L: Lenoble J, ed, 1985:  Radiative Transfer in Absorbing
 *              and Scattering Atmospheres: Standard Computational Procedures, Deepak Publishing, Hampton, Virginia
 *      NT: Nakajima T, Tanaka M, 1988,  Algorithms for Radiative Intensity Calculations in 
 *              Moderately Thick Atmospheres Using a Truncation Approximation, J.Q.S.R.T. 40, 51-69
 *      OS: Ozisik M, Shouman S, 1980,  Source Function Expansion Method for Radiative Transfer in a Two-Layer
 *              Slab, J.Q.S.R.T. 24, 441-449
 *      SS: Stamnes K, Swanson R, 1981,  A New Look at the Discrete Ordinate Method for Radiative Transfer
 *              Calculations in Anisotropically Scattering Atmospheres, J. Atmos. Sci. 38, 387-399
 *      SD: Stamnes K, Dale H, 1981, A New Look at the Discrete Ordinate Method for Radiative Transfer
 *              Calculations in Anisotropically Scattering Atmospheres. II: Intensity Computations,
 *              J. Atmos. Sci. 38, 2696-2706
 *      S1: Stamnes K, 1982, On the Computation of Angular Distributions of Radiation in Planetary
 *              Atmospheres, J.Q.S.R.T. 28, 47-51
 *      S2: Stamnes K, 1982, Reflection and Transmission by a Vertically Inhomogeneous Planetary Atmosphere,
 *              Planet. Space Sci. 30, 727-732
 *      SC: Stamnes K, Conklin P, 1984, A New Multi-Layer Discrete Ordinate Approach to Radiative Transfer
 *              in Vertically Inhomogeneous Atmospheres, J.Q.S.R.T. 31, 273-282
 *      SW: Sweigart A, 1970, Radiative Transfer in Atmospheres Scattering According to the Rayleigh Phase Function
 *              with Absorption, The Astrophysical Journal Supplement Series 22, 1-80
 *    STWJ: Stamnes K, Tsay SC, Wiscombe W, Jayaweera K, 1988, A Numerically Stable Algorithm for
 *              Discrete-Ordinate-Method Radiative Transfer in Multiple Scattering and Emitting Layered Media,
 *              Appl. Opt. 27, 2502-2509
 *    STWL: Stamnes K, Tsay SC, Wiscombe W, Laszlo I: A General-Purpose Numerically Stable Computer
 *              Code for Discrete-Ordinate-Method Radiative Transfer in Scattering and Emitting Layered Media,
 *              DISORT Report v1.1 (2000)
 * VH1,VH2: Van de Hulst, H.C., 1980: Multiple Light Scattering, Tables, Formulas and Applications, Volumes 1 and 2,
 *              Academic Press, New York.
 *       W: Wiscombe, W., 1977:  The Delta-M Method: Rapid Yet Accurate Radiative Flux Calculations, J. Atmos. Sci.
 *              34, 1408-1422
 *-----------------------------------------------------------------------------------------------------------------
 *
 * This file contains the following functions:
 *
 *   c_disort()...................Plane-parallel discrete ordinates radiative transfer program
 *   c_bidir_reflectivity().......Supplies surface bi-directional reflectivity (Fortran name bdref).
 *   c_getmom()...................Calculate phase function Legendre expansion coefficients in various special
 *                                cases.
 *   c_asymmetric_matrix()........Solve eigenfunction problem for real asymmetric matrix known a priori
 *                                to have real eigenvalues (Fortran name asymtx).
 *   c_intensity_components().....Calculate the Fourier intensity components at the quadrature angles for azimuthal
 *                                expansion terms (mazim) in eq. SD(2), STWL(6) (Fortran name cmpint).
 *   c_fluxes()...................Calculate the radiative fluxes, mean intensity, and flux derivative with respect
 *                                to optical depth from the m=0 intensity components (the azimuthally-averaged
 *                                intensity).
 *   c_intensity_correction().....Correct intensity field by using Nakajima-Tanaka (1988) algorithm (Fortran name
 *                                intcor).
 *   c_secondary_scat()...........Calculate secondary scattered intensity of eq. STWL (A7) (Fortran name secsca).
 *   c_new_intensity_correction().Correct intensity field by using Buras-Emde (201X) algorithm (Fortran name
 *                                intcor3).
 *   prep_double_scat_integr()....Prepare integration for double scattering intensity correction (see BDE).
 *   c_new_secondary_scat().......Calculate secondary scattered intensity of eq. BDE (XX) (Fortran name secsca3).
 *   calc_phase_squared().........Calculate double scattering intensity correction (see BDE).
 *   c_disort_set()...............Perform misc. setting-up operations (Fortran name setdis).
 *   c_set_matrix()...............Calculate coefficient matrix for the set of equations obtained from the boundary
 *                                conditions and the continuity-of-intensity-at-layer-interface equations (Fortran
 *                                name setmtx).
 *   c_single_scat()..............Calculates single-scattered intensity from eqs. STWL (65b,d,e) (Fortran name
 *                                sinsca).
 *   c_solve_eigen()..............Solve eigenvalue/vector problem necessary to construct homogeneous part of
 *                                discrete ordinate solution; STWJ(8b), STWL(23f) (Fortran name soleig).
 *   c_solve0()...................Construct right-hand side vector -b- for general boundary conditions STWJ(17) and
 *                                solve system of eqns. obtained from the b.c.s and the
 *                                continuity-of-intensity-at-layer-interface eqns.
 *   c_solve1()...................Construct right-hand side vector -b- for isotropic incidence (only) on either top
 *                                or bottom boundary and solve system of eqns. obtained from the b.c.s and the
 *                                continuity-of-intensity-at-layer-interface eqns.
 *   c_twostr_solve_bc()..........Construct right-hand side vector -b- for general b.c. and solve system of eqns.
 *                                obtained from the b.c.s and the continuity-of-intensity-at-layer-interface eqns.
 *   c_surface_bidir()............Compute user's surface bidirectional properties, STWL(41) (Fortran name surfac).
 *   c_interp_eigenvec()..........Interpolate eigenvectors to user angles; eq SD(8) (Fortran name terpev).
 *   c_interp_source()............Interpolate source functions to user angles, eq. STWL(30) (Fortran name terpso).
 *   c_upbeam()...................Find the incident-beam particular solution of SS(18), STWL(24a).
 *   c_upisot()...................Find the particular solution of thermal radiation of STWL(25).
 *   c_user_intensities().........Compute intensity components at user output angles for azimuthal expansion terms
 *                                in eq. SD(2), STWL(6) (Fortran name usrint).
 *   c_xi_fun()...................Calculates Xi function of eq. STWL (72) (Fortran name xifunc).
 *   c_check_inputs().............Check the input dimensions and variables (Fortran name chekin).
 *   c_dref().....................Flux albedo for given angle of incidence, given a bidirectional reflectivity.
 *   c_legendre_poly()............Compute the normalized associated Legendre polynomial, defined in terms of the
 *                                associated Legendre polynomial (Fortran name lepoly).
 *   c_planck_func1().............Compute Planck function integrated between two wavenumbers (Fortran name plkavg).
 *   c_planck_func2().............Compute Planck function integrated between two wavenumbers, or Planck function at
 *                                a specific wavenumber (Fortran name tplkavg).
 *   c_print_avg_intensities()....Print azimuthally averaged intensities at user angles (Fortran name pravin).
 *   c_print_inputs().............Print values of input variables (Fortran name prtinp).
 *   c_print_intensities()........Print the intensity at user polar and azimuthal angles (Fortran name prtint).
 *   c_gaussian_quadrature()......Compute weights and abscissae for ordinary Gaussian quadrature on the interval
 *                                (0,1); that is, such that sum(i=1 to M) ( GWT(i) f(GMU(i)) ) is a good
 *                                approximation to integral(0 to 1) (f(x) dx) (Fortran name qgausn).
 *   c_ratio()....................Returns ratio a/b with overflow and underflow protection, or 1.+a if b == 0.
 *   c_self_test()................Sets up self test and compares results (Fortran name slftst).
 *   c_albtrans().................DISORT special case to get only albedo and transmissivity of entire medium
 *                                as a function of incident beam angle (Fortran name albtrn).
 *   c_albtrans_intensity().......Computes azimuthally-averaged intensity at top and bottom of medium (related to
 *                                albedo and transmission of medium by reciprocity principles; see Ref S2; Fortran
 *                                name altrin).
 *   c_print_albtrans()...........Print planar albedo and transmissivity of medium as a function of incident beam
 *                                angle (Fortran name praltr).
 *   c_albtrans_spherical().......Calculate spherical albedo and transmissivity for the entire medium from the m=0
 *                                intensity components; a specialized version of fluxes (Fortran name spaltr).
 *   c_errmsg()...................Print out a warning or error message; abort if error.
 *   c_write_bad_var()............Write names of erroneous variable, keep count (Fortran name wrtbad, wrtbad2).
 *   c_write_too_small_dim()......Write name of too-small symbolic dimension and the value to which it should be
 *                                increased (Fortran name wrtdim, wrtdim2).
 *   c_sgbco()....................Factor a band matrix by Gaussian elimination and estimate the condition of the
 *                                matrix.
 *   c_sgbfa()....................Factor a band matrix by elimination.
 *   c_afval()....................Solve the band system A*X = B or transpose(A)*X = B using the factors computed by
 *                                sgbco() or sgbfa().
 *   c_sgeco()....................Factor a matrix by Gaussian elimination and estimate the condition of the matrix.
 *   c_sgefa()....................Factor a matrix by Gaussian elimination.
 *   c_sgesl()....................Solve the system A*X = B or transpose(A)*X = B using the factors computed by
 *                                sgeco() or sgefa().
 *   c_sasum()....................Sum of absolute values of elements in an array.
 *   c_saxpy()....................Compute A*X + Y given scalar A and vectors X and Y.
 *   c_sdot().....................Dot product (inner product) of input vectors X and Y.
 *   c_sscal()....................Multiply given scalar and vector.
 *   c_isamax()...................Return biggest absolute value among the array elements.
 *   c_twostr()...................Solve the radiative transfer equation in the two-stream approximation.
 *                                Based on the general-purpose algorithm DISORT, but both simplified and extended.
 *   c_chapman()..................Calculate the Chapman factor.
 *   c_twostr_check_inputs()......Check input dimensions and variables for two stream code (Fortran name tchekin).
 *   c_twostr_fluxes()............Calculate radiative fluxes, mean intensity, and flux derivative with respect to
 *                                optical depth from the azimuthally-averaged intensity (Fortran name tfluxes).
 *   c_twostr_solns().............Calculates the homogenous and particular solutions to the radiative transfer
 *                                equation in the two-stream approximation, for each layer in the medium (Fortran
 *                                name hopsol).
 *   c_twostr_print_inputs()......Print values of input variables (Fortran name tprtinp).
 *   c_twostr_set()...............Perform miscellaneous setting-up operations (Fortran name settwo).
 *
 * Functions added to C version:
 *
 *   c_fcmp().....................Safe floating-point comparison function (more reliable than ==)
 *   c_disort_state_alloc().......Dynamically allocate memory for disort input arrays, incl. ones the user can ask
 *                                disort() to calculate
 *   c_disort_state_free()........Free memory allocated by disort_state_alloc()
 *   c_disort_out_alloc_()........Dynamically allocate memory for disort output arrays
 *   c_disort_out_free()..........Free memory allocated by disort_out_alloc()
 *   c_twostr_state_alloc().......Dynamically allocate memory for twostr input arrays
 *   c_twostr_state_free(0........Free memory allocated by twostr_state_alloc()
 *   c_twostr_out_alloc().........Dynamically allocate memory for twostr output arrays
 *   c_twostr_out_free()..........Free memory allocated by twostr_out_alloc()
 *   c_dbl_vector()...............Allocate zeroed memory for double-precision vector of given range
 *   c_int_vector()...............Allocate zeroed memory for integer vector of given range
 *   c_free_dbl_vector()..........Free memory allocated by dbl_vector()
 *
 *  C rewrite by Timothy E. Dowling (Univ. of Louisville)
 *  This rewrite includes conversion to double precision, dynamic memory allocation, introduction of C structures
 *  (which may yield beneficial cache-aware memory allocation), and improved readability of subroutine names.
 *  new intensity correction added by Robert Buras (LMU Munich)
 */

#include "cdisort.h"
#include "locate.h"

/*============================= c_disort() ==============================*/

/*-------------------------------------------------------------------------------*
 * Plane-parallel discrete ordinates radiative transfer program                  *
 * C version                                                                     *
 * Fortran ftp site: ftp://climate.gsfc.nasa.gov/pub/wiscombe/Multiple_Scatt/    *
 *-------------------------------------------------------------------------------*

  Calling Tree (omitting calls to c_errmsg, c_dbl_vector, c_int_vector, c_free_dbl_vector):

  c_disort-+-c_self_test-+-c_disort_state_alloc
         |               +-c_disort_out_alloc
         |               +-c_disort_state_free
         |               +-c_disort_out_free
         +-c_check_inputs-+-(c_write_bad_var)
         |                +-c_dref
         +-c_disort_set-+-c_gaussian_quadrature
         +-c_print_inputs
         +-c_albtrans-+-c_legendre_poly
         |            +-c_solve_eigen-+-c_asymmetric_matrix
         |            +-c_interp_eigenvec
         |            +-c_set_matrix
         |            +-(c_sgbco)
         |            +-c_solve1-+-(c_sgbsl)
         |            +-c_atltrin
         |            +-c_albtrans_spherical
         |            +-c_print_albtrans
         +-c_planck_func1
         +-c_legendre_poly
         +-c_surface_bidir-+-c_gaussian_quadrature
         |                 +-c_bidir_reflectivity
         +-c_solve_eigen-+-c_asymmetric_matrix
	 +-c_set_coefficients_beam_source
	 +-c_interp_coefficients_beam_source
         +-c_upbeam_pseudo_spherical-+-(c_sgeco)
         |                           +-(c_sgesl)
         +-c_upbeam-+-(c_sgeco)
         |          +-(c_sgesl)
         +-c_upbeam_general_source-+-(c_sgeco)
         |                         +-(c_sgesl)
         +-c_upisot-+-(c_sgeco)
         |          +-(c_sgesl)
         +-c_interp_eigenvec
         +-c_interp_source
         +-c_set_matrix
         +-c_solve0-+-(c_sgbco)
         |          +-(c_sgbsl)
         +-c_fluxes
         +-c_user_intensities
         +-c_intensity_components
         +-c_print_avg_intensities
         +-c_ratio
         +-c_intensity_correction-+-c_single_scat
         |                        +-c_secondary_scat-+-c_xi_func
         +-c_new_intensity_correction-+-c_single_scat
         |                            +-prep_double_scat_integr
         |                            +-c_new_secondary_scat-+-c_xi_func
         |                            +-calc_phase_squared
         +-c_print_intensities

 +-------------------------------------------------------------------+

  Index conventions (for all loops and all variable descriptions):

  iu       :  for user polar angles
  iq,jq,kq :  for computational polar angles ('quadrature angles')
  iq/2     :  for half the computational polar angles (just the ones in either 0-90 degrees, or 90-180 degrees)
  j        :  for user azimuthal angles
  k,l      :  for Legendre expansion coefficients or, alternatively, subscripts of associated Legendre polynomials
  lu       :  for user levels
  lc       :  for computational layers (each having a different single-scatter albedo and/or phase function)
  lev      :  for computational levels
  mazim    :  for azimuthal components in Fourier cosine expansion of intensity and phase function

 +------------------------------------------------------------------+

               I N T E R N A L    V A R I A B L E S

   AMB(iq/2,iq/2)....First  matrix factor in reduced eigenvalue problem of eqs. SS(12), STWJ(8E), STWL(23f) (used only in solve_eigen);
                     ab[].zero (see cdisort.h)
   APB(iq/2,iq/2)....Second matrix factor in reduced eigenvalue problem of eqs. SS(12), STWJ(8E), STWL(23f) (used only in solve_eigen);
                     ab[].one (see cdisort.h)
   ARRAY(iq,iq)......Scratch matrix for solve_eigen(), upbeam() and upisot()
                     (see each subroutine for definition)
   B()...............Right-hand side vector of eq. SC(5) going into SOLVE0,1;
                     returns as solution vector vector  L, the constants of integration
   BDR(iq/2,0:iq/2)..Bottom-boundary bidirectional reflectivity for a given azimuthal component.  First index always
                     refers to a computational angle.  Second index: if zero, refers to incident beam angle UMU0;
                     if non-zero, refers to a computational angle.
   BEM(iq/2).........Bottom-boundary directional emissivity at computational angles.
   bplanck...........Intensity emitted from bottom boundary
   callnum...........Number of surface calls
   CBAND()...........Matrix of left-hand side of the linear system eq. SC(5), scaled by eq. SC(12);
                     in banded form required by LINPACK solution routines
   CC(iq,iq).........C-sub-IJ in eq. SS(5)
   CH(lc)............The Chapman-factor to correct for pseudo-spherical geometry in the direct beam.
   CHTAU(lc).........The optical depth in spherical geometry.
   CMU(iq)...........Computational polar angles (Gaussian)
   CWT(iq)...........Quadrature weights corresponding to CMU
   corint............When set TRUE, correct intensities for delta-scaling effects (see Nakajima and Tanaka, 1988).
                     When FALSE, intensities are not corrected. In general, CORINT should be set true when beam
                     source is present (FBEAM is not zero) and DELTAM is TRUE in a problem including scattering.
                     However, execution is faster when CORINT is FALSE, and intensities outside the aureole may still be
                     accurate enough.  When CORINT is TRUE, it is important to have a sufficiently high order of
                     Legendre approximation of the phase function. This is because the intensities are corrected by
                     calculating the single-scattered radiation, for which an adequate representation of the phase
                     function is crucial.  In case of a low order Legendre approximation of an otherwise highly
                     anisotropic phase function, the intensities might actually be more accurate when corint is FALSE.
                     When only fluxes are calculated (ds->flag.onlyfl is TRUE), or there is no beam source (FBEAM=0.0), or there
                     is no scattering (SSALB =0. for all layers) corint is set FALSE by the code.
   delm0.............Kronecker delta, delta-sub-M0, where M = MAZIM is the number of the Fourier component in the
                     azimuth cosine expansion
   deltam............TRUE,  use delta-M method ( see Wiscombe, 1977 );
                     FALSE, do not use delta-M method. 
                     In general, for a given number of streams, intensities and fluxes will be more accurate for phase functions
                     with a large forward peak if DELTAM is set true. Intensities close to the forward scattering
                     direction are often less accurate, however, when the delta-M method is applied. The intensity deltam
                     correction of Nakajima and Tanaka is used to improve the accuracy of the intensities.
   dither............Small quantity subtracted from single-scattering albedos of unity, in order to avoid using special
                     case formulas;  prevents an eigenvalue of exactly zero from occurring, which would cause an immediate overflow
   DTAUCPR(lc).......Computational-layer optical depths (delta-M-scaled if DELTAM = TRUE, otherwise equal to DTAUC)
   EMU(iu)...........Bottom-boundary directional emissivity at user angles.
   EVAL(iq)..........Temporary storage for eigenvalues of eq. SS(12)
   EVECC(iq,iq)......Complete eigenvectors of SS(7) on return from solve_eigen; stored permanently in  GC
   EXPBEA(lc)........Transmission of direct beam in delta-M optical depth coordinates
   FLYR(lc)..........Separated fraction in delta-M method
   GL(k,lc)..........Phase function Legendre polynomial expansion coefficients, calculated from PMOM by
                     including single-scattering albedo, factor 2K+1, and (if DELTAM=TRUE) the delta-M scaling
   GC(iq,iq,lc)......Eigenvectors at polar quadrature angles, g in eq. SC(1)
   GU(iu,iq,lc)......Eigenvectors interpolated to user polar angles (g  in eqs. SC(3) and S1(8-9), i.e. g without the l factor)
   IPVT(lc*iq).......Integer vector of pivot indices for LINPACK routines
   KK(iq,lc).........Eigenvalues of coeff. matrix in eq. SS(7)
   kconv.............Counter in azimuth convergence test
   LAYRU(lu).........Computational layer in which user output level UTAU(LU) is located
   LL(iq,lc).........Constants of integration L in eq. SC(1), obtained by solving scaled version of eq. SC(5)
   lyrcut............TRUE, radiation is assumed zero below layer ncut because of almost complete absorption
   naz...............Number of azimuthal components considered
   ncut..............Computational layer number in which absorption optical depth first exceeds ABSCUT
   OPRIM(lc).........Single scattering albedo after delta-M scaling
   pass1.............TRUE on first entry, FALSE thereafter
   PKAG(0:lc)........Integrated Planck function for internal emission
   PRNTU0(l).........logical flag to trigger printing of azimuthally-averaged intensities:
                       l    quantities printed
                      --    ------------------
                       0    azimuthally-averaged intensities at user
                               levels and computational polar angles
                       1    azimuthally-averaged intensities at user
                               levels and user polar angles
   PSI0(iq)..........Sum just after square bracket in  eq. SD(9); psi[].zero (see cdisort.h)
   PSI1(iq)..........Sum in  eq. STWL(31d); psi[].one
   RMU(iu,0:iq)......Bottom-boundary bidirectional reflectivity for a given azimuthal component.  First index always
                     refers to a user angle.  Second index: if zero, refers to incident beam angle UMU0;
                     if non-zero, refers to a computational angle.
   scat_yes..........int, TRUE if scattering, FALSE if not (added to C version)
   TAUC(0:lc)........Cumulative optical depth (un-delta-M-scaled)
   TAUCPR(0:lc)......Cumulative optical depth (delta-M-scaled if DELTAM = TRUE, otherwise equal to TAUC)
   tplanck...........Intensity emitted from top boundary
   UUM(iu,lu)........Expansion coefficients when the intensity (u-super-M) is expanded in Fourier cosine series
                     in azimuth angle
   U0C(iq,lu)........Azimuthally-averaged intensity at quadrature angle
   U0U(iu,lu)........If ds->flag.onlyfl = FALSE, azimuthally-averaged intensity at user angles and user levels
                     If ds->flag.onlyfl = TRUE, azimuthally-averaged intensity at computational
                     (Gaussian quadrature) angles and user levels; the corresponding quadrature angle cosines are
                     returned in UMU.
   UTAUPR(lu)........Optical depths of user output levels in delta-M coordinates; equal to UTAU(LU) if no delta-M
   WK(iq)............Scratch array
   XR0(lc)...........X-sub-zero in expansion of thermal source function preceding eq. SS(14)(has no mu-dependence); b-sub-zero in eq. STWL(24d)
   XR1(lc)...........X-sub-one in expansion of thermal source function; see eqs. SS(14-16); b-sub-one in STWL(24d)
   YLM0(l)...........Normalized associated Legendre polynomial of subscript L at the beam angle (not saved
                     as function of superscipt M)
   YLMC(l,iq)........Normalized associated Legendre polynomial of subscript L at the computational angles
                     (not saved as function of superscipt M)
   YLMU(l,iu)........Normalized associated Legendre polynomial of subscript L at the user angles
                     (not saved as function of superscipt M)
   Z()...............scratch array used in solve0(), albtrans() to solve a linear system for the constants of integration
   Z0(iq)............Solution vectors Z-sub-zero of eq. SS(16); zee[].zero (see cdisort.h)
   Z1(iq)............Solution vectors Z-sub-one  of eq. SS(16); zee[].one
   Z0U(iu,lc)........Z-sub-zero in eq. SS(16) interpolated to user angles from an equation derived from SS(16); zu[].zero (see cdisort.h)
   Z1U(iu,lc)........Z-sub-one  in eq. SS(16) interpolated to user angles from an equation derived from SS(16); zu[].one
   ZBEAM(iu,lc)......Particular solution for beam source
   ZGU(iu,lc)........General source function interpolated to user angles
   ZJ(iq)............Right-hand side vector  X-sub-zero in eq. SS(19), also the solution vector
                     Z-sub-zero after solving that system
   ZJG(iq)...........Right-hand side vector  X-sub-zero in eq. KS(10), also the solution vector
                     Z-sub-zero after solving that system for a general source constant over a layer
   ZZ(iq,lc).........Permanent storage for the beam source vectors ZJ
   ZZG(iq,lc)........Permanent storage for the beam source vectors ZJG
   ZPLK0(iq,lc)......Permanent storage for the thermal source vectors plk[].zero obtained by solving eq. SS(16)
   ZPLK1(iq,lc)......Permanent storage for the thermal source vectors plk[].one  obtained by solving eq. SS(16)

*/ 

void c_disort(disort_state  *ds,
	      disort_output *out)
{
  static int
    self_tested = -1;
  int
    prntu0[2],
    corint,deltam,scat_yes,compare,lyrcut,needdeltam,
    iq,iu,j,kconv,l,lc,lev,lu,mazim,naz,ncol,ncos,ncut,nn;
  static int
    callnum=1;
  int
    ipvt[ds->nstr*ds->nlyr],
    layru[ds->ntau];
  double
    angcos,azerr,azterm,bplanck,cosphi,delm0,
    sgn,tplanck;
  double
    *array,*b,*bdr,*bem,*cband,*cc,*ch,*chtau,
    *cmu,*cwt, *dtaucpr,*emu,*eval,*evecc,*expbea,
    *flyr,*gc,*gl,*gu,*kk,*ll,
    *oprim,*phasa,*phast,*phasm,*phirad,*pkag,
    *rmu,*tauc,*taucpr,*u0c,*utaupr,*uum,
    *wk,*xba,*ylm0,*ylmc,*ylmu,
    *z,*zbeam,
    *zbeama,zbsa=0,*zj,*zjg,*zju,*zgu,*zz,*zzg;
  disort_pair
    *ab,*fl,*plk,*xr,*psi,*xb,*zbeamsp,*zbs,*zee,*zu;
  disort_triplet
    *zbu;
  const double
    dither = 100.*DBL_EPSILON;


  /* Set these here to ensure that memory is correctly allocated. */
  if (!ds->flag.usrtau) {
    ds->ntau = ds->nlyr+1;
  }
  if ( (!ds->flag.usrang || ds->flag.onlyfl)  && ( !ds->flag.ibcnd == SPECIAL_BC)) {
    ds->numu = ds->nstr;
  }
  if (ds->flag.usrang && ds->flag.ibcnd == SPECIAL_BC) {
    ds->numu *= 2;
  }

  if (self_tested == -1) {
    int
      prntu0_test[2] = {FALSE,FALSE};
    disort_state
      ds_test;
    disort_output
      out_test;
    /*
     * Set input values for self-test.
     * Be sure self_test() sets all print flags off.
     */
    self_tested = 0;
    compare     = FALSE;
    c_self_test(compare,prntu0_test,&ds_test,&out_test);
    c_disort(&ds_test,&out_test);
  }

  /*
   * Determine whether there is scattering or not
   */
  scat_yes = FALSE;
  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (SSALB(lc) > 0.) {
      scat_yes = TRUE;
      break;
    }
  }

  /* 
   * Turn on delta-M tranformation
   */
  deltam = TRUE;

  /* delta-M scaling makes only sense if phase function has more
   * moments than streams
   */
  needdeltam = FALSE;
  if( deltam==TRUE ) {
    for (lc=1; lc<=ds->nlyr; lc++)
      if ( PMOM(ds->nstr,lc) != 0.0 )
	needdeltam = TRUE;
    if (needdeltam==FALSE)
      deltam=FALSE;
  }

  /*
   * Turn off intensity correction when only fluxes are calculated, there
   * is no beam source, no scattering, or delta-M transformation is not applied
   */
  corint = ds->flag.intensity_correction;
  if (ds->flag.onlyfl || ds->bc.fbeam == 0. || !scat_yes || !deltam) 
    corint = FALSE;

  prntu0[0] = FALSE;
  prntu0[1] = FALSE;

  /*
   * Allocate zeroed memory
   */
  tauc = c_dbl_vector(0,ds->nlyr,"tauc");

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if(SSALB(lc) == 1.) {
      SSALB(lc) = 1.-dither;
    }
    TAUC(lc) = TAUC(lc-1)+DTAUC(lc);
  }

  /* Check input dimensions and variables */
  c_check_inputs(ds,scat_yes,deltam,corint,tauc,callnum);

  /*-------------------------------------------------------------------------------------------*
   * Special case for getting albedo and transmissivity of medium for many beam angles at once *
   *-------------------------------------------------------------------------------------------*/

  if (ds->flag.ibcnd == SPECIAL_BC) {
    /* 
     * Allocate zeroed memory
     */
    array    = c_dbl_vector(0,ds->nstr*ds->nstr-1,"array");
    b        = c_dbl_vector(0,ds->nstr*ds->nlyr-1,"b");
    bdr      = c_dbl_vector(0,((ds->nstr/2)+1)*(ds->nstr/2)-1,"bdr");
    cband    = c_dbl_vector(0,ds->nstr*ds->nlyr*(9*(ds->nstr/2)-2)-1,"cband");
    ch       = c_dbl_vector(0,ds->nlyr-1,"ch");
    chtau    = c_dbl_vector(0,(2*ds->nlyr+1)-1,"chtau");
    cc       = c_dbl_vector(0,ds->nstr*ds->nstr-1,"cc");
    cmu      = c_dbl_vector(0,ds->nstr-1,"cmu");
    cwt      = c_dbl_vector(0,ds->nstr-1,"cwt");
    dtaucpr  = c_dbl_vector(0,ds->nlyr-1,"dtaucpr");
    eval     = c_dbl_vector(0,(ds->nstr/2)-1,"eval");
    evecc    = c_dbl_vector(0,ds->nstr*ds->nstr-1,"evecc");
    expbea   = c_dbl_vector(0,ds->nlyr,"expbea");
    flyr     = c_dbl_vector(0,ds->nlyr-1,"flyr");
    gc       = c_dbl_vector(0,ds->nlyr*ds->nstr*ds->nstr-1,"gc");
    gl       = c_dbl_vector(0,ds->nlyr*(ds->nstr+1),"gl");
    gu       = c_dbl_vector(0,ds->nlyr*ds->nstr*ds->numu-1,"gu");
    kk       = c_dbl_vector(0,ds->nlyr*ds->nstr-1,"kk");
    ll       = c_dbl_vector(0,ds->nlyr*ds->nstr-1,"ll");
    oprim    = c_dbl_vector(0,ds->nlyr-1,"oprim");
    taucpr   = c_dbl_vector(0,ds->nlyr,"taucpr");
    utaupr   = c_dbl_vector(0,ds->ntau-1,"utaupr");
    wk       = c_dbl_vector(0,ds->nstr-1,"wk");
    ylmc     = c_dbl_vector(0,ds->nstr*(ds->nstr+1)-1,"ylmc");
    ylmu     = c_dbl_vector(0,(ds->numu)*(ds->nstr+1)-1,"ylmu");
    z        = c_dbl_vector(0,ds->nstr*ds->nlyr-1,"z");

    ab       = (disort_pair *)calloc((ds->nstr/2)*(ds->nstr/2),sizeof(disort_pair));
    if (!ab) {
      c_errmsg("disort alloc error for ab", DS_ERROR);
    }
    /*
     * Zero output arrays
     */
    if (!ds->flag.usrtau) {
      memset(ds->utau,0,ds->ntau*sizeof(double));
    }
    if (!ds->flag.usrang || ds->flag.onlyfl) {
      memset(ds->umu,0,(ds->numu+1)*sizeof(double));
    }
    memset(out->rad,   0,ds->ntau*sizeof(disort_radiant));
    memset(out->albmed,0,ds->numu*sizeof(double));
    memset(out->trnmed,0,ds->numu*sizeof(double));
    if (ds->flag.onlyfl == FALSE) {
      memset(out->uu,0,ds->numu*ds->ntau*ds->nphi*sizeof(double));
    }

    /* Perform various setup operations */
    c_disort_set(ds,ch,chtau,cmu,cwt,deltam,dtaucpr,expbea,flyr,gl,layru,&lyrcut,&ncut,&nn,&corint,oprim,tauc,taucpr,utaupr);

    /*  Print input information */
    if(ds->flag.prnt[0]) {
      c_print_inputs(ds,dtaucpr,scat_yes,deltam,corint,flyr,lyrcut,oprim,tauc,taucpr);
    }

    c_albtrans(ds,out,ab,array,b,bdr,cband,cc,cmu,cwt,dtaucpr,eval,evecc,gl,gc,gu,ipvt,kk,ll,nn,taucpr,ylmc,ylmu,z,wk);

    callnum++;

    /*
     * Free allocated memory
     */
    free(array), free(b),    free(bdr), free(cband),  free(cc),    free(ch);
    free(chtau), free(cmu),  free(cwt), free(dtaucpr),free(eval),  free(evecc);
    free(expbea),free(flyr), free(gc),  free(gl),     free(gu),    free(kk);
    free(ll),    free(oprim),free(tauc),free(taucpr), free(utaupr),free(wk);
    free(ylmc),  free(ylmu), free(z),   free(ab);

    return;
  }

  /*--------------*
   * General case *
   *--------------*/

  /* 
   * Allocate zeroed memory
   */
  array   = c_dbl_vector(0,ds->nstr*ds->nstr-1,"array");
  b       = c_dbl_vector(0,ds->nstr*ds->nlyr-1,"b");
  bdr     = c_dbl_vector(0,((ds->nstr/2)+1)*(ds->nstr/2)-1,"bdr");
  bem     = c_dbl_vector(0,(ds->nstr/2)-1,"bem");
  cband   = c_dbl_vector(0,ds->nstr*ds->nlyr*(9*(ds->nstr/2)-2)-1,"cband");
  cc      = c_dbl_vector(0,ds->nstr*ds->nstr-1,"cc");
  ch      = c_dbl_vector(0,ds->nlyr-1,"ch");
  chtau   = c_dbl_vector(0,(2*ds->nlyr+1)-1,"chtau");
  cmu     = c_dbl_vector(0,ds->nstr-1,"cmu");
  cwt     = c_dbl_vector(0,ds->nstr-1,"cwt");
  dtaucpr = c_dbl_vector(0,ds->nlyr-1,"dtaucpr");
  emu     = c_dbl_vector(0,ds->numu-1,"emu");
  eval    = c_dbl_vector(0,(ds->nstr/2)-1,"eval");
  evecc   = c_dbl_vector(0,ds->nstr*ds->nstr-1,"evecc");
  expbea  = c_dbl_vector(0,ds->nlyr,"expbea");
  flyr    = c_dbl_vector(0,ds->nlyr,"flyr");    // We need at least one element
  gc      = c_dbl_vector(0,ds->nlyr*ds->nstr*ds->nstr-1,"gc");
  gl      = c_dbl_vector(0,ds->nlyr*(ds->nstr+1),"gl");
  gu      = c_dbl_vector(0,ds->nlyr*ds->nstr*ds->numu-1,"gu");
  kk      = c_dbl_vector(0,ds->nlyr*ds->nstr-1,"kk");
  ll      = c_dbl_vector(0,ds->nlyr*ds->nstr-1,"ll");
  oprim   = c_dbl_vector(0,ds->nlyr-1,"oprim");
  phasa   = c_dbl_vector(0,ds->nlyr-1,"phasa");
  phast   = c_dbl_vector(0,ds->nlyr-1,"phast");
  phasm   = c_dbl_vector(0,ds->nlyr-1,"phasm");
  if (ds->nphi > 0) {
    phirad = c_dbl_vector(0,ds->nphi-1,"phirad");
  }
  else {
    phirad = NULL;
  }
  pkag   = c_dbl_vector(0,ds->nlyr,"pkag");
  rmu    = c_dbl_vector(0,((ds->nstr/2)+1)*ds->numu-1,"rmu");
  taucpr = c_dbl_vector(0,ds->nlyr,"taucpr");
  u0c    = c_dbl_vector(0,ds->ntau*ds->nstr-1,"u0c");
  utaupr = c_dbl_vector(0,ds->ntau-1,"utaupr");
  uum    = c_dbl_vector(0,ds->ntau*ds->numu-1,"uum");
  wk     = c_dbl_vector(0,ds->nstr-1,"wk");
  xba    = c_dbl_vector(0,ds->nlyr,"xba");
  ylm0   = c_dbl_vector(0,ds->nstr,"ylm0");
  ylmc   = c_dbl_vector(0,ds->nstr*(ds->nstr+1)-1,"ylmc");
  ylmu   = c_dbl_vector(0,ds->numu*(ds->nstr+1)-1,"ylmu");
  z      = c_dbl_vector(0,ds->nstr*ds->nlyr-1,"z");
  zbeam  = c_dbl_vector(0,ds->nlyr*ds->numu-1,"zbeam");
  zbeama = c_dbl_vector(0,ds->nlyr,"zbeama");
  zj     = c_dbl_vector(0,ds->nstr-1,"zj");
  zjg    = c_dbl_vector(0,ds->nstr-1,"zjg");
  zju    = c_dbl_vector(0,ds->numu,"zju");
  zgu    = c_dbl_vector(0,ds->nlyr*ds->numu-1,"zgu");
  zz     = c_dbl_vector(0,ds->nlyr*ds->nstr-1,"zz");
  zzg    = c_dbl_vector(0,ds->nlyr*ds->nstr-1,"zzg");
  /*
   * Using C structures facilitates cache-aware memory allocation, which can reduce
   * cache misses and potentially speed up computer execution.
   */
  fl     = (disort_pair *)calloc(ds->ntau,sizeof(disort_pair));                  if (!fl)      c_errmsg("disort alloc error for fl", DS_ERROR);
  plk    = (disort_pair *)calloc(ds->nlyr*ds->nstr,sizeof(disort_pair));         if (!plk)     c_errmsg("disort alloc error for plk",DS_ERROR);
  ab     = (disort_pair *)calloc((ds->nstr/2)*(ds->nstr/2),sizeof(disort_pair)); if (!ab)      c_errmsg("disort alloc error for ab", DS_ERROR);
  xr     = (disort_pair *)calloc(ds->nlyr,sizeof(disort_pair));                  if (!xr)      c_errmsg("disort alloc error for xr", DS_ERROR);
  psi    = (disort_pair *)calloc(ds->nstr,sizeof(disort_pair));                  if (!psi)     c_errmsg("disort alloc error for psi",DS_ERROR);
  xb     = (disort_pair *)calloc(ds->nlyr*ds->nstr,sizeof(disort_pair));         if (!xb)      c_errmsg("disort alloc error for xb",DS_ERROR);
  zbs    = (disort_pair *)calloc(ds->nstr,sizeof(disort_pair));                  if (!zbs)     c_errmsg("disort alloc error for zbs",DS_ERROR);
  zbeamsp= (disort_pair *)calloc(ds->nlyr*ds->nstr,sizeof(disort_pair));         if (!zbeamsp) c_errmsg("disort alloc error for zbeamsp",DS_ERROR);
  zee    = (disort_pair *)calloc(ds->nstr,sizeof(disort_pair));                  if (!zee)     c_errmsg("disort alloc error for zee",DS_ERROR);
  zu     = (disort_pair *)calloc(ds->nlyr*ds->numu,sizeof(disort_pair));         if (!zu)      c_errmsg("disort alloc error for zu", DS_ERROR);

  zbu    = (disort_triplet *)calloc(ds->nlyr*ds->numu,sizeof(disort_triplet));   if (!zbu)     c_errmsg("disort alloc error for zbu", DS_ERROR);

  /*
   * Zero output arrays
   */ 
  if (!ds->flag.usrtau) {
    memset(ds->utau,0,ds->ntau*sizeof(double));
  }
  if (!ds->flag.usrang || ds->flag.onlyfl) {
    memset(ds->umu,0,(ds->numu)*sizeof(double));
  }
  memset(out->rad,0,ds->ntau*sizeof(disort_radiant));
  if (ds->flag.onlyfl == FALSE) {
    memset(out->uu,0,ds->numu*ds->ntau*ds->nphi*sizeof(double));
  }

  /* Perform various setup operations */
  c_disort_set(ds,ch,chtau,cmu,cwt,deltam,dtaucpr,expbea,flyr,gl,layru,&lyrcut,&ncut,&nn,&corint,oprim,tauc,taucpr,utaupr);


  /*  Print input information */
  if(ds->flag.prnt[0]) {
    c_print_inputs(ds,dtaucpr,scat_yes,deltam,corint,flyr,lyrcut,oprim,tauc,taucpr);
  }

  /*
   * Calculate Planck functions
   */
  if (!ds->flag.planck) {
    bplanck = 0.;
    tplanck = 0.;
  }
  else {
    tplanck = c_planck_func1(ds->wvnmlo,ds->wvnmhi,ds->bc.ttemp)*ds->bc.temis;
    bplanck = c_planck_func1(ds->wvnmlo,ds->wvnmhi,ds->bc.btemp);
    for (lev = 0; lev <= ds->nlyr; lev++) {
      PKAG(lev) = c_planck_func1(ds->wvnmlo,ds->wvnmhi,TEMPER(lev));
    }
  }

  /*
   *--------  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  ---------
   *          (eq STWJ 5, STWL 6)
   */
  kconv = 0;
  naz   = ds->nstr-1;

  /*
   * Azimuth-independent case
   */
  if (ds->bc.fbeam == 0.                         || 
      fabs(1.-ds->bc.umu0) < 1.e-5               || 
      ds->flag.onlyfl                            || 
      (ds->numu == 1 && fabs(1.-UMU(1)) < 1.e-5) || 
      (ds->numu == 1 && fabs(1.+UMU(1)) < 1.e-5) || 
      (ds->numu == 2 && fabs(1.+UMU(1)) < 1.e-5 && fabs(1.-UMU(2)) < 1.e-5)) {
    naz = 0;
  }

  for (mazim = 0; mazim <= naz; mazim++) {
    if (mazim == 0) {
      delm0 = 1.;
    }
    else {
      delm0 = 0.;
    }

    /*
     * Get normalized associated Legendre polynomials for
     *   (a) incident beam angle cosine
     *   (b) computational and user polar angle cosines
     */
    if (ds->bc.fbeam > 0.) {
      ncos   = 1;
      angcos = -ds->bc.umu0;
      c_legendre_poly(ncos,mazim,ds->nstr,ds->nstr-1,&angcos,ylm0);
    }

    if (!ds->flag.onlyfl && ds->flag.usrang) {
      c_legendre_poly(ds->numu,mazim,ds->nstr,ds->nstr-1,ds->umu,ylmu);
    }
    c_legendre_poly(nn,mazim,ds->nstr,ds->nstr-1,cmu,ylmc);

    /*
     * Get normalized associated Legendre polynomials with negative arguments from those with
     * positive arguments; Dave/Armstrong eq. (15), STWL(59)
     */
    sgn = -1.;
    for (l = mazim; l <= ds->nstr-1; l++) {
      sgn *= -1.;
      for (iq = nn+1; iq <= ds->nstr; iq++) {
        YLMC(l,iq) = sgn*YLMC(l,iq-nn);
      }
    }

    /*
     * Specify users bottom reflectivity and emissivity properties
     */
    if (!lyrcut) {
      c_surface_bidir(ds, delm0, cmu, mazim, nn, bdr, emu, bem, rmu,
		      callnum);
    }

    /*--------------  BEGIN LOOP ON COMPUTATIONAL LAYERS  ------------*/
    for (lc = 1; lc <= ncut; lc++) {
      /*
       * Solve eigenfunction problem in eq. STWJ(8B), STWL(23f); return eigenvalues and eigenvectors
       */
      c_solve_eigen(ds,lc,ab,array,cmu,cwt,gl,mazim,nn,ylmc,cc,evecc,eval,kk,gc,wk);
      /*
       * Calculate particular solutions of eq. SS(18), STWL(24a) for incident beam source
       */
      if (ds->bc.fbeam > 0.) {
	if ( ds->flag.spher == TRUE ) {
	  /* Pseudo-spherical approach */
	  c_set_coefficients_beam_source(ds,ch,chtau,cmu,delm0,ds->bc.fbeam,
					 gl,lc,mazim,ds->nstr,
					 taucpr,xba,xb,ylm0,ylmc,zj);
	  
	  if ( ds->flag.usrang == TRUE  ) {
	    /* Get coefficients at umu for pseudo-spherical source */
	    c_interp_coefficients_beam_source(ds,chtau,delm0,ds->bc.fbeam,
					      gl,lc,mazim,ds->nstr,
					      ds->numu,taucpr,zbu,
					      xba,zju,ylm0,ylmu);
	  }
	  c_upbeam_pseudo_spherical(ds,lc,array,cc,cmu,ipvt,nn,wk,
				    xb,xba,zbs,&zbsa,zbeamsp,zbeama);
	}
	else {
	  /* Plane-parallel version */
	  c_upbeam(ds,lc,array,cc,cmu,delm0,gl,ipvt,mazim,nn,wk,ylm0,ylmc,zj,zz);
	}
      }

      /*
       * Calculate particular solutions of eq. SS(18), STWL(24a), KS(5) for 
       * general user specified source.
       */
      if (ds->flag.general_source) {
	c_upbeam_general_source(ds,lc,mazim,array,cc,ipvt,nn,wk,zjg,zzg);
      }

      /*
       * Calculate particular solutions of eq. SS(15), STWL(25) for thermal emission source
       */
      if (ds->flag.planck && mazim == 0) {
        XR1(lc) = 0.;
        if (DTAUCPR(lc) > 1e-4) { /* fix by RPB, caused problems in make check AVHRR CH4/5 */
          XR1(lc) = (PKAG(lc)-PKAG(lc-1))/DTAUCPR(lc);
        }
        XR0(lc) = PKAG(lc-1)-XR1(lc)*TAUCPR(lc-1);
        c_upisot(ds,lc,array,cc,cmu,ipvt,nn,oprim,wk,xr,zee,plk);
      }

      if (!ds->flag.onlyfl && ds->flag.usrang) {
        /*
         * Interpolate eigenvectors to user angles
         */
        c_interp_eigenvec(ds,lc,cwt,evecc,gl,gu,mazim,nn,wk,ylmc,ylmu);
        /*
         * Interpolate source terms to user angles
         */
        c_interp_source(ds,lc,cwt,delm0,gl,mazim,oprim,ylm0,ylmc,ylmu,
			psi,xr,zee,zj,zjg,zbeam,zbu,zbs,zbsa,zgu,zu);
      }
    }
    /*-------------------  END LOOP ON COMPUTATIONAL LAYERS  ----------------*/

    /*
     *
     * Set coefficient matrix of equations combining boundary and layer interface conditions
     */
    c_set_matrix(ds,bdr,cband,cmu,cwt,delm0,dtaucpr,gc,kk,lyrcut,&ncol,ncut,taucpr,wk);

    /*
     * Solve for constants of integration in homogeneous solution (general boundary conditions)
     */
    c_solve0(ds,b,bdr,bem,bplanck,cband,cmu,cwt,expbea,ipvt,ll,lyrcut,
	     mazim,ncol,ncut,nn,tplanck,taucpr,z,zbeamsp,zbeama,zz,zzg,plk);

    /*
     * Compute upward and downward fluxes
     */
    if (mazim == 0) {
      c_fluxes(ds,out,ch,cmu,cwt,gc,kk,layru,ll,lyrcut,ncut,nn,PRNTU0(1),
	       taucpr,utaupr,xr,zbeamsp,zbeama,zz,zzg,plk,fl,u0c);
    }

    if (ds->flag.onlyfl) {
      /*
       * Save azimuthal-avg intensities at quadrature angles
       */
      for (lu = 1; lu <= ds->ntau; lu++) {
        for (iq = 1; iq <= ds->nstr; iq++) {
          U0U(iq,lu) = U0C(iq,lu);
        }
      }
      break;
    }

    memset(uum,0,ds->numu*ds->ntau*sizeof(double));

    if (ds->flag.usrang) {
      /*
       * Compute azimuthal intensity components at user angles
       */
      c_user_intensities(ds,bplanck,cmu,cwt,delm0,dtaucpr,emu,expbea,
			 gc,gu,kk,layru,ll,lyrcut,mazim,
			 ncut,nn,rmu,taucpr,tplanck,utaupr,wk,
			 zbu,zbeam,zbeamsp,
			 zbeama,zgu,zu,zz,zzg,plk,uum);
    }
    else {
      /*
       * Compute azimuthal intensity components at quadrature angles
       */
      c_intensity_components(ds,gc,kk,layru,ll,lyrcut,mazim,ncut,nn,taucpr,utaupr,zz,plk,uum);
    }

    if (mazim == 0) {
      /*
       * Save azimuthally averaged intensities
       */
      for (lu = 1; lu <= ds->ntau; lu++) {
        for (iu = 1; iu <= ds->numu; iu++) {
          U0U(iu,lu) = UUM(iu,lu);
          for (j = 1; j <= ds->nphi; j++) {
            UU(iu,lu,j) = UUM(iu,lu);
          }
        }
      }

      if ( ds->flag.output_uum)
	for (lu = 1; lu <= ds->ntau; lu++) 
	  for (iu = 1; iu <= ds->numu; iu++)
            OUT_UUM(iu,lu,mazim) = UUM(iu,lu);
      
      /*
       * Print azimuthally averaged intensities at user angles
       */
      if (PRNTU0(2)) {
        c_print_avg_intensities(ds,out);
      }

      if (naz > 0) {
        memset(phirad,0,ds->nphi*sizeof(double));
        for (j = 1; j <= ds->nphi; j++) {
          PHIRAD(j) = (PHI(j)-ds->bc.phi0)*DEG;
        }
      }
    }
    else {
      /*
       * Increment intensity by current azimuthal component (Fourier cosine series);  eq SD(2), STWL(6)
       */
      azerr = 0.;
      for (j = 1; j <= ds->nphi; j++) {
        cosphi = cos((double)mazim*PHIRAD(j));
        for (lu = 1; lu <= ds->ntau; lu++) {
          for (iu = 1; iu <= ds->numu; iu++) {
            azterm       = UUM(iu,lu)*cosphi;
            UU(iu,lu,j) += azterm;
            azerr        = MAX(azerr,c_ratio(fabs(azterm),fabs(UU(iu,lu,j))));
          }
        }
      }
      if ( ds->flag.output_uum)
	for (lu = 1; lu <= ds->ntau; lu++) 
	  for (iu = 1; iu <= ds->numu; iu++) 
            OUT_UUM(iu,lu,mazim) = UUM(iu,lu);

      if(azerr <= ds->accur) {
        kconv++;
      }
      if (kconv >= 2) {
        break;
      }
    }
  }
  /*--------------  END LOOP ON AZIMUTHAL COMPONENTS  ----------------*/



  for (iu = 1; iu <= ds->numu; iu++) {
    lu = ds->ntau;
    j =  1;
  }
  if (corint) {
    /*
     * Apply Nakajima/Tanaka intensity corrections
     */
    if (!ds->flag.old_intensity_correction && self_tested == 1) {
      if (ds->flag.quiet==VERBOSE)
	fprintf(stderr,"Using new intensity correction, with phase functions\n");
      c_new_intensity_correction(ds,out,dither,flyr,layru,lyrcut,ncut,oprim,phasa,phast,phasm,phirad,tauc,taucpr,utaupr);
    }
    else {
      if (ds->flag.quiet==VERBOSE) 
	fprintf(stderr,"Using original intensity correction, with phase moments\n");
      c_intensity_correction(ds,out,dither,flyr,layru,lyrcut,ncut,oprim,phasa,phast,phasm,phirad,tauc,taucpr,utaupr);
    }
  }

  
  for (iu = 1; iu <= ds->numu; iu++) {
    lu = ds->ntau;
    j =  1;
  }


  if (ds->flag.prnt[2] && !ds->flag.onlyfl) {
    /*
     * Print intensities
     */
    c_print_intensities(ds,out);
  }

  if (self_tested == 0) {
    /*
     * Compare test case results with correct answers and abort if bad
     */
    compare = TRUE;
    c_self_test(compare,prntu0,ds,out);

    self_tested = 1;
  }

  callnum++;

  /* 
   * Free allocated memory
   */
  free(ab),free(array);
  free(b),free(bdr),free(bem);
  free(cband),free(cc),free(ch),free(chtau),free(cmu),free(cwt);
  free(dtaucpr);
  free(emu),free(eval),free(evecc),free(expbea);
  free(flyr),free(fl);
  free(gc),free(gl),free(gu);
  free(kk);
  free(ll);
  free(oprim);
  free(phasa),free(phast),free(phasm),free(phirad),free(pkag),free(plk),free(psi);
  free(rmu);
  free(tauc),free(taucpr);
  free(u0c),free(utaupr),free(uum);
  free(wk);
  free(xb),free(xba),free(xr);
  free(ylm0),free(ylmc),free(ylmu);
  free(z),free(zbu),free(zbeam),free(zbeamsp),
  free(zbeama),free(zbs),free(zj),free(zjg),free(zju),
  free(zgu),free(zz),free(zzg),free(zee),free(zu);

  return;
}

/*============================= end of c_disort() =======================*/

/*============================= c_bidir_reflectivity() ==================*/

/*
  Supplies surface bi-directional reflectivity.

  NOTE 1: Bidirectional reflectivity in DISORT is defined by eq. 39 in STWL.
  NOTE 2: Both MU and MU0 (cosines of reflection and incidence angles) are positive.

  Translated from fortran to C by Robert Buras; original name BDREF

  INPUT:

    wvnmlo    : Lower wavenumber (inv cm) of spectral interval
    wvnmhi    : Upper wavenumber (inv cm) of spectral interval
    mu        : Cosine of angle of reflection (positive)
    mup       : Cosine of angle of incidence (positive)
    dphi      : Difference of azimuth angles of incidence and reflection
                (radians)
    brdf_type : BRDF type
    brdf      : BRDF input
    callnum   : number of surface calls

  LOCAL VARIABLES:

    ans       :  Return variable
    badmu     :  minimally allowed value for mu1 and mu2
    flxalb    :  
    irmu      :
    rmu       :
    swvnmlo   : value of wvnmlo from last call of this routine
    swvnmhi   : value of wvnmhi from last call of this routine
    srho0     : value of rho0   from last call of this routine
    sk        : value of k      from last call of this routine
    stheta    : value of theta  from last call of this routine
    ssigma    : value of sigma  from last call of this routine
    st1       : value of t1     from last call of this routine
    st2       : value of t2     from last call of this routine
    sscale    : value of scale  from last call of this routine
    siso      : value of iso    from last call of this routine
    svol      : value of vol    from last call of this routine
    sgeo      : value of geo    from last call of this routine
    su10      : value of u10    from last call of this routine
    spcl      : value of pcl    from last call of this routine
    ssal      : value of sal    from last call of this routine

   Called by- c_dref, c_surface_bidir
   Calls- c_dref, c_bidir_reflectivity_hapke,
          c_bidir_reflectivity_rpv, ocean_brdf, ambrals_brdf
-------------------------------------------------------------------------*/

double c_bidir_reflectivity ( double       wvnmlo,
			      double       wvnmhi,
			      double       mu,
			      double       mup,
			      double       dphi,
			      int          brdf_type,
			      disort_brdf *brdf,
			      int          callnum )
{
  int
    irmu;

  double
    ans, rmu, flxalb;

  static double
    badmu, swvnmlo, swvnmhi, srho0, sk,
    stheta, ssigma, st1, st2, sscale;

#if HAVE_BRDF
    static double
    siso, svol, sgeo;
#endif

  ans = 0.0;

  switch (brdf_type) {
  case BRDF_HAPKE:

    ans = c_bidir_reflectivity_hapke ( wvnmlo, wvnmhi, mu, mup, dphi );

    break;
  case BRDF_RPV:
    if ( swvnmlo != wvnmlo      ||
	 swvnmhi != wvnmhi      ||
	 srho0   != brdf->rpv->rho0   ||
	 sk      != brdf->rpv->k      ||	   	   
	 stheta  != brdf->rpv->theta  ||
	 ssigma  != brdf->rpv->sigma  ||
	 st1     != brdf->rpv->t1     ||
	 st2     != brdf->rpv->t2     ||
	 sscale  != brdf->rpv->scale ) {

      swvnmlo = wvnmlo;
      swvnmhi = wvnmhi;
      srho0   = brdf->rpv->rho0;
      sk      = brdf->rpv->k;
      stheta  = brdf->rpv->theta;
      ssigma  = brdf->rpv->sigma;
      st1     = brdf->rpv->t1;
      st2     = brdf->rpv->t2;
      sscale  = brdf->rpv->scale;

      badmu = 0.0;

      for (irmu=100; irmu>=0; irmu--) {

	rmu = ((double)irmu) * 0.01;

	flxalb = c_dref( wvnmlo, wvnmhi, rmu, brdf_type, brdf, callnum );

	if ( flxalb < 0.0 || flxalb > 1.0 ) {
	  badmu = rmu + 0.01;
	  if (badmu > 1.0)
	    badmu = 1.0;
	  fprintf(stderr,"Using %f as limiting mu in RPV \n",badmu);
	  break;
	}
      }
    }

    ans = c_bidir_reflectivity_rpv ( brdf->rpv, mup, mu, dphi, badmu );

    break;
  case BRDF_CAM:

#if HAVE_BRDF
    /* call C tree saving function */
    /*
     * NOTE: Should group brdf->cam input arguments into the single pointer brdf->cam,
     *       in the same manner as brdf->rpv for c_bidir_reflectivity_rpv().
     */
    ans = ocean_brdf ( wvnmlo, wvnmhi, mu, mup, dphi, 
		       brdf->cam->u10, brdf->cam->pcl, brdf->cam->xsal, callnum);

    /* remove BRDFs smaller than 0 */
    if (ans < 0.0)
      ans = 0.0;

    /* check for NaN */
    if ( ans != ans ) {
      fprintf(stderr,"NaN returned from ocean_brdf: %e %e %e %e %e %e %e %e\n",
	      wvnmlo, wvnmhi, mu, mup, dphi, brdf->cam->u10, brdf->cam->pcl, brdf->cam->xsal);
      ans = 1.0;
    }
#else
    c_errmsg("Error, ocean_brdf is not linked with your code!",DS_ERROR);
#endif
    break;
  case BRDF_AMB:

#if HAVE_BRDF
    /* mu = 0 or dmu = 0 cause problems */
    if ( siso != brdf->ambrals->iso ||
	 svol != brdf->ambrals->vol ||
	 sgeo != brdf->ambrals->geo ) {

      siso = brdf->ambrals->iso;
      svol = brdf->ambrals->vol;
      sgeo = brdf->ambrals->geo;

      badmu = 0.0;

      for (irmu=100; irmu>=0; irmu--) {

	rmu = ((double)irmu) * 0.01;

	flxalb = c_dref( wvnmlo, wvnmhi, rmu, brdf_type, brdf, callnum );

	if ( flxalb < 0.0 || flxalb > 1.0 ) {
	  badmu = rmu + 0.01;
	  if (badmu > 1.0)
	    badmu = 1.0;
	  fprintf(stderr,"Using %f as limiting mu in AMBRALS \n",badmu);
	  break;
	}
      }
    }

    /* convert phi to degrees */
    /*    sdphi = dphi;
	  smup  = mup;
	  smu   = mu; probably no longer needed */

    dphi /= DEG;

    if ( badmu > 0.0 ) {
      if ( mu < badmu )
	mu = badmu;
      if ( mup < badmu )
	mup = badmu;
    }

    /*
     * NOTE: Should group brdf->ambrals input arguments into the single pointer brdf->ambrals,
     *       in the same manner as brdf->rpv for c_bidir_reflectivity_rpv().
     */
    ans = ambrals_brdf (brdf->ambrals->iso, brdf->ambrals->vol, brdf->ambrals->geo, mu, mup, dphi);

    /*    dphi = sdphi;
	  mup  = smup;
	  mu   = smu; probably no longer needed */

    /* check for NaN */
    if ( ans != ans ) {
      fprintf(stderr,"NaN returned from ambrals_brdf: %e %e %e %e %e %e %e %e\n",
	      wvnmlo, wvnmhi, mu, mup, dphi, brdf->ambrals->iso, brdf->ambrals->vol, brdf->ambrals->geo);
      ans = 1.0;
    }
#else
    c_errmsg("Error, ambrals_brdf is not linked with your code!",DS_ERROR);
#endif

    break;
  default:
    fprintf(stderr,"bidir_reflectivity--surface BDRF model %d not known",
	    brdf_type);
    c_errmsg("Exiting...",DS_ERROR);
  }

  return ans;
}

/*============================= end of c_bidir_reflectivity() ===========*/

/*============================= c_bidir_reflectivity_hapke() ============*/

/*    
 * Hapke's BRDF model (times Pi/Mu0):
 *   Hapke, B., Theory of reflectance and emittance spectroscopy, Cambridge University Press, 1993, 
 * eq. 8.89 on page 233. Parameters are from Fig. 8.15 on page 231, except for w.
     
  INPUT:

    wvnmlo : Lower wavenumber (inv cm) of spectral interval
    wvnmhi : Upper wavenumber (inv cm) of spectral interval
    mu     : Cosine of angle of reflection (positive)
    mup    : Cosine of angle of incidence (positive)
    dphi   : Difference of azimuth angles of incidence and reflection
                (radians)

  LOCAL VARIABLES:

    iref   : bidirectional reflectance options; 1 - Hapke's BDR model
    b0     : empirical factor to account for the finite size of particles in Hapke's BDR model
    b      : term that accounts for the opposition effect (retroreflectance, hot spot) in Hapke's BDR model
    ctheta : cosine of phase angle in Hapke's BDR model
    gamma  : albedo factor in Hapke's BDR model
    h0     : H(mu0) in Hapke's BDR model
    h      : H(mu) in Hapke's BDR model
    hh     : angular width parameter of opposition effect in Hapke's BDR model
    p      : scattering phase function in Hapke's BDR model
    theta  : phase angle (radians); the angle between incidence and reflection directions in Hapke's BDR model
    w      : single scattering albedo in Hapke's BDR model

   Called by- c_bidir_reflectivity
-------------------------------------------------------------------------*/

double c_bidir_reflectivity_hapke ( double wvnmlo,
				    double wvnmhi,
				    double mu,
				    double mup,
				    double dphi )
{
  double
    b0,b,ctheta,Xgamm,
    h0,h,hh,p,thetah,w;

  ctheta = mu*mup+sqrt((1.-mu*mu)*(1.-mup*mup))*cos(dphi);
  thetah = acos(ctheta);
  p      = 1.+.5*ctheta;
  hh     =  .06;
  b0     = 1.;
  b      = b0*hh/(hh+tan(.5*thetah));
  w      = 0.6;
  Xgamm  = sqrt(1.-w);
  h0     = (1.+2.*mup)/(1.+2.*Xgamm*mup);
  h      = (1.+2.*mu )/(1.+2.*Xgamm*mu );

  return .25*w*((1.+b)*p+h0*h-1.0)/(mu+mup);
}
  
/*============================= end of c_bidir_reflectivity_hapke() =====*/

/*============================= c_bidir_reflectivity_rpv() ==============*/

/*
  Computes the Rahman, Pinty, Verstraete BRDF.  The incident
  and outgoing cosine zenith angles are MU1 and MU2, respectively,
  and the relative azimuthal angle is PHI.  In this case the incident
  direction is where the radiation is coming from, so MU1>0 and 
  the hot spot is MU2=MU1 and PHI=180 (the azimuth convention is
  different from the original Frank Evans code). 
  The reference is:
  Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
  Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
  With NOAA Advanced Very High Resolution Radiometer Data,
  J. Geophys. Res., 98, 20791-20801.

  Translated from fortran to C by Robert Buras; original name RPV_REFLECTION

  INPUT:

    rho0   :  BRDF rpv: rho0
    k      :  BRDF rpv: k
    theta  :  BRDF rpv: theta
    sigma  :  BRDF rpv snow: sigma
    t1     :  BRDF rpv snow: t1
    t2     :  BRDF rpv snow: t2
    scale  :  BRDF rpv: scale
    mu1    :  Cosine of angle of reflection (positive)
    mu2    :  Cosine of angle of incidence (positive)
    phi    :  Difference of azimuth angles of incidence and reflection
                 (radians)
    badmu  :  minimally allowed value for mu1 and mu2

  LOCAL VARIABLES:

    ans    :  Return value

   Called by- c_bidir_reflectivity
-------------------------------------------------------------------------*/

double c_bidir_reflectivity_rpv ( rpv_brdf_spec *brdf,
                                  double         mu1,
				  double         mu2,
				  double         phi,
				  double         badmu )
{
  double
    m, f, h, cosphi, sin1, sin2, cosg, tan1, tan2, capg,
    hspot, t, g;
  double ans;

  /* This function needs more checking; some constraints are 
     required to avoid albedos larger than 1; in particular,
     the BDREF is limited to 5 times the hotspot value to
     avoid extremely large values at low polar angles */


  /* Azimuth convention different from Frank Evans:
     Here PHI=0 means the backward direction while 
     while in DISORT PHI=0 means forward. */
  phi = M_PI - phi;

  /* Don't allow mu's smaller than BADMU because 
     the albedo is larger than 1 for those */
  if ( badmu > 0.0 ) {
    if ( mu1 < badmu )
      mu1 = badmu;
    if ( mu2 < badmu )
      mu2 = badmu;
  }

  /* Hot spot */
  hspot = brdf->rho0 * ( pow ( 2.0 * mu1 * mu1 * mu1 , brdf->k - 1.0 ) * 
		   ( 1.0 - brdf->theta ) / ( 1.0 + brdf->theta ) / ( 1.0 + brdf->theta )
		   *  ( 2.0 - brdf->rho0 ) 
		   + brdf->sigma / mu1 ) * ( brdf->t1 * exp ( M_PI * brdf->t2 ) + 1.0 );

  /* Hot spot region */
  /* is this bug??? phi <= 1e-4 would be more sensible ... RPB */
  if (phi == 1e-4 && mu1 == mu2)
    return hspot * brdf->scale;
      
  m = pow ( mu1 * mu2 * ( mu1 + mu2 ) , brdf->k - 1.0 );
  cosphi = cos(phi);
  sin1 = sqrt ( 1.0 - mu1 * mu1 );
  sin2 = sqrt ( 1.0 - mu2 * mu2 );
  cosg = mu1 * mu2 + sin1 * sin2 * cosphi;
  g = acos ( cosg );
  f = ( 1.0 - brdf->theta * brdf->theta ) /
    pow ( 1.0 + 2.0 * brdf->theta * cosg + brdf->theta * brdf->theta , 1.5);

  tan1 = sin1 / mu1;
  tan2 = sin2 / mu2;
  capg = sqrt( tan1 * tan1 + tan2 * tan2 - 2.0 * tan1 * tan2 * cosphi );
  h = 1.0 + ( 1.0 - brdf->rho0 ) / ( 1.0 + capg );
  t = 1.0 + brdf->t1 * exp ( brdf->t2 * ( M_PI - g ) );

  ans = brdf->rho0 * ( m * f * h + brdf->sigma / mu1 ) * t * brdf->scale;
      
 if (ans < 0.0)
   ans = 0.0;

 return ans;
}

/*============================= end of c_bidir_reflectivity_rpv() =======*/

/*============================= c_getmom() ==============================*/

/*
 * Calculate phase function Legendre expansion coefficients in various special cases.
 *
 *  INPUT
 *    iphas  Phase function options
 *                       ISOTROPIC: Isotropic
 *                        RAYLEIGH: Rayleigh
 *               HENYEY_GREENSTEIN: Henyey-Greenstein with asymmetry factor GG
 *             HAZE_GARCIA_SIEWART: Haze L as specified by Garcia/Siewert
 *            CLOUD_GARCIA_SIEWART: Cloud C.1 as specified by Garcia/Siewert
 *    gg      Asymmetry factor for Henyey-Greenstein case
 *    nmom    Index of highest Legendre coefficient needed (number of streams 'nstr'
 *            chosen for the discrete ordinate method). Set to -1 for no scattering.
 *  OUTPUT
 *    PMOM(k) Legendre expansion coefficients (k = 0 to nmom)
 *
 * Reference:  Garcia, R. and C. Siewert, 1985: Benchmark Results in Radiative Transfer, 
 *               Transp. Theory and Stat. Physics 14, 437-484, Tables 10 And 17
 */

void c_getmom(int    iphas,
              double  gg,
              int     nmom,
              double *pmom)
{
  const double cldmom[299] = {
    2.544,3.883,4.568,5.235,5.887,6.457,7.177,7.859,8.494,9.286,9.856,10.615,11.229,11.851,12.503,
    13.058,13.626,14.209,14.660,15.231,15.641,16.126,16.539,16.934,17.325,17.673,17.999,18.329,18.588,
    18.885,19.103,19.345,19.537,19.721,19.884,20.024,20.145,20.251,20.330,20.401,20.444,20.477,20.489,
    20.483,20.467,20.427,20.382,20.310,20.236,20.136,20.036,19.909,19.785,19.632,19.486,19.311,19.145,
    18.949,18.764,18.551,18.348,18.119,17.901,17.659,17.428,17.174,16.931,16.668,16.415,16.144,15.883,
    15.606,15.338,15.058,14.784,14.501,14.225,13.941,13.662,13.378,13.098,12.816,12.536,12.257,11.978,
    11.703,11.427,11.156,10.884,10.618,10.350,10.090,9.827,9.574,9.318,9.072,8.822,8.584,8.340,8.110,
    7.874,7.652,7.424,7.211,6.990,6.785,6.573,6.377,6.173,5.986,5.790,5.612,5.424,5.255,5.075,4.915,
    4.744,4.592,4.429,4.285,4.130,3.994,3.847,3.719,3.580,3.459,3.327,3.214,3.090,2.983,2.866,2.766,
    2.656,2.562,2.459,2.372,2.274,2.193,2.102,2.025,1.940,1.869,1.790,1.723,1.649,1.588,1.518,1.461,
    1.397,1.344,1.284,1.235,1.179,1.134,1.082,1.040,0.992,0.954,0.909,0.873,0.832,0.799,0.762,0.731,
    0.696,0.668,0.636,0.610,0.581,0.557,0.530,0.508,0.483,0.463,0.440,0.422,0.401,0.384,0.364,0.349,
    0.331,0.317,0.301,0.288,0.273,0.262,0.248,0.238,0.225,0.215,0.204,0.195,0.185,0.177,0.167,0.160,
    0.151,0.145,0.137,0.131,0.124,0.118,0.112,0.107,0.101,0.097,0.091,0.087,0.082,0.079,0.074,0.071,
    0.067,0.064,0.060,0.057,0.054,0.052,0.049,0.047,0.044,0.042,0.039,0.038,0.035,0.034,0.032,0.030,
    0.029,0.027,0.026,0.024,0.023,0.022,0.021,0.020,0.018,0.018,0.017,0.016,0.015,0.014,0.013,0.013,
    0.012,0.011,0.011,0.010,0.009,0.009,0.008,0.008,0.008,0.007,0.007,0.006,0.006,0.006,0.005,0.005,
    0.005,0.005,0.004,0.004,0.004,0.004,0.003,0.003,0.003,0.003,0.003,0.003,0.002,0.002,0.002,0.002,
    0.002,0.002,0.002,0.002,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
    0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  const double hazelm[82] = {
    2.41260,3.23047,3.37296,3.23150,2.89350,2.49594,2.11361,1.74812,1.44692,1.17714,0.96643,0.78237,
    0.64114,0.51966,0.42563,0.34688,0.28351,0.23317,0.18963,0.15788,0.12739,0.10762,0.08597,0.07381,
    0.05828,0.05089,0.03971,0.03524,0.02720,0.02451,0.01874,0.01711,0.01298,0.01198,0.00904,0.00841,
    0.00634,0.00592,0.00446,0.00418,0.00316,0.00296,0.00225,0.00210,0.00160,0.00150,0.00115,0.00107,
    0.00082,0.00077,0.00059,0.00055,0.00043,0.00040,0.00031,0.00029,0.00023,0.00021,0.00017,0.00015,
    0.00012,0.00011,0.00009,0.00008,0.00006,0.00006,0.00005,0.00004,0.00004,0.00003,0.00003,0.00002,
    0.00002,0.00002,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001};
  register int
    k;

  /* 
   * Screen for invalid inputs
   */
  if (iphas < FIRST_IPHAS || iphas > LAST_IPHAS) {
    c_errmsg("getmom--bad input variable iphas",DS_ERROR);
  }
  if (nmom < 2) {
    c_errmsg("getmom--bad input variable nmom",DS_ERROR);
  }

  pmom[0] = 1.;
  for (k = 1; k < nmom; k++) {
    pmom[k] = 0.;
  }

  switch(iphas) {
    case RAYLEIGH:
      pmom[2] = 0.1;
    break;
    case HENYEY_GREENSTEIN:
      if (gg <= -1. || gg >= 1.) {
        c_errmsg("getmom--bad input variable gg",DS_ERROR);
      }
      for(k = 1; k <= nmom; k++) {
        pmom[k] = pow(gg,(double)k);
      }
    break;
    case HAZE_GARCIA_SIEWERT:
      /* Haze-L phase function */
      for (k = 1; k <= IMIN(82,nmom); k++) {
        pmom[k] = hazelm[k-1]/(double)(2*k+1);
      }
    break;
    case CLOUD_GARCIA_SIEWERT:
      /* Cloud C.1 phase function */
      for (k = 1; k <= IMIN(298,nmom); k++) {
        pmom[k] = cldmom[k-1]/(double)(2*k+1);
      }
    break;
  }

  return;
}

/*============================= end of c_getmom() =======================*/

/*============================= c_asymmetric_matrix() ===================*/

/*
  Solves eigenfunction problem for real asymmetric matrix for which it
  is known a priori that the eigenvalues are real. This is an adaptation
  of a subroutine EIGRF in the IMSL library to use real instead of complex
  arithmetic, accounting for the known fact that the eigenvalues and
  eigenvectors in the discrete ordinate solution are real.
  
  EIGRF is based primarily on EISPACK routines.  The matrix is first
  balanced using the Parlett-Reinsch algorithm.  Then the Martin-Wilkinson
  algorithm is applied. There is a statement 'j = wk(i)' that converts a
  double precision variable to an integer variable; this seems dangerous
  to us in principle, but seems to work fine in practice.
  
  References:

  Dongarra, J. and C. Moler, EISPACK -- A Package for Solving Matrix
      Eigenvalue Problems, in Cowell, ed., 1984: Sources and Development of
      Mathematical Software, Prentice-Hall, Englewood Cliffs, NJ
  Parlett and Reinsch, 1969: Balancing a Matrix for Calculation of
      Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
  Wilkinson, J., 1965: The Algebraic Eigenvalue Problem, Clarendon Press,
      Oxford

   I N P U T    V A R I A B L E S:

       aa    :  input asymmetric matrix, destroyed after solved
        m    :  order of aa
       ia    :  first dimension of aa
    ievec    :  first dimension of evec

   O U T P U T    V A R I A B L E S:

       evec  :  (unnormalized) eigenvectors of aa (column j corresponds to EVAL(J))
       eval  :  (unordered) eigenvalues of aa (dimension m)
       ier   :  if != 0, signals that EVAL(ier) failed to converge;
                   in that case eigenvalues ier+1,ier+2,...,m  are
                   correct but eigenvalues 1,...,ier are set to zero.

   S C R A T C H   V A R I A B L E S:

       wk    :  work area (dimension at least 2*m)
       
   Called by- c_solve_eigen
   Calls- c_errmsg
 -------------------------------------------------------------------*/

void c_asymmetric_matrix(double *aa,
                         double *evec,
                         double *eval,
                         int     m,
                         int     ia,
                         int     ievec,
                         int    *ier,
                         double *wk)
{
  const double
   c1 =    .4375,
   c2 =    .5,
   c3 =    .75,
   c4 =    .95,
   c5 =  16.,
   c6 = 256.;
  int
    noconv,notlas,
    i,ii,in,j,k,ka,kkk,l,lb=0,lll,n,n1,n2;
  double
    col,discri,f,g,h,p=0,q=0,r=0,repl,rnorm,row,
    s,scale,sgn,t,tol,uu,vv,w,x,y,z;

  *ier = 0;
  tol = DBL_EPSILON;
  if (m < 1 || ia < m || ievec < m) {
    c_errmsg("asymmetric_matrix--bad input variable(s)",DS_ERROR);
  }

  /*
   * Handle 1x1 and 2x2 special cases
   */
  if (m == 1) {
    EVAL(1)   = AA(1,1);
    EVEC(1,1) = 1.;
    return;
  }
  else if (m == 2) {
    discri = SQR(AA(1,1)-AA(2,2))+4.*AA(1,2)*AA(2,1);
    if(discri < 0.) {
      c_errmsg("asymmetric_matrix--complex evals in 2x2 case",DS_ERROR);
    }
    sgn = 1.;
    if (AA(1,1) < AA(2,2)) {
     sgn = -1.;
    }
    EVAL(1)   = .5*(AA(1,1)+AA(2,2)+sgn*sqrt(discri));
    EVAL(2)   = .5*(AA(1,1)+AA(2,2)-sgn*sqrt(discri));
    EVEC(1,1) = 1.;
    EVEC(2,2) = 1.;
    if (AA(1,1) == AA(2,2) && (AA(2,1) == 0. || AA(1,2) == 0.)) {
      rnorm     = fabs(AA(1,1))+fabs(AA(1,2))+fabs(AA(2,1))+fabs(AA(2,2));
      w         = tol*rnorm;
      EVEC(2,1) =  AA(2,1)/w;
      EVEC(1,2) = -AA(1,2)/w;
    }
    else {
      EVEC(2,1) = AA(2,1)/(EVAL(1)-AA(2,2));
      EVEC(1,2) = AA(1,2)/(EVAL(2)-AA(1,1));
    }
    return;
  }

  /*
   * Initialize output variables
   */
  *ier = 0;
  memset(eval,0,m*sizeof(double));
  memset(evec,0,ievec*ievec*sizeof(double));
  for (i = 1; i <= m; i++) {
    EVEC(i,i) = 1.;
  }

  /*
   * Balance the input matrix and reduce its norm by diagonal similarity transformation stored in wk;
   * then search for rows isolating an eigenvalue and push them down.
   */
  rnorm = 0.;
  l     = 1;
  k     = m;

S50:

  kkk = k;
  for (j = kkk; j >= 1; j--) {
    row = 0.;
    for (i = 1; i <= k; i++) {
      if (i != j) {
        row += fabs(AA(j,i));
      }
    }
    if (row == 0.) {
      WK(k) = (double)j;
      if (j != k) {
        for (i = 1; i <= k; i++) {
          repl    = AA(i,j);
          AA(i,j) = AA(i,k);
          AA(i,k) = repl;
        }
        for (i = l; i <= m; i++) {
          repl    = AA(j,i);
          AA(j,i) = AA(k,i);
          AA(k,i) = repl;
        }
      }
      k--;
      goto S50;
    }
  }

  /*
   * Search for columns isolating an eigenvalue and push them left.
   */

S100:

  lll = l;
  for (j = lll; j <= k; j++) {
    col = 0.;
    for (i = l; i <= k; i++) {
      if (i != j) {
        col += fabs(AA(i,j));
      }
    }
    if (col == 0.) {
      WK(l) = (double)j;
      if (j != l) {
        for (i = 1; i <= k; i++) {
          repl    = AA(i,j);
          AA(i,j) = AA(i,l);
          AA(i,l) = repl;
        }
        for (i = l; i <= m; i++) {
          repl    = AA(j,i);
          AA(j,i) = AA(l,i);
          AA(l,i) = repl;
        }
      }
      l++;
      goto S100;
    }
  }

  /*
   * Balance the submatrix in rows L through K
   */
  for (i = l; i <= k; i++) {
    WK(i) = 1.;
  }

  noconv = TRUE;
  while (noconv) {
    noconv = FALSE;
    for (i = l; i <= k; i++) {
      col = 0.;
      row = 0.;
      for (j = l; j <= k; j++) {
        if (j != i) {
          col += fabs(AA(j,i));
          row += fabs(AA(i,j));
        }
      }

      f = 1.;
      g = row/c5;
      h = col+row;

      while (col < g) {
        f   *= c5;
        col *= c6;
      }

      g = row*c5;

      while (col >= g) {
        f   /= c5;
        col /= c6;
      }

      /*
       * Now balance
       */
      if ((col+row)/f < c4*h) {
        WK(i)  *= f;
        noconv  = TRUE;
        for (j = l; j <= m; j++) {
          AA(i,j) /= f;
        }
        for (j = 1; j <= k; j++) {
          AA(j,i) *= f;
        }
      }
    }
  }

  if (k-1 >= l+1) {
    /*
     * Transfer A to a Hessenberg form.
     */
    for (n = l+1; n <= k-1; n++) {
      h       = 0.;
      WK(n+m) = 0.;
      scale   = 0.;
      /*
       * Scale column
       */
      for (i = n; i <= k; i++) {
        scale += fabs(AA(i,n-1));
      }
      if (scale != 0.) {
        for (i = k; i >= n; i--) {
          WK(i+m)  = AA(i,n-1)/scale;
          h       += SQR(WK(i+m));
        }
        g        = -F77_SIGN(sqrt(h),WK(n+m));
        h       -= WK(n+m)*g;
        WK(n+m) -= g;
        /*
         * Form (I-(U*UT)/H)*A
         */
        for (j = n; j <= m; j++) {
          f = 0.;
          for (i = k; i >= n; i--) {
            f += WK(i+m)*AA(i,j);
          }
          for (i = n; i <= k; i++) {
            AA(i,j) -= WK(i+m)*f/h;
          }
        }
        /*
         * Form (i-(u*ut)/h)*a*(i-(u*ut)/h)
         */
        for (i = 1; i <= k; i++) {
          f = 0.;
          for (j = k; j >= n; j--) {
            f += WK(j+m)*AA(i,j);
          }
          for (j = n; j <= k; j++) {
            AA(i,j) -= WK(j+m)*f/h;
          }
        }
        WK(n+m)   *= scale;
        AA(n,n-1)  = scale*g;
      }
    }

    for (n = k-2; n >= l; n--) {
      n1 = n+1;
      n2 = n+2;
      f = AA(n+1,n);
      if( f != 0.) {
        f *= WK(n+1+m);
        for (i = n+2; i <= k; i++) {
          WK(i+m) = AA(i,n);
        }
        if (n+1 <= k) {
          for (j = 1; j <= m; j++) {
            g = 0.;
            for (i = n+1; i <= k; i++) {
              g += WK(i+m)*EVEC(i,j);
            }
            g /= f;
            for (i = n+1; i <= k; i++) {
              EVEC(i,j) += g*WK(i+m);
            }
          }
        }
      }
    }
  }

  n = 1;
  for (i = 1; i <= m; i++) {
    for (j = n; j <= m; j++) {
      rnorm += fabs(AA(i,j));
    }
    n = i;
    if (i < l || i > k) {
      EVAL(i) = AA(i,i);
    }
  }

  n = k;
  t = 0.;
  /*
   * Search for next eigenvalues
   */

S400:

  if (n < l) {
    goto S550;
  }

  in = 0;
  n1 = n-1;
  n2 = n-2;

  /*
   * Look for single small sub-diagonal element
   */

S410:

  for (i = l; i <= n; i++) {
    lb = n+l-i;
    if (lb == l) {
      break;
    }
    s = fabs(AA(lb-1,lb-1))+fabs(AA(lb,lb));
    if (s == 0.) {
      s = rnorm;
    }
    if (fabs(AA(lb,lb-1)) <= tol*s) {
      break;
    }
  }

  x = AA(n,n);
  if (lb == n) {
    /*
     * One eigenvalue found
     */
    AA(n,n) = x+t;
    EVAL(n) = AA(n,n);
    n       = n1;
    goto S400;
  }

  y = AA(n1,n1);
  w = AA(n,n1)*AA(n1,n);

  if (lb == n1) {
    /*
     * Two eigenvalues found
     */
    p         = (y-x)*c2;
    q         = p*p+w;
    z         = sqrt(fabs(q));
    AA(n,n)   = x+t;
    x         = AA(n,n);
    AA(n1,n1) = y+t;
    /*
     * Real pair
     */
    z        = p+F77_SIGN(z,p);
    EVAL(n1) = x+z;
    EVAL(n)  = EVAL(n1);

    if (z != 0.) {
      EVAL(n) = x-w/z;
    }
    x = AA(n,n1);
    /*
     * Employ scale factor in case X and Z are very small
     */
    r = sqrt(x*x+z*z);
    p = x/r;
    q = z/r;
    /*
     * Row modification
     */
    for (j = n1; j <= m; j++) {
      z        = AA(n1,j);
      AA(n1,j) =  q*z+p*AA(n,j);
      AA(n, j) = -p*z+q*AA(n,j);
    }
    /*
     * Column modification
     */
    for (i = 1; i <= n; i++) {
      z        = AA(i,n1);
      AA(i,n1) =  q*z+p*AA(i,n);
      AA(i,n ) = -p*z+q*AA(i,n);
    }
    /*
     * Accumulate transformations
     */
    for (i = l; i <= k; i++) {
      z          = EVEC(i,n1);
      EVEC(i,n1) =  q*z+p*EVEC(i,n);
      EVEC(i,n ) = -p*z+q*EVEC(i,n);
    }
    n = n2;
    goto S400;
  }

  if (in == 30) {
    /*
     * No convergence after 30 iterations; set error indicator to
     * the index of the current eigenvalue, and return.
     */
    *ier = n;
    return;
  }

  /*
   * Form shift
   */
  if (in == 10 || in == 20) {
    t += x;
    for (i = l; i <= n; i++) {
      AA(i,i) -= x;
    }
    s = fabs(AA(n,n1))+fabs(AA(n1,n2));
    x = c3*s;
    y = x;
    w = -c1*s*s;
  }

  in++;

  /*
   * Look for two consecutive small sub-diagonal elements
   */
  for (j = lb; j <= n2; j++) {
    i  = n2+lb-j;
    z  = AA(i,i);
    r  = x-z;
    s  = y-z;
    p  = (r*s-w)/AA(i+1,i)+AA(i,i+1);
    q  = AA(i+1,i+1)-z-r-s;
    r  = AA(i+2,i+1);
    s  = fabs(p)+fabs(q)+fabs(r);
    p /= s;
    q /= s;
    r /= s;

    if (i == lb) {
      break;
    }

    uu = fabs(AA(i,i-1))*(fabs(q)+fabs(r));
    vv = fabs(p)*(fabs(AA(i-1,i-1))+fabs(z)+fabs(AA(i+1,i+1)));

    if (uu <= tol*vv) {
      break;
    }
  }

  AA(i+2,i) = 0.;
  for (j = i+3; j <= n; j++) {
    AA(j,j-2) = 0.;
    AA(j,j-3) = 0.;
  }

  /*
   * Double QR step involving rows K to N and columns M to N
   */
  for (ka = i; ka <= n1; ka++) {
    notlas = (ka != n1);
    if (ka == i) {
      s = F77_SIGN(sqrt(p*p+q*q+r*r),p);
      if (lb != i) {
        AA(ka,ka-1) *= -1;
      }
    }
    else {
      p = AA(ka,  ka-1);
      q = AA(ka+1,ka-1);
      r = 0.;
      if (notlas) {
        r = AA(ka+2,ka-1);
      }
      x = fabs(p)+fabs(q)+fabs(r);
      if (x == 0.) {
        continue;
      }
      p /= x;
      q /= x;
      r /= x;
      s  = F77_SIGN(sqrt(p*p+q*q+r*r),p);

      AA(ka,ka-1) = -s*x;
    }

    p += s;
    x  = p/s;
    y  = q/s;
    z  = r/s;
    q /= p;
    r /= p;

    /*
     * Row modification
     */
    for (j = ka; j <= m; j++) {
      p = AA(ka,j)+q*AA(ka+1,j);
      if (notlas) {
        p          += r*AA(ka+2,j);
        AA(ka+2,j) -= p*z;
      }
      AA(ka+1,j) -= p*y;
      AA(ka,  j) -= p*x;
    }

    /*
     * Column modification
     */
    for (ii = 1; ii <= IMIN(n,ka+3); ii++) {
      p = x*AA(ii,ka)+y*AA(ii,ka+1);
      if (notlas) {
        p           += z*AA(ii,ka+2);
        AA(ii,ka+2) -= p*r;
      }
      AA(ii,ka+1) -= p*q;
      AA(ii,ka  ) -= p;
    }

    /*
     * Accumulate transformations
     */
    for (ii = l; ii <= k; ii++) {
      p = x*EVEC(ii,ka)+y*EVEC(ii,ka+1);
      if (notlas) {
        p             += z*EVEC(ii,ka+2);
        EVEC(ii,ka+2) -= p*r;
      }
      EVEC(ii,ka+1) -= p*q;
      EVEC(ii,ka  ) -= p;
    }
  }

  goto S410;

  /*
   * All evals found, now backsubstitute real vector
   */

S550:

  if (rnorm != 0.) {
    for (n = m; n >= 1; n--) {
      n2      = n;
      AA(n,n) = 1.;
      for (i = n-1; i >= 1; i--) {
        w = AA(i,i)-EVAL(n);
        if (w == 0.) {
          w = tol*rnorm;
        }
        r = AA(i,n);
        for (j = n2; j <= n-1; j++) {
          r += AA(i,j)*AA(j,n);
        }
        AA(i,n) = -r/w;
        n2      = i;
      }
    }
    /*
     * End backsubstitution vectors of isolated evals
     */
    for (i = 1; i <= m; i++) {
      if (i < l || i > k) {
        for (j = i; j <= m; j++) {
          EVEC(i,j) = AA(i,j);
        }
      }
    }
    /*
     * Multiply by transformation matrix
     */
    if (k != 0) {
      for (j = m; j >= l; j--) {
        for (i = l; i <= k; i++) {
          z = 0.;
          for (n = l; n <= IMIN(j,k); n++) {
            z += EVEC(i,n)*AA(n,j);
          }
          EVEC(i,j) = z;
        }
      }
    }
  }
  for (i = l; i <= k; i++) {
    for (j = 1; j <= m; j++) {
      EVEC(i,j) *= WK(i);
    }
  }

  /*
   * Interchange rows if permutations occurred
   */
  for (i = l-1; i >= 1; i--) {
    j = WK(i);
    if (i != j) {
      for (n = 1; n <= m; n++) {
        repl      = EVEC(i,n);
        EVEC(i,n) = EVEC(j,n);
        EVEC(j,n) = repl;
      }
    }
  }
  for (i = k+1; i <= m; i++) {
    j = WK(i);
    if (i != j) {
      for (n = 1; n <= m; n++) {
        repl      = EVEC(i,n);
        EVEC(i,n) = EVEC(j,n);
        EVEC(j,n) = repl;
      }
    }
  }

  return;
}

/*============================= end of c_asymmetric_matrix() ============*/

/*============================= c_intensity_components() ================*/

/*
    Calculates the Fourier intensity components at the quadrature
    angles for azimuthal expansion terms (mazim) in eq. SD(2),STWL(6)

    I N P U T    V A R I A B L E S:

       ds      :  Disort state variables
       kk      :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       gc      :  Eigenvectors at polar quadrature angles in eq. SC(1)
       ll      :  Constants of integration in eq. SC(1), obtained by solving scaled version of eq. SC(5);
                  exponential term of eq. SC(12) not included
       lyrcut  :  Logical flag for truncation of computational layer
       mazim   :  Order of azimuthal component
       ncut    :  Number of computational layer where absorption optical depth exceeds ABSCUT
       nn      :  Order of double-Gauss quadrature (NSTR/2)
       taucpr  :  Cumulative optical depth (delta-M-scaled)
       utaupr  :  Optical depths of user output levels in delta-M coordinates;  equal to UTAU if no delta-M
       zz      :  Beam source vectors in eq. SS(19), STWL(24b)
       plk     :  Thermal source vectors z0,z1 by solving eq. SS(16), Y-sub-zero, Y-sub-one in STWL(26ab);
                  plk[].zero, plk[].one (see cdisort.h)

    O U T P U T   V A R I A B L E S:

       uum     :  Fourier components of the intensity in eq. SD(12) (at polar quadrature angles)

    I N T E R N A L   V A R I A B L E S:

       fact    :  exp(-utaupr/umu0)
       zint    :  intensity of m=0 case, in eq. SC(1)

   Called by- c_disort
 -------------------------------------------------------------------*/

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
                            double       *uum)
{
  register int
    iq,jq,lu,lyu;
  register double
    zint;

  /*
   * Loop over user levels
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    lyu = LAYRU(lu);
    if (lyrcut && lyu > ncut) {
      continue;
    }
    for (iq = 1; iq <= ds->nstr; iq++) {
      zint = 0.;
      for (jq = 1; jq <= nn; jq++) {
        zint += GC(iq,jq,lyu)*LL(jq,lyu)*exp(-KK(jq,lyu)*(UTAUPR(lu)-TAUCPR(lyu  )));
      }
      for (jq = nn+1; jq <=ds->nstr; jq++) {
        zint += GC(iq,jq,lyu)*LL(jq,lyu)*exp(-KK(jq,lyu)*(UTAUPR(lu)-TAUCPR(lyu-1)));
      }
      UUM(iq,lu) = zint;
      if (ds->bc.fbeam > 0.) {
        UUM(iq,lu) = zint+ZZ(iq,lyu)*exp(-UTAUPR(lu)/ds->bc.umu0);
      }
      if (ds->flag.planck && mazim == 0) {
        UUM(iq,lu) += ZPLK0(iq,lyu)+ZPLK1(iq,lyu)*UTAUPR(lu);
      }
    }
  }

  return;
}

/*============================= end of c_intensity_components() =========*/

/*============================= c_fluxes() ==============================*/

/*
    Calculates the radiative fluxes, mean intensity, and flux
    derivative with respect to optical depth from the m=0 intensity
    components (the azimuthally-averaged intensity)

    I N P U T    V A R I A B L E S:

       ds       :  Disort state variables
       cmu      :  Abscissae for Gauss quadrature over angle cosine
       cwt      :  Weights for Gauss quadrature over angle cosine
       gc       :  Eigenvectors at polar quadrature angles, SC(1)
       kk       :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       layru    :  Layer number of user level UTAU
       ll       :  Constants of integration in eq. SC(1), obtained by solving scaled version of eq. SC(5);
                   exponential term of eq. SC(12) not included
       lyrcut   :  Logical flag for truncation of comput. layer
       ncut     :  Number of computational layer where absorption optical depth exceeds ABSCUT
       nn       :  Order of double-Gauss quadrature (NSTR/2)
       prntu0   :  TRUE, print azimuthally-averaged intensity at quadrature angles
       taucpr   :  Cumulative optical depth (delta-M-scaled)
       utaupr   :  Optical depths of user output levels in delta-M coordinates;  equal to UTAU if no delta-M
       xr       :  Expansion of thermal source function in eq. SS(14,16), STWL(24c); xr[].zero, xr[].one (see cdisort.h)
       zz       :  Beam source vectors in eq. SS(19), STWL(24b)
       zzg      :  Beam source vectors in eq. KS(10)for a general source constant over a layer
       plk      :  Thermal source vectors z0,z1 by solving eq. SS(16), Y0,Y1 in STWL(26b,a);
                   plk[].zero, plk[].one (see cdisort.h)

    O U T P U T    V A R I A B L E S:

       out      : Disort output variables
       u0c      :  Azimuthally averaged intensities (at polar quadrature angles)

    I N T E R N A L    V A R I A B L E S:

       dirint   :  Direct intensity attenuated
       fdntot   :  Total downward flux (direct + diffuse)
       fl       :  fl[].zero: 'fldir' = direct-beam flux (delta-M scaled), fl[].one 'fldn' = diffuse down-flux (delta-M scaled)
       fnet     :  Net flux (total_down-diffuse_up)
       fact     :  EXP(- UTAUPR/UMU0)
       plsorc   :  Planck source function (thermal)
       zint     :  Intensity of m = 0 case, in eq. SC(1)

   Called by- c_disort
 -------------------------------------------------------------------*/

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
              double        *u0c)
{
  register int
    iq,jq,lu,lyu;
  double
    ang1,ang2,dirint,
    fact=0,fdntot,fnet,plsorc,zint;

  if (ds->flag.prnt[1]) {
    fprintf(stdout,"\n\n                     <----------------------- FLUXES ----------------------->\n"
                   "   Optical  Compu    Downward    Downward    Downward      Upward                    Mean      Planck   d(Net Flux)\n"
                   "     Depth  Layer      Direct     Diffuse       Total     Diffuse         Net   Intensity      Source   / d(Op Dep)\n");
  }

  /*
   * Zero DISORT output arrays
   */
  memset(u0c,0,ds->ntau*ds->nstr*sizeof(double));
  memset(fl,0,ds->ntau*sizeof(disort_pair));

  /*
   * Loop over user levels
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    lyu = LAYRU(lu);

    if (lyrcut && lyu > ncut) {
      /*
       * No radiation reaches this level
       */
      fdntot = 0.;
      fnet   = 0.;
      plsorc = 0.;
      if (ds->flag.prnt[1]) {
        fprintf(stdout,"%10.4f%7d%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%14.3e\n",
                        UTAU(lu),lyu,RFLDIR(lu),RFLDN(lu),fdntot,FLUP(lu),fnet,UAVG(lu),plsorc,DFDT(lu));
      }
      continue;
    }

    if (ds->bc.fbeam > 0.) {
      if ( ds->flag.spher == TRUE ) {
	fact         = exp( - UTAUPR(lu) / CH(lyu) );
	RFLDIR( lu ) = fabs(ds->bc.umu0)*ds->bc.fbeam*
	  exp( - UTAU( lu ) / CH(lyu) );
      }
      else {
	fact       = exp(-UTAUPR(lu)/ds->bc.umu0);
	RFLDIR(lu) = ds->bc.umu0*ds->bc.fbeam*exp(-UTAU(lu)/ds->bc.umu0);
      }
      dirint     = ds->bc.fbeam*fact;
      FLDIR(lu)  = ds->bc.umu0*ds->bc.fbeam*fact;
    }
    else {
      dirint     = 0.;
      FLDIR(lu)  = 0.;
      RFLDIR(lu) = 0.;
    }

    for (iq = 1; iq <= nn; iq++) {
      zint = 0.;
      for (jq = 1; jq <= nn; jq++) {
        zint += GC(iq,jq,lyu)*LL(jq,lyu)*exp(-KK(jq,lyu)*(UTAUPR(lu)-TAUCPR(lyu  )));
      }
      for (jq = nn+1; jq <= ds->nstr; jq++) {
        zint += GC(iq,jq,lyu)*LL(jq,lyu)*exp(-KK(jq,lyu)*(UTAUPR(lu)-TAUCPR(lyu-1)));
      }

      U0C(iq,lu) = zint;
      if (ds->bc.fbeam > 0. ) {
	if ( ds->flag.spher == TRUE ) {
	  U0C(iq,lu) += exp(-ZBEAMA(lyu)*UTAUPR(lu))*
	    ( ZBEAM0(iq,lyu)+ZBEAM1(iq,lyu)*UTAUPR(lu) );
	}
	else {
	  U0C(iq,lu) += ZZ(iq,lyu)*fact;
	}
      }
      if ( ds->flag.general_source == TRUE ) {
	U0C(iq,lu) += ZZG(iq,lyu);
      }
      U0C(iq,lu) += ZPLK0(iq,lyu)+ZPLK1(iq,lyu)*UTAUPR(lu);
      UAVG(lu)   += CWT(nn+1-iq)*U0C(iq,lu);
      UAVGDN(lu) += CWT(nn+1-iq)*U0C(iq,lu);
      FLDN(lu)   += CWT(nn+1-iq)*U0C(iq,lu)*CMU(nn+1-iq);
    }

    for (iq = nn+1; iq <= ds->nstr; iq++) {
      zint = 0.;
      for (jq = 1; jq <= nn; jq++) {
        zint += GC(iq,jq,lyu)*LL(jq,lyu)*exp(-KK(jq,lyu)*(UTAUPR(lu)-TAUCPR(lyu  )));
      }
      for (jq = nn+1; jq <= ds->nstr; jq++) {
        zint += GC(iq,jq,lyu)*LL(jq,lyu)*exp(-KK(jq,lyu)*(UTAUPR(lu)-TAUCPR(lyu-1)));
      }

      U0C(iq,lu) = zint;
      if (ds->bc.fbeam > 0.) {
	if ( ds->flag.spher == TRUE ) {
	  U0C(iq,lu) += exp(-ZBEAMA(lyu)*UTAUPR(lu))*
	    ( ZBEAM0(iq,lyu)+ZBEAM1(iq,lyu)*UTAUPR(lu) );
	}
	else {
	  U0C(iq,lu) += ZZ(iq,lyu)*fact;
	}
      }
      if ( ds->flag.general_source == TRUE ) {
	U0C(iq,lu) += ZZG(iq,lyu);
      }
      U0C(iq,lu) += ZPLK0(iq,lyu)+ZPLK1(iq,lyu)*UTAUPR(lu);
      UAVG(lu)   += CWT(iq-nn)*U0C(iq,lu);
      UAVGUP(lu) += CWT(iq-nn)*U0C(iq,lu);
      FLUP(lu)   += CWT(iq-nn)*U0C(iq,lu)*CMU(iq-nn);
    }
    FLUP(lu)  *= 2.*M_PI;
    FLDN(lu)  *= 2.*M_PI;
    fdntot     = FLDN(lu)+FLDIR(lu);
    fnet       = fdntot-FLUP(lu);
    RFLDN(lu)  = fdntot-RFLDIR(lu);
    UAVG(lu)   = (2.*M_PI*UAVG(lu)+dirint)/(4.*M_PI);
    UAVGSO(lu) =  dirint / (4.*M_PI);
    UAVGDN(lu) = (2.*M_PI*UAVGDN(lu) )/(4.*M_PI);
    UAVGUP(lu) = (2.*M_PI*UAVGUP(lu) )/(4.*M_PI);
    plsorc     = XR0(lyu)+XR1(lyu)*UTAUPR(lu);
    DFDT(lu)   = (1.-SSALB(lyu))*4.*M_PI*(UAVG(lu)-plsorc);

    if (ds->flag.prnt[1]) {
      fprintf(stdout,"%10.4f%7d%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%14.3e\n",
                      UTAU(lu),lyu,RFLDIR(lu),RFLDN(lu),fdntot,FLUP(lu),fnet,UAVG(lu),plsorc,DFDT(lu));
    }
  }

  if (prntu0) {
    fprintf(stdout,"\n\n%s\n"," ******** AZIMUTHALLY AVERAGED INTENSITIES ( at polar quadrature angles ) *******");
    for (lu = 1; lu <= ds->ntau; lu++) {
      fprintf(stdout,"\n%s%10.4f\n\n%s\n",
                     " Optical depth =",UTAU(lu),
                     "     Angle (deg)   cos(Angle)     Intensity     Angle (deg)   cos(Angle)     Intensity");
      for (iq = 1; iq <= nn; iq++) {
        ang1 = acos(CMU(2*nn-iq+1))/DEG;
        ang2 = acos(CMU(     iq  ))/DEG;
        fprintf(stdout,"%16.4f%13.5f%14.3e%16.4f%13.5f%14.3e\n",
                        ang1,CMU(2*nn-iq+1),U0C(iq,   lu),
                        ang2,CMU(     iq  ),U0C(iq+nn,lu));
      }
    }
  }

  return;
}

/*============================= end of c_fluxes() =======================*/

/*============================= c_intensity_correction() ================*/

/*
       Corrects intensity field by using Nakajima-Tanaka algorithm
       (1988). For more details, see Section 3.6 of STWL NASA report.
                I N P U T   V A R I A B L E S

       ds      Disort state variables
       dither  small multiple of machine precision
       flyr    separated fraction in delta-M method
       layru   index of UTAU in multi-layered system
       lyrcut  logical flag for truncation of computational layer
       ncut    total number of computational layers considered
       oprim   delta-M-scaled single-scatter albedo
       phirad  azimuthal angles in radians
       tauc    optical thickness at computational levels
       taucpr  delta-M-scaled optical thickness
       utaupr  delta-M-scaled version of UTAU

                O U T P U T   V A R I A B L E S

       out->UU  corrected intensity field; UU(IU,LU,J)
                 iu=1,ds->numu; lu=1,ds->ntau; j=1,ds->nphi

                I N T E R N A L   V A R I A B L E S

       ctheta  cosine of scattering angle
       dtheta  angle (degrees) to define aureole region as
                    direction of beam source +/- DTHETA
       phasa   actual (exact) phase function
       phasm   delta-M-scaled phase function
       phast   phase function used in TMS correction; actual phase
                    function divided by (1-FLYR*SSALB)
       pl      ordinary Legendre polynomial of degree l, P-sub-l
       plm1    ordinary Legendre polynomial of degree l-1, P-sub-(l-1)
       plm2    ordinary Legendre polynomial of degree l-2, P-sub-(l-2)
       theta0  incident zenith angle (degrees)
       thetap  emergent angle (degrees)
       ussndm  single-scattered intensity computed by using exact
                   phase function and scaled optical depth
                   (first term in STWL(68a))
       ussp    single-scattered intensity from delta-M method
                   (second term in STWL(68a))
       duims   intensity correction term from IMS method
                   (delta-I-sub-IMS in STWL(A.19))

   Called by- c_disort
   Calls- c_single_scat, c_secondary_scat
 -------------------------------------------------------------------*/

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
                            double        *utaupr)
{
  register int
    iu,jp,k,lc,ltau,lu;
  double
    ctheta,dtheta,duims,pl,plm1,plm2,
    theta0=0,thetap=0,ussndm,ussp;

  dtheta = 10.;

  /*
   * Start loop over zenith angles
   */
  for (iu = 1; iu <= ds->numu; iu++) {
    if (UMU(iu) < 0.) {
      /*
       * Calculate zenith angles of incident and emerging directions
       */
      theta0 = acos(-ds->bc.umu0)/DEG;
      thetap = acos(UMU(iu))/DEG;
    }
    /*
     * Start loop over azimuth angles
     */
    for (jp = 1; jp <= ds->nphi; jp++) {
      /*
       * Calculate cosine of scattering angle, eq. STWL(4)
       */
      ctheta = -ds->bc.umu0*UMU(iu)+sqrt((1.-SQR(ds->bc.umu0))*(1.-SQR(UMU(iu))))*cos(PHIRAD(jp));
       /*
        * Initialize phase function                              
        */
      for (lc = 1; lc <= ncut; lc++) {
        PHASA(lc) = 1.;
        PHASM(lc) = 1.;
      }
      /*
       * Initialize Legendre poly. recurrence
       */

      plm1 = 1.;
      plm2 = 0.;
      for (k = 1; k <= ds->nmom; k++) {
        /*
         * Calculate Legendre polynomial of P-sub-l by upward recurrence
         */
        pl   = ((double)(2*k-1)*ctheta*plm1-(double)(k-1)*plm2)/k;
        plm2 = plm1;
        plm1 = pl;

        /*
         * Calculate actual phase function
         */
        for (lc = 1; lc <= ncut; lc++) {
          PHASA(lc) += (double)(2*k+1)*pl*PMOM(k,lc);
        }
        /*
         * Calculate delta-M transformed phase function
         */
        if (k <= ds->nstr-1) {
          for (lc = 1; lc <= ncut; lc++) {
            PHASM(lc) += (double)(2*k+1)*pl*(PMOM(k,lc)-FLYR(lc))/(1.-FLYR(lc));
          }
        }
      }
      /*
       * Apply TMS method, eq. STWL(68)
       */
      for (lc = 1; lc <= ncut; lc++) {
        PHAST(lc) = PHASA(lc)/(1.-FLYR(lc)*SSALB(lc));
      }
      for (lu = 1; lu <= ds->ntau; lu++) {
        if (!lyrcut || LAYRU(lu) < ncut) {
          ussndm        = c_single_scat(dither,LAYRU(lu),ncut,phast,ds->ssalb,taucpr,UMU(iu),ds->bc.umu0,UTAUPR(lu),ds->bc.fbeam);
          ussp          = c_single_scat(dither,LAYRU(lu),ncut,phasm,oprim,    taucpr,UMU(iu),ds->bc.umu0,UTAUPR(lu),ds->bc.fbeam);
          UU(iu,lu,jp) += ussndm-ussp;
        }
      }
      if (UMU(iu) < 0. && fabs(theta0-thetap) <= dtheta) {
        /*
         * Emerging direction is in the aureole (theta0 +/- dtheta).
         * Apply IMS method for correction of secondary scattering below top level.
         */
        ltau = 1;
        if (UTAU(1) <= dither) {
          ltau = 2;
        }
        for (lu = ltau; lu <= ds->ntau; lu++) {
          if(!lyrcut || LAYRU(lu) < ncut) {
            duims         = c_secondary_scat(ds,iu,lu,ctheta,flyr,LAYRU(lu),tauc);
	    UU(iu,lu,jp) -= duims;
          }
        }
      } 
    } /* end loop over azimuth angles */
  } /* end loop over zenith angles */

  return;
}

/*============================= end of c_intensity_correction() =========*/

/*============================= c_secondary_scat() ======================*/

/*
   Calculates secondary scattered intensity of eq. STWL (A7)

                I N P U T   V A R I A B L E S

        ds      Disort state variables
        iu      index of user polar angle
        lu      index of user level
        ctheta  cosine of scattering angle
        flyr    separated fraction f in Delta-M method
        layru   index of utau in multi-layered system
        tauc    cumulative optical depth at computational layers

                I N T E R N A L   V A R I A B L E S

        pspike  2*p"-p"*p", where p" is the residual phase function
        wbar    mean value of single scattering albedo
        fbar    mean value of separated fraction f
        dtau    layer optical depth
        stau    sum of layer optical depths between top of atmopshere and layer layru

   Called by- c_intensity_correction
   Calls- c_xi_func
 -------------------------------------------------------------------*/

double c_secondary_scat(disort_state *ds,
                        int           iu,
                        int           lu,
                        double        ctheta,
                        double       *flyr,
                        int           layru,
                        double       *tauc)
{
  register int
    k,lyr;
  const double
    tiny = 1.e-4;
  double
    dtau,fbar,gbar,pl,plm1,plm2,pspike,
    stau,umu0p,wbar;
  register double
    tmp;

  /*
   * Calculate vertically averaged value of single scattering albedo and separated
   * fraction f, eq. STWL (A.15)
   */
  dtau = UTAU(lu)-TAUC(layru-1);
  wbar = SSALB(layru)*dtau;
  fbar = FLYR(layru)*wbar;
  stau = dtau;
  for (lyr = 1; lyr <= layru-1; lyr++) {
    wbar += DTAUC(lyr)*SSALB(lyr);
    fbar += DTAUC(lyr)*SSALB(lyr)*FLYR(lyr);
    stau += DTAUC(lyr);
  }

  if (wbar <= tiny || fbar <= tiny || stau <= tiny || ds->bc.fbeam <= tiny) {
    return 0.;
  }

  fbar /= wbar;
  wbar /= stau;
  /*
   * Calculate pspike = (2p"-p"*p")
   */
  pspike = 1.;
  gbar   = 1.;
  plm1   = 1.;
  plm2   = 0.;
  /*
   * pspike for l <= 2n-1
   */
  for (k = 1; k <= ds->nstr-1; k++) {
    pl      = ((double)(2*k-1)*ctheta*plm1-(double)(k-1)*plm2)/k;
    plm2    = plm1;
    plm1    = pl;
    pspike += gbar*(2.-gbar)*(double)(2*k+1)*pl;
  }
  /*
   * pspike for l > 2n-1
   */
  for (k = ds->nstr; k <= ds->nmom; k++) {
    pl   = ((double)(2*k-1)*ctheta*plm1-(double)(k-1)*plm2)/k;
    plm2 = plm1;
    plm1 = pl;
    dtau = UTAU(lu)-TAUC(layru-1);
    gbar = PMOM(k,layru)*SSALB(layru)*dtau;
    for (lyr = 1; lyr <= layru-1; lyr++) {
      gbar += PMOM(k,lyr)*SSALB(lyr)*DTAUC(lyr);
    }
    tmp = fbar*wbar*stau;
    if (tmp <= tiny) {
      gbar = 0.;
    }
    else {
      gbar /= tmp;
    }
    pspike += gbar*(2.-gbar)*(double)(2*k+1)*pl;
  }
  umu0p = ds->bc.umu0/(1.-fbar*wbar);
  /*
   * Calculate IMS correction term, eq. STWL (A.13)
   */
  return ds->bc.fbeam/(4.*M_PI)*SQR(fbar*wbar)/(1.-fbar*wbar)*pspike*c_xi_func(-UMU(iu),umu0p,UTAU(lu));
}

/*============================= end of c_secondary_scat() ===============*/

/*============================= c_new_intensity_correction() ============*/

/*
       Corrects intensity field by using alternative Buras-Emde algorithm
       (201X).

                I N P U T   V A R I A B L E S

       ds      Disort state variables
       dither  small multiple of machine precision
       flyr    separated fraction in delta-M method
       layru   index of UTAU in multi-layered system
       lyrcut  logical flag for truncation of computational layer
       ncut    total number of computational layers considered
       oprim   delta-M-scaled single-scatter albedo
       phirad  azimuthal angles in radians
       tauc    optical thickness at computational levels
       taucpr  delta-M-scaled optical thickness
       utaupr  delta-M-scaled version of UTAU

                O U T P U T   V A R I A B L E S

       out->UU  corrected intensity field; UU(IU,LU,J)
                 iu=1,ds->numu; lu=1,ds->ntau; j=1,ds->nphi

                I N T E R N A L   V A R I A B L E S

       ctheta    cosine of scattering angle
       dtheta    angle (degrees) to define aureole region as
                      direction of beam source +/- DTHETA
       phasa     actual (exact) phase function
       phasm     delta-M-scaled phase function
       phast     phase function used in TMS correction; actual phase
                      function divided by (1-FLYR*SSALB)
       pl        ordinary Legendre polynomial of degree l, P-sub-l
       plm1      ordinary Legendre polynomial of degree l-1, P-sub-(l-1)
       plm2      ordinary Legendre polynomial of degree l-2, P-sub-(l-2)
       theta0    incident zenith angle (degrees)
       thetap    emergent angle (degrees)
       ussndm    single-scattered intensity computed by using exact
                     phase function and scaled optical depth
                     (first term in STWL(68a))
       ussp      single-scattered intensity from delta-M method
                     (second term in STWL(68a))
       duims     intensity correction term from IMS method
                     (delta-I-sub-IMS in STWL(A.19))
       nf        number of angular phase integration grid point
                     (zenith angle, theta)
       np        number of angular phase integration grid point
                     (azimuth angle, phi)
       nphase    number of angles for which original phase function
                     (ds->phase) is defined
       mu_eq     cos(theta) phase integration grid points,
                     equidistant in abs(f_phas2)
       norm_phas normalization factor for phase integration
       norm      normalization factor for preparation of phas2
       neg_phas  index whether phas2 is negative
       phas2     residual phase function
       phasr     delta-M scaled phase function
       f_phas2   cumulative integrated phase function phas2
       fbar      mean value of separated fraction f

   Called by- c_disort
   Calls- c_single_scat, c__new_secondary_scat,
          prep_double_scat_integr, c_dbl_vector
 -------------------------------------------------------------------*/

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
				double        *utaupr)
{
  register int
    iu,jp,k,lc,ltau,lu;
  double
    ctheta,dtheta,duims,pl,plm1,plm2,
    theta0=0,thetap=0,ussndm,ussp;

  const int
    nf = 100;
  const double
    tiny = 1e-4;
  int it=0, lyr=0;
  int nphase=ds->nphase;

  double *mu_eq=NULL, *norm_phas=NULL, norm=0.0;
  int *neg_phas=NULL;

  double *phas2=NULL, *phasr=NULL;
  double f_phas2=0.0;
  double fbar=0.0;
  int need_secondary_scattering=0;

  dtheta = 10.;

  /* beginning of BDE stuff */

  /* check whether secondary scattering is performed at all */
  for (iu = 1; iu <= ds->numu; iu++) {
    if (UMU(iu) < 0.) {
      /*
       * Calculate zenith angles of incident and emerging directions
       */
      theta0 = acos(-ds->bc.umu0)/DEG;
      thetap = acos(UMU(iu))/DEG;
      if (fabs(theta0-thetap) <= dtheta) {
	need_secondary_scattering=TRUE;
	break;
      }
    }
  }

  if (need_secondary_scattering==TRUE) {
    /* Initialization of new PSPIKE.                                      */

    mu_eq  = c_dbl_vector(0,nf*ds->ntau-1,"mu_eq");
    norm_phas = c_dbl_vector(0,ds->ntau-1,"norm_phas");
    neg_phas  = c_int_vector(0,nf*ds->ntau-1,"neg_phas");
    phas2 = c_dbl_vector(0,ds->nphase*ds->ntau-1,"phas2");
    phasr = c_dbl_vector(0,ds->nlyr-1,"phasr");

    /* Calculate delta-scaled phase function (phasr) */

    for (it=1; it<=ds->nphase; it++) {

      ctheta = ds->MUP(it);

      for (lc=1; lc<=ds->nlyr; lc++)
	PHASR(lc) = 1.0 - FLYR(lc);

      plm1 = 1.0;
      plm2 = 0.0;

      for (k=1; k<=ds->nstr-1; k++) {

	/* ** Calculate Legendre polynomial of */
	/* ** P-sub-l by upward recurrence     */

	pl = ( (2*k-1) * ctheta * plm1 - (k-1) * plm2 ) / k;
	plm2 = plm1;
	plm1 = pl;

	for (lc=1; lc<=ds->nlyr; lc++)
	  PHASR(lc) += (2*k+1) * pl * ( PMOM(k,lc) - FLYR(lc) );

      }

      /* calculate difference between original and delta-scaled phase
	 functions (phas2) */

      for (lu=1; lu<=ds->ntau; lu++) {

	PHAS2(it,lu) = 0.0;

	/* this could be optimized */
	for (lyr=1; lyr<=LAYRU(lu)-1; lyr++)
	  PHAS2(it,lu) += ( DSPHASE(it,lyr) - PHASR(lyr) ) *
	    SSALB(lyr) * DTAUC(lyr);

	lyr = LAYRU(lu);
	PHAS2(it,lu) += ( DSPHASE(it,lyr) - PHASR(lyr) ) *
	  SSALB(lyr) * ( UTAU(lu) - TAUC(lyr-1) );

      }

    } /* end for it<nphas */

    /* normalize by 1/(ssa*beta*f) */

    for (lu=1; lu<=ds->ntau; lu++) {

      lyr = LAYRU(lu);
      fbar = FLYR(lyr) * SSALB(lyr) * ( UTAU(lu) - TAUC(lyr-1) );

      for (lyr=1; lyr<=LAYRU(lu)-1; lyr++)
	fbar += SSALB(lyr) * DTAUC(lyr) * FLYR(lyr);

      if ( fbar <= tiny || ds->bc.fbeam <= tiny )
	for (it=1; it<=ds->nphase; it++)
	  PHAS2(it,lu) = 0.0;
      else {
	fbar = 1. / fbar;
	for (it=1; it<=ds->nphase; it++)
	  PHAS2(it,lu) *= fbar;
      }

      /* normalize phas2 to 2.0 */

      f_phas2 = 0.0;
      for (it=2; it<=ds->nphase; it++)
	f_phas2 +=
	  ( ds->MUP(it) - ds->MUP(it-1) ) * 0.5 *
	  ( PHAS2(it,lu) + PHAS2(it-1,lu) );

      if (f_phas2 != 0.0) {
	norm = 2.0 / f_phas2;
	for (it=1; it<=ds->nphase; it++)
	  PHAS2(it,lu) *= norm;
      }

    } /* end for lu<ntau */

    prep_double_scat_integr (ds->nphase, ds->ntau, nf, ds->mu_phase,
			     phas2, mu_eq, neg_phas, norm_phas);
  } /* end if (need_secondary_scattering) */

  /* end of BDE stuff */

  /*
   * Start loop over zenith angles
   */
  for (iu = 1; iu <= ds->numu; iu++) {
    if (UMU(iu) < 0.) {
      /*
       * Calculate zenith angles of incident and emerging directions
       */
      theta0 = acos(-ds->bc.umu0)/DEG;
      thetap = acos(UMU(iu))/DEG;
    }
    /*
     * Start loop over azimuth angles
     */
    for (jp = 1; jp <= ds->nphi; jp++) {
      /*
       * Calculate cosine of scattering angle, eq. STWL(4)
       */
      ctheta = -ds->bc.umu0*UMU(iu)+sqrt((1.-SQR(ds->bc.umu0))*(1.-SQR(UMU(iu))))*cos(PHIRAD(jp));
      /*
       * Initialize phase function                              
       */
      for (lc = 1; lc <= ncut; lc++) {
        PHASM(lc) = 1.;
      }

      /* BDE ** Interpolate original phase function */
      /* BDE ** to actual phase function            */

      /* !!! +1: locate starts counting from 0! */
      it = locate ( ds->mu_phase, ds->nphase, ctheta ) + 1;

      for (lc=1; lc<=ncut; lc++)
	PHASA(lc) = DSPHASE(it,lc)
	  + ( ctheta - ds->MUP(it) ) /
	  ( ds->MUP(it+1) - ds->MUP(it) ) *
	  ( DSPHASE(it+1,lc) - DSPHASE(it,lc) );
      /*
       * Initialize Legendre poly. recurrence
       */
      plm1 = 1.;
      plm2 = 0.;
      for (k = 1; k <= ds->nstr-1; k++) {
        /*
         * Calculate Legendre polynomial of P-sub-l by upward recurrence
         */
        pl   = ((double)(2*k-1)*ctheta*plm1-(double)(k-1)*plm2)/k;
        plm2 = plm1;
        plm1 = pl;

        /*
         * Calculate delta-M transformed phase function
         */
	for (lc=1; lc <= ncut; lc++) {
	  PHASM(lc) += (double)(2*k+1)*pl*(PMOM(k,lc)-FLYR(lc))/(1.-FLYR(lc));
	}
      }
      /*
       * Apply TMS method, eq. STWL(68)
       */
      for (lc = 1; lc <= ncut; lc++) {
        PHAST(lc) = PHASA(lc)/(1.-FLYR(lc)*SSALB(lc));
      }
      for (lu = 1; lu <= ds->ntau; lu++) {
        if (!lyrcut || LAYRU(lu) < ncut) {
          ussndm        = c_single_scat(dither,LAYRU(lu),ncut,phast,ds->ssalb,taucpr,UMU(iu),ds->bc.umu0,UTAUPR(lu),ds->bc.fbeam);
          ussp          = c_single_scat(dither,LAYRU(lu),ncut,phasm,oprim,    taucpr,UMU(iu),ds->bc.umu0,UTAUPR(lu),ds->bc.fbeam);
          UU(iu,lu,jp) += ussndm-ussp;
        }
      }
      if (UMU(iu) < 0. && fabs(theta0-thetap) <= dtheta) {
        /*
         * Emerging direction is in the aureole (theta0 +/- dtheta).
         * Apply IMS method for correction of secondary scattering below top level.
         */
        ltau = 1;
        if (UTAU(1) <= dither) {
          ltau = 2;
        }
        for (lu = ltau; lu <= ds->ntau; lu++) {
          if(!lyrcut || LAYRU(lu) < ncut) {
            duims = c_new_secondary_scat(ds,iu,lu,it,ctheta,flyr,
					 LAYRU(lu),tauc,
					 nf,
					 phas2, mu_eq, neg_phas,
					 NORM_PHAS(lu));
	    UU(iu,lu,jp) -= duims;
          }
        }
      } 
    } /* end loop over azimuth angles */
  } /* end loop over zenith angles */

  free(mu_eq); free(norm_phas); free(neg_phas);
  free(phas2); free(phasr);

  return;
}

/*============================= end of c_new_intensity_correction() =====*/

/*============================= prep_double_scat_integr () ==============*/

/*
       Prepares double scattering integration according to alternative
       Buras-Emde algorithm(201X).

                I N P U T   V A R I A B L E S

       nphase    number of angles for which original phase function
                     (ds->phase) is defined
       ntau      
       nf        number of angular phase integration grid point
                     (zenith angle, theta)
       mu_phase  cos(theta) grid of phase function
       phas2     residual phase function

                O U T P U T   V A R I A B L E S

       mu_eq     cos(theta) phase integration grid points,
                     equidistant in abs(f_phas2)
       neg_phas  index whether phas2 is negative
       norm_phas normalization factor for phase integration

                I N T E R N A L   V A R I A B L E S

       f_phas2_abs absolute value of integrated phase function
                      phas2
       f_phas2     cumulative integrated phase function phas2
       df          step length for calculating mu_eq

   Called by- c_new_intensity_correction
   Calls- c_dbl_vector, locate
 -------------------------------------------------------------------*/

void prep_double_scat_integr (int nphase, int ntau,
			      int           nf,
			      double       *mu_phase,
			      double       *phas2,
			      double       *mu_eq,
			      int          *neg_phas,
			      double       *norm_phas)
{
  int it=0, i=0, lu=0;
  double *f_phas2_abs=NULL;
  double f_phas2=0.0, df=0.0;

  f_phas2_abs = c_dbl_vector(0,nphase,"f_phas2_abs");

  for (lu=1; lu<=ntau; lu++) {

    /* calculate integral of |phas2| (f_phas2_abs) */

    F_PHAS2_ABS(1) = 0.0;
    for (it=2; it<=nphase; it++)
      F_PHAS2_ABS(it) = F_PHAS2_ABS(it-1) +
	( MUP(it) - MUP(it-1) ) * 0.5 *
	( fabs( PHAS2(it,lu) ) + fabs ( PHAS2(it-1,lu) ) );

    /* define mu grid which is equidistant in f_phas2_abs (mu_eq);
       find areas of negative phas2 (neg_phas);
       define normalization (norm_phas) */

    f_phas2 = 0.0;
    df = F_PHAS2_ABS(nphase) / (nf-1);
    MU_EQ(1,lu) = -1.0;

    if ( PHAS2(1,lu) > 0.0 )
      NEG_PHAS(1,lu) = FALSE;
    else
      NEG_PHAS(1,lu) = TRUE;

    it = 1;
    for (i=2; i<=nf-1; i++) {
      f_phas2 += df;

      while ( F_PHAS2_ABS(it+1) < f_phas2 )
	it++;

      MU_EQ(i,lu) = MUP(it)
	+ ( f_phas2 - F_PHAS2_ABS(it) ) /
	( F_PHAS2_ABS(it+1) - F_PHAS2_ABS(it) ) *
	( MUP(it+1) - MUP(it) );

      if ( PHAS2(it,lu) > 0.0 && PHAS2(it+1,lu) > 0.0 )
	NEG_PHAS(i,lu) = FALSE;
      else {
	if ( PHAS2(it,lu) < 0.0 && PHAS2(it+1,lu) < 0.0 )
	  NEG_PHAS(i,lu) = TRUE;
	else {
	  if ( PHAS2(it,lu) + ( f_phas2 - F_PHAS2_ABS(it) ) /
	       ( F_PHAS2_ABS(it+1) - F_PHAS2_ABS(it) ) *
	       ( PHAS2(it+1,lu) - PHAS2(it,lu) ) > 0.0 )
	    NEG_PHAS(i,lu) = FALSE;
	  else
	    NEG_PHAS(i,lu) = TRUE;
	}
      }

    } /* end for i<nf */

    MU_EQ(nf,lu) = 1.0;
    if ( PHAS2(nphase,lu) > 0.0 )
      NEG_PHAS(nf,lu) = FALSE;
    else
      NEG_PHAS(nf,lu) = TRUE;

    NORM_PHAS(lu) = F_PHAS2_ABS(nphase) / ( (nf-1) * M_PI );

  } /* end for lu<ntau */

  free(f_phas2_abs);
}

/*============================= end of prep_double_scat_integr() ========*/

/*============================= c_new_secondary_scat() ==================*/

/*
   Calculates secondary scattered intensity, new method (see BDE)

                I N P U T   V A R I A B L E S

        ds        Disort state variables
        iu        index of user polar angle
        lu        index of user level
	it	  index where ctheta contained in mu grid of exact
	             phase function
        ctheta    cosine of scattering angle
        flyr      separated fraction f in Delta-M method
        layru     index of utau in multi-layered system
        tauc      cumulative optical depth at computational layers
        nf        number of angular phase integration grid point
                     (zenith angle, theta)
        phas2     residual phase function
        mu_eq     cos(theta) phase integration grid points,
                     equidistant in abs(f_phas2)
        neg_phas  index whether phas2 is negative
        norm_phas normalization factor for phase integration

                I N T E R N A L   V A R I A B L E S

        pspike  2*p"-p"*p", where p" is the residual phase function
        pspike1 2*p", where p" is the residual phase function
        pspike2 p"*p", where p" is the residual phase function
        wbar    mean value of single scattering albedo
        fbar    mean value of separated fraction f
        dtau    layer optical depth
        stau    sum of layer optical depths between top of atmopshere and layer layru
	umu0p
        nphase  number of angles for which original phase function
                   (ds->phase) is defined

   Called by- c_new_intensity_correction
   Calls- calc_phase_squared, c_xi_func
 -------------------------------------------------------------------*/

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
			    double        norm_phas)
{
  register int
    lyr;
  const double
    tiny = 1.e-4;
  double
    dtau,fbar,pspike,
    stau,umu0p,wbar;
  int nphase=ds->nphase;

  double pspike1=0.0, pspike2=0.0;

  /*
   * Calculate vertically averaged value of single scattering albedo and separated
   * fraction f, eq. STWL (A.15)
   */
  dtau = UTAU(lu)-TAUC(layru-1);
  wbar = SSALB(layru)*dtau;
  fbar = FLYR(layru)*wbar;
  stau = dtau;
  for (lyr = 1; lyr <= layru-1; lyr++) {
    wbar += DTAUC(lyr)*SSALB(lyr);
    fbar += DTAUC(lyr)*SSALB(lyr)*FLYR(lyr);
    stau += DTAUC(lyr);
  }

  if (wbar <= tiny || fbar <= tiny || stau <= tiny || ds->bc.fbeam <= tiny) {
    return 0.;
  }

  fbar /= wbar;
  wbar /= stau;

  /* Calculate pspike1=P" */

  pspike1 = PHAS2(it,lu) + ( ctheta - ds->MUP(it) ) /
    ( ds->MUP(it+1) - ds->MUP(it) ) * ( PHAS2(it+1,lu) - PHAS2(it,lu) );

  pspike2 = calc_phase_squared (ds->nphase, lu, ctheta, nf,
				ds->mu_phase, phas2, mu_eq, neg_phas,
				norm_phas);

  pspike = 2.*pspike1 - pspike2;

  umu0p = ds->bc.umu0/(1.-fbar*wbar);

  /*
   * Calculate IMS correction term, eq. STWL (A.13)
   */
  return ds->bc.fbeam/(4.*M_PI)*SQR(fbar*wbar)/(1.-fbar*wbar)*pspike*c_xi_func(-UMU(iu),umu0p,UTAU(lu));
}

/*============================= end of c_new_secondary_scat() ===========*/

/*============================= calc_phase_squared() ====================*/

/*
   Calculates squared phase function (see BDE)

                I N P U T   V A R I A B L E S

        nphase  number of angles for which original phase function
                   (ds->phase) is defined
        lu        index of user level
        ctheta  cosine of scattering angle
        nf        number of angular phase integration grid point
                     (zenith angle, theta)
        mu_phase  cos(theta) grid of phase function
        phas2     residual phase function
        mu_eq     cos(theta) phase integration grid points,
                     equidistant in abs(f_phas2)
        neg_phas  index whether phas2 is negative
        norm_phas normalization factor for phase integration

                I N T E R N A L   V A R I A B L E S

        pspike2  p"*p", where p" is the residual phase function; return value
	mu1arr
	stheta   corresponding sin of ctheta
	smueq    corresponding sin of mu_eq
	phint    phase function integrated over phi
	scr

   Called by- c_new_secondary_scat
 -------------------------------------------------------------------*/

double calc_phase_squared (int           nphase,
			   int           lu,
			   double        ctheta,
			   int           nf,
			   double       *mu_phase,
			   double       *phas2,
			   double       *mu_eq,
			   int          *neg_phas,
			   double        norm_phas)
{
  int j=0, k=0, it=0;

  double pspike2=0.0, stheta=0.0;
  double smueq=0.0, phint=0.0;

  double mumin=0.0, mumax=0.0;
  int imin=0, imax=0;
  double D=0.0, C=0.0, Dp=0.0, Cp=0.0;
  int cutting=FALSE;

  stheta = sqrt( 1.0 - ctheta * ctheta );

  /* calculate pspike2 */

  /* Note: MU_EQ(j.lu) is mu_1; ctheta is mu; MUP(k) is mu_i in BDE(201X) */


  for (j=1;j<=nf;j++) {

    /* special case: second scattering angle does not depend on
       azimuth of first scattering angle */
    if (ctheta==1.0 || MU_EQ(j,lu)==1.0) {
      it = locate ( mu_phase, nphase, MU_EQ(j,lu)*ctheta ) + 1;
      phint = M_PI * ( PHAS2(it,lu)
		       + ( MU_EQ(j,lu)*ctheta - MUP(it) )
		       / ( MUP (it+1) - MUP(it) )
		       * ( PHAS2(it+1,lu) - PHAS2(it,lu) ) );
      if (ctheta==1.0)
	phint /= 2.0;
    }
    else {
      phint = 0.0;

      smueq = sqrt ( 1. - MU_EQ(j,lu)*MU_EQ(j,lu) );

      /* locate integration borders */
      mumin = ctheta *  MU_EQ(j,lu) - stheta * smueq;
      mumax = ctheta *  MU_EQ(j,lu) + stheta * smueq;

      /* cut where mu_1 = mu_2 */
      if (MU_EQ(j,lu) < mumax) {
	mumax = MU_EQ(j,lu);
	cutting=TRUE;
      }
      else
	cutting=FALSE;

      if (mumin<mumax) {
	imin = locate ( mu_phase, nphase, mumin)+1;
	imax = locate ( mu_phase, nphase, mumax)+1;

	k=imin;
	/* assuming SPF is linear in mu */
	D = ( PHAS2(k+1,lu) - PHAS2(k,lu) ) / ( MUP(k+1) - MUP(k) );
	C = PHAS2(k,lu) - MUP(k) * D;

	phint +=  ( D * ctheta * MU_EQ(j,lu) + C ) * M_PI / 2.0;

	for (k=imin+1;k<=imax;k++) {

	  Dp = ( PHAS2(k+1,lu) - PHAS2(k,lu) ) / ( MUP(k+1) - MUP(k) );
	  Cp = PHAS2(k,lu) - MUP(k) * Dp;

	  phint += 
	    ( Dp - D ) * sqrt ( 1.0 - ctheta * ctheta
				- MU_EQ(j,lu) * MU_EQ(j,lu)
				+ 2.0 * ctheta * MU_EQ(j,lu) * MUP(k)
				- MUP(k) * MUP(k) )
	    + ( ( Dp - D )* ctheta * MU_EQ(j,lu) + Cp - C ) *
	    asin ( ( ctheta * MU_EQ(j,lu) - MUP(k) )
		   / ( smueq * stheta ) );

	  D=Dp;
	  C=Cp;
	}

	if (cutting==TRUE)
	  phint += - D * sqrt ( 1.0 - ctheta * ctheta
			      + 2.0 * MU_EQ(j,lu) * MU_EQ(j,lu) *
			      ( ctheta - 1.0 ) )
	    - ( D * ctheta * MU_EQ(j,lu) + C ) *
	    asin ( ( ctheta - 1.0 ) * MU_EQ(j,lu)
		   / ( smueq * stheta ) );
	else
	  phint += ( D * ctheta * MU_EQ(j,lu) + C ) * M_PI / 2.0;
      }
    }

    if (j==1 || j==nf) {
      if ( NEG_PHAS(j,lu) == TRUE )
	pspike2 = pspike2 - 0.5 * phint;
      else
	pspike2 = pspike2 + 0.5 * phint;
    }
    else {
      if ( NEG_PHAS(j,lu) == TRUE )
	pspike2 = pspike2 - phint;
      else
	pspike2 = pspike2 + phint;
    }

  }

  pspike2 *= norm_phas;

  return pspike2;
}

/*============================= end of calc_phase_squared() =============*/

/*============================= c_disort_set() ==========================*/

/*
    Perform miscellaneous setting-up operations

    I N P U T  V A R I A B L E S

       ds         Disort state variables
       deltam
       tauc

    O U T P U T     V A R I A B L E S:

       If ds->flag.usrtau is FALSE
       ds->ntau
       ds->utau 

       If ds->flag.usrang is FALSE
       ds->numu
       ds->umu
  
       cmu,cwt     computational polar angles and corresponding quadrature weights
       dtaucpr
       expbea      transmission of direct beam
       flyr        separated fraction in delta-m method
       gl          phase function legendre coefficients multiplied by (2l+1) and single-scatter albedo
       layru       computational layer in which utau falls
       lyrcut      flag as to whether radiation will be zeroed below layer ncut
       ncut        computational layer where absorption optical depth first exceeds  abscut
       nn          ds->nstr/2
       oprim       delta-m-scaled single-scatter albedo
       taucpr      delta-m-scaled optical depth
       utaupr      delta-m-scaled version of  utau

   Called by- c_disort
   Calls- c_gaussian_quadrature, c_errmsg

 ---------------------------------------------------------------------*/

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
                  double       *utaupr)
{
  register int
    iq,iu,k,lc,lu;
  const double
    abscut = 10.;
  double
    abstau,chtau_tmp,f,taup,zenang;

  if (!ds->flag.usrtau) {
   /*
    * Set output levels at computational layer boundaries
    */
    for (lc = 0;  lc <= ds->ntau-1; lc++) {
      UTAU(lc+1) = TAUC(lc);
    }
  }

  /*
   * Apply delta-M scaling and move description of computational layers to local variables
   */
  TAUCPR(0) = 0.;
  abstau    = 0.;
  for (lc = 1; lc <= ds->nlyr; lc++) {
    PMOM(0,lc)  = 1.;
    if (abstau < abscut) {
      *ncut = lc;
    }
    abstau += (1.-SSALB(lc))*DTAUC(lc);
    if (!deltam) {
      OPRIM(lc)   = SSALB(lc);
      DTAUCPR(lc) = DTAUC(lc);
      TAUCPR(lc)  = TAUC(lc);
      for (k = 0; k <= ds->nstr-1; k++) {
        GL(k,lc)  = (double)(2*k+1)*OPRIM(lc)*PMOM(k,lc);
      }
      f = 0.;
    }
    else {
      /*
       * Do delta-M transformation
       */
      f           = PMOM(ds->nstr,lc);
      OPRIM(lc)   = SSALB(lc)*(1.-f)/(1.-f*SSALB(lc));
      DTAUCPR(lc) = (1.-f*SSALB(lc))*DTAUC(lc);
      TAUCPR(lc)  = TAUCPR(lc-1)+DTAUCPR(lc);
      for (k = 0; k <= ds->nstr-1; k++) {
        GL(k,lc)  = (double)(2*k+1)*OPRIM(lc)*(PMOM(k,lc)-f)/(1.-f);
      }
    }

    FLYR(lc)   = f;
  }

  /*
   * Calculate Chapman function if spherical geometry, set expbea and
   * ch for beam source.
   */
  if( (ds->flag.ibcnd == GENERAL_BC && ds->bc.fbeam > 0.) ||
      (ds->flag.ibcnd == GENERAL_BC && ds->flag.general_source )) {

    CHTAU(0)  = 0.;
    EXPBEA(0) = 1.;
    zenang    = acos(ds->bc.umu0)/DEG;
    
    if( ds->flag.spher == TRUE && ds->bc.umu0 < 0. ) {
      EXPBEA(0) = exp(-c_chapman(1,0.,tauc,ds->nlyr,ds->zd,
				 ds->dtauc,zenang,ds->radius));
    }
    if ( ds->flag.spher == TRUE ) {
      for (lc = 1; lc <= *ncut; lc++) {
        taup        = TAUCPR(lc-1) + DTAUCPR(lc)/2.;
	/* Need Chapman function at top (0.0) and middle (0.5) of layer */
        CHTAU(lc  ) = c_chapman(lc, 0.,   taucpr,ds->nlyr,ds->zd,
				dtaucpr,zenang,ds->radius);
        chtau_tmp   = c_chapman(lc, 0.5,  taucpr,ds->nlyr,ds->zd,
				dtaucpr,zenang,ds->radius);
        CH(lc)      = taup/chtau_tmp;
        EXPBEA(lc)  = exp(-CHTAU(lc));
      }
    }
    else {
      for (lc = 1; lc <= *ncut; lc++) {
        CH(lc)     = ds->bc.umu0;
        EXPBEA(lc) = exp(-TAUCPR(lc)/ds->bc.umu0);
      }
    }
  }
  else {
    for (lc = 1; lc <= *ncut; lc++) {
      EXPBEA(lc) = 0.;
    }
  }

  /*
   * If no thermal emission, cut off medium below absorption optical depth = abscut ( note that
   * delta-M transformation leaves absorption optical depth invariant ).  Not worth the
   * trouble for one-layer problems, though.
   */
  *lyrcut = FALSE;
  if (abstau >= abscut && !ds->flag.planck && ds->flag.ibcnd != SPECIAL_BC && ds->nlyr > 1) {
    *lyrcut = TRUE;
  }
  if(!*lyrcut) *ncut = ds->nlyr;

  /* 
   * Set arrays defining location of user output levels within delta-M-scaled computational mesh
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    for (lc = 1; lc < ds->nlyr; lc++) {
      if (UTAU(lu) >= TAUC(lc-1) && UTAU(lu) <= TAUC(lc)) {
        break;
      }
    }

    UTAUPR(lu) = UTAU(lu);
    if (deltam) {
      UTAUPR(lu) = TAUCPR(lc-1)+(1.-SSALB(lc)*FLYR(lc))*(UTAU(lu)-TAUC(lc-1));
    }
    LAYRU(lu) = lc;
  }

  /*
   * Calculate computational polar angle cosines and associated quadrature weights for Gaussian
   * quadrature on the interval (0,1) (upward)
   */
  *nn = ds->nstr/2;
  c_gaussian_quadrature(*nn,cmu,cwt);

  /*
   * Downward (neg) angles and weights
   */
  for (iq = 1; iq <= *nn; iq++) {
    CMU(iq+*nn) = -CMU(iq);
    CWT(iq+*nn) =  CWT(iq);
  }

  if (ds->flag.ibcnd == GENERAL_BC && ds->bc.fbeam > 0.) {
    /*
     * Compare beam angle to comput. angles
     */
    for (iq = 1; iq <= *nn; iq++) {      
      if (fabs(ds->bc.umu0-CMU(iq))/fabs(ds->bc.umu0) < 1.e-4) {
        c_errmsg("cdisort_set--beam angle=computational angle; change ds.nstr",DS_ERROR);
      }
    }
  }

  if (!ds->flag.usrang || ds->flag.onlyfl) {
    /*
     * Set output polar angles to computational polar angles
     */
    for (iu = 1; iu <= *nn; iu++) {
      UMU(iu) = -CMU(*nn+1-iu);
    }
    for (iu = *nn+1; iu <=ds->nstr; iu++) {
      UMU(iu) =  CMU(iu-*nn);
    }
  }

  if (ds->flag.usrang && ds->flag.ibcnd == SPECIAL_BC) {
    /*
     * Shift positive user angle cosines to upper locations and put negatives in lower locations
     */
    for (iu = 1; iu <= ds->numu/2; iu++) {
      UMU(iu+ds->numu/2) = UMU(iu);
    }
    for (iu = 1; iu <= ds->numu/2; iu++) {
      UMU(iu) = -UMU((ds->numu/2)+1-iu);
    }
  }

  return;
}

/*============================= end of c_disort_set() ===================*/

/*============================= c_set_matrix() ==========================*/

/*
    Calculate coefficient matrix for the set of equations obtained from the
    boundary conditions and the continuity-of-intensity-at-layer-interface equations.

    Store in the special banded-matrix format required by LINPACK routines


    I N P U T      V A R I A B L E S:

       ds       :  Disort state variables
       bdr      :  surface bidirectional reflectivity
       cmu,cwt  :  abscissae, weights for Gauss quadrature over angle cosine
       delm0    :  Kronecker delta, delta-sub-m0
       gc       :  Eigenvectors at polar quadrature angles, SC(1)
       kk       :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       lyrcut   :  Logical flag for truncation of computational layers
       nn       :  Number of streams in a hemisphere (NSTR/2)
       ncut     :  Total number of computational layers considered
       taucpr   :  Cumulative optical depth (delta-M-scaled)

   O U T P U T     V A R I A B L E S:

       cband    :  Left-hand side matrix of linear system eq. SC(5), scaled by eq. SC(12); 
                   in banded form required by LINPACK solution routines
       ncol     :  Number of columns in cband


   I N T E R N A L    V A R I A B L E S:

       irow     :  Points to row in CBAND
       jcol     :  Points to position in layer block
       lda      :  Row dimension of CBAND
       ncd      :  Number of diagonals below or above main diagonal
       nshift   :  For positioning number of rows in band storage
       wk       :  Temporary storage for EXP evaluations


   BAND STORAGE

      LINPACK requires band matrices to be input in a special
      form where the elements of each diagonal are moved up or
      down (in their column) so that each diagonal becomes a row.
      (The column locations of diagonal elements are unchanged.)

      Example:  if the original matrix is

          11 12 13  0  0  0
          21 22 23 24  0  0
           0 32 33 34 35  0
           0  0 43 44 45 46
           0  0  0 54 55 56
           0  0  0  0 65 66

      then its LINPACK input form would be:

           *  *  *  +  +  +  , * = not used
           *  * 13 24 35 46  , + = used for pivoting
           * 12 23 34 45 56
          11 22 33 44 55 66
          21 32 43 54 65  *

      If A is a band matrix, the following program segment
      will convert it to the form (ABD) required by LINPACK
      band-matrix routines:

        n  = (column dimension of a, abd)
        ml = (band width below the diagonal)
        mu = (band width above the diagonal)
        m = ml+mu+1;
        for (j = 1; j <= n; j++) {
          i1 = IMAX(1,j-mu);
          i2 = IMIN(n,j+ml);
          for (i = i1; i <= i2; i++) {
            k = i-j+m;
            ABD(k,j) = A(i,j);
          }
        }

      This uses rows  ml+1 through  2*ml+mu+1  of ABD.
      The total number of rows needed in ABD is 2*ml+mu+1.
      In the example above, n = 6, ml = 1, mu = 2, and the
      row dimension of ABD must be >= 5.

   Called by- c_disort, c_albtrans
 -------------------------------------------------------------------*/

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
                  double       *wk)
{
  int
    mi     = ds->nstr/2,
    mi9m2  = 9*mi-2,
    nnlyri = ds->nstr*ds->nlyr,
    nn     = ds->nstr/2;
  register int
    iq,irow,jcol,jq,k,lc,lda,ncd,nncol,nshift;
  double
    expa,sum;

  memset(cband,0,mi9m2*nnlyri*sizeof(double));

  ncd    = 3*nn-1;
  lda    = 3*ncd+1;
  nshift = lda-2*ds->nstr+1;
  *ncol  = 0;

  /*
   * Use continuity conditions of eq. STWJ(17) to form coefficient matrix in STWJ(20);
   * employ scaling transformation STWJ(22)
   */
  for (lc = 1; lc <= ncut; lc++) {
    for (iq = 1; iq <= nn; iq++) {
      WK(iq) = exp(KK(iq,lc)*DTAUCPR(lc));
    }
    jcol = 0;
    for (iq = 1; iq <= nn; iq++) {
      *ncol += 1;
      irow   = nshift-jcol;
      for (jq = 1; jq <= ds->nstr; jq++) {
        CBAND(irow+ds->nstr,*ncol) =  GC(jq,iq,lc);
        CBAND(irow,         *ncol) = -GC(jq,iq,lc)*WK(iq);
        irow++;
      }
      jcol++;
    }

    for (iq = nn+1; iq <= ds->nstr; iq++) {
      *ncol += 1;
      irow = nshift-jcol;
      for (jq = 1; jq <= ds->nstr; jq++) {
        CBAND(irow+ds->nstr,*ncol) =  GC(jq,iq,lc)*WK(ds->nstr+1-iq);
        CBAND(irow,         *ncol) = -GC(jq,iq,lc);
        irow++;
      }
      jcol++;
    }
  }

  /*
   * Use top boundary condition of STWJ(20a) for first layer
   */
  jcol = 0;
  for (iq = 1; iq <= nn; iq++) {
    expa = exp(KK(iq,1)*TAUCPR(1));
    irow = nshift-jcol+nn;
    for (jq = nn; jq >= 1; jq--) {
      CBAND(irow,jcol+1) = GC(jq,iq,1)*expa;
      irow++;
    }
    jcol++;
  }

  for (iq = nn+1; iq <=ds->nstr; iq++) {
    irow = nshift-jcol+nn;
    for (jq = nn; jq >= 1; jq--) {
      CBAND(irow,jcol+1) = GC(jq,iq,1);
      irow++;
    }
    jcol++;
  }

  /*
   * Use bottom boundary condition of STWJ(20c) for last layer
   */
  nncol = *ncol-ds->nstr;
  jcol  = 0;
  for (iq = 1; iq <= nn; iq++) {
    nncol++;
    irow = nshift-jcol+ds->nstr;
    for (jq = nn+1; jq <= ds->nstr; jq++) {
      if (lyrcut || ( ds->flag.lamber && delm0 == 0. ) ) {
        /*
         * No azimuthal-dependent intensity if Lambert surface; 
         * no intensity component if truncated bottom layer
         */
        CBAND(irow,nncol) = GC(jq,iq,ncut);
      }
      else {
        sum = 0.;
        for (k = 1; k <= nn; k++) {
          sum += CWT(k)*CMU(k)*BDR(jq-nn,k)*GC(nn+1-k,iq,ncut);
        }
        CBAND(irow,nncol) = GC(jq,iq,ncut)-(1.+delm0)*sum;
      }
      irow++;
    }
    jcol++;
  }

  for (iq = nn+1; iq <= ds->nstr; iq++) {
    nncol++;
    irow = nshift-jcol+ds->nstr;
    expa = WK(ds->nstr+1-iq);
    for (jq = nn+1; jq <= ds->nstr; jq++) {
      if (lyrcut || (ds->flag.lamber && delm0 == 0.)) {
        CBAND(irow,nncol) = GC(jq,iq,ncut)*expa;
      }
      else {
        sum = 0.;
        for (k = 1; k <= nn; k++) {
          sum += CWT(k)*CMU(k)*BDR(jq-nn,k)*GC(nn+1-k,iq,ncut);
        }
        CBAND(irow,nncol) = (GC(jq,iq,ncut)-(1.+delm0)*sum)*expa;
      }
      irow++;
    }
    jcol++;
  }

  return;
}

/*============================= end of c_set_matrix() ===================*/

/*============================= c_single_scat() =========================*/

/*
        Calculates single-scattered intensity from eqs. STWL (65b,d,e)

                I N P U T   V A R I A B L E S
        
        dither   small multiple of machine precision
        layru    index of utau in multi-layered system
        nlyr     number of sublayers
        phase    phase functions of sublayers
        omega    single scattering albedos of sublayers
        tau      optical thicknesses of sublayers
        umu      cosine of emergent angle
        umu0     cosine of incident zenith angle
        utau     user defined optical depth for output intensity
        fbeam   incident beam radiation at top


   Called by- c_intensity_correction
 -------------------------------------------------------------------*/

double c_single_scat(double   dither,
                     int      layru,
                     int      nlyr,
                     double  *phase,
                     double  *omega,
                     double  *tau,
                     double   umu,
                     double   umu0,
                     double   utau,
                     double   fbeam)
{
  register int
    lyr;
  double
    ans,exp0,exp1;

  ans  = 0.;
  exp0 = exp(-utau/umu0);

  if (fabs(umu+umu0) <= dither) {
    /*
     * Calculate downward intensity when umu=umu0, eq. STWL (65e)
     */
    for (lyr = 1; lyr <= layru-1; lyr++) {
      ans += OMEGA(lyr)*PHASE(lyr)*(TAU(lyr)-TAU(lyr-1));
    }
    ans = fbeam/(4.*M_PI*umu0)*exp0*(ans+OMEGA(layru)*PHASE(layru)*(utau-TAU(layru-1)));
    return ans;
  }

  if (umu > 0.) {
    /*
     * Upward intensity, eq. STWL (65b)
     */
    for (lyr = layru; lyr <= nlyr; lyr++) {
      exp1  = exp(-((TAU(lyr)-utau)/umu+TAU(lyr)/umu0));
      ans  += OMEGA(lyr)*PHASE(lyr)*(exp0-exp1);
      exp0  = exp1;
    }
  }
  else {
    /*
     * Downward intensity, eq. STWL (65d)
     */
    for (lyr = layru; lyr >= 1; lyr--) {
      exp1  = exp(-((TAU(lyr-1)-utau)/umu+TAU(lyr-1)/umu0));
      ans  += OMEGA(lyr)*PHASE(lyr)*(exp0-exp1);
      exp0  = exp1;
    }
  }
  ans *= fbeam/(4.*M_PI*(1.+umu/umu0));

  return ans;
}

/*============================= end of c_single_scat() ==================*/

/*============================= c_solve_eigen() =========================*/

/*
   Solves eigenvalue/vector problem necessary to construct homogeneous
   part of discrete ordinate solution; STWJ(8b), STWL(23f)
   ** NOTE ** Eigenvalue problem is degenerate when single scattering
              albedo = 1;  present way of doing it seems numerically more
              stable than alternative methods that we tried

   I N P U T     V A R I A B L E S:

       ds     :  Disort state variables
       lc     :
       gl     :  Delta-M scaled Legendre coefficients of phase function
                 (including factors 2l+1 and single-scatter albedo)
       cmu    :  Computational polar angle cosines
       cwt    :  Weights for quadrature over polar angle cosine
       mazim  :  Order of azimuthal component
       nn     :  Half the total number of streams
       ylmc   :  Normalized associated Legendre polynomial
                 at the quadrature angles CMU


   O U T P U T    V A R I A B L E S:

       cc     :  C-sub-ij in eq. SS(5); needed in SS(15&18)
       eval   :  NN eigenvalues of eq. SS(12), STWL(23f) on return
                 from asymmetric_matrix but then square roots taken
       evecc  :  NN eigenvectors  (G+) - (G-)  on return
                 from asymmetric_matrix ( column j corresponds to EVAL(j) )
                 but then  (G+) + (G-)  is calculated from SS(10),
                 G+  and  G-  are separated, and  G+  is stacked on
                 top of  G-  to form NSTR eigenvectors of SS(7)
       gc     :  Permanent storage for all NSTR eigenvectors, but
                 in an order corresponding to KK
       kk     :  Permanent storage for all NSTR eigenvalues of SS(7),
                 but re-ordered with negative values first ( square
                 roots of EVAL taken and negatives added )


   I N T E R N A L   V A R I A B L E S:

       ab            :  Matrices AMB (alpha-beta), APB (alpha+beta) in reduced eigenvalue problem (see cdisort.h)
       array         :  Complete coefficient matrix of reduced eigenvalue
                        problem: (alpha+beta)*(alpha-beta)
       gpplgm        :  (g+) + (g-) (cf. eqs. SS(10-11))
       gpmigm        :  (g+) - (g-) (cf. eqs. SS(10-11))
       wk            :  Scratch array required by asymmetric_matrix

   Called by- c_disort, c_albtrans
   Calls- c_asymmetric_matrix, c_errmsg
 -------------------------------------------------------------------*/

/*
 * NOTE: Here the scratch array ARRAY(,) is half the size in each dimension compared to other subroutines
 */
#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*(ds->nstr/2)]

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
                   double       *wk)
{
  int
    ier;
  register int
    iq,jq,kq,l;
  double
    alpha,beta,gpmigm,gpplgm,sum;

  /*
   * Calculate quantities in eqs. SS(5-6), STWL(8b,15,23f)
   */
  for (iq = 1; iq <= nn; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      sum = 0.;
      for (l = mazim; l <= ds->nstr-1; l++) {
        sum += GL(l,lc)*YLMC(l,iq)*YLMC(l,jq);
      }
      CC(iq,jq) = .5*sum*CWT(jq);
    }
    for (jq = 1; jq <= nn; jq++) {
      /*
       * Fill remainder of array using symmetry relations  C(-mui,muj) = C(mui,-muj) and C(-mui,-muj) = C(mui,muj)
       */
      CC(iq+nn,jq   ) = CC(iq,jq+nn);
      CC(iq+nn,jq+nn) = CC(iq,jq   );
      /*
       * Get factors of coeff. matrix of reduced eigenvalue problem
       */
      alpha      = CC(iq,jq   )/CMU(iq);
      beta       = CC(iq,jq+nn)/CMU(iq);
      AMB(iq,jq) = alpha-beta;
      APB(iq,jq) = alpha+beta;
    }
    AMB(iq,iq) -= 1./CMU(iq);
    APB(iq,iq) -= 1./CMU(iq);
  }
  /*
   * Finish calculation of coefficient matrix of reduced eigenvalue problem: 
   * get matrix product (alpha+beta)*(alpha-beta); SS(12),STWL(23f)
   */
  for (iq = 1; iq <= nn; iq++) {
    for (jq = 1; jq <= nn; jq++) {
      sum = 0.;
      for (kq = 1; kq <= nn; kq++) {
        sum += APB(iq,kq)*AMB(kq,jq);
      }
      ARRAY(iq,jq) = sum;
    }
  }

  /*
   * Find (real) eigenvalues and eigenvectors
   */
  c_asymmetric_matrix(array,evecc,eval,nn,ds->nstr/2,ds->nstr,&ier,wk);

  if (ier > 0) {
    fprintf(stderr,"\n\n asymmetric_matrix--eigenvalue no. %4d didn't converge.  Lower-numbered eigenvalues wrong.\n",ier);
    c_errmsg("asymmetric_matrix--convergence problems",DS_ERROR);
  }

  for (iq = 1; iq <= nn; iq++) {
    EVAL(iq)     = sqrt(fabs(EVAL(iq)));
    KK(iq+nn,lc) = EVAL(iq);
    /*
     * Add negative eigenvalue
     */
    KK(nn+1-iq,lc) = -EVAL(iq);
  }

  /*
   * Find eigenvectors (G+) + (G-) from SS(10) and store temporarily in APB array
   */
  for (jq = 1; jq <= nn; jq++) {
    for (iq = 1; iq <= nn; iq++) {
      sum = 0.;
      for (kq = 1; kq <= nn; kq++) {
        sum += AMB(iq,kq)*EVECC(kq,jq);
      }
      APB(iq,jq) = sum/EVAL(jq);
    }
  }
  for (jq = 1; jq <= nn; jq++) {
    for (iq = 1; iq <= nn; iq++) {
      gpplgm = APB(  iq,jq);
      gpmigm = EVECC(iq,jq);
      /*
       * Recover eigenvectors G+,G- from their sum and difference; stack them to get eigenvectors of full system
       * SS(7) (JQ = eigenvector number)
       */
      EVECC(iq,   jq) = .5*(gpplgm+gpmigm);
      EVECC(iq+nn,jq) = .5*(gpplgm-gpmigm);
      /*
       * Eigenvectors corresponding to negative eigenvalues (corresp. to reversing sign of 'k' in SS(10) )
       */
      gpplgm *= -1; 
      EVECC(iq,   jq+nn)     = .5*(gpplgm+gpmigm);
      EVECC(iq+nn,jq+nn)     = .5*(gpplgm-gpmigm);
      GC(nn+iq,  nn+jq,  lc) = EVECC(iq,   jq   );
      GC(nn-iq+1,nn+jq,  lc) = EVECC(iq+nn,jq   );
      GC(nn+iq,  nn-jq+1,lc) = EVECC(iq,   jq+nn);
      GC(nn-iq+1,nn-jq+1,lc) = EVECC(iq+nn,jq+nn);
    }
  }

  return;
}

/*============================= end of c_solve_eigen() ==================*/

/*============================= c_solve0() ==============================*/

/*
        Construct right-hand side vector B for general boundary
        conditions STWJ(17) and solve system of equations obtained
        from the boundary conditions and the continuity-of-
        intensity-at-layer-interface equations.
        Thermal emission contributes only in azimuthal independence.

    I N P U T      V A R I A B L E S:

       ds       :  Disort input variables
       bdr      :  Surface bidirectional reflectivity
       bem      :  Surface bidirectional emissivity
       bplanck  :  Bottom boundary thermal emission
       cband    :  Left-hand side matrix of linear system eq. SC(5),
                   scaled by eq. SC(12); in banded form required
                   by LINPACK solution routines
       cmu,cwt  :  Abscissae, weights for Gauss quadrature
                   over angle cosine
       expbea   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
       lyrcut   :  Logical flag for truncation of computational layers
       mazim    :  Order of azimuthal component
       ncol     :  Number of columns in CBAND
       nn       :  Order of double-Gauss quadrature (NSTR/2)
       ncut     :  Total number of computational layers considered
       tplanck  :  Top boundary thermal emission
       taucpr   :  Cumulative optical depth (delta-M-scaled)
       zz       :  Beam source vectors in eq. SS(19), STWL(24b)
       zzg      :  Beam source vectors in eq. KS(10)for a general source constant over a layer
       plk      :  Thermal source vectors z0,z1 by solving eq. SS(16), Y0,Y1 in STWL(26b,a);
                   plk[].zero, plk[].one (see cdisort.h)

    O U T P U T     V A R I A B L E S:

       b        :  Right-hand side vector of eq. SC(5) going into
                   sgbsl; returns as solution vector of eq. SC(12),
                   constants of integration without exponential term
      ll        :  Permanent storage for B, but re-ordered

   I N T E R N A L    V A R I A B L E S:

       ipvt     :  Integer vector of pivot indices
       it       :  Pointer for position in  B
       ncd      :  Number of diagonals below or above main diagonal
       rcond    :  Indicator of singularity for cband
       z        :  Scratch array required by sgbco

   Called by- c_disort
   Calls- c_sgbco, c_errmsg, c_sgbsl
 +-------------------------------------------------------------------*/

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
              disort_pair  *plk)
{
  register int
    ipnt,iq,it,jq,lc,ncd;
  double
    rcond,sum,diff;
  
  memset(b,0,ds->nstr*ds->nlyr*sizeof(double));
  
  /*
   * Construct B, STWJ(20a,c) for parallel beam+bottom
   * reflection+thermal emission at top and/or bottom
   */
  if (mazim > 0 && ( ds->bc.fbeam > 0.  || ds->flag.general_source) ) {
    /*
     * Azimuth-dependent case (never called if FBEAM = 0)
     */
    if ( lyrcut == TRUE || ds->flag.lamber == TRUE ) {
      /*
       * No azimuthal-dependent intensity for Lambert surface; no
       * intensity component for truncated bottom layer
       */
      for (iq = 1; iq <= nn; iq++) {
        /*
         * Top boundary
         */
	if ( ds->flag.spher == TRUE ) {
	  B(iq) = - ZBEAM0(nn+1-iq,1);
	}
	else {
	  B(iq) = - ZZ(nn+1-iq,1);
	}
	if ( ds->flag.general_source == TRUE ) {
	  B(iq) -= ZZG(nn+1-iq,1);
	  //aky	  B(iq) = B(iq) - ZZG(nn+1-iq,1);
	}
        /*
         * Bottom boundary
         */
	if ( ds->flag.spher == TRUE ) {
	  B(ncol-nn+iq) = - exp(-ZBEAMA(ncut)*TAUCPR(ncut))*
	    (ZBEAM0(iq+nn,ncut) + ZBEAM1(iq+nn,ncut)*TAUCPR(ncut));
	}
	else {
	  B(ncol-nn+iq) = - ZZ(iq+nn,ncut)*EXPBEA(ncut);
	}
	if ( ds->flag.general_source == TRUE ) {
	  B(ncol-nn+iq) -=  ZZG(iq+nn,ncut);
	  //aky	  B(ncol-nn+iq) = B(ncol-nn+iq)  - ZZG(iq+nn,ncut);
	}
      }
    }
    else {
      for (iq = 1; iq <= nn; iq++) {
	if ( ds->flag.spher == TRUE ) {
	  B(iq) = - ZBEAM0(nn+1-iq,1);
	}
	else {
	  B(iq) = - ZZ(nn+1-iq,1);
	}
	if ( ds->flag.general_source == TRUE ) {
	  B(iq) -= ZZG(nn+1-iq,1);
	  //aky	  B(iq) = B(iq) - ZZG(nn+1-iq,1);
	}
	if ( ds->flag.spher == TRUE ) {
	  c_errmsg("solve0--BDR not implemented for pseudo-spherical geometry",
		   DS_WARNING);
	}
	else {
	  sum   = 0.;
	  for (jq = 1; jq <= nn; jq++) {
	    sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*ZZ(nn+1-jq,ncut)*EXPBEA(ncut);
	  }
	  B(ncol-nn+iq) = sum;
	  if ( ds->flag.general_source == TRUE ) {
	    sum   = 0.;
	    for (jq = 1; jq <= nn; jq++) {
	      sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*ZZG(nn+1-jq,ncut);
	    }
	    B(ncol-nn+iq) += sum;
	  }
	}
        if (ds->bc.fbeam > 0.) {
	  if ( ds->flag.spher == TRUE ) {
	    c_errmsg("solve0--BDR not implemented for pseudo-spherical geometry",
		     DS_WARNING)  ;
	  }
	  else {
	    B(ncol-nn+iq) += (BDR(iq,0)*ds->bc.umu0*ds->bc.fbeam/
			      M_PI-ZZ(iq+nn,ncut))*EXPBEA(ncut);
	  }
        }
	if ( ds->flag.general_source == TRUE ) {
	    B(ncol-nn+iq) += -ZZG(iq+nn,ncut);
	}
      }
    }
    /*
     * Continuity condition for layer interfaces of eq. STWJ(20b)
     */
    it = nn;
    diff = 0;
    for (lc = 1; lc <= ncut-1; lc++) {
      for (iq = 1; iq <= ds->nstr; iq++) {
	if ( ds->flag.general_source == TRUE ) {
	  diff = (ZZG(iq,lc+1)-ZZG(iq,lc));
	}
	if ( ds->flag.spher == TRUE ) {
	  B(++it) = exp(-ZBEAMA(lc+1)*TAUCPR(lc))*
	    (ZBEAM0(iq,lc+1)+ZBEAM1(iq,lc+1)*TAUCPR(lc))
	    -  exp(-ZBEAMA(lc)*TAUCPR(lc))*
	    (ZBEAM0(iq,lc)+ZBEAM1(iq,lc)*TAUCPR(lc))
	    + diff;
	}
	else {
	  B(++it) = (ZZ(iq,lc+1)-ZZ(iq,lc))*EXPBEA(lc)  + diff;
	}
      }
    }
  }
  else {
    /*
     * Azimuth-independent case
     */
    if (ds->bc.fbeam == 0. && ds->flag.general_source == FALSE ) {
      for (iq = 1; iq <= nn; iq++) {
        /*
         * Top boundary
         */
        B(iq) = -ZPLK0(nn+1-iq,1)+ds->bc.fisot+tplanck;
      }
      if ( lyrcut == TRUE ) {
        /*
         * No intensity component for truncated bottom layer
         */
        for (iq = 1; iq <= nn; iq++) {
          /*
           * Bottom boundary
           */
          B(ncol-nn+iq) = -ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut);
        }
      }
      else {
        for (iq = 1; iq <= nn; iq++) {
          sum = 0.;
          for (jq = 1; jq <= nn; jq++) {
            sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*
	      (ZPLK0(nn+1-jq,ncut)+ZPLK1(nn+1-jq,ncut)*TAUCPR(ncut));
          }
          B(ncol-nn+iq) = 2.*sum+BEM(iq)*bplanck-
	    ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut);
        }
      }
      /*
       * Continuity condition for layer interfaces, STWJ(20b)
       */
      it = nn;
      for (lc = 1; lc <= ncut-1; lc++) {
        for (iq = 1; iq <= ds->nstr; iq++) {
          B(++it) = ZPLK0(iq,lc+1)-ZPLK0(iq,lc)+
	    (ZPLK1(iq,lc+1)-ZPLK1(iq,lc))*TAUCPR(lc);
        }
      }
    }
    else {
      if ( ds->flag.spher == TRUE ) {
	for (iq = 1; iq <= nn; iq++) 
	  B(iq) = -ZBEAM0(nn+1-iq,1)-ZPLK0(nn+1-iq,1)+ds->bc.fisot+tplanck;
      }
      else {
	for (iq = 1; iq <= nn; iq++) 
	  B(iq) = -ZZ(nn+1-iq,1)-ZPLK0(nn+1-iq,1)+ds->bc.fisot+tplanck;
      }
      if ( ds->flag.general_source == TRUE ) {
	for (iq = 1; iq <= nn; iq++) 
	  B(iq) -= ZZG(nn+1-iq,1);
	//aky	  B(iq) = B(iq) - ZZG(nn+1-iq,1);
      }
      if (lyrcut) {
	if ( ds->flag.spher == TRUE ) {
	  for (iq = 1; iq <= nn; iq++) {
	    B(ncol-nn+iq) = -exp(-ZBEAMA(ncut)*TAUCPR(ncut))*
	      (ZBEAM0(iq+nn,ncut)+ ZBEAM1(iq+nn,ncut)*TAUCPR(ncut))
	      -ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut);
	  }
	}
	else {
	  for (iq = 1; iq <= nn; iq++) {
	    B(ncol-nn+iq) = -ZZ(iq+nn,ncut)*EXPBEA(ncut)
	      -ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut);
	  }
	}
	if ( ds->flag.general_source == TRUE ) {
	  for (iq = 1; iq <= nn; iq++) 
	    B(ncol-nn+iq) -= ZZG(iq+nn,ncut);
	  //aky	    B(ncol-nn+iq) = B(ncol-nn+iq) - ZZG(iq+nn,ncut);
	}
      }
      else {
	if ( ds->flag.spher == TRUE ) {
	  for (iq = 1; iq <= nn; iq++) {
	    sum = 0.;
	    for (jq = 1; jq <= nn; jq++) {
	      sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*
		( exp(-ZBEAMA(ncut)*TAUCPR(ncut))*
		  (ZBEAM0(nn+1-jq,ncut)+ZBEAM1(nn+1-jq,ncut)*TAUCPR(ncut))
		  + ZZG(nn+1-jq,ncut)
		  + ZPLK0(nn+1-jq,ncut)+ZPLK1(nn+1-jq,ncut)*TAUCPR(ncut));
	    }
	    B(ncol-nn+iq) = 2.0*sum +
	      ( BDR(iq,0)*ds->bc.umu0*ds->bc.fbeam/M_PI) *EXPBEA(ncut)
	      -  exp(-ZBEAMA(ncut)*TAUCPR(ncut))*
	      (ZBEAM0(iq+nn,ncut)+ZBEAM1(iq+nn,ncut)*TAUCPR(ncut))
	      - ZZG(iq+nn,ncut)
	      + BEM(iq)*bplanck
	      -ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut)
	      +ds->bc.fluor;
	  }
	}
	else {
	  for (iq = 1; iq <= nn; iq++) {
	    sum = 0.;
	    for (jq = 1; jq <= nn; jq++) {
	      sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*
		(ZZ(nn+1-jq,ncut)*EXPBEA(ncut)+ZPLK0(nn+1-jq,ncut)
		 + ZZG(nn+1-jq,ncut)
		 +ZPLK1(nn+1-jq,ncut)*TAUCPR(ncut));
	    }
	    B(ncol-nn+iq) = 2.*sum+
	      (BDR(iq,0)*ds->bc.umu0*ds->bc.fbeam/M_PI-ZZ(iq+nn,ncut))
	      *EXPBEA(ncut)
	      - ZZG(iq+nn,ncut)
	      +BEM(iq)*bplanck-ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut)
	      +ds->bc.fluor;
	  }
	}
      }
      it = nn;
      if ( ds->flag.spher == TRUE ) {
	for (lc = 1; lc <= ncut-1; lc++) {
	  for (iq = 1; iq <= ds->nstr; iq++) {
	    B(++it) = exp(-ZBEAMA(lc+1)*TAUCPR(lc))*
	      (ZBEAM0(iq,lc+1)+ZBEAM1(iq,lc+1)*TAUCPR(lc))
	      -exp(-ZBEAMA(lc)*TAUCPR(lc))*
	      (ZBEAM0(iq,lc)+ZBEAM1(iq,lc)*TAUCPR(lc))
	      +ZZG(iq,lc+1)-ZZG(iq,lc)
	      +ZPLK0(iq,lc+1)-ZPLK0(iq,lc)+
	      (ZPLK1(iq,lc+1)-ZPLK1(iq,lc))*TAUCPR(lc);
	  }
	}
      }
      else {
	for (lc = 1; lc <= ncut-1; lc++) {
	  for (iq = 1; iq <= ds->nstr; iq++) {
	    B(++it) = (ZZ(iq,lc+1)-ZZ(iq,lc))*EXPBEA(lc)
	      +ZZG(iq,lc+1)-ZZG(iq,lc)
	      +ZPLK0(iq,lc+1)-ZPLK0(iq,lc)
	      +(ZPLK1(iq,lc+1)-ZPLK1(iq,lc))*TAUCPR(lc);
	  }
	}
      }
    }
  }

  /*
   * Find L-U (lower/upper triangular) decomposition of band matrix
   * CBAND and test if it is nearly singular (note: CBAND is
   * destroyed) (CBAND is in LINPACK packed format)
   */
  rcond = 0.;
  ncd   = 3*nn-1;
  c_sgbco(cband,(9*(ds->nstr/2)-2),ncol,ncd,ncd,ipvt,&rcond,z);

  if (1.+rcond == 1.) {
    c_errmsg("solve0--sgbco says matrix near singular",DS_WARNING);
  }

  /*
   * Solve linear system with coeff matrix CBAND and R.H. side(s) B
   * after CBAND has been L-U decomposed. Solution is returned in B.
   */

  c_sgbsl(cband,(9*(ds->nstr/2)-2),ncol,ncd,ncd,ipvt,b,0);

  /*
   * Zero CBAND (it may contain 'foreign' elements upon returning from
   * LINPACK); necessary to prevent errors
   */
  memset(cband,0,(9*(ds->nstr/2)-2)*(ds->nstr*ds->nlyr)*sizeof(double));

  for (lc = 1; lc <= ncut; lc++) {
    ipnt = lc*ds->nstr-nn;
    for (iq = 1; iq <= nn; iq++) {
      LL(nn-iq+1,lc) = B(ipnt-iq+1);
      LL(nn+iq,  lc) = B(ipnt+iq  );
    }
  }

  return;
}

/*============================= c_surface_bidir() =======================*/

/*
       Computes user's' surface bidirectional properties, STWL(41)

   I N P U T     V A R I A B L E S:

       ds     :  Disort input variables
       cmu    :  Computational polar angle cosines (Gaussian)
       delm0  :  Kronecker delta, delta-sub-m0
       mazim  :  Order of azimuthal component
       nn     :  Order of Double-Gauss quadrature (ds->nstr/2)
       callnum:  number of surface calls

    O U T P U T     V A R I A B L E S:

       bdr :  Fourier expansion coefficient of surface bidirectional
                 reflectivity (computational angles)
       rmu :  Surface bidirectional reflectivity (user angles)
       bem :  Surface directional emissivity (computational angles)
       emu :  Surface directional emissivity (user angles)

    I N T E R N A L     V A R I A B L E S:

       dref   :  Directional reflectivity
       gmu    :  The NMUG angle cosine quadrature points on (0,1)
                 NMUG is set in cdisort.h
       gwt    :  The NMUG angle cosine quadrature weights on (0,1)

   Called by- c_disort
   Calls- c_gaussian_quadrature, c_bidir_reflectivity
+---------------------------------------------------------------------*/

void c_surface_bidir(disort_state *ds,
                     double        delm0,
                     double       *cmu,
                     int           mazim,
                     int           nn,
                     double       *bdr,
                     double       *emu,
                     double       *bem,
                     double       *rmu,
		     int           callnum)
{
  static int
    pass1 = TRUE;
  register int
    iq,iu,jg,jq,k;
  double
    dref,sum;
  static double
    gmu[NMUG],gwt[NMUG];
  
  if (pass1) {
    pass1 = FALSE;
    c_gaussian_quadrature(NMUG/2,gmu,gwt);
    for (k = 1; k <= NMUG/2; k++) {
      GMU(k+NMUG/2) = -GMU(k);
      GWT(k+NMUG/2) =  GWT(k);
    }
  }

  memset(bdr,0,(ds->nstr/2)*((ds->nstr/2)+1)*sizeof(double));
  memset(bem,0,(ds->nstr/2)*sizeof(double));

  /*
   * Compute Fourier expansion coefficient of surface bidirectional reflectance
   * at computational angles eq. STWL (41)
   */
  if (ds->flag.lamber && mazim == 0) {
    for (iq = 1; iq <= nn; iq++) {
      BEM(iq) = 1.-ds->bc.albedo;
      for (jq = 0; jq <= nn; jq++) {
        BDR(iq,jq) = ds->bc.albedo;
      }
    }
  }
  else if (!ds->flag.lamber) {
    for (iq = 1; iq <= nn; iq++) {
      for (jq = 1; jq <= nn; jq++) {
        sum = 0.;
        for (k = 1; k <= NMUG; k++) {
          sum += GWT(k) * 
	    c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, CMU(iq), CMU(jq),
				   M_PI * GMU(k), ds->flag.brdf_type, &ds->brdf, callnum)
	    * cos((double)mazim * M_PI * GMU(k) );
        }
        BDR(iq,jq) = .5*(2.-delm0)*sum;
      }
      if (ds->bc.fbeam > 0.) {
        sum = 0.;
        for(k = 1; k <= NMUG; k++) {
          sum += GWT (k) *
	    c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, CMU(iq), ds->bc.umu0,
				   M_PI * GMU(k), ds->flag.brdf_type, &ds->brdf, callnum )
	    * cos((double)mazim * M_PI * GMU(k) );
        }
        BDR(iq,0) = .5*(2.-delm0)*sum;
      }
    }
    if (mazim == 0) {
      /*
       * Integrate bidirectional reflectivity at reflection polar angle cosines -CMU- and incident angle
       * cosines -GMU- to get directional emissivity at computational angle cosines -CMU-.
       */
      for (iq = 1; iq <= nn; iq++) {
        dref = 0.;
        for (jg = 1; jg <= NMUG; jg++) {
          sum = 0.;
          for (k = 1; k <= NMUG/2; k++) {
            sum += GWT(k) * GMU(k) *
	      c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, CMU(iq), GMU(k),
				     M_PI * GMU(jg), ds->flag.brdf_type, &ds->brdf, callnum );
          }
          dref += GWT(jg)*sum;
        }
        BEM(iq) = 1.-dref;
      }
    }
  }
  /*
   * Compute Fourier expansion coefficient of surface bidirectional reflectance at user angles eq. STWL (41)
   */
  if(!ds->flag.onlyfl && ds->flag.usrang) {
    memset(emu,0,ds->numu*sizeof(double));
    memset(rmu,0,ds->numu*((ds->nstr/2)+1)*sizeof(double));
    for (iu = 1; iu <= ds->numu; iu++) {
      if (UMU(iu) > 0.) {
        if(ds->flag.lamber && mazim == 0) {
          for (iq = 0; iq <= nn; iq++) {
            RMU(iu,iq) = ds->bc.albedo;
          }
          EMU(iu) = 1.-ds->bc.albedo;
        }
        else if (!ds->flag.lamber) {
          for (iq = 1; iq <= nn; iq++) {
            sum = 0.;
            for (k = 1; k <= NMUG; k++) {
              sum += GWT(k) *
		c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, UMU(iu), CMU(iq),
				       M_PI * GMU(k), ds->flag.brdf_type, &ds->brdf, callnum )
		* cos( (double)mazim * M_PI * GMU(k) );
            }
            RMU(iu,iq) = .5*(2.-delm0)*sum;
          }
          if (ds->bc.fbeam > 0.) {
            sum = 0.;
            for (k = 1; k <= NMUG; k++) {
              sum += GWT(k) *
		c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, UMU(iu),
				       ds->bc.umu0, M_PI * GMU(k),
				       ds->flag.brdf_type, &ds->brdf, callnum )
		* cos( (double)mazim * M_PI * GMU(k) );
            }
            RMU(iu,0) = .5*(2.-delm0)*sum;
          }
          if (mazim == 0) {
            /*
             * Integrate bidirectional reflectivity at reflection angle cosines -UMU- and
             * incident angle cosines -GMU- to get directional emissivity at user angle cosines -UMU-.
             */
            dref = 0.;
            for (jg = 1; jg <= NMUG; jg++) {
              sum = 0.;
              for (k = 1; k <= NMUG/2; k++) {
                sum += GWT(k) * GMU(k) *
		  c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, UMU(iu), GMU(k),
					 M_PI*GMU(jg), ds->flag.brdf_type, &ds->brdf, callnum );
              }
              dref += GWT(jg)*sum;
            }
            EMU(iu) = 1.-dref;
          }
        }
      }
    }
  }

  return;
}

/*============================= end of c_surface_bidir() ================*/

/*============================= c_interp_eigenvec() =====================*/

/*
   Interpolate eigenvectors to user angles; eq SD(8)

   Called by- c_disort, c_albtrans
 --------------------------------------------------------------------*/

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
                       double       *ylmu)
{
  register int
    iq,iu,jq,l;
  double
    sum;

  for (iq = 1; iq <= ds->nstr; iq++) {
    for (l = mazim; l <= ds->nstr-1; l++) {
      /*
       * Inner sum in SD(8) times all factors in outer sum but PLM(mu)
       */
      sum = 0.;
      for (jq = 1; jq <= ds->nstr; jq++) {
        sum += CWT(jq)*YLMC(l,jq)*EVECC(jq,iq);
      }
      WK(l+1) = .5*GL(l,lc)*sum;
    }
    /*
     * Finish outer sum in SD(8) and store eigenvectors
     */
    for (iu = 1; iu <= ds->numu; iu++) {
      sum = 0.;
      for (l = mazim; l <=ds->nstr-1; l++) {
        sum += WK(l+1)*YLMU(l,iu);
      }
      if (iq <= nn) {
        GU(iu,nn+iq,lc) = sum;
      }
      else {
        GU(iu,ds->nstr+1-iq,lc) = sum;
      }
    }
  }

  return;
}

/*============================= end of c_interp_eigenvec() ==============*/

/*============================= c_interp_source() =======================*/

/*
    Interpolates source functions to user angles, eq. STWL(30)

    I N P U T      V A R I A B L E S:

       ds     :  Disort state variables
       cwt    :  Weights for Gauss quadrature over angle cosine
       delm0  :  Kronecker delta, delta-sub-m0
       gl     :  Delta-M scaled Legendre coefficients of phase function
                 (including factors 2l+1 and single-scatter albedo)
       mazim  :  Order of azimuthal component
       oprim  :  Single scattering albedo
       xr     :  Expansion of thermal source function, eq. STWL(24d); xr[].zero, xr[].one (see cdisort.h)
       ylm0   :  Normalized associated Legendre polynomial at the beam angle
       ylmc   :  Normalized associated Legendre polynomial at the quadrature angles
       ylmu   :  Normalized associated Legendre polynomial at the user angles
       zbs0   :  Solution vectors z-sub-zero of Eq. KS(10-11), used if pseudo-spherical
       zbs1   :  Solution vectors z-sub-one  of Eq. KS(10-11), used if pseudo-spherical
       zbsa   :  Alfa coefficient in Eq. KS(7), used if pseudo-spherical
       zee    :  Solution vectors Z-sub-zero, Z-sub-one of eq. SS(16), STWL(26a,b)
       zj     :  Solution vector Z-sub-zero after solving eq. SS(19), STWL(24b)
       zjg    :  Right-hand side vector  X-sub-zero in eq. KS(10), also the solution vector
                 Z-sub-zero after solving that system for a general source constant over a layer

    O U T P U T     V A R I A B L E S:

       zbeam  :  Incident-beam source function at user angles
       zu     :  Components 0 and 1 of a linear-in-optical-depth-dependent source (approximating the Planck emission source)
       zgu    :  General source function at user angles

   I N T E R N A L    V A R I A B L E S:

       psi  :   psi[].zero: Sum just after square bracket in eq. SD(9)
                psi[].one:  Sum in eq. STWL(31d)

   Called by- c_disort
 -------------------------------------------------------------------*/

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
                     disort_pair    *zu)
{
  register int
    iq,iu,jq;
  double
    fact,psum,psum0,psum1,sum,sum0,sum1;

  if (ds->bc.fbeam > 0.) {
    /*
     * Beam source terms; eq. SD(9)
     */
    if ( ds->flag.spher == TRUE ) {
      for (iq = mazim; iq <= ds->nstr-1; iq++) {
	psum0 = 0.;
	psum1 = 0.;
	for (jq = 1; jq <= ds->nstr; jq++) {
	  psum0 +=  CWT(jq)*YLMC(iq,jq)*ZBS0(jq);
	  psum1 +=  CWT(jq)*YLMC(iq,jq)*ZBS1(jq);
	}
	PSI0(iq+1) = 0.5*GL(iq,lc)*psum0;
	PSI1(iq+1) = 0.5*GL(iq,lc)*psum1;
      }
      for (iu = 1; iu <= ds->numu; iu++) {
	sum0 = 0.;
	sum1 = 0.;
	for (iq = mazim; iq <= ds->nstr-1; iq++) {
	  sum0 += YLMU(iq,iu)*PSI0(iq+1);
	  sum1 += YLMU(iq,iu)*PSI1(iq+1);
	}
	ZB0U(iu,lc) = sum0 + ZB0U(iu,lc);
	ZB1U(iu,lc) = sum1 + ZB1U(iu,lc);
	ZBAU(iu,lc) = zbsa;
      }
    }
    else {
      for (iq = mazim; iq <= ds->nstr-1; iq++) {
	psum = 0.;
	for (jq = 1; jq <= ds->nstr; jq++) {
	  psum += CWT(jq)*YLMC(iq,jq)*ZJ(jq);
	}
	PSI0(iq+1) = .5*GL(iq,lc)*psum;
      }
      fact = (2.-delm0)*ds->bc.fbeam/(4.*M_PI);
      for (iu = 1; iu <= ds->numu; iu++) {
	sum = 0.;
	for (iq = mazim; iq <= ds->nstr-1; iq++) {
	  sum += YLMU(iq,iu)*(PSI0(iq+1)+fact*GL(iq,lc)*YLM0(iq));
	}
	ZBEAM(iu,lc) = sum;
      }
    }
  }
  if (ds->flag.general_source > 0.) {
    /*
     * General source; eq. SD(9), KS(13)
     */
    for (iq = mazim; iq <= ds->nstr-1; iq++) {
      psum0 = 0.;
      for (jq = 1; jq <= ds->nstr; jq++) {
	psum0 +=  CWT(jq)*YLMC(iq,jq)*ZJG(jq);
      }
      PSI0(iq+1) = 0.5*GL(iq,lc)*psum0;
    }
    for (iu = 1; iu <= ds->numu; iu++) {
      sum0 = 0.;
      for (iq = mazim; iq <= ds->nstr-1; iq++) {
	sum0 += YLMU(iq,iu)*PSI0(iq+1);
      }
      ZGU(iu,lc) = sum0 + GENSRCU(mazim,lc,iu);
    }
  }

  if (ds->flag.planck && mazim == 0) {
    /*
     * Thermal source terms, STWJ(27c), STWL(31c)
     */
    for (iq = mazim; iq <=ds->nstr-1; iq++) {
      psum0 = 0.;
      psum1 = 0.;
      for (jq = 1; jq <= ds->nstr; jq++) {
        psum0 += CWT(jq)*YLMC(iq,jq)*Z0(jq);
        psum1 += CWT(jq)*YLMC(iq,jq)*Z1(jq);
      }
      PSI0(iq+1) = .5*GL(iq,lc)*psum0;
      PSI1(iq+1) = .5*GL(iq,lc)*psum1;
    }
    for (iu = 1; iu <= ds->numu; iu++) {
      sum0 = 0.;
      sum1 = 0.;
      for (iq = mazim; iq <= ds->nstr-1; iq++) {
        sum0 += YLMU(iq,iu)*PSI0(iq+1);
        sum1 += YLMU(iq,iu)*PSI1(iq+1);
      }
      Z0U(iu,lc) = sum0+(1.-OPRIM(lc))*XR0(lc);
      Z1U(iu,lc) = sum1+(1.-OPRIM(lc))*XR1(lc);
    }
  }

  return;
}

/*============================= end of c_interp_source() ================*/


/*============================= c_interp_coefficients_beam_source =======*/

/*
     Find coefficients at user angle, necessary for later use in
     c_interp_source()
*/

/*

    I N P U T      V A R I A B L E S:

       cmu    :   Computational polar angles
       chtau  :   The optical depth in spherical geometry.
       delmo  :   Kronecker delta, delta-sub-m0
       fbeam  :   incident beam radiation at top
       gl     :   Phase function Legendre coefficients multiplied by (2l+1) and single-scatter albedo
       lc:    :   layer index
       mazim  :   order of azimuthal component
       nstr   :   number of streams
       numu   :   number of user angles
       taucpr :   delta-m-scaled optical depth
       xba    :   alfa in eq. KS(7) 
       ylmu   :   Normalized associated Legendre polynomial at the user angles -umu-
       ylm0   :   Normalized associated Legendre polynomial at the beam angle

    O U T P U T     V A R I A B L E S:

       zb0u   :   x-sub-zero in KS(7) at user angles -umu-
       zb1u   :   x-sub-one in KS(7) at user angles -umu-
       zju    :  Solution vector Z-sub-zero after solving eq. SS(19), STWL(24b), at user angles -umu-

   Called by- c_disort

*/

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
				       double         *zju,
				       double         *ylm0,
				       double         *ylmu)
{
  register int 
    iu,k;
  double 
    deltat,sum,q0a,q2a,q0,q2;
  
  /*     Calculate x-sub-zero in STWJ(6d) */
  deltat = TAUCPR(lc) - TAUCPR(lc-1);

  q0a = exp(-CHTAU(lc-1));
  q2a = exp(-CHTAU(lc));
     
  for (iu = 1; iu <= numu; iu++) {
    sum = 0.0;
    for (k = mazim; k <= nstr-1; k++) {
      sum = sum + GL(k,lc)*YLMU(k,iu)*YLM0(k);
    }
    ZJU(iu) = (2.0-delm0)*fbeam*sum/(4.0*M_PI);
  }

  for (iu = 1; iu <= numu; iu++) {
     
    q0 = q0a*ZJU(iu);
    q2 = q2a*ZJU(iu);
     
    /*     x-sub-zero and x-sub-one in Eqs. KS(48-49)   */

    ZB1U(iu,lc)=(1./deltat)*(q2*exp(XBA(lc)*TAUCPR(lc))
			     -q0*exp(XBA(lc)*TAUCPR(lc-1)));
    ZB0U(iu,lc) = q0*exp(XBA(lc)*TAUCPR(lc-1))-ZB1U(iu,lc)*TAUCPR(lc-1);
  }

  return;

}
/*============================= end c_interp_coefficients_beam_source ===*/

/*============================= c_set_coefficients_beam_source() ========*/

/*
       Set coefficients in ks(7) for beam source

    I N P U T      V A R I A B L E S:

       cmu    :   Computational polar angles
       ch     :   The Chapman-factor to correct for pseudo-spherical geometry in the direct beam.
       chtau  :   The optical depth in spherical geometry.
       delmo  :   Kronecker delta, delta-sub-m0
       fbeam  :   incident beam radiation at top
       gl     :   Phase function Legendre coefficients multiplied by (2l+1) and single-scatter albedo
       lc:    :   layer index
       mazim  :   order of azimuthal component
       nstr   :   number of streams
       taucpr :   delta-m-scaled optical depth
       ylmc   :   Normalized associated Legendre polynomial at the quadrature angles -cmu-
       ylm0   :   Normalized associated Legendre polynomial at the beam angle

    O U T P U T     V A R I A B L E S:

       xba    :   alfa in eq. KS(7) 
       xb0    :   x-sub-zero in KS(7)
       xb1    :   x-sub-one in KS(7)
       zj     :  Solution vector Z-sub-zero after solving eq. SS(19), STWL(24b)

   Called by- c_disort
 -------------------------------------------------------------------*/

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
				    double       *zj)
{

  register int 
    iq,k;
  double 
    deltat,sum,q0a,q2a,q0,q2;
  static double
    big;

  big    = sqrt(DBL_MAX)/1.e+10;

  /*     Calculate x-sub-zero in STWJ(6d)   */

  for (iq = 1; iq <= nstr; iq++) {
    sum = 0;
    for (k = mazim; k <= nstr-1; k++) {
      sum += GL(k,lc)*YLMC(k,iq)*YLM0(k);
    }
    ZJ(iq) = (2.-delm0)*fbeam*sum/(4.*M_PI);
  }

  q0a = exp( -CHTAU(lc-1) );
  q2a = exp( -CHTAU(lc) );

  /*     Calculate alfa coefficient  */
     
  deltat = TAUCPR(lc) - TAUCPR(lc-1);

  XBA(lc) = 1./CH(lc);
        
  if ( fabs(XBA(lc)) > big  &&  TAUCPR(lc) > 1.)  XBA(lc) = 0.0;

  if( fabs(XBA(lc)*TAUCPR(lc)) > log(big))	  XBA(lc) = 0.0;

  /*     Dither alfa if it is close to one of the quadrature angles */

  if (  fabs(XBA(lc)) > 0.00001 ) {
    for (iq = 1; iq <= nstr/2; iq++) {
      if (fabs((fabs(XBA(lc))-1.0/CMU(iq))/XBA(lc) ) < 0.05 ) XBA(lc) = XBA(lc) * 1.001;      
    }
  }

  for (iq = 1; iq <= nstr; iq++) {

    q0 = q0a * ZJ(iq);
    q2 = q2a * ZJ(iq);
     
    /*     x-sub-zero and x-sub-one in Eqs. KS(48-49)   */   
       
    XB1(iq,lc) = (1.0/deltat)*(q2*exp(XBA(lc)*TAUCPR(lc)) - q0*exp(XBA(lc)*TAUCPR(lc-1)));    
    XB0(iq,lc) = q0 * exp(XBA(lc)*TAUCPR(lc-1)) - XB1(iq,lc)*TAUCPR(lc-1);

  }
  return;
}
/*============================= end c_set_coefficients_beam_source() ====*/


/*============================= c_upbeam() ==============================*/

/*
   Finds the incident-beam particular solution of SS(18), STWL(24a)

   I N P U T    V A R I A B L E S:

       ds     :  Disort state variables
       cc     :  C-sub-ij in eq. SS(5)
       cmu    :  Abscissae for Gauss quadrature over angle cosine
       delm0  :  Kronecker delta, delta-sub-m0
       gl     :  Delta-M scaled Legendre coefficients of phase function
                 (including factors 2l+1 and single-scatter albedo)
       mazim  :  Order of azimuthal component
       ylm0   :  Normalized associated Legendre polynomial at the beam angle
       ylmc   :  Normalized associated Legendre polynomial at the quadrature angles

   O U T P U T    V A R I A B L E S:

       zj     :  Right-hand side vector X-sub-zero in SS(19),STWL(24b);
                 also the solution vector Z-sub-zero after solving that system
       zz     :  Permanent storage for zj, but re-ordered

   I N T E R N A L    V A R I A B L E S:

       array  :  Coefficient matrix in left-hand side of eq. SS(19), STWL(24b)
       ipvt   :  Integer vector of pivot indices required by LINPACK
       wk     :  Scratch array required by LINPACK

   Called by- c_disort
   Calls- c_sgeco, c_errmsg, c_sgesl
 -------------------------------------------------------------------*/

#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*ds->nstr]

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
              double       *zz)
{
  register int
    iq,jq,k;
  double
    rcond,sum;

  for (iq = 1; iq <= ds->nstr; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      ARRAY(iq,jq) = -CC(iq,jq);
    }
    ARRAY(iq,iq) += 1.+CMU(iq)/ds->bc.umu0;
    sum = 0.;
    for (k = mazim; k <=ds->nstr-1; k++) {
      sum += GL(k,lc)*YLMC(k,iq)*YLM0(k);
    }
    ZJ(iq) = (2.-delm0)*ds->bc.fbeam*sum/(4.*M_PI);
  }
  /*
   * Find L-U (lower/upper triangular) decomposition of ARRAY and see if it is nearly singular
   * (NOTE:  ARRAY is altered)
   */
  rcond = 0.;
  c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);

  if (1.+rcond == 1.) {
    c_errmsg("upbeam--sgeco says matrix near singular",DS_WARNING);
  }

  /*
   * Solve linear system with coeff matrix ARRAY (assumed already L-U decomposed) and R.H. side(s) ZJ;
   * return solution(s) in ZJ
   */
  c_sgesl(array,ds->nstr,ds->nstr,ipvt,zj,0);
  for (iq = 1; iq <= nn; iq++) {
    ZZ(nn+iq,  lc) = ZJ(iq);
    ZZ(nn-iq+1,lc) = ZJ(iq+nn);
  }

  return;
}

/*============================= end of c_upbeam() =======================*/

/*============================= c_upbeam_pseudo_spherical() =============*/

/*

       Finds the particular solution of beam source KS(10-11)

     Routines called:  sgeco, sgesl

   I N P U T     V A R I A B L E S:

       cc     :  capital-c-sub-ij in Eq. SS(5)
       cmu    :  abscissae for gauss quadrature over angle cosine
       xb0    :  EXPansion of beam source function Eq. KS(7)
       xb1    :  EXPansion of beam source function Eq. KS(7)
       xba    :  EXPansion of beam source function Eq. KS(7)
       (remainder are 'disort' input variables)

    O U T P U T    V A R I A B L E S:

       zbs0     :  solution vectors z-sub-zero of Eq. KS(10-11)
       zbs1     :  solution vectors z-sub-one  of Eq. KS(10-11)
       zbsa     :  alfa coefficient in Eq. KS(7)
       zbeam0,  :  permanent storage for -zbs0,zbs1,zbsa-, but rD-ordered
        zbeam1,
        zbeama
 
   I N T E R N A L    V A R I A B L E S:

       array  :  coefficient matrix in left-hand side of Eq. KS(10)
       ipvt   :  integer vector of pivot indices required by *linpack*
       wk     :  scratch array required by *linpack*

   Called by- c_disort
   Calls- c_sgeco, c_errmsg, c_sgesl
 -------------------------------------------------------------------*/

#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*ds->nstr]

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
			       double       *zbeama)
{

  register int
    iq,jq;
  double
    rcond,rmin;


  for (iq = 1; iq <= ds->nstr; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      ARRAY(iq,jq) = -CC(iq,jq);
    }
    ARRAY(iq,iq) += 1.+XBA(lc)*CMU(iq);
    *zbsa     = XBA(lc);
    ZBS1(iq) = XB1(iq,lc);
  }

  /*
   * Find L-U (lower/upper triangular) decomposition of ARRAY and see
   * if it is nearly singular
   * (NOTE: ARRAY is altered)
   */

  rcond = 0.;
  c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);

  if (1.+rcond == 1.) {
    c_errmsg("upbeam_pseudo_spherical--sgeco says matrix near singular",
	     DS_WARNING);
  }
     
  rmin = 1.0e-4;
  if ( rcond < rmin ) {
    /*     Dither alpha if rcond to small   */
    if(XBA(lc) ==0.0)       XBA(lc)=0.000000005;

    XBA(lc) = XBA(lc) * 1.00000005;

    for (iq = 1; iq <= ds->nstr; iq++) {
      for (jq = 1; jq <= ds->nstr; jq++) {
	ARRAY(iq,jq) = -CC(iq,jq);
      }	
      ARRAY(iq,iq) += 1.0+XBA(lc)*CMU(iq);
      *zbsa     = XBA(lc);
      ZBS1(iq) = XB1(iq,lc);
    }
    /*     Solve linear equations KS(10-11) with dithered alpha */
    rcond = 0.;
    c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);               
    if (1.+rcond == 1.) {
      c_errmsg("upbeam_pseudo_spherical--sgeco says matrix near singular",
	       DS_WARNING);
    }
  }

  for (iq = 1; iq <= ds->nstr; iq++)  WK(iq) = ZBS1(iq);
  c_sgesl( array, ds->nstr, ds->nstr, ipvt, wk, 0 );
          
  for (iq = 1; iq <= ds->nstr; iq++) {
    ZBS1(iq) = WK(iq);
    ZBS0(iq) = XB0(iq,lc) + CMU(iq) * ZBS1(iq);
  }

  for (iq = 1; iq <= ds->nstr; iq++)  WK(iq) = ZBS0(iq);
  c_sgesl( array, ds->nstr, ds->nstr, ipvt, wk, 0 );
  for (iq = 1; iq <= ds->nstr; iq++)  ZBS0(iq) = WK(iq);

  /*   ... and now some index gymnastic for the inventive ones...  */

  ZBEAMA(lc)            = *zbsa;
  for (iq = 1; iq <= nn; iq++) {
    ZBEAM0( iq+nn, lc )   = ZBS0( iq );
    ZBEAM1( iq+nn, lc )   = ZBS1( iq );
    ZBEAM0( nn+1-iq, lc ) = ZBS0( iq+nn );
    ZBEAM1( nn+1-iq,lc )  = ZBS1( iq+nn );
  }

 return;

}
  

/*============================= end of c_upbeam_pseudo_spherical() ======*/

/*============================= c_upbeam_general_source() ===============*/

/*
   Finds the incident-beam particular solution of SS(18), STWL(24a)

   I N P U T    V A R I A B L E S:

       ds     :  Disort state variables
       cc     :  C-sub-ij in eq. SS(5)
       cmu    :  Abscissae for Gauss quadrature over angle cosine
       delm0  :  Kronecker delta, delta-sub-m0
       gl     :  Delta-M scaled Legendre coefficients of phase function
                 (including factors 2l+1 and single-scatter albedo)
       mazim  :  Order of azimuthal component
       ylm0   :  Normalized associated Legendre polynomial at the beam angle
       ylmc   :  Normalized associated Legendre polynomial at the quadrature angles

   O U T P U T    V A R I A B L E S:

       zjg    :  Right-hand side vector  X-sub-zero in eq. KS(10), also the solution vector
                 Z-sub-zero after solving that system for a general source constant over a layer
       zzg    :  Permanent storage for zjg, but re-ordered

   I N T E R N A L    V A R I A B L E S:

       array  :  Coefficient matrix in left-hand side of eq. SS(19), STWL(24b)
       ipvt   :  Integer vector of pivot indices required by LINPACK
       wk     :  Scratch array required by LINPACK

   Called by- c_disort
   Calls- c_sgeco, c_errmsg, c_sgesl
 -------------------------------------------------------------------*/

#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*ds->nstr]

void c_upbeam_general_source(disort_state *ds,
			     int           lc,
			     int           maz,
			     double       *array,
			     double       *cc,
			     int          *ipvt,
			     int           nn,
			     double       *wk,
			     double       *zjg,
			     double       *zzg)
{
  register int
    iq,jq;
  double
    rcond;

  for (iq = 1; iq <= nn; iq++) {
    ZJG(iq )    = GENSRC(maz,lc,nn+iq);
    ZJG(nn+iq)  = GENSRC(maz,lc,nn+1-iq);
  }

  for (iq = 1; iq <= ds->nstr; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      ARRAY(iq,jq) = -CC(iq,jq);
    }
    ARRAY(iq,iq) = 1 + ARRAY(iq,iq);
  }

  /*
   * Find L-U (lower/upper triangular) decomposition of ARRAY and see if it is nearly singular
   * (NOTE:  ARRAY is altered)
   */
  rcond = 0.;
  c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);

  if (1.+rcond == 1.) {
    c_errmsg("upbeam_general_source--sgeco says matrix near singular",DS_WARNING);
  }

  /*
   * Solve linear system with coeff matrix ARRAY (assumed already L-U decomposed) and R.H. side(s) ZJG;
   * return solution(s) in ZJG
   */
  c_sgesl(array,ds->nstr,ds->nstr,ipvt,zjg,0);
  for (iq = 1; iq <= nn; iq++) {
    ZZG(nn+iq,  lc) = ZJG(iq);
    ZZG(nn-iq+1,lc) = ZJG(iq+nn);
  }
 
  return;
}

/*============================= end of c_upbeam_general_source() ========*/

/*============================= c_upisot() ==============================*/

/*
    Finds the particular solution of thermal radiation of STWL(25)

    I N P U T     V A R I A B L E S:

       ds     :  Disort state variables
       cc     :  C-sub-ij in eq. SS(5), STWL(8b)
       cmu    :  Abscissae for Gauss quadrature over angle cosine
       oprim  :  Delta-M scaled single scattering albedo
       xr     :  Expansion coefficient b-sub-zero, b-sub-one of thermal source function, eq. STWL(24c)

    O U T P U T    V A R I A B L E S:

       zee    :  Solution vectors Z-sub-zero, Z-sub-one of eq. SS(16), STWL(26a,b)
       plk    :  Permanent storage for zee, but re-ordered

   I N T E R N A L    V A R I A B L E S:

       array  :  Coefficient matrix in left-hand side of eq. SS(16)
       ipvt   :  Integer vector of pivot indices required by LINPACK
       wk     :  Scratch array required by LINPACK

   Called by- c_disort
   Calls- c_sgeco, c_errmsg, c_sgesl
 -------------------------------------------------------------------*/

#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*ds->nstr]

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
              disort_pair  *plk)
{
  register int
    iq,jq;
  double
    rcond;

  for (iq = 1; iq <= ds->nstr; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      ARRAY(iq,jq) = -CC(iq,jq);
    }
    ARRAY(iq,iq) += 1.;
    Z1(iq) = (1.-OPRIM(lc))*XR1(lc);
  }
  /*
   * Solve linear equations: same as in upbeam, except zj replaced by z1 and z0
   */
  rcond = 0.;
  c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);

  if (1.+rcond == 1.) {
    c_errmsg("upisot--sgeco says matrix near singular",DS_WARNING);
  }
  
  for (iq = 1; iq <= ds->nstr; iq++) {
    /* Need to use WK() as a buffer, since Z1 is part of a structure */
    WK(iq) = Z1(iq);
  }
  c_sgesl(array,ds->nstr,ds->nstr,ipvt,wk,0);
  for (iq = 1; iq <= ds->nstr; iq++) {
    Z1(iq) = WK(iq);
  }

  for (iq = 1; iq <= ds->nstr; iq++) {
    Z0(iq) = (1.-OPRIM(lc))*XR0(lc)+CMU(iq)*Z1(iq);
  }

  for (iq = 1; iq <= ds->nstr; iq++) {
    /* Need to use WK() as a buffer, since Z0 is part of a structure */
    WK(iq) = Z0(iq);
  }
  c_sgesl(array,ds->nstr,ds->nstr,ipvt,wk,0);
  for (iq = 1; iq <= ds->nstr; iq++) {
    Z0(iq) = WK(iq);
  }
  for (iq = 1; iq <= nn; iq++) {
    ZPLK0(nn+iq,  lc) = Z0(iq   );
    ZPLK1(nn+iq,  lc) = Z1(iq   );
    ZPLK0(nn-iq+1,lc) = Z0(iq+nn);
    ZPLK1(nn-iq+1,lc) = Z1(iq+nn);
  }

  return;
}

/*============================= end of c_upisot() =======================*/

/*============================= c_user_intensities() ====================*/

/*
   Computes intensity components at user output angles for azimuthal
   expansion terms in eq. SD(2), STWL(6)

   I N P U T    V A R I A B L E S:

       ds     :  Disort state variables
       bplanck:  Integrated Planck function for emission from
                 bottom boundary
       cmu    :  Abscissae for Gauss quadrature over angle cosine
       cwt    :  Weights for Gauss quadrature over angle cosine
       delm0  :  Kronecker delta, delta-sub-M0
       emu    :  Surface directional emissivity (user angles)
       expbea :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
       gc     :  Eigenvectors at polar quadrature angles, SC(1)
       gu     :  Eigenvectors interpolated to user polar angles
                    (i.e., G in eq. SC(1) )
       kk     :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       layru  :  Layer number of user level UTAU
       ll     :  Constants of integration in eq. SC(1), obtained
                 by solving scaled version of eq. SC(5);
                 exponential term of eq. SC(12) not included
       lyrcut :  Logical flag for truncation of computational layer
       mazim  :  Order of azimuthal component
       ncut   :  Total number of computational layers considered
       nn     :  Order of double-Gauss quadrature (NSTR/2)
       rmu    :  Surface bidirectional reflectivity (user angles)
       taucpr :  Cumulative optical depth (delta-M-Scaled)
       tplanck:  Integrated Planck function for emission from
                 top boundary
       utaupr :  Optical depths of user output levels in delta-M
                 coordinates;  equal to UTAU if no delta-M
       zgu    :  General source function at user angles
       zu     :  Z-sub-zero, Z-sub-one in eq. SS(16) interpolated to user angles from an equation derived from SS(16),
                 Y-sub-zero, Y-sub-one on STWL(26b,a); zu[].zero, zu[].one (see cdisort.h)
       zz     :  Beam source vectors in eq. SS(19), STWL(24b)
       zzg    :  Beam source vectors in eq. KS(10)for a general source constant over a layer
       plk    :  Thermal source vectors z0,z1 by solving eq. SS(16),
                 Y-sub-zero,Y-sub-one in STWL(26)
       zbeam  :  Incident-beam source vectors


    O U T P U T    V A R I A B L E S:

       uum    :  Azimuthal components of the intensity in eq. STWJ(5),
                 STWL(6)

    I N T E R N A L    V A R I A B L E S:

       bnddir :  Direct intensity down at the bottom boundary
       bnddfu :  Diffuse intensity down at the bottom boundary
       bndint :  Intensity attenuated at both boundaries, STWJ(25-6)
       dtau   :  Optical depth of a computational layer
       lyrend :  End layer of integration
       lyrstr :  Start layer of integration
       palint :  Intensity component from parallel beam
       plkint :  Intensity component from planck source
       wk     :  Scratch vector for saving exp evaluations

       All the exponential factors (exp1, expn,... etc.)
       come from the substitution of constants of integration in
       eq. SC(12) into eqs. S1(8-9).  They all have negative
       arguments so there should never be overflow problems.

   Called by- c_disort
 -------------------------------------------------------------------*/

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
                        double         *uum)
{
  register int
    negumu,
    iq,iu,jq,lc,lu,lyrend,lyrstr,lyu;
  double
    alfa,bnddfu,bnddir,bndint,
    denom,dfuint,dtau,dtau1,dtau2,
    exp0=0,exp1=0,exp2=0,expn,
    f0n,f1n,fact,genint,
    palint,plkint,sgn;

  /*
   * Incorporate constants of integration into interpolated eigenvectors
   */
  for (lc = 1; lc <= ncut; lc++) {
    for (iq = 1; iq <= ds->nstr; iq++) {
      for (iu = 1; iu <= ds->numu; iu++) {
        GU(iu,iq,lc) *= LL(iq,lc);
      }
    }
  }

  /*
   * Loop over levels at which intensities are desired ('user output levels')
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    if (ds->bc.fbeam > 0.) {
      exp0 = exp(-UTAUPR(lu)/ds->bc.umu0);
    }
    lyu = LAYRU(lu);
    /*
     * Loop over polar angles at which intensities are desired
     */
    for (iu = 1; iu <= ds->numu; iu++) {
      if (lyrcut && lyu > ncut) {
        continue;
      }
      negumu = (UMU(iu) < 0.);
      if (negumu) {
        lyrstr = 1;
        lyrend = lyu-1;
        sgn    = -1.;
      }
      else {
        lyrstr = lyu+1;
        lyrend = ncut;
        sgn    = 1.;
      }

      /*
       * For downward intensity, integrate from top to LYU-1 in eq. S1(8); for upward,
       * integrate from bottom to LYU+1 in S1(9)
       */
      genint = 0.;
      palint = 0.;
      plkint = 0.;
      for (lc = lyrstr; lc <= lyrend; lc++) {
        dtau = DTAUCPR(lc);
        exp1 = exp((UTAUPR(lu)-TAUCPR(lc-1))/UMU(iu));
        exp2 = exp((UTAUPR(lu)-TAUCPR(lc  ))/UMU(iu));

        if (ds->flag.planck && mazim == 0) {
          /*
           * Eqs. STWL(36b,c, 37b,c)
           */
          f0n     = sgn*(exp1-exp2);
          f1n     = sgn*((TAUCPR(lc-1)+UMU(iu))*exp1
                        -(TAUCPR(lc  )+UMU(iu))*exp2);
          plkint += Z0U(iu,lc)*f0n+Z1U(iu,lc)*f1n;
        }

        if (ds->bc.fbeam > 0.) {
	  if ( ds->flag.spher == TRUE ) {
	    denom  =  sgn*1.0/(ZBAU(iu,lc)*UMU(iu)+1.0);
	    palint += (ZB0U(iu,lc)*denom*(exp(-ZBAU(iu,lc)*TAUCPR(lc-1)) *exp1
					  -exp(-ZBAU(iu,lc)*TAUCPR(lc)) *exp2 ) 
		       +ZB1U(iu,lc)*denom*((TAUCPR(lc-1)+sgn*denom*UMU(iu))
					   *exp(-ZBAU(iu,lc)*TAUCPR(lc-1)) *exp1
					   -(TAUCPR(lc)+sgn*denom*UMU(iu) )
					   *exp(-ZBAU(iu,lc)*TAUCPR(lc))*exp2));
	  }
	  else {
	    denom = 1.+UMU(iu)/ds->bc.umu0;
	    if (fabs(denom) < 0.0001) {
	      /*
	       * L'Hospital limit
	       */
	      expn = (dtau/ds->bc.umu0)*exp0;
	    }
	    else {
	      expn = (exp1*EXPBEA(lc-1)
		      -exp2*EXPBEA(lc  ))*sgn/denom;
	    }
	    palint += ZBEAM(iu,lc)*expn;
	  }
        }
	if ( ds->flag.general_source ) {
          genint += ZGU(iu,lc)*sgn*(exp1-exp2);
	}
        /*
         * KK is negative
         */
        for (iq = 1; iq <= nn; iq++) {
          WK(iq) = exp(KK(iq,lc)*dtau);
          denom  = 1.+UMU(iu)*KK(iq,lc);
          if (fabs(denom) < 0.0001) {
            /*
             * L'Hospital limit
             */
            expn = (dtau/UMU(iu))*exp2;
          }
          else {
            expn = sgn*(exp1*WK(iq)-exp2)/denom;
          }
          palint += GU(iu,iq,lc)*expn;
        }

        /*
         * KK is positive
         */
        for (iq = nn+1; iq <= ds->nstr; iq++) {
          denom = 1.+UMU(iu)*KK(iq,lc);
          if (fabs(denom) < 0.0001) {
            /*
             * L'Hospital limit
             */
            expn = -(dtau/UMU(iu))*exp1;
          }
          else {
            expn = sgn*(exp1-exp2*WK(ds->nstr+1-iq))/denom;
          }
          palint += GU(iu,iq,lc)*expn;
        }
      }

      /*
       * Calculate contribution from user output level to next computational level
       */
      dtau1 = UTAUPR(lu)-TAUCPR(lyu-1);
      dtau2 = UTAUPR(lu)-TAUCPR(lyu  );

      if ((fabs(dtau1) >= 1.e-6 || !negumu) && (fabs(dtau2) >= 1.e-6 ||  negumu)) {
        if(negumu) {
          exp1 = exp(dtau1/UMU(iu));
        }
        else {
          exp2 = exp(dtau2/UMU(iu));
        }
        if (ds->bc.fbeam > 0.) {
	  if ( ds->flag.spher == TRUE ) {
	    if ( negumu ) {	     
	      expn = exp1;
	      alfa = ZBAU(iu,lyu);
	      denom = (-1.0/(alfa*UMU(iu)+1.));	        
	      palint += ZB0U(iu,lyu)*denom*(-exp(-alfa*UTAUPR(lu))
					    + expn*exp(-alfa*TAUCPR(lyu-1)))
		+ZB1U(iu,lyu)*denom*( -(UTAUPR(lu)-UMU(iu)*denom)*exp(-alfa*UTAUPR(lu))
				      +(TAUCPR(lyu-1)-UMU(iu)*denom)*expn*exp(-alfa*TAUCPR(lyu-1)));
	    }
	    else {
	      expn = exp2;
	      alfa = ZBAU(iu,lyu);
	      denom = (1.0/(alfa*UMU(iu)+1.0));
	      palint += ZB0U(iu,lyu)*denom*(exp(-alfa*UTAUPR(lu))
					    -exp(-alfa*TAUCPR(lyu))*expn)
		+ZB1U(iu,lyu)*denom*( (UTAUPR(lu) +UMU(iu)*denom)*exp(-alfa*UTAUPR(lu))
				      -(TAUCPR(lyu)+UMU(iu)*denom)*exp(-alfa*TAUCPR(lyu))*expn );	          
	    }
	  }
	  else {
	    denom = 1.+UMU(iu)/ds->bc.umu0;
	    if (fabs(denom) < 0.0001) {
	      expn = (dtau1/ds->bc.umu0)*exp0;
	    }
	    else if (negumu) {
	      expn = (exp0-EXPBEA(lyu-1)*exp1)/denom;
	    }
	    else {
	      expn = (exp0-EXPBEA(lyu  )*exp2)/denom;
	    }
	    palint += ZBEAM(iu,lyu)*expn;
	  }
        }
	if ( ds->flag.general_source ) {
          if (negumu) {
            expn = exp1;
          }
          else {
            expn = exp2;
          }
          genint += ZGU(iu,lyu)*(1.-expn);
	}
        /*
         * KK is negative
         */
        dtau = DTAUCPR(lyu);
        for (iq = 1; iq <= nn; iq++) {
          denom = 1.+UMU(iu)*KK(iq,lyu);
          if (fabs(denom) < 0.0001) {
            expn = -dtau2/UMU(iu)*exp2;
          }
          else if (negumu) {
            expn = (exp(-KK(iq,lyu)*dtau2)
                   -exp( KK(iq,lyu)*dtau )*exp1)/denom;
          }
          else {
            expn = (exp(-KK(iq,lyu)*dtau2)-exp2)/denom;
          }
          palint += GU(iu,iq,lyu)*expn;
        }

        /*
         * KK is positive
         */
        for (iq = nn+1; iq <= ds->nstr; iq++) {
          denom = 1.+UMU(iu)*KK(iq,lyu);
          if (fabs(denom) < 0.0001) {
            expn = -(dtau1/UMU(iu))*exp1;
          }
          else if (negumu) {
            expn = (exp(-KK(iq,lyu)*dtau1)-exp1)/denom;
          }
          else {
            expn = (exp(-KK(iq,lyu)*dtau1)
                   -exp(-KK(iq,lyu)*dtau )*exp2)/denom;
          }
          palint += GU(iu,iq,lyu)*expn;
        }

        if (ds->flag.planck && mazim == 0) {
          /*
           * Eqs. STWL (35-37) with tau-sub-n-1 replaced by tau for upward, and
           * tau-sub-n replaced by tau for downward directions
           */
          if (negumu) {
            expn = exp1;
            fact = TAUCPR(lyu-1)+UMU(iu);
          }
          else {
            expn = exp2;
            fact = TAUCPR(lyu  )+UMU(iu);
          }
          f0n     = 1.-expn;
          f1n     = UTAUPR(lu)+UMU(iu)-fact*expn;
          plkint += Z0U(iu,lyu)*f0n+Z1U(iu,lyu)*f1n;
        }
      }

      /*
       * Calculate intensity components attenuated at both boundaries.
       * NOTE: no azimuthal intensity component for isotropic surface
       */
      bndint = 0.;
      if (negumu && mazim == 0) {
        bndint = (ds->bc.fisot+tplanck)*exp(UTAUPR(lu)/UMU(iu));
      }
      else if (!negumu) {
        if (lyrcut || ( ds->flag.lamber && mazim > 0 ) ) {
          UUM(iu,lu) = palint+plkint;
          continue;
        }

        for (jq = nn+1; jq <= ds->nstr; jq++) {
          WK(jq) = exp(-KK(jq,ds->nlyr)*DTAUCPR(ds->nlyr));
        }
        bnddfu = 0.;
        for (iq = nn; iq >= 1; iq--) {
          dfuint = 0.;
          for (jq = 1; jq <= nn; jq++) {
            dfuint += GC(iq,jq,ds->nlyr)*LL(jq,ds->nlyr);
          }
          for (jq= nn+1; jq <= ds->nstr; jq++) {
            dfuint += GC(iq,jq,ds->nlyr)*LL(jq,ds->nlyr)*WK(jq);
          }
          if (ds->bc.fbeam > 0.) {
	    if ( ds->flag.spher == TRUE ) {
	      dfuint += exp(-ZBEAMA(ds->nlyr)*TAUCPR(ds->nlyr)) *
		(ZBEAM0(iq,ds->nlyr)+ZBEAM1(iq,ds->nlyr)*TAUCPR(ds->nlyr));
	    }
	    else {
	      dfuint += ZZ(iq,ds->nlyr)*EXPBEA(ds->nlyr);
	    }
          }
	  if ( ds->flag.general_source ) {
	    dfuint += ZZG(iq,ds->nlyr);
	  }
          dfuint += delm0*(ZPLK0(iq,ds->nlyr)+ZPLK1(iq,ds->nlyr)*TAUCPR(ds->nlyr));
          bnddfu += (1.+delm0)*RMU(iu,nn+1-iq)*CMU(nn+1-iq)*CWT(nn+1-iq)*dfuint;
        }
        bnddir = 0.;
        if (ds->bc.fbeam > 0. || ds->bc.umu0 >0.) {
          bnddir = ds->bc.umu0*ds->bc.fbeam/M_PI*RMU(iu,0)*EXPBEA(ds->nlyr);
        }
        bndint = (bnddfu+bnddir+delm0*EMU(iu)*bplanck+ds->bc.fluor)*exp((UTAUPR(lu)-TAUCPR(ds->nlyr))/UMU(iu));
      }
      UUM(iu,lu) = palint+plkint+bndint+genint;
    }
  }

  return;
}

/*============================= end of c_user_intensities() =============*/

/*============================= c_xi_func() =============================*/

/*
   Calculates Xi function of eq. STWL (72)

         I N P U T   V A R I A B L E S

   umu1,2    cosine of zenith angle_1, _2
   tau       optical thickness of the layer

   NOTE: Original Fortran version also had argument umu3, but was only
         called for the case umu2 == umu3, so these two arguments are
         fused together here to reduce conditional testing.

   Called by- c_secondary_scat
 -------------------------------------------------------------------*/

double c_xi_func(double umu1,
               double umu2,
               double tau)
{
  double
    exp1,x1;

  x1   = (umu2-umu1)/(umu2*umu1);
  exp1 = exp(-tau/umu1);

  if (x1 != 0.) {
    return ((tau*x1-1.)*exp(-tau/umu2)+exp1)/(x1*x1*umu1*umu2);
  }
  else {
    return tau*tau*exp1/(2.*umu1*umu2);
  }
}

/*============================= end of c_xi_func() ======================*/

/*============================= c_check_inputs() ========================*/

/*
 * Checks the input dimensions and variables
 *
 * Calls- c_write_bad_var, c_dref, c_errmsg
 * Called by- c_disort
 */

void c_check_inputs(disort_state *ds,
		    int           scat_yes,
		    int           deltam,
		    int           corint,
		    double       *tauc,
		    int           callnum)
{
  int
    inperr = FALSE;
  register int
    irmu,iu,j,k,lc,lu, nu;
  double
    flxalb,rmu,umumin;

  if (ds->nstr < 2 || ds->nstr%2 != 0) {
    inperr = c_write_bad_var(VERBOSE,"ds.nstr");
  }
  if (ds->nstr == 2) {
    c_errmsg("check_inputs()--2 streams not recommended;\n\nUse specialized 2-stream code c_twostr() instead",DS_WARNING);
  }
  if (ds->nlyr < 1) {
    inperr = c_write_bad_var(VERBOSE,"ds.nlyr");
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (DTAUC(lc) < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.dtauc");
    }
    if (SSALB(lc) < 0.0 || SSALB(lc) > 1.0) {
      inperr = c_write_bad_var(VERBOSE,"ds.ssalb");
    }
    if (ds->flag.ibcnd == GENERAL_BC) {
      if (ds->flag.planck) {
        if (lc == 1 && TEMPER(0) < 0.) {
          inperr = c_write_bad_var(VERBOSE,"ds.temper");
        }
        if (TEMPER(lc) < 0.) {
          inperr = c_write_bad_var(VERBOSE,"ds.temper");
        }
      }
    }
    else if (ds->flag.ibcnd == SPECIAL_BC) {
      ds->flag.planck = FALSE;
    }
    else {
      c_errmsg("check_inputs---unrecognized ds->flag.ibcnd",DS_ERROR);
    }
  }

  if (ds->nmom < 0 || (scat_yes  && ds->nmom < ds->nstr)) {
    inperr = c_write_bad_var(VERBOSE,"ds.nmom");
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    for (k = 0; k <= ds->nmom; k++) {
      if (PMOM(k,lc) < -1. || PMOM(k,lc) > 1.) {
        inperr = c_write_bad_var(VERBOSE,"PMOM(k,lc)");
      }
    }
  }

  if( ds->flag.spher == TRUE ) {
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (ds->ZD(lc) > ds->ZD(lc-1)) {
        inperr     = c_write_bad_var(ds->flag.quiet,"zd");
      }
    }
  }
 
  if (ds->flag.ibcnd == GENERAL_BC) {
    if (ds->flag.usrtau) {
      if (ds->ntau < 1) {
        inperr = c_write_bad_var(VERBOSE,"ds.ntau");
      }
      for (lu = 1; lu <= ds->ntau; lu++) {
	/* Do a relative check to see if we are just beyond the bottom boundary */
	/* This might happen due to numerical rounding off problems.  ak20110224*/
        if (fabs(UTAU(lu)-TAUC(ds->nlyr)) <= 1.e-6*TAUC(ds->nlyr)) {
          UTAU(lu) = TAUC(ds->nlyr);
        }
        if(UTAU(lu) < 0. || UTAU(lu) > TAUC(ds->nlyr)) {
          inperr = c_write_bad_var(VERBOSE,"ds.utau");
        }
      }
    }
  }

  if (ds->flag.usrang) {
    if (ds->numu < 0) {
      inperr = c_write_bad_var(VERBOSE,"ds.numu");
    }
    if (!ds->flag.onlyfl && ds->numu == 0) {
      inperr = c_write_bad_var(VERBOSE,"ds.numu");
    }
    nu = ds->numu;
    if (ds->flag.ibcnd == SPECIAL_BC ) nu = ds->numu/2;
    for (iu = 1; iu <= nu; iu++) {
      if (UMU(iu) < -1. || UMU(iu) > 1. || UMU(iu) == 0.) {
        inperr = c_write_bad_var(VERBOSE,"ds.umu");
      }
      if (ds->flag.ibcnd == SPECIAL_BC && UMU(iu) < 0.) {
        inperr = c_write_bad_var(VERBOSE,"ds.umu");
      }
      if (iu > 1) {
        if (UMU(iu) < UMU(iu-1)) {
          inperr = c_write_bad_var(VERBOSE,"ds.umu");
        }
      }
    }
  }

  if (!ds->flag.onlyfl && ds->flag.ibcnd != SPECIAL_BC) {
    if (ds->nphi <= 0) {
      inperr = c_write_bad_var(VERBOSE,"ds.nphi");
    }
    for (j=1; j <=ds->nphi; j++) {
      if (PHI(j) < 0. || PHI(j) > 360.) {
        inperr = c_write_bad_var(VERBOSE,"ds.phi");
      }
    }
  }

  if (ds->flag.ibcnd != GENERAL_BC && ds->flag.ibcnd != SPECIAL_BC) {
    inperr = c_write_bad_var(VERBOSE,"ds.flag.ibcnd");
  }

  if (ds->flag.ibcnd == GENERAL_BC) {
    if (ds->bc.fbeam < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.fbeam");
    }
    else if (ds->bc.fbeam > 0.) {
      umumin = 0.;
      if( ds->flag.spher == TRUE ) {
	umumin = -1.;
      }
      if (ds->bc.umu0 <= umumin || ds->bc.umu0 > 1.) {
        inperr = c_write_bad_var(VERBOSE,"ds.bc.umu0");
      }
      if (ds->bc.phi0 < 0. || ds->bc.phi0 > 360.) {
        inperr = c_write_bad_var(VERBOSE,"ds.bc.phi0");
      }
    }

    if (ds->bc.fisot < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.fisot");
    }

    if (ds->flag.lamber) {
      if (ds->bc.albedo < 0. || ds->bc.albedo > 1.) {
        inperr = c_write_bad_var(VERBOSE,"ds.bc.albedo");
      }
    }
    else {
      /*
       * Make sure flux albedo at dense mesh of incident angles does not assume unphysical values
       */
      for (irmu = 0; irmu <= 100; irmu++) {
        rmu    = (double)irmu*0.01;
        flxalb = c_dref(ds->wvnmlo, ds->wvnmhi, rmu, ds->flag.brdf_type, &ds->brdf, callnum);
        if (flxalb < 0. || flxalb > 1.) {
          inperr = c_write_bad_var(VERBOSE,"bidir_reflectivity()");
        }
      }
    }
  }
  else {
    if (ds->bc.albedo < 0. || ds->bc.albedo > 1.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.albedo");
    }
  }

  if (ds->flag.planck && ds->flag.ibcnd != SPECIAL_BC) {
    if (ds->wvnmlo < 0. || ds->wvnmhi <= ds->wvnmlo) {
      inperr = c_write_bad_var(VERBOSE,"ds.wvnmlo,hi");
    }
    if (ds->bc.temis < 0. || ds->bc.temis > 1.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.temis");
    }
    if (ds->bc.btemp < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.btemp");
    }
    if (ds->bc.ttemp < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.ttemp");
    }
  }

  if (ds->accur < 0. || ds->accur > 1.e-2) {
    inperr = c_write_bad_var(VERBOSE,"ds.accur");
  }

  if (inperr) {
    c_errmsg("DISORT--input and/or dimension errors",DS_ERROR);
  }

  if (ds->flag.planck && ds->flag.quiet == VERBOSE) {
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (fabs(TEMPER(lc)-TEMPER(lc-1)) > 10.) {
        c_errmsg("check_inputs--vertical temperature step may be too large for good accuracy",DS_WARNING);
      }
    }
  }
  if(!corint && (!ds->flag.onlyfl && ds->bc.fbeam > 0. && scat_yes && deltam)) {
    c_errmsg("check_inputs--intensity correction is off;\nintensities may be less accurate",DS_WARNING);
  }

  return;
}

/*============================= end of c_check_inputs() =================*/

/*============================= c_dref() ================================*/

/*
  Flux albedo for given angle of incidence, given a bidirectional reflectivity.

  INPUTS
    wvnmlo    :  Lower wavenumber (inv-cm) of spectral interval
    wvnmhi    :  Upper wavenumber (inv-cm) of spectral interval
    mu        :  Cosine of incidence angle
    brdf_type :  BRDF type
    brdf      :  pointer to disort_brdf structure
    callnum   :  number of surface calls

  INTERNAL VARIABLES

    gmu    : The NMUG angle cosine quadrature points on (0,1)
             NMUG is set in cdisort.h
    gwt    : The NMUG angle cosine quadrature weights on (0,1)

   Called by- c_check_inputs
   Calls- c_gaussian_quadrature, c_errmsg, c_bidir_reflectivity
 --------------------------------------------------------------------*/

double c_dref(double       wvnmlo,
              double       wvnmhi,
              double       mu,
	      int          brdf_type,
	      disort_brdf *brdf,
	      int          callnum )
{
  static int
    pass1 = TRUE;
  register int
    jg,k;
  double
    ans,sum;
  static double
    gmu[NMUG],gwt[NMUG];

  if (pass1) {
    pass1 = FALSE;
    c_gaussian_quadrature(NMUG/2,gmu,gwt);
    for (k = 1; k <= NMUG/2; k++) {
      GMU(k+NMUG/2) = -GMU(k);
      GWT(k+NMUG/2) =  GWT(k);
    }
  }

  if (fabs(mu) > 1.) {
    c_errmsg("dref--input argument error(s)",DS_ERROR);
  }

  ans = 0.;
  /*
   * Loop over azimuth angle difference
   */
  for (jg = 1; jg <= NMUG; jg++) {
    /*
     * Loop over angle of reflection
     */
    sum = 0.;
    for (k = 1; k <= NMUG/2; k++) {
      sum += GWT(k) * GMU(k) *
	c_bidir_reflectivity ( wvnmlo, wvnmhi, GMU(k), mu, M_PI*GMU(jg), brdf_type, brdf, callnum );
    }
    ans += GWT(jg)*sum;
  }
  if (ans < 0. || ans > 1.) {
    c_errmsg("DREF--albedo value not in [0,1]",DS_WARNING);
  }

  return ans;
}

/*============================= end of c_dref() =========================*/

/*============================= c_legendre_poly() =======================*/

/*
       Computes the normalized associated Legendre polynomial, defined
       in terms of the associated Legendre polynomial Plm = P-sub-l-super-m as

          Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)

       for fixed order m and all degrees from l = m to TWONM1.
       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
       from a prior call to the routine.

       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of High-Order
                    Associated Legendre Polynomials, J. Quant. Spectrosc. Radiat. Transfer 10,
                    557-562, 1970. (hereafter D/A)

       METHOD: Varying degree recurrence relationship.

       NOTES:
       (1) The D/A formulas are transformed by setting m=n-1; l=k-1.
       (2) Assumes that routine is called first with  m = 0, then with
           m = 1, etc. up to  m = twonm1.


  I N P U T     V A R I A B L E S:

       nmu    :  Number of arguments of YLM
       m      :  Order of YLM
       maxmu  : 
       twonm1 :  Max degree of YLM
       MU(i)  :  Arguments of YLM (i = 1 to nmu)

       If m > 0, YLM(m-1,i) for i = 1 to nmu is assumed to exist from a prior call.


  O U T P U T     V A R I A B L E:

       YLM(l,i) :  l = m to twonm1, normalized associated Legendre polynomials
                   evaluated at argument MU(i)

   Called by- c_disort, c_albtrans
 -------------------------------------------------------------------*/

void c_legendre_poly(int     nmu,
                     int     m,
                     int     maxmu,
                     int     twonm1,
                     double *mu,
                     double *ylm)
{
  register int
    i,l;
  register double
    tmp1,tmp2;

  if (m == 0) {
    /*
     * Upward recurrence for ordinary Legendre polynomials
     */
    for (i = 1; i <= nmu; i++) {
      YLM(0,i) = 1.;
      YLM(1,i) = MU(i);
    }
    for (l = 2; l <= twonm1; l++) {
      for (i = 1; i <= nmu; i++) {
        YLM(l,i) = ((double)(2*l-1)*MU(i)*YLM(l-1,i)-(double)(l-1)*YLM(l-2,i))/l;
      }
    }
  }
  else {
    for (i = 1; i <= nmu; i++) {
      /*
       * Y-sub-m-super-m; derived from D/A eqs. (11,12), STWL(58c)
       */
      YLM(m,i) = -sqrt((1.-1./(2*m))*(1.-SQR(MU(i))))*YLM(m-1,i);

      /*
       * Y-sub-(m+1)-super-m; derived from D/A eqs.(13,14) using eqs.(11,12), STWL(58f)
       */
      YLM(m+1,i) = sqrt(2.*m+1.)*MU(i)*YLM(m,i);
    }
    /*
     * Upward recurrence; D/A eq.(10), STWL(58a)
     */
    for (l = m+2; l <= twonm1; l++) {
      tmp1 = sqrt((l-m  )*(l+m  ));
      tmp2 = sqrt((l-m-1)*(l+m-1));
      for (i = 1; i <= nmu; i++) {
        YLM(l,i) = ((double)(2*l-1)*MU(i)*YLM(l-1,i)-tmp2*YLM(l-2,i))/tmp1;
      }
    }
  }

  return;
}

/*============================= end of c_legendre_poly() ================*/

/*============================= c_planck_func1() ========================*/

/*
  Computes Planck function integrated between two wavenumbers

  INPUT :  wnumlo : Lower wavenumber (inv cm) of spectral interval
           wnumhi : Upper wavenumber
           t      : Temperature (K)

  OUTPUT : ans    : Integrated Planck function ( Watts/sq m ) = Integral (WNUMLO to WNUMHI) of
                    2h c*c nu*nu*nu / ( exp(hc nu/kT) - 1) (where h=Plancks constant, c=speed of
                    light, nu=wavenumber, T=temperature, and k = Boltzmann constant)

  Reference : Specifications of the Physical World: New Value of the Fundamental Constants,
                Dimensions/N.B.S., Jan. 1974

  Method :  For wnumlo close to wnumhi, a Simpson-rule quadrature is done to avoid
            ill-conditioning; otherwise

            (1)  For wnumlo or wnumhi small, integral(0 to WNUMLO/HI) is calculated by expanding
                 the integrand in a power series and integrating term by term;

            (2)  Otherwise, integral (wnumlo/hi to infinity) is calculated by expanding the
                 denominator of the integrand in powers of the exponential and integrating
                 term by term.

  Accuracy :  At least 6 significant digits, assuming the physical constants are infinitely accurate

  ERRORS WHICH ARE NOT TRAPPED:

      * Power or exponential series may underflow, giving no significant digits.  
        This may or may not be of concern, depending on the application.

      * Simpson-rule special case is skipped when denominator of integrand will cause overflow.
        In that case the normal procedure is used, which may be inaccurate if the wavenumber limits
        (wnumlo, wnumhi) are close together.

  LOCAL VARIABLES

        a1,2,... :  Power series coefficients
        c2       :  h * c / k, in units cm*K (h = Plancks constant,
                      c = speed of light, k = Boltzmann constant)
        D(I)     :  Exponential series expansion of integral of
                       Planck function from wnumlo (i=1) or wnumhi
                       (i=2) to infinity
        ex       :  exp( - V(I) )
        exm      :  ex**m
        mmax     :  No. of terms to take in exponential series
        mv       :  Multiples of V(I)
        P(I)     :  Power series expansion of integral of
                       Planck function from zero to WNUMLO (I=1) or
                       WNUMHI (I=2)
        sigma    :  Stefan-Boltzmann constant (W/m**2/K**4)
        sigdpi   :  sigma/pi
        smallv   :  Number of times the power series is used (0,1,2)
        V(i)     :  c2 * (wnumlo(i=1) or wnumhi(i=2))/temperature
        vcut     :  Power-series cutoff point
        vcp      :  Exponential series cutoff points
        vmax     :  Largest allowable argument of EXP function

   Called by- c_disort
   Calls- c_errmsg
 ----------------------------------------------------------------------*/

#define PLKF(x) ({const double _x = (x); _x*_x*_x/(exp(_x)-1.);})
#define A1    (1./3.)
#define A2    (-1./8.)
#define A3    (1./60.)
#define A4    (-1./5040.)
#define A5    (1./272160.)
#define A6    (-1./13305600.)
#define C2    (1.438786)
#define SIGMA (5.67032E-8)
#define VCUT  (1.5)

double c_planck_func1(double wnumlo,
                      double wnumhi,
                      double t)
{
  register int
    i,k,m,mmax,n,smallv;
  int
    converged;
  static int
    initialized = FALSE;
  const double
    vcp[7] = {10.25,5.7,3.9,2.9,2.3,1.9,0.0};
  double
    del,ex,exm,hh,mv,oldval,
    val,val0,vsq,d[2],p[2],v[2],
    ans;
  static double
    vmax,sigdpi,conc;

  if (!initialized) {
    vmax   = log(DBL_MAX);
    sigdpi = SIGMA/M_PI;
    conc   = 15./pow(M_PI,4.);

    initialized = TRUE;
  }

  if (t < 0. || wnumhi <= wnumlo || wnumlo < 0.) {
    c_errmsg("planck_func1--temperature or wavenums. wrong",DS_ERROR);
  }

  if (t < 1.e-4) {
    return 0.;
  }

  v[0] = C2*wnumlo/t;
  v[1] = C2*wnumhi/t;

  if (v[0] > DBL_EPSILON && v[1] < vmax && (wnumhi-wnumlo)/wnumhi < 1.e-2) {
    /*
     * Wavenumbers are very close.  Get integral by iterating Simpson rule to convergence.
     */
    hh     = v[1]-v[0];
    oldval = 0.;
    val0   = PLKF(v[0])+PLKF(v[1]);
    converged = FALSE;
    for (n = 1; n <= 10; n++) {
      del = hh/(2*n);
      val = val0;
      for (k = 1; k <= 2*n-1; k++) {
        val += (double)(2*(1+k%2))*PLKF(v[0]+(double)k*del);
      }
      val *= del*A1;
      if (fabs((val-oldval)/val) <= 1.e-6) {
        /* convergence */
        converged = TRUE;
        break;
      }
      oldval = val;
    }
    if (!converged) {
      c_errmsg("planck_func1--Simpson rule didn't converge",DS_WARNING);
    }

    return sigdpi*pow(t,4.0)*conc*val;
  }

  /*
   * General case
   */
  smallv = 0;
  for (i = 1; i <= 2; i++) {
    if(v[i-1] < VCUT) {
      /*
       * Use power series
       */
      smallv++;
      vsq    = SQR(v[i-1]);
      p[i-1] = conc*vsq*v[i-1]*(A1+v[i-1]*(A2+v[i-1]*(A3+vsq*(A4+vsq*(A5+vsq*A6)))));
    }
    else {
      /*
       * Use exponential series
       *
       * Find upper limit of series
       */
      mmax = 1;
      while (v[i-1] < vcp[mmax-1]) {
        mmax++;
      }

      ex     = exp(-v[i-1]);
      exm    = 1.;
      d[i-1] = 0.;
      for (m = 1; m <= mmax; m++) {
        mv      = (double)m*v[i-1];
        exm     = ex*exm;
        d[i-1] += exm*(6.+mv*(6.+mv*(3.+mv)))/SQR(m*m);
      }
      d[i-1] *= conc;
    }
  }
  /*
   * Handle ill-conditioning
   */
  if (smallv == 2) {
    /*
     * wnumlo and wnumhi both small
     */
    ans = p[1]-p[0];
  }
  else if (smallv == 1) {
    /*
     * wnumlo small, wnumhi large
     */
    ans = 1.-p[0]-d[1];
  }
  else {
    /*
     * wnumlo and wnumhi both large
     */
    ans = d[0]-d[1];
  }
  ans *= sigdpi*pow(t,4.0);
  if (ans == 0.) {
    c_errmsg("planck_func1--returns zero; possible underflow",DS_WARNING);
  }

  return ans;
}

#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6
#undef C2
#undef SIGMA
#undef VCUT
#undef PLKF

/*============================= end of c_planck_func1() =================*/

/*============================= c_print_avg_intensities() ===============*/

/*
   Print azimuthally averaged intensities at user angles

   Called by- c_disort
 -------------------------------------------------------------------*/

void c_print_avg_intensities(disort_state *ds,
			     disort_output *out)
{
  register int
    iu,iumax,iumin,
    lenfmt,lu,np,npass;

  if(ds->numu < 1) {
    return;
  }

  fprintf(stdout,"\n\n *******  AZIMUTHALLY AVERAGED INTENSITIES (at user polar angles)  ********\n");
  lenfmt = 8;
  npass  = 1+(ds->numu-1)/lenfmt;
  fprintf(stdout,"\n   Optical   Polar Angle Cosines"
                 "\n     Depth");

  for (np = 1; np <= npass; np++) {
    iumin = 1+lenfmt*(np-1);
    iumax = IMIN(lenfmt*np,ds->numu);
    fprintf(stdout,"\n          "); 
    for (iu = iumin; iu <= iumax; iu++) {
      fprintf(stdout,"%14.5f",UMU(iu));
    }
    fprintf(stdout,"\n");

    for (lu = 1; lu <= ds->ntau; lu++) {
      fprintf(stdout,"%10.4f",UTAU(lu));
      for (iu = iumin; iu <= iumax; iu++) {
        fprintf(stdout,"%14.4e",U0U(iu,lu));
      }
      fprintf(stdout,"\n");
    }
  }

  return;
}

/*============================= end of c_print_avg_intensities() ========*/

/*============================= c_print_inputs() ========================*/

/*
   Print values of input variables

   Called by- c_disort
 --------------------------------------------------------------------*/

void c_print_inputs(disort_state *ds,
                    double       *dtaucpr,
                    int           scat_yes,
                    int           deltam,
                    int           corint,
                    double       *flyr,
                    int           lyrcut,
                    double       *oprim,
                    double       *tauc,
                    double       *taucpr)
{
  register int
    iq,iu,j,k,lc,lu;

  fprintf(stdout,"\n\n"
                 " ****************************************************************************************************\n"
                 " DISORT: %s\n"
                 " ****************************************************************************************************\n",
                 ds->header);

  fprintf(stdout,"\n No. streams =%4d     No. computational layers =%4d\n",ds->nstr,ds->nlyr);

  if (ds->flag.ibcnd != SPECIAL_BC) {
    fprintf(stdout,"%4d User optical depths :",ds->ntau);
    for (lu = 1; lu <= ds->ntau; lu++) {
      fprintf(stdout,"%10.4f",UTAU(lu));
      if (lu%10 == 0) {
        fprintf(stdout,"\n                          ");
      }
    }
    fprintf(stdout,"\n");
  }

  if (!ds->flag.onlyfl) {
    fprintf(stdout,"%4d User polar angle cosines :",ds->numu);
    for (iu = 1; iu <= ds->numu; iu++) {
      fprintf(stdout,"%9.5f",UMU(iu));
      if (iu%10 == 0) {
        fprintf(stdout,"\n                               ");
      }
    }
    fprintf(stdout,"\n");
  }

  if (!ds->flag.onlyfl && ds->flag.ibcnd != SPECIAL_BC) {
    fprintf(stdout,"%4d User azimuthal angles :",ds->nphi);
    for (j = 1; j <= ds->nphi; j++) {
      fprintf(stdout,"%9.2f",PHI(j));
      if (j%10 == 0) {
        fprintf(stdout,"n                            ");
      }
    }
    fprintf(stdout,"\n");
  }

  if (!ds->flag.planck || ds->flag.ibcnd == SPECIAL_BC) {
    fprintf(stdout," No thermal emission\n");
  }

  if (ds->flag.spher == TRUE) {
    fprintf(stdout," Pseudo-spherical geometry invoked\n");
  }

  if (ds->flag.general_source == TRUE) {
    fprintf(stdout," Calculation with general source term\n");
  }

  if (ds->flag.ibcnd == GENERAL_BC) {
    fprintf(stdout," Boundary condition flag: ds.flag.ibcnd = GENERAL_BC\n");
    fprintf(stdout,"    Incident beam with intensity =%11.3e and polar angle cosine = %8.5f  and azimuth angle =%7.2f\n",
                   ds->bc.fbeam,ds->bc.umu0,ds->bc.phi0);
    fprintf(stdout,"    plus isotropic incident intensity =%11.3e\n",ds->bc.fisot);

    if (ds->bc.fluor > 0.0 ) {
      fprintf(stdout,"    Bottom isotropic exiting intensity =%11.3e\n",ds->bc.fluor);
    }
    if (ds->flag.lamber) {
      fprintf(stdout,"    Bottom albedo (Lambertian) =%8.4f\n",ds->bc.albedo);
    }
    else {
      fprintf(stdout,"    Bidirectional reflectivity at bottom\n");
    }

    if(ds->flag.planck) {
      fprintf(stdout,"    Thermal emission in wavenumber interval :%14.4f%14.4f\n",ds->wvnmlo,ds->wvnmhi);
      fprintf(stdout,"    Bottom temperature =%10.2f    Top temperature =%10.2f    Top emissivity =%8.4f\n",
                     ds->bc.btemp,ds->bc.ttemp,ds->bc.temis);
    }
  }
  else if (ds->flag.ibcnd == SPECIAL_BC) {
    fprintf(stdout," Boundary condition flag: ds.flag.ibcnd = SPECIAL_BC\n");
    fprintf(stdout,"    Isotropic illumination from top and bottom\n");
    fprintf(stdout,"    Bottom albedo (Lambertian) =%8.4f\n",ds->bc.albedo);
  }
  else {
    c_errmsg("Unrecognized ds.flag.ibcnd",DS_WARNING);
  }

  if (deltam) {
    fprintf(stdout," Uses delta-M method\n");
  }
  else {
    fprintf(stdout," Does not use delta-M method\n");
  }

  if (corint) {
    fprintf(stdout," Uses TMS/IMS method\n");
  }
  else {
    fprintf(stdout," Does not use TMS/IMS method\n");
  }

  if (ds->flag.ibcnd == SPECIAL_BC) {
    fprintf(stdout," Calculate albedo and transmissivity of medium vs. incident beam angle\n");
  }
  else if (ds->flag.onlyfl) {
    fprintf(stdout," Calculate fluxes only\n");
  }
  else {
    fprintf(stdout," Calculate fluxes and intensities\n");
  }

  fprintf(stdout," Relative convergence criterion for azimuth series =%11.2e\n",ds->accur);

  if (lyrcut) {
    fprintf(stdout," Sets radiation = 0 below absorption optical depth 10\n");
  }

  /*
   * Print layer variables (to read, skip every other line)
   */
  if(ds->flag.planck) {
    fprintf(stdout,"\n                                     <------------- Delta-M --------------->");
    fprintf(stdout,"\n                   Total    Single                           Total    Single");
    fprintf(stdout,"\n       Optical   Optical   Scatter   Separated   Optical   Optical   Scatter    Asymm");
    fprintf(stdout,"\n         Depth     Depth    Albedo    Fraction     Depth     Depth    Albedo   Factor   Temperature\n");
  }
  else {
    fprintf(stdout,"\n                                     <------------- Delta-M --------------->");
    fprintf(stdout,"\n                   Total    Single                           Total    Single");
    fprintf(stdout,"\n       Optical   Optical   Scatter   Separated   Optical   Optical   Scatter    Asymm");
    fprintf(stdout,"\n         Depth     Depth    Albedo    Fraction     Depth     Depth    Albedo   Factor\n");
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    
    if (ds->flag.planck) {
      fprintf(stdout,"%4d%10.4f%10.4f%10.5f%12.5f%10.4f%10.4f%10.5f%9.4f%14.3f\n",
                     lc,DTAUC(lc),TAUC(lc),SSALB(lc),FLYR(lc),DTAUCPR(lc),TAUCPR(lc),OPRIM(lc),PMOM(1,lc),TEMPER(lc-1));
    }
    else {
      fprintf(stdout,"%4d%10.4f%10.4f%10.5f%12.5f%10.4f%10.4f%10.5f%9.4f\n",
                     lc,DTAUC(lc),TAUC(lc),SSALB(lc),FLYR(lc),DTAUCPR(lc),TAUCPR(lc),OPRIM(lc),PMOM(1,lc));
    }
  }
  if (ds->flag.planck) {
    fprintf(stdout,"                                                                                     %14.3f\n",
            TEMPER(ds->nlyr));
  }

  if (ds->flag.prnt[4] && scat_yes) {
    fprintf(stdout,"\n Number of Phase Function Moments = %5d\n",ds->nmom+1);
    fprintf(stdout," Layer   Phase Function Moments\n");
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (SSALB(lc) > 0.) {
        fprintf(stdout,"%6d",lc);
        for (k = 0; k <= ds->nmom; k++) {
          fprintf(stdout,"%11.6f",PMOM(k,lc));
          if ((k+1)%10 == 0) {
            fprintf(stdout,"\n      ");
          } 
        }
        fprintf(stdout,"\n");
      }
    }
  }

  if (ds->flag.general_source == TRUE) {
    fprintf(stdout," Calculation with general source term\n");
    j = 0;
    for (lc = 1; lc <= ds->nlyr; lc++) {
      fprintf(stdout,"%4d%10.4f",lc,DTAUC(lc));
      for (iq = 1; iq <= ds->nstr; iq++) {	
	fprintf(stdout,"%13.6e",GENSRC(j,lc,iq));
      }
      fprintf(stdout,"\n");
    }
  }


  return;
}

/*============================= end of c_print_inputs() =================*/

/*============================= c_print_intensities() ===================*/

/*
   Prints the intensity at user polar and azimuthal angles
   All arguments are disort state or output variables

   Called by- c_disort
 -------------------------------------------------------------------*/

void c_print_intensities(disort_state  *ds,
                         disort_output *out)
{
  register int
    iu,j,jmax,jmin,lenfmt,lu,np,npass;

  if (ds->nphi < 1) {
    return;
  }

  fprintf(stdout,"\n\n *********  I N T E N S I T I E S  *********\n");
  lenfmt = 10;
  npass  = 1+(ds->nphi-1)/lenfmt;
  fprintf(stdout,"\n             Polar   Azimuth angles (degrees)");
  fprintf(stdout,"\n   Optical   Angle");
  fprintf(stdout,"\n    Depth   Cosine\n");
  for (lu = 1; lu <= ds->ntau; lu++) {
    for (np = 1; np <= npass; np++) {
      jmin = 1+lenfmt*(np-1);
      jmax = IMIN(lenfmt*np,ds->nphi);
      fprintf(stdout,"\n                  ");
      for (j = jmin; j <= jmax; j++) {
        fprintf(stdout,"%11.2f",PHI(j));
      }
      fprintf(stdout,"\n");
      if (np == 1) {
        fprintf(stdout,"%10.4f%8.4f",UTAU(lu),UMU(1));
        for (j = jmin; j <= jmax; j++) {
          fprintf(stdout,"%11.3e",UU(1,lu,j));
        }
        fprintf(stdout,"\n");
      }
      else {
        fprintf(stdout,"          %8.4f",UMU(1));
        for (j = jmin; j <= jmax; j++) {
          fprintf(stdout,"%11.3e",UU(1,lu,j));
        }
        fprintf(stdout,"\n");
      }
      for (iu = 2; iu <= ds->numu; iu++) {
        fprintf(stdout,"          %8.4f",UMU(iu));
        for (j = jmin; j <= jmax; j++) {
          fprintf(stdout,"%11.3e",UU(iu,lu,j));
        }
        fprintf(stdout,"\n");
      }
    }
  }

  return;
}

/*============================= end of c_print_intensities() ============*/

/*============================= c_gaussian_quadrature() =================*/

/*
   Compute weights and abscissae for ordinary Gaussian quadrature
   on the interval (0,1);  that is, such that
       sum(i=1 to M) ( GWT(i) f(GMU(i)) )
   is a good approximation to integral(0 to 1) ( f(x) dx )

   INPUT :     m        order of quadrature rule

   OUTPUT :    GMU(I)   array of abscissae (I = 1 TO M)
               GWT(I)   array of weights (I = 1 TO M)

   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
                 Integration, Academic Press, New York, pp. 87, 1975

   METHOD:     Compute the abscissae as roots of the Legendre polynomial P-sub-M using a cubically convergent
               refinement of Newton's method.  Compute the weights from eq. 2.7.3.8 of Davis/Rabinowitz.  Note
               that Newton's method can very easily diverge; only a very good initial guess can guarantee convergence.
               The initial guess used here has never led to divergence even for M up to 1000.

   ACCURACY:   Relative error no better than TOL or computer precision (DBL_EPSILON), whichever is larger

   INTERNAL VARIABLES:
    iter      : Number of Newton Method iterations
    pm2,pm1,p : 3 successive Legendre polynomials
    ppr       : Derivative of Legendre polynomial
    p2pri     : 2nd derivative of Legendre polynomial
    tol       : Convergence criterion for Legendre poly root iteration
    x,xi      : Successive iterates in cubically-convergent version of Newtons Method (seeking roots of Legendre poly.)

   Called by- c_dref, c_disort_set, c_surface_bidir
   Calls- c_errmsg
 -------------------------------------------------------------------*/

/* Maximum allowed iterations of Newton Method */
#define MAXIT 1000

void c_gaussian_quadrature(int    m,
                           double *gmu,
                           double *gwt)
{
  static int
    initialized = FALSE;
  register int
    iter,k,lim,nn,np1;
  double
    cona,t,en,nnp1,p=0,p2pri,pm1,pm2,ppr,
    prod,tmp,x,xi;
  static double
    tol;

  if (!initialized) {
    tol         = 10.*DBL_EPSILON;
    initialized = TRUE;
  }

  if (m < 1) {
    c_errmsg("gaussian_quadrature--Bad value of m",DS_ERROR);
  }

  if (m == 1) {
    GMU(1) = 0.5;
    GWT(1) = 1.0;
    return;
  }

  en   = (double)m;
  np1  = m+1;
  nnp1 = m*np1;
  cona = (double)(m-1)/(8*m*m*m);
  lim  = m/2;
  for (k = 1; k <= lim; k++) {
    /*
     * Initial guess for k-th root of Legendre polynomial, from Davis/Rabinowitz (2.7.3.3a)
     */
    t = (double)(4*k-1)*M_PI/(4*m+2);
    x = cos(t+cona/tan(t));

    /*
     * Upward recurrence for Legendre polynomials
     */
    for (iter = 1; iter <= MAXIT+1; iter++) {
      if (iter > MAXIT) {
        c_errmsg("gaussian_quadrature--max iteration count",DS_ERROR);
      }
      pm2 = 1.;
      pm1 = x;
      for (nn = 2; nn <= m; nn++) {
        p   = ((double)(2*nn-1)*x*pm1-(double)(nn-1)*pm2)/nn;
        pm2 = pm1;
        pm1 = p;
      }
      /*
       * Newton Method
       */
      tmp   = 1./(1.-x*x);
      ppr   = en*(pm2-x*p)*tmp;
      p2pri = (2.*x*ppr-nnp1*p)*tmp;
      xi    = x-p/ppr*(1.+p/ppr*p2pri/(2.*ppr));
      /*
       * Check for convergence
       */
      if (fabs(xi-x) <= tol) {
        break;
      }
      else {
        x = xi;
      }
    }

    /*
     * Iteration finished--calculate weights, abscissae for (-1,1)
     */
    GMU(k)     = -x;
    GWT(k)     = 2./(tmp*SQR(en*pm2));
    GMU(np1-k) = -GMU(k);
    GWT(np1-k) =  GWT(k);
  }

  /*
   * Set middle abscissa and weight for rules of odd order
   */
  if (m%2 != 0) {
    GMU(lim+1) = 0.;
    prod       = 1.;
    for (k = 3; k <= m; k+=2) {
      prod *= (double)k/(k-1);
    }
    GWT(lim+1) = 2./SQR(prod);
  }
  /*
   * Convert from (-1,1) to (0,1)
   */
  for (k = 1; k <= m; k++) {
    GMU(k) = 0.5*GMU(k)+0.5;
    GWT(k) = 0.5*GWT(k);
  }

  return;
}

#undef MAXIT

/*============================= end of c_gaussian_quadrature() ==========*/

/*============================= c_ratio() ===============================*/

/*
 * Calculate ratio a/b with overflow and underflow protection
 * (thanks to Prof. Jeff Dozier for some suggestions here).
 * 
 * Modification in this C version: in the case b == 0., returns 1.+a.
 *
 * Called by: c_disort
 */

double c_ratio(double a,
             double b)
{
  static int
    initialized = FALSE;
  static double
    tiny,huge,powmax,powmin;
  double
    ans,absa,absb,powa,powb;

  if(!initialized) {
    tiny   = DBL_MIN;
    huge   = DBL_MAX;
    powmax = log10(huge);
    powmin = log10(tiny);
   
    initialized = TRUE;
  }

  if (c_fcmp(b,0.) == 0) {
    ans = 1.+a;
  }
  else if (c_fcmp(a,0.) == 0) {
    ans = 0.;
  }
  else {
    absa = fabs(a);
    absb = fabs(b);
    powa = log10(absa);
    powb = log10(absb);
    if (c_fcmp(absa,tiny) < 0 && c_fcmp(absb,tiny) < 0) {
      ans = 1.;
    }
    else if (c_fcmp(powa-powb,powmax) >= 0) {
      ans = huge;
    }
    else if(c_fcmp(powa-powb,powmin) <= 0) {
      ans = tiny;
    }
    else {
      ans = absa/absb;
    }

   /*
    * NOTE: Don't use old trick of determining sign from a*b because a*b
    *       may overflow or underflow.
    */
    if ( (a > 0. && b < 0.) || (a < 0. && b > 0.) ) {
      ans *= -1;
    }
  }

  return ans;
}

/*============================= end of c_ratio() ========================*/

/*============================= c_fcmp() ================================*/

/*
 * Derived from fcmp(), version 1.2.2, 
 * Copyright (c) 1998-2000 Theodore C. Belding
 * University of Michigan Center for the Study of Complex Systems
 * <mailto:Ted.Belding@umich.edu>
 * <http://fcmp.sourceforge.net>
 *
 * The major modification we have made is to remove the "epsilon" argument
 * and set epsilon inside the fcmp() function.
 *
 * Description:
 *   It is generally not wise to compare two floating-point values for
 *   exact equality, for example using the C == operator.  The function
 *   fcmp() implements Knuth's suggestions for safer floating-point
 *   comparison operators, from:
 *   Knuth, D. E. (1998). The Art of Computer Programming.
 *   Volume 2: Seminumerical Algorithms. 3rd ed. Addison-Wesley.
 *   Section 4.2.2, p. 233. ISBN 0-201-89684-2.
 *
 * Input parameters:
 *   x1, x2: numbers to be compared
 *
 * Returns:
 *   -1 if x1 < x2
 *    0 if x1 == x2
 *    1 if x1 > x2		
 */

int c_fcmp(double x1,
           double x2) {
  int 
    exponent;
  double
    delta,
    difference;
  const double
    epsilon = DBL_EPSILON;
  
  /* 
   * Get exponent(max(fabs(x1),fabs(x2))) and store it in exponent. 
   *
   * If neither x1 nor x2 is 0,
   * this is equivalent to max(exponent(x1),exponent(x2)).
   *
   * If either x1 or x2 is 0, its exponent returned by frexp would be 0,
   * which is much larger than the exponents of numbers close to 0 in
   * magnitude. But the exponent of 0 should be less than any number
   * whose magnitude is greater than 0.
   *
   * So we only want to set exponent to 0 if both x1 and x2 are 0. 
   * Hence, the following works for all x1 and x2. 
   */
  frexp(fabs(x1) > fabs(x2) ? x1 : x2,&exponent);

  /* 
   * Do the comparison.
   *
   * delta = epsilon*pow(2,exponent)
   *
   * Form a neighborhood around x2 of size delta in either direction.
   * If x1 is within this delta neighborhood of x2, x1 == x2.
   * Otherwise x1 > x2 or x1 < x2, depending on which side of
   * the neighborhood x1 is on.
   */
  delta      = ldexp(epsilon,exponent); 
  difference = x1-x2;

  if (difference > delta) {
    /* x1 > x2 */
    return 1;
  }
  else if (difference < -delta) {
    /* x1 < x2 */
    return -1;
  }
  else  {
    /* -delta <= difference <= delta */
    return 0;  /* x1 == x2 */
  }
}

/*============================= end of c_fcmp() =========================*/

/*============================= c_self_test() ===========================*/

/*
 * If  compare is FALSE, set up self-test disort_state ds_test.
 * If  compare is TRUE, compare self-test results with correct
 * answers, abort if error, and free self-test memory.
 *
 * (See file 'DISORT.txt' for variable definitions.)
 *
 *    I N T E R N A L    V A R I A B L E S:
 *
 *       acc        Relative accuracy required for passing self-test
 *       error      Relative errors in DISORT output variables
 *       ok         Logical variable for determining failure of self-test
 *
 * Called by- c_disort
 * Calls- c_errmsg
 */

void c_self_test(int            compare,
                 int           *prntu0,
                 disort_state  *ds,
                 disort_output *out)
{
  const double
    acc = 1.e-4;
  int
    i,ok;    
  double
    error;

  if(compare == FALSE) {
    for (i = 0; i < 5; i++) {
      ds->flag.prnt[i] = FALSE;
    }
    ds->flag.ibcnd     = GENERAL_BC;
    ds->flag.usrang    = TRUE;
    ds->flag.usrtau    = TRUE;
    ds->flag.lamber    = TRUE;
    ds->flag.onlyfl    = FALSE;
    ds->flag.planck    = TRUE;
    ds->flag.quiet     = QUIET;
    ds->flag.spher     = FALSE;
    ds->flag.general_source = FALSE;
    ds->flag.brdf_type = BRDF_NONE;
    ds->flag.intensity_correction     = TRUE;
    ds->flag.old_intensity_correction = TRUE;
    ds->flag.output_uum=FALSE;

    ds->nstr = 4;
    ds->nlyr = 1;
    ds->nmom = 4;
    ds->numu = 1;
    ds->ntau = 1;
    ds->nphi = 1;

    /* Allocate memory for self test */
    c_disort_state_alloc(ds);
    c_disort_out_alloc(ds,out);

    ds->accur  = 1.e-4;
    ds->wvnmlo =     0.;
    ds->wvnmhi = 50000.;

    ds->bc.fbeam  =  M_PI;
    ds->bc.umu0   =   .866;
    ds->bc.phi0   =   0.;
    ds->bc.fisot  =   1.;
    ds->bc.fluor  =   0.;
    ds->bc.albedo =    .7;
    ds->bc.ttemp  = 100.;
    ds->bc.btemp  = 300.;
    ds->bc.temis  =    .8;

    TEMPER(0) = 210.;
    TEMPER(1) = 200.;

    DTAUC(1)  = 1.;
    SSALB(1)  =  .9;

    /* Haze L moments */
    PMOM(0,1) = 1.;
    PMOM(1,1) =  .8042;
    PMOM(2,1) =  .646094;
    PMOM(3,1) =  .481851;
    PMOM(4,1) =  .359056;

    UMU(1)  =  0.5;
    UTAU(1) =  0.5;
    PHI(1)  = 90.0;

    return;
  }
  else if (compare == TRUE) {
    /*
     * Compare test case results with correct answers and abort if bad
     */
    ok = TRUE;

    error = (out->uu[0]-47.865571)/47.865571;
    if (fabs(error) > acc) {
      ok = FALSE;
      fprintf(stderr,"Output variable uu differed by %g percent from correct value.\n",100.*error);
    }

    error = (out->rad[0].rfldir-1.527286)/1.527286;
    if (fabs(error) > acc) {
      ok = FALSE;
      fprintf(stderr,"Output variable rfldir differed by %g percent from correct value.\n",100.*error);
    }

    error = (out->rad[0].rfldn-28.372225)/28.372225;
    if (fabs(error) > acc) {
      ok = FALSE;
      fprintf(stderr,"Output variable rfldn differed by %g percent from correct value.\n",100.*error);
    }

    error = (out->rad[0].flup-152.585284)/152.585284;
    if (fabs(error) > acc) {
      ok = FALSE;
      fprintf(stderr,"Output variable flup differed by %g percent from correct value.\n",100.*error);
    }

    /* Free allocated memory for self test */
    c_disort_out_free(ds,out);
    c_disort_state_free(ds);

    if (!ok) {
      c_errmsg("DISORT--self-test failed",DS_ERROR);
    }

    return;
  }
  else {
    fprintf(stderr,"**error--self_test(): compare=%d not recognized\n",compare);
    exit(1);
  }
}

/*============================= end of c_self_test() =====================*/

/******************************************************************
 ********** end of DISORT service routines ************************
 ******************************************************************

 ******************************************************************
 ********** ds->flag.ibcnd = SPECIAL_BC routines ******************
 ******************************************************************/

/*============================= c_albtrans() =============================*/

/*
   DISORT special case to get only albedo and transmissivity of entire medium as a function of incident beam angle
   (many simplifications because boundary condition is just isotropic illumination, there are no thermal sources, and
   particular solutions do not need to be computed).  See Ref. S2 and references therein for details.
   The basic idea is as follows.  The reciprocity principle leads to the following relationships for a plane-parallel,
   vertically inhomogeneous medium lacking thermal (or other internal) sources:
  
      albedo(theta) = u_0(theta) for unit-intensity isotropic
                       illumination at *top* boundary
       trans(theta) =  u_0(theta) for unit-intensity isotropic
                       illumination at *bottom* boundary
    where

       albedo(theta) = albedo for beam incidence at angle theta
       trans(theta)  = transmissivity for beam incidence at angle theta
       u_0(theta)    = upward azim-avg intensity at top boundary
                       at angle theta

   O U T P U T    V A R I A B L E S:
  
       ALBMED(iu)   Albedo of the medium as a function of incident
                    beam angle cosine UMU(IU)
  
       TRNMED(iu)   Transmissivity of the medium as a function of
                    incident beam angle cosine UMU(IU)

    I N T E R N A L   V A R I A B L E S:

       ncd         number of diagonals below/above main diagonal
       rcond       estimate of the reciprocal condition of matrix CBAND; for system  CBAND*X = B, relative
                   perturbations in CBAND and B of size epsilon may cause relative perturbations in X of size
                   epsilon/RCOND.  If RCOND is so small that
                          1.0 + RCOND .eq. 1.0
                   is true, then CBAND may be singular to working precision.
       cband       Left-hand side matrix of linear system eq. SC(5), scaled by eq. SC(12);
                   in banded form required by LINPACK solution routines
       ncol        number of columns in CBAND matrix
       ipvt        INTEGER vector of pivot indices (most others documented in DISORT)
  
   Called by- c_disort
   Calls- c_legendre_poly, c_sgbco, c_solve_eigen, c_interp_eigenvec, c_set_matrix, c_solve1,
          c_albtrans_intensity, c_albtrans_spherical, c_print_albtrans
 --------------------------------------------------------------------------------------*/

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
                double        *wk)
{
  int
    lyrcut,ncol;
  register int
    iq,iu,l,lc,mazim,ncd,ncut;
  double
    delm0,rcond,sgn,sphalb,sphtrn;

  mazim = 0;
  delm0 = 1.;
  /*
   * Set DISORT variables that are ignored in this special case but are needed below in argument
   * lists of subroutines shared with general case
   */
  ncut            = ds->nlyr;
  lyrcut          = FALSE;
  ds->bc.fisot    = 1.;
  ds->bc.fluor    = 0.;
  ds->flag.lamber = TRUE;

  /*
   * Get Legendre polynomials for computational and user polar angle cosines
   */
  c_legendre_poly(ds->numu,mazim,ds->nstr,ds->nstr-1,ds->umu,ylmu);
  c_legendre_poly(nn,      mazim,ds->nstr,ds->nstr-1,cmu,    ylmc);

  /*
   * Evaluate Legendre polynomials with negative arguments from those with positive arguments;
   * Dave/Armstrong eq. (15), STWL(59)
   */
  sgn = -1.0;
  for (l = mazim; l <= ds->nstr-1; l++) {
    sgn *= -1;
    for (iq = nn+1; iq <= ds->nstr; iq++) {
      YLMC(l,iq) = sgn*YLMC(l,iq-nn);
    }
  }

  /*
   * Zero out bottom reflectivity (ALBEDO is used only in analytic formulae involving ALBEDO = 0
   * solutions; eqs 16-17 of Ref S2)
   */
  memset(bdr,0,(ds->nstr/2)*((ds->nstr/2)+1)*sizeof(double));

  /*-------------------  BEGIN LOOP ON COMPUTATIONAL LAYERS  -------------*/
  for (lc = 1; lc <= ds->nlyr; lc++) {
    /*
     * Solve eigenfunction problem in eq. STWJ(8b), STWL(23f)
     */
    c_solve_eigen(ds,lc,ab,array,cmu,cwt,gl,mazim,nn,ylmc,cc,evecc,eval,kk,gc,wk);
    /*
     * Interpolate eigenvectors to user angles
     */
    c_interp_eigenvec(ds,lc,cwt,evecc,gl,gu,mazim,nn,wk,ylmc,ylmu);
  }
  /*------------------  END LOOP ON COMPUTATIONAL LAYERS  ---------------*/

  /*
   * Set coefficient matrix (CBAND) of equations
   * combining boundary and layer interface
   * conditions (in band-storage mode required by
   * LINPACK routines)
   */
  c_set_matrix(ds,bdr,cband,cmu,cwt,delm0,dtaucpr,gc,kk,lyrcut,&ncol,ncut,taucpr,wk);

  /*
   * LU-decompose the coeff. matrix (LINPACK)
   */
  ncd = 3*nn-1;
  c_sgbco(cband,(9*(ds->nstr/2)-2),ncol,ncd,ncd,ipvt,&rcond,z);
  if (1.+rcond == 1.) {
    c_errmsg("albtrans--sgbco says matrix near singular",DS_WARNING);
  }

  /*
   * First, illuminate from top; if only one layer, this will give us everything
   * Solve for constants of integration in homogeneous solution
   */
  c_solve1(ds,cband,TOP_ILLUM,ipvt,ncol,ncut,nn,b,ll);

  /*
   * Compute azimuthally-averaged intensity at user angles; gives albedo if multi-layer (eq. 9 of Ref S2);
   * gives both albedo and transmissivity if single layer (eqs. 3-4 of Ref S2)
   */
  c_albtrans_intensity(ds,out,gu,kk,ll,nn,taucpr,wk);

  /*
   * Get beam-incidence albedos from reciprocity principle
   */

  for (iu = 1; iu <= ds->numu/2; iu++) {
    ALBMED(iu) = U0U(iu+ds->numu/2,1);
  }
  if (ds->nlyr == 1) {
    for (iu = 1; iu <= ds->numu/2; iu++) {
      /*
       * Get beam-incidence transmissivities from reciprocity principle (1 layer);
       * flip them end over end to correspond to positive UMU instead of negative
       */
      TRNMED(iu) = U0U(ds->numu/2+1-iu,2)+exp(-TAUCPR(ds->nlyr)/UMU(iu+ds->numu/2));
    }
  }
  else {
    /*
     * Second, illuminate from bottom (if multiple layers)
     */
    c_solve1(ds,cband,BOT_ILLUM,ipvt,ncol,ncut,nn,b,ll);
    c_albtrans_intensity(ds,out,gu,kk,ll,nn,taucpr,wk);
    /*
     * Get beam-incidence transmissivities from reciprocity principle
     */
    for (iu = 1; iu <= ds->numu/2; iu++) {
      TRNMED(iu) = U0U(iu+ds->numu/2,1)+exp(-TAUCPR(ds->nlyr)/UMU(iu+ds->numu/2));
    }
  }

  if (ds->bc.albedo > 0.) {
    /*
     * Get spherical albedo and transmissivity
     */
    if (ds->nlyr == 1) {
      c_albtrans_spherical(ds,cmu,cwt,gc,kk,ll,nn,taucpr,&sphalb,&sphtrn);
    }
    else {
      c_albtrans_spherical(ds,cmu,cwt,gc,kk,ll,nn,taucpr,&sphtrn,&sphalb);
    }
    /*
     * Ref. S2, eqs. 16-17 (these eqs. have a simple physical interpretation
     * like that of adding-doubling eqs.)
     */
    for (iu = 1; iu <= ds->numu; iu++) {
      
      ALBMED(iu) += ds->bc.albedo/(1.-ds->bc.albedo*sphalb)*sphtrn*TRNMED(iu);
      TRNMED(iu) += ds->bc.albedo/(1.-ds->bc.albedo*sphalb)*sphalb*TRNMED(iu);
    }
  }
  /*
   * Return UMU to all positive values, to agree with ordering in ALBMED, TRNMED
   */
  ds->numu /= 2;
  for (iu = 1; iu <= ds->numu; iu++) {
    UMU(iu) = UMU(iu+ds->numu);
  }
  if (ds->flag.prnt[3]) {
    c_print_albtrans(ds,out);
  }
  
  /* CE: I want to output the the spherical albedo and transmittance, and use the */
  /* variables ALBMED and TRNMED for this. They are not used so far otherwise in uvspec */
  /* If somebody needs these variables I will include new variables for sphtrn and sphalb*/
  ALBMED(1)=sphalb;
  TRNMED(1)=sphtrn;
  
  return;
}

/*============================= end of c_albtrans() ======================*/

/*============================= c_albtrans_intensity() ===================*/

/*
   Computes azimuthally-averaged intensity at top and bottom of medium
   (related to albedo and transmission of medium by reciprocity principles;
   see Ref S2).  User polar angles are used as incident beam angles.
   (This is a very specializedversion of user_intensities)

   ** NOTE **  User input values of UMU (assumed positive) are temporarily in
               upper locations of  UMU  and corresponding negatives are in
               lower locations (this makes GU come out right); the contents
               of the temporary UMU array are:
                   -UMU(ds->numu),..., -UMU(1), UMU(1),..., UMU(ds->numu)

   I N P U T    V A R I A B L E S:

       ds     :  Disort state variables
       gu     :  Eigenvectors interpolated to user polar angles (i.e., g in eq. SC(1), STWL(31ab))
       kk     :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       ll     :  Constants of integration in eq. SC(1), obtained by solving scaled version of eq. SC(5);
                 exponential term of eq. SC(12) not included
       nn     :  Order of double-Gauss quadrature (NSTR/2)
       taucpr :  Cumulative optical depth (delta-M-scaled)

   O U T P U T    V A R I A B L E:

       out->u0u : Diffuse azimuthally-averaged intensity at top and bottom of medium (directly transmitted component,
                  corresponding to bndint in user_intensities, is omitted).

   I N T E R N A L    V A R I A B L E S:

       dtau   :  Optical depth of a computational layer
       palint :  Non-boundary-forced intensity component
       utaupr :  Optical depths of user output levels (delta-M scaled)
       wk     :  Scratch vector for saving 'EXP' evaluations
       All the exponential factors (i.e., exp1, expn,... etc.)
       come from the substitution of constants of integration in
       eq. SC(12) into eqs. S1(8-9).  All have negative arguments.

   Called by- c_albtrans
 -------------------------------------------------------------------*/

void c_albtrans_intensity(disort_state *ds,
			  disort_output *out,
                          double       *gu,
                          double       *kk,
                          double       *ll,
                          int           nn,
                          double       *taucpr,
                          double       *wk)
{
  register int
    iq,iu,iumax,iumin,lc,lu;
  double
    denom,dtau,exp1,exp2,expn,mu,palint,sgn,utaupr[2];

  UTAUPR(1) = 0.;
  UTAUPR(2) = TAUCPR(ds->nlyr);

  for (lu = 1; lu <= 2; lu++) {
    if (lu == 1) {
      iumin = ds->numu/2+1;
      iumax = ds->numu;
      sgn   = 1.;
    }
    else {
      iumin = 1;
      iumax = ds->numu/2;
      sgn   = -1.;
    }

    /*
     * Loop over polar angles at which albedos/transmissivities desired
     * ( upward angles at top boundary, downward angles at bottom )
     */
    for (iu = iumin; iu <= iumax; iu++) {
      mu = UMU(iu);
      /*
       * Integrate from top to bottom computational layer
       */
      palint = 0.;
      for (lc = 1; lc <= ds->nlyr; lc++) {
        dtau = TAUCPR(lc)-TAUCPR(lc-1);
        exp1 = exp((UTAUPR(lu)-TAUCPR(lc-1))/mu);
        exp2 = exp((UTAUPR(lu)-TAUCPR(lc  ))/mu);
        /*
         * KK is negative
         */
        for (iq = 1; iq <= nn; iq++) {
          WK(iq) = exp(KK(iq,lc)*dtau);
          denom  = 1.+mu*KK(iq,lc);
          if (fabs(denom) < 0.0001) {
            /*
             * L'Hospital limit
             */
            expn = dtau/mu*exp2;
          }
          else {
            expn = (exp1*WK(iq)-exp2)*sgn/denom;
          }
          palint += GU(iu,iq,lc)*LL(iq,lc)*expn;
        }
        /*
         * KK is positive
         */
        for (iq = nn+1; iq <= ds->nstr; iq++) {
          denom = 1.+mu*KK(iq,lc);
          if (fabs(denom) < 0.0001) {
            expn = -dtau/mu*exp1;
          }
          else {
            expn = (exp1-exp2*WK(ds->nstr+1-iq))*sgn/denom;
          }
          palint += GU(iu,iq,lc)*LL(iq,lc)*expn;
        }
      }
      U0U(iu,lu) = palint;
    }
  }

  return;
}

/*============================= end of c_albtrans_intensity() ============*/

/*============================= c_print_albtrans() =======================*/

/*
   Print planar albedo and transmissivity of medium as a function of
   incident beam angle

   Called by- c_albtrans
 --------------------------------------------------------------------*/

void c_print_albtrans(disort_state  *ds,
                      disort_output *out)
{
  register int
    iu;

  fprintf(stdout,"\n\n\n *******  Flux Albedo and/or Transmissivity of entire medium  ********\n");
  fprintf(stdout," Beam Zen Ang   cos(Beam Zen Ang)      Albedo   Transmissivity\n");
  for (iu = 1; iu <= ds->numu; iu++) {
    fprintf(stdout,"%13.4f%20.6f%12.5f%17.4e\n",acos(UMU(iu))/DEG,UMU(iu),ALBMED(iu),TRNMED(iu));
  }

  return;
}

/*============================= end of c_print_albtrans() ================*/

/*============================= c_solve1() ===============================*/

/*
     Construct right-hand side vector -b- for isotropic incidence
     (only) on either top or bottom boundary and solve system
     of equations obtained from the boundary conditions and the
     continuity-of-intensity-at-layer-interface equations

     I N P U T      V A R I A B L E S:

       ds       :  Disort state variables
       cband    :  Left-hand side matrix of banded linear system
                   eq. SC(5), scaled by eq. SC(12); assumed already
                   in LU-decomposed form, ready for LINPACK solver
       ihom     :  Direction-of-illumination flag (TOP_ILLUM, top; BOT_ILLUM, bottom)
       ipvt     :
       ncol     :  Number of columns in CBAND
       ncut     :
       nn       :  Order of double-Gauss quadrature (NSTR/2)

    O U T P U T     V A R I A B L E S:

       b        :  Right-hand side vector of eq. SC(5) going into
                   sgbsl; returns as solution vector of eq.
                   SC(12), constants of integration without
                   exponential term
       ll       :  permanent storage for -b-, but re-ordered


    I N T E R N A L    V A R I A B L E S:

       ipvt     :  INTEGER vector of pivot indices
       ncd      :  Number of diagonals below or above main diagonal

   Called by- c_albtrans
   Calls- c_sgbsl
 +-------------------------------------------------------------------+
*/

void c_solve1(disort_state *ds,
              double       *cband,
              int           ihom,
              int          *ipvt,
              int           ncol,
              int           ncut,
              int           nn,
              double       *b,
              double       *ll)
{
  register int
    i,ipnt,iq,lc,ncd;

  memset(b,0,ds->nstr*ds->nlyr*sizeof(double));

  if (ihom == TOP_ILLUM) {
    /*
     * Because there are no beam or emission sources, remainder of B array is zero
     */
    for (i = 1; i <= nn; i++) {
      B(i)         = ds->bc.fisot;
      B(ncol-nn+i) = 0.;
    }
  }
  else if (ihom == BOT_ILLUM) {
    for (i = 1; i <= nn; i++) {
      B(i)         = 0.;
      B(ncol-nn+i) = ds->bc.fisot;
    }
  }
  else {
    c_errmsg("solve1---unrecognized ihom",DS_ERROR);
  }

  ncd = 3*nn-1;
  c_sgbsl(cband,(9*(ds->nstr/2)-2),ncol,ncd,ncd,ipvt,b,0);
  for (lc = 1; lc <= ncut; lc++) {
    ipnt = lc*ds->nstr-nn;
    for (iq = 1; iq <= nn; iq++) {
      LL(nn-iq+1,lc) = B(ipnt-iq+1);
      LL(nn+iq,  lc) = B(ipnt+iq  );
    }
  }

  return;
}

/*============================= end of c_solve1() ========================*/

/*============================= c_albtrans_spherical() ===================*/

/*
    Calculates spherical albedo and transmissivity for the entire medium
    from the m=0 intensity components (this is a specialized version of fluxes)

    I N P U T    V A R I A B L E S:

       ds      :  Disort state variables
       cmu,cwt :  Abscissae, weights for Gaussian quadrature over angle cosine
       kk      :  Eigenvalues of coeff. matrix in eq. SS(7)
       gc      :  Eigenvectors at polar quadrature angles, SC(1)
       ll      :  Constants of integration in eq. SC(1), obtained by solving
                  scaled version of eq. SC(5); exponential term of eq. SC(12) not incl.
       nn      :  Order of double-Gauss quadrature (NSTR/2)

    O U T P U T   V A R I A B L E S:

       sflup   :  Up-flux at top (equivalent to spherical albedo due to
                  reciprocity).  For illumination from below it gives
                  spherical transmissivity

       sfldn   :  Down-flux at bottom (for single layer, equivalent to
                  spherical transmissivity due to reciprocity)

    I N T E R N A L   V A R I A B L E S:

       zint    :  Intensity of m=0 case, in eq. SC(1)

   Called by- c_albtrans
 --------------------------------------------------------------------*/

void c_albtrans_spherical(disort_state *ds,
                          double       *cmu,
                          double       *cwt,
                          double       *gc,
                          double       *kk,
                          double       *ll,
                          int           nn,
                          double       *taucpr,
                          double       *sflup,
                          double       *sfldn)
{
  register int
    iq,jq;
  double
    zint;

  *sflup = 0.;
  for (iq = nn+1; iq <= ds->nstr; iq++) {
    zint = 0.;
    for (jq = 1; jq <= nn; jq++) {
      zint += GC(iq,jq,1)*LL(jq,1)*exp(KK(jq,1)*TAUCPR(1));
    }
    for (jq = nn+1; jq <= ds->nstr; jq++) {
      zint += GC(iq,jq,1)*LL(jq,1);
    }
    *sflup += CWT(iq-nn)*CMU(iq-nn)*zint;
  }

  *sfldn = 0.;
  for (iq = 1; iq <= nn; iq++) {
    zint = 0.;
    for (jq = 1; jq <= nn; jq++) {
      zint += GC(iq,jq,ds->nlyr)*LL(jq,ds->nlyr);
    }
    for (jq = nn+1; jq <=ds->nstr; jq++) {
      zint += GC(iq,jq,ds->nlyr)*LL(jq,ds->nlyr)*exp(-KK(jq,ds->nlyr)*(TAUCPR(ds->nlyr)-TAUCPR(ds->nlyr-1)));
    }
    *sfldn += CWT(nn+1-iq)*CMU(nn+1-iq)*zint;
  }

  *sflup *= 2.;
  *sfldn *= 2.;

  return;
}

/*============================= end of c_albtrans_spherical() ============*/

/******************************************************************
 ********** End of ds->flag.ibcnd = SPECIAL_BC routines ***********
 ******************************************************************/

/*============================= c_errmsg() ===============================*/

/*
 * Print out a warning or error message;  abort if type == DS_ERROR
 */

#ifdef ENABLE_ORIGINAL_OUTPUT_HANDLING

#define MAX_WARNINGS 100

void c_errmsg(const char *messag,
              int   type)
{
  static int
    warning_limit = FALSE,
    num_warnings  = 0;

  if (type == DS_ERROR) {
    fprintf(stderr,"\n ******* ERROR >>>>>>  %s\n",messag);
    exit(1);
  }

  if (warning_limit) return;

  if (++num_warnings <= MAX_WARNINGS) {
    fprintf(stderr,"\n ******* WARNING >>>>>>  %s\n",messag);
  }
  else {
    fprintf(stderr,"\n\n >>>>>>  TOO MANY WARNING MESSAGES --  ','They will no longer be printed  <<<<<<<\n\n");
    warning_limit = TRUE;
  }

  return;
}

#undef MAX_WARNINGS

/*============================= end of c_errmsg() ========================*/

/*============================= c_write_bad_var() ========================*/

/*
   Write name of erroneous variable and return TRUE; count and abort
   if too many errors.

   Input : quiet  = VERBOSE or QUIET
           varnam = name of erroneous variable to be written
 ----------------------------------------------------------------------*/

int c_write_bad_var(int   quiet,
                    const char *varnam)
{
  const int
    maxmsg = 50;
  static int
    nummsg = 0;

  nummsg++;
  if (quiet != QUIET) {
    fprintf(stderr,"\n ****  Input variable %s in error  ****\n",varnam);
    if (nummsg == maxmsg) {
      c_errmsg("Too many input errors.  Aborting...",DS_ERROR);
    }
  }

  return TRUE;
}

/*============================= end of c_write_bad_var() =================*/

/*============================= c_write_too_small_dim() ==================*/

/*
    Write name of too-small symbolic dimension and the value it should be
    increased to;  return TRUE

    Input :  quiet  = VERBOSE or QUIET
             dimnam = name of symbolic dimension which is too small
             minval = value to which that dimension should be increased
 ----------------------------------------------------------------------*/

int c_write_too_small_dim(int   quiet,
                          const char *dimnam,
                          int   minval)
{
  if (quiet != QUIET) {
    fprintf(stderr," ****  Symbolic dimension %s should be increased to at least %d  ****\n",
            dimnam,minval);
  }

  return TRUE;
}

/*============================= end of c_write_too_small_dim =============*/

#endif /* ENABLE_ORIGINAL_OUTPUT_HANDLING */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Call tree:

   c_sgbco
       c_sasum
       c_sdot
       c_saxpy
       c_sgbfa
           c_isamax
           c_saxpy
           c_sscal
       c_sscal
   c_sgbsl
       c_sdot
       c_saxpy
   c_sgeco
       c_sasum
       c_sdot
       c_saxpy
       c_sgefa
           c_isamax
           c_saxpy
           c_sscal
       c_sscal
   c_sgesl
       c_sdot
       c_saxpy
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*============================= c_sgbco() ================================*/

/*
     Factors a real band matrix by Gaussian elimination and estimates the
     condition of the matrix.
     Revision date:  8/1/82
     Author:  Moler, C.B. (Univ. of New Mexico)

     If  RCOND  is not needed, sgbfa is slightly faster.
     To solve  A*X = B , follow sgbco by sgbsl.

     Inputs:
        abd     double(LDA,N), contains the matrix in band storage.
                The columns of the matrix are stored in the columns of abd
                and the diagonals of the matrix are stored in rows
                ml+1 through 2*ml+mu+1 of  abd.
                See the comments below for details.
        lda     int, the leading dimension of the array abd.
                lda must be >= 2*ml+mu+1.
        n       int,the order of the original matrix.
        ml      int, number of diagonals below the main diagonal.
                0 <= ml < n.
        mu      int, number of diagonals above the main diagonal.
                0 <= mu < n.
                more efficient if  ml <= mu.

     Outputs:
        abd     an upper triangular matrix in band storage and
                the multipliers which were used to obtain it.
                The factorization can be written  A = L*U  where
                L  is a product of permutation and unit lower
                triangular matrices and  U  is upper triangular.
        ipvt    int[n], an integer vector of pivot indices.
        rcond   double, an estimate of the reciprocal condition of A.
                For the system  A*X = B, relative perturbations
                in A and B of size epsilon may cause relative
                perturbations in  X  of size  epsilon/rcond.
                If rcond  is so small that the logical expression
                   1.+RCOND == 1.
                is true, then  A  may be singular to working
                precision.  In particular, rcond is zero if exact
                singularity is detected or the estimate underflows.
        z       double[n], a work vector whose contents are usually
                unimportant. If A is close to a singular matrix, then
                z is an approximate null vector in the sense that
                norm(a*z) = rcond*norm(a)*norm(z).

     Band storage:
           If A is a band matrix, the following program segment
           will set up the input (with unit-offset arrays):
                   ml = (band width below the diagonal)
                   mu = (band width above the diagonal)
                   m = ml+mu+1
                   for (j = 1; j <= n; j++) {
                     i1 = IMAX(1,j-mu);
                     i2 = IMIN(n,j+ml);
                     for (i = i1; i <= i2; i++) {
                       k = i-j+m;
                       ABD(K,J) = A(I,J);
                     }
                   }
           This uses rows ml+1 through 2*ml+mu+1 of abd.
           In addition, the first ml rows in abd are used for
           elements generated during the triangularization.
           The total number of rows needed in abd is 2*ml+mu+1.
           The ml+mu by ml+mu upper left triangle and the
           ml by ml lower right triangle are not referenced.

     Example:  if the original matrix is

           11 12 13  0  0  0
           21 22 23 24  0  0
            0 32 33 34 35  0
            0  0 43 44 45 46
            0  0  0 54 55 56
            0  0  0  0 65 66

      then  n = 6, ml = 1, mu = 2, lda >= 5  and abd should contain
            *  *  *  +  +  +  , * = not used
            *  * 13 24 35 46  , + = used for pivoting
            * 12 23 34 45 56
           11 22 33 44 55 66
           21 32 43 54 65  *

 --------------------------------------------------------------------*/

void c_sgbco(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             double *rcond,
             double *z)
{
  int
    info;
  register int
    is,j,ju,k,kb,kp1,l,la,lm,lz,m,mm;
  double
    anorm,ek,s,sm,t,wk,wkm,ynorm;

  /*
   * compute 1-norm of A
   */
  anorm = 0.;
  l  = ml+1;
  is = l+mu;
  for (j = 1; j <= n; j++) {
    anorm = MAX(anorm,c_sasum(l,&ABD(is,j)));
    if (is > ml+1) {
      is--;
    }
    if (j <= mu) {
      l++;
    }
    if (j >= n-ml) {
      l--;
    }
  }
  /*
   * factor
   */
  c_sgbfa(abd,lda,n,ml,mu,ipvt,&info);
  /*
   * rcond = 1/(norm(A)*(estimate of norm(inverse(A)))) .
   * estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E.
   * trans(A) is the transpose of A.  The components of E are
   * chosen to cause maximum local growth in the elements of W where
   * trans(U)*W = E. The vectors are frequently rescaled to avoid overflow.
   * solve trans(U)*W = E
   */
  ek = 1.;

  memset(z,0,n*sizeof(double));

  m  = ml+mu+1;
  ju = 0;
  for (k = 1; k <= n; k++) {
    if (Z(k) != 0.) {
      ek = F77_SIGN(ek,-Z(k));
    }
    if (fabs(ek-Z(k)) > fabs(ABD(m,k))) {
      s = fabs(ABD(m,k))/fabs(ek-Z(k));
      c_sscal(n,s,z);
      ek *= s;
    }
    wk  =  ek-Z(k);
    wkm = -ek-Z(k);
    s   = fabs(wk);
    sm  = fabs(wkm);
    if (ABD(m,k) != 0.) {
      wk  /= ABD(m,k);
      wkm /= ABD(m,k);
    }
    else {
      wk  = 1.;
      wkm = 1.;
    }
    kp1 = k+1;
    ju  = IMIN(IMAX(ju,mu+IPVT(k)),n);
    mm  = m;
    if (kp1 <= ju) {
      for (j = kp1; j <= ju; j++) {
        mm--;
        sm   += fabs(Z(j)+wkm*ABD(mm,j));
        Z(j) += wk*ABD(mm,j);
        s    += fabs(Z(j));
      }
      if (s < sm) {
        t  = wkm-wk;
        wk = wkm;
        mm = m;
        for (j = kp1; j <= ju; j++) {
          mm--;
          Z(j) += t*ABD(mm,j);
        }
      }
    }
    Z(k) = wk;
  }

  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);

  /*
   * solve trans(L)*Y = W
   */
  for (kb = 1; kb <= n; kb++) {
    k  = n+1-kb;
    lm = IMIN(ml,n-k);
    if (k < n) {
      Z(k) += c_sdot(lm,&ABD(m+1,k),&Z(k+1));
    }
    if (fabs(Z(k)) > 1.) {
      s = 1./fabs(Z(k));
      c_sscal(n,s,z);
    }

    l    = IPVT(k);
    t    = Z(l);
    Z(l) = Z(k);
    Z(k) = t;
  }

  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);

  ynorm = 1.;
  /*
   * solve L*V = Y
   */
  for (k = 1; k <= n; k++) {
    l    = IPVT(k);
    t    = Z(l);
    Z(l) = Z(k);
    Z(k) = t;
    lm   = IMIN(ml,n-k);
    if (k < n) {
      c_saxpy(lm,t,&ABD(m+1,k),&Z(k+1));
    }
    if (fabs(Z(k)) > 1.) {
      s = 1./fabs(Z(k));
      c_sscal(n,s,z);
      ynorm *= s;
    }
  }

  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);

  ynorm *= s;
  /*
   * solve  U*Z = W
   */
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (fabs(Z(k)) > fabs(ABD(m,k))) {
      s = fabs(ABD(m,k))/fabs(Z(k));
      c_sscal(n,s,z);
      ynorm *= s;
    }
    if (ABD(m,k) != 0.) {
      Z(k) /= ABD(m,k);
    }
    else {
      Z(k) = 1.;
    }
    lm = IMIN(k,m)-1;
    la = m-lm;
    lz = k-lm;
    t  = -z[k-1];
    c_saxpy(lm,t,&ABD(la,k),&Z(lz));
  }

  /*
   * make znorm = 1.
   */
  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);

  ynorm *= s;
  if(anorm != 0.) {
    *rcond = ynorm/anorm;
  }
  else {
    *rcond = 0.;
  }

  return;
}

/*============================= end of c_sgbco() =========================*/

/*============================= c_sgbfa() ================================*/

/*
    Factors a real band matrix by elimination.
    Revision date:  8/1/82
    Author:  Moler, C. B. (U. of New Mexico)
    c_sgbfa is usually called by c_sgbco, but it can be called
    directly with a saving in time if rcond is not needed.

    Inputs:  same as c_sgbco
    Outputs:
        abd,ipvt    same as c_sgbco
        info    int,
                = 0  normal value.
                = k  if  u(k,k) == 0.  This is not an error
                     condition for this subroutine, but it does
                     indicate that sgbsl will divide by zero if
                     called.  Use  rcond  in c_sgbco for a reliable
                     indication of singularity.
    (see c_sgbco for description of band storage mode)

    NOTE: using memset() to zero columns in abd
 ----------------------------------------------------------------*/

void c_sgbfa(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             int    *info)
{
  register int
    i0,j,j0,j1,ju,jz,k,kp1,l,lm,m,mm,nm1;
  double
    t;

  m     = ml+mu+1;
  *info = 0;
  /*
   * zero initial fill-in columns
   */
  j0 = mu+2;
  j1 = IMIN(n,m)-1;
  for (jz = j0; jz <= j1; jz++) {
    i0 = m+1-jz;
    memset(&ABD(i0,jz),0,(ml-i0+1)*sizeof(double));
  }
  jz = j1;
  ju = 0;

  /*
   * Gaussian elimination with partial pivoting
   */
  nm1 = n-1;
  for (k = 1; k <= nm1; k++) {
    kp1 = k+1;
   /*
    * zero next fill-in column
    */
    jz++;
    if (jz <= n) {
      memset(&ABD(1,jz),0,ml*sizeof(double));
    }
    /*
     * find L = pivot index
     */
    lm      = IMIN(ml,n-k);
    l       = c_isamax(lm+1,&ABD(m,k))+m-1;
    IPVT(k) = l+k-m;
    if (ABD(l,k) == 0.) {
     /*
      * zero pivot implies this column already triangularized
      */
      *info = k;
    }
    else {
      /*
       * interchange if necessary
       */
      if (l != m) {
        t        = ABD(l,k);
        ABD(l,k) = ABD(m,k);
        ABD(m,k) = t;
      }
      /*
       * compute multipliers
       */
      t = -1./ABD(m,k);
      c_sscal(lm,t,&ABD(m+1,k));
      /*
       * row elimination with column indexing
       */
      ju = IMIN(IMAX(ju,mu+IPVT(k)),n);
      mm = m;
      for (j = kp1; j <= ju; j++) {
        l--;
        mm--;
        t = ABD(l,j);
        if (l != mm) {
          ABD(l,j)  = ABD(mm,j);
          ABD(mm,j) = t;
        }
        c_saxpy(lm,t,&ABD(m+1,k),&ABD(mm+1,j));
      }
    }
  }
  IPVT(n) = n;
  if (ABD(m,n) == 0.) {
    *info = n;
  }

  return;
}

/*============================= end of c_sgbfa() =========================*/

/*============================= c_sgbsl() ================================*/

/*
    Solves the real band system
       A * X = B  or  transpose(A) * X = B
    using the factors computed by sgbco or sgbfa.
    Revision date:  8/1/82
    Author:  Moler, C. B. (Univ. of New Mexico)

    Inputs:
        abd     double(lda, n), the output from sgbco or sgbfa.
        lda     int, the leading dimension of the array abd.
        n       int, the order of the original matrix.
        ml      int, number of diagonals below the main diagonal.
        mu      int, number of diagonals above the main diagonal.
        ipvt    int(n), the pivot vector from sgbco or sgbfa.
        b       double(n), the right hand side vector.
        job     int,
                = 0         to solve  A*X = B ,
                = nonzero   to solve  transpose(A)*X = B

     Outputs:
        b       the solution vector  X

     Error condition:
        A division by zero will occur if the input factor contains a
        zero on the diagonal.  Technically, this indicates singularity,
        but it is often caused by improper arguments or improper
        setting of lda.  It will not occur if the subroutines are
        called correctly and if c_sgbco has set rcond > 0.0
        or sgbfa has set info = 0 .
     To compute  inverse(a)*c  where c is a matrix
     with p columns
      c_sgbco(abd,lda,n,ml,mu,ipvt,&rcond,z)
      if (rcond is too small) ...
        for (j = 1; j <= p; j++) {
          c_sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
        }
 --------------------------------------------------------*/

void c_sgbsl(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             double *b,
             int     job)
{
  register int
    k,kb,l,la,lb,lm,m,nm1;
  double
    t;

  m   = mu+ml+1;
  nm1 = n-1;
  if (job == 0) {
   /*
    * solve  A*X = B;  first solve L*Y = B
    */
    if (ml != 0) {
      for (k = 1; k <= nm1; k++) {
        lm = IMIN(ml,n-k);
        l  = IPVT(k);
        t  = B(l);
        if (l != k) {
          B(l) = B(k);
          B(k) = t;
        }
        c_saxpy(lm,t,&ABD(m+1,k),&B(k+1));
      }
    }
    /*
     * now solve  U*X = Y
     */
    for (kb = 1; kb <= n; kb++) {
      k     = n+1-kb;
      B(k) /= ABD(m,k);
      lm    = IMIN(k,m)-1;
      la    = m-lm;
      lb    = k-lm;
      t     = -B(k);
      c_saxpy(lm,t,&ABD(la,k),&B(lb));
    }
  }
  else {
    /*
     * solve  trans(A)*X = B; first solve trans(U)*Y = B
     */
    for (k = 1; k <= n; k++) {
      lm   = IMIN(k,m)-1;
      la   = m-lm;
      lb   = k-lm;
      t    = c_sdot(lm,&ABD(la,k),&B(lb));
      B(k) = (B(k)-t)/ABD(m,k);
    }
    /*
     * now solve trans(L)*X = Y
     */
    if (ml != 0) {
      for (kb = 1; kb <= nm1; kb++) {
        k     = n-kb;
        lm    = IMIN(ml,n-k);
        B(k) += c_sdot(lm,&ABD(m+1,k),&B(k+1));
        l     = IPVT(k);
        if (l != k) {
          t    = B(l);
          B(l) = B(k);
          B(k) = t;
        }
      }
    }
  }

  return;
}

/*============================= end of c_sgbsl() =========================*/

/*============================= c_sgeco() ================================*/

/*
   Factors a real matrix by Gaussian elimination
   and estimates the condition of the matrix.
   Revision date:  8/1/82
   Author:  Moler, C. B. (Univ. of New Mexico)
   If rcond is not needed, sgefa is slightly faster.
   To solve  A*X = B, follow sgeco by sgesl.

     Inputs:
        a       double(lda, n), the matrix to be factored.
        lda     int, the leading dimension of the array a.
        n       int, the order of the matrix a.

     Outputs:
        a       an upper triangular matrix and the multipliers
                which were used to obtain it.
                The factorization can be written  A = L*U , where
                L  is a product of permutation and unit lower
                triangular matrices and U is upper triangular.
        ipvt    int(n), an integer vector of pivot indices.
        rcond   double, an estimate of the reciprocal condition of a.
                For the system A*X = B, relative perturbations
                in A and B of size epsilon may cause relative
                perturbations in X of size epsilon/rcond.
                If rcond is so small that the logical expression
                  1.+rcond == 1.
                is true, then A may be singular to working precision.
                In particular, rcond is zero if exact singularity
                is detected or the estimate underflows.
        z       double(n), a work vector whose contents are usually
                unimportant. If A is close to a singular matrix, then z 
                is an approximate null vector in the sense that
                norm(A*Z) = rcond*norm(A)*norm(Z) .
 ------------------------------------------------------------------*/

void c_sgeco(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             double *rcond,
             double *z)
{
  int
    info;
  register int
    j,k,kb,kp1,l;
  double
    anorm,ek,s,sm,t,wk,wkm,ynorm;

  /*
   * compute 1-norm of A
   */
  anorm = 0.;
  for (j = 1; j <= n; j++) {
    anorm = MAX(anorm,c_sasum(n,&A(1,j)));
  }

  /*
   * factor
   */
  c_sgefa(a,lda,n,ipvt,&info);

  /*
   * rcond = 1/(norm(A)*(estimate of norm(inverse(A)))).
   * estimate = norm(Z)/norm(Y) where A*Z = Y and trans(A)*Y = E.
   * trans(A) is the transpose of A. The components of E are
   * chosen to cause maximum local growth in the elements of W where
   * trans(U)*W = E.  The vectors are frequently rescaled to avoid overflow.
   * solve trans(U)*W = E
   */
  ek = 1.;
  memset(z,0,n*sizeof(double));

  for (k = 1; k <= n; k++) {
    if (Z(k) != 0.) {
      ek = F77_SIGN(ek,-Z(k));
    }
    if (fabs(ek-Z(k)) > fabs(A(k,k))) {
      s = fabs(A(k,k))/fabs(ek-Z(k));
      c_sscal(n,s,z);
      ek *= s;
    }
    wk  =  ek-Z(k);
    wkm = -ek-Z(k);
    s   = fabs(wk);
    sm  = fabs(wkm);
    if (A(k,k) != 0.) {
      wk  /= A(k,k);
      wkm /= A(k,k);
    }
    else {
      wk  = 1.;
      wkm = 1.;
    }
    kp1 = k+1;
    if (kp1 <= n) {
      for (j = kp1; j <= n; j++) {
        sm   += fabs(Z(j)+wkm*A(k,j));
        Z(j) += wk*A(k,j);
        s    += fabs(Z(j));
      }
      if (s < sm) {
        t  = wkm-wk;
        wk = wkm;
        for (j = kp1; j <= n; j++) {
          Z(j) += t*A(k,j);
        }
      }
    }
    Z(k) = wk;
  }

  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);
  /*
   * solve trans(L)*Y = W
   */
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (k < n) {
      Z(k) += c_sdot(n-k,&A(k+1,k),&Z(k+1));
    }
    if (fabs(Z(k)) > 1.) {
      s = 1./fabs(Z(k));
      c_sscal(n,s,z);
    }
    l    = IPVT(k);
    t    = Z(l);
    Z(l) = Z(k);
    Z(k) = t;
  }
  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);
  /*
   * solve L*V = Y
   */
  ynorm = 1.;
  for (k = 1; k <= n; k++) {
    l    = IPVT(k);
    t    = Z(l);
    Z(l) = Z(k);
    Z(k) = t;
    if (k < n) {
      c_saxpy(n-k,t,&A(k+1,k),&Z(k+1));
    }
    if (fabs(Z(k)) > 1.) {
      s = 1./fabs(Z(k));
      c_sscal(n,s,z);
      ynorm *= s;
    }
  }
  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);
  /*
   * solve U*Z = V
   */
  ynorm *= s;
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (fabs(Z(k)) > fabs(A(k,k))) {
      s = fabs(A(k,k))/fabs(Z(k));
      c_sscal(n,s,z);
      ynorm *= s;
    }
    if (A(k,k) != 0.) {
      Z(k) /= A(k,k);
    }
    else {
      Z(k) = 1.;
    }
    t = -Z(k);
    c_saxpy(k-1,t,&A(1,k),&Z(1));
  }
  /*
   * make znorm = 1.0
   */
  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);
  ynorm *= s;
  if (anorm != 0.) {
    *rcond = ynorm/anorm;
  }
  else {
    *rcond = 0.;
  }

  return;
}

/*============================= end of c_sgeco() =========================*/

/*============================= c_sgefa() ================================*/

/*
   Factors a real matrix by Gaussian elimination.
   Revision date:  8/1/82
   Author:  Moler, C. B. (Univ. of New Mexico)
   c_sgefa is usually called by c_sgeco, but it can be called directly with a
   saving in time if rcond is not needed.
   (time for c_sgeco) = (1+9/n)*(time for c_sgefa).

   Inputs:  same as c_sgeco

   Outputs:
        a,ipvt  same as c_sgeco
        info    int,
                = 0  normal value.
                = k  if  u(k,k) = 0.  This is not an error condition for
                     this subroutine, but it does indicate that c_sgesl or
                     c_sgedi will divide by zero if called.  Use rcond in
                     c_sgeco for a reliable indication of singularity.
 ---------------------------------------------------------------------*/

void c_sgefa(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             int    *info)
{
  register int
    j,k,kp1,l,nm1;
  double
    t;

  /*
   * Gaussian elimination with partial pivoting
   */
  *info = 0;
  nm1   = n-1;
  for (k = 1; k <= nm1; k++) {
    kp1 = k+1;
    /*
     * find L = pivot index
     */
    l       = c_isamax(n-k+1,&A(k,k))+k-1;
    IPVT(k) = l;
    if (A(l,k) == 0.) {
      /*
       * zero pivot implies this column already triangularized
       */
      *info = k;
    }
    else {
      /*
       * interchange if necessary
       */
      if (l != k) {
        t      = A(l,k);
        A(l,k) = A(k,k);
        A(k,k) = t;
      }
      /*
       * compute multipliers
       */
      t = -1./A(k,k);
      c_sscal(n-k,t,&A(k+1,k));
      /*
       * row elimination with column indexing
       */
      for (j = kp1; j <= n; j++) {
        t = A(l,j);
        if (l != k) {
          A(l,j) = A(k,j);
          A(k,j) = t;
        }
        c_saxpy(n-k,t,&A(k+1,k),&A(k+1,j));
      }
    }
  }
  IPVT(n) = n;
  if (A(n,n) == 0.) {
    *info = n;
  }

  return;
}

/*============================= end of c_sgefa() =========================*/

/*============================= c_sgesl() ================================*/

/*
  Solves the real system
     A*X = B  or  transpose(A)*X = B
  using the factors computed by sgeco or sgefa.
  Revision date:  8/1/82
  Author:  Moler, C. B. (Univ. of New Mexico)

     Inputs:
        a       double(lda, n), the output from sgeco or sgefa.
        lda     int, the leading dimension of the array  A
        n       int, the order of the matrix  A
        ipvt    int(n), the pivot vector from sgeco or sgefa.
        b       double(n), the right hand side vector.
        job     int, 
                = 0         to solve  A*X = B ,
                = nonzero   to solve  transpose(A)*X = B

     Outputs:
        b       the solution vector x

     Error condition:
        A division by zero will occur if the input factor contains a
        zero on the diagonal. Technically, this indicates singularity,
        but it is often caused by improper arguments or improper setting
        of lda. It will not occur if the subroutines are called correctly
        and if sgeco has set rcond > 0. or sgefa has set info = 0 .
     To compute  inverse(a)*c where c is a matrix with p columns
           c_sgeco(a,lda,n,ipvt,rcond,z);
           if (rcond is too small) ...
           for (j = 1; j <= p; j++) {
             c_sgesl(a,lda,n,ipvt,c(1,j),0);
           }
 ---------------------------------------------------------------------*/

void c_sgesl(double *a,
             int     lda,
             int     n,
             int    *ipvt,
             double *b,
             int     job)
{
  register int
    k,kb,l,nm1;
  double
    t;

  nm1 = n-1;
  if (job == 0) {
    /*
     * solve  A*X = B; first solve L*Y = B
     */
    for (k = 1; k <= nm1; k++) {
      l = IPVT(k);
      t = B(l);
      if (l != k) {
        B(l) = B(k);
        B(k) = t;
      }
      c_saxpy(n-k,t,&A(k+1,k),&B(k+1));
    }
    /*
     * now solve  U*X = Y
     */
    for (kb = 1; kb <= n; kb++) {
      k     = n+1-kb;
      B(k) /= A(k,k);
      t     = -B(k);
      c_saxpy(k-1,t,&A(1,k),&B(1));
    }
  }
  else {
    /*
     * solve trans(A)*X = B; first solve trans(U)*Y = B
     */
    for (k = 1; k <= n; k++) {
      t    = c_sdot(k-1,&A(1,k),&B(1));
      B(k) = (B(k)-t)/A(k,k);
    }
    /*
     * now solve  trans(l)*x = y
     */
    for (kb = 1; kb <= nm1; kb++) {
      k = n-kb;
      B(k) += c_sdot(n-k,&A(k+1,k),&B(k+1));
      l     = IPVT(k);
      if (l != k) {
        t    = B(l);
        B(l) = B(k);
        B(k) = t;
      }
    }
  }

  return;
}

/*============================= end of c_sgesl() =========================*/

/*============================= c_sasum() ================================*/

/*
  Input--   n     Number of elements in vector to be summed
            sx    array, length n, containing vector

  OUTPUT--  ans   Sum from i = 1 to n of fabs(SX(i))

  NOTE: Fortran input incx removed because it is not used by
        disort or twostr
 ----------------------------------------------------------*/

double c_sasum(int     n,
             double *sx)
{
  register int
    i,m;
  double
    ans;

  ans = 0.;
  if (n <= 0) {
    return ans;
  }

  m = n%4;
  if (m != 0) {
    /*
     * clean-up loop so remaining vector length is a multiple of 4.
     */
    for (i = 1; i <= m; i++) {
      ans += fabs(SX(i));
    }
  }
  /*
   * unroll loop for speed
   */
  for (i = m+1; i <= n; i+=4) {
    ans += fabs(SX(i  ))
          +fabs(SX(i+1))
          +fabs(SX(i+2))
          +fabs(SX(i+3));
  }

  return ans;
}

/*============================= end of c_sasum() =========================*/

/*============================= c_saxpy() ================================*/

/*
  y = a*x + y  (x, y = vectors, a = scalar)

  INPUT--
        n   Number of elements in input vectors x and y
       sa   Scalar multiplier a
       sx   Array containing vector x
       sy   Array containing vector Y

 OUTPUT--
       sy   For i = 1 to n, overwrite  SY(i) with sa*SX(i)+SY(i)

  NOTE: Fortran inputs incx, incy removed because they are not used
        by disort or twostr
 ------------------------------------------------------------*/

void c_saxpy(int     n,
             double  sa,
             double *sx,
             double *sy)
{
  register int
    i,m;

  if (n <= 0 || sa == 0.) {
    return;
  }

  m = n%4;
  if (m != 0) {
    /*
     * clean-up loop so remaining vector length is a multiple of 4.
     */
    for (i = 1; i <= m; i++) {
      SY(i) += sa*SX(i);
    }
  }
  /*
   * unroll loop for speed
   */
  for (i = m+1; i <= n; i+=4) {
    SY(i  ) += sa*SX(i  );
    SY(i+1) += sa*SX(i+1);
    SY(i+2) += sa*SX(i+2);
    SY(i+3) += sa*SX(i+3);
  }

  return;
}

/*============================= c_saxpy() ================================*/

/*============================= c_sdot() =================================*/

/*
  Dot product of vectors x and y

  INPUT--
        n  Number of elements in input vectors x and y
       sx  Array containing vector x
       sy  Array containing vector y

 OUTPUT--
      ans  Sum for i = 1 to n of  SX(i)*SY(i),

  NOTE: Fortran input arguments incx, incy removed because they
        are not used in disort or twostr
 ------------------------------------------------------------------*/

double c_sdot(int     n,
              double *sx,
              double *sy)
{
  register int
    i,m;
  double
    ans;

  ans = 0.;
  if (n <= 0) {
    return ans;
  }

  m = n%4;
  if (m != 0) {
    /*
     * clean-up loop so remaining vector length is a multiple of 4.
     */
    for (i = 1; i <= m; i++) {
      ans += SX(i)*SY(i);
    }
  }
  /*
   * unroll loop for speed
   */
  for (i = m+1; i <= n; i+=4) {
    ans += SX(i  )*SY(i  )
          +SX(i+1)*SY(i+1)
          +SX(i+2)*SY(i+2)
          +SX(i+3)*SY(i+3);
  }

  return ans;
}

/*============================= end of c_sdot() ==========================*/

/*============================= c_sscal() ================================*/

/*
  Multiply vector sx by scalar sa

  INPUT--  n  Number of elements in vector
          sa  Scale factor
          sx  Array, length n, containing vector

 OUTPUT-- sx  Replace SX(i) with sa*SX(i) for i = i to n

 NOTE: Fortran input argument incx removed since it is not used
       in disort or twostr

 ---------------------------------------------------------------------*/

void c_sscal(int    n,
             double  sa,
             double *sx)
{
  register int
    i,m;

  if (n <= 0) {
    return;
  }
  m = n%4;
  if (m != 0) {
    /*
     * clean-up loop so remaining vector length is a multiple of 4.
     */
    for (i = 1; i <= m; i++) {
      SX(i) *= sa;
    }
  }
  /*
   * unroll loop for speed
   */
  for (i = m+1; i <= n; i+=4) {
    SX(i  ) *= sa;
    SX(i+1) *= sa;
    SX(i+2) *= sa;
    SX(i+3) *= sa;
  }

  return;
}

/*============================= end of c_sscal() =========================*/

/*============================= c_isamax() ===============================*/

/*
 INPUT--  n        Number of elements in vector of interest
          sx       Array, length n, containing vector

 OUTPUT-- ans      First i, i = 1 to n, to maximize fabs(SX(i))

 NOTE: Fortran input incx removed because it is not used by
       disort or twostr
 ---------------------------------------------------------------------*/

int c_isamax(int     n,
             double *sx)
{
  register int
    ans=0,i;
  double
   smax,xmag;

  if (n <= 0) {
    ans = 0;
  }
  else if (n == 1) {
    ans = 1;
  }
  else {
    smax = 0.;
    for (i = 1; i <= n; i++) {
      xmag = fabs(SX(i));
      if (smax < xmag) {
        smax = xmag;
        ans  = i;
      }
    }
  } 

  return ans;
}

/*============================= end of c_isamax() ========================*/

/*
 * Shift macros that are different than for disort
 */
#undef  KK
#define KK(lyu) kk[lyu-1]

/*============================= c_twostr() ===============================*/

/*
 Copyright (C) 1993, 1994, 1995 Arve Kylling

 C rewrite by Timothy E. Dowling (Univ. of Louisville)

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 1, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 To obtain a copy of the GNU General Public License write to the
 Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
 USA.

+---------------------------------------------------------------------+

     AUTHOR :  Arve Kylling (July 1993)
               Arve.Kylling@itek.norut.no

     REFERENCES (cited in the programs using the acronyms shown):

     DS: Dahlback, A. and K. Stamnes 1991: A new spherical
      model for computing the radiation field available
      for photolysis and heating at twilight, Planet.
      Space Sci. 39, 671-683.

     KS: Kylling, A., and K. Stamnes 1992: Efficient yet accurate
      solution of the linear transport  equation in the 
      presence of internal sources: the exponential-linear
      approximation, J. Comp. Phys. 102, 265-276.

    KST: Kylling, A., K. Stamnes and S.-C. Tsay 1995: A reliable
      and efficient two-stream algorithm for radiative
      transfer; Documentation of accuracy in realistic
      layered media, in print, Journal of Atmospheric
      Chemistry 21, 115-150.

    STWJ: Stamnes, K., S.-C. Tsay, W. Wiscombe and K. Jayaweera
      1988: Numerically stable algorithm for discrete-
      ordinate-method radiative transfer in multiple
      scattering and emitting layered media, Appl.
      Optics., 27, 2502.

    WW: Wiscombe, W., 1977:  The Delta-M Method: Rapid Yet
      Accurate Radiative Flux Calculations, J. Atmos. Sci.
      34, 1408-1422

+---------------------------------------------------------------------+

    I n t r o d u c t o r y    n o t e

    (References are given as author-last-name strings, e.g., KST.)

    twostr() solves the radiative transfer equation in an absorbing,
    emitting and multiple scattering, layered pseudo-spherical
    medium in the two-stream approximation. For a discussion of the
    theory behind the present implementation see (KST).

    twostr() is based on the general n-stream algorithm DISORT 
    described in Stamnes et al. (1988, STWJ), and incorporates 
    all the advanced features of that algorithm. Furthermore it
    has been extended to include spherical geometry using the 
    perturbation approach of Dahlback and Stamnes (1991). Relative
    to DISORT, it is both simplified and extended as follows:

     1) The only quantities calculated are mean intensities and fluxes.

     2) The medium may be taken to be pseudo-spherical (flag.spher is TRUE)

     3) Only Lambertian reflection at the bottom boundary is allowed


    General remarks about the structure of the input/output parameters

    The list of input variables is more easily comprehended if
    the following simple facts are borne in mind :

    * there is one vertical coordinate, measured in optical depth units;

    * the layers necessary for computational purposes are entirely
      decoupled from the levels at which the user desires results.

    The computational layering is usually constrained by the problem,
    in the sense that each computational layer must be reasonably
    homogeneous and not have a temperature variation of more than
    about 20 K across it (if thermal sources are considered).
    For example, a clear layer overlain by a cloud overlain by a
    dusty layer would suggest three computational layers.

    However the radiant quantities can be returned to the user at ANY
    level.  For example, the user may have picked 3 computational 
    layers, but he can then request intensities from e.g. only the
    middle of the 2nd layer.

+---------------------------------------------------------------------+

    I n p u t    v a r i a b l e s

    Note on units:

       The radiant output units are determined by the sources of
    radiation driving the problem.  Lacking thermal emission, the
    radiant output units are the units of the sources ds.bc.fbeam and
    ds.bc.fisot.
       If thermal emission of any kind is included, subprogram planck_func2()
    determines the units.  The default planck_func2() has mks units [w/sq m].
    ds.bc.fbeam and ds.bc.fisot must have the same units as planck_func2() when
    thermal emission is present.


    ********  Computational layer structure  ********

        ===========================================================
        == Note:  Layers are numbered from the top boundary down ==
        ===========================================================

    ds.nlyr     Number of computational layers

    DTAUC(lc)   lc = 1 to ds.nlyr,
                optical depths of computational layers

    SSALB(lc)   lc = 1 to ds.nlyr,
                single-scatter albedos of computational layers

    GG(lc)      lc = 1 to ds.nlyr,
                asymmetry factor of computational layers
                Should be <= 1.0 (complete forward scattering) and
                >= -1.0 (complete backward scattering).
                NOTE. GG is changed by twostr() if deltam = TRUE.

    TEMPER(lev) lev = 0 to ds.nlyr, temperatures [K] of levels.
                (Note that temperature is specified at levels
                rather than for layers.)  Don't forget to put top
                temperature in 'TEMPER(0)', not 'TEMPER(1)'.  Top and
                bottom values do not need to agree with top and
                bottom boundary temperatures ds.bc.ttemp and ds.bc.btemp 
                (i.e. slips are allowed).
                Needed only if ds.flag.planck is TRUE.

    ZD(lev)     lev = 0 to ds.nlyr, altitude of level above
                the ground, i.e. ZD(nlyr) = 0., the surface of 
                the planet. Typically in units of (km)
                Must have same units as -radius-. 
                Used to calculate the Chapman function when
                spherical geometry is needed. 
                Needed only if flag.spher is TRUE.

    ds.wvnmlo,  Wavenumbers (inv cm) of spectral interval
      ds.wvnmhi ( used only for calculating Planck function )
                needed only if ds.flag.planck is true.
                If ds.wvnmlo < ds.wvnmhi the Planck function is
                integrated over this interval. If ds.wvnmlo ==  ds.wvnmhi
                the Planck function at wvnmlo is returned.


    ********  User level organization  ********

    ds.flag.usrtau = FALSE, radiant quantities are to be returned
                     at boundary of every computational layer.

                   = TRUE,  radiant quantities are to be returned
                     at user-specified optical depths, as follows:

    ds.ntau        Number of optical depths

    UTAU(lu)       lu = 1 to ds.ntau, user optical depths, in increasing order.
                   UTAU(ntau) must be no greater than the total optical depth of the medium.

     ******** Top and bottom boundary conditions  *********

    ds.bc.fbeam : Intensity of incident parallel beam at top boundary.
                  (units w/sq m if thermal sources active, otherwise
                  arbitrary units).  Corresponding incident flux
                  is  'umu0'  times 'fbeam'.  Note that this is an
                  infinitely wide beam, not a searchlight beam.

    ds.bc.umu0  : Polar angle cosine of incident beam.

    ds.bc.fisot : Intensity of top-boundary isotropic illumination.
                  (units w/sq m if thermal sources active, otherwise
                  arbitrary units).  Corresponding incident flux
                  is  pi (M_PI = 3.14159...)  times 'fisot'.

    ds.bc.albedo: Bottom-boundary albedo, bottom boundary is
                  assumed to be Lambert reflecting.

    ds.bc.btemp : Temperature of bottom boundary (K)  (bottom
                  emissivity is calculated from -albedo-,
                  so it need not be specified).
                  Needed only if -planck- is true.

    ds.bc.ttemp : Temperature of top boundary (K)
                  Needed only if -planck- is true.

    ds.bc.temis : Emissivity of top boundary
                  Needed only if -planck- is true.

    radius      : Distance from center of planet to the planets 
                  surface (km) 

    **********  Control flags  **************

    ds.flag.planck  = TRUE, include thermal emission
                      FALSE, ignore all thermal emission (saves computer time)
                     ( If ds.flag.planck = FALSE, it is not necessary to set any of
                      the variables having to do with thermal emission )

    ds.flag.prnt[0] = TRUE, print input variables 
    ds.flag.prnt[1] = TRUE, print fluxes, mean intensities and flux divergence.

    deltam  = TRUE,  use delta-m method ( see Wiscombe, 1977 )
            = FALSE, don't use delta-m method
            In general intensities and fluxes will be more accurate
            for phase functions with a large forward peak (i.e.
            an asymmetry factor close to 1.) if 'deltam' is set true.

    ds.flag.spher = TRUE, spherical geometry accounted for. In this case
                    -radius- and -zd- must be set also. NOTE: this option
                    increases the execution time, hence use it only when
                    necessary if speed is of concern.
                  = FALSE, plane-parallel atmosphere assumed

    ds.header   : A 127- (or less) character header for prints


+---------------------------------------------------------------------+
               O u t p u t    v a r i a b l e s
+---------------------------------------------------------------------+

    == Note on units == If thermal sources are specified, fluxes come
                        out in [w/sq m] and intensities in [w/sq m/steradian].  
                        Otherwise, the flux and intensity units are determined
                        by the units of -fbeam- and -fisot-.

    If ds.flag.usrtau = FALSE :

         ds.ntau      Number of optical depths at which radiant
                      quantities are evaluated ( = nlyr+1 )

         UTAU(lu)     lu = 1 to ntau, optical depths, in increasing
                      order, corresponding to boundaries of
                      computational layers (see -dtauc-)

    RFLDIR(lu)    :   Direct-beam flux (without delta-m scaling)

    RFLDN(lu)     :   Diffuse down-flux (total minus direct-beam)
                      (without delta-m scaling)

    FLUP(lu)      :   Diffuse up-flux

    DFDT(lu)      :   Flux divergence  d(net flux)/d(optical depth),
                      where 'net flux' includes the direct beam
                      (an exact result;  not from differencing fluxes)

    UAVG(lu)      :   Mean intensity (including the direct beam)

    IERROR(i)     :   Error flag array, if IERROR(i) is zero everything
                      is ok, otherwise twostr() found a fatal error, in this
                      case, twostr return immediately and reports the error in IERROR.

                      i =  1 : ds.nlyr <  1

                      i =  3 : dtauc   <  0.
                      i =  4 : ssalb   <  0. || ssalb > 1.
                      i =  5 : temper  <  0. 
                      i =  6 : gg      < -1. || gg > 1.
                      i =  7 : ZD(lc)  >  ZD(lc-1)
                      i =  8 : ds.ntau <  1

                      i = 10 : UTAU(lu) < 0. || UTAU(lu) > TAUC(nlyr) 

                      i = 12 : fbeam    < 0.
                      i = 13 : if flag.spher = FALSE
                                    umu0  < 0. || umu0 > 1.
                               if flag.spher = TRUE
                                    umu0  < 0. || umu0 > 1.
                      i = 14 : fisot   < 0.
                      i = 15 : albedo  < 0. || albedo > 1.
                      i = 16 : wvnmlo  < 0. || wvnmhi < wvnmlo
                      i = 17 : temis   < 0. || temis  > 1.
                      i = 18 : btemp   < 0.
                      i = 19 : ttemp   < 0.

                      i = 22 : !ds->flag.usrtau && ds->ntau < ds->nlyr+1

                     NOTE: i = 2, 9, 11, 20, and 21 are eliminated in the C version by the
                           change from static to dynamic memory allocation

+---------------------------------------------------------------------+

                 I/O variable specifications

+---------------------------------------------------------------------+
      Routines called (in order): c_twostr_check_inputs, c_twostr_set, c_twostr_print_inputs,
                                  c_twostr_solns, c_set_matrix, c_twostr_solve_bc, c_twostr_fluxes
+---------------------------------------------------------------------+

  Index conventions (for all loops and all variable descriptions):

     iq     :  For quadrature angles
     lu     :  For user levels
     lc     :  For computational layers (each having a different single-scatter albedo and/or phase function)
     lev    :  For computational levels
     ls     :  Runs from 0 to 2*ds->nlyr+1, ls = 1,2,3 refers to top, center and bottom of layer 1, 
               ls = 3,4,5 refers to top, center and bottom of layer 2, etc.

+---------------------------------------------------------------------+

               I n t e r n a l    v a r i a b l e s

   B()...........Right-hand side vector of eqs. KST(38-41), set in twostr_solve_bc()
   bplanck.......Intensity emitted from bottom boundary
   CBAND().......Matrix of left-hand side of the linear system eqs. KST(38-41); in tridiagonal form
   CH(lc)........The Chapman-factor to correct for pseudo-spherical geometry in the direct beam.
   CHTAU(lc).....The optical depth in spherical geometry.
   cmu...........Computational polar angle, single or double Gaussian quadrature rule used, see twostr_set()
   EXPBEA(lc)....Transmission of direct beam in delta-m optical depth coordinates
   FLDIR(lu).....Direct beam flux (delta-m scaled); fl[].zero (see cdisort.h)
   FLDN(lu)......Diffuse down flux (delta-m scaled); fl[].one (see cdisort.h)
   FLYR(lc)......Truncated fraction in delta-m method
   KK(lc)........Eigenvalues in eq. KST(20)
   LAYRU(lu).....Computational layer in which user output level UTAU(lu) is located
   LL(iq,lc).....Constants of integration C-tilde in eqs. KST(42-43) obtained by solving eqs. KST(38-41)
   lyrcut........True, radiation is assumed zero below layer -ncut- because of almost complete absorption
   ncut..........Computational layer number in which absorption optical depth first exceeds abscut
   OPRIM(lc).....Single scattering albedo after delta-m scaling
   pass1.........TRUE on first entry, FALSE thereafter
   PKAG(0:lc)....Integrated Planck function for internal emission at layer boundaries
   PKAGC(lc).....Integrated Planck function for internal emission at layer center
   RR(lc)........Eigenvectors at polar quadrature angles.
   TAUC(0:lc)....Cumulative optical depth (un-delta-m-scaled)
   TAUCPR(0:lc)..Cumulative optical depth (delta-m-scaled if deltam = TRUE, otherwise equal to TAUC)
   tplanck.......Intensity emitted from top boundary
   U0C(iq,lu)....Azimuthally-averaged intensity
   UTAUPR(lu)....Optical depths of user output levels in delta-m coordinates;  equal to UTAU(lu) if no delta-m

   The following are members of the structure twostr_xyz:
   XB_0D(lc).....x-sub-zero-sup-minus in expansion of pseudo-spherical beam source, eq. KST(22)
   XB_0U(lc).....x-sub-zero-sup-plus  in expansion of pseudo-spherical beam source, eq. KST(22)
   XB_1D(lc).....x-sub-one-sup-minus  in expansion of pseudo-spherical beam source, eq. KST(22)
   XB_1U(lc).....x-sub-one-sup-plus   in expansion of pseudo-spherical beam source, eq. KST(22)
   XP_0(lc)......x-sub-zero in expansion of thermal source function; see eq. KST(22) (has no (mu) dependence)
   XP_1(lc)......x-sub-one  in expansion of thermal source function; see eq. KST(22) (has no (mu) dependence)
   YB_0D(lc).....y-sub-zero-sup-minus in eq. KST(23), solution for pseudo-spherical beam source
   YB_0U(lc).....y-sub-zero-sup-plus  in eq. KST(23), solution for pseudo-spherical beam source
   YB_1D(lc).....y-sub-one-sup-minus  in eq. KST(23), solution for pseudo-spherical beam source
   YB_1U(lc).....y-sub-one-sup-plus   in eq. KST(23), solution for pseudo-spherical beam source
   YP_0D(lc).....y-sub-zero-sup-minus in eq. KST(23), solution for thermal source
   YP_0U(lc).....y-sub-zero-sup-plus  in eq. KST(23), solution for thermal source
   YP_1D(lc).....y-sub-one-sup-minus  in eq. KST(23), solution for thermal source
   YP_1U(lc).....y-sub-one-sup-plus   in eq. KST(23), solution for thermal source
   ZB_A(lc)......Alpha coefficient in eq. KST(22) for pseudo-spherical beam source
   ZP_A(lc)......Alpha coefficient in eq. KST(22) for thermal source
*/

void c_twostr(disort_state  *ds,
              disort_output *out,
              int            deltam,
              double        *gg,
              int           *ierror,
              double         radius)
{
  register int
    lc,ierr;
  int
    lyrcut,iret,ncut,nn;
  int
    ipvt[ds->nstr*ds->nlyr],
    layru[ds->ntau];
  double
    cmu,bplanck,tplanck;
  double
    *b,*cband,*ch,*chtau,*dtaucpr,*expbea,*flyr,*ggprim,
    *kk,*ll,*oprim,*pkag,*pkagc,*rr,*tauc,*taucpr,*u0c,*utaupr;
  disort_pair
    *fl;
  twostr_xyz
    *ts;
  twostr_diag
    *diag;
  const double
    dither = 100.*DBL_EPSILON;

  /*
   * Allocate zeroed memory
   */
  b       = c_dbl_vector(0,ds->nstr*ds->nlyr-1,"b");
  cband   = c_dbl_vector(0,ds->nstr*ds->nlyr*(9*(ds->nstr/2)-2)-1,"cband");
  ch      = c_dbl_vector(0,ds->nlyr-1,"ch");
  chtau   = c_dbl_vector(0,(2*ds->nlyr+1)-1,"chtau");
  dtaucpr = c_dbl_vector(0,ds->nlyr-1,"dtaucpr");
  expbea  = c_dbl_vector(0,ds->nlyr,"expbea");
  flyr    = c_dbl_vector(0,ds->nlyr-1,"flyr");
  ggprim  = c_dbl_vector(0,ds->nlyr-1,"ggprim");
  kk      = c_dbl_vector(0,ds->nlyr-1,"kk");
  ll      = c_dbl_vector(0,ds->nlyr*ds->nstr-1,"ll");
  oprim   = c_dbl_vector(0,ds->nlyr-1,"oprim");
  pkag    = c_dbl_vector(0,ds->nlyr,"pkag");
  pkagc   = c_dbl_vector(0,ds->nlyr-1,"pkagc");
  rr      = c_dbl_vector(0,ds->nlyr-1,"rr");
  tauc    = c_dbl_vector(0,ds->nlyr,"tauc");
  taucpr  = c_dbl_vector(0,ds->nlyr,"taucpr");
  u0c     = c_dbl_vector(0,ds->ntau*ds->nstr-1,"u0c");
  utaupr  = c_dbl_vector(0,ds->ntau-1,"utaupr");
  /*
   * Using C structures to facilitate cache-aware memory allocation, which tends to
   * reduce cache misses and speed up computer execution.
   */
  fl   = (disort_pair *)calloc(ds->ntau,  sizeof(disort_pair)); if (!fl)   c_errmsg("twostr alloc error for fl",  DS_ERROR);
  ts   = (twostr_xyz  *)calloc(ds->nlyr,  sizeof(twostr_xyz )); if (!ts)   c_errmsg("twostr alloc error for ts",  DS_ERROR);
  diag = (twostr_diag *)calloc(2*ds->nlyr,sizeof(twostr_diag)); if (!diag) c_errmsg("twostr alloc error for diag",DS_ERROR);

  if(ds->flag.prnt[0]) {
    fprintf(stdout,"\n\n\n\n"
            " ************************************************************************************************************************\n"
            "                         Two stream method radiative transfer program, version 1.13\n"
            " ************************************************************************************************************************\n");
  }

  memset(ierror,0,TWOSTR_NERR*sizeof(int));

  /*
   * Calculate cumulative optical depth and dither single-scatter albedo to improve numerical behavior of
   * eigenvalue/vector computation
   */

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if(SSALB(lc) == 1.) {
      SSALB(lc) = 1.-dither;
    }
    TAUC(lc) = TAUC(lc-1)+DTAUC(lc);
  }

  /*
   * Check input dimensions and variables
   */
  c_twostr_check_inputs(ds,gg,ierror,tauc);

  iret = 0;
  for (ierr = 1; ierr <= TWOSTR_NERR; ierr++) {
    if (IERROR(ierr) != 0) {
      iret = 1;
      if (ds->flag.quiet==VERBOSE) {
        fprintf(stderr,"\ntwostr reports fatal error: %d\n",ierr);
      }
    }
  }
  if (iret == 1) {
    goto free_local_memory_and_return;
  }

 /*
  * Perform various setup operations
  */
  c_twostr_set(ds,&bplanck,ch,chtau,&cmu,deltam,dtaucpr,expbea,flyr,gg,ggprim,layru,&lyrcut,
             &ncut,&nn,oprim,pkag,pkagc,radius,tauc,taucpr,&tplanck,utaupr);

  /*
   * Print input information
   */
  if (ds->flag.prnt[0]) {
    c_twostr_print_inputs(ds,deltam,flyr,gg,lyrcut,oprim,tauc,taucpr);
  }

  /*
   * Calculate the homogenous and particular solutions
   */
  c_twostr_solns(ds,ch,chtau,cmu,ncut,oprim,pkag,pkagc,taucpr,ggprim,kk,rr,ts);

  /*
   * Solve for constants of integration in homogeneous solution (general boundary conditions)
   */
  c_twostr_solve_bc(ds,ts,bplanck,cband,cmu,expbea,lyrcut,nn,ncut,tplanck,taucpr,kk,rr,ipvt,b,ll,diag);

  /*
   * Compute upward and downward fluxes, mean intensities and flux divergences.
   */
  c_twostr_fluxes(ds,ts,ch,cmu,kk,layru,ll,lyrcut,ncut,oprim,rr,taucpr,utaupr,out,u0c,fl);

  /*
   * Free allocated memory
   */
 free_local_memory_and_return:
  free(b),     free(cband), free(ch),  free(chtau),free(dtaucpr),free(expbea),free(flyr),free(fl);
  free(ggprim),free(kk),    free(ll),  free(oprim),free(pkag),   free(pkagc), free(rr),  free(tauc),free(taucpr);
  free(u0c),   free(utaupr),free(diag),free(ts);
  
  return;
}

/*============================= end of c_twostr() ========================*/

/*============================= c_chapman() ==============================*/

/*
 Calculates the Chapman factor.

 I n p u t       v a r i a b l e s:

      lc        : Computational layer
      taup      :
      tauc      :
      nlyr      : Number of layers in atmospheric model
      zd(lc)    : lc = 0, nlyr. zd(lc) is distance from bottom
                  surface to top of layer lc. zd(nlyr) = 0. km
      dtau_c    : Optical thickness of layer lc (un-delta-m-scaled)
      zenang    : Solar zenith angle as seen from bottom surface
      r         : Radial parameter, see Velinow & Kostov (2001). NOTE: Use the same dimension as zd,
                  for instance both in km.

 O u t p u t      v a r i a b l e s:

      ch        : Chapman-factor. In a pseudo-spherical atmosphere, replace exp(-tau/umu0) by exp(-ch(lc)) in the
                  beam source in

 I n t e r n a l     v a r i a b l e s:

      dhj       : delta-h-sub-j in eq. B2 (DS)
      dsj       : delta-s-sub-j in eq. B2 (DS)
      fact      : =1 for first  sum in eq. B2 (DS)
                  =2 for second sum in eq. B2 (DS)
      rj        : r-sub-j   in eq. B1 (DS)
      rjp1      : r-sub-j+1 in eq. B1 (DS)
      xpsinz    : The length of the line OG in Fig. 1, (DS)

 
 NOTE: Assumes a spherical planet. One might consider generalizing following
       Velinow YPI, Kostov VI, 2001, Generalization on Chapman Function for the Atmosphere of an Oblate Rotating Planet, 
         Comptes Rendus de l'Academie Bulgare des Sciences 54, 29-34.
*/

double c_chapman(int     lc,
                 double  taup,
                 double *tauc,
                 int     nlyr,
                 double *zd,
                 double *dtau_c,
                 double  zenang,
                 double  r)
{
  register int
    id,j;
  double
    zenrad,xp,xpsinz,
    sum,fact,fact2,rj,rjp1,dhj,dsj;

  zenrad = zenang*DEG;
  xp     = r+ZD(lc)+(ZD(lc-1)-ZD(lc))*taup;
  xpsinz = xp*sin(zenrad);

  if (zenang > 90. && xpsinz < r) {
    return 1.e+20;
  }

  /*
   * Find index of layer in which the screening height lies
   */
  id = lc;
  if (zenang > 90.) {
    for (j= lc; j <= nlyr; j++) {
      if (xpsinz < (ZD(j-1)+r) && (xpsinz >= ZD(j)+r)) {
        id = j;
      }
    }
  }

  sum = 0.;
  for (j = 1; j <= id; j++) {
    fact  = 1.;
    fact2 = 1.;
    /*
     * Include factor of 2 for zenang > 90., second sum in eq. B2 (DS)
     */
    if (j > lc) {
      fact = 2.;
    }
    else if (j == lc && lc == id && zenang > 90.) {
      fact2 = -1.;
    }

    rj   = r+ZD(j-1);
    rjp1 = r+ZD(j  );
    if (j == lc && id == lc) {
      rjp1 = xp;
    }

    dhj = ZD(j-1)-ZD(j);
    if (id > lc && j == id) {
      dsj = sqrt(rj*rj-xpsinz*xpsinz);
    }
    else {
      dsj = sqrt(rj*rj-xpsinz*xpsinz)-fact2*sqrt(rjp1*rjp1-xpsinz*xpsinz);
    }
    sum += DTAU_C(j)*fact*dsj/dhj;
  }
  /*
   * Add third term in eq. B2 (DS)
   */
  if (id > lc) {
    dhj  = ZD(lc-1)-ZD(lc);
    dsj  = sqrt(xp*xp-xpsinz*xpsinz)-sqrt(SQR(ZD(lc)+r)-xpsinz*xpsinz);
    sum += DTAU_C(lc)*dsj/dhj;
  }

  return sum;
}

/*============================= end of c_chapman() =======================*/

double c_chapman_simpler(int     lc,
                 double  taup,
                 int     nlyr,
                 double *zd,
                 double *dtau_c,
                 double  zenang,
                 double  r)
{
  register int
    id,j;
  double
    zenrad,xp,xpsinz,
    sum,fact,fact2,rj,rjp1,dhj,dsj;

  zenrad = zenang*DEG;
  xp     = r+ZD(lc)+(ZD(lc-1)-ZD(lc))*taup;
  xpsinz = xp*sin(zenrad);

  if (zenang > 90. && xpsinz < r) {
    return 1.e+20;
  }

  /*
   * Find index of layer in which the screening height lies
   */
  id = lc;
  if (zenang > 90.) {
    for (j= lc; j <= nlyr; j++) {
      if (xpsinz < (ZD(j-1)+r) && (xpsinz >= ZD(j)+r)) {
        id = j;
      }
    }
  }

  sum = 0.;
  for (j = 1; j <= id; j++) {
    fact  = 1.;
    fact2 = 1.;
    /*
     * Include factor of 2 for zenang > 90., second sum in eq. B2 (DS)
     */
    if (j > lc) {
      fact = 2.;
    }
    else if (j == lc && lc == id && zenang > 90.) {
      fact2 = -1.;
    }

    rj   = r+ZD(j-1);
    rjp1 = r+ZD(j  );
    if (j == lc && id == lc) {
      rjp1 = xp;
    }

    dhj = ZD(j-1)-ZD(j);
    if (id > lc && j == id) {
      dsj = sqrt(rj*rj-xpsinz*xpsinz);
    }
    else {
      dsj = sqrt(rj*rj-xpsinz*xpsinz)-fact2*sqrt(rjp1*rjp1-xpsinz*xpsinz);
    }
    sum += DTAU_C(j)*fact*dsj/dhj;
  }
  /*
   * Add third term in eq. B2 (DS)
   */
  if (id > lc) {
    dhj  = ZD(lc-1)-ZD(lc);
    dsj  = sqrt(xp*xp-xpsinz*xpsinz)-sqrt(SQR(ZD(lc)+r)-xpsinz*xpsinz);
    sum += DTAU_C(lc)*dsj/dhj;
  }

  return sum;
}

/*============================= end of c_chapman_simpler() ===============*/

/*============================= c_twostr_check_inputs() ==================*/

/*
 * Checks the twostr input dimensions and variables
 */

void c_twostr_check_inputs(disort_state *ds,
                           double       *gg,
                           int          *ierror,
                           double       *tauc)
{
  int
    inperr,lc,lu;
  double
    umumin;

  inperr = FALSE;

  if (ds->nlyr < 1) {
    inperr    = c_write_bad_var(ds->flag.quiet,"nlyr");
    IERROR(1) = 1;
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (DTAUC(lc) < 0.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"dtauc");
      IERROR(3) += 1;
    }
    if (SSALB(lc) < 0. || SSALB(lc) > 1.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"ssalb");
      IERROR(4) += 1;
    }
    if (ds->flag.planck) {
      if (lc == 1 && TEMPER(0) < 0.) {
        inperr     = c_write_bad_var(ds->flag.quiet,"temper");
        IERROR(5) += 1;
      }
      if (TEMPER(lc) < 0.) {
        inperr     = c_write_bad_var(ds->flag.quiet,"temper");
        IERROR(5) += 1;
      }
    }
    if (GG(lc) < -1. || GG(lc) > 1.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"gg");
      IERROR(6) += 1;
    }
  }

  if(ds->flag.spher==TRUE) {
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (ds->ZD(lc) > ds->ZD(lc-1)) {
        inperr     = c_write_bad_var(ds->flag.quiet,"zd");
        IERROR(7) += 1;
      }
    }
  }

  if (ds->flag.usrtau) {
    if (ds->ntau < 1) {
      inperr    = c_write_bad_var(ds->flag.quiet,"ntau");
      IERROR(8) = 1;
    }
    for (lu = 1; lu <= ds->ntau; lu++) {
      if (fabs(UTAU(lu)-TAUC(ds->nlyr)) <= 1.e-6*TAUC(ds->nlyr)) { /* relative check copied from c_check_inputs() */
        UTAU(lu)= TAUC(ds->nlyr);
      }
      if (UTAU(lu) < 0. || UTAU(lu) > TAUC(ds->nlyr)) {
        inperr      = c_write_bad_var(ds->flag.quiet,"utau");
        IERROR(10) += 1;
      }
    }
  }

  if (ds->bc.fbeam < 0.) {
    inperr     = c_write_bad_var(ds->flag.quiet,"fbeam");
    IERROR(12) = 1;
  }

  umumin = 0.;
  if(ds->flag.spher==TRUE) {
    umumin = -1.;
  }

  if (ds->bc.fbeam > 0. && (ds->bc.umu0 <= umumin || ds->bc.umu0 > 1.)) {
    inperr     = c_write_bad_var(ds->flag.quiet,"umu0");
    IERROR(13) = 1;
  }
  if (ds->bc.fisot < 0.) {
    inperr     = c_write_bad_var(ds->flag.quiet,"fisot");
    IERROR(14) = 1;
  }
  if (ds->bc.albedo < 0. || ds->bc.albedo > 1.) {
    inperr     = c_write_bad_var(ds->flag.quiet,"albedo");
    IERROR(15) = 1;
  }

  if(ds->flag.planck) {
    if (ds->wvnmlo < 0. || ds->wvnmhi < ds->wvnmlo) {
      inperr     = c_write_bad_var(ds->flag.quiet,"wvnmlo,hi");
      IERROR(16) = 1;
    }
    if (ds->bc.temis < 0. || ds->bc.temis > 1.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"temis");
      IERROR(17) = 1;
    }
    if (ds->bc.btemp < 0.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"btemp");
      IERROR(18) = 1;
    }
    if (ds->bc.ttemp < 0.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"ttemp");
      IERROR(19) = 1;
    }
  }

  if (!ds->flag.usrtau && ds->ntau < ds->nlyr+1) {
    inperr = c_write_too_small_dim(ds->flag.quiet,"ds.ntau",ds->nlyr+1);
    IERROR(22) = 1;
  }

  if (ds->bc.fluor < 0.) {
    inperr     = c_write_bad_var(ds->flag.quiet,"fluor");
    IERROR(23) = 1;
  }

  if (inperr) {
    c_errmsg("twostr_check_inputs--input and/or dimension errors",DS_ERROR);
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (ds->flag.planck && fabs(TEMPER(lc)-TEMPER(lc-1)) > 50. && ds->flag.quiet==VERBOSE) {
      c_errmsg("twostr_check_inputs--vertical temperature step may be too large for good accuracy",DS_WARNING);
    }
  }

  return;
}

/*============================= end of c_twostr_check_inputs() ===========*/

/*============================= c_twostr_fluxes() ========================*/

/*
 Calculates the radiative fluxes, mean intensity, and flux derivative
 with respect to optical depth from the azimuthally-averaged intensity

 I n p u t     v a r i a b l e s:

   ds         :  'Disort' state variables
   ts         :  twostr_xyz structure variables (xp_0, yb_0d, zb_a...; see cdisort.h)
   ch         :  Chapman factor
   cmu        :  Abscissa for gauss quadrature over angle cosine
   kk         :  Eigenvalues
   layru      :  Layer numbers of user levels -utau-
   ll         :  Constants of integration in eqs. KST(42-43), obtaine by solving eqs. KST(38-41)
   lyrcut     :  Logical flag for truncation of comput. layer
   ncut       :  Number of computational layer where absorption optical depth exceeds -abscut-
   oprim      :  Delta-m scaled single scattering albedo
   rr         :  Eigenvectors at polar quadrature angles
   flag.spher :  TRUE turns on pseudo-spherical effects
   taucpr     :  Cumulative optical depth (delta-m-scaled)
   utaupr     :  Optical depths of user output levels in delta-m coordinates; equal to  -utau- if no delta-m

 O u t p u t     v a r i a b l e s:

   out      :  'Disort' output variables
   u0c      :  Azimuthally averaged intensities at polar quadrature angle cmu

 I n t e r n a l       v a r i a b l e s:

   dirint   :  direct intensity attenuated
   fdntot   :  total downward flux (direct + diffuse)
   fldir    :  fl[].zero, direct-beam flux (delta-m scaled)
   fldn     :  fl[].one, diffuse down-flux (delta-m scaled)
   fnet     :  net flux (total-down - diffuse-up)
   fact     :  EXP( - utaupr / ch ), where ch is the Chapman factor
   plsorc   :  Planck source function (thermal)
 ---------------------------------------------------------------------*/

void c_twostr_fluxes(disort_state  *ds,
                     twostr_xyz    *ts,
                     double        *ch,
                     double         cmu,
                     double        *kk,
                     int           *layru,
                     double        *ll,
                     int            lyrcut,
                     int            ncut,
                     double        *oprim,
                     double        *rr,
                     double        *taucpr,
                     double        *utaupr,
                     disort_output *out,
                     double        *u0c,
                     disort_pair   *fl)
{
  register int
    lu,lyu;
  double
    fdntot,fnet,plsorc,dirint;
  register double
    fact1,fact2;

  if (ds->flag.prnt[1]) {
    fprintf(stdout,"\n\n                     <----------------------- Fluxes ----------------------->\n"
                   "   optical  compu    downward    downward    downward       upward                    mean      Planck   d(net flux)\n"
                   "     depth  layer      direct     diffuse       total      diffuse         net   intensity      source   / d(op dep)\n");
  }

  memset(out->rad,0,ds->ntau*sizeof(disort_radiant));

  /*
   * Loop over user levels
   */
  if (ds->flag.planck) {
    for (lu = 1; lu <= ds->ntau; lu++) {
      lyu        = LAYRU(lu);
      fact1      = exp(-ZP_A(lyu)*UTAUPR(lu));
      U0C(1,lu) += fact1*(YP_0D(lyu)+YP_1D(lyu)*UTAUPR(lu));
      U0C(2,lu) += fact1*(YP_0U(lyu)+YP_1U(lyu)*UTAUPR(lu));
    }
  }
  for (lu = 1; lu <= ds->ntau; lu++) {
    lyu = LAYRU(lu);
    if (lyrcut && lyu > ncut) {
      /*
       * No radiation reaches this level
       */
      fdntot = 0.;
      fnet   = 0.;
      plsorc = 0.;
    }
    else {
      if (ds->bc.fbeam > 0.) {
        fact1      = exp(-ZB_A(lyu)*UTAUPR(lu));
        U0C(1,lu) += fact1*(YB_0D(lyu)+YB_1D(lyu)*UTAUPR(lu));
        U0C(2,lu) += fact1*(YB_0U(lyu)+YB_1U(lyu)*UTAUPR(lu));
        if (ds->bc.umu0 > 0. || ds->flag.spher) {
          fact1      = ds->bc.fbeam*exp(-UTAUPR(lu)/CH(lyu));
          dirint     = fact1;
          FLDIR(lu)  = fabs(ds->bc.umu0)*fact1;
          RFLDIR(lu) = fabs(ds->bc.umu0)*ds->bc.fbeam*exp(-UTAU(lu)/CH(lyu));
        }
        else {
          dirint     = 0.;
          FLDIR(lu)  = 0.;
          RFLDIR(lu) = 0.;
        }
      }
      else {
        dirint     = 0.;
        FLDIR(lu)  = 0.;
        RFLDIR(lu) = 0.;
      }
      fact1      = LL(1,lyu)*exp( KK(lyu)*(UTAUPR(lu)-TAUCPR(lyu  )));
      fact2      = LL(2,lyu)*exp(-KK(lyu)*(UTAUPR(lu)-TAUCPR(lyu-1)));
      U0C(1,lu) += fact2+RR(lyu)*fact1;
      U0C(2,lu) += fact1+RR(lyu)*fact2;
      /*
       * Calculate fluxes and mean intensities; downward and upward fluxes from eq. KST(9)
       */
      fact1     = 2.*M_PI*cmu;
      FLDN(lu)  = fact1*U0C(1,lu);
      FLUP(lu)  = fact1*U0C(2,lu);
      fdntot    = FLDN(lu)+FLDIR(lu);
      fnet      = fdntot-FLUP(lu);
      RFLDN(lu) = fdntot-RFLDIR(lu);
      /*
       * Mean intensity from eq. KST(10)
       */
      UAVG(lu) = U0C(1,lu)+U0C(2,lu);
      UAVG(lu) = (2.*M_PI*UAVG(lu)+dirint)/(4.*M_PI);

      /*
       * Flux divergence from eqs. KST(11-12)
       */
      plsorc   = 1./(1.-OPRIM(lyu))*exp(-ZP_A(lyu)*UTAUPR(lu))*(XP_0(lyu)+XP_1(lyu)*UTAUPR(lu));
      DFDT(lu) = (1.-SSALB(lyu))*4.*M_PI*(UAVG(lu)-plsorc);
    }
    if (ds->flag.prnt[1]) {
      fprintf(stdout,"%10.4f%7d%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%14.3e\n",
                     UTAU(lu),lyu,RFLDIR(lu),RFLDN(lu),fdntot,FLUP(lu),fnet,UAVG(lu),plsorc,DFDT(lu));
    }
  }

  return;
}

/*============================= end of c_twostr_fluxes() =================*/

/*============================= c_twostr_solns() =========================*/

/*
    Calculates the homogenous and particular solutions to the
    radiative transfer equation in the two-stream approximation,
    for each layer in the medium.

    I n p u t     v a r i a b l e s:

      ds         : 'Disort' state variables
      ch         : Chapman correction factor
      chtau      :
      cmu        : Abscissa for gauss quadrature over angle cosine
      ncut       : Number of computational layer where absorption optical depth exceeds -abscut-
      oprim      : Delta-m scaled single scattering albedo
      pkag,c     : Planck function in each layer
      flag.spher : spher = true => spherical geometry invoked
      taucpr     : Cumulative optical depth (delta-m-scaled)
      ggprim     :

   O u t p u t     v a r i a b l e s:

      kk         :  Eigenvalues
      rr         :  Eigenvectors at polar quadrature angles
      ts         :  twostr_xyz structure variables (see cdisort.h)
  ----------------------------------------------------------------------*/

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
                    twostr_xyz   *ts)
{
  register int
    lc;
  static int
    initialized = FALSE;
  static double
    big,large,small,little;
  double
    q_1,q_2,qq,q0a,q0,q1a,q2a,q1,q2,
    deltat,denomb,z0p,z0m,arg,sgn,fact3,denomp,
    beta,fact1,fact2;

  if (!initialized) {
    /*
     * The calculation of the particular solutions require some care; small,little,
       big, and large have been set so that no problems should occur in double precision.
     */
    small  = 1.e+30*DBL_MIN;
    little = 1.e+20*DBL_MIN;
    big    = sqrt(DBL_MAX)/1.e+10;
    large  = log(DBL_MAX)-20.;

    initialized = TRUE;
  }

  /*----------------  Begin loop on computational layers  ---------------------*/

  for (lc = 1; lc <= ncut; lc++) {
    /*
     * Calculate eigenvalues -kk- and eigenvector -rr-, eqs. KST(20-21)
     */
    beta   = 0.5*(1.-3.*GGPRIM(lc)*cmu*cmu);
    fact1  = 1.-OPRIM(lc);
    fact2  = 1.-OPRIM(lc)+2.*OPRIM(lc)*beta;
    KK(lc) = (1./cmu)*sqrt(fact1*fact2);
    RR(lc) = (sqrt(fact2)-sqrt(fact1))/(sqrt(fact2)+sqrt(fact1));

    if (ds->bc.fbeam > 0.) {
      /*
       * Set coefficients in KST(22) for beam source
       */
      q_1 = ds->bc.fbeam/(4.*M_PI)*OPRIM(lc)*(1.-3.*GGPRIM(lc)*cmu*ds->bc.umu0);
      q_2 = ds->bc.fbeam/(4.*M_PI)*OPRIM(lc)*(1.+3.*GGPRIM(lc)*cmu*ds->bc.umu0);

      if (ds->bc.umu0 >= 0.) {
        qq = q_2;
      }
      else {
        qq = q_1;
      }

      if (ds->flag.spher) {
        q0a = exp(-CHTAU(lc-1));
        q0  = q0a*qq;
        if (q0 <= small) {
          q1a = 0.;
          q2a = 0.;
        }
        else {
          q1a = exp(-CHTAU(lc-1  ));
          q2a = exp(-CHTAU(lc));
        }
      }
      else {
        q0a = exp(-TAUCPR(lc-1)/ds->bc.umu0);
        q0  = q0a*qq;
        if (q0 <= small) {
          q1a = 0.;
          q2a = 0.;
        }
        else {
          q1a = exp(-(TAUCPR(lc-1)+TAUCPR(lc))/(2.*ds->bc.umu0));
          q2a = exp(-TAUCPR(lc)/ds->bc.umu0);
        }
      }
      q1 = q1a*qq;
      q2 = q2a*qq;

      /*
       * Calculate alpha coefficient
       */
      deltat     = TAUCPR(lc)-TAUCPR(lc-1);
      ZB_A(lc)   = 1./CH(lc);
      if (fabs(ZB_A(lc)*TAUCPR(lc-1)) > large || fabs(ZB_A(lc)*TAUCPR(lc)) > large) {
        ZB_A(lc) = 0.;
      }

      /*
       * Dither alpha if it is close to an eigenvalue
       */
      denomb = fact1*fact2-SQR(ZB_A(lc)*cmu);
      if (denomb < 1.e-03) {
        ZB_A(lc) = 1.02*ZB_A(lc);
      }
      q0 = q0a*q_1;
      q2 = q2a*q_1;

      /*
       * Set constants in eq. KST(22)
       */
      if (deltat < 1.e-07) {
        XB_1D(lc) = 0.;
      }
      else {
        XB_1D(lc) = 1./deltat*(q2*exp(ZB_A(lc)*TAUCPR(lc))-q0*exp(ZB_A(lc)*TAUCPR(lc-1)));
      }
      XB_0D(lc) = q0*exp(ZB_A(lc)*TAUCPR(lc-1))-XB_1D(lc)*TAUCPR(lc-1);
      q0        = q0a*q_2;
      q2        = q2a*q_2;

      if (deltat < 1.e-07) {
        XB_1U(lc) = 0.;
      }
      else {
        XB_1U(lc) = 1./deltat*(q2*exp(ZB_A(lc)*TAUCPR(lc))-q0*exp(ZB_A(lc)*TAUCPR(lc-1)));
      }
      XB_0U(lc) = q0*exp(ZB_A(lc)*TAUCPR(lc-1))-XB_1U(lc)*TAUCPR(lc-1);

      /*
       * Calculate particular solutions for incident beam source in pseudo-spherical geometry, eqs. KST(24-25)
       */
      denomb    = fact1*fact2-SQR(ZB_A(lc)*cmu);
      YB_1D(lc) = (OPRIM(lc)*beta*XB_1D(lc)+(1.-OPRIM(lc)*(1.-beta)+ZB_A(lc)*cmu)*XB_1U(lc))/denomb;
      YB_1U(lc) = (OPRIM(lc)*beta*XB_1U(lc)+(1.-OPRIM(lc)*(1.-beta)-ZB_A(lc)*cmu)*XB_1D(lc))/denomb;
      z0p       = XB_0U(lc)-cmu*YB_1D(lc);
      z0m       = XB_0D(lc)+cmu*YB_1U(lc);
      YB_0D(lc) = (OPRIM(lc)*beta*z0m+(1.-OPRIM(lc)*(1.-beta)+ZB_A(lc)*cmu)*z0p)/denomb;
      YB_0U(lc) = (OPRIM(lc)*beta*z0p+(1.-OPRIM(lc)*(1.-beta)-ZB_A(lc)*cmu)*z0m)/denomb;
    }

    if(ds->flag.planck) {
      /*
       * Set coefficients in KST(22) for thermal source
       * Calculate alpha coefficient
       */
      q0     = (1.-OPRIM(lc))*PKAG(lc-1);
      q1     = (1.-OPRIM(lc))*PKAGC(lc);
      q2     = (1.-OPRIM(lc))*PKAG(lc);
      deltat = TAUCPR(lc)-TAUCPR(lc-1);

      if ((q2 < q0*1.e-02 || q2 <= little) && q1 > little && q0 > little) {
        /*
         * Case 1: source small at bottom layer; alpha eq. KS(50)
         */
        ZP_A(lc) = MIN(2./deltat*log(q0/q1),big);
        if (ZP_A(lc)*TAUCPR(lc-1) >= log(big)) {
          XP_0(lc) = big;
        }
        else {
          XP_0(lc) = q0;
        }
        XP_1(lc) = 0.;
      }
      else if ((q2 <= q1*1.e-02 || q2 <= little) && (q1 <= q0*1.e-02 || q1 <= little) && q0 > little) {
        /*
         * Case 2: Source small at center and bottom of layer
         */
        ZP_A(lc) = big/TAUCPR(ncut);
        XP_0(lc) = q0;
        XP_1(lc) = 0.;
      }
      else if (q2 <= little && q1 <= little && q0 <= little) {
        /*
         * Case 3: All sources zero
         */
        ZP_A(lc) = 0.;
        XP_0(lc) = 0.;
        XP_1(lc) = 0.;
      }
      else if ( ( fabs((q2-q0)/q2) < 1.e-04 && fabs((q2-q1)/q2) < 1.e-04 ) || deltat < 1.e-04) {
        /*
         * Case 4: Sources same at center, bottom and top of layer or layer optically very thin
         */
        ZP_A(lc) = 0.;
        XP_0(lc) = q0;
        XP_1(lc) = 0.;
      }
      else {
        /*
         *  Case 5: Normal case
         */
        arg = MAX(SQR(q1/q2)-q0/q2,0.);
        /*
         * alpha eq. (44). For source that has its maximum at the top of the layer, use negative solution
         */
        sgn = 1.;
        if (PKAG(lc-1) > PKAG(lc)) {
         sgn = -1.;
        }
        fact3 = log(q1/q2+sgn*sqrt(arg));

        /* Be careful with log of numbers close to one */
        if (fabs(fact3) <= 0.005) {
          /* numbers close to one */
          q1    = 0.99*q1; 
          fact3 = log(q1/q2+sgn*sqrt(arg));
        }

        ZP_A(lc) = 2./deltat*fact3;
        if (fabs(ZP_A(lc)*TAUCPR(lc)) > log(DBL_MAX)-log(q0*100.)) {
          ZP_A(lc) = 0.;
        }

        /*
         * Dither alpha if it is close to an eigenvalue
         */
        denomp = fact1*fact2-SQR(ZP_A(lc)*cmu);
        if (denomp < 1.e-03) {
          ZP_A(lc) *= 1.01;
        }

        /*
         * Set constants in eqs. KST(22)
         */
        if(deltat < 1.e-07) {
          XP_1(lc) = 0.;
        }
        else {
          XP_1(lc) = 1./deltat*(q2*exp(ZP_A(lc)*TAUCPR(lc))-q0*exp(ZP_A(lc)*TAUCPR(lc-1)));
        }
        XP_0(lc) = q0*exp(ZP_A(lc)*TAUCPR(lc-1))-XP_1(lc)*TAUCPR(lc-1);
      }

      /*
       * Calculate particular solutions eqs. KST(24-25) for internal thermal so
       */
      denomp    = fact1*fact2-SQR(ZP_A(lc)*cmu);
      YP_1D(lc) = (OPRIM(lc)*beta*XP_1(lc)+(1.-OPRIM(lc)*(1.-beta)+ZP_A(lc)*cmu)*XP_1(lc))/denomp;
      YP_1U(lc) = (OPRIM(lc)*beta*XP_1(lc)+(1.-OPRIM(lc)*(1.-beta)-ZP_A(lc)*cmu)*XP_1(lc))/denomp;
      z0p       = XP_0(lc)-cmu*YP_1D(lc);
      z0m       = XP_0(lc)+cmu*YP_1U(lc);
      YP_0D(lc) = (OPRIM(lc)*beta*z0m+(1.-OPRIM(lc)*(1.-beta)+ZP_A(lc)*cmu)*z0p)/denomp;
      YP_0U(lc) = (OPRIM(lc)*beta*z0p+(1.-OPRIM(lc)*(1.-beta)-ZP_A(lc)*cmu)*z0m)/denomp;
    }
  }

  return;
}

/*============================= end of c_twostr_solns() ==================*/

/*============================= c_twostr_print_inputs() ==================*/

/*
 * Print values of twostream input variables
 */
void c_twostr_print_inputs(disort_state *ds,
                           int           deltam,
                           double       *flyr,
                           double       *gg,
                           int           lyrcut,
                           double       *oprim,
                           double       *tauc,
                           double       *taucpr)
{
  register int
    lu,lc;

  fprintf(stdout,"\n\n"
                 " ****************************************************************************************************\n"
                 " %s\n"
                 " ****************************************************************************************************\n",
                 ds->header);

  fprintf(stdout,"\n No. streams = %4d     No. computational layers =%4d\n",ds->nstr,ds->nlyr);
  fprintf(stdout,"%4d User optical depths :",ds->ntau);
  for (lu = 1; lu <= ds->ntau; lu++) {
    fprintf(stdout,"%10.4f",UTAU(lu));
    if (lu%10 == 0) {
      fprintf(stdout,"\n                          ");
    }
  }
  fprintf(stdout,"\n");

  if (ds->flag.spher) {
    fprintf(stdout," Pseudo-spherical geometry invoked\n");
  }

  if(!ds->flag.planck) {
    fprintf(stdout," No thermal emission\n");
  }

  fprintf(stdout,"    Incident beam with intensity =%11.3e and polar angle cosine = %8.5f\n"
                 "    plus isotropic incident intensity =%11.3e\n",
                 ds->bc.fbeam,ds->bc.umu0,ds->bc.fisot);

  fprintf(stdout,"    Bottom albedo (lambertian) =%8.4f\n",ds->bc.albedo);

  if(ds->flag.planck) {
    fprintf(stdout,"    Thermal emission in wavenumber interval :%14.4f%14.4f\n"
                   "    bottom temperature =%10.2f     top temperature =%10.2f    top emissivity =%8.4f\n",
                   ds->wvnmlo,ds->wvnmhi,ds->bc.btemp,ds->bc.ttemp,ds->bc.temis);
  }

  if(deltam) {
    fprintf(stdout," Uses delta-m method\n");
  }
  else {
    fprintf(stdout," Does not use delta-m method\n");
  }

  if(lyrcut) {
    fprintf(stdout," Sets radiation = 0 below absorption optical depth 10\n");
  }

  if(ds->flag.planck) {
    fprintf(stdout,"\n                                     <------------- delta-m --------------->"
                   "\n                   total    single                           total    single"
                   "\n       optical   optical   scatter   truncated   optical   optical   scatter    asymm"
                   "\n         depth     depth    albedo    fraction     depth     depth    albedo   factor   temperature\n");
  }
  else {
    fprintf(stdout,"\n                                     <------------- delta-m --------------->"
                   "\n                   total    single                           total    single"
                   "\n       optical   optical   scatter   truncated   optical   optical   scatter    asymm"
                   "\n         depth     depth    albedo    fraction     depth     depth    albedo   factor\n");
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (ds->flag.planck) {
      fprintf(stdout,"%4d%10.4f%10.4f%10.5f%12.5f%10.4f%10.4f%10.5f%9.4f%14.3f\n",
                     lc,DTAUC(lc),TAUC(lc),SSALB(lc),FLYR(lc),TAUCPR(lc)-TAUCPR(lc-1),TAUCPR(lc),OPRIM(lc),GG(lc),TEMPER(lc-1));
    }
    else {
      fprintf(stdout,"%4d%10.4f%10.4f%10.5f%12.5f%10.4f%10.4f%10.5f%9.4f\n",
                     lc,DTAUC(lc),TAUC(lc),SSALB(lc),FLYR(lc),TAUCPR(lc)-TAUCPR(lc-1),TAUCPR(lc),OPRIM(lc),GG(lc));
    }
  }

  if(ds->flag.planck) {
    fprintf(stdout,"                                                                                     %14.3f\n",TEMPER(ds->nlyr));
  }

  return;
}

/*============================= end of c_twostr_print_inputs() ===========*/

/*============================= c_twostr_set() ===========================*/

/*
 Perform miscellaneous setting-up operations

 Routines called: c_errmsg

 Input :  ds         'Disort' input variables

 Output:  ntau,utau  If ds->flag.usrtau = FALSE
          bplanck    Intensity emitted from bottom boundary
          ch         The Chapman factor
          cmu        Computational polar angle
          expbea     Transmission of direct beam
          flyr       Truncated fraction in delta-m method
          layru      Computational layer in which utau falls
          lyrcut     Flag as to whether radiation will be zeroed below layer ncut
          ncut       Computational layer where absorption optical depth first exceeds abscut
          nn         nstr/2 = 1
          nstr       No.of streams (=2)
          oprim      Delta-m-scaled single-scatter albedo
          pkag,c     Planck function in each layer
          taucpr     Delta-m-scaled optical depth
          tplanck    Intensity emitted from top boundary
          utaupr     Delta-m-scaled version of utau

 Internal Variables
          abscut     Absorption optical depth, medium is cut off below this depth
          tempc      Temperature at center of layer, assumed to be average of
                     layer boundary temperatures
  ---------------------------------------------------------------------*/

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
                  double       *utaupr)
{
  static int
    firstpass = TRUE;
  register int
    lc,lu,lev;
  double
    zenang,abstau,chtau_tmp,f,tempc,taup,
    abscut = 10.;

  if (firstpass) {
    firstpass = FALSE;
    ds->nstr  = 2;
    *nn       = ds->nstr/2;
  }

  if (!ds->flag.usrtau) {
    /*
     * Set output levels at computational layer boundaries
     */
    ds->ntau = ds->nlyr+1;
    for (lc = 0; lc <= ds->ntau-1; lc++) {
      UTAU(lc+1) = TAUC(lc);
    }
  }
  /*
   * Apply delta-m scaling and move description of computational layers to local variables
   */

  /*
   * NOTE: If not using calloc() to dynamically allocate memory, then need to zero-out
   *       taucpr, expbea, flyr, oprim here.
   */

  abstau = 0.;
  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (abstau < abscut) {
      *ncut = lc;
    }
    abstau += (1.-SSALB(lc))*DTAUC(lc);
    if (!deltam) {
      OPRIM(lc)   = SSALB(lc);
      TAUCPR(lc)  = TAUC(lc);
      f           = 0.;
      GGPRIM(lc)  = GG(lc);
      DTAUCPR(lc) = DTAUC(lc);
    }
    else {
     /*
      * Do delta-m transformation eqs. WW(20a,20b,14)
      */
      f           = SQR(GG(lc));
      TAUCPR(lc)  = TAUCPR(lc-1)+(1.-f*SSALB(lc))*DTAUC(lc);
      OPRIM(lc)   = SSALB(lc)*(1.-f)/(1.-f*SSALB(lc));
      GGPRIM(lc)  = (GG(lc)-f)/(1.-f);
      DTAUCPR(lc) = TAUCPR(lc)-TAUCPR(lc-1);
    }
    FLYR(lc) = f;
  }
  /*
   * If no thermal emission, cut off medium below absorption optical
   * depth = abscut (note that delta-m transformation leaves absorption
   * optical depth invariant). Not worth the trouble for one-layer problems, though.
   */
  *lyrcut = FALSE;
  if (abstau >= abscut && !ds->flag.planck && ds->nlyr > 1) {
    *lyrcut = TRUE;
  }
  if (!*lyrcut) {
    *ncut = ds->nlyr;
  }
  /*
   * Calculate Chapman function if spherical geometry, set expbea and ch for beam source.
   */
  if (ds->bc.fbeam > 0.) {
    CHTAU(0) = 0.;
    EXPBEA(0) = 1.;
    zenang    = acos(ds->bc.umu0)/DEG;
    
    if(ds->flag.spher == TRUE && ds->bc.umu0 < 0.) {
      EXPBEA(0) = exp(-c_chapman(1,0.,tauc,ds->nlyr,ds->zd,ds->dtauc,zenang,radius));
    }
    if (ds->flag.spher == TRUE) {
      for (lc = 1; lc <= *ncut; lc++) {
        taup        = TAUCPR(lc-1)+DTAUCPR(lc)/2.;
        CHTAU(lc  ) = c_chapman(lc, 0.0,      taucpr,ds->nlyr,ds->zd,dtaucpr,zenang,radius);
        chtau_tmp   = c_chapman(lc, 0.5,taucpr,ds->nlyr,ds->zd,dtaucpr,zenang,radius);
        CH(lc)      = taup/chtau_tmp;
        EXPBEA(lc)  = exp(-CHTAU(lc));
      }
    }
    else {
      for (lc = 1; lc <= *ncut; lc++) {
        CH(lc)     = ds->bc.umu0;
        EXPBEA(lc) = exp(-TAUCPR(lc)/ds->bc.umu0);
      }
    }
  }
  /*
   * Set arrays defining location of user output levels within delta-m-scaled computational mesh
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    for (lc = 1; lc <= ds->nlyr-1; lc++) {
      if (UTAU(lu) >= TAUC(lc-1) && UTAU(lu) <= TAUC(lc)) {
        break;
      }
    }
    UTAUPR(lu) = UTAU(lu);
    if (deltam) {
      UTAUPR(lu) = TAUCPR(lc-1)+(1.-SSALB(lc)*FLYR(lc))*(UTAU(lu)-TAUC(lc-1));
    }
    LAYRU(lu) = lc;
  }

  /*
   * Set computational polar angle cosine for double gaussian
   * quadrature; cmu = 0.5, or  single gaussian quadrature; cmu = 1./sqrt(3
   * See KST for discussion of which is better for your specific applicatio
   */
  if(ds->flag.planck && ds->bc.fbeam == 0.) {
    *cmu = 0.5;
  }
  else {
    *cmu = sqrt(1./3.);
  }
  /*
   * Calculate planck functions
   */
  if (!ds->flag.planck) {
    *bplanck = 0.;
    *tplanck = 0.;
    /*
     * NOTE: If not using calloc() for dynamic memory allocation, need to zero-out
     *       pkag and pkagc here.
     */
  }
  else {
    *tplanck = c_planck_func2(ds->wvnmlo,ds->wvnmhi,ds->bc.ttemp)*ds->bc.temis;
    *bplanck = c_planck_func2(ds->wvnmlo,ds->wvnmhi,ds->bc.btemp);
    for (lev = 0; lev <= ds->nlyr; lev++) {
      PKAG(lev) = c_planck_func2(ds->wvnmlo,ds->wvnmhi,TEMPER(lev));
    }
    for (lc = 1; lc <=ds->nlyr; lc++) {
      tempc     = .5*(TEMPER(lc-1)+TEMPER(lc));
      PKAGC(lc) = c_planck_func2(ds->wvnmlo,ds->wvnmhi,tempc);
    }
  }

  return;
}

/*============================= end of c_twostr_set() ====================*/

/*============================= c_twostr_solve_bc() ======================*/

/*
 Construct right-hand side vector -b- for general boundary conditions
 and solve system of equations obtained from the boundary conditions
 and the continuity-of-intensity-at-layer-interface equations.

 Routines called: c_sgbfa, c_sgbsl

 I n p u t      v a r i a b l e s:

       ds       : 'Disort' state variables
       ts       :  twostr_xyz structure variables (see cdisort.h)
       bplanck  :  Bottom boundary thermal emission
       cband    :  Left-hand side matrix of linear system eqs. KST(38-41)
                   in banded form required by linpack solution routines
       cmu      :  Abscissa for gauss quadrature over angle cosine
       expbea   :  Transmission of incident beam, EXP(-taucpr/ch)
       lyrcut   :  Logical flag for truncation of comput. layer
       ncol     :  Counts of columns in -cband-
       nn       :  Order of double-gauss quadrature (nstr/2)
       ncut     :  Total number of computational layers considered
       tplanck  :  Top boundary thermal emission
       taucpr   :  Cumulative optical depth (delta-m-scaled)
       kk       :
       rr       :
       ipvt     :

 O u t p u t     v a r i a b l e s:

       b        :  Right-hand side vector of eqs. KST(38-41) going into
                   sgbsl; returns as solution vector of eqs. KST(38-41)
                   constants of integration
       ll       :  Permanent storage for -b-, but re-ordered

 I n t e r n a l    v a r i a b l e s:

       diag     : diag[].super, diag[].on, diag[].sub

 ---------------------------------------------------------------------*/

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
                       twostr_diag  *diag)
{
  int
    info;
  register int
    irow,lc,nloop,nrow,job;
  double
    wk0,wk1,wk,rpp1_m,rp_m,rpp1_p,rp_p,sum,refflx;
  register double
    fact1,fact2,fact3,fact4;

  /*
   * First top row, top boundary condition
   */
  irow = 1;
  lc   = 1;
  /*
   * SUBD(irow) is undefined
   */
  DIAG(irow)   = RR(lc)*exp(-KK(lc)*TAUCPR(lc));
  SUPERD(irow) = 1.;
  /*
   * next from layer no. 2 to nlyr-1
   */
  nloop = ncut-1;
  for (lc = 1; lc <= nloop; lc++) {
    irow++;
    wk0          = exp(-KK(lc  )*(TAUCPR(lc  )-TAUCPR(lc-1)));
    wk1          = exp(-KK(lc+1)*(TAUCPR(lc+1)-TAUCPR(lc  )));
    SUBD(irow)   = 1.-RR(lc)*RR(lc+1);
    DIAG(irow)   = (RR(lc)-RR(lc+1))*wk0;
    SUPERD(irow) = -(1.-SQR(RR(lc+1)))*wk1;
    irow++;
    SUBD(irow)   = (1.-SQR(RR(lc)))*wk0;
    DIAG(irow)   = (RR(lc)-RR(lc+1))*wk1;
    SUPERD(irow) = -(1.-RR(lc+1)*RR(lc));
  }
  /*
   * bottom layer
   */
  irow++;
  lc = ncut;
  /*
   * SUPERD(irow) = undefined
   */
  wk = exp(-KK(lc)*(TAUCPR(lc)-TAUCPR(lc-1)));
  if (lyrcut) {
    SUBD(irow) = 1.;
    DIAG(irow) = RR(lc)*wk;
  }
  else {
    SUBD(irow) = 1.-2.*ds->bc.albedo*cmu*RR(lc);
    DIAG(irow) = (RR(lc)-2.*ds->bc.albedo*cmu)*wk;
  }

  /*
   * NOTE: If not allocating memory with calloc(), need to zero out b here.
   */

  /*
   * Construct -b-, for parallel beam + bottom reflection + thermal emission at top and/or bottom
   * 
   * Top boundary, right-hand-side of eq. KST(28)
   */
  lc   = 1;
  irow = 1;
  B(irow) = -YB_0D(lc)-YP_0D(lc)+ds->bc.fisot+tplanck;
  /*
   * Continuity condition for layer interfaces, right-hand-side of eq. KST(29)
   */
  for (lc = 1; lc <= nloop; lc++) {
    fact1     = exp(-ZB_A(lc+1)*TAUCPR(lc));
    fact2     = exp(-ZP_A(lc+1)*TAUCPR(lc));
    fact3     = exp(-ZB_A(lc  )*TAUCPR(lc));
    fact4     = exp(-ZP_A(lc  )*TAUCPR(lc));
    rpp1_m    = fact1*(YB_0D(lc+1)+YB_1D(lc+1)*TAUCPR(lc))+fact2*(YP_0D(lc+1)+YP_1D(lc+1)*TAUCPR(lc));
    rp_m      = fact3*(YB_0D(lc  )+YB_1D(lc  )*TAUCPR(lc))+fact4*(YP_0D(lc  )+YP_1D(lc  )*TAUCPR(lc));
    rpp1_p    = fact1*(YB_0U(lc+1)+YB_1U(lc+1)*TAUCPR(lc))+fact2*(YP_0U(lc+1)+YP_1U(lc+1)*TAUCPR(lc));
    rp_p      = fact3*(YB_0U(lc  )+YB_1U(lc  )*TAUCPR(lc))+fact4*(YP_0U(lc  )+YP_1U(lc  )*TAUCPR(lc));
    B(++irow) = rpp1_p-rp_p-RR(lc+1)*(rpp1_m-rp_m);
    B(++irow) = rpp1_m-rp_m-RR(lc  )*(rpp1_p-rp_p);
  }
  /*
   * Bottom boundary
   */
  lc = ncut;
  if (lyrcut) {
    /*
     * Right-hand-side of eq. KST(30)
     */
    B(++irow) = -exp(-ZB_A(ncut)*TAUCPR(ncut))*(YB_0U(ncut)+YB_1U(ncut)*TAUCPR(ncut))
                -exp(-ZP_A(ncut)*TAUCPR(ncut))*(YP_0U(ncut)+YP_1U(ncut)*TAUCPR(ncut));
  }
  else {
    sum = cmu*ds->bc.albedo*(exp(-ZB_A(ncut)*TAUCPR(ncut))*(YB_0D(ncut)+YB_1D(ncut)*TAUCPR(ncut))
                            +exp(-ZP_A(ncut)*TAUCPR(ncut))*(YP_0D(ncut)+YP_1D(ncut)*TAUCPR(ncut)));
   if (ds->bc.umu0 <= 0.) {
     refflx = 0.;
   }
   else {
     refflx = 1.;
   }
   B(++irow) = 2.*sum+ds->bc.albedo*ds->bc.umu0*ds->bc.fbeam/M_PI*refflx*EXPBEA(ncut)+(1.-ds->bc.albedo)*bplanck
               -exp(-ZB_A(ncut)*TAUCPR(ncut))*(YB_0U(ncut)+YB_1U(ncut)*TAUCPR(ncut))
               -exp(-ZP_A(ncut)*TAUCPR(ncut))*(YP_0U(ncut)+YP_1U(ncut)*TAUCPR(ncut));

 }
 /*
  * solve for constants of integration by inverting matrix KST(38-41)
  */
  nrow = irow;

  /*
   * NOTE: If not allocating memory with calloc(), need to zero out cband here.
   */

  for (irow = 1; irow <= nrow; irow++) {
    CBAND(1,irow) = 0.;
    CBAND(3,irow) = DIAG(irow);
  }
  for (irow = 1; irow <= nrow-1; irow++) {
    CBAND(2,irow+1) = SUPERD(irow);
  }
  for (irow = 2; irow <= nrow; irow++) {
    CBAND(4,irow-1) = SUBD(irow);
  }

  c_sgbfa(cband,(9*(ds->nstr/2)-2),nrow,1,1,ipvt,&info);
  job = 0;
  c_sgbsl(cband,(9*(ds->nstr/2)-2),nrow,1,1,ipvt,b,job);

  /*
   * unpack
   */
  irow = 0;
  for (lc = 1; lc <= ncut; lc++) {
    /* downward direction */
    LL(1,lc) = B(++irow);

    /* upward direction */
    LL(2,lc) = B(++irow);
  }

  return;
}

/*============================= end of c_twostr_solve_bc() ===============*/

/*============================= c_planck_func2() =========================*/

/*
  Computes Planck function integrated between two wavenumbers,
  except if wnmulo = wnmuhi, then the Planck function at wnumlo is returned

  I N P U T :  wnumlo : Lower wavenumber [inv cm] of spectral interval
               wnumhi : Upper wavenumber
               t      : Temperature [K]

  O U T P U T :  ans  : Integrated Planck function [Watts/sq m]
                         = integral (wnumlo to wnumhi) of 2h c*c nu*nu*nu/(exp(hc nu/(kT))-1), 
                         where h = Plancks constant, c = speed of light, nu = wavenumber,
                         T=temperature,and k = Boltzmann constant

  REFERENCE : Specifications of the physical world: New value of the fundamental constants,
                Dimensions/N.B.S., Jan. 1974

  METHOD :  For  -wnumlo-  close to  -wnumhi-, a Simpson-rule quadrature is done
            to avoid ill-conditioning; otherwise

            (1)  For wavenumber (wnumlo or wnumhi) small, integral(0 to wnum) is calculated by expanding
                 the integrand in a power series and integrating term by term;

            (2)  Otherwise, integral(wnumlo/hi to infinity) is calculated by expanding the denominator of the
                 integrand in powers of the exponential and integrating term by term.

  ACCURACY :  At least 6 significant digits, assuming the physical constants are infinitely accurate

  ERRORS that are not trapped:

      * Power or exponential series may underflow, giving no significant digits.
        This may or may not be of concern, depending on the application.

      * Simpson-rule special case is skipped when denominator of integrand will cause overflow.  
        In that case the normal procedure is used, which may be inaccurate if the wavenumber limits
        (wnumlo, wnumhi) are close together.
 ----------------------------------------------------------------------

        LOCAL VARIABLES

        a1,2,... :  Power series coefficients
        c2       :  h*c/k, in units cm*k (h = Planck's constant, c = speed of light, k = Boltzmann constant)
        D(I)     :  Exponential series expansion of integral of Planck function from wnumlo (i=1) 
                    or wnumhi (i=2) to infinity
        ex       :  exp(-V(I))
        exm      :  pow(ex,m)
        mmax     :  No. of terms to take in exponential series
        mv       :  multiples of 'V(i)'
        P(I)     :  Power series expansion of integral of Planck function from zero to wnumlo (i=1) or wnumhi (i=2)
        sigma    :  Stefan-Boltzmann constant (W m-2 K-4)
        sigdpi   :  sigma/pi
        smallv   :  Number of times the power series is used (0,1,2)
        V(I)     :  c2*(wnumlo(i=1) or wnumhi(i=2))/temperature
        vcut     :  Power-series cutoff point
        vcp      :  Exponential series cutoff points
        vmax     :  Largest allowable argument of 'exp' function
  ----------------------------------------------------------------------*/

#define A1    (1./3.)
#define A2    (-1./8.)
#define A3    (1./60.)
#define A4    (-1./5040.)
#define A5    (1./272160.)
#define A6    (-1./13305600.)
#define C2    (1.438786)
#define SIGMA (5.67032e-8)
#define VCUT  (1.5)
#define PLKF(x) ({const double _x = (x); _x*_x*_x/(exp(_x)-1.);})

double c_planck_func2(double wnumlo,
                      double wnumhi,
                      double t)
{
  register int
    m,n,smallv,k,i,mmax;
  static int
    initialized = FALSE;
  double
    ans,del,val,val0,oldval,exm,
    ex,mv,vsq,wvn,arg,hh,
    d[2],p[2],v[2];
  const double
    vcp[7] = {10.25,5.7,3.9,2.9,2.3,1.9,0.0};
  static double
    sigdpi,vmax,conc,c1;

  if (!initialized) {
    sigdpi = SIGMA/M_PI;
    vmax   = log(DBL_MAX);
    conc   = 15./pow(M_PI,4.);
    c1     = 1.1911e-18;  

    initialized = TRUE;
  }
  if (t < 0. || wnumhi < wnumlo || wnumlo < 0.) {
    c_errmsg("planck_func2--temperature or wavenumbers wrong",DS_ERROR);
  }
  if (t < 1.e-4) {
    return 0.;
  }
  if (wnumhi == wnumlo) {
    wvn    = wnumhi;
    arg    = exp(-C2*wvn/t);
    return c1*wvn*wvn*wvn*arg/(1.-arg);
  }

  v[0] = C2*wnumlo/t;
  v[1] = C2*wnumhi/t;

  if (v[0] > DBL_EPSILON && v[1] < vmax && (wnumhi-wnumlo)/wnumhi < 1.e-2) {
    /*
     * Wavenumbers are very close. Get integral by iterating Simpson rule to convergence.
     */
    hh     = v[1]-v[0];
    oldval = 0.;
    val0   = PLKF(v[0])+PLKF(v[1]);
    for (n = 1; n <= 10; n++) {
      del = hh/(2*n);
      val = val0;
      for (k = 1; k <=2*n-1; k++) {
        val += (double)(2*(1+k%2))*PLKF(v[0]+(double)k*del);
      }
      val *= del*A1;
      if (fabs((val-oldval)/val) <= 1.e-6) {
        return sigdpi*SQR(t*t)*conc*val;
      }
      oldval = val;
    }
    c_errmsg("planck_func2--Simpson rule did not converge",DS_WARNING);
    return sigdpi*SQR(t*t)*conc*val;
  }

  smallv = 0;
  for (i = 0; i <= 1; i++) {
    if(v[i] < VCUT) {
      /*
       * Use power series
       */
      smallv++;
      vsq  = v[i]*v[i];
      p[i] = conc*vsq*v[i]*(A1+v[i]*(A2+v[i]*(A3+vsq*(A4+vsq*(A5+vsq*A6)))));
    }
    else {
      /*
       * Use exponential series
       *
       * Find upper limit of series
       */ 
      mmax = 1;
      while (v[i] < vcp[mmax-1]) {
        mmax++;
      }
      
      ex   = exp(-v[i]);
      exm  = 1.;
      d[i] = 0.;

      for (m = 1; m <= mmax; m++) {
        mv    = (double)m*v[i];
        exm  *= ex;
        d[i] += exm*(6.+mv*(6.+mv*(3.+mv)))/SQR(m*m);
      }
      d[i] *= conc;
    }
  }

  if (smallv == 2) {
    /*
     * wnumlo and wnumhi both small
     */
    ans = p[1]-p[0];
  }
  else if (smallv == 1) {
    /*
     * wnumlo small, wnumhi large
     */
    ans = 1.-p[0]-d[1];
  }
  else {
    /*
     * wnumlo and wnumhi both large
     */
    ans = d[0]-d[1];
  }
  ans *= sigdpi*SQR(t*t);
  if (ans == 0.) {
    c_errmsg("planck_func2--returns zero; possible underflow",DS_WARNING);
  }

  return ans;
}

#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6
#undef C2
#undef SIGMA
#undef VCUT
#undef PLKF

/*============================= end of c_planck_func2() =================*/

/*============================= c_disort_state_alloc() ==================*/

/*
 * Dynamically allocate memory for disort input arrays, including
 * ones that the user can optionally ask disort() to calculate.
 */
void c_disort_state_alloc(disort_state *ds)
{
  int
    nu=0;

  ds->dtauc = c_dbl_vector(0,ds->nlyr,"ds->dtauc");  
  ds->ssalb = c_dbl_vector(0,ds->nlyr,"ds->ssalb");
  /*
   * NOTE: PMOM is used in the code even when ds->nmom is not set by the user
   *       (such as when there is no scattering). Its first dimension needs to be
   *       at least 0:ds->nstr, hence we introduce ds->nmom_nstr in the C version.
   */
  ds->nmom_nstr = IMAX(ds->nmom,ds->nstr);
  ds->pmom      = c_dbl_vector(0,(ds->nmom_nstr+1)*ds->nlyr-1,"ds->pmom");

  if (ds->flag.ibcnd == SPECIAL_BC) {
    ds->flag.planck = FALSE;
    ds->flag.lamber = TRUE;
    ds->flag.usrtau = FALSE;
  }

  /* range 0 to nlyr */
  if (ds->flag.planck == TRUE) {
    ds->temper = c_dbl_vector(0,ds->nlyr,"ds->temper");
  }
  else {
    ds->temper = NULL;
  }

  if (ds->flag.general_source == TRUE) {
    ds->gensrc  = c_dbl_vector(0,ds->nstr*ds->nlyr*ds->nstr,"ds->gensrc");
    ds->gensrcu = c_dbl_vector(0,ds->nstr*ds->nlyr*ds->numu,"ds->gensrcu");
  }
  else {
    ds->gensrc  = NULL;
    ds->gensrcu = NULL;
  }

  if (ds->flag.usrtau == FALSE) {
    ds->ntau = ds->nlyr+1;
  }
  ds->utau = c_dbl_vector(0,ds->ntau-1,"ds->utau");

  //20130723ak Treat as in c_twostr_state_alloc. See comment there.
  //20130723ak Thanks to Tim for reporting this one.
  ds->zd        = c_dbl_vector(0,ds->nlyr+1,"ds->zd");

  /* range starts at 0 */
  nu = ds->numu;
  if ( (!ds->flag.usrang || ds->flag.onlyfl)) nu = ds->nstr;

  if (ds->flag.ibcnd == SPECIAL_BC) 
    ds->umu = c_dbl_vector(0,2*nu,"ds->umu");
  else
    ds->umu = c_dbl_vector(0,nu,"ds->umu");

  if (ds->nphi >= 1) {
    ds->phi = c_dbl_vector(0,ds->nphi-1,"ds->nphi");
  }
  else {
    ds->phi = NULL;
  }

  if (!ds->flag.old_intensity_correction) {
    if (ds->nphase >= 1) {
      ds->mu_phase = c_dbl_vector(0,ds->nphase-1,"ds->mu_phase");
      ds->phase = c_dbl_vector(0,ds->nlyr*ds->nphase-1,"ds->phase");
    }
    else {
      ds->mu_phase = NULL;
      ds->phase = NULL;
    }
  }

  switch(ds->flag.brdf_type) {
    case BRDF_RPV:
      ds->brdf.rpv = (rpv_brdf_spec *)calloc(1,sizeof(rpv_brdf_spec));
      if (!ds->brdf.rpv) {
        c_errmsg("calloc error for ds->brdf.rpv",DS_ERROR);
      }
    break;
#if HAVE_BRDF
    case BRDF_AMB:
      ds->brdf.ambrals = (ambrals_brdf_spec *)calloc(1,sizeof(ambrals_brdf_spec));
      if (!ds->brdf.ambrals) {
        c_errmsg("calloc error for ds->brdf.ambrals",DS_ERROR);
      }
    break;
    case BRDF_CAM:
      ds->brdf.cam = (cam_brdf_spec *)calloc(1,sizeof(cam_brdf_spec));
      if (!ds->brdf.cam) {
        c_errmsg("calloc error for ds->brdf.cam",DS_ERROR);
      }
    break;
#endif
    default:
      ;
    break;
  }

  return;
}

/*============================= end of c_disort_state_alloc() ============*/

/*============================= c_disort_state_free() ===================*/

/*
 *  Free memory allocated by disort_state_alloc()
 */
void c_disort_state_free(disort_state *ds)
{
  if (ds->phi)    free(ds->phi);
  if (ds->umu)    free(ds->umu);
  if (ds->utau)   free(ds->utau);
  if (ds->temper) free(ds->temper);
  if (ds->pmom)   free(ds->pmom);
  if (ds->ssalb)  free(ds->ssalb);
  if (ds->dtauc)  free(ds->dtauc);
  if (ds->zd)     free(ds->zd);
  if (ds->flag.general_source == TRUE) {
    if (ds->gensrc)     free(ds->gensrc);
    if (ds->gensrcu)    free(ds->gensrcu);
  }
  if (!ds->flag.old_intensity_correction) {
    if (ds->nphase >= 1) {
      free(ds->mu_phase);
      free(ds->phase);
    }
  }

  switch(ds->flag.brdf_type) {
    case BRDF_RPV:
      free(ds->brdf.rpv);
    break;
#if HAVE_BRDF
    case BRDF_AMB:
      free(ds->brdf.ambrals);
    break;
    case BRDF_CAM:
      free(ds->brdf.cam);
    break;
#endif
    default:
      ;
    break;
  }

  return;
}

/*============================= end of c_disort_state_free() ============*/

/*============================= c_disort_out_alloc() ====================*/

/*
 *   Dynamically allocate memory for disort output arrays
 */
void c_disort_out_alloc(disort_state  *ds,
                        disort_output *out)
{

  int
    nu;

  out->rad = (disort_radiant *)calloc(ds->ntau,sizeof(disort_radiant));

  if (!out->rad) {
    c_errmsg("disort_out_alloc---error allocating out->rad array",DS_ERROR);
  }
  nu = ds->numu;
  if ( (!ds->flag.usrang || ds->flag.onlyfl)) {
    nu = ds->nstr;
  }
  out->uu = c_dbl_vector(0,ds->nphi*nu*ds->ntau,"out->uu");

  out->u0u = c_dbl_vector(0,ds->ntau*nu,"out->u0u");

  if ( ds->flag.output_uum )
    out->uum = c_dbl_vector(0,ds->nstr*nu*ds->ntau,"out->uum");

  if (ds->flag.ibcnd == SPECIAL_BC) {
    out->albmed = c_dbl_vector(0,ds->numu,"out->albmed");
    out->trnmed = c_dbl_vector(0,ds->numu,"out->trnmed");
  }
  else {
    out->albmed = NULL;
    out->trnmed = NULL;
  }

  return;
}

/*============================= end of c_disort_out_alloc() =============*/

/*============================= c_disort_out_free() =====================*/

/*
 * Free memory allocated by disort_out_alloc()
 */
void c_disort_out_free(disort_state  *ds,
                       disort_output *out)
{

  if (out->trnmed) free(out->trnmed);
  if (out->albmed) free(out->albmed);
  if (out->u0u)    free(out->u0u);
  if (out->uu)     free(out->uu);
  if (out->rad)    free(out->rad);
  if ( ds->flag.output_uum ) 
    if (out->uum) free (out->uum);

  return;
}

/*============================= end of c_disort_out_free() ==============*/

/*============================= c_twostr_state_alloc() ==================*/

/*
 * Dynamically allocate memory for twostr input arrays.
 */
void c_twostr_state_alloc(disort_state *ds)
{
  /* Set to two streams */
  ds->nstr = 2;

  /* Set flags not controlled by user */
  ds->flag.prnt[2] = FALSE;
  ds->flag.prnt[3] = FALSE;
  ds->flag.prnt[4] = FALSE;
  ds->flag.onlyfl  = TRUE;

  ds->dtauc = c_dbl_vector(0,ds->nlyr-1,"ds->dtauc");
  ds->ssalb = c_dbl_vector(0,ds->nlyr-1,"ds->ssalb");

  /* range 0 to nlyr */
  if (ds->flag.planck == TRUE) {
    ds->temper = c_dbl_vector(0,ds->nlyr,"ds->temper");
  }
  else {
    ds->temper = NULL;
  }

  if (ds->flag.usrtau == FALSE) {
    ds->ntau = ds->nlyr+1;
  }
  ds->utau = c_dbl_vector(0,ds->ntau-1,"ds->utau");

  //20120820ak Tim says: if spher is false
  //20120820ak during allocation, it has a seg fault later if you turn spher
  //20120820ak on and try to use that functionality.  Probably would be
  //20120820ak better to just allocate this little guy regardless of the status of spher.
  //20120820ak So I commented the following.
  //20120820ak  if (ds->flag.spher == TRUE) {
  ds->zd        = c_dbl_vector(0,ds->nlyr+1,"ds->zd");	
  //20120820ak  }

  return;
}

/*============================= endof c_twostr_state_alloc() ============*/

/*============================= c_twostr_state_free() ===================*/

/*
 *  Free memory allocated by twostr_state_alloc()
 */
void c_twostr_state_free(disort_state *ds)
{
  if (ds->utau)   free(ds->utau);
  if (ds->temper) free(ds->temper);
  if (ds->ssalb ) free(ds->ssalb);
  if (ds->dtauc ) free(ds->dtauc);
  if (ds->zd )    free(ds->zd);
  return;
}

/*============================= end of c_twostr_state_free() ============*/

/*============================= c_twostr_out_alloc() ====================*/

/*
 *   Dynamically allocate memory for twostr output arrays
 */
void c_twostr_out_alloc(disort_state  *ds,
                        disort_output *out)
{
  out->rad = (disort_radiant *)calloc(ds->ntau,sizeof(disort_radiant));
  if (!out->rad) {
    c_errmsg("disort_out_alloc---error allocating out->rad array",DS_ERROR);
  }

  return;
}

/*============================= end of c_twostr_out_alloc() =============*/

/*============================= c_twostr_out_free() =====================*/

/*
 * Free memory allocated by twostr_out_alloc()
 */
void c_twostr_out_free(disort_state  *ds,
                       disort_output *out)
{
  if (out->rad) free(out->rad);

  return;
}

/*============================= end of c_twostr_out_free() ==============*/

/*============================= c_dbl_vector() ==========================*/

/*
 * Allocates memory for a 1D double-precision array with range [nl..nh].
 *
 * NOTE: calloc() zeros the memory it allocates.
 */ 

double *c_dbl_vector(int  nl, 
		     int  nh,
		     char *name)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  double         
    *m;

  if (nh < nl) {
    fprintf(stderr,"\n\n**error:%s, variable %s, range (%d,%d)\n","dbl_vector",name,nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (double *)calloc(len_safe,sizeof(double));

  if (!m) {
    c_errmsg("dbl_vector---alloc error",DS_ERROR);
  }
  m -= nl_safe;

  return m;
}

/*============================= end of c_dbl_vector() ===================*/

/*============================= c_int_vector() ==========================*/

/*
 * Allocates memory for a 1D integer array with range [nl..nh].
 *
 * NOTE: calloc() zeros the memory it allocates.
 */ 

int *c_int_vector(int  nl, 
		  int  nh,
		  char *name)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  int         
    *m;

  if (nh < nl) {
    fprintf(stderr,"\n\n**error:%s, variable %s, range (%d,%d)\n","int_vector",name,nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (int *)calloc(len_safe,sizeof(int));

  if (!m) {
    c_errmsg("int_vector---alloc error",DS_ERROR);
  }
  m -= nl_safe;

  return m;
}

/*============================= end of c_int_vector() ===================*/

/*============================= c_free_dbl_vector() =====================*/

/*
 * Frees memory allocated by dbl_vector().
 *
 * NOTE: If the array is zero-offset, can just use free().
 * NOTE: Argument nh is not used, but kept to match dbl_vector().
 */

void c_free_dbl_vector(double *m, 
                       int     nl, 
                       int     nh)
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m      += nl_safe;
  free(m);

  return;
}

/*============================= end of c_free_dbl_vector() ==============*/

/*============================= c_setout() ==============================*/

/*-------------------------------------------------------------------
 * Copyright (C) 1994 Arve Kylling
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 1, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * To obtain a copy of the GNU General Public License write to the
 * Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
 * USA.
 *-------------------------------------------------------------------
 * Linearly interpolate to get approximate tau corresponding to 
 * altitude zout
 *
 * Input/output variables described in phodis.f
 *
 *
 * This code was translated to c from fortran by Robert Buras
 *
 */

int c_setout( float *sdtauc,
	      int    nlyr,
	      int    ntau,
	      float *sutau,
	      float *z,
	      float *zout )
{
  int itau=0, lc=0, itype=0;
  double hh=0.0;
  double *tauint=NULL;

  tauint = c_dbl_vector(0,nlyr+1,"tauint");

  if (tauint==NULL) {
    fprintf(stderr,"Error allocating tauint!\n");
    return -1;
  }

  /* */     

  TAUINT (1) = 0.0;
  for (lc=1; lc<=nlyr; lc++)
    TAUINT (lc+1) = TAUINT (lc) + SDTAUC (lc);

  itype = 2;

  for (itau=1; itau<=ntau; itau++) 
    SUTAU (itau) = c_inter( nlyr+1, itype, ZOUT (itau),
			    z, tauint, &hh );

  free(tauint);

  return 0;
}

/*============================= end of c_setout() =======================*/

/*============================= c_inter() ===============================*/

/*-------------------------------------------------------------------
 * Copyright (C) 1994 Arve Kylling
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 1, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * To obtain a copy of the GNU General Public License write to the
 * Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
 * USA.
 *-------------------------------------------------------------------
 *
 *     Interpolates at the x-point arg from x-value array xarr and
 *     y-value array yarr. xarr and yarr are expected to have
 *     descending arguments, i.e. for atmospheric applications
 *     xarr typically holds the altitude and inter expects
 *     xarr(1) = altitude at top of atmosphere.
 *
 *     Input variables:
 *     dim       Array dimension of xarr and yarr
 *     npoints   No. points in arrays xarr and yarr
 *     itype     Interpolation type
 *     arg       Interpolation argument
 *     xarr      array of x values
 *     yarr      array of y values
 *
 *     Output variables:
 *     ynew      Interpolated function value at arg
 *     hh        gradient or scale height value  
 *
 * This code was translated to c from fortran by Robert Buras
 *
 */

double c_inter( int     npoints,
		int     itype,
		double  arg,
		float  *xarr,
		double *yarr,
		double *hh )
{
  int iq=0, ip=0;
  double ynew=0.0;

  if ( arg <= XARR (1) && arg >= XARR (npoints) ) {
    for (iq=1;iq<=npoints-1;iq++)
      if ( arg <= XARR (iq) && arg >= XARR (iq+1) )
	ip=iq;
    if ( arg == XARR (npoints) )
      ip = npoints - 1;
  }
  else {
    if ( arg > XARR (1) )
      ip = 1;
    else {
      if ( arg < XARR (npoints) )
	ip = npoints - 1;
    }
  }

  /* Interpolate function value at arg from data points ip to ip+1 */

  switch(itype) {
  case 1:
    /*     exponential interpolation */
    if ( YARR (ip+1) == YARR (ip) ) {
      *hh = 0.0;
      ynew = YARR (ip);
    }
    else {
      *hh = -( XARR (ip+1) - XARR (ip) ) / 
	log( YARR (ip+1) / YARR (ip));
      ynew = YARR (ip) * exp(- ( arg - XARR (ip) ) / *hh );
    }
    break;
  case 2:
    /*     linear interpolation */
    *hh = ( YARR (ip+1) - YARR (ip) ) / ( XARR (ip+1) - XARR (ip) );
    ynew = YARR (ip) + *hh * ( arg - XARR (ip) );
    break;
  default:
    fprintf (stderr, "Error, unknown itype %d (line %d, function '%s' in '%s')\n",
	     itype, __LINE__, __func__, __FILE__);
    return -999.0;
  }

  return ynew;
}

/*============================= end of c_inter() ========================*/

/*============================= c_gaussian_quadrature_test() ============*/

int c_gaussian_quadrature_test(int nstr, float *sza, double umu0)
{

  /* Test if the solar zenith angle coincides with one of
     the computational angles */

  int nn=0, iq=0, result=0;
  double umu0s=0.0, *cmu=NULL, *cwt=NULL;

  cmu = c_dbl_vector(0,nstr,"cmu");
  if (cmu==NULL) {
    fprintf(stderr,"Error allocating cmu!\n");
    return -1;
  }

  cwt = c_dbl_vector(0,nstr,"cwt");
  if (cwt==NULL) {
    fprintf(stderr,"Error allocating cwt!\n");
    return -1;
  }
    
  nn = nstr / 2.0;

  c_gaussian_quadrature ( nn, cmu, cwt );

  for (iq=1; iq<=nn; iq++) {
    if( fabs( (umu0 - CMU (iq)) / umu0 ) < 1.0e-4 ) {
      umu0s = umu0;
      if ( umu0 < CMU (iq) )
	umu0  = CMU (iq) * (1. - 1.1e-4);
      else
	umu0  = CMU (iq) * (1. + 1.1e-4);

      *sza   = acos (umu0)/DEG;
      fprintf(stderr,"%s %s %s %f %s %f\n",
	      "******* WARNING >>>>>> \n",
	      "SETDIS--beam angle=computational angle;\n",
	      "******* changing cosine of solar zenith angle, umu0, from ",
	      umu0s, "to", umu0 );
      result=-1;
    }
  }

  free(cwt);
  free(cmu);      
  return result;
}

/*============================= end of c_gaussian_quadrature_test() =====*/

/* * * * * * * * * * * * * * end of cdisort.c * * * * * * * * * * * */

