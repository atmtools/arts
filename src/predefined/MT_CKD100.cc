#include <propagationmatrix.h>

#include <vector>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "matpack_concepts.h"

namespace Absorption::PredefinedModel::MT_CKD100 {
Numeric XINT_FUN(const Numeric V1A,
                 const Numeric /* V2A */,
                 const Numeric DVA,
                 const matpack::vector_like auto& A,
                 const Index nA,
                 const Numeric VI) {
  // ----------------------------------------------------------------------
  //     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED
  //     FROM V1A TO V2A IN INCREMENTS OF DVA INTO XINT
  // ----------------------------------------------------------------------

  const Numeric ONEPL = 1.001;  // original value given in F77 code
  // FIXME const Numeric ONEMI  = 0.999;     // original value given in F77 code

  //const Numeric ONEPL  = 0.001;  // modified value for C/C++ code

  Numeric RECDVA = 1.00e0 / DVA;

  int J = (int)((VI - V1A) * RECDVA + ONEPL);
  Numeric VJ = V1A + DVA * (Numeric)(J - 1);
  Numeric P = RECDVA * (VI - VJ);
  Numeric C = (3.00e0 - 2.00e0 * P) * P * P;
  Numeric B = 0.500e0 * P * (1.00e0 - P);
  Numeric B1 = B * (1.00e0 - P);
  Numeric B2 = B * P;

  Numeric xint = 0.;
  if (J - 1 > 0 && J + 2 < nA) {
    xint = -A[J - 1] * B1 + A[J] * (1.00e0 - C + B2) + A[J + 1] * (C + B1) -
           A[J + 2] * B2;
  }

  /*
  cout << (J-1) << " <-> " << (J+2)
       << ",  V=" << VI << ", VJ=" << VJ << "\n";
  cout << "xint=" << xint  << " " << A[J-1] << " " << A[J] << " " << A[J+1] << " " << A[J+2] << "\n";
  */

  return xint;
}

Numeric XINT_FUN(const Numeric V1A,
                 const Numeric /* V2A */,
                 const Numeric DVA,
                 const matpack::vector_like auto& A,
                 const Numeric VI) {
  // ----------------------------------------------------------------------
  //     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED
  //     FROM V1A TO V2A IN INCREMENTS OF DVA INTO XINT
  // ----------------------------------------------------------------------

  const Numeric ONEPL = 1.001;  // original value given in F77 code
  // FIXME const Numeric ONEMI  = 0.999;     // original value given in F77 code

  //const Numeric ONEPL  = 0.001;  // modified value for C/C++ code

  Numeric RECDVA = 1.00e0 / DVA;

  int J = (int)((VI - V1A) * RECDVA + ONEPL);
  Numeric VJ = V1A + DVA * (Numeric)(J - 1);
  Numeric P = RECDVA * (VI - VJ);
  Numeric C = (3.00e0 - 2.00e0 * P) * P * P;
  Numeric B = 0.500e0 * P * (1.00e0 - P);
  Numeric B1 = B * (1.00e0 - P);
  Numeric B2 = B * P;

  Numeric xint = 0.;
  if (J - 1 > 0 && J + 2 < A.nelem()) {
    xint = -A[J - 1] * B1 + A[J] * (1.00e0 - C + B2) + A[J + 1] * (C + B1) -
           A[J + 2] * B2;
  }

  /*
  cout << (J-1) << " <-> " << (J+2)
       << ",  V=" << VI << ", VJ=" << VJ << "\n";
  cout << "xint=" << xint  << " " << A[J-1] << " " << A[J] << " " << A[J+1] << " " << A[J+2] << "\n";
  */

  return xint;
}

Numeric RADFN_FUN(const Numeric VI, const Numeric XKT) {
  // ---------------------------------------------------------------------- B18060
  //              LAST MODIFICATION:    12 AUGUST 1991                      B17940
  //                                                                        B17950
  //                 IMPLEMENTATION:    R.D. WORSHAM                        B17960
  //                                                                        B17970
  //            ALGORITHM REVISIONS:    S.A. CLOUGH                         B17980
  //                                    R.D. WORSHAM                        B17990
  //                                    J.L. MONCET                         B18000
  //                                                                        B18010
  //                                                                        B18020
  //                    ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18030
  //                    840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18040
  //                                                                        B18050
  //                                                                        B18070
  //              WORK SUPPORTED BY:    THE ARM PROGRAM                     B18080
  //                                    OFFICE OF ENERGY RESEARCH           B18090
  //                                    DEPARTMENT OF ENERGY                B18100
  //                                                                        B18110
  //                                                                        B18120
  //     SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18130
  //                                                                        B18140
  //                                            FASCOD3                     B18150
  //                                                                        B18160
  // ---------------------------------------------------------------------- B18060
  //                                                                        B18170
  //  IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED

  Numeric XVI = VI;
  Numeric RADFN = 0.00e0;

  if (XKT > 0.0) {
    Numeric XVIOKT = XVI / XKT;

    if (XVIOKT <= 0.01e0) {
      RADFN = 0.500e0 * XVIOKT * XVI;
    } else if (XVIOKT <= 10.0e0) {
      Numeric EXPVKT = exp(-XVIOKT);
      RADFN = XVI * (1.00e0 - EXPVKT) / (1.00e0 + EXPVKT);
    } else {
      RADFN = XVI;
    }
  } else {
    RADFN = XVI;
  }

  return RADFN;
}

//! Ported from legacy continua.  Original documentation
//! CKD version MT 1.00 O2-O2 collision induced absorption (fundamental band)
/*!
  Model reference:
  F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman,
  J.-M. Hartmann, Ch. Boulet,
  "Infrared collision-induced absorption by O2 near 6.4 microns for
  atmospheric applications: measurements and emprirical modeling",
  Appl. Optics, 35, 5911-5917, (1996).

   \param[out] pxsec        cross section (absorption/volume mixing ratio) of
                            O2-O2 CIA fundamental band according to CKD_MT 1.00   [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set
                            (Cin)<br>
                            or choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid            [Hz]
   \param    abs_p          predefined pressure grid             [Pa]
   \param    abs_t          predefined temperature grid          [K]
   \param    vmr            O2 volume mixing ratio profile       [1]

   \remark   F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman,
             J.-M. Hartmann, Ch. Boulet,<br>
             Infrared collision-induced absorption by O2 near 6.4 microns for
             atmospheric applications: measurements and emprirical modeling,<br>
             Appl. Optics, 35, 5911-5917, (1996).

   \note     This absorption model is taken from the FORTRAN77 code of
             CKD_MT version 1.00 written by<br>
             Atmospheric and Environmental Research Inc. (AER),<br>
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br>
             http://www.rtweb.aer.com/continuum_frame.html

   \author Thomas Kuhn
   \date 2002-28-08
 */
//! New implementation
void oxygen_cia(PropagationMatrix& propmat_clearsky,
                const Vector& f_grid,
                const Numeric p_pa,
                const Numeric t,
                const Numeric vmr) {
  constexpr Numeric O2O2_O2F_ckd_mt_100_v1 = 1340.000;
  constexpr Numeric O2O2_O2F_ckd_mt_100_v2 = 1850.000;
  constexpr Numeric O2O2_O2F_ckd_mt_100_dv = 5.000;
  constexpr int O2O2_O2F_ckd_mt_100_npt = 103;
  static constexpr std::array O2O2_O2Fo_ckd_mt_100{
      0.000E+00, 0.000E+00, 9.744E-09, 2.256E-08, 3.538E-08, 4.820E-08,
      6.100E-08, 7.400E-08, 8.400E-08, 9.600E-08, 1.200E-07, 1.620E-07,
      2.080E-07, 2.460E-07, 2.850E-07, 3.140E-07, 3.800E-07, 4.440E-07,
      5.000E-07, 5.710E-07, 6.730E-07, 7.680E-07, 8.530E-07, 9.660E-07,
      1.100E-06, 1.210E-06, 1.330E-06, 1.470E-06, 1.590E-06, 1.690E-06,
      1.800E-06, 1.920E-06, 2.040E-06, 2.150E-06, 2.260E-06, 2.370E-06,
      2.510E-06, 2.670E-06, 2.850E-06, 3.070E-06, 3.420E-06, 3.830E-06,
      4.200E-06, 4.450E-06, 4.600E-06, 4.530E-06, 4.280E-06, 3.960E-06,
      3.680E-06, 3.480E-06, 3.350E-06, 3.290E-06, 3.250E-06, 3.230E-06,
      3.230E-06, 3.210E-06, 3.190E-06, 3.110E-06, 3.030E-06, 2.910E-06,
      2.800E-06, 2.650E-06, 2.510E-06, 2.320E-06, 2.130E-06, 1.930E-06,
      1.760E-06, 1.590E-06, 1.420E-06, 1.250E-06, 1.110E-06, 9.900E-07,
      8.880E-07, 7.910E-07, 6.780E-07, 5.870E-07, 5.240E-07, 4.640E-07,
      4.030E-07, 3.570E-07, 3.200E-07, 2.900E-07, 2.670E-07, 2.420E-07,
      2.150E-07, 1.820E-07, 1.600E-07, 1.460E-07, 1.280E-07, 1.030E-07,
      8.700E-08, 8.100E-08, 7.100E-08, 6.400E-08, 5.807E-08, 5.139E-08,
      4.496E-08, 3.854E-08, 3.212E-08, 2.569E-08, 1.927E-08, 1.285E-08,
      6.423E-09, 0.000E+00};
  static constexpr std::array O2O2_O2Ft_ckd_mt_100{
      0.000E+00,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,
      4.670E+02,  4.000E+02,  3.150E+02,  3.790E+02,  3.680E+02,  4.750E+02,
      5.210E+02,  5.310E+02,  5.120E+02,  4.420E+02,  4.440E+02,  4.300E+02,
      3.810E+02,  3.350E+02,  3.240E+02,  2.960E+02,  2.480E+02,  2.150E+02,
      1.930E+02,  1.580E+02,  1.270E+02,  1.010E+02,  7.100E+01,  3.100E+01,
      -6.000E+00, -2.600E+01, -4.700E+01, -6.300E+01, -7.900E+01, -8.800E+01,
      -8.800E+01, -8.700E+01, -9.000E+01, -9.800E+01, -9.900E+01, -1.090E+02,
      -1.340E+02, -1.600E+02, -1.670E+02, -1.640E+02, -1.580E+02, -1.530E+02,
      -1.510E+02, -1.560E+02, -1.660E+02, -1.680E+02, -1.730E+02, -1.700E+02,
      -1.610E+02, -1.450E+02, -1.260E+02, -1.080E+02, -8.400E+01, -5.900E+01,
      -2.900E+01, 4.000E+00,  4.100E+01,  7.300E+01,  9.700E+01,  1.230E+02,
      1.590E+02,  1.980E+02,  2.200E+02,  2.420E+02,  2.560E+02,  2.810E+02,
      3.110E+02,  3.340E+02,  3.190E+02,  3.130E+02,  3.210E+02,  3.230E+02,
      3.100E+02,  3.150E+02,  3.200E+02,  3.350E+02,  3.610E+02,  3.780E+02,
      3.730E+02,  3.380E+02,  3.190E+02,  3.460E+02,  3.220E+02,  2.910E+02,
      2.900E+02,  3.500E+02,  3.710E+02,  5.040E+02,  4.000E+02,  4.000E+02,
      4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,
      4.000E+02,  4.000E+02};

  // ************************** CKD stuff ************************************

  constexpr Numeric xLosmt = 2.686763e19;  // Loschmidt Number [molecules/cm^3]
  constexpr Numeric T1 = 273.0e0;
  constexpr Numeric TO = 296.0e0;
  constexpr Numeric PO = 1013.0e0;

  // It is assumed here that f_grid is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_grid. n_f_new < n_f
  Numeric V1ABS = f_grid[0] / (Constant::c * 1.00e2);  // [cm^-1]
  Numeric V2ABS =
      f_grid[f_grid.nelem() - 1] / (Constant::c * 1.00e2);  // [cm^-1]

  // ------------------- subroutine O2_VER_1 ----------------------------

  // retrieve the appropriate array sequence of the CKD model array.
  constexpr Numeric DVC = O2O2_O2F_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;

  int I1 = (int)((V1C - O2O2_O2F_ckd_mt_100_v1) / O2O2_O2F_ckd_mt_100_dv);
  if (V1C < O2O2_O2F_ckd_mt_100_v1) I1 = -1;
  V1C = O2O2_O2F_ckd_mt_100_v1 + (O2O2_O2F_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int)((V2C - O2O2_O2F_ckd_mt_100_v1) / O2O2_O2F_ckd_mt_100_dv);

  int NPTC = I2 - I1 + 3;
  if (NPTC > O2O2_O2F_ckd_mt_100_npt) NPTC = O2O2_O2F_ckd_mt_100_npt + 1;

  V2C = V1C + O2O2_O2F_ckd_mt_100_dv * (Numeric)(NPTC - 1);

  if (NPTC < 1) {
    return;
  }

  Vector xo2(NPTC + 1, 0.);
  Vector xo2t(NPTC + 1, 0.);

  for (Index J = 1; J <= NPTC; ++J) {
    Index I = I1 + J;
    if ((I > 0) && (I <= O2O2_O2F_ckd_mt_100_npt)) {
      xo2[J] = O2O2_O2Fo_ckd_mt_100[I];
      xo2t[J] = O2O2_O2Ft_ckd_mt_100[I];
    }
  }

  // ------------------- subroutine O2_VER_1 ----------------------------

  // Loop pressure/temperature:
  Numeric Tave = t;                  // [K]
  Numeric Pave = (p_pa * 1.000e-2);  // [hPa]
  // FIXME Numeric vmro2  = vmr[i];                                 // [1]
  Numeric WTOT = xLosmt * (Pave / PO) * (T1 / Tave);  // [molecules/cm^2]
  Numeric tau_fac = WTOT * (Pave / PO) * (T1 / Tave);

  Numeric XKT = Tave / 1.4387752;  // = (T*k_B) / (h*c)

  Numeric xktfac = (1.000e0 / TO) - (1.000e0 / Tave);  // [1/K]
  Numeric factor = (1.000e0 / xLosmt);

  // Molecular cross section calculated by CKD.
  // The cross sectionis calculated on the predefined
  // CKD wavenumber grid.
  Vector k(NPTC + 1 + 1, 0.);  // [1/cm]
  for (Index J = 1; J <= NPTC; ++J) {
    Numeric VJ = V1C + (DVC * (Numeric)(J - 1));
    Numeric SO2 = 0.0e0;
    if (xo2[J] > 0.0e0) {
      Numeric C0 = factor * xo2[J] * exp(xo2t[J] * xktfac) / VJ;
      SO2 = tau_fac * C0;
    }

    // CKD cross section without radiative field
    k[J] = SO2 * RADFN_FUN(VJ, XKT);  // [1]
  }

  // Loop input frequency array. The previously calculated cross section
  // has therefore to be interpolated on the input frequencies.
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // calculate the associated wave number (= 1/wavelength)
    Numeric V = f_grid[s] / (Constant::c * 1.00e2);  // [cm^-1]
    if ((V > O2O2_O2F_ckd_mt_100_v1) && (V < O2O2_O2F_ckd_mt_100_v2)) {
      // arts cross section [1/m]
      // interpolate the k vector on the f_grid grid
      propmat_clearsky.Kjj()[s] +=
          vmr * 1.000e2 * XINT_FUN(V1C, V2C, DVC, k, V);
    }
  }
}

//! Ported from legacy continua.  Original documentation
//! CKD version MT 1.00 O2 v0<-v0 band absorption
/*!
  Model reference:
  CKD_MT 1.00 implementation of oxygen collision induced fundamental model of
  O2 continuum formulated by
  Mate et al. over the spectral region 7550-8486 cm-1:
  B. Mate, C. Lugez, G.T. Fraser, W.J. Lafferty,
  "Absolute Intensities for the O2 1.27 micron
  continuum absorption",
  J. Geophys. Res., 104, 30,585-30,590, 1999.

  The units of these continua coefficients are  1 / (amagat_O2*amagat_air)

  Also, refer to the paper "Observed  Atmospheric
  Collision Induced Absorption in Near Infrared Oxygen Bands",
  Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
  Journal of Geophysical Research (1997).

   \param[out] pxsec        cross section (absorption/volume mixing ratio) of
                            O2 v0<-v0 band according to CKD_MT 1.00  [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set
                            (Cin)<br>
                            or choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid            [Hz]
   \param    abs_p          predefined pressure grid             [Pa]
   \param    abs_t          predefined temperature grid          [K]
   \param    vmr            O2 volume mixing ratio profile       [1]
   \param    abs_n2         N2 volume mixing ratio profile       [1]

   \remark   B. Mate, C. Lugez, G.T. Fraser, W.J. Lafferty,<br>
             Absolute Intensities for the O2 1.27 micron continuum absorption,<br>
             J. Geophys. Res., 104, 30,585-30,590, 1999.

   \note     This absorption model is taken from the FORTRAN77 code of
             CKD_MT version 1.00 written by<br>
             Atmospheric and Environmental Research Inc. (AER),<br>
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br>
             http://www.rtweb.aer.com/continuum_frame.html<br>
             <br>
       Oxygen band absorption model for the \f$a^1\Delta_g\f$
             \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly
             \f$X^3\Sigma^-_g\f$ band system considering the
             \f$\nu=0\f$
             \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly
             \f$\nu=0\f$
             transitions.

   \author   Thomas Kuhn
   \date     2002-28-08
 */
//! New implementation
void oxygen_v0v0(PropagationMatrix& propmat_clearsky,
                 const Vector& f_grid,
                 const Numeric p_pa,
                 const Numeric t,
                 const Numeric vmr,
                 const Numeric n2) {
  constexpr Numeric O2_00_ckd_mt_100_v1 = 7536.000e0;
  constexpr Numeric O2_00_ckd_mt_100_v2 = 8500.000e0;
  constexpr Numeric O2_00_ckd_mt_100_dv = 2.000e0;
  constexpr int O2_00_ckd_mt_100_npt = 483;
  static constexpr std::array O2_00_ckd_mt_100{
      0.000E+00, 0.000E+00, 4.355E-11, 8.709E-11, 1.742E-10, 3.484E-10,
      6.968E-10, 1.394E-09, 2.787E-09, 3.561E-09, 3.314E-09, 3.368E-09,
      3.435E-09, 2.855E-09, 3.244E-09, 3.447E-09, 3.891E-09, 4.355E-09,
      3.709E-09, 4.265E-09, 4.772E-09, 4.541E-09, 4.557E-09, 4.915E-09,
      4.688E-09, 5.282E-09, 5.755E-09, 5.096E-09, 5.027E-09, 4.860E-09,
      4.724E-09, 5.048E-09, 5.248E-09, 5.473E-09, 4.852E-09, 5.362E-09,
      6.157E-09, 6.150E-09, 6.347E-09, 6.388E-09, 6.213E-09, 6.521E-09,
      8.470E-09, 8.236E-09, 8.269E-09, 8.776E-09, 9.122E-09, 9.189E-09,
      9.778E-09, 8.433E-09, 9.964E-09, 9.827E-09, 1.064E-08, 1.063E-08,
      1.031E-08, 1.098E-08, 1.156E-08, 1.295E-08, 1.326E-08, 1.467E-08,
      1.427E-08, 1.452E-08, 1.456E-08, 1.554E-08, 1.605E-08, 1.659E-08,
      1.754E-08, 1.757E-08, 1.876E-08, 1.903E-08, 1.876E-08, 1.869E-08,
      2.036E-08, 2.203E-08, 2.221E-08, 2.284E-08, 2.288E-08, 2.394E-08,
      2.509E-08, 2.663E-08, 2.720E-08, 2.839E-08, 2.923E-08, 2.893E-08,
      2.949E-08, 2.962E-08, 3.057E-08, 3.056E-08, 3.364E-08, 3.563E-08,
      3.743E-08, 3.813E-08, 3.946E-08, 4.082E-08, 4.201E-08, 4.297E-08,
      4.528E-08, 4.587E-08, 4.704E-08, 4.962E-08, 5.115E-08, 5.341E-08,
      5.365E-08, 5.557E-08, 5.891E-08, 6.084E-08, 6.270E-08, 6.448E-08,
      6.622E-08, 6.939E-08, 7.233E-08, 7.498E-08, 7.749E-08, 8.027E-08,
      8.387E-08, 8.605E-08, 8.888E-08, 9.277E-08, 9.523E-08, 9.880E-08,
      1.037E-07, 1.076E-07, 1.114E-07, 1.151E-07, 1.203E-07, 1.246E-07,
      1.285E-07, 1.345E-07, 1.408E-07, 1.465E-07, 1.519E-07, 1.578E-07,
      1.628E-07, 1.685E-07, 1.760E-07, 1.847E-07, 1.929E-07, 2.002E-07,
      2.070E-07, 2.177E-07, 2.262E-07, 2.365E-07, 2.482E-07, 2.587E-07,
      2.655E-07, 2.789E-07, 2.925E-07, 3.023E-07, 3.153E-07, 3.296E-07,
      3.409E-07, 3.532E-07, 3.680E-07, 3.859E-07, 3.951E-07, 4.074E-07,
      4.210E-07, 4.381E-07, 4.588E-07, 4.792E-07, 4.958E-07, 5.104E-07,
      5.271E-07, 5.501E-07, 5.674E-07, 5.913E-07, 6.243E-07, 6.471E-07,
      6.622E-07, 6.831E-07, 6.987E-07, 7.159E-07, 7.412E-07, 7.698E-07,
      7.599E-07, 7.600E-07, 7.918E-07, 8.026E-07, 8.051E-07, 8.049E-07,
      7.914E-07, 7.968E-07, 7.945E-07, 7.861E-07, 7.864E-07, 7.741E-07,
      7.675E-07, 7.592E-07, 7.400E-07, 7.362E-07, 7.285E-07, 7.173E-07,
      6.966E-07, 6.744E-07, 6.597E-07, 6.413E-07, 6.265E-07, 6.110E-07,
      5.929E-07, 5.717E-07, 5.592E-07, 5.411E-07, 5.235E-07, 5.061E-07,
      4.845E-07, 4.732E-07, 4.593E-07, 4.467E-07, 4.328E-07, 4.161E-07,
      4.035E-07, 3.922E-07, 3.820E-07, 3.707E-07, 3.585E-07, 3.475E-07,
      3.407E-07, 3.317E-07, 3.226E-07, 3.134E-07, 3.016E-07, 2.969E-07,
      2.894E-07, 2.814E-07, 2.749E-07, 2.657E-07, 2.610E-07, 2.536E-07,
      2.467E-07, 2.394E-07, 2.337E-07, 2.302E-07, 2.241E-07, 2.191E-07,
      2.140E-07, 2.093E-07, 2.052E-07, 1.998E-07, 1.963E-07, 1.920E-07,
      1.862E-07, 1.834E-07, 1.795E-07, 1.745E-07, 1.723E-07, 1.686E-07,
      1.658E-07, 1.629E-07, 1.595E-07, 1.558E-07, 1.523E-07, 1.498E-07,
      1.466E-07, 1.452E-07, 1.431E-07, 1.408E-07, 1.381E-07, 1.362E-07,
      1.320E-07, 1.298E-07, 1.262E-07, 1.247E-07, 1.234E-07, 1.221E-07,
      1.197E-07, 1.176E-07, 1.142E-07, 1.121E-07, 1.099E-07, 1.081E-07,
      1.073E-07, 1.061E-07, 1.041E-07, 1.019E-07, 9.969E-08, 9.727E-08,
      9.642E-08, 9.487E-08, 9.318E-08, 9.116E-08, 9.046E-08, 8.827E-08,
      8.689E-08, 8.433E-08, 8.324E-08, 8.204E-08, 8.036E-08, 7.951E-08,
      7.804E-08, 7.524E-08, 7.392E-08, 7.227E-08, 7.176E-08, 6.975E-08,
      6.914E-08, 6.859E-08, 6.664E-08, 6.506E-08, 6.368E-08, 6.262E-08,
      6.026E-08, 6.002E-08, 5.866E-08, 5.867E-08, 5.641E-08, 5.589E-08,
      5.499E-08, 5.309E-08, 5.188E-08, 5.139E-08, 4.991E-08, 4.951E-08,
      4.833E-08, 4.640E-08, 4.524E-08, 4.479E-08, 4.304E-08, 4.228E-08,
      4.251E-08, 4.130E-08, 3.984E-08, 3.894E-08, 3.815E-08, 3.732E-08,
      3.664E-08, 3.512E-08, 3.463E-08, 3.503E-08, 3.218E-08, 3.253E-08,
      3.107E-08, 2.964E-08, 2.920E-08, 2.888E-08, 2.981E-08, 2.830E-08,
      2.750E-08, 2.580E-08, 2.528E-08, 2.444E-08, 2.378E-08, 2.413E-08,
      2.234E-08, 2.316E-08, 2.199E-08, 2.088E-08, 1.998E-08, 1.920E-08,
      1.942E-08, 1.859E-08, 1.954E-08, 1.955E-08, 1.749E-08, 1.720E-08,
      1.702E-08, 1.521E-08, 1.589E-08, 1.469E-08, 1.471E-08, 1.543E-08,
      1.433E-08, 1.298E-08, 1.274E-08, 1.226E-08, 1.204E-08, 1.201E-08,
      1.298E-08, 1.220E-08, 1.220E-08, 1.096E-08, 1.080E-08, 9.868E-09,
      9.701E-09, 1.130E-08, 9.874E-09, 9.754E-09, 9.651E-09, 9.725E-09,
      8.413E-09, 7.705E-09, 7.846E-09, 8.037E-09, 9.163E-09, 8.098E-09,
      8.160E-09, 7.511E-09, 7.011E-09, 6.281E-09, 6.502E-09, 7.323E-09,
      7.569E-09, 5.941E-09, 5.867E-09, 5.676E-09, 4.840E-09, 5.063E-09,
      5.207E-09, 4.917E-09, 5.033E-09, 5.356E-09, 3.795E-09, 4.983E-09,
      4.600E-09, 3.635E-09, 3.099E-09, 2.502E-09, 3.823E-09, 3.464E-09,
      4.332E-09, 3.612E-09, 3.682E-09, 3.709E-09, 3.043E-09, 3.593E-09,
      3.995E-09, 4.460E-09, 3.583E-09, 3.290E-09, 3.132E-09, 2.812E-09,
      3.109E-09, 3.874E-09, 3.802E-09, 4.024E-09, 3.901E-09, 2.370E-09,
      1.821E-09, 2.519E-09, 4.701E-09, 3.855E-09, 4.685E-09, 5.170E-09,
      4.387E-09, 4.148E-09, 4.043E-09, 3.545E-09, 3.392E-09, 3.609E-09,
      4.635E-09, 3.467E-09, 2.558E-09, 3.389E-09, 2.672E-09, 2.468E-09,
      1.989E-09, 2.816E-09, 4.023E-09, 2.664E-09, 2.219E-09, 3.169E-09,
      1.654E-09, 3.189E-09, 2.535E-09, 2.618E-09, 3.265E-09, 2.138E-09,
      1.822E-09, 2.920E-09, 2.002E-09, 1.300E-09, 3.764E-09, 3.212E-09,
      3.222E-09, 2.961E-09, 2.108E-09, 1.708E-09, 2.636E-09, 2.937E-09,
      2.939E-09, 2.732E-09, 2.218E-09, 1.046E-09, 6.419E-10, 1.842E-09,
      1.112E-09, 1.265E-09, 4.087E-09, 2.044E-09, 1.022E-09, 5.109E-10,
      2.554E-10, 1.277E-10, 6.386E-11, 0.000E+00};

  // ************************** CKD stuff ************************************

  // FIXME const Numeric xLosmt    = 2.686763e19; // Loschmidt Number [molecules/cm^3]
  constexpr Numeric T1 = 273.0e0;
  // FIXME const Numeric TO        =  296.0e0;
  constexpr Numeric PO = 1013.0e0;

  // It is assumed here that f_grid is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_grid. n_f_new < n_f
  Numeric V1ABS = f_grid[0] / (Constant::c * 1.00e2);  // [cm^-1]
  Numeric V2ABS =
      f_grid[f_grid.nelem() - 1] / (Constant::c * 1.00e2);  // [cm^-1]

  // ------------------- subroutine O2INF1 ----------------------------

  // retrieve the appropriate array sequence of the CKD model array.
  Numeric DVC = O2_00_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;

  int I1 = (int)((V1C - O2_00_ckd_mt_100_v1) / O2_00_ckd_mt_100_dv);
  if (V1C < O2_00_ckd_mt_100_v1) I1 = I1 - 1;
  V1C = O2_00_ckd_mt_100_v1 + (O2_00_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int)((V2C - O2_00_ckd_mt_100_v1) / O2_00_ckd_mt_100_dv);

  int NPTC = I2 - I1 + 3;

  V2C = V1C + O2_00_ckd_mt_100_dv * (Numeric)(NPTC - 1);

  if (NPTC < 1) {
    return;
  }

  std::vector<Numeric> CO(NPTC + 1);

  for (Index J = 1; J <= NPTC; ++J) {
    CO[J] = 0.000e0;
    Index I = I1 + J;
    if ((I > 0) && (I <= O2_00_ckd_mt_100_npt)) {
      Numeric VJ = V1C + (DVC * (Numeric)(J - 1));
      CO[J] = O2_00_ckd_mt_100[I] / VJ;
    }
  }

  // ------------------- subroutine O2INF1 ----------------------------

  // Loop pressure/temperature:
  Numeric Tave = t;                  // [K]
  Numeric Pave = (p_pa * 1.000e-2);  // [hPa]
  Numeric vmro2 = vmr;               // [1]
  Numeric vmrn2 = n2;                // [1]
  Numeric ADJWO2 = (vmro2 + 0.300e0 * vmrn2) / 0.446e0 * (Pave / PO) *
                   (Pave / PO) * (T1 / Tave) * (T1 / Tave);
  Numeric XKT = Tave / 1.4387752e0;  // = (T*k_B) / (h*c)

  // Molecular cross section calculated by CKD.
  // The cross sectionis calculated on the predefined
  // CKD wavenumber grid. The abs. coeff. is then the
  // cross section times the number density.
  std::vector<Numeric> k(NPTC + 1);  // [1/cm]
  k[0] = 0.00e0;                     // not used array field
  for (Index J = 1; J <= NPTC; ++J) {
    Numeric VJ = V1C + (DVC * (Numeric)(J - 1));
    Numeric SO2 = 0.0e0;
    if (CO[J] > 0.0e0) {
      SO2 = ADJWO2 * CO[J];
    }

    // CKD (cross section * number density) with radiative field
    k[J] = SO2 * RADFN_FUN(VJ, XKT);  // [1/cm]
  }

  // Loop input frequency array. The previously calculated cross section
  // has therefore to be interpolated on the input frequencies.
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // calculate the associated wave number (= 1/wavelength)
    Numeric V = f_grid[s] / (Constant::c * 1.00e2);  // [cm^-1]
    if ((V > O2_00_ckd_mt_100_v1) && (V < O2_00_ckd_mt_100_v2)) {
      // arts cross section [1/m]
      // interpolate the k vector on the f_grid grid
      propmat_clearsky.Kjj()[s] +=
          vmr * 1.000e2 * XINT_FUN(V1C, V2C, DVC, k, NPTC + 1, V);
    }
  }
}

//! Ported from legacy continua.  Original documentation
//! CKD version MT 1.00 O2 v1<-v0 band absorption
/*!
  Model reference:
  CKD_MT 1.00 implementation of oxygen v1<-v0 band model of
  Mlawer, Clough, Brown, Stephen, Landry, Goldman, Murcray,
  "Observed  Atmospheric Collision Induced Absorption in Near Infrared Oxygen Bands",
  Journal of Geophysical Research, vol 103, no. D4, pp. 3859-3863, 1998.

   \param[out] pxsec        cross section (absorption/volume mixing ratio) of
                            O2 v1<-v0 band according to CKD_MT 1.00  [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set
                            (Cin)<br>
                            or choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid            [Hz]
   \param    abs_p          predefined pressure grid             [Pa]
   \param    abs_t          predefined temperature grid          [K]
   \param    vmr            O2 volume mixing ratio profile       [1]

   \remark   Mlawer, Clough, Brown, Stephen, Landry, Goldman, Murcray,<br>
             Observed  Atmospheric Collision Induced Absorption in Near Infrared Oxygen Bands,<br>
             J. Geophys. Res., 103, D4, 3859-3863, 1998.

   \note     This absorption model is taken from the FORTRAN77 code of
             CKD_MT version 1.00 written by<br>
             Atmospheric and Environmental Research Inc. (AER),<br>
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br>
       http://www.rtweb.aer.com/continuum_frame.html<br>
             <br>
       Oxygen band absorption model for the \f$a^1\Delta_g\f$
             \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly
             \f$X^3\Sigma^-_g\f$ band system considering the
             \f$\nu=0\f$
             \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly
             \f$\nu=1\f$
             transitions.

   \author Thomas Kuhn
   \date 2002-28-08
 */
//! New implementation
void oxygen_v0v1(PropagationMatrix& propmat_clearsky,
                 const Vector& f_grid,
                 const Numeric p_pa,
                 const Numeric t,
                 const Numeric vmr) {
  using Math::pow2;

  constexpr Numeric O2_10_ckd_mt_100_v1 = 9100.000e0;
  constexpr Numeric O2_10_ckd_mt_100_v2 = 11000.000e0;
  constexpr Numeric O2_10_ckd_mt_100_dv = 2.000e0;

  // ************************** CKD stuff ************************************

  constexpr Numeric xLosmt = 2.686763e19;  // Loschmidt Number [molecules/cm^3]
  constexpr Numeric T1 = 273.0e0;
  constexpr Numeric TO = 296.0e0;
  constexpr Numeric PO = 1013.0e0;
  // FIXME const Numeric vmr_argon = 9.000e-3;    // VMR of argon is assumed to be const.

  // CKD_MT 1.00 implementation of oxygen v1<-v0 band model of
  // Mlawer, Clough, Brown, Stephen, Landry, Goldman, Murcray,
  // "Observed  Atmospheric Collision Induced Absorption in Near Infrared Oxygen Bands",
  // Journal of Geophysical Research, vol 103, no. D4, pp. 3859-3863, 1998.
  const Numeric V1S = O2_10_ckd_mt_100_v1;
  const Numeric V2S = O2_10_ckd_mt_100_v2;
  const Numeric DVS = O2_10_ckd_mt_100_dv;
  const Numeric V1_osc = 9375.000e0;
  const Numeric HW1 = 58.960e0;
  const Numeric V2_osc = 9439.000e0;
  const Numeric HW2 = 45.040e0;
  const Numeric S1 = 1.166e-4;
  const Numeric S2 = 3.086e-5;

  // It is assumed here that f_grid is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_grid. n_f_new < n_f
  Numeric V1ABS = f_grid[0] / (Constant::c * 1.00e2);  // [cm^-1]
  Numeric V2ABS =
      f_grid[f_grid.nelem() - 1] / (Constant::c * 1.00e2);  // [cm^-1]

  // ------------------- subroutine O2INF2 ----------------------------

  // retrieve the appropriate array sequence of the CKD model array.
  Numeric DVC = DVS;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;

  int NPTC = (int)(((V2C - V1C) / DVC) + 3);

  V2C = V1C + (DVC * (Numeric)(NPTC - 1));

  if (NPTC < 1) {
    return;
  }

  std::vector<Numeric> C(NPTC + 1);
  C[0] = 0.000e0;  // not used field of array

  for (Index J = 1; J <= NPTC; ++J) {
    C[J] = 0.000e0;
    Numeric VJ = V1C + (DVC * (Numeric)(J - 1));

    if ((VJ > V1S) && (VJ < V2S)) {
      Numeric DV1 = VJ - V1_osc;
      Numeric DV2 = VJ - V2_osc;

      Numeric DAMP1 = 1.00e0;
      Numeric DAMP2 = 1.00e0;

      if (DV1 < 0.00e0) {
        DAMP1 = exp(DV1 / 176.100e0);
      }

      if (DV2 < 0.00e0) {
        DAMP2 = exp(DV2 / 176.100e0);
      }

      Numeric O2INF = 0.31831e0 *
                      (((S1 * DAMP1 / HW1) / (1.000e0 + pow2(DV1 / HW1))) +
                       ((S2 * DAMP2 / HW2) / (1.000e0 + pow2(DV2 / HW2)))) *
                      1.054e0;
      C[J] = O2INF / VJ;
    }
  }

  // ------------------- subroutine O2INF2 ----------------------------

  // Loop pressure/temperature:
  Numeric Tave = t;                  // [K]
  Numeric Pave = (p_pa * 1.000e-2);  // [hPa]
  Numeric vmro2 = vmr;               // [1]
  Numeric WTOT =
      1.000e-20 * xLosmt * (Pave / PO) * (T1 / Tave);  // [molecules/cm^2]
  Numeric ADJWO2 = (vmro2 / 0.209e0) * WTOT * (Pave / PO) * (TO / Tave);
  Numeric XKT = Tave / 1.4387752;  // = (T*k_B) / (h*c)

  // Molecular cross section calculated by CKD.
  // The cross sectionis calculated on the predefined
  // CKD wavenumber grid.
  std::vector<Numeric> k(NPTC + 1);  // [1/cm]
  k[0] = 0.00e0;                     // not used array field
  for (Index J = 1; J <= NPTC; ++J) {
    Numeric VJ = V1C + (DVC * (Numeric)(J - 1));
    Numeric SO2 = 0.0e0;
    if (C[J] > 0.0e0) {
      SO2 = ADJWO2 * C[J];
    }

    // CKD cross section without radiative field
    k[J] = SO2 * RADFN_FUN(VJ, XKT);  // [1]
  }

  // Loop input frequency array. The previously calculated cross section
  // has therefore to be interpolated on the input frequencies.
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // calculate the associated wave number (= 1/wavelength)
    Numeric V = f_grid[s] / (Constant::c * 1.00e2);  // [cm^-1]
    if ((V > V1S) && (V < V2S)) {
      // arts cross section [1/m]
      // interpolate the k vector on the f_grid grid
      propmat_clearsky.Kjj()[s] +=
          vmr * 1.000e2 * XINT_FUN(V1C, V2C, DVC, k, NPTC + 1, V);
    }
  }
}
}  // namespace Absorption::PredefinedModel::MT_CKD100
