#include <arts_constexpr_math.h>
#include <propagationmatrix.h>

#include <array>
#include <numeric>

#include "arts_constants.h"

namespace Absorption::PredefinedModel::PWR98 {
//! Ported from legacy continua.  Original documentation
//! PWR98H2OAbsModel
/*!
   \param[out] pxsec        cross section (absorption/volume mixing ratio) of
                            H2O (lines+continuum) according to P. W. Rosenkranz, 1998 [1/m]
   \param    CCin           scaling factor for the H2O-continuum  [1]
   \param    CLin           scaling factor for the line strengths [1]
   \param    CWin           scaling factor for the line widths    [1]
   \param    model          allows user defined input parameter set
                            (CCin, CLin, and CWin)<br> or choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid       [Hz]
   \param    abs_p          predefined pressure grid       [Pa]
   \param    abs_t          predefined temperature grid     [K]
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, and CWin
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', 'RosenkranzLines', 'RosenkranzContinuum',
             and 'user'. See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and
             Radio Science, Vol. 34(4), 1025, 1999.

   \author Thomas Kuhn
   \date 2001-11-05
 */
//! New implementation
void water(PropagationMatrix& propmat_clearsky,
           const Vector& f_grid,
           const Numeric p_pa,
           const Numeric t,
           const Numeric vmr) noexcept {
  using Math::pow2;
  using Math::pow3;

  //   REFERENCES:
  //   LINE INTENSITIES FROM HITRAN92 (SELECTION THRESHOLD=
  //     HALF OF CONTINUUM ABSORPTION AT 1000 MB).
  //   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED:
  //     H.J.LIEBE AND T.A.DILLON, J.CHEM.PHYS. V.50, PP.727-732 (1969) &
  //     H.J.LIEBE ET AL., JQSRT V.9, PP. 31-47 (1969)  (22GHz);
  //     A.BAUER ET AL., JQSRT V.37, PP.531-539 (1987) &
  //     ASA WORKSHOP (SEPT. 1989) (380GHz);
  //     AND A.BAUER ET AL., JQSRT V.41, PP.49-54 (1989) (OTHER LINES).
  //   AIR-BROADENED CONTINUUM BASED ON LIEBE & LAYTON, NTIA
  //     REPORT 87-224 (1987); SELF-BROADENED CONTINUUM BASED ON
  //     LIEBE ET AL, AGARD CONF. PROC. 542 (MAY 1993),
  //     BUT READJUSTED FOR LINE SHAPE OF
  //     CLOUGH et al, ATMOS. RESEARCH V.23, PP.229-241 (1989).
  //
  // Coefficients are from P. W. Rosenkranz., Radio Science, 33(4), 919, 1998
  // line frequencies [GHz]
  static constexpr std::array PWRfl{22.2350800,
                                    183.3101170,
                                    321.2256400,
                                    325.1529190,
                                    380.1973720,
                                    439.1508120,
                                    443.0182950,
                                    448.0010750,
                                    470.8889470,
                                    474.6891270,
                                    488.4911330,
                                    556.9360020,
                                    620.7008070,
                                    752.0332270,
                                    916.1715820};

  // line intensities at 300K [Hz * cm2] (see Janssen Appendix to Chap.2 for this)
  static constexpr std::array PWRs1{1.31e-14,
                                    2.273e-12,
                                    8.036e-14,
                                    2.694e-12,
                                    2.438e-11,
                                    2.179e-12,
                                    4.624e-13,
                                    2.562e-11,
                                    8.369e-13,
                                    3.263e-12,
                                    6.659e-13,
                                    1.531e-9,
                                    1.707e-11,
                                    1.011e-9,
                                    4.227e-11};

  // T coeff. of intensities [1]
  static constexpr std::array PWRb2{2.144,
                                    0.668,
                                    6.179,
                                    1.541,
                                    1.048,
                                    3.595,
                                    5.048,
                                    1.405,
                                    3.597,
                                    2.379,
                                    2.852,
                                    0.159,
                                    2.391,
                                    0.396,
                                    1.441};

  // air-broadened width parameters at 300K [GHz/hPa]
  static constexpr std::array PWRw3{0.00281,
                                    0.00281,
                                    0.00230,
                                    0.00278,
                                    0.00287,
                                    0.00210,
                                    0.00186,
                                    0.00263,
                                    0.00215,
                                    0.00236,
                                    0.00260,
                                    0.00321,
                                    0.00244,
                                    0.00306,
                                    0.00267};

  // T-exponent of air-broadening [1]
  static constexpr std::array PWRx{0.69,
                                   0.64,
                                   0.67,
                                   0.68,
                                   0.54,
                                   0.63,
                                   0.60,
                                   0.66,
                                   0.66,
                                   0.65,
                                   0.69,
                                   0.69,
                                   0.71,
                                   0.68,
                                   0.70};

  // self-broadened width parameters at 300K [GHz/hPa]
  static constexpr std::array PWRws{0.01349,
                                    0.01491,
                                    0.01080,
                                    0.01350,
                                    0.01541,
                                    0.00900,
                                    0.00788,
                                    0.01275,
                                    0.00983,
                                    0.01095,
                                    0.01313,
                                    0.01320,
                                    0.01140,
                                    0.01253,
                                    0.01275};

  // T-exponent of self-broadening [1]
  static constexpr std::array PWRxs{0.61,
                                    0.85,
                                    0.54,
                                    0.74,
                                    0.89,
                                    0.52,
                                    0.50,
                                    0.67,
                                    0.65,
                                    0.64,
                                    0.72,
                                    1.00,
                                    0.68,
                                    0.84,
                                    0.78};

  // Loop pressure/temperature:
  // here the total pressure is not multiplied by the H2O vmr for the
  // P_H2O calculation because we calculate pxsec and not abs: abs = vmr * pxsec
  const Numeric pvap_dummy = 1e-2 * p_pa;
  // water vapor partial pressure [hPa]
  const Numeric pvap = 1e-2 * p_pa * vmr;
  // dry air partial pressure [hPa]
  const Numeric pda = (1e-2 * p_pa) - pvap;
  // Rosenkranz number density  (Rosenkranz H2O mass density in [g/m³])
  // [g/m³]    =  [g*K / Pa*m³]  *  [Pa/K]
  // rho       =   (M_H2O / R)   *  (P_H2O / T)
  // rho       =      2.1667     *  abs_p * vmr / abs_t
  // den       = 3.335e16 * rho
  // FIXME Numeric den        = 3.335e16 * (2.1667 * abs_p[i] * vmr[i] / abs_t[i]);
  const Numeric den_dummy = 3.335e16 * (2.1667 * p_pa / t);
  // inverse relative temperature [1]
  const Numeric ti = (300.0 / t);
  const Numeric ti2 = pow(ti, (Numeric)2.5);

  // continuum term [Np/km/GHz2]
  const Numeric con = pvap_dummy * pow3(ti) * 1.000e-9 *
                      ((0.543 * pda) + (17.96 * pvap * pow(ti, (Numeric)4.5)));

  // Loop over input frequency
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // input frequency in [GHz]
    const Numeric ff = f_grid[s] * 1e-9;
    // line contribution at position f
    Numeric sum = 0.000;

    // Loop over spectral lines

    for (Index l = 0; l < 15; l++) {
      const Numeric width = (PWRw3[l] * pda * pow(ti, PWRx[l])) +
                            (PWRws[l] * pvap * pow(ti, PWRxs[l]));
      //        Numeric width    = CW * ( PWRw3[l] * pda  * pow(ti, PWRx[l]) +
      //          PWRws[l] * pvap * pow(ti, PWRxs[l]) );
      const Numeric wsq = width * width;
      const Numeric strength = PWRs1[l] * ti2 * exp(PWRb2[l] * (1.0 - ti));
      // frequency differences
      const Numeric df0 = ff - PWRfl[l];
      const Numeric df1 = ff + PWRfl[l];
      // use Clough's definition of local line contribution
      const Numeric base = width / (wsq + 562500.000);
      // positive and negative resonances
      Numeric res = 0.000;
      if (fabs(df0) < 750.0) res += width / (df0 * df0 + wsq) - base;
      if (fabs(df1) < 750.0) res += width / (df1 * df1 + wsq) - base;
      sum += strength * res * pow2(ff / PWRfl[l]);
    }

    // line term [Np/km]
    const Numeric absl = 0.3183e-4 * den_dummy * sum;
    // pxsec = abs/vmr [1/m] (Rosenkranz model in [Np/km])
    // 4.1907e-5 = 0.230259 * 0.1820 * 1.0e-3    (1/(10*log(e)) = 0.230259)
    propmat_clearsky.Kjj()[s] += vmr * 1.000e-3 * (absl + (con * ff * ff));
  }
}

//! Ported from legacy continua.  Original documentation
//! Oxygen complex at 60 GHz plus mm O2 lines plus O2 continuum
/*!
  REFERENCES FOR EQUATIONS AND COEFFICIENTS:
  P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
  BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
  H.J. Liebe et al, JQSRT V.48, PP.629-643 (1992).
  M.J. Schwartz, Ph.D. thesis, M.I.T. (1997).
  SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
  This version differs from Liebe's MPM92 in two significant respects:
  1. It uses the modification of the 1- line width temperature dependence
  recommended by Schwartz: (1/T).
  2. It uses the same temperature dependence (X) for submillimeter
  line widths as in the 60 GHz band: (1/T)**0.8

  history:
  05-01-95  P. Rosenkranz
  11-05-97  P. Rosenkranz - 1- line modification.
  12-16-98  pwr - updated submm freq's and intensities from HITRAN96

   \param[out] pxsec        cross section (absorption/volume mixing ratio) of
                            O2 according to the P. W. Rosenkranz, 1993 [1/m]
   \param    CCin           O2-continuum scale factor  [1]
   \param    CLin           O2 line strength scale factor [1]
   \param    CWin           O2 line broadening scale factor [1]
   \param    COin           O2 line coupling scale factor [1]
   \param    model          allows user defined input parameter set
                            (CCin, CLin, CWin, and COin)<br> or choice of
                            pre-defined parameters of specific models (see note below).
   \param    version        determines model version: 1988, 1993, 1998
   \param    f_grid         predefined frequency grid        [Hz]
   \param    abs_p          predefined pressure              [Pa]
   \param    abs_t          predefined temperature grid      [K]
   \param    vmrh2o         H2O volume mixing ratio profile  [1]
   \param    vmr            O2 volume mixing ratio profile   [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, CWin, and COin
             are neglected (model dominates over parameters).<br>
             Allowed models:<br>
             'Rosenkranz', 'RosenkranzLines', 'RosenkranzContinuum',
             'RosenkranzNoCoupling', and 'user'. <br>
       For the parameter  version the following three string values are allowed:
       'PWR88', 'PWR93', 'PWR98'.<br>
             See the user guide for detailed explanations.

   \remark   Reference:  P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.

   \author Thomas Kuhn
   \date 2001-11-05
 */
//! New implementation
void oxygen(PropagationMatrix& propmat_clearsky,
            const Vector& f_grid,
            const Numeric p_pa,
            const Numeric t,
            const Numeric vmr,
            const Numeric h2o) {
  using Math::pow2;
  using Math::pow3;
  constexpr Numeric VMRCalcLimit = 1.000e-25;

  // widths in MHz/mbar for the O2 continuum
  constexpr Numeric WB300 = 0.56;  // [MHz/mbar]=[MHz/hPa]
  constexpr Numeric X = 0.80;      // [1]

  // intensities in the submm range are updated according to HITRAN96
  static constexpr std::array F{
      118.7503, 56.2648,  62.4863,  58.4466,  60.3061, 59.5910, 59.1642,
      60.4348,  58.3239,  61.1506,  57.6125,  61.8002, 56.9682, 62.4112,
      56.3634,  62.9980,  55.7838,  63.5685,  55.2214, 64.1278, 54.6712,
      64.6789,  54.1300,  65.2241,  53.5957,  65.7648, 53.0669, 66.3021,
      52.5424,  66.8368,  52.0214,  67.3696,  51.5034, 67.9009, 368.4984,
      424.7632, 487.2494, 715.3931, 773.8397, 834.1458};

  // intensities in the submm range are updated according to HITRAN96
  static constexpr std::array S300{
      0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14, 0.3351E-14, 0.3292E-14,
      0.3721E-14, 0.3891E-14, 0.3640E-14, 0.4005E-14, 0.3227E-14, 0.3715E-14,
      0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14, 0.1391E-14, 0.1808E-14,
      0.9124E-15, 0.1230E-14, 0.5603E-15, 0.7842E-15, 0.3228E-15, 0.4689E-15,
      0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15, 0.4264E-16, 0.6899E-16,
      0.1924E-16, 0.3229E-16, 0.8191E-17, 0.1423E-16, 0.6494E-15, 0.7083E-14,
      0.3025E-14, 0.1835E-14, 0.1158E-13, 0.3993E-14};

  // y parameter for the calculation of Y [1/bar]
  static constexpr std::array Y300{
      -0.0233, 0.2408,  -0.3486, 0.5227,  -0.5430, 0.5877,  -0.3970, 0.3237,
      -0.1348, 0.0311,  0.0725,  -0.1663, 0.2832,  -0.3629, 0.3970,  -0.4599,
      0.4695,  -0.5199, 0.5187,  -0.5597, 0.5903,  -0.6246, 0.6656,  -0.6942,
      0.7086,  -0.7325, 0.7348,  -0.7546, 0.7702,  -0.7864, 0.8083,  -0.8210,
      0.8439,  -0.8529, 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000};

  // line width parameter [GHz/bar]
  static constexpr std::array W300{
      1.630, 1.646, 1.468, 1.449, 1.382, 1.360, 1.319, 1.297, 1.266, 1.248,
      1.221, 1.207, 1.181, 1.171, 1.144, 1.139, 1.110, 1.108, 1.079, 1.078,
      1.050, 1.050, 1.020, 1.020, 1.000, 1.000, 0.970, 0.970, 0.940, 0.940,
      0.920, 0.920, 0.890, 0.890, 1.920, 1.920, 1.920, 1.810, 1.810, 1.810};

  // temperature exponent of the line strength in [1]
  static constexpr std::array BE{
      0.009, 0.015, 0.083, 0.084, 0.212, 0.212, 0.391, 0.391, 0.626, 0.626,
      0.915, 0.915, 1.260, 1.260, 1.660, 1.665, 2.119, 2.115, 2.624, 2.625,
      3.194, 3.194, 3.814, 3.814, 4.484, 4.484, 5.224, 5.224, 6.004, 6.004,
      6.844, 6.844, 7.744, 7.744, 0.048, 0.044, 0.049, 0.145, 0.141, 0.145};

  // v parameter for the calculation of Y [1/bar]
  static constexpr std::array V{
      0.0079, -0.0978, 0.0844, -0.1273, 0.0699, -0.0776, 0.2309, -0.2825,
      0.0436, -0.0584, 0.6056, -0.6619, 0.6451, -0.6759, 0.6547, -0.6675,
      0.6135, -0.6139, 0.2952, -0.2895, 0.2654, -0.2590, 0.3750, -0.3680,
      0.5085, -0.5002, 0.6206, -0.6091, 0.6526, -0.6393, 0.6640, -0.6475,
      0.6729, -0.6545, 0.0000, 0.0000,  0.0000, 0.0000,  0.0000, 0.0000};

  // Loop pressure/temperature:
  // check if O2-VMR is exactly zero (caused by zeropadding), then return 0.
  if (vmr == 0.) {
    return;
  }

  // check if O2-VMR will cause an underflow due to division by zero:
  ARTS_USER_ERROR_IF(
      vmr < VMRCalcLimit,
      "ERROR: PWR98 O2 full absorption model has detected a O2 volume mixing ratio of ",
      vmr,
      " which is below the threshold of ",
      VMRCalcLimit,
      ".\n"
      "Therefore no calculation is performed.\n")

  // relative inverse temperature [1]
  const Numeric TH = 3.0000e2 / t;
  const Numeric TH1 = (TH - 1.000e0);
  const Numeric B = pow(TH, X);
  // partial pressure of H2O and dry air [hPa]
  const Numeric PRESWV = 1e-2 * (p_pa * h2o);
  const Numeric PRESDA = 1e-2 * (p_pa * (1.000e0 - h2o));
  const Numeric DEN = 0.001 * (PRESDA * B + 1.1 * PRESWV * TH);  // [hPa]
  const Numeric DENS = 0.001 * (PRESDA + 1.1 * PRESWV) * TH;     // [hPa]
  const Numeric DFNR = WB300 * DEN;                              // [GHz]

  // continuum absorption [1/m/GHz]
  const Numeric CCONT = 1.23e-10 * pow2(TH) * p_pa;

  // Loop over input frequency
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // initial O2 line absorption at position ff
    // Numeric O2ABS  = 0.000e0;cd safff

    // input frequency in [GHz]
    const Numeric ff = 1e-9 * f_grid[s];

    // continuum absorption [Neper/km]
    const Numeric CONT = CCONT * (ff * ff * DFNR / (ff * ff + DFNR * DFNR));

    // Loop over Rosnekranz '93 spectral line frequency:
    Numeric SUM = 0.000e0;
    for (std::size_t l = 0; l < W300.size(); ++l) {
      const Numeric DF =
          W300[l] * ((fabs((F[l] - 118.75)) < 0.10) ? DENS : DEN);  // [hPa]
      // 118 line update according to M. J. Schwartz, MIT, 1997

      const Numeric Y = 0.001 * 0.01 * p_pa * B * (Y300[l] + V[l] * TH1);
      const Numeric STR = S300[l] * exp(-BE[l] * TH1);
      const Numeric SF1 =
          (DF + (ff - F[l]) * Y) / ((ff - F[l]) * (ff - F[l]) + DF * DF);
      const Numeric SF2 =
          (DF - (ff + F[l]) * Y) / ((ff + F[l]) * (ff + F[l]) + DF * DF);
      SUM += STR * (SF1 + SF2) * (ff / F[l]) * (ff / F[l]);
    }

    // O2 absorption [Neper/km]
    // Rosenkranz uses the factor 0.5034e12 in the calculation of the abs coeff.
    // This factor is the product of several terms:
    // 0.5034e12 = ISORATIO *   VMR   * (Hz/GHz) * (k_B*300K)^-1
    //           = 0.995262 * 0.20946 *   10^-9  * 2.414322e21(hPa*cm^2*km)^-1
    //             |---- 0.2085 ----|   |---- 2.414322e12(hPa*cm^2*km)^-1 ---|
    //             |---- 0.2085 ----|   |---- 2.414322e10( Pa*cm^2*km)^-1 ---|
    // O2ABS = 2.4143e12 * SUM * PRESDA * pow(TH, 3.0) / PI;
    // O2ABS = CONT + (2.414322e10 * SUM * abs_p[i] * pow(TH, 3.0) / PI);
    // unit conversion x Nepers/km = y 1/m  --->  y = x * 1.000e-3
    // therefore 2.414322e10 --> 2.414322e7
    // pxsec [1/m]
    propmat_clearsky.Kjj()[s] +=
        vmr * (CONT + (2.414322e7 * SUM * p_pa * pow3(TH) / Constant::pi));
  }
}
}  // namespace Absorption::PredefinedModel::PWR98
