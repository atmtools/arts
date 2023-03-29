#include <propagationmatrix.h>

#include <array>
#include <cstddef>

#include "arts_constants.h"
#include "debug.h"

namespace Absorption::PredefinedModel::TRE05 {
//
// #################################################################################
//
/**

   \retval   MPMLineShapeO2Function  O2-line shape function value         [1]
   \param    gamma                   O2-line width                        [Hz]
   \param    fl                      H2O-line central frequency of the    [Hz]
   \param    f                       frequency position of calculation    [Hz]
   \param    delta                   O2-line mixing parameter             [1]

   \note     This function calculates the line shape function of Van Vleck and Weisskopf
             for O2 with line mixing.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */

constexpr Numeric MPMLineShapeO2Function(const Numeric gamma,
                                         const Numeric fl,
                                         const Numeric f,
                                         const Numeric delta) noexcept {
  /*
    this routine calculates the line shape function of Van Vleck and Weisskopf
    for O2 with line mixing.

    creation  TKS, 14.07.01

    input:   gamma   [GHz]    line width of line L
             fl      [GHz]    central frequency of line L
             f       [GHz]    frequency position of calculation
             delta   [1]     line mixing parameter

    output:  value   [1]     line shape function value at f for the line parameters
                             of line L

   */

  double f_minus, f_plus; /* internal variables */
  double value;           /* return value       */

  // line at fl
  f_minus = (gamma - delta * (fl - f)) / ((fl - f) * (fl - f) + gamma * gamma);

  // mirror line at -fl
  f_plus = (gamma - delta * (fl + f)) / ((fl + f) * (fl + f) + gamma * gamma);

  // VVW line shape function value
  value = f * (f_minus + f_plus);

  return value;
}

// #################################################################################
//! TRE05O2AbsModel
/*!
 *  \param[out] pxsec        cross section (absorption/volume mixing ratio) of
 *                           O2 according to TRE05 [1/m]
 *  \param    CCin           scaling factor for the O2-continuum   [1]
 *  \param    CLin           scaling factor for the O2-line strengths [1]
 *  \param    CWin           scaling factor for the O2-line widths    [1]
 *  \param    COin           scaling factor for the O2-line coupling  [1]
 *  \param    model          allows user defined input parameter set
 *                           (CCin, CLin, CWin, and COin)<br> or choice of
 *                           pre-defined parameters of specific models (see note below).
 *  \param    f_grid         predefined frequency grid           [Hz]
 *  \param    abs_p          predefined pressure                 [Pa]
 *  \param    abs_t          predefined temperature grid         [K]
 *  \param    abs_h2o        H2O volume mixing ratio profile    [1]
 *  \param    vmr            O2 volume mixing ratio profile     [1]
 * 
 *  \note     Except for  model 'user' the input parameters CCin, CLin, CWin, and COin
 *            are neglected (model dominates over parameters).<br>
 *            Allowed models: 'TRE05', 'TRE05Lines', 'TRE05Continuum', 'TRE05NoCoupling',
 *            'TRE05NoCutoff', and 'user'. See the user guide for detailed explanations.
 * 
 *  \remark   References: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
 *            <i>Propagation modeling of moist air and suspended water/ice
 *            particles at frequencies below 1000 GHz</i>,<br>
 *            AGARD 52nd Specialists Meeting of the Electromagnetic Wave
 *            Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21
 *  
 *            M.Yu. Tretyakov, M.A. Koshelev, V.V. Dorovskikh,
 *            D.S. Makarov, P.W. Rosenkranz; 60-GHz oxygen band: precise broadening and central frequencies
 *            of fine-structure lines, absolute absorption profile
 *            at atmospheric pressure, and revision of mixing coefficients
 *            doi:10.1016/j.jms.2004.11.011
 * 
 * \remark    This is a copy of MPM93O2AbsModel with an exception of having new values from
 *            the Tretyakov etal. 2005 paper.
 * 
 *  \author Richard Larsson
 *  \date 2013-09-20
 */

void oxygen(PropagationMatrix& propmat_clearsky,
            const Vector& f_grid,
            const Numeric p_pa,
            const Numeric t,
            const Numeric oxygen_vmr,
            const Numeric water_vmr) {
  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0             1           2         3         4      5        6
  //         f0           a1           a2       a3        a4      a5       a6
  //        [GHz]      [kHz/hPa]      [1]    [MHz/hPa]    [1]       [10³/hPa]
  static constexpr std::array tre05{
      std::array{
          50.474214, 0.975 / 10, 9.651, 0.669, 0.0, 0.2566, 0.685},  // 37-
      std::array{
          50.987745, 2.529 / 10, 8.653, 0.717, 0.0, 0.2246, 0.680},  // 35-
      std::array{
          51.503360, 6.193 / 10, 7.709, 0.764, 0.0, 0.1947, 0.6729},  // 33-
      std::array{
          52.021429, 14.32 / 10, 6.819, 0.811, 0.0, 0.1667, 0.6640},  // 31-
      std::array{
          52.542418, 31.24 / 10, 5.983, 0.858, 0.0, 0.1388, 0.6526},  // 29-
      std::array{
          53.066934, 64.29 / 10, 5.201, 0.906, 0.0, 0.1349, 0.6206},  // 27-
      std::array{
          53.595775, 124.6 / 10, 4.474, 0.955, 0.0, 0.2227, 0.5085},  // 25-
      std::array{
          54.130025, 227.3 / 10, 3.800, 0.996, 0.0, 0.3170, 0.3750},  // 23-
      std::array{
          54.671180, 389.7 / 10, 3.182, 1.037, 0.0, 0.3558, 0.2654},  // 21-
      std::array{
          55.221384, 627.1 / 10, 2.618, 1.089, 0.0, 0.2560, 0.2952},  // 19-
      std::array{
          55.783815, 945.3 / 10, 2.109, 1.134, 0.0, -0.1172, 0.6135},  // 17-
      std::array{
          56.264774, 543.4 / 10, 0.014, 1.703, 0.0, 0.3525, -0.0978},  //  1+
      std::array{
          56.363399, 1331.8 / 10, 1.654, 1.189, 0.0, -0.2378, 0.6547},  // 15-
      std::array{
          56.968211, 1746.6 / 10, 1.255, 1.223, 0.0, -0.3545, 0.6451},  // 13-
      std::array{
          57.612486, 2120.1 / 10, 0.910, 1.262, 0.0, -0.5416, 0.6056},  // 11-
      std::array{
          58.323877, 2363.7 / 10, 0.621, 1.295, 0.0, -0.1932, 0.0436},  //  9-
      std::array{
          58.446588, 1442.1 / 10, 0.083, 1.491, 0.0, 0.6768, -0.1273},  //  3+
      std::array{
          59.164204, 2379.9 / 10, 0.387, 1.353, 0.0, -0.6561, 0.2309},  //  7-
      std::array{
          59.590983, 2090.7 / 10, 0.207, 1.408, 0.0, 0.6957, -0.0776},  //  5+
      std::array{
          60.306056, 2103.4 / 10, 0.207, 1.415, 0.0, -0.6395, 0.0699},  //  5-
      std::array{
          60.434778, 2438.0 / 10, 0.386, 1.339, 0.0, 0.6342, -0.2825},  //  7+
      std::array{
          61.150562, 2479.5 / 10, 0.621, 1.292, 0.0, 0.1014, -0.0584},  //  9+
      std::array{
          61.800158, 2275.9 / 10, 0.910, 1.263, 0.0, 0.5014, -0.6619},  // 11+
      std::array{
          62.411220, 1915.4 / 10, 1.255, 1.217, 0.0, 0.3029, -0.6759},  // 13+
      std::array{
          62.486253, 1503.0 / 10, 0.083, 1.513, 0.0, -0.4499, 0.0844},  //  3-
      std::array{
          62.997984, 1490.2 / 10, 1.654, 1.174, 0.0, 0.1856, -0.6675},  // 15+
      std::array{
          63.568526, 1078.0 / 10, 2.108, 1.134, 0.0, 0.0658, -0.6139},  // 17+
      std::array{
          64.127775, 728.7 / 10, 2.617, 1.088, 0.0, -0.3036, -0.2895},  // 19+
      std::array{
          64.678910, 461.3 / 10, 3.181, 1.038, 0.0, -0.3968, -0.2590},  // 21+
      std::array{
          65.224078, 274.0 / 10, 3.800, 0.996, 0.0, -0.3528, -0.3680},  // 23+
      std::array{
          65.764779, 153.0 / 10, 4.473, 0.955, 0.0, -0.2548, -0.5002},  // 25+
      std::array{
          66.302096, 80.40 / 10, 5.200, 0.906, 0.0, -0.1660, -0.6091},  // 27+
      std::array{
          66.836834, 39.80 / 10, 5.982, 0.858, 0.0, -0.1680, -0.6393},  // 29+
      std::array{
          67.369601, 18.56 / 10, 6.818, 0.811, 0.0, -0.1956, -0.6475},  // 31+
      std::array{
          67.900868, 8.172 / 10, 7.708, 0.764, 0.0, -0.2216, -0.6545},  // 33+
      std::array{
          68.431006, 3.397 / 10, 8.652, 0.717, 0.0, -0.2492, -0.660},  // 35+
      std::array{
          68.960312, 1.334 / 10, 9.650, 0.669, 0.0, -0.2773, -0.665},  // 37+
      std::array{
          118.750334, 940.3 / 10, 0.010, 1.664, 0.0, -0.0439, 0.0079},  //  1-
      std::array{368.498246, 67.4 / 10, 0.048, 1.64, 0.0, 0.0, 0.0},    // QN1
      std::array{424.763020, 637.7 / 10, 0.044, 1.64, 0.0, 0.0, 0.0},   // QN2
      std::array{487.249273, 237.4 / 10, 0.049, 1.60, 0.0, 0.0, 0.0},   // QN3
      std::array{715.392902, 98.1 / 10, 0.145, 1.60, 0.0, 0.0, 0.0},    // QN4
      std::array{773.839490, 572.3 / 10, 0.141, 1.62, 0.0, 0.0, 0.0},   // QN5
      std::array{834.145546, 183.1 / 10, 0.145, 1.47, 0.0, 0.0, 0.0}};  // QN6

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we
  // have top devide the line strength by this value since arts multiplies pxsec
  // by these variables later in abs_coefCalc.
  constexpr Numeric VMRISO = 0.2085;

  // O2 continuum parameters of TRE05:
  constexpr Numeric S0 =
      6.140e-5;  // line strength                        [ppm]
  constexpr Numeric G0 =
      0.560e-3;  // line width                           [GHz/hPa]
  constexpr Numeric X0 = 0.800;  // temperature dependence of line width [1]

  if (oxygen_vmr == 0.) return;

  constexpr Numeric VMRCalcLimit = 1.000e-25;
  // check if O2-VMR will cause an underflow due to division by zero:
  ARTS_USER_ERROR_IF(
      oxygen_vmr < VMRCalcLimit,
      "ERROR: TRE05 O2 full absorption model has detected a O2 volume mixing ratio of ",
      oxygen_vmr,
      " which is below the threshold of ",
      VMRCalcLimit,
      ".\n"
      "Therefore no calculation is performed.\n")
  constexpr Numeric dB_km_to_1_m =
      (1.00000e-3 / (10.0 * Constant::log10_euler));
  constexpr Numeric Hz_to_GHz = 1.000000e-9;  // [GHz/Hz]
  constexpr Numeric Pa_to_hPa = 1.000000e-2;  // [Pa/hPa]

  // relative inverse temperature [1]
  const Numeric theta = (300.0 / t);
  // H2O partial pressure [hPa]
  const Numeric pwv = Pa_to_hPa * p_pa * water_vmr;
  // dry air partial pressure [hPa]
  const Numeric pda = (Pa_to_hPa * p_pa) - pwv;
  // here the total pressure is devided by the O2 vmr for the
  // P_dry calculation because we calculate pxsec and not abs: abs = vmr * pxsec
  // old version without VMRISO: Numeric pda_dummy = pda / oxygen_vmr;
  const Numeric pda_dummy = pda;
  // O2 continuum strength [ppm]
  const Numeric strength_cont = S0 * pda_dummy * pow(theta, (Numeric)2.);
  // O2 continuum pseudo line broadening [GHz]
  const Numeric gam_cont = G0 * (pwv + pda) * pow(theta, X0);  // GHz

  // Loop over input frequency
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // input frequency in [GHz]
    const Numeric ff = f_grid[s] * Hz_to_GHz;
    // O2 continuum absorption [1/m]
    // cross section: pxsec = absorption / var
    // the vmr of O2 will be multiplied at the stage of absorption calculation:
    const Numeric Nppc = strength_cont * ff * gam_cont /
                         (pow(ff, (Numeric)2.) + pow(gam_cont, (Numeric)2.));

    // Loop over TRE05 O2 spectral lines:
    Numeric Nppl = 0.0;
    for (const auto& l : tre05) {
      // line strength [ppm]   S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
      const Numeric strength = 1.000e-6 * pda_dummy * l[1] / l[0] *
                               pow(theta, (Numeric)3.) *
                               exp(l[2] * (1.0 - theta));
      // line broadening parameter [GHz]
      const Numeric gam =
          (l[3] * 0.001 *
           ((pda * pow(theta, ((Numeric)0.8 - l[4]))) + (1.10 * pwv * theta)));
      // line mixing parameter [1]
      //      if (l < 11) CD = 1.1000;
      const Numeric delta = ((l[5] + l[6] * theta) * (pda + pwv) *
                             pow(theta, (Numeric)0.8) * (Numeric)0.001);
      // absorption [dB/km] like in the original TRE05
      Nppl += strength * MPMLineShapeO2Function(gam, l[0], ff, delta);
    }

    // in TRE05 there is a cutoff for O2 line absorption if abs_l < 0
    // absorption cannot be less than 0 according to TRE05 philosophy.
    if (Nppl < 0.000) Nppl = 0.0000;  // <---!!IMPORTANT FEATURE!!

    //
    // O2 line absorption [1/m]
    // cross section: pxsec = absorption / var
    // the vmr of O2 will be multiplied at the stage of absorption calculation:
    propmat_clearsky.Kjj()[s] +=
        oxygen_vmr * dB_km_to_1_m * 0.1820 * ff * (Nppl + Nppc) / VMRISO;
  }
}
}  // namespace Absorption::PredefinedModel::TRE05
