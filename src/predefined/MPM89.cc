#include <arts_constexpr_math.h>
#include <propagationmatrix.h>

#include <array>
#include <numeric>

#include "arts_constants.h"

namespace Absorption::PredefinedModel::MPM89 {
constexpr Numeric MPMLineShapeFunction(const Numeric gamma,
                                       const Numeric fl,
                                       const Numeric f) noexcept {
  /*
    this routine calculates the line shape function of Van Vleck and Weisskopf
    with the factor (f/f_o)¹. for the MPM pseudo continuum line.

    creation  TKS, 4.11.00

    input:   gamma   [Hz]    line width of line L
             fl      [Hz]    central frequency of line L
             f       [Hz]    frequency position of calculation

    output:  value   [1/Hz]  line shape function value at f for the line parameters
                              of line L

   */

  double f_minus, f_plus; /* internal variables */
  double value;           /* return value       */

  // line at fl
  f_minus = 1.000 / ((f - fl) * (f - fl) + gamma * gamma);

  // mirror line at -fl
  f_plus = 1.000 / ((f + fl) * (f + fl) + gamma * gamma);

  // VVW line shape function value
  value = fabs(f / fl) * gamma * (f_minus + f_plus);

  return value;
}

void water(PropagationMatrix& propmat_clearsky,
           const Vector& f_grid,
           const Numeric p_pa,
           const Numeric t,
           const Numeric vmr) noexcept {
  using Math::pow3;
  constexpr Numeric dB_km_to_1_m = (1e-3 / (10.0 * Constant::log10_euler));

  //
  // Coefficients are from Liebe, Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631
  //         0           1        2       3        4      5      6
  //         f0          b1       b2      b3       b4     b5     b6
  //        [GHz]     [kHz/kPa]   [1]   [MHz/kPa]  [1]    [1]    [1]
  static constexpr std::array<std::array<Numeric, 7>, 30> mpm89 = {
      std::array{22.235080, 0.1090, 2.143, 28.11, 0.69, 4.80, 1.00},
      std::array{67.813960, 0.0011, 8.735, 28.58, 0.69, 4.93, 0.82},
      std::array{119.995940, 0.0007, 8.356, 29.48, 0.70, 4.78, 0.79},
      std::array{183.310074, 2.3000, 0.668, 28.13, 0.64, 5.30, 0.85},
      std::array{321.225644, 0.0464, 6.181, 23.03, 0.67, 4.69, 0.54},
      std::array{325.152919, 1.5400, 1.540, 27.83, 0.68, 4.85, 0.74},
      std::array{336.187000, 0.0010, 9.829, 26.93, 0.69, 4.74, 0.61},
      std::array{380.197372, 11.9000, 1.048, 28.73, 0.69, 5.38, 0.84},
      std::array{390.134508, 0.0044, 7.350, 21.52, 0.63, 4.81, 0.55},
      std::array{437.346667, 0.0637, 5.050, 18.45, 0.60, 4.23, 0.48},
      std::array{439.150812, 0.9210, 3.596, 21.00, 0.63, 4.29, 0.52},
      std::array{443.018295, 0.1940, 5.050, 18.60, 0.60, 4.23, 0.50},
      std::array{448.001075, 10.6000, 1.405, 26.32, 0.66, 4.84, 0.67},
      std::array{470.888947, 0.3300, 3.599, 21.52, 0.66, 4.57, 0.65},
      std::array{474.689127, 1.2800, 2.381, 23.55, 0.65, 4.65, 0.64},
      std::array{488.491133, 0.2530, 2.853, 26.02, 0.69, 5.04, 0.72},
      std::array{503.568532, 0.0374, 6.733, 16.12, 0.61, 3.98, 0.43},
      std::array{504.482692, 0.0125, 6.733, 16.12, 0.61, 4.01, 0.45},
      std::array{556.936002, 510.0000, 0.159, 32.10, 0.69, 4.11, 1.00},
      std::array{620.700807, 5.0900, 2.200, 24.38, 0.71, 4.68, 0.68},
      std::array{658.006500, 0.2740, 7.820, 32.10, 0.69, 4.14, 1.00},
      std::array{752.033227, 250.0000, 0.396, 30.60, 0.68, 4.09, 0.84},
      std::array{841.073593, 0.0130, 8.180, 15.90, 0.33, 5.76, 0.45},
      std::array{859.865000, 0.1330, 7.989, 30.60, 0.68, 4.09, 0.84},
      std::array{899.407000, 0.0550, 7.917, 29.85, 0.68, 4.53, 0.90},
      std::array{902.555000, 0.0380, 8.432, 28.65, 0.70, 5.10, 0.95},
      std::array{906.205524, 0.1830, 5.111, 24.08, 0.70, 4.70, 0.53},
      std::array{916.171582, 8.5600, 1.442, 26.70, 0.70, 4.78, 0.78},
      std::array{970.315022, 9.1600, 1.920, 25.50, 0.64, 4.94, 0.67},
      std::array{987.926764, 138.0000, 0.258, 29.85, 0.68, 4.55, 0.90}};

  // here the total pressure is not multiplied by the H2O vmr for the
  // P_H2O calculation because we calculate pxsec and not abs: abs = vmr * pxsec
  const Numeric pwv_dummy = 1e-3 * p_pa;
  // relative inverse temperature [1]
  const Numeric theta = (300.0 / t);
  // H2O partial pressure [kPa]
  const Numeric pwv = pwv_dummy * vmr;
  // dry air partial pressure [kPa]
  const Numeric pda = pwv_dummy - pwv;
  // H2O continuum absorption [dB/km/GHz^2] like in the original MPM89
  const Numeric Nppc = pwv_dummy * pow3(theta) * 1.000e-5 *
                       ((0.113 * pda) + (3.57 * pwv * pow(theta, 7.5)));

  // Loop over input frequency
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // input frequency in [GHz]
    const Numeric ff = f_grid[s] * 1e-9;
    // H2O line contribution at position f
    const Numeric Nppl = std::transform_reduce(
        mpm89.begin(),
        mpm89.end(),
        0.0,
        std::plus{},
        [pwv_dummy, theta, pwv, pda, ff](auto& l) {
          const Numeric strength =
              pwv_dummy * l[1] * pow(theta, 3.5) * exp(l[2] * (1.000 - theta));
          // line broadening parameter [GHz]
          const Numeric gam =
              l[3] * 0.001 *
              (l[5] * pwv * pow(theta, l[6]) + pda * pow(theta, l[4]));
          return strength * MPMLineShapeFunction(gam, l[0], ff);
        });

    //
    // H2O line absorption [1/m]
    propmat_clearsky.Kjj()[s] +=
        vmr * dB_km_to_1_m * 0.1820 * ff * (Nppl + (Nppc * ff));
  }
}

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

void oxygen(PropagationMatrix& propmat_clearsky,
            const Vector& f_grid,
            const Numeric p_pa,
            const Numeric t,
            const Numeric vmr,
            const Numeric h2o) {
  using Math::pow2;
  using Math::pow3;
  constexpr Numeric VMRCalcLimit = 1.000e-25;
  constexpr Numeric dB_km_to_1_m = (1e-3 / (10.0 * Constant::log10_euler));

  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0            1           2        3          4        5        6
  //         f0           a1          a2       a3         a4      a5        a6
  //        [GHz]      [kHz/hPa]     [1]    [MHz/hPa]    [1]         [1/kPa]
  static constexpr std::array<std::array<Numeric, 7>, 44> mpm89 = {
      std::array{50.474238, 0.94, 9.694, 8.60, 0.0, 1.600, 5.520},        // 0
      std::array{50.987749, 2.46, 8.694, 8.70, 0.0, 1.400, 5.520},        // 1
      std::array{51.503350, 6.08, 7.744, 8.90, 0.0, 1.165, 5.520},        // 2
      std::array{52.021410, 14.14, 6.844, 9.20, 0.0, 0.883, 5.520},       // 3
      std::array{52.542394, 31.02, 6.004, 9.40, 0.0, 0.579, 5.520},       // 4
      std::array{53.066907, 64.10, 5.224, 9.70, 0.0, 0.252, 5.520},       // 5
      std::array{53.595749, 124.70, 4.484, 10.00, 0.0, -0.066, 5.520},    // 6
      std::array{54.130000, 228.00, 3.814, 10.20, 0.0, -0.314, 5.520},    // 7
      std::array{54.671159, 391.80, 3.194, 10.50, 0.0, -0.706, 5.520},    // 8
      std::array{55.221367, 631.60, 2.624, 10.79, 0.0, -1.151, 5.514},    // 9
      std::array{55.783802, 953.50, 2.119, 11.10, 0.0, -0.920, 5.025},    // 10
      std::array{56.264775, 548.90, 0.015, 16.46, 0.0, 2.881, -0.069},    // 11
      std::array{56.363389, 1344.00, 1.660, 11.44, 0.0, -0.596, 4.750},   // 12
      std::array{56.968206, 1763.00, 1.260, 11.81, 0.0, -0.556, 4.104},   // 13
      std::array{57.612484, 2141.00, 0.915, 12.21, 0.0, -2.414, 3.536},   // 14
      std::array{58.323877, 2386.00, 0.626, 12.66, 0.0, -2.635, 2.686},   // 15
      std::array{58.446590, 1457.00, 0.084, 14.49, 0.0, 6.848, -0.647},   // 16
      std::array{59.164207, 2404.00, 0.391, 13.19, 0.0, -6.032, 1.858},   // 17
      std::array{59.590983, 2112.00, 0.212, 13.60, 0.0, 8.266, -1.413},   // 18
      std::array{60.306061, 2124.00, 0.212, 13.82, 0.0, -7.170, 0.916},   // 19
      std::array{60.434776, 2461.00, 0.391, 12.97, 0.0, 5.664, -2.323},   // 20
      std::array{61.150560, 2504.00, 0.626, 12.48, 0.0, 1.731, -3.039},   // 21
      std::array{61.800154, 2298.00, 0.915, 12.07, 0.0, 1.738, -3.797},   // 22
      std::array{62.411215, 1933.00, 1.260, 11.71, 0.0, -0.048, -4.277},  // 23
      std::array{62.486260, 1517.00, 0.083, 14.68, 0.0, -4.290, 0.238},   // 24
      std::array{62.997977, 1503.00, 1.665, 11.39, 0.0, 0.134, -4.860},   // 25
      std::array{63.568518, 1087.00, 2.115, 11.08, 0.0, 0.541, -5.079},   // 26
      std::array{64.127767, 733.50, 2.620, 10.78, 0.0, 0.814, -5.525},    // 27
      std::array{64.678903, 463.50, 3.195, 10.50, 0.0, 0.415, -5.520},    // 28
      std::array{65.224071, 274.80, 3.815, 10.20, 0.0, 0.069, -5.520},    // 29
      std::array{65.764772, 153.00, 4.485, 10.00, 0.0, -0.143, -5.520},   // 30
      std::array{66.302091, 80.09, 5.225, 9.70, 0.0, -0.428, -5.520},     // 31
      std::array{66.836830, 39.46, 6.005, 9.40, 0.0, -0.726, -5.520},     // 32
      std::array{67.369598, 18.32, 6.845, 9.20, 0.0, -1.002, -5.520},     // 33
      std::array{67.900867, 8.01, 7.745, 8.90, 0.0, -1.255, -5.520},      // 34
      std::array{68.431005, 3.30, 8.695, 8.70, 0.0, -1.500, -5.520},      // 35
      std::array{68.960311, 1.28, 9.695, 8.60, 0.0, -1.700, -5.520},      // 36
      std::array{118.750343, 945.00, 0.009, 16.30, 0.0, -0.247, 0.003},   // 37
      std::array{368.498350, 67.90, 0.049, 19.20, 0.6, 0.000, 0.000},     // 38
      std::array{424.763124, 638.00, 0.044, 19.16, 0.6, 0.000, 0.000},    // 39
      std::array{487.249370, 235.00, 0.049, 19.20, 0.6, 0.000, 0.000},    // 40
      std::array{715.393150, 99.60, 0.145, 18.10, 0.6, 0.000, 0.000},     // 41
      std::array{773.839675, 671.00, 0.130, 18.10, 0.6, 0.000, 0.000},    // 42
      std::array{834.145330, 180.00, 0.147, 18.10, 0.6, 0.000, 0.000}     // 43
  };

  // O2 continuum parameters of MPM92:
  const Numeric S0 = 6.140e-4;  // line strength                        [ppm]
  const Numeric G0 = 5.60e-3;  // line width                           [GHz/kPa]
  const Numeric X0 = 0.800;    // temperature dependence of line width [1]

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we
  // have top devide the line strength by this value since arts multiplies pxsec
  // by these variables later in abs_coefCalc.
  const Numeric VMRISO = 0.2085;

  // check if O2-VMR is exactly zero (caused by zeropadding), then return 0.
  if (vmr == 0.) {
    return;
  }

  // check if O2-VMR will cause an underflow due to division by zero:
  ARTS_USER_ERROR_IF(
      vmr < VMRCalcLimit,
      "ERROR: MPM89 O2 full absorption model has detected a O2 volume mixing ratio of ",
      vmr,
      " which is below the threshold of ",
      VMRCalcLimit,
      ".\n"
      "Therefore no calculation is performed.\n")

  // relative inverse temperature [1]
  const Numeric theta = (300.0 / t);
  // H2O partial pressure [kPa]
  const Numeric pwv = 1e-3 * p_pa * h2o;
  // dry air partial pressure [kPa]
  const Numeric pda = (1e-3 * p_pa) - pwv;
  // here the total pressure is devided by the O2 vmr for the
  // P_dry calculation because we calculate pxsec and not abs: abs = vmr * pxsec
  const Numeric pda_dummy = pda;
  // O2 continuum strength [ppm]
  const Numeric strength_cont = S0 * pda_dummy * pow2(theta);
  // O2 continuum pseudo line broadening [GHz]
  const Numeric gam_cont = G0 * (pwv + pda) * pow(theta, X0);  // GHz

  // Loop over input frequency
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // input frequency in [GHz]
    const Numeric ff = f_grid[s] * 1e-9;
    // O2 continuum absorption [1/m]
    // cross section: pxsec = absorption / var
    // the vmr of O2 will be multiplied at the stage of absorption calculation:
    const Numeric Nppc =
        strength_cont * ff * gam_cont / (pow2(ff) + pow2(gam_cont));

    // Loop over MPM89 O2 spectral lines:
    const Numeric Nppl = std::transform_reduce(
        mpm89.begin(),
        mpm89.end(),
        0.0,
        std::plus{},
        [pda_dummy, theta, pda, pwv, ff](auto& l) {
          // line strength [ppm]   S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
          const Numeric strength = l[1] * 1.000e-6 * pda_dummy * pow3(theta) *
                                   exp(l[2] * (1.000 - theta)) / l[0];
          // line broadening parameter [GHz]
          const Numeric gam = (l[3] * 1.000e-3 *
                               ((pda * pow(theta, ((Numeric)0.80 - l[4]))) +
                                (1.10 * pwv * theta)));
          // line mixing parameter [1]
          const Numeric delta = ((l[5] + l[6] * theta) * 1.000e-3 * pda *
                                 pow(theta, (Numeric)0.8));
          // absorption [dB/km] like in the original MPM92
          return strength * MPMLineShapeO2Function(gam, l[0], ff, delta);
        });

    //
    // O2 line absorption [1/m]
    propmat_clearsky.Kjj()[s] += vmr * dB_km_to_1_m * 0.1820 * ff *
                                 (((Nppl < 0.000) ? 0.0 : Nppl) + Nppc) /
                                 VMRISO;
  }
}
}  // namespace Absorption::PredefinedModel::MPM89
