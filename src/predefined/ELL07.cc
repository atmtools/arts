#include <propagationmatrix.h>

#include <algorithm>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "compare.h"
#include "debug.h"

namespace Absorption::PredefinedModel::ELL07 {
//! Ported from legacy continua.  Original documentation
//! ELL07WaterDropletAbs
/*!
   \param[out] pxsec        cross section (absorption/volume mixing ratio) of
                            water clouds according to ELL07 [1/m]
   \param    model          allows choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid       [Hz]
   \param    abs_p          predefined pressure grid        [Pa]
   \param    abs_t          predefined temperature grid     [K]
   \param    vmr            suspended water droplet density profile [kg/m³]
                            (valid range (Tab.1 of MPM93): 0-0.005)

   \note     Allowed models: 'ELL07'.
             See the user guide for detailed explanations.

   \remark   Reference: W. J. Ellison, <br>
            <i>Permittivity of Pure Water, at Standard Atmospheric Pressure, over the
             Frequency Range 0-25 THz and Temperature Range 0-100C</i>,<br>
             J. Phys. Chem. Ref. Data, Vol. 36, No. 1, 2007

   \author Stuart Fox
   \date 2015-06-03
 */
//! New implementation
void compute(PropagationMatrix& propmat_clearsky,
             const Vector& f_grid,
             const Numeric t,
             const Numeric lwc) {
  using Cmp::gt;
  using Constant::pi;
  using Constant::two_pi;
  using Math::pow2;
  using Math::pow3;
  constexpr Numeric LIQUID_AND_ICE_TREAT_AS_ZERO = 1e-10;
  constexpr Numeric dB_km_to_1_m = (1e-3 / (10.0 * Constant::log10_euler));

  if (lwc < LIQUID_AND_ICE_TREAT_AS_ZERO) {
    return;
  }

  constexpr Numeric m =
      1.00e3;  // specific weight of the droplet,  fixed value:  1.00e3 kg/m3
  constexpr Numeric low_lim_den =
      -LIQUID_AND_ICE_TREAT_AS_ZERO;  // lower limit of suspended droplet particle density vector [kg/m3]
  constexpr Numeric high_lim_den =
      5.00e-3;  // upper limit of suspended droplet particle density vector [kg/m3]

  // ELL07 model parameters - table 2 in Ellison (2007)
  constexpr Numeric a1 = 79.23882;
  constexpr Numeric a2 = 3.815866;
  constexpr Numeric a3 = 1.634967;
  constexpr Numeric tc = 133.1383;
  constexpr Numeric b1 = 0.004300598;
  constexpr Numeric b2 = 0.01117295;
  constexpr Numeric b3 = 0.006841548;
  constexpr Numeric c1 = 1.382264e-13;
  constexpr Numeric c2 = 3.510354e-16;
  constexpr Numeric c3 = 6.30035e-15;
  constexpr Numeric d1 = 652.7648;
  constexpr Numeric d2 = 1249.533;
  constexpr Numeric d3 = 405.5169;
  constexpr Numeric p0 = 0.8379692;
  constexpr Numeric p1 = -0.006118594;
  constexpr Numeric p2 = -0.000012936798;
  constexpr Numeric p3 = 4235901000000.0;
  constexpr Numeric p4 = -14260880000.0;
  constexpr Numeric p5 = 273815700.0;
  constexpr Numeric p6 = -1246943.0;
  constexpr Numeric p7 = 9.618642e-14;
  constexpr Numeric p8 = 1.795786e-16;
  constexpr Numeric p9 = -9.310017E-18;
  constexpr Numeric p10 = 1.655473e-19;
  constexpr Numeric p11 = 0.6165532;
  constexpr Numeric p12 = 0.007238532;
  constexpr Numeric p13 = -0.00009523366;
  constexpr Numeric p14 = 15983170000000.0;
  constexpr Numeric p15 = -74413570000.0;
  constexpr Numeric p16 = 497448000.0;
  constexpr Numeric p17 = 2.882476e-14;
  constexpr Numeric p18 = -3.142118e-16;
  constexpr Numeric p19 = 3.528051e-18;

  // Check limits of suspended water droplet density ("vmr") [kg/m³]
  ARTS_USER_ERROR_IF(((lwc < low_lim_den) || (lwc > high_lim_den)),
                     "ERROR in ELL07WaterDropletAbs:\n"
                     "Valid range is ",
                     low_lim_den,
                     "-",
                     high_lim_den,
                     "kg/m3,\n"
                     "but found a value = ",
                     lwc)

  ARTS_USER_ERROR_IF(std::any_of(f_grid.begin(), f_grid.end(), gt(25e12)),
                     "Liquid cloud absorption model ELL07 only valid at\n"
                     "frequencies up to 25THz. Yours are above.")

  ARTS_USER_ERROR_IF(
      t < 210 or t > 373,
      "Only valid for temperatures 210-373 K, not for your value of ",
      t,
      " K")

  // Temperature in celsius
  const Numeric t_cels = t - 273.15;
  // static permittivity
  const Numeric epsilon_s = 87.9144 - 0.404399 * t_cels -
                            9.58726e-4 * pow2(t_cels) -
                            1.32802e-6 * pow3(t_cels);
  // Model parameters
  const Numeric delta1 = a1 * exp(-b1 * t_cels);
  const Numeric delta2 = a2 * exp(-b2 * t_cels);
  const Numeric delta3 = a3 * exp(-b3 * t_cels);
  const Numeric tau1 = c1 * exp(d1 / (t_cels + tc));
  const Numeric tau2 = c2 * exp(d2 / (t_cels + tc));
  const Numeric tau3 = c3 * exp(d3 / (t_cels + tc));
  const Numeric delta4 = p0 + p1 * t_cels + p2 * pow2(t_cels);
  const Numeric f0 = p3 + p4 * t_cels + p5 * pow2(t_cels) + p6 * pow3(t_cels);
  const Numeric tau4 =
      p7 + p8 * t_cels + p9 * pow2(t_cels) + p10 * pow3(t_cels);
  const Numeric delta5 = p11 + p12 * t_cels + p13 * pow2(t_cels);
  const Numeric f1 = p14 + p15 * t_cels + p16 * pow2(t_cels);
  const Numeric tau5 = p17 + p18 * t_cels + p19 * pow2(t_cels);

  // Loop frequency:
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // real part of the complex permittivity of water (triple-debye + 2 resonances)
    const Numeric Reepsilon =
        epsilon_s -
        pow2((two_pi * f_grid[s])) *
            (pow2(tau1) * delta1 / (1. + pow2(two_pi * f_grid[s] * tau1)) +
             pow2(tau2) * delta2 / (1. + pow2(two_pi * f_grid[s] * tau2)) +
             pow2(tau3) * delta3 / (1. + pow2(two_pi * f_grid[s] * tau3))) -
        pow2(two_pi * tau4) * delta4 / 2. *
            (f_grid[s] * (f0 + f_grid[s]) /
                 (1. + pow2(two_pi * tau4 * (f0 + f_grid[s]))) -
             f_grid[s] * (f0 - f_grid[s]) /
                 (1. + pow2(two_pi * tau4 * (f0 - f_grid[s])))) -
        pow2(two_pi * tau5) * delta5 / 2. *
            (f_grid[s] * (f1 + f_grid[s]) /
                 (1. + pow2(two_pi * tau5 * (f1 + f_grid[s]))) -
             f_grid[s] * (f1 - f_grid[s]) /
                 (1. + pow2(two_pi * tau5 * (f1 - f_grid[s]))));
    // imaginary part of the complex permittivity of water (triple-debye + 2 resonances)
    const Numeric Imepsilon =
        two_pi * f_grid[s] *
            (tau1 * delta1 / (1. + pow2(two_pi * f_grid[s] * tau1)) +
             tau2 * delta2 / (1. + pow2(two_pi * f_grid[s] * tau2)) +
             tau3 * delta3 / (1. + pow2(two_pi * f_grid[s] * tau3))) +
        pi * f_grid[s] * tau4 * delta4 *
            (1. / (1. + pow2(two_pi * tau4 * (f0 + f_grid[s]))) +
             1. / (1. + pow2(two_pi * tau4 * (f0 - f_grid[s])))) +
        pi * f_grid[s] * tau5 * delta5 *
            (1. / (1. + pow2(two_pi * tau5 * (f1 + f_grid[s]))) +
             1. / (1. + pow2(two_pi * tau5 * (f1 - f_grid[s]))));

    // the imaginary part of the complex refractivity of suspended liquid water particle [ppm]
    // In MPM93 w is in g/m³ and m is in g/cm³. Because of the units used in arts,
    // a factor of 1.000e6 must be multiplied with the ratio (w/m):
    // MPM93: (w/m)_MPM93  in   (g/m³)/(g/cm³)
    // arts:  (w/m)_arts   in  (kg/m³)/(kg/m³)
    // =====> (w/m)_MPM93   =   1.0e6 * (w/m)_arts
    // the factor of 1.0e6 is included below in pxsec calculation.
    const Numeric ImNw =
        1.500 / m *
        (3.000 * Imepsilon / (pow2((Reepsilon + 2.000)) + pow2(Imepsilon)));
    // liquid water particle absorption cross section [1/m]
    // The vmr of H2O will be multiplied at the stage of absorption
    // calculation: abs = vmr * pxsec.
    // pxsec = abs/vmr [1/m] but MPM93 is in [dB/km] --> conversion necessary
    propmat_clearsky.Kjj()[s] +=
        lwc * 1.000e6 * dB_km_to_1_m * 0.1820 * (f_grid[s] * 1e-9) * ImNw;
  }
}
}  // namespace Absorption::PredefinedModel::ELL07
