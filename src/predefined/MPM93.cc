#include <rtepack.h>

#include "arts_constants.h"

namespace Absorption::PredefinedModel::MPM93 {
//! Ported from legacy continua. 
/*!
  see publication side of National Telecommunications and Information Administration
  http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
  and ftp side for downloading the MPM93 original source code:
  ftp://ftp.its.bldrdoc.gov/pub/mpm93/

   \param[out] propmat_clearksy absorption coefficient of
                                N2-continuum according to MPM93 [1/m]
   \param    f_grid             predefined frequency grid            [Hz]
   \param    p_pa               pressure                             [Pa]
   \param    t                  temperature                          [K]
   \param    n2                 N2 volume mixing ratio               [1]
   \param    h2o                H2O volume mixing ratio profile      [1]

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Original author Thomas Kuhn, new implementation by Stuart Fox
   \date 2001-11-05
 */
//! New implementation
void nitrogen(PropmatVector& propmat_clearsky,
              const Vector& f_grid,
              const Numeric p_pa,
              const Numeric t,
              const Numeric n2,
              const Numeric h2o) {
  using std::pow;
  using Constant::pi, Constant::speed_of_light;

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 N2 continuum model
  // (AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  constexpr Numeric xT = 3.500;  // temperature exponent [1]
  constexpr Numeric xf = 1.500;  // frequency exponent [1]
  constexpr Numeric gxf =
      9.000 * xf;  // needed for the unit conversion of G
  constexpr Numeric S = 2.296e-31;  // line strength  [1/Pa² * 1/Hz]
  static const Numeric G =
      1.930e-5 *
      pow(10.000, -gxf);  // frequency factor [1/Hz^xf]
  // ---------------------------------------------------------------------------------------

  constexpr Numeric fac = 4.0 * pi / speed_of_light;  //  = 4 * pi / c

  const Numeric th = 300.0 / t;
  const Numeric strength =
        S * pow((p_pa * (1.0000 - h2o)), 2.0) *
        pow(th, xT);

    // Loop frequency:
    for (Index s = 0; s < f_grid.nelem(); ++s) {
      propmat_clearsky[s].A() += n2 * 
                     fac * 
                     strength *               // strength
                     pow(f_grid[s], 2.0) /  (1.000 + G * pow(f_grid[s], xf)) * // frequency dependence
                     n2;  // N2 vmr
    }
}
}  // namespace Absorption::PredefinedModel::MPM93
