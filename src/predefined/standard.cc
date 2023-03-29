#include <propagationmatrix.h>

#include "arts_constexpr_math.h"

namespace Absorption::PredefinedModel::Standard {
//! Ported from legacy continua.  Original documentation//! Standard_O2_continuum
/*!
   \param[out] pxsec        cross section (absorption/volume mixing ratio) of
                            O2-continuum according to Rosenkranz 1993 [1/m]
   \param    Cin            O2-continuum coefficient                  [1/(Hz*Pa*m)]
   \param    G0in           line width                                [Hz/Pa]
   \param    G0Ain          dry air broadening parameter              [1]
   \param    G0Bin          water vapor broadening parameter          [1]
   \param    XG0din         line strength temperature exponent        [1]
   \param    XG0win         line widths temperature exponent          [1]
   \param    model          allows user defined input parameter set
                            (S0in, G0in, XS0in, and XG0in)<br> or choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid                 [Hz]
   \param    abs_p          predefined pressure grid                  [Pa]
   \param    abs_t          predefined temperature grid               [K]
   \param    abs_h2o        H2O volume mixing ratio profile           [1]
   \param    vmr            O2 volume mixing ratio profile            [1]

   \note     Except for  model 'user' the input parameters S0in, G0in, XS0in, and XG0in
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz' and 'user'.
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.<br>
             <br>
             Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */
//! New implementation
void oxygen(PropagationMatrix& propmat_clearsky,
            const Vector& f_grid,
            const Numeric p_pa,
            const Numeric t,
            const Numeric o2,
            const Numeric h2o) {
  using Math::pow2;

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // P. W. Rosenkranz, Chapter 2, in M. A. Janssen,
  // Atmospheric Remote Sensing by Microwave Radiometry, John Wiley & Sons, Inc., 1993
  // ftp://mesa.mit.edu/phil/lbl_rt
  constexpr Numeric C =
      (1.108e-14 / pow2(3.0e2));  // [1/(Hz*Pa*m)] line strength
  const Numeric G0 = 5600.000;    // line width [Hz/Pa]
  const Numeric G0A = 1.000;      // line width [1]
  const Numeric G0B = 1.100;      // line width [1]
  const Numeric XG0d = 0.800;     // temperature dependence of line width [1]
  const Numeric XG0w = 1.000;     // temperature dependence of line width [1]

  const Numeric TH = 3.0e2 / t;  // relative temperature  [1]

  const Numeric ph2o = p_pa * h2o;   // water vapor partial pressure [Pa]
  const Numeric pdry = p_pa - ph2o;  // dry air partial pressure     [Pa]

  // pseudo broadening term [Hz]
  const Numeric gamma =
      G0 * (G0A * pdry * pow(TH, XG0d) + G0B * ph2o * pow(TH, XG0w));

  // Loop over frequency grid:
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // division by vmr of O2 is necessary because of the absorption calculation
    // abs = vmr * pxsec.
    propmat_clearsky.Kjj()[s] +=
        o2 * C * p_pa * pow2(TH) *
        (gamma * pow2(f_grid[s]) / (pow2(f_grid[s]) + pow2(gamma)));
  }
}

//! Ported from legacy continua.  Original documentation//! 4) N2-N2
/*!
  P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen,
  "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993

   \param[out] pxsec        cross section (absorption/volume mixing ratio) of
                            N2-continuum according to Rosenkranz, 1993 [1/m]
   \param    Cin            continuum strength [1/m * 1/(Hz*Pa)²]
   \param    xfin           continuum frequency exponent [1]
   \param    xtin           continuum strength temperature exponent [1]
   \param    xpin           continuum strength pressure exponent [1]
   \param    model          allows user defined input parameter set
                            (Cin, xfin, xtin, and xpin)<br> or choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid      [Hz]
   \param    abs_p          predefined pressure grid       [Pa]
   \param    abs_t          predefined temperature grid    [K]
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters Cin, xfin, xtin, and xpin
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', and 'user'.
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.

   \author Thomas Kuhn
   \date 2001-11-05
 */
//! New implementation
void nitrogen(PropagationMatrix& propmat_clearsky,
              const Vector& f_grid,
              const Numeric p_pa,
              const Numeric t,
              const Numeric n2) {
  using std::pow;

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model, Chapter 2, pp 74, in M. A. Janssen,
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
  constexpr Numeric C = 1.05e-38;  // [1/(Pa²*Hz²*m)]
  constexpr Numeric xf = 2.00;     // [1]
  constexpr Numeric xt = 3.55;     // [1]
  constexpr Numeric xp = 2.00;     // [1]
  // ---------------------------------------------------------------------------------------

  // Loop over frequency grid:
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    // The second N2-VMR will be multiplied at the stage of absorption
    // calculation: abs = vmr * pxsec.
    propmat_clearsky.Kjj()[s] +=
        n2 * C *               // strength [1/(m*Hz²Pa²)]
        pow(300.00 / t, xt) *  // T dependence        [1]
        pow(f_grid[s], xf) *   // f dependence    [Hz^xt]
        pow(p_pa, xp) *        // p dependence    [Pa^xp]
        pow(n2, xp - 1);       // last N2-VMR at the stage
                               // of absorption calculation
  }
}

//! Ported from legacy continua.  Original documentation
//! Standard_H2O_foreign_continuum
/*!
   \param[out] pxsec        cross section (absorption/volume mixing ratio) of the
                            H2O-dry air continuum [1/m]
   \param    Cin            constant absorption strength [1/m / (Hz*Pa)²]
   \param    xin            temperature exponent         [1]
   \param    model          allows user defined input parameter set
                            (C and x)<br> or choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid    [Hz]
   \param    abs_t          predefined temperature grid  [K]
   \param    abs_p          predefined pressure          [Pa]
   \param    vmr            H2O volume mixing ratio     [1]

   \note     Except for  model 'user' the input parameters C and x
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', 'CruzPol', 'MPM89', 'MPM87', and 'user'.
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and
             Radio Science, Vol. 34(4), 1025, 1999.

   \author Thomas Kuhn
   \date 2001-08-03
 */
//! New implementation
void water_foreign(PropagationMatrix& propmat_clearsky,
                   const Vector& f_grid,
                   const Numeric p_pa,
                   const Numeric t,
                   const Numeric h2o) {
  using Math::pow2;

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model (Radio Science, 33(4), 919, 1998):
  constexpr Numeric C = 5.43e-35;  //  [1/m / (Hz²*Pa²)]
  constexpr Numeric x = 0.0;       //  [1]

  // Dry air partial pressure: p_dry := p_tot - p_h2o.
  const Numeric pdry = p_pa * (1.000e0 - h2o);
  // Dummy scalar holds everything except the quadratic frequency dependence.
  // The vmr of H2O will be multiplied at the stage of absorption
  // calculation: abs = vmr * pxsec.
  const Numeric dummy = C * pow(300. / t, x + 3) * p_pa * pdry;

  // Loop frequency:
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    propmat_clearsky.Kjj()[s] += h2o * dummy * pow2(f_grid[s]);
  }
}

//! Ported from legacy continua.  Original documentation//! Standard_H2O_self_continuum
/*!
   \param[out] pxsec        cross section (absorption/volume mixing ratio) of the
                            H2O-H2O continuum [1/m]
   \param    Cin            constant absorption strength     [1/m / (Hz*Pa)²]
   \param    xin            temperature exponent of (300/T)  [1]
   \param    model          allows user defined input parameter set
                            (C and x)<br> or choice of
                            pre-defined parameters of specific models (see note below).
   \param    f_grid         predefined frequency grid        [Hz]
   \param    abs_p          predefined pressure grid         [Pa]
   \param    abs_t          predefined temperature grid      [K]
   \param    vmr            H2O volume mixing ratio          [1]

   \note     Except for  model 'user' the input parameters C and x
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', 'CruzPol', 'MPM89', 'MPM87', and 'user'.
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and
             Radio Science, Vol. 34(4), 1025, 1999.
   \author Thomas Kuhn
   \date 2001-11-05
 */
//! New implementation
void water_self(PropagationMatrix& propmat_clearsky,
                const Vector& f_grid,
                const Numeric p_pa,
                const Numeric t,
                const Numeric h2o) {
  using Math::pow2;

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model (Radio Science, 33(4), 919, 1998):
  constexpr Numeric C = 1.796e-33;  //  [1/m / (Hz²*Pa²)]
  constexpr Numeric x = 4.5;        //  [1]
  // Dummy scalar holds everything except the quadratic frequency dependence.
  // The second vmr of H2O will be multiplied at the stage of absorption
  // calculation: abs = vmr * pxsec.
  const Numeric dummy = C * pow(300. / t, x + 3) * pow2(p_pa) * h2o;

  // Loop over frequency grid:
  for (Index s = 0; s < f_grid.nelem(); ++s) {
    propmat_clearsky.Kjj()[s] += h2o * dummy * pow2(f_grid[s]);
    //    cout << "pxsec(" << s << "," << i << "): " << pxsec(s,i) << "\n";
  }
}
}  // namespace Absorption::PredefinedModel::Standard