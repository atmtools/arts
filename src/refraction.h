/*===========================================================================
  ===  File description 
  ===========================================================================*/

/*!
  \file   refraction.h
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2003-01-17
  
  \brief  Refraction functions.
  
   This file contains the definition of the functions in refraction.cc.
*/

#ifndef refraction_h
#define refraction_h

#include <workspace.h>

void complex_n_water_liebe93(Matrix& complex_n,
                             const Vector& f_grid,
                             const Numeric& t);

void complex_n_ice_matzler06(Matrix& complex_n,
                             const Vector& f_grid,
                             const Numeric& t);

/** Refractive index of water and steam for the optical and near infrared \n
 *
 * From:
 *  Revised formulation for the Refractive
 *  Index of Water and Steam as a Function of
 *  Wavelength, Temperature and Density
 *  \n
 *  Journal of Physical and Chemical Reference Data 27, 761 (1998);
 *  https://doi.org/10.1063/1.556029 27, 761
 *  \n
 *  see also
 *  http://www.iapws.org/release.html
 *
 * @param[out]n refractive index
 * @param[in] only_valid_range flag if true refractive index is calculated only
 *                             within range of validity. If false no check is made,
 *                             so you at your own risk.
 * @param[in] frequency frequency
 * @param[in] temperature temperature
 * @param[in] density density of water or steam
  */
void refractive_index_water_and_steam_VisNIR(Numeric& n,
                                             const Index& only_valid_range,
                                             const Numeric& frequency,
                                             const Numeric& temperature,
                                             const Numeric& density);

#endif  // refraction_h
