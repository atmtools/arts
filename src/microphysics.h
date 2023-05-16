/*!
  \file   microphysics.h
  \author Jana Mendrok, Daniel Kreyling, Manfred Brath, Patrick Eriksson
  \date   2017-07-10 
  
  \brief  Internal functions for microphysics calculations (size distributions etc.)
*/

#ifndef microphysics_h
#define microphysics_h

#include "array.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "matpack_data.h"
#include "messages.h"
#include "optproperties.h"
#include "ppath.h"


//! asymmetry_parameter
/*! Calculates the asymmetry parameter, TRO data only

    \param sa_grid   Scattering angle grid   
    \param pfun      Phase function of TRO-type

    \return Assymetry parameter

    \author Patrick Eriksson
    \date 2022-03-06
*/
Numeric asymmetry_parameter(ConstVectorView za_grid,
                            ConstVectorView pfun);

/*! Derives a and b for relationship mass = a * x^b

    The parameters a and b are derived by a fit including all data inside the
    size range [x_fit_start,x_fit_end].

    The vector x must have been checked to have at least 2 elements.

    An error is thrown if less than two data points are found inside
    [x_fit_start,x_fit_end].

    \return a           Derived a parameter.
    \return b           Derived b parameter.
    \param  x           Size grid
    \param  mass        Particle masses
    \param  x_fit_start Start point of x-range to use for fitting
    \param  x_fit_end   Endpoint of x-range to use for fitting
  
  \author Jana Mendrok, Patrick Eriksson
  \date 2017-10-18

*/
void derive_scat_species_a_and_b(Numeric& a,
                                 Numeric& b,
                                 const Vector& x,
                                 const Vector& mass,
                                 const Numeric& x_fit_start,
                                 const Numeric& x_fit_end);

#endif  //microphysics_h
