/* Copyright (C) 2005-2012 Cory Davis <cory@met.ed.ac.uk>
                            
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   mc_antenna.cc
  \author Cory Davis <cdavis@staffmail.ed.ac.uk>
  \date   2005-12-01 

  \brief  Monte Carlo Antenna implementation
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <sstream>

#include "mc_antenna.h"


//! ran_gaussian

/*!
  Draw a random normal (Gaussian) deviate.
  This has been copied with minor changes from the GSL function gsl_ran_gaussian.
  Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122 
  
  \param rng  Rng random number generator instance
  \param sigma standard deviation parameter for gaussian distribution

  Returns the gaussian random deviate.

  \author Cory Davis
  \date   2003-12-01
  

*/

Numeric ran_gaussian (
                      Rng& rng, 
                      const Numeric sigma)
{
  Numeric x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * rng.draw();
      y = -1 + 2 * rng.draw();

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

//! makes the antenna pattern a pencil beam
void MCAntenna::set_pencil_beam (void){
  atype=ANTENNA_TYPE_PENCIL_BEAM;
}

//! makes the antenna pattern a 2D gaussian specified by za and aa standard deviations
/*!
Givees the MCAntenna object a 2D gaussian response function
\param za_sigma The std. dev. parameter for zenith angle
\param aa_sigma The std. dev. parameter for azimuthal angle.
\author Cory Davis
\date 2005-12-02
 */
void MCAntenna::set_gaussian (const Numeric& za_sigma,
                              const Numeric& aa_sigma)
{
  atype=ANTENNA_TYPE_GAUSSIAN;
  sigma_za=za_sigma;
  sigma_aa=aa_sigma;
}

//! makes the antenna pattern a 2D gaussian specified by za and aa FWHM
/*!
Givees the MCAntenna object a 2D gaussian response function
\param za_fwhm The full width half maximum zenith angle
\param aa_fwhm The full width half maximum azimuthal angle.
\author Cory Davis
\date 2005-12-02
 */
void MCAntenna::set_gaussian_fwhm (const Numeric& za_fwhm,
                                   const Numeric& aa_fwhm)
{
  atype=ANTENNA_TYPE_GAUSSIAN;
  sigma_za=za_fwhm/2.3548;
  sigma_aa=aa_fwhm/2.3548;
}

//! makes the antenna pattern use a 2D lookup table to define the antenna response
/*!
The lookup antenna type is not yet implemented

\param za_grid_ zenith angle grid for the antenna response lookup table
\param aa_grid_ azimuthal angle grid for the antenna response lookup table
\param G_lookup_ the lookup table data
\author Cory Davis
\date 2005-12-02
 */
void MCAntenna::set_lookup (ConstVectorView& za_grid_,
                            ConstVectorView& aa_grid_,
                            ConstMatrixView& G_lookup_)
{
  atype=ANTENNA_TYPE_LOOKUP;
  za_grid=za_grid_;
  aa_grid=aa_grid_;
  G_lookup=G_lookup_;
}

//! returns the antenna type
/*!

\author Cory Davis
\date 2006-6-16
 */
AntennaType MCAntenna::get_type(void) const
{
  return atype;
}




//! draws a line of sight by sampling the antenna response function
/*!

\param sampled_rte_los Output: The sampled line of sight
\param rng a random number generator
\param bore_sight_los the bore sight LOS
\author Cory Davis
\date 2005-12-02
 */
void MCAntenna::draw_los(VectorView& sampled_rte_los,
                         Rng& rng,
                         ConstVectorView bore_sight_los) const
{

  switch ( atype )
    {
    case ANTENNA_TYPE_PENCIL_BEAM:
      sampled_rte_los=bore_sight_los;
      break;
    case ANTENNA_TYPE_GAUSSIAN:
      sampled_rte_los[0]=bore_sight_los[0]+ran_gaussian(rng,sigma_za);
      sampled_rte_los[1]=bore_sight_los[1]+ran_gaussian(rng,sigma_aa);
      if ( sampled_rte_los[1]>180 )
        {
          sampled_rte_los[1]-=360;
        }
      break;
      //    case ANTENNA_TYPE_LOOKUP:
      //ostringstream os;
      //os << "Antenna type ANTENNA_TYPE_LOOKUP not yet implemented.";
      //throw runtime_error( os.str() );
      //break;
    default:
      ostringstream os;
      os << "invalid Antenna type.";
      throw runtime_error( os.str() );
  }
  
}

ostream& operator<< (ostream& os, const MCAntenna&)
{
  os << "MCAntenna: Output operator not implemented";
  return os;
}
