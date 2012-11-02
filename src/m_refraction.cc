/* Copyright (C) 2003-2012 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
  ===  File description 
  ===========================================================================*/

/*!
  \file   m_refraction.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2003-01-09

  \brief  Workspace methods releated to refraction.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "absorption.h"
#include "arts.h"
#include "check_input.h"
#include "math_funcs.h"
#include "matpackI.h"
#include "messages.h"
#include "refraction.h"
#include "special_interp.h"
#include "abs_species_tags.h"

extern const Numeric ELECTRON_CHARGE;
extern const Numeric ELECTRON_MASS;
extern const Numeric PI;
extern const Numeric VACUUM_PERMITTIVITY;
extern const Numeric TORR2PA;

/*===========================================================================
  === WSMs
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void refr_indexFreeElectrons(
          Numeric&   refr_index,
          Numeric&   refr_index_group,
    const Vector&    f_grid,
    const Numeric&   rte_edensity,
    const Verbosity&)
{
  // The expression used in found in many textbooks, e.g. Rybicki and Lightman
  // (1979). Note that the refractive index corresponds to the phase velocity.

  static const Numeric k = ELECTRON_CHARGE * ELECTRON_CHARGE / 
                         ( VACUUM_PERMITTIVITY * ELECTRON_MASS * 4 * PI * PI );

  if( rte_edensity > 0 )
    {
      const Numeric f = ( f_grid[0] + last(f_grid) ) / 2.0;
      const Numeric a = rte_edensity*k/(f*f);
      const Numeric n = sqrt( 1 - a );

      if( a > 0.25 ) 
        {
          ostringstream os;
          os << "The frequency must at least be twice the plasma frequency.\n"
             << "For this particular point, the plasma frequency is: " 
             << sqrt(rte_edensity*k)/1e6 << " MHz.";
          throw runtime_error( os.str() );
        }

      refr_index       += n - 1;
      refr_index_group += 1/n - 1;
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refr_indexIR(
          Numeric&   refr_index,
          Numeric&   refr_index_group,
    const Numeric&   rte_pressure,
    const Numeric&   rte_temperature,
    const Verbosity&)
{
  static const Numeric bn0  = 1.000272620045304;
  static const Numeric bn02 = bn0*bn0;
  static const Numeric bk   = 288.16 * (bn02-1.0) / (1013.25*(bn02+2.0));

  // Pa -> hPa
  const Numeric n = sqrt( (2.0*bk*rte_pressure/100.0+rte_temperature) / 
                          ( rte_temperature-bk*rte_pressure/100.0) ) - 1; 

  refr_index       += n;
  refr_index_group += n;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refr_indexThayer(
          Numeric&   refr_index,
          Numeric&   refr_index_group,
    const Numeric&   rte_pressure,
    const Numeric&   rte_temperature,
    const Vector&    rte_vmr_list,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Verbosity& )
{
  if( abs_species.nelem() != rte_vmr_list.nelem() )
    throw runtime_error( "The number of tag groups differ between "
                                         "*rte_vmr_list* and *abs_species*." );

  Index   firstH2O = find_first_species_tg( abs_species,
                                      species_index_from_species_name("H2O") );

  if( firstH2O < 0 )
    throw runtime_error(
       "Water vapour is a required (must be a tag group in *abs_species*)." );

  const Numeric   e = rte_pressure * rte_vmr_list[firstH2O];

  const Numeric n = ( 77.6e-8 * ( rte_pressure - e ) + 
             ( 64.8e-8 + 3.776e-3 / rte_temperature ) * e ) / rte_temperature;

  refr_index       += n;
  refr_index_group += n;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refr_indexMWgeneral(
          Numeric&   refr_index,
          Numeric&   refr_index_group,
    const Numeric&   rte_pressure,
    const Numeric&   rte_temperature,
    const Vector&    rte_vmr_list,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Verbosity& )
{
//FIXME: Shall n be rescaled for sum(VMW)=1? Doing so now, but is it correct?
//       Test sensitivity to applying or not applying rescaling.

/*
   for now, hard-coding the reference refindices and refT/p. could make that
   some re-setabel parameters (like iso ratios... also regarding storing them in
   file/data struct per species)
*/
  const Numeric p0 = 760.*TORR2PA;
  const Numeric T0 = 273.15;

  // Number of refractive species:
  const Index nrs = 5;

  // This is hardwired here and quite primitive, but should do the job.
  // Set refractive index species names.
  ArrayOfString ref_spec_names(nrs);
  ref_spec_names[0]  = "N2";
  ref_spec_names[1]  = "O2";
  ref_spec_names[2] = "CO2";
  ref_spec_names[3]  = "H2";
  ref_spec_names[4]  = "He";
  // Set reference refractive indices
  // Values from Newell and Baird, 1965
  Vector ref_n(nrs);
  ref_n[0] = 293.81e-6;
  ref_n[1] = 266.95e-6;
  ref_n[2] = 495.16e-6;
  ref_n[3] = 135.77e-6;
  ref_n[4] =  34.51e-6;

// Checks
  if( abs_species.nelem() != rte_vmr_list.nelem() )
    throw runtime_error( "The number of tag groups differ between "
                                         "*rte_vmr_list* and *abs_species*." );

/*
   further checks:
   ? non-neg T
   ? VMRs (or refacting VMRs) have to at least add up to a threshold (see
   xsec_species ff)
   ?
*/


// Data management
  // Find the location of all refractive species in abs_species. Set to -1 if
  // not found. The length of array ref_spec_locations is the number of
  // considered refractive species (in this method: N2, O2, CO2, H2, He).
  // The value means:
  // -1 = not in abs_species
  // N  = species is number N in abs_species

  //Can't use this one as it inside gets the broadening species names and
  //number. Also, we would have to make a workaround for this_species. So,
  //instead we use a modified version of this function directly included here.
  /*find_broad_spec_locations(ref_spec_locations,
                             abs_species,
                             this_species);*/

  ArrayOfIndex ref_spec_locations(nrs);
  
  // Loop over all broadening species and see if we can find them in abs_species.
  for (Index i=0; i<nrs; ++i) {
    // Find associated internal species index (we do the lookup by index, not by name). 
    const Index isi = species_index_from_species_name(ref_spec_names[i]); 
    
    // Find position of broadening species isi in abs_species. The called
    // function returns -1 if not found, which is already the correct
    // treatment for this case that we also want here.
    ref_spec_locations[i] = find_first_species_tg(abs_species,isi);
  }


// The actual calculation
  // N_tot = sum (Nref_i *     p_i/p_0 * T0/T)
  //       = sum (Nref_i * vmr_i*p/p_0 * T0/T)
  //       = p/p_0 * T0/T *  sum (  Nref_i  * vmr_i)
  //       = p/p_0 * T0/T *  sum ( (1+Np_i) * vmr_i)
  //       = p/p_0 * T0/T * (sum (vmr_i) + sum(Np_i * vmr_i))
  // VMR rescale: /=sum(vmr_i)
  //       = p/p_0 * T0/T * (1 + sum(Np_i * vmr_i)/sum(vmr_i) )

  const Numeric ratioT = T0/rte_temperature;
  const Numeric ratiop = rte_pressure/p0;

  Numeric ref_spec_vmr_sum = 0.;
  Numeric n = 0.;

  // Add up refractive species, where available:
  for (Index i=0; i<nrs; ++i) {
      if ( ref_spec_locations[i] >= 0 ) {

          // Add to VMR sum:
          ref_spec_vmr_sum += rte_vmr_list[ref_spec_locations[i]];

          // refraction contribution (excluding the constant factor p/p_0 * T0/T):
          n += ref_n[i] * rte_vmr_list[ref_spec_locations[i]];
      }
  }

  if ( abs(ref_spec_vmr_sum-1) > 0.1 )
      {
        ostringstream os;
        os << "Error: The total VMR of all your defined refractive\n"
             << "species is " << ref_spec_vmr_sum
             << ", more than 10% " << "different from 1.\n";
        throw runtime_error(os.str());
      }
    
  // normalize refractive index with the considered total VMR and add offset
  // (n=1) part:
  n /= ref_spec_vmr_sum;
  n += 1.;
  // as above, but with out normalization:
  //n += ref_spec_vmr_sum;

  // now applying the constant factor p/p_0 * T0/T and
  n *= (ratioT*ratiop);

  refr_index       += n;
  refr_index_group += n;
}

