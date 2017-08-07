/* Copyright (C) 2011-2017 Jana Mendrok <jana.mendrok@gmail.com>
                      
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
   USA. 
*/

/*!
  \file   microphysics.cc
  \author Jana Mendrok, Daniel Kreyling, Manfred Brath, Patrick Eriksson
  \date   2017-08-01
  
  \brief  Internal functions for microphysics calculations (size distributions etc.)
*/

#include "microphysics.h"

extern const Numeric PI;
extern const Numeric DENSITY_OF_ICE;
extern const Numeric DENSITY_OF_WATER;

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>
#include <stdexcept>

#include "arts.h"
#include "check_input.h"
#include "cloudbox.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "mc_antenna.h"
#include "messages.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rng.h"
#include "sorting.h"



/*! Calculates the particle number density field for McFarquhar and Heymsfield
    (1997) size distribution. To be used for cloud ice.

    Negative IWC trigger an error (unless robust=1, where IWC=0 is applied,
    hence pnd_field=0 is returned. Negative temperatures trigger an error,
    temperatures >280K are only accepted if robust=1. For temperatures >273.15K
    (=0C), the distribution is evaluated assuming T=273.15K.

    \return pnd_field   Particle number density field
    \param IWC_field    mass content (cloud ice) field [kg/m3]
    \param t_field      atmospheric temperature [K]
    \param limits       pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta    scattering meta data for all scattering elements
    \param scat_species array index of scattering species handled by this distribution
    \param part_string  scat_species tag for profile/distribution handled here
    \param delim        delimiter string of *scat_species* elements
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2012-04-03

*/
void pnd_fieldMH97 (Tensor4View pnd_field,
                    const Tensor3& IWC_field,
                    const Tensor3& t_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfArrayOfScatteringMetaData& scat_meta,
                    const Index& scat_species,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity)
{
  const String psdname="MH97";
  const Index N_se = scat_meta[scat_species].nelem();
  const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
  ArrayOfIndex intarr;
  Vector mass_unsorted ( N_se, 0.0 );
  Vector mass ( N_se, 0.0 );
  Vector diameter_mass_equivalent ( N_se, 0.0 );
  Vector pnd ( N_se, 0.0 );
  Vector dNdD ( N_se, 0.0 );

  String psd_param;
  String partfield_name;
  ArrayOfString psd_options;

  //split String and copy to ArrayOfString
  parse_psd_param( psd_param, part_string, delim);
  parse_partfield_name( partfield_name, part_string, delim);

  bool noisy = (psd_param == "MH97n");
  bool robust = false;
  parse_psd_options( psd_options, part_string, delim);
  for ( Index i=0; i<psd_options.nelem(); i++ )
  {
    robust = (robust || (psd_options[i]=="robust") );
  }

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].mass) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element mass.\n"
             << "But mass is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      mass_unsorted[i] = ( scat_meta[scat_species][i].mass );
    }
  get_sorted_indexes(intarr, mass_unsorted);

  // extract scattering meta data
  for ( Index i=0; i< N_se; i++ )
  {
    mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
    if ( isnan(scat_meta[scat_species][i].diameter_volume_equ) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element volume equ. diameter.\n"
             << "But volume equ. diameter is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
    diameter_mass_equivalent[i] =
      scat_meta[scat_species][intarr[i]].diameter_volume_equ;
  }
  
  if (mass.nelem() > 0)
  // mass.nelem()=0 implies no selected scattering element for the respective
  // scattering species field. should not occur.
  {
      assert( is_increasing(diameter_mass_equivalent) );

      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // a valid IWC value encountered. Calculating dNdD.
            if (IWC_field ( p, lat, lon ) > 0.)
              {
                Numeric T = t_field ( p, lat, lon );

                // abort if T is negative
                if ( T<0. )
                  {
                    ostringstream os;
                    os << "Negative temperatures not allowed.\n"
                       << "Yours is " << T << "K.";
                    throw runtime_error ( os.str() );
                  }
                // abort if T is too high
                if ( !robust && T>280. )
                  {
                    ostringstream os;
                    os << "Temperatures above 280K not allowed by MH97"
                       << " (to allow: run with robust option).\n"
                       << "Yours is " << T << "K.";
                    throw runtime_error ( os.str() );
                  }
                // allow some margin on T>0C (but use T=0C for PSD calc)
                T = min( T, 273.15 );

                psd_cloudice_MH97 ( dNdD, diameter_mass_equivalent,
                              IWC_field(p,lat,lon), T, noisy );

                // ensure that any particles where produced
                if ( dNdD.sum() == 0.0 )
                  { // no particles at all were produced. means, there's some
                    // issue with the setup (none of the included particles
                    // produces numerically considerable amount of particles
                    // with this PSD, likely due to all the particles being
                    // either too large or too small)
                    ostringstream os;
                    os <<  psdname << " PSD did not produce any non-zero pnd values "
                       << "(likely results\n"
                       << "from considered particles being too large or too "
                       << "small)\n"
                       << "The problem occured for profile '"<< partfield_name
                       << "' at: " << "p = " << p << ", lat = " << lat
                       << ", lon = " << lon << ".\n";
                    throw runtime_error ( os.str() );
                  }

                // derive pnds from scaling dNdD with bin width
                if (diameter_mass_equivalent.nelem() > 1)
                  bin_integral( pnd, diameter_mass_equivalent, dNdD );
                else
                  pnd = dNdD;
	    
                // calculate error of pnd sum and real XWC
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), mass,
                             p, lat, lon, partfield_name, verbosity );
	    
                // writing pnd vector to WSV pnd_field
                for ( Index i = 0; i< N_se; i++ )
                  {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                  }
              }

            // MH97 requires mass density. If not set, abort calculation.
            else if ( isnan(IWC_field ( p, lat, lon )) )
              {
                ostringstream os;
                os << "Size distribution " << psdname
                   << " requires knowledge of mass "
                   << "density of atmospheric ice.\n"
                   << "At grid point (" << p << ", " << lat << ", " << lon
                   << ") in (p,lat,lon) a NaN value is encountered, "
                   << "i.e. mass density is unknown.";
                throw runtime_error( os.str() );
              }

            // MH97 PSD is parameterized in IWC and can not handle negative
            // numbers, hence abort.
            else if (!robust && IWC_field ( p, lat, lon ) < 0.)
              {
                ostringstream os;
                os << "Size distribution " << psdname
                   << " is parametrized in ice mass content.\n"
                   << "It can not handle negative values like IWC="
                   << IWC_field (p,lat,lon) << " kg/m3\n"
                   << "found at grid point ("
                   << p << ", " << lat << ", " << lon << ")";
                throw runtime_error( os.str() );
              }

            // for IWC==0 or robust&&IWC<0. we just set pnd_field=0
            else
                for ( Index i = 0; i< N_se; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = 0.;
                }
          }
        }
      }
  }
}


/*! Calculates the particle number density field for Heymsfield (2011, personal
    comm.) size distribution. To be used for atmospheric ice, particularly cloud
    ice and snow.

    \return pnd_field   Particle number density field
    \param IWC_field    mass content (cloud ice or snow) field [kg/m3]
    \param t_field      atmospheric temperature [K]
    \param limits       pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta    scattering meta data for all scattering elements
    \param scat_species array index of scattering species handled by this distribution
    \param part_string  scat_species tag for profile/distribution handled here
    \param delim        delimiter string of *scat_species* elements
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2012-04-05

*/
void pnd_fieldH11 (Tensor4View pnd_field,
                   const Tensor3& IWC_field,
                   const Tensor3& t_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
  const String psdname="H11";
  const Index N_se = scat_meta[scat_species].nelem();
  const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
  ArrayOfIndex intarr;
  Vector diameter_max_unsorted ( N_se, 0.0 );
  Vector diameter_max ( N_se, 0.0 );
  Vector mass ( N_se, 0.0 );
  Vector pnd ( N_se, 0.0 );
  Vector dNdD ( N_se, 0.0 );
  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].diameter_max) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element maximum diameter.\n"
             << "But maximum diameter is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      diameter_max_unsorted[i] = ( scat_meta[scat_species][i].diameter_max );
    }
  get_sorted_indexes(intarr, diameter_max_unsorted);
      
  // extract scattering meta data
  for ( Index i=0; i< N_se; i++ )
  {
      diameter_max[i] = scat_meta[scat_species][intarr[i]].diameter_max; // [m]

      if ( isnan(scat_meta[scat_species][intarr[i]].mass) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element mass.\n"
             << "But mass is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
  }

  if (diameter_max.nelem() > 0)
  // diameter_max.nelem()=0 implies no selected scattering element for the respective
  // scattering species field. should not occur anymore.
  {
      assert( is_increasing(diameter_max) );

      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // H11 requires mass density. If not set, abort calculation.
            if ( isnan(IWC_field ( p, lat, lon )) )
              {
                ostringstream os;
                os << "Size distribution " << psdname << " requires knowledge of mass "
                   << "density of atmospheric ice.\n"
                   << "At grid point (" << p << ", " << lat << ", " << lon
                   << ") in (p,lat,lon) a NaN value is encountered, "
                   << "i.e. mass density is unknown.";
                throw runtime_error( os.str() );
              }

            // A valid IWC value encountered (here, we can also handle negative
            // values!). Calculating dNdD.
            else if (IWC_field ( p, lat, lon ) != 0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<diameter_max.nelem(); i++ ) //loop over number of scattering elements
                {
                    // calculate particle size distribution for H11
                    // [# m^-3 m^-1]
                    dNdD[i] = IWCtopnd_H11 ( diameter_max[i], t_field ( p, lat, lon ) );
                }

                // ensure that any particles where produced
                if ( dNdD.sum() == 0.0 )
                  { // no particles at all were produced. means, there's some
                    // issue with the setup (none of the included particles
                    // produces numerically considerable amount of particles
                    // with this PSD, likely due to all the particles being
                    // either too large or too small)
                    ostringstream os;
                    os <<  psdname << " PSD did not produce any non-zero pnd values "
                       << "(likely results\n"
                       << "from considered particles being too large or too "
                       << "small)\n"
                       << "The problem occured for profile '"<< partfield_name
                       << "' at: " << "p = " << p << ", lat = " << lat
                       << ", lon = " << lon << ".\n";
                    throw runtime_error ( os.str() );
                  }

                // scale pnds by scale width
                if (diameter_max.nelem() > 1)
                    bin_integral( pnd, diameter_max, dNdD ); //[# m^-3]
                else
                    pnd = dNdD;

                // calculate proper scaling of pnd sum from real IWC and apply
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), mass,
                             p, lat, lon, partfield_name, verbosity );

                // writing pnd vector to wsv pnd_field
                for ( Index i =0; i< N_se; i++ )
                {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                }
            }

            // for IWC==0, we just set pnd_field=0
            else
            {
                for ( Index i = 0; i< N_se; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = 0.;
                }
            }
          }
        }
      }
  }
}


/*! Calculates the particle number density field for Heymsfield (2013, personal
    comm.) size distribution. To be used for atmospheric ice, particularly cloud
    ice and snow.

    \return pnd_field   Particle number density field
    \param IWC_field    mass content (cloud ice or snow) field [kg/m3]
    \param t_field      atmospheric temperature [K]
    \param limits       pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta    scattering meta data for all scattering elements
    \param scat_species array index of scattering species handled by this distribution
    \param part_string  scat_species tag for profile/distribution handled here
    \param delim        delimiter string of *scat_species* elements
  
  \author Johan Strandgren, Daniel Kreyling
  \date 2013-08-26

*/
void pnd_fieldH13 (Tensor4View pnd_field,
                   const Tensor3& IWC_field,
                   const Tensor3& t_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
  const String psdname="H13";
  const Index N_se = scat_meta[scat_species].nelem();
  const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
  ArrayOfIndex intarr;
  Vector diameter_max_unsorted ( N_se, 0.0 );
  Vector diameter_max ( N_se, 0.0 );
  Vector mass ( N_se, 0.0 );
  Vector pnd ( N_se, 0.0 );
  Vector dNdD ( N_se, 0.0 );
  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].diameter_max) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element maximum diameter.\n"
             << "But maximum diameter is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      diameter_max_unsorted[i] = ( scat_meta[scat_species][i].diameter_max );
    }
  get_sorted_indexes(intarr, diameter_max_unsorted);
      
  // extract scattering meta data
  for ( Index i=0; i< N_se; i++ )
  {
      diameter_max[i] = scat_meta[scat_species][intarr[i]].diameter_max; // [m]

      if ( isnan(scat_meta[scat_species][intarr[i]].mass) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element mass.\n"
             << "But mass is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
  }

  if (diameter_max.nelem() > 0)
  // diameter_max.nelem()=0 implies no selected scattering elements for the respective
  // scattering species field. should not occur anymore.
  {
      assert( is_increasing(diameter_max) );

      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // H13 requires mass density. If not set, abort calculation.
            if ( isnan(IWC_field ( p, lat, lon )) )
              {
                ostringstream os;
                os << "Size distribution " << psdname << " requires knowledge of mass "
                   << "density of atmospheric ice.\n"
                   << "At grid point (" << p << ", " << lat << ", " << lon
                   << ") in (p,lat,lon) a NaN value is encountered, "
                   << "i.e. mass density is unknown.";
                throw runtime_error( os.str() );
              }

            // A valid IWC value encountered (here, we can also handle negative
            // values!). Calculating dNdD.
            else if (IWC_field ( p, lat, lon ) != 0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<diameter_max.nelem(); i++ ) //loop over number of scattering elements
                {
                    // calculate particle size distribution for H13
                    // [# m^-3 m^-1]
                    dNdD[i] = IWCtopnd_H13 ( diameter_max[i], t_field ( p, lat, lon ) );
                }

                // ensure that any particles where produced
                if ( dNdD.sum() == 0.0 )
                  { // no particles at all were produced. means, there's some
                    // issue with the setup (none of the included particles
                    // produces numerically considerable amount of particles
                    // with this PSD, likely due to all the particles being
                    // either too large or too small)
                    ostringstream os;
                    os <<  psdname << " PSD did not produce any non-zero pnd values "
                       << "(likely results\n"
                       << "from considered particles being too large or too "
                       << "small)\n"
                       << "The problem occured for profile '"<< partfield_name
                       << "' at: " << "p = " << p << ", lat = " << lat
                       << ", lon = " << lon << ".\n";
                    throw runtime_error ( os.str() );
                  }

                // scale pnds by scale width
                if (diameter_max.nelem() > 1)
                    bin_integral( pnd, diameter_max, dNdD ); //[# m^-3]
                else
                    pnd = dNdD;

                // calculate proper scaling of pnd sum from real IWC and apply
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), mass,
                             p, lat, lon, partfield_name, verbosity );

                // writing pnd vector to wsv pnd_field
                for ( Index i =0; i< N_se; i++ )
                {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                }
            }

            // for IWC==0, we just set pnd_field=0
            else
            {
                for ( Index i = 0; i< N_se; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = 0.;
                }
            }
          }
        }
      }
  }
}

/*! Calculates the particle number density field for Heymsfield (2013, personal
    comm.) size distribution. To be used for atmospheric ice, particularly cloud
    ice and snow.

    \return pnd_field   Particle number density field
    \param IWC_field    mass content (cloud ice or snow) field [kg/m3]
    \param t_field      atmospheric temperature [K]
    \param limits       pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta    scattering meta data for all scattering elements
    \param scat_species array index of scattering species handled by this distribution
    \param part_string  scat_species tag for profile/distribution handled here
    \param delim        delimiter string of *scat_species* elements
  
  \author Johan Strandgren
  \date 2013-08-26

*/
void pnd_fieldH13Shape (Tensor4View pnd_field,
                        const Tensor3& IWC_field,
                        const Tensor3& t_field,
                        const ArrayOfIndex& limits,
                        const ArrayOfArrayOfScatteringMetaData& scat_meta,
                        const Index& scat_species,
                        const String& part_string,
                        const String& delim,
                        const Verbosity& verbosity)
{ 
  CREATE_OUT1;  
    
  const String psdname="H13shape";
  const Index N_se = scat_meta[scat_species].nelem();
  const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
  ArrayOfIndex intarr;
  Vector diameter_max_unsorted ( N_se, 0.0 );
  Vector diameter_max ( N_se, 0.0 );
  Vector mass ( N_se, 0.0 );
  Vector pnd ( N_se, 0.0 );
  Vector diameter_area_equivalent ( N_se, 0.0 );
  Vector Rarea ( N_se, 0.0 ); // Area ratio = diameter_area_equivalent^2/diameter_max^2
  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].diameter_max) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element maximum diameter.\n"
             << "But maximum diameter is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      diameter_max_unsorted[i] = ( scat_meta[scat_species][i].diameter_max );
    }
  get_sorted_indexes(intarr, diameter_max_unsorted);
  
  // extract scattering meta data
  for ( Index i=0; i< N_se; i++ )
  {
      diameter_max[i] = scat_meta[scat_species][intarr[i]].diameter_max; // [m]

      if ( isnan(scat_meta[scat_species][intarr[i]].diameter_area_equ_aerodynamical) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element area equivalent diameter.\n"
             << "But area equivalent diameter is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      diameter_area_equivalent[i] = scat_meta[scat_species][intarr[i]].diameter_area_equ_aerodynamical; // [m]

      if ( isnan(scat_meta[scat_species][intarr[i]].mass) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element mass.\n"
             << "But mass is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]

      Rarea[i] = (diameter_area_equivalent[i]*diameter_area_equivalent[i]) / 
                 (diameter_max[i]*diameter_max[i]); // [m2/m2]
  }
    // Collect all unique maximum diameters
    vector<Numeric> diameter_max_in;
    for (Iterator1D it = diameter_max.begin(); it != diameter_max.end(); ++it)
        if (find(diameter_max_in.begin(), diameter_max_in.end(), *it) == diameter_max_in.end())
            diameter_max_in.push_back(*it);
    
    // Collect all unique area ratios
    vector<Numeric> Rarea_in;
    for (Iterator1D it = Rarea.begin(); it != Rarea.end(); ++it)
        if (find(Rarea_in.begin(), Rarea_in.end(), *it) == Rarea_in.end())
            Rarea_in.push_back(*it);
                
    Vector diameter_max_input;
    Vector Rarea_input;
    diameter_max_input=diameter_max_in;
    Rarea_input=Rarea_in;
    
    // Check size and shape limits
    if (diameter_max[0]<7.7*1e-5)
    {
        ostringstream os;
        os << "The " << psdname << " parametrization is only valid for particles with\n"
           << "a maximum diameter >= to 77 um, but the smallest value of\n"
           << "*diameter_max* in this simulation is " << diameter_max[0] << " um.\n"
           << "Set a new *diameter_max_grid*!\n";
        throw runtime_error(os.str());
    }

    if (Rarea_input.nelem()==1)
    { 
        out1 << "WARNING! Only one unique area ratio is used.\n"
             << "Using parametrization H13 will generate the same results "
             << "as " << psdname << "\n"
             << "but with less computational effort, hence in a shorter time.\n";
    }
    
    Vector dNdD ( diameter_max_input.nelem(), 0.0 );
    Vector Ar ( diameter_max_input.nelem(), 0.0 );
    Vector pnd_temp ( diameter_max_input.nelem(), 0.0 );

    
  //const bool suppress=true;
  //const Verbosity temp_verb(0,0,0);

  if (diameter_max_input.nelem() > 0)
  // diameter_max.nelem()=0 implies no selected scattering elements for the respective
  // scattering species field. should not occur anymore.
  {
      assert( is_increasing(diameter_max_input) );

      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // H13 requires mass density. If not set, abort calculation.
            if ( isnan(IWC_field ( p, lat, lon )) )
              {
                ostringstream os;
                os << "Size distribution " << psdname << " requires knowledge of mass "
                   << "density of atmospheric ice.\n"
                   << "At grid point (" << p << ", " << lat << ", " << lon
                   << ") in (p,lat,lon) a NaN value is encountered, "
                   << "i.e. mass density is unknown.";
                throw runtime_error( os.str() );
              }

            // A valid IWC value encountered (here, we can also handle negative
            // values!). Calculating dNdD.
            else if (IWC_field ( p, lat, lon ) != 0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<diameter_max_input.nelem(); i++ ) //loop over number of scattering elements
                {
                    // calculate particle size distribution for H13Shape
                    // [# m^-3 m^-1]
                    dNdD[i] = IWCtopnd_H13Shape ( diameter_max_input[i], t_field ( p, lat, lon ) );
                    
                    // calculate Area ratio distribution for H13Shape
                    Ar[i] = area_ratioH13 (diameter_max_input[i], t_field (p, lat, lon ) );
                }

                // ensure that any particles where produced
                if ( dNdD.sum() == 0.0 )
                  { // no particles at all were produced. means, there's some
                    // issue with the setup (none of the included particles
                    // produces numerically considerable amount of particles
                    // with this PSD, likely due to all the particles being
                    // either too large or too small)
                    ostringstream os;
                    os <<  psdname << " PSD did not produce any non-zero pnd values "
                       << "(likely results\n"
                       << "from considered particles being too large or too "
                       << "small)\n"
                       << "The problem occured for profile '"<< partfield_name
                       << "' at: " << "p = " << p << ", lat = " << lat
                       << ", lon = " << lon << ".\n";
                    throw runtime_error ( os.str() );
                  }

                // scale pnds by scale width
                if (diameter_max_input.nelem() > 1)
                    bin_integral( pnd_temp, diameter_max_input, dNdD ); //[# m^-3]
                else
                    pnd_temp = dNdD;

                // Check which element in Rarea is closest to Ar and assign
                // the PND for that size to that scattering element and zeros to the rest
                Index l;
                l=Rarea_input.nelem();

                Vector diff;
                
                for ( Index k=0, j=0; j<pnd_temp.nelem(); k+=l,j++ )
                {   
                    diff = Rarea[Range(k,l)];
                    
                    diff -= Ar[j];
                    
                    Numeric minval = std::numeric_limits<Numeric>::infinity();
                    Index   minpos = -1;
                    
                    for (Index i = 0; i < diff.nelem(); i++)
                    {
                        if (abs(diff[i]) < minval)
                        {
                            minval = abs(diff[i]);
                            minpos = i;
                        }
                    }
                        pnd[Range(k,l)]=0;
                        pnd[minpos+k]=pnd_temp[j];
                     
                    
                }

                // calculate proper scaling of pnd sum from real IWC and apply
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), mass,
                             p, lat, lon, partfield_name, verbosity );
               
                // writing pnd vector to wsv pnd_field
                for ( Index i =0; i< N_se; i++ )
                {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                }
            }

            // for IWC==0, we just set pnd_field=0
            else
            {
                for ( Index i = 0; i< N_se; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = 0.;
                }
            }
          }
        }
      }
  }
}


/*! Calculates the particle number density field for Field (2007) size
 *  distribution, both for tropics or midlatitudes.
 *
 *  For the estimation of the second moment the mass dimension
 *  relationship is estimated by regression from the meta data.
 *
 *  To be used for snow, to be given in terms of mass content.

    Negative SWC trigger an error (unless robust=1, where SWC=0 is applied,
    hence pnd_field=0 is returned. Negative temperatures trigger an error.
    No further temperature validity limits are implemented (although the
    parametrization itself has been derived from measurements limited to
    -60C<=T<=0C, i.e. shall strictly only be applied then).
 
 \param pnd_field Particle number density field
 \param SWC_field (snow) mass content field [kg/m3]
 \param t_field atmospheric temperature [K]
 \\param limits pnd_field boundaries (indices in p, lat, lon)
 \param scat_meta_array particle meta data for particles
 \param scat_data_start start index for particles handled by this distribution
 \param npart number of particles handled by this distribution
 \param part_string part_species tag for profile/distribution handled here
 \param delim Delimiter string of *part_species* elements
 
 \author Manfred Brath, Jana Mendrok, D. Kreyling
 \date 2014-12-02
 
 */
void pnd_fieldF07 (Tensor4View pnd_field,
                   const Tensor3& SWC_field,
                   const Tensor3& t_field,
                   const String& regime,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
  assert( (regime=="TR") || (regime=="ML") );
  const String psdname=string("F07")+regime;

  String partfield_name;
  ArrayOfString psd_options;

  const Index N_se = scat_meta[scat_species].nelem();
  const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
  ArrayOfIndex intarr;

  Vector diameter_max_unsorted ( N_se, 0.0 );
  Vector diameter_max ( N_se, 0.0 );
  Vector mass ( N_se, 0.0 );
  Vector pnd ( N_se, 0.0 );
  Vector dNdD ( N_se, 0.0 );

  Numeric alpha;
  Numeric beta;
  Vector log_m( N_se, 0.0 );
  Vector log_D( N_se, 0.0 );
  Vector q;
    
  // diameter_max.nelem()=0 implies no selected scattering element for the respective
  // scattering species field. should not occur anymore.
  if ( diameter_max.nelem() == 0 )
    return;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);
    
  bool robust = false;
  parse_psd_options( psd_options, part_string, delim);
  for ( Index i=0; i<psd_options.nelem(); i++ )
    {
      robust = (robust || (psd_options[i]=="robust") );
    }

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].diameter_max) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element maximum diameter.\n"
            << "But maximum diameter is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
        diameter_max_unsorted[i] = ( scat_meta[scat_species][i].diameter_max );
    }
  get_sorted_indexes(intarr, diameter_max_unsorted);
    
  // extract scattering meta data
  for ( Index i=0; i< N_se; i++ )
    {
      diameter_max[i] = scat_meta[scat_species][intarr[i]].diameter_max; // [m]
        
      if ( isnan(scat_meta[scat_species][intarr[i]].mass) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element mass.\n"
            << "But mass is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
      mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
        
      // logarithm of Dmax, needed for estimating mass-dimension-relationship
      log_D[i]=log(diameter_max[i]);
        
      // logarithm of the mass, even though it is  little weird to have
      // a logarithm of something with a unit...
      log_m[i]=log(mass[i]);
    }
    
  if ( N_se>1 )
    {
      //estimate mass-dimension relationship from meta data by linear regression
      // Assumption of a power law for the mass dimension relationship
      // Approach: log(m) = log(alpha)+beta*log(dmax/D0)
      linreg(q,log_D, log_m);

      alpha=exp(q[0]);
      beta=q[1];
    }
  else
    {
      // for a monodispersion we can't estimate the m-D relation (1 relation, 2
      // unknowns), hence we fix one of the parameters and calculate the other
      // such that we have them consistent. but shouldn't make any difference on
      // the end result whatever we choose here (all ice has to end up with this
      // scattering anyways)
      beta=2;
      alpha=mass[0]/pow(diameter_max[0],beta);
    }
    
  CREATE_OUT2;
  out2 << "Mass-dimension relationship m=alpha*(dmax)^beta:\n"
       << "alpha = " << alpha << " kg \n"
       << "beta = " << beta << "\n";

  assert( is_increasing(diameter_max) );
  // iteration over all atm. levels
  for ( Index p=limits[0]; p<limits[1]; p++ )
    {
      for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            if (SWC_field ( p, lat, lon ) > 0.)
            {
              Numeric T = t_field ( p, lat, lon );

              // abort if T is negative
              if ( T<0. )
                {
                  ostringstream os;
                  os << "Negative temperatures not allowed.\n"
                     << "Yours is " << T << "K.";
                  throw runtime_error ( os.str() );
                }

              psd_snow_F07( dNdD, diameter_max,
                            SWC_field ( p, lat, lon ), T,
                            alpha, beta, regime);

              // ensure that any particles where produced
              if ( dNdD.sum() == 0.0 )
              { // no particles at all were produced. means, there's some
                // issue with the setup (none of the included particles
                // produces numerically considerable amount of particles
                // with this PSD, likely due to all the particles being
                // either too large or too small)
                ostringstream os;
                os <<  psdname << " PSD did not produce any non-zero pnd values "
                  << "(likely results\n"
                  << "from considered particles being too large or too small)\n"
                  << "The problem occured for profile '" << partfield_name
                  << "' at: " << "p = " << p << ", lat = " << lat
                  << ", lon = " << lon << ".\n";
                throw runtime_error ( os.str() );
              }

              // scale pnds by scale width
              if (diameter_max.nelem() > 1)
                bin_integral( pnd, diameter_max, dNdD ); //[# m^-3]
              else
                pnd = dNdD;
                        
              // calculate proper scaling of pnd sum from real IWC and apply
              chk_pndsum ( pnd, SWC_field ( p,lat,lon ), mass,
                           p, lat, lon, partfield_name, verbosity );
                        
              // writing pnd vector to wsv pnd_field
              for ( Index i =0; i< N_se; i++ )
              {
                pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                            lat-limits[2], lon-limits[4] ) = pnd[i];
              }
            }

            // F07 requires mass density. If not set, abort calculation.
            else if ( isnan(SWC_field ( p, lat, lon )) )
            {
              ostringstream os;
              os << "Size distribution " << psdname << " requires knowledge"
                 << " of mass density of atmospheric ice.\n"
                 << "At grid point (" << p << ", " << lat << ", " << lon
                 << ") in (p,lat,lon) a NaN value is encountered, "
                 << "i.e. mass density is unknown.";
              throw runtime_error( os.str() );
            }
                    
            // F07 can not handle negative SWC.
            else if (!robust && SWC_field ( p, lat, lon ) < 0.)
            {
              ostringstream os;
              os << "Size distribution " << psdname
                 << " is parametrized in snow mass content.\n"
                 << "It does not handle negative values like SWC="
                 << SWC_field (p,lat,lon) << " kg/m3\n"
                 << "found at grid point ("
                 << p << ", " << lat << ", " << lon << ")";
              throw runtime_error( os.str() );
            }

            // for SWC==0 or robust&&SWC<0. we just set pnd_field=0
            else
            {
              for ( Index i = 0; i< N_se; i++ )
              {
                pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                            lat-limits[2], lon-limits[4] ) = 0.;
              }
            }
          }
        }
    }
}


/*! Calculates the particle number density field according
 *  to the two moment scheme of Axel Seifert, that is used the ICON model.
 
 \param pnd_field Particle number density field
 \param WC_field  mass content field [kg/m3]
 \param N_field Total Number density field [1/m^3]
 \\param limits pnd_field boundaries (indices in p, lat, lon)
 \param scat_meta_array particle meta data for particles
 \param scat_data_start start index for particles handled by this distribution
 \param npart number of particles handled by this distribution
 \param part_string part_species tag for profile/distribution handled here
 \param delim Delimiter string of *part_species* elements
 \param psd_type string with a tag defining the (hydrometeor) scheme
 
 \author Manfred Brath (parts of the function is based on pnd_H11 (of J. Mendrok & D. Kreyling))
 \date 2015-01-19
 
 */
void pnd_fieldS2M (Tensor4View pnd_field,
                     const Tensor3& WC_field,
                     const Tensor3& N_field,
                     const ArrayOfIndex& limits,
                     const ArrayOfArrayOfScatteringMetaData& scat_meta,
                     const Index& scat_species,
                     const String& part_string,
                     const String& delim,
                     const Verbosity& verbosity)
{
    const String psdname="S2M";
    const Index N_se = scat_meta[scat_species].nelem();
    const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
    ArrayOfIndex intarr;
    Vector mass_unsorted ( N_se, 0.0 );
    Vector mass ( N_se, 0.0 );
    Vector pnd ( N_se, 0.0 );
    Vector dNdD ( N_se, 0.0 );
    String partfield_name;
    String psd_str;
    String psd_param;
    bool logic_M;
    Numeric N_tot;
    ArrayOfString substrings;
    
    
    //split String and copy to ArrayOfString
    parse_partfield_name( partfield_name, part_string, delim);
    parse_psd_param( psd_param, part_string, delim);
    
    // Check if N_field should be handled as number density field
    // or as mean particle mass field
    psd_param.split(substrings,"_");
    
    if (substrings.size()==3)
    {
        // case: N_field is mean particle mass
        if (substrings[2] == "M")
        {
            psd_str=psd_param.substr(0,psd_param.length()-2);

            logic_M=true;
        }
        else
        {
            ostringstream os;
            os << "You use a wrong tag! The last substring of the tag\n"
            << "has to be an uppercase M";
            throw runtime_error( os.str() );
        }
    }
    // case: N_field is total number density
    else
    {
        logic_M=false;
        psd_str=psd_param;
    }

    for ( Index i=0; i < N_se; i++ )
    {
        if ( isnan(scat_meta[scat_species][i].mass) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element mass.\n"
            << "But mass is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
        mass_unsorted[i] = ( scat_meta[scat_species][i].mass );
    }
    get_sorted_indexes(intarr, mass_unsorted);
    
    // extract scattering meta data
    for ( Index i=0; i< N_se; i++ )
    {
        mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
    }

    if (mass.nelem() > 0)
        // diameter_max.nelem()=0 implies no selected scattering element for the respective
        // scattering species field. should not occur anymore.
    {
        assert( is_increasing(mass) );

        // iteration over all atm. levels
        for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // S2M requires mass density. If not set, abort calculation.
                    if ( isnan(WC_field ( p, lat, lon )) || isnan(N_field ( p, lat, lon )))
                    {
                        if (logic_M)
                        {
                            ostringstream os;
                            os << "Size distribution " << psdname << " requires knowledge of mass "
                            << "concentration of hydrometeors and mean particle mass.\n"
                            << "At grid point (" << p << ", " << lat << ", " << lon
                            << ") in (p,lat,lon) a NaN value is encountered, "
                            << "i.e. mass concentration amd/or mean particle mass is unknown.";
                            throw runtime_error( os.str() );
                        }
                        else
                        {
                            ostringstream os;
                            os << "Size distribution " << psdname << " requires knowledge of mass "
                            << "concentration of hydrometeors and number density.\n"
                            << "At grid point (" << p << ", " << lat << ", " << lon
                            << ") in (p,lat,lon) a NaN value is encountered, "
                            << "i.e. mass concentration amd/or mean particle mass is unknown.";
                            throw runtime_error( os.str() );
                        }
                    }
                    
                    // A valid IWC value encountered (here, we can also handle negative
                    // values!). Calculating dNdD.
                    else if (WC_field ( p, lat, lon ) != 0. ) // hier n0 abfrage weg
                    {
                        
                        if (logic_M)
                        {
                            if (N_field ( p, lat, lon ) != 0.)
                            {
                                //case: N_field is mean particle mass
                                N_tot=WC_field ( p, lat, lon )/N_field ( p, lat, lon );
                            }
                            else//if mean particle mass is zero, then the number density is
                                //set to zero. The number density will be adjusted to
                                //a non-zero within psd_S2M due to the limits
                                //of the scheme
                            {
                                N_tot=0;
                            }
                        }
                        else
                        {
                            //case: N_fied is total number density
                            N_tot=N_field ( p, lat, lon );
                        }
                            
                        
                        
                        // iteration over all given size bins
                        for ( Index i=0; i<mass.nelem(); i++ ) //loop over number of scattering elements
                        {
                            // calculate particle size distribution for H11
                            // [# m^-3 m^-1]
                            dNdD[i] = psd_S2M( mass[i], N_tot,
                                                     WC_field ( p, lat, lon ),
                                                     psd_str);
                        }
                        
                        // sometimes there is possibility if WC_field is very small
                        // but still greater than zero and N_tot is large that dNdD
                        // could be zero
                        if (dNdD.sum()==0)
                        {
                            
                            for ( Index i = 0; i< N_se; i++ )
                            {
                                pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                           lat-limits[2], lon-limits[4] ) = 0.;
                            }
                        }
                        else
                        {
                        
                            // scale pnds by scale width
                            if (mass.nelem() > 1)
                                bin_integral( pnd, mass, dNdD ); //[# m^-3]
                            else
                                pnd = dNdD;
                            
                           
                            
                            // calculate proper scaling of pnd sum from real IWC and apply
                            chk_pndsum ( pnd, WC_field ( p,lat,lon ), mass,
                                        p, lat, lon, partfield_name, verbosity );
                            
                            // writing pnd vector to wsv pnd_field
                            for ( Index i =0; i< N_se; i++ )
                            {
                                pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                           lat-limits[2], lon-limits[4] ) = pnd[i];
                            }
                        }
                    }
                    
                    // if either WC_field or N_field is zero, we just set pnd_field=0
                    else
                    {
                        for ( Index i = 0; i< N_se; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = 0.;
                        }
                    }
                }
            }
        }
    }
}


/*! Calculates the particle number density field according
 *  to the Milbrandt and Yau two moment scheme, which is used in the GEM model.
 *  See also Milbrandt and Yau, 2005.
 
 \param pnd_field Particle number density field
 \param WC_field  mass content field [kg/m3]
 \param N_field Total Number density field [1/m^3]
 \\param limits pnd_field boundaries (indices in p, lat, lon)
 \param scat_meta_array particle meta data for particles
 \param scat_data_start start index for particles handled by this distribution
 \param npart number of particles handled by this distribution
 \param part_string part_species tag for profile/distribution handled here
 \param delim Delimiter string of *part_species* elements
 \param psd_type string with a tag defining the (hydrometeor) scheme
 
 \author Manfred Brath (parts of the function is based on pnd_H11 (of J. Mendrok & D. Kreyling))
 \date 2017-08-01
 
 */
void pnd_fieldMY2 (Tensor4View pnd_field,
                   const Tensor3& WC_field,
                   const Tensor3& N_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
    const String psdname="MY2";
    const Index N_se = scat_meta[scat_species].nelem();
    const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
    ArrayOfIndex intarr;
    Vector diameter_max_unsorted ( N_se, 0.0 );
    Vector diameter_max ( N_se, 0.0 );
    Vector mass ( N_se, 0.0 );
    Vector pnd ( N_se, 0.0 );
    Vector dNdD ( N_se, 0.0 );
    String partfield_name;
    String psd_str;
    String psd_param;
    bool logic_M;
    Numeric N_tot;
    ArrayOfString substrings;
    
    
    //split String and copy to ArrayOfString
    parse_partfield_name( partfield_name, part_string, delim);
    parse_psd_param( psd_param, part_string, delim);
    
    // Check if N_field should be handled as number density field
    // or as mean particle mass field
    psd_param.split(substrings,"_");
    
    if (substrings.size()==3)
    {
        // case: N_field is mean particle mass
        if (substrings[2] == "M")
        {
            psd_str=psd_param.substr(0,psd_param.length()-2);
            
            logic_M=true;
        }
        else
        {
            ostringstream os;
            os << "You use a wrong tag! The last substring of the tag\n"
            << "has to be an uppercase M";
            throw runtime_error( os.str() );
        }
    }
    // case: N_field is total number density
    else
    {
        logic_M=false;
        psd_str=psd_param;
    }
    
    for ( Index i=0; i < N_se; i++ )
    {
        if ( isnan(scat_meta[scat_species][i].mass) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element mass.\n"
            << "But mass is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
        diameter_max_unsorted[i] = ( scat_meta[scat_species][i].diameter_max );
    }
    get_sorted_indexes(intarr, diameter_max_unsorted);
    
    // extract scattering meta data
    for ( Index i=0; i< N_se; i++ )
    {
        if ( isnan(scat_meta[scat_species][intarr[i]].diameter_max) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element diameter_max.\n"
            << "But mass is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
       
        
        diameter_max[i] = scat_meta[scat_species][intarr[i]].diameter_max; // [m]
        
                mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
        
    }
    
    if (diameter_max.nelem() > 0)
        // diameter_max.nelem()=0 implies no selected scattering element for the respective
        // scattering species field. should not occur anymore.
    {
        assert( is_increasing(diameter_max) );
        
        // iteration over all atm. levels
        for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // MY2 requires mass density. If not set, abort calculation.
                    if ( isnan(WC_field ( p, lat, lon )) || isnan(N_field ( p, lat, lon )))
                    {
                        if (logic_M)
                        {
                            ostringstream os;
                            os << "Size distribution " << psdname << " requires knowledge of mass "
                            << "concentration of hydrometeors and mean particle mass.\n"
                            << "At grid point (" << p << ", " << lat << ", " << lon
                            << ") in (p,lat,lon) a NaN value is encountered, "
                            << "i.e. mass concentration amd/or mean particle mass is unknown.";
                            throw runtime_error( os.str() );
                        }
                        else
                        {
                            ostringstream os;
                            os << "Size distribution " << psdname << " requires knowledge of mass "
                            << "concentration of hydrometeors and number density.\n"
                            << "At grid point (" << p << ", " << lat << ", " << lon
                            << ") in (p,lat,lon) a NaN value is encountered, "
                            << "i.e. mass concentration amd/or mean particle mass is unknown.";
                            throw runtime_error( os.str() );
                        }
                    }                                        
                    
                    
                    // A valid IWC value encountered (here, we can also handle negative
                    // values!). Calculating dNdD.
                    else if (WC_field ( p, lat, lon ) != 0. && N_field ( p, lat, lon ) != 0.)
                    {
                        
                        if (logic_M)
                        {
                            //case: N_field is mean particle mass
                            N_tot=WC_field ( p, lat, lon )/N_field ( p, lat, lon );
                        }
                        else
                        {
                            //case: N_fied is total number density
                            N_tot=N_field ( p, lat, lon );
                        }
                        
                        
                        
                        // iteration over all given size bins
                        for ( Index i=0; i<diameter_max.nelem(); i++ ) //loop over number of scattering elements
                        {
                            // calculate particle size distribution for H11
                            // [# m^-3 m^-1]
                            dNdD[i] = psd_MY2( diameter_max[i], N_tot,
                                                  WC_field ( p, lat, lon ),
                                                  psd_str);
                        }
                        
                        // sometimes there is possibility if WC_field is very small
                        // but still greater than zero and N_tot is large that dNdD
                        // could be zero
                        if (dNdD.sum()==0)
                        {
                            
                            for ( Index i = 0; i< N_se; i++ )
                            {
                                pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                           lat-limits[2], lon-limits[4] ) = 0.;
                            }
                        }
                        else
                        {
                            
                            // scale pnds by scale width
                            if (mass.nelem() > 1)
                                bin_integral( pnd, diameter_max, dNdD ); //[# m^-3]
                            else
                                pnd = dNdD;
                            
                            
                            
                            // calculate proper scaling of pnd sum from real IWC and apply
                            chk_pndsum ( pnd, WC_field ( p,lat,lon ), mass,
                                        p, lat, lon, partfield_name, verbosity );
                            
                            // writing pnd vector to wsv pnd_field
                            for ( Index i =0; i< N_se; i++ )
                            {
                                pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                           lat-limits[2], lon-limits[4] ) = pnd[i];
                            }
                        }
                    }
                    
                    // if either WC_field or N_field is zero, we just set pnd_field=0
                    else
                    {
                        for ( Index i = 0; i< N_se; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = 0.;
                        }
                    }
                }
            }
        }
    }
}





/*! Calculates the particle number density field for Cloud liquid water according 
 * to the modified gamma distribution for cloud water inside Geer and Baordo (2014),
 * see table A1
 * Assumptions are: density of particles is constant and particle shape is sphere.
 
 \param pnd_field Particle number density field
 \param LWC_field (liquid cloud water) mass content field [kg/m3]
 \\param limits pnd_field boundaries (indices in p, lat, lon)
 \param scat_meta_array particle meta data for particles
 \param scat_data_start start index for particles handled by this distribution
 \param npart number of particles handled by this distribution
 \param part_string part_species tag for profile/distribution handled here
 \param delim Delimiter string of *part_species* elements
 
 \author Manfred Brath (parts of the function is based on pnd_H11 (of J. Mendrok & D. Kreyling))
 \date 2014-12-02
 
 */
void pnd_fieldMGD_LWC (Tensor4View pnd_field,
                     const Tensor3& LWC_field,
                     const ArrayOfIndex& limits,
                     const ArrayOfArrayOfScatteringMetaData& scat_meta,
                     const Index& scat_species,
                     const String& part_string,
                     const String& delim,
                     const Verbosity& verbosity)
{
    const String psdname="MGD_LWC";
    const Index N_se = scat_meta[scat_species].nelem();
    const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
    ArrayOfIndex intarr;
    Vector diameter_volume_equ_unsorted ( N_se, 0.0 );
    Vector diameter_volume_equ ( N_se, 0.0 );
    Vector mass ( N_se, 0.0 );
    Vector rho ( N_se, 0.0 );
    Numeric mean_rho;
    Numeric max_rho;
    Numeric delta_rho;
    Vector pnd ( N_se, 0.0 );
    Vector dNdD ( N_se, 0.0 );
    String partfield_name;
    
    
    
    //split String and copy to ArrayOfString
    parse_partfield_name( partfield_name, part_string, delim);
    
    for ( Index i=0; i < N_se; i++ )
    {
        if ( isnan(scat_meta[scat_species][i].diameter_volume_equ) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element volume equivalent diameter.\n"
            << "But volume equvalent diameter is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
        diameter_volume_equ_unsorted[i] = ( scat_meta[scat_species][i].diameter_volume_equ );
    }
    get_sorted_indexes(intarr, diameter_volume_equ_unsorted);
    
    // extract scattering meta data
    for ( Index i=0; i< N_se; i++ )
    {
        diameter_volume_equ[i] = scat_meta[scat_species][intarr[i]].diameter_volume_equ; // [m]
        
        if ( isnan(scat_meta[scat_species][intarr[i]].mass) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element mass.\n"
            << "But mass is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
        mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
        rho[i] = mass[i]/(pow(diameter_volume_equ[i],3) *PI/6);
    }
    mean_rho=rho.sum()/(Numeric)rho.nelem();
    
    // checking if the particles have the same density
    max_rho=max(rho);
    delta_rho=(max_rho-min(rho))/max_rho;
    
    if (delta_rho >= 0.1)
    {
        ostringstream os;
        os << "MGD_LWC is valid only for particles with the same or\n"
        "at least almost the same density. The difference between\n"
        "maximum and minimum density must be lower than 10 percent.\n"
        "Your difference is  " << delta_rho << ".\n"
        "Check your scattering particles";
        throw runtime_error( os.str() );
    }

    if (diameter_volume_equ.nelem() > 0)
        // diameter_max.nelem()=0 implies no selected scattering element for the respective
        // scattering species field. should not occur anymore.
    {
        assert( is_increasing(diameter_volume_equ) );

        // iteration over all atm. levels
        for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // MGD_LWC requires mass density. If not set, abort calculation.
                    if ( isnan(LWC_field ( p, lat, lon )) )
                    {
                        ostringstream os;
                        os << "Size distribution " << psdname << " requires knowledge of mass "
                        << "density of cloud liquid water.\n"
                        << "At grid point (" << p << ", " << lat << ", " << lon
                        << ") in (p,lat,lon) a NaN value is encountered, "
                        << "i.e. mass density is unknown.";
                        throw runtime_error( os.str() );
                    }
                    
                    // A valid IWC value encountered (here, we can also handle negative
                    // values!). Calculating dNdD.
                    else if (LWC_field ( p, lat, lon ) != 0.)
                    {
                        // iteration over all given size bins
                        for ( Index i=0; i<diameter_volume_equ.nelem(); i++ ) //loop over number of scattering elements
                        {
                            // calculate particle size distribution for H11
                            // [# m^-3 m^-1]
                            dNdD[i] = LWCtopnd_MGD_LWC( diameter_volume_equ[i],
                                                        mean_rho,
                                                        LWC_field(p, lat, lon) );
                        }

                        // ensure that any particles where produced
                        if ( dNdD.sum() == 0.0 )
                        { // no particles at all were produced. means, there's some
                          // issue with the setup (none of the included particles
                          // produces numerically considerable amount of particles
                          // with this PSD, likely due to all the particles being
                          // either too large or too small)
                          ostringstream os;
                          os <<  psdname << " PSD did not produce any non-zero pnd values "
                             << "(likely results\n"
                             << "from considered particles being too large or too "
                             << "small)\n"
                             << "The problem occured for profile '"<< partfield_name
                             << "' at: " << "p = " << p << ", lat = " << lat
                             << ", lon = " << lon << ".\n";
                          throw runtime_error ( os.str() );
                        }

                        // scale pnds by scale width
                        if (diameter_volume_equ.nelem() > 1)
                            bin_integral( pnd, diameter_volume_equ, dNdD ); //[# m^-3]
                        else
                            pnd = dNdD;
                        
                        // calculate proper scaling of pnd sum from real IWC and apply
                        chk_pndsum ( pnd, LWC_field ( p,lat,lon ), mass,
                                    p, lat, lon, partfield_name, verbosity );
                        
                        // writing pnd vector to wsv pnd_field
                        for ( Index i =0; i< N_se; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = pnd[i];
                        }
                    }
                    
                    // for IWC==0, we just set pnd_field=0
                    else
                    {
                        for ( Index i = 0; i< N_se; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = 0.;
                        }
                    }
                }
            }
        }
    }
}

/*! Calculates the particle number density field for Cloud ice according
 * to the modified gamma distribution for cloud water inside Geer and Baordo (2014),
 * see table A1
 * Assumptions are: density of particles is constant and particle shape is sphere.
 
 \param pnd_field Particle number density field
 \param IWC_field (cloud ice) mass content field [kg/m3]
 \\param limits pnd_field boundaries (indices in p, lat, lon)
 \param scat_meta_array particle meta data for particles
 \param scat_data_start start index for particles handled by this distribution
 \param npart number of particles handled by this distribution
 \param part_string part_species tag for profile/distribution handled here
 \param delim Delimiter string of *part_species* elements
 
 \author Manfred Brath (parts of the function is based on pnd_H11 (of J. Mendrok & D. Kreyling))
 \date 2014-12-02
 
 */
void pnd_fieldMGD_IWC (Tensor4View pnd_field,
                       const Tensor3& IWC_field,
                       const ArrayOfIndex& limits,
                       const ArrayOfArrayOfScatteringMetaData& scat_meta,
                       const Index& scat_species,
                       const String& part_string,
                       const String& delim,
                       const Verbosity& verbosity)
{
    const String psdname="MGD_IWC";
    const Index N_se = scat_meta[scat_species].nelem();
    const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
    ArrayOfIndex intarr;
    Vector diameter_volume_equ_unsorted( N_se, 0.0 );
    Vector diameter_volume_equ( N_se, 0.0 );
    Vector mass( N_se, 0.0 );
    Vector rho( N_se, 0.0 );
    Numeric mean_rho;
    Numeric max_rho;
    Numeric delta_rho;
    Vector pnd( N_se, 0.0 );
    Vector dNdD( N_se, 0.0 );
    String partfield_name;
    
    
    
    //split String and copy to ArrayOfString
    parse_partfield_name( partfield_name, part_string, delim);
    
    for ( Index i=0; i < N_se; i++ )
    {
        if ( isnan(scat_meta[scat_species][i].diameter_max) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element volume equivalent diameter.\n"
            << "But volume equivalent diameter is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
        diameter_volume_equ_unsorted[i] = ( scat_meta[scat_species][i].diameter_volume_equ );
    }
    get_sorted_indexes(intarr, diameter_volume_equ_unsorted);
    
    // extract scattering meta data
    for ( Index i=0; i< N_se; i++ )
    {
        diameter_volume_equ[i] = scat_meta[scat_species][intarr[i]].diameter_volume_equ; // [m]
        
        if ( isnan(scat_meta[scat_species][intarr[i]].mass) )
        {
            ostringstream os;
            os << "Use of size distribution " << psdname << " (as requested for\n"
            << "scattering species #" << scat_species << ")\n"
            << "requires knowledge of scattering element mass.\n"
            << "But mass is not given for scattering elements #"
            << i << "!";
            throw runtime_error( os.str() );
        }
        mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
        rho[i] = mass[i]/(pow(diameter_volume_equ[i],3) *PI/6);
    }
    mean_rho=rho.sum()/(Numeric)rho.nelem();
    
    // checking if the particles have the same density
    max_rho=max(rho);
    delta_rho=(max_rho-min(rho))/max_rho;
    
    if (delta_rho >= 0.1)
    {
        ostringstream os;
        os << "MGD_IWC is valid only for particles with the same\n"
        "at least almost the same density. The difference between\n"
        " maximum and minimum density must be lower than 10 percent.\n"
        "Your difference is  " << delta_rho << ".\n"
        "Check your scattering particles";
        throw runtime_error( os.str() );
    }
    
    
    if (diameter_volume_equ.nelem() > 0)
        // diameter_max.nelem()=0 implies no selected scattering element for the respective
        // scattering species field. should not occur anymore.
    {
        assert( is_increasing(diameter_volume_equ) );

        // iteration over all atm. levels
        for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // MGD requires mass density. If not set, abort calculation.
                    if ( isnan(IWC_field ( p, lat, lon )) )
                    {
                        ostringstream os;
                        os << "Size distribution " << psdname << " requires knowledge of mass "
                        << "density of cloud liquid water.\n"
                        << "At grid point (" << p << ", " << lat << ", " << lon
                        << ") in (p,lat,lon) a NaN value is encountered, "
                        << "i.e. mass density is unknown.";
                        throw runtime_error( os.str() );
                    }
                    
                    // A valid IWC value encountered (here, we can also handle negative
                    // values!). Calculating dNdD.
                    else if (IWC_field ( p, lat, lon ) != 0.)
                    {
                        // iteration over all given size bins
                        for ( Index i=0; i<diameter_volume_equ.nelem(); i++ ) //loop over number of scattering elements
                        {
                            // calculate particle size distribution for H11
                            // [# m^-3 m^-1]
                            dNdD[i] = IWCtopnd_MGD_IWC( diameter_volume_equ[i],
                                                        mean_rho,
                                                        IWC_field(p, lat, lon) );
                        }

                        // ensure that any particles where produced
                        if ( dNdD.sum() == 0.0 )
                        { // no particles at all were produced. means, there's some
                          // issue with the setup (none of the included particles
                          // produces numerically considerable amount of particles
                          // with this PSD, likely due to all the particles being
                          // either too large or too small)
                          ostringstream os;
                          os <<  psdname << " PSD did not produce any non-zero pnd values "
                             << "(likely results\n"
                             << "from considered particles being too large or too "
                             << "small)\n"
                             << "The problem occured for profile '"<< partfield_name
                             << "' at: " << "p = " << p << ", lat = " << lat
                             << ", lon = " << lon << ".\n";
                          throw runtime_error ( os.str() );
                        }

                        // scale pnds by scale width
                        if (diameter_volume_equ.nelem() > 1)
                            bin_integral( pnd, diameter_volume_equ, dNdD ); //[# m^-3]
                        else
                            pnd = dNdD;
                        
                        // calculate proper scaling of pnd sum from real IWC and apply
                        chk_pndsum ( pnd, IWC_field ( p,lat,lon ), mass,
                                    p, lat, lon, partfield_name, verbosity );
                        
                        // writing pnd vector to wsv pnd_field
                        for ( Index i =0; i< N_se; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = pnd[i];
                        }
                    }
                    
                    // for IWC==0, we just set pnd_field=0
                    else
                    {
                        for ( Index i = 0; i< N_se; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = 0.;
                        }
                    }
                }
            }
        }
    }
}


/*! Calculates the particle number density field for Marshall and Palmer (1948)
    size distribution. To be used for precipitation, particularly rain.

    \return pnd_field   Particle number density field
    \param PR_field     precipitation rate (mass flux) field [kg/(m2*s)]
    \param limits       pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta    scattering meta data for all scattering elements
    \param scat_species array index of scattering species handled by this distribution
    \param part_string  scat_species tag for profile/distribution handled here
    \param delim        delimiter string of *scat_species* elements
  
  \author Jana Mendrok
  \date 2012-04-04

*/
void pnd_fieldMP48 (Tensor4View pnd_field,
                    const Tensor3& PR_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfArrayOfScatteringMetaData& scat_meta,
                    const Index& scat_species,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity)
{
  const String psdname="MP48";
  const Index N_se = scat_meta[scat_species].nelem();
  const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
  ArrayOfIndex intarr;
  Vector mass_unsorted ( N_se, 0.0 );
  Vector mass ( N_se, 0.0 );
  Vector diameter_melted_equivalent ( N_se, 0.0 );
  Vector vol ( N_se, 0.0 );
  Vector pnd ( N_se, 0.0 );
  Vector dNdD ( N_se, 0.0 );

  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].mass) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element mass.\n"
             << "But mass is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      mass_unsorted[i] = ( scat_meta[scat_species][i].mass );
    }
  get_sorted_indexes(intarr, mass_unsorted);
	
  // extract scattering meta data
  for ( Index i=0; i< N_se; i++ )
  {
      mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
      diameter_melted_equivalent[i] = pow(6.*mass[i]/PI/DENSITY_OF_WATER,1./3.); // [m]

      if ( isnan(scat_meta[scat_species][intarr[i]].mass) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element volume.\n"
             << "But volume is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      vol[i] = PI/6. *
        pow(scat_meta[scat_species][intarr[i]].diameter_volume_equ,3.); // [m3]
  }

  // conversion factor from PR [kg/(m^2*s)] to PR[mm/hr]
  /* for the conversion we need to know the mean density of distribution:
     rho_mean = mass_total/vol_total
              = int(mass[D]*dNdD[D])dD/int(vol[D]*dNdD[D])dD

     However, we do not yet know dNdD[D], which in turn is dependent on rho[D].
      proper solution requires iterative approach (does it?), but is it really
      worth to implement this? at least rain drops density should not really
      vary with particle size. might be different for snow/hail/graupel.
     So, here we set conversion factor to pseudo-PR[mm/hr]*[kg/m^3] and divide by
      density later on.
  */
  const Numeric convfac=3.6e6; // [PR [kg/(m^2*s)] to PR[mm/hr]*[kg/m3]
  Numeric fac, rho_mean, vol_total, mass_total, tPR;
      
  // set parameterisation constants here instead of in PRtopnd_MP48, since we
  // also need them for PR to PWC conversion
  const Numeric N0 = 0.08*1e8; // [#/cm^3/cm] converted to [#/m^3/m]
  const Numeric lambda_fac = 41.*1e2; // [cm^-1] converted to [m^-1] to fit d[m]
  const Numeric lambda_exp = -0.21;
  Numeric PWC, lambda = NAN;

  if (diameter_melted_equivalent.nelem() > 0)
  // diameter_melted_equivalent.nelem()=0 implies no selected scattering elements for the respective
  // scattering species field. should not occur.
  {
      assert( is_increasing(diameter_melted_equivalent) );

      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // A valid PR value encountered (here, we can also handle negative
            // values!). Calculating dNdD.
            if (PR_field ( p, lat, lon ) > 0.)
            {
              Index n_it = 0;

              // initializing for proper start of while loop.
              rho_mean = 0.; //this has to be off the real value in order to go
                             //into while loop
              mass_total = mass.sum();
              vol_total = vol.sum();
              while (abs( rho_mean/(mass_total/vol_total)-1.)>1e-3)
              // did bulk mean density change?
              {
                if (n_it>10)
                {
                    ostringstream os;
                    os<< "ERROR: Bulk mean density for " << psdname << " distribution is"
                      << " not converging.\n"
                      << "fractional change: "
                      << abs( rho_mean/(mass_total/vol_total)-1.) << "\n"
                      << "at p: " << p << " lat: " << lat << " lon" << lon
                      << " for PR=" << PR_field ( p, lat, lon ) << "\n";
                    throw runtime_error ( os.str() );
                }
                rho_mean = mass_total/vol_total;
                fac = convfac/rho_mean;

                // do PR [kg/(m^2*s)] to PR [mm/hr] conversion
                tPR = PR_field ( p, lat, lon ) * fac;
                // get slope of distribution [m^-1]
                lambda = lambda_fac * pow(tPR,lambda_exp);

                // derive particle number density for all given sizes
                for ( Index i=0; i<diameter_melted_equivalent.nelem(); i++ )
                {
                    // calculate particle size distribution with MP48
                    // output: [# m^-3 m^-1]
                    dNdD[i] = PRtopnd_MP48 ( tPR, diameter_melted_equivalent[i]);
                }

                // ensure that any particles where produced
                if ( dNdD.sum() == 0.0 )
                  { // no particles at all were produced. means, there's some
                    // issue with the setup (none of the included particles
                    // produces numerically considerable amount of particles
                    // with this PSD, likely due to all the particles being
                    // either too large or too small)
                    ostringstream os;
                    os <<  psdname << " PSD did not produce any non-zero pnd values "
                       << "(likely results\n"
                       << "from considered particles being too large or too "
                       << "small)\n"
                       << "The problem occured for profile '"<< partfield_name
                       << "' at: " << "p = " << p << ", lat = " << lat
                       << ", lon = " << lon << ".\n";
                    throw runtime_error ( os.str() );
                  }

                // scale pnds by bin width
                if (diameter_melted_equivalent.nelem() > 1)
                    bin_integral( pnd, diameter_melted_equivalent, dNdD );
                else
                    pnd = dNdD;

                // derive mass and volume over whole size distribution for
                // updated mean density
                mass_total = vol_total = 0.;
                for ( Index i=0; i<diameter_melted_equivalent.nelem(); i++ )
                {
                    mass_total += mass[i]*pnd[i];
                    vol_total += vol[i]*pnd[i];
                }
                n_it++;
              }

              // Make sure lambda was initialized in the while loop
              assert(!isnan(lambda));

              // calculate error of pnd sum and real XWC
              PWC = rho_mean*PI*N0 / pow(lambda,4.);
              chk_pndsum ( pnd, PWC, mass,
                           p, lat, lon, partfield_name, verbosity );
	    
              // writing pnd vector to wsv pnd_field
              for ( Index i = 0; i< N_se; i++ )
              {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = pnd[i];
              }
            }

            // MP48 requires mass flux (actually, it's precip rate. but we can
            // convert these). If not set, abort calculation.
            else if ( isnan(PR_field ( p, lat, lon )) )
              {
                ostringstream os;
                os << "Size distribution " << psdname << " requires knowledge of mass "
                   << "flux of atmospheric (liquid/solid) water.\n"
                   << "At grid point (" << p << ", " << lat << ", " << lon
                   << ") in (p,lat,lon) a NaN value is encountered, "
                   << "i.e. mass flux is unknown.";
                throw runtime_error( os.str() );
              }

            // MP48 PSD is parameterized in PR and can not handle negative
            // numbers, hence abort.
            else if (PR_field ( p, lat, lon ) < 0.)
              {
                ostringstream os;
                os << "Size distribution " << psdname << " is parametrized in "
                   << "precipitation rate (mass flux).\n"
                   << "It can not handle negative values like PR="
                   << PR_field (p,lat,lon) << " kg/m2/s\n"
                   << "found at grid point ("
                   << p << ", " << lat << ", " << lon << ")";
                throw runtime_error( os.str() );
              }

            // for PR==0, we just set pnd_field=0
            else
                for ( Index i = 0; i< N_se; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = 0.;
                }
          }
        }
      }
  }
}


/*! Calculates the particle number density field for Wang et al. (2016) size
    distribution. To be used for rain.

    \return pnd_field   Particle number density field
    \param RWC_field    mass content (rain) field [kg/m3]
    \param limits       pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta    scattering meta data for all scattering elements
    \param scat_species array index of scattering species handled by this distribution
    \param part_string  scat_species tag for profile/distribution handled here
    \param delim        delimiter string of *scat_species* elements
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2017-06-13

*/
void pnd_fieldW16 (Tensor4View pnd_field,
                    const Tensor3& RWC_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfArrayOfScatteringMetaData& scat_meta,
                    const Index& scat_species,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity)
{
  const String psdname="W16";
  const Index N_se = scat_meta[scat_species].nelem();
  const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
  ArrayOfIndex intarr;
  Vector mass_unsorted ( N_se, 0.0 );
  Vector mass ( N_se, 0.0 );
  Vector diameter_mass_equivalent ( N_se, 0.0 );
  Vector pnd ( N_se, 0.0 );
  Vector dNdD ( N_se, 0.0 );

  String psd_param;
  String partfield_name;
  ArrayOfString psd_options;

  //split String and copy to ArrayOfString
  parse_psd_param( psd_param, part_string, delim);
  parse_partfield_name( partfield_name, part_string, delim);

  bool robust = false;
  parse_psd_options( psd_options, part_string, delim);
  for ( Index i=0; i<psd_options.nelem(); i++ )
    robust = (robust || (psd_options[i]=="robust") );

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].mass) ||
           isnan(scat_meta[scat_species][i].diameter_volume_equ) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element mass.\n"
             << "But mass is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      mass_unsorted[i] = ( scat_meta[scat_species][i].mass );
    }
  get_sorted_indexes(intarr, mass_unsorted);

  // extract scattering meta data
  for ( Index i=0; i< N_se; i++ )
  {
      mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
      if ( isnan(scat_meta[scat_species][intarr[i]].diameter_volume_equ) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element volume equivalent diameter.\n"
             << "But volume equivalent diameter is not given for scattering elements #"
             << intarr[i] << "!";
          throw runtime_error( os.str() );
        }
      diameter_mass_equivalent[i] = scat_meta[scat_species][intarr[i]].diameter_volume_equ;
  }
  
  if (mass.nelem() > 0)
  // mass.nelem()=0 implies no selected scattering element for the respective
  // scattering species field. should not occur.
  {
      assert( is_increasing(diameter_mass_equivalent) );
      
      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // a valid RWC value encountered. Calculating dNdD.
            if (RWC_field ( p, lat, lon ) > 0.)
              {
                psd_rain_W16 ( dNdD, diameter_mass_equivalent,
                               RWC_field(p,lat,lon) );

                // ensure that any particles where produced
                if ( dNdD.sum() == 0.0 )
                  { // no particles at all were produced. means, there's some
                    // issue with the setup (none of the included particles
                    // produces numerically considerable amount of particles
                    // with this PSD, likely due to all the particles being
                    // either too large or too small)
                    ostringstream os;
                    os <<  psdname << " PSD did not produce any non-zero pnd values "
                       << "(likely results\n"
                       << "from considered particles being too large or too "
                       << "small)\n"
                       << "The problem occured for profile '"<< partfield_name
                       << "' at: " << "p = " << p << ", lat = " << lat
                       << ", lon = " << lon << ".\n";
                    throw runtime_error ( os.str() );
                  }

                // derive pnds from scaling dNdD with bin width
                if (diameter_mass_equivalent.nelem() > 1)
                  bin_integral( pnd, diameter_mass_equivalent, dNdD );
                else
                  pnd = dNdD;
	    
                // calculate error of pnd sum and real XWC
                chk_pndsum ( pnd, RWC_field ( p,lat,lon ), mass,
                             p, lat, lon, partfield_name, verbosity );
	    
                // writing pnd vector to WSV pnd_field
                for ( Index i = 0; i< N_se; i++ )
                  {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                  }
              }

            // W16 requires mass density. If not set, abort calculation.
            else if ( isnan(RWC_field ( p, lat, lon )) )
              {
                ostringstream os;
                os << "Size distribution " << psdname << " requires knowledge of mass "
                   << "density of rain.\n"
                   << "At grid point (" << p << ", " << lat << ", " << lon
                   << ") in (p,lat,lon) a NaN value is encountered, "
                   << "i.e. mass density is unknown.";
                throw runtime_error( os.str() );
              }

            // W16 is parameterized in RWC and does not handle negative
            // numbers, hence abort.
            else if (RWC_field ( p, lat, lon ) < 0.)
              {
                ostringstream os;
                os << "Size distribution " << psdname
                   << " is parametrized in rain mass content.\n"
                   << "It can not handle negative values like RWC="
                   << RWC_field (p,lat,lon) << " kg/m3\n"
                   << "found at grid point ("
                   << p << ", " << lat << ", " << lon << ")";
                throw runtime_error( os.str() );
              }

            // for RWC==0, we just set pnd_field=0
            else
                for ( Index i = 0; i< N_se; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = 0.;
                }
          }
        }
      }
  }
}


/*! Calculates the particle number density field for Hess et al. (1998)
    size distribution, namely the continental stratus case. To be used for
    liquid clouds.

    \return pnd_field   Particle number density field
    \param LWC_field    mass content (liquid water) field [kg/m3]
    \param limits       pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta    scattering meta data for all scattering elements
    \param scat_species array index of scattering species handled by this distribution
    \param part_string  scat_species tag for profile/distribution handled here
    \param delim        delimiter string of *scat_species* elements
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2012-04-04

*/
void pnd_fieldH98 (Tensor4View pnd_field,
                   const Tensor3& LWC_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
  const String psdname="H98";
  const Index N_se = scat_meta[scat_species].nelem();
  const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
  ArrayOfIndex intarr;
  Vector diameter_volume_equivalent_unsorted ( N_se, 0.0 );
  Vector diameter_volume_equivalent ( N_se, 0.0 );
  Vector mass ( N_se, 0.0 );
  Vector radius ( N_se, 0.0 );
  Vector pnd ( N_se, 0.0 );
  Vector dNdr ( N_se, 0.0 );

  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].diameter_volume_equ) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element volume equivalent diameter.\n"
             << "But volume equivalent diameter is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      diameter_volume_equivalent_unsorted[i] = ( scat_meta[scat_species][i].diameter_volume_equ );
    }
  get_sorted_indexes(intarr, diameter_volume_equivalent_unsorted);
      
  // extract scattering meta data
  for ( Index i=0; i< N_se; i++ )
  {
      diameter_volume_equivalent[i]= scat_meta[scat_species][intarr[i]].diameter_volume_equ; // [m]
      radius[i] = 0.5 * diameter_volume_equivalent[i]; // [m]

      if ( isnan(scat_meta[scat_species][intarr[i]].mass) )
        {
          ostringstream os;
          os << "Use of size distribution " << psdname << " (as requested for\n"
             << "scattering species #" << scat_species << ")\n"
             << "requires knowledge of scattering element mass.\n"
             << "But mass is not given for scattering elements #"
             << i << "!";
          throw runtime_error( os.str() );
        }
      mass[i] = scat_meta[scat_species][intarr[i]].mass; // [kg]
  }

  if (radius.nelem() > 0)
  // radius.nelem()=0 implies no selected scattering elements for the respective
  // scattering species field. should not occur anymore
  {
      assert( is_increasing(radius) );

      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // H98 requires mass density. If not set, abort calculation.
            if ( isnan(LWC_field ( p, lat, lon )) )
              {
                ostringstream os;
                os << "Size distribution " << psdname << " requires knowledge of mass "
                   << "density of atmospheric liquid water.\n"
                   << "At grid point (" << p << ", " << lat << ", " << lon
                   << ") in (p,lat,lon) a NaN value is encountered, "
                   << "i.e. mass density is unknown.";
                throw runtime_error( os.str() );
              }

            // A valid LWC value encountered (here, we can also handle negative
            // values!). Calculating dNdD.
            if (LWC_field ( p,lat,lon ) != 0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<radius.nelem(); i++ ) //loop over number of scattering elements
                {
                    // calculate particle size distribution for liquid
                    // [# m^-3 m^-1]
                    dNdr[i] = LWCtopnd ( LWC_field ( p,lat,lon ), radius[i] );
                }

                // scale pnds by scale width. output: [# m^-3]
                if (radius.nelem() > 1)
                    bin_integral( pnd, radius, dNdr );
                else
                    pnd = dNdr;
	    
                // calculate error of pnd sum and real XWC
                chk_pndsum ( pnd, LWC_field ( p,lat,lon ), mass,
                             p, lat, lon, partfield_name, verbosity );

                // writing pnd vector to wsv pnd_field
                for ( Index i =0; i< N_se; i++ )
                {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                    //dlwc[q] = pnd2[q]*vol[q]*rho[q];
                }
            }

            // for LWC==0, we just set pnd_field=0
            else
            {
                for ( Index i = 0; i< N_se; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = 0.;
                }
            }
          }
        }
      }
  }
}


/*! Calculates particle size distribution according to a general modified gamma
 *  distribution
 *  
 *  Uses all four free parameters (N0, mu, lambda, gamma) to calculate
 *    psd(D) = N0 * D^mu * exp( -lambda * D^gamma )
 *  
 *  Reference: Petty & Huang, JAS, (2011).
 *  Ported from Atmlab.

    \param psd     particle number density per size interval [#/m3*m]
    \param diameter  size of the scattering elements (volume equivalent
                       diameter) [m]
    \param N0      Intercept parameter. See above. [ ? ]
    \param mu      See above. [ ? ]
    \param lambda  See above. [ ? ]
    \param gamma   See above. [ ? ]
  
  \author Jana Mendrok, Patrick Eriksson
  \date 2017-06-07

*/
void psd_general_MGD ( Vector& psd,
                   const Vector& diameter,
                   const Numeric& N0,
                   const Numeric& mu,
                   const Numeric& lambda,
                   const Numeric& gamma )
{
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

  // ensure numerical stability
  if( (mu+1)/gamma <= 0. )
    {
      ostringstream os;
      os << "(mu+1) / gamma must be > 0.";
      throw runtime_error( os.str() );
    }

  // skip calculation if N0 is 0.0
  if ( N0 == 0.0 )
  {
    return;
  }
  assert (N0>0.);

  if( gamma == 1 )
    {
      if( mu == 0 )
        // Exponential distribution
        for( Index iD=0; iD<nD; iD++ )
          {
            psd[iD] = N0 *
                      exp( -lambda*diameter[iD] );
          }
      else
        // Gamma distribution
        for( Index iD=0; iD<nD; iD++ )
          {
            psd[iD] = N0 * pow( diameter[iD], mu ) *
                      exp( -lambda*diameter[iD] );
          }
    }
  else
    {
      // Complete MGD
      for( Index iD=0; iD<nD; iD++ )
        {
          psd[iD] = N0 * pow( diameter[iD], mu ) * 
                    exp( -lambda*pow( diameter[iD], gamma ) );
        }
    }

  //for( Index iD=0; iD<nD; iD++ )
  //  if( isnan(psd[iD]) )
  //    psd[iD] = 0.0;
}


/*! Calculates particle size distribution of cloud ice using MH97 parametrization.
 *  
 *  Handles a vector of sizes at a time. Implicitly assumes particles of water
 *  ice. Strictly requires IWC and T to be positive, i.e. calling method needs
 *  to ensure this.
 *  
 *  Adapted from the 'old' IWCtopnd_MH97.

    \param psd     particle number density per size interval [#/m3*m]
    \param diameter  size of the scattering elements (supposed to be mass (aka
                       volume) equivalent diameter of pure ice particle) [m]
    \param iwc     atmospheric ice water content [kg/m3]
    \param t       atmospheric temperature [K]
    \param noisy   flag whether to add noise onto PSD parameters according to
                     their reported error statistics
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2017-06-07

*/
void psd_cloudice_MH97 ( Vector& psd,
                   const Vector& diameter,
                   const Numeric& iwc,
                   const Numeric& t,
                   const bool noisy )
{
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

  assert (t>0.);

  // skip calculation if IWC is 0.0
  if ( iwc == 0.0 )
  {
    return;
  }
  assert (iwc>0.);

  // convert m to microns
  Vector d_um(nD);
  for ( Index iD=0; iD<nD; iD++ )
    d_um[iD] = 1e6 * diameter[iD];
  //convert T from Kelvin to Celsius
  Numeric Tc = t-273.15;

  //[kg/m3] -> [g/m3] as used by parameterisation
  Numeric ciwc = iwc*1e3;
  Numeric cdensity = DENSITY_OF_ICE*1e3;

  Numeric sig_a=0., sig_b1=0.;
  Numeric sig_b2=0., sig_m=0.;
  Numeric sig_aamu=0., sig_bamu=0.;
  Numeric sig_abmu=0., sig_bbmu=0.;
  Numeric sig_aasigma=0., sig_basigma=0;
  Numeric sig_absigma=0., sig_bbsigma=0.;

  if ( noisy )
  {
    Rng rng;                      //Random Number generator
    Index mc_seed;
    mc_seed = (Index)time(NULL);
    rng.seed(mc_seed, Verbosity());

    sig_a=0.068, sig_b1=0.054;
    sig_b2=5.5e-3, sig_m=0.0029;
    sig_aamu=0.02, sig_bamu=0.0005;
    sig_abmu=0.023, sig_bbmu=0.5e-3;
    sig_aasigma=0.02, sig_basigma=0.5e-3;
    sig_absigma=0.023, sig_bbsigma=4.7e-4;

    sig_a = ran_gaussian(rng,sig_a);
    sig_b1 = ran_gaussian(rng,sig_b1);
    sig_b2 = ran_gaussian(rng,sig_b2);
    sig_m = ran_gaussian(rng,sig_m);
    sig_aamu = ran_gaussian(rng,sig_aamu);
    sig_bamu = ran_gaussian(rng,sig_bamu);
    sig_abmu = ran_gaussian(rng,sig_abmu);
    sig_bbmu = ran_gaussian(rng,sig_bbmu);
    sig_aasigma = ran_gaussian(rng,sig_aasigma);
    sig_basigma = ran_gaussian(rng,sig_basigma);
    sig_absigma = ran_gaussian(rng,sig_absigma);
    sig_bbsigma = ran_gaussian(rng,sig_bbsigma);
  }

  // Split IWC in small and large size modes

  // Calculate IWC in each mode
  Numeric a=0.252+sig_a; //g/m^3
  Numeric b1=0.837+sig_b1;
  Numeric IWCs100=min ( ciwc,a*pow ( ciwc,b1 ) );
  Numeric IWCl100=ciwc-IWCs100;


  // Gamma distribution component (small mode)

  Numeric b2=-4.99e-3+sig_b2; //micron^-1
  Numeric m=0.0494+sig_m; //micron^-1
  Numeric alphas100=b2-m*log10 ( IWCs100 ); //micron^-1

  // alpha, and hence dNdD1, becomes NaN if IWC>0.
  // this should be ensured to not happen before.
  //
  // alpha, and hence dNdD1, becomes negative for large IWC.
  // towards this limit, particles anyway get larger than 100um, i.e., outside
  // the size region the small-particle mode gamma distrib is intended for.
  // hence, leave dNdD1 at 0 for those cases.
  Vector dNdD1(nD, 0.);
  if (alphas100>0.)
  {
    Numeric Ns100 = 6*IWCs100 * pow ( alphas100,5. ) /
                    ( PI*cdensity*gamma_func ( 5. ) );//micron^-5
    for ( Index iD=0; iD<nD; iD++ )
      dNdD1[iD] = 1e18 * Ns100*d_um[iD] *
                  exp ( -alphas100*d_um[iD] ); //micron^-4 -> m^-3 micron^-1
  }


  // Log normal distribution component (large mode)

  // for small IWC, IWCtotal==IWC<100 & IWC>100=0.
  // this will give dNdD2=NaN. avoid that by explicitly setting to 0
  Vector dNdD2(nD, 0.);
  if (IWCl100>0.)
  {
    //FIXME: Do we need to ensure mul100>0 and sigmal100>0?

    Numeric aamu=5.20+sig_aamu;
    Numeric bamu=0.0013+sig_bamu;
    Numeric abmu=0.026+sig_abmu;
    Numeric bbmu=-1.2e-3+sig_bbmu;
    Numeric amu=aamu+bamu*Tc;
    Numeric bmu=abmu+bbmu*Tc;
    Numeric mul100=amu+bmu*log10 ( IWCl100 );

    Numeric aasigma=0.47+sig_aasigma;
    Numeric basigma=2.1e-3+sig_basigma;
    Numeric absigma=0.018+sig_absigma;
    Numeric bbsigma=-2.1e-4+sig_bbsigma;
    Numeric asigma=aasigma+basigma*Tc;
    Numeric bsigma=absigma+bbsigma*Tc;
    Numeric sigmal100=asigma+bsigma*log10 ( IWCl100 );

    if ( (mul100>0.) & (sigmal100>0.) )
    {
      Numeric a1 = 6*IWCl100; //g/m^3
      Numeric a2_fac = pow ( PI,3./2. ) * cdensity*sqrt(2) *
                   exp( 3*mul100+9./2. * pow ( sigmal100,2 ) ) * 
                   sigmal100 * pow ( 1.,3 );
      //a2 = a2_fac * d_um; //g/m^3/micron^4
      for ( Index iD=0; iD<nD; iD++ )
        dNdD2[iD] = 1e18 * a1 / (a2_fac*d_um[iD]) *
                    exp ( -0.5 * pow ( ( log ( d_um[iD] )-mul100 ) /sigmal100,2 ) );
                    //micron^-4 -> m^-3 micron^-1
    }
  }

  for( Index iD=0; iD<nD; iD++ )
    {
      // FIXME: Do we still need this check here? Non-NaN of each mode should
      // now be ensure by the checks/if-loops for each of the modes separately.
      // I, JM, think (and hope).
      //if ( !isnan(dNdD1[iD]) && !isnan(dNdD2[iD]) )
        psd[iD] = ( dNdD1[iD]+dNdD2[iD] ) *1e6; // m^-3 m^-1
    }
}


/*! Calculates particle size distribution of tropical or midlatitude snow using
 *  F07 parametrization.
 *
 *  Handles a vector of sizes at a time.
 *  Strictly requires SWC and T to be positive and regime to be either "TR" or
 *  "ML", i.e. calling methods need to ensure these.
 *  No further limitations on the allowed temperatures here. Strictly valid it's
 *  only within -60<=t<=0C, the measured t-range the parametrization is based
 *  on. However, this is left to be handled by the calling methods.
 *
 *  Adapted from the 'old' IWCtopnd_F07TR/ML.
 
    \param psd     particle number density per size interval [#/m3*m]
    \param diameter  size of the scattering elements (supposed to be maximum
                       diameter of the ice particles) [m]
    \param swc     atmospheric snow water content [kg/m^3]
    \param t       atmospheric temperature [K]
    \param alpha   mass-dimension relationship scaling factor
                     (m=alpha*(Dmax/D0)^beta) [kg]
    \param beta    mass-dimension relationship exponent [-]
    \param regime  parametrization regime to apply (TR=tropical, ML=midlatitude)
 
 \author Manfred Brath, Jana Mendrok
 \date 2017-06-13
 
 */
void psd_snow_F07 ( Vector& psd,
                   const Vector& diameter,
                   const Numeric& swc,
                   const Numeric& t,
                   const Numeric alpha,
                   const Numeric beta,
                   const String& regime )
{
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

  assert (t>0.);

  // skip calculation if SWC is 0.0
  if ( swc == 0.0 )
  {
    return;
  }
  assert (swc>0.);

  Numeric An, Bn, Cn;
  Numeric M2, Mn, M2Mn;
  Numeric base, pp;
  Numeric x, phi23;
  //Numeric dN;

  Vector q(5);

  assert( (regime=="TR") || (regime=="ML") );

  //factors of phi23
  if( regime=="TR" )
    q={152., -12.4, 3.28, -0.78, -1.94};
  else // if( regime=="ML" )
    q={141., -16.8, 102.,  2.07, -4.82};

  //Factors of factors of the moment estimation parametrization
  Vector Aq{13.6,-7.76,0.479};
  Vector Bq{-0.0361,0.0151,0.00149};
  Vector Cq{0.807,0.00581,0.0457};

  //convert T from Kelvin to Celsius
  Numeric Tc = t-273.15;

  // estimate second moment
  M2 = swc/alpha;
  if (beta!=2)
    {
      // calculate factors of the moment estimation parametrization
      An = exp(Aq[0]+Aq[1]*beta+Aq[2]*pow(beta,2));
      Bn = Bq[0]+Bq[1]*beta+Bq[2]*pow(beta,2);
      Cn = Cq[0]+Cq[1]*beta+Cq[2]*pow(beta,2);

      base = M2*exp(-Bn*Tc)/An;
      pp = 1./(Cn);

      M2 = pow(base,pp);
    }

  // order of the moment parametrization
  const Numeric n=3;

  // calculate factors of the moment estimation parametrization
  An = exp(Aq[0]+Aq[1]*n+Aq[2]*pow(n,2));
  Bn = Bq[0]+Bq[1]*n+Bq[2]*pow(n,2);
  Cn = Cq[0]+Cq[1]*n+Cq[2]*pow(n,2);

  // moment parametrization
  Mn = An*exp(Bn*Tc)*pow(M2,Cn);

  M2Mn = pow(M2,4.)/pow(Mn,3.);

  for( Index iD=0; iD<nD; iD++ )
    {
      // define x
      x = diameter[iD]*M2/Mn;

      // characteristic function
      phi23 = q[0]*exp(q[1]*x)+q[2]*pow(x,q[3])*exp(q[4]*x);
    
      // distribution function
      //dN = phi23*M2Mn;

      //if ( !isnan(psd[dN]) )
      //  psd[iD] = dN;

      // set psd directly. Non-NaN should be (and is, hopefully) ensured by
      // checks above (if we encounter further NaN, that should e handled above.
      // which intermediate quantities make problems? at which parametrisation
      // values, ie WC, T, alpha, beta?).
      psd[iD] = phi23*M2Mn;
    }
}


/*! Calculates particle size distribution of (stratiform) rain using Wang16
 *  parametrization.
 *  
 *  Uses rain water water content. PSD follows an exponential distribution.
 *  Handles a vector of sizes at a time.
 *  
 *  Reference: Wang et al., 2016, "Investigation of liquid cloud microphysical
 *  properties of deep convective systems: 1. Parameterization raindrop size
 *  distribution and its application for stratiform rain estimation".
 *  Ported from CloudArts matlab implementation.

    \param psd     particle number density per size interval [#/m3*m]
    \param diameter  size of the scattering elements (volume equivalent
                       diameter) [m]
    \param rwc     atmospheric rain water content [kg/m3]
  
  \author Jana Mendrok, Bengt Rydberg
  \date 2017-06-07

*/
void psd_rain_W16 ( Vector& psd,
                   const Vector& diameter,
                   const Numeric& rwc )
{
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

  // skip calculation if RWC is 0.0
  if ( rwc == 0.0 )
  {
    return;
  }
  assert (rwc>0.);

  // a and b relates N0 to lambda N0 = a*lambda^b
  Numeric a = 0.000141;
  Numeric b = 1.49;
  Numeric c1 = DENSITY_OF_WATER * PI / 6;
  Numeric base = c1 / rwc * a * tgamma(4);
  Numeric exponent = 1. / (4 - b);
  Numeric lambda = 1e2 * pow( base, exponent );
  Numeric N0 = 1e8 * a * pow( lambda, b );

  //psd_general_MGD( psd, N0, 0, lambda, 1 );
  for( Index iD=0; iD<nD; iD++ )
    {
      psd[iD] = N0 * exp( -lambda*diameter[iD] );
    }
}



/*! Calculates particle size distribution using H11 parametrization.
 *  Each diameter of the scattering elements is a node in the distribution.
 *  One call of this function calculates number density for one scattering
 *  element.  

    \return dNdD particle number density per diameter interval [#/m3/m]
          
    \param diameter_max  maximum diameter of scattering scattering element [m]
    \param t     atmospheric temperature [K]
  
  \author Daniel Kreyling
  \date 2011-10-28

*/
Numeric IWCtopnd_H11 ( const Numeric diameter_max,
                       const Numeric t)
{  
  Numeric dNdD;
  Numeric la;
  Numeric mu;

  // convert m to cm
  Numeric dmax = diameter_max * 1e2;
  //convert T from Kelvin to Celsius
  Numeric T = t-273.15;
  //choose parametrization depending on T
  if ( T >= -56. )
  {
    la= 12.13 * exp( -0.055*T );
  }
  else
  {
    la= 0.83 * exp( -0.103*T );
  }
  if ( T >= -68. )
  {
    mu= -0.57 - 0.028*T;
  }
  else
  {
    mu= -30.93 - 0.472*T;
  }
 
  //Distribution function H11

  dNdD = pow( dmax, mu ) * exp ( -la * dmax );

  if (isnan(dNdD)) dNdD = 0.0;
  return dNdD;
}



/*! Calculates particle size distribution using H13 parametrization.
 *  Each diameter of the scattering elements is a node in the distribution.
 *  One call of this function calculates number density for one scattering
 *  element.  

    \return dNdD particle number density per diameter interval [#/m3/m]
          
    \param diameter_max  maximum diameter of scattering scattering element [m]
    \param t     atmospheric temperature [K]
  
  \author Johan Strandgren, Daniel Kreyling  
  \date 2013-08-26

*/
Numeric IWCtopnd_H13 ( const Numeric diameter_max,
                       const Numeric t)
{  
  Numeric dNdD;
  Numeric la;
  Numeric mu;

  // convert m to cm
  Numeric dmax = diameter_max * 1e2;
  //convert T from Kelvin to Celsius
  Numeric T = t-273.15;
  //choose parametrization depending on T
  if ( T >= -58. )
  {
    la= 9.88 * exp( -0.060*T );
  }
  else
  {
    la= 0.75 * exp( -0.1057*T );
  }
  if ( T >= -61. )
  {
    mu= -0.59 - 0.030*T;
  }
  else
  {
    mu= -14.09 - 0.248*T;
  }
 
  //Distribution function H13

  dNdD = pow( dmax, mu ) * exp ( -la * dmax );

  if (isnan(dNdD)) dNdD = 0.0;
  return dNdD;
}



/*! Calculates particle size and shape distribution using H13 parametrization.
 *  Each diameter of the scattering elements is a node in the distribution.
 *  One call of this function calculates number density for one scattering
 *  element.  

    \return dNdD particle number density per diameter interval [#/m3/m]
          
    \param diameter_max  maximum diameter of scattering scattering element [m]
    \param t     atmospheric temperature [K]
  
  \author Johan Strandgren  
  \date 2013-08-26

*/
Numeric IWCtopnd_H13Shape ( const Numeric diameter_max,
                            const Numeric t)
{  
  Numeric dNdD;
  Numeric la;
  Numeric mu;
  // convert m to cm
  

  Numeric dmax = diameter_max * 1e2;
  //convert T from Kelvin to Celsius
  Numeric T = t-273.15;
  //choose parametrization depending on T
  if ( T >= -58. )
  {
    la= 9.88 * exp( -0.060*T );
  }
  else
  {
    la= 0.75 * exp( -0.1057*T );
  }
  if ( T >= -61. )
  {
    mu= -0.59 - 0.030*T;
  }
  else
  {
    mu= -14.09 - 0.248*T;
  }
  
  //Distribution function H13Shape

  dNdD = pow( dmax, mu ) * exp ( -la * dmax );

  if (isnan(dNdD)) dNdD = 0.0;
  return dNdD;
}



/*! Calculates area ratio using H13 shape parametrization.
 *  Each scattering element is a node in the aspect ratio distribution.
 *  One call of this function calculates one aspect ratio.  

    \return dNdD particle number density per diameter interval [#/m3/m]
          
    \param diameter_max  maximum diameter of scattering scattering element [m]
    \param t     atmospheric temperature [K]
  
  \author Johan Strandgren  
  \date 2013-08-26

*/
Numeric area_ratioH13 ( const Numeric diameter_max,
                        const Numeric t)
{  
  Numeric Ar;
  Numeric alpha;
  Numeric beta;

  // convert m to cm
  Numeric dmax = diameter_max * 1e2;
  //convert T from Kelvin to Celsius
  Numeric T = t-273.15;
  //Parameterize for all temperatures
  
  alpha = 0.25*exp(0.0161*T);
  
  beta = -0.25+0.0033*T;
  
  // Area ratio function depending on temperature

  Ar = alpha*pow(dmax,beta);

  if (isnan(Ar)) Ar = 0.0;
  return Ar;
}



/*! Calculates particle size distribution for liquid water clouds using a gamma
 *  parametrization by Hess et al., 1998 (continental stratus).
 *  Each scattering element is a node in the distribution.
 *  One call of this function calculates one particle number density.
 *  Implicitly assumes particles of liquid water.

	\return n  particle number density per radius interval [#/m3*m]
         
	\param lwc atmospheric liquid water content [kg/m3]
	\param r radius of scattering scattering element [m]
  
  \author Daniel Kreyling
  \date 2010-12-16

*/
Numeric LWCtopnd (const Numeric lwc, //[kg/m^3]
                  const Numeric radius // [m]
		  )
{ 
  // skip calculation if LWC is 0.0
  if ( lwc == 0.0 )
  {
    return 0.0;
  }
	Numeric rc = 4.7*1e-6; //[um]->[m]
	Numeric alpha = 5.0;
	Numeric gam = 1.05;
	
	Numeric a4g = (alpha+4.)/gam;
	Numeric B = (alpha/gam) / pow(rc,gam); 
	Numeric A = 0.75/PI * lwc/DENSITY_OF_WATER * gam * pow(B,a4g) /
                    gamma_func(a4g);
	Numeric dNdr = A * (pow(radius,alpha) * exp(-B*pow(radius,gam))); // [#/m3/m]
	
/* alternative implementation
  Numeric rc = 4.7; //micron
	Numeric alpha = 5.0;
	Numeric gam = 1.05;
	
	Numeric B=(alpha/gam)/pow(rc,gam); 
	Numeric A=gam*pow(B,((alpha+1)/gam))/gamma_func((alpha+1)/gam);
	Numeric dNdr=A*(pow(radius*1e6,alpha)*exp(-B*pow(radius*1e6,gam)))*1e6; // [#/m3/m]
*/

  if (isnan(dNdr)) dNdr = 0.0;
	return dNdr;
}

/*! Calculates the particle number density field according
 *  to the two moment scheme of Axel Seifert, that is used the ICON model.
 *  One call of this function calculates one particle number density.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param mass   Mass of scattering particle [kg]
 \param N_tot  Total number of particles (0th moment) [#/m3/m/kg^mu]
 \param M      Total mass concentration of Particles (1st moment) [kg/m^3]
 \param psd_type string with a tag defining the (hydrometeor) scheme

 
 \author Manfred Brath
 \date 2015-01-19
 
 */
Numeric psd_S2M (const Numeric mass,
                     const Numeric N_tot,
                     const Numeric M1,
                     const String psd_type)
{
    Numeric dN;
    Numeric N0;
    Numeric Lambda;
    Numeric arg1;
    Numeric arg2;
    Numeric temp;
    Numeric mu;
    Numeric gamma;
    Numeric xmin;
    Numeric xmax;
    Numeric M0min;
    Numeric M0max;
    Numeric M0;
    
    
    
    // Get the coefficients for the right hydrometeor
    if ( psd_type == "S2M_IWC" ) //Cloud ice water
    {
        mu=0.;
        gamma=1./3.;
        xmin=1e-12;
        xmax=1e-5;
    }
    else if ( psd_type == "S2M_RWC" ) //Rain
    {
        mu=0.;
        gamma=1./3.;
        xmin=2.6e-10;
        xmax=3e-6;
    }
    else if ( psd_type == "S2M_SWC" ) //Snow
    {
        mu=0.;
        gamma=1./2.;
        xmin=1e-10;
        xmax=2e-5;
    }
    else if ( psd_type == "S2M_GWC" ) //Graupel
    {
        mu=1.;
        gamma=1./3.;
        xmin=1e-9;
        xmax=5e-4;
    }
    else if ( psd_type == "S2M_HWC" ) //Hail
    {
        mu=1.;
        gamma=1./3.;
        xmin=2.6e-10;
        xmax=5e-4;
    }
    else if ( psd_type == "S2M_LWC" ) //Cloud liquid water
    {
        mu=1;
        gamma=1;
        xmin=4.2e-15;
        xmax=2.6e-10;
    }
    else
    {
        ostringstream os;
        os << "You use a wrong tag! ";
        throw runtime_error( os.str() );
    }

    M0=N_tot;
    
    // lower and upper limit check is taken from the ICON code of the two moment
    //scheme
    
    M0max=M1/xmax;
    M0min=M1/xmin;
    
    //check lower limit of the scheme
    if (M0>M0min)
    {
        M0=M0min;
    }
    
    //check upper limit of the scheme
    if (M0<M0max)
    {
        M0=M0max;
    }
    
    
    //arguments for Gamma function
    arg2=(mu+2)/gamma;
    arg1=(mu+1)/gamma;
    
    
    temp=M0/M1*gamma_func(arg2)/gamma_func(arg1);
    
    //Lambda (parameter for modified gamma distribution)
    Lambda=pow(temp, gamma);
    
    //N0
    N0=M0*gamma/gamma_func(arg1)*pow(Lambda, arg1);
    
    //Distribution function
    dN=mod_gamma_dist(mass, N0,Lambda, mu, gamma);
    
    if (isnan(dN)) dN = 0.0;
    
    
    return dN;
}


/*! Calculates the particle number density field according
 *  to the Milbrandt and Yau two moment scheme, which is used in the GEM model.
 *  See also milbrandt and yau, 2005.
 *  One call of this function calculates one particle number density.

 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param mass   Mass of scattering particle [kg]
 \param N_tot  Total number of particles (0th moment) [#/m3/m/kg^mu]
 \param M      Total mass concentration of Particles (1st moment) [kg/m^3]
 \param psd_type string with a tag defining the (hydrometeor) scheme
 
 
 \author Manfred Brath
 \date 2015-08-01
 
 */
Numeric psd_MY2 (const Numeric diameter_max,
                     const Numeric N_tot,
                     const Numeric M,
                     const String psd_type)
{
    Numeric dN;
    Numeric N0;
    Numeric Lambda;
    Numeric arg1;
    Numeric arg2;
    Numeric temp;
    Numeric mu;
    Numeric gamma;
    Numeric alpha;
    Numeric beta;
    
    
    // Get the coefficients for the right hydrometeor
    if ( psd_type == "MY2_IWC" ) //Cloud ice water
    {
        mu=0.;
        gamma=1.;
        alpha=440.; //[kg]
        beta=3;
    }
    else if ( psd_type == "MY2_RWC" ) //Rain
    {
        mu=0.;
        gamma=1;
        alpha=523.5988; //[kg]
        beta=3;
    }
    else if ( psd_type == "MY2_SWC" ) //Snow
    {
        mu=0.;
        gamma=1;
        alpha=52.35988; //[kg]
        beta=3;
    }
    else if ( psd_type == "MY2_GWC" ) //Graupel
    {
        mu=0.;
        gamma=1;
        alpha=209.4395; //[kg]
        beta=3;
    }
    else if ( psd_type == "MY2_HWC" ) //Hail
    {
        mu=0.;
        gamma=1;
        alpha=471.2389; //[kg]
        beta=3;
    }
    else if ( psd_type == "MY2_LWC" ) //Cloud liquid water
    {
        mu=1;
        gamma=1;
        alpha=523.5988; //[kg]
        beta=3;
    }
    else
    {
        ostringstream os;
        os << "You use a wrong tag! ";
        throw runtime_error( os.str() );
    }
    
    
    
    
    //Calculate Number density only, if mass is between xmin and xmax

    //arguments for Gamma function
    arg2=(mu+beta+1)/gamma;
    arg1=(mu+1)/gamma;
    
    
    temp=alpha*N_tot/M*gamma_func(arg2)/gamma_func(arg1);
    
    //Lambda (parameter for modified gamma distribution)
    Lambda=pow(temp, gamma/beta);
    
    //N0
    N0=N_tot*gamma/gamma_func(arg1)*pow(Lambda, arg1);
    
    //Distribution function
    dN=mod_gamma_dist(diameter_max, N0,Lambda, mu, gamma);
    
    if (isnan(dN)) dN = 0.0;
    
    
    return dN;
}



/*! Calculates the particle number density field for Cloud liquid water according
 *  the modified gamma distribution for cloud water inside Geer and Baordo (2014),
 *  see table A1
 *  One call of this function calculates one particle number density.
 *  Assumptions are: density of particles is constant and particle shape is sphere.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param d volume equivalent diameter of scattering particle [m]
 \param lwc liquid water content [kg/m^3]
 
 \author Manfred Brath
 \date 2014-11-28
 
 */
Numeric LWCtopnd_MGD_LWC ( const Numeric d, const Numeric rho, const Numeric lwc
                          )
{
    Numeric dN;
    Numeric N0;
    
    // coefficients of modified gamma distri
    const Numeric mu=2; //
    const Numeric lambda=2.13e5; //
    const Numeric gamma=1;
    const Numeric fc=1.4863e30;
    
   
    // N0
    N0=fc*lwc/rho; //[m^-4]
    
    //Distribution function
    dN=N0*pow(d ,mu)*exp(-lambda*pow(d,gamma));
    
    if (isnan(dN)) dN = 0.0;
    return dN;
}


/*! Calculates the particle number density field for Cloud ice according
 *  the modified gamma distribution for cloud ice inside Geer and Baordo (2014),
 *  see table A1
 *  One call of this function calculates one particle number density.
 *  Assumptions are: density of particles is constant and particle shape is sphere.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param d volume equivalent diameter of scattering particle [m]
 \param iwc ice water content [kg/m^3]
 
 \author Manfred Brath
 \date 2014-11-28
 
 */
Numeric IWCtopnd_MGD_IWC ( const Numeric d, const Numeric rho, const Numeric iwc)
{
    Numeric dN;
    Numeric N0;
    
    // coefficients of modified gamma distri
    const Numeric mu=2; //
    const Numeric lambda=2.05e5; //
    const Numeric gamma=1;
    const Numeric fc=1.1813e30;

    
    // N0
    N0=fc*iwc/rho; //[m^-4]
    
    //Distribution function
    dN=N0*pow(d ,mu)*exp(-lambda*pow(d,gamma));
    
    if (isnan(dN)) dN = 0.0;
    return dN;
}




/*! Calculates particle size distribution using MP48 parametrization for
 *  precipitating hydrometeors (rain, snow, ...).
 *  Each scattering element is a node in the distribution.
 *  One call of this function calculates one particle number density.  

    \return dNdD particle number density per diameter interval [#/m3*m]
          
    \param PR    precipitation rate [mm/hr]
    \param diameter_melted_equivalent melted equivalent diameter of scattering scattering element [m]
 
  \author Jana Mendrok
  \date 2012-04-04

*/
Numeric PRtopnd_MP48 (const Numeric PR,
                      const Numeric diameter_melted_equivalent)
{
  // skip calculation if PR is 0.0
  if ( PR == 0.0 )
  {
    return 0.0;
  }

  Numeric N0 = 0.08*1e8; // [#/cm3/cm] converted to [#/m3/m]
  Numeric lambda = 41.*1e2*pow(PR,-0.21); // [1/cm] converted to [1/m] to fit
                                          // diameter_melted_equivalent [m]

  Numeric n = N0*exp(-lambda*diameter_melted_equivalent);
  return n;
}




