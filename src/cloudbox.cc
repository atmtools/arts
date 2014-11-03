/* Copyright (C) 2002-2012 Claudia Emde <claudia.emde@dlr.de>
                      
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
  \file   cloudbox.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Thu May  23 10:59:55 2002
  
  \brief  Internal functions for scattering calculations.
*/

#include "cloudbox.h"

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Numeric PI;

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <limits>

#include "arts.h"
#include "messages.h"
#include "make_vector.h"
#include "logic.h"
#include "ppath.h"
#include "physics_funcs.h"
#include "math_funcs.h"
#include "check_input.h"
#include "rng.h"
#include <ctime>
#include "mc_antenna.h"
#include "sorting.h"



//! Check whether hydromet grid size is equal to atmospheric grid size 
//! and if hydromet profile is zero (no cloud) in each grid point.
/*!
  \param dim atmosphere dimension 
  \param hydromet hydrometeor grid in p, lat or lon dimension
  \param p_grid p grid of current atmosphere
  \param lat_grid lat grid of current atmosphere
  \param lon_grid lon grid of current atmosphere

  \author Daniel Kreyling
  \date   2011-01-27
*/  
void chk_massdensity_field( bool& empty_flag,
			 const Index&  dim,	
			 const Tensor3& massdensity, 
			 const Vector& p_grid,
			 const Vector& lat_grid,
			 const Vector& lon_grid
		       )
{
  
  // check p
  if ( massdensity.npages() != p_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *p_grid* ("
         << p_grid.nelem()
         << ") is not equal to the number of pages of *massdensity* ("
         << massdensity.npages()
         <<").";
      throw runtime_error(os.str() );
  }
  
  // check lat
  if(dim >= 2 )
  {
    if ( massdensity.nrows() != lat_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *lat_grid* ("
         << lat_grid.nelem()
         << ") is not equal to the number of rows of *massdensity* ("
         << massdensity
         << ").";
      throw runtime_error(os.str() );
      
    }
  }
  
  // check lon
  if(dim == 3 )
  {
    if ( massdensity.ncols() != lon_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *lon_grid* ("
         << lon_grid.nelem()
         << ") is not equal to the number of columns of *massdensity*"
         << massdensity
         << ").";
      throw runtime_error(os.str() );
      
    }
  }
  
  empty_flag = false;
  // set empty_flag to true if a single value of hydromet_field is unequal zero    
    for (Index j=0; j<massdensity.npages(); j++) {
      for (Index k=0; k<massdensity.nrows(); k++) {
	    for (Index l=0; l<massdensity.ncols(); l++) {
	      if ( massdensity(j,k,l) != 0.0) empty_flag = true;
	}
      }
    }  
}



//! Check whether particle number density is zero at a specified pressure level
/*!
  \param i_p Pressure index
  \param pnd_field_raw Particle number density data
  \param pnd_field_file pnd field filename

  \author Claudia Emde
  \date   2005-05-09
*/ 
void chk_if_pnd_zero_p(const Index& i_p,
                       const GriddedField3& pnd_field_raw,
                       const String& pnd_field_file,
                       const Verbosity& verbosity)
  
{
  const ConstVectorView pfr_p_grid = pnd_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView pfr_lat_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView pfr_lon_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LON_GRID);

  for (Index i = 0; i < pfr_lat_grid.nelem(); i++ ) 
    {
      for (Index j = 0; j < pfr_lon_grid.nelem(); j++ )
        {
          if ( pnd_field_raw.data(i_p, i, j) != 0. )
            {
              CREATE_OUT1;
              ostringstream os;
              os << "Warning: \n"
                 << "The particle number density field contained in the file '"
                 << pnd_field_file << "'\nis non-zero outside the cloudbox "
                 << "or close the cloudbox boundary at the \n"
                 << "following position:\n"
                 << "pressure = " << pfr_p_grid[i_p] 
                 << ", p_index = " << i_p << "\n"
                 << "latitude = " << pfr_lat_grid[i] 
                 << ", lat_index = " << i << "\n"
                 << "longitude = " << pfr_lon_grid[j] 
                 << ", lon_index = " << j << "\n"
                 << "\n";
              // throw runtime_error( os.str() );
              out1 << os.str();
            }
        }
    }
}



//! Check whether particle number density is zero at a specified latitude
/*!
  \param i_lat          Latitude index
  \param pnd_field_raw  Particle number density data
  \param pnd_field_file pnd field filename

  \author Claudia Emde
  \date   2005-05-09
*/           
void chk_if_pnd_zero_lat(const Index& i_lat,
                         const GriddedField3& pnd_field_raw,
                         const String& pnd_field_file,
                         const Verbosity& verbosity)
  
{
  const ConstVectorView pfr_p_grid = pnd_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView pfr_lat_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView pfr_lon_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LON_GRID);

  for (Index i = 0; i < pfr_p_grid.nelem(); i++ ) 
    {
      for (Index j = 0; j < pfr_lon_grid.nelem(); j++ )
        {
          if ( pnd_field_raw.data(i, i_lat, j) != 0. )
            {
              CREATE_OUT1;
              ostringstream os;
              os << "Warning: \n" 
                 << "The particle number density field contained in the file '"
                 << pnd_field_file << "'\nis non-zero outside the cloudbox "
                 << "or close the cloudbox boundary at the \n"
                 << "following position:\n"
                 << "pressure = " << pfr_p_grid[i] << ", p_index = "
                 << i << "\n"
                 << "latitude = " << pfr_lat_grid[i_lat] 
                 << ", lat_index = "<<i_lat<< "\n"
                 << "longitude = " << pfr_lon_grid[j] 
                 << ", lon_index = " << j << "\n"
                 << "\n";
                //throw runtime_error( os.str() );
              out1 << os.str();
            }
        }
    }
}



//! Check whether particle number density is zero at a specified longitude
/*!
  \param i_lon          Latitude index
  \param pnd_field_raw  Particle number density data
  \param pnd_field_file pnd field filename

  \author Claudia Emde
  \date   2005-05-09
*/           
void chk_if_pnd_zero_lon(const Index& i_lon,
                         const GriddedField3& pnd_field_raw,
                         const String& pnd_field_file,
                         const Verbosity& verbosity)
  
{
  const ConstVectorView pfr_p_grid = pnd_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView pfr_lat_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView pfr_lon_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LON_GRID);

  for (Index i = 0; i < pfr_p_grid.nelem(); i++ ) 
    {
      for (Index j = 0; j < pfr_lat_grid.nelem(); j++ )
        {
          if ( pnd_field_raw.data(i, j, i_lon) != 0. )
            {
              CREATE_OUT1;
              ostringstream os;
              os << "Warning: \n" 
                 << "The particle number density field contained in the file '"
                 << pnd_field_file << "'\nis non-zero outside the cloudbox "
                 << "or close the cloudbox boundary at the \n"
                 << "following position:\n"
                 << "pressure = " << pfr_p_grid[i] 
                 << ", p_index = " << i << "\n"
                 << "latitude = " << pfr_lat_grid[j] 
                 << ", lat_index = " << j << "\n"
                 << "longitude = " << pfr_lon_grid[i_lon] 
                 << ", lon_index = "
                 << i_lon << "\n"
                 << "\n";
              //throw runtime_error( os.str() );
              out1 << os.str();
            }
        }
    }
}



//! Check particle number density files
/*!
  This function checks, whether the particle number density file
  has the right atmospheric dimension (check for no non-zero pnd values outside
  cloudbox removed. done in pnd_fieldCalc). 

  \param pnd_field_raw   pnd field data
  \param pnd_field_file  pnd field filename
  \param atmosphere_dim  Atmospheric dimension

  \author Claudia Emde
  \date   2005-04-05
*/ 
void chk_pnd_data(
                  const GriddedField3& pnd_field_raw,
                  const String& pnd_field_file,
                  const Index& atmosphere_dim,
                  const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  const ConstVectorView pfr_p_grid = pnd_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView pfr_lat_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView pfr_lon_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LON_GRID);

  // The consistency of the dimensions is checked in the reading routine. 
  // Here we have to check whether the atmospheric dimension is correct and whether 
  // the particle number density is 0 on the cloudbox boundary and outside the cloudbox.
  
  out3 << "Check particle number density file " << pnd_field_file 
       << "\n"; 
 
  if (atmosphere_dim == 1 && (pfr_lat_grid.nelem() != 1 
                              || pfr_lon_grid.nelem() != 1) )
    {
      ostringstream os; 
      os << "The atmospheric dimension is 1D but the particle "
         << "number density file * " << pnd_field_file 
         << " is for a 3D atmosphere. \n";
      throw runtime_error( os.str() );
    }
      
  
  else if( atmosphere_dim == 3) 
    {
      if(pfr_lat_grid.nelem() == 1 
         || pfr_lon_grid.nelem() == 1)
        {
          ostringstream os; 
          os << "The atmospheric dimension is 3D but the particle "
             << "number density file * " << pnd_field_file 
             << " is for a 1D or a 2D atmosphere. \n";
          throw runtime_error( os.str() );
        }
    } 
  
  out3 << "Particle number density data is o.k. \n";  
}



//! Check particle number density files (pnd_field_raw)
/*!
  
  \param pnd_field_raw   pnd field raw data (array for all particle types)
  \param pnd_field_file  pnd field filename
  \param atmosphere_dim  Atmospheric dimension
 
  \author Claudia Emde
  \date   2005-04-05
*/ 
void chk_pnd_raw_data(
                      const ArrayOfGriddedField3& pnd_field_raw,
                      const String& pnd_field_file,
                      const Index& atmosphere_dim,
                      const Verbosity& verbosity
                      )
{
  CREATE_OUT3;
  
  for( Index i = 0; i < pnd_field_raw.nelem(); i++)
    {
      out3 << "Element in pnd_field_raw_file:" << i << "\n";
      chk_pnd_data(pnd_field_raw[i],
                   pnd_field_file, atmosphere_dim,
                   verbosity);
    }
}



//! chk_pnd_field_raw_only_in_cloudbox
/*! 
    Checks whether the pnd_field is zero outside the cloudbox.
    This is of a higher level than chk_pnd_data because it does
    not require any filename and because it works on all pnd_field_raw
    rather than just one element. Otherwise, it is mostly a new
    implementation of the same functionality.

    \param    dim                The atmospheric dimensionality.
    \param    pnd_field_raw      All pnd_field_raw data.
    \param    p_grid             Pressure grid.
    \param    lat_grid           Latitude grid.
    \param    lon_grid           Longitude grid.
    \param    cloudbox_limits    The edges of the cloudbox.

    \author Gerrit Holl
    \date   2011-03-24
*/
void chk_pnd_field_raw_only_in_cloudbox(
        const Index&                 dim,
        const ArrayOfGriddedField3&  pnd_field_raw,  
        ConstVectorView              p_grid,
        ConstVectorView              lat_grid,
        ConstVectorView              lon_grid,
        const ArrayOfIndex&          cloudbox_limits )
{
    Numeric p, lat, lon, v;
    Index n, p_i, lat_i, lon_i;
    // For any non-zero point, verify we're outside the cloudbox
    for (n=0; n < pnd_field_raw.nelem(); n++) {
        for (p_i=0; p_i < pnd_field_raw[n].data.npages(); p_i++) {
            for (lat_i=0; lat_i < pnd_field_raw[n].data.nrows(); lat_i++) {
                for (lon_i=0; lon_i < pnd_field_raw[n].data.ncols(); lon_i++) {
                    v = pnd_field_raw[n].data(p_i, lat_i, lon_i);
                    if (v != 0) {
                        // Verify pressure is between cloudbox limits
                        p = pnd_field_raw[n].get_numeric_grid(GFIELD3_P_GRID)[p_i];
//                        if (!((p <= p_grid[cloudbox_limits[0]]) &
//                              (p >= p_grid[cloudbox_limits[1]]))) {
                        if ( (p <= p_grid[cloudbox_limits[1]]) ||
                             ( (p >= p_grid[cloudbox_limits[0]]) &&
                              (cloudbox_limits[0]!=0)) ) {
                            ostringstream os;
                            os << "Found non-zero pnd outside cloudbox. "
                               << "Cloudbox extends from p="
                               << p_grid[cloudbox_limits[0]]
                               << " Pa to p="
                               << p_grid[cloudbox_limits[1]]
                               << " Pa, but found pnd=" << v
                               << "/m³ at p=" << p << " Pa for particle #" << n
                               << ".";
                            throw runtime_error(os.str());
                        }
                        // Verify latitude is too
                        if (dim > 1) {
                            lat = pnd_field_raw[n].get_numeric_grid(GFIELD3_LAT_GRID)[lat_i];
                            if (!((lat > lat_grid[cloudbox_limits[2]]) &
                                  (lat < lat_grid[cloudbox_limits[3]]))) {
                                ostringstream os;
                                os << "Found non-zero pnd outside cloudbox. "
                                   << "Cloudbox extends from lat="
                                   << lat_grid[cloudbox_limits[2]]
                                   << "° to lat="
                                   << lat_grid[cloudbox_limits[3]]
                                   << "°, but found pnd=" << v
                                   << "/m³ at lat=" << lat << "° for particle #" << n
                                   << ".";
                                throw runtime_error(os.str());
                            }
                        }
                        // Etc. for longitude
                        if (dim > 2) {
                            lon = pnd_field_raw[n].get_numeric_grid(GFIELD3_LON_GRID)[lon_i];
                            if (!((lon > lon_grid[cloudbox_limits[4]]) &
                                  (lon < lon_grid[cloudbox_limits[5]]))) {
                                ostringstream os;
                                os << "Found non-zero pnd outside cloudbox. "
                                   << "Cloudbox extends from lon="
                                   << lon_grid[cloudbox_limits[4]]
                                   << "° to lat="
                                   << lon_grid[cloudbox_limits[5]]
                                   << "°, but found pnd=" << v
                                   << "/m³ at lon=" << lon << "° for particle #" << n
                                   << ".";
                                throw runtime_error(os.str());
                            }
                        }
                    }
                }
            }
        }
    }
}



//!  Check validity of part_species setting
/*!
  This function checks, whether number of elements in each particle string is
  ok, and whether the entries for size limits are indeed numbers (or '*').

	\param part_species Array of particle species tags.
  \param delim Delimiter string of *part_species* elements.

  \author Jana Mendrok
  \date 2012-10-25

*/
void chk_part_species (
                      const ArrayOfString& part_species,
                      const String& delim
                      )
{
  ArrayOfString strarr;
  Numeric sizeval;

  for ( Index k=0; k<part_species.nelem(); k++ )
    {
      part_species[k].split ( strarr, delim );
      if ( strarr.nelem() > 5 )
        {     
          ostringstream os;
          os << "Individual strings in part_species can only contain up to 5\n"
             << "elements, but entry #" << k << " contains the following "
             << strarr.nelem() << ":\n" << strarr << "\n";
          throw runtime_error ( os.str() );
        }

      if ( strarr.nelem() > 3 && strarr[3] != "*" )
        {
          istringstream os1 ( strarr[3] );
          os1 >> sizeval;
          if (os1.fail())
            {
              ostringstream os;
              os << "Sizemin specification can only be a number or wildcard ('*')"
                 << ", but is '" << strarr[3] << "'\n";
              throw runtime_error ( os.str() );
            }
        }

      if ( strarr.nelem() > 4 && strarr[4] != "*" )
        {
          istringstream os1 ( strarr[4] );
          os1 >> sizeval;
          if (os1.fail())
            {
              ostringstream os;
              os << "Sizemax specification can only be a number or wildcard ('*')"
                 << ", but is '" << strarr[4] << "'\n";
              throw runtime_error ( os.str() );
            }
        }
    }
}



//! Check scattering data general
/*!
  FIXME
  
  \param scat_data_array Array of single scattering data
  \param scat_meta_array Array of scattering meta data

  \author Daniel Kreyling
  \date 2010-12-02
*/

void chk_scattering_data(const ArrayOfSingleScatteringData& scat_data_array,
                         const ArrayOfScatteringMetaData& scat_meta_array,
                         const Verbosity&)
{
  if (scat_data_array.nelem() != scat_meta_array.nelem())
  {
    ostringstream os;
    os << "The number of elments in *scat_data_array*\n"
    << "and *scat_meta_array* do not match.\n"
    << "Each SingleScattering file must correspond\n"
    << "to one ScatteringMeta data file.";
    throw runtime_error( os.str());
  }

}

//! Check scattering data meta files
/*!
  FIXME
  
  \param scat_meta scattering meta data
  \param scat_meta_file filename of the data to be checked

  \author Daniel Kreyling
  \date 2010-12-02
*/
void chk_scattering_meta_data(const ScatteringMetaData& scat_meta _U_,
                              const String& scat_meta_file,
                              const Verbosity& verbosity)
{
  CREATE_OUT2;
  out2 << "  Check scattering meta properties file "<< scat_meta_file << "\n";
  
/* this check is outdated. type now is free from!
   however, we might want to have other things checked here!?
   - which parameters at least are needed? -> radius, ...?
   - ...
  if  (scat_meta.type != "Ice" && scat_meta.type != "Water" && scat_meta.type != "Aerosol")
  {
	  ostringstream os; 
	  os << "Type in " << scat_meta_file << " must be 'Ice', 'Water' or 'Aerosol'\n";     
	  throw runtime_error( os.str() );
	}
*/
  //(more) checks need to be included
}



//! Check single scattering data files
/*!
  This function checks, whether a datafile containing the single scattering 
  properties of a particle type includes the required frequencies and 
  temperatures and whether the angular grids are defined correctly.
  Furthermore it checks the self consistency of the data by checking the 
  dimensions of pha_mat, ext_mat and abs_vec depending on the particle_type case.
  
  \param scat_data_array Single scattering data
  \param scat_data_file Filename of the data to be checked.
  \param f_grid        Frequency grid
  
  \author Claudia Emde
  \date   2005-04-04
*/
void chk_scat_data(const SingleScatteringData& scat_data_array,
                                const String& scat_data_file,
                                ConstVectorView f_grid,
                                const Verbosity& verbosity)
{
  CREATE_OUT2;
  out2 << "  Check single scattering properties file "<< scat_data_file 
       << "\n";

  if (scat_data_array.particle_type != 10 &&
      scat_data_array.particle_type != 20 &&
      scat_data_array.particle_type != 30)
    {
      ostringstream os; 
      os << "Particle type value in file" << scat_data_file << "is wrong."
         << "It must be \n"
         << "10 - arbitrary oriented particles \n"
         << "20 - randomly oriented particles or \n"
         << "30 - horizontally aligned particles.\n";
      throw runtime_error( os.str() );
    }
    
    chk_interpolation_grids("scat_data_array.f_grid to f_grid",
			    scat_data_array.f_grid,
			    f_grid);
  
/*  if (!(scat_data_array.f_grid[0] <= f_grid[0] &&
        last(f_grid) <= 
        last(scat_data_array.f_grid) ))
    {
      ostringstream os;
      os << "The range of frequency grid in the single scattering"
         << " properties datafile " 
         << scat_data_file << " does not contain all values of"
         << "*f_grid*.";
      throw runtime_error( os.str() );
    }*/

  // Here we only check whether the temperature grid is of the unit K, not 
  // whether it corresponds to the required values it T_field. The second 
  // option is not trivial since here one has to look whether the pnd_field 
  // is none zero for the corresponding temperature. This check done in the 
  // functions where the multiplication with the particle number density is 
  // done. 

  if (!(0. < scat_data_array.T_grid[0] && last(scat_data_array.T_grid) < 1001.))
    {
      ostringstream os;
      os << "The temperature values in " <<  scat_data_file 
         << " are negative or very large. Check whether you have used the "
         << "right unit [Kelvin].";
      throw runtime_error( os.str() );
    }
  
  if (scat_data_array.za_grid[0] != 0.)
    {
      ostringstream os;
      os << "The first value of the za grid in the single" 
         << " scattering properties datafile " 
         << scat_data_file << " must be 0.";
        throw runtime_error( os.str() );
    } 

  if (last(scat_data_array.za_grid) != 180.)
    {
      ostringstream os;
      os << "The last value of the za grid in the single"
         << " scattering properties datafile " 
         << scat_data_file << " must be 180.";
      throw runtime_error( os.str() );
    } 
  
  if (scat_data_array.particle_type == 10 && scat_data_array.aa_grid[0] != -180.)
     {
       ostringstream os;
       os << "For particle_type = 10 (general orientation) the first value"
          << " of the aa grid in the single scattering"
          << " properties datafile " 
          << scat_data_file << "must be -180.";
         throw runtime_error( os.str() );
     } 
  
  if (scat_data_array.particle_type == 30 && scat_data_array.aa_grid[0] != 0.)
    {
      ostringstream os;
      os << "For particle_type = 30 (horizontal orientation)"
         << " the first value"
         << " of the aa grid in the single scattering"
         << " properties datafile " 
         << scat_data_file << "must be 0.";
        throw runtime_error( os.str() );
    }   
  
  if (scat_data_array.particle_type == 30 && last(scat_data_array.aa_grid) != 180.)
    {
      ostringstream os;
      os << "The last value of the aa grid in the single"
         << " scattering properties datafile " 
         << scat_data_file << " must be 180.";
        throw runtime_error( os.str() );
    }   

  ostringstream os_pha_mat;
  os_pha_mat << "pha_mat* in the file *" << scat_data_file;
  ostringstream os_ext_mat;
  os_ext_mat << "ext_mat* in the file * " << scat_data_file;
  ostringstream os_abs_vec;
  os_abs_vec << "abs_vec* in the file * " << scat_data_file;
  
  switch (scat_data_array.particle_type){
    
  case PARTICLE_TYPE_GENERAL:
    
    out2 << "  Datafile is for arbitrarily orientated particles. \n";
    
    chk_size(os_pha_mat.str(), scat_data_array.pha_mat_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             scat_data_array.za_grid.nelem(), scat_data_array.aa_grid.nelem(),
             scat_data_array.za_grid.nelem(), scat_data_array.aa_grid.nelem(), 
              16); 
    
    chk_size(os_ext_mat.str(), scat_data_array.ext_mat_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             scat_data_array.za_grid.nelem(), scat_data_array.aa_grid.nelem(),
             7);
    
    chk_size(os_abs_vec.str(), scat_data_array.abs_vec_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             scat_data_array.za_grid.nelem(), scat_data_array.aa_grid.nelem(),
             4);
    break;
    
  case PARTICLE_TYPE_MACROS_ISO: 
    
    out2 << "  Datafile is for randomly oriented particles, i.e., "
         << "macroscopically isotropic and mirror-symmetric scattering "
         << "media. \n";
    
    chk_size(os_pha_mat.str(), scat_data_array.pha_mat_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             scat_data_array.za_grid.nelem(), 1, 1, 1, 6);
    
    chk_size(os_ext_mat.str(), scat_data_array.ext_mat_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             1, 1, 1);
    
    chk_size(os_abs_vec.str(), scat_data_array.abs_vec_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             1, 1, 1);
    break; 
    
  case PARTICLE_TYPE_HORIZ_AL:
    
    out2 << "  Datafile is for horizontally aligned particles. \n"; 
    
    chk_size(os_pha_mat.str(), scat_data_array.pha_mat_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             scat_data_array.za_grid.nelem(), scat_data_array.aa_grid.nelem(),
             scat_data_array.za_grid.nelem()/2+1, 1, 
             16); 

    chk_size(os_ext_mat.str(), scat_data_array.ext_mat_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             scat_data_array.za_grid.nelem()/2+1, 1, 
             3);
    
    chk_size(os_abs_vec.str(), scat_data_array.abs_vec_data,
             scat_data_array.f_grid.nelem(), scat_data_array.T_grid.nelem(),
             scat_data_array.za_grid.nelem()/2+1, 1, 
             2);
    break;

  case PARTICLE_TYPE_SPHERICAL:
    throw runtime_error(
                        "Special case for spherical particles not"
                        "implemented."
                        "Use p20 (randomly oriented particles) instead."
                        );
    
  }

}





/*! Checks, whether a gridpoint is inside the cloudbox.

    \return true is returned if the point is inside the 
          cloudbox.
          
  \param gp_p  pressure GridPos
  \param gp_lat latitude GridPos
  \param gp_lon longitude GridPos
  \param cloudbox_limits The limits of the cloudbox.
  \param include_boundaries boolean: determines whther or not points on the 
  boundary are considered to be inside the cloudbox.

  \author Claudia Emde (rewritten by Cory Davis 2005-07-03)
  \date 2003-06-06

*/
bool is_gp_inside_cloudbox(
   const GridPos&      gp_p,
   const GridPos&      gp_lat,
   const GridPos&      gp_lon,
   const ArrayOfIndex& cloudbox_limits,
   const bool&         include_boundaries,
   const Index&        atmosphere_dim )
                        
{
  if (include_boundaries)
    {
      // Pressure dimension
      double ipos = fractional_gp( gp_p );
      if( ipos < double( cloudbox_limits[0] )  ||
          ipos > double( cloudbox_limits[1] ) )
        { return false; }
    
      else if( atmosphere_dim >= 2 )
        {
          // Latitude dimension
          ipos = fractional_gp( gp_lat );
          if( ipos < double( cloudbox_limits[2] )  || 
              ipos > double( cloudbox_limits[3] ) )
            { return false; }
      
          else if( atmosphere_dim == 3 )
            {
              // Longitude dimension
              ipos = fractional_gp( gp_lon );
              if( ipos < double( cloudbox_limits[4] )  || 
                  ipos > double( cloudbox_limits[5] ) )
                { return false; } 
            } 
        }
      return true;
    }
  else
    {
      // Pressure dimension
      double ipos = fractional_gp( gp_p );
      if( ipos <= double( cloudbox_limits[0] )  ||
          ipos >= double( cloudbox_limits[1] ) )
        { return false; }
      
      else if( atmosphere_dim >= 2 )
        {
          // Latitude dimension
          ipos = fractional_gp( gp_lat );
          if( ipos <= double( cloudbox_limits[2] )  || 
              ipos >= double( cloudbox_limits[3] ) )
            { return false; }
        
          else if( atmosphere_dim == 3 )
            {
              // Longitude dimension
              ipos = fractional_gp( gp_lon );
              if( ipos <= double( cloudbox_limits[4] )  || 
                  ipos >= double( cloudbox_limits[5] ) )
                { return false; }           
            }
        }
      return true;
    }
}



/*! Checks, whether the last point of a propagation path 
  is inside the cloudbox.

  Works only for 3D !!!

    \return true is returned if the point is inside the 
          cloudbox.
          
  \param ppath_step Propagation path step.
  \param cloudbox_limits The limits of the cloudbox.
  \param include_boundaries boolean: determines whther or not points on the 
  boundary are considered to be inside the cloudbox.

  \author Claudia Emde (rewritten by Cory Davis 2005-07-03)
  \date 2003-06-06

*/
bool is_inside_cloudbox(const Ppath& ppath_step,
                        const ArrayOfIndex& cloudbox_limits,
                        const bool include_boundaries)
                        
{
  assert( cloudbox_limits.nelem() == 6 );
  const Index np=ppath_step.np;
  
  return is_gp_inside_cloudbox(ppath_step.gp_p[np-1],ppath_step.gp_lat[np-1],
                               ppath_step.gp_lon[np-1],cloudbox_limits,include_boundaries);
  
}



/*! Calculates the particle number density field for McFarquhar and Heymsfield
    (1997) size distribution. To be used for cloud ice.

    \return pnd_field Particle number density field
    \param IWC_field mass content (cloud ice) field [kg/m3]
    \param t_field atmospheric temperature [K]
    \param limits pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta_array particle meta data for particles
    \param scat_data_start start index for particles handled by this distribution
    \param npart number of particles handled by this distribution
    \param part_string part_species tag for profile/distribution handled here
    \param delim Delimiter string of *part_species* elements
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2012-04-03

*/
void pnd_fieldMH97 (Tensor4View pnd_field,
                    const Tensor3& IWC_field,
                    const Tensor3& t_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfScatteringMetaData& scat_meta_array,
                    const Index& scat_data_start,
                    const Index& npart,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity)
{
  ArrayOfIndex intarr;
  Index pos;
  Vector vol_unsorted ( npart, 0.0 );
  Vector vol ( npart, 0.0 );
  Vector dm ( npart, 0.0 );
  Vector rho ( npart, 0.0 );
  Vector pnd ( npart, 0.0 );
  Vector dN ( npart, 0.0 );

  String psd_param;
  String partfield_name;

  //split String and copy to ArrayOfString
  parse_psd_param( psd_param, part_string, delim);
  parse_partfield_name( partfield_name, part_string, delim);

  bool noisy = (psd_param == "MH97n");

  for ( Index i=0; i < npart; i++ )
      vol_unsorted[i] = ( scat_meta_array[i+scat_data_start].volume );
  get_sorted_indexes(intarr, vol_unsorted);

  // extract scattering meta data
  for ( Index i=0; i< npart; i++ )
  {
      pos = intarr[i]+scat_data_start;

      vol[i] = ( scat_meta_array[pos].volume ); //m^3
      // calculate melted diameter from volume [m]
      dm[i] = pow ( 
             ( scat_meta_array[pos].volume *6./PI ),
             ( 1./3. ) );
      // get density from meta data [kg/m^3]
      rho[i] = scat_meta_array[pos].density;
      // get aspect ratio from meta data [ ]
  }
  
  if (dm.nelem() > 0)
  // dm.nelem()=0 implies no selected particles for the respective particle
  // field. should not occur.
  {
      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // we only need to go through PSD calc if there is any material
            if (IWC_field ( p, lat, lon ) > 0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<dm.nelem(); i++ )
                {
                  // calculate particle size distribution with MH97
                  // [# m^-3 m^-1]
                  dN[i] = IWCtopnd_MH97 ( IWC_field ( p, lat, lon ),
                                          dm[i], t_field ( p, lat, lon ),
                                          rho[i], noisy );
                  //dN2[i] = dN[i] * vol[i] * rho[i];
                }
            
                // scale pnds by bin width
                if (dm.nelem() > 1)
                  scale_pnd( pnd, dm, dN );
                else
                  pnd = dN;
	    
                // calculate error of pnd sum and real XWC
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), vol, rho,
                             p, lat, lon, partfield_name, verbosity );
	    
                // writing pnd vector to wsv pnd_field
                for ( Index i = 0; i< npart; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = pnd[i];
                }
            }
            else
                for ( Index i = 0; i< npart; i++ )
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

    \param pnd_field Particle number density field
    \param IWC_field mass content (cloud ice or snow) field [kg/m3]
    \param t_field atmospheric temperature [K]
    \\param limits pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta_array particle meta data for particles
    \param scat_data_start start index for particles handled by this distribution
    \param npart number of particles handled by this distribution
    \param part_string part_species tag for profile/distribution handled here
    \param delim Delimiter string of *part_species* elements
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2012-04-05

*/
void pnd_fieldH11 (Tensor4View pnd_field,
                   const Tensor3& IWC_field,
                   const Tensor3& t_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfScatteringMetaData& scat_meta_array,
                   const Index& scat_data_start,
                   const Index& npart,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
  ArrayOfIndex intarr;
  Index pos;
  Vector dmax_unsorted ( npart, 0.0 );
  Vector vol ( npart, 0.0 );
  Vector dm ( npart, 0.0 );
  Vector rho ( npart, 0.0 );
  Vector pnd ( npart, 0.0 );
  Vector dN ( npart, 0.0 );
  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < npart; i++ )
      dmax_unsorted[i] = ( scat_meta_array[i+scat_data_start].diameter_max );
  get_sorted_indexes(intarr, dmax_unsorted);
      
  // extract scattering meta data
  for ( Index i=0; i< npart; i++ )
  {
      pos = intarr[i]+scat_data_start;

      vol[i]= scat_meta_array[pos].volume; //[m^3]
      // get maximum diameter from meta data [m]
      dm[i] = scat_meta_array[pos].diameter_max;
      // get density from meta data [kg/m^3]
      rho[i] = scat_meta_array[pos].density;
      // get aspect ratio from meta data [ ]
  }

  if (dm.nelem() > 0)
  // dm.nelem()=0 implies no selected particles for the respective particle
  // field. should not occur anymore.
  {
      // itertation over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // we only need to go through PSD calc if there is any material
            if (IWC_field ( p, lat, lon ) > 0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<dm.nelem(); i++ ) //loop over number of particles
                {
                    // calculate particle size distribution for H11
                    // [# m^-3 m^-1]
                    dN[i] = IWCtopnd_H11 ( dm[i], t_field ( p, lat, lon ) );
                }
                // scale pnds by scale width
                if (dm.nelem() > 1)
                    scale_pnd( pnd, dm, dN ); //[# m^-3]
                else
                    pnd = dN;

                // scale H11 distribution (which is independent of Ice or 
                // Snow massdensity) to current massdensity.
                /* JM120411: we don't need this - it's doing exactly, what
                             chk_pndsum is doing. instead we redefine verbosity
                             for checksum here to suppress the 'scaling is off'
                             warnings and let chk_pndsum do the proper scaling.
                scale_H11 ( pnd, IWC_field ( p,lat,lon ), vol, rho ); */

                // calculate proper scaling of pnd sum from real IWC and apply
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), vol, rho,
                             p, lat, lon, partfield_name, verbosity );

                // writing pnd vector to wsv pnd_field
                for ( Index i =0; i< npart; i++ )
                {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                }
            }
            else
            {
                for ( Index i = 0; i< npart; i++ )
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

    \param pnd_field Particle number density field
    \param IWC_field mass content (cloud ice or snow) field [kg/m3]
    \param t_field atmospheric temperature [K]
    \\param limits pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta_array particle meta data for particles
    \param scat_data_start start index for particles handled by this distribution
    \param npart number of particles handled by this distribution
    \param part_string part_species tag for profile/distribution handled here
    \param delim Delimiter string of *part_species* elements
  
  \author Johan Strandgren, Daniel Kreyling
  \date 2013-08-26

*/
void pnd_fieldH13 (Tensor4View pnd_field,
                   const Tensor3& IWC_field,
                   const Tensor3& t_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfScatteringMetaData& scat_meta_array,
                   const Index& scat_data_start,
                   const Index& npart,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
  ArrayOfIndex intarr;
  Index pos;
  Vector dmax_unsorted ( npart, 0.0 );
  Vector vol ( npart, 0.0 );
  Vector dm ( npart, 0.0 );
  Vector rho ( npart, 0.0 );
  Vector pnd ( npart, 0.0 );
  Vector dN ( npart, 0.0 );
//  Vector ar ( npart, 0.0 );
  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < npart; i++ )
      dmax_unsorted[i] = ( scat_meta_array[i+scat_data_start].diameter_max );
  get_sorted_indexes(intarr, dmax_unsorted);
      
  // extract scattering meta data
  for ( Index i=0; i< npart; i++ )
  {
      pos = intarr[i]+scat_data_start;

      vol[i]= scat_meta_array[pos].volume; //[m^3]
      // get maximum diameter from meta data [m]
      dm[i] = scat_meta_array[pos].diameter_max;
      // get density from meta data [kg/m^3]
      rho[i] = scat_meta_array[pos].density;
      // get aspect ratio from meta data [ ]
//      ar[i] = scat_meta_array[pos].aspect_ratio;
  }

/*
    // Collect all unique aspect ratios and check if the are more than one
    vector<Numeric> ar_in;
    for (Iterator1D it = ar.begin(); it != ar.end(); ++it)
        if (find(ar_in.begin(), ar_in.end(), *it) == ar_in.end())
            ar_in.push_back(*it);
    
    if (ar_in.size()>1)
    {    
        ostringstream os;
        os << "There are " << ar_in.size() << " unique aspect ratios in *scat_meta_array*.\n"
        "This parametrization is only valid for one single\n"
        "aspect ratio\n";
        throw runtime_error(os.str());
    }
*/
    
  if (dm.nelem() > 0)
  // dm.nelem()=0 implies no selected particles for the respective particle
  // field. should not occur anymore.
  {
      // itertation over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // we only need to go through PSD calc if there is any material
            if (IWC_field ( p, lat, lon ) > 0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<dm.nelem(); i++ ) //loop over number of particles
                {
                    // calculate particle size distribution for H13
                    // [# m^-3 m^-1]
                    dN[i] = IWCtopnd_H13 ( dm[i], t_field ( p, lat, lon ) );
                }
                // scale pnds by scale width
                if (dm.nelem() > 1)
                    scale_pnd( pnd, dm, dN ); //[# m^-3]
                else
                    pnd = dN;

                // scale H13 distribution (which is independent of Ice or 
                // Snow massdensity) to current massdensity.
                /* JM120411: we don't need this - it's doing exactly, what
                             chk_pndsum is doing. instead we redefine verbosity
                             for checksum here to suppress the 'scaling is off'
                             warnings and let chk_pndsum do the proper scaling.
                scale_H13 ( pnd, IWC_field ( p,lat,lon ), vol, rho );*/

                // calculate proper scaling of pnd sum from real IWC and apply
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), vol, rho,
                             p, lat, lon, partfield_name, verbosity );

                // writing pnd vector to wsv pnd_field
                for ( Index i =0; i< npart; i++ )
                {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                }
            }
            else
            {
                for ( Index i = 0; i< npart; i++ )
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

    \param pnd_field Particle number density field
    \param IWC_field mass content (cloud ice or snow) field [kg/m3]
    \param t_field atmospheric temperature [K]
    \\param limits pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta_array particle meta data for particles
    \param scat_data_start start index for particles handled by this distribution
    \param npart number of particles handled by this distribution
    \param part_string part_species tag for profile/distribution handled here
    \param delim Delimiter string of *part_species* elements
  
  \author Johan Strandgren
  \date 2013-08-26

*/
void pnd_fieldH13Shape (Tensor4View pnd_field,
                   const Tensor3& IWC_field,
                   const Tensor3& t_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfScatteringMetaData& scat_meta_array,
                   const Index& scat_data_start,
                   const Index& npart,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{ 
  CREATE_OUT1;  
    
  ArrayOfIndex intarr;
  Index pos;
  Vector dmax_unsorted ( npart, 0.0 );
  Vector vol ( npart, 0.0 );
  Vector dm ( npart, 0.0 );
  Vector rho ( npart, 0.0 );
  Vector pnd ( npart, 0.0 );
  Vector ar ( npart, 0.0 ); // Aspect ratio set for the T-matrix calculations
  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < npart; i++ )
      dmax_unsorted[i] = ( scat_meta_array[i+scat_data_start].diameter_max );
  get_sorted_indexes(intarr, dmax_unsorted);
  
  // extract scattering meta data
  for ( Index i=0; i< npart; i++ )
  {
      pos = intarr[i]+scat_data_start;

      vol[i]= scat_meta_array[pos].volume; //[m^3]
      // get maximum diameter from meta data [m]
      dm[i] = scat_meta_array[pos].diameter_max;
      // get density from meta data [kg/m^3]
      rho[i] = scat_meta_array[pos].density;
      // get aspect ratio from meta data [ ]
      ar[i] = scat_meta_array[pos].aspect_ratio;
  }
    // Collect all unique maximum diameters
    vector<Numeric> dm_in;
    for (Iterator1D it = dm.begin(); it != dm.end(); ++it)
        if (find(dm_in.begin(), dm_in.end(), *it) == dm_in.end())
            dm_in.push_back(*it);
    
    // Collect all unique aspect ratios
    vector<Numeric> ar_in;
    for (Iterator1D it = ar.begin(); it != ar.end(); ++it)
        if (find(ar_in.begin(), ar_in.end(), *it) == ar_in.end())
            ar_in.push_back(*it);
                
    Vector dm_input;
    Vector ar_input;
    dm_input=dm_in;
    ar_input=ar_in;
    
    // Check size and shape limits
    if (dm[0]<7.7*1e-5)
    {
        ostringstream os;
        os << "The *H13Shape* parametrization is only valid for particles with\n"
        "a maximum diameter >= to 77 um, the smallest value of *diameter_max* in\n"
        "this simulation is " << dm[0] << " um and thus to large. Set a new *diameter_max_grid*!\n";
        throw runtime_error(os.str());
    }

    if (ar_input.nelem()==1)
    { 
        out1 << "WARNING! Only one unique aspect ratio is used. The parametrization\n"
             << "*H13* will generate the same result but with less computations\n"
             << "and thus on a sorter time\n";
    }
    
    if (ar_input[ar_input.nelem()-1] >= 1)
    {    
        ostringstream os;
        os << "H13Shape is only valid for prolate speheroids\n"
        "and cylinders at the moment, i.e for aspect ratios smaller\n"
        "than one. The maximum aspect ratio chosen is " << ar_input[ar_input.nelem()-1] << ".\n"
        "Set a new *aspect_ratio_grid";
        throw runtime_error( os.str() );
     }
    
    Vector dN ( dm_input.nelem(), 0.0 );
    Vector Ar ( dm_input.nelem(), 0.0 );
    Vector arH13 ( dm_input.nelem(), 0.0 );
    Vector pnd_temp ( dm_input.nelem(), 0.0 );

    
  //const bool suppress=true;
  //const Verbosity temp_verb(0,0,0);

  if (dm_input.nelem() > 0)
  // dm.nelem()=0 implies no selected particles for the respective particle
  // field. should not occur anymore.
  {
      // itertation over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // we only need to go through PSD calc if there is any material
            if (IWC_field ( p, lat, lon ) > 0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<dm_input.nelem(); i++ ) //loop over number of particles
                {
                    // calculate particle size distribution for H13Shape
                    // [# m^-3 m^-1]
                    dN[i] = IWCtopnd_H13Shape ( dm_input[i], t_field ( p, lat, lon ) );
                    
                    // calculate Area ratio distribution for H13Shape
                    Ar[i] = area_ratioH13 (dm_input[i], t_field (p, lat, lon ) );
                    
                    // Aspect ratio equals area ratio (for prolate particles)
                    arH13[i] = Ar[i];
                }

                // scale pnds by scale width
                if (dm_input.nelem() > 1)
                    scale_pnd( pnd_temp, dm_input, dN ); //[# m^-3]
                else
                    pnd_temp = dN;

                // Check which element in arthat is closest to arH13 and assign
                // the PND for that size to that particle and zeros to the rest
                Index l;
                l=ar_input.nelem();

                Vector diff;
                
                for ( Index k=0, j=0; j<pnd_temp.nelem(); k+=l,j++ )
                {   
                    diff = ar[Range(k,l)];
                    
                    diff -= arH13[j];
                    
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

                // scale H13 distribution (which is independent of Ice or 
                // Snow massdensity) to current massdensity.
                /* JM120411: we don't need this - it's doing exactly, what
                             chk_pndsum is doing. instead we redefine verbosity
                             for checksum here to suppress the 'scaling is off'
                             warnings and let chk_pndsum do the proper scaling.
                scale_H13 ( pnd, IWC_field ( p,lat,lon ), vol, rho );*/

                // calculate proper scaling of pnd sum from real IWC and apply
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), vol, rho,
                             p, lat, lon, partfield_name, verbosity );
//                             p, lat, lon, partfield_name, temp_verb );
               
                // writing pnd vector to wsv pnd_field
                for ( Index i =0; i< npart; i++ )
                {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                }
            }
            else
            {
                for ( Index i = 0; i< npart; i++ )
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
 *  distribution for tropics. 
 *  For the estimation of the second moment the mass dimension
 *  relationship of Cotton et al., QJRMS, 2013 is used.
 To be used for snow. For this distribution the snow has to be as mass content.

\param pnd_field Particle number density field
\param SWC_field (snow) mass content field [kg/m3]
\param t_field atmospheric temperature [K]
\\param limits pnd_field boundaries (indices in p, lat, lon)
\param scat_meta_array particle meta data for particles
\param scat_data_start start index for particles handled by this distribution
\param npart number of particles handled by this distribution
\param part_string part_species tag for profile/distribution handled here
\param delim Delimiter string of *part_species* elements

\author Manfred Brath (most of the function is based on pnd_H11 (of J. Mendrok & D. Kreyling))
\date 2014-10-22

*/
void pnd_fieldF07TR (Tensor4View pnd_field,
                      const Tensor3& SWC_field,
                      const Tensor3& t_field,
                      const ArrayOfIndex& limits,
                      const ArrayOfScatteringMetaData& scat_meta_array,
                      const Index& scat_data_start,
                      const Index& npart,
                      const String& part_string,
                      const String& delim,
                      const Verbosity& verbosity)
{
  ArrayOfIndex intarr;
  Index pos;
  Vector dmax_unsorted ( npart, 0.0 );
  Vector vol ( npart, 0.0 );
  Vector dmax ( npart, 0.0 );
  Vector rho ( npart, 0.0 );
  Vector pnd ( npart, 0.0 );
  Vector dN ( npart, 0.0 );
  String partfield_name;
  
  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);
  
  for ( Index i=0; i < npart; i++ )
    dmax_unsorted[i] = ( scat_meta_array[i+scat_data_start].diameter_max );
  get_sorted_indexes(intarr, dmax_unsorted);
  
  // extract scattering meta data
  for ( Index i=0; i< npart; i++ )
  {
    pos = intarr[i]+scat_data_start;
    
    vol[i]= scat_meta_array[pos].volume; //[m^3]
    // get maximum diameter from meta data [m]
    dmax[i] = scat_meta_array[pos].diameter_max;
    // get density from meta data [kg/m^3]
    rho[i] = scat_meta_array[pos].density;
    // get aspect ratio from meta data [ ]
  }
  
  if (dmax.nelem() > 0)
    // dm.nelem()=0 implies no selected particles for the respective particle
    // field. should not occur anymore.
  {
    // itertation over all atm. levels
    for ( Index p=limits[0]; p<limits[1]; p++ )
    {
      for ( Index lat=limits[2]; lat<limits[3]; lat++ )
      {
        for ( Index lon=limits[4]; lon<limits[5]; lon++ )
        {
          // we only need to go through PSD calc if there is any material
          if (SWC_field ( p, lat, lon ) > 0.)
          {
            // iteration over all given size bins
            for ( Index i=0; i<dmax.nelem(); i++ ) //loop over number of particles
            {
              // calculate particle size distribution for H11
              // [# m^-3 m^-1]
              dN[i] = IWCtopnd_F07TR( dmax[i], t_field ( p, lat, lon ), SWC_field ( p, lat, lon ));
            }
            // scale pnds by scale width
            if (dmax.nelem() > 1)
              scale_pnd( pnd, dmax, dN ); //[# m^-3]
            else
              pnd = dN;
            
            
            // calculate proper scaling of pnd sum from real IWC and apply
            chk_pndsum ( pnd, SWC_field ( p,lat,lon ), vol, rho,
                        p, lat, lon, partfield_name, verbosity );
            
            // writing pnd vector to wsv pnd_field
            for ( Index i =0; i< npart; i++ )
            {
              pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                         lat-limits[2], lon-limits[4] ) = pnd[i];
            }
          }
          else
          {
            for ( Index i = 0; i< npart; i++ )
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

/*! Calculates the particle number density field for Gunn and Marshall (1958)
 size distribution. To be used for snow fall rates.
 
 \param pnd_field Particle number density field
 \param PR_field precipitation rate field [kg/(m2*s)]
 \\param limits pnd_field boundaries (indices in p, lat, lon)
 \param scat_meta_array particle meta data for particles
 \param scat_data_start start index for particles handled by this distribution
 \param npart number of particles handled by this distribution
 \param part_string part_species tag for profile/distribution handled here
 \param delim Delimiter string of *part_species* elements
 
 \author Manfred Brath (most of the function was taken by me from pnd_mp48 (of j. mendrok))
 \date 2014-10-20
 
 
 
 */
void pnd_fieldGM58 (Tensor4View pnd_field,
                    const Tensor3& PR_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfScatteringMetaData& scat_meta_array,
                    const Index& scat_data_start,
                    const Index& npart,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity)
{
    ArrayOfIndex intarr;
    Index pos;
    Vector vol_unsorted ( npart, 0.0 );
    Vector vol ( npart, 0.0 );
    Vector vol_me ( npart, 0.0 );
    Vector dve ( npart, 0.0 );
    Vector dme ( npart, 0.0 );
    Vector rho_snow ( npart, 0.0 );
    Vector pnd ( npart, 0.0 );
    Vector dN ( npart, 0.0 );
    
    String partfield_name;
    
    // Density of water, needed for melted Diameter
    const Numeric rho_water=1000.;// [kg/m^3]
    
    //split String and copy to ArrayOfString
    parse_partfield_name( partfield_name, part_string, delim);
    
    for ( Index i=0; i < npart; i++ )
        vol_unsorted[i] = ( scat_meta_array[i+scat_data_start].volume );
    get_sorted_indexes(intarr, vol_unsorted);
    

    
    // extract scattering meta data
    for ( Index i=0; i< npart; i++ )
    {
        pos = intarr[i]+scat_data_start;
        
        vol[i] = ( scat_meta_array[pos].volume ); //m^3
        
        // calculate volume equivalent diameter [m]
        dve[i] = pow(
                     ( scat_meta_array[pos].volume *6./PI ),
                     ( 1./3. ) );
        // get density from meta data [kg/m^3]
        rho_snow[i] = scat_meta_array[pos].density;
        
        // Calculate the melted Diameter
        dme[i] =pow( ( rho_snow[i]/rho_water ),( 1./3. ) )*dve[i];
        
        // calculate the mass equivalent(melted) Volume
        vol_me=( rho_snow[i]/rho_water )*vol[i];
        
        
    }
    
    
    Numeric fac, tPR, N0;
    Numeric PWC, lambda = NAN;
    
    
    // conversion factor from PR [kg/(m^2*s)] to PR[mm/hr]
    const Numeric convfac=3.6e6; // [PR [kg/(m^2*s)] to PR[mm/hr]*[kg/m3]
    
    
    
    // set parameterisation constants
    const Numeric lambda_fac = 25.5*1e2; // [cm^-1] converted to [m^-1] to fit d[m]
    const Numeric lambda_exp = -0.48;
    const Numeric N0_exp= -0.87;
    const Numeric N0_fac=0.038*1e8;
    
    
    
    if (dve.nelem() > 0)
        // dve.nelem()=0 implies no selected particles for the respective particle
        // field. should not occur.
    {
        // iteration over all atm. levels
        for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // we only need to go through PSD calc if there is any precipitation
                    if (PR_field ( p, lat, lon ) > 0.)
                    {
                        
                        fac = convfac/rho_water;
                        
                        // do PR [kg/(m^2*s)] to PR [mm/hr] conversion
                        tPR = PR_field ( p, lat, lon ) * fac;
                        
                        
                        //calculate Intercept
                        N0 = pow(tPR,N0_exp)*N0_fac; // [#/cm^3/cm] converted to [#/m^3/m]
                        
                        // get slope of distribution [m^-1]
                        lambda = lambda_fac * pow(tPR,lambda_exp);
                        
                        // derive particle number density for all given sizes
                        for ( Index i=0; i<dve.nelem(); i++ )
                        {
                            
                            dN[i] = N0 * exp(-lambda*dme[i]);
                            
                        }
                        
                        // scale pnds by bin width
                        if (dme.nelem() > 1)
                            scale_pnd( pnd, dme, dN );
                        else
                            pnd = dN;
                        
                      
                        
                        
                        // Make sure lambda was initialized in the while loop
                        assert(!isnan(lambda));
                        
                        // calculate error of pnd sum and real XWC
                        PWC = rho_water*PI*N0 / pow(lambda,4.);
                        chk_pndsum ( pnd, PWC, vol, rho_snow,
                                    p, lat, lon, partfield_name, verbosity );

                        
                        // writing pnd vector to wsv pnd_field
                        for ( Index i = 0; i< npart; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = pnd[i];
                        }
                    }
                    else
                        for ( Index i = 0; i< npart; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = 0.;
                        }
                }
            }
        }
    }
}



/*! Calculates the particle number density field for Sekhon and Srivastava (1970)
 size distribution. To be used for snow fall rates.
 
 \param pnd_field Particle number density field
 \param PR_field precipitation rate field [kg/(m2*s)]
 \\param limits pnd_field boundaries (indices in p, lat, lon)
 \param scat_meta_array particle meta data for particles
 \param scat_data_start start index for particles handled by this distribution
 \param npart number of particles handled by this distribution
 \param part_string part_species tag for profile/distribution handled here
 \param delim Delimiter string of *part_species* elements
 
 \author Manfred Brath (most of the function was taken from pnd_mp48 (of j. mendrok))
 \date 2014-10-22
 
 
 
 */
void pnd_fieldSS70 (Tensor4View pnd_field,
                    const Tensor3& PR_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfScatteringMetaData& scat_meta_array,
                    const Index& scat_data_start,
                    const Index& npart,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity)
{
    ArrayOfIndex intarr;
    Index pos;
    Vector vol_unsorted ( npart, 0.0 );
    Vector vol ( npart, 0.0 );
    //Vector vol_me ( npart, 0.0 );
    Vector dve ( npart, 0.0 );
    Vector dme ( npart, 0.0 );
    Vector rho_snow ( npart, 0.0 );
    Vector pnd ( npart, 0.0 );
    Vector pnd_old ( npart, 0.0 );
    Vector dN ( npart, 0.0 );
    
    Numeric fac, tPR, N0;//, temp,temp2,temp3;
    Numeric PWC, lambda = NAN;
    
    String partfield_name;
    
    // Density of water, needed for melted Diameter
    const Numeric rho_water=1000.;// [kg/m^3]
    
    //split String and copy to ArrayOfString
    parse_partfield_name( partfield_name, part_string, delim);
    
    for ( Index i=0; i < npart; i++ )
        vol_unsorted[i] = ( scat_meta_array[i+scat_data_start].volume );
    get_sorted_indexes(intarr, vol_unsorted);
    
    
    
    // extract scattering meta data
    for ( Index i=0; i< npart; i++ )
    {
        pos = intarr[i]+scat_data_start;
        
        vol[i] = ( scat_meta_array[pos].volume ); //m^3
        
        // calculate volume equivalent diameter [m]
        dve[i] = pow(
                     ( scat_meta_array[pos].volume *6./PI ),
                     ( 1./3. ) );
        
        // get density from meta data [kg/m^3]
        rho_snow[i] = scat_meta_array[pos].density;
        
        // Calculate the melted Diameter
        dme[i] =pow( ( rho_snow[i]/rho_water ),( 1./3. ) )*dve[i];
        
        // calculate the mass equivalent(melted) Volume
        //vol_me=( rho_snow[i]/rho_water )*vol[i];
        
        
    }
    
    
    
    // conversion factor from PR [kg/(m^2*s)] to PR[mm/hr]
    const Numeric convfac=3.6e6; // [PR [kg/(m^2*s)] to PR[mm/hr]*[kg/m3]
    
    
    
    // set parameterisation constants
    const Numeric lambda_fac = 22.9*1e2; // [cm^-1] converted to [m^-1] to fit d[m]
    const Numeric lambda_exp = -0.45;
    const Numeric N0_exp= -0.94;
    const Numeric N0_fac=0.025*1e8;
    
    
    
    if (dve.nelem() > 0)
        // dve.nelem()=0 implies no selected particles for the respective particle
        // field. should not occur.
    {
        // iteration over all atm. levels
        for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // we only need to go through PSD calc if there is any precipitation
                    if (PR_field ( p, lat, lon ) > 0.)
                    {
                        
                        fac = convfac/rho_water;
                        
                        // do PR [kg/(m^2*s)] to PR [mm/hr] conversion
                        tPR = PR_field ( p, lat, lon ) * fac;
                        
                        
                        
                        //calculate Intercept
                        N0 = pow(tPR,N0_exp)*N0_fac; // [#/cm^3/cm] converted to [#/m^3/m]
                        
                        // get slope of distribution [m^-1]
                        lambda = lambda_fac * pow(tPR,lambda_exp);
                        
                        // derive particle number density for all given sizes
                        for ( Index i=0; i<dve.nelem(); i++ )
                        {
                            dN[i] = N0 * exp(-lambda*dme[i]);
                        }
                        
                        
                        // scale pnds by bin width
                        if (dme.nelem() > 1)
                            scale_pnd( pnd, dme, dN );
                        else
                            pnd = dN;
                        
                        
                        
                        // Make sure lambda was initialized in the while loop
                        assert(!isnan(lambda));
                        
//                        pnd_old=pnd;
                        
                        // calculate error of pnd sum and real XWC
                        PWC = rho_water*PI*N0 / pow(lambda,4.);
                        chk_pndsum ( pnd, PWC, vol, rho_snow,
                                    p, lat, lon, partfield_name, verbosity );
                        
                        
                        // writing pnd vector to wsv pnd_field
                        for ( Index i = 0; i< npart; i++ )
                        {
//                            temp=pnd[i];
//                            temp2=dN[i];
//                            temp3=pnd_old[i];
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = pnd[i];
                            
                        }
                    }
                    else
                        for ( Index i = 0; i< npart; i++ )
                        {
                            pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                       lat-limits[2], lon-limits[4] ) = 0.;
                        }
                }
            }
        }
    }
}

/*! Calculates the particle number density field for Marshall and Palmer (1948)
    size distribution. To be used for precipitation, particularly rain.

    \param pnd_field Particle number density field
    \param PR_field precipitation rate field [kg/(m2*s)]
    \\param limits pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta_array particle meta data for particles
    \param scat_data_start start index for particles handled by this distribution
    \param npart number of particles handled by this distribution
    \param part_string part_species tag for profile/distribution handled here
    \param delim Delimiter string of *part_species* elements
  
  \author Jana Mendrok
  \date 2012-04-04

*/
void pnd_fieldMP48 (Tensor4View pnd_field,
                    const Tensor3& PR_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfScatteringMetaData& scat_meta_array,
                    const Index& scat_data_start,
                    const Index& npart,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity)
{
  ArrayOfIndex intarr;
  Index pos;
  Vector vol_unsorted ( npart, 0.0 );
  Vector vol ( npart, 0.0 );
  Vector d ( npart, 0.0 );
  Vector rho ( npart, 0.0 );
  Vector pnd ( npart, 0.0 );
  Vector pnd_old ( npart, 0.0 );
  Vector dN ( npart, 0.0 );

  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < npart; i++ )
      vol_unsorted[i] = ( scat_meta_array[i+scat_data_start].volume );
  get_sorted_indexes(intarr, vol_unsorted);
	
  // extract scattering meta data
  for ( Index i=0; i< npart; i++ )
  {
      pos = intarr[i]+scat_data_start;

      vol[i] = ( scat_meta_array[pos].volume ); //m^3
      // calculate volume equivalent diameter [m]
      d[i] = pow ( 
             ( scat_meta_array[pos].volume *6./PI ),
             ( 1./3. ) );
      // get density from meta data [kg/m^3]
      rho[i] = scat_meta_array[pos].density;
  }

  // conversion factor from PR [kg/(m^2*s)] to PR[mm/hr]
  /* for the conversion we need to know the mean density of distribution:
     rho_mean = mass_total/vol_total
              = int(vol[D]*rho[D]*dN[D])dD/int(vol[D]*dN[D])dD

     however, we do not yet know dN[D], which in turn is dependent on rho[D].
      proper solution requires iterative approach (does it?), but is it really
      worth to implement this? at least rain drops density should not really
      vary with particle size. might be different for snow/hail/graupel.
     so, here we set conversion factor to pseudo-PR[mm/hr]*[kg/m^3] and divide by
      density later on 
  */
  const Numeric convfac=3.6e6; // [PR [kg/(m^2*s)] to PR[mm/hr]*[kg/m3]
    Numeric fac, rho_mean, vol_total, mass_total, tPR;
      
  // set parameterisation constants here instead of in PRtopnd_MP48, since we
  // also need them for PR to PWC conversion
  const Numeric N0 = 0.08*1e8; // [#/cm^3/cm] converted to [#/m^3/m]
  const Numeric lambda_fac = 41.*1e2; // [cm^-1] converted to [m^-1] to fit d[m]
  const Numeric lambda_exp = -0.21;
  Numeric PWC, lambda = NAN;

  if (d.nelem() > 0)
  // d.nelem()=0 implies no selected particles for the respective particle
  // field. should not occur.
  {
      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            // we only need to go through PSD calc if there is any precipitation
            if (PR_field ( p, lat, lon ) > 0.)
            {
              Index n_it = 0;

              // initializing for proper start of while loop
              lambda = 0.;
              rho_mean = 0.;
              mass_total = rho.sum()/Numeric(npart);
              vol_total = 1.;
              while (abs( rho_mean/(mass_total/vol_total)-1.)>1e-2)
              // did bulk mean density change?
              {
                if (n_it>10)
                {
                    ostringstream os;
                    os<< "ERROR: Bulk mean density for MP48 distribution is"
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
                for ( Index i=0; i<d.nelem(); i++ )
                {
                    // calculate particle size distribution with MP48
                    // output: [# m^-3 m^-1]
                    //dN[i] = PRtopnd_MP48 ( tPR, d[i]);
                    // too much a hassle to have a separate function. so we do
                    // the calculation directly here.
                    dN[i] = N0 * exp(-lambda*d[i]);
                    //dN2[i] = dN[i] * vol[i] * rho[i];
                }

                // scale pnds by bin width
                if (d.nelem() > 1)
                    scale_pnd( pnd, d, dN );
                else
                    pnd = dN;

                // derive mass and volume over whole size distribution for
                // updated mean density
                mass_total = vol_total = 0.;
                for ( Index i=0; i<d.nelem(); i++ )
                {
                    mass_total += vol[i]*rho[i]*pnd[i];
                    vol_total += vol[i]*pnd[i];
                }
                //cout << n_it << ". iteration changing bulk density from "
                //     << rho_mean << " to " << mass_total/vol_total << "\n";
                n_it++;
              }
              //cout << "fine at p: " << p << " lat: " << lat << " lon" << lon
              //     << " for PR=" << PR_field ( p, lat, lon )*fac << "mm/hr.\n";

              // Make sure lambda was initialized in the while loop
              assert(!isnan(lambda));
                

              // calculate error of pnd sum and real XWC
              PWC = rho_mean*PI*N0 / pow(lambda,4.);
              chk_pndsum ( pnd, PWC, vol, rho,
                           p, lat, lon, partfield_name, verbosity );
	    
              // writing pnd vector to wsv pnd_field
              for ( Index i = 0; i< npart; i++ )
              {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = pnd[i];
              }
            }
            else
                for ( Index i = 0; i< npart; i++ )
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

    \param pnd_field Particle number density field
    \param LWC_field mass content (liquid water) field [kg/m3]
    \\param limits pnd_field boundaries (indices in p, lat, lon)
    \param scat_meta_array particle meta data for particles
    \param scat_data_start start index for particles handled by this distribution
    \param npart number of particles handled by this distribution
    \param part_string part_species tag for profile/distribution handled here
    \param delim Delimiter string of *part_species* elements
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2012-04-04

*/
void pnd_fieldH98 (Tensor4View pnd_field,
                   const Tensor3& LWC_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfScatteringMetaData& scat_meta_array,
                   const Index& scat_data_start,
                   const Index& npart,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
  ArrayOfIndex intarr;
  Index pos;
  Vector vol_unsorted ( npart, 0.0 );
  Vector vol ( npart, 0.0 );
  Vector r ( npart, 0.0 );
  Vector rho ( npart, 0.0 );
  Vector pnd ( npart, 0.0 );
  Vector dN ( npart, 0.0 );

  String partfield_name;

  //split String and copy to ArrayOfString
  parse_partfield_name( partfield_name, part_string, delim);

  for ( Index i=0; i < npart; i++ )
      vol_unsorted[i] = ( scat_meta_array[i+scat_data_start].volume );
  get_sorted_indexes(intarr, vol_unsorted);
      
  // extract scattering meta data
  for ( Index i=0; i< npart; i++ )
  {
      pos = intarr[i]+scat_data_start;

      vol[i]= scat_meta_array[pos].volume; // [m^3]
      // calculate radius from volume [m]
      r[i] = 0.5 * pow (
               ( 6*scat_meta_array[pos].volume/PI ), ( 1./3. ) );
      // get density from meta data [kg/m^3]
      rho[i] = scat_meta_array[pos].density;
  }

  if (r.nelem() > 0)
  // r.nelem()=0 implies no selected particles for the respective particle
  // field. should not occur anymore
  {
      // iteration over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
      {
        for ( Index lat=limits[2]; lat<limits[3]; lat++ )
        {
          for ( Index lon=limits[4]; lon<limits[5]; lon++ )
          {
            if (LWC_field ( p,lat,lon )>0.)
            {
                // iteration over all given size bins
                for ( Index i=0; i<r.nelem(); i++ ) //loop over number of particles
                {
                    // calculate particle size distribution for liquid
                    // [# m^-3 m^-1]
                    dN[i] = LWCtopnd ( LWC_field ( p,lat,lon ), rho[i], r[i] );
                    //dN2[i] = LWCtopnd2 ( r[i] );  // [# m^-3 m^-1]
                    //dN2[i] = dN[i] * vol[i] * rho[i];
                }

                // scale pnds by scale width. output: [# m^-3]
                if (r.nelem() > 1)
                    scale_pnd( pnd, r, dN );
                else
                    pnd = dN;
	    
                // calculate error of pnd sum and real XWC
                chk_pndsum ( pnd, LWC_field ( p,lat,lon ), vol, rho,
                             p, lat, lon, partfield_name, verbosity );

                // writing pnd vector to wsv pnd_field
                for ( Index i =0; i< npart; i++ )
                {
                    pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                                lat-limits[2], lon-limits[4] ) = pnd[i];
                    //dlwc[q] = pnd2[q]*vol[q]*rho[q];
                }
            }
            else
            {
                for ( Index i = 0; i< npart; i++ )
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


 
 


/*! Calculates particle size distribution using MH97 parametrization.
 *  Each diameter of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density.

    \return dN particle number density per diameter interval [#/m3*m]
 
    \param iwc atmospheric ice water content [kg/m3]
    \param dm melted diameter of scattering particle [m]
    \param t atmospheric temperature [K]
    \param density of the scattering particle [kg/m3]
    \param perturb PSD parameters according to their error statistics?
  
  \author Daniel Kreyling
  \date 2010-12-06

*/
Numeric IWCtopnd_MH97 (	const Numeric iwc,
			const Numeric dm,
			const Numeric t,
			const Numeric density,
			const bool noisy )
{
  // skip calculation if IWC is 0.0
  if ( iwc == 0.0 )
  {
    return 0.0;
  }

  //[kg/m3] -> [g/m3] as used by parameterisation
  Numeric ciwc = iwc*1e3;
  Numeric cdensity = density*1e3;

  Numeric sig_a=0.068, sig_b1=0.054;
  Numeric sig_b2=5.5e-3, sig_m=0.0029;
  Numeric sig_aamu=0.02, sig_bamu=0.0005;
  Numeric sig_abmu=0.023, sig_bbmu=0.5e-3;
  Numeric sig_aasigma=0.02, sig_basigma=0.5e-3;
  Numeric sig_absigma=0.023, sig_bbsigma=4.7e-4;

  if ( noisy )
  {
    Rng rng;                      //Random Number generator
    Index mc_seed;
    mc_seed = (Index)time(NULL);
    rng.seed(mc_seed, Verbosity());

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
  else
  {
    sig_a=0., sig_b1=0.;
    sig_b2=0., sig_m=0.;
    sig_aamu=0., sig_bamu=0.;
    sig_abmu=0., sig_bbmu=0.;
    sig_aasigma=0., sig_basigma=0;
    sig_absigma=0., sig_bbsigma=0.;
  }

  Numeric dN1;
  Numeric dN2;
  Numeric dN;

  // convert m to microns
  Numeric Dm = dm * 1e6;
  //convert T from Kelvin to Celsius
  Numeric T = t-273.15;
  //split IWC in IWCs100 and IWCl100

  // determine iwc<100um and iwc>100um
  Numeric a=0.252+sig_a; //g/m^3
  Numeric b1=0.837+sig_b1;
  Numeric IWCs100=min ( ciwc,a*pow ( ciwc,b1 ) );
  Numeric IWCl100=ciwc-IWCs100;


  //Gamma distribution component

  Numeric b2=-4.99e-3+sig_b2; //micron^-1
  Numeric m=0.0494+sig_m; //micron^-1
  Numeric alphas100=b2-m*log10 ( IWCs100 ); //micron^-1
  // for large IWC alpha, hence dN1, goes negative. avoid that.
  // towards this limit, particles anyway get larger 100um, i.e.,
  // outside the size region gamma distrib is intended for
  if (alphas100>0.)
  {
    Numeric Ns100 = 6*IWCs100 * pow ( alphas100,5. ) /
                    ( PI*cdensity*gamma_func ( 5. ) );//micron^-5
    Numeric Nm1 = Ns100*Dm*exp ( -alphas100*Dm ); //micron^-4
    dN1 = Nm1*1e18; // m^-3 micron^-1
  }
  else
  {
    dN1 = 0.;
  }



  //Log normal distribution component

  // for small IWC, IWCtotal==IWC<100 & IWC>100=0.
  // this will give dN2=NaN. avoid that by explicitly setting to 0
  if (IWCl100>0.)
  {
    //FIXME: Do we need to ensure mul100>0 and sigmal100>0?

    Numeric aamu=5.20+sig_aamu;
    Numeric bamu=0.0013+sig_bamu;
    Numeric abmu=0.026+sig_abmu;
    Numeric bbmu=-1.2e-3+sig_bbmu;
    Numeric amu=aamu+bamu*T;
    Numeric bmu=abmu+bbmu*T;
    Numeric mul100=amu+bmu*log10 ( IWCl100 );

    Numeric aasigma=0.47+sig_aasigma;
    Numeric basigma=2.1e-3+sig_basigma;
    Numeric absigma=0.018+sig_absigma;
    Numeric bbsigma=-2.1e-4+sig_bbsigma;
    Numeric asigma=aasigma+basigma*T;
    Numeric bsigma=absigma+bbsigma*T;
    Numeric sigmal100=asigma+bsigma*log10 ( IWCl100 );

    if ( (mul100>0.) & (sigmal100>0.) )
    {
      Numeric a1 = 6*IWCl100; //g/m^3
      Numeric a2 = pow ( PI,3./2. ) * cdensity*sqrt(2) *
                   exp( 3*mul100+9./2. * pow ( sigmal100,2 ) ) * 
                   sigmal100 * pow ( 1.,3 ) * Dm; //g/m^3/micron^4
      Numeric Nm2 = a1/a2 *
                    exp ( -0.5 * pow ( ( log ( Dm )-mul100 ) /sigmal100,2 ) );
                    //micron^-4
      dN2 = Nm2*1e18; // m^-3 micron^-1
    }
    else
    {
      dN2 = 0.;
/*      ostringstream os;
      os<< "ERROR: not a valid MH97 size distribution:\n"
        << "mu>100=" << mul100 << " resulting from amu="<< amu << " and bmu="<< bmu << "\n"
        << "(aamu="<<aamu<<", bamu="<<bamu<<", abmu="<<abmu<<", bbmu="<<bbmu<<")\n"
        << "sigma>100=" << sigmal100 << " resulting from asigma="<< asigma << " and bsigma="<< bsigma << "\n"
        << "(aasigma="<<aasigma<<", basigma="<<basigma<<", absigma="<<absigma<<", bbsigma="<<bbsigma<<")\n";
      throw runtime_error ( os.str() );      */
    }
  }
  else
  {
    dN2 = 0.;
  }



  dN = ( dN1+dN2 ) *1e6; // m^-3 m^-1
  if (isnan(dN)) dN = 0.0;
  return dN;
}



/*! Calculates particle size distribution using H11 parametrization.
 *  Each diameter of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density.  

    \return dN particle number density per diameter interval [#/m3/m]
          
    \param d maximum diameter of scattering particle [m]
    \param t atmospheric temperature [K]
  
  \author Daniel Kreyling
  \date 2011-10-28

*/
Numeric IWCtopnd_H11 ( const Numeric d,
                       const Numeric t)
{  
  Numeric dN;
  Numeric la;
  Numeric mu;
  // convert m to cm
  Numeric D = d * 1e2;
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

  dN = pow( D, mu ) * exp ( -la * D );

  if (isnan(dN)) dN = 0.0;
  return dN;
}



/*! Calculates particle size distribution using H13 parametrization.
 *  Each diameter of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density.  

    \return dN particle number density per diameter interval [#/m3/m]
          
    \param d maximum diameter of scattering particle [m]
    \param t atmospheric temperature [K]
  
  \author Johan Strandgren, Daniel Kreyling  
  \date 2013-08-26

*/
Numeric IWCtopnd_H13 ( const Numeric d,
                       const Numeric t)
{  
  Numeric dN;
  Numeric la;
  Numeric mu;
  // convert m to cm
  

  Numeric D = d * 1e2;
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

  dN = pow( D, mu ) * exp ( -la * D );

  if (isnan(dN)) dN = 0.0;
  return dN;
}



/*! Calculates particle size distribution using H13 parametrization.
 *  Each diameter of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density.  

    \return dN particle number density per diameter interval [#/m3/m]
          
    \param d maximum diameter of scattering particle [m]
    \param t atmospheric temperature [K]
  
  \author Johan Strandgren  
  \date 2013-08-26

*/
Numeric IWCtopnd_H13Shape ( const Numeric d,
                       const Numeric t)
{  
  Numeric dN;
  Numeric la;
  Numeric mu;
  // convert m to cm
  

  Numeric D = d * 1e2;
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

  dN = pow( D, mu ) * exp ( -la * D );

  if (isnan(dN)) dN = 0.0;
  return dN;
}



/*! Calculates area ratio from the temperature and maximum diameter using H13 parametrization.
 *  Each diameter of the scattering particles is a node in the aspect ratio distribution.
 *  One call of this function calculates one aspect ratio.  

    \return dN particle number density per diameter interval [#/m3/m]
          
    \param d maximum diameter of scattering particle [m]
    \param t atmospheric temperature [K]
  
  \author Johan Strandgren  
  \date 2013-08-26

*/
Numeric area_ratioH13 ( const Numeric d,
                        const Numeric t)
{  
  Numeric Ar;
  Numeric alpha;
  Numeric beta;
  // convert m to cm
  

  Numeric D = d * 1e2;
  //convert T from Kelvin to Celsius
  Numeric T = t-273.15;
  //Parameterize for all temperatures
  
  alpha = 0.25*exp(0.0161*T);
  
  beta = -0.25+0.0033*T;
  
  // Area ratio function depending on temperature

  Ar = alpha*pow(D,beta);

  if (isnan(Ar)) Ar = 0.0;
  return Ar;
}

/*! Calculates particle size distribution using F07 parametrization for
 *  tropics. For the estimation of the second moment the mass dimension 
 *  relationship of Cotton et al., QJRMS, 2013 is used.
 *  Each diameter of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param d maximum diameter of scattering particle [m]
 \param t atmospheric temperature [K]
 \param swc snow water content [kg/m^3]
 
 \author Manfred Brath
 \date 2014-10-22
 
 */
Numeric IWCtopnd_F07TR ( const Numeric d, const Numeric t,
                         const Numeric swc)
{
    Numeric dN;
    Numeric An;
    Numeric Bn;
    Numeric Cn;
    Numeric M2;
    Numeric Mn;
    Numeric x;

    
    // Factor for the Mass-Dimension relationship taken from Cotton et al., QJRMS, 2013, m=alpha*(Dmax/D0)^2 [kg]
    const Numeric alpha=0.0257;
    
    //factors of phi23
    Vector q=MakeVector(152,-12.4,3.28,-0.78,-1.94);
    
    //Factors of factors of the moment estimation parametrization
    Vector Aq=MakeVector(13.6,-7.76,0.479);
    Vector Bq=MakeVector(-0.0361,0.0151,0.00149);
    Vector Cq=MakeVector(0.807,0.00581,0.0457);
    
    
    //convert T from Kelvin to Celsius
    Numeric T = t-273.15;
    
    // estimate second moment
    M2=swc/alpha;
    
    // order of the moment parametrization
    const Numeric n=3;
    
    // calculate factors of the moment estimation parametrization
    An=exp(Aq[0]+Aq[1]*n+Aq[2]*pow(n,2));
    Bn=Bq[0]+Bq[1]*n+Bq[2]*pow(n,2);
    Cn=Cq[0]+Cq[1]*n+Cq[2]*pow(n,2);
    
    // moment parametrization
    Mn=An*exp(Bn*T)*pow(M2,Cn);
    
    //Define x
    x=d*M2/Mn;
    
    //Distribution function
    dN = q[0]*exp(q[1]*x)+q[2]*pow(x,q[3])*exp(q[4]*x);
    
    if (isnan(dN)) dN = 0.0;
    return dN;
}


/*! Calculates particle size distribution for liquid water clouds using a gamma
 *  parametrization by Hess et al., 1998 (continental stratus).
 *  Each radius of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density. 

	\return n particle number density per radius interval [#/m3*m]
         
	\param lwc atmospheric liquid water content [kg/m3]
	\param r radius of scattering particle [m]
  
  \author Daniel Kreyling
  \date 2010-12-16

*/
Numeric LWCtopnd (const Numeric lwc, //[kg/m^3]
		  const Numeric density, //[kg/m^3]
		  const Numeric r // [m]
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
	Numeric A = 0.75/PI * lwc/density * gam * pow(B,a4g) /
                    gamma_func(a4g);
	Numeric n = A * (pow(r,alpha) * exp(-B*pow(r,gam)));
	// n in [# m^-3 m^-1]
	
	if (isnan(n)) n = 0.0;
	return n;
}

// ONLY FOR TESTING PURPOSES
Numeric LWCtopnd2 (
		   const Numeric r // [m]
		  )
{ 	
	Numeric rc = 4.7; //micron
	Numeric alpha = 5.0;
	Numeric gam = 1.05;
	
	Numeric B=(alpha/gam)/pow(rc,gam); 
	Numeric A=gam*pow(B,((alpha+1)/gam))/gamma_func((alpha+1)/gam);
	Numeric n=A*(pow(r*1e6,alpha)*exp(-B*pow(r*1e6,gam)))*1e6;
	// [# m^-3 m^-1]
	
	if (isnan(n)) n = 0.0;
	return n;
}



/*! Calculates particle size distribution using MP48 parametrization for rain.
 *  Each diameter of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density.  

    \return dN particle number density per diameter interval [#/m3*m]
          
    \param R precipitation rate [mm/hr]
    \param D volume equivalent diameter of scattering particle [m]
 
  \author Jana Mendrok
  \date 2012-04-04

*/
Numeric PRtopnd_MP48 (	const Numeric R,
			const Numeric D)
{
  // skip calculation if PR is 0.0
  if ( R == 0.0 )
  {
    return 0.0;
  }

  Numeric N0 = 0.08*1e-2; // [#/cm^3/cm] converted to [#/m^3/um]
  Numeric lambda = 41.*1e2*pow(R,-0.21); // [cm^-1] converted to [m^-1] to fit D[m]

  Numeric n = N0*exp(-lambda*D);
  return n;
}



/*! Scaling pnd values by width of size bin. 
 * Bin width is determined from preceeding and following particle size.
 * Vector y and x must be equal in size. Vector w holds the weights.
 * Derived from trapezoid integration rule.
         
    \param w weights
    \param x e.g. particle radius [m]
    \param y e.g. particle number density per radies interval [#/m3*m]
  
  \author Daniel Kreyling
  \date 2010-12-15

*/
void scale_pnd  (  Vector& w,
		   const Vector& x,
		   const Vector& y) 
{
    // check if vectors have same size
    if (x.nelem() != y.nelem()) 
    {
        throw std::logic_error("x and y must be the same size");
    }
    
    if (x.nelem()>1) // calc. integration weights (using trapezoid integration)
    {
        for (Index i = 0; i<x.nelem(); i++)
        {
            if (i == 0) // first value
            {
                w[i] = 0.5*(x[i+1]-x[i])*y[i]; // m^-3
            }
            else if (i == x.nelem()-1) //last value
            {
                w[i] = 0.5*(x[i]-x[i-1])*y[i]; // m^-3
            }
            else // all values inbetween
            {
                w[i] = 0.5*(x[i+1]-x[i-1])*y[i]; // m^-3
            }
        }
    }
    else // for monodisperse pnd=dN
    {
      w[0] = y[0];
    }
}



/*! Check sum of pnd vector against total mass density value.
 *  Deviation is calculated and used to adjust the output of vector pnd.
         
	\param pnd particle number density [#/m3]
	\param xwc atmospheric massdensity [kg/m3]
	\param density scattering particle density [kg/m3]
	\param vol scattering particle volume [m3]
  
  \author Daniel Kreyling
  \date 2010-12-15

*/
void chk_pndsum (Vector& pnd,
                 const Numeric xwc,
                 const Vector& vol,
                 const Vector& density,
                 const Index& p,
                 const Index& lat,
                 const Index& lon,
                 const String& partfield_name,
                 const Verbosity& verbosity)

{
  CREATE_OUT2;
  
  // set vector x to pnd size
  Vector x ( pnd.nelem(), 0.0 );
  Numeric error;

  //cout << "p = " << p << ", pnd.nelem:" << pnd.nelem() << ", xwc: " << xwc << "\n";
  for ( Index i = 0; i<pnd.nelem(); i++ )
  {
    // convert from particles/m^3 to kg/m^3
    x[i] = pnd[i]*density[i]*vol[i];
    /*cout<< "p = " << p << ", i: " << i << "\n"
        << "pnd[i]: " << pnd[i] << ", density[i]: " << density[i] << ", vol[i]: " << vol[i] << "\n"
        << "x[i]: " << x[i] << "\n";*/
  }

  //cout<<"at p = "<<p<<", lat = "<<lat<<", lon = "<<lon
  //<<" given mass density: "<<xwc<<", calc mass density: "<<x.sum() << "\n";
  if ( x.sum() == 0.0 )
    if ( xwc == 0.0 )
    {
      // set error and all pnd values to zero, IF there is 
      // no scattering particles at this atmospheric level.
      error = 0.0;
      pnd = 0.0;
    }
    else
    { // when x.sum()==0, but xwc!=0, obviously something went wrong in pnd calc
      ostringstream os;
      os<< "ERROR: in WSM chk_pndsum in pnd_fieldSetup!\n" 
      << "Given mass density != 0, but calculated mass density == 0.\n"
      << "Seems, something went wrong in pnd_fieldSetup. Check!\n"
      << "The problem occured for profile '"<< partfield_name <<"' at: "
      << "p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<".\n";
     throw runtime_error ( os.str() );
    }
  else
  {
    error = xwc/x.sum();
  // correct all pnd values with error
    pnd *= error;
  // give warning if deviations are more than 10%
    if ( error > 1.10 || error < 0.90 )
    {
      CREATE_OUT1;
      //ostringstream os;
      out1<< "WARNING: in WSM chk_pndsum in pnd_fieldSetup!\n" 
      << "The deviation of the sum of nodes in the particle size distribution\n"
      << "to the initial input mass density (IWC/LWC) is larger than 10%!\n"
      << "The deviation of: "<< error-1.0<<" occured in the atmospheric level: "
      << "p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<".\n";
      //cerr<<os;
    }
  }
  //cout<<", error: "<<error<<"corrected to: "<<x.sum()*error<<"\n";

  out2 << "PND scaling factor in atm. level "
       << "(p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<"): "<< error <<"\n";
    //cout<<"\npnd_scaled\t"<<pnd<<endl;
    //cout<<"\nPND\t"<<pnd.sum()<<"\nXWC\t"<<xwc<<"\nerror\t"<<error<<endl;
    //cout<<"\n"<<x.sum()<<"\n"<<xwc<<"\n"<<error<<endl;
}



/*! The H11 PSD is scaled to the initial 'ice' or 'snow' massdensity, after
 *  the distribution has been evaluated. This function applies the scaling.
         
	\param xwc atmospheric massdensity [kg/m3]
	\param density scattering particle density [kg/m3]
	\param vol scattering particle volume [m3]
  
  \author Daniel Kreyling
  \date 2011-10-31

*/
void scale_H11 (
		 Vector& pnd,
		 const Numeric xwc,
		 const Vector& density,
		 const Vector& vol)
		 // const Index& p,
		 // const Index& lat,
		 // const Index& lon,
		 // const Verbosity& verbosity)

{
  // set vector x to pnd size
  Vector x (pnd.nelem(), 0.0);

  for ( Index i = 0; i<pnd.nelem(); i++ )
  {
    // convert from particles/m^3 to kg/m^3
    x[i] = pnd[i]*density[i]*vol[i];
    //out0<<x[i]<<"\n"<< pnd[i]<< "\n";
  }

  if ( x.sum() == 0.0 )
  {
    // set ratio and all pnd values to zero, IF there are 
    // no scattering particles at this atmospheric level.
    pnd = 0.0;
  }
  else
  {
      // calculate the ratio of initial massdensity (xwc) to sum of pnds
    const Numeric ratio = xwc/x.sum();
    // scale each pnd to represent the total massdensity
    pnd *= ratio;
  }

}

/*! The H13 PSD is scaled to the initial 'ice' or 'snow' massdensity, after
 *  the distribution has been evaluated. This function applies the scaling.
         
    \param xwc atmospheric massdensity [kg/m3]
    \param density scattering particle density [kg/m3]
    \param vol scattering particle volume [m3]
  
  \author Johan Strandgren, Daniel Kreyling
  \date 2013-09-03

*/
void scale_H13 (
         Vector& pnd,
         const Numeric xwc,
         const Vector& density,
         const Vector& vol)
         // const Index& p,
         // const Index& lat,
         // const Index& lon,
         // const Verbosity& verbosity)

{
  // set vector x to pnd size
  Vector x (pnd.nelem(), 0.0);

  for ( Index i = 0; i<pnd.nelem(); i++ )
  {
    // convert from particles/m^3 to kg/m^3
    x[i] = pnd[i]*density[i]*vol[i];
    //out0<<x[i]<<"\n"<< pnd[i]<< "\n";
  }

  if ( x.sum() == 0.0 )
  {
    // set ratio and all pnd values to zero, IF there are 
    // no scattering particles at this atmospheric level.
    pnd = 0.0;
  }
  else
  {
      // calculate the ratio of initial massdensity (xwc) to sum of pnds
    const Numeric ratio = xwc/x.sum();
    // scale each pnd to represent the total massdensity
    pnd *= ratio;
  }

}


/*! Splitting part_species string and parse type of massdensity_field

	\param  partfield_name name of atmospheric particle field
	\param  part_string containing infos about scattering particle calculations
  \param delim Delimiter string of *part_species* elements

  \author Daniel Kreyling
  \date 2011-02-21

*/
void parse_partfield_name (//WS Output:
                      String& partfield_name,
                      // WS Input:
                      const String& part_string,
                      const String& delim)
{
  ArrayOfString strarr;

  // split part_species string at delim and write to ArrayOfString
  part_string.split ( strarr, delim );

  //first entry is particle field name (e.g. "IWC", "LWC" etc.)
  if (strarr.size()>0 && part_string[0]!=delim[0])
  {
      partfield_name = strarr[0];
  }
  else
  {
      ostringstream os;
      os << "No information on particle field name in '" << part_string << "'\n";
      throw runtime_error ( os.str() );

  }

  /* no restrictions on profile naming anymore
  if (  partfield_name != "IWC" && 
	partfield_name != "Snow" &&
	partfield_name != "LWC" &&
	partfield_name != "Rain" )
  {
    ostringstream os;
    os << "First substring in " << part_string
       << " must be rather 'LWC', 'IWC', 'Rain' or 'Snow'\n"
       <<"Check input in *part_species*!\n";
    throw runtime_error ( os.str() );
  }*/
}



/*! Splitting part_species string and parse psd_param
	\param psd_param particle size distribution parametrization
	\param part_string containing infos about scattering particle calculations
  \param delim Delimiter string of *part_species* elements
  
  \author Daniel Kreyling
  \date 2011-02-21

*/
void parse_psd_param (//WS Output:
                      String& psd_param,
                      // WS Input:
                      const String& part_string,
                      const String& delim)
{
  ArrayOfString strarr;

  // split part_species string at delim and write to ArrayOfString
  part_string.split ( strarr, delim );

  // second entry is particle size distribution parametrisation  ( e.g."MH97")
  // check, whether we have a second entry
  if (strarr.size()>1)
      psd_param = strarr[1];
  else
      psd_param = "";

/*  // jm120218: FIX! <- DONE (120401)
  // this should not be checked here, but in pnd_fieldSetup (e.g., via
  // a case switch default)
  if ( psd_param.substr(0,4) != "MH97" && psd_param != "liquid" &&
       psd_param != "H11" )
  {
    ostringstream os;
    os <<"The chosen PSD parametrisation in " << part_string
       << " can not be handeled in the moment.\n"
       <<"Choose either 'MH97', 'H11' or 'liquid'!\n" ;
    throw runtime_error ( os.str() );
  }*/
}



/*! Splitting part_species string and parse type of particle (from ScattData array)

	\param  part_type particle type (material, phase). 
	\param  part_string containing infos about scattering particle calculations
  \param delim Delimiter string of *part_species* elements

  \author Jana Mendrok
  \date 2012-04-03

*/
void parse_part_material (//WS Output:
                      String& part_material,
                      // WS Input:
                      const String& part_string,
                      const String& delim)
{
  ArrayOfString strarr;

  // split part_species string at delim and write to ArrayOfString
  part_string.split ( strarr, delim );

  // third entry is requested particle (material) type ( e.g."Water", "Ice")
  // check, whether we have a third entry
  if (strarr.size()>2)
  {
      part_material = strarr[2];
  }
  else
  {
      ostringstream os;
      os << "No information on particle type in '" << part_string << "'\n";
      throw runtime_error ( os.str() );
  }
}



/*! Splitting part_species string and parse min and max particle radius
	\param sizemin min scattering particle radius
	\param sizemax max scattering particle radius
	\param part_string containing infos about scattering particle calculations
  \param delim Delimiter string of *part_species* elements
  
  \author Daniel Kreyling
  \date 2011-02-21

*/
void parse_part_size (//WS Output:
                      Numeric& sizemin,
                      Numeric& sizemax,
                      // WS Input:
                      const String& part_string,
                      const String& delim)
{
  ArrayOfString strarr;

  // split part_species string at delim and write to ArrayOfString
  part_string.split ( strarr, delim );
  
  // convert String for size range, into Numeric
  // 1. third entry is minimum particle radius
  if ( strarr.nelem() < 4 || strarr[3] == "*" )
  {
    sizemin = 0.;
  }
  else
  {
    istringstream os1 ( strarr[3] );
    os1 >> sizemin;
  }
  // 2. fourth entry is maximum particle radius
  if ( strarr.nelem() < 5 || strarr[4] == "*" )
  {
    sizemax = -1.;
  }
  else
  {
    istringstream os2 ( strarr[4] );
    os2 >> sizemax;
  }
}

