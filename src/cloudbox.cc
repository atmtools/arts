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
extern const Numeric DENSITY_OF_ICE;
extern const Numeric DENSITY_OF_WATER;

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



//! Check whether field of a specific scattering species zero everywhere.
/*!
  \return empty_flag        flag whether all field entries are zero
  \param scat_species_field scattering species field (e.g. mass density,
                             mass flux, total number density)
  \param fieldname          name of scattering species field (just for info)
  \param dim                the atmosphere dimension 
  \param p_grid             pressure grid of current atmosphere
  \param lat_grid           latitude grid of current atmosphere
  \param lon_grid           longitude grid of current atmosphere

  \author Daniel Kreyling
  \date   2011-01-27
*/  
void chk_scat_species_field(bool& empty_flag,
                            const Tensor3& scat_species_field, 
                            const String& fieldname,
                            const Index&  dim,	
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid)
{
  
  // check p
  if ( scat_species_field.npages() != p_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *p_grid* (" << p_grid.nelem()
         << ") is unequal the number of pages of *" << fieldname << "* ("
         << scat_species_field.npages() << ").";
      throw runtime_error(os.str() );
  }
  
  // check lat
  if(dim >= 2 )
  {
    if ( scat_species_field.nrows() != lat_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *lat_grid* (" << lat_grid.nelem()
         << ") is unequal the number of rows of *" << fieldname << "* ("
         << scat_species_field.nrows() << ").";
      throw runtime_error(os.str() );
      
    }
  }
  
  // check lon
  if(dim == 3 )
  {
    if ( scat_species_field.ncols() != lon_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *lon_grid* (" << lon_grid.nelem()
         << ") is unequal the number of columns of *" << fieldname << "* ("
         << scat_species_field.ncols() << ").";
      throw runtime_error(os.str() );
      
    }
  }
  
  empty_flag = false;
  // set empty_flag to true if a single value of hydromet_field is unequal zero    
    for (Index j=0; j<scat_species_field.npages(); j++) {
      for (Index k=0; k<scat_species_field.nrows(); k++) {
	    for (Index l=0; l<scat_species_field.ncols(); l++) {
	      if ( scat_species_field(j,k,l) != 0.0) empty_flag = true;
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
  
  \param pnd_field_raw   pnd field raw data (array for all scattering elements)
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
                               << "/m³ at p=" << p << " Pa for scattering "
                               << "element #" << n << ".";
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
                                   << "/m³ at lat=" << lat << "° for scattering "
                                   << "element #" << n << ".";
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
                                   << "/m³ at lon=" << lon << "° for scattering "
                                   << "element #" << n << ".";
                                throw runtime_error(os.str());
                            }
                        }
                    }
                }
            }
        }
    }
}



//!  Check validity of scat_species setting
/*!
  This function checks, whether number of elements in each scattering species
  string is ok, and whether the entries for size limits are indeed numbers (or
  '*').

	\param scat_species Array of scattering species tags.
  \param delim        delimiter string of *scat_species* elements.

  \author Jana Mendrok
  \date 2012-10-25

*/
void chk_scat_species (
                      const ArrayOfString& scat_species,
                      const String& delim
                      )
{
  ArrayOfString strarr;
  Numeric sizeval;

  for ( Index k=0; k<scat_species.nelem(); k++ )
    {
      scat_species[k].split ( strarr, delim );
      if ( strarr.nelem() > 4 )
        {     
          ostringstream os;
          os << "Individual strings in scat_species can only contain up to 4\n"
             << "elements, but entry #" << k << " contains the following "
             << strarr.nelem() << ":\n" << strarr << "\n";
          throw runtime_error ( os.str() );
        }

      if ( strarr.nelem() > 2 && strarr[2] != "*" )
        {
          istringstream os1 ( strarr[3] );
          os1 >> sizeval;
          if (os1.fail())
            {
              ostringstream os;
              os << "Sizemin specification can only be a number or wildcard ('*')"
                 << ", but is '" << strarr[2] << "'\n";
              throw runtime_error ( os.str() );
            }
        }

      if ( strarr.nelem() > 3 && strarr[3] != "*" )
        {
          istringstream os1 ( strarr[3] );
          os1 >> sizeval;
          if (os1.fail())
            {
              ostringstream os;
              os << "Sizemax specification can only be a number or wildcard ('*')"
                 << ", but is '" << strarr[3] << "'\n";
              throw runtime_error ( os.str() );
            }
        }
    }
}



//! Check scattering data general
/*!
  FIXME
  
  \param scat_data Array of single scattering data
  \param scat_meta Array of scattering meta data

  \author Daniel Kreyling
  \date 2010-12-02
*/

void chk_scattering_data(const ArrayOfSingleScatteringData& scat_data,
                         const ArrayOfScatteringMetaData& scat_meta,
                         const Verbosity&)
{
  if (scat_data.nelem() != scat_meta.nelem())
  {
    ostringstream os;
    os << "The number of elements in *scat_data* and *scat_meta* do not match.\n"
       << "Each SingleScattering file must correspond to one ScatteringMeta"
       << " data file.";
    throw runtime_error( os.str());
  }

}

//! Check scattering data meta files
/*!
  FIXME
  
  \param scat_meta_single scattering meta data
  \param scat_meta_file filename of the data to be checked

  \author Daniel Kreyling
  \date 2010-12-02
*/
void chk_scattering_meta_data(const ScatteringMetaData& scat_meta_single _U_,
                              const String& scat_meta_file,
                              const Verbosity& verbosity)
{
  CREATE_OUT3;
  out3 << "  Check scattering meta data file "<< scat_meta_file 
       << "\n";

  /* this check is outdated. type now is free from!
   however, we might want to have other things checked here!?
   - which parameters at least are needed? -> radius, ...?
   - ...
  if  (scat_meta_single.type != "Ice" && scat_meta_single.type != "Water" && scat_meta_single.type != "Aerosol")
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
  properties of a scattering element includes the required frequencies and 
  temperatures and whether the angular grids are defined correctly.
  Furthermore it checks the self consistency of the data by checking the 
  dimensions of pha_mat, ext_mat and abs_vec depending on the ptype case.
  
  \param scat_data Single scattering data
  \param scat_data_file Filename of the data to be checked.
  \param f_grid        Frequency grid
  
  \author Claudia Emde
  \date   2005-04-04
*/
void chk_scat_data(const SingleScatteringData& scat_data,
                   const String& scat_data_file,
                   ConstVectorView f_grid,
                   const Verbosity& verbosity)
{
  CREATE_OUT2;
  out2 << "  Check single scattering properties file "<< scat_data_file 
       << "\n";

  assert(scat_data.ptype == PTYPE_GENERAL ||
         scat_data.ptype == PTYPE_MACROS_ISO ||
         scat_data.ptype == PTYPE_HORIZ_AL);

  chk_interpolation_grids("scat_data.f_grid to f_grid",
			    scat_data.f_grid,
			    f_grid);
  
/*  if (!(scat_data.f_grid[0] <= f_grid[0] &&
        last(f_grid) <= 
        last(scat_data.f_grid) ))
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

  if (!(0. < scat_data.T_grid[0] && last(scat_data.T_grid) < 1001.))
    {
      ostringstream os;
      os << "The temperature values in " <<  scat_data_file 
         << " are negative or very large. Check whether you have used the "
         << "right unit [Kelvin].";
      throw runtime_error( os.str() );
    }
  
  if (scat_data.za_grid[0] != 0.)
    {
      ostringstream os;
      os << "The first value of the za grid in the single" 
         << " scattering properties datafile " 
         << scat_data_file << " must be 0.";
        throw runtime_error( os.str() );
    } 

  if (last(scat_data.za_grid) != 180.)
    {
      ostringstream os;
      os << "The last value of the za grid in the single"
         << " scattering properties datafile " 
         << scat_data_file << " must be 180.";
      throw runtime_error( os.str() );
    } 
  
  if (scat_data.ptype == PTYPE_GENERAL && scat_data.aa_grid[0] != -180.)
     {
       ostringstream os;
       os << "For ptype = \"general\" the first value"
          << " of the aa grid in the single scattering"
          << " properties datafile " 
          << scat_data_file << "must be -180.";
         throw runtime_error( os.str() );
     } 
  
  if (scat_data.ptype == PTYPE_HORIZ_AL && scat_data.aa_grid[0] != 0.)
    {
      ostringstream os;
      os << "For ptype = \"horizontally_aligned\""
         << " the first value"
         << " of the aa grid in the single scattering"
         << " properties datafile " 
         << scat_data_file << "must be 0.";
        throw runtime_error( os.str() );
    }   
  
  if (scat_data.ptype == PTYPE_HORIZ_AL && last(scat_data.aa_grid) != 180.)
    {
      ostringstream os;
      os << "For ptype = \"horizontally_aligned\""
         << " the last value of the aa grid in the single"
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
  
  switch (scat_data.ptype){
    
  case PTYPE_GENERAL:
    
    out2 << "  Datafile is for arbitrarily orientated particles. \n";
    
    chk_size(os_pha_mat.str(), scat_data.pha_mat_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             scat_data.za_grid.nelem(), scat_data.aa_grid.nelem(),
             scat_data.za_grid.nelem(), scat_data.aa_grid.nelem(), 
              16); 
    
    chk_size(os_ext_mat.str(), scat_data.ext_mat_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             scat_data.za_grid.nelem(), scat_data.aa_grid.nelem(),
             7);
    
    chk_size(os_abs_vec.str(), scat_data.abs_vec_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             scat_data.za_grid.nelem(), scat_data.aa_grid.nelem(),
             4);
    break;
    
  case PTYPE_MACROS_ISO: 
    
    out2 << "  Datafile is for randomly oriented particles, i.e., "
         << "macroscopically isotropic and mirror-symmetric scattering "
         << "media. \n";
    
    chk_size(os_pha_mat.str(), scat_data.pha_mat_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             scat_data.za_grid.nelem(), 1, 1, 1, 6);
    
    chk_size(os_ext_mat.str(), scat_data.ext_mat_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             1, 1, 1);
    
    chk_size(os_abs_vec.str(), scat_data.abs_vec_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             1, 1, 1);
    break; 
    
  case PTYPE_HORIZ_AL:
    
    out2 << "  Datafile is for horizontally aligned particles. \n"; 
    
    chk_size(os_pha_mat.str(), scat_data.pha_mat_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             scat_data.za_grid.nelem(), scat_data.aa_grid.nelem(),
             scat_data.za_grid.nelem()/2+1, 1, 
             16); 

    chk_size(os_ext_mat.str(), scat_data.ext_mat_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             scat_data.za_grid.nelem()/2+1, 1, 
             3);
    
    chk_size(os_abs_vec.str(), scat_data.abs_vec_data,
             scat_data.f_grid.nelem(), scat_data.T_grid.nelem(),
             scat_data.za_grid.nelem()/2+1, 1, 
             2);
    break;

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

  //split String and copy to ArrayOfString
  parse_psd_param( psd_param, part_string, delim);
  parse_partfield_name( partfield_name, part_string, delim);

  bool noisy = (psd_param == "MH97n");

  for ( Index i=0; i < N_se; i++ )
    {
      if ( isnan(scat_meta[scat_species][i].mass) )
        {
          ostringstream os;
          os << "Use of size distribution MH97 (as requested for\n"
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
      diameter_mass_equivalent[i] = pow(6.*mass[i]/PI/DENSITY_OF_ICE,1./3.); // [m]
  }
  
  if (mass.nelem() > 0)
  // mass.nelem()=0 implies no selected scattering element for the respective
  // scattering species field. should not occur.
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
                for ( Index i=0; i<diameter_mass_equivalent.nelem(); i++ )
                {
                  // calculate particle size distribution with MH97
                  // [# m^-3 m^-1]
                  dNdD[i] = IWCtopnd_MH97 ( IWC_field ( p, lat, lon ),
                                            diameter_mass_equivalent[i],
                                            t_field ( p, lat, lon ),
                                            noisy );
                }
            
                // scale pnds by bin width
                if (diameter_mass_equivalent.nelem() > 1)
                  scale_pnd( pnd, diameter_mass_equivalent, dNdD );
                else
                  pnd = dNdD;
	    
                // calculate error of pnd sum and real XWC
                chk_pndsum ( pnd, IWC_field ( p,lat,lon ), mass,
                             p, lat, lon, partfield_name, verbosity );
	    
                // writing pnd vector to wsv pnd_field
                for ( Index i = 0; i< N_se; i++ )
                {
                  pnd_field ( intarr[i]+scat_data_start, p-limits[0],
                              lat-limits[2], lon-limits[4] ) = pnd[i];
                }
            }
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
          os << "Use of size distribution H11 (as requested for\n"
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
          os << "Use of size distribution H11 (as requested for\n"
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
                for ( Index i=0; i<diameter_max.nelem(); i++ ) //loop over number of scattering elements
                {
                    // calculate particle size distribution for H11
                    // [# m^-3 m^-1]
                    dNdD[i] = IWCtopnd_H11 ( diameter_max[i], t_field ( p, lat, lon ) );
                }
                // scale pnds by scale width
                if (diameter_max.nelem() > 1)
                    scale_pnd( pnd, diameter_max, dNdD ); //[# m^-3]
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
          os << "Use of size distribution H13 (as requested for\n"
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
          os << "Use of size distribution H13 (as requested for\n"
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
                for ( Index i=0; i<diameter_max.nelem(); i++ ) //loop over number of scattering elements
                {
                    // calculate particle size distribution for H13
                    // [# m^-3 m^-1]
                    dNdD[i] = IWCtopnd_H13 ( diameter_max[i], t_field ( p, lat, lon ) );
                }
                // scale pnds by scale width
                if (diameter_max.nelem() > 1)
                    scale_pnd( pnd, diameter_max, dNdD ); //[# m^-3]
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
          os << "Use of size distribution H13Shape (as requested for\n"
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
          os << "Use of size distribution H13Shape (as requested for\n"
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
          os << "Use of size distribution H13Shape (as requested for\n"
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
        os << "The *H13Shape* parametrization is only valid for particles with\n"
           << "a maximum diameter >= to 77 um, but the smallest value of\n"
           << "*diameter_max* in this simulation is " << diameter_max[0] << " um.\n"
           << "Set a new *diameter_max_grid*!\n";
        throw runtime_error(os.str());
    }

    if (Rarea_input.nelem()==1)
    { 
        out1 << "WARNING! Only one unique area ratio is used.\n"
             << "Using parametrization *H13* will generate the same results "
             << "as *H13shape*\n"
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
                for ( Index i=0; i<diameter_max_input.nelem(); i++ ) //loop over number of scattering elements
                {
                    // calculate particle size distribution for H13Shape
                    // [# m^-3 m^-1]
                    dNdD[i] = IWCtopnd_H13Shape ( diameter_max_input[i], t_field ( p, lat, lon ) );
                    
                    // calculate Area ratio distribution for H13Shape
                    Ar[i] = area_ratioH13 (diameter_max_input[i], t_field (p, lat, lon ) );
                }

                // scale pnds by scale width
                if (diameter_max_input.nelem() > 1)
                    scale_pnd( pnd_temp, diameter_max_input, dNdD ); //[# m^-3]
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
    \param PR_field     precipitation rate field [kg/(m2*s)]
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
          os << "Use of size distribution MP48 (as requested for\n"
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
          os << "Use of size distribution MP48 (as requested for\n"
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
              mass_total = mass.sum();
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
                for ( Index i=0; i<diameter_melted_equivalent.nelem(); i++ )
                {
                    // calculate particle size distribution with MP48
                    // output: [# m^-3 m^-1]
                    //dNdD[i] = PRtopnd_MP48 ( tPR, diameter_melted_equivalent[i]);
                    // too much a hassle to have a separate function. so we do
                    // the calculation directly here.
                    dNdD[i] = N0 * exp(-lambda*diameter_melted_equivalent[i]);
                    //dNdD2[i] = dNdD[i] * vol[i] * rho[i];
                }

                // scale pnds by bin width
                if (diameter_melted_equivalent.nelem() > 1)
                    scale_pnd( pnd, diameter_melted_equivalent, dNdD );
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
          os << "Use of size distribution H98 (as requested for\n"
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
          os << "Use of size distribution H98 (as requested for\n"
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
                for ( Index i=0; i<radius.nelem(); i++ ) //loop over number of scattering elements
                {
                    // calculate particle size distribution for liquid
                    // [# m^-3 m^-1]
                    dNdr[i] = LWCtopnd ( LWC_field ( p,lat,lon ), radius[i] );
                }

                // scale pnds by scale width. output: [# m^-3]
                if (radius.nelem() > 1)
                    scale_pnd( pnd, radius, dNdr );
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



/*! Calculates particle size distribution using MH97 parametrization.
 *  Each diameter of the scattering elements is a node in the distribution.
 *  One call of this function calculates number density for one scattering
 *  element. Implicitly assumes particles of water ice.

    \return dNdD particle number density per diameter interval [#/m3*m]
          
    \param iwc     atmospheric ice water content [kg/m3]
    \param diameter_mass_equivalent  mass equivalent diameter of scattering element [m]
    \param t       atmospheric temperature [K]
    \param perturb flag whether to add noise onto PSD parameters according to
                   their reported error statistics
  
  \author Daniel Kreyling
  \date 2010-12-06

*/
Numeric IWCtopnd_MH97 ( const Numeric iwc,
                        const Numeric diameter_mass_equivalent,
                        const Numeric t,
                        const bool noisy )
{
  // skip calculation if IWC is 0.0
  if ( iwc == 0.0 )
  {
    return 0.0;
  }

  //[kg/m3] -> [g/m3] as used by parameterisation
  Numeric ciwc = iwc*1e3;
  Numeric cdensity = DENSITY_OF_ICE*1e3;

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

  Numeric dNdD1;
  Numeric dNdD2;
  Numeric dNdD;

  // convert m to microns
  Numeric dmass = diameter_mass_equivalent * 1e6;
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
  // for large IWC alpha, hence dNdD1, goes negative. avoid that.
  // towards this limit, particles anyway get larger 100um, i.e.,
  // outside the size region gamma distrib is intended for
  if (alphas100>0.)
  {
    Numeric Ns100 = 6*IWCs100 * pow ( alphas100,5. ) /
                    ( PI*cdensity*gamma_func ( 5. ) );//micron^-5
    Numeric Nm1 = Ns100*dmass*exp ( -alphas100*dmass ); //micron^-4
    dNdD1 = Nm1*1e18; // m^-3 micron^-1
  }
  else
  {
    dNdD1 = 0.;
  }



  //Log normal distribution component

  // for small IWC, IWCtotal==IWC<100 & IWC>100=0.
  // this will give dNdD2=NaN. avoid that by explicitly setting to 0
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
                   sigmal100 * pow ( 1.,3 ) * dmass; //g/m^3/micron^4
      Numeric Nm2 = a1/a2 *
                    exp ( -0.5 * pow ( ( log ( dmass )-mul100 ) /sigmal100,2 ) );
                    //micron^-4
      dNdD2 = Nm2*1e18; // m^-3 micron^-1
    }
    else
    {
      dNdD2 = 0.;
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
    dNdD2 = 0.;
  }



  dNdD = ( dNdD1+dNdD2 ) *1e6; // m^-3 m^-1
  if (isnan(dNdD)) dNdD = 0.0;
  return dNdD;
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

  Numeric N0 = 0.08*1e-2; // [#/cm3/cm] converted to [#/m3/um]
  Numeric lambda = 41.*1e2*pow(PR,-0.21); // [1/cm] converted to [1/m] to fit diameter_melted_equivalent[m]

  Numeric n = N0*exp(-lambda*diameter_melted_equivalent);
  return n;
}



/*! Scaling pnd values by width of size bin. 
 * Bin width is determined from preceeding and following scattering element size.
 * Vector y and x must be equal in size. Vector w holds the weights.
 * Derived from trapezoid integration rule.
         
    \param w weights
    \param x e.g. scattering element radius [m]
    \param y e.g. particle number density per radius interval [#/m3*m]
  
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
    else // for monodisperse pnd=dNdD
    {
      w[0] = y[0];
    }
}



/*! Check sum of pnd vector against total mass density value.
 *  Deviation is calculated and used to adjust the output of vector pnd.
         
	\param pnd   particle number density [#/m3]
	\param xwc   scattering species massdensity [kg/m3]
	\param mass  scattering element mass [kg]
  
  \author Daniel Kreyling
  \date 2010-12-15

*/
void chk_pndsum (Vector& pnd,
                 const Numeric xwc,
                 const Vector& mass,
                 const Index& p,
                 const Index& lat,
                 const Index& lon,
                 const String& partfield_name,
                 const Verbosity& verbosity)

{
  CREATE_OUT2;
  
  if ( xwc == 0.0 )
    {
      // set all pnd values to zero, IF particle mass density/flux is
      // zero at this atmospheric level.
      pnd = 0.0;
      return;
    }

  // set vector x to pnd size
  Vector x ( pnd.nelem(), 0.0 );
  Numeric error;

  //cout << "p = " << p << ", pnd.nelem:" << pnd.nelem() << ", xwc: " << xwc << "\n";
  for ( Index i = 0; i<pnd.nelem(); i++ )
  {
    // convert from particles/m^3 to kg/m^3
    x[i] = pnd[i]*mass[i];
    /*cout<< "p = " << p << ", i: " << i << "\n"
        << "pnd[i]: " << pnd[i] << ", density[i]: " << density[i] << ", vol[i]: " << vol[i] << "\n"
        << "x[i]: " << x[i] << "\n";*/
  }

  if ( x.sum() == 0.0 )
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
          out1<< "WARNING: in WSM chk_pndsum in pnd_fieldSetup!\n" 
              << "The deviation of the sum of nodes in the particle size distribution\n"
              << "to the initial input mass density (IWC/LWC) is larger than 10%!\n"
              << "The deviation of: "<< error-1.0<<" occured in the atmospheric level: "
              << "p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<".\n";
        }
    }

  out2 << "PND scaling factor in atm. level "
       << "(p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<"): "<< error <<"\n";
}



/*! Splitting scat_species string and parse type of scattering species field

  \param  partfield_name name of atmospheric scattering species field
  \param  part_string    scattering species tag from *scat_species*
  \param  delim          delimiter string of *scat_species* elements

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

  // split scat_species string at delim and write to ArrayOfString
  part_string.split ( strarr, delim );

  //first entry is scattering species field name (e.g. "IWC", "LWC" etc.)
  if (strarr.size()>0 && part_string[0]!=delim[0])
  {
      partfield_name = strarr[0];
  }
  else
  {
      ostringstream os;
      os << "No information on scattering species field name in '"
         << part_string << "'\n";
      throw runtime_error ( os.str() );

  }
}



/*! Splitting scat_species string and parse psd_param
	\param  psd_param   particle size distribution parametrization
	\param  part_string scattering species tag from *scat_species*
  \param  delim       delimiter string of *scat_species* elements
  
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

  // split scat_species string at delim and write to ArrayOfString
  part_string.split ( strarr, delim );

  // second entry is particle size distribution parametrisation  ( e.g."MH97")
  // check, whether we have a second entry
  if (strarr.size()>1)
      psd_param = strarr[1];
  else
      psd_param = "";
}



/*! Splitting scat_species string and parse min and max particle radius
	\param  sizemin     minimum radius of particles to consider
	\param  sizemax     maximum radius of particles to consider
	\param  part_string scattering species tag from *scat_species*
  \param  delim       delimiter string of *scat_species* elements
  
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

  // split scat_species string at delim and write to ArrayOfString
  part_string.split ( strarr, delim );
  
  // convert String for size range, into Numeric
  // 1. third entry is minimum particle size
  if ( strarr.nelem() < 3 || strarr[2] == "*" )
  {
    sizemin = 0.;
  }
  else
  {
    istringstream os1 ( strarr[2] );
    os1 >> sizemin;
  }
  // 2. fourth entry is maximum particle size
  if ( strarr.nelem() < 4 || strarr[3] == "*" )
  {
    sizemax = -1.;
  }
  else
  {
    istringstream os2 ( strarr[3] );
    os2 >> sizemax;
  }
}

