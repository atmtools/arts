/* Copyright (C) 2002-2008 Claudia Emde <claudia.emde@dlr.de>
                      
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

#include "arts.h"
#include "messages.h"
#include "make_vector.h"
#include "logic.h"
#include "ppath.h"
#include "physics_funcs.h"
#include "math_funcs.h"
#include "check_input.h"


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
void chk_if_pnd_zero_p(
                       const Index& i_p,
                       const GriddedField3& pnd_field_raw,
                       const String& pnd_field_file
                       )
  
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
void chk_if_pnd_zero_lat(
                       const Index& i_lat,
                       const GriddedField3& pnd_field_raw,
                       const String& pnd_field_file
                       )
  
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
void chk_if_pnd_zero_lon(
                       const Index& i_lon,
                       const GriddedField3& pnd_field_raw,
                       const String& pnd_field_file
                       )
  
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
  has the right atmospheric dimension and whether the cloudbox includes
  all points where the particle number density is non-zero. 

  \param pnd_field_raw   pnd field data
  \param pnd_field_file  pnd field filename
  \param atmosphere_dim  Atmospheric dimension
  \param p_grid          Pressure grid
  \param lat_grid        Latitude grid
  \param lon_grid        Longitude grid
  \param cloudbox_limits Cloudbox limits

  \author Claudia Emde
  \date   2005-04-05
*/ 
void chk_pnd_data(
                  const GriddedField3& pnd_field_raw,
                  const String& pnd_field_file,
                  const Index& atmosphere_dim,
                  ConstVectorView p_grid,
                  ConstVectorView lat_grid,
                  ConstVectorView lon_grid,
                  const ArrayOfIndex& cloudbox_limits
                  )
{
  const ConstVectorView pfr_p_grid = pnd_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView pfr_lat_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView pfr_lon_grid = pnd_field_raw.get_numeric_grid(GFIELD3_LON_GRID);

  // The consistency of the dimensions is checked in the reading routine. 
  // Here we have to check whether the atmospheric dimension is correct and whether 
  // the particle number density is 0 on the cloudbox boundary and outside the cloudbox.
  
  out3 << "Check particle number density file " << pnd_field_file 
       << "\n"; 
 
  Index i_p;
 
  // Lower pressure limit
  for (i_p = 0; pfr_p_grid[i_p] > p_grid[cloudbox_limits[0]]; i_p++)
    { chk_if_pnd_zero_p(i_p, pnd_field_raw, pnd_field_file); }
  // The first point inside the cloudbox also needs to be zero !!
  //chk_if_pnd_zero_p(i_p, pnd_field_raw, pnd_field_file);
  
  //Upper pressure limit 
  for (i_p = pfr_p_grid.nelem()-1;
       pfr_p_grid[i_p] < p_grid[cloudbox_limits[1]]; i_p--)
    { chk_if_pnd_zero_p(i_p, pnd_field_raw, pnd_field_file); }
  //chk_if_pnd_zero_p(i_p, pnd_field_raw, pnd_field_file);
  
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

      Index i_lat;
      Index i_lon;

      // Lower latitude limit
      for (i_lat = 0; pfr_lat_grid[i_lat] > 
                      lat_grid[cloudbox_limits[2]]; i_lat++)
        { chk_if_pnd_zero_lat(i_lat, pnd_field_raw, pnd_field_file); }

      // The first point inside the cloudbox also needs to be zero !!
      // chk_if_pnd_zero_lat(i_lat+1, pnd_field_raw, pnd_field_file);

      //Upper latitude limit 
      for (i_lat = pfr_lat_grid.nelem()-1;
           pfr_lat_grid[i_lat] < lat_grid[cloudbox_limits[3]]; 
           i_lat--)
        { chk_if_pnd_zero_lat(i_lat, pnd_field_raw, pnd_field_file); }
      //chk_if_pnd_zero_lat(i_lat-1, pnd_field_raw, pnd_field_file);
      
      // Lower longitude limit
      for (i_lon = 0; pfr_lon_grid[i_lon] > 
           lon_grid[cloudbox_limits[4]]; i_lon++)
        { chk_if_pnd_zero_lon(i_lon, pnd_field_raw, pnd_field_file); }
      // The first point inside the cloudbox also needs to be zero !!
      // chk_if_pnd_zero_lon(i_lon+1, pnd_field_raw, pnd_field_file);
      
      //Upper longitude limit 
      for (i_lon = pfr_lon_grid.nelem()-1;
           pfr_lon_grid[i_lon] < lon_grid[cloudbox_limits[5]]; 
           i_lon--)
        { chk_if_pnd_zero_lon(i_lon, pnd_field_raw, pnd_field_file); }
      //chk_if_pnd_zero_lon(i_lon-1, pnd_field_raw, pnd_field_file);
    } 
  
  out3 << "Particle number density data is o.k. \n";
  
}

//! Check particle number density files (pnd_field_raw)
/*!
  
  \param pnd_field_raw   pnd field raw data (array for all particle types)
  \param pnd_field_file  pnd field filename
  \param atmosphere_dim  Atmospheric dimension
  \param p_grid          Pressure grid
  \param lat_grid        Latitude grid
  \param lon_grid        Longitude grid
  \param cloudbox_limits Cloudbox limits
 
  \author Claudia Emde
  \date   2005-04-05
*/ 
void chk_pnd_raw_data(
                      const ArrayOfGriddedField3& pnd_field_raw,
                      const String& pnd_field_file,
                      const Index& atmosphere_dim,
                      ConstVectorView p_grid,
                      ConstVectorView lat_grid,
                      ConstVectorView lon_grid,
                      const ArrayOfIndex& cloudbox_limits
                      )
{
  for( Index i = 0; i < pnd_field_raw.nelem(); i++)
    {
      out3 << "Element in pnd_field_raw_file:" << i << "\n";
      chk_pnd_data(pnd_field_raw[i],
                   pnd_field_file, atmosphere_dim,
                   p_grid, lat_grid, lon_grid, cloudbox_limits);
    }
}

//! Check scattering data general
/*!
  FIXME
  
  \param scat_data_raw Array of single scattering data
  \param scat_data_meta_array Array of scattering meta data

  \author Daniel Kreyling
  \date 2010-12-02
*/

void chk_scattering_data(
			 const ArrayOfSingleScatteringData& scat_data_raw,
			 const ArrayOfScatteringMetaData& scat_data_meta_array
			 )
{
  if (scat_data_raw.nelem() != scat_data_meta_array.nelem())
    {
      ostringstream os;
      os << "The number of elments in *scat_data_raw*\n"
	 << "and *scat_data_meta_array* do not match.\n"
	 << "Each SingleScattering file must correspond\n"
	 << "to one ScatteringMeta data file.";
	throw runtime_error( os.str());
    }

}

//! Check scattering data meta files
/*!
  FIXME
  
  \param scat_data_meta scattering meta data
  \param scat_data_meta_file filename of the data to be checked

  \author Daniel Kreyling
  \date 2010-12-02
*/

void chk_scattering_meta_data(
			      const ScatteringMetaData& scat_data_meta,
			      const String& scat_data_meta_file
			      )
{
  out2 << "  Check scattering meta properties file "<< scat_data_meta_file 
       << "\n";
  
       if  (scat_data_meta.type != "Ice" && scat_data_meta.type != "Water" && scat_data_meta.type != "Aerosol")
       {
	  ostringstream os; 
	  os << "Type in " << scat_data_meta_file << " must be 'Ice', 'Water' or 'Aerosol'\n";     
	  throw runtime_error( os.str() );
	}
    //(more) checks need to be included
}


//! Check single scattering data files
/*!
  This function checks, whether a datafile containing the single scattering 
  properties of a particle type includes the required frequencies and 
  temperatures and whether the angular grids are defined correctly.
  Furthermore it checks the self consistency of the data by checking the 
  dimensions of pha_mat, ext_mat and abs_vec depening on the ptype case. 
  
  \param scat_data_raw Single scattering data
  \param scat_data_file Filename of the data to be checked.
  \param f_grid        Frequency grid
  
  \author Claudia Emde
  \date   2005-04-04
*/
void chk_single_scattering_data(
                                const SingleScatteringData& scat_data_raw,
                                const String& scat_data_file,
                                ConstVectorView f_grid
                                )
{
  out2 << "  Check single scattering properties file "<< scat_data_file 
       << "\n";

  if (scat_data_raw.ptype != 10 && 
      scat_data_raw.ptype != 20 &&
      scat_data_raw.ptype != 30)
    {
      ostringstream os; 
      os << "Ptype value in file" << scat_data_file << "is wrong."
         << "It must be \n"
         << "10 - arbitrary oriented particles \n"
         << "20 - randomly oriented particles or \n"
         << "30 - horizontally aligned particles.\n";
      throw runtime_error( os.str() );
    }
    
    chk_interpolation_grids("scat_data_raw.f_grid to f_grid",
			    scat_data_raw.f_grid,
			    f_grid);
  
/*  if (!(scat_data_raw.f_grid[0] <= f_grid[0] &&
        last(f_grid) <= 
        last(scat_data_raw.f_grid) ))
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

  if (!(0. < scat_data_raw.T_grid[0] && last(scat_data_raw.T_grid) < 1001.))
    {
      ostringstream os;
      os << "The temperature values in " <<  scat_data_file 
         << " are negative or very large. Check whether you have used the "
         << "right unit [Kelvin].";
      throw runtime_error( os.str() );
    }
  
  if (scat_data_raw.za_grid[0] != 0.)
    {
      ostringstream os;
      os << "The first value of the za grid in the single" 
         << " scattering properties datafile " 
         << scat_data_file << " must be 0.";
        throw runtime_error( os.str() );
    } 

  if (last(scat_data_raw.za_grid) != 180.)
    {
      ostringstream os;
      os << "The last value of the za grid in the single"
         << " scattering properties datafile " 
         << scat_data_file << " must be 180.";
      throw runtime_error( os.str() );
    } 
  
  if (scat_data_raw.ptype == 10 && scat_data_raw.aa_grid[0] != -180.)
     {
       ostringstream os;
       os << "For ptype = 10 (general orientation) the first value"
          << " of the aa grid in the single scattering"
          << " properties datafile " 
          << scat_data_file << "must be -180.";
         throw runtime_error( os.str() );
     } 
  
  if (scat_data_raw.ptype == 30 && scat_data_raw.aa_grid[0] != 0.)
    {
      ostringstream os;
      os << "For ptype = 30 (horizontal orientation)"
         << " the first value"
         << " of the aa grid in the single scattering"
         << " properties datafile " 
         << scat_data_file << "must be 0.";
        throw runtime_error( os.str() );
    }   
  
  if (scat_data_raw.ptype == 30 && last(scat_data_raw.aa_grid) != 180.)
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
  
  switch (scat_data_raw.ptype){
    
  case PARTICLE_TYPE_GENERAL:
    
    out2 << "  Datafile is for arbitrarily orientated particles. \n";
    
    chk_size(os_pha_mat.str(), scat_data_raw.pha_mat_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             scat_data_raw.za_grid.nelem(), scat_data_raw.aa_grid.nelem(),
             scat_data_raw.za_grid.nelem(), scat_data_raw.aa_grid.nelem(), 
              16); 
    
    chk_size(os_ext_mat.str(), scat_data_raw.ext_mat_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             scat_data_raw.za_grid.nelem(), scat_data_raw.aa_grid.nelem(),
             7);
    
    chk_size(os_abs_vec.str(), scat_data_raw.abs_vec_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             scat_data_raw.za_grid.nelem(), scat_data_raw.aa_grid.nelem(),
             4);
    break;
    
  case PARTICLE_TYPE_MACROS_ISO: 
    
    out2 << "  Datafile is for randomly oriented particles, i.e., "
         << "macroscopically isotropic and mirror-symmetric scattering "
         << "media. \n";
    
    chk_size(os_pha_mat.str(), scat_data_raw.pha_mat_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             scat_data_raw.za_grid.nelem(), 1, 1, 1, 6);
    
    chk_size(os_ext_mat.str(), scat_data_raw.ext_mat_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             1, 1, 1);
    
    chk_size(os_abs_vec.str(), scat_data_raw.abs_vec_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             1, 1, 1);
    break; 
    
  case PARTICLE_TYPE_HORIZ_AL:
    
    out2 << "  Datafile is for horizontally aligned particles. \n"; 
    
    chk_size(os_pha_mat.str(), scat_data_raw.pha_mat_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             scat_data_raw.za_grid.nelem(), scat_data_raw.aa_grid.nelem(),
             scat_data_raw.za_grid.nelem()/2+1, 1, 
             16); 

    chk_size(os_ext_mat.str(), scat_data_raw.ext_mat_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             scat_data_raw.za_grid.nelem()/2+1, 1, 
             3);
    
    chk_size(os_abs_vec.str(), scat_data_raw.abs_vec_data,
             scat_data_raw.f_grid.nelem(), scat_data_raw.T_grid.nelem(),
             scat_data_raw.za_grid.nelem()/2+1, 1, 
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
bool is_gp_inside_cloudbox(const GridPos& gp_p,
                           const GridPos& gp_lat,
                           const GridPos& gp_lon,
                           const ArrayOfIndex& cloudbox_limits,
                           const bool include_boundaries)
                        
{
  bool result=true;
  // Pressure dimension
  double ipos = fractional_gp( gp_p );
  if (include_boundaries){
    if( ipos < double( cloudbox_limits[0] )  ||
        ipos > double( cloudbox_limits[1] ) )
      { result=false; }
    
    else {
      // Latitude dimension
      ipos = fractional_gp( gp_lat );
      if( ipos < double( cloudbox_limits[2] )  || 
          ipos > double( cloudbox_limits[3] ) )
        { result=false; }
      
      else
        {
          // Longitude dimension
          ipos = fractional_gp( gp_lon );
          if( ipos < double( cloudbox_limits[4] )  || 
              ipos > double( cloudbox_limits[5] ) )
            { result=false; } 
        }
    }
  }
  else
    {
      if( ipos <= double( cloudbox_limits[0] )  ||
          ipos >= double( cloudbox_limits[1] ) )
        { result=false; }
      
      else {
        // Latitude dimension
        ipos = fractional_gp( gp_lat );
        if( ipos <= double( cloudbox_limits[2] )  || 
            ipos >= double( cloudbox_limits[3] ) )
          { result=false; }
        
        else
          {
            // Longitude dimension
            ipos = fractional_gp( gp_lon );
            if( ipos <= double( cloudbox_limits[4] )  || 
                ipos >= double( cloudbox_limits[5] ) )
              { result=false; }           
          }
      }
    }
  return result;
  
}


/*! Checks, whether the last point of a propagation path 
  is inside the cloudbox.

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
  const Index np=ppath_step.np;
  
  return is_gp_inside_cloudbox(ppath_step.gp_p[np-1],ppath_step.gp_lat[np-1],
                               ppath_step.gp_lon[np-1],cloudbox_limits,include_boundaries);
  
}

/*! barometric heightformula for isothermal earth atmosphere 
    \return p1 pressure in displacement level [Pa]
          
    \param p atmospheric pressure at starting level [Pa]
    \param dh vertical displacement to starting pressure level [m]
    
  
  \author Daniel Kreyling
  \date 2011-01-20
*/
Numeric barometric_heightformula ( //output is p1
				   //input
				   const Numeric& p, 
				   const Numeric& dh
				   )

{
 /* taken from: Seite „Barometrische Höhenformel“. In: Wikipedia, 
 * Die freie Enzyklopädie. Bearbeitungsstand: 3. April 2011, 20:28 UTC.
 * URL: http://de.wikipedia.org/w/index.php?title=Barometrische_H%C3%B6henformel&oldid=87257486 
 * (Abgerufen: 15. April 2011, 15:41 UTC) 
 */

  
  //barometirc height formula
  Numeric M = 0.02896; //mean molar mass of air [kg mol^-1]
  Numeric g = 9.807; //earth acceleration [kg m s^-1]
  Numeric R = 8.314; //universal gas constant [J K^−1 mol^−1]
  Numeric T = 253; //median tropospheric reference temperature [K]
    
  // calculation
  Numeric p1 = p * exp(-(-dh)/(R*T/(M*g)));
  
  return p1;
  
}

/*! Calculates particle size distribution using MH97 parametrization.
 *  Each diameter of the scattering particles is a node in the ditribution.
 *  One call of this function, calculates one particle number density.  

    \return dN particle number density per diameter interval [#/m3*m]
          
    \param iwc atmospheric ice water content [g/m3]
    \param dm melted diameter of scattering particle [m]
    \param t atmospheric temperature [K]
    \param density of the scattering particle [g/m3]
  
  \author Daniel Kreyling
  \date 2010-12-06

*/
Numeric IWCtopnd_MH97 (	const Numeric iwc,
			const Numeric dm,
			const Numeric t,
			const Numeric density)
{
  // skip calculation if IWC is 0.0
  if ( iwc == 0.0 )
  {
    return 0.0;
  }
  Numeric dN1;
  Numeric dN2;
  Numeric dN;
  // convert m to microns
  Numeric Dm = dm * 1e6;
  //convert T from Kelvin to Celsius
  Numeric T = t-273.15;
  //split IWC in IWCs100 and IWCl100
  Numeric a=0.252; //g/m^3
  Numeric b1=0.837;
  Numeric IWC0=1; //g/m^3
  Numeric IWCs100=min ( iwc,a*pow ( iwc/IWC0,b1 ) );
  Numeric IWCl100=iwc-IWCs100;


  //Gamma distribution component

  Numeric b2=-4.99*1e-3; //micron^-1
  Numeric m=0.0494; //micron^-1
  Numeric alfas100=b2-m*log10 ( IWCs100/IWC0 ); //miron^-1
  Numeric Ns100=6*IWCs100*pow ( alfas100,5. ) / ( PI*density*gamma_func ( 5. ) );//micron^-5
  Numeric Nm1=Ns100*Dm*exp ( -alfas100*Dm ); //micron^-4
  dN1 = Nm1*1e18; // m^-3 micron^-1



  //Log normal distribution component

  Numeric aamu=5.20;
  Numeric bamu=0.0013;
  Numeric abmu=0.026;
  Numeric bbmu=-1.2*1e-3;
  Numeric amu=aamu+bamu*T;
  Numeric bmu=abmu+bbmu*T;
  Numeric mul100=amu+bmu*log10 ( IWCl100/IWC0 );

  Numeric aasigma=0.47;
  Numeric basigma=2.1*1e-3;
  Numeric absigma=0.018;
  Numeric bbsigma=-2.1*1e-4;
  Numeric asigma=aasigma+basigma*T;
  Numeric bsigma=absigma+bbsigma*T;
  Numeric sigmal100=asigma+bsigma*log10 ( IWCl100/IWC0 );

  Numeric D0=1.0; //micron
  Numeric a1=6*IWCl100; //g/m^3
  Numeric a2=pow ( PI,3./2. ) *density*sqrt ( 2 ) *exp ( 3*mul100+9./2.*pow ( sigmal100,2 ) ) *sigmal100*pow ( D0,3 ) *Dm; //g/m^3/micron^4
  Numeric Nm2=a1/a2*exp ( -0.5*pow ( ( log ( Dm/D0 )-mul100 ) /sigmal100,2 ) ); //micron^-4
  dN2 = Nm2*1e18; // m^-3 micron^-1



  dN = ( dN1+dN2 ) *1e6; // m^-3 m^-1
  if (isnan(dN)) dN = 0.0;
  return dN;
}


/*! Calculates particle size distribution for liquid particles using gamma parametrization.
 *  Each radius of the scattering particles is a node in the ditribution.
 *  One call of this function, calculates one particle number density. 

	\return n particle number density per radius interval [#/m3*m]
         
	\param lwc atmospheric liquid water content [g/m3]
	\param r radius of scattering particle [m]
  
  \author Daniel Kreyling
  \date 2010-12-16

*/
Numeric LWCtopnd (const Numeric lwc, //[g/m^3]
		  //const Numeric density,
		  const Numeric r // [m]
		  )
{ 
  // skip calculation if LWC is 0.0
  if ( lwc == 0.0 )
  {
    return 0.0;
  }
	Numeric rc = 4.7; //[micron]
	Numeric alpha = 5.0;
	Numeric gam = 1.05;
	
	Numeric B=(alpha/gam)/pow(rc,gam); 
	// factor 1e12 is density of water [1 g/cm^3] in units of [g/micron^3]
	Numeric A=(((3*lwc*gam*pow(B,((alpha+4)/gam))*1e12)/4)/PI)/gamma_func((alpha+4)/gam); 
	Numeric n=A*(pow(r*1e6,alpha)*exp(-B*pow(r*1e6,gam)))*1e6; //
	// n in [# m^-3 m^-1]
	
	if (isnan(n)) n = 0.0;
	return n;
}

// ONLY FOR TESTING PURPOSES
Numeric LWCtopnd2 (//const Numeric vol, //[g/m^3]
		   //const Numeric density,
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


/*! Scaling pnd values by width of size bin. 
 * Bin width is detemined from preceding and following particle size.
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
    
    // calc. weights (derived from trapezoid integration)  
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
		else { // all values inbetween
		  w[i] = 0.5*(x[i+1]-x[i-1])*y[i]; // m^-3   
		}
	}	

}

/*! Check sum of pnd vector against total mass density value.
 *  Deviation is calculated and used to adjust the output of vector pnd.
         
	\param pnd particle number density [g/m3]
	\param xwc atmospheric massdensity [g/m3]
	\param density scattering particle density [g/m3]
	\param vol scattering particle volume [m3]
  
  \author Daniel Kreyling
  \date 2010-12-15

*/
void chk_pndsum (Vector& pnd,
		 const Numeric xwc,
		 const Vector& density,
		 const Vector& vol,
		 const Index& p,
		 const Index& lat,
		 const Index& lon
		)

{
  // set vector x to pnd size
  Vector x ( pnd.nelem(), 0.0 );
  Numeric error;

  for ( Index i = 0; i<pnd.nelem(); i++ )
  {
    // convert from particles/m^3 to g/m^3
    x[i] = pnd[i]*density[i]*vol[i];
    //out0<<x[i]<<"\n"<< pnd[i]<< "\n";
  }

  if ( x.sum() == 0.0 )
  {
    // set error and all pnd values to zero, IF there is 
    // no scattering particles at this atmospheric level.
    error = 0.0;
    pnd = 0.0;
  }
  else
  {
    error = xwc/x.sum();
  // correct all pnd values with error
    pnd *= error;
  // give warning if deviations are more than 10%
    if ( error > 1.10 || error < 0.90 )
    {
      //ostringstream os;
      out1<< "WARNING: in WSM chk_pndsum in pnd_fieldSetup!\n" 
         << "The deviation of the sum of nodes in the particle size distribution\n"
	 << "to the initial input mass density (IWC/LWC) is larger than 10%!\n"
	 << "The deviation of: "<< error-1.0<<" occured in the atmospheric level: p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<".\n";
      //cerr<<os;
    }
  }

  out2 << "PND scaling factor in atm. level (p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<"): "<< error <<"\n";
    //cout<<"\npnd_scaled\t"<<pnd<<endl;
    //cout<<"\nPND\t"<<pnd.sum()<<"\nXWC\t"<<xwc<<"\nerror\t"<<error<<endl;
    //cout<<pnd.sum()<<"\n"<<xwc<<"\n"<<error<<endl;
}


// ONLY FOR TESTING PURPOSES
void chk_pndsum2 (Vector& pnd,
		 const Numeric xwc)

{
  // set vector x to pnd size
  Vector x=pnd;
  Numeric error;

  if ( xwc == 0.0 )
  {
    // set error and all pnd values to zero
    error = 0.0;
    pnd = 0.0;
  }
  else
  {
    error = xwc/x.sum();

    // correct all pnd values with error
    if ( error > 1.05 || error < 0.95 )
    {
      pnd *= error
             ;
    }

  }
    //cout<<"\nPND2\t"<<pnd.sum()<<"\nXWC2\t"<<xwc<<"\nerror2\t"<<error<<endl;
    //cout<<pnd.sum()<<"\n"<<xwc<<"\n"<<error<<endl;
    //cout<<"\npnd_scaled2\t"<<pnd<<endl;
}

/*! Splitting part_species string and parse type of massdensity_field

	\param  part_type type of atmospheric scattering particle profile 
	\param  part_string containing infos about scattering particle calculations

  \author Daniel Kreyling
  \date 2011-02-21

*/
void parse_part_type ( //WS Output:
  String& part_type,
  // WS Input:
  const String& part_string )
{

  ArrayOfString strarr;

  // split part_species string at "-" and write to ArrayOfString
  part_string.split ( strarr, "-" );

  //first entry is hydrometeor type (e.g. "IWC", "LWC" etc.)
  part_type = strarr[0];

  if (  part_type != "IWC" && 
	part_type != "Snow" &&
	part_type != "LWC" &&
	part_type != "Rain" )
  {
    ostringstream os;
    os << "First substring in " << part_string << " must be rather 'LWC', 'IWC', 'Rain' or 'Snow'\n"
    <<"Ckeck input in *part_species*!\n";
    throw runtime_error ( os.str() );
  }
}


/*! Splitting part_species string and parse psd_param
	\param psd_param particle size distribution parametrization
	\param part_string containing infos about scattering particle calculations
  
  \author Daniel Kreyling
  \date 2011-02-21

*/
void parse_psd_param ( //WS Output:
  String& psd_param,
  // WS Input:
  const String& part_string )
{

  ArrayOfString strarr;

  // split part_species string at "-" and write to ArrayOfString
  part_string.split ( strarr, "-" );
  // second entry is particle size distribution parametrisation  ( e.g."MH97")
  psd_param = strarr[1];

  if ( psd_param != "MH97" && psd_param != "liquid" )
  {
    ostringstream os;
    os <<"The chosen PSD parametrisation in " << part_string << " can not be handeled in the moment.\n"
    <<"Choose either 'MH97' or 'liquid'!\n" ;
    throw runtime_error ( os.str() );
  }
}

/*! Splitting part_species string and parse min and max particle radius
	\param sizemin min scattering particle radius
	\param sizemax max scattering particle radius
	\param part_string containing infos about scattering particle calculations
  
  \author Daniel Kreyling
  \date 2011-02-21

*/
void parse_part_size ( //WS Output:
  Numeric& sizemin,
  Numeric& sizemax,
  // WS Input:
  const String& part_string )
{

  ArrayOfString strarr;

  // split part_species string at "-" and write to ArrayOfString
  part_string.split ( strarr, "-" );
  
  // convert String for size range, into Numeric
  // 1. third entry is minimum particle radius
  if ( strarr.nelem() < 3 || strarr[2] == "*" )
  {
    sizemin = 0.;
  }
  else
  {
    istringstream os1 ( strarr[2] );
    os1 >> sizemin;
  }
  // 2. fourth entry is maximum particle radius
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
