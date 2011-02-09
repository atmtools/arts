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
bool chk_hydromet_field(
			 const Index&  dim,	
			 const Tensor3& hydromet,			 
			 const Vector& p_grid,
			 const Vector& lat_grid,
			 const Vector& lon_grid
		       )
{
  bool x = true;
  
  // check p
  if (hydromet.npages() != p_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *" << p_grid <<"* is not equal to size of *" << hydromet <<"*.";
            throw runtime_error(os.str() );
  }
  
  // check lat
  if(dim >= 2 )
  {
    if (hydromet.nrows() != lat_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *" <<lat_grid <<"* is not equal to size of *" << hydromet <<"*.";
            throw runtime_error(os.str() );
      
    }
  }
  
  // check lon
  if(dim == 3 )
  {
    if (hydromet.ncols() != lon_grid.nelem()) {
    
      ostringstream os;
      os << "The size of *" <<lon_grid <<"* is not equal to size of *" << hydromet <<"*.";
            throw runtime_error(os.str() );
      
    }
  }
  // set x to false if a single value of hydromet_field is unequal zero    
    for (Index j=0; j<hydromet.npages(); j++) {
      for (Index k=0; k<hydromet.nrows(); k++) {
	for (Index l=0; l<hydromet.ncols(); l++) {
	  if (hydromet(j,k,l) != 0.0) x = false;
	}
      }
    }  
  return x;
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
      os << "The number of elments in *"<< scat_data_raw << "*\n"
	 << "and *" << scat_data_meta_array << "* do not match.\n"
	 << "Each SingleScattering file must correspond to one\n"
	 << "to one scattering meta data file.";
	throw runtime_error( os.str());
    }

}

//! Check scattering data meta files
/*!
  FIXME
  
  \param scat_data_meta scattering meta data
  \param scat_data_meta_file Filename of the data to be checked

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
  
  if (!(scat_data_raw.f_grid[0] <= f_grid[0] &&
        last(f_grid) <= 
        last(scat_data_raw.f_grid) ))
    {
      ostringstream os;
      os << "The range of frequency grid in the single scattering"
         << " properties datafile " 
         << scat_data_file << " does not contain all values of"
         << "*f_grid*.";
      throw runtime_error( os.str() );
    }

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
    
  case PTYPE_GENERAL:
    
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
    
  case PTYPE_MACROS_ISO: 
    
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
    
  case PTYPE_HORIZ_AL:
    
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

  case PTYPE_SPHERICAL:
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

/*! barometric heightformula for NONE isothermal atmosphere
    \return Numeric p1.
          
    \param Numeric p,
    \param Numeric dh,
    \param Numeric T0
  
  \author Daniel Kreyling
  \date 2011-01-20
*/
Numeric barometric_heightformula ( //output is p1
				   //input
				   const Numeric p,
				   const Numeric dh,
				   const Numeric T0
				 )
{
				  
  //barometirc height formula
  Numeric M = 0.02896; //mean molar mass kg mol^-1
  Numeric g = 9.807; //earth accelaration kg m s^-1
  Numeric R = 8.314; //universal gas constant J K^−1 mol^−1
  Numeric a = 9.79e-3; //dry adiabatic temp.gradient K m^-1
    
  // calculation
  Numeric p1 = p * pow((1-a*dh/T0),(M*g)/(R*a));
  
  return p1;
  
}

/*! Calculates dN/dDm on one atmospheric grid point per particle type using MH97

    \return Numeric dN.
          
    \param Numeric iwc,
    \param Numeric dm,
    \param Numeric t,
    \param Numeric density
  
  \author Daniel Kreyling
  \date 2010-12-06

*/
Numeric IWCtopnd_MH97 (	const Numeric iwc,
			Numeric dm,
			const Numeric t,
			const Numeric density)
{
	      Numeric dN;
	      // convert m to microns
	      dm *= 1e6;
	      //convert T from Kelvin to Celsius
	      Numeric T = t-273.15;
	      //split IWC in IWCs100 and IWCl100
	      Numeric a=0.252; //g/m^3
	      Numeric b1=0.837;
	      Numeric IWC0=1; //g/m^3
	      Numeric IWCs100=min(iwc,a*pow((iwc/IWC0),(b1)));
	      Numeric IWCl100=iwc-IWCs100;
	      
	      if (dm <= 100) 
		  {
		    //Gamma distribution component

		    Numeric b2=-4.99*1e-3; //micron^-1
		    Numeric m=0.0494; //micron^-1
		    Numeric alfas100=(b2-m*log10(IWCs100/IWC0)); //miron^-1
		    Numeric Ns100=6*IWCs100*(pow(alfas100,(5.)))/(PI*density*gamma((Numeric)5.));//micron^-5
		    Numeric Nm1=Ns100*dm*exp(-alfas100*dm); //micron^-4
		    dN = Nm1*1e18; // m^-3 micron^-1
		     
		    
		   }
	      else {
		    //Log normal distribution component

		    Numeric aamu=5.20;
		    Numeric bamu=0.0013;
		    Numeric abmu=0.026;
		    Numeric bbmu=-1.2*1e-3;
		    Numeric amu=aamu+(bamu*T); 
		    Numeric bmu=abmu+(bbmu*T);
		    Numeric mul100=amu+(bmu*log10(IWCl100/IWC0));
		    
		    Numeric aasigma=0.47;
		    Numeric basigma=2.1*1e-3;
		    Numeric absigma=0.018;
		    Numeric bbsigma=-2.1*1e-4;
		    Numeric asigma=aasigma+(basigma*T);
		    Numeric bsigma=absigma+(bbsigma*T);
		    Numeric sigmal100=asigma+(bsigma*log10(IWCl100/IWC0));
		    
		    Numeric D0=1.0; //micron
		    Numeric a1=6*IWCl100; //g/m^3
		    Numeric a2=(pow(PI,(3./2.)))*density*sqrt(2)*exp(3*mul100+(9./2.)*pow(sigmal100,2))*sigmal100*(pow(D0,3))*dm; //g/m^3/micron^4
		    Numeric Nm2=(a1/a2)*exp(-(1./2.)*pow(((log(dm/D0)-mul100)/sigmal100),2)); //micron^-4
		    dN = Nm2*1e18; // m^-3 micron^-1
		     
		       
		  }
	      if (isnan(dN)) dN = 0.0;
	      return dN;
}


/*! Calculates dN/dr on one atmospheric grid point per liquid water particle type. 

	\return Numeric n.
         
    \param const Numeric lwc,
    \param const Numeric r
  
  \author Daniel Kreyling
  \date 2010-12-16

*/
Numeric LWCtopnd (const Numeric lwc, //[g/m^3]
		  //const Numeric density,
		  const Numeric r // [m]
		  )
{ 	
	Numeric rc = 4.7; //micron
	Numeric alpha = 5.0;
	Numeric gam = 1.05;
	
	Numeric B=(alpha/gam)/pow(rc,gam); 
	Numeric A=(((3*lwc*gam*pow(B,((alpha+4)/gam)))/4)/PI)/gamma((alpha+4)/gam);
	Numeric n=A*(pow(r*1e6,alpha)*exp(-B*pow(r*1e6,gam))); //
	n *= 1e18; // [# m^-3 m^-1]
	//out0<<A;
	//out0<<"\n";
	//out0<<n;
	//out0<<"\n";
	
	if (isnan(n)) n = 0.0;
	return n;
}

/*! Calculates dN/dr on one atmospheric grid point per liquid water particle type. 
         
  \param const Numeric r
  
  \author Daniel Kreyling
  \date 2010-12-16

*/
Numeric LWCtopnd2 (//const Numeric vol, //[g/m^3]
		   //const Numeric density,
		   const Numeric r // [m]
		  )
{ 	
	Numeric rc = 4.7; //micron
	Numeric alpha = 5.0;
	Numeric gam = 1.05;
	
	Numeric B=(alpha/gam)/pow(rc,gam); 
	Numeric A=gam*pow(B,((alpha+1)/gam))/gamma((alpha+1)/gam);
	Numeric n=A*(pow(r*1e6,alpha)*exp(-B*pow(r*1e6,gam)))*1e6;
	// [# m^-3 m^-1]
	
	
	//out0<<A;
	//out0<<"\n";
	//out0<<n;
	//out0<<"\n";
	
	if (isnan(n)) n = 0.0;
	return n;
}


/*! does the trapezoid integration for the function y = f(x). 
 * Vector y and x must be equal in size. Vector w holds the weights.
         
    \param Vector& w,
    \param const Vector& x,
    \param const Vector& y
  
  \author Daniel Kreyling
  \date 2010-12-15

*/
void trapezoid_integrate(  Vector& w,
			   const Vector& x,
			   const Vector& y) 
{
    // check if vectors have same size
    if (x.nelem() != y.nelem()) 
    {
        throw std::logic_error("x and y must be the same size");
    }
    
    // start trapezoid integration  
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
		else { // all vallues in between
		  w[i] = 0.5*(x[i+1]-x[i-1])*y[i]; // m^-3   
		}
	}	

}

/*! Check sum of vector pnd against total XWC.
 * relative error is calculated and used to adjust the output of vector pnd.
         
	\param	 Vector& pnd,
	\param	 const Numeric iwc,
	\param	 const Vector density,
	\param	 const Vector vol
  
  \author Daniel Kreyling
  \date 2010-12-15

*/
void chk_pndsum (Vector& pnd,
		 const Numeric xwc,
		 const Vector density,
		 const Vector vol)
{
    // set vector x to pnd size
    Vector x (pnd.nelem(), 0.0);
    Numeric err; //relerr, abserr;
    
    
    for (Index i = 0; i<pnd.nelem(); i++)
    {
	// convert from particles/m^3 to g/m^3 
	x[i] = pnd[i]*density[i]*vol[i]; 
	//out0<<x[i];
	//out0<<"\n";
	//out0<< pnd[i];
	//out0<< "\n";
    }
    if (xwc == 0.0) 
    {
      // set error and all pnd values to zero
      err = 0.0;
      pnd[Range(joker)] = 0.0;
    }
    else {
      err = xwc/x.sum();
      //abserr = xwc-x.sum();
      //relerr = abs(abserr/xwc);
      //out0<<abserr;
      //out0<<"\n";
      //out0<<err;
      //out0<<"\n";
    }
    // correct all pnd values with error 
    if (err > 1.05 && err < 0.95)
	{
	  pnd *= err;
	}
    //out0<<pnd;
    //out0<<"\n";
    
}
