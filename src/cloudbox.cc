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
          if ( pnd_field_raw(i_p, i, j) != 0. )
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
          if ( pnd_field_raw(i, i_lat, j) != 0. )
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
          if ( pnd_field_raw(i, j, i_lon) != 0. )
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

  if (!(0. < scat_data_raw.T_grid[0] && last(scat_data_raw.T_grid) < 321.))
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
