/* Copyright (C) 2002-2007 Claudia Emde <claudia.emde@dlr.de>
                      
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
  for (Index i = 0; i < pnd_field_raw.lat_grid.nelem(); i++ ) 
    {
      for (Index j = 0; j < pnd_field_raw.lon_grid.nelem(); j++ )
        {
          if ( pnd_field_raw.data(i_p, i, j) != 0. )
            {
              ostringstream os;
              os << "Warning: \n"
                 << "The particle number density field contained in the file '"
                 << pnd_field_file << "'\nis non-zero outside the cloudbox "
                 << "or close the cloudbox boundary at the \n"
                 << "following position:\n"
                 << "pressure = " << pnd_field_raw.p_grid[i_p] 
                 << ", p_index = " << i_p << "\n"
                 << "latitude = " << pnd_field_raw.lat_grid[i] 
                 << ", lat_index = " << i << "\n"
                 << "longitude = " << pnd_field_raw.lon_grid[j] 
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
  for (Index i = 0; i < pnd_field_raw.p_grid.nelem(); i++ ) 
    {
      for (Index j = 0; j < pnd_field_raw.lon_grid.nelem(); j++ )
        {
          if ( pnd_field_raw.data(i, i_lat, j) != 0. )
            {
              ostringstream os;
              os << "Warning: \n" 
                 << "The particle number density field contained in the file '"
                 << pnd_field_file << "'\nis non-zero outside the cloudbox "
                 << "or close the cloudbox boundary at the \n"
                 << "following position:\n"
                 << "pressure = " << pnd_field_raw.p_grid[i] << ", p_index = "
                 << i << "\n"
                 << "latitude = " << pnd_field_raw.lat_grid[i_lat] 
                 << ", lat_index = "<<i_lat<< "\n"
                 << "longitude = " << pnd_field_raw.lon_grid[j] 
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
  for (Index i = 0; i < pnd_field_raw.p_grid.nelem(); i++ ) 
    {
      for (Index j = 0; j < pnd_field_raw.lat_grid.nelem(); j++ )
        {
          if ( pnd_field_raw.data(i, j, i_lon) != 0. )
            {
              ostringstream os;
              os << "Warning: \n" 
                 << "The particle number density field contained in the file '"
                 << pnd_field_file << "'\nis non-zero outside the cloudbox "
                 << "or close the cloudbox boundary at the \n"
                 << "following position:\n"
                 << "pressure = " << pnd_field_raw.p_grid[i] 
                 << ", p_index = " << i << "\n"
                 << "latitude = " << pnd_field_raw.lat_grid[j] 
                 << ", lat_index = " << j << "\n"
                 << "longitude = " << pnd_field_raw.lon_grid[i_lon] 
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
  // The consistency of the dimensions is checked in the reading routine. 
  // Here we have to check whether the atmospheric dimension is correct and whether 
  // the particle number density is 0 on the cloudbox boundary and outside the cloudbox.
  
  out3 << "Check particle number density file " << pnd_field_file 
       << "\n"; 
 
  Index i_p;
 
  // Lower pressure limit
  for (i_p = 0; pnd_field_raw.p_grid[i_p] > p_grid[cloudbox_limits[0]]; i_p++)
    { chk_if_pnd_zero_p(i_p, pnd_field_raw, pnd_field_file); }
  // The first point inside the cloudbox also needs to be zero !!
  //chk_if_pnd_zero_p(i_p, pnd_field_raw, pnd_field_file);
  
  //Upper pressure limit 
  for (i_p = pnd_field_raw.p_grid.nelem()-1;
       pnd_field_raw.p_grid[i_p] < p_grid[cloudbox_limits[1]]; i_p--)
    { chk_if_pnd_zero_p(i_p, pnd_field_raw, pnd_field_file); }
  //chk_if_pnd_zero_p(i_p, pnd_field_raw, pnd_field_file);
  
  if (atmosphere_dim == 1 && (pnd_field_raw.lat_grid.nelem() != 1 
                              || pnd_field_raw.lon_grid.nelem() != 1) )
    {
      ostringstream os; 
      os << "The atmospheric dimension is 1D but the particle "
         << "number density file * " << pnd_field_file 
         << " is for a 3D atmosphere. \n";
      throw runtime_error( os.str() );
    }
      
  
  else if( atmosphere_dim == 3) 
    {
      if(pnd_field_raw.lat_grid.nelem() == 1 
         || pnd_field_raw.lon_grid.nelem() == 1)
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
      for (i_lat = 0; pnd_field_raw.lat_grid[i_lat] > 
                      lat_grid[cloudbox_limits[2]]; i_lat++)
        { chk_if_pnd_zero_lat(i_lat, pnd_field_raw, pnd_field_file); }

      // The first point inside the cloudbox also needs to be zero !!
      // chk_if_pnd_zero_lat(i_lat+1, pnd_field_raw, pnd_field_file);

      //Upper latitude limit 
      for (i_lat = pnd_field_raw.lat_grid.nelem()-1;
           pnd_field_raw.lat_grid[i_lat] < lat_grid[cloudbox_limits[3]]; 
           i_lat--)
        { chk_if_pnd_zero_lat(i_lat, pnd_field_raw, pnd_field_file); }
      //chk_if_pnd_zero_lat(i_lat-1, pnd_field_raw, pnd_field_file);
      
      // Lower longitude limit
      for (i_lon = 0; pnd_field_raw.lon_grid[i_lon] > 
           lon_grid[cloudbox_limits[4]]; i_lon++)
        { chk_if_pnd_zero_lon(i_lon, pnd_field_raw, pnd_field_file); }
      // The first point inside the cloudbox also needs to be zero !!
      // chk_if_pnd_zero_lon(i_lon+1, pnd_field_raw, pnd_field_file);
      
      //Upper longitude limit 
      for (i_lon = pnd_field_raw.lon_grid.nelem()-1;
           pnd_field_raw.lon_grid[i_lon] < lon_grid[cloudbox_limits[5]]; 
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
  out2 << "Check single scattering properties file "<< scat_data_file 
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
    
    out2 << "Datafile is for arbitrarily orientated particles. \n";
    
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
    
    out2 << "Datafile is for randomly oriented particles, i.e., "
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
    
    out2 << "Datafile is for horizontally aligned particles. \n"; 
    
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



//! Interpolation of cloud box intensity field
/* 
   See WSM *iyInterpCloudboxField*.

   \param iy                Out: As the WSV with same name.
   \param scat_i_p          In: As the WSV with same name.
   \param scat_i_lat        In: As the WSV with same name.
   \param scat_i_lon        In: As the WSV with same name.
   \param doit_i_field1D_spectrum In: As the WSV with same name.
   \param rte_gp_p          In: As the WSV with same name.
   \param rte_gp_lat        In: As the WSV with same name.
   \param rte_gp_lon        In: As the WSV with same name.
   \param rte_los           In: As the WSV with same name.
   \param cloudbox_on       In: As the WSV with same name.
   \param cloudbox_limits   In: As the WSV with same name.
   \param atmosphere_dim    In: As the WSV with same name.
   \param stokes_dim        In: As the WSV with same name.
   \param scat_za_grid      In: As the WSV with same name. 
   \param scat_aa_grid      In: As the WSV with same name. 
   \param f_grid            In: As the WSV with same name.
   \param p_grid            In: As the WSV with same name.
   \param interpmeth        Interpolation method. Can be "linear" or 
   "polynomial".
 
   \author Claudia Emde and Patrick Eriksson
   \date 2004-09-29
*/
void iy_interp_cloudbox_field(
                              Matrix&               iy,
                              const Tensor7&        scat_i_p,
                              const Tensor7&        scat_i_lat,
                              const Tensor7&        scat_i_lon,
                              const Tensor4&        doit_i_field1D_spectrum, 
                              const GridPos&        rte_gp_p,
                              const GridPos&        rte_gp_lat,
                              const GridPos&        rte_gp_lon,
                              const Vector&         rte_los,
                              const Index&          cloudbox_on,
                              const ArrayOfIndex&   cloudbox_limits,
                              const Index&          atmosphere_dim,
                              const Index&          stokes_dim,
                              const Vector&         scat_za_grid,
                              const Vector&         scat_aa_grid,
                              const Vector&         f_grid,
                              const String&         interpmeth )
{
  //--- Check input -----------------------------------------------------------
  if( !(atmosphere_dim == 1  ||  atmosphere_dim == 3) )
    throw runtime_error( "The atmospheric dimensionality must be 1 or 3.");
  if( !cloudbox_on )
    throw runtime_error( "The cloud box is not activated and no outgoing "
                         "field can be returned." );
  if ( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error(
       "*cloudbox_limits* is a vector which contains the upper and lower\n"
       "limit of the cloud for all atmospheric dimensions.\n"
       "So its length must be 2 x *atmosphere_dim*" ); 
  if( scat_za_grid.nelem() == 0 )
    throw runtime_error( "The variable *scat_za_grid* is empty. Are dummy "
                         "values from *cloudboxOff used?" );
  if( !( interpmeth == "linear"  ||  interpmeth == "polynomial" ) )
    throw runtime_error( "Unknown interpolation method. Possible choices are "
                         "\"linear\" and \"polynomial\"." );
  if( interpmeth == "polynomial"  &&  atmosphere_dim != 1  )
    throw runtime_error( "Polynomial interpolation method is only available"
                         "for *atmosphere_dim* = 1." );
  // Size of scat_i_p etc is checked below
  //---------------------------------------------------------------------------


  //--- Determine if at border or inside of cloudbox (or outside!)
  //
  // Let us introduce a number coding for cloud box borders.
  // Borders have the same number as position in *cloudbox_limits*.
  // Innside cloud box is coded as 99, and outside as > 100.
  Index  border  = 999;
  //
  //- Check if at any border
  if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[0] ) )
    { border = 0; }
  else if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[1] ) )
    { border = 1; }
  if( atmosphere_dim == 3  &&  border > 100 )
    {
      if( is_gridpos_at_index_i( rte_gp_lat, cloudbox_limits[2] ) )
        { border = 2; }
      else if( is_gridpos_at_index_i( rte_gp_lat, cloudbox_limits[3] ) )
        { border = 3; }
      else if( is_gridpos_at_index_i( rte_gp_lon, cloudbox_limits[4] ) )
        { border = 4; }
      else if( is_gridpos_at_index_i( rte_gp_lon, cloudbox_limits[5] ) )
        { border = 5; }
    }

  //
  //- Check if inside
  if( border > 100 )
    {
      // Assume inside as it is easiest to detect if outside (can be detected
      // check in one dimension at the time)
      bool inside = true;
      Numeric fgp;

      // Check in pressure dimension
      fgp = fractional_gp( rte_gp_p );
      if( fgp < Numeric(cloudbox_limits[0])  || 
          fgp > Numeric(cloudbox_limits[1]) )
        { inside = false; }

      // Check in lat and lon dimensions
     if( atmosphere_dim == 3  &&  inside )
       {
         fgp = fractional_gp( rte_gp_lat );
         if( fgp < Numeric(cloudbox_limits[2])  || 
             fgp > Numeric(cloudbox_limits[3]) )
           { inside = false; }
         fgp = fractional_gp( rte_gp_lon );
         if( fgp < Numeric(cloudbox_limits[4])  || 
             fgp > Numeric(cloudbox_limits[5]) )
           { inside = false; }
       }

     if( inside )
       { border = 99; }
    }

  // If outside, something is wrong
  if( border > 100 )
    {
      
      throw runtime_error( 
                 "Given position has been found to be outside the cloudbox." );
    }

  //- Sizes
  const Index   nf  = f_grid.nelem();
  DEBUG_ONLY (const Index   np  = cloudbox_limits[1] - cloudbox_limits[0] + 1);
  const Index   nza  = scat_za_grid.nelem();

  //- Resize *iy*
  iy.resize( nf, stokes_dim );

  
  // Sensor points inside the cloudbox
  if( border == 99 )
    {
      if (atmosphere_dim == 3)
        {
          throw runtime_error(
                              "3D DOIT calculations are not implemented\n"
                              "for observations from inside the cloudbox.\n"
                              );
        }
      else
        {
          assert(atmosphere_dim == 1);
          
          // *doit_i_field1D_spectra* is normally calculated internally:
          assert( is_size(doit_i_field1D_spectrum, nf, np, nza, stokes_dim) );
          
          out3 << "    Interpolating outgoing field:\n";
          out3 << "       zenith_angle: " << rte_los[0] << "\n";
          out3 << " Sensor inside cloudbox at position:  " << 
            rte_gp_p << "\n";
          
          // Grid position in *scat_za_grid*
          GridPos gp_za;
          gridpos( gp_za, scat_za_grid, rte_los[0] );
          
          // Corresponding interpolation weights
          Vector itw_za(2);
          interpweights( itw_za, gp_za );
          
          // Grid position in *p_grid* (only cloudbox part because 
          // doit_i_field1D_spectra is only defined inside the cloudbox
          GridPos gp_p;
          gp_p = rte_gp_p;
          gp_p.idx = rte_gp_p.idx - cloudbox_limits[0]; 

          Vector itw_p(2);
          interpweights( itw_p, gp_p );

          Vector iy_p(nza);

          if( interpmeth == "linear" )
            {
              for(Index is = 0; is < stokes_dim; is++ )
                {
                  for(Index iv = 0; iv < nf; iv++ )
                    {
                      for (Index i_za = 0; i_za < nza; i_za++)
                        {
                          iy_p[i_za] = interp
                          (itw_p, doit_i_field1D_spectrum(iv, joker, i_za, is),
                             gp_p);
                        }
                      iy(iv,is) = interp( itw_za, iy_p, gp_za);
                    }
                }
            }
          else
            {   
              for(Index is = 0; is < stokes_dim; is++ )
                {
                  for(Index iv = 0; iv < nf; iv++ )
                    {
                      for (Index i_za = 0; i_za < nza; i_za++)
                        {
                          iy_p[i_za] = interp
                          (itw_p, doit_i_field1D_spectrum(iv, joker, i_za, is),
                             gp_p);
                        }
                      iy(iv,is) =  interp_poly( scat_za_grid, iy_p, rte_los[0],
                                                gp_za );
                    }
                }
            }
          
        }
      
    }
  
  // Sensor outside the cloudbox

  // --- 1D ------------------------------------------------------------------
  else if( atmosphere_dim == 1 )
    {
      // Use assert to check *scat_i_p* as this variable should to 99% be
      // calculated internally, and thus make it possible to avoid this check.
      assert( is_size( scat_i_p, nf,2,1,1,scat_za_grid.nelem(),1,stokes_dim ));

      out3 << "    Interpolating outgoing field:\n";
      out3 << "       zenith_angle: " << rte_los[0] << "\n";
      if( border )
        out3 << "       top side\n";
      else
        out3 << "       bottom side\n";
      
      // Grid position in *scat_za_grid*
      GridPos gp;
      gridpos( gp, scat_za_grid, rte_los[0] );

      // Corresponding interpolation weights
      Vector itw(2);
      interpweights( itw, gp );

      if( interpmeth == "linear" )
        {
          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {

                  iy(iv,is) = interp( itw, 
                                    scat_i_p( iv, border, 0, 0, joker, 0, is ),
                                      gp );
                }
            }
        }
      else
        {
          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp_poly( scat_za_grid, 
                       scat_i_p( iv, border, 0, 0, joker, 0, is ) , rte_los[0],
                                            gp );
                }
            }
        }
    }
  
  // --- 3D ------------------------------------------------------------------
  else
    {
      // Use asserts to check *scat_i_XXX* as these variables should to 99% be
      // calculated internally, and thus make it possible to avoid this check.
      assert ( is_size( scat_i_p, nf, 2, scat_i_p.nshelves(), 
                        scat_i_p.nbooks(), scat_za_grid.nelem(), 
                        scat_aa_grid.nelem(), stokes_dim ));

      assert ( is_size( scat_i_lat, nf, scat_i_lat.nvitrines(), 2, 
                        scat_i_p.nbooks(), scat_za_grid.nelem(), 
                        scat_aa_grid.nelem(), stokes_dim ));
      assert ( is_size( scat_i_lon, nf, scat_i_lat.nvitrines(), 
                        scat_i_p.nshelves(), 2, scat_za_grid.nelem(), 
                        scat_aa_grid.nelem(), stokes_dim ));

      out3 << "    Interpolating outgoing field:\n";
      out3 << "       zenith angle : " << rte_los[0] << "\n";
      out3 << "       azimuth angle: " << rte_los[1]+180. << "\n";

      
      // Scattering angle grid positions
      GridPos gp_za, gp_aa;
      gridpos( gp_za, scat_za_grid, rte_los[0] );
      gridpos( gp_aa, scat_aa_grid, rte_los[1]+180. );

      // Interpolation weights (for 4D "red" interpolation)
      Vector   itw(16);

      // Outgoing from pressure surface
      if( border <= 1 )
        {
          // Lat and lon grid positions with respect to cloud box 
          GridPos cb_gp_lat, cb_gp_lon;
          cb_gp_lat      = rte_gp_lat;
          cb_gp_lon      = rte_gp_lon;
          cb_gp_lat.idx -= cloudbox_limits[2];
          cb_gp_lon.idx -= cloudbox_limits[4]; 
          
          interpweights( itw, cb_gp_lat, cb_gp_lon, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                        scat_i_p( iv, border, joker, joker, joker, joker, is ),
                                      cb_gp_lat, cb_gp_lon, gp_za, gp_aa );
                }
            }
        }

      // Outgoing from latitude surface
      else if( border <= 3 )
        {
          // Pressure and lon grid positions with respect to cloud box 
          GridPos cb_gp_p, cb_gp_lon;
          cb_gp_p        = rte_gp_p;
          cb_gp_lon      = rte_gp_lon;
          cb_gp_p.idx   -= cloudbox_limits[0];
          cb_gp_lon.idx -= cloudbox_limits[4]; 
          
          interpweights( itw, cb_gp_p, cb_gp_lon, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                    scat_i_lat( iv, joker, border-2, joker, joker, joker, is ),
                                      cb_gp_p, cb_gp_lon, gp_za, gp_aa );
                }
            }
        }

      // Outgoing from longitude surface
      else
        {
          // Pressure and lat grid positions with respect to cloud box 
          GridPos cb_gp_p, cb_gp_lat;
          cb_gp_p        = rte_gp_p;
          cb_gp_lat      = rte_gp_lat;
          cb_gp_p.idx   -= cloudbox_limits[0]; 
          cb_gp_lat.idx -= cloudbox_limits[2];
          
          interpweights( itw, cb_gp_p, cb_gp_lat, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                     scat_i_lon( iv, joker, joker, border-4, joker, joker, is ),
                                      cb_gp_p, cb_gp_lat, gp_za, gp_aa );
                }
            }
        }
    }
}

