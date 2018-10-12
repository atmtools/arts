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
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "mc_antenna.h"
#include "messages.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rng.h"
#include "sorting.h"



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
  cloudbox removed. done in pnd_fieldCalcFrompnd_field_raw). 

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
  Index nelem=2;

  for ( Index k=0; k<scat_species.nelem(); k++ )
    {
      scat_species[k].split ( strarr, delim );
      if ( strarr.nelem() < nelem )
        {     
          ostringstream os;
          os << "Individual strings in scat_species must contain at least "
             << nelem << " elements,\n"
             << "but entry #" << k << " contains only the following "
             << strarr.nelem() << ":\n" << strarr << "\n";
          throw runtime_error ( os.str() );
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
    os << "The number of elements in in current scat_species'  *scat_data* and "
       << "*scat_meta* do not match.\n"
       << "Each *scat_data* entry must correspond to one entry in *scat_meta*.";
    throw runtime_error( os.str());
  }

}

//! Check scattering data meta
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


//! Check single scattering data
/*!
  This function checks the self consistency of the data by checking the
  dimensions of pha_mat, ext_mat and abs_vec depending on the ptype case.
  It furthermore checks whether the angular grids are defined correctly
  depending on ptype and the sanity of the temperature grid.
  
  \param scat_data_single[in]  Single scattering data of a single scattering element

  \author Claudia Emde
  \date   2005-04-04
*/
void chk_scat_data(const SingleScatteringData& scat_data_single,
                   const Verbosity& verbosity)
{
  CREATE_OUT3;

  assert(scat_data_single.ptype == PTYPE_GENERAL ||
         scat_data_single.ptype == PTYPE_TOTAL_RND ||
         scat_data_single.ptype == PTYPE_AZIMUTH_RND);

  if (scat_data_single.za_grid[0] != 0.)
    {
      ostringstream os;
      os << "The first value of the zenith angle grid in the single" 
         << " scattering properties data must be 0.";
        throw runtime_error( os.str() );
    } 

  if (last(scat_data_single.za_grid) != 180.)
    {
      ostringstream os;
      os << "The last value of the zenith angle grid in the single"
         << " scattering properties data must be 180.";
      throw runtime_error( os.str() );
    } 
  
  if (scat_data_single.ptype == PTYPE_GENERAL && scat_data_single.aa_grid[0] != -180.)
     {
       ostringstream os;
       os << "For ptype = \"general\" the first value"
          << " of the azimuth angle grid in the single scattering"
          << " properties data must be -180.";
         throw runtime_error( os.str() );
     } 
  
  if (scat_data_single.ptype == PTYPE_AZIMUTH_RND && scat_data_single.aa_grid[0] != 0.)
    {
      ostringstream os;
      os << "For ptype = \"azimuthally_random\""
         << " the first value"
         << " of the azimuth angle grid in the single scattering"
         << " properties data must be 0.";
        throw runtime_error( os.str() );
    }   
  
  if (scat_data_single.ptype != PTYPE_TOTAL_RND && last(scat_data_single.aa_grid) != 180.)
    {
      ostringstream os;
      os << "For ptypes = \"azimuthally_random\" and \"general\""
         << " the last value of the azimuth angle grid in the single"
         << " scattering properties data must be 180.";
        throw runtime_error( os.str() );
    }   

  ostringstream os_pha_mat;
  os_pha_mat << "pha_mat ";
  ostringstream os_ext_mat;
  os_ext_mat << "ext_mat ";
  ostringstream os_abs_vec;
  os_abs_vec << "abs_vec ";
  
  switch (scat_data_single.ptype){
    
  case PTYPE_GENERAL:
    
    out3 << "  Data is for arbitrarily orientated particles. \n";
    
    chk_size(os_pha_mat.str(), scat_data_single.pha_mat_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             scat_data_single.za_grid.nelem(), scat_data_single.aa_grid.nelem(),
             scat_data_single.za_grid.nelem(), scat_data_single.aa_grid.nelem(), 
              16); 
    
    chk_size(os_ext_mat.str(), scat_data_single.ext_mat_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             scat_data_single.za_grid.nelem(), scat_data_single.aa_grid.nelem(),
             7);
    
    chk_size(os_abs_vec.str(), scat_data_single.abs_vec_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             scat_data_single.za_grid.nelem(), scat_data_single.aa_grid.nelem(),
             4);
    break;
    
  case PTYPE_TOTAL_RND:
    
    out3 << "  Data is for macroscopically isotropic and mirror-symmetric "
         << "scattering media, i.e. for totally randomly oriented particles "
         << "with at least one plane of symmetry. \n";
    
    chk_size(os_pha_mat.str(), scat_data_single.pha_mat_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             scat_data_single.za_grid.nelem(), 1, 1, 1, 6);
    
    chk_size(os_ext_mat.str(), scat_data_single.ext_mat_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             1, 1, 1);
    
    chk_size(os_abs_vec.str(), scat_data_single.abs_vec_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             1, 1, 1);
    break; 
    
  case PTYPE_AZIMUTH_RND:
    
    out3 << "  Data is for azimuthally randomly oriented particles. \n";
    
    chk_size(os_pha_mat.str(), scat_data_single.pha_mat_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             scat_data_single.za_grid.nelem(), scat_data_single.aa_grid.nelem(),
             scat_data_single.za_grid.nelem(), 1, 
             16); 

    chk_size(os_ext_mat.str(), scat_data_single.ext_mat_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             scat_data_single.za_grid.nelem(), 1, 
             3);
    
    chk_size(os_abs_vec.str(), scat_data_single.abs_vec_data,
             scat_data_single.f_grid.nelem(), scat_data_single.T_grid.nelem(),
             scat_data_single.za_grid.nelem(), 1, 
             2);
    break;

  }

  // Here we only check whether the temperature grid is of the unit K, not 
  // whether it corresponds to the required values in t_field. The second 
  // option is not trivial since here one has to look whether the pnd_field 
  // is non-zero for the corresponding temperature. This check is done in the 
  // functions where the multiplication with the particle number density is 
  // done. 
  if (scat_data_single.T_grid[0]<0. || last(scat_data_single.T_grid)>1001.)
    {
      ostringstream os;
      os << "The temperature values in the single scattering data" 
         << " are negative or very large. Check whether you use the "
         << "right unit [Kelvin].";
      throw runtime_error( os.str() );
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



/*! Derives weights of a bin-type quadrature for arbitrary wide bins.
 *
 * Note: Rectangular and trapezoidal rule essentially give the same weights
 * when considering the same nodes x (not the mid-points between the nodes as
 * rectangular is sometimes using) and limiting the quadrature range by the
 * first and last node).
 * Keyword order=0 calculates rectangular bins (ie bins extend beyond the first
 * and last nodes), order=1 to trapezoidal bins (ie bins end exactly at nodes).
         
 \param w      resulting weights at ordinates.
 \param x      ordinates.
 \param order  order of quadrature (see above for details).
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2017-06-16

*/
void bin_quadweights( Vector& w,
                const Vector& x,
                const Index& order )
{
  Index nx = x.nelem();

  assert( nx>1 );
  assert( is_increasing(x) );

  if( order==0 )
    {
      w[0] = min(x[1]-x[0],0.5*(x[1]+x[0])); // the latter is the half distance
                                             // from x0 to x1 plus the distance
                                             // to 0, ie 0.5(x1-x0)+x0.
      w[nx-1] = x[nx-1]-x[nx-2];
    }
  else
    {
      w[0] = 0.5*(x[1]-x[0]);
      w[nx-1] = 0.5*(x[nx-1]-x[nx-2]);
    }
  for (Index i = 1; i<nx-1; i++)
    {
      w[i] = 0.5*(x[i+1]-x[i-1]);
    }
}



/*! Derives bin-wise integral of y over x
 *
 * Used eg for Scaling pnd values by width of size bin. 
 * Bin width is determined from preceeding and following ordinates of a node.
 * Vector y and x must be equal in size. Vector w holds the local (ie over the
 * individual bins) quadrature results.
 *
 * Derived from rectengular/trapezoid integration rule (both provide the same
 * results).
         
  \param w   local quadrature result of y over x bins
  \param x   function nodes, eg scattering element size in [m]
  \param y   function values, eg particle number density per size interval [#/m3/m]
  
  \author Daniel Kreyling, Jana Mendrok
  \date 2010-12-15

*/
void bin_integral( Vector& w,
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
      bin_quadweights(w, x);
      for (Index i = 0; i<x.nelem(); i++)
        w[i] *= y[i];
    }
    else // for monodisperse pnd=dNdD
    {
      w[0] = y[0];
    }
}



/*! Check sum of pnd vector against total mass density value.
 *  Deviation is calculated and used to adjust the output of vector pnd.
         
	\param pnd   particle number density [#/m3]
	\param xwc   scattering species mass density [kg/m3]
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
  CREATE_OUT3;
  
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
    /*
    cout<< "p = " << p << ", i: " << i
        << " - pnd[i]: " << pnd[i] << ", mass[i]: " << mass[i]
        << " => x[i]: " << x[i] << "\n";
    */
  }

  if ( x.sum() == 0.0 )
    { // when x.sum()==0, but xwc!=0, obviously something went wrong in pnd calc
      ostringstream os;
      os<< "ERROR: in WSM chk_pndsum in pnd_fieldCalcFromscat_speciesFields!\n" 
      << "Given mass density != 0, but calculated mass density == 0.\n"
      << "Seems, something went wrong in pnd_fieldCalcFromscat_speciesFields. Check!\n"
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
          //CREATE_OUT1;
          out3 << "WARNING: in WSM chk_pndsum in pnd_fieldCalcFromscat_speciesFields!\n" 
               << "The deviation of the sum of nodes in the particle size distribution\n"
               << "to the initial input mass density of '"<< partfield_name
               <<"' is larger than 10%!\n"
               << "The deviation of: "<< error-1.0<<" occured in the atmospheric level: "
               << "p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<".\n";
        }
    }

  out3 << "PND scaling factor in atm. level "
       << "(p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<"): "<< error <<"\n";
}



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
	      if ( scat_species_field(j,k,l) != 0.0 &&
            !std::isnan(scat_species_field(j,k,l)) ) empty_flag = true;
//	      if ( scat_species_field(j,k,l) != 0.0 ) empty_flag = true;
	}
      }
    }  
}


//! Adjust uppermost and lowermost cloudy level for one scat_species_*_*_field.
/*!

  lower and upper levels have to be preinitialized before calling this function.
 
  \param[in,out] lower              lowermost level containing scattering particles
  \param[in,out] upper              uppermost level containing scattering particles
  \param[in]     scat_species_field scattering species field (e.g. mass density, mass
                                    flux, total number density)
  \param[in]     atmosphere_dim     the atmosphere dimension
  \param[in]     cloudbox_margin    flag whether to determine lowermost level or set to
                                    surface

  \author Daniel Kreyling, Jana Mendrok
  \date   2015-02-09
*/  
void find_cloudlimits(Index&          lower,
                      Index&          upper,
                      const Tensor3&  scat_species_field,
                      const Index&    atmosphere_dim,
                      const Numeric&  cloudbox_margin)
{
  if ( atmosphere_dim == 1 )
  {
    // scattering species profiles
    ConstVectorView ss_prof = scat_species_field ( joker, 0 , 0 );

    Index i = 0;

    // find lower cloudbox_limit to surface if margin != -1 (cloudbox not
    // forced to reach down to surface)
    if ( cloudbox_margin != -1 )
    {
      // find index of first pressure level where hydromet_field is
      // unequal 0, starting from the surface
      for ( i=0; i<lower; i++ )
      {
        //cout << "for lower limit checking level #" << i << "\n";

        // if any of the scat species fields contains a non-zero, non-NaN
        // value at this atm level we found a potential lower limit value
        if ( ss_prof[i] != 0.0 && !std::isnan(ss_prof[i]) )
        {
          //cout << "found particles\n";

          // check if lower is the lowest index in all selected
          // scattering species fields
          if ( lower > i )
          {
            lower = i;
            //cout << "new lower limit at level #" << lower << "\n";
          }
          break;
        }
      }

    }

    // find index of highest pressure level, where scat_species_mass_density_field is
    // unequal 0, starting from top of the atmosphere
    for ( Index j=scat_species_field.npages()-1; j>=max(i,upper); j-- )
    {
      //cout << "for upper limit checking level #" << j << "\n";

      // if any of the scat species fields contains a non-zero, non-NaN
      // value at this atm level we found a potential lower limit value
      if ( ss_prof[j] != 0.0 && !std::isnan(ss_prof[j]) )
      {
        //cout << "found particles\n";

        // check if upper is the highest index in all selected
        // scattering species fields
        if ( upper < j )
        {
          upper = j;
          //cout << "new upper limit at level #" << upper << "\n";
        }
        break;
      }
    }
  }

  else
  {
    ostringstream os;
    os << "Not yet available for 2D and 3D cases.";
    throw runtime_error( os.str() );
  }

    /*  //NOT WORKING YET
      // Latitude limits
      else if ( atmosphere_dim == 2 )
      {
        MatrixView hydro_lat = hydromet_field ( nhyd, joker, joker, 0 );

        for ( i=0; i<hydro_lat.nrows(); i++ )
        {
          for ( j=0; j<hydro_lat.ncols(); j++ )
          {
            if ( hydro_lat[i,j] != 0.0 )
            {

              if ( lat1 <= j ) lat1 =j;
              //cloudbox_limits[2] = lat1;
              //break;
            }

          }
          if ( lower <= i )    lower = i;
        }

        for ( k=hydro_lat.nelem()-1; k>=i; k-- )
        {
          if ( hydro_lat[k] != 0.0 )
          {
            lat2 = k;
            cloudbox_limits[3] = lat2;
            break;

          }

        }
      }

      // Longitude limits
      if ( atmosphere_dim == 3 )
      {
        Tensor3View hydro_lon = hydromet_field ( nhyd, joker, joker, joker );

        for ( i=0; i<hydro_lon.nelem(); i++ )
        {
          if ( hydro_lon[i] != 0.0 )
          {
            lon1 = i;
            cloudbox_limits[4] = lon1;
            break;
          }

        }
        for ( j=hydro_lon.nelem()-1; j>=i; j-- )
        {
          if ( hydro_lon[j] != 0.0 )
          {
            lon2 = j;
            cloudbox_limits[5] = lon2;
            break;

          }
}*/
  
}
  
 

/*! Parse atm_field_compact fieldname for species type

  \param  species_type  species indentifier (first part of field_name)
  \param  field_name    fieldname of atm_field_compact data entry
  \param  delim         delimiter string of field_name

  \author Jana Mendrok
  \date 2014-11-20

*/
void parse_atmcompact_speciestype (//WS Output:
                                   String& species_type,
                                   // WS Input:
                                   const String& field_name,
                                   const String& delim)
{
  ArrayOfString strarr;

  // split field_name string at '-' and write to ArrayOfString
  field_name.split ( strarr, delim );

  // first entry is species type
  // (i.e. "abs_species" or "scat_species". or "T" or "z", which are ignored.)
  if (strarr.size()>0 && field_name[0]!='-')
  {
      species_type = strarr[0];
  }
  else
  {
      ostringstream os;
      os << "No information on field species type found in '"
         << field_name << "'\n";
      throw runtime_error ( os.str() );

  }
}


/*! Parse atm_field_compact fieldname for species name

  \param  species_type  species name (second part of field_name)
  \param  field_name    fieldname of atm_field_compact data entry
  \param  delim         delimiter string of field_name

  \author Jana Mendrok
  \date 2014-11-20

*/
void parse_atmcompact_speciesname (//WS Output:
                                   String& species_name,
                                   // WS Input:
                                   const String& field_name,
                                   const String& delim)
{
  ArrayOfString strarr;

  // split field_name string at '-' and write to ArrayOfString
  field_name.split ( strarr, delim );

  // second entry is species name
  // (e.g. "H2O, "O3" etc. for abs_species or "IWC", "LWC" etc. for scat_species)
  if (strarr.size()>1)
  {
      species_name = strarr[1];
  }
  else
  {
      ostringstream os;
      os << "No information on field species name found in '"
         << field_name << "'\n";
      throw runtime_error ( os.str() );

  }
}


/*! Parse atm_field_compact fieldname for type of scat_species field

  \param  scat_type     species name (second part of field_name)
  \param  field_name    fieldname of atm_field_compact data entry
  \param  delim         delimiter string of field_name

  \author Jana Mendrok
  \date 2014-11-20

*/
void parse_atmcompact_scattype (//WS Output:
                                String& scat_type,
                                // WS Input:
                                const String& field_name,
                                const String& delim)
{
  ArrayOfString strarr;

  // split field_name string at '-' and write to ArrayOfString
  field_name.split ( strarr, delim );

  // third entry is type of scat_species field
  // (e.g. "mass_density", "mass_flux", "number_density")
  if (strarr.size()>2)
  {
      scat_type = strarr[2];
  }
  else
  {
      ostringstream os;
      os << "No information on type of scat_species field found in '"
         << field_name << "'\n";
      throw runtime_error ( os.str() );

  }
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

/*! Splitting scat_species string and parse additional options
	\param  psd_options contents of part_string positions >2
	\param  part_string scattering species tag from *scat_species*
  \param  delim       delimiter string of *scat_species* elements
  
  \author Jana Mendrok
  \date 2016-06-03

*/
void parse_psd_options (//WS Output:
                        ArrayOfString& psd_options,
                        // WS Input:
                        const String& part_string,
                        const String& delim)
{
  ArrayOfString strarr;

  // split scat_species string at delim and write to ArrayOfString
  part_string.split ( strarr, delim );

  // everything beyond second entry can hold psd-specific options
  if (strarr.size()>2)
    {
        psd_options.resize(strarr.nelem()-2);
        std::copy(strarr.begin()+2, strarr.end(), psd_options.begin());
    }
    else
        psd_options.resize(0);
}

