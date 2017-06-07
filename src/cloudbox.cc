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
#include "lin_alg.h"



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
       << "Each scat_data entry must correspond to one entry in scat_meta.";
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
  This function checks, whether the single scattering 
  properties of a scattering element includes the required frequencies.

  \param scat_data_single[in]  Single scattering data of a single scattering element
  \param f_grid[in]            Frequency grid
  
  \author Claudia Emde
  \date   2005-04-04
*/
void chk_scat_data_fgrid(const SingleScatteringData& scat_data_single,
                         ConstVectorView f_grid,
                         const String& infostring)
{
  chk_interpolation_grids( infostring, scat_data_single.f_grid, f_grid);
  
/*  if (!(scat_data_single.f_grid[0] <= f_grid[0] &&
        last(f_grid) <= 
        last(scat_data_single.f_grid) ))
    {
      ostringstream os;
      os << "The range of frequency grid in the single scattering"
         << " properties data does not contain all values of"
         << "*f_grid*.";
      throw runtime_error( os.str() );
    }*/
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
  CREATE_OUT2;

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
    
    out2 << "  Data is for arbitrarily orientated particles. \n";
    
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
    
    out2 << "  Data is for macroscopically isotropic and mirror-symmetric "
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
    
    out2 << "  Data is for azimuthally randomly oriented particles. \n";
    
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
    robust = (robust || (psd_options[i]=="robust") );

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
            // a valid IWC value encountered. Calculating dNdD.
            if (IWC_field ( p, lat, lon ) > 0.)
              {
                Numeric T = t_field ( p, lat, lon );

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

                psdFromMH97 ( dNdD, diameter_mass_equivalent,
                              IWC_field(p,lat,lon), T, noisy );
/*
                // iteration over all given size bins
                for ( Index i=0; i<diameter_mass_equivalent.nelem(); i++ )
                  {

                    // calculate particle size distribution with MH97
                    // [# m^-3 m^-1]
                    dNdD[i] = IWCtopnd_MH97 ( IWC_field ( p, lat, lon ),
                                              diameter_mass_equivalent[i],
                                              t_field ( p, lat, lon ),
                                              noisy, robust );
                  }
*/

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
                  scale_pnd( pnd, diameter_mass_equivalent, dNdD );
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
                os << "Size distribution " << psdname << " requires knowledge of mass "
                   << "density of atmospheric ice.\n"
                   << "At grid point (" << p << ", " << lat << ", " << lon
                   << ") in (p,lat,lon) a NaN value is encountered, "
                   << "i.e. mass density is unknown.";
                throw runtime_error( os.str() );
              }

            // MH97 PSD is parameterized in IWC and can not handle negative
            // numbers, hence abort.
            else if (IWC_field ( p, lat, lon ) < 0.)
              {
                ostringstream os;
                os << "Size distribution " << psdname << " is parametrized in ice mass"
                   << "content.\n"
                   << "It can not handle negative values like IWC="
                   << IWC_field (p,lat,lon) << " kg/m3\n"
                   << "found at grid point ("
                   << p << ", " << lat << ", " << lon << ")";
                throw runtime_error( os.str() );
              }

            // for IWC==0, we just set pnd_field=0
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
      // itertation over all atm. levels
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
      // itertation over all atm. levels
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
      // itertation over all atm. levels
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
 *  distribution for tropics.
 *  For the estimation of the second moment the mass dimension
 *  relationship is estimated by regression from the meta data.
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
 
 \author Manfred Brath (parts of the function is based on pnd_H11 (of J. Mendrok & D. Kreyling))
 \date 2014-1202
 
 */
void pnd_fieldF07TR (Tensor4View pnd_field,
                   const Tensor3& SWC_field,
                   const Tensor3& t_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity)
{
    const String psdname="F07TR";
    const Index N_se = scat_meta[scat_species].nelem();
    const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
    ArrayOfIndex intarr;
    Vector diameter_max_unsorted ( N_se, 0.0 );
    Vector diameter_max ( N_se, 0.0 );
    Vector mass ( N_se, 0.0 );
    Vector pnd ( N_se, 0.0 );
    Vector dNdD ( N_se, 0.0 );
    String partfield_name;
    Numeric alpha;
    Numeric beta;
    Vector log_m( N_se, 0.0 );
    Vector log_D( N_se, 0.0 );
    Vector q;
    
    if ( diameter_max.nelem() > 0 )
    // diameter_max.nelem()=0 implies no selected scattering element for the respective
    // scattering species field. should not occur anymore.
    {
      //unit conversion
      const Numeric D0=1; //[m]
    
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
        
        // logarithm of Dmax, needed for estimating mass-dimension-relationship
        log_D[i]=log(diameter_max[i]/D0);
        
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
        alpha=mass[0]/(diameter_max[0]*diameter_max[0]);
      }
    
      CREATE_OUT2;
      out2 << "Mass-dimension relationship m=alpha*(dmax/D0)^beta:\n"
           << "alpha = " << alpha << " kg \n"
           << "beta = " << beta << "\n";

      // itertation over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // F07 requires mass density. If not set, abort calculation.
                    if ( isnan(SWC_field ( p, lat, lon )) )
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
                    else if (SWC_field ( p, lat, lon ) != 0.)
                    {
                        // iteration over all given size bins
                        for ( Index i=0; i<diameter_max.nelem(); i++ ) //loop over number of scattering elements
                        {
                            // calculate particle size distribution for H11
                            // [# m^-3 m^-1]
                            dNdD[i] = IWCtopnd_F07TR( diameter_max[i], t_field ( p, lat, lon ),
                                                   SWC_field ( p, lat, lon ), alpha, beta);
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
                            scale_pnd( pnd, diameter_max, dNdD ); //[# m^-3]
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
 *  distribution for mid latitude.
 *  For the estimation of the second moment the mass dimension
 *  relationship is estimated by regression from the meta data.
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
 
 \author Manfred Brath (parts of the function is based on pnd_H11 (of J. Mendrok & D. Kreyling))
 \date 2014-12-02
 
 */
void pnd_fieldF07ML (Tensor4View pnd_field,
                     const Tensor3& SWC_field,
                     const Tensor3& t_field,
                     const ArrayOfIndex& limits,
                     const ArrayOfArrayOfScatteringMetaData& scat_meta,
                     const Index& scat_species,
                     const String& part_string,
                     const String& delim,
                     const Verbosity& verbosity)
{
    const String psdname="F07ML";
    const Index N_se = scat_meta[scat_species].nelem();
    const Index scat_data_start = FlattenedIndex(scat_meta, scat_species);
    ArrayOfIndex intarr;
    Vector diameter_max_unsorted ( N_se, 0.0 );
    Vector diameter_max ( N_se, 0.0 );
    Vector mass ( N_se, 0.0 );
    Vector pnd ( N_se, 0.0 );
    Vector dNdD ( N_se, 0.0 );
    String partfield_name;
    Numeric alpha;
    Numeric beta;
    Vector log_m( N_se, 0.0 );
    Vector log_D( N_se, 0.0 );
    Vector q;
    
    if ( diameter_max.nelem() > 0 )
    // diameter_max.nelem()=0 implies no selected scattering element for the respective
    // scattering species field. should not occur anymore.
    {
      //unit conversion
      const Numeric D0=1; //[m]
    
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
        
        // logarithm of Dmax, needed for estimating mass-dimension-relationship
        log_D[i]=log(diameter_max[i]/D0);
        
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
        alpha=mass[0]/(diameter_max[0]*diameter_max[0]);
      }
    
      CREATE_OUT2;
      out2 << "Mass-dimension relationship m=alpha*(dmax/D0)^beta:\n"
           << "alpha = " << alpha << " kg \n"
           << "beta = " << beta << "\n";
    
      // itertation over all atm. levels
      for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // F07 requires mass density. If not set, abort calculation.
                    if ( isnan(SWC_field ( p, lat, lon )) )
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
                    else if (SWC_field ( p, lat, lon ) != 0.)
                    {
                        // iteration over all given size bins
                        for ( Index i=0; i<diameter_max.nelem(); i++ ) //loop over number of scattering elements
                        {
                            // calculate particle size distribution for H11
                            // [# m^-3 m^-1]
                            dNdD[i] = IWCtopnd_F07ML( diameter_max[i], t_field ( p, lat, lon ),
                                                     SWC_field ( p, lat, lon ), alpha, beta);
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
                            scale_pnd( pnd, diameter_max, dNdD ); //[# m^-3]
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


/*! Calculates the particle number density field according
 *  to the two moment scheme of Axel Seifert, that is used the ICON model.
 *  Important, depending on the hydrometeor type, there is a lower and a upper
 *  boundary for valid particle mass. For masses below and above these boundary
 *  the number density is set to zero but the mass is conserved.
 
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
        // itertation over all atm. levels
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
                        for ( Index i=0; i<mass.nelem(); i++ ) //loop over number of scattering elements
                        {
                            // calculate particle size distribution for H11
                            // [# m^-3 m^-1]
                            dNdD[i] = WCtopnd_S2M( mass[i], N_tot,
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
                                scale_pnd( pnd, mass, dNdD ); //[# m^-3]
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
        // itertation over all atm. levels
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
                            scale_pnd( pnd, diameter_volume_equ, dNdD ); //[# m^-3]
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
        // itertation over all atm. levels
        for ( Index p=limits[0]; p<limits[1]; p++ )
        {
            for ( Index lat=limits[2]; lat<limits[3]; lat++ )
            {
                for ( Index lon=limits[4]; lon<limits[5]; lon++ )
                {
                    // F07 requires mass density. If not set, abort calculation.
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
                            scale_pnd( pnd, diameter_volume_equ, dNdD ); //[# m^-3]
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


/*! Calculates particle size distribution of cloud ice using MH97 parametrization.
 *  
 *  Handles a vector of sizes at a time. Implicitly assumes particles of water
 *  ice. Strictly requires temperature to be <=273.15K (<=0C), i.e. calling
 *  method needs to ensure this.
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
void psdFromMH97 ( Vector& psd,
                   const Vector& diameter,
                   const Numeric& iwc,
                   const Numeric& t,
                   const bool noisy )
{
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

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
  Numeric T = t-273.15;

  assert (T<=0.);

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
  Vector dNdD1(nD, 0.);
  if (alphas100>0.)
  {
    Numeric Ns100 = 6*IWCs100 * pow ( alphas100,5. ) /
                    ( PI*cdensity*gamma_func ( 5. ) );//micron^-5
    for ( Index iD=0; iD<nD; iD++ )
      dNdD1[iD] = 1e18 * Ns100*d_um[iD] *
                  exp ( -alphas100*d_um[iD] ); //micron^-4 -> m^-3 micron^-1
  }


  //Log normal distribution component

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
      if ( !isnan(dNdD1[iD]) && !isnan(dNdD2[iD]) )
        psd[iD] = ( dNdD1[iD]+dNdD2[iD] ) *1e6; // m^-3 m^-1
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

/*! Calculates particle size distribution using F07 parametrization for
 *  tropics.
 *  Each diameter of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param d maximum diameter of scattering particle [m]
 \param t atmospheric temperature [K]
 \param swc snow water content [kg/m^3]
 \param alpha Factor for the mass-dimension (m=alpha*(Dmax/D0)^beta) relationship [kg]
 \param beta Exponent for the mass-dimension relationship [pure number]
 
 \author Manfred Brath
 \date 2014-11-28
 
 */
Numeric IWCtopnd_F07TR ( const Numeric d, const Numeric t,
                        const Numeric swc,const Numeric alpha,
                        const Numeric beta)
{
    Numeric dN;
    Numeric phi23;
    Numeric An;
    Numeric Bn;
    Numeric Cn;
    Numeric M2;
    Numeric Mn;
    
    Numeric x;
    
    Numeric Mbeta;
    Numeric base;
    Numeric pp;
    
    
    
    //factors of phi23
    Vector q=MakeVector(152,-12.4,3.28,-0.78,-1.94);
    
    //Factors of factors of the moment estimation parametrization
    Vector Aq=MakeVector(13.6,-7.76,0.479);
    Vector Bq=MakeVector(-0.0361,0.0151,0.00149);
    Vector Cq=MakeVector(0.807,0.00581,0.0457);
    
    
    //convert T from Kelvin to Celsius
    Numeric Tc = t-273.15;
    
    // estimate second moment
    if (beta==2)
        M2=swc/alpha;
    else
    {
        Mbeta=swc/alpha;
        
        // calculate factors of the moment estimation parametrization
        An=exp(Aq[0]+Aq[1]*beta+Aq[2]*pow(beta,2));
        Bn=Bq[0]+Bq[1]*beta+Bq[2]*pow(beta,2);
        Cn=Cq[0]+Cq[1]*beta+Cq[2]*pow(beta,2);
        
        base=Mbeta*exp(-Bn*Tc)/An;
        pp=1./(Cn);
        
        M2=pow(base,pp);
    }
    
    // order of the moment parametrization
    const Numeric n=3;
    
    // calculate factors of the moment estimation parametrization
    An=exp(Aq[0]+Aq[1]*n+Aq[2]*pow(n,2));
    Bn=Bq[0]+Bq[1]*n+Bq[2]*pow(n,2);
    Cn=Cq[0]+Cq[1]*n+Cq[2]*pow(n,2);
    
    // moment parametrization
    Mn=An*exp(Bn*Tc)*pow(M2,Cn);
    
    //Define x
    x=d*M2/Mn;
    
    //Characteristic function
    phi23 = q[0]*exp(q[1]*x)+q[2]*pow(x,q[3])*exp(q[4]*x);
    
    //Distribution function
    dN=phi23*pow(M2,4.)/pow(Mn,3.);
    
    if (isnan(dN)) dN = 0.0;
    return dN;
}

/*! Calculates particle size distribution using F07 parametrization for
 *  mid latitude.
 *  Each diameter of the scattering particles is a node in the distribution.
 *  One call of this function calculates one particle number density.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param d maximum diameter of scattering particle [m]
 \param t atmospheric temperature [K]
 \param swc snow water content [kg/m^3]
 \param alpha Factor for the mass-dimension (m=alpha*(Dmax/D0)^beta) relationship [kg]
 \param beta Exponent for the mass-dimension relationship [pure number]
 
 \author Manfred Brath
 \date 2014-11-14
 
 */
Numeric IWCtopnd_F07ML ( const Numeric d, const Numeric t,
                        const Numeric swc,const Numeric alpha,
                        const Numeric beta)
{
    Numeric dN;
    Numeric phi23;
    Numeric An;
    Numeric Bn;
    Numeric Cn;
    Numeric M2;
    Numeric Mn;
    
    Numeric x;
    
    Numeric Mbeta;
    Numeric base;
    Numeric pp;
    
    //factors of phi23
    Vector q=MakeVector(141,-16.8,102,2.07,-4.82);
    
    //Factors of factors of the moment estimation parametrization
    Vector Aq=MakeVector(13.6,-7.76,0.479);
    Vector Bq=MakeVector(-0.0361,0.0151,0.00149);
    Vector Cq=MakeVector(0.807,0.00581,0.0457);
    
    
    //convert T from Kelvin to Celsius
    Numeric Tc = t-273.15;
    
    // estimate second moment
    if (beta==2)
        M2=swc/alpha;
    else
    {
        Mbeta=swc/alpha;
        
        // calculate factors of the moment estimation parametrization
        An=exp(Aq[0]+Aq[1]*beta+Aq[2]*pow(beta,2));
        Bn=Bq[0]+Bq[1]*beta+Bq[2]*pow(beta,2);
        Cn=Cq[0]+Cq[1]*beta+Cq[2]*pow(beta,2);
        
        base=Mbeta*exp(-Bn*Tc)/An;
        pp=1./Cn;
        
        M2=pow(base,pp);
    }
    
    // order of the moment parametrization
    const Numeric n=3;
    
    // calculate factors of the moment estimation parametrization
    An=exp(Aq[0]+Aq[1]*n+Aq[2]*pow(n,2));
    Bn=Bq[0]+Bq[1]*n+Bq[2]*pow(n,2);
    Cn=Cq[0]+Cq[1]*n+Cq[2]*pow(n,2);
    
    // moment parametrization
    Mn=An*exp(Bn*Tc)*pow(M2,Cn);
    
    //Define x
    x=d*M2/Mn;
    
    //Characteristic function
    phi23 = q[0]*exp(q[1]*x)+q[2]*pow(x,q[3])*exp(q[4]*x);
    
    //Distribution function
    dN=phi23*pow(M2,4.)/pow(Mn,3.);
    
    if (isnan(dN)) dN = 0.0;
    return dN;
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
 *  Important, depending on the hydrometeor type, there is a lower and a upper
 *  boundary for valid particle mass. For masses below and above these boundary
 *  the number density is set to zero. Within this function the mass is NOT 
 *  conserved.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param mass   Mass of scattering particle [kg]
 \param N_tot  Total number of particles (0th moment) [#/m3/m/kg^mu]
 \param M      Total mass concentration of Particles (1st moment) [kg/m^3]
 \param psd_type string with a tag defining the (hydrometeor) scheme

 
 \author Manfred Brath
 \date 2015-01-19
 
 */
Numeric WCtopnd_S2M (const Numeric mass,
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
    Numeric xmin;
    Numeric xmax;
    
    
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
        mu=0;
        gamma=1./2.;
        xmin=1e-10;
        xmax=2e-5;
    }
    else if ( psd_type == "S2M_GWC" ) //Graupel
    {
        mu=1;
        gamma=1./3.;
        xmin=1e-9;
        xmax=5e-4;
    }
    else if ( psd_type == "S2M_HWC" ) //Hail
    {
        mu=1;
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

    
    
    
    //Calculate Number density only, if mass is between xmin and xmax
    if (xmin<mass && mass<xmax)
    {
    
        //arguments for Gamma function
        arg2=(mu+2)/gamma;
        arg1=(mu+1)/gamma;
        
        
        temp=N_tot/M*gamma_func(arg2)/gamma_func(arg1);
        
        //Lambda (parameter for modified gamma distribution)
        Lambda=pow(temp, gamma);
        
        //N0
        N0=N_tot*gamma/gamma_func(arg1)*pow(Lambda, arg1);
        
        //Distribution function
        dN=mod_gamma_dist(mass, N0,Lambda, mu, gamma);
        
        if (isnan(dN)) dN = 0.0;
    }
    else
    {
        dN=0;
    }
    
    
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
          out2<< "WARNING: in WSM chk_pndsum in pnd_fieldCalcFromscat_speciesFields!\n" 
              << "The deviation of the sum of nodes in the particle size distribution\n"
              << "to the initial input mass density of '"<< partfield_name
              <<"' is larger than 10%!\n"
              << "The deviation of: "<< error-1.0<<" occured in the atmospheric level: "
              << "p = "<<p<<", lat = "<<lat<<", lon = "<<lon<<".\n";
        }
    }

  out2 << "PND scaling factor in atm. level "
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
            !isnan(scat_species_field(j,k,l)) ) empty_flag = true;
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
        if ( ss_prof[i] != 0.0 && !isnan(ss_prof[i]) )
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
      if ( ss_prof[j] != 0.0 && !isnan(ss_prof[j]) )
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

