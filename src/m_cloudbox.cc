/* Copyright (C) 2002-2012
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
   Claudia Emde     <claudia.emde@dlr.de>
   Cory Davis       <cory.davis@metservice.com>	   
                         
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
  === File description 
  ===========================================================================*/

/*!
  \file   m_cloudbox.cc
  \author Patrick Eriksson, Claudia Emde and Sreerekha T. R.
  \date   2002-05-08 

  \brief  Workspace functions related to the definintion of the cloud box.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <cstdlib>
#include <cmath>

#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "xml_io.h"
#include "messages.h"
#include "gridded_fields.h"
#include "logic.h"
#include "rte.h"
#include "interpolation.h"
#include "special_interp.h"
#include "cloudbox.h"
#include "optproperties.h"
#include "math_funcs.h"
#include "physics_funcs.h"
#include "sorting.h"
#include "file.h"
#include "parameters.h"

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Numeric PI;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxOff (
         Index&                       cloudbox_on,
         ArrayOfIndex&                cloudbox_limits,
         Agenda&                      iy_cloudbox_agenda,
         Tensor4&                     pnd_field,
         ArrayOfArrayOfSingleScatteringData& scat_data,
         Matrix&                      particle_masses,
   const Verbosity&)
{
  cloudbox_on = 0;
  cloudbox_limits.resize ( 0 );
  iy_cloudbox_agenda = Agenda();
  iy_cloudbox_agenda.set_name ( "iy_cloudbox_agenda" );
  pnd_field.resize(0,0,0,0);
  scat_data.resize(0);
  particle_masses.resize(0,0);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetAutomatically (// WS Output:
                               //Workspace& /* ws */,
                               Index&          cloudbox_on,
                               ArrayOfIndex&   cloudbox_limits,
                               //Agenda&  iy_cloudbox_agenda,
                               // WS Input:
                               const Index&    atmosphere_dim,
                               const Vector&   p_grid,
                               const Vector&   lat_grid,
                               const Vector&   lon_grid,
                               const Tensor4&  scat_species_mass_density_field,
                               const Tensor4&  scat_species_mass_flux_field,
                               const Tensor4&  scat_species_number_density_field,
                               // Control Parameters
                               const Numeric&  cloudbox_margin,
                               const Verbosity& verbosity)
{
  // Check existing WSV
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );
  // includes p_grid chk_if_decresasing
  chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid ); 
  // Set cloudbox_on
  cloudbox_on = 1;

  if ( atmosphere_dim > 1 )
    {
      ostringstream os;
      os << "cloudboxSetAutomatically not yet available for 2D and 3D cases.";
      throw runtime_error( os.str() );
    }

  // Allocate cloudbox_limits
  cloudbox_limits.resize ( atmosphere_dim*2 );

  // Variables
  Index p1;
  if ( cloudbox_margin == -1 )
    {
      cloudbox_limits[0] = 0;
      p1 = 0;
    }
  else p1 = scat_species_mass_density_field.npages()-1;
  Index p2 = 0;

  // OLE: Commented out until code that uses it at the end of this function is commented back in
//  if ( atmosphere_dim > 1 )
//    {
//      Index lat1 = scat_species_mass_density_field.nrows()-1;
//      Index lat2 = 0;
//    }
//  if ( atmosphere_dim > 2 )
//    {
//      Index lon1 = scat_species_mass_density_field.ncols()-1;
//      Index lon2 = 0;
//    }

  bool not_empty_any=false;
  bool not_empty_md=true;
  bool not_empty_mf=true;
  bool not_empty_nd=true;

  Index nss=0;

  if (scat_species_mass_density_field.nbooks() == 0)
    not_empty_md = 0;
  else
    nss = scat_species_mass_density_field.nbooks();

  if (scat_species_mass_flux_field.nbooks() == 0)
    not_empty_mf = 0;
  else if (nss!=0)
    {
      if (nss!=scat_species_mass_flux_field.nbooks())
        {
          ostringstream os;
          os << "Inconsistent number of scattering elements in\n"
             << "scat_species_mass_density_field and "
             << "scat_species_mass_flux_field.";
          throw runtime_error( os.str() );
        }
     }
  else
    nss = scat_species_mass_flux_field.nbooks();

  if (scat_species_number_density_field.nbooks() == 0)
    not_empty_nd = 0;
  else if (nss!=0)
    {
      if (nss!=scat_species_number_density_field.nbooks())
        {
          ostringstream os;
          os << "Inconsistent number of scattering elements in\n"
             << "scat_species_number_density_field and "
             << "scat_species_mass_density/flux_field.";
          throw runtime_error( os.str() );
        }
     }
  else
    nss = scat_species_mass_flux_field.nbooks();


  //--------- Start loop over scattering species ------------------------------
  for ( Index l=0; l<nss; l++ )
  {
    //cout << "for scatt species #" << l << ":\n";
    bool not_empty;

    //not_empty is set to true, if a single value of scat_species_XX_field
    //is unequal zero (and not NaN), i.e. if we actually have some amount of
    //these scattering species in the atmosphere.
    if (not_empty_md)
    {
      chk_scat_species_field ( not_empty,
                               scat_species_mass_density_field ( l, joker, joker, joker ),
                               "scat_species_mass_density_field",
                               atmosphere_dim, p_grid, lat_grid, lon_grid );
      //if particles found, enter detailed search
      if (not_empty)
      {
        not_empty_any=true;
        find_cloudlimits(p1, p2, scat_species_mass_density_field ( l, joker, joker, joker ),
                         atmosphere_dim, cloudbox_margin);
      }
      //cout << "particles in mass density field: " << not_empty << "\n";
    }

    if (not_empty_mf)
    {
      chk_scat_species_field ( not_empty,
                               scat_species_mass_flux_field ( l, joker, joker, joker ),
                               "scat_species_mass_flux_field",
                               atmosphere_dim, p_grid, lat_grid, lon_grid );
      if (not_empty)
      {
        not_empty_any=true;
        find_cloudlimits(p1, p2, scat_species_mass_flux_field ( l, joker, joker, joker ),
                         atmosphere_dim, cloudbox_margin);
      }
      //cout << "particles in mass flux field: " << not_empty << "\n";
    }

    if (not_empty_nd)
    {
      chk_scat_species_field ( not_empty,
                               scat_species_number_density_field ( l, joker, joker, joker ),
                               "scat_species_number_density_field",
                               atmosphere_dim, p_grid, lat_grid, lon_grid );
      if (not_empty)
      {
        not_empty_any=true;
        find_cloudlimits(p1, p2, scat_species_number_density_field ( l, joker, joker, joker ),
                         atmosphere_dim, cloudbox_margin);
      }
      //cout << "particles in number density field: " << not_empty << "\n";
    }
    //cout << "particles in any field: " << not_empty_any << "\n";
  }

  // decrease lower cb limit by one to ensure that linear interpolation of 
  // particle number densities is possible.
  Index p0 = 0; //only for the use of function *max*
  p1 = max(p1-1, p0);

  Numeric p_margin1;

  // alter lower cloudbox_limit by cloudbox_margin, using barometric
  // height formula
  p_margin1 = barometric_heightformula ( p_grid[p1], cloudbox_margin );
  Index k = 0;
  while ( p_grid[k+1] >= p_margin1 && k+1 < p_grid.nelem() ) k++;
  cloudbox_limits[0]= k;

  // increase upper cb limit by one to ensure that linear interpolation of 
  // particle number densities is possible.
  p2 = min(p2+1, scat_species_mass_density_field.npages()-1);
  // set upper cloudbox_limit
  // if cloudbox reaches to the upper most pressure level
  if ( p2 >= scat_species_mass_density_field.npages()-1)
  {
    CREATE_OUT2;
    out2<<"The cloud reaches to TOA!\n"
    <<"Check scat_species_mass_density_field data, if realistic!\n";
  }
  cloudbox_limits[1] = p2;

  //out0<<"\n"<<p2<<"\n"<<p_grid[p2]<<"\n";

  // check if all selected scattering species fields are zero at each level,
  // than switch cloudbox off, skipping scattering calculations
  if ( !not_empty_any )
  {
    CREATE_OUT0;
    //cloudboxOff ( cloudbox_on, cloudbox_limits, iy_cloudbox_agenda );
    cloudbox_on = 0;
    out0<<"Cloudbox is switched off!\n";

    return;
  }


  // assert keyword arguments

  // The pressure in *p1* must be bigger than the pressure in *p2*.
  assert ( p_grid[p1] > p_grid[p2] );
  // The pressure in *p1* must be larger than the last value in *p_grid*.
  assert ( p_grid[p1] > p_grid[p_grid.nelem()-1] );
  // The pressure in *p2* must be smaller than the first value in *p_grid*."
  assert ( p_grid[p2] < p_grid[0] );

  /*
  if ( atmosphere_dim >= 2 )
  {
    // The latitude in *lat2* must be bigger than the latitude in *lat1*.
    assert ( lat_grid[lat2] > lat_grid[lat1] );
    // The latitude in *lat1* must be >= the second value in *lat_grid*.
    assert ( lat_grid[lat1] >= lat_grid[1] );
    // The latitude in *lat2* must be <= the next to last value in *lat_grid*.
    assert ( lat_grid[lat2] <= lat_grid[lat_grid.nelem()-2] );
  }
  if ( atmosphere_dim == 3 )
  {
    // The longitude in *lon2* must be bigger than the longitude in *lon1*.
    assert ( lon_grid[lon2] > lon_grid[lon1] );
    // The longitude in *lon1* must be >= the second value in *lon_grid*.
    assert ( lon_grid[lon1] >= lon_grid[1] );
    // The longitude in *lon2* must be <= the next to last value in *lon_grid*.
    assert ( lon_grid[lon2] <= lon_grid[lon_grid.nelem()-2] );
  }
  */
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetManually(// WS Output:
                         Index&         cloudbox_on,
                         ArrayOfIndex&  cloudbox_limits,
                         // WS Input:
                         const Index&   atmosphere_dim,
                         const Vector&  p_grid,
                         const Vector&  lat_grid,
                         const Vector&  lon_grid,
                         // Control Parameters:
                         const Numeric& p1,
                         const Numeric& p2,
                         const Numeric& lat1,
                         const Numeric& lat2,
                         const Numeric& lon1,
                         const Numeric& lon2,
                         const Verbosity&)
{
  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Check keyword arguments
  if( p1 <= p2 )
    throw runtime_error( "The pressure in *p1* must be bigger than the "
                         "pressure in *p2*." );
  if( p1 <= p_grid[p_grid.nelem()-1] )
    throw runtime_error( "The pressure in *p1* must be larger than the "
                         "last value in *p_grid*." );
  if( p2 >= p_grid[0] )
    throw runtime_error( "The pressure in *p2* must be smaller than the "
                         "first value in *p_grid*." );
  if( atmosphere_dim >= 2 )
    {
      if( lat2 <= lat1 )
        throw runtime_error( "The latitude in *lat2* must be bigger than the "
                             "latitude in *lat1*.");
      if( lat1 < lat_grid[1] )
        throw runtime_error( "The latitude in *lat1* must be >= the "
                             "second value in *lat_grid*." );
      if( lat2 > lat_grid[lat_grid.nelem()-2] )
        throw runtime_error( "The latitude in *lat2* must be <= the "
                             "next to last value in *lat_grid*." );
    }
  if( atmosphere_dim == 3 )
    {
      if( lon2 <= lon1 )
        throw runtime_error( "The longitude in *lon2* must be bigger than the "
                             "longitude in *lon1*.");
      if( lon1 < lon_grid[1] )
        throw runtime_error( "The longitude in *lon1* must be >= the "
                             "second value in *lon_grid*." );
      if( lon2 > lon_grid[lon_grid.nelem()-2] )
        throw runtime_error( "The longitude in *lon2* must be <= the "
                             "next to last value in *lon_grid*." );
    }

  // Set cloudbox_on
  cloudbox_on = 1;

  // Allocate cloudbox_limits
  cloudbox_limits.resize( atmosphere_dim*2 );

  // Pressure limits
  if( p1 > p_grid[1] )
    {
      cloudbox_limits[0] = 0;
    }
  else
    {
      for( cloudbox_limits[0]=1; p_grid[cloudbox_limits[0]+1]>=p1; 
                                                     cloudbox_limits[0]++ ) {}
    }
  if( p2 < p_grid[p_grid.nelem()-2] )
    {
      cloudbox_limits[1] = p_grid.nelem() - 1;
    }
  else
    {
      for( cloudbox_limits[1]=p_grid.nelem()-2; 
                    p_grid[cloudbox_limits[1]-1]<=p2; cloudbox_limits[1]-- ) {}
    }

  // Latitude limits
  if( atmosphere_dim >= 2 )
    {
      for( cloudbox_limits[2]=1; lat_grid[cloudbox_limits[2]+1]<=lat1; 
                                                     cloudbox_limits[2]++ ) {}
      for( cloudbox_limits[3]=lat_grid.nelem()-2; 
                lat_grid[cloudbox_limits[3]-1]>=lat2; cloudbox_limits[3]-- ) {}
    }

  // Longitude limits
  if( atmosphere_dim == 3 )
    {
      for( cloudbox_limits[4]=1; lon_grid[cloudbox_limits[4]+1]<=lon1; 
                                                     cloudbox_limits[4]++ ) {}
      for( cloudbox_limits[5]=lon_grid.nelem()-2; 
                lon_grid[cloudbox_limits[5]-1]>=lon2; cloudbox_limits[5]-- ) {}
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetManuallyAltitude(// WS Output:
                                 Index&         cloudbox_on,
                                 ArrayOfIndex&  cloudbox_limits,
                                 // WS Input:
                                 const Index&   atmosphere_dim,
                                 const Tensor3& z_field,
                                 const Vector&  lat_grid,
                                 const Vector&  lon_grid,
                                 // Control Parameters:
                                 const Numeric& z1,
                                 const Numeric& z2,
                                 const Numeric& lat1,
                                 const Numeric& lat2,
                                 const Numeric& lon1,
                                 const Numeric& lon2,
                                 const Verbosity&)
{
  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  // Check keyword arguments
  if( z1 >= z2 )
    throw runtime_error( "The altitude in *z1* must be smaller than the "
                         "altitude in *z2*." );
  /* These cases are in fact handled
  if( z1 <= z_field(0, 0, 0) )
    throw runtime_error( "The altitude in *z1* must be larger than the "
                         "first value in *z_field*." );
  if( z2 >= z_field(z_field.npages()-1, 0, 0) )
    throw runtime_error( "The altitude in *z2* must be smaller than the "
                         "last value in *z_field*." );
  */
  if( atmosphere_dim == 3 )
    {
      if( lat2 <= lat1 )
        throw runtime_error( "The latitude in *lat2* must be bigger than the "
                             " latitude in *lat1*.");
      if( lat1 < lat_grid[1] )
        throw runtime_error( "The latitude in *lat1* must be >= the "
                             "second value in *lat_grid*." );
      if( lat2 > lat_grid[lat_grid.nelem()-2] )
        throw runtime_error( "The latitude in *lat2* must be <= the "
                             "next to last value in *lat_grid*." );
      if( lon2 <= lon1 )
        throw runtime_error( "The longitude in *lon2* must be bigger than the "
                             "longitude in *lon1*.");
      if( lon1 < lon_grid[1] )
        throw runtime_error( "The longitude in *lon1* must be >= the "
                             "second value in *lon_grid*." );
      if( lon2 > lon_grid[lon_grid.nelem()-2] )
        throw runtime_error( "The longitude in *lon2* must be <= the "
                             "next to last value in *lon_grid*." );
    }
  
  // Set cloudbox_on
  cloudbox_on = 1;

  // Allocate cloudbox_limits
  cloudbox_limits.resize( atmosphere_dim*2 );

  // Pressure/altitude limits
  if( z1 < z_field(1, 0, 0) )
    {
      cloudbox_limits[0] = 0;
    }
  else
    {
      for( cloudbox_limits[0]=1; z_field(cloudbox_limits[0]+1, 0, 0) <= z1; 
                                                     cloudbox_limits[0]++ ) {}
    }
  if( z2 > z_field(z_field.npages()-2, 0, 0) )
    {
      cloudbox_limits[1] = z_field.npages() - 1;
    }
  else
    {
      for( cloudbox_limits[1]=z_field.npages()- 2; 
           z_field(cloudbox_limits[1]-1, 0, 0) >= z2; cloudbox_limits[1]-- ) {}
    }

  // Latitude limits
  if( atmosphere_dim >= 2 )
    {
      for( cloudbox_limits[2]=1; lat_grid[cloudbox_limits[2]+1]<=lat1; 
                                                     cloudbox_limits[2]++ ) {}
      for( cloudbox_limits[3]=lat_grid.nelem()-2; 
                lat_grid[cloudbox_limits[3]-1]>=lat2; cloudbox_limits[3]-- ) {}
    }

  // Longitude limits
  if( atmosphere_dim == 3 )
    {
      for( cloudbox_limits[4]=1; lon_grid[cloudbox_limits[4]+1]<=lon1; 
                                                     cloudbox_limits[4]++ ) {}
      for( cloudbox_limits[5]=lon_grid.nelem()-2; 
                lon_grid[cloudbox_limits[5]-1]>=lon2; cloudbox_limits[5]-- ) {}
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void scat_species_fieldCleanup (//WS Output:
                          Tensor4& scat_species_field_out,
                          //WS Input:
                          const Tensor4& scat_species_field_in,
                          const Numeric& threshold,
                          const Verbosity&)
{
  // First make a copy of scat_species_field_in as it might be identical with
  // scat_species_field_out. And resize scat_species_field_out appropriately.
  Tensor4 infield = scat_species_field_in;

  // Check that scat_species_field contains realistic values
  //(values smaller than the threshold will be set to 0)
  for ( Index i=0; i<scat_species_field_in.nbooks(); i++ )
  {
    for ( Index j=0; j<scat_species_field_in.npages(); j++ )
    {
      for ( Index k=0; k<scat_species_field_in.nrows(); k++ )
      {
        for ( Index l=0; l<scat_species_field_in.ncols(); l++ )
        {
          if ( scat_species_field_in ( i,j,k,l ) < threshold ) 
          {
            infield ( i,j,k,l ) = 0.0;
          }
        }
      }
    }
  }
  scat_species_field_out = infield;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void scat_speciesInit (ArrayOfString&  scat_species,
                       const Verbosity&)
{
  scat_species.resize ( 0 );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void scat_speciesSet (// WS Generic Output:
                      ArrayOfString&  scat_species,
                      // Control Parameters:
                      const ArrayOfString& scatspecies_tags,
                      const String& delim,
                      const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  scat_species.resize ( scatspecies_tags.nelem() );
  //assign input strings to scat_species
  scat_species = scatspecies_tags;

  chk_scat_species (scat_species, delim);

  // Print list of scattering species settings to the most verbose output stream:
  out3 << "  Defined scattering species settings: ";
  for ( Index i=0; i<scat_species.nelem(); ++i )
  {
    out3 << "\n  " << i << ": "<<scat_species[i];

  }
  out3 << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeInit (//WS Output:
                       ArrayOfArrayOfSingleScatteringData& scat_data,
                       ArrayOfGriddedField3& pnd_field_raw,
                       const Verbosity&)
{
  scat_data.resize(0);
  pnd_field_raw.resize(0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeAdd( //WS Output:
                     ArrayOfArrayOfSingleScatteringData& scat_data,
                     ArrayOfGriddedField3&  pnd_field_raw,
                     // WS Input (needed for checking the datafiles):
                     const Index& atmosphere_dim,
                     const Vector& f_grid,
//                     const Vector& p_grid,
//                     const Vector& lat_grid,
//                     const Vector& lon_grid,
                     // Keywords:
                     const String& scat_data_file,
                     const String& pnd_field_file,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  //--- Check input ---------------------------------------------------------
  
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  //chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Frequency grid
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  

  //--- Reading the data ---------------------------------------------------

  Index last_species = scat_data.nelem()-1;
  if (last_species == -1)
  {
      scat_data.resize(1);
      last_species = 0;
  }

  // Append *scat_data* and *pnd_field_raw* with empty Arrays of Tensors.
  SingleScatteringData scat_data_single;
  scat_data[last_species].push_back(scat_data_single);
  
  GriddedField3 pnd_field_data;
  pnd_field_raw.push_back(pnd_field_data);
  
  out2 << "  Read single scattering data\n";
  xml_read_from_file(scat_data_file,
                     scat_data[last_species][scat_data[last_species].nelem()-1],
                     verbosity);

  chk_scat_data(scat_data[last_species][scat_data[last_species].nelem()-1],
                             scat_data_file, f_grid, verbosity);       
  
  out2 << "  Read particle number density field\n";
  if (pnd_field_file.nelem() < 1)
  {
    CREATE_OUT1;
    out1 << "Warning: No pnd_field_file specified. Ignored. \n";
  }
  else
    {
      xml_read_from_file(pnd_field_file, pnd_field_raw[pnd_field_raw.nelem()-1],
                         verbosity);
      
      chk_pnd_data(pnd_field_raw[pnd_field_raw.nelem()-1],
                   pnd_field_file, atmosphere_dim, verbosity);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeAddAll (//WS Output:
                         ArrayOfArrayOfSingleScatteringData& scat_data,
                         ArrayOfGriddedField3&  pnd_field_raw,
                         // WS Input(needed for checking the datafiles):
                         const Index& atmosphere_dim,
                         const Vector& f_grid,
//                         const Vector& p_grid,
//                         const Vector& lat_grid,
//                         const Vector& lon_grid,
                         // Keywords:
                         const String& filelist_scat_data,
                         const String& pnd_fieldarray_file,
                         const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  //--- Check input ---------------------------------------------------------

  // Atmosphere
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );
  //chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Frequency grid
  if ( f_grid.nelem() == 0 )
    throw runtime_error ( "The frequency grid is empty." );
  chk_if_increasing ( "f_grid", f_grid );


  //--- Reading the data ---------------------------------------------------
  ArrayOfString data_files;
  xml_read_from_file ( filelist_scat_data, data_files, verbosity );
  scat_data.resize(1);
  scat_data[0].resize ( data_files.nelem() );

  for ( Index i = 0; i<data_files.nelem(); i++ )
  {

    out2 << "  Read single scattering data\n";
    xml_read_from_file ( data_files[i], scat_data[0][i], verbosity );

    chk_scat_data ( scat_data[0][i],
                                 data_files[i], f_grid,
                                 verbosity );

  }

  out2 << "  Read particle number density data \n";
  xml_read_from_file ( pnd_fieldarray_file, pnd_field_raw, verbosity );

  chk_pnd_raw_data ( pnd_field_raw,
                     pnd_fieldarray_file, atmosphere_dim, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleType2abs_speciesAdd( //WS Output:
                     ArrayOfArrayOfSingleScatteringData& scat_data,
                     ArrayOfGriddedField3& vmr_field_raw,
                     ArrayOfArrayOfSpeciesTag& abs_species,
                     Index& propmat_clearsky_agenda_checked,
                     Index& abs_xsec_agenda_checked,
                     // WS Input (needed for checking the datafiles):
                     const Index& atmosphere_dim,
                     const Vector& f_grid,
//                     const Vector& p_grid,
//                     const Vector& lat_grid,
//                     const Vector& lon_grid,
                     // Keywords:
                     const String& scat_data_file,
                     const String& pnd_field_file,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  //--- Check input ---------------------------------------------------------
  
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  //chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Frequency grid
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  

  //--- Reading the data ---------------------------------------------------

  Index last_species = scat_data.nelem()-1;
  if (last_species == -1)
  {
      scat_data.resize(1);
      last_species = 0;
  }

  // Append *scat_data* and *pnd_field_raw* with empty Arrays of Tensors.
  SingleScatteringData scat_data_single;
  scat_data[last_species].push_back(scat_data_single);
  
  out2 << "  Read single scattering data\n";
  xml_read_from_file(scat_data_file,
                     scat_data[last_species][scat_data[last_species].nelem()-1],
                     verbosity);

  chk_scat_data(scat_data[last_species][scat_data[last_species].nelem()-1],
                             scat_data_file, f_grid, verbosity);       
  
  out2 << "  Read particle number density field\n";
  if (pnd_field_file.nelem() < 1)
  {
    CREATE_OUT1;
    out1 << "Warning: No pnd_field_file specified. Ignored here,\n"
         << "but user HAS TO add that later on!\n";
  }
  else
    {
      GriddedField3 pnd_field_data;
      vmr_field_raw.push_back(pnd_field_data);
  
      xml_read_from_file(pnd_field_file, vmr_field_raw[vmr_field_raw.nelem()-1],
                         verbosity);
      
      chk_pnd_data(vmr_field_raw[vmr_field_raw.nelem()-1],
                   pnd_field_file, atmosphere_dim, verbosity);
    }

  out2 << "  Append 'particle' field to abs_species\n";
  ArrayOfString species;
  species.push_back("particles");

  abs_speciesAdd( abs_species,
                  propmat_clearsky_agenda_checked, abs_xsec_agenda_checked,
                  species, verbosity );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatteringParticleTypeAndMetaRead (//WS Output:
                                        ArrayOfArrayOfSingleScatteringData& scat_data,
                                        ArrayOfArrayOfScatteringMetaData& scat_meta,
                                        const Vector& f_grid,
                                        // Keywords:
                                        const ArrayOfString& scat_data_files,
                                        const Verbosity& verbosity)
{
  CREATE_OUT3;

  //--- Reading the data ---------------------------------------------------
  ArrayOfSingleScatteringData arr_ssd;
  ArrayOfScatteringMetaData arr_smd;

  arr_ssd.resize ( scat_data_files.nelem() );
  arr_smd.resize ( scat_data_files.nelem() );

  for ( Index i = 0; i<scat_data_files.nelem(); i++ )
    {
      out3 << "  Read single scattering data\n";
      xml_read_from_file ( scat_data_files[i], arr_ssd[i], verbosity );

      chk_scat_data ( arr_ssd[i],
                      scat_data_files[i], f_grid,
                      verbosity );

      // make meta data name from scat data name
      ArrayOfString strarr;
      scat_data_files[i].split ( strarr, ".xml" );
      String scat_meta_file = strarr[0]+".meta.xml";
      cout << "looking for " << scat_meta_file << "\n";

      out3 << "  Read scattering meta data\n";
      xml_read_from_file ( scat_meta_file, arr_smd[i], verbosity );
            
      //FIXME: currently nothing is done in chk_scattering_meta_data!
      chk_scattering_meta_data ( arr_smd[i],
                                 scat_meta_file, verbosity );
            
    }

  // check if arrays have same size
  chk_scattering_data ( arr_ssd,
                        arr_smd, verbosity );

  // append as new scattering species
  scat_data.push_back(arr_ssd);
  scat_meta.push_back(arr_smd);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatteringParticlesSelect (//WS Output:
                                ArrayOfArrayOfSingleScatteringData& scat_data,
                                ArrayOfArrayOfScatteringMetaData& scat_meta,
                                // WS Input:
                                const ArrayOfString& scat_species,
                                const String& delim,
                                const Verbosity& verbosity)
{ 
  CREATE_OUT1;
  CREATE_OUT3;
  //--- Adjusting data to user specified input (scat_species)-------------------

  // first check that sizes of scat_species and scat_data/scat_meta agree
  Index nspecies = scat_species.nelem();
  if ( nspecies != scat_data.nelem() || nspecies != scat_meta.nelem() )
    {
      ostringstream os;
      os << "Number of scattering species specified by scat_species does\n"
         << "not agree with number of scattering species in\n"
         << "scat_data or scat_meta:\n"
         << "scat_species has " << nspecies << " entries, while scat_data has "
         << scat_data.nelem() << " and scat_meta has " << scat_meta.nelem()
         << ".";
      throw runtime_error ( os.str() );
    }

  // create temporary containers for selected elements
  ArrayOfArrayOfSingleScatteringData scat_data_tmp;
  ArrayOfArrayOfScatteringMetaData scat_meta_tmp;
  scat_data_tmp.resize(scat_species.nelem());
  scat_meta_tmp.resize(scat_species.nelem());

  // loop over array of scat_species--------------------------------------------
  // no more sorting by material tag. only left to select scattering elements
  // within specified size range (in terms of volume equiv diameter).
  for ( Index i_ss=0; i_ss<nspecies; i_ss++ )
  {
   
    String partfield_name;
    Numeric sizemin;
    Numeric sizemax;
    const Numeric tolerance = 1e-6;

    //split scat_species string and copy values to parameter
    parse_partfield_name( partfield_name, scat_species[i_ss], delim);
    parse_part_size(sizemin, sizemax, scat_species[i_ss], delim);

  // choosing the specified SingleScatteringData and ScatteringMetaData
    for ( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
      {
        // scattering element volume equivalent diameter is extracted from the
        // scattering element's meta data and checked whether it's within size
        // selected range (sizemax < 0 check follows from wildcard usage and
        // means consider all sizes on the upper end)
        if ( scat_meta[i_ss][i_se].diameter_volume_equ*1e6 > sizemin-sizemin*tolerance &&
             ( sizemax+sizemax*tolerance > scat_meta[i_ss][i_se].diameter_volume_equ*1e6 ||
               sizemax < 0. ) )
          {
            // copy selected scattering element to temp arrays
            scat_data_tmp[i_ss].push_back(scat_data[i_ss][i_se]);
            scat_meta_tmp[i_ss].push_back(scat_meta[i_ss][i_se]);
          }
      }

    // To use a particle species field without associated scattering element
    // data poses a high risk of accidentially neglecting these species. That's
    // unlikely what the user intends. Hence throw error.
    if (scat_meta_tmp[i_ss].nelem()<1)
      {
        ostringstream os;
        os << "For scattering species " << partfield_name << " no scattering "
           << "element matching the requested size range found.\n"
           << "Check scat_data and scat_meta input as well as the scat_species "
           << "definition!";
        throw runtime_error ( os.str() );
      }
  }

  // check if array is empty. should only apply when scat_species is empty,
  // hence just post a warning, don't throw error.
  if ( !TotalNumberOfElements(scat_meta_tmp) )
  {
    out1 << "WARNING! No scattering elements selected.\n"
         << "Continuing without any selected scattering elements.\n";
  }

  scat_meta = scat_meta_tmp;
  scat_data = scat_data_tmp;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesFromMetaDataSingleCategory(
        Matrix&                    particle_masses,
  const ArrayOfArrayOfScatteringMetaData& scat_meta,
  const Verbosity&)
{
  const Index np_total = TotalNumberOfElements(scat_meta);

  particle_masses.resize(np_total,1);

  Index i_se_flat = 0;
  for( Index i_ss=0; i_ss<scat_meta.nelem(); i_ss++ )
  {
      for( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
      {
          if( isnan(scat_meta[i_ss][i_se].mass) ||
              scat_meta[i_ss][i_se].mass <= 0 ||
              scat_meta[i_ss][i_se].mass > 1. )
          {
              ostringstream os;
              os << "A presumably incorrect value found for "
              << "scat_meta[" << i_ss << "][" << i_se << "].mass.\n"
              << "The value is " << scat_meta[i_ss][i_se].mass;
              throw std::runtime_error(os.str());
          }

          if( scat_meta[i_ss][i_se].diameter_volume_equ <= 0 ||
             scat_meta[i_ss][i_se].diameter_volume_equ > 0.5 )
          {
              ostringstream os;
              os << "A presumably incorrect value found for "
              << "scat_meta[" << i_ss << "][" << i_se << "].diameter_volume_equ.\n"
              << "The value is " << scat_meta[i_ss][i_se].diameter_volume_equ;
              throw std::runtime_error(os.str());
          }

          particle_masses(i_se_flat,0) = scat_meta[i_ss][i_se].diameter_volume_equ;

          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesFromMetaData
                        (//WS Output:
                         Matrix& particle_masses,
                         // WS Input:
                         const ArrayOfArrayOfScatteringMetaData& scat_meta,
                         const Verbosity& )
{
  // resize particle_masses to required diemsions and properly initialize values
  particle_masses.resize ( TotalNumberOfElements(scat_meta), scat_meta.nelem() );
  particle_masses = 0.;

  // calculate and set particle_masses
  Index i_se_flat = 0;
  for ( Index i_ss=0; i_ss<scat_meta.nelem(); i_ss++ )
  {
    for ( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
    {
      particle_masses (i_se_flat, i_ss) = scat_meta[i_ss][i_se].diameter_volume_equ;
      i_se_flat++;
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFrompnd_field_raw(//WS Output:
                   Tensor4& pnd_field,
                   //WS Input
                   const Vector& p_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const ArrayOfGriddedField3& pnd_field_raw,
                   const Index& atmosphere_dim,
                   const ArrayOfIndex& cloudbox_limits,
                   const Index& zeropadding,
                   const Verbosity& verbosity)
{
  // Basic checks of input variables
  //
  // Particle number density data
  // 
  if (pnd_field_raw.nelem() == 0)
  {
    ostringstream os;
    os << "No particle number density data given. Please use WSMs\n"
       << "*ParticleTypeInit* and *ParticleTypeAdd(All)* for reading\n"
       << "scattering element data.\n";
    throw runtime_error(os.str());
  }
  
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  ArrayOfIndex cloudbox_limits_tmp;
  /*if ( cloudbox_limits.nelem()== 0 )
    {
      //If no limits set, the cloud(box) is supposed to cover the
      //complete atmosphere. This particularly to facilitate use of
      //scat_data_single&pnd_fields for non-scatt, greybody cloud calculations.
      //Hence, set the limits to first & last elements of the different grids.
      //Note: no margin left at lat/lon_grid edges!.
      cloudbox_limits_tmp.resize(2*atmosphere_dim);

      // Pressure limits
      cloudbox_limits_tmp[0] = 0;
      cloudbox_limits_tmp[1] = p_grid.nelem() - 1;
      // Latitude limits
      if( atmosphere_dim >= 2 )
        {
          cloudbox_limits_tmp[2] = 0;
          cloudbox_limits_tmp[3] = lat_grid.nelem() - 1;
        }
      // Longitude limits
      if( atmosphere_dim == 3 )
        {
          cloudbox_limits_tmp[4] = 0;
          cloudbox_limits_tmp[5] = lon_grid.nelem() - 1;
        }
    }
  else */
  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");
  else
    cloudbox_limits_tmp = cloudbox_limits;

  // Check that pnd_field_raw has at least 2 grid-points in each dimension.
  // Otherwise, interpolation further down will fail with assertion.
  
  for (Index d = 0; d < atmosphere_dim; d++)
    {
      for (Index i = 0; i < pnd_field_raw.nelem(); i++)
        {
          if (pnd_field_raw[i].get_grid_size(d) < 2)
            {
              ostringstream os;
              os << "Error in pnd_field_raw data. ";
              os << "Dimension " << d << " (name: \"";
              os << pnd_field_raw[i].get_grid_name(d);
              os << "\") has only ";
              os << pnd_field_raw[i].get_grid_size(d);
              os << " element";
              os << ((pnd_field_raw[i].get_grid_size(d)==1) ? "" : "s");
              os << ". Must be at least 2.";
              throw runtime_error(os.str());
            }
        }
    }
  const Index Np_cloud = cloudbox_limits_tmp[1]-cloudbox_limits_tmp[0]+1;
  
  ConstVectorView p_grid_cloud = p_grid[Range(cloudbox_limits_tmp[0], Np_cloud)];

  // Check that no scatterers exist outside the cloudbox
  chk_pnd_field_raw_only_in_cloudbox(atmosphere_dim, pnd_field_raw,
                                     p_grid, lat_grid, lon_grid,
                                     cloudbox_limits_tmp);

  //==========================================================================
  if ( atmosphere_dim == 1)
    {
        ArrayOfGriddedField3 pnd_field_tmp;

        GriddedFieldPRegrid(pnd_field_tmp, p_grid_cloud, pnd_field_raw,
                            1, zeropadding, verbosity);

        FieldFromGriddedField(pnd_field,
                              p_grid_cloud,
                              pnd_field_tmp[0].get_numeric_grid(1),
                              pnd_field_tmp[0].get_numeric_grid(2),
                              pnd_field_tmp, verbosity);
    }
  else if(atmosphere_dim == 2)
    {
      const Index Nlat_cloud = cloudbox_limits_tmp[3]-cloudbox_limits_tmp[2]+1;

      ConstVectorView lat_grid_cloud = 
        lat_grid[Range(cloudbox_limits_tmp[2],Nlat_cloud)];           

      if (zeropadding)
        {
          // FIXME: OLE: Implement this
          CREATE_OUT0;
          out0 << "WARNING: zeropadding currently only supported for 1D.";
        }

      //Resize variables
      pnd_field.resize( pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 1 );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      
      // Interpolate pnd_field. 
      // Loop over the scattering element:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos( gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID),
                     p_grid_cloud);
          gridpos( gp_lat, pnd_field_raw[i].get_numeric_grid(GFIELD3_LAT_GRID),
                   lat_grid_cloud);
          
          // Interpolation weights:
          Tensor3 itw( Np_cloud, Nlat_cloud, 4 );
          interpweights( itw, gp_p, gp_lat );
          
          // Interpolate:
          interp( pnd_field(i,joker,joker,0), itw, 
                  pnd_field_raw[i].data(joker,joker,0), gp_p, gp_lat );
        }
    }
  else
    {
      const Index Nlat_cloud = cloudbox_limits_tmp[3]-cloudbox_limits_tmp[2]+1;
      const Index Nlon_cloud = cloudbox_limits_tmp[5]-cloudbox_limits_tmp[4]+1;

      if (zeropadding)
        {
          // FIXME: OLE: Implement this
          CREATE_OUT0;
          out0 << "WARNING: zeropadding currently only supported for 1D.";
        }

      ConstVectorView lat_grid_cloud =
        lat_grid[Range(cloudbox_limits_tmp[2],Nlat_cloud)];           
      ConstVectorView lon_grid_cloud = 
        lon_grid[Range(cloudbox_limits_tmp[4],Nlon_cloud)];
      
      //Resize variables
      pnd_field.resize( pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 
                        Nlon_cloud );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      ArrayOfGridPos gp_lon(Nlon_cloud);
      
      // Interpolate pnd_field. 
      // Loop over the scattering element types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos( gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID),
                     p_grid_cloud);
          gridpos( gp_lat, pnd_field_raw[i].get_numeric_grid(GFIELD3_LAT_GRID),
                   lat_grid_cloud);
          gridpos( gp_lon, pnd_field_raw[i].get_numeric_grid(GFIELD3_LON_GRID),
                   lon_grid_cloud);
          
          // Interpolation weights:
          Tensor4 itw( Np_cloud, Nlat_cloud, Nlon_cloud, 8 );
          interpweights( itw, gp_p, gp_lat, gp_lon );
          
          // Interpolate:
          interp( pnd_field(i,joker,joker,joker), itw, 
                  pnd_field_raw[i].data, gp_p, gp_lat, gp_lon );
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldExpand1D(Tensor4&        pnd_field,
                       const Index&    atmosphere_dim,
                       const Index&    cloudbox_on,    
                       const ArrayOfIndex&   cloudbox_limits,
                       const Index&    nzero,
                       const Verbosity&)
{
/*  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );
*/
  if( atmosphere_dim == 1 )
    { throw runtime_error( "No use in calling this method for 1D." ); }
  if( !cloudbox_on )
    { throw runtime_error( 
                "No use in calling this method with cloudbox off." ); }

  if( nzero < 1 )
    { throw runtime_error( "The argument *nzero must be > 0." ); }

  // Sizes
  const Index   npart = pnd_field.nbooks();
  const Index   np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index   nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
        Index   nlon = 1;
  if( atmosphere_dim == 3 )  
    { nlon = cloudbox_limits[5] - cloudbox_limits[4] + 1; }

  if( pnd_field.npages() != np  ||  pnd_field.nrows() != 1  ||  
      pnd_field.ncols() != 1 )
    { throw runtime_error( "The input *pnd_field* is either not 1D or does not "
                           "match pressure size of cloudbox." );}

  // Temporary container
  Tensor4 pnd_temp = pnd_field;

  // Resize and fill
  pnd_field.resize( npart, np, nlat, nlon );
  pnd_field = 0;
  //
  for( Index ilon=nzero; ilon<nlon-nzero; ilon++ )
    {
      for( Index ilat=nzero; ilat<nlat-nzero; ilat++ )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              for( Index is=0; is<npart; is++ )
                { pnd_field(is,ip,ilat,ilon) = pnd_temp(is,ip,0,0); }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldZero(//WS Output:
                   Tensor4& pnd_field,
                   ArrayOfArrayOfSingleScatteringData& scat_data,
                   //WS Input:
                   const Vector& p_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const Verbosity&)
{
  // 3D  atmosphere
  if (lat_grid.nelem()>1)
    {
      //Resize pnd_field and set it to 0:
      pnd_field.resize(1, p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem());
      pnd_field = 0.;
    }
  else // 1D atmosphere
     {
      //Resize pnd_field and set it to 0:
      pnd_field.resize(1, p_grid.nelem(), 1, 1);
      pnd_field = 0.;
     }
  
  //Resize scat_data and set it to 0:
  // Number of scattering elements
  scat_data.resize(1);
  scat_data[0].resize(1);
  scat_data[0][0].ptype = PTYPE_MACROS_ISO;
  scat_data[0][0].description = " ";
  // Grids which contain full ranges which one wants to calculate
  nlinspace(scat_data[0][0].f_grid, 1e9, 3.848043e+13 , 5);
  nlinspace(scat_data[0][0].T_grid, 0, 400, 5);
  nlinspace(scat_data[0][0].za_grid, 0, 180, 5);
  nlinspace(scat_data[0][0].aa_grid, 0, 360, 5);
  // Resize the data arrays
  scat_data[0][0].pha_mat_data.resize(5,5,5,1,1,1,6);
  scat_data[0][0].pha_mat_data = 0.;
  scat_data[0][0].ext_mat_data.resize(5,5,1,1,1);
  scat_data[0][0].ext_mat_data = 0.;
  scat_data[0][0].abs_vec_data.resize(5,5,1,1,1);
  scat_data[0][0].abs_vec_data = 0.;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFromscat_speciesFields (//WS Output:
                     Tensor4& pnd_field,
                     //WS Input:
                     const Index& atmosphere_dim,
                     const Index& cloudbox_on,
                     const ArrayOfIndex& cloudbox_limits,
                     const Tensor4& scat_species_mass_density_field,
                     const Tensor4& scat_species_mass_flux_field,
                     const Tensor4& scat_species_number_density_field,
                     const Tensor3& t_field,
                     const ArrayOfArrayOfScatteringMetaData& scat_meta,
                     const ArrayOfString& scat_species,
                     const String& delim,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;

  // Cloudbox on/off?
  if ( !cloudbox_on )
  {
    /* Must initialise pnd_field anyway; but empty */
    pnd_field.resize(0, 0, 0, 0);

    CREATE_OUT0;
    out0 << "  Cloudbox is off, pnd_field is set to empty.\n";
    return;
  }

  // Check consistency of scat_species and scat_data
  if ( scat_species.nelem() != scat_meta.nelem() )
    {
      ostringstream os;
      os << "Number of scattering species specified by scat_species does\n"
         << "not agree with number of scattering species in scat_meta:\n"
         << "scat_species has " << scat_species.nelem()
         << " entries, while scat_meta has " << scat_meta.nelem() << ".";
      throw runtime_error ( os.str() );
    }

  // ------- set pnd_field boundaries to cloudbox boundaries -------------------
  //set pnd_field boundaries
  ArrayOfIndex limits(6);
  //pressure
  limits[0] = cloudbox_limits[0];
  limits[1] = cloudbox_limits[1]+1;
  //latitude
  if ( atmosphere_dim > 1 )
  {
      limits[2] = cloudbox_limits[2];
      limits[3] = cloudbox_limits[3]+1;
  }
  else
  {
      limits[2] = 0;
      limits[3] = 1;
  }
  if ( atmosphere_dim > 2 )
  {
      limits[4] = cloudbox_limits[4];
      limits[5] = cloudbox_limits[5]+1;
  }
  else
  {
      limits[4] = 0;
      limits[5] = 1;
  }

  /* Do some checks. Not foolproof, but catches at least some. */

  if ((limits[1] > scat_species_mass_density_field.npages()) ||
      (limits[1] > t_field.npages()) ||
      (limits[3] > scat_species_mass_density_field.nrows()) ||
      (limits[3] > t_field.nrows()) ||
      (limits[5] > scat_species_mass_density_field.ncols()) ||
      (limits[5] > t_field.ncols()))
  {
    ostringstream os;
    os << "Cloudbox out of bounds compared to fields. "
       << "Upper limits: (p, lat, lon): "
       << "(" << limits[1] << ", " << limits[3] << ", " << limits[5] << "). "
       << "*scat_species_mass_density_field*: "
       << "(" << scat_species_mass_density_field.npages() << ", "
       << scat_species_mass_density_field.nrows() << ", "
       << scat_species_mass_density_field.ncols() << "). "
       << "*t_field*: "
       << "(" << t_field.npages() << ", "
       << t_field.nrows() << ", "
       << t_field.ncols() << ").";
    throw runtime_error(os.str());
  }

  //resize pnd_field to required atmospheric dimension and scattering elements
  pnd_field.resize ( TotalNumberOfElements(scat_meta), limits[1]-limits[0],
                     limits[3]-limits[2], limits[5]-limits[4] );
  ArrayOfIndex intarr;

  //-------- Start pnd_field calculations---------------------------------------

  // loop over nelem of scat_species
  for ( Index i_ss=0; i_ss<scat_species.nelem(); i_ss++ )
  {

    //cout << "for scat species #" << i_ss << " using PSD ";
    String psd_param;
    String partfield_name;
    String psd;

    //split String and copy to ArrayOfString
    parse_psd_param( psd_param, scat_species[i_ss], delim);
    parse_partfield_name( partfield_name, scat_species[i_ss], delim);

  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   * NOTE: when adding further distributions here, document them (particularly
   * the tags) in the WSM online documentation in methods.cc. Also, create a
   * wrapper WSM to the dNdD core method.
   * See ARTS Developer Guide for details.
   * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    //---- pnd_field calculations for MH97 -------------------------------
    if ( psd_param.substr(0,4) == "MH97" )
    {
        psd = "MH97";
        //cout << psd << "\n";

        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
                 << psd << " should should only be applied to cloud"
                 << " ice.\n";
        }

        pnd_fieldMH97( pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       t_field, limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }

    //---- pnd_field calculations for H11 ----------------------------
    else if ( psd_param == "H11" )
    {
        psd = "H11";
        //cout << psd << "\n";

        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
                 << psd << " should only be applied to cloud or precipitating"
                 << " ice.\n";
        }

        pnd_fieldH11 (pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       t_field, limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }
    
    //---- pnd_field calculations for H13 ----------------------------
    else if ( psd_param == "H13" )
    {
        psd = "H13";
        //cout << psd << "\n";

        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC"  )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
                 << psd << " should only be applied to cloud or precipitating"
                 << " ice.\n";
        }

        pnd_fieldH13 (pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       t_field, limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }
    
//---- pnd_field calculations for H13Shape ----------------------------
    else if ( psd_param == "H13Shape" )
    {
        psd = "H13Shape";
        //cout << psd << "\n";

        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC"  )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
                 << psd << " should only be applied to cloud or precipitating"
                 << " ice.\n";
        }

        pnd_fieldH13Shape (pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       t_field, limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }
      
    //---- pnd_field calculations for F07TR ----------------------------
    else if ( psd_param == "F07TR" )
    {
        psd = "F07TR";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to cloud or precipitating"
            << " ice.\n";
        }
        
        pnd_fieldF07TR (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      t_field, limits,
                      scat_meta, i_ss,
                      scat_species[i_ss], delim,
                      verbosity);
    }
      
      //---- pnd_field calculations for F07ML ----------------------------
    else if ( psd_param == "F07ML" )
    {
        psd = "F07ML";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW"
            && partfield_name != "Snow" && partfield_name != "SWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to cloud or precipitating"
            << " ice.\n";
        }
        
        pnd_fieldF07ML (pnd_field,
                        scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                        t_field, limits,
                        scat_meta, i_ss,
                        scat_species[i_ss], delim,
                        verbosity);
    }
      
      //---- pnd_field calculations for S2M_LWC ----------------------------
    else if ( psd_param == "S2M_LWC" )
    {
        psd = "S2M_LWC";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "LWC" && partfield_name != "CLW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to liquid cloud water.\n";
        }
        
        pnd_fieldS2M (pnd_field,
                        scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                        scat_species_number_density_field ( i_ss, joker, joker, joker ),
                        limits,
                        scat_meta,
                        i_ss,
                        scat_species[i_ss],
                        delim,
                        psd,
                        verbosity);
    }
      
      //---- pnd_field calculations for S2M_IWC ----------------------------
    else if ( psd_param == "S2M_IWC" )
    {
        psd = "S2M_IWC";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to cloud ice.\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      psd,
                      verbosity);
    }
    
      //---- pnd_field calculations for S2M_RWC ----------------------------
    else if ( psd_param == "S2M_RWC" )
    {
        psd = "S2M_RWC";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "Rain" && partfield_name != "RWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to precipitating liquid water\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      psd,
                      verbosity);
    }
      
      
      //---- pnd_field calculations for S2M_SWC ----------------------------
    else if ( psd_param == "S2M_SWC" )
    {
        psd = "S2M_SWC";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "Snow" && partfield_name != "SWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to snow\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      psd,
                      verbosity);
    }
      
      
      //---- pnd_field calculations for S2M_GWC ----------------------------
    else if ( psd_param == "S2M_GWC" )
    {
        psd = "S2M_GWC";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "Graupel" && partfield_name != "GWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to graupel\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      psd,
                      verbosity);
    }
      
      //---- pnd_field calculations for S2M_HWC ----------------------------
    else if ( psd_param == "S2M_HWC" )
    {
        psd = "S2M_HWC";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "Hail" && partfield_name != "HWC")
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to hail\n";
        }
        
        pnd_fieldS2M (pnd_field,
                      scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                      scat_species_number_density_field ( i_ss, joker, joker, joker ),
                      limits,
                      scat_meta,
                      i_ss,
                      scat_species[i_ss],
                      delim,
                      psd,
                      verbosity);
    }
      
      //---- pnd_field calculations for MGD_LWC ----------------------------
    else if ( psd_param == "MGD_LWC" )
    {
        psd = "MGD_LWC";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "LWC" && partfield_name != "CLW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to liquid cloud water.\n";
        }
        
        pnd_fieldMGD_LWC (pnd_field,
                        scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                        limits,
                        scat_meta, i_ss,
                        scat_species[i_ss], delim,
                        verbosity);
    }
      
      //---- pnd_field calculations for MGD_LWC ----------------------------
    else if ( psd_param == "MGD_IWC" )
    {
        psd = "MGD_IWC";
        //cout << psd << "\n";
        
        //check for expected scattering species field name
        if ( partfield_name != "IWC" && partfield_name != "CIW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
            << psd << " should only be applied to cloud ice.\n";
        }
        
        pnd_fieldMGD_IWC (pnd_field,
                          scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                          limits,
                          scat_meta, i_ss,
                          scat_species[i_ss], delim,
                          verbosity);
    }
      
    //---- pnd_field calculations for MP48 -------------------------------
    else if ( psd_param == "MP48" )
    {
        psd = "MP48";
        //cout << psd << "\n";

        //check for expected scattering species field name
        if ( partfield_name != "Rain" && partfield_name != "Snow"
            && partfield_name != "RR" && partfield_name != "SR" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
                 << psd << " should only be applied to liquid or frozen"
                 << " precipitation.\n";
        }

        pnd_fieldMP48( pnd_field,
                       scat_species_mass_flux_field ( i_ss, joker, joker, joker ),
                       limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }

    // ---- pnd_field calculations for liquid water clouds ---------------------
    //      (parameters from Hess98, Stratus continental)
    else if ( psd_param == "H98_STCO" || psd_param == "STCO" || psd_param == "liquid" )
    {
        psd = "H98_STCO";
        //cout << psd << "\n";

        //check for expected scattering species field name
        if ( partfield_name != "LWC" && partfield_name != "CLW" )
        {
            out2 << "WARNING! Unexpected scattering species field name ("
                 << partfield_name << ") for distribution " << psd << ".\n"
                 << psd << " should should only be applied to liquid or frozen"
                 << " precipitation.\n";
        }

        pnd_fieldH98 (pnd_field,
                       scat_species_mass_density_field ( i_ss, joker, joker, joker ),
                       limits,
                       scat_meta, i_ss,
                       scat_species[i_ss], delim,
                       verbosity);
    }

    else
    {
        ostringstream os;
        os << "Size distribution function '" << psd_param << "' is unknown!\n";
        throw runtime_error( os.str() );
    }

  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   * NOTE: when adding further distributions here, document them (particularly
   * the tags) in the WSM online documentation in methods.cc. Also, create a
   * wrapper WSM to the dNdD core method.
   * See ARTS Developer Guide for details.
   * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MH97 (//WS Output:
              Vector& dNdD,
              //WS Input:
              const Vector& diameter_mass_equivalent,
              const Numeric& IWC,
              const Numeric& t,
              const Index& noisy,
              const Verbosity&)
{
  Index n_se = diameter_mass_equivalent.nelem();
  dNdD.resize(n_se);

  for ( Index i=0; i<n_se; i++ )
    {
      // calculate particle size distribution with MH97
      // [# m^-3 m^-1]
      dNdD[i] = IWCtopnd_MH97 ( IWC, diameter_mass_equivalent[i], t, noisy );
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_H11 (//WS Output:
             Vector& dNdD,
             //WS Input:
             const Vector& Dmax,
             const Numeric& t,
             const Verbosity&)
{
  Index n_se = Dmax.nelem();
  dNdD.resize(n_se);

  for ( Index i=0; i<n_se; i++ ) //loop over number of scattering elementss
    {
      // calculate particle size distribution for H11
      // [# m^-3 m^-1]
      dNdD[i] = IWCtopnd_H11 ( Dmax[i], t );
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_Ar_H13 (//WS Output:
                Vector& dNdD,
                Vector& Ar,
                //WS Input:
                const Vector& Dmax,
                const Numeric& t,
                const Verbosity&)
{
  Index n_se = Dmax.nelem();
  dNdD.resize(n_se);
  Ar.resize(n_se);

  for ( Index i=0; i<n_se; i++ ) //loop over number of scattering elementss
    {
      // calculate particle size distribution for H13Shape
      // [# m^-3 m^-1]
      dNdD[i] = IWCtopnd_H13Shape ( Dmax[i], t );
      // calculate Area ratio distribution for H13Shape
      Ar[i] = area_ratioH13 ( Dmax[i], t );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_F07TR (//WS Output:
                Vector& dNdD,
                //WS Input:
                const Vector& diameter_max,
                const Numeric& SWC,
                const Numeric& t,
                const Numeric& alpha,
                const Numeric& beta,
                const Verbosity&)
{
    Index n_se = diameter_max.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with MH97
        // [# m^-3 m^-1]
        dNdD[i] = IWCtopnd_F07TR(diameter_max[i], t, SWC, alpha, beta) ;
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_F07ML (//WS Output:
                 Vector& dNdD,
                 //WS Input:
                 const Vector& diameter_max,
                 const Numeric& SWC,
                 const Numeric& t,
                 const Numeric& alpha,
                 const Numeric& beta,
                 const Verbosity&)
{
    Index n_se = diameter_max.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with MH97
        // [# m^-3 m^-1]
        dNdD[i] = IWCtopnd_F07ML(diameter_max[i], t, SWC, alpha, beta) ;
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_S2M (//WS Output:
                 Vector& dNdD,
                 //WS Input:
               const Vector& mass,
               const Numeric& N_tot,
               const Numeric& M,
               const String& psd_type,
               const Verbosity&)
{
    Index n_se = mass.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with S2M
        // [# m^-3 kg^-1]
        dNdD[i] = WCtopnd_S2M(mass[i], N_tot, M, psd_type) ;
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MGD_LWC (//WS Output:
                 Vector& dNdD,
                 //WS Input:
                 const Vector& diameter_volume_equ,
                 const Numeric& rho,
                 const Numeric& LWC,
                 const Verbosity&)
{
    Index n_se = diameter_volume_equ.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with MH97
        // [# m^-3 m^-1]
        dNdD[i] = LWCtopnd_MGD_LWC( diameter_volume_equ[i],rho ,LWC );
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MGD_IWC (//WS Output:
                   Vector& dNdD,
                   //WS Input:
                   const Vector& diameter_volume_equ,
                   const Numeric& rho,
                   const Numeric& IWC,
                   const Verbosity&)
{
    Index n_se = diameter_volume_equ.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with MH97
        // [# m^-3 m^-1]
        dNdD[i] = IWCtopnd_MGD_IWC( diameter_volume_equ[i],rho,IWC );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_H98 (//WS Output:
             Vector& dNdD,
             //WS Input:
             const Vector& diameter_volume_equivalent,
             const Numeric& LWC,
             const Verbosity&)
{
  Index n_se = diameter_volume_equivalent.nelem();
  dNdD.resize(n_se);
  const Numeric dDdR = 2.; 

  for ( Index i=0; i<n_se; i++ )
    {
      // calculate particle size distribution for liquid
      // and compensate for LWCtopnd providing dNdR
      // [# m^-3 m^-1]
      dNdD[i] = LWCtopnd ( LWC, diameter_volume_equivalent[i]/2. ) / dDdR;
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MP48 (//WS Output:
              Vector& dNdD,
              //WS Input:
              const Vector& diameter_melted_equivalent,
              const Numeric& PR,
              const Verbosity&)
{
  Index n_se = diameter_melted_equivalent.nelem();
  dNdD.resize(n_se);

  // derive particle number density for all given sizes
  for ( Index i=0; i<n_se; i++ )
    {
      // calculate particle size distribution with MP48
      // output: [# m^-3 m^-1]
      dNdD[i] = PRtopnd_MP48 ( PR, diameter_melted_equivalent[i]);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromdNdD (//WS Output:
              Vector& pnd,
              //WS Input:
              const Vector& dNdD,
              const Vector& diameter,
              const Numeric& total_content,
              const Vector& scatelem_content,
              const Verbosity& verbosity)
{
  pnd.resize( diameter.nelem() );

  // scale dNdD by size bin width
  if (diameter.nelem() > 1)
      scale_pnd( pnd, diameter, dNdD ); //[# m^-3]
  else
      pnd = dNdD;

  // scaling pnd to real mass/number density/flux (some PSDs have implicit
  // scaling - then this is only a check -, other don't
  chk_pndsum ( pnd, total_content, scatelem_content,
               0, 0, 0, "your field", verbosity );
}
