/* Copyright (C) 2002-2008
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

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Numeric PI;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxOff (
  // WS Output:
  Index&           cloudbox_on,
  ArrayOfIndex&    cloudbox_limits,
  Agenda&          iy_cloudbox_agenda )
{
  cloudbox_on = 0;
  cloudbox_limits.resize ( 0 );
  iy_cloudbox_agenda = Agenda();
  iy_cloudbox_agenda.set_name ( "iy_cloudbox_agenda" );
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetAutomatically (
  // WS Output:
  Index&          cloudbox_on,
  ArrayOfIndex&   cloudbox_limits,
  Agenda&  iy_cloudbox_agenda,
  // WS Input:
  const Index&    atmosphere_dim,
  const ArrayOfString&  part_species,
  const Vector&   p_grid,
  const Vector&   lat_grid,
  const Vector&   lon_grid,
  const Tensor4&  massdensity_field,
  // Control Parameters
  const Numeric&  cloudbox_margin
)
{
  // Variables
  Numeric p1 = massdensity_field.npages()-1;
  Numeric p2 = 0;
  Numeric lat1 = massdensity_field.nrows()-1;
  Numeric lat2 = 0;
  Numeric lon1 = massdensity_field.ncols()-1;
  Numeric lon2 = 0;
  Numeric p_margin1;

  Index nhyd=0, i=0, j=0, k=0, l=0;
  bool x = false; //initialize flag, telling if all selected *massdensity_fields* are zero(false) or not(true)


  // Check existing WSV
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid ); // includes p_grid chk_if_decresasing

  // Set cloudbox_on
  cloudbox_on = 1;

  // Allocate cloudbox_limits
  cloudbox_limits.resize ( atmosphere_dim*2 );

  //--------- Start loop over particles ---------------------------------------------------------------
  for ( l=0; l<part_species.nelem(); l++ )
  {
    String hydromet_type;

    //split String and copy to ArrayOfString
    // split part_species string at "-" and write to ArrayOfString
    parse_part_type(hydromet_type, part_species[l]);
    
    // select hydrometeor type according to user input
    // Index nhyd, describes the column index in *massdensity_field*
    if ( hydromet_type == "LWC" ) nhyd = 0;
    else if ( hydromet_type == "IWC" ) nhyd = 1;
    else if ( hydromet_type == "Rain" ) nhyd = 2;
    else if ( hydromet_type == "Snow" ) nhyd = 3;
   
    bool y; // empty_flag
    //y is set to true, if a single value of massdensity_field is unequal zero
    chk_massdensity_field ( y,
                         atmosphere_dim,
                         massdensity_field ( nhyd, joker, joker, joker ),
                         p_grid,
                         lat_grid,
                         lon_grid );

    //-----------Start setting cloudbox limits------------------------------------------------------
    if ( y )
    {
      //massdensity_field unequal zero -> x is true
      x = true;

      if ( atmosphere_dim == 1 )
      {
        // Pressure limits
        ConstVectorView hydro_p = massdensity_field ( nhyd, joker, 0 , 0 );


        // set lower cloudbox_limit to surface if margin = -1
        if ( cloudbox_margin == -1 )
        {
          cloudbox_limits[0] = 0;
          i = 0;
        }
        else
        {
          // find index of first pressure lvl where hydromet_field is unequal 0, starting from surface
          for ( i=0; i<hydro_p.nelem(); i++ )
          {
            if ( hydro_p[i] != 0.0 )
            {
              if ( p1 > i ) //check if p1 is the lowest index in all selected massdensity fields
              {
                p1 = i;
              }
              break;
            }
          }

        }
        // find index of highest pressure lvl where massdensity_field is unequal 0, starting from top of atm.
        for ( j=hydro_p.nelem()-1; j>=i; j-- )
        {
          if ( hydro_p[j] != 0.0 )
          {
            if ( p2 < j ) //check if p2 is the highest index in all selected massdensity fields
            {
	      // j+1 to ensure that linear interpolation of particle number densities is possible
              p2 = j+1;
            }
            break;
          }
        }
        //out0<<"\n"<<p1<<"\n"<<p_grid[p1]<<"\n";

        // alter lower cloudbox_limit by cloudbox_margin, using barometric height formula
        p_margin1 = barometric_heightformula ( p_grid[p1], cloudbox_margin );
        while ( p_grid[k] >= p_margin1 ) k++;
        cloudbox_limits[0]= k;

        //out0<<"\n"<<cloudbox_limits[0]<<"\n"<<p_margin1<<"\n";

        // set upper cloudbox_limit
	// if cloudbox reaches to the upper most pressure level it needs to be reduced 
	if (p2 >= hydro_p.nelem())
	{
	  
	  cloudbox_limits[1] = hydro_p.nelem()-1;
	  
	  out2<<"The cloudbox reaches to TOA!\n"
	      <<"The upper limit is lowered by one.\n"
	      <<"Check hydromet_field data, if realistic!\n";
	      
	}
	else cloudbox_limits[1] = p2;

        //out0<<"\n"<<p2<<"\n"<<p_grid[p2]<<"\n";
      }
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
          if ( p1 <= i )    p1 = i;
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
  // check if all selected massdensity fields are zero at each level switch cloudbox off,
  // skipping scattering calculations
  if ( !x )
  {
    cloudboxOff ( cloudbox_on, cloudbox_limits, iy_cloudbox_agenda );

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

  if ( atmosphere_dim >= 2 )
  {
    // The latitude in *lat2* must be bigger than the latitude in *lat1*.
    assert ( lat_grid[lat2] > lat_grid[lat1] );
    // The latitude in *lat1* must be >= than the second value in *lat_grid*.
    assert ( lat_grid[lat1] >= lat_grid[1] );
    // The latitude in *lat2* must be <= than the next to last value in *lat_grid*.
    assert ( lat_grid[lat2] <= lat_grid[lat_grid.nelem()-2] );
  }
  if ( atmosphere_dim == 3 )
  {
    // The longitude in *lon2* must be bigger than the longitude in *lon1*.
    assert ( lon_grid[lon2] > lon_grid[lon1] );
    // The longitude in *lon1* must be >= than the second value in *lon_grid*.
    assert ( lon_grid[lon1] >= lon_grid[1] );
    // The longitude in *lon2* must be <= than the next to last value in *lon_grid*.
    assert ( lon_grid[lon2] <= lon_grid[lon_grid.nelem()-2] );
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetManually(
        // WS Output:
        Index&          cloudbox_on,
        ArrayOfIndex&   cloudbox_limits,
        // WS Input:
        const Index&    atmosphere_dim,
        const Vector&   p_grid,
        const Vector&   lat_grid,
        const Vector&   lon_grid,
        // Control Parameters:
        const Numeric& p1,
        const Numeric& p2,
        const Numeric& lat1,
        const Numeric& lat2,
        const Numeric& lon1,
        const Numeric& lon2 )
{
  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Check keyword arguments
  if( p1 <= p2 )
    throw runtime_error( 
            "The pressure in *p1* must be bigger than the pressure in *p2*." );
  if( p1 <= p_grid[p_grid.nelem()-1] )
    throw runtime_error( "The pressure in *p1* must be larger than the "
                                                   "last value in *p_grid*." );
  if( p2 >= p_grid[0] )
    throw runtime_error( "The pressure in *p2* must be smaller than the "
                                                  "first value in *p_grid*." );
  if( atmosphere_dim >= 2 )
    {
      if( lat2 <= lat1 )
        throw runtime_error( 
         "The latitude in *lat2* must be bigger than the latitude in *lat1*.");
      if( lat1 < lat_grid[1] )
        throw runtime_error( "The latitude in *lat1* must be >= than the "
                                               "second value in *lat_grid*." );
      if( lat2 > lat_grid[lat_grid.nelem()-2] )
        throw runtime_error( "The latitude in *lat2* must be <= than the "
                                         "next to last value in *lat_grid*." );
    }
  if( atmosphere_dim == 3 )
    {
      if( lon2 <= lon1 )
        throw runtime_error( 
       "The longitude in *lon2* must be bigger than the longitude in *lon1*.");
      if( lon1 < lon_grid[1] )
        throw runtime_error( "The longitude in *lon1* must be >= than the "
                                               "second value in *lon_grid*." );
      if( lon2 > lon_grid[lon_grid.nelem()-2] )
        throw runtime_error( "The longitude in *lon2* must be <= than the "
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
void cloudboxSetManuallyAltitude(
        // WS Output:
        Index&          cloudbox_on,
        ArrayOfIndex&   cloudbox_limits,
        // WS Input:
        const Index&    atmosphere_dim,
        const Tensor3&  z_field,
        const Vector&   lat_grid,
        const Vector&   lon_grid,
        // Control Parameters:
        const Numeric& z1,
        const Numeric& z2,
        const Numeric& lat1,
        const Numeric& lat2,
        const Numeric& lon1,
        const Numeric& lon2 )
{
  // Check existing WSV
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  // Check keyword arguments
  if( z1 >= z2 )
    throw runtime_error( 
                        "The altitude in *z1* must be smaller than the altitude in *z2*." );
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
        throw runtime_error( 
                            "The latitude in *lat2* must be bigger than the latitude in *lat1*.");
      if( lat1 < lat_grid[1] )
        throw runtime_error( "The latitude in *lat1* must be >= than the "
                             "second value in *lat_grid*." );
      if( lat2 > lat_grid[lat_grid.nelem()-2] )
        throw runtime_error( "The latitude in *lat2* must be <= than the "
                             "next to last value in *lat_grid*." );
      if( lon2 <= lon1 )
        throw runtime_error
          ("The longitude in *lon2* must be bigger than the longitude in *lon1*.");
      if( lon1 < lon_grid[1] )
        throw runtime_error( "The longitude in *lon1* must be >= than the "
                             "second value in *lon_grid*." );
      if( lon2 > lon_grid[lon_grid.nelem()-2] )
        throw runtime_error( "The longitude in *lon2* must be <= than the "
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
void Massdensity_cleanup ( //WS Output:
  Tensor4& massdensity_field,
  //WS Input:
  const Numeric& massdensity_threshold
)
{
  // Check that hydromet_fields contain realistic values
  //(values smaller than hydromet_threshold will be set to 0)
  for ( Index i=0; i<massdensity_field.nbooks(); i++ )
  {
    for ( Index j=0; j<massdensity_field.npages(); j++ )
    {
      for ( Index k=0; k<massdensity_field.nrows(); k++ )
      {
        for ( Index l=0; l<massdensity_field.ncols(); l++ )
        {
          if ( massdensity_field ( i,j,k,l ) < massdensity_threshold ) massdensity_field ( i,j,k,l ) = 0.0;
        }
      }
    }
 }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleSpeciesInit ( ArrayOfString&  part_species )
{
  part_species.resize ( 0 );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleSpeciesSet ( // WS Generic Output:
  ArrayOfString&  part_species,
  // Control Parameters:
  const ArrayOfString& names
)
{
  part_species.resize ( names.nelem() );
  //assign input strings to part_species
  part_species = names;

  // Print list of particle settings to the most verbose output stream:
  out3 << "  Defined particle settings: ";
  for ( Index i=0; i<part_species.nelem(); ++i )
  {
    out3 << "\n  " << i << ": "<<part_species[i];

  }
  out3 << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeInit ( //WS Output:
  ArrayOfSingleScatteringData& scat_data_raw,
  ArrayOfGriddedField3& pnd_field_raw
)
{
  scat_data_raw.reserve ( 20 );
  pnd_field_raw.reserve ( 20 );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeAddAll ( //WS Output:
  ArrayOfSingleScatteringData& scat_data_raw,
  ArrayOfGriddedField3&  pnd_field_raw,
  // WS Input(needed for checking the datafiles):
  const Index& atmosphere_dim,
  const Vector& f_grid,
  const Vector& p_grid,
  const Vector& lat_grid,
  const Vector& lon_grid,
  const ArrayOfIndex& cloudbox_limits,
  // Keywords:
  const String& filename_scat_data,
  const String& pnd_field_file )
{
  //--- Check input ---------------------------------------------------------

  // Atmosphere
  chk_if_in_range ( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Cloudbox limits
  if ( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error (
      "*cloudbox_limits* is a vector which contains"
      "the upper and lower\n"
      "limit of the cloud for all atmospheric dimensions.\n"
      "So its length must be 2 x *atmosphere_dim*" );
  // Frequency grid
  if ( f_grid.nelem() == 0 )
    throw runtime_error ( "The frequency grid is empty." );
  chk_if_increasing ( "f_grid", f_grid );


  //--- Reading the data ---------------------------------------------------
  ArrayOfString data_files;
  xml_read_from_file ( filename_scat_data, data_files );
  scat_data_raw.resize ( data_files.nelem() );

  for ( Index i = 0; i<data_files.nelem(); i++ )
  {

    out2 << "  Read single scattering data\n";
    xml_read_from_file ( data_files[i],
                         scat_data_raw[i] );

    chk_single_scattering_data ( scat_data_raw[i],
                                 data_files[i], f_grid );

  }

  out2 << "  Read particle number density data \n";
  xml_read_from_file ( pnd_field_file, pnd_field_raw );

  chk_pnd_raw_data ( pnd_field_raw,
                     pnd_field_file, atmosphere_dim, p_grid, lat_grid,
                     lon_grid, cloudbox_limits );
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatteringParticleTypeAndMetaRead ( //WS Output:
  ArrayOfSingleScatteringData& scat_data_raw,
  ArrayOfScatteringMetaData& scat_data_meta_array,
  const Vector& f_grid,
  // Keywords:
  const String& filename_scat_data,
  const String& filename_scat_meta_data )
{
  //--- Reading the data ---------------------------------------------------
  ArrayOfString data_files;
  ArrayOfString meta_data_files;
  
  // single scattering data read to temporary ArrayOfSingleScatteringData
  xml_read_from_file ( filename_scat_data, data_files );
  scat_data_raw.resize ( data_files.nelem() );

  for ( Index i = 0; i<data_files.nelem(); i++ )
  {
    out3 << "  Read single scattering data\n";
    xml_read_from_file ( data_files[i],
                         scat_data_raw[i] );

    chk_single_scattering_data ( scat_data_raw[i],
                                 data_files[i], f_grid );

  }

  // scattering meta data read to temporary ArrayOfScatteringMetaData
  xml_read_from_file ( filename_scat_meta_data, meta_data_files );
  scat_data_meta_array.resize ( meta_data_files.nelem() );

  for ( Index i = 0; i<meta_data_files.nelem(); i++ )
  {

    out3 << "  Read scattering meta data\n";
    xml_read_from_file ( meta_data_files[i],
                         scat_data_meta_array[i] );

    chk_scattering_meta_data ( scat_data_meta_array[i],
                               meta_data_files[i] );

  }
  
  // check if arrays have same size
  chk_scattering_data ( scat_data_raw,
                        scat_data_meta_array );
  
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatteringParticlesSelect (//WS Output:
  ArrayOfSingleScatteringData& scat_data_raw,
  ArrayOfScatteringMetaData& scat_data_meta_array,
  ArrayOfIndex& scat_data_nelem,
  // WS Input:
  const ArrayOfString& part_species)
{
  //--- Adjusting data to user specified input (part_species)---------------------------------
  // 
  String type;
  Numeric intarr_total = 0;
  ArrayOfIndex intarr;
  
  //make temporary copy
  ArrayOfSingleScatteringData scat_data_raw_tmp = scat_data_raw;
  ArrayOfScatteringMetaData scat_data_meta_array_tmp = scat_data_meta_array;
  
  scat_data_nelem.resize( part_species.nelem() );
  
  // loop over array of part_species
  for ( Index k=0; k<part_species.nelem(); k++ )
  {
   
    String hydromet_type;
    Numeric sizemin;
    Numeric sizemax;

    //split part_species string and copy values to parameter
    parse_part_type(hydromet_type, part_species[k]);
    // set type according to *part_species* input
    if ( hydromet_type == "IWC" || hydromet_type== "Snow" ) type = "Ice";
    else if ( hydromet_type== "LWC" || hydromet_type == "Rain" ) type = "Water";
    
    //split part_species string and copy values to parameter
    parse_part_size(sizemin, sizemax, part_species[k]);
    

    // choosing the specified SingleScatteringData and ScatteringMetaData
    for ( Index j=0; j<scat_data_meta_array_tmp.nelem(); j++ )
    {
      // check for particle phase type (e.g. "Ice", "Water",...)
      if ( scat_data_meta_array_tmp[j].type == type ) 
      {       
        //particle radius is calculated from particle volume given in scattering meta data
        Numeric r_particle = pow (    3./4. * scat_data_meta_array_tmp[j].V *1e18 /PI ,  1./3.  );
	
	//check if particle is in size range
        if ( r_particle  >= sizemin && sizemax >= r_particle )
	{
	  // fill ArrayOfIndex with indices of selected scattering data
          intarr.push_back ( j );
        }
      }
    }
    // WSV scat_data_nelem gets the number of elements of scattering data
    // connected to each selection String in *part_species*   
    scat_data_nelem[k] = intarr.nelem() - intarr_total;
    intarr_total = intarr.nelem();
  }
  // check if array is empty
  if ( !intarr.nelem() )
  {
    ostringstream os;
    os<<"The selection in " << part_species << " is NOT choosing any of the given Scattering Data.\n"
    << "--> Does the selection in *part_species* fit any of the Single Scattering Data input? \n";
    throw runtime_error ( os.str() );
  }


  // resize WSVs to size of intarr
  scat_data_raw.resize ( intarr.nelem() );
  scat_data_meta_array.resize ( intarr.nelem() );

  for ( Index j=0; j<intarr.nelem(); j++ )
  {
    //append to WSV Arrays
    scat_data_meta_array[j] = scat_data_meta_array_tmp[intarr[j]] ;
    scat_data_raw[j] = scat_data_raw_tmp[intarr[j]] ;
  }


}


/* Workspace method: Doxygen documentation will be auto-generated */
void ParticleTypeAdd( //WS Output:
                 ArrayOfSingleScatteringData& scat_data_raw,
                 ArrayOfGriddedField3&  pnd_field_raw,
                 // WS Input (needed for checking the datafiles):
                 const Index& atmosphere_dim,
                 const Vector& f_grid,
                 const Vector& p_grid,
                 const Vector& lat_grid,
                 const Vector& lon_grid,
                 const ArrayOfIndex& cloudbox_limits,
                 // Keywords:
                 const String& scat_data_file,
                 const String& pnd_field_file
          )
{
  //--- Check input ---------------------------------------------------------
  
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Cloudbox limits
  if ( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains"
                        "the upper and lower\n"
                        "limit of the cloud for all atmospheric dimensions.\n"
                        "So its length must be 2 x *atmosphere_dim*" ); 
  // Frequency grid
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  

  //--- Reading the data ---------------------------------------------------

  // Append *scat_data_raw* and *pnd_field_raw* with empty Arrays of Tensors. 
  SingleScatteringData scat_data;
  scat_data_raw.push_back(scat_data);
  
  GriddedField3 pnd_field_data;
  pnd_field_raw.push_back(pnd_field_data);
  
  out2 << "  Read single scattering data\n";
  xml_read_from_file(scat_data_file, scat_data_raw[scat_data_raw.nelem()-1]);

  chk_single_scattering_data(scat_data_raw[scat_data_raw.nelem()-1],
                             scat_data_file, f_grid);       
  
  out2 << "  Read particle number density field\n";
  if (pnd_field_file.nelem() < 1)
    out1 << "Warning: No pnd_field_file specified. Ignored. \n";
  else
    {
      xml_read_from_file(pnd_field_file,
                         pnd_field_raw[pnd_field_raw.nelem()-1]);
      
      chk_pnd_data(pnd_field_raw[pnd_field_raw.nelem()-1],
                   pnd_field_file, atmosphere_dim, p_grid, lat_grid,
                   lon_grid, cloudbox_limits);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalc(//WS Output:
                   Tensor4& pnd_field,
                   //WS Input
                   const Vector& p_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const ArrayOfGriddedField3& pnd_field_raw,
                   const Index& atmosphere_dim,
                   const ArrayOfIndex& cloudbox_limits
                   )
{
  // Basic checks of input variables
  //
  // Particle number density data
  // 
  if (pnd_field_raw.nelem() == 0)
    throw runtime_error(
                        "No particle number density data given. Please\n"
                        "use WSMs *ParticleTypeInit* and \n"
                        "*ParticleTypeAdd(All)* for reading cloud particle\n"
                        "data.\n"
                        );
  
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");

  const Index Np_cloud = cloudbox_limits[1]-cloudbox_limits[0]+1;
  
  ConstVectorView p_grid_cloud = p_grid[Range(cloudbox_limits[0], Np_cloud)];
      
  //==========================================================================
  if ( atmosphere_dim == 1)
    {
      //Resize variables
      pnd_field.resize(pnd_field_raw.nelem(), Np_cloud, 1, 1 );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
         
      // Interpolate pnd_field. 
      // Loop over the particle types:
      for (Index i = 0; i < pnd_field_raw.nelem(); ++ i)
        {
          // Calculate grid positions:
          p2gridpos( gp_p, pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID), 
                     p_grid_cloud );

          // Interpolation weights:
          Matrix itw(Np_cloud, 2);

          interpweights( itw, gp_p);
          // Interpolate:
          interp( pnd_field(i,joker,0,0), itw, 
                  pnd_field_raw[i].data(joker,0,0), gp_p );
        }
    }

  else if(atmosphere_dim == 2)
    {
      const Index Nlat_cloud = cloudbox_limits[3]-cloudbox_limits[2]+1;

      ConstVectorView lat_grid_cloud = 
        lat_grid[Range(cloudbox_limits[2],Nlat_cloud)];           
      
      //Resize variables
      pnd_field.resize( pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 1 );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      
      // Interpolate pnd_field. 
      // Loop over the particle types:
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
      const Index Nlat_cloud = cloudbox_limits[3]-cloudbox_limits[2]+1;
      const Index Nlon_cloud = cloudbox_limits[5]-cloudbox_limits[4]+1;

      ConstVectorView lat_grid_cloud = 
        lat_grid[Range(cloudbox_limits[2],Nlat_cloud)];           
      ConstVectorView lon_grid_cloud = 
        lon_grid[Range(cloudbox_limits[4],Nlon_cloud)];
      
      //Resize variables
      pnd_field.resize( pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 
                        Nlon_cloud );
      
      // Gridpositions:
      ArrayOfGridPos gp_p(Np_cloud);
      ArrayOfGridPos gp_lat(Nlat_cloud);
      ArrayOfGridPos gp_lon(Nlon_cloud);
      
      // Interpolate pnd_field. 
      // Loop over the particle types:
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
void pnd_fieldExpand1D(
            Tensor4&        pnd_field,
      const Vector&         p_grid,
      const Vector&         lat_grid,
      const Vector&         lon_grid,
      const Index&          atmosphere_dim,
      const Index&          cloudbox_on,    
      const ArrayOfIndex&   cloudbox_limits,
      const Index&          nzero )
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_cloudbox( atmosphere_dim, p_grid, lat_grid, lon_grid,
                                                cloudbox_on, cloudbox_limits );

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
    { nlon = cloudbox_limits[5] - cloudbox_limits[5] + 1; }

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
                      ArrayOfSingleScatteringData& scat_data_raw,
                      //WS Input:
                      const Vector& p_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid)
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
  
  //Resize scat_data_raw and set it to 0:
  // Number iof particle types
  scat_data_raw.resize(1);
  scat_data_raw[0].ptype = PARTICLE_TYPE_MACROS_ISO;
  scat_data_raw[0].description = " ";
  // Grids which contain full ranges which one wants to calculate
  nlinspace(scat_data_raw[0].f_grid, 1e9, 3.848043e+13 , 5);  
  nlinspace(scat_data_raw[0].T_grid, 0, 400, 5);
  nlinspace(scat_data_raw[0].za_grid, 0, 180, 5);
  nlinspace(scat_data_raw[0].aa_grid, 0, 360, 5);
  // Resize the data arrays
  scat_data_raw[0].pha_mat_data.resize(5,5,5,1,1,1,6);
  scat_data_raw[0].pha_mat_data = 0.;
  scat_data_raw[0].ext_mat_data.resize(5,5,1,1,1);
  scat_data_raw[0].ext_mat_data = 0.;
  scat_data_raw[0].abs_vec_data.resize(5,5,1,1,1);
  scat_data_raw[0].abs_vec_data = 0.;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldSetup ( //WS Output:
  Tensor4& pnd_field,
  //WS Input:
  const Index& atmosphere_dim,
  const Index& cloudbox_on,
  const ArrayOfIndex& cloudbox_limits,
  const Tensor4& massdensity_field,
  const Tensor3& t_field,
  const ArrayOfScatteringMetaData& scat_data_meta_array,
  const ArrayOfString& part_species,
  const ArrayOfIndex& scat_data_nelem
)
{
  // Cloudbox on/off?
  if ( !cloudbox_on ) return;



  //set iteration limits to cloudbox boundaries according to atmosphere_dim
  //initialize iteration boundaries
  Index p_cbstart = 0;
  Index p_cbend = 1;
  Index lat_cbstart = 0;
  Index lat_cbend = 1;
  Index lon_cbstart = 0;
  Index lon_cbend = 1;
  //pressure
  p_cbstart = cloudbox_limits[0];
  p_cbend = cloudbox_limits[1]+1;

  //latitude
  if ( atmosphere_dim >= 2 )
  {
    lat_cbstart = cloudbox_limits[2];
    lat_cbend = cloudbox_limits[3]+1;
  }
  //longitude
  if ( atmosphere_dim == 3 )
  {
    lon_cbstart = cloudbox_limits[4];
    lon_cbend = cloudbox_limits[5]+1;

  }

//resize pnd_field to number of size bins and atm. dimensions
  pnd_field.resize ( scat_data_meta_array.nelem(), p_cbend-p_cbstart, lat_cbend-lat_cbstart, lon_cbend-lon_cbstart );
  Index scat_data_start = 0;
  ArrayOfIndex intarr;

  //-------- Start pnd_field calculations----------------------------------------------------------

  // loop over nelem of part_species
  for ( Index k=0; k<part_species.nelem(); k++ )
  {

    String psd_param;

    //split String and copy to ArrayOfString
    parse_psd_param( psd_param, part_species[k]);

    // initialize control parameters
    Vector vol_unsorted ( scat_data_nelem[k], 0.0 );
    Vector vol ( scat_data_nelem[k], 0.0 );
    Vector dm ( scat_data_nelem[k], 0.0 );
    Vector r ( scat_data_nelem[k], 0.0 );
    Vector rho ( scat_data_nelem[k], 0.0 );
    Vector pnd ( scat_data_nelem[k], 0.0 );
    Vector pnd2 ( scat_data_nelem[k], 0.0 ); //temporary
    Vector dN ( scat_data_nelem[k], 0.0 );
    Vector dN2 ( scat_data_nelem[k], 0.0 ); //temporary
    Vector dlwc ( scat_data_nelem[k], 0.0 ); //temporary


// start pnd_field calculations for MH97
    if ( psd_param == "MH97" )
    {
      
      for ( Index i=0; i < scat_data_nelem[k]; i++ )
      {
	vol_unsorted[i] = ( scat_data_meta_array[i+scat_data_start].V ); //m^3
      }
      get_sorted_indexes(intarr, vol_unsorted);
      cout<<"intarr\t"<<intarr<<endl;
	
      
      // extract IWC_field and convert from kg/m^3 to g/m^3
      Tensor3 IWC_field = massdensity_field ( 1, joker, joker, joker );
      IWC_field*=1000; //IWC [g/m^3]
      //out0<<"\n"<<IWC_field<<"\n";

      // extract scattering meta data
      //Index j = 0; //counter for parameter assignment inside loop
      for ( Index i=0; i< scat_data_nelem[k]; i++ )
      {
        vol[i] = ( scat_data_meta_array[intarr[i]+scat_data_start].V ); //m^3
        // calculate melted diameter from volume [m]
        dm[i] = pow ( ( ( 6*scat_data_meta_array[intarr[i]+scat_data_start].V ) /PI ), ( 1./3. ) );
	// get density from meta data [g/m^3]
        rho[i] = scat_data_meta_array[intarr[i]+scat_data_start].density * 1000;

        //check for correct particle type
        //if (scat_data_meta_array[i].type != strarr[0])
        if ( scat_data_meta_array[intarr[i]+scat_data_start].type != "Ice" )
        {
          throw runtime_error ( "The particle phase is unequal 'Ice'.\n"
                                "MH97 can only be applied to ice paricles.\n"
                                "Check ScatteringMetaData!" );
        }
        //j++; // increment counter
      }
      // test vector
      //linspace(dm, 2, 2000, 1000);
      cout<<"dm\t"<<dm<<endl;
      

      // itertation over all atm. levels
      for ( Index p=p_cbstart; p<p_cbend; p++ )
      {
        for ( Index lat=lat_cbstart; lat<lat_cbend; lat++ )
        {
          for ( Index lon=lon_cbstart; lon<lon_cbend; lon++ )
          {
            // iteration over all given size bins
            for ( Index i=0; i<dm.nelem(); i++ )
            {
              // calculate particle size distribution with MH97
              dN[i] = IWCtopnd_MH97 ( IWC_field ( p, lat, lon ), dm[i], t_field ( p, lat, lon ), rho[i] );// [# m^-3 m^-1]
	      //dN2[i] = dN[i] * vol[i] * rho[i];
            }

            //out0<<dN<<"\n";
            
            
            // scale pnds by bin width
            scale_pnd( pnd, dm, dN );
	    //scale_pnd( pnd2, dm, dN2 );
            // calculate error of pnd sum and real XWC
            chk_pndsum ( pnd, IWC_field ( p,lat,lon ), vol, rho );
            //chk_pndsum2 (pnd2, IWC_field ( p,lat,lon ));

            // writing pnd vector to wsv pnd_field
            //Index q = 0; //counter for parameter assignment inside loop
            for ( Index i = 0; i< scat_data_nelem[k]; i++ )
            {
              pnd_field ( intarr[i]+scat_data_start, p-p_cbstart, lat-lat_cbstart, lon-lon_cbstart ) = pnd[i];

             // q++; // increment counter
            }

          }
        }
      }

    }
    // start pnd_field calculations for liquid
    else if ( psd_param == "liquid" )
    {
      for ( Index i=0; i < scat_data_nelem[k]; i++ )
      {
	vol_unsorted[i] = ( scat_data_meta_array[i+scat_data_start].V ); //m^3
      }
      get_sorted_indexes(intarr, vol_unsorted);
      cout<<"intarr\t"<<intarr<<endl;
      
      // extract LWC_field and convert from kg/m^3 to g/m^3
      Tensor3 LWC_field = massdensity_field ( 0, joker, joker, joker );
      LWC_field *= 1000; //LWC [g/m^3]

      // extract scattering meta data
      //Index j = 0; //counter for parameter assignment inside loop
      for ( Index i=0; i< scat_data_nelem[k]; i++ )
      {
        vol[i]= scat_data_meta_array[intarr[i]+scat_data_start].V; //m^3
        // calculate diameter from volume [m]
        dm[i] = pow ( ( 6*scat_data_meta_array[intarr[i]+scat_data_start].V/PI ), ( 1./3. ) );
        // diameter to radius
        r[i] = dm[i]/2; // [m]
        // get density from meta data [g/m^3]
        rho[i] = scat_data_meta_array[intarr[i]+scat_data_start].density * 1000; // get density from meta data [g/m^3]

        //check for correct particle type
        if ( scat_data_meta_array[intarr[i]+scat_data_start].type != "Water" )
        {
          throw runtime_error ( "The particle phase is unequal 'Water'.\n"
                                "All particles must be of liquid phase to apply this PSD.\n"
                                "Check ScatteringMetaData!" );
        }
        //j++; // increment counter
      }
      cout<<"\nr\t"<<r<<endl;

      // itertation over all atm. levels
      for ( Index p=p_cbstart; p<p_cbend; p++ )
      {
        for ( Index lat=lat_cbstart; lat<lat_cbend; lat++ )
        {
          for ( Index lon=lon_cbstart; lon<lon_cbend; lon++ )
          {
            // iteration over all given size bins
            for ( Index i=0; i<r.nelem(); i++ ) //loop over number of particles
            {
              // calculate particle size distribution for liquid
              dN[i] = LWCtopnd ( LWC_field ( p,lat,lon ), r[i] ); // [# m^-3 m^-1]
              //dN2[i] = LWCtopnd2 ( r[i] );  // [# m^-3 m^-1]
	      //dN2[i] = dN[i] * vol[i] * rho[i];
            }

            //dlwc *= LWC_field(p, lat, lon)/dlwc.sum();

            //out0<<"\n"<<dN<<"\n"<<dN2<<"\n";

            // scale pnds by scale width
            scale_pnd( pnd, r, dN ); //[# m^-3]
	    //scale_pnd( pnd2, r, dN2 );
            //trapezoid_integrate ( pnd2, r, dN2 );//[# m^-3]
            //out0<<"\n"<<"HIER!"<<"\n"<<pnd<<"\n"<<pnd2<<"\n";
            // calculate error of pnd sum and real XWC
            chk_pndsum ( pnd, LWC_field ( p,lat,lon ), vol, rho );
	    //chk_pndsum2 (pnd2, LWC_field ( p,lat,lon ));
            //chk_pndsum (pnd, testiwc[p], vol, rho);


            // writing pnd vector to wsv pnd_field
            //Index q = 0; //counter for parameter assignment inside loop
            for ( Index i =0; i< scat_data_nelem[k]; i++ )
            {
              pnd_field ( intarr[i]+scat_data_start, p-p_cbstart, lat-lat_cbstart, lon-lon_cbstart ) = pnd[i];
              //dlwc[q] = pnd2[q]*vol[q]*rho[q];

              //q++; // increment counter
            }

            //out0<<"\n"<<LWC_field(p, lat, lon)<<"\n"<< dlwc.sum()<<"\n";

            //pnd2 *= ( LWC_field ( p, lat, lon ) /dlwc.sum() );

            //out0<<"\n"<<pnd<<"\n"<<pnd2<<"\n";

          }
        }
      }
    }

    // alter starting index of current scattering data array to starting index of next iteration step
    scat_data_start = scat_data_start + scat_data_nelem[k];

  }
}


