/* Copyright (C) 2002,2003
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   Stefan Buehler   <sbuehler@uni-bremen.de>
   Claudia Emde     <claudia@sat.physik.uni-bremen.de>
                            
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

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/



//! cloudboxOff
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void cloudboxOff(
        // WS Output:
        Index&           cloudbox_on,
        ArrayOfIndex&    cloudbox_limits,
        Tensor7&         scat_i_p,
        Tensor7&         scat_i_lat,
        Tensor7&         scat_i_lon,
        Vector&          scat_za_grid,
        Vector&          scat_aa_grid )
{
  cloudbox_on = 0;
  cloudbox_limits.resize(0);
  scat_i_p.resize(0,0,0,0,0,0,0);
  scat_i_lat.resize(0,0,0,0,0,0,0);
  scat_i_lon.resize(0,0,0,0,0,0,0);
  scat_za_grid.resize(0);
  scat_aa_grid.resize(0);
}

//! cloudboxSetIncomingForTauCalc1D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-10-15
*/
void CloudboxSetIncomingForTauCalc1D(// WS Output:
                                     Tensor7& scat_i_p,
                                     Tensor7& scat_i_lat,
                                     Tensor7& scat_i_lon,
                                     // WS Input:
                                     const ArrayOfIndex&  cloudbox_limits,
                                     const Index& atmosphere_dim,
                                     const Index& stokes_dim,
                                     const Vector&  scat_za_grid,
                                     const Vector&  f_grid )
{
  // Check the atmosphere dimensionality and that only the first Stokes
  // component is used(?)
  if( atmosphere_dim != 1 )
    throw runtime_error( "This function is only defined for 1D "
      "atmosphere.\n");
  if( stokes_dim != 1 )
    throw runtime_error( "Only the first Stokes component must be set.\n");

  out2 <<"Setting clearsky field on cloudbox boundary for transmission "
       <<"calculation.\n" ;

  // Get sizes and lengths
  Index Nf = f_grid.nelem();
  Index Nza = scat_za_grid.nelem();

  // Resize the output variables
  scat_i_p.resize(Nf, 2, 1, 1, Nza, 1, 1);

  // Calculate the clearsky field on the lower side. Calculations are made for
  // zenith angles over 90 degrees.
  for( Index k=0; k<Nza; k++ ) {
    if( scat_za_grid[k]>90 ) {
      scat_i_p(Range(joker),0,0,0,k,0,0)=1/-cos(scat_za_grid[k]*DEG2RAD);
    } else {
      scat_i_p(Range(joker),0,0,0,k,0,0)=0.0;
    }
  }

  // Set the clearsky field on the upper side to zero.
  scat_i_p( Range(joker), 1, 0,0, Range(joker),0, 0) = 0.0;

  // Set incoming radiation at vertical cloudbox limits to empty vectors.
  scat_i_lat.resize(0,0,0,0,0,0,0);
  scat_i_lon.resize(0,0,0,0,0,0,0);

}

//! cloudboxSetManually
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-19
*/
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

  if( atmosphere_dim == 2 )
    { throw runtime_error( "The cloud box is not defined for 2D." ); }

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


//! This method interpolates clear sky field on the cloudbox boundary 
//on all grid points inside the cloud box. 

/*! 
  This method uses a linear 3D interpolation scheme to obtain the 
  radiation field on all grid points inside the cloud box form the 
  clear sky field on the cloud bod boundary.
  
  \param i_field Output : Intensity field
  \param scat_i_p Input : Intensity field on cloudbox boundary 
  (equal pressure surfaces)
  \param scat_i_lat Input : Intensity field on cloudbox boundary 
  (equal latitude surfaces)
  \param scat_i_lon Input : Intensity field on cloudbox boundary
  (equal longitude surfaces)
  \param f_grid Input : frequency grid
  \param f_index Input : the frequency index for scattering calculation
  \param p_grid Input : the pressure grid
  \param lat_grid Input : the latitude grid
  \param lon_grid Input : the longitude grid
  \param cloudbox_limits Input : Limits of the cloud box
  \param atmosphere_dim Input : dimension of atmosphere
  
  \author Sreerekha T.R.
   \date   2002-07-24
*/
void i_fieldSetClearsky(Tensor6& i_field,
                const Tensor7& scat_i_p,
                const Tensor7& scat_i_lat,
                const Tensor7& scat_i_lon,
                const Vector& f_grid,
                const Index& f_index,
                const Vector& p_grid,
                const Vector& lat_grid,
                const Vector& lon_grid,
                const ArrayOfIndex& cloudbox_limits,
                const Index& atmosphere_dim )
{
  
  out2 << "Interpolate boundary clearsky field to obtain the initial field.\n";

  if(atmosphere_dim == 1)
    {
      Index  N_f = scat_i_p.nlibraries();
      if (f_grid.nelem() != N_f){

        throw runtime_error(" scat_i_p should have same frequency  "
                            " dimension as f_grid");
      }
     
      if(scat_i_p.nvitrines() != 2){
        throw runtime_error("scat_i_p should have only two elements "
                            "in pressure grid which corresponds "
                            "to the two pressure surfaces");
      }
      
     
      Index N_za = scat_i_p.npages() ;
   
      Index N_aa = scat_i_p.nrows();
   
      Index N_i = scat_i_p.ncols();
         
      //1. interpolation - pressure grid
      
      
      i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1,
                     1,
                     1,
                     N_za,
                     N_aa,
                     N_i);

      i_field = 0.;

      

      /*the old grid is having only two elements, corresponding to the 
        cloudbox_limits and the new grid have elements corresponding to
        all grid points inside the cloudbox plus the cloud_box_limits*/

      ArrayOfGridPos p_gp((cloudbox_limits[1]- cloudbox_limits[0])+1);
      
      p2gridpos(p_gp,
              p_grid[Range(cloudbox_limits[0], 
                           2,
                           (cloudbox_limits[1]- cloudbox_limits[0]))],
              p_grid[Range(cloudbox_limits[0], 
                           (cloudbox_limits[1]- cloudbox_limits[0])+1)]);
      
      Matrix itw((cloudbox_limits[1]- cloudbox_limits[0])+1, 2);
      interpweights ( itw, p_gp );

     
 
      for (Index za_index = 0; za_index < N_za ; ++ za_index)
        {
          for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
            {
              for (Index i = 0 ; i < N_i ; ++ i)
                {
                  
                  VectorView target_field = i_field(Range(joker),
                                                    0,
                                                    0,
                                                    za_index,
                                                    aa_index,
                                                    i);
                  
                  ConstVectorView source_field = scat_i_p(f_index,
                                                          Range(joker),    
                                                          0,
                                                          0,
                                                          za_index,
                                                          aa_index,
                                                          i);
                  
                  interp(target_field,
                         itw,
                         source_field,
                         p_gp);
                }
              
            }
        }

    }
  if(atmosphere_dim == 3)
    {
      Index  N_f = scat_i_p.nlibraries();
      if (scat_i_lat.nlibraries() != N_f || 
          scat_i_lon.nlibraries() != N_f){
        
        throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
                            "same frequency dimension");
      }
      Index N_p = cloudbox_limits[1] - cloudbox_limits[0] + 1;
      if(scat_i_lon.nvitrines() != N_p ||
         scat_i_lat.nvitrines() != N_p ){
        throw runtime_error("scat_i_lat and scat_i_lon should have  "
                            "same pressure grid dimension as p_grid");
      }
      
      Index N_lat =  cloudbox_limits[3] - cloudbox_limits[2] + 1;
      
      if(scat_i_lon.nshelves() != N_lat ||
         scat_i_p.nshelves()   != N_lat){
        throw runtime_error("scat_i_p and scat_i_lon should have  "
                            "same latitude grid dimension as lat_grid");
      }
  
      Index N_lon = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      if(scat_i_lat.nbooks() != N_lon ||
         scat_i_p.nbooks()   != N_lon ){
        throw runtime_error("scat_i_p and scat_i_lat should have  "
                            "same longitude grid dimension as lon_grid");
      }
      if(scat_i_p.nvitrines() != 2){
        throw runtime_error("scat_i_p should have only two elements "
                            "in pressure grid which corresponds "
                            "to the two pressure surfaces");
      }
      
      if(scat_i_lat.nshelves() != 2){
        throw runtime_error("scat_i_lat should have only two elements "
                            "in latitude grid which corresponds "
                            "to the two latitude surfaces");
        
      }
      if(scat_i_lon.nbooks() != 2){
        throw runtime_error("scat_i_lon should have only two elements "
                            "in longitude grid which corresponds "
                            "to the two longitude surfaces");
        
      }
      Index N_za = scat_i_p.npages() ;
      if (scat_i_lat.npages() != N_za || 
          scat_i_lon.npages() != N_za){
        
        throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
                            "same dimension for zenith angles");
      }
      Index N_aa = scat_i_p.nrows();
      if (scat_i_lat.nrows() != N_aa || 
          scat_i_lon.nrows() != N_aa){
        
        throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
                            "same dimension for azimuth angles");
      }
      Index N_i = scat_i_p.ncols();
      if (scat_i_lat.ncols() != N_i || 
          scat_i_lon.ncols() != N_i){
        
        throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
                            "same value for stokes_dim and can take only"
                            "values 1,2,3 or 4");
      }
      
      //1. interpolation - pressure grid, latitude grid and longitude grid
      
   
      //i_field
      i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1, 
                     (cloudbox_limits[3]- cloudbox_limits[2])+1,
                     (cloudbox_limits[5]- cloudbox_limits[4])+1,
                     N_za, 
                     N_aa,
                     N_i);
      

      
      ArrayOfGridPos p_gp((cloudbox_limits[1]- cloudbox_limits[0])+1);
      ArrayOfGridPos lat_gp((cloudbox_limits[3]- cloudbox_limits[2])+1);
      ArrayOfGridPos lon_gp((cloudbox_limits[5]- cloudbox_limits[4])+1);

      /*the old grid is having only two elements, corresponding to the 
        cloudbox_limits and the new grid have elements corresponding to
        all grid points inside the cloudbox plus the cloud_box_limits*/
      
      p2gridpos(p_gp,
              p_grid[Range(cloudbox_limits[0], 
                           2,
                           (cloudbox_limits[1]- cloudbox_limits[0]))],
              p_grid[Range(cloudbox_limits[0], 
                           (cloudbox_limits[1]- cloudbox_limits[0])+1)]);
      gridpos(lat_gp,
              lat_grid[Range(cloudbox_limits[2], 
                             2,
                             (cloudbox_limits[3]- cloudbox_limits[2]))],
              lat_grid[Range(cloudbox_limits[2], 
                             (cloudbox_limits[3]- cloudbox_limits[2])+1)]);
      gridpos(lon_gp,
              lon_grid[Range(cloudbox_limits[4], 
                             2,
                             (cloudbox_limits[5]- cloudbox_limits[4]))],
              lon_grid[Range(cloudbox_limits[4], 
                             (cloudbox_limits[5]- cloudbox_limits[4])+1)]);

 
      //interpolation weights corresponding to pressure, latitude and 
      //longitude grids.

      Matrix itw_p((cloudbox_limits[1]- cloudbox_limits[0])+1, 2);
      Matrix itw_lat((cloudbox_limits[3]- cloudbox_limits[2])+1, 2);
      Matrix itw_lon((cloudbox_limits[5]- cloudbox_limits[4])+1, 2);

      interpweights ( itw_p, p_gp );
      interpweights ( itw_lat, lat_gp );
      interpweights ( itw_lon, lon_gp );

      // interpolation - pressure grid
      for (Index lat_index = 0; 
           lat_index <= (cloudbox_limits[3]-cloudbox_limits[2]); ++ lat_index)
        {
          for (Index lon_index = 0; 
               lon_index <= (cloudbox_limits[5]-cloudbox_limits[4]);
                            ++ lon_index)
            {
              for (Index za_index = 0; za_index < N_za ; ++ za_index)
                {
                  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
                    {
                      for (Index i = 0 ; i < N_i ; ++ i)
                        {
                          
                          VectorView target_field = i_field(Range(joker),
                                                            lat_index,
                                                            lon_index,
                                                            za_index,
                                                            aa_index,
                                                            i);
                          
                          ConstVectorView source_field = scat_i_p(f_index,
                                                                  Range(joker),    
                                                                  lat_index,
                                                                  lon_index,
                                                                  za_index,
                                                                  aa_index,
                                                                  i);
                          
                          interp(target_field,
                                 itw_p,
                                 source_field,
                                 p_gp);
                        }
                    }
                }
            } 
        }
      //interpolation latitude
      for (Index p_index = 0; 
           p_index <= (cloudbox_limits[1]-cloudbox_limits[0]) ; ++ p_index)
        {
          for (Index lon_index = 0; 
               lon_index <= (cloudbox_limits[5]-cloudbox_limits[4]) ;
               ++ lon_index)
            {
              for (Index za_index = 0; za_index < N_za ; ++ za_index)
                {
                  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
                    {
                      for (Index i = 0 ; i < N_i ; ++ i)
                        {
                          
                          VectorView target_field = i_field(p_index,
                                                            Range(joker),
                                                            lon_index,
                                                            za_index,
                                                            aa_index,
                                                            i);
                          
                          ConstVectorView source_field = scat_i_lat(f_index,
                                                                  p_index,
                                                                  Range(joker),    
                                                                  lon_index,
                                                                  za_index,
                                                                  aa_index,
                                                                  i);
                          
                          interp(target_field,
                                 itw_lat,
                                 source_field,
                                 lat_gp);
                        }
                    }
                }
            } 
        }
      //interpolation -longitude
      for (Index p_index = 0; 
           p_index <= (cloudbox_limits[1]-cloudbox_limits[0]); ++ p_index)
        {
          for (Index lat_index = 0; 
               lat_index <= (cloudbox_limits[3]-cloudbox_limits[2]);
               ++ lat_index)
            {
              for (Index za_index = 0; za_index < N_za ; ++ za_index)
                {
                  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
                    {
                      for (Index i = 0 ; i < N_i ; ++ i)
                        {
                          
                          VectorView target_field = i_field(p_index,
                                                            lat_index,
                                                            Range(joker),
                                                            za_index,
                                                            aa_index,
                                                            i);
                          
                          ConstVectorView source_field = scat_i_lon(f_index,
                                                                  p_index,    
                                                                  lat_index,
                                                                  Range(joker),
                                                                  za_index,
                                                                  aa_index,
                                                                  i);
                          
                          interp(target_field,
                                 itw_lon,
                                 source_field,
                                 lon_gp);
                        }
                    }
                }
            } 
        }
      //end of interpolation
    }//ends atmosphere_dim = 3
}


      
//! Set constant initial field

/*! 
  This method sets the initial field inside the cloudbox to a constant value.

  The user can specify a value for each Stokes dimension in the control file
  by the variable i_field_value, which is a vector containing 4 elements, the
  value of the initial field for each Stokes dimension.
  
  \param i_field Output : Intensity field
  \param scat_i_p Input : Intensity field on cloudbox boundary 
  (equal pressure surfaces)
  \param scat_i_lat Input : Intensity field on cloudbox boundary 
  (equal latitude surfaces)
  \param scat_i_lon Input : Intensity field on cloudbox boundary
  (equal longitude surfaces)
  \param f_grid Input : frequency grid
  \param f_index Input : the frequency index for scattering calculation
  \param p_grid Input : the pressure grid
  \param lat_grid Input : the latitude grid
  \param lon_grid Input : the longitude grid
  \param cloudbox_limits Input : Limits of the cloud box
  \param atmosphere_dim Input : dimension of atmosphere
  \param stokes_dim     FIXME: Add documentation.
  \param i_field_values Keyword : value of the constant initial field 

  \author Claudia Emde
  \date 2002-08-26
  
*/
void i_fieldSetConst(//WS Output:
                        Tensor6& i_field,
                        //WS Input:
                        const Tensor7& scat_i_p,
                        const Tensor7& scat_i_lat,
                        const Tensor7& scat_i_lon,
                        const Vector& f_grid,
                        const Index& f_index,
                        const Vector& p_grid,
                        const Vector& lat_grid,
                        const Vector& lon_grid,
                        const ArrayOfIndex& cloudbox_limits,
                        const Index& atmosphere_dim,
                        const Index& stokes_dim,
                        // Keyword       
                        const Vector& i_field_values)
{
  out2 << "Set initial field to constant values: " << i_field_values << "\n"; 

  // In the 1D case the atmospheric layers are defined by p_grid and the
  // required interface is scat_i_p.
  Index N_za = scat_i_p.npages();
  Index N_aa = scat_i_p.nrows();
  Index N_i = stokes_dim; 
  Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;

  // Check the input:
  assert( f_index < f_grid.nelem());
 
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  // Grids have to be adapted to atmosphere_dim.
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  
   // Check the input:
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*"); 


 
  if(atmosphere_dim == 1)
    {
      cout << "atm_dim = 1" << endl; 
      
    // Define the size of i_field.
    i_field.resize((cloudbox_limits[1] - cloudbox_limits[0])+1, 1, 1,  N_za,
                   1, N_i);
    i_field = 0.;

    // Loop over all zenith angle directions.
    for (Index za_index = 0; za_index < N_za; za_index++)
      {
        for (Index i = 0; i < stokes_dim; i++)
          { 
            //set the value for the upper boundary
            i_field(cloudbox_limits[1]-cloudbox_limits[0], 0, 0, za_index,
                    0, i) = 
	      scat_i_p(0, 1, 0, 0, za_index, 0, i);
            //set the value for the lower boundary 
            i_field(0, 0, 0, za_index, 0, i) =  
	      scat_i_p(0, 0, 0, 0, za_index, 0, i);
            for (Index scat_p_index = 1; scat_p_index < cloudbox_limits[1] - 
                   cloudbox_limits[0]; scat_p_index++ )
              // The field inside the cloudbox is set to some arbitrary value.
              i_field(scat_p_index, 0, 0, za_index, 0, i) =  i_field_values[i];
          }    
      }
    }
  else 
    {
      if ( !is_size(scat_i_p, f_grid.nelem(), 2, Nlat_cloud, 
                Nlon_cloud, N_za, N_aa, stokes_dim)  
           || !is_size(scat_i_lat, f_grid.nelem(), Np_cloud, 2, 
                       Nlon_cloud, N_za, N_aa, stokes_dim)  
           || !is_size(scat_i_lon, f_grid.nelem(), Np_cloud,  
                       Nlat_cloud, 2, N_za, N_aa, stokes_dim) )
        throw runtime_error(
                            "One of the interface variables (*scat_i_p*, "
                            "*scat_i_lat* or *scat_i_lon*) does not have "
                            "the right dimensions.  \n Probably you have "
                            "calculated them before for another value of "
                            "*stokes_dim*."
                            );
      

      
      cout << "atm_dim = 3" << endl;      
      i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1, 
                     (cloudbox_limits[3]- cloudbox_limits[2])+1,
                     (cloudbox_limits[5]- cloudbox_limits[4])+1,
                     N_za, 
                     N_aa,
                     N_i);
      
      i_field = 0.;
      

      // Loop over all directions:
      for (Index za_index = 0; za_index < N_za; za_index++)
        {
          for (Index aa_index = 0; aa_index < N_aa; aa_index++)
            {
              for (Index i = 0; i < stokes_dim; i++)
                { 
                  // pressure boundaries
                  // 
                   for (Index lat_index = cloudbox_limits[2]; 
                         lat_index <= cloudbox_limits[3]; lat_index++)
                      {
                        for (Index lon_index = cloudbox_limits[4]; 
                             lon_index <= cloudbox_limits[5]; lon_index++)
                          {
                            //set the value for the upper pressure boundary
                            i_field(cloudbox_limits[1]-cloudbox_limits[0], 
                                    lat_index-cloudbox_limits[2],
                                    lon_index-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_p(0, 1, lat_index-cloudbox_limits[2],
                                       lon_index-cloudbox_limits[4],
                                       za_index, aa_index, i);
                            //set the value for the lower pressure boundary 
                            i_field(0, lat_index-cloudbox_limits[2],
                                    lon_index-cloudbox_limits[4],
                                    za_index, aa_index, i) =  
                              scat_i_p(0, 0, lat_index-cloudbox_limits[2],
                                       lon_index-cloudbox_limits[4],
                                       za_index, aa_index, i);
                          }
                      }
                   
                   for (Index p_index = cloudbox_limits[0]; 
                        p_index <= cloudbox_limits[1]; p_index++)
                     {
                       // latitude boundaries
                       //
                        for (Index lon_index = cloudbox_limits[4]; 
                             lon_index <= cloudbox_limits[5]; lon_index++)
                          {
                            // first boundary
                            i_field(p_index-cloudbox_limits[0], 
                                    cloudbox_limits[3]-cloudbox_limits[2],
                                    lon_index-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_lat(0, p_index-cloudbox_limits[0],
                                         1, lon_index-cloudbox_limits[4],
                                         za_index, aa_index, i);
                            // second boundary
                            i_field(p_index-cloudbox_limits[0], 0, 
                                    lon_index-cloudbox_limits[4], 
                                    za_index, aa_index, i) =  
                              scat_i_lat(0, p_index-cloudbox_limits[0], 0,
                                         lon_index-cloudbox_limits[4],
                                         za_index, aa_index, i);
                            
                          }                                                    
                        // longitude boundaries
                   for (Index lat_index = cloudbox_limits[2]; 
                             lat_index <= cloudbox_limits[3]; lat_index++)
                          {
                            // first boundary
                            i_field(p_index-cloudbox_limits[0],
                                    lat_index-cloudbox_limits[2],
                                    cloudbox_limits[5]-cloudbox_limits[4],
                                    za_index, aa_index, i) = 
                              scat_i_lon(0, p_index-cloudbox_limits[0],
                                         lat_index-cloudbox_limits[2], 1,
                                         za_index, aa_index, i);
                            // second boundary
                            i_field(p_index-cloudbox_limits[0],  
                                    lat_index-cloudbox_limits[2],
                                    0, 
                                    za_index, aa_index, i) =  
                              scat_i_lon(0, p_index-cloudbox_limits[0],
                                         lat_index-cloudbox_limits[2], 0,
                                         za_index, aa_index, i);
                          } //lat_grid loop
                     } //p_grid loop
                   //
                   // Set the initial field to a constant value inside the 
                   // cloudbox:
                   // 
                   for( Index p_index = (cloudbox_limits[0]+1); 
                                 p_index <  cloudbox_limits[1] ;
                                 p_index ++)
                     {
                       for (Index lat_index = (cloudbox_limits[2]+1); 
                            lat_index < cloudbox_limits[3]; 
                            lat_index++)
                         {
                           for (Index lon_index = (cloudbox_limits[4]+1); 
                                lon_index < cloudbox_limits[5];
                                lon_index++)
                             {
                               i_field(p_index-cloudbox_limits[0],
                                       lat_index-cloudbox_limits[2],
                                       lon_index-cloudbox_limits[4],
                                       za_index, aa_index, i) =  
                                 i_field_values[i];
                             }
                         }
                     }
                } // stokes loop
            } // aa_grid loop
        } // za_grid loop
       
    } // atmosphere dim = 3
}


//! Initialize variables containing information about the particles.
/*! 
 
  WS Output:
  \param scat_data_raw  Single scattering data.
  \param pnd_field_raw  Particle number density field data.

*/
void ParticleTypeInit( //WS Output:
                      ArrayOfSingleScatteringData& scat_data_raw,
                      ArrayOfArrayOfTensor3& pnd_field_raw
                      )
{
  scat_data_raw.reserve(20);
  pnd_field_raw.reserve(20); 
}


//! Initialize variables containing information about the particles.
/*! 
 
  WS Output:
  \param amp_mat_raw    Amplitude matrix data.
  \param pnd_field_raw  Particle number density field data.

*/
void ParticleTypeInitAmpl( //WS Output:
                      ArrayOfArrayOfTensor6& amp_mat_raw,
                      ArrayOfArrayOfTensor3& pnd_field_raw
                      )
{
  amp_mat_raw.reserve(20);
  pnd_field_raw.reserve(20); 
}


//! Read single scattering data and particle number density field from 
//  data base.
/*! 
  This method allows the user to chose hydro-meteor species and particle
  number density fields. 
  There is one database for particle number density fields ( ....),
  which includes the following particle types:

  Another database (....) contains the single scattering properties for 
  hydro-meteor species.
  
  \param scat_data_raw Single scattering data.
  \param pnd_field_raw Particle number density field data.
  \param scat_data_file Filename for scattering data.
  \param pnd_field_file Filename for pnd field data.
*/
void ParticleTypeAdd( //WS Output:
                 ArrayOfSingleScatteringData& scat_data_raw,
                 ArrayOfArrayOfTensor3& pnd_field_raw,
                 // Keyword:
                 const String& scat_data_file,
                 const String& pnd_field_file)
{

  // Append *amp_mat_raw* and *pnd_field_raw* with empty Arrays of Tensors. 
  SingleScatteringData scat_data;
  scat_data_raw.push_back(scat_data);
  
  //ArrayOfTensor3 pnd_field_data;
  //pnd_field_raw.push_back(pnd_field_data);
    
  // Read amplitude matrix data:
  out2 << "Read single scattering data\n";
  xml_read_from_file( scat_data_file, scat_data_raw[scat_data_raw.nelem()-1]);
  //read_gridded_tensor3 (pnd_field_file, pnd_field_raw[pnd_field_raw.nelem()-1]);
       
}



//! Read amplitute matrix and particle number density field from data base.
/*! 
  This method allows the user to chose particle types and particle number 
  density fields. 
  There is one database for particle number density fields ( ....),
  which includes the following particle types:

  Another database (....) containes the amplitude matrices for those particle
  types from which all optical properties can be derived.
  
  \param amp_mat_raw Amplitude matrix data.
  \param pnd_field_raw Particle number density field data.
  \param amp_mat_file Filename for amplitude matrix data.
  \param pnd_field_file Filename for pnd field data.
*/
void ParticleTypeAddAmpl( //WS Output:
                 ArrayOfArrayOfTensor6& amp_mat_raw,
                 ArrayOfArrayOfTensor3& pnd_field_raw,
                 // Keyword:
                 const String& amp_mat_file,
                 const String& pnd_field_file)
{

  // Append *amp_mat_raw* and *pnd_field_raw* with empty Arrays of Tensors. 
  ArrayOfTensor6 amp_mat_data;
  amp_mat_raw.push_back(amp_mat_data);
  
  //ArrayOfTensor3 pnd_field_data;
  //pnd_field_raw.push_back(pnd_field_data);
    
  // Read amplitude matrix data:
  out2 << "Read amplitude matrix_data\n";
  read_gridded_tensor6( amp_mat_file, amp_mat_raw[amp_mat_raw.nelem()-1]);
  //read_gridded_tensor3 (pnd_field_file, pnd_field_raw[pnd_field_raw.nelem()-1]);
       
}


//! Method for the communication between cloudbox and clearsky.
/*
  This method puts the scattered radiation field into the interface variables 
  between the cloudbox and the clearsky, which are *scat_i_p*, *scat_i_lat*,
  *scat_i_lon*. As i_field is only stored for one frequency given by
  *f_index* this method has to be executed after each scattering 
  calculation to store the scattered field on the boundary of the cloudbox.
  
  The best way to calculate spectra including the influence of 
  scattering is to set up the *scat_mono_agenda* where this method 
  can be included.
 
  
 \param scat_i_p i_field on pressure boundaries
 \param scat_i_lat i_field on latitude boundaries
 \param scat_i_lon i_field on longitude boundaries
 \param i_field radiation field inside cloudbox
 \param f_grid frequency grid
 \param f_index index for scattering calculation
 \param p_grid pressure grid
 \param lat_grid latitude grid
 \param lon_grid longitude grid
 \param scat_za_grid zenith angle grid
 \param scat_aa_grid azimuth angle grid
 \param stokes_dim stokes dimension
 \param atmosphere_dim atmospheric dimension
 \param cloudbox_limits limits of the cloudbox

 \author Claudia Emde
 \date 2002-09-09
 
*/
void scat_iPut(//WS Output:
               Tensor7&  scat_i_p,
               Tensor7& scat_i_lat,
               Tensor7& scat_i_lon,
               //WS Input:
               const Tensor6& i_field,
               const Vector& f_grid,
               const Index& f_index,
               const Vector& p_grid,
               const Vector& lat_grid,
               const Vector& lon_grid,
               const Vector& scat_za_grid,
               const Vector& scat_aa_grid,
               const Index& stokes_dim,
               const Index& atmosphere_dim,
               const ArrayOfIndex& cloudbox_limits)
{
  // Some sizes:
  Index N_f = f_grid.nelem();
  Index N_za = scat_za_grid.nelem();
  Index N_aa = scat_aa_grid.nelem();
  Index N_p = cloudbox_limits[1] - cloudbox_limits[0] +1;

  // Some checks:
  
  assert( f_index < f_grid.nelem() );
 
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  // Grids have to be adapted to atmosphere_dim.
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  
   // Check the input:
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*"); 

 

 

  // End of checks.


  // Put the i_field at the cloudbox boundary into the interface variable 
  // scat_i_p.
  if(atmosphere_dim == 1)
    {
      // Check size of i_field.

       assert ( is_size( i_field, 
                    (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                    1, 
                    1,
                    N_za, 
                    1,
                    stokes_dim));

       assert ( is_size( scat_i_p,
                          N_f, 2, 1, 1, N_za, 1, stokes_dim ));

      
       for (Index za = 0; za < N_za; za++)
            {
              for (Index i = 0; i < stokes_dim; i++)
                {  
                  
                  //i_field at lower boundary
                  scat_i_p(f_index, 0, 0, 0,
                           za, 0, i) = 
                    i_field(0, 0, 0, za, 0, i);
                  //i_field at upper boundary
                  scat_i_p(f_index, 1, 0, 0,
                          za, 0, i) = 
                    i_field(cloudbox_limits[1] - cloudbox_limits[0],
                            0, 0, za, 0, i); 

                }//end stokes_dim
            }//end za loop
      

    }//end atmosphere_dim = 1
                                                        
      

      
  if(atmosphere_dim == 3)
    {
      // Some sizes relevant for 3D atmosphere
      Index N_lat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      Index N_lon = cloudbox_limits[5] - cloudbox_limits[4] + 1;


      // Check size of i_field.
      
      assert ( is_size( i_field, 
                        cloudbox_limits[1] - cloudbox_limits[0] + 1,
                        cloudbox_limits[3] - cloudbox_limits[2] + 1, 
                        cloudbox_limits[5] - cloudbox_limits[4] + 1, 
                        N_za, 
                        N_aa,
                        stokes_dim));

      // Resize interface variables:
      scat_i_p.resize(N_f, 2, N_lat, N_lon, N_za, N_aa, stokes_dim);
      scat_i_lat.resize(N_f, N_p, 2, N_lon, N_za, N_aa, stokes_dim);
      scat_i_lon.resize(N_f, N_p, N_lat, 2, N_za, N_aa, stokes_dim);


      //  assert ( is_size( scat_i_p,
//                         N_f, 2, N_lat, N_lon, N_za, N_aa, stokes_dim ));
      
//       assert ( is_size( scat_i_lat,
//                         N_f, N_p, 2, N_lon, N_za, N_aa, stokes_dim ));

//       assert ( is_size( scat_i_lon,
//                         N_f, N_p, N_lat, 2, N_za, N_aa, stokes_dim ));
      
 
      for (Index za = 0; za < N_za; za++)
        {
          for (Index aa = 0; aa < N_aa; aa++)
            {
              for (Index i = 0; i < stokes_dim; i++)
                {  
                  //
                  // Put i_field in scat_i_p:
                  //
                  for (Index lat = 0; lat < N_lat; lat++)
                    {
                      for (Index lon = 0; lon < N_lon; lon++)
                        {
                          //i_field at lower boundary
                          scat_i_p(f_index, 0, lat, lon,
                                   za, aa, i) = 
                            i_field(0, lat, lon, za, aa, i);
                          //i_field at upper boundary
                          scat_i_p(f_index, 1, lat, lon,
                                   za, aa, i) = 
                            i_field(cloudbox_limits[1]-cloudbox_limits[0],
                                    lat, lon, za, aa, i);
                        }
                    }
                  // 
                  // Put i_field in scat_i_lat:
                  //
                  for (Index p = 0; p < N_p; p++)
                    {
                      for (Index lon = 0; lon < N_lon; lon++)
                        {
                          //i_field at lower boundary
                          scat_i_lat(f_index, p, 0, lon,
                                     za, aa, i) = 
                            i_field(p, 0, lon, za, aa, i);
                          //i_field at upper boundary
                          scat_i_lat(f_index, p, 1, lon,
                                     za, aa, i) = 
                            i_field(p, cloudbox_limits[3]-
                                    cloudbox_limits[2],
                                    lon, za, aa, i);
                        }
                      //
                      // Put i_field in scat_i_lon:
                      for (Index lat = 0; lat < N_lat; lat++)
                        {
                          //i_field at lower boundary
                          scat_i_lon(f_index, p, lat, 0,
                                     za, aa, i) = 
                            i_field(p, lat, 0, za, aa, i);
                          //i_field at upper boundary
                          scat_i_lon(f_index, p, lat, 1,
                                     za, aa, i) = 
                            i_field(p, lat, cloudbox_limits[5]-
                                    cloudbox_limits[4], za, aa, i);
                        } 
                    }
                }
            }
        }
    }
}


//! Scattered radiance on the cloudbox boundary.
/* 
 This method returns the radiances for a given direction and position on the 
 boundary of the cloudbox. It interpolates from *scat_za_grid* on the 
 requested direction. The variable *i_out* is a matrix with the 
 dimensions [f_grid, stokes_dim].
  
  \param i_out Scattered radiance.
  \param i_out_name name for the outgoing radiance.
  \param scat_i_p i_field on pressure boundaries.
  \param scat_i_lat i_field on latitude boundaries.
  \param scat_i_lon i_field on longitude boundaries.
  \param rte_gp_p grid poition of the point on the boundary.
  \param rte_gp_lat grid position.
  \param rte_gp_lon grid position.
  \param rte_los direction.
  \param cloudbox_on is needed internally.
  \param atmosphere_dim
  \param stokes_dim Stokes dimension.
  \param scat_za_grid 
  \param scat_aa_grid
  \param f_grid

  \author Claudia Emde
  \date 2002-09-10

 */    
void CloudboxGetOutgoing(// WS Generic Output:
                         Matrix&   i_out,
                         // WS Generic Output Names:
                         const String&   i_out_name,
                         //WS Specific Input:
                         const Tensor7& scat_i_p,
                         const Tensor7& scat_i_lat,
                         const Tensor7& scat_i_lon,
                         const GridPos& rte_gp_p,
                         const GridPos& rte_gp_lat,
                         const GridPos& rte_gp_lon,
                         const Vector& rte_los,
                         const Index& cloudbox_on,
                         const ArrayOfIndex& cloudbox_limits,
                         const Index& atmosphere_dim,
                         const Index& stokes_dim,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Vector& f_grid)
{
  // Check the input:
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");

  if( !cloudbox_on )
    throw runtime_error( "The cloud box is not activated and no outgoing "
                                                    "field can be returned." );

  if( scat_za_grid.nelem() == 0 )
    throw runtime_error( "The variable *scat_za_grid* is empty. Are dummy "

                                            "values from *cloudboxOff used?" );

 if(atmosphere_dim == 1)
   {
     assert ( is_size( scat_i_p,
                       f_grid.nelem(), 2, 1, 1,
                       scat_za_grid.nelem(), 1,
                       stokes_dim ));
     cout << "rte_gp_p.idx" << rte_gp_p.idx << "\n";
     cout << "cloudbox_limits[0] " << cloudbox_limits[0] << "\n";
      cout << "cloudbox_limits[1]" << cloudbox_limits[1] << "\n";

     //Check, if grid_positions correspond to cloudbox boundary
     if (rte_gp_p.idx != cloudbox_limits[0] &&
	   rte_gp_p.idx != cloudbox_limits[1])
       throw runtime_error(
			   "Gridpositions have to be on the boundary of the "
			   "cloudbox defined by *cloudbox_limits*."
			   );

     //Define a vector to interpolate the outgoing radiance which is
     //defined on scat_za_grid on the requested zenith angle in
     //*cloudbox_los*.
       Vector zenith_angle(1);

       zenith_angle[0] = rte_los[0];

       //Array to store grid positions
       ArrayOfGridPos gp(1);
       gridpos(gp, scat_za_grid, zenith_angle);

       //Matrix to store interpolation weights
       //Matrix itw(scat_za_grid.nelem(),2);
       Matrix itw(gp.nelem(),2);
       interpweights(itw, gp);
       
       for(Index i = 0; i < stokes_dim; i++)
	 {
	   for(Index f_index = 0; f_index < f_grid.nelem(); f_index ++)
	     {
	       //This variable holds the radiation for a specified frequency.
	       //It is neccessairy because the interpolation is done for 
	       //each frequency separately.
	       Vector i_out_f(scat_za_grid.nelem());
	       
	       //lower boundary
	       if(rte_gp_p.idx == cloudbox_limits[0])
		 {
		   ConstVectorView i_f = scat_i_p(f_index, 0, 0, 0, 
						  Range(joker), 0, i);
		   i_out_f = i_f;
		 }
	       //upper boundary
	       else if(rte_gp_p.idx == cloudbox_limits[1])
		 {
		   ConstVectorView i_f = scat_i_p(f_index, 1, 0, 0,
						  Range(joker), 0, i);
		   i_out_f = i_f;
		 }
	       //Define vector for the interpolated radiance.
	       Vector i_out_los(1);
	       
	       //Do the interpolation:
	       interp(i_out_los, itw, i_out_f, gp);
	       
	       //Put the value into the matrix:
	       i_out(f_index, i) = i_out_los[0];
	     }//end frequency loop
	 }//end stokes_dim loop
            
   }// end atmosphere_dim 1
 

 if(atmosphere_dim == 3)
    {
      Index N_lat = scat_i_p.nshelves();
      Index N_lon = scat_i_p.nbooks();
      Index N_p = scat_i_lat.nvitrines();

      //Check consistency of input.

      assert ( is_size( scat_i_p,
                        f_grid.nelem(), 2, N_lat, N_lon, 
                        scat_za_grid.nelem(), scat_aa_grid.nelem(),
                        stokes_dim ));

      assert ( is_size( scat_i_lat,
                        f_grid.nelem(), N_p, 2, N_lon, 
                        scat_za_grid.nelem(), scat_aa_grid.nelem(),
                        stokes_dim ));
      assert ( is_size( scat_i_lon,
                        f_grid.nelem(), N_p, N_lat, 2, 
                        scat_za_grid.nelem(), scat_aa_grid.nelem(),
                        stokes_dim ));

      //Check, if grid_positions correspond to cloudbox boundary
     if ( (rte_gp_p.idx != cloudbox_limits[0] &&
           rte_gp_p.idx != cloudbox_limits[1]) &&
          (rte_gp_lat.idx != cloudbox_limits[2] &&
           rte_gp_lat.idx != cloudbox_limits[3]) &&
          (rte_gp_lon.idx != cloudbox_limits[4] &&
           rte_gp_lon.idx != cloudbox_limits[5] ) )
       throw runtime_error(
                           "Gridpositions have to be on the boundary of the "
                           "cloudbox defined by *cloudbox_limits*."
                           );

     //Define vectors to interpolate the outgoing radiance which is 
     //defined on scat_za_grid and scat_aa_grid on the requested zenith and
     // azimuth angle in 
     //*cloudbox_los*.
     Vector zenith_angle(1), azimuth_angle(1);
     zenith_angle[0] = rte_los[0];
     azimuth_angle[0] = rte_los[1];

     //Arrays to store grid positions
     ArrayOfGridPos gp_za(1), gp_aa(1);
     gridpos(gp_za, scat_za_grid, zenith_angle);
     gridpos(gp_aa, scat_aa_grid, azimuth_angle);

     //Matrices to store interpolation weights
     Tensor3 itw(1, 1, 4);
     interpweights(itw, gp_za, gp_aa);

     for(Index i = 0; i < stokes_dim; i++)
       {
         for(Index f_index = 0; f_index < f_grid.nelem(); f_index ++)
           {
             //This variable holds the radiation for a specified frequency.
             //It is neccessairy because the interpolation is done for 
             //each frequency separately.
             Matrix i_out_f(scat_za_grid.nelem(), scat_aa_grid.nelem());

             // The requested point lies on one of the pressure
             // boundaries.
             if ((rte_gp_p.idx == cloudbox_limits[0] ||
                  rte_gp_p.idx == cloudbox_limits[1] ) &&
                 ((rte_gp_lat.idx >= cloudbox_limits[2]) &&
                  rte_gp_lat.idx <= cloudbox_limits[3]) &&
                 ((rte_gp_lon.idx >= cloudbox_limits[4]) &&
                  rte_gp_lon.idx <= cloudbox_limits[5]) )
               {
                 
                 for(Index lat_index = 0; 
                     lat_index < cloudbox_limits[3]-cloudbox_limits[2]+1; lat_index++)
                   {
                     for(Index lon_index = 0;
                         lon_index < cloudbox_limits[5]-cloudbox_limits[4]+1;
                         lon_index ++)
                       {
                         if(rte_gp_p.idx == cloudbox_limits[0])
                           {
                             i_out_f = scat_i_p(f_index, 0, lat_index, 
                                                lon_index, joker, joker, i);
                           }
                         else if(rte_gp_p.idx == cloudbox_limits[1])
                           {
                             i_out_f = scat_i_p(f_index, 1, lat_index, 
                                                lon_index, joker, joker, i);
                           }
                       }
                   }
               }
             //
             // The requested point lies on one of the latitude
             // boundaries.
             else if ((rte_gp_lat.idx == cloudbox_limits[2] ||
                  rte_gp_lat.idx == cloudbox_limits[3] ) &&
                 ((rte_gp_p.idx >= cloudbox_limits[0]) &&
                  rte_gp_p.idx <= cloudbox_limits[1]) &&
                 ((rte_gp_lon.idx >= cloudbox_limits[4]) &&
                  rte_gp_lon.idx <= cloudbox_limits[5]) )
               
               {
                 
                 for(Index p_index = 0; 
                     p_index < cloudbox_limits[1]-cloudbox_limits[0]+1; p_index++)
                   {
                     for(Index lon_index = 0;
                         lon_index < cloudbox_limits[5]-cloudbox_limits[4]+1;
                         lon_index ++)
                       {
                         if(rte_gp_lat.idx == cloudbox_limits[2])
                           {
                             i_out_f = scat_i_lat(f_index, p_index, 0, 
                                                lon_index, joker, joker, i);
                           }
                         else if(rte_gp_lat.idx == cloudbox_limits[3])
                           {
                             i_out_f = scat_i_lat(f_index, p_index, 1, 
                                                  lon_index, joker, joker, i);
                           }
                       }
                   }
               }
             //
             // The requested point lies on one of the longitude
             // boundaries.
             else if ((rte_gp_lon.idx == cloudbox_limits[4] ||
                  rte_gp_lon.idx == cloudbox_limits[5] ) &&
                 ((rte_gp_p.idx >= cloudbox_limits[0]) &&
                  rte_gp_p.idx <= cloudbox_limits[1]) &&
                 ((rte_gp_lat.idx >= cloudbox_limits[2]) &&
                  rte_gp_lat.idx <= cloudbox_limits[3]) )
               
               {
                 
                 for(Index p_index = 0; 
                     p_index < cloudbox_limits[1]-cloudbox_limits[0]+1; p_index++)
                   {
                     for(Index lat_index = 0;
                         lat_index < cloudbox_limits[3]-cloudbox_limits[2]+1;
                         lat_index ++)
                       {
                         if(rte_gp_lon.idx == cloudbox_limits[4])
                           {
                             i_out_f = scat_i_lon(f_index, p_index, lat_index, 
                                                  0, joker, joker, i);
                           }
                         else if(rte_gp_lon.idx == cloudbox_limits[5])
                           {
                             i_out_f = scat_i_lon(f_index, p_index, lat_index, 
                                                  1, joker, joker, i);
                           }
                       }
                   }
               }
             
             else 
               {
                 throw runtime_error(
                                     "Error in CloudboxGetOutgoing:\n"
                                     "The point where you want to "
                                     "calculcate the radiation coming out \n"
                                     "of the cloud is not on the cloudbox \n"
                                     "boundary."
                                     );
               }
              
             // Outgoing radiances has to be interpolated on the
             // required direction.
             
             Matrix i_out_los(1,1);
             
             // Do the interpolation.
             interp(i_out_los, itw, i_out_f, gp_za, gp_aa);
             
             //Put the value into the matrix:
             i_out(f_index, i) = i_out_los(0,0);
             
           }//end frequency loop
       }//end stokes_dim loop
    }// end atmosphere_dim 3
   
}

//! This method gives the radiation field at the cloudbox boundary. 
/* 
   For each grid point on the cloudbox boundary in each propagation direction
   the radiation field is calculated using the method *RteCalc*, which 
   performs a clearsky radiative transfer calculation.
   
  \param scat_i_p i_field on pressure boundaries.
  \param scat_i_lat i_field on latitude boundaries.
  \param scat_i_lon i_field on longitude boundaries.
  \param i_in radiation coming into the cloudbox.
  \param cloudbox_pos Position on the cloudbox boundary.
  \param cloudbox_los Direction of radiation.
  \param cloudbox_limits Cloudbox limits.
  \param atmosphere_dim Atmospheric dimension.
  \param stokes_dim Stokes dimension.

  \author Sreerekha T.R., Claudia Emde
  \date 2002-10-07

 */    
void CloudboxGetIncoming(// WS Output:
                         Tensor7& scat_i_p,
                         Tensor7& scat_i_lat,
                         Tensor7& scat_i_lon,
                         Ppath& ppath,
                         Ppath& ppath_step,
                         Matrix& i_rte,
                         Matrix& i_space,
                         Matrix& ground_emission,
                         Matrix& ground_los,
                         Tensor4& ground_refl_coeffs,
                         Vector& rte_los,
                         Vector& rte_pos,
                         GridPos& rte_gp_p,
                         GridPos& rte_gp_lat,
                         GridPos& rte_gp_lon,
                         //WS Specific Input:
                         const ArrayOfIndex& cloudbox_limits,
                         const Index& atmosphere_dim,
                         const Index& stokes_dim,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Vector& f_grid,
                         const Agenda& ppath_step_agenda,
                         const Agenda& rte_agenda,
                         const Agenda& i_space_agenda,
                         const Agenda& ground_refl_agenda,
                         const Vector& p_grid,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Tensor3& z_field,
                         const Tensor3& t_field,
                         const Matrix& r_geoid,
                         const Matrix& z_ground
                         )

{

  out2 <<"Function: CloudboxGetIncoming \n"
       <<"Get clearsky field on cloudbox boundary. \n"
       <<"---------------------------------------- \n" ;

  Index Nf = f_grid.nelem();
  Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;

  Index Nza = scat_za_grid.nelem();

  Index Ni = stokes_dim;

  // Assign dummies for variables associated with sensor.
  bool     apply_sensor = false;
  Vector   mblock_za_grid_dummy(1);
           mblock_za_grid_dummy[0] = 0;
  Vector   mblock_aa_grid_dummy(0), sensor_rot_dummy(0);
  Matrix   sensor_pol_dummy;
  Index    antenna_dim_dummy = 1;
  Sparse   sensor_response_dummy;

  // Dummy variable for flag cloudbox_on. It has to be 0 for clearsky
  // calculations. We want to calculate the clearsky field on the boundary
  // so we have to set the variable to 0.
  Index cloudbox_on_dummy = 0;

  // Variable to avoid duplication of input checks in rte_calc
  bool   check_input = true;

  // Dummy for measurement vector
  Vector   y_dummy(0);

  if(atmosphere_dim == 1)
    {
      // In the 1D case this has to be resized again, as some grids are empty,
      // then Nlat, Nlon, Naa is 0.
      scat_i_p.resize(Nf, 2, 1, 1, Nza, 1, Ni);

      // Empty dummies for interface variables scat_i_p, scat_i_lat,
      // scat_i_lon.
      Tensor7 scat_i_p_dummy;
      Tensor7 scat_i_lat_dummy;
      Tensor7 scat_i_lon_dummy;

      //Define the variables for position and direction.
      Matrix sensor_los(1,1);
      Matrix sensor_pos(1,1);

      // Get scat_i_p at lower boundary
      sensor_pos(0,0) = r_geoid(0,0) + z_field(cloudbox_limits[0], 0, 0);

      for (Index scat_za_index = 0; scat_za_index < Nza; scat_za_index ++)
        {
          sensor_los(0,0) =  scat_za_grid[scat_za_index];

          rte_calc( y_dummy, ppath, ppath_step, i_rte,
               rte_pos, rte_los, rte_gp_p,
               rte_gp_lat, rte_gp_lon,
               i_space, ground_emission, ground_los,
               ground_refl_coeffs, ppath_step_agenda, rte_agenda, 
               i_space_agenda, ground_refl_agenda,
               atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, t_field,
               r_geoid,
               z_ground, cloudbox_on_dummy, cloudbox_limits, scat_i_p_dummy,
               scat_i_lat_dummy, scat_i_lon_dummy, scat_za_grid, scat_aa_grid,
               sensor_response_dummy, sensor_pos, sensor_los,
               sensor_pol_dummy, sensor_rot_dummy,
               f_grid, stokes_dim, antenna_dim_dummy,
               mblock_za_grid_dummy, mblock_aa_grid_dummy,
               check_input, apply_sensor, 1 );

          check_input = false;

          scat_i_p( Range(joker), 0, 0, 0,
                    scat_za_index,0,
                    Range(joker)) = i_rte;

                }

      // Get scat_i_p at upper boundary
       sensor_pos(0, 0) = r_geoid(0,0)+z_field(cloudbox_limits[1],0,0);

      for (Index scat_za_index = 0; scat_za_index < Nza;  scat_za_index ++)
        {
          sensor_los(0, 0) =  scat_za_grid[scat_za_index];

          rte_calc( y_dummy, ppath, ppath_step, i_rte,
                 rte_pos, rte_los, rte_gp_p,
                 rte_gp_lat, rte_gp_lon,
                 i_space, ground_emission, ground_los,
                 ground_refl_coeffs, ppath_step_agenda, rte_agenda,
                 i_space_agenda, ground_refl_agenda, 
                 atmosphere_dim, p_grid, lat_grid, lon_grid, z_field,
                 t_field, r_geoid,
                 z_ground, cloudbox_on_dummy, cloudbox_limits, scat_i_p_dummy,
                 scat_i_lat_dummy,
                 scat_i_lon_dummy, scat_za_grid, scat_aa_grid,
                 sensor_response_dummy, sensor_pos,
                 sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                 f_grid, stokes_dim, antenna_dim_dummy,
                 mblock_za_grid_dummy, mblock_aa_grid_dummy,
                 check_input, apply_sensor, 1 );

          
          scat_i_p( Range(joker), 1, 0,0,
                    scat_za_index,0,
                    Range(joker)) = i_rte;

           }

    }

  // atmosphere_dim = 3:
  else
    {
      Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      Index Naa = scat_aa_grid.nelem();
       // Convert scat_za_grid to "sensor coordinates"
      //(-180° < azimuth angle < 180°)
      //
      Vector aa_grid(Naa);
      for(Index i = 0; i<Naa; i++)
        aa_grid[i] = scat_aa_grid[i] - 180;

      // Resize interface variables:
      scat_i_p.resize(Nf, 2, Nlat_cloud, Nlon_cloud, Nza, Naa, Ni);
      scat_i_lat.resize(Nf, Np_cloud, 2, Nlon_cloud, Nza, Naa, Ni);
      scat_i_lon.resize(Nf, Np_cloud, Nlat_cloud, 2, Nza, Naa, Ni);

      // Empty dummies for interface variables scat_i_p, scat_i_lat,
      // scat_i_lon.
      Tensor7 scat_i_p_dummy;
      Tensor7 scat_i_lat_dummy;
      Tensor7 scat_i_lon_dummy;

      // Define the variables for position and direction.
      // LOS for 1 measurement block defined by zenith angle and azimuth angle.
      Matrix sensor_los(1,2);

      // Position defined by pressure, latitude, longitude.
      Matrix sensor_pos(1,3);


      // Get scat_i_p at lower boundary

      for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
        {
          for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
            {
              sensor_pos(0,0) = r_geoid(lat_index + cloudbox_limits[2],
                                        lon_index + cloudbox_limits[4]) 
                + z_field(cloudbox_limits[0],
                          lat_index + cloudbox_limits[2],
                          lon_index + cloudbox_limits[4]);
              sensor_pos(0,1) = lat_grid[lat_index + cloudbox_limits[2]];
              sensor_pos(0,2) = lon_grid[lon_index + cloudbox_limits[4]];

              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      sensor_los(0,0) = scat_za_grid[scat_za_index];
                      sensor_los(0,1) = aa_grid[scat_aa_index];
                      
                      rte_calc( y_dummy, ppath, ppath_step, i_rte,
                               rte_pos, rte_los, rte_gp_p, 
                               rte_gp_lat, rte_gp_lon,
                               i_space, ground_emission, ground_los, 
                               ground_refl_coeffs, ppath_step_agenda,
                               rte_agenda, 
                               i_space_agenda, ground_refl_agenda, 
                               atmosphere_dim, p_grid, lat_grid, lon_grid,
                               z_field,
                               t_field, r_geoid, z_ground, cloudbox_on_dummy,
                               cloudbox_limits, scat_i_p_dummy,
                               scat_i_lat_dummy, scat_i_lon_dummy,
                               scat_za_grid,
                               aa_grid, sensor_response_dummy, sensor_pos,
                               sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                               f_grid, stokes_dim, 
                               antenna_dim_dummy, 
                               mblock_za_grid_dummy, mblock_aa_grid_dummy, 
                               check_input, apply_sensor, 1 );

                      check_input = false;

                      scat_i_p( Range(joker), 0, lat_index, lon_index, 
                                scat_za_index, scat_aa_index,
                                Range(joker)) 
                        = i_rte;
                      
                    }
                }
            }
        }
      

      // Get scat_i_p at upper boundary

       for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
        {
          for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
            {
              sensor_pos(0,0) = r_geoid(lat_index + cloudbox_limits[2],
                                        lon_index + cloudbox_limits[4]) 
                + z_field(cloudbox_limits[1],
                          lat_index + cloudbox_limits[2],
                          lon_index + cloudbox_limits[4]);
              sensor_pos(0,1) = lat_grid[lat_index + cloudbox_limits[2]];
              sensor_pos(0,2) = lon_grid[lon_index + cloudbox_limits[4]];

              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      sensor_los(0,0) = scat_za_grid[scat_za_index];
                      sensor_los(0,1) = aa_grid[scat_aa_index];
                      
                      rte_calc( y_dummy, ppath, ppath_step, i_rte,
                               rte_pos, rte_los, rte_gp_p, 
                               rte_gp_lat, rte_gp_lon,
                               i_space, ground_emission, ground_los, 
                               ground_refl_coeffs, ppath_step_agenda,
                               rte_agenda, 
                               i_space_agenda, ground_refl_agenda, 
                               atmosphere_dim, p_grid, lat_grid, lon_grid,
                               z_field,
                               t_field, r_geoid, z_ground, cloudbox_on_dummy,
                               cloudbox_limits, scat_i_p_dummy,
                               scat_i_lat_dummy, scat_i_lon_dummy,
                               scat_za_grid,
                               aa_grid, sensor_response_dummy, sensor_pos,
                               sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                               f_grid, stokes_dim, antenna_dim_dummy, 
                               mblock_za_grid_dummy, mblock_aa_grid_dummy,
                               check_input, apply_sensor, 1 );
          
                      scat_i_p( Range(joker), 1, lat_index, lon_index, 
                                scat_za_index, scat_aa_index,
                                Range(joker)) 
                        = i_rte;
                      
                    }
                }
            }
        }
              
       
      // Get scat_i_lat (1st boundary):

      for (Index p_index = 0; p_index < Np_cloud; p_index++ )
        {
          for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
            {
              sensor_pos(0,0) = r_geoid(cloudbox_limits[2],
                                        lon_index + cloudbox_limits[4]) +
                z_field(p_index + cloudbox_limits[0],
                        cloudbox_limits[2],
                        lon_index + cloudbox_limits[4]);
              sensor_pos(0,1) = lat_grid[cloudbox_limits[2]];
              sensor_pos(0,2) = lon_grid[lon_index + cloudbox_limits[4]];

              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      sensor_los(0,0) = scat_za_grid[scat_za_index];
                      sensor_los(0,1) = aa_grid[scat_aa_index];
                      
                      rte_calc( y_dummy, ppath, ppath_step, i_rte,
                               rte_pos, rte_los, rte_gp_p, 
                               rte_gp_lat, rte_gp_lon,
                               i_space, ground_emission, ground_los, 
                               ground_refl_coeffs, ppath_step_agenda,
                               rte_agenda, 
                               i_space_agenda, ground_refl_agenda, 
                               atmosphere_dim, p_grid, lat_grid, lon_grid,
                               z_field,
                               t_field, r_geoid, z_ground, cloudbox_on_dummy,
                               cloudbox_limits, scat_i_p_dummy,
                               scat_i_lat_dummy, scat_i_lon_dummy,
                               scat_za_grid,
                               aa_grid, sensor_response_dummy, sensor_pos,
                               sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                               f_grid, stokes_dim, antenna_dim_dummy, 
                               mblock_za_grid_dummy, mblock_aa_grid_dummy,
                               check_input, apply_sensor, 1 );
          
                      scat_i_lat( Range(joker), p_index, 0, lon_index, 
                                scat_za_index, scat_aa_index,
                                Range(joker)) 
                        = i_rte;
                      
                    }
                }
            }
        }

      
      // Get scat_i_lat (2nd boundary)

      for (Index p_index = 0; p_index < Np_cloud; p_index++ )
        {
          for (Index lon_index = 0; lon_index < Nlon_cloud; lon_index++ )
            {
              sensor_pos(0,0) = r_geoid(cloudbox_limits[3],
                                        lon_index + cloudbox_limits[4]) +
                 z_field(p_index + cloudbox_limits[0],
                         cloudbox_limits[3],
                         lon_index + cloudbox_limits[4]);
               sensor_pos(0,1) = lat_grid[cloudbox_limits[3]];
               sensor_pos(0,2) = lon_grid[lon_index + cloudbox_limits[4]];
               
               for (Index scat_za_index = 0; scat_za_index < Nza;
                    scat_za_index ++)
                 {
                   for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      sensor_los(0,0) = scat_za_grid[scat_za_index];
                      sensor_los(0,1) = aa_grid[scat_aa_index];
                      
                      rte_calc( y_dummy, ppath, ppath_step, i_rte,
                               rte_pos, rte_los, rte_gp_p, 
                               rte_gp_lat, rte_gp_lon,
                               i_space, ground_emission, ground_los, 
                               ground_refl_coeffs, ppath_step_agenda,
                               rte_agenda, 
                               i_space_agenda, ground_refl_agenda, 
                               atmosphere_dim, p_grid, lat_grid, lon_grid,
                               z_field,
                               t_field, r_geoid, z_ground, cloudbox_on_dummy,
                               cloudbox_limits, scat_i_p_dummy,
                               scat_i_lat_dummy, scat_i_lon_dummy,
                               scat_za_grid,
                               aa_grid, sensor_response_dummy, sensor_pos,
                               sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                               f_grid, stokes_dim, antenna_dim_dummy, 
                               mblock_za_grid_dummy, mblock_aa_grid_dummy,
                               check_input, apply_sensor, 1 );
          
                      scat_i_lat( Range(joker), p_index, 1, lon_index, 
                                scat_za_index, scat_aa_index,
                                Range(joker)) 
                        = i_rte;
                      
                    }
                }
            }
        }    


       // Get scat_i_lon (1st boundary):

      for (Index p_index = 0; p_index < Np_cloud; p_index++ )
        {
          for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
            {
              sensor_pos(0,0) = r_geoid(lat_index + cloudbox_limits[2],
                                        cloudbox_limits[4]) + 
                                 z_field(p_index + cloudbox_limits[0],
                                         lat_index + cloudbox_limits[2],
                                         cloudbox_limits[4]);
              sensor_pos(0,1) = lat_grid[lat_index + cloudbox_limits[2]];
              sensor_pos(0,2) = lon_grid[cloudbox_limits[4]];

              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      sensor_los(0,0) = scat_za_grid[scat_za_index];
                      sensor_los(0,1) = aa_grid[scat_aa_index];
                      
                      rte_calc( y_dummy, ppath, ppath_step, i_rte,
                               rte_pos, rte_los, rte_gp_p, 
                               rte_gp_lat, rte_gp_lon,
                               i_space, ground_emission, ground_los, 
                               ground_refl_coeffs, ppath_step_agenda,
                               rte_agenda, 
                               i_space_agenda, ground_refl_agenda, 
                               atmosphere_dim, p_grid, lat_grid, lon_grid,
                               z_field,
                               t_field, r_geoid, z_ground, cloudbox_on_dummy,
                               cloudbox_limits, scat_i_p_dummy,
                               scat_i_lat_dummy, scat_i_lon_dummy,
                               scat_za_grid,
                               aa_grid, sensor_response_dummy, sensor_pos,
                               sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                               f_grid, stokes_dim, antenna_dim_dummy, 
                               mblock_za_grid_dummy, mblock_aa_grid_dummy,
                               check_input, apply_sensor, 1 );
          
                      scat_i_lon( Range(joker), p_index, lat_index, 0, 
                                scat_za_index, scat_aa_index,
                                Range(joker)) 
                        = i_rte;
                      
                    }
                }
            }
        }

      
      // Get scat_i_lon (2nd boundary)

      for (Index p_index = 0; p_index < Np_cloud; p_index++ )
        {
          for (Index lat_index = 0; lat_index < Nlat_cloud; lat_index++ )
            {
              sensor_pos(0,0) = r_geoid(lat_index + cloudbox_limits[2],
                                        cloudbox_limits[5]) + 
                                 z_field(p_index + cloudbox_limits[0],
                                         lat_index + cloudbox_limits[2],
                                         cloudbox_limits[5]);
              sensor_pos(0,1) = lat_grid[lat_index + cloudbox_limits[2]];
              sensor_pos(0,2) = lon_grid[cloudbox_limits[5]];
              
              for (Index scat_za_index = 0; scat_za_index < Nza;
                   scat_za_index ++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < Naa; 
                       scat_aa_index ++)
                    {
                      sensor_los(0,0) = scat_za_grid[scat_za_index];
                      sensor_los(0,1) = aa_grid[scat_aa_index];
                      
                      rte_calc( y_dummy, ppath, ppath_step, i_rte,
                               rte_pos, rte_los, rte_gp_p, 
                               rte_gp_lat, rte_gp_lon,
                               i_space, ground_emission, ground_los, 
                               ground_refl_coeffs, ppath_step_agenda,
                               rte_agenda, 
                               i_space_agenda, ground_refl_agenda, 
                               atmosphere_dim, p_grid, lat_grid, lon_grid,
                               z_field,
                               t_field, r_geoid, z_ground, cloudbox_on_dummy,
                               cloudbox_limits, scat_i_p_dummy,
                               scat_i_lat_dummy, scat_i_lon_dummy,
                               scat_za_grid,
                               aa_grid, sensor_response_dummy, sensor_pos,
                               sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                               f_grid, stokes_dim, antenna_dim_dummy, 
                               mblock_za_grid_dummy, mblock_aa_grid_dummy,
                               check_input, apply_sensor, 1 );
          
                      scat_i_lon( Range(joker), p_index, lat_index, 1, 
                                scat_za_index, scat_aa_index,
                                Range(joker)) 
                        = i_rte;
                      
                    }
                }
            }
        }

    }// End atmosphere_dim = 3.
  //
  out3 << "Finished calculation of incoming field on cloudbox boundary.\n";
}


//! This method gives the radiation field on the cloudbox boundary for a 
// spherically symmetric clearsky atmosphere. 

/* 
   For each grid point on the cloudbox boundary in each propagation direction
   the radiation field is calculated using the method *RteCalc*, which 
   performs a clearsky radiative transfer calculation.
   
  \param scat_i_p i_field on pressure boundaries.
  \param scat_i_lat i_field on latitude boundaries.
  \param scat_i_lon i_field on longitude boundaries.
  \param i_in radiation coming into the cloudbox.
  \param cloudbox_pos Position on the cloudbox boundary.
  \param cloudbox_los Direction of radiation.
  \param cloudbox_limits Cloudbox limits.
  \param atmosphere_dim Atmospheric dimension.
  \param stokes_dim Stokes dimension.

  \author Sreerekha T.R., Claudia Emde
  \date 2002-10-07

 */    
void CloudboxGetIncoming1DAtm(// WS Output:
                         Tensor7& scat_i_p,
                         Tensor7& scat_i_lat,
                         Tensor7& scat_i_lon,
                         Ppath& ppath,
                         Ppath& ppath_step,
                         Matrix& i_rte, 
                         Matrix& i_space,
                         Matrix& ground_emission,
                         Matrix& ground_los,
                         Tensor4& ground_refl_coeffs,
                         Vector& rte_los,
                         Vector& rte_pos,
                         GridPos& rte_gp_p,
                         GridPos& rte_gp_lat,
                         GridPos& rte_gp_lon,
                         //WS Specific Input:
                         const ArrayOfIndex& cloudbox_limits,
                         const Index& atmosphere_dim,
                         const Index& stokes_dim,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Vector& f_grid,
                         const Agenda& ppath_step_agenda,
                         const Agenda& rte_agenda,
                         const Agenda& i_space_agenda,
                         const Agenda& ground_refl_agenda,
                         const Vector& p_grid,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Tensor3& z_field,
                         const Tensor3& t_field,
                         const Matrix& r_geoid,
                         const Matrix& z_ground
                         )

{

  out2 <<"Function: CloudboxGetIncoming \n"
       <<"Get clearsky field on cloudbox boundary. \n"
       <<"---------------------------------------- \n" ;

  Index Nf = f_grid.nelem();
  Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;

  Index Nza = scat_za_grid.nelem();
  Index Naa = scat_aa_grid.nelem();
  Index Ni = stokes_dim;

  // Assign dummies for variables associated with sensor.
  bool     apply_sensor = false;
  Vector   mblock_za_grid_dummy(1);
           mblock_za_grid_dummy[0] = 0;
  Vector   mblock_aa_grid_dummy(0), sensor_rot_dummy(0);
  Matrix   sensor_pol_dummy;
  Index    antenna_dim_dummy = 1; 
  Sparse   sensor_response_dummy;

  // Dummy variable for flag cloudbox_on. It has to be 0 for clearsky
  // calculations. We want to calculate the clearsky field on the boundary
  // so we have to set the variable to 0.
  Index cloudbox_on_dummy = 0;

  // Variable to avoid duplication of input checks in rte_calc
  bool   check_input = true;

  // Dummy for measurement vector
  Vector   y_dummy(0);

  Vector aa_grid(Naa);
  for(Index i = 0; i<Naa; i++)
    aa_grid[i] = scat_aa_grid[i] - 180;
  
  Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1; 
  Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
  
  
  // As the atmosphere is spherically symmetric we only have to calculate 
  // one azimuth angle.
  Index scat_aa_index = 0;

  // Resize interface variables:
  scat_i_p.resize(Nf, 2, Nlat_cloud, Nlon_cloud, Nza, Naa, Ni);
  scat_i_lat.resize(Nf, Np_cloud, 2, Nlon_cloud, Nza, Naa, Ni);
  scat_i_lon.resize(Nf, Np_cloud, Nlat_cloud, 2, Nza, Naa, Ni);
  
  // Empty dummies for interface variables scat_i_p, scat_i_lat,
  // scat_i_lon.
  Tensor7 scat_i_p_dummy;
  Tensor7 scat_i_lat_dummy;
  Tensor7 scat_i_lon_dummy;

  // Define the variables for position and direction.
  // LOS for 1 measurement block defined by zenith angle and azimuth angle.
  Matrix sensor_los(1,2); 
  
  // Position defined by pressure, latitude, longitude.
  Matrix sensor_pos(1,3);

  Index lat_index = 0;
  Index lon_index = 0;

  // These variables are constant for all calculations:

  sensor_pos(0,1) = lat_grid[cloudbox_limits[2]];
  sensor_pos(0,2) = lon_grid[lon_index + cloudbox_limits[4]];

  sensor_los(0,1) = scat_aa_grid[scat_aa_index];

  // Get scat_i_p at lower boundary
  sensor_pos(0,0) = r_geoid(lat_index + cloudbox_limits[2],
                            lon_index + cloudbox_limits[4]) 
    + z_field(cloudbox_limits[0],
              lat_index + cloudbox_limits[2],
              lon_index + cloudbox_limits[4]);

  for (Index scat_za_index = 0; scat_za_index < Nza;
       scat_za_index ++)
    {
      sensor_los(0,0) = scat_za_grid[scat_za_index];
      
      rte_calc( y_dummy, ppath, ppath_step, i_rte,
                rte_pos, rte_los, rte_gp_p, 
                rte_gp_lat, rte_gp_lon,
                i_space, ground_emission, ground_los, 
                ground_refl_coeffs, ppath_step_agenda,
                rte_agenda, 
                i_space_agenda, ground_refl_agenda, 
                atmosphere_dim, p_grid, lat_grid, lon_grid,
                z_field,
                t_field, r_geoid, z_ground, cloudbox_on_dummy,
                cloudbox_limits, scat_i_p_dummy,
                scat_i_lat_dummy, scat_i_lon_dummy,
                scat_za_grid,
                aa_grid, sensor_response_dummy, sensor_pos,
                sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                f_grid, stokes_dim, 
                antenna_dim_dummy, 
                mblock_za_grid_dummy, mblock_aa_grid_dummy, 
                check_input, apply_sensor, 1 );
      
      check_input = false;

      for (Index lat = 0; lat < Nlat_cloud; lat ++)
        {
          for (Index lon = 0; lon < Nlon_cloud; lon ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_p( Range(joker), 0, lat, lon, 
                            scat_za_index, aa,
                            Range(joker)) 
                    = i_rte;
                }
            }
        }
    }
      

  // Get scat_i_p at upper boundary

  sensor_pos(0,0) = r_geoid(lat_index + cloudbox_limits[2],
                            lon_index + cloudbox_limits[4]) 
    + z_field(cloudbox_limits[1],
              lat_index + cloudbox_limits[2],
              lon_index + cloudbox_limits[4]);
  
  for (Index scat_za_index = 0; scat_za_index < Nza;
       scat_za_index ++)
    {
      sensor_los(0,0) = scat_za_grid[scat_za_index];
      
      rte_calc( y_dummy, ppath, ppath_step, i_rte,
                rte_pos, rte_los, rte_gp_p, 
                rte_gp_lat, rte_gp_lon,
                i_space, ground_emission, ground_los, 
                ground_refl_coeffs, ppath_step_agenda,
                rte_agenda, 
                i_space_agenda, ground_refl_agenda, 
                atmosphere_dim, p_grid, lat_grid, lon_grid,
                z_field,
                t_field, r_geoid, z_ground, cloudbox_on_dummy,
                cloudbox_limits, scat_i_p_dummy,
                scat_i_lat_dummy, scat_i_lon_dummy,
                scat_za_grid,
                aa_grid, sensor_response_dummy, sensor_pos,
                sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                f_grid, stokes_dim, antenna_dim_dummy, 
                mblock_za_grid_dummy, mblock_aa_grid_dummy,
                check_input, apply_sensor, 1 );
      

      for (Index lat = 0; lat < Nlat_cloud; lat ++)
        {
          for (Index lon = 0; lon < Nlon_cloud; lon ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_p( Range(joker), 1, lat, lon, 
                            scat_za_index, aa,
                            Range(joker)) 
                    = i_rte;
                }
            }
        }
    }

              
       
  // Get scat_i_lat (1st boundary):

  
  for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      sensor_pos(0,0) = r_geoid(cloudbox_limits[2],
                                lon_index + cloudbox_limits[4]) +
        z_field(p_index + cloudbox_limits[0],
                cloudbox_limits[2],
                lon_index + cloudbox_limits[4]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza;
           scat_za_index ++)
        {
          sensor_los(0,0) = scat_za_grid[scat_za_index];
          
          rte_calc( y_dummy, ppath, ppath_step, i_rte,
                    rte_pos, rte_los, rte_gp_p, 
                    rte_gp_lat, rte_gp_lon,
                    i_space, ground_emission, ground_los, 
                    ground_refl_coeffs, ppath_step_agenda,
                    rte_agenda, 
                    i_space_agenda, ground_refl_agenda, 
                    atmosphere_dim, p_grid, lat_grid, lon_grid,
                    z_field,
                    t_field, r_geoid, z_ground, cloudbox_on_dummy,
                    cloudbox_limits, scat_i_p_dummy,
                    scat_i_lat_dummy, scat_i_lon_dummy,
                    scat_za_grid,
                    aa_grid, sensor_response_dummy, sensor_pos,
                    sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                    f_grid, stokes_dim, antenna_dim_dummy, 
                    mblock_za_grid_dummy, mblock_aa_grid_dummy,
                    check_input, apply_sensor, 1 );

          
          for (Index lon = 0; lon < Nlon_cloud; lon ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_lat( Range(joker), p_index, 0, lon, 
                                scat_za_index, aa,
                                Range(joker)) 
                        = i_rte;
                }
            }
        }
    }
    

  
  // Get scat_i_lat (2nd boundary)
  
  for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      sensor_pos(0,0) = r_geoid(cloudbox_limits[3],
                                lon_index + cloudbox_limits[4]) +
        z_field(p_index + cloudbox_limits[0],
                cloudbox_limits[3],
                lon_index + cloudbox_limits[4]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza;
           scat_za_index ++)
        {
          
          sensor_los(0,0) = scat_za_grid[scat_za_index];
          
          rte_calc( y_dummy, ppath, ppath_step, i_rte,
                    rte_pos, rte_los, rte_gp_p, 
                    rte_gp_lat, rte_gp_lon,
                    i_space, ground_emission, ground_los, 
                    ground_refl_coeffs, ppath_step_agenda,
                    rte_agenda, 
                    i_space_agenda, ground_refl_agenda, 
                    atmosphere_dim, p_grid, lat_grid, lon_grid,
                    z_field,
                    t_field, r_geoid, z_ground, cloudbox_on_dummy,
                    cloudbox_limits, scat_i_p_dummy,
                    scat_i_lat_dummy, scat_i_lon_dummy,
                    scat_za_grid,
                    aa_grid, sensor_response_dummy, sensor_pos,
                    sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                    f_grid, stokes_dim, antenna_dim_dummy, 
                    mblock_za_grid_dummy, mblock_aa_grid_dummy,
                    check_input, apply_sensor, 1 );
          
          for (Index lon = 0; lon < Nlon_cloud; lon ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_lat( Range(joker), p_index, 1, lon, 
                            scat_za_index, aa,
                            Range(joker)) 
                    = i_rte;
                }
            }
        }
    }    


       // Get scat_i_lon (1st boundary):

  for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      sensor_pos(0,0) = r_geoid(lat_index + cloudbox_limits[2],
                                cloudbox_limits[4]) + 
        z_field(p_index + cloudbox_limits[0],
                lat_index + cloudbox_limits[2],
                cloudbox_limits[4]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza;
           scat_za_index ++)
        {
          sensor_los(0,0) = scat_za_grid[scat_za_index];
          
          rte_calc( y_dummy, ppath, ppath_step, i_rte,
                    rte_pos, rte_los, rte_gp_p, 
                    rte_gp_lat, rte_gp_lon,
                    i_space, ground_emission, ground_los, 
                    ground_refl_coeffs, ppath_step_agenda,
                    rte_agenda, 
                    i_space_agenda, ground_refl_agenda, 
                    atmosphere_dim, p_grid, lat_grid, lon_grid,
                    z_field,
                    t_field, r_geoid, z_ground, cloudbox_on_dummy,
                    cloudbox_limits, scat_i_p_dummy,
                    scat_i_lat_dummy, scat_i_lon_dummy,
                    scat_za_grid,
                    aa_grid, sensor_response_dummy, sensor_pos,
                    sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                    f_grid, stokes_dim, antenna_dim_dummy, 
                    mblock_za_grid_dummy, mblock_aa_grid_dummy,
                    check_input, apply_sensor, 1 );

          for (Index lat = 0; lat < Nlat_cloud; lat ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_lon( Range(joker), p_index, lat, 0, 
                            scat_za_index, aa,
                            Range(joker)) 
                    = i_rte;
                }
            }
        }
    }
  
  // Get scat_i_lon (2nd boundary)
  
  for (Index p_index = 0; p_index < Np_cloud; p_index++ )
    {
      sensor_pos(0,0) = r_geoid(lat_index + cloudbox_limits[2],
                                cloudbox_limits[5]) + 
        z_field(p_index + cloudbox_limits[0],
                lat_index + cloudbox_limits[2],
                cloudbox_limits[5]);
      
      for (Index scat_za_index = 0; scat_za_index < Nza;
           scat_za_index ++)
        {
          sensor_los(0,0) = scat_za_grid[scat_za_index];
          
          rte_calc( y_dummy, ppath, ppath_step, i_rte,
                    rte_pos, rte_los, rte_gp_p, 
                    rte_gp_lat, rte_gp_lon,
                    i_space, ground_emission, ground_los, 
                    ground_refl_coeffs, ppath_step_agenda,
                    rte_agenda, 
                    i_space_agenda, ground_refl_agenda, 
                    atmosphere_dim, p_grid, lat_grid, lon_grid,
                    z_field,
                    t_field, r_geoid, z_ground, cloudbox_on_dummy,
                    cloudbox_limits, scat_i_p_dummy,
                    scat_i_lat_dummy, scat_i_lon_dummy,
                    scat_za_grid,
                    aa_grid, sensor_response_dummy, sensor_pos,
                    sensor_los, sensor_pol_dummy, sensor_rot_dummy,
                    f_grid, stokes_dim, antenna_dim_dummy, 
                    mblock_za_grid_dummy, mblock_aa_grid_dummy,
                    check_input, apply_sensor, 1 );

          
          for (Index lat = 0; lat < Nlat_cloud; lat ++)
            {
              for (Index aa = 0; aa < Naa; aa ++)
                {
                  scat_i_lon( Range(joker), p_index, lat, 1, 
                            scat_za_index, aa,
                            Range(joker)) 
                    = i_rte;
                }
            }
        }
    }
  
  //
  out3 << "Finished calculation of incoming field on cloudbox boundary.\n";
}


//! The function does a batch calculation for metoffice fields.
/*!
  This method is used for simulating ARTS for metoffice model field
  The method reads the amsu data file stored in the variable 
  *met_amsu_data*.  This variable holds the latitude, longitude, 
  satellite zenith angle, and amsu-b corrected and uncorrected BTs.  
  The metoffice profiles are extracted for each of these lat and lon.
  The profiles are for pressure, temperature, altitude, humidity  and
  ice water content. In ARTS, the temperature data is stored in t_field_raw,
  which is an AraayOfTensor3, vmrs in vmr_field_raw which is an 
  ArrayOfArrayOfTensor3 and IWC is converted to corresponding particle 
  number density and stored in pnd_field_raw which is an ArrayOfArrayOfTensor3.
  Corresponding to each lat-lon pixel in the AMSU data, you have files
  storing t_field_raw, vmr_field_raw, pnd_field_raw.  The batch method loops
  over the number of lat-lon points, reads the data, sets the cloudbox and 
  perform the radiative transfer calculation.  The output *ybatch* holds the 
  spectra corresponding to all lat-lon points.
  
  \param ybatch Spectra for a batch of metoffice profiles
  \param y Agenda output: Spectra
  \param t_field_raw Agenda input: Temperature field
  \param z_field_raw Agenda input: Altitude field
  \param vmr_field_raw Agenda input: VMR field
  \param pnd_field_raw Agenda input: Raw particle number density field
  \param pnd_field Agenda input: Particle number density field
  \param p_grid Agenda input: Pressure grid
  \param sensor_los Agenda input: Sensor line of sight
  \param cloudbox_on Agenda input: Flag to activate cloudbox
  \param cloudbox_limits Agenda input: Limits of the cloudbox
  \param z_ground Agenda input: Ground height
  \param gas_species Tag group absorption
  \param met_profile_path Path of metoffice data
  \param met_profile_calc_agenda Agenda for absorption calculation and RT methods
  \param f_grid Frequency grid
  \param met_amsu_data Amsu data set
  \param sensor_pos Sensor position 
  \param r_geoid Geoid radius
  \param lat_grid Latitude grid
  \param lon_grid Longitude grid
  \param atmosphere_dim Atmospheric dimensionality
  \param scat_data_raw Single scattering data
  \param nelem_p_grid Keyword: Number of elements in pressure grid
  \param met_profile_path Keyword: Path of temperature, altitude, humidity fields
  \param met_profile_pnd_path Keyword: Path of pnd_field
  \author Sreerekha T.R.
  \date 2003-04-17
*/
void ybatchMetProfiles(//Output
		       Matrix& ybatch,
		       //Agenda communication variables
		       //Output of met_profile_calc_agenda
		       Vector& y,
		       //Input to met_profile_calc_agenda
		       ArrayOfTensor3& t_field_raw,
		       ArrayOfTensor3& z_field_raw,
		       ArrayOfArrayOfTensor3& vmr_field_raw,
		       ArrayOfArrayOfTensor3& pnd_field_raw,
		       Vector& p_grid,
		       Matrix& sensor_los,
		       Index& cloudbox_on,
		       ArrayOfIndex& cloudbox_limits,
		       Matrix& z_ground,
		       //Input
		       const ArrayOfArrayOfSpeciesTag& gas_species,
		       const Agenda& met_profile_calc_agenda,
		       const Vector& f_grid,
		       const Matrix& met_amsu_data,
		       const Matrix& sensor_pos,
		       const Matrix& r_geoid,
		       const Vector& lat_grid,
		       const Vector& lon_grid,
		       const Index& atmosphere_dim,
		       const ArrayOfSingleScatteringData& scat_data_raw,
		       //Keyword
		       const Index& nelem_p_grid,
		       const String& met_profile_path,
		       const String& met_profile_pnd_path)
{
  Index no_profiles = met_amsu_data.nrows();
  
  //*vmr_field_raw* is an ArrayOfArrayOfTensor3 where the first array
      //holds the gaseous species. 
      //Resize *vmr_field_raw* according to the number of gaseous species
      //elements
      vmr_field_raw.resize(gas_species.nelem());
  
  //The second array holds the gridded fields, one array for 
  //pressure, latitude, longitude and the data each.
  for (Index i = 0; i < gas_species.nelem(); ++ i)
    {
      vmr_field_raw[i].resize(4);
    }

  //pnd_field_raw is an ArrayOfArrayOfTensor3 where the first array
  //holds particle species.
  // Number of particle types:
  const Index N_pt = scat_data_raw.nelem();
  pnd_field_raw.resize(N_pt);
  
  // The satellite zenith angle is read in from the amsu data
  // and converted to arts sensor_los
  ConstVectorView sat_za_from_data = met_amsu_data(Range(joker),3);
  
  sensor_los.resize(1,1);
  
  // The lat and lon are extracted to get the proper file names of 
  // profiles
  ConstVectorView lat = met_amsu_data(Range(joker),0);
  ConstVectorView lon = met_amsu_data(Range(joker),1);
  
  z_ground.resize(1,1);
  
  // The spectra .
  y.resize(f_grid.nelem());
  
  // The batch spectra.
  ybatch.resize(no_profiles, f_grid.nelem());
  
  // Loop over the number of profiles.
  for (Index i = 0; i < no_profiles; ++ i)
    {
      ostringstream lat_os, lon_os;
      
      Index lat_prec = 3;
      if(lat[i] < 0) lat_prec--;
      if(abs(lat[i])>=10 )
	{
	  lat_prec--;
	  if(abs(lat[i])>=100 ) lat_prec--;
	}
      
      lat_os.setf (ios::showpoint | ios::fixed);
      lat_os << setprecision(lat_prec) << lat[i];
      
      Index lon_prec = 4;
      if(lon[i] < 0) lon_prec--;
      if(abs(lon[i])>=10 )
	{
	  lon_prec--;
	  if(abs(lon[i])>=100 ) lon_prec--;
	}
      lon_os.setf (ios::showpoint | ios::fixed);
      lon_os << setprecision(lon_prec) << lon[i];
      
      sensor_los(0,0) = 
      	180.0 - (asin(r_geoid(0,0) * sin(sat_za_from_data[i] * DEG2RAD) /sensor_pos(0,0)))* RAD2DEG;
      
      //Reads the t_field_raw from file
      xml_read_from_file(met_profile_path +"profile.lat_"+lat_os.str()+".lon_"+lon_os.str() + ".t.xml",
			 t_field_raw);
      
      //Reads the z_field_raw from file
      xml_read_from_file(met_profile_path +"profile.lat_"+lat_os.str()+".lon_"+lon_os.str()  + ".z.xml",
			 z_field_raw);
      
      //Reads the humidity from file - it is only an ArrayofTensor3
      xml_read_from_file(met_profile_path +"profile.lat_"+lat_os.str()+".lon_"+lon_os.str() + ".H2O.xml", 
       			 vmr_field_raw[0]);
      
      //Reads the pnd_field_raw for one particle
      //xml_read_from_file("/rinax/storage/users/rekha/uk_data/profiles/new_obs/newest_forecastfields/reff100/profile.lat_"+lat_os.str()+".lon_"+lon_os.str() + ".pnd100.xml",  pnd_field_raw[0]);
      
      xml_read_from_file(met_profile_pnd_path +"profile.lat_"+lat_os.str()+".lon_"+lon_os.str() + ".pnd100.xml",  pnd_field_raw[0]);
      
      //Write the profile number into a file.
      out2<<"profile_number.xml"<< i;
      
      // Set z_ground from lowest level of z_field 
      z_ground(0,0) = z_field_raw[3](0,0,0);
      
      /* The vmr_field_raw is an ArrayofArrayofTensor3 where the outer 
	 array is for species.
	 
	 The oxygen and nitrogen VMRs are set to constant values of 0.209
	 and 0.782, respectively and are used along with humidity field 
	 to generate *vmr_field_raw*.*/
            
      /*The second element of the species.  The first 3 Tensors in the
	array are the same .  They are pressure grid, latitude grid and
	longitude grid.  The third tensor which is the vmr is set to a 
	constant value of 0.782, corresponding to N2.*/
      vmr_field_raw[1][3].resize(vmr_field_raw[0][0].npages(),
				 vmr_field_raw[0][1].nrows(),
				 vmr_field_raw[0][2].ncols());
      vmr_field_raw[1][0] = vmr_field_raw[0][0]; //pressure grid for 1st element
      vmr_field_raw[1][1] = vmr_field_raw[0][1]; //latitude grid for 1st element
      vmr_field_raw[1][2] = vmr_field_raw[0][2]; //longitude grid for 1st element
      vmr_field_raw[1][3] = 0.782;//vmr of N2
      
      /*the third element of the species.  the first 3 Tensors in the
	array are the same .  They are pressure grid, latitude grid and
	longitude grid.  The third tensor which is the vmr is set to a 
	constant value of 0.209, corresponding to O2.*/
      vmr_field_raw[2][3].resize(vmr_field_raw[0][0].npages(),
				 vmr_field_raw[0][1].nrows(),
				 vmr_field_raw[0][2].ncols());
      vmr_field_raw[2][0] = vmr_field_raw[0][0];//pressure grid for 2nd element
      vmr_field_raw[2][1] = vmr_field_raw[0][1];//latitude grid for 2nd element
      vmr_field_raw[2][2] = vmr_field_raw[0][2];//longitude grid for 2nd element
      vmr_field_raw[2][3] =  0.209;//vmr of O2
      
      // N_p is the number of elements in the pressure grid
      Index N_p = t_field_raw[0].npages();
      
      //Making a p_grid with the first and last element taken from the profile.
      VectorNLogSpace(p_grid, 
		      "p_grid", 
		      t_field_raw[0](0,0,0), 
		      t_field_raw[0](N_p -1,0,0), 
		      nelem_p_grid);
      
      /*To set the cloudbox limits, the lower and upper cloudbox limits
	are to be set.  The lower cloudbox limit is set to the lowest
	pressure level.  The upper level is the highest level where the 
	ice water content is non-zero.*/
      Numeric cl_grid_min, cl_grid_max;
      
      //Lower limit = lowest pressure level of the original grid.
      //Could it be the interpolated p_grid? FIXME STR
      cl_grid_min = t_field_raw[0](0, 0, 0);
      
      // A counter for non-zero ice content
      Index level_counter = 0;
      
      // Loop over all pressure levels
      for (Index ip = 0; ip< N_p; ++ip)
 	{
	  //Checking for non-zero ice content. 0.001 is a threshold for
	  //ice water content.	  
	  if(pnd_field_raw[0][3](ip, 0, 0) > 0.001) 
	    {
	      ++level_counter;
	      //if non-zero ice content is found, it is set to upper 
	      //cloudbox limit. Moreover, we take one level higher 
	      // than the upper limit because we want the upper limit
	      //to have 0 pnd.
	      cl_grid_max = t_field_raw[0](ip +1, 0, 0);
	    }
	}
      
      //cloudbox limits have dimensions 2*atmosphere_dim
      cloudbox_limits.resize( atmosphere_dim*2 );
      
      //if there is no cloud in the considered profile, still we
      //need to set the upper limit. I here set the first level 
      //for the upper cloudbox limit.
      if(level_counter == 0)
	{
	  cl_grid_max = p_grid[1];
	}
      
      //Cloudbox is set.
      cloudboxSetManually(cloudbox_on, 
			  cloudbox_limits,
			  atmosphere_dim,
			  p_grid,
			  lat_grid,
			  lon_grid,
			  cl_grid_min,
			  cl_grid_max,
			  0,0,0,0);
      
      /*executing the met_profile_calc_agenda
	Agenda communication variables are
	Output of met_profile_calc_agenda : y
	Input to met_profile_calc_agenda  : t_field_raw,
	z_field_raw, vmr_field_raw, pnd_field_raw, p_grid,
	sensor_los, cloudbox_on, cloudbox_limits, z_ground, */
		       
      met_profile_calc_agenda.execute();
      
      //putting in the spectra *y* for each profile, thus assigning y
      //to the ith row of ybatch
      ybatch(i, Range(joker)) = y;
      
    }// closing the loop over profile basenames
}

//! The function does a batch calculation for metoffice fields.
/*!
  This method is used for simulating ARTS for metoffice model field
  This method loops over *met_profile_basenames* which contains the
  basenames of the metoffice profile files as an ArrayOfString.
  Corresponding to each basename we have temperature field, altitude
  field, humidity field and particle number density field.  The
  temperature field and altitude field are stored in the same dimensions
  as *t_field_raw* and *z_field_raw*.  The oxygen and nitrogen VMRs are
  set to constant values of 0.209 and 0.782, respectively and are used
  along with humidity field to generate *vmr_field_raw*.  
  
  The three fields *t_field_raw*, *z_field_raw*, and *vmr_field_raw* are
  given as input to *met_profile_calc_agenda* which is called in this
  method.  See documentation of WSM *met_profile_calc_agenda* for more
  information on this agenda
  
  \param ybatch spectra for a batch of metoffice profiles
  \param t_field_raw temperature field
  \param z_field_raw altitude field
  \param vmr_field_raw VMR field
  \param y spectra
  \param p_grid pressure grid
  \param sensor_los sensor line of sight
  \param z_ground ground height
  \param gas_species species under consideration
  \param met_profile_path Path to the MO profiles
  \param met_profile_calc_agenda agenda for absorption calculation and RT methods
  \param f_grid frequency grid
  \param met_amsu_data the latlon information in amsu obs
  \param sensor_pos sensor position
  \param r_geoid geoid radius
  \param nelem_p_grid number of levels in the pressure grid.
  
  \author Sreerekha T.R.
  \date 2003-04-17
*/
void ybatchMetProfilesClear(//Output
			    Matrix& ybatch,
			    ArrayOfTensor3& t_field_raw,
			    ArrayOfTensor3& z_field_raw,
			    ArrayOfArrayOfTensor3& vmr_field_raw,
			    Vector& y,
			    Vector& p_grid,
			    Matrix& sensor_los,
			    Matrix& z_ground,
			    //Input
			    const ArrayOfArrayOfSpeciesTag& gas_species,
			    const Agenda& met_profile_calc_agenda,
			    const Vector& f_grid,
			    const Matrix& met_amsu_data,
			    const Matrix& sensor_pos,
			    const Matrix& r_geoid,
			    //Keyword
			    const Index& nelem_p_grid,
			    const String& met_profile_path)
{
  Index no_profiles = met_amsu_data.nrows();
  //Index no_profiles = met_profile_basenames.nelem();
  // The humidity data is stored as  an ArrayOfTensor3 whereas
  // vmr_field_raw is an ArrayOfArrayOfTensor3
  ArrayOfTensor3 vmr_field_raw_h2o;
  
  vmr_field_raw.resize(gas_species.nelem());
  
  for (Index i = 0; i < gas_species.nelem(); ++ i)
    {
      vmr_field_raw[i].resize(4);
    }
  
  y.resize(f_grid.nelem());
  ybatch.resize(no_profiles, f_grid.nelem());
  
  Vector sat_za_from_profile;
  sat_za_from_profile = met_amsu_data(Range(joker),3);
  Numeric sat_za;
  
  sensor_los.resize(1,1);
    
  Vector lat, lon;
  lat = met_amsu_data(Range(joker),0);
  lon = met_amsu_data(Range(joker),1);

  Vector oro_height;
  oro_height = met_amsu_data(Range(joker),5);
  
  z_ground.resize(1,1);
  for (Index i = 0; i < no_profiles; ++ i)
    {
      ostringstream lat_os, lon_os;

      Index lat_prec = 3;
      if(lat[i] < 0) lat_prec--;
      if(abs(lat[i])>=10 )
	{
	  lat_prec--;
	  if(abs(lat[i])>=100 ) lat_prec--;
	}

      lat_os.setf (ios::showpoint | ios::fixed);
      lat_os << setprecision(lat_prec) << lat[i];
      
      Index lon_prec = 4;
      if(lon[i] < 0) lon_prec--;
      if(abs(lon[i])>=10 )
	{
	  lon_prec--;
	  if(abs(lon[i])>=100 ) lon_prec--;
	}
      lon_os.setf (ios::showpoint | ios::fixed);
      lon_os << setprecision(lon_prec) << lon[i];
      cout<<lat_os.str()<<endl;
      cout<<lon_os.str()<<endl;

      
      sat_za = sat_za_from_profile[i];
      
      //sensor_los(Range(joker),0) = 
      //	180.0 - (asin(r_geoid(0,0) * sin(sat_za * PI/180.) /sensor_pos(0,0)))*180./PI;
      sensor_los(Range(joker),0) = 
      	180.0 - (asin(r_geoid(0,0) * sin(sat_za * PI/180.) /sensor_pos(0,0)))*180./PI;
      cout<<"sensor_los"<<sat_za_from_profile[i]<<endl;
      cout<<"sensor_los"<<sat_za<<endl;
      cout<<"sensor_los"<<sensor_los<<endl;
      //Reads the t_field_raw from file
      
      xml_read_from_file(met_profile_path +"profile.lat_"+lat_os.str()+".lon_"+lon_os.str() + ".t.xml",
			 t_field_raw);
      //Reads the z_field_raw from file
      xml_read_from_file(met_profile_path +"profile.lat_"+lat_os.str()+".lon_"+lon_os.str()  + ".z.xml",
			 z_field_raw);
      
      //Reads the humidity from file - it is only an ArrayofTensor3
      // The vmr_field_raw is an ArrayofArrayofTensor3 where the outer 
      // array is for species
      xml_read_from_file(met_profile_path +"profile.lat_"+lat_os.str()+".lon_"+lon_os.str() + ".H2O.xml", 
      		 vmr_field_raw_h2o);
      //xml_read_from_file("/home/home01/rekha/uk/profiles/sat_vmr/profile.lat_"+lat_os.str()//+".lon_"+lon_os.str() + ".H2O_es.xml", vmr_field_raw_h2o);
      
      cout << "--------------------------------------------------------------------------"<<endl;
      cout << "The file" << met_profile_path +"profile.lat_"+lat_os.str()+".lon_"+lon_os.str()<< "is executed now"<<endl;
      cout << "--------------------------------------------------------------------------"<<endl; 
      xml_write_to_file("profile_number.xml",  i);
      // the first element of the species is water vapour. 
      
      // N_p is the number of elements in the pressure grid
      //z_ground(0,0) = oro_height[i]+ 0.01;
      z_ground(0,0) = z_field_raw[3](0,0,0);
      cout<<"z_ground"<<z_ground<<endl;
      Index N_p = t_field_raw[0].npages();
      
      vmr_field_raw[0] = vmr_field_raw_h2o;
      
      // the second element of the species.  the first 3 Tensors in the
      //array are the same .  They are pressure grid, latitude grid and
      // longitude grid.  The third tensor which is the vmr is set to a 
      // constant value of 0.782.
      vmr_field_raw[1][3].resize(vmr_field_raw[0][0].npages(),
				 vmr_field_raw[0][1].nrows(),
				 vmr_field_raw[0][2].ncols());
      vmr_field_raw[1][0] = vmr_field_raw[0][0];
      vmr_field_raw[1][1] = vmr_field_raw[0][1];
      vmr_field_raw[1][2] = vmr_field_raw[0][2];
      vmr_field_raw[1][3](joker, joker, joker) = 0.782;
      
      // the second element of the species.  the first 3 Tensors in the
      //array are the same .  They are pressure grid, latitude grid and
      // longitude grid.  The third tensor which is the vmr is set to a 
      // constant value of 0.209.
      vmr_field_raw[2][3].resize(vmr_field_raw[0][0].npages(),
				 vmr_field_raw[0][1].nrows(),
				 vmr_field_raw[0][2].ncols());
      vmr_field_raw[2][0] = vmr_field_raw[0][0];
      vmr_field_raw[2][1] = vmr_field_raw[0][1];
      vmr_field_raw[2][2] = vmr_field_raw[0][2];
      vmr_field_raw[2][3] (joker, joker, joker) = 0.209;
      
      //xml_write_to_file(met_profile_basenames[i]+ ".N2.xml", vmr_field_raw[1]);
      //xml_write_to_file(met_profile_basenames[i]+ ".O2.xml", vmr_field_raw[2]);
     
      //Making a p_grid with the first and last element taken from the profile.
      // this is because of the extrapolation problem.
      
      VectorNLogSpace(p_grid, 
		      "p_grid", 
		      t_field_raw[0](0,0,0), 
		      t_field_raw[0](N_p -1,0,0), 
		      nelem_p_grid);
      cout<<"t_field_raw[0](0,0,0)"<<t_field_raw[0](0,0,0)<<endl;
      cout<<"t_field_raw[0](N_p -1,0,0)"<<t_field_raw[0](N_p -1,0,0)<<endl;
      xml_write_to_file("p_grid.xml", p_grid);

      // executing the met_profile_calc_agenda
      met_profile_calc_agenda.execute();
      
      //putting in the spectra *y* for each profile
      ybatch(i, Range(joker)) = y;
      
    }// closing the loop over profile basenames
}
