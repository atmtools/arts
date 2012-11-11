#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "matpackVII.h"
#include "logic.h"
#include "ppath.h"
#include "agenda_class.h"
#include "physics_funcs.h"
#include "lin_alg.h"
#include "math_funcs.h"
#include "messages.h"
#include "xml_io.h"
#include "rte.h"
#include "special_interp.h"
#include "doit.h"
#include "m_general.h"
#include "wsv_aux.h"
#include "geodetic.h"



/* Workspace method: Doxygen documentation will be auto-generated */
void CloudboxGetIncoming2(
         Workspace&      ws,
         Tensor7&        doit_i_field,
   const Agenda&         iy_main_agenda,
   const Index&          atmosphere_dim,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Tensor3&        z_field,
   const Tensor3&        t_field,
   const Tensor4&        vmr_field,
   const Matrix&         z_surface,
   const Index&          cloudbox_on,
   const ArrayOfIndex&   cloudbox_limits,
   const Index&          basics_checked,
   const Index&          cloudbox_checked,
   const Vector&         f_grid,
   const Index&          stokes_dim,
   const String&         iy_unit,
   const Agenda&         blackbody_radiation_agenda,
   const Vector&         scat_za_grid,
   const Vector&         scat_aa_grid,
   const Index&          fill_interior,
   const Verbosity&)
{
  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) return;
  
  // Basics and cloudbox OK?
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  // Main sizes
  const Index  Nf   = f_grid.nelem();
  const Index  Np   = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index  Nza  = scat_za_grid.nelem();
  const Index  Ni   = stokes_dim;
        Index  Nlat=1, Nlon=1, Naa=1;

  //--- Special checks -------------------------------------------------------
  if( !(atmosphere_dim == 1  ||  atmosphere_dim == 3) )
    throw runtime_error( "The atmospheric dimensionality must be 1 or 3.");
  // DOIT requires frequency based radiance:
  if( iy_unit != "1"  || 
      !chk_if_std_blackbody_agenda( ws, blackbody_radiation_agenda ) )
    {
      ostringstream os;
      os << "It is assumed that you use this method together with DOIT.\n"
         << "Usage of this method then demands that the *iy_main_agenda*\n"
         << "returns frequency based radiance (ie. [W/m2/Hz/sr]).\n"
         << "This requires that *iy_unit* is set to \"1\" and that\n"
         << "*blackbody_radiation_agenda uses *blackbody_radiationPlanck*\n"
         << "or a corresponding WSM.\n"
         << "At least one of these requirements is not met.\n";
      throw runtime_error( os.str() );
    }
  if( scat_za_grid[0] != 0. || scat_za_grid[Nza-1] != 180. )
        throw runtime_error(
               "*scat_za_grid* must include 0 and 180 degrees as endpoints." );
  if( atmosphere_dim == 3 )
    {
      Naa = scat_aa_grid.nelem();
      if( scat_aa_grid[0] != 0. || scat_aa_grid[Naa-1] != 360. )
        throw runtime_error(
               "*scat_aa_grid* must include 0 and 360 degrees as endpoints." );
      Nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      Nlon = cloudbox_limits[5] - cloudbox_limits[4] + 1;
    }
  //--------------------------------------------------------------------------

  // Size doit_i_field and set to dummy value
  doit_i_field.resize( Nf, Np, Nlat, Nlon, Nza, Naa, Ni );
  doit_i_field = -999e9;

  // Convert scat_aa_grid to "sensor coordinates"
  // (-180° < azimuth angle < 180°)
  Vector aa_grid(Naa);
  for(Index i = 0; i<Naa; i++)
    { aa_grid[i] = scat_aa_grid[i] - 180; }

  // Define the variables for position and direction.
  Vector pos( atmosphere_dim );
  Vector los( max( Index(1), atmosphere_dim-1 ) );

  
  for( Index o=0; o<Nlon; o++ ) {
    for( Index a=0; a<Nlat; a++ ) {

      if( atmosphere_dim == 3 )
        {
          pos[1] = lat_grid[ a + cloudbox_limits[2] ];
          pos[2] = lon_grid[ o + cloudbox_limits[4] ];
        }

      for( Index p=Np-1; p>=0; p-- ) {

        if( fill_interior || p==0 || p==Np-1 || ( atmosphere_dim==3 && 
                              ( a==0 || a==Nlat-1 || o==0 || o==Nlon-1 ) ) ) {

          pos[0] = z_field( p+cloudbox_limits[0], 0, 0 );

          if( pos[0] <= z_surface(a,o) )
            throw runtime_error( "The surface is (not yet) allowed to be "
                                 "inside the cloudbox." );

          for( Index za=0; za<Nza; za++ ) {

            los[0] = scat_za_grid[za];
            Matrix iy;

            for( Index aa=0; aa<Naa; aa++ ) {
              // For end points of scat_za_grid, we need only to
              // perform calculations for first aa
              if( ( za != 0  &&  za != (Nza-1) )  ||  aa == 0 )
                {
                  if( atmosphere_dim == 3 )
                    { los[1] = aa_grid[aa]; }

                  get_iy( ws, iy, t_field, z_field, vmr_field, 0, f_grid, 
                          pos, los, Vector(0), iy_main_agenda );
                }
              
              doit_i_field( joker, p, a, o, za, aa, joker ) = iy;
            }
          } 
        } 
      } 
    } 
  }
}




/* Workspace method: Doxygen documentation will be auto-generated */
/*
void iyInterpCloudboxField2(
         Matrix&         iy,
   const Tensor7&        doit_i_field,
   const Vector&         rte_pos,
   const Vector&         rte_los,
   const Index&          jacobian_do,
   const Index&          cloudbox_on,
   const ArrayOfIndex&   cloudbox_limits,
   const Index&          basics_checked,
   const Index&          cloudbox_checked,
   const Index&          atmosphere_dim,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Tensor3&        z_field,
   const Index&          stokes_dim,
   const Vector&         scat_za_grid,
   const Vector&         scat_aa_grid,
   const Vector&         f_grid,
   const String&         interpmeth,
   const Verbosity& )
{
  //--- Check input -----------------------------------------------------------
  if( !(atmosphere_dim == 1  ||  atmosphere_dim == 3) )
    throw runtime_error( "The atmospheric dimensionality must be 1 or 3.");
  // Basics and cloudbox OK?
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );
  if( !cloudbox_on )
    throw runtime_error( "The cloud box is not activated and no outgoing "
                         "field can be returned." );
  if( scat_za_grid.nelem() == 0 )
    throw runtime_error( "The variable *scat_za_grid* is empty. Are dummy "
                         "values from *cloudboxOff used?" );
  if( !( interpmeth == "linear"  ||  interpmeth == "polynomial" ) )
    throw runtime_error( "Unknown interpolation method. Possible choices are "
                         "\"linear\" and \"polynomial\"." );
  if( interpmeth == "polynomial"  &&  atmosphere_dim != 1  )
    throw runtime_error( "Polynomial interpolation method is only available "
                         "for *atmosphere_dim* = 1." );
  if( doit_i_field.nlibraries() != f_grid.nelem() )
    throw runtime_error( "Inconsistency in size between *f_grid* and "
                         "*doit_i_field*." );
  if( jacobian_do )
    throw runtime_error( 
        "This method does not provide any jacobians (jacobian_do must be 0)" );
  //---------------------------------------------------------------------------


  // Convert rte_pos to grid positions
  GridPos gp_p, gp_lat, gp_lon;
  rte_pos2gridpos( gp_p, gp_lat, gp_lon, atmosphere_dim, 
                   p_grid, lat_grid, lon_grid, z_field, rte_pos );

  // Check pressure dimension
  Numeric fg = fractional_gp( gp_p );
  if( fg < Numeric(cloudbox_limits[0]) || fg > Numeric(cloudbox_limits[1]) )
    throw runtime_error( "The given *rte_pos* is outside the cloudbox "
                         "(in pressure dimension)!" );

  // Same for lat and lon
  if( atmosphere_dim == 3 )
    {
      fg = fractional_gp( gp_lat );
      if( fg < Numeric(cloudbox_limits[2]) || fg > Numeric(cloudbox_limits[3]) )
        throw runtime_error( "The given *rte_pos* is outside the cloudbox "
                             "(in latitude dimension)!" );
      fg = fractional_gp( gp_lon );
      if( fg < Numeric(cloudbox_limits[4]) || fg > Numeric(cloudbox_limits[5]) )
        throw runtime_error( "The given *rte_pos* is outside the cloudbox "
                             "(in longitude dimension)!" );
    }
}

*/
