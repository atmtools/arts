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
   const Index&          sensor_checked,
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
  if( !sensor_checked )
    throw runtime_error( "The sensor variables must be flagged to have passed"
                         "a consistency check (sensor_checked=1)." );

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

  // A matrix to flag if surface done for lat/lon
  Matrix  surface_done( Nlat, Nlon, -1 );

  // Convert scat_aa_grid to "sensor coordinates"
  // (-180° < azimuth angle < 180°)
  Vector aa_grid(Naa);
  for(Index i = 0; i<Naa; i++)
    { aa_grid[i] = scat_aa_grid[i] - 180; }

  // Define the variables for position and direction.
  Vector pos( atmosphere_dim );
  Vector los( max( Index(1), atmosphere_dim-1 ) );

  for( Index o=0; o<Nlon; o++ ) 
    {
      for( Index a=0; a<Nlat; a++ ) 
        {

          if( atmosphere_dim == 3 )
            {
              pos[1] = lat_grid[ a + cloudbox_limits[2] ];
              pos[2] = lon_grid[ o + cloudbox_limits[4] ];
            }

          for( Index p=Np-1; p>=0; p-- ) // Note looping downwards
            {

              if( fill_interior || p==0 || p==Np-1 || ( atmosphere_dim==3 && 
                                 ( a==0 || a==Nlat-1 || o==0 || o==Nlon-1 ) ) ) 
                {

                  pos[0] = z_field( p+cloudbox_limits[0], 0, 0 );

                  // Above the surface
                  if( pos[0] > z_surface(a,o) )
                    {
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

                              get_iy( ws, iy, t_field, z_field, vmr_field, 0, 
                                      f_grid, pos, los, Vector(0), 
                                      iy_main_agenda );
                            }
              
                          doit_i_field( joker, p, a, o, za, aa, joker ) = iy;
                        }
                      } 
                    }
                  // Surface or below
                  else
                    {
                      // If surface already done for lat/lon, just copy
                      // from above
                      if( surface_done( a, o ) > 0 )
                        {
                          for( Index za=0; za<Nza; za++ ) {
                            for( Index aa=0; aa<Naa; aa++ ) {
                              for( Index iv=0; iv<Nf; iv++ ) {
                                for( Index is=0; is<stokes_dim; is++ ) {
                                  doit_i_field(iv,p,a,o,za,aa,is) =
                                    doit_i_field(iv,p+1,a,o,za,aa,is);
                            } } } }
                        }
                      // Calculate for surface altitude
                      else
                        {
                          pos[0] = z_surface(a,o);
                          for( Index za=0; za<Nza; za++ ) {
                            //
                            los[0] = scat_za_grid[za];
                            Matrix iy;
                            //
                            for( Index aa=0; aa<Naa; aa++ ) {
                              // For end points of scat_za_grid, we need only to
                              // perform calculations for first aa
                              if( ( za != 0  &&  za != (Nza-1) )  ||  aa == 0 )
                                {
                                  if( atmosphere_dim == 3 )
                                    { los[1] = aa_grid[aa]; }

                                  get_iy( ws, iy, t_field, z_field, vmr_field, 
                                          0, f_grid, pos, los, Vector(0), 
                                          iy_main_agenda );
                                }
                              doit_i_field(joker,p,a,o,za,aa,joker) = iy;
                            }
                          }
                         surface_done( a, o ) = 1;
                        }
                    }
                }
            } 
        } 
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
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
   const Index&          za_order,
   const Verbosity& )
{
  // Basics and cloudbox OK?
  if( !(atmosphere_dim == 1  ||  atmosphere_dim == 3) )
    throw runtime_error( "The atmospheric dimensionality must be 1 or 3.");
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );
  if( !cloudbox_on )
    throw runtime_error( "The cloud box is not activated and no outgoing "
                         "field can be returned." );

  // Main sizes
  const Index nf  = f_grid.nelem();
  const Index nza = scat_za_grid.nelem();
  const Index np  = cloudbox_limits[1]-cloudbox_limits[0]+1;
        Index naa  = 1;
        Index nlat = 1;
        Index nlon = 1;
  if( atmosphere_dim == 3 )
    {
      naa  = scat_aa_grid.nelem();
      nlat = cloudbox_limits[3]-cloudbox_limits[2]+1;
      nlon = cloudbox_limits[5]-cloudbox_limits[4]+1;
    }

  //--- Check input -----------------------------------------------------------
  // rte_pos checked as part rte_pos2gridpos
  chk_rte_los( atmosphere_dim, rte_los );
  if( nza == 0 )
    throw runtime_error( "The variable *scat_za_grid* is empty. Are dummy "
                         "values from *cloudboxOff used?" );
  if( za_order < 1  ||  za_order > 5 )
    throw runtime_error( "The argument *za_order* must be in the range [1,5].");
  if( doit_i_field.ncols() != stokes_dim )
    throw runtime_error( "Inconsistency in size between *stokes_dim* and "
                         "*doit_i_field*." );
  if( doit_i_field.nrows() != naa )
    throw runtime_error( "Inconsistency in size between *scat_aa_grid* and "
                         "*doit_i_field*." );      
  if( doit_i_field.npages() != nza )
    throw runtime_error( "Inconsistency in size between *scat_za_grid* and "
                         "*doit_i_field*." );
  if( doit_i_field.nshelves() != nlat )
    throw runtime_error( "Inconsistency in size between *doit_i_field* and "
                         "latitude part of *cloudbox_limits*." );
  if( doit_i_field.nbooks() != nlon )
    throw runtime_error( "Inconsistency in size between *doit_i_field* and "
                         "longitude part of *cloudbox_limits*." );
  if( doit_i_field.nvitrines() != np )
    throw runtime_error( "Inconsistency in size between *doit_i_field* and "
                         "pressure part of *cloudbox_limits*." );
  if( doit_i_field.nlibraries() != nf )
    throw runtime_error( "Inconsistency in size between *f_grid* and "
                         "*doit_i_field*." );
  if( min( doit_i_field ) <= -999  )
    throw runtime_error( "A very large negative value found in *doit_i_field*. "
                         "Has the radiation field been calculated properly?" );
  if( jacobian_do )
    throw runtime_error( 
        "This method does not provide any jacobians (jacobian_do must be 0)" );
  //---------------------------------------------------------------------------

  // The approach is:
  //   1. if 3D shrink this to 1D.
  //   2. Interpolate in pressure (surface comes in here)
  //   3. Interpolate in zenith angle

  // Convert rte_pos to grid positions
  GridPos gp_p, gp_lat, gp_lon;
  rte_pos2gridpos( gp_p, gp_lat, gp_lon, atmosphere_dim, p_grid, lat_grid, 
                                                  lon_grid, z_field, rte_pos );

  // Remove lat, lon and azimuth angle dimensions
  //
  Tensor4 I4;
  //
  if( atmosphere_dim == 1 )
    { I4 = doit_i_field( joker, joker, 0, 0, joker, 0, joker ); }
  else
    {
      // Not ready!!!

      // Check and shift lat and lon grid positions
      Numeric fg = fractional_gp( gp_lat );
      if( fg < Numeric(cloudbox_limits[2]) || 
          fg > Numeric(cloudbox_limits[3]) )
        throw runtime_error( "The given *rte_pos* is outside the cloudbox "
                             "(in latitude dimension)!" );
      gp_lat.idx -= cloudbox_limits[2];
      gridpos_upperend_check( gp_lat, nlat );
      //
      fg = fractional_gp( gp_lon );
      if( fg < Numeric(cloudbox_limits[4]) || 
          fg > Numeric(cloudbox_limits[5]) )
        throw runtime_error( "The given *rte_pos* is outside the cloudbox "
                             "(in longitude dimension)!" );
      gp_lon.idx -= cloudbox_limits[4];
      gridpos_upperend_check( gp_lon, nlon );
    }      

  // Check pressure dimension
  Numeric fg = fractional_gp( gp_p );
  if( fg < Numeric(cloudbox_limits[0])-1e-6  ||   // 1-e6 to have some margin 
      fg > Numeric(cloudbox_limits[1])+1e-6 )     // for numerical issues
    throw runtime_error( "The given *rte_pos* is outside the cloudbox "
                         "(in pressure dimension)!" );

  // Pressure (no interpolation if at top or bottom)
  //
  Tensor3 I3;
  //
  if( fg <= Numeric(cloudbox_limits[0]) )
    { I3 = I4(joker,0,joker,joker); }
  else if( fg >= Numeric(cloudbox_limits[1]) )
    { I3 = I4(joker,np-1,joker,joker); }
  else
    { // Interpolation needed!:
      // Shift grid position
      gp_p.idx -= cloudbox_limits[0];
      gridpos_upperend_check( gp_p, np );

      // Interpolate in pressure
      Vector itw(2);
      interpweights( itw, gp_p );
      //
      I3.resize(nf,nza,stokes_dim);
      //
      for( Index iv=0; iv<nf; iv++ ) {
        for( Index iza=0; iza<nza; iza++ ) {
          for( Index is=0; is<stokes_dim; is++ ) {
            I3(iv,iza,is) = interp( itw, I4(iv,joker,iza,is), gp_p );
        } } }
    }

  iy.resize(nf,stokes_dim);

  // Interpolate in zenith angle
  GridPosPoly gp;
  gridpos_poly( gp, scat_za_grid, rte_los[0], za_order );
  Vector itw( za_order+1 );
  interpweights( itw, gp );
  for( Index iv=0; iv<nf; iv++ ) 
    {
      for( Index is=0; is<stokes_dim; is++ ) 
        {
          iy(iv,is) = interp( itw, I3(iv,joker,is), gp ); 
        } 
    }
}
