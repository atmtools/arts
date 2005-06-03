/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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
  \file   m_rte.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-11 

  \brief  Workspace functions for solving clear sky radiative transfer.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "jacobian.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


//! RteCalc
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-09-17
*/
void RteCalc(
         Vector&                     y,
         Ppath&                      ppath,
         Ppath&                      ppath_step,
         Vector&                     ppath_p,
         Vector&                     ppath_t,
         Matrix&                     ppath_vmr,
         Matrix&                     iy,
         Vector&                     rte_pos,
         GridPos&                    rte_gp_p,
         GridPos&                    rte_gp_lat,
         GridPos&                    rte_gp_lon,
         Vector&                     rte_los,
         Sparse&                     jacobian,
         ArrayOfIndex&               rte_do_vmr_jacs,
         Tensor4&                    diy_dvmr,
         Index&                      rte_do_t_jacs,
         Tensor3&                    diy_dt,
   const Agenda&                     ppath_step_agenda,
   const Agenda&                     rte_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     iy_surface_agenda,
   const Agenda&                     iy_cloudbox_agenda,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   gas_species,
   const Matrix&                     r_geoid,
   const Matrix&                     z_surface,
   const Index&                      cloudbox_on, 
   const ArrayOfIndex&               cloudbox_limits,
   const Sparse&                     sensor_response,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Vector&                     f_grid,
   const Index&                      stokes_dim,
   const Index&                      antenna_dim,
   const Vector&                     mblock_za_grid,
   const Vector&                     mblock_aa_grid,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices )

{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index nmblock = sensor_pos.nrows();
  const Index nza     = mblock_za_grid.nelem();

  // Number of azimuthal direction for pencil beam calculations
  Index naa = mblock_aa_grid.nelem();
  if( antenna_dim == 1 )
    { naa = 1; }


  //--- Check input -----------------------------------------------------------
  //---------------------------------------------------------------------------

  // Agendas (agendas not always used are checked elsewhere when used)
  //
  chk_not_empty( "ppath_step_agenda", ppath_step_agenda );
  chk_not_empty( "rte_agenda", rte_agenda );

  // Stokes
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );

  // Basic checks of atmospheric, geoid and surface variables
  //  
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "z_field", z_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, 
                                                                    lon_grid );
  chk_atm_surface( "z_surface", z_surface, atmosphere_dim, lat_grid, 
                                                                    lon_grid );

  // Check that z_field has strictly increasing pages.
  //
  for( Index row=0; row<z_field.nrows(); row++ )
    {
      for( Index col=0; col<z_field.ncols(); col++ )
        {
          ostringstream os;
          os << "z_field (for latitude nr " << row << " and longitude nr " 
             << col << ")";
          chk_if_increasing( os.str(), z_field(joker,row,col) ); 
        }
    }

  // Check that there is no gap between the surface and lowest pressure 
  // surface
  //
  for( Index row=0; row<z_surface.nrows(); row++ )
    {
      for( Index col=0; col<z_surface.ncols(); col++ )
        {
          if( z_surface(row,col)<z_field(0,row,col) ||
                   z_surface(row,col)>=z_field(z_field.npages()-1,row,col) )
            {
              ostringstream os;
              os << "The surface altitude (*z_surface*) cannot be outside "
                 << "of the altitudes in *z_field*.";
              if( atmosphere_dim > 1 )
                os << "\nThis was found to be the case for:\n"
                   << "latitude " << lat_grid[row];
              if( atmosphere_dim > 2 )
                os << "\nlongitude " << lon_grid[col];
              throw runtime_error( os.str() );
            }
        }
    }

  // Cloud box
  //  
  chk_cloudbox( atmosphere_dim, p_grid, lat_grid, lon_grid,
                                                cloudbox_on, cloudbox_limits );

  // Frequency grid
  //
  if( nf == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );

  // Antenna
  //
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
  if( nza == 0 )
    throw runtime_error( "The measurement block zenith angle grid is empty." );
  chk_if_increasing( "mblock_za_grid", mblock_za_grid );
  if( antenna_dim == 1 )
    {
      if( mblock_aa_grid.nelem() != 0 )
        throw runtime_error( 
          "For antenna_dim = 1, the azimuthal angle grid must be empty." );
    }
  else
    {
      if( atmosphere_dim < 3 )
        throw runtime_error( "2D antennas (antenna_dim=2) can only be "
                                                 "used with 3D atmospheres." );
      if( mblock_aa_grid.nelem() == 0 )
        throw runtime_error(
                      "The measurement block azimuthal angle grid is empty." );
      chk_if_increasing( "mblock_aa_grid", mblock_aa_grid );
    }

  // Sensor
  //
  if( sensor_response.ncols() != nf * nza * naa * stokes_dim ) 
    {
      ostringstream os;
      os << "The *sensor_response* matrix does not have the right size, \n"
         << "either the method *sensor_responseInit* has not been run \n"
         << "prior to the call to *RteCalc* or some of the other sensor\n"
         << "response methods has not been correctly configured.";
      throw runtime_error( os.str() );
    }

  // Sensor position and LOS.
  //
  // That the angles are inside OK ranges are checked inside ppathCalc.
  //
  if( sensor_pos.ncols() != atmosphere_dim )
    throw runtime_error( "The number of columns of sensor_pos must be "
                              "equal to the atmospheric dimensionality." );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( 
                      "For 1D and 2D, sensor_los shall have one column." );
  if( atmosphere_dim == 3  &&  sensor_los.ncols() != 2 )
    throw runtime_error( "For 3D, sensor_los shall have two columns." );
  if( sensor_los.nrows() != nmblock )
    {
      ostringstream os;
      os << "The number of rows of sensor_pos and sensor_los must be "
         << "identical, but sensor_pos has " << nmblock << " rows,\n"
         << "while sensor_los has " << sensor_los.nrows() << " rows.";
      throw runtime_error( os.str() );
    }
  //---------------------------------------------------------------------------
  //--- End: Check input ------------------------------------------------------


  //--- Init y and ib ---------------------------------------------------------

  // Number of elements of *y* for one mblock
  Index    nblock = sensor_response.nrows();
  
  // Resize *y* to have correct length.
  y.resize( nmblock*nblock );

  // Create vector for MPB radiances for 1 measurement block.
  Vector ib( nf*nza*naa*stokes_dim );


  //--- Init Jacobian part ----------------------------------------------------

  rte_do_vmr_jacs.resize(0);
  rte_do_t_jacs   = 0;

  ArrayOfIndex    jqi_vmr(0);         // Index in jacobian_quantities of VMRs
  ArrayOfIndex    ji0_vmr(0);         // Start index in jacobian for anal. VMRs
  ArrayOfIndex    jin_vmr(0);         // Length of x for anal. VMRs
  Index           jqi_t = 0;          // As above, but for temperature
  Index           ji0_t = 0;        
  Index           jin_t = 0;
  ArrayOfMatrix   ib_vmr_jacs(0);
  Matrix          ib_t_jacs(0,0);

  for( Index i=0; i<jacobian_quantities.nelem(); i++ )
    {
      if ( jacobian_quantities[i].MainTag() == "Gas species"  &&  
                                          jacobian_quantities[i].Analytical() )
        { 
          jqi_vmr.push_back( i );
          // Find index in *gas_species* of jacobian species
          ArrayOfSpeciesTag   tags;
          array_species_tag_from_string( tags, 
                                             jacobian_quantities[i].Subtag() );
          Index si = chk_contains( "gas_species", gas_species, tags );
          rte_do_vmr_jacs.push_back( si ); 
          // Set size of MPB matrix
          ArrayOfIndex ji = jacobian_indices[i];
          const Index  nx = ji[1]-ji[0]+1;
          ji0_vmr.push_back( ji[0] );
          jin_vmr.push_back( nx );
          ib_vmr_jacs.push_back( Matrix(ib.nelem(),nx,0.0) );
        }
      if ( jacobian_quantities[i].MainTag() == "Temperature"  &&  
                                          jacobian_quantities[i].Analytical() )
        { 
          jqi_t         = i;
          rte_do_t_jacs = 1; 
          // Set size of MPB matrix
          ArrayOfIndex ji = jacobian_indices[i];
          const Index  nx = ji[1]-ji[0]+1;
          ji0_t = ji[0];
          jin_t = nx;
          ib_t_jacs = Matrix(ib.nelem(),nx,0.0);
        }
    }


  //--- Loop:  measurement block / zenith angle / azimuthal angle
  //
  Index    nydone = 0;                 // Number of positions in y done
  Index    nbdone;                     // Number of positions in ib done
  Vector   los( sensor_los.ncols() );  // LOS of interest
  //
  for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {
      nbdone = 0;

      for( Index iza=0; iza<nza; iza++ )
        {
          for( Index iaa=0; iaa<naa; iaa++ )
            {
              //--- Argument for verbosity of agendas
              bool  ag_verb = ( (iaa + iza + mblock_index) != 0 );

              //--- LOS of interest
              los     = sensor_los( mblock_index, joker );
              los[0] += mblock_za_grid[iza];
              if( antenna_dim == 2 )
                { los[1] += mblock_aa_grid[iaa]; }

              // Set the basic jacobian variables to be empty.
              // If they are coming back empty, then the jacobains can be
              // assumed to be zero.
              diy_dvmr.resize(0,0,0,0);
              diy_dt.resize(0,0,0);

              //--- Calculate *iy*
              iy_calc( iy, ppath, ppath_step, ppath_p, ppath_t, ppath_vmr,
                 rte_pos, rte_gp_p, rte_gp_lat, rte_gp_lon, rte_los, 
                 ppath_step_agenda, rte_agenda, iy_space_agenda, 
                 iy_surface_agenda, iy_cloudbox_agenda, atmosphere_dim, 
                 p_grid, lat_grid, lon_grid, z_field, t_field, vmr_field,
                 r_geoid, z_surface, 
                 cloudbox_on,  cloudbox_limits, sensor_pos(mblock_index,joker),
                 los, f_grid, stokes_dim, ag_verb );

              //--- Copy *iy* to *ib*
              for( Index is=0; is<stokes_dim; is++ )
                { ib[Range(nbdone+is,nf,stokes_dim)] = iy(joker,is); }


              //--- Jacobian part: --------------------------------------------

              //--- Gas species ---
              if( diy_dvmr.ncols() )  // See comment above
                {
                  for( Index ig=0; ig<rte_do_vmr_jacs.nelem(); ig++ )
                    {
                      //- Scale to other species retrieval modes
                      const String mode = 
                                       jacobian_quantities[jqi_vmr[ig]].Mode();
                      if( mode == "vmr" )
                        {}
                      else if( mode == "rel" )
                        {
                          for( Index ip=0; ip<ppath.np; ip++ )
                            { diy_dvmr(joker,joker,ip,ig) *= 
                                             ppath_vmr(rte_do_vmr_jacs[ig],ip);
                            }
                        }
                      else if( mode == "nd" )
                        {
                          for( Index ip=0; ip<ppath.np; ip++ )
                            { 
                              diy_dvmr(joker,joker,ip,ig) /= 
                                       number_density(ppath_p[ip],ppath_t[ip]);
                            }
                        }                  
                      else
                        { assert(0); }  // Should have been catched before

                      //- Map from ppath to retrieval quantities
                      Tensor3   diy_dx;
                      jacobian_from_path_to_rgrids( 
                                       diy_dx, diy_dvmr(joker,joker,joker,ig),
                                       atmosphere_dim, ppath, ppath_p, 
                                       jacobian_quantities[jqi_vmr[ig]] );

                      //- Copy obtained values to *ib_vmr_jacs*
                      for( Index is=0; is<stokes_dim; is++ )
                        { 
                         ib_vmr_jacs[ig](Range(nbdone+is,nf,stokes_dim),joker)=
                                                        diy_dx(joker,is,joker);
                        }
                    }
                }

              //--- Temperature ---
              if( diy_dt.ncols() && rte_do_t_jacs )
                {
                  //- Map from ppath to retrieval quantities
                  Tensor3   diy_dx;
                  jacobian_from_path_to_rgrids( 
                           diy_dx, diy_dt(joker,joker,joker), atmosphere_dim, 
                           ppath, ppath_p, jacobian_quantities[jqi_t] );

                  //- Copy obtained values to *ib_vmr_jacs*
                  for( Index is=0; is<stokes_dim; is++ )
                    { 
                      ib_t_jacs(Range(nbdone+is,nf,stokes_dim),joker)=
                                                        diy_dx(joker,is,joker);
                    }
                }

              //--- End of jacobian part --------------------------------------


              // Increase nbdone
              nbdone += nf*stokes_dim;

            } // iaa loop
        } // iza loop


      //--- Apply sensor response matrix on ib
      mult( y[Range(nydone,nblock)], sensor_response, ib );


      //--- Apply sensor response matrix on jacobians
      for( Index ig=0; ig<rte_do_vmr_jacs.nelem(); ig++ )
        {
          Matrix  K(nblock,jin_vmr[ig]);
          mult( K, sensor_response, ib_vmr_jacs[ig] );
          for( Index col=0; col<jin_vmr[ig]; col++ )
            {
              for( Index row=0; row<nblock; row++ )
                {
                  if( K(row,col) != 0 )
                    { jacobian.rw(nydone+row,ji0_vmr[ig]+col) = K(row,col); }
                }
            }
        }
      if( rte_do_t_jacs )
        {
          Matrix  K(nblock,jin_t);
          mult( K, sensor_response, ib_t_jacs );
          for( Index col=0; col<jin_t; col++ )
            {
              for( Index row=0; row<nblock; row++ )
                {
                  if( K(row,col) != 0 )
                    { jacobian.rw(nydone+row,ji0_t+col) = K(row,col); }
                }
            }
        }


      //--- Increase nydone
      nydone += nblock;
    }
}



//! RteCalcNoJacobian
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2005-05-30
*/
void RteCalcNoJacobian(
         Vector&                     y,
         Ppath&                      ppath,
         Ppath&                      ppath_step,
         Vector&                     ppath_p,
         Vector&                     ppath_t,
         Matrix&                     ppath_vmr,
         Matrix&                     iy,
         Vector&                     rte_pos,
         GridPos&                    rte_gp_p,
         GridPos&                    rte_gp_lat,
         GridPos&                    rte_gp_lon,
         Vector&                     rte_los,
   const Agenda&                     ppath_step_agenda,
   const Agenda&                     rte_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     iy_surface_agenda,
   const Agenda&                     iy_cloudbox_agenda,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const Matrix&                     r_geoid,
   const Matrix&                     z_surface,
   const Index&                      cloudbox_on, 
   const ArrayOfIndex&               cloudbox_limits,
   const Sparse&                     sensor_response,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Vector&                     f_grid,
   const Index&                      stokes_dim,
   const Index&                      antenna_dim,
   const Vector&                     mblock_za_grid,
   const Vector&                     mblock_aa_grid )
{
  Sparse                     jacobian;
  ArrayOfRetrievalQuantity   jacobian_quantities;
  ArrayOfArrayOfIndex        jacobian_indices;
  ArrayOfIndex               rte_do_vmr_jacs;
  Tensor4                    diy_dvmr(0,0,0,0);
  Index                      rte_do_t_jacs;
  Tensor3                    diy_dt(0,0,0);
  ArrayOfArrayOfSpeciesTag   gas_species(0);
  
  jacobianOff( jacobian, jacobian_quantities, jacobian_indices, 
               rte_do_vmr_jacs, rte_do_t_jacs );

  RteCalc( y, ppath, ppath_step, ppath_p, ppath_t, ppath_vmr, iy, rte_pos, 
           rte_gp_p, rte_gp_lat, rte_gp_lon, rte_los, jacobian, 
           rte_do_vmr_jacs, diy_dvmr, rte_do_t_jacs, diy_dt, 
           ppath_step_agenda, rte_agenda, iy_space_agenda, iy_surface_agenda,
           iy_cloudbox_agenda, atmosphere_dim, p_grid, lat_grid, lon_grid, 
           z_field, t_field, vmr_field, gas_species, r_geoid, z_surface, 
           cloudbox_on,  cloudbox_limits, sensor_response, sensor_pos, 
           sensor_los, f_grid, stokes_dim, antenna_dim, mblock_za_grid, 
           mblock_aa_grid, jacobian_quantities, jacobian_indices );
}



//! RteStd
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Claudia Emde and Patrick Eriksson
   \date   2003-01-07
*/
void RteStd(
      // WS Output:
             Matrix&         iy,
             Vector&         emission,
             Matrix&         abs_scalar_gas,
             Numeric&        rte_pressure,
             Numeric&        rte_temperature,
             Vector&         rte_vmr_list,
             Index&          f_index,
             Index&          ppath_index,
             Tensor4&        diy_dvmr,
             Tensor3&        diy_dt,
       // WS Input:
       const Ppath&          ppath,
       const Vector&         ppath_p,
       const Vector&         ppath_t,
       const Matrix&         ppath_vmr,
       const Vector&         f_grid,
       const Index&          stokes_dim,
       const Agenda&         emission_agenda,
       const Agenda&         scalar_gas_absorption_agenda,
       const ArrayOfIndex&   rte_do_gas_jacs,
       const Index&          rte_do_t_jacs )
{
  Tensor4 dummy(0,0,0,0);

  rte_std( iy, emission, abs_scalar_gas, rte_pressure, 
           rte_temperature, rte_vmr_list, f_index, ppath_index, 
           dummy, diy_dvmr, diy_dt,
           ppath, ppath_p, ppath_t, ppath_vmr, f_grid, stokes_dim, 
           emission_agenda, scalar_gas_absorption_agenda,
           rte_do_gas_jacs, rte_do_t_jacs, false );
}



//! RteStdWithTransmissions
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2005-05-19
*/
void RteStdWithTransmissions(
      // WS Output:
             Matrix&         iy,
             Vector&         emission,
             Matrix&         abs_scalar_gas,
             Numeric&        rte_pressure,
             Numeric&        rte_temperature,
             Vector&         rte_vmr_list,
             Index&          f_index,
             Index&          ppath_index,
             Tensor4&        ppath_transmissions,
             Tensor4&        diy_dvmr,
             Tensor3&        diy_dt,
       // WS Input:
       const Ppath&          ppath,
       const Vector&         ppath_p,
       const Vector&         ppath_t,
       const Matrix&         ppath_vmr,
       const Vector&         f_grid,
       const Index&          stokes_dim,
       const Agenda&         emission_agenda,
       const Agenda&         scalar_gas_absorption_agenda,
       const ArrayOfIndex&   rte_do_gas_jacs,
       const Index&          rte_do_t_jacs )
{
  rte_std( iy, emission, abs_scalar_gas, rte_pressure, 
           rte_temperature, rte_vmr_list, f_index, ppath_index, 
           ppath_transmissions, diy_dvmr, diy_dt,
           ppath, ppath_p, ppath_t, ppath_vmr, f_grid, stokes_dim, 
           emission_agenda, scalar_gas_absorption_agenda,
           rte_do_gas_jacs, rte_do_t_jacs, true );
}


