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
         Matrix&                     iy,
         Matrix&                     jacobian,
         Index&                      ppath_array_do,
         ArrayOfPpath&               ppath_array,
         Index&                      ppath_array_index,
   const Agenda&                     ppath_step_agenda,
   const Agenda&                     rte_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     surface_prop_agenda,
   const Agenda&                     iy_cloudbox_agenda,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
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

  ArrayOfIndex     rte_do_vmr_jacs (0);
  Index            rte_do_t_jacs = 0;
  ArrayOfTensor4   diy_dvmr;
  ArrayOfTensor4   diy_dt;

  ArrayOfIndex     jqi_vmr(0);        // Index in jacobian_quantities of VMRs
  ArrayOfIndex     ji0_vmr(0);        // Start index in jacobian for anal. VMRs
  ArrayOfIndex     jin_vmr(0);        // Length of x for anal. VMRs
  Index            jqi_t = 0;         // As above, but for temperature
  Index            ji0_t = 0;        
  Index            jin_t = 0;
  ArrayOfMatrix    ib_vmr_jacs(0);    // Correspondance to *ib* for VMR jac.
  Matrix           ib_t_jacs(0,0);    // Correspondance to *ib* for t jac.

  for( Index i=0; i<jacobian_quantities.nelem(); i++ )
    {
      if ( jacobian_quantities[i].MainTag() == "Gas species"  &&  
                                          jacobian_quantities[i].Analytical() )
        { 
          ppath_array_do = 1;
          jqi_vmr.push_back( i );
          // Find index in *abs_species* of jacobian species
          ArrayOfSpeciesTag   tags;
          array_species_tag_from_string( tags, 
                                             jacobian_quantities[i].Subtag() );
          Index si = chk_contains( "abs_species", abs_species, tags );
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
          ppath_array_do = 1;
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

              //--- Set *ppath_array* and *diy_dX*-variables to be empty
              ppath_array_index = -1;
              ppath_array.resize(0);
              diy_dvmr.resize(0);
              diy_dt.resize(0);

              //--- Calculate *iy*
              iy_calc( iy, ppath, ppath_step,
                       ppath_array_index,
                       ppath_array, diy_dvmr, diy_dt,
                       ppath_step_agenda, rte_agenda, iy_space_agenda, 
                       surface_prop_agenda, iy_cloudbox_agenda, atmosphere_dim,
                       p_grid, lat_grid, lon_grid, z_field, t_field, vmr_field,
                       r_geoid, z_surface, cloudbox_on, cloudbox_limits, 
                       sensor_pos(mblock_index,joker), los, f_grid, stokes_dim,
                       ppath_array_do,
                       rte_do_vmr_jacs, rte_do_t_jacs, ag_verb );

              //--- Copy *iy* to *ib*
              for( Index is=0; is<stokes_dim; is++ )
                { ib[Range(nbdone+is,nf,stokes_dim)] = iy(joker,is); }


              //--- Jacobian part: --------------------------------------------

              //--- Absorption species ---
              for( Index ig=0; ig<rte_do_vmr_jacs.nelem(); ig++ )
                {
                  //- Scale to other species retrieval modes
                  const String mode = jacobian_quantities[jqi_vmr[ig]].Mode();
                  if( mode == "vmr" )
                    {}
                  else if( mode == "rel" )
                    {
                      for( Index ia=0; ia<ppath_array.nelem(); ia++ )
                        {
                          if( ppath_array[ia].np > 1 )
                            {
                              for( Index ip=0; ip<ppath_array[ia].np; ip++ )
                                { diy_dvmr[ia](ig,ip,joker,joker) *= 
                                   ppath_array[ia].vmr(rte_do_vmr_jacs[ig],ip);
                                }
                            }
                        }
                    }
                  else if( mode == "nd" )
                    {
                      for( Index ia=0; ia<ppath_array.nelem(); ia++ )
                        {
                          if( ppath_array[ia].np > 1 )
                            {
                              for( Index ip=0; ip<ppath_array[ia].np; ip++ )
                                { 
                                  diy_dvmr[ia](ig,ip,joker,joker) /= 
                                       number_density( ppath_array[ia].p[ip],
                                                       ppath_array[ia].t[ip] );
                                }
                            }
                        }    
                    }              
                  else
                    { assert(0); }  // Should have been catched before

                  //- Map from ppath to retrieval quantities
                  {
                    jacobian_from_path_to_rgrids( ib_vmr_jacs[ig], nbdone,
                                     diy_dvmr, ig, atmosphere_dim, ppath_array,
                                     jacobian_quantities[jqi_vmr[ig]] );
                  }
                }

              //--- Temperature ---
              if( rte_do_t_jacs )
                {
                  //- Map from ppath to retrieval quantities
                  jacobian_from_path_to_rgrids( ib_t_jacs, nbdone, diy_dt, 0,
                                                atmosphere_dim, ppath_array, 
                                                jacobian_quantities[jqi_t] );
                }

              //--- End of jacobian part --------------------------------------


              // Increase nbdone
              nbdone += nf*stokes_dim;

            } // iaa loop
        } // iza loop


      //--- Apply sensor response matrix on ib
      mult( y[Range(nydone,nblock)], sensor_response, ib );


      //--- Apply sensor response matrix on jacobians, and store
      for( Index ig=0; ig<rte_do_vmr_jacs.nelem(); ig++ )
        {
          mult( jacobian(Range(nydone,nblock),Range(ji0_vmr[ig],jin_vmr[ig])),
                                            sensor_response, ib_vmr_jacs[ig] );
          /*
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
          */
        }
      if( rte_do_t_jacs )
        {
          mult( jacobian(Range(nydone,nblock),Range(ji0_t,jin_t)), 
                                                  sensor_response, ib_t_jacs );
          /*
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
          */
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
         Matrix&                     iy,
   const Agenda&                     ppath_step_agenda,
   const Agenda&                     rte_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     surface_prop_agenda,
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
  Matrix                     jacobian;
  ArrayOfRetrievalQuantity   jacobian_quantities;
  ArrayOfArrayOfIndex        jacobian_indices;
  ArrayOfArrayOfSpeciesTag   abs_species(0);
  Index                      ppath_array_do;
  ArrayOfPpath               ppath_array;
  Index                      ppath_array_index;


  jacobianOff( jacobian, jacobian_quantities, jacobian_indices );

  RteCalc( y, ppath, ppath_step, iy, jacobian, 
           ppath_array_do, ppath_array, ppath_array_index,
           ppath_step_agenda, rte_agenda, iy_space_agenda, surface_prop_agenda,
           iy_cloudbox_agenda, atmosphere_dim, p_grid, lat_grid, lon_grid, 
           z_field, t_field, vmr_field, abs_species, r_geoid, z_surface, 
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
             Matrix&                  iy,
             Vector&                  emission,
             Matrix&                  abs_scalar_gas,
             ArrayOfTensor4&          diy_dvmr,
             ArrayOfTensor4&          diy_dt,
       // WS Input:
       const Ppath&                   ppath,
       const ArrayOfPpath&            ppath_array, 
       const Index&                   ppath_array_index,
       const Vector&                  f_grid,
       const Index&                   stokes_dim,
       const Agenda&                  emission_agenda,
       const Agenda&                  abs_scalar_gas_agenda,
       const ArrayOfIndex&            rte_do_gas_jacs,
       const Index&                   rte_do_t_jacs )
{
  Tensor4 dummy(0,0,0,0);

  rte_std( iy, emission, abs_scalar_gas,
           dummy, diy_dvmr, diy_dt,
           ppath, ppath_array, ppath_array_index, f_grid, stokes_dim, 
           emission_agenda, abs_scalar_gas_agenda,
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
             Matrix&                  iy,
             Vector&                  emission,
             Matrix&                  abs_scalar_gas,
             Tensor4&                 ppath_transmissions,
             ArrayOfTensor4&          diy_dvmr,
             ArrayOfTensor4&          diy_dt,
       // WS Input:
       const Ppath&                   ppath,
       const ArrayOfPpath&            ppath_array, 
       const Index&                   ppath_array_index,
       const Vector&                  f_grid,
       const Index&                   stokes_dim,
       const Agenda&                  emission_agenda,
       const Agenda&                  abs_scalar_gas_agenda,
       const ArrayOfIndex&            rte_do_gas_jacs,
       const Index&                   rte_do_t_jacs )
{
  rte_std( iy, emission, abs_scalar_gas,
           ppath_transmissions, diy_dvmr, diy_dt,
           ppath, ppath_array, ppath_array_index, f_grid, stokes_dim, 
           emission_agenda, abs_scalar_gas_agenda,
           rte_do_gas_jacs, rte_do_t_jacs, true );
}


