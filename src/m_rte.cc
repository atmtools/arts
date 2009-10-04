/* Copyright (C) 2002-2008
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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
#include "arts_omp.h"
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

extern const String ABSSPECIES_MAINTAG;
extern const String TEMPERATURE_MAINTAG;





/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void RteCalc(
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
         Matrix&                     jacobian,
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
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Vector&                     sensor_response_za,
   const Vector&                     sensor_response_aa,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Vector&                     f_grid,
   const Index&                      stokes_dim,
   const Index&                      antenna_dim,
   const Vector&                     mblock_za_grid,
   const Vector&                     mblock_aa_grid,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const String&                     y_unit,
   const String&                     jacobian_unit )
{
  // Consistency checks of input. 
  // Also returning some basic sizes (nblock = length(y) for one mblock)
  //
  Index nf=0, nmblock=0, nza=0, naa=0, nblock=0;  
  //
  rtecalc_check_input( nf, nmblock, nza, naa, nblock, stokes_dim, f_grid, 
                       atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                       t_field, r_geoid, z_surface, sensor_response, 
                       sensor_response_f, sensor_response_pol,
                       sensor_response_za, sensor_pos, sensor_los, antenna_dim, 
                       mblock_za_grid, mblock_aa_grid, cloudbox_on, 
                       cloudbox_limits, y_unit, jacobian_unit );

  // Agendas not checked elsewhere
  //
  chk_not_empty( "ppath_step_agenda", ppath_step_agenda );
  chk_not_empty( "rte_agenda", rte_agenda );


  //--- Init y and ib ---------------------------------------------------------

  // Resize *y* and auxiliary to have correct length.
  y.resize( nmblock*nblock );
  y_f.resize( nmblock*nblock );
  y_pol.resize( nmblock*nblock );
  y_pos.resize( nmblock*nblock, sensor_pos.ncols() );
  y_los.resize( nmblock*nblock, sensor_los.ncols() );

  // Create vector for MPB radiances for 1 measurement block.
  Vector ib( nf*nza*naa*stokes_dim );


  //--- Init Jacobian part ----------------------------------------------------
  //
  ArrayOfIndex     rte_do_vmr_jacs(0);
  Index            rte_do_t_jacs = 0;
  //
  ArrayOfIndex     jqi_vmr(0);      // Index in jacobian_quantities of VMRs
  ArrayOfIndex     ji0_vmr(0);      // Start index in jacobian for analyt. VMRs
  ArrayOfIndex     jin_vmr(0);      // Length of x for anal. VMRs
  Index            jqi_t = 0;       // As above, but for temperature
  Index            ji0_t = 0;       
  Index            jin_t = 0;
  ArrayOfMatrix    ib_vmr_jacs(0);  // Correspondance to *ib* for VMR jac.
  Matrix           ib_t_jacs(0,0);  // Correspondance to *ib* for t jac.
  //
  Index            ppath_array_do = 0;
  //
  String j_unit = jacobian_unit;
  if ( jacobian_unit == "-" )
    { j_unit = y_unit; }
  //
  for( Index i=0; i<jacobian_quantities.nelem(); i++ )
    {
      if ( jacobian_quantities[i].MainTag() == ABSSPECIES_MAINTAG  &&  
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
      if ( jacobian_quantities[i].MainTag() == TEMPERATURE_MAINTAG  &&  
                                          jacobian_quantities[i].Analytical() )
        { 
          ppath_array_do = 1;
          jqi_t          = i;
          rte_do_t_jacs  = 1; 
          // Set size of MPB matrix
          ArrayOfIndex ji = jacobian_indices[i];
          const Index  nx = ji[1]-ji[0]+1;
          ji0_t = ji[0];
          jin_t = nx;
          ib_t_jacs = Matrix(ib.nelem(),nx,0.0);
        }
    }

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_iy_space_agenda (iy_space_agenda);
  Agenda l_ppath_step_agenda (ppath_step_agenda);
  Agenda l_rte_agenda (rte_agenda);
  Agenda l_surface_prop_agenda (surface_prop_agenda);
  Agenda l_iy_cloudbox_agenda (iy_cloudbox_agenda);


  //--- Loop:  measurement block / zenith angle / azimuthal angle
  //
  Index    nydone = 0;                 // Number of positions in y done
  //
  for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {

#pragma omp parallel for                                        \
  if(!arts_omp_in_parallel() && nza>1)                          \
  default(none)                                                 \
  shared(nza, naa, nf, joker, mblock_index, rte_do_t_jacs,      \
         rte_do_vmr_jacs, ppath_array_do, ib, jqi_vmr,          \
         ib_vmr_jacs, j_unit, jqi_t, ib_t_jacs)                 \
  firstprivate(l_ws, l_iy_space_agenda,         \
               l_ppath_step_agenda, l_rte_agenda,               \
               l_surface_prop_agenda, l_iy_cloudbox_agenda)
      for( Index iza=0; iza<nza; iza++ )
        {
          // The try block here is necessary to correctly handle
          // exceptions inside the parallel region. 
          try
            {
              //--- Define *iy* and ppath variables
              Matrix           iy;
              Ppath            ppath;
              ArrayOfPpath     ppath_array;
              Index            ppath_array_index;
              ArrayOfTensor4   diy_dvmr;
              ArrayOfTensor4   diy_dt;

              for( Index iaa=0; iaa<naa; iaa++ )
                {
                  //--- Start index in *ib* for data to include 
                  const Index   nbdone = ( iza*naa + iaa ) * nf * stokes_dim;

                  //--- LOS of interest
                  //
                  Vector los( sensor_los.ncols() );
                  //
                  los     = sensor_los( mblock_index, joker );
                  los[0] += mblock_za_grid[iza];
                  //
                  if( antenna_dim == 2 )
                    {
                      throw runtime_error(
    "2D antennas are not yet correctly handled. Contact Patrick for details." );
                      los[1] += mblock_aa_grid[iaa]; 
                    }

                  //--- Set *ppath_array* and *diy_dX*-variables to be empty
                  //
                  ppath_array_index = -1;
                  ppath_array.resize(0);
                  //
                  diy_dvmr.resize(0);
                  diy_dt.resize(0);

                  //--- Calculate *iy*
                  iy_calc( l_ws, iy, ppath, ppath_array_index, ppath_array, 
                           diy_dvmr, diy_dt,
                           l_ppath_step_agenda, l_rte_agenda, 
                           l_iy_space_agenda, l_surface_prop_agenda, 
                           l_iy_cloudbox_agenda, atmosphere_dim, p_grid, 
                           lat_grid, lon_grid, z_field, t_field, vmr_field,
                           r_geoid, z_surface, cloudbox_on, cloudbox_limits, 
                           sensor_pos(mblock_index,joker), los, f_grid, 
                           stokes_dim, ppath_array_do, rte_do_vmr_jacs, 
                           rte_do_t_jacs );

                  //--- Unit conversions
                  apply_y_unit( iy, y_unit, f_grid );

                  //--- Copy *iy* to *ib*
                  for( Index is=0; is<stokes_dim; is++ )
                    { ib[Range(nbdone+is,nf,stokes_dim)] = iy(joker,is); }


                  //--- Jacobian part: -----------------------------------------

                  //--- Absorption species ---
                  for( Index ig=0; ig<rte_do_vmr_jacs.nelem(); ig++ )
                    {
                      //- Scale to other species retrieval modes
                      const String mode = 
                                        jacobian_quantities[jqi_vmr[ig]].Mode();
                      if( mode == "vmr" )
                        {}
                      else if( mode == "rel" )
                        {
                          for( Index ia=0; ia<ppath_array.nelem(); ia++ )
                            {
                              if( ppath_array[ia].np > 1 )
                                {
                                  for( Index ip=0; ip<ppath_array[ia].np; ip++ )
                                    diy_dvmr[ia](ig,ip,joker,joker) *= 
                                    ppath_array[ia].vmr(rte_do_vmr_jacs[ig],ip);
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
                                    diy_dvmr[ia](ig,ip,joker,joker) /= 
                                        number_density( ppath_array[ia].p[ip],
                                                        ppath_array[ia].t[ip] );
                                }
                            }    
                        }              
                      else
                        { assert(0); }  // Should have been catched before

                      //- Map from ppath to retrieval quantities
                      jacobian_from_path_to_rgrids( ib_vmr_jacs[ig], nbdone,
                                     diy_dvmr, ig, atmosphere_dim, ppath_array,
                                     jacobian_quantities[jqi_vmr[ig]] );

                      //--- Unit conversions
                      apply_y_unit( 
                            ib_vmr_jacs[ig](Range(nbdone,nf*stokes_dim),joker), 
                                                               j_unit, f_grid );
                    }

                  //--- Temperature ---
                  if( rte_do_t_jacs )
                    {
                      //- Map from ppath to retrieval quantities
                      jacobian_from_path_to_rgrids( ib_t_jacs, nbdone, diy_dt,
                                                 0, atmosphere_dim, ppath_array, 
                                                 jacobian_quantities[jqi_t] );

                      //--- Unit conversions
                      apply_y_unit(ib_t_jacs(Range(nbdone,nf*stokes_dim),joker), 
                                                               j_unit, f_grid );
                    }

                  //--- End of jacobian part -----------------------------------

                } // iaa loop

            } // end try block
          catch (runtime_error e)
            {
              exit_or_rethrow(e.what());
            }
        } // iza loop

 
      //--- Apply sensor response matrix on ib
      mult( y[Range(nydone,nblock)], sensor_response, ib );

      //--- Auxiliary variables
      for( Index ii=0; ii<nblock; ii++ )
        { 
          y_f[nydone+ii]         = sensor_response_f[ii];
          y_pol[nydone+ii]       = sensor_response_pol[ii]; 
          y_pos(nydone+ii,joker) = sensor_pos(mblock_index,joker);
          y_los(nydone+ii,0)     = sensor_los(mblock_index,0) +
                                   sensor_response_za[ii];
          if( sensor_response_aa.nelem() )
            { 
              y_los(nydone+ii,1) = sensor_los(mblock_index,0) +
                                   sensor_response_aa[ii]; 
            }
        }

      //--- Apply sensor response matrix on jacobians, and store
      for( Index ig=0; ig<rte_do_vmr_jacs.nelem(); ig++ )
        {
          mult( jacobian(Range(nydone,nblock),Range(ji0_vmr[ig],jin_vmr[ig])),
                                            sensor_response, ib_vmr_jacs[ig] );
        }
      if( rte_do_t_jacs )
        {
          mult( jacobian(Range(nydone,nblock),Range(ji0_t,jin_t)), 
                                                  sensor_response, ib_t_jacs );
        }


      //--- Increase nydone
      nydone += nblock;
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void RteCalcNoJacobian(
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
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
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Vector&                     sensor_response_za,
   const Vector&                     sensor_response_aa,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Vector&                     f_grid,
   const Index&                      stokes_dim,
   const Index&                      antenna_dim,
   const Vector&                     mblock_za_grid,
   const Vector&                     mblock_aa_grid,
   const String&                     y_unit )
{
  Matrix                     jacobian;
  ArrayOfRetrievalQuantity   jacobian_quantities;
  ArrayOfArrayOfIndex        jacobian_indices;
  String                     jacobian_unit;
  ArrayOfArrayOfSpeciesTag   abs_species(0);

  Index dummy;
  Agenda adummy;

  jacobianOff( dummy, adummy, jacobian, jacobian_quantities, jacobian_indices, 
                                                               jacobian_unit );

  RteCalc( ws, y, y_f, y_pol, y_pos, y_los, jacobian, 
           ppath_step_agenda, rte_agenda, iy_space_agenda, surface_prop_agenda,
           iy_cloudbox_agenda, atmosphere_dim, p_grid, lat_grid, lon_grid, 
           z_field, t_field, vmr_field, abs_species, r_geoid, z_surface, 
           cloudbox_on,  cloudbox_limits, sensor_response, sensor_response_f,
           sensor_response_pol, sensor_response_za, sensor_response_aa,
           sensor_pos, sensor_los, f_grid, stokes_dim, antenna_dim, 
           mblock_za_grid, mblock_aa_grid,
           jacobian_quantities, jacobian_indices,
           y_unit, jacobian_unit );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void RteCalcMC(
         Workspace&                     ws,
         Vector&                        y,
         Vector&                        y_f,
         ArrayOfIndex&                  y_pol,
         Matrix&                        y_pos,
         Matrix&                        y_los,
         Vector&                        mc_error,
   const Agenda&                        iy_space_agenda,
   const Agenda&                        surface_prop_agenda,
   const Agenda&                        opt_prop_gas_agenda,
   const Agenda&                        abs_scalar_gas_agenda, 
   const Index&                         atmosphere_dim,
   const Vector&                        p_grid,
   const Vector&                        lat_grid,
   const Vector&                        lon_grid,
   const Tensor3&                       z_field,
   const Tensor3&                       t_field,
   const Tensor4&                       vmr_field,
   const Matrix&                        r_geoid,
   const Matrix&                        z_surface,
   const Index&                         cloudbox_on, 
   const ArrayOfIndex&                  cloudbox_limits,
   const Tensor4&                       pnd_field,
   const ArrayOfSingleScatteringData&   scat_data_raw,
   const Sparse&                        sensor_response,
   const Vector&                        sensor_response_f,
   const ArrayOfIndex&                  sensor_response_pol,
   const Vector&                        sensor_response_za,
   const Vector&                        sensor_response_aa,
   const Matrix&                        sensor_pos,
   const Matrix&                        sensor_los,
   const Vector&                        f_grid,
   const Index&                         stokes_dim,
   const Index&                         antenna_dim,
   const Vector&                        mblock_za_grid,
   const Vector&                        mblock_aa_grid,
   const String&                        y_unit,
   const Numeric&                       mc_std_err,
   const Index&                         mc_max_time,
   const Index&                         mc_max_iter,
   const Index&                         mc_z_field_is_1D )
{

  // Consistency checks of input. Also returning some basic sizes
  //
  Index nf=0, nmblock=0, nza=0, naa=0, nblock=0;
  //
  if( atmosphere_dim != 3 )
        throw runtime_error( 
          "Monte Carlos calculations require that *atmosphere_dim* is 3." );
  //
  rtecalc_check_input( nf, nmblock, nza, naa, nblock, stokes_dim, f_grid, 
                       atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                       t_field, r_geoid, z_surface, sensor_response, 
                       sensor_response_f, sensor_response_pol,
                       sensor_response_za, sensor_pos, sensor_los, antenna_dim, 
                       mblock_za_grid, mblock_aa_grid, cloudbox_on, 
                       cloudbox_limits, y_unit, "-" );

  // Some MC variables are only local here
  Tensor3  mc_points;
  //
  MCAntenna mc_antenna;
  mc_antenna.set_pencil_beam();


  //--- Init y and ib ---------------------------------------------------------

  // Resize *y*, its aux variables and *mc_error* to have correct length.
  y.resize( nmblock*nblock );
  y_f.resize( nmblock*nblock );
  y_pol.resize( nmblock*nblock );
  y_pos.resize( nmblock*nblock, sensor_pos.ncols() );
  y_los.resize( nmblock*nblock, sensor_los.ncols() );
  mc_error.resize( nmblock*nblock );
  mc_error = 0.0;                     // Needed as values are accumulated below

  // Create vectors for MPB radiances for 1 measurement block.
  Vector ib( nf*nza*naa*stokes_dim );
  Vector ib_error( nf*nza*naa*stokes_dim );


  //--- Loop:  measurement block / zenith angle / azimuthal angle
  //
  Index    nydone = 0;                 // Number of positions in y done
  //

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_iy_space_agenda (iy_space_agenda);
  Agenda l_surface_prop_agenda (surface_prop_agenda);
  Agenda l_opt_prop_gas_agenda (opt_prop_gas_agenda);
  Agenda l_abs_scalar_gas_agenda (abs_scalar_gas_agenda);

  for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {
#pragma omp parallel for                                        \
  if(!arts_omp_in_parallel() && nza>1)                          \
  default(none)                                                 \
  shared(nza, naa, joker, mblock_index, nf, mc_antenna,         \
         mc_points, ib, ib_error)                               \
  firstprivate(l_ws, l_iy_space_agenda,         \
               l_surface_prop_agenda,l_opt_prop_gas_agenda,     \
               l_abs_scalar_gas_agenda)
      for( Index iza=0; iza<nza; iza++ )
        {
          // The try block here is necessary to correctly handle
          // exceptions inside the parallel region. 
          try
            {
              Matrix   los(1,sensor_los.ncols());  // LOS of interest
              Matrix   pos(1,sensor_pos.ncols());  // POS of interest

              // Vectors for a monochromatic single pencil beam calculation
              Vector iyf;
              Vector iyf_error;

              // Some MC variables
              Index mc_iteration_count;
              Index mc_seed;

              // Loop azimuth angles
              for( Index iaa=0; iaa<naa; iaa++ )
                {
                  //--- POS of interest
                  pos(0,joker)  = sensor_pos( mblock_index, joker );

                  //--- LOS of interest
                  los(0,joker)  = sensor_los( mblock_index, joker );
                  los(0,0)     += mblock_za_grid[iza];
                  if( antenna_dim == 2 )
                    { los(0,1) += mblock_aa_grid[iaa]; }

                  for( Index f_index=0; f_index<nf; f_index++ )
                    {
                      ArrayOfSingleScatteringData   scat_data_mono;

                      scat_data_monoCalc( scat_data_mono, scat_data_raw, 
                                          f_grid, f_index );

                      // Seed reset for each loop. If not done, the errors 
                      // appear to be highly correlated.
                      MCSetSeedFromTime( mc_seed );
                  
                      MCGeneral( l_ws,
                                 iyf, mc_iteration_count, iyf_error, mc_points, 
                                 mc_antenna, f_grid, f_index, pos, los,
                                 stokes_dim, l_iy_space_agenda,
                                 l_surface_prop_agenda, l_opt_prop_gas_agenda,
                                 l_abs_scalar_gas_agenda,
                                 p_grid, lat_grid, lon_grid, 
                                 z_field, r_geoid, z_surface,
                                 t_field, vmr_field, 
                                 cloudbox_limits, pnd_field,
                                 scat_data_mono, mc_seed, 
                                 y_unit, mc_std_err, mc_max_time, mc_max_iter, 
                                 mc_z_field_is_1D ); 
                  
                      //--- Start index in *ib* for data to include 
                      const Index   nbdone = ( ( iza*naa + iaa ) * nf + 
                                               f_index ) * stokes_dim;

                      //--- Copy *iyf* to *ib*
                      for( Index is=0; is<stokes_dim; is++ )
                        { 
                          ib[nbdone+is]       = iyf[is]; 
                          ib_error[nbdone+is] = iyf_error[is]; 
                        } // is loop

                    } // f_index loop
                } // iaa loop

            } // end try block
          catch (runtime_error e)
            {
              exit_or_rethrow(e.what());
            }
        } // iza loop


      //--- Apply sensor response matrix on ib and ib_error
      //
      mult( y[Range(nydone,nblock)], sensor_response, ib );
      //
      for( Index irow=0; irow<nblock; irow++ )
        {
          for( Index icol=0; icol<sensor_response.ncols(); icol++ )
            { mc_error[nydone+irow] += 
               pow( sensor_response(irow,icol)*ib_error[icol], (Numeric)2.0 ); 
            }
        }

      //--- Auxiliary variables
      for( Index ii=0; ii<nblock; ii++ )
        { 
          y_f[nydone+ii]         = sensor_response_f[ii];
          y_pol[nydone+ii]       = sensor_response_pol[ii]; 
          y_pos(nydone+ii,joker) = sensor_pos(mblock_index,joker);
          y_los(nydone+ii,0)     = sensor_los(mblock_index,0) +
                                   sensor_response_za[ii];
          if( sensor_response_aa.nelem() )
            { 
              y_los(nydone+ii,1) = sensor_los(mblock_index,0) +
                                   sensor_response_aa[ii]; 
            }
        }

      //--- Increase nydone
      nydone += nblock;
    }

  // Convert *mc_error* from variances to std. devs.
  transform( mc_error, sqrt, mc_error );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void RteStd( Workspace&               ws,
      // WS Output:
             Matrix&                  iy,
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

  rte_std( ws, iy, dummy, diy_dvmr, diy_dt,
           ppath, ppath_array, ppath_array_index, f_grid, stokes_dim, 
           emission_agenda, abs_scalar_gas_agenda,
           rte_do_gas_jacs, rte_do_t_jacs, false );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void RteStdWithTransmissions(
             Workspace&               ws,
      // WS Output:
             Matrix&                  iy,
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
  rte_std( ws, iy, ppath_transmissions, diy_dvmr, diy_dt,
           ppath, ppath_array, ppath_array_index, f_grid, stokes_dim, 
           emission_agenda, abs_scalar_gas_agenda,
           rte_do_gas_jacs, rte_do_t_jacs, true );
}




/* Workspace method: Doxygen documentation will be auto-generated */
void yUnit(
              Vector&   y,
        const String&   y_unit,
        const Vector&   y_f )
{
  const Index n = y.nelem();

  if( y_f.nelem() != n )
    {
      throw runtime_error( "The length of *y* and *y_f* must be the same" );
    }

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
    {
      for( Index i=0; i<n; i++ )
        { 
          y[i] = invrayjean( y[i], y_f[i] ); 
        }
    }

  else if( y_unit == "PlanckBT" )
    {
      for( Index i=0; i<n; i++ )
        { 
          y[i] = invplanck( y[i], y_f[i] ); 
        }
    }
  else
    {
      ostringstream os;
      os << "Unknown option: y_unit = \"" << y_unit << "\"\n"
         << "Recognised choices are: \"1\", \"RJBT\" and \"PlanckBT\"";
      throw runtime_error( os.str() );      
    }
}



// --------
// New stuff
// --------

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;

//! zaaa2cart
/*! 
   Converts zenith and azimuth angles to a cartesian unit vector.

   This function and the sister function cart2zaaa handles
   transformation of line-of-sights. This in contrast to the sph/poslos
   functions that handles positions, or combinations of positions and
   line-of-sight.

   The cartesian coordinate system used for these two functions can 
   be defined as
    z : za = 0
    x : za=90, aa=0
    y : za=90, aa=90

   \param   dx    Out: x-part of LOS unit vector.
   \param   dy    Out: y-part of LOS unit vector.
   \param   dz    Out: z-part of LOS unit vector.
   \param   za    LOS zenith angle at observation position.
   \param   aa    LOS azimuth angle at observation position.

   \author Patrick Eriksson
   \date   2009-10-02
*/
void zaaa2cart(
             double&   dx,
             double&   dy,
             double&   dz,
       const double&   za,
       const double&   aa )
{
  const double   zarad  = DEG2RAD * za;
  const double   aarad  = DEG2RAD * aa;

  dz = cos( zarad );
  dx = sin( zarad );
  dy = sin( aarad ) * dx;
  dx = cos( aarad ) * dx;
}



//! cart2zaaa
/*! 
   Converts a cartesian directional vector to zenith and azimuth

   This function and the sister function cart2zaaa handles
   transformation of line-of-sights. This in contrast to the sph/poslos
   functions that handles positions, or combinations of positions and
   line-of-sight.

   The cartesian coordinate system used for these two functions can 
   be defined as
    z : za = 0
    x : za=90, aa=0
    y : za=90, aa=90

   \param   za    Out: LOS zenith angle at observation position.
   \param   aa    Out: LOS azimuth angle at observation position.
   \param   dx    x-part of LOS unit vector.
   \param   dy    y-part of LOS unit vector.
   \param   dz    z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2009-10-02
*/
void cart2zaaa(
             double&   za,
             double&   aa,
       const double&   dx,
       const double&   dy,
       const double&   dz )
{
  const double r = sqrt( dx*dx + dy*dy + dz*dz );

  assert( r > 0 );

  za = RAD2DEG * acos( dz / r );
  aa = RAD2DEG * atan2( dy, dx );
}



//! rotationmat3D
/*! 
   Creates a 3D rotation matrix

   Creates a rotation matrix such that R * x, operates on x by rotating 
   x around the origin a radians around line connecting the origin to the 
   point vrot.

   The function is based on rotationmat3D.m, by Belechi (the function 
   is added to atmlab).

   \param   R     Out: Rotation matrix
   \param   vrot  Rotation axis
   \param   a     Rotation angle

   \author Bileschi and Patrick Eriksson
   \date   2009-10-02
*/
void rotationmat3D( 
           Matrix&   R, 
   ConstVectorView   vrot, 
    const Numeric&   a )
{
  assert( R.ncols() == 3 );
  assert( R.nrows() == 3 );
  assert( vrot.nelem() == 3 );

  const double u = vrot[0];
  const double v = vrot[1];
  const double w = vrot[2];

  const double u2 = u * u;
  const double v2 = v * v;
  const double w2 = w * w;

  assert( sqrt( u2 + v2 + w2 ) );

  const double c = cos( DEG2RAD * a );
  const double s = sin( DEG2RAD * a );
  
  // Fill R
  R(0,0) = u2 + (v2 + w2)*c;
  R(0,1) = u*v*(1-c) - w*s;
  R(0,2) = u*w*(1-c) + v*s;
  R(1,0) = u*v*(1-c) + w*s;
  R(1,1) = v2 + (u2+w2)*c;
  R(1,2) = v*w*(1-c) - u*s;
  R(2,0) = u*w*(1-c) - v*s;
  R(2,1) = v*w*(1-c)+u*s;
  R(2,2) = w2 + (u2+v2)*c;
}



/*! 
   Maps MBLOCK_AA_GRID values to correct ZA and AA

   Sensor LOS azimuth angles and mblock_aa_grid values can not be added in a
   straightforward way due to properties of the polar coordinate system used to
   define line-of-sights. This function performs a "mapping" ensuring that the
   pencil beam directions specified by mblock_za_grid and mblock_aa_grid form
   a rectangular grid (on the unit sphere) for any za.

   za0 and aa0 match the angles of the ARTS WSV sensor_los.
   aa_grid shall hold values "close" to 0. The limit is here set to 5 degrees.

   \param   za         Out: Zenith angle matching aa0+aa_grid
   \param   aa         Out: Azimuth angles matching aa0+aa_grid
   \param   za0        Zenith angle
   \param   aa0        Centre azimuth angle
   \param   aa_grid    MBLOCK_AA_GRID values

   \author Patrick Eriksson
   \date   2009-10-02
*/
void map_daa(
             double&   za,
             double&   aa,
       const double&   za0,
       const double&   aa0,
       const double&   aa_grid )
{
  assert( abs( aa_grid ) <= 5 );

  Vector  xyz(3);
  Vector  vrot(3);
  Vector  u(3);

  // Unit vector towards aa0 at za=90
  //
  zaaa2cart( xyz[0], xyz[1], xyz[2], 90, aa0 );
    
  // Find vector around which rotation shall be performed
  // 
  // We can write this as cross([0 0 1],xyz). It turns out that the result 
  // of this operation is just [-y,x,0].
  //
  vrot[0] = -xyz[1];
  vrot[1] = xyz[0];
  vrot[2] = 0;

  // Unit vector towards aa0+aa at za=90
  //
  zaaa2cart( xyz[0], xyz[1], xyz[2], 90, aa0+aa_grid );

  // Apply rotation
  //
  Matrix R(3,3);
  rotationmat3D( R, vrot, za0-90 );
  mult( u, R, xyz );

  // Calculate za and aa for rotated u
  //
  cart2zaaa( za, aa, u[0], u[1], u[2] );
}



void atm_checkedCalc(
         Index&          atm_checked,
   const Index&          atmosphere_dim,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Tensor3&        z_field,
   const Tensor3&        t_field,
   const Matrix&         r_geoid,
   const Matrix&         z_surface,
   const Index&          cloudbox_on, 
   const ArrayOfIndex&   cloudbox_limits )
{
  atm_checked = 1;

  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  // Consistency between dim girds and tmospheric fields/surfaces
  //
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "z_field", z_field, atmosphere_dim, p_grid, lat_grid, 
                                                                     lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid, 
                                                                     lon_grid );
  // vmr_field excluded, could maybe be empty 
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_surface", z_surface, atmosphere_dim, lat_grid, lon_grid );

  // Check that z_field has strictly increasing pages.
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
}



//! apply_y_unit
/*!
    Performs conversion from radiance to other units.

    The function takes radiances in the form of a Tensor3, to cover
    all needs for standard RT calculations. Tensor3 is needed for
    handling conversion of jacobians. Though, conversion of jacobian
    data shall go be made through apply_j_unit (that is mainly an interface
    to this function).

    If the radiance data come as a matrix or a vector, use "telescoping" 
    (see further AUG):
    
      Matrix r
      apply_y_unit( Tensor3View(r), y_unit, f_grid )

    \param   iy       In/Out: Tensor3 with data to be converted, where each 
                      column corresponds to a frequency.
    \param   y_unit   As the WSV.
    \param   f_grid   As the WSV.

    \author Patrick Eriksson 
    \date   2007-10-31
*/
void apply_y_unit( 
      Tensor3View   iy, 
    const String&   y_unit, 
    const Vector&   f_grid )
{
  assert( f_grid.nelem() == iy.ncols() );

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
    {
      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          const Numeric scfac = invrayjean( 1, f_grid[iv] );
          for( Index irow=0; irow<iy.nrows(); irow++ )
            {
              for( Index ipage=0; ipage<iy.npages(); ipage++ )
                {
                  iy(ipage,irow,iv) *= scfac;
                }
            }
        }
    }

  else if( y_unit == "PlanckBT" )
    {
      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          for( Index irow=0; irow<iy.nrows(); irow++ )
            {
              for( Index ipage=0; ipage<iy.npages(); ipage++ )
                {
                  iy(ipage,irow,iv) = invplanck( iy(ipage,irow,iv), f_grid[iv] );
                }
            }
        }

      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          for( Index icol=0; icol<iy.ncols(); icol++ )
            {
            }
        }
    }
  else
    {
      ostringstream os;
      os << "Unknown option: y_unit = \"" << y_unit << "\"\n" 
         << "Recognised choices are: \"1\", \"RJBT\" and \"PlanckBT\"";
      throw runtime_error( os.str() );      
    }

  // If adding more options here, also add in yUnit and jacobianUnit.
}



//! apply_j_unit
/*!
    As apply_y_unit but takes jacobian_unit as input.

    \param   iy       In/Out: Tensor3 with data to be converted, where each 
                      column corresponds to a frequency.
    \param   j_unit   As the WSV.
    \param   f_grid   As the WSV.

    \author Patrick Eriksson 
    \date   2007-10-31
*/
void apply_j_unit( 
      Tensor3View   iy, 
    const String&   jacobian_unit, 
    const Vector&   f_grid )
{
  if( !( jacobian_unit=="1"  ||  jacobian_unit=="RJBT"  ||  
         jacobian_unit=="-" ) )
    {
      ostringstream os;
      os << "Allowed options for *jacobian_unit* are: ""1"", ""RJBT"", and "
         << """-"".";
      throw runtime_error( os.str() );
    }

  apply_y_unit( iy, jacobian_unit, f_grid );
}




// jacobian now sized here.
// n_max_aux. No unit conversion for y_aux, but sensor_response applied
// iy is now "transposed"

void yCalc(
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
         Matrix&                     y_aux,
         Matrix&                     jacobian,
   const Index&                      atm_checked,
   const Index&                      atmosphere_dim,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const Index&                      cloudbox_on,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Vector&                     mblock_za_grid,
   const Vector&                     mblock_aa_grid,
   const Index&                      antenna_dim,
   const Sparse&                     sensor_response,
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Vector&                     sensor_response_za,
   const Vector&                     sensor_response_aa,
   const Agenda&                     iy_agenda,
   const String&                     y_unit,
   const Agenda&                     jacobian_agenda,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const String&                     jacobian_unit )
{
  // Some sizes
  const Index   nf      = f_grid.nelem();
  const Index   nza     = mblock_za_grid.nelem();
        Index   naa     = mblock_aa_grid.nelem();   // Can be set to 1 below
  const Index   n1y     = sensor_response.nrows();
  const Index   nmblock = sensor_pos.nrows();
  const Index   nib     = nf * nza * naa * stokes_dim;


  //---------------------------------------------------------------------------
  // Input checks
  //---------------------------------------------------------------------------

  // Atmosphere OK?
  //
  if( !atm_checked )
    throw runtime_error( "The atmosphere must be flagged to have passed a "
                         "consistency check (atm_checked=1)." );

  // Basic dimensionalities
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );

  // Sensor position and LOS.
  //
  if( sensor_pos.ncols() != atmosphere_dim )
    throw runtime_error( "The number of columns of sensor_pos must be "
                              "equal to the atmospheric dimensionality." );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( "For 1D and 2D, sensor_los shall have one column." );
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
  if( max( sensor_los(joker,0) ) > 180 )
    throw runtime_error( 
       "First column of *sensor_los* is not allowed to have values above 180." );
  if( atmosphere_dim == 2 )
    {
      if( min( sensor_los(joker,0) ) < -180 )
          throw runtime_error( "For atmosphere_dim = 2, first column of "
                      "*sensor_los* is not allowed to have values below -180." );
    }     
  else
    {
      if( min( sensor_los(joker,0)  ) < 0 )
          throw runtime_error( "For atmosphere_dim != 2, first column of "
                         "*sensor_los* is not allowed to have values below 0." );
    }    
  if( atmosphere_dim == 3  &&  max( sensor_los(joker,1) ) > 180 )
    throw runtime_error( 
      "Second column of *sensor_los* is not allowed to have values above 180." );

  // Antenna
  //
  if( nza == 0 )
    throw runtime_error( "The measurement block zenith angle grid is empty." );
  chk_if_increasing( "mblock_za_grid", mblock_za_grid );
  //
  if( antenna_dim == 1 )
    {
      if( naa != 0 )
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
  if( sensor_response.ncols() != nib ) 
    {
      ostringstream os;
      os << "The *sensor_response* matrix does not have the right size,\n"
         << "either the method *sensor_responseInit* has not been run or some\n"
         << "of the other sensor response methods has not been correctly\n"
         << "configured.";
      throw runtime_error( os.str() );
    }

  // Sensor aux variables
  //
  if( n1y != sensor_response_f.nelem()  || n1y != sensor_response_pol.nelem() ||
      n1y != sensor_response_za.nelem() || n1y != sensor_response_za.nelem() )
    {
      ostringstream os;
      os << "Sensor auxiliary variables do not have the correct size.\n"
         << "The following variables should all have same size:\n"
         << "length(y) for one block      : " << n1y << "\n"
         << "sensor_response_f.nelem()    : " << sensor_response_f.nelem()
         << "\nsensor_response_pol.nelem(): " << sensor_response_pol.nelem()
         << "\nsensor_response_za.nelem() : " << sensor_response_za.nelem() 
         << "\nsensor_response_aa.nelem() : " << sensor_response_za.nelem() 
         << "\n";
      throw runtime_error( os.str() );
    }



  //---------------------------------------------------------------------------
  // Allocations and resizing
  //---------------------------------------------------------------------------

  // Effective length of mblock_aa_grid
  //
  if( antenna_dim == 1 )
    { naa = 1; }

  // Resize *y* and *y_XXX*
  //
  y.resize( nmblock*n1y );
  y_f.resize( nmblock*n1y );
  y_pol.resize( nmblock*n1y );
  y_pos.resize( nmblock*n1y, sensor_pos.ncols() );
  y_los.resize( nmblock*n1y, sensor_los.ncols() );
  y_aux.resize( 0, 0 );        // Size can only be determined later

  // Create containers for data of 1 measurement block.
  //
  // We do not know the number of aux variables. The parallelisation makes it
  // unsafe to do the allocation of ib_aux inside the for-loops. We must use a 
  // hard-coded maximum number.
  //
  const Index n_aux_max = 3;
        Index n_aux     = -1;
  //
  Vector ib( nib );
  Matrix ib_aux( nib, n_aux_max );

  // Jacobian variables
  //
  Index  j_analytical_do = 0;
  Matrix jib( 0, 0 );
  //
  if( jacobian_do )
    {
      jacobian.resize( nmblock*n1y, 
                              jacobian_indices[jacobian_indices.nelem()-1][1] );
      //
      Index nx_analytical = 0;
      //
      for( Index i=0; i<jacobian_quantities.nelem(); i++ )
        {
          if( jacobian_quantities[i].Analytical() )
            { 
              j_analytical_do  = 1; 
                nx_analytical += jacobian_indices[i][1]-jacobian_indices[i][0];
            }
        }
      //
      if( j_analytical_do )
        { jib.resize( nib, nx_analytical ); }
    }



  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  for( Index imblock=0; imblock<nmblock; imblock++ )
    {

// #pragma ...

      for( Index iza=0; iza<nza; iza++ )
        {
          // The try block here is necessary to correctly handle
          // exceptions inside the parallel region. 
          try
            {
              for( Index iaa=0; iaa<naa; iaa++ )
                {

                  //--- LOS of interest
                  //
                  Vector los( sensor_los.ncols() );
                  //
                  los     = sensor_los( imblock, joker );
                  los[0] += mblock_za_grid[iza];

                  // Handle za/aa_grid "out-of-bounds" and mapping effects
                  //
                  if( antenna_dim == 2 )
                    {
                      map_daa( los[0], los[1], los[0], los[1], 
                                                          mblock_aa_grid[iaa] );
                    }
                  else if( atmosphere_dim == 1  && abs(los[0]-90) > 90 )
                    {
                      if( los[0] < 0 )          { los[0] = -los[0]; }
                      else if( los[0] > 180 )   { los[0] = 360 - los[0]; }
                    }
                  else if( atmosphere_dim == 2  && abs(los[0]) > 180 )
                    {
                      los[0] = -sign(los[0])*360 + los[0];
                    }
                  else if( atmosphere_dim == 3  &&  abs(los[0]-90) > 90 )
                    {
                      map_daa( los[0], los[1], los[0], los[1], 0 );
                    }

                  cout << los[0] << " " << los[1] << endl;

                  // Calculate iy
                  //
                  Matrix         iy;
                  Tensor3        iy_aux;
                  ArrayOfTensor3 diy_dx;
                  //
                  iy_agendaExecute( ws, iy, iy_aux, diy_dx,
                                    sensor_pos(imblock,joker), los, cloudbox_on, 
                                    jacobian_do, f_grid, t_field, vmr_field,
                                    iy_agenda );

                  // Check sizes
                  //
                  assert( iy.ncols() == nf );
                  assert( iy.nrows() == stokes_dim );
                  //
                  if( n_aux < 0 )
                    { 
                      n_aux = iy_aux.npages(); 
                      if( n_aux > n_aux_max )
                        {
                          ostringstream os;
                          os << "The number of auxilary variables (columns of "
                             << "iy_aux) is hard coded.\nIt is presently set to "
                             << n_aux_max << " variables.";
                          throw runtime_error( os.str() );      
                        }
                    }
                  //
                  assert( iy_aux.ncols() == nf );
                  assert( iy_aux.nrows() == stokes_dim );
                  assert( iy_aux.npages() == n_aux );
                  //
                  if( j_analytical_do )
                    {
                      for( Index i=0; i<jacobian_quantities.nelem(); i++ )
                        {
                          if( jacobian_quantities[i].Analytical() )
                            { 
                              assert( diy_dx[i].ncols() == stokes_dim );
                              assert( diy_dx[i].nrows() == nf );
                              assert( diy_dx[i].npages() == 
                                                        jacobian_indices[i][1] -
                                                        jacobian_indices[i][0] );
                            }
                          else
                            { assert( diy_dx[i].ncols() == 0 ); }
                        }
                    }
                  
                  // iy    : unit conversion and copy to ib
                  // iy_aux: copy to ib_aux
                  //
                  apply_y_unit( Tensor3View(iy), y_unit, f_grid );
                  //
                  const Index row0 =( iza*naa + iaa ) * nf * stokes_dim;
                  //
                  for( Index is=0; is<stokes_dim; is++ )
                    { 
                      ib[Range(row0+is,nf,stokes_dim)] = iy(is,joker); 
                      //
                      for( Index iaux=0; iaux<n_aux; iaux++ )
                        { 
                          ib_aux(Range(row0+is,nf,stokes_dim),iaux) = 
                                                         iy_aux(iaux,is,joker);
                        }
                    }

                  // diy_dx : unit conversion and copy to dib_dx
                  if( j_analytical_do )
                    {
                      // Dummy code !!!!!
                      apply_y_unit( diy_dx[0], jacobian_unit, f_grid );
                    }                  

                }  // End aa loop
            }  // End try

          catch (runtime_error e)
            {
              exit_or_rethrow(e.what());
            }
        }  // End za loop

      // Apply sensor response matrix on ib
      //
      const Index row0 = imblock * n1y;
      //
      mult( y[Range(row0,n1y)], sensor_response, ib );

      // Information and auxilary variables
      //
      for( Index i=0; i<n1y; i++ )
        { 
          y_f[row0+i]         = sensor_response_f[i];
          y_pol[row0+i]       = sensor_response_pol[i]; 
          y_pos(row0+i,joker) = sensor_pos(imblock,joker);
          y_los(row0+i,0)     = sensor_los(imblock,0) + sensor_response_za[i];
          if( sensor_response_aa.nelem() )
            { 
              y_los(row0+i,1) = sensor_los(imblock,0) + sensor_response_aa[i]; 
            }
        }

      // Auxiliary variables
      //
      if( n_aux > 0 )
        {
          if( y_aux.nrows() == 0 )
            { y_aux.resize( nmblock*n1y, n_aux ); }
          //
          mult( y_aux(Range(row0,n1y),joker), sensor_response, 
                                                 ib_aux(joker,Range(0,n_aux)) );
        }

      // Run jacobian_agenda (can be empty)
      if( jacobian_do  &&  jacobian_agenda.nelem() > 0 )
        { jacobian_agendaExecute( ws, jacobian, jacobian_agenda ); }

    }  // End mblock loop
}

