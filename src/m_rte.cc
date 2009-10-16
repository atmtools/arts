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
                  apply_y_unit( Tensor3View(iy), y_unit, f_grid );

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
                      apply_j_unit( 
                            ib_vmr_jacs[ig](Range(nbdone,nf*stokes_dim),joker), 
                                                        j_unit, y_unit, f_grid );
                    }

                  //--- Temperature ---
                  if( rte_do_t_jacs )
                    {
                      //- Map from ppath to retrieval quantities
                      jacobian_from_path_to_rgrids( ib_t_jacs, nbdone, diy_dt,
                                                 0, atmosphere_dim, ppath_array, 
                                                 jacobian_quantities[jqi_t] );

                      //--- Unit conversions
                      apply_j_unit(ib_t_jacs(Range(nbdone,nf*stokes_dim),joker), 
                                                        j_unit, y_unit, f_grid );
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

  jacobianOff( dummy, adummy, jacobian_quantities, jacobian_indices, 
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
    { throw runtime_error( "The length of *y* and *y_f* must be the same" ); }

  apply_y_unit( Tensor3View( MatrixView( y ) ), y_unit, y_f );
}



// --------
// New stuff
// --------

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;


#define FOR_ANALYTICAL_JACOBIANS_DO(what_to_do) \
  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ ) \
    { \
      if( jacobian_quantities[iq].Analytical() ) \
        { what_to_do } \
    } 



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



/*! Maps MBLOCK_AA_GRID values to correct ZA and AA

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



/* Workspace method: Doxygen documentation will be auto-generated */
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



//! iy_transmission_for_scalar_tau
/*!
    Sets *iy_transmission* to match scalar optical thicknesses.

    *iy_transmission* is sized by the function.

    \param   iy_transmission   Out: As the WSV.
    \param   stokes_dim        As the WSV.
    \param   tau               Vector with the optical thickness for each 
                               frequency.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void iy_transmission_for_scalar_tau( 
       Tensor3&     iy_transmission,
  const Index&      stokes_dim,
  ConstVectorView   tau )
{
  iy_transmission.resize( tau.nelem(), stokes_dim, stokes_dim );
  iy_transmission = 0;
  for( Index i=0; i<tau.nelem(); i++ )
    { 
      const Numeric t = exp( -tau[i] );
      for( Index is=0; is<stokes_dim; is++ )
        { 
          iy_transmission(i,is,is) = t;
        }
    }
}



//! iy_transmission_mult_scalar_tau
/*!
    Multiplicates iy_transmission with scalar optical thicknesses.

    The functions can incorporate the transmission of a clear-sky
    path. That is, the transmission can be described by a single value
    The transmission of this path is gives as the optical depth for
    each frequency.

    The "new path" is assumed to be further away from the sensor than 
    the propagtion path already included in iy_transmission. That is,
    the operation can be written as:
    
       Ttotal = Told * Tnew

    where Told is the transmission corresponding to *iy_transmission*
    and Tnew corresponds to *tau*.

    *iy_trans_new* is sized by the function.

    \param   iy_trans_new      Out: Updated version of *iy_transmission*
    \param   iy_transmission   As the WSV.
    \param   tau               Vector with the optical thickness for each 
                               frequency.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void iy_transmission_mult_scalar_tau( 
       Tensor3&      iy_trans_new,
  ConstTensor3View   iy_transmission,
  ConstVectorView    tau )
{
  const Index nf = iy_transmission.npages();
  const Index ns = iy_transmission.ncols();

  assert( ns == iy_transmission.nrows() );
  assert( nf == tau.nelem() );

  iy_trans_new.resize( nf, ns, ns );

  Matrix  tm( ns, ns, 0.0 );

  for( Index iv=0; iv<nf; iv++ )
    {
      const Numeric t = exp( -tau[iv] );
      for( Index is=0; is<ns; is++ )
        { tm(is,is) = t; }
      mult( iy_trans_new(iv,joker,joker), iy_transmission(iv,joker,joker), tm );
    } 
}



//! iy_transmission_mult
/*!
    Multiplicates iy_transmission with (vector) transmissions.

    That is, a multiplication of *iy_transmission* with another
    variable having same structure and holding transmission values.

    The "new path" is assumed to be further away from the sensor than 
    the propagtion path already included in iy_transmission. That is,
    the operation can be written as:
    
       Ttotal = Told * Tnew

    where Told is the transmission corresponding to *iy_transmission*
    and Tnew corresponds to *tau*.

    *iy_trans_new* is sized by the function.

    \param   iy_trans_new      Out: Updated version of *iy_transmission*
    \param   iy_transmission   As the WSV.
    \param   trans             A variable matching *iy_transmission.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void iy_transmission_mult( 
       Tensor3&      iy_trans_new,
  ConstTensor3View   iy_transmission,
  ConstTensor3View   trans )
{
  const Index nf = iy_transmission.npages();
  const Index ns = iy_transmission.ncols();

  assert( ns == iy_transmission.nrows() );
  assert( nf == trans.npages() );
  assert( ns == trans.nrows() );
  assert( ns == trans.ncols() );

  iy_trans_new.resize( nf, ns, ns );

  for( Index iv=0; iv<nf; iv++ )
    {
      mult( iy_trans_new(iv,joker,joker), iy_transmission(iv,joker,joker), 
                                                    trans(iv,joker,joker) );
    } 
}



//! get_ptvmr_for_ppath
/*!
    Determines pressure, temperature and VMR for each propgataion path point.

    The output variables are sized inside the function. For VMR the
    dimensions are [ species, propagation path point ].

    \param   ppath_p     Out: Pressure for each ppath point.
    \param   ppath_t     Out: Temperature for each ppath point.
    \param   ppath_vmr   Out: VMR values for each ppath point.
    \param   ppath             As the WSV.
    \param   atmosphere_dim    As the WSV.
    \param   p_grid            As the WSV.
    \param   lat_grid          As the WSV.
    \param   lon_grid          As the WSV.
    \param   t_field           As the WSV.
    \param   vmr_field         As the WSV.
    \param   f_grid            As the WSV.

    \author Patrick Eriksson 
    \date   2009-10-05
*/
void get_ptvmr_for_ppath( 
        Vector&      ppath_p, 
        Vector&      ppath_t, 
        Matrix&      ppath_vmr, 
  const Ppath&       ppath,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstVectorView    lat_grid,
  ConstVectorView    lon_grid,
  ConstTensor3View   t_field,
  ConstTensor4View   vmr_field )
{
  const Index   np  = ppath.np;

  // Pressure:
  ppath_p.resize(np);
  Matrix itw_p(np,2);
  interpweights( itw_p, ppath.gp_p );      
  itw2p( ppath_p, p_grid, ppath.gp_p, itw_p );
  
  // Temperature:
  ppath_t.resize(np);
  Matrix   itw_field;
  interp_atmfield_gp2itw( itw_field, atmosphere_dim, p_grid, lat_grid, 
                          lon_grid, ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
  interp_atmfield_by_itw( ppath_t,  atmosphere_dim, p_grid, lat_grid, 
                          lon_grid, t_field, ppath.gp_p, 
                          ppath.gp_lat, ppath.gp_lon, itw_field );

  //  VMR fields:
  const Index ns = vmr_field.nbooks();
  ppath_vmr.resize(ns,np);
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( ppath_vmr(is, joker), atmosphere_dim,
             p_grid, lat_grid, lon_grid, vmr_field( is, joker, joker,  joker ), 
             ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
}



//! get_step_vars_for_standardRT
/*!
    Determines variables for each step of "standard" RT integration

    The output variables are sized inside the function. The dimension order is 
       [ frequency, absorption species, ppath point ]

    \param   ws                Out: The workspace
    \param   ppath_abs_sclar   Out: Mean absorption coefficient 
    \param   ppath_emission    Out: Mean emission
    \param   ppath_tau         Out: Optical thickness of each ppath step 
    \param   total_tau         Out: Total optical thickness of path
    \param   abs_scalar_agenda As the WSV.    
    \param   emission_agenda   As the WSV.    
    \param   ppath_p           Pressure for each ppath point.
    \param   ppath_t           Temperature for each ppath point.
    \param   ppath_vmr         VMR values for each ppath point.
    \param   nf                Number of frequencies, length of f_grid
    \param   emission_do       Flag for calculation of emission. Should be
                               set to 0 for pure transmission calculations.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void get_step_vars_for_standardRT( 
        Workspace&   ws,
        Tensor3&     ppath_abs_scalar, 
        Matrix&      ppath_emission, 
        Matrix&      ppath_tau,
        Vector&      total_tau,
  const Agenda&      abs_scalar_agenda,
  const Agenda&      emission_agenda,
  const Ppath&       ppath,
  ConstVectorView    ppath_p, 
  ConstVectorView    ppath_t, 
  ConstMatrixView    ppath_vmr, 
  const Index&       nf,
  const Index&       emission_do )
{
  // Sizes
  const Index   np   = ppath.np;
  const Index   nabs = ppath_vmr.nrows();

  // Init variables
  ppath_abs_scalar.resize( nf, nabs, np );
  ppath_tau.resize( nf, np-1 );
  total_tau.resize( nf );
  total_tau = 0;
  if( emission_do )
    { ppath_emission.resize( nf, np-1 ); }

  // Log of the pressure
  Vector ppath_logp( np );
  transform( ppath_logp, log, ppath_p  );

  for( Index ip=0; ip<np-1; ip++ )
    {
      // Mean of p, t and VMRs for the step
      const Numeric   p_mean = exp( 0.5*( ppath_logp[ip+1]+ppath_logp[ip] ) );
      const Numeric   t_mean = ( ppath_t[ip+1] + ppath_t[ip] ) / 2.0;
            Vector    vmr_mean( nabs );
      for( Index ia=0; ia<nabs; ia++ )
        { vmr_mean[ia] = 0.5 * ( ppath_vmr(ia,ip+1) + ppath_vmr(ia,ip) ); }

      // Calculate emission and absorption terms
      //
      // We must use temporary vectors as the agenda input must be
      // free to be resized
      Vector   evector;
      Matrix   sgmatrix;
      //
      abs_scalar_gas_agendaExecute( ws, sgmatrix, -1, p_mean, t_mean, 
                                                  vmr_mean, abs_scalar_agenda );
      ppath_abs_scalar(joker,joker,ip) = sgmatrix;
      //
      if( emission_do )
        {
          emission_agendaExecute( ws, evector, t_mean, emission_agenda );
          ppath_emission(joker,ip) = evector;
        }

      // Partial and total tau
      //
      for( Index iv=0; iv<nf; iv++ )
        { 
          ppath_tau(iv,ip)  = ppath.l_step[ip] * sgmatrix(iv,joker).sum(); 
          total_tau[iv]    += ppath_tau(iv,ip);
        }
    }
}



//! Help function for analytical jacobian calculations
/*!
    The function determines which terms in jacobian_quantities that are 
    analytical absorption species and temperature jacobians. 

    *abs_species_i* and *is_t* shall be sized to have the same length
    as *jacobian_quantities*. For analytical absorption species
    jacobians, *abs_species_i* is set to the matching index in
    *abs_species*. Otherwise, to -1. For analytical temperature
    jacobians, *is_t* is set to 1. Otherwise to 0.

    \param   abs_species_i         Out: Matching index in abs_species 
    \param   is_t                  Out: Flag for: Is a temperature jacobian?
    \param   jacobian_quantities   As the WSV.
    \param   abs_species           As the WSV.


    \author Patrick Eriksson 
    \date   2009-10-07
*/
void get_pointers_for_analytical_jacobians( 
        ArrayOfIndex&               abs_species_i, 
        ArrayOfIndex&               is_t,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfSpeciesTag&   abs_species )
{
  FOR_ANALYTICAL_JACOBIANS_DO( 
    //
    if( jacobian_quantities[iq].MainTag() == TEMPERATURE_MAINTAG )
      { is_t[iq] = 1; }
    else
      { is_t[iq] = 0; }
    //
    if( jacobian_quantities[iq].MainTag() == ABSSPECIES_MAINTAG )
      {
        ArrayOfSpeciesTag  atag;
        array_species_tag_from_string( atag, jacobian_quantities[iq].Subtag() );
        abs_species_i[iq] = chk_contains( "abs_species", abs_species, atag );
      }
    else
      { abs_species_i[iq] = -1; }
  )
}



//! diy_from_path_to_rgrids
/*!
    Maps jacobian data for points along the propagation path, to
    jacobian retrieval grid data.

    \param   diy_dx              Out: Jacobians for selected retrieval grids.
    \param   jacobian_quantity   As the WSV.
    \param   diy_dpath           Jacobians along the propagation path.
    \param   atmosphere_dim      As the WSV.
    \param   ppath               As the WSV.
    \param   ppath_p             The pressure at each ppath point.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
// A small help function, to make the code below cleaner
void from_dpath_to_dx(
         MatrixView   diy_dx,
    ConstMatrixView   diy_dq,
    const Numeric&    w )
{
  for( Index irow=0; irow<diy_dx.nrows(); irow++ )
    { 
      for( Index icol=0; icol<diy_dx.ncols(); icol++ )
        { diy_dx(irow,icol) += w * diy_dq(irow,icol); }
    }
}
//
void diy_from_path_to_rgrids(
         Tensor3View          diy_dx,
   const RetrievalQuantity&   jacobian_quantity,
   ConstTensor3View           diy_dpath,
   const Index&               atmosphere_dim,
   const Ppath&               ppath,
   ConstVectorView            ppath_p )
{
  // We want here an extrapolation to infinity -> 
  //                                        extremly high extrapolation factor
  const Numeric   extpolfac = 1.0e99;

  // Retrieval grid of interest
  Vector r_grid;

  if( ppath.np > 1 )  // Otherwise nothing to do here
    {
      // Pressure
      r_grid = jacobian_quantity.Grids()[0];
      Index            nr1 = r_grid.nelem();
      ArrayOfGridPos   gp_p(ppath.np);
      p2gridpos( gp_p, r_grid, ppath_p, extpolfac );

      // Latitude
      Index            nr2 = 1;
      ArrayOfGridPos   gp_lat;
      if( atmosphere_dim > 1 )
        {
          gp_lat.resize(ppath.np);
          r_grid = jacobian_quantity.Grids()[1];
          nr2    = r_grid.nelem();
          gridpos( gp_lat, r_grid, ppath.pos(joker,1), extpolfac );
        }

      // Longitude
      Index            nr3 = 1;
      ArrayOfGridPos   gp_lon;
      if( atmosphere_dim > 2 )
        {
          gp_lon.resize(ppath.np);
          r_grid = jacobian_quantity.Grids()[2];
          nr3    = r_grid.nelem();
          gridpos( gp_lon, r_grid, ppath.pos(joker,2), extpolfac );
        }
      
      //- 1D
      if( atmosphere_dim == 1 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              if( gp_p[ip].fd[0] < 1 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[1] );
                }
              if( gp_p[ip].fd[0] > 0 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx+1,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[0] );
                }
            }
        }

      //- 2D
      else if( atmosphere_dim == 2 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[1]*gp_p[ip].fd[1] );
              // Low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[1]*gp_p[ip].fd[0] );
              // High lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[0]*gp_p[ip].fd[1] );
              // High lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[0]*gp_p[ip].fd[0] );
            }
        }

      //- 3D
      else if( atmosphere_dim == 3 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr2*nr1*gp_lon[ip].idx +
                           nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lon, low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // Low lon, low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // Low lon, high lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // Low lon, high lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);

              // Increase *ix* (to be valid for high lon level)
              ix += nr2*nr1;

              // High lon, low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // High lon, low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // High lon, high lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // High lon, high lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);
            }
        }
    }
}



//! vmrunitscf
/*!
    Scale factor for conversion between gas species units.

    The function finds the factor with which the total absorption of a
    gas species shall be multiplicated to match the selected
    (jacobian) unit. 

    \param   x      Out: scale factor
    \param   unit   Unit selected.
    \param   vmr    VMR value.
    \param   p      Pressure
    \param   t      Temperature.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
void vmrunitscf(  
        Numeric&   x, 
  const String&    unit, 
  const Numeric&   vmr,
  const Numeric&   p,
  const Numeric&   t )
{
  if( unit == "rel" )
    { x = 1; }
  else if( unit == "vmr" )
    { x = 1 / vmr; }
  else if( unit == "nd" )
    { x = 1 / ( vmr * number_density( p, t ) ); }
  else
    {
      throw runtime_error( "Allowed options for gas species jacobians are "
                           "\"rel\", \"vmr\" and \"nd\"." );
    }
}



//! get_iy_of_background
/*!
    Determines iy of the "background" of a propgation path.

    The task is to determine *iy* and related variables for the
    background, or to continue the raditiave calculations
    "backwards". The details here depends on the method selected for
    *iy_clearsky_agenda*. 

    Each background is handled by an agenda. Several of these agandes
    can involve recursive calls of *iy_clearsky_agenda*.

    \param   ws                    Out: The workspace
    \param   iy                    Out: As the WSV.
    \param   iy_aux                Out: As the WSV.
    \param   diy_aux               Out: As the WSV.
    \param   iy_transmission       As the WSV.
    \param   iy_aux_do             As the WSV.
    \param   jacobian_do           As the WSV.
    \param   ppath                 As the WSV.
    \param   atmosphere_dim        As the WSV.
    \param   t_field               As the WSV.
    \param   vmr_field             As the WSV.
    \param   cloudbox_on           As the WSV.
    \param   stokes_dim            As the WSV.
    \param   f_grid                As the WSV.
    \param   iy_clearsky_agenda    As the WSV.
    \param   iy_space_agenda       As the WSV.
    \param   surface_prop_agenda   As the WSV.
    \param   iy_cloudbox_agenda    As the WSV.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
void get_iy_of_background(
        Workspace&        ws,
        Matrix&           iy,
        Tensor3&          iy_aux,
        ArrayOfTensor3&   diy_dx,
  ConstTensor3View        iy_transmission,
  const Index&            iy_aux_do,
  const Index&            jacobian_do,
  const Ppath&            ppath,
  const Index&            atmosphere_dim,
  ConstTensor3View        t_field,
  ConstTensor4View        vmr_field,
  const Index&            cloudbox_on,
  const Index&            stokes_dim,
  ConstVectorView         f_grid,
  const Agenda&           iy_clearsky_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           surface_prop_agenda,
  const Agenda&           iy_cloudbox_agenda )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index np      = ppath.np;

  // Set rte_pos, rte_gp_XXX and rte_los to match the last point in ppath.
  // The agendas below use different combinations of these variables.
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  Vector rte_pos;
  Vector rte_los;
  GridPos rte_gp_p;
  GridPos rte_gp_lat;
  GridPos rte_gp_lon;
  rte_pos.resize( atmosphere_dim );
  rte_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  rte_los.resize( ppath.los.ncols() );
  rte_los = ppath.los(np-1,joker);
  gridpos_copy( rte_gp_p, ppath.gp_p[np-1] );
  if( atmosphere_dim > 1 )
    { gridpos_copy( rte_gp_lat, ppath.gp_lat[np-1] ); }
  if( atmosphere_dim > 2 )
    { gridpos_copy( rte_gp_lon, ppath.gp_lon[np-1] ); }

  out3 << "Radiative background: " << ppath.background << "\n";


  // Handle the different background cases
  //
  switch( ppath_what_background( ppath ) )
    {
      // All "interfaces" still work with transposed iy !!!

    case 1:   //--- Space ---------------------------------------------------- 
      {
        iy_space_agendaExecute( ws, iy, rte_pos, rte_los, iy_space_agenda );
        
        if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
          {
            ostringstream os;
            os << "expected size = [" << stokes_dim << "," << nf << "]\n"
               << "size of iy    = [" << iy.nrows() << "," << iy.ncols()<< "]\n"
               << "The size of *iy* returned from *iy_space_agenda* is "
               << "not correct.";
            throw runtime_error( os.str() );      
          }
      }
      break;

    case 2:   //--- The surface -----------------------------------------------
      {
        // Call *surface_prop_agenda*
        //
        Matrix    surface_los;
        Tensor4   surface_rmatrix;
        Matrix    surface_emission;
        //
        surface_prop_agendaExecute( ws, surface_emission, surface_los, 
                                    surface_rmatrix, rte_pos, rte_los, 
                                    rte_gp_p, rte_gp_lat, rte_gp_lon,
                                    surface_prop_agenda );

        // Check output of *surface_prop_agenda*
        //
        const Index   nlos = surface_los.nrows();
        if( nlos )   // if 0, blackbody ground and some checks are not needed
          {
            if( surface_los.ncols() != rte_los.nelem() )
              throw runtime_error( 
                        "Number of columns in *surface_los* is not correct." );
            if( nlos != surface_rmatrix.nbooks() )
              throw runtime_error( 
                  "Mismatch in size of *surface_los* and *surface_rmatrix*." );
            if( surface_rmatrix.npages() != nf )
              throw runtime_error( 
                       "Mismatch in size of *surface_rmatrix* and *f_grid*." );
            if( surface_rmatrix.nrows() != stokes_dim  ||  
                surface_rmatrix.ncols() != stokes_dim ) 
              throw runtime_error( 
              "Mismatch between size of *surface_rmatrix* and *stokes_dim*." );
          }
        if( surface_emission.ncols() != stokes_dim )
          throw runtime_error( 
             "Mismatch between size of *surface_emission* and *stokes_dim*." );
        if( surface_emission.nrows() != nf )
          throw runtime_error( 
                       "Mismatch in size of *surface_emission* and f_grid*." );
        //---------------------------------------------------------------------

        // Variable to hold down-welling radiation
        Tensor3   I( nlos, nf, stokes_dim );
 
        // Loop *surface_los*-es. If no such LOS, we are ready.
        if( nlos > 0 )
          {
            for( Index ilos=0; ilos<nlos; ilos++ )
              {
                // Include surface reflection matrix in *iy_transmission*
                Tensor3 iy_trans_new;
                //
                iy_transmission_mult(  iy_trans_new, iy_transmission, 
                                       surface_rmatrix(ilos,joker,joker,joker) );

                // Calculate downwelling radiation for LOS ilos 
                iy_clearsky_agendaExecute( ws, iy, iy_aux, diy_dx,
                             0, rte_pos, surface_los(ilos,joker), iy_trans_new, 
                             cloudbox_on, jacobian_do, iy_aux_do, 
                             f_grid, t_field, vmr_field, iy_clearsky_agenda );

                I(ilos,joker,joker) = iy;
              }
          }

        // Add up
        surface_calc( iy, I, surface_los, surface_rmatrix, surface_emission );
      }
      break;


    case 3:   //--- Cloudbox surface or interior --------------------------------
    case 4:
      {
        // Pass a local copy of ppath to the cloudbox agenda
        Ppath   ppath_local;
        ppath_init_structure( ppath_local, atmosphere_dim, ppath.np );
        ppath_copy( ppath_local, ppath );

        // The cloudbox agenda is not yet handling iy_aux and iy_transmission.
        // And is returning *iy* transposed
        //
        iy_cloudbox_agendaExecute( ws, iy, ppath_local,
                                   rte_pos, rte_los, rte_gp_p,
                                   rte_gp_lat, rte_gp_lon,
                                   iy_cloudbox_agenda );

        if( iy.nrows() != nf  ||  iy.ncols() != stokes_dim )
          {
            out1 << "expected size = [" << nf << "," << stokes_dim << "]\n";
            out1 << "iy size = [" << iy.nrows() << "," << iy.ncols()<< "]\n";
            throw runtime_error( "The size of *iy* returned from "
                                      "*iy_cloudbox_agenda* is not correct." );
          }
      }
      break;

    default:  //--- ????? ----------------------------------------------------
      // Are we here, the coding is wrong somewhere
      assert( false );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyMC(
        Workspace&                  ws,
        Matrix&                     iy,
        Tensor3&                    iy_aux,
        ArrayOfTensor3&             diy_dx,
  const Index&                      iy_agenda_call1,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      iy_aux_do,
  const Index&                      jacobian_do,
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
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     abs_scalar_gas_agenda, 
  const Agenda&                     opt_prop_gas_agenda,
  const Tensor4&                    pnd_field,
  const Tensor3&                    iy_transmission,
  const String&                     y_unit,
  const Numeric&                    mc_std_err,
  const Index&                      mc_max_time,
  const Index&                      mc_max_iter,
  const Index&                      mc_z_field_is_1D )
{
  // Throw error if unsupported features are requested
  if( atmosphere_dim != 3 )
    throw runtime_error( 
                  "Only 3D atmospheres are allowed (atmosphere_dim must be 3)" );
  if( !cloudbox_on )
    throw runtime_error( 
                      "The cloudbox must be activated (cloudbox_on must be 1)" );
  if( iy_aux_do )
    throw runtime_error( 
        "This method does not provide any auxilary data (iy_aux_do must be 0)" );
  if( jacobian_do )
    throw runtime_error( 
          "This method does not provide any jacobians (jacobian_do must be 0)" );
  if( !iy_agenda_call1 )
    throw runtime_error( 
                    "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );

  // Size output variables
  //
  const Index   nf  = f_grid.nelem();
  //
  iy.resize( nf, stokes_dim );
  //
  iy_aux.resize( 0, 0, 0 );
  diy_dx.resize(0);

  // Some MC variables are only local here
  Tensor3  mc_points;
  Index    mc_iteration_count, mc_seed;
  Vector   mc_error;
  //
  MCAntenna mc_antenna;
  mc_antenna.set_pencil_beam();

  // Pos and los must be matrices 
  Matrix pos(1,3), los(1,2);
  //
  pos(0,joker) = rte_pos;
  los(0,joker) = rte_los;

  for( Index f_index=0; f_index<nf; f_index++ )
    {
      ArrayOfSingleScatteringData   scat_data_mono;
      
      scat_data_monoCalc( scat_data_mono, scat_data_raw, 
                                          f_grid, f_index );

      // Seed reset for each loop. If not done, the errors 
      // appear to be highly correlated.
      MCSetSeedFromTime( mc_seed );

      Vector y;
                  
      MCGeneral( ws, y, mc_iteration_count, mc_error, mc_points, 
                 mc_antenna, f_grid, f_index, pos, los, stokes_dim, 
                 iy_space_agenda, surface_prop_agenda, opt_prop_gas_agenda,
                 abs_scalar_gas_agenda, p_grid, lat_grid, lon_grid, z_field, 
                 r_geoid, z_surface, t_field, vmr_field, cloudbox_limits, 
                 pnd_field, scat_data_mono, mc_seed, y_unit, mc_std_err, 
                 mc_max_time, mc_max_iter, mc_z_field_is_1D ); 

      assert( y.nelem() == stokes_dim );
 
      // Data returned can not be in Tb
      if ( y_unit == "RJBT" )
        { 
          const Numeric scfac = invrayjean( 1, f_grid[f_index] );
          y /= scfac;
        }
     
      iy(f_index,joker) = y;
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandardClearsky(
        Workspace&                  ws,
        Matrix&                     iy,
        Tensor3&                    iy_aux,
        ArrayOfTensor3&             diy_dx,
  const Index&                      iy_agenda_call1,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      iy_aux_do,
  const Index&                      jacobian_do,
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
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Agenda&                     ppath_step_agenda,
  const Agenda&                     emission_agenda,
  const Agenda&                     abs_scalar_agenda,
  const Agenda&                     iy_clearsky_agenda,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     iy_cloudbox_agenda,
  const Tensor3&                    iy_transmission,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices )
{
  // A minimal amount of checks, for efficiency reasons
  assert( rte_pos.nelem() == atmosphere_dim );
  assert( ( atmosphere_dim < 3   &&  rte_los.nelem() == 1 )  ||
          ( atmosphere_dim == 3  &&  rte_los.nelem() == 2 ) );

  // Determine if there are any jacobians to handle
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }

  // If primary call, initilise iy_aux and diy_dx
  //
  const Index   nf  = f_grid.nelem();
  //
  if( iy_agenda_call1 )
    {
      if( iy_aux_do ) { iy_aux.resize( 3, nf, stokes_dim ); iy_aux = 0; }
      else            { iy_aux.resize( 0, 0, 0 ); }
      //
      if( j_analytical_do ) 
        {
          diy_dx.resize( jacobian_indices.nelem() ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dx[iq].resize( jacobian_indices[iq][1] - 
                               jacobian_indices[iq][0] + 1, nf, stokes_dim ); 
            diy_dx[iq] = 0.0;
          ) 
        }
      else { diy_dx.resize( 0 ); }
    }

  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
              lat_grid, lon_grid, z_field, r_geoid, z_surface,
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 1 );
  //
  const Index   i_background = ppath_what_background( ppath );


  // Get quantities required for each ppath point, and update iy_transmission
  //
  // If np = 1, we only need to determine the radiative background
  //
  const Index     np  = ppath.np;
        Vector    ppath_p, ppath_t, total_tau;
        Matrix    ppath_vmr, ppath_emission, ppath_tau;
        Tensor3   ppath_abs_scalar, iy_trans_new;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      //
      get_ptvmr_for_ppath( ppath_p, ppath_t, ppath_vmr, ppath, atmosphere_dim, 
                           p_grid, lat_grid, lon_grid, t_field, vmr_field );

      // Get emission, absorption, optical thickness for each step, and total
      // optical thickness
      //
      get_step_vars_for_standardRT( ws, ppath_abs_scalar, ppath_emission, 
                                    ppath_tau, total_tau, 
                                    abs_scalar_agenda, emission_agenda,
                                    ppath, ppath_p, ppath_t, ppath_vmr, nf, 1 );
    }

  // Handle iy_transmission (the variable not needed when background is space)
  //

  if( i_background > 1 || iy_aux_do )
    {
      if( iy_agenda_call1 )
        { iy_transmission_for_scalar_tau( iy_trans_new, stokes_dim, 
                                                                  total_tau ); }
      else
        { iy_transmission_mult_scalar_tau( iy_trans_new, iy_transmission, 
                                                                  total_tau ); }
    }

  // Radiative background
  //
  get_iy_of_background( ws, iy, iy_aux, diy_dx, iy_trans_new,  
                      iy_aux_do, jacobian_do, ppath, atmosphere_dim, t_field,
                      vmr_field, cloudbox_on, stokes_dim, f_grid, iy_clearsky_agenda,
                      iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda );
  
  // Set iy_aux 
  //
  if( iy_aux_do ) 
    {
      // Fill transmission if background is TOA or surface
      if( i_background <= 2 )
        {
          for( Index iv=0; iv<nf; iv++ )
            { 
              for( Index is=0; is<stokes_dim; is++ )
                { iy_aux( 0, iv, is ) = iy_trans_new( iv, is, is ); }
            }
        }
      //
      // Set background if primary call
      if( iy_agenda_call1 )
        { iy_aux( 1, joker, joker ) = (Numeric)i_background; }
      // 
      // Flag intersection with cloudbox
      if( i_background >= 3 )
        { iy_aux( 2, joker, joker ) = 1.0; }
    }    

  // Do RT calculations
  //
  if( np > 1 )
    {
      // Create container for the derivatives with respect to changes
      // at each ppath point. And find matching abs_species-index and 
      // "temperature flag" (for analytical jacobians).
      //
      ArrayOfTensor3  diy_dpath; 
      ArrayOfIndex    abs_species_i, is_t; 
      //
      const Numeric   dt = 0.1;
            Matrix    ppath_e2;
            Tensor3   ppath_as2;
      //
      if( j_analytical_do )
        { 
          diy_dpath.resize( jacobian_indices.nelem() ); 
          abs_species_i.resize( jacobian_indices.nelem() ); 
          is_t.resize( jacobian_indices.nelem() ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dpath[iq].resize( np, nf, stokes_dim ); 
            diy_dpath[iq] = 0.0;
          )
          get_pointers_for_analytical_jacobians( abs_species_i, is_t, 
                                              jacobian_quantities, abs_species );
          //
          // Get emission and absorption for disturbed temperature
          Index do_t=0;
          for( Index iq=0; iq<is_t.nelem() && do_t==0; iq++ )
            { if( is_t[iq] ) { do_t = 1; } }
          if( do_t )
            {
              Matrix mdummy; Vector vdummy, t2=ppath_t;
              for( Index i=0; i<t2.nelem(); i++ ) { t2[i] += dt; }
              get_step_vars_for_standardRT( ws, ppath_as2, ppath_e2, mdummy, 
                                   vdummy, abs_scalar_agenda, emission_agenda,
                                   ppath, ppath_p, t2, ppath_vmr, nf, 1 );
            }
        }

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Loop frequencies
          for( Index iv=0; iv<nf; iv++ )
            { 
              // Jacobian quantities
              if( j_analytical_do )
                {
                  const Numeric A = ppath.l_step[ip] * exp( -total_tau[iv] );
                  const Numeric B = 0.5 * A * ( ppath_emission(iv,ip)-iy(iv,0) );
                  
                  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ ) 
                    {
                      // Absorbing species
                      const Index isp = abs_species_i[iq]; 
                      if( isp >= 0 )
                        {
                          // Scaling factors to handle retrieval unit
                          Numeric unitscf1, unitscf2;
                          vmrunitscf( unitscf1, jacobian_quantities[iq].Mode(), 
                                   ppath_vmr(isp,ip), ppath_p[ip], ppath_t[ip] );
                          vmrunitscf( unitscf2, jacobian_quantities[iq].Mode(), 
                             ppath_vmr(isp,ip+1), ppath_p[ip+1], ppath_t[ip+1] );
                          //
                          // Stokes component 1
                          const Numeric w0 = B * ppath_abs_scalar(iv,isp,ip);
                          diy_dpath[iq](ip  ,iv,0) += unitscf1 * w0;
                          diy_dpath[iq](ip+1,iv,0) += unitscf2 * w0;
                          //
                          // Higher stokes components
                          for( Index is=1; is<stokes_dim; is++ )
                            { 
                              const Numeric wi = -0.5 * A * iy(iv,is) *
                                                     ppath_abs_scalar(iv,isp,ip);
                              diy_dpath[iq](ip  ,iv,is) += unitscf1 * wi;
                              diy_dpath[iq](ip+1,iv,is) += unitscf2 * wi;
                            }
                        }

                      // Temperature
                      else if( is_t[iq] )
                        {
                          const Numeric dkdt = ( ppath_as2(iv,joker,ip).sum() -
                                      ppath_abs_scalar(iv,joker,ip).sum() ) / dt;
                          const Numeric dbdt = ( ppath_e2(iv,ip) -
                                                    ppath_emission(iv,ip) ) / dt;
                          //
                          // Stokes component 1
                          const Numeric w0 = B * dkdt + ( 
                                        exp(-(total_tau[iv]-ppath_tau(iv,ip))) - 
                                        exp( -total_tau[iv]) ) * dbdt;
                          diy_dpath[iq](ip  ,iv,0) += w0;
                          diy_dpath[iq](ip+1,iv,0) += w0;
                          //
                          // Higher stokes components
                          for( Index is=1; is<stokes_dim; is++ )
                            { 
                              const Numeric wi = -0.5 * A * iy(iv,is) * dkdt;
                              diy_dpath[iq](ip  ,iv,is) += wi;
                              diy_dpath[iq](ip+1,iv,is) += wi;
                            }
                        }
                    }
                  
                  // Update total_tau (not used for spectrum calculation)
                  total_tau[iv] -= ppath_tau(iv,ip);
                }

              // Spectrum
              //
              const Numeric step_tr = exp( -ppath_tau(iv,ip) );
              //
              iy(iv,0)  = iy(iv,0)*step_tr + ppath_emission(iv,ip)*(1-step_tr); 
              //
              for( Index is=1; is<stokes_dim; is++ )
                { iy(iv,is) *= step_tr; }
            }
        } 

      // Map jacobians from ppath to retrieval grids
      if( j_analytical_do )
        { 
          // Weight with iy_transmission
          if( !iy_agenda_call1 )
            {
              FOR_ANALYTICAL_JACOBIANS_DO( 
                for( Index iv=0; iv<nf; iv++ )
                  { 
                    Matrix X;
                    X = diy_dpath[iq](joker,iv,joker);
                    mult( diy_dpath[iq](joker,iv,joker), 
                                            iy_transmission(iv,joker,joker), X );
                  }
               )
            }

          // Map to retrieval grids
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                                diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
          )
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */

/*
void iyTransmissionStandardClearsky(
        Workspace&                  ws,
        Matrix&                     iy,
        Tensor3&                    iy_aux,
        ArrayOfTensor3&             diy_dx,
  const Index&                      iy_agenda_call1,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      iy_aux_do,
  const Index&                      jacobian_do,
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
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Agenda&                     ppath_step_agenda,
  const Agenda&                     abs_scalar_agenda,
  const Agenda&                     iy_clearsky_agenda,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     iy_cloudbox_agenda,
  const Tensor3&                    iy_transmission,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices )
{
  // A minimal amount of checks, for efficiency reasons
  assert( rte_pos.nelem() == atmosphere_dim );
  assert( ( atmosphere_dim < 3   &&  rte_los.nelem() == 1 )  ||
          ( atmosphere_dim == 3  &&  rte_los.nelem() == 2 ) );

  // Determine if there are any jacobians to handle
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }

  // If primary call, initilise iy_aux and diy_dx
  //
  const Index     nf  = f_grid.nelem();
  //
  if( iy_aux.nrows() == 0 )
    { 
      if( !iy_aux_do ) { iy_aux.resize( 0, 0, 0 ); }
      else             { iy_aux.resize( 3, stokes_dim, nf ); iy_aux = 0; }
    }
  //
  if( iy_agenda_call1 )
    { 
      if( !j_analytical_do ) { diy_dx.resize( 0 ); }
      else
        {
          diy_dx.resize( jacobian_indices.nelem() ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dx[iq].resize( jacobian_indices[iq][1]-jacobian_indices[iq][0], 
                                                              stokes_dim, nf ); 
            diy_dx[iq] = 0.0;
          ) 
        }
    }

  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
              lat_grid, lon_grid, z_field, r_geoid, z_surface,
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 1 );

  // Get quantities required for each ppath point, and update iy_transmission
  //
  // If np = 1, we only need to determine the radiative background
  //
  const Index     np  = ppath.np;
        Vector    ppath_p, ppath_t, total_tau;
        Matrix    ppath_vmr, dummy, ppath_tau;
        Tensor3   ppath_abs_scalar, iy_trans_new;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      //
      get_ptvmr_for_ppath( ppath_p, ppath_t, ppath_vmr, ppath, atmosphere_dim, 
                           p_grid, lat_grid, lon_grid, t_field, vmr_field );

      // Get absorption, optical thickness for each step, and total
      // optical thickness 
      //
      get_step_vars_for_standardRT( ws, ppath_abs_scalar, dummy, 
                                    ppath_tau, total_tau, 
                                    abs_scalar_agenda, Agenda(),
                                    ppath, ppath_p, ppath_t, ppath_vmr, nf, 0 );

      // New iy_transmission
      iy_transmission_mult_scalar_tau( iy_trans_new, iy_transmission, total_tau);
    }
  else
    {
      // If np=1, *iy_trans_new* can be identical to *iy_transmission*
      //
      iy_trans_new = iy_transmission;   // How to avoid copying of values?
    }

  // Radiative background
  //
  Index i_background = 1;
  //
  get_iy_of_background( ws, iy, iy_aux, iy_transmission, diy_dx, 
                      iy_aux_do, jacobian_do, ppath, atmosphere_dim, t_field,
                      vmr_field, cloudbox_on, stokes_dim, f_grid, iy_clearsky_agenda,
                      iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda );
  
  // Set iy_aux 
  //
  if( iy_aux_do ) 
    {
      // Set background if primary call
      if( iy_agenda_call1 )
        { iy_aux( 0, joker, joker ) = (Numeric)i_background; }
      // 
      // Flag intersection with cloudbox
      if( i_background >= 3 )
        { iy_aux( 1, joker, joker ) = 1.0; }
    }    

  // Do RT calculations
  //
  if( np > 1 )
    {
      for( Index iv=0; iv<nf; iv++ )
        { 
          const Numeric path_tr = exp( -total_tau[iv] );
          //
          for( Index is=0; is<stokes_dim; is++ )
            { iy(is,iv) *= path_tr; }
        }      
    }
}
*/



//! ib_calc
/*!
    ...


    \author Patrick Eriksson 
    \date   2009-10-16
*/
void ib_calc(
        Workspace&                  ws,
        Vector&                     ib,
        Matrix&                     ib_aux,
        Index&                      n_aux,
        ArrayOfMatrix&              dib_dx,
  const Index&                      imblock,
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
  const Agenda&                     iy_clearsky_agenda,
  const Index&                      iy_aux_do,
  const String&                     y_unit,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const String&                     jacobian_unit )
{
  // Sizes
  const Index   nf      = f_grid.nelem();
  const Index   nza     = mblock_za_grid.nelem();
        Index   naa     = mblock_aa_grid.nelem();   
  if( antenna_dim == 1 )  
    { naa = 1; }
  const Index   nib     = nf * nza * naa * stokes_dim;

  // Create containers for data of 1 measurement block.
  //
  // We do not know the number of aux variables. The parallelisation makes it
  // unsafe to do the allocation of ib_aux inside the for-loops. We must use a 
  // hard-coded maximum number.
  //
  const Index n_aux_max = 3;
  n_aux     = -1;   // -1 flags value yet unknown
  //
  ib.resize( nib );
  ib_aux.resize( nib, n_aux_max );
  //
  if( j_analytical_do )
    {
      dib_dx.resize( jacobian_indices.nelem() );
      FOR_ANALYTICAL_JACOBIANS_DO(
        dib_dx[iq].resize( nib, jacobian_indices[iq][1] -
                                jacobian_indices[iq][0] + 1 );
      )
    }
  else
    { dib_dx.resize( 0 ); }


  // Start of actual calculations

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
                { map_daa( los[0], los[1], los[0], los[1], 
                                                    mblock_aa_grid[iaa] ); }
              else if( atmosphere_dim == 1  && abs(los[0]-90) > 90 )
                { if( los[0] < 0 )          { los[0] = -los[0]; }
                  else if( los[0] > 180 )   { los[0] = 360 - los[0]; } }
              else if( atmosphere_dim == 2  && abs(los[0]) > 180 )
                { los[0] = -sign(los[0])*360 + los[0]; }
              else if( atmosphere_dim == 3  &&  abs(los[0]-90) > 90 )
                { map_daa( los[0], los[1], los[0], los[1], 0 ); }


              // Calculate iy
              //
              Matrix         iy;
              Tensor3        iy_aux, iy_transmission;
              ArrayOfTensor3 diy_dx; 
              //
              iy_clearsky_agendaExecute( ws, iy, iy_aux, diy_dx,
                        1, sensor_pos(imblock,joker), los, iy_transmission, 
                        cloudbox_on, j_analytical_do, iy_aux_do, 
                        f_grid, t_field, vmr_field, iy_clearsky_agenda );

              // Check sizes
              //
              assert( iy.ncols() == stokes_dim );
              assert( iy.nrows() == nf );
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
              if( n_aux )
                { 
                  assert( iy_aux.ncols() == stokes_dim );
                  assert( iy_aux.nrows() == nf );
                  assert( iy_aux.npages() == n_aux );
                }
              //
              if( j_analytical_do )
                {
                  FOR_ANALYTICAL_JACOBIANS_DO(
                    assert( diy_dx[iq].ncols() == stokes_dim );
                    assert( diy_dx[iq].nrows() == nf );
                    assert( diy_dx[iq].npages() == jacobian_indices[iq][1] -
                                               jacobian_indices[iq][0] + 1 );
                  )
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
                  ib[Range(row0+is,nf,stokes_dim)] = iy(joker,is); 
                  //
                  for( Index iaux=0; iaux<n_aux; iaux++ )
                    { 
                      ib_aux(Range(row0+is,nf,stokes_dim),iaux) = 
                                                       iy_aux(iaux,joker,is);
                    }
                }

              // Jacobian part
              if( j_analytical_do )
                {
                  // Basically copy calculations for *iy*
                  FOR_ANALYTICAL_JACOBIANS_DO(
                    //
                    apply_j_unit( diy_dx[iq], jacobian_unit, y_unit, f_grid );
                    //
                    for( Index ip=0; ip<jacobian_indices[iq][1] -
                                        jacobian_indices[iq][0]+1; ip++ )
                      {
                        for( Index is=0; is<stokes_dim; is++ )
                          { 
                            dib_dx[iq](Range(row0+is,nf,stokes_dim),ip)=
                                                     diy_dx[iq](ip,joker,is); 
                          }
                      }                              
                  )
                }
            }  // End aa loop
        }  // End try

      catch (runtime_error e)
        {
          exit_or_rethrow(e.what());
        }
    }  // End za loop


 
}

/* Workspace method: Doxygen documentation will be auto-generated */
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
  const Agenda&                     iy_clearsky_agenda,
  const Index&                      iy_aux_do,
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
        Index   naa     = mblock_aa_grid.nelem();   
  if( antenna_dim == 1 )  
    { naa = 1; }
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

  // Resize *y* and *y_XXX*
  //
  y.resize( nmblock*n1y );
  y_f.resize( nmblock*n1y );
  y_pol.resize( nmblock*n1y );
  y_pos.resize( nmblock*n1y, sensor_pos.ncols() );
  y_los.resize( nmblock*n1y, sensor_los.ncols() );
  y_aux.resize( 0, 0 );        // Size can only be determined later

  // Jacobian variables
  //
  Index  j_analytical_do = 0;
  //
  if( jacobian_do )
    {
      jacobian.resize( nmblock*n1y, 
                            jacobian_indices[jacobian_indices.nelem()-1][1]+1 );
      //
      FOR_ANALYTICAL_JACOBIANS_DO(
        j_analytical_do  = 1; 
      )
    }


  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  for( Index imblock=0; imblock<nmblock; imblock++ )
    {
      // Calculate monochromatic pencil beam data for 1 measurement block
      //
      Vector          ib, yb(n1y);
      Matrix          ib_aux;
      Index           n_aux;
      ArrayOfMatrix   dib_dx;      
      //
      ib_calc( ws, ib, ib_aux, n_aux, dib_dx, imblock, atmosphere_dim, 
               t_field, vmr_field, cloudbox_on, stokes_dim, f_grid, 
               sensor_pos, sensor_los, mblock_za_grid, mblock_aa_grid, 
               antenna_dim, iy_clearsky_agenda, iy_aux_do, y_unit, 
               j_analytical_do, jacobian_quantities, jacobian_indices, 
               jacobian_unit );

      // Apply sensor response matrix on ib
      //
      mult( yb, sensor_response, ib );
      //
      const Index row0 = imblock * n1y;
      //
      y[Range(row0,n1y)] = yb;

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

      // dib_dx part of *jacobian*
      //
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO(
            mult( jacobian(Range(row0,n1y), Range(jacobian_indices[iq][0],
                          jacobian_indices[iq][1]-jacobian_indices[iq][0]+1)),
                                                   sensor_response, dib_dx[iq] );
          )
        }

      // Rest of *jacobian*: run jacobian_agenda (can be empty)
      //
      if( jacobian_do  &&  jacobian_agenda.nelem() > 0 )
        { jacobian_agendaExecute( ws, jacobian, jacobian_agenda ); }

    }  // End mblock loop
}

