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
  \file   rte.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-29

  \brief  Functions to solve radiative transfer tasks.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "auto_md.h"
#include "check_input.h"
#include "logic.h"
#include "math_funcs.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "lin_alg.h"




/*===========================================================================
  === The functions in alphabetical order
  ===========================================================================*/

//! apply_y_unit
/*!
    Performs conversion from radiance to other units, following the keyword
    argument *y_unit* used in the set of RteCalc functions.

    \param   iy       In/Out: Matrix with data to be converted, where each 
                      row corresponds to a frequency.
    \param   y_unit   As the keyword argument for *RteCalc*.
    \param   f_grid   Frequency grid.

    \author Patrick Eriksson 
    \date   2007-10-31
*/
void apply_y_unit( 
       MatrixView   iy, 
    const String&   y_unit, 
    const Vector&   f_grid )
{
  assert( f_grid.nelem() == iy.nrows() );

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
    {
      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          const Numeric scfac = invrayjean( 1, f_grid[iv] );
          for( Index icol=0; icol<iy.ncols(); icol++ )
            {
              iy(iv,icol) *= scfac;
            }
        }
    }

  else if( y_unit == "PlanckBT" )
    {
      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          for( Index icol=0; icol<iy.ncols(); icol++ )
            {
              iy(iv,icol) = invplanck( iy(iv,icol), f_grid[iv] );
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



//! apply_y_unit_single
/*!
    A version of apply_y_unit handle monochormatic input. Just an interface
    to apply_y_unit

    \param   i        In/Out: Vector with data to be converted, where each 
                      position corresponds to a Stokes dimension.
    \param   y_unit   As the keyword argument for *RteCalc*.
    \param   f        Frequency value.

    \author Patrick Eriksson 
    \date   2007-10-31
*/
void apply_y_unit_single( 
          Vector&   i, 
    const String&   y_unit, 
    const Numeric&  f )
{
  // Create frequency grid vector
  Vector f_grid(1,f);

  // Create matrix matching WSV iy
  Matrix iy(1,i.nelem());
  iy(0,joker) = i;

  apply_y_unit( iy, y_unit, f_grid );

  i = iy(0,joker);
}



//! get_radiative_background
/*!
  Sets *iy* to the radiative background for a propagation path.

  The function uses *ppath* to determine the radiative background
  for a propagation path and calls the relevant agenda. Coding of
  backgrounds is described in the header of the function
  ppath_set_background (in ppath.cc).

  The main purpose of the function is set *iy*. It is NOT needed to set 
  *iy* to the correct size before calling the function.

  \param[in,out] ws Current workspace
  \param[out] iy
  \param[out] ppath
  \param[out] ppath_array_index
  \param[out] ppath_array
  \param[out] diy_dvmr
  \param[out] diy_dt
  \param[in]  ppath_step_agenda
  \param[in]  rte_agenda
  \param[in]  surface_prop_agenda
  \param[in]  iy_space_agenda
  \param[in]  iy_cloudbox_agenda
  \param[in]  atmosphere_dim
  \param[in]  p_grid
  \param[in]  lat_grid
  \param[in]  lon_grid
  \param[in]  z_field
  \param[in]  t_field
  \param[in]  vmr_field
  \param[in]  r_geoid
  \param[in]  z_surface
  \param[in]  cloudbox_on
  \param[in]  cloudbox_limits
  \param[in]  f_grid
  \param[in]  stokes_dim
  \param[in]  ppath_array_do
  \param[in]  rte_do_vmr_jacs
  \param[in]  rte_do_t_jacs

  \author Patrick Eriksson 
  \date   2002-09-17
*/
void get_radiative_background(
              Workspace&               ws,
              Matrix&                  iy,
              Ppath&                   ppath,
              Index&                   ppath_array_index,
              ArrayOfPpath&            ppath_array,
              ArrayOfTensor4&          diy_dvmr,
              ArrayOfTensor4&          diy_dt,
        const Agenda&                  ppath_step_agenda,
        const Agenda&                  rte_agenda,
        const Agenda&                  surface_prop_agenda,
        const Agenda&                  iy_space_agenda,
        const Agenda&                  iy_cloudbox_agenda,
        const Index&                   atmosphere_dim,
        const Vector&                  p_grid,
        const Vector&                  lat_grid,
        const Vector&                  lon_grid,
        const Tensor3&                 z_field,
        const Tensor3&                 t_field,
        const Tensor4&                 vmr_field,
        const Matrix&                  r_geoid,
        const Matrix&                  z_surface,
        const Index&                   cloudbox_on, 
        const ArrayOfIndex&            cloudbox_limits,
        const Vector&                  f_grid,
        const Index&                   stokes_dim,
        const Index&                   ppath_array_do,
        const ArrayOfIndex&            rte_do_vmr_jacs,
        const Index&                   rte_do_t_jacs )
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


  // Initialize iy to the radiative background
  switch ( ppath_what_background( ppath ) )
    {

    case 1:   //--- Space ---------------------------------------------------- 
      {
        chk_not_empty( "iy_space_agenda", iy_space_agenda );

        iy_space_agendaExecute( ws, iy, rte_pos, rte_los, 
                                iy_space_agenda );
     
        if( iy.nrows() != nf  ||  iy.ncols() != stokes_dim )
          {
            out1 << "expected size = [" << nf << "," << stokes_dim << "]\n";
            out1 << "iy size = [" << iy.nrows() << "," << iy.ncols()<< "]\n";
            throw runtime_error( "The size of *iy* returned from "
                                          "*iy_space_agenda* is not correct.");
          }
      }
      break;


    case 2:   //--- The surface -----------------------------------------------
      {
        // Call *surface_prop_agenda*
        //
        chk_not_empty( "surface_prop_agenda", surface_prop_agenda );
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
            if( surface_rmatrix.ncols() != stokes_dim  ||  
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
            // Make local version of *ppath* that later must be restored 
            Ppath   pp_copy;
            ppath_init_structure( pp_copy, atmosphere_dim, ppath.np );
            ppath_copy( pp_copy, ppath );
            //
            const Index  pai = ppath_array_index;

            for( Index ilos=0; ilos<nlos; ilos++ )
              {
                // Calculate downwelling radiation for LOS ilos 
                iy_calc( ws, iy, ppath, ppath_array_index, ppath_array,
                   diy_dvmr, diy_dt, ppath_step_agenda, rte_agenda, 
                   iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda, 
                   atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                   t_field, vmr_field, r_geoid, z_surface, cloudbox_on, 
                   cloudbox_limits, rte_pos, surface_los(ilos,joker), f_grid, 
                   stokes_dim, ppath_array_do, rte_do_vmr_jacs, rte_do_t_jacs );

                I(ilos,joker,joker) = iy;

                // Include reflection matrix in jacobian quantities
                //
                // Assume polarisation effects in surface_rmatrix
                // (not of any impoartnce if stokes_dim=1)
                const bool pol_r = true;
                //
                if( rte_do_vmr_jacs.nelem() )
                  {
                    for( Index iv=0; iv<nf; iv++ )
                      {
                        include_trans_in_diy_dq( diy_dvmr, iv, pol_r,
                                          surface_rmatrix(ilos,iv,joker,joker),
                                          ppath_array, ppath_array_index );
                      }
                  }
                if( rte_do_t_jacs )
                  {
                    for( Index iv=0; iv<nf; iv++ )
                      {
                        include_trans_in_diy_dq( diy_dt, iv, pol_r,
                                          surface_rmatrix(ilos,iv,joker,joker),
                                          ppath_array, ppath_array_index );
                      }
                  }

                // Reset *ppath_array_index*
                ppath_array_index = pai;
              }

            // Copy data back to *ppath*.
            ppath_init_structure( ppath, atmosphere_dim, pp_copy.np );
            ppath_copy( ppath, pp_copy );
          }

        // Add up
        surface_calc( iy, I, surface_los, surface_rmatrix, surface_emission );
      }
      break;


    case 3:   //--- Cloudbox surface -----------------------------------------
      {
        // Pass a local copy of ppath to the cloudbox agenda, to preserve
        // the original ppath when going back to *rte_calc*
        Ppath   ppath_local;
        ppath_init_structure( ppath_local, atmosphere_dim, ppath.np );
        ppath_copy( ppath_local, ppath );

        chk_not_empty( "iy_cloudbox_agenda", iy_cloudbox_agenda );

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


    case 4: // inside the cloudbox
      {
        chk_not_empty( "iy_cloudbox_agenda", iy_cloudbox_agenda );

        iy_cloudbox_agendaExecute( ws, iy, ppath, rte_pos, rte_los, rte_gp_p, 
                     rte_gp_lat, rte_gp_lon, iy_cloudbox_agenda );

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



//! include_cumtrans_in_diy_dq
/*!
    Multiplicates a diy-variable with the transmission to the path's end,
    on the same time as the total transmission is calculated.

    The function handles only a single frequency.

    \param   diy_dq              Input/Output: Corresponds to diy_vmr or
                                               diy_dt.
    \param   trans               Output: The complete transmission of considered
                                        path.
    \param   iv                  Input: Frequency index.
    \param   any_trans_polarised Input: Boolean to indicate if any transmission
                                        matrix in *ppath_transmission* gives
                                        any polarisation effect. If set to 
                                        false, diagonal matrices are assumed.
    \param   ppath_transmission  Input: As the WSV

    \author Patrick Eriksson, 
    \date   2005-06-10
*/
void include_cumtrans_in_diy_dq( 
           Tensor4&   diy_dq,
           Matrix&    trans,
     const Index&     iv,
     const bool&      any_trans_polarised,
     const Tensor3&   ppath_transmissions )
{
  const Index    stokes_dim = diy_dq.ncols();
  const Index    np         = diy_dq.npages();
  const Index    ng         = diy_dq.nbooks();

  if( any_trans_polarised )
    {
      Matrix  mtmp( stokes_dim, stokes_dim );
      Vector  vtmp( stokes_dim );

      // Transmission of 1
      id_mat( trans );

      for( Index ip=0; ip<np-1; ip++ )
        {
          mtmp = trans;
          mult( trans, ppath_transmissions(ip,joker,joker), mtmp );
          for( Index ig=0; ig<ng; ig++ )
            {
              vtmp = diy_dq(ig,ip+1,iv,joker); 
              mult( diy_dq(ig,ip+1,iv,joker), trans, vtmp );
            }
        }
    }
  else
    {
      // Transmission of 1
      id_mat( trans );

      for( Index ip=0; ip<np-1; ip++ )
        {
          for( Index is=0; is<stokes_dim; is++ )
            {
              trans(is,is) *= ppath_transmissions(ip,is,is);
              for( Index ig=0; ig<ng; ig++ )
                { diy_dq(ig,ip+1,iv,is) *= trans(is,is); }
            }
        }
    }
}



//! include_trans_in_diy_dq
/*!
    Multiplicates a diy-variable with a transmission matrix.

    The transmission is included for the propagation path part indicated by
    *ppath_array_index*, and backwards from this ppath part.

    The function handles only a single frequency.

    \param   diy_dq              Input/Output: Corresponds to diy_vmr or
                                               diy_dt.
    \param   iv                  Input: Frequency index.
    \param   pol_trans           Input: Boolean to indicate if *trans* gives
                                        any polarisation effect. If set to 
                                        false, a diagonal *trans* is assumed.
    \param   trans               Input: A transmission, or reflection, matrix.
    \param   ppath_array         Input: As the WSV
    \param   ppath_array_index   Input: As the WSV

    \author Patrick Eriksson, 
    \date   2005-06-10
*/
void include_trans_in_diy_dq( 
            ArrayOfTensor4&   diy_dq,  
      const Index&            iv,
            bool              pol_trans,
      ConstMatrixView         trans,
      const ArrayOfPpath&     ppath_array, 
      const Index&            ppath_array_index )
{
  const Index   stokes_dim = diy_dq[0].ncols();

  if( stokes_dim == 1 )
    { pol_trans = false; }

  Index   pai = ppath_array_index;

  if( pol_trans )
    {
      for( Index ig=0; ig<diy_dq[pai].nbooks(); ig++ )
        {
          for( Index ip=0; ip<ppath_array[pai].np; ip++ )
            {
              const Vector  vtmp = diy_dq[pai](ig,ip,iv,joker);
              mult( diy_dq[pai](ig,ip,iv,joker), trans, vtmp );
            }
        }
    }
  else
    {
      for( Index ig=0; ig<diy_dq[pai].nbooks(); ig++ )
        {
          for( Index ip=0; ip<ppath_array[pai].np; ip++ )
            {
              for( Index is=0; is<stokes_dim; is++ )
                { diy_dq[pai](ig,ip,iv,is) *= trans(is,is); }
            }
        }
    }

  for( Index ia=0; ia<ppath_array[ppath_array_index].next_parts.nelem(); ia++ )
    {
      pai = ppath_array[ppath_array_index].next_parts[ia];
      include_trans_in_diy_dq( diy_dq, iv, pol_trans, trans, ppath_array, pai );
    }
}



//! iy_calc
/*!
  Solves the monochromatic pencil beam RTE.

  The function performs three basic tasks: <br>
    1. Determines the propagation path (by call of ppath_calc). <br>
    2. Determines the radiative background. <br>
    3. Executes *rte_agenda*. <br> <br>

  The start position and LOS shall be put into *los* and *pos* (and
  not *rte_pos* and *rte_los*).

  It is NOT needed to set *iy* to the correct size before calling
  the function.

  \param[in,out] ws Current Workspace
  \param[out] iy
  \param[out] ppath
  \param[out] ppath_array_index
  \param[out] ppath_array
  \param[out] diy_dvmr
  \param[out] diy_dt
  \param[in]  ppath_step_agenda
  \param[in]  rte_agenda
  \param[in]  surface_prop_agenda
  \param[in]  iy_space_agenda
  \param[in]  iy_cloudbox_agenda
  \param[in]  atmosphere_dim
  \param[in]  p_grid
  \param[in]  lat_grid
  \param[in]  lon_grid
  \param[in]  z_field
  \param[in]  t_field
  \param[in]  vmr_field
  \param[in]  r_geoid
  \param[in]  z_surface
  \param[in]  cloudbox_on
  \param[in]  cloudbox_limits
  \param[in]  pos
  \param[in]  los
  \param[in]  f_grid
  \param[in]  stokes_dim
  \param[in]  ppath_array_do
  \param[in]  rte_do_vmr_jacs
  \param[in]  rte_do_t_jacs

  \author Patrick Eriksson 
  \date   2002-09-17
*/
void iy_calc( Workspace&               ws,
              Matrix&                  iy,
              Ppath&                   ppath,
              Index&                   ppath_array_index,
              ArrayOfPpath&            ppath_array,
              ArrayOfTensor4&          diy_dvmr,
              ArrayOfTensor4&          diy_dt,
        const Agenda&                  ppath_step_agenda,
        const Agenda&                  rte_agenda,
        const Agenda&                  iy_space_agenda,
        const Agenda&                  surface_prop_agenda,
        const Agenda&                  iy_cloudbox_agenda,
        const Index&                   atmosphere_dim,
        const Vector&                  p_grid,
        const Vector&                  lat_grid,
        const Vector&                  lon_grid,
        const Tensor3&                 z_field,
        const Tensor3&                 t_field,
        const Tensor4&                 vmr_field,
        const Matrix&                  r_geoid,
        const Matrix&                  z_surface,
        const Index&                   cloudbox_on, 
        const ArrayOfIndex&            cloudbox_limits,
        const Vector&                  pos,
        const Vector&                  los,
        const Vector&                  f_grid,
        const Index&                   stokes_dim,
        const Index&                   ppath_array_do,
        const ArrayOfIndex&            rte_do_vmr_jacs,
        const Index&                   rte_do_t_jacs )
{
  //- Determine propagation path
  const bool  outside_cloudbox = true;
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
              lat_grid, lon_grid, z_field, r_geoid, z_surface,
              cloudbox_on, cloudbox_limits, pos, los, outside_cloudbox );
  //
  const Index   np = ppath.np;


  //- Determine atmospheric fields at each ppath point --------------------
  if( np > 1 )
    {
      // Pressure:
      ppath.p.resize(np);
      Matrix itw_p(np,2);
      interpweights( itw_p, ppath.gp_p );      
      itw2p( ppath.p, p_grid, ppath.gp_p, itw_p );
      
      // Temperature:
      ppath.t.resize(np);
      Matrix   itw_field;
      interp_atmfield_gp2itw( itw_field, atmosphere_dim, p_grid, lat_grid, 
                            lon_grid, ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
      interp_atmfield_by_itw( ppath.t,  atmosphere_dim, p_grid, lat_grid, 
                              lon_grid, t_field, ppath.gp_p, 
                              ppath.gp_lat, ppath.gp_lon, itw_field );

      //  VMR fields:
      const Index ns = vmr_field.nbooks();
      ppath.vmr.resize(ns,np);
      for( Index is=0; is<ns; is++ )
        {
          interp_atmfield_by_itw( ppath.vmr(is, joker), atmosphere_dim,
            p_grid, lat_grid, lon_grid, vmr_field( is, joker, joker,  joker ), 
            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
        }
    }


  // Handle *ppath_array*
  if( ppath_array_do )
    {
      // Include link from previous ppath part to calculated ppath part
      // (if not first part, which is indicated by -1)
      if( ppath_array_index >= 0 )
        { ppath_array[ppath_array_index].next_parts.push_back(
                                                       ppath_array.nelem() ); }
      
      // Add this ppath part to *ppath_array*
      ppath_array_index = ppath_array.nelem();
      ppath_array.push_back( ppath );

      // Extend jacobian quantities
      //
      Index   n = np;
      if( np == 1 )
        { n = 0; }
      //
      if( rte_do_vmr_jacs.nelem() )
        {
          diy_dvmr.push_back( 
            Tensor4(rte_do_vmr_jacs.nelem(),n,f_grid.nelem(),stokes_dim,0.0) );
        }
      if( rte_do_t_jacs )
        {
          diy_dt.push_back( Tensor4(1,n,f_grid.nelem(),stokes_dim,0.0) );
        }
    }
  

  // Determine the radiative background
  //
  iy.resize(f_grid.nelem(),stokes_dim);
  //
  get_radiative_background( ws, iy, ppath, ppath_array_index,
              ppath_array, diy_dvmr, diy_dt, ppath_step_agenda, rte_agenda, 
              surface_prop_agenda, iy_space_agenda, iy_cloudbox_agenda, 
              atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, t_field, 
              vmr_field, r_geoid, z_surface, cloudbox_on, cloudbox_limits, 
              f_grid, stokes_dim, ppath_array_do, rte_do_vmr_jacs, 
              rte_do_t_jacs );


  // If the number of propagation path points is 0 or 1, we are already ready,
  // the observed spectrum equals then the radiative background.
  if( np > 1 )
    {
      rte_agendaExecute( ws, iy, diy_dvmr, diy_dt, ppath,
                         ppath_array, ppath_array_index,
                         rte_do_vmr_jacs, rte_do_t_jacs,
                         stokes_dim, f_grid, rte_agenda );
    }
}



//! iy_calc_no_jacobian
/*!
    Interface to *iy_calc* where jacobian variables can be left out.
*/
void iy_calc_no_jacobian(
              Workspace&      ws,
              Matrix&         iy,
              Ppath&          ppath,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         iy_space_agenda,
        const Agenda&         surface_prop_agenda,
        const Agenda&         iy_cloudbox_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        const Tensor4&        vmr_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         pos,
        const Vector&         los,
        const Vector&         f_grid,
        const Index&          stokes_dim )
{
  ArrayOfIndex               rte_do_vmr_jacs(0);
  ArrayOfTensor4             diy_dvmr(0);
  Index                      rte_do_t_jacs=0;
  ArrayOfTensor4             diy_dt(0);
  Index                      ppath_array_do = 0;
  ArrayOfPpath               ppath_array(0);
  Index                      ppath_array_index=-1;

  iy_calc( ws, iy, ppath, 
           ppath_array_index, ppath_array, diy_dvmr, diy_dt, 
           ppath_step_agenda, rte_agenda, iy_space_agenda, surface_prop_agenda,
           iy_cloudbox_agenda, atmosphere_dim, p_grid, lat_grid, lon_grid, 
           z_field, t_field, vmr_field, r_geoid, z_surface, cloudbox_on, 
           cloudbox_limits, pos, los, f_grid, stokes_dim, 
           ppath_array_do,
           rte_do_vmr_jacs, rte_do_t_jacs );
}



//! rte_step_std
/*!
    Solves monochromatic VRTE for an atmospheric slab with constant 
    conditions.

    The function can be used for clearsky and cloudbox calculations.

    The function is best explained by considering a homogenous layer. That is,
    the physical conditions inside the layer are constant. In reality they
    are not constant, so in practical all coefficients have to be averaged 
    before calling this function. 
    Total extinction and absorption inside the layer are described by
    *ext_mat_av* and *abs_vec_av* respectively,
    the blackbdody radiation of the layer is given by *rte_planck_value*
    and the propagation path length through the layer is *l_step*.

    There is an additional scattering source term in the 
    VRTE, the scattering integral term. For this function a constant
    scattering term is assumed. The radiative transfer step is only a part 
    the iterative solution of the scattering problem, for more 
    information consider AUG. In the clearsky case this variable has to be
    set to 0.

    When calling the function, the vector *stokes_vec* shall contain the
    Stokes vector for the incoming radiation. The function returns this
    vector, then containing the outgoing radiation on the other side of the 
    layer.

    The function performs the calculations differently depending on the
    conditions to improve the speed. There are three cases: <br>
       1. Scalar absorption (stokes_dim = 1). <br>
       2. The matrix ext_mat_gas is diagonal (unpolarised absorption). <br>
       3. The total general case.

    \param   stokes_vec         Input/Output: A Stokes vector.
    \param   trans_mat          Input/Output: Transmission matrix of slab.
    \param   ext_mat_av         Input: Averaged extinction matrix.
    \param   abs_vec_av         Input: Averaged absorption vector.
    \param   sca_vec_av         Input: averaged scattering vector.
    \param   l_step             Input: The length of the RTE step.
    \param   rte_planck_value   Input: Blackbody radiation.

    \author Claudia Emde and Patrick Eriksson, 
    \date   2002-11-22
*/
void rte_step_std(//Output and Input:
              VectorView stokes_vec,
              MatrixView trans_mat,
              //Input
              ConstMatrixView ext_mat_av,
              ConstVectorView abs_vec_av,
              ConstVectorView sca_vec_av,
              const Numeric& l_step,
              const Numeric& rte_planck_value )
{
  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Check inputs:
  assert(is_size(trans_mat, stokes_dim, stokes_dim));
  assert(is_size(ext_mat_av, stokes_dim, stokes_dim));
  assert(is_size(abs_vec_av, stokes_dim));
  assert(is_size(sca_vec_av, stokes_dim));
  assert( rte_planck_value >= 0 );
  assert( l_step >= 0 );
  assert (!is_singular( ext_mat_av ));


  // Check, if only the first component of abs_vec is non-zero:
  bool unpol_abs_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (abs_vec_av[i] != 0)
      unpol_abs_vec = false;

  bool unpol_sca_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (sca_vec_av[i] != 0)
      unpol_sca_vec = false;


  //--- Scalar case: ---------------------------------------------------------
  if( stokes_dim == 1 )
    {
      trans_mat(0,0) = exp(-ext_mat_av(0,0) * l_step);
      stokes_vec[0]  = stokes_vec[0] * trans_mat(0,0) +
                       (abs_vec_av[0] * rte_planck_value + sca_vec_av[0]) /
        ext_mat_av(0,0) * (1 - trans_mat(0,0));
    }


  //--- Vector case: ---------------------------------------------------------

  // We have here two cases, diagonal or non-diagonal ext_mat_gas
  // For diagonal ext_mat_gas, we expect abs_vec_gas to only have a
  // non-zero value in position 1.

  //- Unpolarised
  else if( is_diagonal(ext_mat_av) && unpol_abs_vec && unpol_sca_vec )
    {
      trans_mat      = 0;
      trans_mat(0,0) = exp(-ext_mat_av(0,0) * l_step);

      // Stokes dim 1
      //   assert( ext_mat_av(0,0) == abs_vec_av[0] );
      //   Numeric transm = exp( -l_step * abs_vec_av[0] );
      stokes_vec[0] = stokes_vec[0] * trans_mat(0,0) +
                      (abs_vec_av[0] * rte_planck_value + sca_vec_av[0]) /
                      ext_mat_av(0,0) * (1 - trans_mat(0,0));

      // Stokes dims > 1
      for( Index i=1; i<stokes_dim; i++ )
        {
          //      assert( abs_vec_av[i] == 0.);
          trans_mat(i,i) = trans_mat(0,0);
          stokes_vec[i]  = stokes_vec[i] * trans_mat(i,i) +
                       sca_vec_av[i] / ext_mat_av(i,i)  * (1 - trans_mat(i,i));
        }
    }


  //- General case
  else
    {
      //Initialize internal variables:

      // Matrix LU used for LU decompostion and as dummy variable:
      Matrix LU(stokes_dim, stokes_dim);
      ArrayOfIndex indx(stokes_dim); // index for pivoting information
      Vector b(stokes_dim); // dummy variable
      Vector x(stokes_dim); // solution vector for K^(-1)*b
      Matrix I(stokes_dim, stokes_dim);

      Vector B_abs_vec(stokes_dim);
      B_abs_vec = abs_vec_av;
      B_abs_vec *= rte_planck_value;

      for (Index i=0; i<stokes_dim; i++)
        b[i] = B_abs_vec[i] + sca_vec_av[i];  // b = abs_vec * B + sca_vec

      // solve K^(-1)*b = x
      ludcmp(LU, indx, ext_mat_av);
      lubacksub(x, LU, b, indx);

      Matrix ext_mat_ds(stokes_dim, stokes_dim);
      ext_mat_ds = ext_mat_av;
      ext_mat_ds *= -l_step; // ext_mat_ds = -ext_mat*ds

      Index q = 10;  // index for the precision of the matrix exp function
      //Matrix exp_ext_mat(stokes_dim, stokes_dim);
      //matrix_exp(exp_ext_mat, ext_mat_ds, q);
      matrix_exp( trans_mat, ext_mat_ds, q);

      Vector term1(stokes_dim);
      Vector term2(stokes_dim);

      id_mat(I);
      for(Index i=0; i<stokes_dim; i++)
        {
          for(Index j=0; j<stokes_dim; j++)
            LU(i,j) = I(i,j) - trans_mat(i,j); // take LU as dummy variable
          // LU(i,j) = I(i,j) - exp_ext_mat(i,j); // take LU as dummy variable
        }
      mult(term2, LU, x); // term2: second term of the solution of the RTE with
                          //fixed scattered field

      // term1: first term of solution of the RTE with fixed scattered field
      //mult(term1, exp_ext_mat, stokes_vec);
      mult( term1, trans_mat, stokes_vec );

      for (Index i=0; i<stokes_dim; i++)
        stokes_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector
    }
}



//! rte_step_std_clearsky
/*!
    Solves monochromatic clear sky VRTE for an atmospheric slab with constant
    temperature and absorption. 

    \param   stokes_vec         Input/Output: A Stokes vector.
    \param   trans_mat          Input/Output: Transmission matrix of slab.
    \param   absorption         Input: Gas absorption for each species.
    \param   l_step             Input: The length of the RTE step.
    \param   emission           Input: Blackbody radiation.

    \author Patrick Eriksson, 
    \date   2005-06-01
*/
void rte_step_std_clearsky(
              //Output and Input:
                    VectorView   stokes_vec,
                    MatrixView   trans_mat,
              //Input
              const Vector&      absorption,
              const Numeric&     l_step,
              const Numeric&     emission )
{
  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Check inputs:
  assert( is_size( trans_mat, stokes_dim, stokes_dim ) ); 
  assert( emission >= 0 );
  assert( l_step >= 0 );

  // Temporary solution
  const bool abs_polarised = false;


  //--- Polarised absorption -------------------------------------------------
  if( abs_polarised )
    { 
      // We do not allow scalar case here
      assert( stokes_dim > 1 );

      throw runtime_error( "Polarised absorption not yet handled.");
    }


  //--- Unpolarised absorption -----------------------------------------------
  else
    {
      // Init transmission matrix to zero
      trans_mat      = 0;

      // Stokes dim 1
      trans_mat(0,0) = exp( -absorption.sum() * l_step );
      stokes_vec[0]  = stokes_vec[0] * trans_mat(0,0) + 
                                             emission * ( 1 - trans_mat(0,0) );

      // Higher Stokes dims
      for( Index i=1; i<stokes_dim; i++ )
        {
          trans_mat(i,i) = trans_mat(0,0); 
          stokes_vec[i]  *= trans_mat(i,i);
        }
      
    }
}



//! rte_std
/*! 
  Core function for the different versions of WSM RteStd. 

  See the the online help (arts -d RteStd)

  \param[in,out] ws Current Workspace
  \param[out] iy
  \param[out] ppath_transmissions
  \param[out] diy_dvmr
  \param[out] diy_dt
  \param[in]  ppath
  \param[in]  ppath_array
  \param[in]  ppath_array_index
  \param[in]  f_grid
  \param[in]  stokes_dim
  \param[in]  emission_agenda
  \param[in]  abs_scalar_gas_agenda
  \param[in]  rte_do_vmr_jacs
  \param[in]  rte_do_t_jacs
  \param[in]  do_transmissions   Boolean to fill *ppath_transmissions* or not.

  \author Patrick Eriksson
  \date   2005-05-19
*/
void rte_std(Workspace&               ws,
             Matrix&                  iy,
             Tensor4&                 ppath_transmissions,
             ArrayOfTensor4&          diy_dvmr,
             ArrayOfTensor4&          diy_dt,
       const Ppath&                   ppath,
       const ArrayOfPpath&            ppath_array, 
       const Index&                   ppath_array_index,
       const Vector&                  f_grid,
       const Index&                   stokes_dim,
       const Agenda&                  emission_agenda,
       const Agenda&                  abs_scalar_gas_agenda,
       const ArrayOfIndex&            rte_do_vmr_jacs,
       const Index&                   rte_do_t_jacs,
       const bool&                    do_transmissions )
{
  // Relevant checks are assumed to be done in RteCalc

  // Some sizes
  const Index   nf = f_grid.nelem();
  const Index   np = ppath.np;
  const Index   ns = ppath.vmr.nrows();        // Number of species
  const Index   ng = rte_do_vmr_jacs.nelem();  // Number of jacobian species
        Vector  rte_vmr_list(ns);

  // Log of pressure
  Vector   logp_ppath(np);
  transform( logp_ppath, log, ppath.p  );

  // If f_index < 0, scalar gas absorption is calculated for 
  // all frequencies in f_grid.
  Index f_index = -1;

  // Transmission variables
  Matrix    trans(stokes_dim,stokes_dim);
  bool      save_transmissions = false;
  bool      any_abs_polarised = false;
  //
  if( do_transmissions  ||  ng  || rte_do_t_jacs ) 
    {
      save_transmissions = true;
      ppath_transmissions.resize(np-1,nf,stokes_dim,stokes_dim); 
    }

  // Local declaration of corresponding WSVs
  //
  Vector emission;
  Matrix abs_scalar_gas;

  // Loop the propagation path steps
  //
  // The number of path steps is np-1.
  // The path points are stored in such way that index 0 corresponds to
  // the point closest to the sensor.
  //
  for( Index ip=np-1; ip>0; ip-- )
    {
      // Calculate mean of atmospheric parameters
      Numeric rte_pressure = exp( 0.5 * ( logp_ppath[ip] + logp_ppath[ip-1] ) );
      Numeric rte_temperature = 0.5 * ( ppath.t[ip] + ppath.t[ip-1] );
      for( Index is=0; is<ns; is++ )
        { rte_vmr_list[is] = 0.5 * ( ppath.vmr(is,ip) + ppath.vmr(is, ip-1) );}
      
      // Call agendas for RT properties
      emission_agendaExecute( ws, emission, rte_temperature, emission_agenda );
      abs_scalar_gas_agendaExecute( ws, abs_scalar_gas, f_index,
          rte_pressure, rte_temperature, rte_vmr_list, abs_scalar_gas_agenda );

      // Calculate "disturbed" emission and absorption if temperature 
      // jacobians shall be calculated
      //
      Numeric dt = 0.1;
      Vector emission2(0);
      Matrix abs_scalar_gas2(0,0);
      //
      if( rte_do_t_jacs )
        {
          emission_agendaExecute( ws, emission2, rte_temperature+dt, 
                                                              emission_agenda );
          abs_scalar_gas_agendaExecute( ws, abs_scalar_gas2, f_index,
                                        rte_pressure, rte_temperature+dt, 
                                        rte_vmr_list, abs_scalar_gas_agenda );
        }

      // Polarised absorption?
      // What check to do here?
      bool abs_polarised = false;
      if( abs_polarised )
        { any_abs_polarised = true; }


      //--- Loop frequencies
      for( Index iv=0; iv<nf; iv++ )
        {
          //--- Jacobians -----------------------------------------------------
          if( ng  || rte_do_t_jacs )  
            {
              if( abs_polarised )
                {}  // To be coded when format for polarised absorption is set
              else
                {
                  const Numeric tau = abs_scalar_gas(iv,joker).sum() * 
                                                             ppath.l_step[ip-1];
                  const Numeric tr = exp( -tau );

                  // Gas species
                  if( ng )
                    {
                      const Numeric cf  = 0.5 * ppath.l_step[ip-1] * tr;
                      const Numeric cf0 = cf * ( emission[iv] - iy(iv,0) );

                      for( Index ig=0; ig<ng; ig++ )
                        {
                          const Index   is = rte_do_vmr_jacs[ig];
                          const Numeric k  = abs_scalar_gas(iv,is) / 
                                                               rte_vmr_list[is];
                          Numeric w  = cf0 * k;
                          //         
                          diy_dvmr[ppath_array_index](ig,ip,iv,0)   += w;
                          diy_dvmr[ppath_array_index](ig,ip-1,iv,0) += w;
                      
                          // Higher Stokes components
                          for( Index s=1; s<stokes_dim; s++ )
                            {
                              w = -cf * k * iy(iv,s);
                              //
                              diy_dvmr[ppath_array_index](ig,ip,iv,s)   += w;
                              diy_dvmr[ppath_array_index](ig,ip-1,iv,s) += w;
                            }
                        }
                    }

                  //- Temperature
                  if( rte_do_t_jacs )
                    {
                      const Numeric dtaudT = (abs_scalar_gas2(iv,joker).sum() *
                                               ppath.l_step[ip-1] - tau ) / 
                                                                   ( 2.0 * dt );
                            Numeric w = ( 1.0 -tr ) *
                               ( emission2[iv] - emission[iv] ) / ( 2.0 * dt ) +
                                      tr * ( emission[iv] - iy(iv,0) ) * dtaudT;
                      //         
                      diy_dt[ppath_array_index](0,ip,iv,0)   += w;
                      diy_dt[ppath_array_index](0,ip-1,iv,0) += w;

                      // Higher Stokes components
                      for( Index s=1; s<stokes_dim; s++ )
                        {
                          w = -tr * iy(iv,s) * dtaudT;
                          //
                          diy_dvmr[ppath_array_index](0,ip,iv,s)   += w;
                          diy_dvmr[ppath_array_index](0,ip-1,iv,s) += w;
                        }
                    }
                }
            }
          //--- End jacobians -------------------------------------------------


          // Perform the RTE step.
          rte_step_std_clearsky( iy(iv,joker), trans, abs_scalar_gas(iv,joker),
                                            ppath.l_step[ip-1], emission[iv] );

          if( save_transmissions )
            { ppath_transmissions(ip-1,iv,joker,joker) = trans; }
        }
    }


  //--- Postprocessing of Jacobians
  if( ng  ||  rte_do_t_jacs )
    {
      for( Index iv=0; iv<nf; iv++ )
        {      
          if( ng )
            include_cumtrans_in_diy_dq( diy_dvmr[ppath_array_index], trans, 
                                        iv, any_abs_polarised, 
                                    ppath_transmissions(joker,iv,joker,joker) );
          if( rte_do_t_jacs )
            include_cumtrans_in_diy_dq( diy_dt[ppath_array_index], trans, 
                                        iv, any_abs_polarised, 
                                    ppath_transmissions(joker,iv,joker,joker) );
          
          for( Index ia=0; 
                   ia<ppath_array[ppath_array_index].next_parts.nelem(); ia++ )
            {
              const Index  pai = ppath_array[ppath_array_index].next_parts[ia];
              if( ng )
                include_trans_in_diy_dq( diy_dvmr, iv, any_abs_polarised,
                                         trans, ppath_array, pai );
              if( rte_do_t_jacs )
                include_trans_in_diy_dq( diy_dt, iv, any_abs_polarised,
                                         trans, ppath_array, pai );
            }
        }
    }
}



//! rtecalc_check_input
/*! 
   Common sub-function for RteCalc functions to check consistency of input.

   \param[out] nf        Number of frequencies in f_grid.
   \param[out] nmblock   Number of measurement blocks.
   \param[out] nza       Number of zenith angles / measurement block.
   \param[out] naa       Number of azimuth angles / measurement block.
   \param[out] nblock    Length of for one measurement block.
   \param[in]  atmosphere_dim
   \param[in]  p_grid
   \param[in]  lat_grid
   \param[in]  lon_grid
   \param[in]  z_field
   \param[in]  t_field
   \param[in]  r_geoid
   \param[in]  z_surface
   \param[in]  cloudbox_on
   \param[in]  cloudbox_limits
   \param[in]  sensor_response
   \param[in]  sensor_pos
   \param[in]  sensor_los
   \param[in]  f_grid
   \param[in]  stokes_dim
   \param[in]  antenna_dim
   \param[in]  mblock_za_grid
   \param[in]  mblock_aa_grid
   \param[in]  y_unit
   \param[in]  jacobian_unit

   \author Patrick Eriksson
   \date   2007-09-19
*/
void rtecalc_check_input(
         Index&                      nf,
         Index&                      nmblock,
         Index&                      nza,
         Index&                      naa,
         Index&                      nblock,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
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
   const String&                     y_unit,
   const String&                     jacobian_unit )
{
  // Some sizes
  nf      = f_grid.nelem();
  nmblock = sensor_pos.nrows();
  nza     = mblock_za_grid.nelem();

  // Number of azimuthal direction for pencil beam calculations
  naa = mblock_aa_grid.nelem();
  if( antenna_dim == 1 )
    { naa = 1; }

  // Number of elements of *y* for one mblock
  nblock = sensor_response.nrows();
  
  // Stokes
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );

  // Basic checks of atmospheric, geoid and surface variables
  //  
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "z_field", z_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  // vmr_field excluded, could maybe be empty 
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

  // *y_unit*
  //
  if( !( y_unit=="1"  ||  y_unit=="RJBT"  ||  y_unit=="PlanckBT" ) )
    {
      ostringstream os;
      os << "Allowed options for *y_unit* are: ""1"", ""RJBT"", and "
         << """PlanckBT"".";
      throw runtime_error( os.str() );
    }

  // *jacobian_unit*
  //
  if( !( jacobian_unit=="1"  ||  jacobian_unit=="RJBT"  ||  
         jacobian_unit=="-" ) )
    {
      ostringstream os;
      os << "Allowed options for *jacobian_unit* are: ""1"", ""RJBT"", and "
         << """-"".";
      throw runtime_error( os.str() );
    }
}


//! surface_calc
/*!
    Weights together downwelling radiation and surface emission.

    *iy* must have correct size when function is called.

    \param   iy                 In/Out: Radiation matrix, amtching 
                                        the WSV with the same name.
    \param   I                  Input: Downwelling radiation, with dimensions
                                       matching: 
                                       (surface_los, f_grid, sokes_dim)
    \param   surface_los        Input: As the WSV with the same name.
    \param   surface_rmatrix    Input: As the WSV with the same name.
    \param   surface_emission   Input: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2005-04-07
*/
void surface_calc(
              Matrix&         iy,
        const Tensor3&        I,
        const Matrix&         surface_los,
        const Tensor4&        surface_rmatrix,
        const Matrix&         surface_emission )
{
  // Some sizes
  const Index   nf         = I.nrows();
  const Index   stokes_dim = I.ncols();
  const Index   nlos       = surface_los.nrows();

  iy = surface_emission;
  
  /*
  cout << surface_emission << "\n";
  cout << iy << "\n";
  cout << I << "\n";
  */

  // Loop *surface_los*-es. If no such LOS, we are ready.
  if( nlos > 0 )
    {
      for( Index ilos=0; ilos<nlos; ilos++ )
        {
          Vector rtmp(stokes_dim);  // Reflected Stokes vector for 1 frequency

          for( Index iv=0; iv<nf; iv++ )
            {
          mult( rtmp, surface_rmatrix(ilos,iv,joker,joker), I(ilos,iv,joker) );
          iy(iv,joker) += rtmp;
            }
        }
    }
}






