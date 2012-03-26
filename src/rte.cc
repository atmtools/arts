/* Copyright (C) 2002-2008
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
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
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
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
#include "geodetic.h"
#include "logic.h"
#include "math_funcs.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "lin_alg.h"

extern const Numeric BOLTZMAN_CONST;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



/*===========================================================================
  === The functions in alphabetical order
  ===========================================================================*/

//! apply_y_unit
/*!
    Performs conversion from radiance to other units.

    Use *apply_y_unit_jac* for conversion of jacobian data.

    \param   iy       In/Out: Tensor3 with data to be converted, where 
                      column dimension corresponds to Stokes dimensionality
                      and row dimension corresponds to frequency.
    \param   y_unit   As the WSV.
    \param   f_grid   As the WSV.
    \param   i_pol    Polarisation indexes. See documentation of y_pol.

    \author Patrick Eriksson 
    \date   2010-04-07
*/
void apply_y_unit( 
            MatrixView   iy, 
         const String&   y_unit, 
       ConstVectorView   f_grid,
   const ArrayOfIndex&   i_pol )
{
  // The code is largely identical between the two apply_y_unit functions.
  // If any change here, remember to update the other function.

  const Index nf = iy.nrows();
  const Index ns = iy.ncols();

  assert( f_grid.nelem() == nf );
  assert( i_pol.nelem() == ns );

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = invrayjean( 1, f_grid[iv] );
          for( Index is=0; is<ns; is++ )
            {
              if( i_pol[is] < 5 )           // Stokes components
                { iy(iv,is) *= scfac; }
              else                          // Measuement single pols
                { iy(iv,is) *= 2*scfac; }
            }
        }
    }

  else if( y_unit == "PlanckBT" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          for( Index is=ns-1; is>=0; is-- ) // Order must here be reversed
            {
              if( i_pol[is] == 1 )
                { iy(iv,is) = invplanck( iy(iv,is), f_grid[iv] ); }
              else if( i_pol[is] < 5 )
                { 
                  assert( i_pol[0] == 1 );
                  iy(iv,is) = 
                    invplanck( 0.5*(iy(iv,0)+iy(iv,is)), f_grid[iv] ) -
                    invplanck( 0.5*(iy(iv,0)-iy(iv,is)), f_grid[iv] );
                }
              else
                { iy(iv,is) = invplanck( 2*iy(iv,is), f_grid[iv] ); }
            }
        }
    }
  
  else if ( y_unit == "W/(m^2 m sr)" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = f_grid[iv] * ( f_grid[iv] / SPEED_OF_LIGHT );
          for( Index is=0; is<ns; is++ )
            { iy(iv,is) *= scfac; }
        }
    }
  
  else if ( y_unit == "W/(m^2 m-1 sr)" )
    {
      iy *= SPEED_OF_LIGHT;
    }

  else
    {
      ostringstream os;
      os << "Unknown option: y_unit = \"" << y_unit << "\"\n" 
         << "Recognised choices are: \"1\", \"RJBT\", \"PlanckBT\""
         << "\"W/(m^2 m sr)\" and \"W/(m^2 m-1 sr)\""; 
      
      throw runtime_error( os.str() );      
    }
}



//! apply_y_unit2
/*!
    Largely as *apply_y_unit* but operates on jacobian data.

    The associated spectrum data *iy* must be in radiance. That is, the
    spectrum can only be converted to Tb after the jacobian data. 

    *iy* must be a single spectrum, and is accordingly here a matrix (and not a
    *Tensor3 as for apply_y_unit).

    \param   J        In/Out: Tensor3 with data to be converted, where 
                      column dimension corresponds to Stokes dimensionality
                      and row dimension corresponds to frequency.
    \param   iy       Associated radiance data.
    \param   y_unit   As the WSV.
    \param   f_grid   As the WSV.

    \author Patrick Eriksson 
    \date   2010-04-10
*/
void apply_y_unit2( 
   Tensor3View           J,
   ConstMatrixView       iy, 
   const String&         y_unit, 
   ConstVectorView       f_grid,
   const ArrayOfIndex&   i_pol )
{
  // The code is largely identical between the two apply_y_unit functions.
  // If any change here, remember to update the other function.

  const Index nf = iy.nrows();
  const Index ns = iy.ncols();
  const Index np = J.npages();

  assert( J.nrows() == nf );
  assert( J.ncols() == ns );
  assert( f_grid.nelem() == nf );
  assert( i_pol.nelem() == ns );

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = invrayjean( 1, f_grid[iv] );
          for( Index is=0; is<ns; is++ )
            {
              if( i_pol[is] < 5 )           // Stokes componenets
                {
                  for( Index ip=0; ip<np; ip++ )
                    { J(ip,iv,is) *= scfac; }
                }
              else                          // Measuement single pols
                {
                  for( Index ip=0; ip<np; ip++ )
                    { J(ip,iv,is) *= 2*scfac; }
                }
            }
        }
    }

  else if( y_unit == "PlanckBT" )
    {
      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          for( Index is=ns-1; is>=0; is-- )
            {
              Numeric scfac = 1;
              if( i_pol[is] == 1 )
                { scfac = dinvplanckdI( iy(iv,is), f_grid[iv] ); }
              else if( i_pol[is] < 5 )
                {
                  assert( i_pol[0] == 1 );
                  scfac = 
                    dinvplanckdI( 0.5*(iy(iv,0)+iy(iv,is)), f_grid[iv] ) +
                    dinvplanckdI( 0.5*(iy(iv,0)-iy(iv,is)), f_grid[iv] );
                }
              else
                { scfac = dinvplanckdI( 2*iy(iv,is), f_grid[iv] ); }
              //
              for( Index ip=0; ip<np; ip++ )
                { J(ip,iv,is) *= scfac; }
            }
        }
    }

  else if ( y_unit == "W/(m^2 m sr)" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = f_grid[iv] * ( f_grid[iv] / SPEED_OF_LIGHT );
          for( Index ip=0; ip<np; ip++ )
            {
              for( Index is=0; is<ns; is++ )
                { J(ip,iv,is) *= scfac; }
            }
        }
    }
  
  else if ( y_unit == "W/(m^2 m-1 sr)" )
    {
      J *= SPEED_OF_LIGHT;
    }
  
  else
    {
      ostringstream os;
      os << "Unknown option: y_unit = \"" << y_unit << "\"\n" 
         << "Recognised choices are: \"1\", \"RJBT\", \"PlanckBT\""
         << "\"W/(m^2 m sr)\" and \"W/(m^2 m-1 sr)\""; 
      
      throw runtime_error( os.str() );      
    }
  
}



//! ext2trans
/*!
    Converts an extinction matrix to a transmission matrix

    The function performs the calculations differently depending on the
    conditions to improve the speed. There are three cases: <br>
       1. Scalar absorption. <br>
       2. The matrix ext_mat_av is diagonal. <br>
       3. The total general case.

    \param   trans_mat          Input/Output: Transmission matrix of slab.
    \param   ext_mat            Input: Averaged extinction matrix.
    \param   lstep             Input: The length of the RTE step.

    \author Claudia Emde and Patrick Eriksson, 
    \date   2010-10-15
*/
void ext2trans(//Output and Input:
              MatrixView trans_mat,
              //Input
              ConstMatrixView ext_mat_av,
              const Numeric& lstep )
{
  //Stokes dimension:
  Index stokes_dim = ext_mat_av.nrows();

  //Check inputs:
  assert( is_size(trans_mat, stokes_dim, stokes_dim) );
  assert( is_size(ext_mat_av, stokes_dim, stokes_dim) );
  assert( lstep >= 0 );
  assert( !is_singular( ext_mat_av ) );

  // Any changes here should also be implemented in rte_step_std.

  //--- Scalar case: ---------------------------------------------------------
  if( stokes_dim == 1 )
    {
      trans_mat(0,0) = exp(-ext_mat_av(0,0) * lstep);
    }


  //--- Vector case: ---------------------------------------------------------

  //- Unpolarised
  else if( is_diagonal(ext_mat_av) )
    {
      const Numeric tv = exp(-ext_mat_av(0,0) * lstep);

      trans_mat = 0;

      for( Index i=0; i<stokes_dim; i++ )
        {
          trans_mat(i,i)  = tv;
        }
    }


  //- General case
  else
    {
      Matrix ext_mat_ds = ext_mat_av;
      ext_mat_ds *= -lstep; 

      Index q = 10;  // index for the precision of the matrix exp function

      matrix_exp( trans_mat, ext_mat_ds, q );
    }
}



//! get_iy_of_background
/*!
    Determines iy of the "background" of a propgation path.

    The task is to determine *iy* and related variables for the
    background, or to continue the raditiave calculations
    "backwards". The details here depends on the method selected for
    the agendas.

    Each background is handled by an agenda. Several of these agandes
    can involve recursive calls of *iy_clearsky_agenda*. It is also
    allowed to input *iy_clearsky_basic_agenda* instead of
    *iy_clearsky_agenda*.

    \param   ws                    Out: The workspace
    \param   iy                    Out: As the WSV.
    \param   iy_error              Out: As the WSV.
    \param   iy_error_type         Out: As the WSV.
    \param   iy_aux                Out: As the WSV.
    \param   diy_aux               Out: As the WSV.
    \param   iy_transmission       As the WSV.
    \param   jacobian_do           As the WSV.
    \param   ppath                 As the WSV.
    \param   atmosphere_dim        As the WSV.
    \param   p_grid                As the WSV.
    \param   lat_grid              As the WSV.
    \param   lon_grid              As the WSV.
    \param   t_field               As the WSV.
    \param   z_field               As the WSV.
    \param   vmr_field             As the WSV.
    \param   cloudbox_on           As the WSV.
    \param   stokes_dim            As the WSV.
    \param   f_grid                As the WSV.
    \param   iy_clearsky_agenda    As the WSV or iy_clearsky_basic_agenda.
    \param   iy_space_agenda       As the WSV.
    \param   surface_prop_agenda   As the WSV.
    \param   iy_cloudbox_agenda    As the WSV.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
void get_iy_of_background(
        Workspace&        ws,
        Matrix&           iy,
        Matrix&           iy_error,
        Index&            iy_error_type,
        Matrix&           iy_aux,
        ArrayOfTensor3&   diy_dx,
  ConstTensor3View        iy_transmission,
  const Index&            jacobian_do,
  const Ppath&            ppath,
  const Index&            atmosphere_dim,
  ConstVectorView         lat_grid,
  ConstVectorView         lon_grid,
  ConstTensor3View        t_field,
  ConstTensor3View        z_field,
  ConstTensor4View        vmr_field,
  const Vector&           refellipsoid,
  const Matrix&           z_surface,
  const Index&            cloudbox_on,
  const Index&            stokes_dim,
  ConstVectorView         f_grid,
  const Agenda&           iy_clearsky_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           surface_prop_agenda,
  const Agenda&           iy_cloudbox_agenda,
  const Verbosity&        verbosity)
{
  CREATE_OUT3
  
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
        //
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
                Vector los = surface_los(ilos,joker);

                // Include surface reflection matrix in *iy_transmission*
                // If iy_transmission is empty, this is interpreted as the
                // variable is not needed.
                //
                Tensor3 iy_trans_new;
                //
                if( iy_transmission.npages() )
                  {
                    iy_transmission_mult(  iy_trans_new, iy_transmission, 
                                     surface_rmatrix(ilos,joker,joker,joker) );
                  }

                // Calculate angular tilt of the surface
                // (sign(los[0]) to handle negative za for 2D)
                //
                Numeric atilt = 0.0;
                //
                if( atmosphere_dim == 2 )
                  {
                    const double rslope = plevel_slope_2d( lat_grid, 
                        refellipsoid, z_surface(joker,0), rte_gp_lat, los[0] );
                    Vector itw(2); interpweights( itw, rte_gp_lat );
                    const Numeric rv_surface = 
                              refell2r( refellipsoid, rte_pos[1] ) +
                              interp( itw, z_surface(joker,0), rte_gp_lat );
                    atilt = plevel_angletilt( rv_surface, rslope);
                  }
                else if ( atmosphere_dim == 3 )
                  {
                    const double rslope = plevel_slope_3d( lat_grid, lon_grid, 
                     refellipsoid, z_surface, rte_gp_lat, rte_gp_lon, los[1] );
                    Vector itw(4); interpweights( itw, rte_gp_lat, rte_gp_lon );
                    const Numeric rv_surface = 
                              refell2r( refellipsoid, rte_pos[1] ) +
                              interp( itw, z_surface, rte_gp_lat, rte_gp_lon );
                    atilt = plevel_angletilt( rv_surface, rslope);
                  }

                const Numeric zamax = 90 - sign(los[0])*atilt;

                // I considered to add a check that surface_los is above the
                // horizon, but that would force e.g. surface_specular_los to
                // actually calculate the surface tilt, which causes
                // unnecessary calculation overhead. The angles are now moved
                // to be just above the horizon, which should be acceptable.

                // Make sure that we have some margin to the "horizon"
                // (otherwise numerical problems can create an infinite loop)
                if( atmosphere_dim == 2  && los[0]<0 ) //2D with za<0
                  { los[0] = max( -zamax+0.1, los[0] ); }
                else
                  { los[0] = min( zamax-0.1, los[0] ); }

                // Calculate downwelling radiation for LOS ilos 
                //
                // The variable iy_clearsky_agenda can in fact be 
                // iy_clearsky_BASIC_agenda
                //
                if( iy_clearsky_agenda.name() == "iy_clearsky_basic_agenda" )
                  {
                    iy_clearsky_basic_agendaExecute( ws, iy, rte_pos, los,
                                              cloudbox_on, iy_clearsky_agenda);
                  }
                else
                  {
                    iy_clearsky_agendaExecute( ws, iy, iy_error, iy_error_type,
                                  iy_aux, diy_dx, 0, iy_trans_new, rte_pos, 
                                  los, cloudbox_on, jacobian_do, t_field,
                                  z_field, vmr_field, -1, iy_clearsky_agenda );
                  }

                I(ilos,joker,joker) = iy;
              }
          }

        // Add up
        surface_calc( iy, I, surface_los, surface_rmatrix, surface_emission );
      }
      break;


    case 3:   //--- Cloudbox boundary or interior ------------------------------
    case 4:
      {
        iy_cloudbox_agendaExecute( ws, iy, iy_error, iy_error_type,
                                   iy_aux, diy_dx, iy_transmission,
                                   rte_pos, rte_los, rte_gp_p, rte_gp_lat, 
                                   rte_gp_lon, iy_cloudbox_agenda );

        if( iy.nrows() != nf  ||  iy.ncols() != stokes_dim )
          {
            CREATE_OUT1
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



//! get_ppath_atmvars
/*!
    Determines pressure, temperature, VMR and winds for each propgataion path
    point.

    The output variables are sized inside the function. For VMR the
    dimensions are [ species, propagation path point ].

    \param   ppath_p           Out: Pressure for each ppath point.
    \param   ppath_t           Out: Temperature for each ppath point.
    \param   ppath_vmr         Out: VMR values for each ppath point.
    \param   ppath_wind_u      Out: U-wind for each ppath point.
    \param   ppath_wind_v      Out: V-wind for each ppath point.
    \param   ppath_wind_w      Out: W-wind for each ppath point.
    \param   ppath             As the WSV.
    \param   atmosphere_dim    As the WSV.
    \param   p_grid            As the WSV.
    \param   lat_grid          As the WSV.
    \param   lon_grid          As the WSV.
    \param   t_field           As the WSV.
    \param   vmr_field         As the WSV.
    \param   wind_u_field      As the WSV.
    \param   wind_v_field      As the WSV.
    \param   wind_w_field      As the WSV.

    \author Patrick Eriksson 
    \date   2009-10-05
*/
void get_ppath_atmvars( 
        Vector&      ppath_p, 
        Vector&      ppath_t, 
        Matrix&      ppath_vmr, 
        Vector&      ppath_wind_u, 
        Vector&      ppath_wind_v, 
        Vector&      ppath_wind_w, 
  const Ppath&       ppath,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstTensor3View   t_field,
  ConstTensor4View   vmr_field,
  ConstTensor3View   wind_u_field,
  ConstTensor3View   wind_v_field,
  ConstTensor3View   wind_w_field )
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
  interp_atmfield_gp2itw( itw_field, atmosphere_dim, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
  interp_atmfield_by_itw( ppath_t,  atmosphere_dim, t_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );

  // VMR fields:
  const Index ns = vmr_field.nbooks();
  ppath_vmr.resize(ns,np);
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( ppath_vmr(is, joker), atmosphere_dim,
                              vmr_field( is, joker, joker,  joker ), 
                              ppath.gp_p, ppath.gp_lat, ppath.gp_lon, 
                              itw_field );
    }

  // Winds:
  //
  ppath_wind_w.resize(np);
  if( wind_w_field.npages() > 0 ) 
    { 
      interp_atmfield_by_itw( ppath_wind_w,  atmosphere_dim, wind_w_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  else
    { ppath_wind_w = 0; }
  //
  ppath_wind_v.resize(np);
  if( wind_v_field.npages() > 0 ) 
    { 
      interp_atmfield_by_itw( ppath_wind_v,  atmosphere_dim, wind_v_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  else
    { ppath_wind_v = 0; }
  //
  ppath_wind_u.resize(np);
  if( atmosphere_dim > 2 )
    {
      if( wind_u_field.npages() > 0 ) 
        { 
          interp_atmfield_by_itw( ppath_wind_u,  atmosphere_dim, wind_u_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
        }
      else
        { ppath_wind_u = 0; }
    }
  else
    { ppath_wind_u = 0; }
}



//! get_ppath_pnd
/*!
    Determines particle number densities for each propgataion path
    point.

    The output variable is sized inside the function as
    [ particle type, propagation path point ].

    \param   ppath_pnd         Out: PND for each ppath point.
    \param   ppath             As the WSV.
    \param   atmosphere_dim    As the WSV.
    \param   pnd_field         As the WSV.

    \author Patrick Eriksson 
    \date   2011-07-14
*/
void get_ppath_pnd( 
        Matrix&         ppath_pnd, 
  const Ppath&          ppath,
  const Index&          atmosphere_dim,
  const ArrayOfIndex&   cloudbox_limits,
  ConstTensor4View      pnd_field )
{
  ppath_pnd.resize( pnd_field.nbooks(), ppath.np );

  Matrix   itw_field;
  ArrayOfGridPos   gpc_p, gpc_lat, gpc_lon;
  //
  interp_cloudfield_gp2itw( itw_field, gpc_p, gpc_lat, gpc_lon, 
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, 
                            atmosphere_dim, cloudbox_limits );
  for( Index ip=0; ip<pnd_field.nbooks(); ip++ )
    {
      interp_atmfield_by_itw( ppath_pnd(ip,joker), atmosphere_dim,
                              pnd_field(ip,joker,joker,joker), 
                              gpc_p, gpc_lat, gpc_lon, itw_field );
    }
}


//! get_ppath_rtvars
/*!
    Determines variables for each step of "standard" RT integration

    The output variables are sized inside the function. The dimension order is 
       [ frequency, absorption species, ppath point ]

    \param   ws                Out: The workspace
    \param   ppath_abs_sclar   Out: Absorption coefficients at each ppath point 
    \param   ppath_tau         Out: Optical thickness of each ppath step 
    \param   total_tau         Out: Total optical thickness of path
    \param   ppath_emission    Out: Emission source term at each ppath point 
    \param   abs_scalar_gas_agenda As the WSV.    
    \param   emission_agenda   As the WSV.    
    \param   ppath_p           Pressure for each ppath point.
    \param   ppath_t           Temperature for each ppath point.
    \param   ppath_vmr         VMR values for each ppath point.
    \param   ppath_wind_u      U-wind for each ppath point.
    \param   ppath_wind_v      V-wind for each ppath point.
    \param   ppath_wind_w      W-wind for each ppath point.
    \param   f_grid            As the WSV.    
    \param   atmosphere_dim    As the WSV.    
    \param   emission_do       Flag for calculation of emission. Should be
                               set to 0 for pure transmission calculations.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void get_ppath_rtvars( 
        Workspace&   ws,
        Tensor3&     ppath_abs_scalar, 
        Matrix&      ppath_tau,
        Vector&      total_tau,
        Matrix&      ppath_emission, 
  const Agenda&      abs_scalar_gas_agenda,
  const Agenda&      emission_agenda,
  const Ppath&       ppath,
  ConstVectorView    ppath_p, 
  ConstVectorView    ppath_t, 
  ConstMatrixView    ppath_vmr, 
  ConstVectorView    ppath_wind_u, 
  ConstVectorView    ppath_wind_v, 
  ConstVectorView    ppath_wind_w, 
  ConstVectorView    f_grid, 
  const Index&       atmosphere_dim,
  const Index&       emission_do )
{
  // Sizes
  const Index   np   = ppath.np;
  const Index   nabs = ppath_vmr.nrows();
  const Index   nf   = f_grid.nelem();

  // Init variables
  ppath_abs_scalar.resize( nf, nabs, np );
  ppath_tau.resize( nf, np-1 );
  total_tau.resize( nf );
  total_tau = 0;
  if( emission_do )
    { ppath_emission.resize( nf, np ); }
  else
    { ppath_emission.resize( 0, 0 ); }

  // Mean of extreme frequencies
  const Numeric f0 = (f_grid[0]+f_grid[nf-1])/2.0;

  for( Index ip=0; ip<np; ip++ )
    {
      // Doppler shift
      //
      Numeric rte_doppler = 0;
      //
      if( ppath_wind_v[ip]!=0 || ppath_wind_u[ip]!=0 || ppath_wind_w[ip]!=0 )
        {
          // za nd aa of winds and LOS
          const Numeric v = sqrt( pow(ppath_wind_u[ip],2) + 
                           pow(ppath_wind_v[ip],2) + pow(ppath_wind_w[ip],2) );
          Numeric aa_w = 0, aa_p = 0;
          if( atmosphere_dim < 3 )
            {
              if( ppath_wind_v[ip]<0 ) 
                { aa_w = PI; }  // Negative v-winds for 1 and 2D
              if( atmosphere_dim == 2  &&  ppath.los(ip,0) < 0 )
                { aa_p = PI; }  // Negative sensor za for 2D
            }
          else //3D 
            { 
              aa_w = atan2( ppath_wind_u[ip], ppath_wind_v[ip] ); 
              aa_p = DEG2RAD * ppath.los(ip,1);
            }
          const Numeric za_w = acos( ppath_wind_w[ip] / v );
          const Numeric za_p = DEG2RAD * fabs( ppath.los(ip,0) );
          //
          // Actual shift
          const Numeric costerm = cos(za_w) * cos(za_p) + 
                                  sin(za_w) * sin(za_p)*cos(aa_w-aa_p); 

          rte_doppler = -v * f0 * costerm / SPEED_OF_LIGHT;
        }

      // We must use temporary variables as the agenda input must be
      // free to be resized
      Vector   evector;
      Matrix   sgmatrix;
      //
      abs_scalar_gas_agendaExecute( ws, sgmatrix, -1, rte_doppler, ppath_p[ip], 
                                    ppath_t[ip], ppath_vmr(joker,ip), 
                                    abs_scalar_gas_agenda );
      ppath_abs_scalar(joker,joker,ip) = sgmatrix;
      //
      if( emission_do )
        {
          emission_agendaExecute( ws, evector, ppath_t[ip], emission_agenda );
          ppath_emission(joker,ip) = evector;
        }

      // Partial and total tau
      //
      if( ip > 0 )
        {
          for( Index iv=0; iv<nf; iv++ )
            { 
              ppath_tau(iv,ip-1) = 0.5 * ppath.lstep[ip-1] * (
                                         ppath_abs_scalar(iv,joker,ip-1).sum() +
                                         ppath_abs_scalar(iv,joker,ip).sum() );
              total_tau[iv] += ppath_tau(iv,ip-1);
            }
        }
    }
}



//! get_ppath_cloudrtvars
/*!
    Determines variables for each step of "standard" RT integration

    See code for details about dimensions etc.

    \param   ws                  Out: The workspace
    \param   ppath_asp_abs_vec   Out: Absorption vectors for absorption species
    \param   ppath_asp_ext_vec   Out: Extinction matrices for absorption species
    \param   ppath_pnd_abs_vec   Out: Absorption vectors for particles
    \param   ppath_pnd_ext_vec   Out: Extinction matrices for particles
    \param   ppath_transmission  Out: Transmissions of each ppath step 
    \param   total_transmission  Out: Total transmission of path
    \param   ppath_emission      Out: Emission source term at each ppath point 
    \param   scat_data           Out: Extracted scattering data. Length of
                                      array affected by *use_mean_scat_data*.
    \param   abs_scalar_gas_agenda As the WSV.    
    \param   emission_agenda     As the WSV.    
    \param   opt_prop_gas_agenda As the WSV.    
    \param   ppath_p             Pressure for each ppath point.
    \param   ppath_t             Temperature for each ppath point.
    \param   ppath_vmr           VMR values for each ppath point.
    \param   ppath_wind_u        U-wind for each ppath point.
    \param   ppath_wind_v        V-wind for each ppath point.
    \param   ppath_wind_w        W-wind for each ppath point.
    \param   ppath_pnd           Particle number densities for each ppath point.
    \param   use_mean_scat_data  As the WSV.    
    \param   scat_data_raw       As the WSV.    
    \param   stokes_dim          As the WSV.    
    \param   f_grid              As the WSV.    
    \param   atmosphere_dim      As the WSV.    
    \param   emission_do         Flag for calculation of emission. Should be
                                 set to 0 for pure transmission calculations.

    \author Patrick Eriksson 
    \date   2011-07-14
*/
void get_ppath_cloudrtvars( 
        Workspace&                     ws,
        Tensor3&                       ppath_asp_abs_vec, 
        Tensor4&                       ppath_asp_ext_mat, 
        Tensor3&                       ppath_pnd_abs_vec, 
        Tensor4&                       ppath_pnd_ext_mat, 
        Tensor4&                       ppath_transmission,
        Tensor3&                       total_transmission,
        Matrix&                        ppath_emission, 
  Array<ArrayOfSingleScatteringData>&  scat_data,
  const Agenda&                        abs_scalar_gas_agenda,
  const Agenda&                        emission_agenda,
  const Agenda&                        opt_prop_gas_agenda,
  const Ppath&                         ppath,
  ConstVectorView                      ppath_p, 
  ConstVectorView                      ppath_t, 
  ConstMatrixView                      ppath_vmr, 
  ConstVectorView                      ppath_wind_u, 
  ConstVectorView                      ppath_wind_v, 
  ConstVectorView                      ppath_wind_w, 
  ConstMatrixView                      ppath_pnd, 
  const Index&                         use_mean_scat_data,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Index&                         stokes_dim,
  ConstVectorView                      f_grid, 
  const Index&                         atmosphere_dim,
  const Index&                         emission_do,
  const Verbosity&                     verbosity)
{
  // Sizes
  const Index   np   = ppath.np;
  const Index   nf   = f_grid.nelem();

  // Init variables
  ppath_asp_abs_vec.resize( nf, stokes_dim, np );
  ppath_asp_ext_mat.resize( nf, stokes_dim, stokes_dim, np );
  ppath_pnd_abs_vec.resize( nf, stokes_dim, np );
  ppath_pnd_ext_mat.resize( nf, stokes_dim, stokes_dim, np );
  ppath_transmission.resize( nf, stokes_dim, stokes_dim, np-1 );
  total_transmission.resize( nf, stokes_dim, stokes_dim );
  //
  if( emission_do )
    { ppath_emission.resize( nf, np ); }
  else
    { ppath_emission.resize( 0, 0 ); }

  // Mean of extreme frequencies
  const Numeric f0 = (f_grid[0]+f_grid[nf-1])/2.0;


  // Particle single scattering properties (are independent of position)
  //
  if( use_mean_scat_data )
    {
      scat_data.resize( 1 );
      scat_data_monoCalc( scat_data[0], scat_data_raw, Vector(1,f0), 0, 
                          verbosity );
    }
  else
    {
      scat_data.resize( nf );
      for( Index iv=0; iv<nf; iv++ )
        { scat_data_monoCalc( scat_data[iv], scat_data_raw, f_grid, iv, 
                              verbosity ); 
        }
    }


  // Loop ppath points
  //
  for( Index ip=0; ip<np; ip++ )
    {
      // Doppler shift 
      // 
      assert( ppath_wind_u[ip] == 0 );
      assert( ppath_wind_v[ip] == 0 );
      assert( ppath_wind_w[ip] == 0 );
      //
      const Numeric rte_doppler = 0;

      // Emission source term
      //
      if( emission_do )
        {
          Vector   evector;   // Agenda must be free to resize
          emission_agendaExecute( ws, evector, ppath_t[ip], emission_agenda );
          ppath_emission(joker,ip) = evector;
        }

      // Absorption species properties 
        {
          Matrix   asgmatrix;   // Agendas must be free to resize
          Matrix   abs_vec( nf, stokes_dim, 0 );
          Tensor3  ext_mat( nf, stokes_dim, stokes_dim, 0 );
          //
          abs_scalar_gas_agendaExecute( ws, asgmatrix, -1, rte_doppler, 
                                        ppath_p[ip], ppath_t[ip], 
                                        ppath_vmr(joker,ip), 
                                        abs_scalar_gas_agenda );
          opt_prop_gas_agendaExecute( ws, ext_mat, abs_vec, -1, asgmatrix, 
                                      opt_prop_gas_agenda );
          ppath_asp_ext_mat(joker,joker,joker,ip) = ext_mat;
          ppath_asp_abs_vec(joker,joker,ip)       = abs_vec;
        }

      // Particle properties 
        {
          // Direction of outgoing scattered radiation (which is reversed to
          // LOS). Note that rte_los2 is only used for extracting scattering
          // properties.
          Vector rte_los2;
          mirror_los( rte_los2, ppath.los(ip,joker), atmosphere_dim );

          // Extinction and absorption
          if( use_mean_scat_data )
            {
              Vector   abs_vec( stokes_dim, 0 );
              Matrix   ext_mat( stokes_dim, stokes_dim, 0 );
              opt_propCalc( ext_mat, abs_vec,
                            rte_los2[0], rte_los2[1], scat_data[0], 
                            stokes_dim, ppath_pnd(joker,ip), ppath_t[ip],
                            verbosity);
              for( Index iv=0; iv<nf; iv++ )
                { 
                  ppath_pnd_ext_mat(iv,joker,joker,ip) = ext_mat;
                  ppath_pnd_abs_vec(iv,joker,ip)       = abs_vec;
                }
            }
          else
            {
              for( Index iv=0; iv<nf; iv++ )
                { 
                  Vector   abs_vec( stokes_dim, 0 );
                  Matrix   ext_mat( stokes_dim, stokes_dim, 0 );
                  opt_propCalc( ext_mat, abs_vec,
                                rte_los2[0], rte_los2[1], scat_data[iv], 
                                stokes_dim, ppath_pnd(joker,ip), ppath_t[ip],
                                verbosity );
                  ppath_pnd_ext_mat(iv,joker,joker,ip) = ext_mat;
                  ppath_pnd_abs_vec(iv,joker,ip)       = abs_vec;
                }
            }
        }
       
      // Transmission
      //
      if( ip == 0 )   // Set total_transmission to identity matrices
        {
          for( Index iv=0; iv<nf; iv++ )
            { id_mat( total_transmission(iv,joker,joker) ); } 
        }
      else
        {
          for( Index iv=0; iv<nf; iv++ )
            { 
              // Average extinction matrix
              Matrix  ext_mat_av( stokes_dim, stokes_dim,0 );
              for( Index is1=0; is1<stokes_dim; is1++ )
                { 
                  for( Index is2=0; is2<stokes_dim; is2++ )
                    {
                      ext_mat_av(is1,is2) = 0.5 * (
                                          ppath_asp_ext_mat(iv,is1,is2,ip-1) +
                                          ppath_asp_ext_mat(iv,is1,is2,ip)   +
                                          ppath_pnd_ext_mat(iv,is1,is2,ip-1) +
                                          ppath_pnd_ext_mat(iv,is1,is2,ip) );
                    }
                }
              // Transmission
              ext2trans( ppath_transmission(iv,joker,joker,ip-1),
                         ext_mat_av, ppath.lstep[ip-1] );  
              const Matrix tmp = total_transmission(iv,joker,joker);
              mult( total_transmission(iv,joker,joker), tmp,
                    ppath_transmission(iv,joker,joker,ip-1) );
            }
        }
    }
}



//! get_rowindex_for_mblock
/*!
    Returns the "range" of *y* corresponding to a measurement block

    \return  The range.
    \param   sensor_response    As the WSV.
    \param   mblock_index            Index of the measurement block.

    \author Patrick Eriksson 
    \date   2009-10-16
*/
Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    mblock_index )
{
  const Index   n1y = sensor_response.nrows();
  return Range( n1y*mblock_index, n1y );
}



//! iyb_calc
/*!
    Calculation of pencil beam monochromatic spectra for 1 measurement block.

    All in- and output variables as the WSV with the same name.

    \author Patrick Eriksson 
    \date   2009-10-16
*/
void iyb_calc(
        Workspace&                  ws,
        Vector&                     iyb,
        Vector&                     iyb_error,
        Index&                      iy_error_type,
        Vector&                     iyb_aux,
        ArrayOfMatrix&              diyb_dx,
  const Index&                      mblock_index,
  const Index&                      atmosphere_dim,
  ConstTensor3View                  t_field,
  ConstTensor3View                  z_field,
  ConstTensor4View                  vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  ConstVectorView                   f_grid,
  ConstMatrixView                   sensor_pos,
  ConstMatrixView                   sensor_los,
  ConstVectorView                   mblock_za_grid,
  ConstVectorView                   mblock_aa_grid,
  const Index&                      antenna_dim,
  const Agenda&                     iy_clearsky_agenda,
  const String&                     y_unit,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const Verbosity&                  verbosity )
{
  // Sizes
  const Index   nf   = f_grid.nelem();
  const Index   nza  = mblock_za_grid.nelem();
        Index   naa  = mblock_aa_grid.nelem();   
  if( antenna_dim == 1 )  
    { naa = 1; }
  const Index   niyb = nf * nza * naa * stokes_dim;

  // Set up size of containers for data of 1 measurement block.
  // (can not be made below due to parallalisation)
  iyb.resize( niyb );
  iyb_error.resize( niyb );
  iyb_aux.resize( niyb );
  Index aux_set = 0;
  //
  iyb_error = 0;
  //
  if( j_analytical_do )
    {
      diyb_dx.resize( jacobian_indices.nelem() );
      FOR_ANALYTICAL_JACOBIANS_DO(
        diyb_dx[iq].resize( niyb, jacobian_indices[iq][1] -
                                  jacobian_indices[iq][0] + 1 );
      )
    }
  else
    { diyb_dx.resize( 0 ); }

  // Polarisation index variable
  // If some kind of modified Stokes vector will be introduced, this must be
  // made to a WSV
  ArrayOfIndex i_pol(stokes_dim);
  for( Index is=0; is<stokes_dim; is++ )
    { i_pol[is] = is + 1; }

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_iy_clearsky_agenda (iy_clearsky_agenda);
  
  // Start of actual calculations
/*#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && nza>1)                            \
  default(none)                                                   \
  firstprivate(l_ws, l_iy_clearsky_agenda)                        \
  shared(sensor_los, mblock_za_grid, mblock_aa_grid, vmr_field,   \
         t_field, f_grid, sensor_pos, \
         joker, naa) */
#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && nza>1)                            \
  firstprivate(l_ws, l_iy_clearsky_agenda)
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
              los     = sensor_los( mblock_index, joker );
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

              // Calculate iy and associated variables
              //
              Matrix         iy, iy_error, iy_aux;
              ArrayOfTensor3 diy_dx;
              //
              iyCalc( l_ws, iy, iy_aux, iy_error, iy_error_type, diy_dx, 
                      1, t_field, z_field, vmr_field, cloudbox_on, 1, 
                      sensor_pos(mblock_index,joker), los, j_analytical_do, 
                      mblock_index, l_iy_clearsky_agenda, verbosity );

              // Start row in iyb etc. for present LOS
              //
              const Index row0 = ( iza*naa + iaa ) * nf * stokes_dim;

              // Jacobian part (must be converted to Tb before iy for PlanckBT)
              // 
              if( j_analytical_do )
                {
                  FOR_ANALYTICAL_JACOBIANS_DO(
                    //
                    apply_y_unit2( diy_dx[iq], iy, y_unit, f_grid, i_pol );
                    //
                    for( Index ip=0; ip<jacobian_indices[iq][1] -
                                        jacobian_indices[iq][0]+1; ip++ )
                      {
                        for( Index is=0; is<stokes_dim; is++ )
                          { 
                            diyb_dx[iq](Range(row0+is,nf,stokes_dim),ip)=
                                                     diy_dx[iq](ip,joker,is); 
                          }
                      }                              
                  )
                }

              // iy       : unit conversion and copy to iyb
              // iy_error : unit conversion and copy to iyb_error
              // iy_aux   : copy to iyb_aux (if iy_aux filled)
              //
              if( iy_error_type > 0 )
                { apply_y_unit2( iy_error, iy, y_unit, f_grid, i_pol ); }
              //
              apply_y_unit( iy, y_unit, f_grid, i_pol );
              //
              for( Index is=0; is<stokes_dim; is++ )
                { 
                  iyb[Range(row0+is,nf,stokes_dim)] = iy(joker,is); 
                  //
                  if( iy_error_type > 0 )
                    { 
                      iyb_error[Range(row0+is,nf,stokes_dim)] = 
                                                            iy_error(joker,is); 
                    }
                  //
                  if( iy_aux.nrows() )
                    {
                      iyb_aux[Range(row0+is,nf,stokes_dim)] = iy_aux(joker,is);
                      aux_set = 1;
                    }
                }
            }  // End aa loop
        }  // End try

      catch (runtime_error e)
        {
          CREATE_OUT0
          exit_or_rethrow(e.what(), out0);
        }
    }  // End za loop

  // If no aux, set to size 0 to flag this
  if( !aux_set )
    { iyb_aux.resize(0); }
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

  assert( iy_transmission.ncols() == iy_transmission.nrows() );
  assert( nf == tau.nelem() );

  iy_trans_new = iy_transmission;

  for( Index iv=0; iv<nf; iv++ )
    { iy_trans_new(iv,joker,joker) *= exp( -tau[iv] ); } 
}



//! mirror_los
/*!
    Determines the backward direction for a given line-of-sight.

    This function can be used to get the LOS to apply for extracting single
    scattering properties, if the propagation path LOS is given.

    A viewing direction of aa=0 is assumed for 1D. This corresponds to 
    positive za for 2D.

    \param   los_mirrored      Out: The line-of-sight for reversed direction.
    \param   los               A line-of-sight
    \param   atmosphere_dim    As the WSV.

    \author Patrick Eriksson 
    \date   2011-07-15
*/
void mirror_los(
        Vector&     los_mirrored,
  ConstVectorView   los, 
  const Index&      atmosphere_dim )
{
  los_mirrored.resize(2);
  //
  if( atmosphere_dim == 1 )
    { 
      los_mirrored[0] = 180 - los[0]; 
      los_mirrored[1] = 180; 
    }
  else if( atmosphere_dim == 2 )
    {
      los_mirrored[0] = 180 - fabs( los[0] ); 
      if( los[0] >= 0 )
        { los_mirrored[1] = 180; }
      else
        { los_mirrored[1] = 0; }
    }
  else if( atmosphere_dim == 3 )
    { 
      los_mirrored[0] = 180 - los[0]; 
      los_mirrored[1] = los[1] + 180; 
      if( los_mirrored[1] > 180 )
        { los_mirrored[1] -= 360; }
    }
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
    and the propagation path length through the layer is *lstep*.

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
    \param   lstep             Input: The length of the RTE step.
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
              const Numeric& lstep,
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
  assert( lstep >= 0 );
  assert (!is_singular( ext_mat_av ));

  // Any changes here associated with the extinction matrix should also be
  // implemented in ext2mat.

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
      trans_mat(0,0) = exp(-ext_mat_av(0,0) * lstep);
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
      trans_mat(0,0) = exp(-ext_mat_av(0,0) * lstep);

      // Stokes dim 1
      //   assert( ext_mat_av(0,0) == abs_vec_av[0] );
      //   Numeric transm = exp( -lstep * abs_vec_av[0] );
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
      ext_mat_ds *= -lstep; // ext_mat_ds = -ext_mat*ds

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



//! surface_calc
/*!
    Weights together downwelling radiation and surface emission.

    *iy* must have correct size when function is called.

    \param   iy                 In/Out: Radiation matrix, amtching 
                                        the WSV with the same name.
    \param   I                  Input: Downwelling radiation, with dimensions
                                       matching: 
                                       (surface_los, f_grid, stokes_dim)
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



//! trans_step_std
/*!
    Solves monochromatic (vector) Beer-Lambert for one step

    The function can be used for clearsky and cloudbox calculations.

    (As rte_step_std, but without emission and scattering source term)

    When calling the function, the vector *stokes_vec* shall contain the
    Stokes vector for the incoming radiation. The function returns this
    vector, then containing the outgoing radiation on the other side of the 
    layer.

    Transmission calculated by *ext2trans*.

    \param   stokes_vec         Input/Output: A Stokes vector.
    \param   trans_mat          Input/Output: Transmission matrix of slab.
    \param   ext_mat            Input: Averaged extinction matrix.
    \param   lstep             Input: The length of the RTE step.

    \author Claudia Emde and Patrick Eriksson, 
    \date   2010-10-15
*/
void trans_step_std(//Output and Input:
              VectorView stokes_vec,
              MatrixView trans_mat,
              //Input
              ConstMatrixView ext_mat_av,
              const Numeric& lstep )
{
  // Checks made in *ext2trans*.
  ext2trans( trans_mat, ext_mat_av, lstep );  

  Vector tmp = stokes_vec;

  mult( stokes_vec, trans_mat, tmp );
}

