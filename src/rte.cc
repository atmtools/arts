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

extern const Numeric BOLTZMAN_CONST;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;



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

    \author Patrick Eriksson 
    \date   2010-04-07
*/
void apply_y_unit( 
      MatrixView      iy, 
    const String&     y_unit, 
    ConstVectorView   f_grid )
{
  assert( f_grid.nelem() == iy.nrows() );

  // The code is largely identical between the two apply_y_unit functions.
  // If any change here, remember to update the other function.

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
    {
      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          const Numeric scfac = invrayjean( 1, f_grid[iv] );
          for( Index is=0; is<iy.ncols(); is++ )
            {
              iy(iv,is) *= scfac;
            }
        }
    }

  else if( y_unit == "PlanckBT" )
    {
      // Use always double to avoid numerical problem (see invrayjean)
      static const double a = PLANCK_CONST / BOLTZMAN_CONST;
      static const double b = 2*PLANCK_CONST / (SPEED_OF_LIGHT*SPEED_OF_LIGHT);

      const Index ns = iy.ncols();

      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          const double c = a * f_grid[iv]; 
          const double d = b * pow(f_grid[iv],3); 

          if( ns == 1 )
            {
              iy(iv,0) = c / log( d/iy(iv,0) + 1.0 ); 
            }
          else
            {
              const Numeric i = iy(iv,0);
              iy(iv,0) = c / log( d / i + 1.0 ); 
              static const double e = pow(iy(iv,0),2)/(c*i*(1+i/d));
              for( Index is=1; is<ns; is++ )
                {
                  iy(iv,is) *= e;
                }
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
    Tensor3View       J,
    ConstMatrixView   iy, 
    const String&     y_unit, 
    ConstVectorView   f_grid )
{
  assert( J.ncols() == iy.ncols() );
  assert( f_grid.nelem() == iy.nrows() );
  assert( f_grid.nelem() == J.nrows() );
  assert( max(iy) < 1e-3 );   // If fails, iy is already in Tb

  // The code is largely identical between the two apply_y_unit functions.
  // If any change here, remember to update the other function.

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
    {
      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          const Numeric scfac = invrayjean( 1, f_grid[iv] );
          for( Index is=0; is<J.ncols(); is++ )
            {
              for( Index ip=0; ip<J.npages(); ip++ )
                {
                  J(ip,iv,is) *= scfac;
                }
            }
        }
    }

  else if( y_unit == "PlanckBT" )
    {
      // Use always double to avoid numerical problem (see invrayjean)
      static const double a = PLANCK_CONST / BOLTZMAN_CONST;
      static const double b = 2*PLANCK_CONST / (SPEED_OF_LIGHT*SPEED_OF_LIGHT);

      const Index ns = J.ncols();

      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          const double  c = a * f_grid[iv]; 
          const double  d = b * pow(f_grid[iv],3); 
          const Numeric i = iy(iv,0);
          const Numeric y = c / log( d / i + 1.0 ); 
          const Numeric e = pow(y,2) / ( c * i * ( 1.0 + i/d ) );

          for( Index ip=0; ip<J.npages(); ip++ )
            {
              for( Index is=0; is<ns; is++ )
                {
                  J(ip,iv,is) *= e;
                }
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
  interp_atmfield_gp2itw( itw_field, atmosphere_dim, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
  interp_atmfield_by_itw( ppath_t,  atmosphere_dim, t_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );

  //  VMR fields:
  const Index ns = vmr_field.nbooks();
  ppath_vmr.resize(ns,np);
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( ppath_vmr(is, joker), atmosphere_dim,
                              vmr_field( is, joker, joker,  joker ), 
                              ppath.gp_p, ppath.gp_lat, ppath.gp_lon, 
                              itw_field );
    }
}



//! get_rowindex_for_mblock
/*!
    Returns the "range" of *y* corresponding to a measurement block

    \return  The range.
    \param   sensor_response    As the WSV.
    \param   imblock            Index of the measurement block.

    \author Patrick Eriksson 
    \date   2009-10-16
*/
Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    imblock )
{
  const Index   n1y = sensor_response.nrows();
  return Range( n1y*imblock, n1y );
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

