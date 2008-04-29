/* Copyright (C) 2005-2008 Cory Davis <cory@met.ed.ac.uk>
                            
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
  \file   mc_interp.cc
  \author Cory Davis <cory@met.ed.ac.uk>
  \date   2005-02-28 

  \brief  Interpolation classes and functions created for use within Monte 
  Carlo scattering simulations 

*/
/*===========================================================================
  === External declarations
  ===========================================================================*/
#include "mc_interp.h"
#include "logic.h"
#include "montecarlo.h"


//! Perform sequential interpolation
/*!
  \param x1 desired x1 value
  \param x2 desired x2 value

  \return interpolated y value at x1,x2

  \author Cory Davis <cdavis@staffmail.ed.ac.uk>
  \date 2005-02-28
*/

Numeric SLIData2::interpolate(Numeric x1, Numeric x2) const
{
  GridPos gp1,gpl,gpr;
  Vector itw1(2),itwl(2),itwr(2);
  Numeric yl,yr;
  //Get interpolation weights for x1 interpolation
  gridpos(gp1,this->x1a,x1);
  interpweights(itw1,gp1);
  //Get y values on bounding x1 grid points for desired x2
  gridpos(gpl,this->x2a[gp1.idx],x2);
  interpweights(itwl,gpl);
  gridpos(gpr,this->x2a[gp1.idx+1],x2);
  interpweights(itwr,gpr);
  yl=interp(itwl,this->ya[gp1.idx],gpl);
  yr=interp(itwr,this->ya[gp1.idx+1],gpr);
  //interpolate these two y values useing x1 interpolation weights
  return itw1[0]*yl+itw1[1]*yr;
}

//void SLIData2::check() const
//{
//  Index nx1=this->x1a.nelem();
//  assert(nx1>0);  
//}


ostream& operator<< (ostream& os, const SLIData2& /* sli */)
{
  os << "SLIData2    : Output operator not implemented";
  return os;
}









                
//! Red 1D Interpolate.
/*! 
  This is a slight modifiaction of Stefan's code to do 1_D interpolation
  to get a Matrix from an array of Matrices

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate.(ArrayOfMatrix)
  \param tc   The grid position for the column dimension.

  \return Interpolated value.

  \author Cory Davis (modified original code by Stefan Buehler)
  \date   2003-06-19
*/

Matrix interp( ConstVectorView itw,
                ArrayOfMatrix a,    
                const GridPos&  tc )
{
  DEBUG_ONLY (const Numeric sum_check_epsilon = 1e-6);

  assert(is_size(itw,2));       // We need 2 interpolation
                                // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  assert( is_same_within_epsilon( itw.sum(),
                                  1,
                                  sum_check_epsilon ) );
  
  // To store interpolated value:
  Matrix tia(a[0].nrows(),a[0].ncols(),0.0);
  Matrix b(a[0].nrows(),a[0].ncols(),0.0);
  Index iti = 0;
  for ( Index c=0; c<2; ++c )
    {
      b=a[tc.idx+c];
      b*=itw[iti];   
      tia += b;
      ++iti;
    }

  return tia;
}               


//! Red 1D Interpolate.
/*! 
  This is a slight modifiaction of Stefan's code to do 1_D interpolation
  to get a Vector from an array of Vectors

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate. (ArrayOfVector)
  \param tc   The grid position for the column dimension.

  \return Interpolated value.
  \author Cory Davis (modified original code by Stefan Buehler)
  \date   2003-06-19

*/

Vector interp( ConstVectorView itw,
                ArrayOfVector a,    
                const GridPos&  tc )
{
  DEBUG_ONLY (const Numeric sum_check_epsilon = 1e-6);
  
  assert(is_size(itw,2));       // We need 2 interpolation
                                // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  assert( is_same_within_epsilon( itw.sum(),
                                  1,
                                  sum_check_epsilon ) );
  
  // To store interpolated value:
  Vector tia(a[0].nelem(),0.0);
  Vector b(a[0].nelem(),0.0);
  Index iti = 0;
  for ( Index c=0; c<2; ++c )
    {
      b=a[tc.idx+c];
      b*=itw[iti];   
      tia += b;
      ++iti;
    }

  return tia;
}               

void interp_scat_angle_temperature(//Output:
                                   VectorView pha_mat_int,
                                   Numeric& theta_rad,
                                   //Input:
                                   const SingleScatteringData& scat_data,
                                   const Numeric& za_sca,
                                   const Numeric& aa_sca,
                                   const Numeric& za_inc,
                                   const Numeric& aa_inc,
                                   const Numeric& rte_temperature
                                   )
{
  Numeric ANG_TOL=1e-7;

  //Calculate scattering angle from incident and scattered directions.
  //The two special cases are implemented here to avoid NaNs that can 
  //sometimes occur in in the acos... formula in forward and backscatter
  //cases. CPD 5/10/03.
  
  if(abs(aa_sca-aa_inc)<ANG_TOL)
    {
      theta_rad=DEG2RAD*abs(za_sca-za_inc);
    }
  else if (abs(abs(aa_sca-aa_inc)-180)<ANG_TOL)
    {
      theta_rad=DEG2RAD*(za_sca+za_inc);
      if (theta_rad>PI){theta_rad=2*PI-theta_rad;}
    }
  else
    {
      const Numeric za_sca_rad = za_sca * DEG2RAD;
      const Numeric za_inc_rad = za_inc * DEG2RAD;
      const Numeric aa_sca_rad = aa_sca * DEG2RAD;
      const Numeric aa_inc_rad = aa_inc * DEG2RAD;
      
      // cout << "Interpolation on scattering angle" << endl;
      assert (scat_data.pha_mat_data.ncols() == 6);
      // Calculation of the scattering angle:
      theta_rad = acos(cos(za_sca_rad) * cos(za_inc_rad) + 
                       sin(za_sca_rad) * sin(za_inc_rad) * 
                       cos(aa_sca_rad - aa_inc_rad));
   }
      const Numeric theta = RAD2DEG * theta_rad;
      
  // Interpolation of the data on the scattering angle:
 
  GridPos thet_gp;
  gridpos(thet_gp, scat_data.za_grid, theta);
  GridPos t_gp;
  gridpos(t_gp, scat_data.T_grid, rte_temperature);
          
  Vector itw(4);
  interpweights(itw, t_gp,thet_gp);

  for (Index i = 0; i < 6; i++)
    {
      pha_mat_int[i] = interp(itw, scat_data.pha_mat_data(0,joker,joker, 0, 0, 0, i), 
                              t_gp,thet_gp);
    }
} 



//! interpTarray

/*!
   Interpolates several arrays calculated by TarrayCalc to give values at a 
   given pathlength

   \param[out]  T             transmittance matrix ( I may have made this term up ).
   \param[out]  K_abs         absorption coefficient vector
   \param[out]  temperature   temperature
   \param[out]  K             extinction matrix at interpolation point
   \param[out]  rte_pos       position at pathlength along ppath
   \param[out]  rte_los       LOS at pathlength along ppath
   \param[in]   pnd_vec       pnd vector
   \param[in]   TArray        array of transmittance matrices
   \param[in]   ext_matArray  array of extinction matrices
   \param[in]   abs_vecArray  array of absorption coefficients
   \param[in]   t_ppath       array of temperatures
   \param[in]   pnd_ppath     array of pressures
   \param[in]   cum_l_step    vector of cumulative pathlengths
   \param[in]   pathlength    pathlength at which to calculate above values
   \param[in]   stokes_dim    length of Stokes vector
   \param[in]   ppath         the Ppath


   \author Cory Davis
   \date   2003-06-19
*/


//interpolates TArray and PPath to give T and rte_pos(los) at a given pathlength
void interpTArray(Matrix& T,
                  Vector& K_abs,
                  Numeric& temperature,
                  MatrixView& K,
                  Vector& rte_pos,
                  Vector& rte_los,
                  VectorView& pnd_vec,
                  const ArrayOfMatrix& TArray,
                  const ArrayOfMatrix& ext_matArray,
                  const ArrayOfVector& abs_vecArray,
                  const Vector& t_ppath,
                  const Matrix& pnd_ppath,
                  const Vector& cum_l_step,
                  const Numeric& pathlength,
                  const Index& stokes_dim,
                  const Ppath& ppath
                  )
{
  //Internal Declarations
  Matrix incT(stokes_dim,stokes_dim,0.0);
  Matrix opt_depth_mat(stokes_dim,stokes_dim);
  Vector itw(2);
  Numeric delta_s;
  Index N_pt=pnd_vec.nelem();
  ArrayOfGridPos gp(1);
                  
  //interpolate transmittance matrix
  gridpos(gp, cum_l_step, pathlength);
  interpweights(itw,gp[0]);
  K = interp(itw,ext_matArray,gp[0]);
  delta_s = pathlength - cum_l_step[gp[0].idx];
  opt_depth_mat = K;
  opt_depth_mat+=ext_matArray[gp[0].idx];
  opt_depth_mat*=-delta_s/2;
  if ( stokes_dim == 1 )
    {
      incT(0,0)=exp(opt_depth_mat(0,0));
    }
  else if ( is_diagonal( opt_depth_mat))
    {
      for ( Index i=0;i<stokes_dim;i++)
        {
          incT(i,i)=exp(opt_depth_mat(i,i));
        }
    }
  else
    {
      //matrix_exp(incT,opt_depth_mat,10);
      matrix_exp_p30(incT,opt_depth_mat);
    }
  mult(T,TArray[gp[0].idx],incT);
  
  K_abs = interp(itw, abs_vecArray,gp[0]);
 
  temperature=interp(itw,t_ppath,gp[0]);

  for (Index i=0;i<N_pt;i++)
    {
      pnd_vec[i]=interp(itw,pnd_ppath(i,Range(joker)),gp[0]);
    }

  for (Index i=0; i<2; i++)
    {
      rte_pos[i] = interp(itw,ppath.pos(Range(joker),i),gp[0]);
      rte_los[i] = interp(itw,ppath.los(Range(joker),i),gp[0]);
    }
  rte_pos[2] = interp(itw,ppath.pos(Range(joker),2),gp[0]);
}
