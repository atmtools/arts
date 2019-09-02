/* Copyright (C) 2005-2012 Cory Davis <cory@met.ed.ac.uk>
                            
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

Numeric SLIData2::interpolate(Numeric x1, Numeric x2) const {
  GridPos gp1, gpl, gpr;
  Vector itw1(2), itwl(2), itwr(2);
  Numeric yl, yr;
  //Get interpolation weights for x1 interpolation
  gridpos(gp1, this->x1a, x1);
  interpweights(itw1, gp1);
  //Get y values on bounding x1 grid points for desired x2
  gridpos(gpl, this->x2a[gp1.idx], x2);
  interpweights(itwl, gpl);
  gridpos(gpr, this->x2a[gp1.idx + 1], x2);
  interpweights(itwr, gpr);
  yl = interp(itwl, this->ya[gp1.idx], gpl);
  yr = interp(itwr, this->ya[gp1.idx + 1], gpr);
  //interpolate these two y values useing x1 interpolation weights
  return itw1[0] * yl + itw1[1] * yr;
}

//void SLIData2::check() const
//{
//  Index nx1=this->x1a.nelem();
//  assert(nx1>0);
//}

ostream& operator<<(ostream& os, const SLIData2& /* sli */) {
  os << "SLIData2    : Output operator not implemented";
  return os;
}

//! Red 1D Interpolate.
/*! 
  This is a slight modifiaction of Stefan's code to do 1_D interpolation
  to get a Matrix from an array of Matrices

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param[out] tia  Interpolated value.
  \param[in]  itw  Interpolation weights.
  \param[in]  a    The field to interpolate.(ArrayOfMatrix)
  \param[in]  tc   The grid position for the column dimension.


  \author Cory Davis (modified original code by Stefan Buehler)
  \date   2003-06-19
*/

void interp(MatrixView tia,
            ConstVectorView itw,
            const ArrayOfMatrix& a,
            const GridPos& tc) {
  DEBUG_ONLY(const Numeric sum_check_epsilon = 1e-6);

  assert(is_size(itw, 2));  // We need 2 interpolation
                            // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  assert(is_same_within_epsilon(itw.sum(), 1, sum_check_epsilon));

  Index anr = a[0].nrows();
  Index anc = a[0].ncols();

  assert(tia.nrows() == anr);
  assert(tia.ncols() == anc);

  for (Index inr = 0; inr < anr; inr++)
    for (Index inc = 0; inc < anc; inc++) {
      tia(inr, inc) =
          a[tc.idx](inr, inc) * itw[0] + a[tc.idx + 1](inr, inc) * itw[1];
    }
}

//! Red 1D Interpolate.
/*! 
  This is a slight modifiaction of Stefan's code to do 1_D interpolation
  to get a Vector from an array of Vectors

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param[out] tia  Interpolated value.
  \param[in]  itw  Interpolation weights.
  \param[in]  a    The field to interpolate. (ArrayOfVector)
  \param[in]  tc   The grid position for the column dimension.

  \author Cory Davis (modified original code by Stefan Buehler)
  \date   2003-06-19

*/
void interp(VectorView tia,
            ConstVectorView itw,
            const ArrayOfVector& a,
            const GridPos& tc) {
  DEBUG_ONLY(const Numeric sum_check_epsilon = 1e-6);
  assert(is_size(itw, 2));  // We need 2 interpolation
                            // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  assert(is_same_within_epsilon(itw.sum(), 1, sum_check_epsilon));

  Index an = a[0].nelem();

  assert(tia.nelem() == an);

  for (Index i = 0; i < an; ++i) {
    tia[i] = a[tc.idx][i] * itw[0] + a[tc.idx + 1][i] * itw[1];
  }
}

void interp_scat_angle_temperature(  //Output:
    VectorView pha_mat_int,
    Numeric& theta_rad,
    //Input:
    const SingleScatteringData& scat_data_single,
    const Numeric& za_sca,
    const Numeric& aa_sca,
    const Numeric& za_inc,
    const Numeric& aa_inc,
    const Numeric& rtp_temperature) {
  Numeric ANG_TOL = 1e-7;

  //Calculate scattering angle from incident and scattered directions.
  //The two special cases are implemented here to avoid NaNs that can
  //sometimes occur in in the acos... formula in forward and backscatter
  //cases. CPD 5/10/03.

  if (abs(aa_sca - aa_inc) < ANG_TOL) {
    theta_rad = DEG2RAD * abs(za_sca - za_inc);
  } else if (abs(abs(aa_sca - aa_inc) - 180) < ANG_TOL) {
    theta_rad = DEG2RAD * (za_sca + za_inc);
    if (theta_rad > PI) {
      theta_rad = 2 * PI - theta_rad;
    }
  } else {
    const Numeric za_sca_rad = za_sca * DEG2RAD;
    const Numeric za_inc_rad = za_inc * DEG2RAD;
    const Numeric aa_sca_rad = aa_sca * DEG2RAD;
    const Numeric aa_inc_rad = aa_inc * DEG2RAD;

    // cout << "Interpolation on scattering angle" << endl;
    assert(scat_data_single.pha_mat_data.ncols() == 6);
    // Calculation of the scattering angle:
    theta_rad =
        acos(cos(za_sca_rad) * cos(za_inc_rad) +
             sin(za_sca_rad) * sin(za_inc_rad) * cos(aa_sca_rad - aa_inc_rad));
  }
  const Numeric theta = RAD2DEG * theta_rad;

  // Interpolation of the data on the scattering angle:

  GridPos thet_gp;
  gridpos(thet_gp, scat_data_single.za_grid, theta);
  GridPos t_gp;

  if (scat_data_single.T_grid.nelem() == 1) {
    Vector itw(2);
    interpweights(itw, thet_gp);

    for (Index i = 0; i < 6; i++) {
      pha_mat_int[i] = interp(
          itw, scat_data_single.pha_mat_data(0, 0, joker, 0, 0, 0, i), thet_gp);
    }
  } else {
    gridpos(t_gp, scat_data_single.T_grid, rtp_temperature);

    Vector itw(4);
    interpweights(itw, t_gp, thet_gp);

    for (Index i = 0; i < 6; i++) {
      pha_mat_int[i] =
          interp(itw,
                 scat_data_single.pha_mat_data(0, joker, joker, 0, 0, 0, i),
                 t_gp,
                 thet_gp);
    }
  }
}
