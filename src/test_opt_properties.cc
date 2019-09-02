/* Copyright (C) 2002-2012 Sreerekha Ravi <rekha@sat.physik.uni-bremen.de>

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. *

/*!
  \file   test_opt_properties.cc
  \author Sreerekha T.R. 
  \date   Mon May 13 11:36:11 2002
*/

#include <iostream>
#include "scatproperties.h"
extern const Numeric DEG2RAD;
int main(void) {
  Vector amp(8);
  Vector za(18), aa(36);
  Matrix l_cross(4, 4);
  Numeric freq = 325 * 1e+9;
  // elements of amplitude matrix
  amp[0] = 1;
  amp[1] = 0;
  amp[2] = 0;
  amp[3] = 0;
  amp[4] = 0;
  amp[5] = 0;
  amp[6] = 0;
  amp[7] = 0;

  cout << "\n Extinction Matrix  test \n";
  cout << "Amplitude matrix: \n";
  cout << amp;

  cout << "\n";
  za[0] = 0;
  aa[0] = 0;
  amp2ext(l_cross, amp, freq);
  cout << "Extinction Matrix: \n";
  for (Index i = 0; i < 4; i++) {
    for (Index j = 0; j < 4; j++)
      //ext_matrix(l_cross,amp,freq);
      cout << " " << l_cross(i, j);
  }
  cout << "\n";
  /*===========================================================
        checking integration
        ==========================================================*/
  Vector A(361);
  Vector A1(361);
  Numeric X;
  Numeric llim = 0.;
  Numeric ulim = 180.;
  Numeric N = 40;
  Numeric h = (ulim - llim) / N;
  cout << "Lower Limit :" << llim << "\n";
  cout << "Upper Limit :" << ulim << "\n";
  cout << "Number of points :" << N << "\n";
  cout << "Step length :" << h << "\n";
  for (Numeric theta = 0; theta < ulim + 1; theta = theta + 1) {
    A[theta] = sin(theta * DEG2RAD);
    //A[theta]=theta;
    A1[theta] = DEG2RAD * A[theta];
    //cout<<theta<<"  "<<A[theta] <<"\n";
  }
  single_trapez(X, A1, llim, ulim, h);
  cout << "Integral of the Expression "
       << " " << X << "\n";
  /*=========================================================
        testing for double integration============================
        =======================================================*/
  Matrix B(181, 361);
  Matrix B1(181, 361);
  //Vector B1(5);
  //Vector X2(3);
  //Numeric X1;
  Numeric Y2;
  //Index N1=4;
  Numeric llim2 = 0;
  Numeric ulim2 = 180;
  Numeric llim1 = 0.;
  Numeric ulim1 = 180;
  Numeric N1 = 180;
  Numeric h1 = (ulim1 - llim1) / N1;
  Numeric N2 = 180;
  Numeric h2 = (ulim2 - llim2) / N2;
  cout << "Number of points inner integral :" << N1 << "\n";
  cout << "Step length - inner integral:" << h1 << "\n";
  cout << "Number of points outer integral:" << N2 << "\n";
  cout << "Step length - outer integral:" << h2 << "\n";
  for (Numeric i = llim2; i < ulim2 + 1; i = i + 1) {
    for (Numeric j = llim1; j < ulim1 + 1; j = j + 1) {
      B(i, j) = sin(j * DEG2RAD);
      B1(i, j) = DEG2RAD * B(i, j);
    }
  }
  double_trapez(Y2, B1, llim1, llim2, ulim1, ulim2, h1, h2);
  cout << " Double Integral of the Expression "
       << " " << Y2 << "\n";
  /*=========================================================
        testing for double integration through successive single integration
        =======================================================*/
  Matrix C(181, 361);  //Defining the matrix dimension
  Vector CC1(181);
  Numeric CCC, CC;
  for (Numeric i = llim2; i < ulim2 + 1; i = i + 1) {
    for (Numeric j = llim1; j < ulim1 + 1; j = j + 1) {
      C(i, j) = sin(j * DEG2RAD);  //Matrix definition;simple sin(theta)
    }
  }
  for (Numeric i = llim2; i < ulim2 + 1; i = i + 1) {
    Vector C1 = C(i, Range(joker));
    C1 *= DEG2RAD;
    //trapezoidal integration for the columns corresponding to each rows
    single_trapez(CC, C1, llim1, ulim1, h1);
    CC1[i] = CC;
  }
  //trapezoidal integration for the  rows
  single_trapez(CCC, CC1, llim2, ulim2, h2);
  cout
      << " Double Integral of the Expression through successive single integrations "
      << " " << CCC << "\n";
  return 0;
}
