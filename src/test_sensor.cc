/* Copyright (C) 2003 Mattias Ekström <ekstrom@rss.chalmers.se>

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

#include "sensor.h"
#include <cmath>
#include <iostream>
#include "matpackI.h"
#include "xml_io.h"
//#include "array.h"
//#include "make_array.h"
//#include "mystring.h"
//#include "make_vector.h"
//#include "math_funcs.h"
//#include "describe.h"
//#include "matpackII.h"

extern const Numeric DEG2RAD;
extern const Numeric PI;

//! antenna_diagram_gaussian
/*! 
  This function has moved here from sensor.cc where it became obsolete.
  
  Sets up a matrix containing a standardised Gaussian antenna diagram,
  described for a certain frequency. The function is called with the
  half-power beam width that determines the shape of the curve for the
  reference frequency, and a matrix that must contain the grid that sets
  up the antenna diagram as the first column. The antenna diagram values
  will then be put in the second column.

  \param   srm     The antenna diagram matrix.
  \param   theta   The antenna average width.
  
  \author Mattias Ekström
  \date   2003-06-18
*/
void antenna_diagram_gaussian(
           MatrixView   srm,
       const Numeric&   theta )
{
  //Assert that srm has the rigth size
  assert( srm.ncols()==2 );
  
  //Initialise variables
  Numeric ln2 = log(2.0);
  
  //Loop over grid points to calculate antenna diagram
  for( Index i=0; i<srm.nrows(); i++ )
    srm(i,1)=exp(-4*ln2*pow(srm(i,0)*DEG2RAD/theta,2));
}

void test1()
{
  //Test antenna_transfer_matrix with sparse matrix
  cout << "\nTest 1:\n";

  Sparse H(3,15);
  Vector m_za(-8,5,4);
  Matrix srm(21,2);
  Vector a_grid(-10,21,1);
  srm(Range(joker),0) = a_grid;
  Vector f_grid(2,3,2);

  antenna_diagram_gaussian(srm, 2);

  antenna_transfer_matrix( H, m_za, srm, f_grid );

  cout << "H:\n" << H << "\n";
  cout << "m_za:\n" << m_za << "\n";
  cout << "srm:\n" << srm << "\n";
  cout << "f_grid:\n" << f_grid << "\n";
}

void test2()
{
  //Test spectrometer_transfer_matrix
  cout << "\nTest 2:\n";

  Sparse H(4,5);
  Matrix srm(7,2);
  Vector x_srm(-3,7,1);
//  Vector values_srm(x_srm.nelem());
//  antenna_diagram_gaussian(values_srm, x_srm, 2);
  srm(Range(joker),0) = x_srm;
  antenna_diagram_gaussian(srm, 2);
//  srm(Range(joker),2) = values_srm;
  Vector f_grid(1,5,2);
  Vector cf_grid(2,4,2);

  spectrometer_transfer_matrix( H, srm, cf_grid, f_grid);

  cout << "H:\n" << H << "\n";
  cout << "srm:\n" << srm << "\n";
  cout << "cf_grid:\n" << cf_grid << "\n";
  cout << "f_grid:\n" << f_grid << "\n";
}

void test3()
{
  //Actually a test for Sparse resize function
  cout << "\nTest 3:\n";

  Sparse H(3,3);
  for( Index i=0; i<H.nrows(); i++) {
    H.rw(i,i) = i+1;
  }

  cout << "H initial:\n" << H << "\n";

  H.resize(4,4);

  cout << "H resized:\n" << H << "\n";

  for( Index i=0; i<H.nrows(); i++ ) {
   H.rw(i,i) = 4-i;
  }

  cout << "H refilled:\n" << H << "\n";
}

void test4()
{
  //Test string comparison
  cout << "\nTest 4:\n";

  String s1="test";

  if( s1=="test" ) {
    cout << "ok!\n";
  } else {
    cout << "not equal!\n";
  }
}

void test5()
{
  //Test mixer_transfer_matrix
  cout << "\nTest 5:\n";

  Sparse H(1,1);
  Vector f_mixer(1);
  Vector f_grid(7);
  f_grid[0]=2;
  f_grid[1]=4;
  f_grid[2]=6;
  f_grid[3]=8;
  f_grid[4]=13;
  f_grid[5]=15;
  f_grid[6]=17;
  bool is_upper=true;
  Numeric lo=10;
  Matrix sfrm(2,2);
  sfrm(0,0)=0;
  sfrm(1,0)=20;
  sfrm(0,1)=0;
  sfrm(1,1)=2;

  mixer_transfer_matrix( H, f_mixer, f_grid, is_upper, lo, sfrm);

  cout << "H:\n" << H << "\n";

  cout << "f_mixer:\n" << f_mixer << "\n";
}

void test6()
{
  //Test backend_transfer_matrix
  cout << "\nTest 6:\n";

  Vector f_mixer(11,7,1);
  Vector f_backend(12,5,1);
  /*
  f_backend[0]=1;
  f_backend[1]=2;
  f_backend[2]=4;
  f_backend[3]=5;
  f_backend[4]=7;
  */

  Sparse H(f_backend.nelem() ,f_mixer.nelem());
  Matrix srm(3,2,0.0);
  srm(0,0)=-4;
  srm(2,0)=4;
  srm(1,1)=0.25;

  spectrometer_transfer_matrix( H, srm, f_backend, f_mixer);

  cout << "f_mixer:\n" << f_mixer << "\n";
  cout << "f_backend:\n" << f_backend << "\n";
  cout << "srm:\n" << srm << "\n";
  cout << "H:\n" << H << "\n";
}

void test7()
{
  //Test sensor_integration_vector
  cout << "\nTest 7:\n";
  
  //Calculate a gaussian response
  Matrix srm(17,2);
  for( Index i=0; i<srm.nrows(); i++ ) {
    srm(i,0)=3+i*0.25;
	srm(i,1)=1/(0.5*sqrt(2*PI))*exp(-pow(srm(i,0)-5,2.0)/(2*pow(0.5,2.0)));
  }
  
  Vector h(10), f_grid(1,10,1);

  sensor_integration_vector(h, srm(joker,1), srm(joker,0), f_grid);

  cout << "srm:\n" << srm << "\n";

  cout << "h:\n" << h << "\n";

  h*=f_grid;
  
  cout << "h*g:\n" << h.sum() << "\n";
}

int main()
{
  test1();
//  test2();
//  test3();
//  test4();
//  test5();
//  test6();
//  test7();

  return 0;
}
