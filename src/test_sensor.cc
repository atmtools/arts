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
  //Test antenna_transfer_matrix with ArrayOfArrayOfMatrix
  cout << "\nTest 1:\n";

  Vector ant_za(-0.03,2,0.06);
  cout << "ant_za:\n" <<ant_za<<endl;
  Vector f(3.175e11,8,1e9);
  cout << "f:\n"<<f<<endl;
  Vector m_za(-0.1,9,0.025);
  cout << "m_za:\n" << m_za << endl;
  Vector a_grid(-0.05,11,0.01);
  cout << "a_grid:\n" << a_grid << endl;
  ArrayOfMatrix adiag1, adiag2, adiag3;
  ArrayOfArrayOfMatrix aadiag;
  Index n_pol = 2;

  Matrix diag(a_grid.nelem(), 2);
  diag(joker,0) = a_grid;
  antenna_diagram_gaussian(diag(joker,Range(0,2,1)),0.09);

  Matrix diag1(a_grid.nelem(),6);
  diag1(joker,0) = a_grid;
  antenna_diagram_gaussian(diag1(joker,Range(0,2,1)),1.8);
  antenna_diagram_gaussian(diag1(joker,Range(0,2,2)),1.9);
  antenna_diagram_gaussian(diag1(joker,Range(0,2,3)),2.0);
  antenna_diagram_gaussian(diag1(joker,Range(0,2,4)),2.1);
  antenna_diagram_gaussian(diag1(joker,Range(0,2,5)),2.2);

  Matrix diag2(a_grid.nelem(),6);
  diag2(joker,0) = a_grid;
  antenna_diagram_gaussian(diag2(joker,Range(0,2,1)),2.3);
  antenna_diagram_gaussian(diag2(joker,Range(0,2,2)),2.4);
  antenna_diagram_gaussian(diag2(joker,Range(0,2,3)),2.5);
  antenna_diagram_gaussian(diag2(joker,Range(0,2,4)),2.6);
  antenna_diagram_gaussian(diag2(joker,Range(0,2,5)),2.7);

  Matrix diag3(a_grid.nelem(),6);
  diag3(joker,0) = a_grid;
  antenna_diagram_gaussian(diag3(joker,Range(0,2,1)),2.8);
  antenna_diagram_gaussian(diag3(joker,Range(0,2,2)),2.9);
  antenna_diagram_gaussian(diag3(joker,Range(0,2,3)),3.0);
  antenna_diagram_gaussian(diag3(joker,Range(0,2,4)),3.1);
  antenna_diagram_gaussian(diag3(joker,Range(0,2,5)),3.2);

  adiag1.push_back(diag);
  //adiag1.push_back(diag1);
  adiag2.push_back(diag2);
  adiag3.push_back(diag3);
  aadiag.push_back(adiag1);
  //aadiag.push_back(adiag1);
  //aadiag.push_back(adiag2);
  //aadiag.push_back(adiag3);

  cout << " Initialising H...";
  Sparse H(f.nelem()*ant_za.nelem()*n_pol,m_za.nelem()*f.nelem()*n_pol);
  cout << "done.\n";

  cout << " Calculating antenna transfer matrix H...";
  antenna_transfer_matrix(H,m_za,aadiag,f,ant_za,n_pol);
  cout << "done.\n";

  cout << "H:["<<H.nrows()<<","<<H.ncols()<<"]:\n"<<H<<"\n";

  cout << " Initialising identity matrix I...";
  Sparse I(H.ncols(),H.ncols());
  for (Index i=0;i<I.ncols(); i++)
    I.rw(i,i)=1;
  cout << "done.\n";

  cout << " Initialising R...";
  Sparse R;
  cout << "done.\n";

  cout << " Resizing R...";
  R.resize(H.nrows(),H.ncols());
  cout << "done.\n";

  cout << "R:["<<R.nrows()<<","<<R.ncols()<<"]:\n"<<R<<"\n";
  
  cout << " Multiplying H and I, store product in R...";
  mult(R,H,I);
  cout << "done.\n";

  cout << "R:["<<R.nrows()<<","<<R.ncols()<<"]:\n"<<R<<"\n";

}

void test2()
{
  //Test spectrometer_transfer_matrix
  cout << "\nTest 2:\n";

  Vector sensor_f(10,5,10);
  Vector ch_f(15,3,15);
  Vector a_grid(-5,11,1);
  Index n_za = 3;
  Index n_pol = 2;
  ArrayOfMatrix aresp;

  Matrix ch_response1(11,4);
  ch_response1(joker,0) = a_grid;
  antenna_diagram_gaussian(ch_response1(joker,Range(0,2,1)),0.09);
  antenna_diagram_gaussian(ch_response1(joker,Range(0,2,2)),0.1);
  antenna_diagram_gaussian(ch_response1(joker,Range(0,2,3)),0.11);

  Matrix ch_response2(11,4);
  ch_response2(joker,0) = a_grid;
  antenna_diagram_gaussian(ch_response2(joker,Range(0,2,1)),0.14);
  antenna_diagram_gaussian(ch_response2(joker,Range(0,2,2)),0.15);
  antenna_diagram_gaussian(ch_response2(joker,Range(0,2,3)),0.16);

  aresp.push_back(ch_response1);
  aresp.push_back(ch_response2);

  Sparse H(ch_f.nelem()*n_za*n_pol,sensor_f.nelem()*n_za*n_pol);

  spectrometer_transfer_matrix(H,aresp,ch_f,sensor_f,n_za,n_pol);

  cout << "H:\n" << H << "\n";
//  cout << "ch_response1:\n" << ch_response1 << "\n";
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
  Numeric lo=10;
  Matrix filter(2,2);
  filter(0,0)=0;
  filter(1,0)=20;
  filter(0,1)=0;
  filter(1,1)=2;

  mixer_transfer_matrix( H, f_mixer, f_grid, lo, filter, 2, 2);

  cout << "H:\n" << H << "\n";

  cout << "f_mixer:\n" << f_mixer << "\n";
}

void test6()
{
  // Test spectrometer_transfer_matrix with one polarisation and
  // one viewing angle
  cout << "\nTest 6:\n";

  Vector sensor_f(10,5,10);
  Vector ch_f(15,3,15);
  Vector a_grid(-5,11,1);
  Index n_za = 1;
  Index n_pol = 1;
  ArrayOfMatrix aresp;

  Matrix ch_response(11,2);
  ch_response(joker,0) = a_grid;
//  antenna_diagram_gaussian(ch_response(joker,Range(0,2,2)),1.9);
  antenna_diagram_gaussian(ch_response(joker,Range(0,2,1)),2.0);
//  antenna_diagram_gaussian(ch_response(joker,Range(0,2,4)),2.1);

  aresp.push_back(ch_response);

  Sparse H(ch_f.nelem()*n_za*n_pol,sensor_f.nelem()*n_za*n_pol);

  spectrometer_transfer_matrix(H,aresp,ch_f,sensor_f,n_za,n_pol);

  cout << "H["<<H.nrows()<<","<<H.ncols()<<"]:\n" << H << "\n";
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
//  test1();
//  test2();
//  test3();
//  test4();
  test5();
//  test6();
//  test7();

  return 0;
}
