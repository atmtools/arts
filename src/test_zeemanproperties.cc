/*!
  \file   test_zeemanproperties.cc
  \author Nikolay Koulev <nkoulev@uni-bremen.de>
  \date   Fri Jan 10 16:46:32 2003
  
  \brief  This is the file for testing the Zeeman effect related routines in ARTS.
  
  
*/
#include <cmath>
#include <iostream>
#include <sstream>
#include "math_funcs.h"
#include "matpackIII.h"
#include "zeemanproperties.h"
#include "xml_io.h"

using namespace std;

int main(void)
{//Input & Output
  Vector f_grid;         // Frequency grid.

  //Output
  Tensor3 ext_mat_zee;   // Tensor3 of the Extinction Matrix for the specified frequency grid. 
  Matrix abs_vec_zee;    // Matrix of the Absorption Vector for the specified frequency grid.
  Matrix xi_mat;         // Matrix of the relative intensities of the Zeeman components.
  Matrix f_z_mat;        // Matrix of the frequencies of the Zeeman components.
 
  //Input
  Matrix N;          // Matrix variable for the content of the zeeman_intensity_coeff.xml. 
  Numeric N_r;       // Rotational quantum number indentification of the Zeeman line.
  Numeric BN_r;      // Number of the split components for a given polarization.
  Numeric AN_r;      // Maximum value of the magnetic quantum number for a given polarization.
  Numeric f_c;       // Central frequency of the unsplit line [GHz].
  Numeric a1;        // First intensity coeficient of the unsplit line.
  Numeric a2;        // Second intensity coeficient of the unsplit line.
  Numeric a3;        // Third intensity coeficient of the unsplit line.





  // Reading the content of zeeman_intensity_coeff.xml. 
  xml_read_from_file ("zeeman_intensity_coeff.xml", N);
  for (Index i=N.nrows()-9;i<N.nrows()-8;i++) // Dummy case through the 5- line.
    {
      N_r = N(i,0); 
      a1 = N(i,2);
      a2 = N(i,3);
      f_c = N(i,1)*1e9; // Reading in and converitng the central frequency to SI unit [Hz].
      cout.precision(10);
      cout << "central frequency = " << f_c << " Hz"<<endl;
      
      // Frequency grid range around the central frequency of the unsplit line.
      Numeric range; 
      range = 6e+6; // dummy value of 6 MHz around the central line
      
      // Starting point of the frequency grid.
      Numeric f_start;
      f_start = f_c - range;
      
      // Number of the frequency grid steps.
      int n_f;
      n_f = 1000; // dummy value
      
      // Size of of the frequency grid step.
      Numeric f_step;
      f_step = 2*range/n_f;
      
      // Initializing the frequency grid.
      Vector f_grid(f_start, n_f, f_step);

      // Calling the function defined in zeemanproperties.cc .
      Zeeman(f_grid , ext_mat_zee, abs_vec_zee, xi_mat, f_z_mat, N_r, BN_r, AN_r, f_c, a1, a2, a3);
 

      // Defining a base string for the names of the output files.
      ostringstream os;
      if (N_r>0)
	{
	  os << abs(N_r) << "+" << "_" << n_f << ".xml"; 
	}
      else if (N_r<0)
	{
	  os << abs(N_r) << "-" << "_" << n_f << ".xml"; 
	}
      else
	{
	  cout << "Rotational number is wrong!" << endl;
	}
      string basestring=os.str();
      
      // Writing out the output files for the respective variables.
      xml_write_to_file(string("f_grid")+basestring,f_grid);
      xml_write_to_file(string("ext_mat_zee")+basestring,ext_mat_zee);
      xml_write_to_file(string("abs_vec_zee")+basestring,abs_vec_zee);
      xml_write_to_file(string("f_z_mat")+basestring,f_z_mat);
      xml_write_to_file(string("xi_mat")+basestring,xi_mat);



    }    
  
  
}
 
