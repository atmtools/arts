/*!
  \file   test_zeeman.cc
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
{
  Tensor3 ext_mat_zee;
  Matrix abs_vec_zee;
  Vector f_grid(1000);
  Matrix N;
  Matrix xi_mat;
  Matrix f_z_mat;
  Numeric N_r;
  Numeric BN_r;
  Numeric AN_r;
  Numeric f_c;
  Numeric a1;
  Numeric a2;


  xml_read_from_file ("zeeman_intensity_coeff.xml", N);
  for (int i=N.nrows()-9;i<N.nrows()-8;i++)
    {
      N_r = N(i,0); 
      a1 = N(i,2);
      a2 = N(i,3);
      f_c = N(i,1);
      cout.precision(10);
      cout << "central frequency = " << f_c << " GHz"<<endl;
      

      for (Index j=0; j<f_grid.nelem(); j++)
	{ 
	  f_grid[j]= f_c - 0.006 + (0.012/f_grid.nelem())*Numeric(j); //Numeric(j) - converts j in Numeric 
	  //f_grid[j]= f_c - 0.002 + 0.0000004*Numeric(j); //Numeric(j) - converts j in Numeric 
	}
      
      Zeeman(f_grid, ext_mat_zee, abs_vec_zee, xi_mat, f_z_mat, N_r, AN_r, BN_r, f_c, a1, a2);
      
      ostringstream os1;
      ostringstream os2;
      ostringstream os3;
      ostringstream os4;
      ostringstream os5;
      
      if (N_r<0)
	{
	  os1 << "f_grid" << "_" << abs(N_r) << "-" << "_" << f_grid.nelem() << ".xml";
	  string namestring1 = os1.str();
	  xml_write_to_file(namestring1,f_grid);
	  
	  
	  os2 << "ext_mat_zee" << "_" << abs(N_r) << "-" << "_" << f_grid.nelem() << ".xml";
	  string namestring2 = os2.str();
	  xml_write_to_file(namestring2,ext_mat_zee);
	  
	  os3 << "abs_vec_zee" << "_" << abs(N_r) << "-" << "_" << f_grid.nelem() << ".xml";
	  string namestring3 = os3.str();
	  xml_write_to_file(namestring3,abs_vec_zee);
	  
	  os4 << "f_z_mat" << "_" << abs(N_r) << "-" << "_" << f_grid.nelem() << ".xml";
	  string namestring4 = os4.str();
	  xml_write_to_file(namestring4,f_z_mat);
	  
	  os5 << "xi_mat" << "_" << abs(N_r) << "-" << "_" << f_grid.nelem() << ".xml";
	  string namestring5 = os5.str();
	  xml_write_to_file(namestring5,xi_mat);
	}
      else if (N_r>0)
	{
	  os1 << "f_grid" << "_" << abs(N_r) << "+" << "_" << f_grid.nelem() << ".xml";
	  string namestring1 = os1.str();
	  xml_write_to_file(namestring1,f_grid);
	  
	  
	  os2 << "ext_mat_zee" << "_" << abs(N_r) << "+" << "_" << f_grid.nelem() << ".xml";
	  string namestring2 = os2.str();
	  xml_write_to_file(namestring2,ext_mat_zee);
	  
	  os3 << "abs_vec_zee" << "_" << abs(N_r) << "+" << "_" << f_grid.nelem() << ".xml";
	  string namestring3 = os3.str();
	  xml_write_to_file(namestring3,abs_vec_zee);
	  
	  os4 << "f_z_mat" << "_" << abs(N_r) << "+" <<  "_" << f_grid.nelem() <<".xml";
	  string namestring4 = os4.str();
	  xml_write_to_file(namestring4,f_z_mat);
	  
	  os5 << "xi_mat" << "_" << abs(N_r) << "+" << "_" << f_grid.nelem() << ".xml";
	  string namestring5 = os5.str();
	  xml_write_to_file(namestring5,xi_mat);
	}
      
    }    
  //  xml_write_to_file("f_grid.xml",f_grid);
  //    xml_write_to_file("ext_mat_zee.xml", ext_mat_zee);
  //    xml_write_to_file("abs_vec_zee.xml", abs_vec_zee);
  
  //    xml_write_to_file("f_z_mat.xml",f_z_mat);
  //    xml_write_to_file("xi_mat.xml",xi_mat);
  
}
