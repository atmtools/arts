/*!
  \file   test_zeeman.cc
  \author Nikolay Koulev <nkoulev@uni-bremen.de>
  \date   Fri Jan 10 16:46:32 2003
  
  \brief  This is the file for testing the Zeeman effect related routines in ARTS.
  
  
*/
#include <cmath>
#include <iostream>
#include "math_funcs.h"
#include "matpackIII.h"
#include "zeemanproperties.h"
#include "xml_io.h"

int main(void)
{
  Tensor3 ext_mat_zee;
  Matrix abs_vec_zee;
  Vector f_grid(100000);
  Matrix N;
  Matrix xi_mat;
  Matrix f_z_mat;
   

  for (Index i=0; i<f_grid.nelem(); i++)
    { 
      f_grid[i]=/*!118.749343*/61.149570 + 0.00000002*Numeric(i); //Numeric(i) - converts i in Numeric 
    }
  Zeeman(f_grid, ext_mat_zee, abs_vec_zee, xi_mat, f_z_mat);

 
  xml_write_to_file("f_grid.xml",f_grid);
  xml_write_to_file("ext_mat_zee.xml", ext_mat_zee);
  xml_write_to_file("abs_vec_zee.xml", abs_vec_zee);

  xml_write_to_file("f_z_mat.xml",f_z_mat);
  xml_write_to_file("xi_mat.xml",xi_mat);
  
}
