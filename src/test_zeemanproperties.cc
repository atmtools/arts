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
  Tensor3 extZ;
  Matrix absvZ;
  Vector f_grid(1000);
 
  
  for (Index i=0; i<f_grid.nelem(); i++)
    { 
      f_grid[i]=117.5 + 1000.*Numeric(i); //Numeric(i) - converts i in Numeric
    }
  Zeeman(f_grid, extZ, absvZ );

 

  xml_write_to_file("extZ.xml", extZ);
  xml_write_to_file("absvZ.xml", absvZ);
}
