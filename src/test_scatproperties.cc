/*!
  \file   test_scatproperties.cc
  \author Sreerekha T.R. 
  \date   Mon May 13 11:36:11 2002
*/

#include <iostream>
#include "scatproperties.h"
#include "math_funcs.h"

extern const Numeric DEG2RAD;
int main(void)
{
  Vector amp(8);
  Vector za(18),aa(36);
  Matrix l_cross(4,4);
  Numeric freq =  5.01*1e+11; 
  // elements of amplitude matrix
  amp[0] = 1;
  amp[1] = 1;
  amp[2] = 1;
  amp[3] = 1;
  amp[4] = 1;
  amp[5] = 1;
  amp[6] = 1;
  amp[7] = 1;
  
  cout << "\n Extinction Matrix  test \n";
  cout << "Amplitude matrix: \n";
  cout << amp;  
  
  cout << "\n";
  za[0]=0;
  aa[0]=0;
  amp2ext(l_cross,amp,freq);
  cout << "Extinction Matrix: \n";
  for( Index i = 0; i<4; i++)
    {
      for (Index j = 0; j<4; j++)
	//ext_matrix(l_cross,amp,freq);
	cout << " " << l_cross(i,j);
    }
      cout << "\n";
    
      /*testing the new code for matrix integration*/
      Matrix D(181,181);
      Matrix D1(181,181);
      Vector za_g(181);
      Vector aa_g(181);
      for (Index i = 0;i < za_g.nelem(); ++i)
	{
	  za_g[i] = i;
	  
	  for (Index j = 0;j < aa_g.nelem(); ++j)
	    {
	      aa_g[j]=j;
	      D(i,j) =  sin(aa_g[j]*DEG2RAD);
	      
	      D(i,j) *= DEG2RAD ;
	    }
	}
      AngIntegrate_trapezoid(D, za_g, aa_g);
      //cout<<"zaa"<<"  "<<D(0,30)<<"\n";
      return 0;
}
