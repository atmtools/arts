/*!
  \file   geomag_calc.cc
  \author Nikolay Koulev <nkoulev@uni-bremen.de>
  \date   Mon Jul 28 11:38:22 2003
  
  \brief  Routine for calculating the geomagnetic field.
  
  
*/
#include <iostream>
#include <cmath>
#include "arts.h"
#include "matpackIII.h"
#include "xml_io.h"
#include "legendre.h"


extern const Numeric PI;
extern const Numeric DEG2RAD;

void magfield_nk(  // Output
	    Numeric& B_r, // radial component of the geomagnetic field in [nT].
	    Numeric& B_th, // latitudinal component of the geomagnetic field in [nT].
	    Numeric& B_ph, // longitudinal component of the geomagnetic field in [nT].
	    
	    // Input
	    const Numeric r, // radial distance to the point in [km]
	    const Numeric theta, // geocentric colatitude of the point in [deg]
	    const Numeric phi, // longitude of the point in [deg].
	    // All coordinates - geocentric!

	    const Index Ny // number of elapsed years after an epoch year, J - [0,4]
	    )
  
{
  
  Numeric a; // mean radius of the Earth in [km].
  Numeric Phi = phi * DEG2RAD; // Longitude angle in radian.
  Numeric Theta = PI/2 - theta * DEG2RAD; // Colatitude angle in radian.
    
  // Initializing values of the magnetic field components.
  B_r = 0;
  B_th = 0;
  B_ph = 0;

  Matrix M;
  
  xml_read_from_file ("geomag_coefficients.xml", M);
  
  // M(i,0) and M(i,1) - the vectors with the values of the first and second coefficients 
  // of the IGRF model.
  // M(i,2) and M(i,3) - the vectors with the values of the anual rate of change of the 
  // first and second coefficient of of the IGRF model. 
  
  
  // Loop over the degree number l of the Legendre polynommes.
  for (Index l = 1; l <= 10; l++)
    {
      // Loop over the order number m of the Legendre polynommes.
      for (Index m = 0; m <= l; m++)
	{
	  
	  // Relating the row index in M to the coresponding 
	  // degree number n and order number l.
 	  Index	  j = l * (l + 1) / 2 + m - 1; 

	  // Calculating the associated Schmidt quasi-normalized Legendre 
	  // polynomial for a degree number l and order number m.
	  Numeric P_lm = 
	    legendre_poly_norm_schmidt (l, m, cos(Theta));

	  // Calculating the derivative of the associated Schmidt quasi-normalized 
	  // Legendre polynomial for a degree number l and order number m.
	  Numeric dP_lm = 
	    legendre_poly_norm_schmidt_deriv (l, m, cos(Theta));
	  
	  // Calculating the radial (upward) component of the magnetic field.
	  B_r += pow(l + 2, a / r) * (l + 1) * 
	    ((M(j,0) + Ny * M(j,2)) * cos(m * Phi) 
	     + (M(j,1) + Ny * M(j,3)) * sin(m * Phi)) 
	       * P_lm;
	  
	  // Calculating the latitudinal (southward) component of the magnetic field. 
	  B_th += pow(l + 2, a / r) * 
	    ((M(j,0) + Ny * M(j,2)) * cos(m * Phi) 
	     + (M(j,1) + Ny * M(j,3)) * sin(m * Phi)) *
	      dP_lm * sin(Theta);
      
      
          // Calculating the longitudinal (eastward) component of the magnetic field.
          B_ph += pow(l + 2, a / r) * m *
	   ((M(j,0) + Ny * M(j,2)) * sin(m * Phi) 
	    - (M(j,1) + Ny * M(j,3)) * cos(m * Phi)) *
	     P_lm / sin(Theta);
      
	  
	}
    }
  
}
