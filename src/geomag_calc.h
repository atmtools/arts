/*!
  \file   geomag_calc.h
  \author Nikolay Koulev <nkoulev@uni-bremen.de>
  \date   Mon Jul 28 11:35:24 2003
  
 \brief The header file for the functions in geomag_calc.cc


*/


#ifndef geomag_calc_h
#define geomag_calc_h

#include "arts.h"
#include "matpackI.h"
#include "matpackIV.h"

void magfield_nk(  // Output
	    Numeric& B_r, // radial component of the geomagnetic field
	    Numeric& B_th, // colatitudinal component of the geomagnetic field
	    Numeric& B_ph, // longitudinal component of the geomagnetic field
	    
	    // Input
	    const Numeric r, // radial distance to the point 
	    const Numeric theta, // geocentric colatitude of the point
	    const Numeric phi, // longitude of the point
	    // All coordinates - geocentric!

	    const Index Ny // number of elapsed years after an epoch year, J - [0,4]
	    );
#endif 
