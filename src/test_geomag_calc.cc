/*!
  \file   test_geomag_calc.cc
  \author Nikolay Koulev <nkoulev@uni-bremen.de>
  \date   Thu Aug 21 12:12:36 2003
  
  \brief  Test file of the geomagnetic calculations.
  
  
*/
#include <cmath>
#include <iostream>
#include <sstream>
#include "exceptions.h"
#include "messages.h"
#include "math_funcs.h"
#include "matpackIII.h"
#include "geomag_calc.h"
#include "xml_io.h"

extern const Numeric EARTH_RADIUS;

int main(void) 
{// Output
  
  Numeric B_r; // radial component of the geomagnetic field in [nT].
  Numeric B_th; // latitudinal component of the geomagnetic field in [nT].
  Numeric B_ph; // longitudinal component of the geomagnetic field in [nT].
  Numeric B_tot; // Absolute value of the magnetic field in [nT].

  // Input
  Numeric z; // altitutde in [km]
  Numeric theta; // geocentric colatitude of the point
  Numeric phi; // longitude of the point
  // All coordinates - geocentric!
  
  Index Ny; // number of elapsed years after an epoch year, J - [0,4]
  
  
  
  extern Messages messages;

    {
      // Reporting was not specified, set default. (Only the
      // important stuff to the screen, everything to the file.)
      messages.screen = 2;
      messages.file   = 0;
    }


  // Feed in altitutde above the mean radius of the Earth in [km].
  cout << "Altitude in km" << endl;
  cin  >> z;
  
  // Feed in latitute in degrees.
  cout << "Latitude in degrees" << endl;
  cin  >> theta;
  
  // Feed in longitude in degrees.
  cout << "Longitude in degrees" << endl;
  cin >> phi ;
  
  // Feed in number of elapsed years after the epoch year.
  cout <<  "Ny" <<  endl;
  cin >>  Ny;


  // Defining the geocetric radius to the point.
  const Numeric r = EARTH_RADIUS + z;

  try
    {
      magfield_nk(B_r, B_th, B_ph, r, theta, phi, Ny);
      
      // Calculating of the total field.
      B_tot = sqrt(B_r * B_r + B_th * B_th + B_ph * B_ph);

    }
  catch (runtime_error e)
    {
      cerr << e.what ();
      exit(1);
    }

  // Output of the radial component of the geomagnetic field in [nT].
  cout << "B_r = " << B_r << " nT" << endl;
  
  // Output of the latitudinal component of the geomagnetic field in [nT].
  cout << "B_th = " << B_th << " nT" << endl;
  
  // Output of the longitudinal component of the geomagnetic field in [nT].
  cout << "B_ph = " << B_ph << " nT" << endl;
  
  // Output of the total geomagnetic field in [nT].
  cout << "B_tot = " <<  B_tot << " nT" << endl;
  
  
  
  
  
  
  
  
  
}
