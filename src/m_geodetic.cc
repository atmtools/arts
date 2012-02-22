/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
                            
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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_geodetic.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2012-02-06

  \brief  Workspace functions of geodetic character.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/




/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "check_input.h"
#include "geodetic.h"
#include "matpackI.h"
#include "messages.h"


extern const Numeric DEG2RAD;
extern const Numeric EARTH_RADIUS;


// Ref. 1:
// Seidelmann, P. Kenneth; Archinal, B. A.; A'hearn, M. F. et al (2007).
// "Report of the IAU/IAG Working Group on cartographic coordinates and
// rotational elements: 2006". Celestial Mechanics and Dynamical Astronomy 98
// (3): 155–180. Bibcode 2007CeMDA..98..155S. doi:10.1007/s10569-007-9072-y


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidEarth(            
         Vector&    refellipsoid,
   const String&    model,
   const Verbosity& )
{
  refellipsoid.resize(2);

  if( model == "Sphere" )
    {
      refellipsoid[0] = EARTH_RADIUS;
      refellipsoid[1] = 0;
    }

  else if( model == "WGS84" )
    { // Values taken from atmlab's ellipsoidmodels.m
      refellipsoid[0] = 6378137; 
      refellipsoid[1] = 0.081819190842621;
    }

  else
    throw runtime_error( "Unknown selection for input argument *model*." );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidForAzimuth( 
         Vector&     refellipsoid,
   const Numeric&    latitude,
   const Numeric&    azimuth,
   const Verbosity& )
{
  if( refellipsoid.nelem() != 2 )
    throw runtime_error( "Input *refellispoid must be a vector of length 2*." );

  if( refellipsoid[1] > 0 )
    {
      const Numeric e2  = refellipsoid[1] * refellipsoid[1];
      const Numeric a   = 1 - e2 * pow( sin( DEG2RAD*latitude ), 2.0 );

      const Numeric rn = 1 / sqrt( a );
      const Numeric rm = ( 1 - e2 ) * ( rn / a );

      const Numeric v = DEG2RAD * azimuth;

      refellipsoid[0] = refellipsoid[0] / ( pow(cos(v),2.0)/rm +
                                            pow(sin(v),2.0)/rn ); 
      refellipsoid[1] = 0;
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidJupiter(            
         Vector&    refellipsoid,
   const String&    model,
   const Verbosity& )
{
  refellipsoid.resize(2);

  if( model == "Sphere" )
    { 
      refellipsoid[0] = 69911e3;   // From Ref. 1 (see above)
      refellipsoid[1] = 0;
    }

  else if( model == "Ellipsoid" )
    {
      refellipsoid[0] = 71492e3;   // From Ref. 1
      refellipsoid[1] = 0.3543;    // Based on Ref. 1
    }

  else
    throw runtime_error( "Unknown selection for input argument *model*." );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidMars(            
         Vector&    refellipsoid,
   const String&    model,
   const Verbosity& )
{
  refellipsoid.resize(2);

  if( model == "Sphere" )
    { 
      refellipsoid[0] = 3389.5e3;   // From Ref. 1 (see above)
      refellipsoid[1] = 0;
    }

  else if( model == "Ellipsoid" )
    {
      refellipsoid[0] = 3396.19e3;   // From Ref. 1
      refellipsoid[1] = 0.1083;      // Based on Ref. 1
    }

  else
    throw runtime_error( "Unknown selection for input argument *model*." );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidMoon(            
         Vector&    refellipsoid,
   const String&    model,
   const Verbosity& )
{
  refellipsoid.resize(2);

  if( model == "Sphere" )
    { 
      refellipsoid[0] = 1737.4e3;  // From Ref. 1 (see above)
      refellipsoid[1] = 0;
    }

  else if( model == "Ellipsoid" )
    { // Values taken from Wikipedia, with reference to:
      // Williams, Dr. David R. (2 February 2006). "Moon Fact Sheet". 
      // NASA (National Space Science Data Center). Retrieved 31 December 2008.
      refellipsoid[0] = 1738.14e3; 
      refellipsoid[1] = 0.0500;
    }

  else
    throw runtime_error( "Unknown selection for input argument *model*." );
}




/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidOrbitPlane( 
         Vector&     refellipsoid,
   const Numeric&    orbitinc,
   const Verbosity& )
{
  if( refellipsoid.nelem() != 2 )
    throw runtime_error( "Input *refellispoid must be a vector of length 2*." );
  chk_if_in_range( "orbitinc", orbitinc,  0, 180 );    

  // Radius at maximum latitude
  const Numeric rp = refell2r( refellipsoid, orbitinc );

  // New eccentricity
  refellipsoid[1] = sqrt( 1 - pow( rp/refellipsoid[0], 2.0 ) );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidSet(            
         Vector&    refellipsoid,
   const Numeric&   re,
   const Numeric&   e,
   const Verbosity& )
{
  refellipsoid.resize(2);

  refellipsoid[0] = re;
  refellipsoid[1] = e;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidVenus(            
         Vector&    refellipsoid,
   const String&    model,
   const Verbosity& )
{
  refellipsoid.resize(2);

  if( model == "Sphere" )
    { 
      refellipsoid[0] = 6051.8e3;   // From Ref. 1 (see above)
      refellipsoid[1] = 0;
    }

  else
    throw runtime_error( "Unknown selection for input argument *model*." );
}


