/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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



/*****************************************************************************
 ***  File description 
 *****************************************************************************/

/*!
  \file   m_general.cc
  \author Patrick Eriksson
  \date   2002-05-08 

  \brief  Workspace functions of a general and overall character.

  This file is for general functions that do not fit in any other "m_"-file.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*****************************************************************************
 *** External declarations
 *****************************************************************************/

#include "array.h"
#include "arts.h"
#include "check_input.h"
#include "messages.h"
#include "mystring.h"



/*****************************************************************************
 *** The functions (in alphabetical order)
 *****************************************************************************/

//! Exit
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void Exit()
{
  out1 << "  Forced exit.\n";
  exit(0);
}



//! PrintArrayOfIndex
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void PrintArrayOfIndex(
        // WS Generic Input:
        const ArrayOfIndex&   x,
        // WS Generic Input Names:
        const String&         x_name )
{
  cout << "  *" << x_name <<"* =";
  for( Index i=0; i<x.nelem(); i++ )
    cout << " " << x[i];
  cout << "\n";
}



//! PrintArrayOfString
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void PrintArrayOfString(
        // WS Generic Input:
        const ArrayOfString&   x,
        // WS Generic Input Names:
        const String&          x_name )
{
  cout << "  *" << x_name <<"*:\n";
  for( Index i=0; i<x.nelem(); i++ )
    cout << "     " << x[i] << "\n";
}



//! PrintIndex
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void PrintIndex(
        // WS Generic Input:
        const Index&   x,
        // WS Generic Input Names:
        const String&  x_name )
{
  cout << "  *" << x_name <<"* = " << x << "\n";;
}



//! PrintMatrix
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-16
*/
void PrintMatrix(
        // WS Generic Input:
        const Matrix&   x,
        // WS Generic Input Names:
        const String&   x_name )
{
  cout << "  *" << x_name <<"*:\n";
  for( Index i=0; i<x.nrows(); i++ )
    {
      cout << "     ";
      for( Index j=0; j<x.ncols(); j++ )
        cout << x(i,j) << " ";
      cout << "\n";
    }
}



//! PrintNumeric
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void PrintNumeric(
        // WS Generic Input:
        const Numeric&   x,
        // WS Generic Input Names:
        const String&    x_name )
{
  cout << "  *" << x_name <<"* = " << x << "\n";;
}



//! PrintString
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void PrintString(
        // WS Generic Input:
        const String&  x,
        // WS Generic Input Names:
        const String&  x_name )
{
  cout << "  *" << x_name <<"* = " << x << "\n";;
}



//! PrintVector
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-16
*/
void PrintVector(
        // WS Generic Input:
        const Vector&   x,
        // WS Generic Input Names:
        const String&   x_name )
{
  cout << "  *" << x_name <<"*:\n";
  for( Index i=0; i<x.nelem(); i++ )
    cout << "     " << x[i] << "\n";
}


//! SetAtmosphere1D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void SetAtmosphere1D(
        // WS Output:
              Index&    atmosphere_dim,
              Vector&   lat_grid,
              Vector&   lon_grid )
{
  atmosphere_dim = 1;
  lat_grid.resize(0);
  lon_grid.resize(0);
}



//! SetAtmosphere2D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void SetAtmosphere2D(
        // WS Output:
              Index&    atmosphere_dim,
              Vector&   lon_grid,
              Numeric&  latitude_1d,
              Numeric&  azimuth_angle_1d )
{
  atmosphere_dim = 2;
  lon_grid.resize(0);
  latitude_1d = -999;
  azimuth_angle_1d = -999;
}



//! SetAtmosphere3D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void SetAtmosphere3D(
        // WS Output:
              Index&    atmosphere_dim,
              Numeric&  latitude_1d,
              Numeric&  azimuth_angle_1d )
{
  atmosphere_dim = 3;
  latitude_1d = -999;
  azimuth_angle_1d = -999;
}



//! Test
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-05-15
*/
void Test( )
{
  // This function can be used to test stuff.

}
