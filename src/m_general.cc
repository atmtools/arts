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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_general.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-08 

  \brief  Workspace functions of a general and overall character.

  This file is for general functions that do not fit in any other "m_"-file.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "check_input.h"
#include "messages.h"
#include "mystring.h"

#include "math_funcs.h"
#include "make_vector.h"
#include "sensor.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! AntennaSet1D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void AntennaSet1D(
        // WS Output:
              Index&    antenna_dim,
              Vector&   mblock_aa_grid )
{
  out2 << "  Sets the antenna dimensionality to 1.\n";
  out3 << "    antenna_dim = 1\n";
  out3 << "    mblock_aa_grid is set to be an empty vector\n";
  antenna_dim = 1;
  mblock_aa_grid.resize(0);
}



//! AntennaSet2D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void AntennaSet2D(
        // WS Output:
              Index&   antenna_dim,
        // WS Input:
        const Index&   atmosphere_dim )
{
  if( atmosphere_dim != 3 )
    throw runtime_error("Antenna dimensionality 2 is only allowed when the "
                                          "atmospheric dimensionality is 3." );
  out2 << "  Sets the antenna dimensionality to 1.\n";
  out3 << "    antenna_dim = 2\n";
  antenna_dim = 2;
}



//! ArrayOfIndexPrint
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void ArrayOfIndexPrint(
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



//! ArrayOfStringPrint
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void ArrayOfStringPrint(
        // WS Generic Input:
        const ArrayOfString&   x,
        // WS Generic Input Names:
        const String&          x_name )
{
  cout << "  *" << x_name <<"*:\n";
  for( Index i=0; i<x.nelem(); i++ )
    cout << "     " << x[i] << "\n";
}



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



//! IndexPrint
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void IndexPrint(
        // WS Generic Input:
        const Index&   x,
        // WS Generic Input Names:
        const String&  x_name )
{
  cout << "  *" << x_name <<"* = " << x << "\n";;
}



//! MatrixPrint
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-16
*/
void MatrixPrint(
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



//! NumericPrint
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void NumericPrint(
        // WS Generic Input:
        const Numeric&   x,
        // WS Generic Input Names:
        const String&    x_name )
{
  cout << "  *" << x_name <<"* = " << x << "\n";;
}



//! StringPrint
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-18
*/
void StringPrint(
        // WS Generic Input:
        const String&  x,
        // WS Generic Input Names:
        const String&  x_name )
{
  cout << "  *" << x_name <<"* = " << x << "\n";;
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



//! VectorPrint
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-16
*/
void VectorPrint(
        // WS Generic Input:
        const Vector&   x,
        // WS Generic Input Names:
        const String&   x_name )
{
  cout << "  *" << x_name <<"*:\n";
  for( Index i=0; i<x.nelem(); i++ )
    cout << "     " << x[i] << "\n";
}



