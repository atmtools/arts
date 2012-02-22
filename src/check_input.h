/* Copyright (C) 2002-2008
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

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
  ===  File description 
  ===========================================================================*/

/*!
   \file   check_input.h
   \author Patrick Eriksson <patrick.eriksson@chalmers.se>
   \date   2002-04-15 

   This file contains the declaration of functions in check_input.cc.
*/



#ifndef checkinput_h
#define checkinput_h



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "agenda_class.h"
#include "exceptions.h"
#include "matpackVII.h"
#include "mystring.h"
#include "gridded_fields.h"


/*===========================================================================
  === Functions in check_input.cc
  ===========================================================================*/

void chk_if_bool( 
        const String&   x_name,
        const Index&    x );

void chk_if_in_range( 
        const String&   x_name,
        const Index&    x, 
        const Index&    x_low, 
        const Index&    x_high );

void chk_if_increasing( 
        const String&       x_name,
        const ArrayOfIndex& x ); 

void chk_not_negative( 
        const String&    x_name,
        const Numeric&   x );

void chk_if_in_range( 
        const String&    x_name,
        const Numeric&   x, 
        const Numeric&   x_low, 
        const Numeric&   x_high );

void chk_vector_length( 
        const String&      x_name,
        ConstVectorView    x,
        const Index&       l );

void chk_vector_length( 
        const String&      x1_name,
        const String&      x2_name,
        ConstVectorView    x1, 
        ConstVectorView    x2 );

void chk_if_increasing( 
        const String&      x_name,
        ConstVectorView    x );

void chk_if_decreasing( 
        const String&      x_name,
        ConstVectorView    x );

void chk_if_equal(
        const String&   x1_name,
        const String&   x2_name,
        ConstVectorView v1,
        ConstVectorView v2,
        Numeric         margin=1e-6);

void chk_interpolation_grids(const String&   which_interpolation,
                             ConstVectorView old_grid,
                             ConstVectorView new_grid,
                             const Index     order=1,                             
                             const Numeric&  extpolfac=0.5 );

void chk_interpolation_grids(const String&   which_interpolation,
                             ConstVectorView old_grid,
                             const Numeric&  new_grid,
                             const Index     order=1,
                             const Numeric&  extpolfac=0.5 );

void chk_matrix_ncols( 
        const String&      x_name,
        ConstMatrixView    x,
        const Index&       l );

void chk_matrix_nrows( 
        const String&      x_name,
        ConstMatrixView    x,
        const Index&       l );

void chk_atm_grids( 
        const Index&      dim,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid );

void chk_atm_field( 
        const String&     x_name,
        ConstTensor3View  x, 
        const Index&      dim,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid );

void chk_atm_field( 
        const String&   x_name,
        ConstTensor4View  x, 
        const Index&    dim,
        const Index&    nspecies,
        ConstVectorView p_grid,
        ConstVectorView lat_grid,
        ConstVectorView lon_grid );

void chk_atm_surface( 
        const String&     x_name,
        const Matrix&     x, 
        const Index&      dim,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid );

void chk_not_empty( 
        const String&      x_name,
        const Agenda&      x );

void chk_pnd_field_raw_only_in_cloudbox(
        const Index&                 dim,
        const ArrayOfGriddedField3&  pnd_field_raw,
        ConstVectorView              p_grid,
        ConstVectorView              lat_grid,
        ConstVectorView              lon_grid,
        const ArrayOfIndex&          cloudbox_limits);

/*===========================================================================
  === Template Functions for Arrays
  ===========================================================================*/

//! Check if an array contains a value.
/*!
  This makes sure that the array *x* contains the element with
  value *what* exactly once.

  As a bonus, it returns the index of *what* in *x*.

  This template function can be used for arrays of anything, provided
  that the "==" operator is defined.

  \return The index of the thing we looked for.
  \param x_name Name of the array to check
  \param x The array to check
  \param what The value to look for.

  \author Stefan Buehler
  \date   2002-11-28
*/
template <class T>
Index chk_contains( const String&   x_name,
                    const Array<T>& x,
                    const T&        what )
{
  // To generate error messages:
  ostringstream os;

  // To store the positions:
  ArrayOfIndex pos;

  // Find all positions of what in x and store in pos:
  find_all( pos, x, what );

  switch ( pos.nelem() ){

  case 0:
    // Not found.
    os << "The array *" << x_name
      <<  "* must contain the element " << what << ",\n"
      << "but it does not.";
    throw runtime_error( os.str() );
    break;

  case 1:
    // Found once, this is what we want!
    return pos[0];

  default:
    // Found more than once.
    os << "The array *" << x_name
      <<  "* must contain the element " << what << "\n"
      << "exactly once, but it does contain it "
      << pos.nelem() << " times.";
    throw runtime_error( os.str() );
    break;
  }

  return -1;
}

//! Check the size of an array.
/*! 
    Checks the size of an Array. Cloned from
    Patricks similar function for Vector.

    The function throws a runtime_error if the size is not correct.  

    This is a template function that works for any array type.

    \param    x_name   The name of the variable.
    \param    x        A variable of type ArrayOfIndex.
    \param    c        The size to match

    \author Stefan Buehler
    \date   2007-05-18
*/
template <class T>
void chk_size( const String&   x_name,
               const Array<T>& x,
               const Index&    c ) 
{
  if ( x.nelem() != c )
    {
      ostringstream os;
      os << "The array *" << x_name << "*\n"
         << "does not have the right size.\n"
         << "The size should be: " << c << "\n" 
         << "but it is:          " << x.nelem();
      throw runtime_error( os.str() );
    }
}



/*===========================================================================
  === Functions for Tensors
  ===========================================================================*/

void chk_size( const String&    x_name,
               ConstVectorView  x,
               const Index&     c );

void chk_size( const String&    x_name,
               ConstMatrixView  x,
               const Index&     r,
               const Index&     c );

void chk_size( const String&    x_name,
               ConstTensor3View x,
               const Index&     p,
               const Index&     r,
               const Index&     c );

void chk_size( const String&    x_name,
               ConstTensor4View x,
               const Index&     b,
               const Index&     p,
               const Index&     r,
               const Index&     c );

void chk_size( const String&    x_name,
               ConstTensor5View x,
               const Index&     s,
               const Index&     b,
               const Index&     p,
               const Index&     r,
               const Index&     c );

void chk_size( const String&    x_name,
               ConstTensor6View x,
               const Index&     v,
               const Index&     s,
               const Index&     b,
               const Index&     p,
               const Index&     r,
               const Index&     c );

void chk_size( const String&    x_name,
               ConstTensor7View x,
               const Index&     l,
               const Index&     v,
               const Index&     s,
               const Index&     b,
               const Index&     p,
               const Index&     r,
               const Index&     c );

#endif  // checkinput_h
