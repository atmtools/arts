/* Copyright (C) 2002 Patrick Eriksson <patrick@rss.chalmers.se>
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
  \file   check_input.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date 2002-04-15 

  General functions to check the size and logic of input to functions.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include "check_input.h"
#include "array.h"
#include "logic.h"



/*===========================================================================
  === Functions for Index
  ===========================================================================*/

//! chk_if_bool 
/*! 
    Checks that a variable of type Index has the value 0 or 1.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Index.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_bool( 
        const String&   x_name,
        const Index&    x )
{
  if ( !is_bool(x) )
    {
      ostringstream os;
      os << "The variable *" << x_name <<  "* must be a boolean (0 or 1).\n" 
         << "The present value of *"<< x_name <<  "* is " << x << ".";
      throw runtime_error( os.str() );
    }
}



//! chk_if_in_range
/*! 
    Checks that a variable of type Index has a value inside the specified
    range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Index.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_in_range( 
        const String&   x_name,
        const Index&    x, 
        const Index&    x_low, 
        const Index&    x_high )
{
  if ( (x<x_low) || (x>x_high) )
    {
      ostringstream os;
      os << "The variable *" << x_name <<  "* must fulfill:\n"
         << "   " << x_low << " <= " << x_name << " <= " << x_high << "\n" 
         << "The present value of *"<< x_name <<  "* is " << x << ".";
      throw runtime_error( os.str() );
    }
}


/*===========================================================================
  === Functions for Numeric
  ===========================================================================*/

//! chk_if_over_0 
/*! 
    Checks that a variable of type Numeric is 0 or is positive.
    range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_over_0( 
        const String&    x_name,
        const Numeric&   x ) 
{
  if ( x <= 0 )
    {
      ostringstream os;
      os << "The variable *" << x_name <<  "* must exceed 0.\n"
         << "The present value of *"<< x_name <<  "* is " << x << ".";
      throw runtime_error( os.str() );
    }
}



//! chk_if_in_range
/*! 
    Checks that a variable of type Numeric has a value inside the specified
    range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_in_range( 
        const String&    x_name,
        const Numeric&   x, 
        const Numeric&   x_low, 
        const Numeric&   x_high )
{
  if ( (x<x_low) || (x>x_high) )
    {
      ostringstream os;
      os << "The variable *" << x_name <<  "* must fulfill:\n"
         << "   " << x_low << " <= " << x_name << " <= " << x_high << "\n" 
         << "The present value of *"<< x_name <<  "* is " << x << ".";
      throw runtime_error( os.str() );
    }
}



/*===========================================================================
  === Functions for Vector
  ===========================================================================*/

//! chk_vector_length
/*! 
    Checks that a vector has the specified length.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.
    \param    l        The expected length of x.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_vector_length( 
        const String&      x_name,
        ConstVectorView    x,
        const Index&       l ) 
{
  if ( x.nelem() != l )
    {
      ostringstream os;
      os << "The vector *" << x_name <<  "* must have the length " << l 
         << ".\n" 
         << "The present length of *"<< x_name <<  "* is " << x.nelem() << ".";
      throw runtime_error( os.str() );
    }
}



//! chk_vector_length
/*! 
    Checks if two vectors have the same length.

    The function gives an error message if this is not the case.

    \param    x1_name   The name of the first vector
    \param    x2_name   The name of the second vector
    \param    x1        The first vector.
    \param    x2        The second vector.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_vector_length( 
        const String&      x1_name,
        const String&      x2_name,
        ConstVectorView    x1, 
        ConstVectorView    x2 ) 
{
  if ( x1.nelem() != x2.nelem() )
    {
      ostringstream os;
      os << "The vectors *" << x1_name <<  "* and *" << x1_name 
         <<  "* must have the same length.\n"
         << "The length of *"<< x1_name <<  "* is " << x1.nelem() << ".\n"
         << "The length of *"<< x2_name <<  "* is " << x2.nelem() << ".";
      throw runtime_error( os.str() );
    }
}



//! chk_if_increasing
/*! 
    Checks if a vector is strictly increasing.

    Duplicated values are not allowed.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_increasing( 
        const String&      x_name,
        ConstVectorView    x ) 
{
  if ( !is_increasing(x) )
    {
      ostringstream os;
      os << "The vector *" << x_name <<  "* must have strictly\n"
         << "increasing values, but this is not the case.\n";
      os << "x = " << x << "\n";
      throw runtime_error( os.str() );
    }
}



//! chk_if_decreasing
/*! 
    Checks if a vector is strictly decreasing.

    Duplicated values are not allowed.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_decreasing( 
        const String&      x_name,
        ConstVectorView    x ) 
{
  if ( !is_decreasing(x) )
    {
      ostringstream os;
      os << "The vector *" << x_name <<  "* must have strictly\ndecreasing "
         << "values, but this is not the case.\n";
      throw runtime_error( os.str() );
    }
}



/*===========================================================================
  === Functions for Matrix
  ===========================================================================*/

//! chk_matrix_ncols
/*! 
    Checks that a matrix has the specified number of columns.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A matrix.
    \param    l        The expected length of x.

    \author Patrick Eriksson 
    \date   2002-05-16
*/
void chk_matrix_ncols( 
        const String&      x_name,
        ConstMatrixView    x,
        const Index&       l ) 
{
  if ( x.ncols() != l )
    {
      ostringstream os;
      os << "The matrix *" << x_name <<  "* must have " << l << " columns,\n"
         << "but the number of columns is " << x.ncols() << ".";
      throw runtime_error( os.str() );
    }
}



//! chk_matrix_nrows
/*! 
    Checks that a matrix has the specified number of rows.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A matrix.
    \param    l        The expected length of x.

    \author Patrick Eriksson 
    \date   2002-05-16
*/
void chk_matrix_nrows( 
        const String&      x_name,
        ConstMatrixView    x,
        const Index&       l ) 
{
  if ( x.nrows() != l )
    {
      ostringstream os;
      os << "The matrix *" << x_name <<  "* must have " << l << " rows,\n"
         << "but the number of rows is " << x.nrows() << ".";
      throw runtime_error( os.str() );
    }
}



/*===========================================================================
  === Functions related to atmospheric grids, fields and surfaces.
  ===========================================================================*/

//! chk_atm_grids 
/*! 
    Checks if the atmospheric grids and the specified atmospheric 
    dimensionality match, and if the grids are ordered correctly.

    The function gives an error message if this is not the case.

    \param    dim          The atmospheric dimensionality.
    \param    p_grid       The pressure grid.
    \param    lat_grid     The latitude grid.
    \param    lon_grid     The longitude grid.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_atm_grids( 
        const Index&      dim,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid )
{
  // p_grid
  if( p_grid.nelem() < 2 )
    throw runtime_error( "The length of *p_grid* must be >= 2." );
  chk_if_decreasing( "p_grid", p_grid );

  // lat_grid
  if( dim == 1 )
    {
      if( lat_grid.nelem() != 0 )
        throw runtime_error(
                          "For dim=1, the length of *lat_grid* must be 0." );
    }
  else
    {
      if( lat_grid.nelem() < 2 )
        throw runtime_error(
                         "For dim>1, the length of *lat_grid* must be >= 2.");
      chk_if_increasing( "lat_grid", lat_grid );
    }

  // lon_grid
  if( dim < 3 )
    { 
      if( lon_grid.nelem() != 0 )
        throw runtime_error(
                           "For dim<3, the length of *lon_grid* must be 0." );
    }
  else
    {
      if( lon_grid.nelem() < 2 )
        throw runtime_error(
                          "For dim=3, the length of *lon_grid* must be >= 2.");
      chk_if_increasing( "lon_grid", lon_grid );
    }

  // Check that latitude and longitude grids are inside OK ranges for 3D
  if( dim == 3 )
    {
      if( lat_grid[0] < -90 )
        throw runtime_error( 
                  "The latitude grid cannot extend below -90 degrees for 3D" );
      if( lat_grid[lat_grid.nelem() - 1] > 90 )
        throw runtime_error( 
                  "The latitude grid cannot extend above 90 degrees for 3D" );
      if( lon_grid[0] < -180 )
        throw runtime_error( 
                "The longitude grid cannot extend below -180 degrees for 3D" );
      if( lon_grid[lon_grid.nelem() - 1] > 180 )
        throw runtime_error( 
                "The longitude grid cannot extend above 180 degrees for 3D" );
    }
}



//! chk_atm_field (simple fields)
/*! 
    Checks if an atmospheric field matches the dimensionality and the grids.

    The function gives an error message if this is not the case.

    \param    x_name       The name of the atmospheric field.
    \param    x            A variable holding an atmospheric field.
    \param    dim          The atmospheric dimensionality.
    \param    p_grid       The pressure grid.
    \param    lat_grid     The latitude grid.
    \param    lon_grid     The longitude grid.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_atm_field( 
        const String&     x_name,
        const Tensor3&    x, 
        const Index&      dim,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid )
{
  Index npages=p_grid.nelem(), nrows=1, ncols=1;
  if( dim > 1 )
    nrows = lat_grid.nelem();
  if( dim > 2 )
    ncols = lon_grid.nelem();
  if( x.ncols()!=ncols || x.nrows()!=nrows || x.npages()!=npages ) 
    {
      ostringstream os;
      os << "The atmospheric field *" << x_name <<  "* has wrong size.\n"
         << "Expected size is " << npages << " x " << nrows << " x " 
         << ncols << ", while actual size is " << x.npages() << " x " 
         << x.nrows() << " x " << x.ncols() << ".";
      throw runtime_error( os.str() );
    }
}

//! chk_atm_field (fields with one more dimension)
/*! 
    Checks if an atmospheric field matches the dimensionality and the
    grids. This is the version for fields like vmr_field, which are a
    Tensor4, not a Tensor3. (First dimension is the gas species.)

    The function gives an error message if this is not the case.

    \param    x_name       The name of the atmospheric field.
    \param    x            A variable holding an atmospheric field.
    \param    dim          The atmospheric dimensionality.
    \param    p_grid       The pressure grid.
    \param    lat_grid     The latitude grid.
    \param    lon_grid     The longitude grid.

    \author Stefan Buehler, cloned from Patrick Eriksson 
    \date   2002-12-20
*/
void chk_atm_field( 
        const String&   x_name,
        const Tensor4&  x, 
        const Index&    dim,
        const Index&    nspecies,
        ConstVectorView p_grid,
        ConstVectorView lat_grid,
        ConstVectorView lon_grid )
{
  Index npages=p_grid.nelem(), nrows=1, ncols=1;
  if( dim > 1 )
    nrows = lat_grid.nelem();
  if( dim > 2 )
    ncols = lon_grid.nelem();

  const Index nbooks=nspecies;

  if( x.ncols()!=ncols || x.nrows()!=nrows || x.npages()!=npages ||
      x.nbooks()!=nbooks ) 
    {
      ostringstream os;
      os << "The atmospheric field *" << x_name <<  "* has wrong size.\n"
         << "Expected size is "
         << nbooks << " x " << npages << " x "
         << nrows << " x " << ncols << ",\n"
         << "while actual size is "
         << x.nbooks() << " x " << x.npages() << " x "
         << x.nrows() << " x " << x.ncols() << ".";
      throw runtime_error( os.str() );
    }
}


//! chk_atm_surface
/*! 
    Checks if an atmospheric surface matches the dimensionality and the grids.

    An example on an atmospheric surface is *z_ground*.

    The function gives an error message if this is not the case.

    \param    x_name       The name of the atmospheric surface.
    \param    x            A variable holding an atmospheric surface.
    \param    dim          The atmospheric dimensionality.
    \param    lat_grid     The latitude grid.
    \param    lon_grid     The longitude grid.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_atm_surface( 
        const String&     x_name,
        const Matrix&     x, 
        const Index&      dim,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid )
{
  Index ncols=1, nrows=1;
  if( dim > 1 )
    nrows = lat_grid.nelem();
  if( dim > 2 )
    ncols = lon_grid.nelem();
  if( x.ncols()!=ncols || x.nrows()!=nrows ) 
    {
      ostringstream os;
      os << "The atmospheric surface *" << x_name <<  "* has wrong size.\n"
         << "Expected size is " << nrows << " x " << ncols << ","
         << " while actual size is " << x.nrows() << " x " << x.ncols() << ".";
      throw runtime_error( os.str() );
    }
}



/*===========================================================================
  === Function(s) releated with the cloud box.
  ===========================================================================*/

//! chk_cloudbox
/*! 
    Checks the consistency of the cloud box workspace variables. 

    The consistency is checked both internally and with respect to the grids.

    The function gives an error message if a consistency failure is found.

    \param    dim                The atmospheric dimensionality.
    \param    p_grid             The pressure grid.
    \param    lat_grid           The latitude grid.
    \param    lon_grid           The longitude grid.
    \param    cloudbox_on        Flag to activate the cloud box.
    \param    cloudbox_limits    Index limits of the cloud box.

    \author Patrick Eriksson 
    \date   2002-05-11
*/
void chk_cloudbox(
        const Index&          dim,
        ConstVectorView       p_grid,
        ConstVectorView       lat_grid,
        ConstVectorView       lon_grid,
        const Index&          cloudbox_on,    
        const ArrayOfIndex&   cloudbox_limits )
{
  chk_if_bool(  "cloudbox_on", cloudbox_on );

  if( cloudbox_on )
    {
      if( cloudbox_limits.nelem() != dim*2 )
        {
          ostringstream os;
          os << "The array *cloudbox_limits* has incorrect length.\n"
             << "For dim = " << dim << " the length shall be " << dim*2
             << " but it is " << cloudbox_limits.nelem() << ".";
          throw runtime_error( os.str() );
        }
       if( cloudbox_limits[1]<=cloudbox_limits[0] || cloudbox_limits[0]<0 ||
                                           cloudbox_limits[1]>=p_grid.nelem() )
        {
          ostringstream os;
          os << "Incorrect value(s) for cloud box pressure limit(s) found."
             << "\nValues are either out of range or upper limit is not "
             << "greater than lower limit.\nWith present length of "
             << "*p_grid*, OK values are 0 - " << p_grid.nelem()-1
             << ".\nThe pressure index limits are set to " 
             << cloudbox_limits[0] << " - " << cloudbox_limits[1] << ".";
          throw runtime_error( os.str() );
        }
      if( dim >= 2 )
        {
          if( cloudbox_limits[3]<=cloudbox_limits[2] || cloudbox_limits[2]<1 ||
                                cloudbox_limits[3]>=lat_grid.nelem()-1 )
            {
              ostringstream os;
              os << "Incorrect value(s) for cloud box latitude limit(s) found."
                 << "\nValues are either out of range or upper limit is not "
                 << "greater than lower limit.\nWith present length of "
                 << "*lat_grid*, OK values are 1 - " << lat_grid.nelem()-2
                 << ".\nThe latitude index limits are set to " 
                 << cloudbox_limits[2] << " - " << cloudbox_limits[3] << ".";
              throw runtime_error( os.str() );
            }
        }
      if( dim >= 3 )
        {
          if( cloudbox_limits[5]<=cloudbox_limits[4] || cloudbox_limits[4]<1 ||
                                cloudbox_limits[5]>=lon_grid.nelem()-1 )
            {
              ostringstream os;
              os << "Incorrect value(s) for cloud box longitude limit(s) found"
                 << ".\nValues are either out of range or upper limit is not "
                 << "greater than lower limit.\nWith present length of "
                 << "*lon_grid*, OK values are 1 - " << lon_grid.nelem()-2
                 << ".\nThe longitude limits are set to " 
                 << cloudbox_limits[4] << " - " << cloudbox_limits[5] << ".";
              throw runtime_error( os.str() );
            }
        }
    }
}



/*===========================================================================
  === Functions for Agendas
  ===========================================================================*/

//! chk_not_empty
/*! 
    Checks that an agenda is not empty.

    The function gives an error message if the agenda is empty.

    \param    x_name   The name of the agenda.
    \param    x        A variable of type Agenda.

    \author Patrick Eriksson 
    \date   2002-08-20
*/
void chk_not_empty( 
        const String&      x_name,
        const Agenda&      x ) 
{
  if( x.nelem() == 0 )
    {
      ostringstream os;
      os << "The agenda *" << x_name << "* is empty.\nIt is not allowed \n"
         << "that an agenda that is actually used to be empty.\n"
         << "Empty agendas are only created of methods setting dummy values \n"
         << "to variables.";
      throw runtime_error( os.str() );
    }
}

/*===========================================================================
  === Functions for Tensors
  ===========================================================================*/

//! Runtime check for size of Vector.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Agenda.
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size( const String&    x_name,
               ConstVectorView  x,
               const Index&     c ) 
{
  if ( !is_size(x,c) )
    {
      ostringstream os;
      os << "The object *" << x_name
         << "* does not have the right size.\n"
         << "Dimensions should be:"
         << " " << c 
         << ",\nbut they are:         "
         << " " << x.nelem()      
         << ".";
      throw runtime_error( os.str() );
    }
}

//! Runtime check for size of Matrix.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Agenda.
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size( const String&    x_name,
               ConstMatrixView  x,
               const Index&     r,
               const Index&     c ) 
{
  if ( !is_size(x,r,c) )
    {
      ostringstream os;
      os << "The object *" << x_name
         << "* does not have the right size.\n"
         << "Dimensions should be:"
         << " " << r 
         << " " << c 
         << ",\nbut they are:         "
         << " " << x.nrows()      
         << " " << x.ncols()      
         << ".";
      throw runtime_error( os.str() );
    }
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Agenda.
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size( const String&    x_name,
               ConstTensor3View x,
               const Index&     p,
               const Index&     r,
               const Index&     c ) 
{
  if ( !is_size(x,p,r,c) )
    {
      ostringstream os;
      os << "The object *" << x_name
         << "* does not have the right size.\n"
         << "Dimensions should be:"
         << " " << p 
         << " " << r 
         << " " << c 
         << ",\nbut they are:         "
         << " " << x.npages()     
         << " " << x.nrows()      
         << " " << x.ncols()      
         << ".";
      throw runtime_error( os.str() );
    }
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Agenda.
  \param    b        Required number of books
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size( const String&    x_name,
               ConstTensor4View x,
               const Index&     b,
               const Index&     p,
               const Index&     r,
               const Index&     c ) 
{
  if ( !is_size(x,b,p,r,c) )
    {
      ostringstream os;
      os << "The object *" << x_name
         << "* does not have the right size.\n"
         << "Dimensions should be:"
         << " " << b 
         << " " << p 
         << " " << r 
         << " " << c 
         << ",\nbut they are:         "
         << " " << x.nbooks()     
         << " " << x.npages()     
         << " " << x.nrows()      
         << " " << x.ncols()      
         << ".";
      throw runtime_error( os.str() );
    }
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Agenda.
  \param    s        Required number of shelves
  \param    b        Required number of books
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size( const String&    x_name,
               ConstTensor5View x,
               const Index&     s,
               const Index&     b,
               const Index&     p,
               const Index&     r,
               const Index&     c ) 
{
  if ( !is_size(x,s,b,p,r,c) )
    {
      ostringstream os;
      os << "The object *" << x_name
         << "* does not have the right size.\n"
         << "Dimensions should be:"
         << " " << s 
         << " " << b 
         << " " << p 
         << " " << r 
         << " " << c 
         << ",\nbut they are:         "
         << " " << x.nshelves()   
         << " " << x.nbooks()     
         << " " << x.npages()     
         << " " << x.nrows()      
         << " " << x.ncols()      
         << ".";
      throw runtime_error( os.str() );
    }
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Agenda.
  \param    v        Required number of vitrines
  \param    s        Required number of shelves
  \param    b        Required number of books
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size( const String&    x_name,
               ConstTensor6View x,
               const Index&     v,
               const Index&     s,
               const Index&     b,
               const Index&     p,
               const Index&     r,
               const Index&     c ) 
{
  if ( !is_size(x,v,s,b,p,r,c) )
    {
      ostringstream os;
      os << "The object *" << x_name
         << "* does not have the right size.\n"
         << "Dimensions should be:"
         << " " << v 
         << " " << s 
         << " " << b 
         << " " << p 
         << " " << r 
         << " " << c 
         << ",\nbut they are:         "
         << " " << x.nvitrines()  
         << " " << x.nshelves()   
         << " " << x.nbooks()     
         << " " << x.npages()     
         << " " << x.nrows()      
         << " " << x.ncols()      
         << ".";
      throw runtime_error( os.str() );
    }
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Agenda.
  \param    l        Required number of libraries
  \param    v        Required number of vitrines
  \param    s        Required number of shelves
  \param    b        Required number of books
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size( const String&    x_name,
               ConstTensor7View x,
               const Index&     l,
               const Index&     v,
               const Index&     s,
               const Index&     b,
               const Index&     p,
               const Index&     r,
               const Index&     c ) 
{
  if ( !is_size(x,l,v,s,b,p,r,c) )
    {
      ostringstream os;
      os << "The object *" << x_name
         << "* does not have the right size.\n"
         << "Dimensions should be:"
         << " " << l 
         << " " << v 
         << " " << s 
         << " " << b 
         << " " << p 
         << " " << r 
         << " " << c 
         << ",\nbut they are:         "
         << " " << x.nlibraries() 
         << " " << x.nvitrines()  
         << " " << x.nshelves()   
         << " " << x.nbooks()     
         << " " << x.npages()     
         << " " << x.nrows()      
         << " " << x.ncols()      
         << ".";
      throw runtime_error( os.str() );
    }
}

