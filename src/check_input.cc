/* Copyright (C) 2000, 2001 Stefan Buehler  <sbuehler@uni-bremen.de>
                            Axel von Engeln <engeln@uni-bremen.de>

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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////

/**
  \file   check_input.cc

  General functions to check the size and logic of input to functions.

  \author Patrick Eriksson
  \date 2002-04-15 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "check_input.h"
#include "logic.h"



////////////////////////////////////////////////////////////////////////////
//   Functions to assert the size of a vector, matrix or tensor.
////////////////////////////////////////////////////////////////////////////

//// assert_size (Vector) /////////////////////////////////////////////////////
/** 
    Asserts that a vector has the expected length.

    \param    x   A variable of type Vector.
    \param    l   The expected length of x.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void assert_size( 
	ConstVectorView   x,
	const Index&      l ) 
{
  assert( x.nelem() == l );
}



//// assert_size (Matrix) /////////////////////////////////////////////////////
/** 
    Asserts that a matrix has the expected size.

    \param    x       A variable of type Matrix.
    \param    nrows   The expected number of rows of x.
    \param    ncols   The expected number of columns of x.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void assert_size( 
	ConstMatrixView   x,
	const Index&      nrows,
	const Index&      ncols ) 
{
  assert( x.ncols() == ncols );
  assert( x.nrows() == nrows );
}



//// assert_size (Tensor3) ////////////////////////////////////////////////////
/** 
    Asserts that a tensor of order 3 has the expected size.

    \param    x       A variable of type Tensor3.
    \param    npages  The expected number of pages of x.
    \param    nrows   The expected number of rows of x.
    \param    ncols   The expected number of columns of x.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void assert_size( 
	const Tensor3&    x,
	const Index&      npages,
	const Index&      nrows,
	const Index&      ncols ) 
{
  assert( x.ncols() == ncols );
  assert( x.nrows() == nrows );
  assert( x.npages() == npages );
}



//// assert_maxdim_of_sensor //////////////////////////////////////////////////
/** 
    Asserts that the effective dimension of a tensor of order 3 does not
    exceed the expected value.

    For dimensions not used, the size shall be 1.

    \param    x     A variable of type Tensor3.
    \param    dim   The expected effective dimension of x.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void assert_maxdim_of_tensor(
	const Tensor3&   x,
	const Index&     dim )
{
  assert( dim >= 1 );
  assert( dim <= 3 );

  if( dim == 1 )
    {
      assert( x.nrows() == 1 );
      assert( x.ncols() == 1 );
    }
  else if( dim == 2 )
    assert( x.ncols() == 1 );
}



////////////////////////////////////////////////////////////////////////////
//   Functions for Index
////////////////////////////////////////////////////////////////////////////

//// chk_if_bool //////////////////////////////////////////////////////////////
/** 
    Checks that a variable of type Index has the value 0 or 1.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Index.

    \author Patrick Eriksson 
    \date   2004-04-15
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



//// chk_if_in_range //////////////////////////////////////////////////////////
/** 
    Checks that a variable of type Index has a value inside the specified
    range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Index.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Patrick Eriksson 
    \date   2004-04-15
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



////////////////////////////////////////////////////////////////////////////
//   Functions for Numeric
////////////////////////////////////////////////////////////////////////////

//// chk_if_over_0 ////////////////////////////////////////////////////////////
/** 
    Checks that a variable of type Numeric is 0 or is positive.
    range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.

    \author Patrick Eriksson 
    \date   2004-04-15
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



//// chk_if_in_range //////////////////////////////////////////////////////////
/** 
    Checks that a variable of type Numeric has a value inside the specified
    range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Patrick Eriksson 
    \date   2004-04-15
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



////////////////////////////////////////////////////////////////////////////
//   Functions for Vector
////////////////////////////////////////////////////////////////////////////

//// chk_vector_length ////////////////////////////////////////////////////////
/** 
    Checks that a vector has the specified length.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.
    \param    l        The expected length of x.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void chk_vector_length( 
	const String&      x_name,
        ConstVectorView&   x,
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



//// chk_vector_length ////////////////////////////////////////////////////////
/** 
    Checks if two vectors have the same length.

    The function gives an error message if this is not the case.

    \param    x1_name   The name of the first vector
    \param    x2_name   The name of the second vector
    \param    x1        The first vector.
    \param    x2        The second vector.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void chk_vector_length( 
	const String&      x1_name,
	const String&      x2_name,
        ConstVectorView&   x1, 
        ConstVectorView&   x2 ) 
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



//// chk_if_increasing ////////////////////////////////////////////////////////
/** 
    Checks if a vector is strictly increasing.

    Duplicated values are not allowed.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void chk_if_increasing( 
	const String&      x_name,
        ConstVectorView&   x ) 
{
  if ( !is_increasing(x) )
    {
      ostringstream os;
      os << "The vector *" << x_name <<  "* must have strictly\nincreasing "
         << "values, but this is not the case.";
      throw runtime_error( os.str() );
    }
}



//// chk_if_decreasing ////////////////////////////////////////////////////////
/** 
    Checks if a vector is strictly decreasing.

    Duplicated values are not allowed.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void chk_if_decreasing( 
	const String&      x_name,
        ConstVectorView&   x ) 
{
  if ( !is_decreasing(x) )
    {
      ostringstream os;
      os << "The vector *" << x_name <<  "* must have strictly\ndecreasing "
         << "values, but this is not the case.";
      throw runtime_error( os.str() );
    }
}



////////////////////////////////////////////////////////////////////////////
//   Functions related to atmospheric fields, grids etc.
////////////////////////////////////////////////////////////////////////////

//// chk_atm_grids ////////////////////////////////////////////////////////////
/** 
    Checks if the atmospheric grids and the specified atmospheric 
    dimensionality match, and if the grids are ordered correctly.

    The function gives an error message if this is not the case.

    \param    dim          The atmospheric dimensionality.
    \param    p_grid       The pressure grid.
    \param    alpha_grid   The latitude grid.
    \param    beta_grid    The longitude grid.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void chk_atm_grids( 
	const Index&      dim,
	ConstVectorView   p_grid,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid )
{
  if( p_grid.nelem() < 2 )
    throw runtime_error( "The length of *p_grid* must be >= 2." );
  chk_if_decreasing( "p_grid", p_grid );

  if( dim == 1 )
    {
      if( alpha_grid.nelem() != 0 )
	throw runtime_error(
                          "For dim=1, the length of *alpha_grid* must be 0." );
    }
  else
    {
      if( alpha_grid.nelem() < 2 )
	throw runtime_error(
                         "For dim>1, the length of *alpha_grid* must be >=2.");
      chk_if_increasing( "alpha_grid", alpha_grid );
    }

  if( dim < 3 )
    { 
      if( beta_grid.nelem() != 0 )
	throw runtime_error(
                           "For dim<3, the length of *beta_grid* must be 0." );
    }
  else
    {
      if( beta_grid.nelem() < 2 )
	throw runtime_error(
                          "For dim=3, the length of *beta_grid* must be >=2.");
      chk_if_increasing( "beta_grid", beta_grid );
    }
}



//// chk_atm_field ////////////////////////////////////////////////////////////
/** 
    Checks if an atmospheric field matches the dimensionality and the grids.

    The function gives an error message if this is not the case.

    \param    x_name       The name of the atmospheric field.
    \param    x            A variable holding an atmospheric field.
    \param    dim          The atmospheric dimensionality.
    \param    p_grid       The pressure grid.
    \param    alpha_grid   The latitude grid.
    \param    beta_grid    The longitude grid.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void chk_atm_field( 
	const String&     x_name,
        const Tensor3&    x, 
	const Index&      dim,
	ConstVectorView   p_grid,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid )
{
  Index npages=p_grid.nelem(), nrows=1, ncols=1;
  if( dim > 1 )
    nrows = alpha_grid.nelem();
  if( dim > 2 )
    ncols = beta_grid.nelem();
  if( x.ncols()!=ncols || x.nrows()!=nrows || x.npages()!=npages ) 
    {
      ostringstream os;
      os << "The atmospheric field *" << x_name <<  "* has wrong size.\n"
         << "Expected size is " << npages << " x " << nrows << " x " 
         << ncols << ",while actual size is " << x.npages() << " x " 
         << x.nrows() << " x " << x.ncols() << ".";
      throw runtime_error( os.str() );
    }
}



//// chk_atm_surface ////////////////////////////////////////////////////////
/** 
    Checks if an atmospheric surface matches the dimensionality and the grids.

    An example on an atmospheric surface is *z_ground*.

    The function gives an error message if this is not the case.

    \param    x_name       The name of the atmospheric surface.
    \param    x            A variable holding an atmospheric surface.
    \param    dim          The atmospheric dimensionality.
    \param    alpha_grid   The latitude grid.
    \param    beta_grid    The longitude grid.

    \author Patrick Eriksson 
    \date   2004-04-15
*/
void chk_atm_surface( 
	const String&     x_name,
        const Matrix&     x, 
	const Index&      dim,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid )
{
  Index ncols=1, nrows=1;
  if( dim > 1 )
    nrows = alpha_grid.nelem();
  if( dim > 2 )
    ncols = beta_grid.nelem();
  if( x.ncols()!=ncols || x.nrows()!=nrows ) 
    {
      ostringstream os;
      os << "The atmospheric surface *" << x_name <<  "* has wrong size.\n"
         << "Expected size is " << nrows << " x " << ncols << ","
         << "while actual size is " << x.nrows() << " x " << x.ncols() << ".";
      throw runtime_error( os.str() );
    }
}
