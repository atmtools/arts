/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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

/*!
  \file   wsv_aux.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Mar 14 10:28:23 2002
  
  \brief  Implementation of member functions for the workspace.
*/

#include "auto_wsv.h"

//! Default constructor.
WorkSpace::WorkSpace()
  : moccupied(N_WSV,0)
{
  // Nothing to do here.
}

//! Checks, whether a variable is occupied.
/*! 
  \param i Index of WSV to check.

  \return True if occupied.
*/
bool WorkSpace::is_occupied(Index i) const
{
  assert( 0 <= i );
  assert( i <  N_WSV );
  return 0 != moccupied[i];
}

//! Sets a variable as occupied.
/*! 
  \param i Index of WSV to set.
*/
void WorkSpace::set(Index i)
{
  assert( 0 <= i );
  assert( i <  N_WSV );
  moccupied[i] = 1;
}

//! Sets a variable as unoccupied.
/*! 
  \param i Index of WSV to set.
*/
void WorkSpace::free(Index i)
{
  assert( 0 <= i );
  assert( i <  N_WSV );
  moccupied[i] = 0;
}
