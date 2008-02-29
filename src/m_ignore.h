/* Copyright (C) 2002-2007 Stefan Buehler <sbuehler@ltu.se>

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
  \file   m_ignore.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 14 17:09:05 2002
  
  \brief  Implementation of Ignore.
  
  This file contains the implementation of the supergeneric method
  Ignore.
*/

#ifndef m_ignore_h
#define m_ignore_h

#include "exceptions.h"
#include "messages.h"
#include "mystring.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template< class T >
void Ignore(// WS Generic Input:
            const T& in _U_)
{
}


/* Workspace method: Doxygen documentation will be auto-generated */
template< class T >
void DoNothing(
            // WS Generic Output:
                  T&      out,
            // WS Generic Input:
                  T&      in)
{
  if( &out != &in )
    throw runtime_error(
                   "*DoNothing* requires that input and output is same WSV." );
}

#endif // m_ignore_h
