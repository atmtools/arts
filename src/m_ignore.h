/* Copyright (C) 2002
   Stefan Buehler <sbuehler@uni-bremen.de>

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
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri Jun 14 17:09:05 2002
  
  \brief  Implementation of Copy.
  
  This file contains the implementation of the supergeneric method
  Ignore.
*/

#ifndef m_ignore_h
#define m_ignore_h

//! Supergeneric Copy.
/*! 
  This is the implementation of the supergeneric Ignore method. See
  arts -d Ignore for a description what the method does.

  \param in Source WSV.
  \param inname Name of source WSV.
*/
template< class T >
void Ignore(// WS Generic Input:
            const T& in,
            // WS Generic Input Names:
            const String& inname)
{
  // Nothing to do here.
  out2 << "  Ignoring " << inname << ".\n";
}

#endif // m_ignore_h
