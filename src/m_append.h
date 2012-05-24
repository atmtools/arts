/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

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
  \file   m_append.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 14 17:09:05 2002
  
  \brief  Implementation of Append.
  
  This file contains the implementation of the supergeneric method
  Append.
*/

#ifndef m_append_h
#define m_append_h

#include "array.h"
#include "exceptions.h"
#include "matpackI.h"


/* Implementations for supported types follow. */

/* Implementation for array types */
template< class T >
void Append(// WS Generic Output:
            Array<T>& out,
            // WS Generic Input:
            const Array<T>& in,
            const Verbosity&)
{
  // Reserve memory in advance to avoid reallocations:
  out.reserve(out.nelem()+in.nelem());
  // Append in to end of out:
  for (Index i=0; i<in.nelem(); ++i)
    out.push_back(in[i]);
}

/* Implementation for array types to append single element */
template< class T >
void Append(// WS Generic Output:
            Array<T>& out,
            // WS Generic Input:
            const T& in,
            const Verbosity&)
{
  // Append in to end of out:
  out.push_back(in);
}

/* Implementation for Vector */
void Append(// WS Generic Output:
            Vector& out,
            // WS Generic Input:
            const Vector& in,
            const Verbosity&)
{
  // Get backup of out:
  Vector dummy = out;

  // Make out the right size:
  out.resize(dummy.nelem()+in.nelem());

  // Copy dummy to first part of out:
  out[Range(0,dummy.nelem())] = dummy;
  
  // Copy in to last part of out:
  out[Range(dummy.nelem(),in.nelem())] = in;
}

/* Implementation for String */
void Append(// WS Generic Output:
            String& out,
            // WS Generic Input:
            const String& in,
            const Verbosity&)
{
  // String stream for easy string operations:
  ostringstream os;
   
  os << out << in;

  out = os.str();
}

#endif // m_append_h
