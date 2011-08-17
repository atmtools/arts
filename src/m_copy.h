/* Copyright (C) 2002-2008 Stefan Buehler <sbuehler@ltu.se>

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
  \file   m_copy.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 14 17:09:05 2002
  
  \brief  Implementation of Copy.
  
  This file contains the implementation of the supergeneric method
  Copy.
*/

#ifndef m_copy_h
#define m_copy_h

#include "messages.h"
#include "mystring.h"
#include "workspace_ng.h"
#include "agenda_class.h"


/* Workspace method: Doxygen documentation will be auto-generated */
template< class T >
void Copy(// WS Generic Output:
          T& out,
          // WS Generic Input:
          const T& in,
          const Verbosity&)
{
  // The use of CloneSize should not be necessary anymore, thanks to
  // the new copy semantics.  
  // CloneSize( out, outname, in, inname );
  out = in;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Copy(Workspace& ws _U_,
          // WS Generic Output:
          Agenda& out,
          // WS Generic Input:
          const Agenda& in,
          const Verbosity&)
{
  // The use of CloneSize should not be necessary anymore, thanks to
  // the new copy semantics.  
  // CloneSize( out, outname, in, inname );
  out = in;
}

#endif // m_copy_h
