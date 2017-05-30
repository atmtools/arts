/* Copyright (C) 2002-2012 Oliver Lemke <olemke@core-dump.info>

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
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-08-27
  
  \brief  Implementation of Select.
  
  This file contains the implementation of the supergeneric method Select.
*/

#ifndef m_select_h
#define m_select_h

#include "messages.h"
#include "mystring.h"
#include "workspace_ng.h"
#include "agenda_class.h"
#include "matpackII.h"


/* Workspace method: Doxygen documentation will be auto-generated */
template< class T >
void Select(// WS Generic Output:
            Array<T>& needles,
            // WS Generic Input:
            const Array<T>& haystack,
            const ArrayOfIndex& needleind,
            const Verbosity&)
{
  // We construct the output in this dummy variable, so that the
  // method also works properly if needles and haystack are the same
  // variable.
  Array<T> dummy( needleind.nelem() );

  // If needleind only contains -1 as the only element, copy the whole thing
  if (needleind.nelem() == 1 && needleind[0] == -1)
    {
      needles = haystack;
      return;
    }

  for( Index i = 0; i < needleind.nelem(); i++)
    {
      if (haystack.nelem() <= needleind[i])
        {
          ostringstream os;
          os << "The input vector only has " << haystack.nelem()
            << " elements. But one of the needle indexes is "
            << needleind[i] << "." << endl;
          os << "The indexes must be between 0 and " << haystack.nelem() - 1;
          throw runtime_error (os.str());
        }
      else if (needleind[i] < 0)
        {
          ostringstream os;
          os << "One of the needle indexes is " << needleind[i] << "." << endl;
          os << "The indexes must be between 0 and " << haystack.nelem() - 1;
          throw runtime_error (os.str());
        }
      else
        dummy[i] = haystack[needleind[i]];
    }

  needles = dummy;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Select(Workspace& /* ws */,
            // WS Generic Output:
            ArrayOfAgenda& needles,
            // WS Generic Input:
            const ArrayOfAgenda& haystack,
            const ArrayOfIndex& needleind,
            const Verbosity& verbosity)
{
  Select(needles, haystack, needleind, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Select(// WS Generic Output:
            Vector& needles,
            // WS Generic Input:
            const Vector& haystack,
            const ArrayOfIndex& needleind,
            const Verbosity&)
{
  // We construct the output in this dummy variable, so that the
  // method also works properly if needles and haystack are the same
  // variable.
  Vector dummy( needleind.nelem() );

  // If needleind only contains -1 as the only element, copy the whole thing
  if (needleind.nelem() == 1 && needleind[0] == -1)
    {
      needles = haystack;
      return;
    }

  for( Index i = 0; i < needleind.nelem(); i++)
    {
      if (haystack.nelem() <= needleind[i])
        {
          ostringstream os;
          os << "The input vector only has " << haystack.nelem()
            << " elements. But one of the needle indexes is "
            << needleind[i] << "." << endl;
          os << "The indexes must be between 0 and " << haystack.nelem() - 1;
          throw runtime_error (os.str());
        }
      else if (needleind[i] < 0)
        {
          ostringstream os;
          os << "One of the needle indexes is " << needleind[i] << "." << endl;
          os << "The indexes must be between 0 and " << haystack.nelem() - 1;
          throw runtime_error (os.str());
        }
      else
        dummy[i] = haystack[needleind[i]];
    }

  needles = dummy;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Select(// WS Generic Output:
            Matrix& needles,
            // WS Generic Input:
            const Matrix& haystack,
            const ArrayOfIndex& needleind,
            const Verbosity&)
{
  // We construct the output in this dummy variable, so that the
  // method also works properly if needles and haystack are the same
  // variable.
  Matrix dummy( needleind.nelem(), haystack.ncols() );

  // If needleind only contains -1 as the only element, copy the whole thing
  if (needleind.nelem() == 1 && needleind[0] == -1)
    {
      needles = haystack;
      return;
    }

  for( Index i = 0; i < needleind.nelem(); i++)
    {
      if (haystack.nrows() <= needleind[i])
        {
          ostringstream os;
          os << "The input matrix only has " << haystack.nrows()
            << " rows. But one of the needle indexes is "
            << needleind[i] << "." << endl;
          os << "The indexes must be between 0 and " << haystack.nrows() - 1;
          throw runtime_error (os.str());
        }
      else if (needleind[i] < 0)
        {
          ostringstream os;
          os << "One of the needle indexes is " << needleind[i] << "." << endl;
          os << "The indexes must be between 0 and " << haystack.nrows() - 1;
          throw runtime_error (os.str());
        }
      else
        dummy(i, joker) = haystack(needleind[i], joker);
    }

  needles = dummy;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Select(// WS Generic Output:
            Sparse& needles,
            // WS Generic Input:
            const Sparse& haystack,
            const ArrayOfIndex& needleind,
            const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // We construct the output in this dummy variable, so that the
  // method also works properly if needles and haystack are the same
  // variable.
  Sparse dummy( needleind.nelem(), haystack.ncols() );

  // If needleind only contains -1 as the only element, copy the whole thing
  if (needleind.nelem() == 1 && needleind[0] == -1)
    {
      needles = haystack;
      return;
    }

  for( Index i = 0; i < needleind.nelem(); i++)
    {
      if (haystack.nrows() <= needleind[i])
        {
          ostringstream os;
          os << "The input matrix only has " << haystack.nrows()
            << " rows. But one of the needle indexes is "
            << needleind[i] << "." << endl;
          os << "The indexes must be between 0 and " << haystack.nrows() - 1;
          throw runtime_error (os.str());
        }
      else if (needleind[i] < 0)
        {
          ostringstream os;
          os << "One of the needle indexes is " << needleind[i] << "." << endl;
          os << "The indexes must be between 0 and " << haystack.nrows() - 1;
          throw runtime_error (os.str());
        }
      else
        {
          // Copy this row of the sparse matrix.
          // This code is inefficient for Sparse, but I leave it like
          // this to be consistent with the other data types for which
          // Select is implemented.
          for ( Index j=0; j<haystack.ncols(); ++j)
            {
              Numeric value = haystack(needleind[i],j);
              if (0 != value)
                dummy.rw(i,j) = value;
            }
        }
    }

  if (dummy.nnz()==haystack.nnz())
    {
      // No data was actually removed.
      out3 << "  Number of nonzero elements has stayed the same.\n";
    }
  else
    {
      out3 << "  Number of nonzero elements reduced from "
           << haystack.nnz() << " to " << dummy.nnz() << ".\n";
    }

  needles = dummy;
}


#endif // m_select_h
