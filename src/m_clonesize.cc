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
  \file   m_clonesize.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri Jun 14 17:09:05 2002
  
  \brief  Implementation of CloneSize.
  
  This file contains the implementation of the supergeneric method
  CloneSize. In this case we do not use a C++ template function,
  because the different groups behave too differently. We have to make
  explicit handlers for all groups.

  If you are uncertain what CloneSize should do for a particular
  group, remember this: Copy calls CloneSize, then the assignment
  operator "=". So, CloneSize should prepare out in such a way, that
  the subsequent assignment will do the right thing.
*/

#include "arts.h"
#include "mystring.h"
#include "messages.h"
#include "matpackVII.h"
#include "array.h"
#include "absorption.h"
#include "ppath.h"
#include "optproperties.h"
#include "agenda_class.h"
#include "auto_md.h"

void CloneSize(// WS Generic Output:
               Index& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Index& in,
               // WS Generic Input Names:
               const String& inname)
{
  out3 << "  You asked me to make " << outname
       << " the same size as " << inname
       << " (" << in
       << ").\n"
       << "  These are scalars, so there is nothing to do.\n";
  out = 0;
}

void CloneSize(// WS Generic Output:
               Numeric& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Numeric& in,
               // WS Generic Input Names:
               const String& inname)
{
  out3 << "  You asked me to make " << outname
       << " the same size as " << inname
       << " (" << in
       << ").\n"
       << "  These are scalars, so there is nothing to do.\n";
  out = 0.0;
}

void CloneSize(// WS Generic Output:
               String& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const String& in,
               // WS Generic Input Names:
               const String& inname)
{
  out3 << "  You asked me to make " << outname
       << " the same size as " << inname
       << " (" << in
       << ").\n"
       << "  These are scalars, so there is nothing to do.\n";
  out = "";
}

void CloneSize(// WS Generic Output:
               Vector& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Vector& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out = 0;
}

void CloneSize(// WS Generic Output:
               Matrix& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Matrix& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nrows(),
              in.ncols() );
  out = 0;
}

void CloneSize(// WS Generic Output:
               Tensor3& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Tensor3& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.npages(),
              in.nrows(),
              in.ncols() );
  out = 0;
}

void CloneSize(// WS Generic Output:
               Tensor4& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Tensor4& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nbooks(),
              in.npages(),
              in.nrows(),
              in.ncols() );
  out = 0;
}

void CloneSize(// WS Generic Output:
               Tensor5& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Tensor5& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nshelves(),
              in.nbooks(),
              in.npages(),
              in.nrows(),
              in.ncols() );
  out = 0;
}

void CloneSize(// WS Generic Output:
               Tensor6& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Tensor6& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nvitrines(),
              in.nshelves(),
              in.nbooks(),
              in.npages(),
              in.nrows(),
              in.ncols() );
  out = 0;
}

void CloneSize(// WS Generic Output:
               Tensor7& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Tensor7& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nlibraries(),
              in.nvitrines(),
              in.nshelves(),
              in.nbooks(),
              in.npages(),
              in.nrows(),
              in.ncols() );
  out = 0;
}

void CloneSize(// WS Generic Output:
               ArrayOfIndex& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfIndex& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out = 0;
}

void CloneSize(// WS Generic Output:
               ArrayOfString& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfString& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out = "";
}

void CloneSize(// WS Generic Output:
               ArrayOfVector& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfVector& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out2 << "  Created an Array of empty Vectors.\n";
}

void CloneSize(// WS Generic Output:
               ArrayOfMatrix& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfMatrix& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out2 << "  Created an Array of empty Matrices.\n";
}


void CloneSize(// WS Generic Output:
               ArrayOfTensor3& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfTensor3& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out2 << "  Created an Array of empty Tensor3.\n";
}

void CloneSize(// WS Generic Output:
               ArrayOfArrayOfTensor3& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfArrayOfTensor3& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out2 << "  Created an Array of Arrays of empty Tensor3.\n";
}

void CloneSize(// WS Generic Output:
               ArrayOfTensor6& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfTensor6& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out2 << "  Created an Array of empty Tensor6.\n";
}

void CloneSize(// WS Generic Output:
               ArrayOfArrayOfTensor6& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfArrayOfTensor6& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out2 << "  Created an Array of Arrays of empty Tensor6.\n";
}


void CloneSize(// WS Generic Output:
               ArrayOfLineRecord& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfLineRecord& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
}

void CloneSize(// WS Generic Output:
               ArrayOfArrayOfLineRecord& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfArrayOfLineRecord& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out2 << "  Created an Array of empty ArrayOfLineRecords.\n";
}

void CloneSize(// WS Generic Output:
               ArrayOfArrayOfSpeciesTag& out,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfArrayOfSpeciesTag& in,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out.resize( in.nelem() );
  out2 << "  Created an Array of empty ArrayOfArrayOfSpeciesTag.\n";
}

void CloneSize(// WS Generic Output:
               Ppath&,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Ppath&,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out2 << "  Nothing to do for this group.\n";
}

void CloneSize(// WS Generic Output:
               Agenda&,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const Agenda&,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out2 << "  Nothing to do for this group.\n";
}

void CloneSize(// WS Generic Output:
               GridPos&,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const GridPos&,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out2 << "  Nothing to do for this group.\n";
}

void CloneSize(// WS Generic Output:
               GasAbsLookup&,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const GasAbsLookup&,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out2 << "  Nothing to do for this group.\n";
}

void CloneSize(// WS Generic Output:
               SingleScatteringData&,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const SingleScatteringData&,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out2 << "  Nothing to do for this group.\n";
}

void CloneSize(// WS Generic Output:
               ArrayOfSingleScatteringData&,
               // WS Generic Output Names:
               const String& outname,
               // WS Generic Input:
               const ArrayOfSingleScatteringData&,
               // WS Generic Input Names:
               const String& inname)
{
  out2 << "  Making " << outname << " the same size as " << inname << ".\n";
  out2 << "  Nothing to do for this group.\n";
}
