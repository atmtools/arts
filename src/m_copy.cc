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
  \file   m_copy.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri Jun 14 17:09:05 2002
  
  \brief  Implementation of Copy.
  
  This file contains the implementation of the supergeneric method
  Copy.
*/

#include "arts.h"
#include "mystring.h"
#include "matpackVII.h"
#include "array.h"
#include "absorption.h"
#include "ppath.h"
#include "agenda_class.h"
#include "auto_md.h"


void Copy_sg_Index(// WS Generic Output:
          Index& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Index& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Index( out, outname, in, inname );
  out = in;
}

void Copy_sg_Numeric(// WS Generic Output:
          Numeric& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Numeric& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Numeric( out, outname, in, inname );
  out = in;
}

void Copy_sg_String(// WS Generic Output:
          String& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const String& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_String( out, outname, in, inname );
  out = in;
}

void Copy_sg_Vector(// WS Generic Output:
          Vector& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Vector& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Vector( out, outname, in, inname );
  out = in;
}

void Copy_sg_Matrix(// WS Generic Output:
          Matrix& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Matrix& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Matrix( out, outname, in, inname );
  out = in;
}

void Copy_sg_Tensor3(// WS Generic Output:
          Tensor3& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Tensor3& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Tensor3( out, outname, in, inname );
  out = in;
}

void Copy_sg_Tensor4(// WS Generic Output:
          Tensor4& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Tensor4& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Tensor4( out, outname, in, inname );
  out = in;
}

void Copy_sg_Tensor5(// WS Generic Output:
          Tensor5& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Tensor5& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Tensor5( out, outname, in, inname );
  out = in;
}

void Copy_sg_Tensor6(// WS Generic Output:
          Tensor6& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Tensor6& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Tensor6( out, outname, in, inname );
  out = in;
}

void Copy_sg_Tensor7(// WS Generic Output:
          Tensor7& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Tensor7& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Tensor7( out, outname, in, inname );
  out = in;
}

void Copy_sg_ArrayOfIndex(// WS Generic Output:
          ArrayOfIndex& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const ArrayOfIndex& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_ArrayOfIndex( out, outname, in, inname );
  out = in;
}

void Copy_sg_ArrayOfString(// WS Generic Output:
          ArrayOfString& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const ArrayOfString& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_ArrayOfString( out, outname, in, inname );
  out = in;
}

void Copy_sg_ArrayOfVector(// WS Generic Output:
          ArrayOfVector& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const ArrayOfVector& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_ArrayOfVector( out, outname, in, inname );
  out = in;
}

void Copy_sg_ArrayOfMatrix(// WS Generic Output:
          ArrayOfMatrix& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const ArrayOfMatrix& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_ArrayOfMatrix( out, outname, in, inname );
  out = in;
}

void Copy_sg_ArrayOfLineRecord(// WS Generic Output:
          ArrayOfLineRecord& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const ArrayOfLineRecord& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_ArrayOfLineRecord( out, outname, in, inname );
  out = in;
}

void Copy_sg_ArrayOfArrayOfLineRecord(// WS Generic Output:
          ArrayOfArrayOfLineRecord& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const ArrayOfArrayOfLineRecord& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_ArrayOfArrayOfLineRecord( out, outname, in, inname );
  out = in;
}

void Copy_sg_TagGroups(// WS Generic Output:
          TagGroups& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const TagGroups& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_TagGroups( out, outname, in, inname );
  out = in;
}

void Copy_sg_Ppath(// WS Generic Output:
          Ppath& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Ppath& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Ppath( out, outname, in, inname );
  out = in;
}

void Copy_sg_Agenda(// WS Generic Output:
          Agenda& out,
          // WS Generic Output Names:
          const String& outname,
          // WS Generic Input:
          const Agenda& in,
          // WS Generic Input Names:
          const String& inname)
{
  CloneSize_sg_Agenda( out, outname, in, inname );
  out = in;
}

