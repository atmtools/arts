/* Copyright (C) 2008
   Oliver Lemke <olemke@core-dump.info>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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
  ===  File description
  ===========================================================================*/

/*!
  \file   m_create.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-07-18

  \brief  Workspace functions for creating variables

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"
#include "array.h"
#include "matpackI.h"
#include "matpackII.h"
#include "matpackIII.h"
#include "matpackIV.h"
#include "matpackV.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "mystring.h"
#include "absorption.h"
#include "gridded_fields.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfGField1Create(// WS Generic Output:
                          ArrayOfGField1& aogf1 )
{
  aogf1 = ArrayOfGField1();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfIndexCreate(// WS Generic Output:
                        ArrayOfIndex& aoi )
{
  aoi = ArrayOfIndex();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfLineRecordCreate(// WS Generic Output:
                             ArrayOfLineRecord& aolr )
{
  aolr = ArrayOfLineRecord();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfLineshapeSpecCreate(// WS Generic Output:
                                ArrayOfLineshapeSpec& aolsp )
{
  aolsp = ArrayOfLineshapeSpec();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfMatrixCreate(// WS Generic Output:
                         ArrayOfMatrix& aom )
{
  aom = ArrayOfMatrix();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfStringCreate(// WS Generic Output:
                         ArrayOfString& aos )
{
  aos = ArrayOfString();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfVectorCreate(// WS Generic Output:
                         ArrayOfVector& aov )
{
  aov = ArrayOfVector();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GField1Create(// WS Generic Output:
                 GField1& g )
{
  g = GField1();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GField2Create(// WS Generic Output:
                 GField2& g )
{
  g = GField2();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GField3Create(// WS Generic Output:
                 GField3& g )
{
  g = GField3();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GField4Create(// WS Generic Output:
                 GField4& g )
{
  g = GField4();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IndexCreate(// WS Generic Output:
                 Index& i )
{
  i = 0;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixCreate(// WS Generic Output:
                  Matrix& m )
{
  m = Matrix();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericCreate(// WS Generic Output:
                   Numeric& n )
{
  n = 0.;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void SparseCreate(// WS Generic Output:
                  Sparse& s )
{
  s = Sparse();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void StringCreate(// WS Generic Output:
                  String& s )
{
  s = "";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3Create(// WS Generic Output:
                   Tensor3& t )
{
  t = Tensor3();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor4Create(// WS Generic Output:
                   Tensor4& t )
{
  t = Tensor4();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor5Create(// WS Generic Output:
                   Tensor5& t _U_ )
{
  t = Tensor5();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor6Create(// WS Generic Output:
                   Tensor6& t _U_ )
{
  t = Tensor6();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor7Create(// WS Generic Output:
                   Tensor7& t )
{
  t = Tensor7();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorCreate(// WS Generic Output:
                  Vector& v )
{
  v = Vector();
}

