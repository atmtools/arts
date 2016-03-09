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
            const String& direction _U_,
            const Verbosity&)
{
    const Array<T>* in_pnt;
    Array<T> in_copy;

    if (&in == &out)
    {
        in_copy = in;
        in_pnt = &in_copy;
    }
    else
        in_pnt = &in;

    const Array<T>& in_ref = *in_pnt;

    // Reserve memory in advance to avoid reallocations:
    out.reserve(out.nelem() + in_ref.nelem());
    // Append in to end of out:
    for (Index i = 0; i < in_ref.nelem(); ++i)
        out.push_back(in_ref[i]);
}


/* Implementation for array types to append single element */
template< class T >
void Append(// WS Generic Output:
            Array<T>& out,
            // WS Generic Input:
            const T& in,
            const String& direction _U_,
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
            const String& direction _U_,
            const Verbosity&)
{
    const Vector* in_pnt;
    Vector in_copy;

    if (&in == &out)
    {
        in_copy = in;
        in_pnt = &in_copy;
    }
    else
        in_pnt = &in;

    const Vector &in_ref = *in_pnt;

    // Get backup of out:
    Vector dummy = out;

    // Make out the right size:
    out.resize(dummy.nelem() + in_ref.nelem());

    // Copy dummy to first part of out:
    if (dummy.nelem())
        out[Range(0, dummy.nelem())] = dummy;

    // Copy in to last part of out:
    if (in_ref.nelem())
        out[Range(dummy.nelem(), in_ref.nelem())] = in_ref;
}


/* Implementation for Matrix */
void Append(// WS Generic Output:
            Matrix& out,
            // WS Generic Input:
            const Matrix& in,
            const String& direction,
            const Verbosity&)
{
    const Matrix* in_pnt;
    Matrix in_copy;

    if (&in == &out)
    {
        in_copy = in;
        in_pnt = &in_copy;
    }
    else
        in_pnt = &in;

    const Matrix &in_ref = *in_pnt;

    // Get backup of out:
    Matrix dummy = out;

    if (!out.nrows() || !out.ncols())
    {
        out = in_ref;
    }
    else if (direction == "leading")
    {
        if (out.ncols() != in_ref.ncols())
            throw runtime_error("Input and output matrix must have the same number of columns.");

        out.resize(dummy.nrows() + in_ref.nrows(), dummy.ncols());

        if (dummy.nrows() && dummy.ncols())
            out(Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
        if (dummy.nrows() && in_ref.nrows() && in_ref.ncols())
            out(Range(dummy.nrows(), in_ref.nrows()), Range(0, in_ref.ncols())) = in_ref;
    }
    else if (direction == "trailing")
    {
        if (out.nrows() != in_ref.nrows())
            throw runtime_error("Input and output matrix must have the same number of rows.");

        out.resize(dummy.nrows(), dummy.ncols() + in_ref.ncols());

        if (dummy.nrows() && dummy.ncols())
            out(Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
        if (dummy.ncols() && in_ref.nrows() && in_ref.ncols())
            out(Range(0, in_ref.nrows()), Range(dummy.ncols(), in_ref.ncols())) = in_ref;
    }
    else throw runtime_error("Dimension must be either \"leading\" or \"trailing\".");
}


/* Implementation for Matrix/Vector */
void Append(// WS Generic Output:
            Matrix& out,
            // WS Generic Input:
            const Vector& in,
            const String& direction,
            const Verbosity&)
{
    // Get backup of out:
    Matrix dummy = out;

    if (direction == "leading")
    {
        if (!out.nrows() || !out.ncols())
        {
            out = in;
        }
        else
        {
            if (out.ncols() != in.nelem())
                throw runtime_error("Number of elements in the input Vector has to match "
                                    "the number of columns in the output Matrix.");

            out.resize(dummy.nrows() + 1, dummy.ncols());
            out(Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
            out(Range(dummy.nrows(), 1), Range(0, in.nelem())) = transpose(in);
        }
    }
    else if (direction == "trailing")
    {
        if (!out.nrows() || !out.ncols())
        {
            out = transpose(in);
        }
        else if (in.nelem())
        {
            if (out.nrows() != in.nelem() && out.nrows() && out.ncols())
                throw runtime_error("Number of elements in the input Vector has to match "
                                    "the number of rows in the output Matrix.");

            out.resize(dummy.nrows(), dummy.ncols() + 1);
            out(Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
            out(Range(0, in.nelem()), Range(dummy.ncols(), 1)) = in;
        }
    }
    else throw runtime_error("Dimension must be either \"leading\" or \"trailing\".");
}


/* Implementation for Vector/Numeric */
void Append(// WS Generic Output:
            Vector& out,
            // WS Generic Input:
            const Numeric& in,
            const String& direction _U_,
            const Verbosity&)
{
    // Get backup of out:
    Vector dummy = out;

    // Make out the right size:
    out.resize(dummy.nelem() + 1);

    // Copy dummy to first part of out:
    if (dummy.nelem())
        out[Range(0, dummy.nelem())] = dummy;

    // Copy in to last part of out:
    out[Range(dummy.nelem(), 1)] = in;
}


/* Implementation for Tensor4 */
void Append(// WS Generic Output:
            Tensor4& out,
            // WS Generic Input:
            const Tensor4& in,
//            const String& direction,
            const String& direction _U_,
            const Verbosity&)
{
    const Tensor4* in_pnt;
    Tensor4 in_copy;

    if (&in == &out)
    {
        in_copy = in;
        in_pnt = &in_copy;
    }
    else
        in_pnt = &in;

    const Tensor4 &in_ref = *in_pnt;

    // Get backup of out:
    Tensor4 dummy = out;

    if (out.npages() != in_ref.npages() || out.nrows() != in_ref.nrows() ||
        out.ncols() != in_ref.ncols())
        throw runtime_error("Input and output Tensor4 must have the same number"
                            "of books.");

    out.resize(dummy.nbooks() + in_ref.nbooks(), dummy.npages(),
               dummy.nrows(), dummy.ncols());

    if (dummy.nbooks() && dummy.npages() && dummy.nrows() && dummy.ncols())
        out(Range(0, dummy.nbooks()), Range(0, dummy.npages()),
            Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
    if (dummy.nbooks() && in_ref.nbooks() && in_ref.npages() &&
        in_ref.nrows() && in_ref.ncols())
        out(Range(dummy.nbooks(), in_ref.nbooks()), Range(0, in_ref.npages()),
            Range(0, in_ref.nrows()), Range(0, in_ref.ncols())) = in_ref;
}


/* Implementation for String */
void Append(// WS Generic Output:
            String& out,
            // WS Generic Input:
            const String& in,
            const String& direction _U_,
            const Verbosity&)
{
  // String stream for easy string operations:
  ostringstream os;
   
  os << out << in;

  out = os.str();
}

#endif // m_append_h
