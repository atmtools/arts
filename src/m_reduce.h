/* Copyright (C) 2014 Richard Larsson <ric.larsson@gmail.com>

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
  
  \file   m_reduce.h
  \author Richard Larsson <ric.larsson@gmail.com>
  \date   Wed Feb 25 2014
  
  \brief  Implementation of Reduce.
  
  This file contains the implementation of the supergeneric method
  Reduce.
*/

#ifndef m_reduce_h
#define m_reduce_h

#include <cstring>
#include "array.h"
#include "exceptions.h"
#include "matpackVII.h"

//////// Helper function /////////

Index num_elem_from_dim_sizes(const ArrayOfIndex& dim_sizes) {
  Index m = 1;
  for (ArrayOfIndex::const_iterator it = dim_sizes.begin();
       it != dim_sizes.end();
       it++)
    m *= *it;

  return m;
}

// Need this for each matpack type: Vector
void select_dims_by_size(ArrayOfIndex& dim_sizes,
                         const Index min_num_elem,
                         const Vector& type) {
  dim_sizes.resize(0);
  if (type.nelem() > min_num_elem) dim_sizes.push_back(type.nelem());
}

// Need this for each matpack type: Matrix
void select_dims_by_size(ArrayOfIndex& dim_sizes,
                         const Index min_num_elem,
                         const Matrix& type) {
  dim_sizes.resize(0);
  if (type.nrows() > min_num_elem) dim_sizes.push_back(type.nrows());
  if (type.ncols() > min_num_elem) dim_sizes.push_back(type.ncols());
}

// Need this for each matpack type: Tensor3
void select_dims_by_size(ArrayOfIndex& dim_sizes,
                         const Index min_num_elem,
                         const Tensor3& type) {
  dim_sizes.resize(0);
  if (type.npages() > min_num_elem) dim_sizes.push_back(type.npages());
  if (type.nrows() > min_num_elem) dim_sizes.push_back(type.nrows());
  if (type.ncols() > min_num_elem) dim_sizes.push_back(type.ncols());
}

// Need this for each matpack type: Tensor4
void select_dims_by_size(ArrayOfIndex& dim_sizes,
                         const Index min_num_elem,
                         const Tensor4& type) {
  dim_sizes.resize(0);
  if (type.nbooks() > min_num_elem) dim_sizes.push_back(type.nbooks());
  if (type.npages() > min_num_elem) dim_sizes.push_back(type.npages());
  if (type.nrows() > min_num_elem) dim_sizes.push_back(type.nrows());
  if (type.ncols() > min_num_elem) dim_sizes.push_back(type.ncols());
}

// Need this for each matpack type: Tensor5
void select_dims_by_size(ArrayOfIndex& dim_sizes,
                         const Index min_num_elem,
                         const Tensor5& type) {
  dim_sizes.resize(0);
  if (type.nshelves() > min_num_elem) dim_sizes.push_back(type.nshelves());
  if (type.nbooks() > min_num_elem) dim_sizes.push_back(type.nbooks());
  if (type.npages() > min_num_elem) dim_sizes.push_back(type.npages());
  if (type.nrows() > min_num_elem) dim_sizes.push_back(type.nrows());
  if (type.ncols() > min_num_elem) dim_sizes.push_back(type.ncols());
}

// Need this for each matpack type: Tensor6
void select_dims_by_size(ArrayOfIndex& dim_sizes,
                         const Index min_num_elem,
                         const Tensor6& type) {
  dim_sizes.resize(0);
  if (type.nvitrines() > min_num_elem) dim_sizes.push_back(type.nvitrines());
  if (type.nshelves() > min_num_elem) dim_sizes.push_back(type.nshelves());
  if (type.nbooks() > min_num_elem) dim_sizes.push_back(type.nbooks());
  if (type.npages() > min_num_elem) dim_sizes.push_back(type.npages());
  if (type.nrows() > min_num_elem) dim_sizes.push_back(type.nrows());
  if (type.ncols() > min_num_elem) dim_sizes.push_back(type.ncols());
}

// Need this for each matpack type: Tensor7
void select_dims_by_size(ArrayOfIndex& dim_sizes,
                         const Index min_num_elem,
                         const Tensor7& type) {
  dim_sizes.resize(0);
  if (type.nlibraries() > min_num_elem) dim_sizes.push_back(type.nlibraries());
  if (type.nvitrines() > min_num_elem) dim_sizes.push_back(type.nvitrines());
  if (type.nshelves() > min_num_elem) dim_sizes.push_back(type.nshelves());
  if (type.nbooks() > min_num_elem) dim_sizes.push_back(type.nbooks());
  if (type.npages() > min_num_elem) dim_sizes.push_back(type.npages());
  if (type.nrows() > min_num_elem) dim_sizes.push_back(type.nrows());
  if (type.ncols() > min_num_elem) dim_sizes.push_back(type.ncols());
}

//////// End helper functions /////////

// To Numeric

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Numeric& o,
    // WS Input:
    // WS Generic Input:
    const Vector& i,
    const Verbosity&) {
  if (i.nelem() == 1)
    o = i[0];
  else {
    ostringstream os;
    os << "The Vector is not also a Numeric";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Numeric& o,
    // WS Input:
    // WS Generic Input:
    const Matrix& i,
    const Verbosity&) {
  if (i.ncols() == 1 && i.nrows() == 1)
    o = i(0, 0);
  else {
    ostringstream os;
    os << "The Matrix is not also a Numeric";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Numeric& o,
    // WS Input:
    // WS Generic Input:
    const Tensor3& i,
    const Verbosity&) {
  if (i.ncols() == 1 && i.nrows() == 1 && i.npages() == 1)
    o = i(0, 0, 0);
  else {
    ostringstream os;
    os << "The Tensor3 is not also a Numeric";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Numeric& o,
    // WS Input:
    // WS Generic Input:
    const Tensor4& i,
    const Verbosity&) {
  if (i.ncols() == 1 && i.nrows() == 1 && i.npages() == 1 && i.nbooks() == 1)
    o = i(0, 0, 0, 0);
  else {
    ostringstream os;
    os << "The Tensor4 is not also a Numeric";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Numeric& o,
    // WS Input:
    // WS Generic Input:
    const Tensor5& i,
    const Verbosity&) {
  if (i.ncols() == 1 && i.nrows() == 1 && i.npages() == 1 && i.nbooks() == 1 &&
      i.nshelves() == 1)
    o = i(0, 0, 0, 0, 0);
  else {
    ostringstream os;
    os << "The Tensor5 is not also a Numeric";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Numeric& o,
    // WS Input:
    // WS Generic Input:
    const Tensor6& i,
    const Verbosity&) {
  if (i.ncols() == 1 && i.nrows() == 1 && i.npages() == 1 && i.nbooks() == 1 &&
      i.nshelves() == 1 && i.nvitrines() == 1)
    o = i(0, 0, 0, 0, 0, 0);
  else {
    ostringstream os;
    os << "The Tensor6 is not also a Numeric";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Numeric& o,
    // WS Input:
    // WS Generic Input:
    const Tensor7& i,
    const Verbosity&) {
  if (i.ncols() == 1 && i.nrows() == 1 && i.npages() == 1 && i.nbooks() == 1 &&
      i.nshelves() == 1 && i.nvitrines() == 1 && i.nlibraries() == 1)
    o = i(0, 0, 0, 0, 0, 0, 0);
  else {
    ostringstream os;
    os << "The Tensor7 is not also a Numeric";
    throw std::runtime_error(os.str());
  }
}

//To Vector

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Vector& o,
    // WS Input:
    // WS Generic Input:
    const Matrix& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 1;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Matrix of size (" << dim_sizes << ") \n"
       << "does not fit a Vector";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Vector& o,
    // WS Input:
    // WS Generic Input:
    const Tensor3& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 1;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor3 of size (" << dim_sizes << ") \n"
       << "does not fit a Vector";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Vector& o,
    // WS Input:
    // WS Generic Input:
    const Tensor4& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 1;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor4 of size (" << dim_sizes << ") \n"
       << "does not fit a Vector";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Vector& o,
    // WS Input:
    // WS Generic Input:
    const Tensor5& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 1;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor5 of size (" << dim_sizes << ") \n"
       << "does not fit a Vector";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Vector& o,
    // WS Input:
    // WS Generic Input:
    const Tensor6& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 1;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor6 of size (" << dim_sizes << ") \n"
       << "does not fit a Vector";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Vector& o,
    // WS Input:
    // WS Generic Input:
    const Tensor7& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 1;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor7 of size (" << dim_sizes << ") \n"
       << "does not fit a Vector";
    throw std::runtime_error(os.str());
  }
}

// To Matrix

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Matrix& o,
    // WS Input:
    // WS Generic Input:
    const Tensor3& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 2;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor3 of size (" << dim_sizes << ") \n"
       << "does not fit a Matrix";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Matrix& o,
    // WS Input:
    // WS Generic Input:
    const Tensor4& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 2;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor4 of size (" << dim_sizes << ") \n"
       << "does not fit a Matrix";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Matrix& o,
    // WS Input:
    // WS Generic Input:
    const Tensor5& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 2;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor5 of size (" << dim_sizes << ") \n"
       << "does not fit a Matrix";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Matrix& o,
    // WS Input:
    // WS Generic Input:
    const Tensor6& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 2;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor6 of size (" << dim_sizes << ") \n"
       << "does not fit a Matrix";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Matrix& o,
    // WS Input:
    // WS Generic Input:
    const Tensor7& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 2;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor7 of size (" << dim_sizes << ") \n"
       << "does not fit a Matrix";
    throw std::runtime_error(os.str());
  }
}

// To Tensor 3

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor3& o,
    // WS Input:
    // WS Generic Input:
    const Tensor4& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 3;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1], dim_sizes[2]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor4 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor3";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor3& o,
    // WS Input:
    // WS Generic Input:
    const Tensor5& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 3;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1], dim_sizes[2]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor5 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor3";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor3& o,
    // WS Input:
    // WS Generic Input:
    const Tensor6& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 3;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1], dim_sizes[2]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor6 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor3";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor3& o,
    // WS Input:
    // WS Generic Input:
    const Tensor7& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 3;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1], dim_sizes[2]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor7 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor3";
    throw std::runtime_error(os.str());
  }
}

// To Tensor4

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor4& o,
    // WS Input:
    // WS Generic Input:
    const Tensor5& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 4;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1], dim_sizes[2], dim_sizes[3]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor5 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor4";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor4& o,
    // WS Input:
    // WS Generic Input:
    const Tensor6& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 4;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1], dim_sizes[2], dim_sizes[3]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor6 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor4";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor4& o,
    // WS Input:
    // WS Generic Input:
    const Tensor7& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 4;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0], dim_sizes[1], dim_sizes[2], dim_sizes[3]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor7 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor4";
    throw std::runtime_error(os.str());
  }
}

// To Tensor5

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor5& o,
    // WS Input:
    // WS Generic Input:
    const Tensor6& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 5;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(
        dim_sizes[0], dim_sizes[1], dim_sizes[2], dim_sizes[3], dim_sizes[4]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor6 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor5";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor5& o,
    // WS Input:
    // WS Generic Input:
    const Tensor7& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 5;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(
        dim_sizes[0], dim_sizes[1], dim_sizes[2], dim_sizes[3], dim_sizes[4]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor7 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor5";
    throw std::runtime_error(os.str());
  }
}

// To Tensor6

/* Workspace method: Doxygen documentation will be auto-generated */
void Reduce(
    // WS Generic Output:
    Tensor6& o,
    // WS Input:
    // WS Generic Input:
    const Tensor7& i,
    const Verbosity&) {
  ArrayOfIndex dim_sizes;
  Index test = 6;

  select_dims_by_size(dim_sizes, 1, i);
  Index num = dim_sizes.nelem();
  if (num == test) {
    o.resize(dim_sizes[0],
             dim_sizes[1],
             dim_sizes[2],
             dim_sizes[3],
             dim_sizes[4],
             dim_sizes[5]);
    memcpy(o.get_c_array(),
           i.get_c_array(),
           sizeof(Numeric) * num_elem_from_dim_sizes(dim_sizes));
  } else {
    ostringstream os;
    os << "The Tensor7 of size (" << dim_sizes << ") \n"
       << "does not fit a Tensor6";
    throw std::runtime_error(os.str());
  }
}

#endif /* m_reduce_h */
