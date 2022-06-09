/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

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

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   bifstream.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-01-23

  \brief This file contains the class declaration of bifstream.

*/

#ifndef BIFSTREAM_H_INCLUDED
#define BIFSTREAM_H_INCLUDED

#include <fstream>

#include "binio.h"
#include "debug.h"

//! Binary output file stream class
/*!
  Handles writing to an output file stream in binary format. It makes it
  possible to use the operator<< for binary output.
*/
class bifstream final : public binistream, public ifstream {
 public:
  bifstream() : ifstream() {}

  explicit bifstream(const char* name,
                     ios::openmode mode = ios::in | ios::binary)
      : ifstream(name, mode) {
    // Open a second file descriptor for fast array reading
    if (!(this->mfilep = fopen(name, "rb"))) {
      ARTS_USER_ERROR("Failed to open ", name);
    }
  }

   ~bifstream() final {
    if (mfilep) {
      fclose(mfilep);
    }
  }

  void seek(long spos, Offset offs) final;
  streampos pos() final;

  bifstream::Byte getByte() final;
  void getRaw(char* c, streamsize n) final {
    if (n <= 8) {
      this->read(c, n);
    } else {
      fseek(mfilep, this->tellg(), SEEK_SET);
      size_t nread = fread(c, sizeof(char), n, mfilep);
      ARTS_USER_ERROR_IF((streamsize)nread != n,
                         "Unexpectedly reached end of binary input file.");
      seek(nread, Add);
    }
  }

private : FILE* mfilep{nullptr};
};

/* Overloaded input operators */
bifstream& operator>>(bifstream& bif, double& n);

bifstream& operator>>(bifstream& bif, float& n);

bifstream& operator>>(bifstream& bif, long& n);

bifstream& operator>>(bifstream& bif, int& n);

#endif
