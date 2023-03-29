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
  \file   bofstream.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-01-21

  \brief This file contains the class declaration of bofstream.

*/

#ifndef BOFSTREAM_H_INCLUDED
#define BOFSTREAM_H_INCLUDED

#include <fstream>

#include "binio.h"

//! Binary output file stream class
/*!
  Handles writing to an output file stream in binary format. It makes it
  possible to use the operator<< for binary output.
*/
class bofstream : public binostream, public ofstream {
 public:
  bofstream() : ofstream() {}

  explicit bofstream(const char* name,
                     ios::openmode mode = ios::out | ios::trunc | ios::binary)
      : ofstream(name, mode) {}

  void seek(long spos, Offset offs) override final;
  streampos pos() override final;

  void putByte(bofstream::Byte b) override final;
  void putRaw(const char* c, streamsize n) override final { this->write(c, n); }
};

/* Overloaded output operators */
bofstream& operator<<(bofstream& bof, double n);

bofstream& operator<<(bofstream& bof, float n);

bofstream& operator<<(bofstream& bof, std::int64_t n);

bofstream& operator<<(bofstream& bof, int n);

#endif
