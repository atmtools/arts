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
  \file   bofstream.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-01-21

  \brief This file contains the class implementation of bofstream.

*/

#include "bofstream.h"
#include <fstream>
#include <stdexcept>

void bofstream::seek(long spos, Offset offs) {
  if (!in) {
    err = NotOpen;
    return;
  }

  switch (offs) {
    case Set:
      this->seekp(spos, ios::beg);
      break;
    case Add:
      this->seekp(spos, ios::cur);
      break;
    case End:
      this->seekp(spos, ios::end);
      break;
  }
}

streampos bofstream::pos() {
  if (!in) {
    err = NotOpen;
    return 0;
  }
  return streamoff(this->tellp());
}

void bofstream::putByte(bofstream::Byte b) {
  if (!this->good()) {
    err |= NotOpen;
    ARTS_USER_ERROR ("Cannot open binary file for writing");
    return;
  }

  this->put(b);
  if (this->bad()) {
    err |= Fatal;
    ARTS_USER_ERROR ("Writing to binary file failed");
  }
}

/* Overloaded output operators */
bofstream& operator<<(bofstream& bof, double n) {
  bof.writeFloat(n, binio::Double);
  return (bof);
}

bofstream& operator<<(bofstream& bof, float n) {
  bof.writeFloat(n, binio::Double);
  return (bof);
}

bofstream& operator<<(bofstream& bof, long n) {
  bof.writeInt(n, 4);
  return (bof);
}

bofstream& operator<<(bofstream& bof, int n) {
  bof.writeInt(n, 4);
  return (bof);
}
