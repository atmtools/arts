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
  \file   bifstream.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-01-23

  \brief This file contains the class implementation of bifstream.

*/

#include "bifstream.h"
#include <fstream>
#include <stdexcept>

void bifstream::seek(long spos, Offset offs) {
  if (!in) {
    err = NotOpen;
    return;
  }

  switch (offs) {
    case Set:
      this->seekg(spos, ios::beg);
      break;
    case Add:
      this->seekg(spos, ios::cur);
      break;
    case End:
      this->seekg(spos, ios::end);
      break;
  }
}

streampos bifstream::pos() {
  if (!in) {
    err = NotOpen;
    return 0;
  }
  return streampos(this->tellg());
}

bifstream::Byte bifstream::getByte() {
  if (this->good()) {
    int iread;
    iread = this->get();
    if (iread == EOF) err |= Eof;
    return (Byte)iread;
  } else {
    err |= NotOpen;
    throw runtime_error("Reading from binary file failed");
    return 0;
  }
}

/* Overloaded input operators */
bifstream& operator>>(bifstream& bif, double& n) {
  n = (double)bif.readFloat(binio::Double);
  return (bif);
}

bifstream& operator>>(bifstream& bif, float& n) {
  n = (float)bif.readFloat(binio::Double);
  return (bif);
}

bifstream& operator>>(bifstream& bif, std::int64_t& n) {
  n = static_cast<std::int64_t>(bif.readInt(4));
  return (bif);
}

bifstream& operator>>(bifstream& bif, int& n) {
  n = (int)bif.readInt(4);
  return (bif);
}
