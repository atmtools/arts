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

bofstream& operator<<(bofstream& bof, std::int64_t n) {
  bof.writeInt(n, 4);
  return (bof);
}

bofstream& operator<<(bofstream& bof, int n) {
  bof.writeInt(n, 4);
  return (bof);
}
