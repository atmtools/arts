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
class bofstream : public binostream, public std::ofstream {
 public:
  bofstream() : std::ofstream() {}

  explicit bofstream(const char* name,
                     std::ios::openmode mode = std::ios::out | std::ios::trunc | std::ios::binary)
      : std::ofstream(name, mode) {}

  void seek(long spos, Offset offs) final;
  std::streampos pos() final;

  void putByte(bofstream::Byte b) final;
  void putRaw(const char* c, std::streamsize n) final { this->write(c, n); }
};

/* Overloaded output operators */
bofstream& operator<<(bofstream& bof, double n);

bofstream& operator<<(bofstream& bof, float n);

bofstream& operator<<(bofstream& bof, std::int64_t n);

bofstream& operator<<(bofstream& bof, int n);

#endif
