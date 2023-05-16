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
