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

#include <debug.h>

#include <stdexcept>

bifstream::bifstream(const char* name, std::ios::openmode mode)
    : std::ifstream(name, mode) {
  // Open a second file descriptor for fast array reading
  if (!(this->mfilep = fopen(name, "rb"))) {
    ARTS_USER_ERROR("Failed to open {}", name);
  }
}

void bifstream::getRaw(char* c, std::streamsize n) {
  if (n <= 8) {
    this->read(c, n);
  } else {
    fseek(mfilep, this->tellg(), SEEK_SET);
    size_t nread = fread(c, sizeof(char), n, mfilep);
    ARTS_USER_ERROR_IF((std::streamsize)nread != n,
                       "Unexpectedly reached end of binary input file.");
    seek(nread, Add);
  }
}

void bifstream::seek(long spos, Offset offs) {
  if (!in) {
    err = NotOpen;
    return;
  }

  switch (offs) {
    case Set: this->seekg(spos, std::ios::beg); break;
    case Add: this->seekg(spos, std::ios::cur); break;
    case End: this->seekg(spos, std::ios::end); break;
  }
}

std::streampos bifstream::pos() {
  if (!in) {
    err = NotOpen;
    return 0;
  }
  return std::streampos(this->tellg());
}

bifstream::Byte bifstream::getByte() {
  if (this->good()) {
    int iread;
    iread = this->get();
    if (iread == EOF) err |= Eof;
    return (Byte)iread;
  }

  err |= NotOpen;
  throw std::runtime_error("Reading from binary file failed");
  return 0;
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
