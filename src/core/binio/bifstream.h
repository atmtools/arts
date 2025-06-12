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

#include <cstdint>
#include <fstream>

#include "binio.h"

//! Binary output file stream class
/*!
  Handles writing to an output file stream in binary format. It makes it
  possible to use the operator<< for binary output.
*/
class bifstream final : public binistream, public std::ifstream {
 public:
  bifstream() : std::ifstream() {}

  explicit bifstream(const char* name,
                     std::ios::openmode mode = std::ios::in | std::ios::binary);

  ~bifstream() final {
    if (mfilep) {
      fclose(mfilep);
    }
  }

  void seek(long spos, Offset offs) final;
  std::streampos pos() final;

  bifstream::Byte getByte() final;
  void getRaw(char* c, std::streamsize n) final;

 private:
  FILE* mfilep{nullptr};
};

/* Overloaded input operators */
bifstream& operator>>(bifstream& bif, double& n);

bifstream& operator>>(bifstream& bif, float& n);

bifstream& operator>>(bifstream& bif, std::int64_t& n);

bifstream& operator>>(bifstream& bif, int& n);

#endif
