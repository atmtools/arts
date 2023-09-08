/* -*-C++-*-
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * binio.h - Binary stream I/O classes
 * Copyright (C) 2002, 2003 Simon Peter <dn.tlp@gmx.net>
 */

#ifndef H_BINIO_BINIO
#define H_BINIO_BINIO

#include <cstdint>
#include <ios>

/***** Configuration *****/

// BINIO_ENABLE_STRING - Build std::string supporting methods
//
// Set to 1 to build std::string supporting methods. You need the STL to
// do this.
#define BINIO_ENABLE_STRING 1

// BINIO_ISO_STDLIB - Build with ISO C++ standard library compliance
//
// Set to 1 to build for the ISO standard C++ library (i.e. namespaces, STL and
// templatized iostream). Set to 0 to build for the traditional C++ library.
#define BINIO_ISO_STDLIB 1

/***** Implementation *****/

// Disable annoying multiple inheritance compiler warning on MSVC6
#ifdef _MSC_VER
#pragma warning(disable : 4250)
#endif

#if BINIO_ENABLE_STRING
#include <string>
#endif

class binio {
 public:
  using Flag = enum { BigEndian = 1 << 0, FloatIEEE = 1 << 1 };

  using ErrorCode = enum {
    NoError = 0,
    Fatal = 1 << 0,
    Unsupported = 1 << 1,
    NotOpen = 1 << 2,
    Denied = 1 << 3,
    NotFound = 1 << 4,
    Eof = 1 << 5
  };

  using Offset = enum { Set, Add, End };
  using FType = enum { Single, Double };
  using Error = int;

  binio();

  void setFlag(Flag f, bool set = true);
  bool getFlag(Flag f);

  Error error();
  bool eof();

  virtual void seek(long, Offset = Set) = 0;
  virtual std::streampos pos() = 0;

 protected:
  using Int = std::int64_t;
  using Float = double;
  using Byte = unsigned char;  // has to be unsigned!

  using Flags = int;

  Flags my_flags{system_flags};
  static const Flags system_flags;
  Error err{NoError};

 private:
  static Flags detect_system_flags();
};

class binistream : virtual public binio {
 public:
  Int readInt(unsigned int size);
  Float readFloat(FType ft);

  void readDoubleArray(double *d, unsigned long size);

  unsigned long readString(char *str, unsigned long amount);
  unsigned long readString(char *str, unsigned long maxlen, const char delim);


#if BINIO_ENABLE_STRING
  std::string readString(const char delim = '\0');
#endif

  Int peekInt(unsigned int size);
  Float peekFloat(FType ft);

  bool ateof();
  void ignore(unsigned long amount = 1);

 protected:
  virtual Byte getByte() = 0;
  virtual void getRaw(char *c, std::streamsize n) = 0;
};

class binostream : virtual public binio {
 public:
  void writeInt(Int val, unsigned int size);
  void writeFloat(Float f, FType ft);
  unsigned long writeString(const char *str, unsigned long amount = 0);
#if BINIO_ENABLE_STRING
  unsigned long writeString(const std::string &str);
#endif

 protected:
  virtual void putByte(Byte) = 0;
  virtual void putRaw(const char *c, std::streamsize n) = 0;
};

#endif
