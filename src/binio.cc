/*
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
 * binio.cpp - Binary stream I/O classes
 * Copyright (C) 2002, 2003 Simon Peter <dn.tlp@gmx.net>
 */

#include <cstring>
#include <stdexcept>

#include "binio.h"
#include "debug.h"

/***** Defines *****/

#if BINIO_ENABLE_STRING
// String buffer size for std::string readString() method
#define STRINGBUFSIZE 256
#endif

/***** binio *****/

const binio::Flags binio::system_flags = binio::detect_system_flags();

binio::Flags binio::detect_system_flags() {
  Flags f = 0;

  // Endian test
  union {
    int word;
    Byte byte;
  } endian_test;

  endian_test.word = 1;
  if (endian_test.byte != 1) f |= BigEndian;

  // IEEE-754 floating-point test
  float fl = 6.5f;
  Byte *dat = (Byte *)&fl;

  if (sizeof(float) == 4 && sizeof(double) == 8) {
    if (f & BigEndian) {
      if (dat[0] == 0x40 && dat[1] == 0xD0 && !dat[2] && !dat[3])
        f |= FloatIEEE;
    } else {
      if (dat[3] == 0x40 && dat[2] == 0xD0 && !dat[1] && !dat[0])
        f |= FloatIEEE;
    }
  }

  return f;
}

binio::binio() {
  if (!(system_flags & FloatIEEE)) {
    ARTS_USER_ERROR(
        "IEEE floating-point numbers are not supported on "
        "this system.");
  }

  // Set Little Endian mode, with IEEE-754 floats.
  this->setFlag(binio::BigEndian, false);  // remove flag
  this->setFlag(binio::FloatIEEE);         // set flag
}

void binio::setFlag(Flag f, bool set) {
  if (set)
    my_flags |= f;
  else
    my_flags &= !f;
}

bool binio::getFlag(Flag f) { return (my_flags & f ? true : false); }

binio::Error binio::error() {
  Error e = err;

  err = NoError;
  return e;
}

bool binio::eof() { return (err & Eof ? true : false); }

/***** binistream *****/
binistream::Int binistream::readInt(unsigned int size) {
  unsigned int i;
  Int val = 0, in;

  // Check if 'size' doesn't exceed our system's biggest type.
  if (size > sizeof(Int)) {
    err |= Unsupported;
    throw runtime_error(
        "The size of the integer to be read exceeds our system's biggest type");
    return 0;
  }

  for (i = 0; i < size; i++) {
    in = getByte();
    if (getFlag(BigEndian))
      val <<= 8;
    else
      in <<= i * 8;
    val |= in;
  }

  return val;
}

binistream::Float binistream::readFloat(FType ft) {
  unsigned int i = 0;
  unsigned int size = 0;
  Byte in[8];
  bool swap;

  // Determine appropriate size for given type.
  switch (ft) {
    case Single:
      size = 4;
      break;  // 32 bits
    case Double:
      size = 8;
      break;  // 64 bits
  }

  // Determine byte ordering, depending on what we do next
  swap = getFlag(BigEndian) ^ (system_flags & BigEndian);

  if (!swap && ((size == sizeof(float)) || (size == sizeof(double)))) {
    if (size == 4) {
      float f;
      getRaw((char *)&f, size);
      return (Float)f;
    }

    double d;
    getRaw((char *)&d, size);
    return (Float)d;
  }

  // Read the float byte by byte, converting endianess
  for (i = 0; i < size; i++)
    if (swap)
      in[size - i - 1] = getByte();
    else
      in[i] = getByte();

  // Compatible system, let the hardware do the conversion
  switch (ft) {
    case Single:
      return *(float *)in;
    case Double:
      return *(double *)in;
  }

  // User tried to read a (yet) unsupported floating-point type. Bail out.
  ARTS_USER_ERROR("Unsupported floating-point type");
}

void swap_endian(double *d) {
  auto *val = reinterpret_cast<uint64_t *>(d);
  // Using this is much faster than swapping bytes manually
  *val = ((((*val) >> 56) & 0x00000000000000FF) |
          (((*val) >> 40) & 0x000000000000FF00) |
          (((*val) >> 24) & 0x0000000000FF0000) |
          (((*val) >> 8) & 0x00000000FF000000) |
          (((*val) << 8) & 0x000000FF00000000) |
          (((*val) << 24) & 0x0000FF0000000000) |
          (((*val) << 40) & 0x00FF000000000000) |
          (((*val) << 56) & 0xFF00000000000000));
}

void binistream::readDoubleArray(double *d, unsigned long size) {
  getRaw((char *)d, sizeof(double) * size);

  const bool swap = getFlag(BigEndian) ^ (system_flags & BigEndian);
  if (swap)
    for (double *dptr = d; dptr < d + size; dptr++) swap_endian(dptr);
}

unsigned long binistream::readString(char *str, unsigned long maxlen) {
  unsigned long i;

  for (i = 0; i < maxlen; i++) {
    str[i] = (char)getByte();
    if (err) {
      str[i] = '\0';
      return i;
    }
  }

  return maxlen;
}

unsigned long binistream::readString(char *str,
                                     unsigned long maxlen,
                                     const char delim) {
  unsigned long i;

  for (i = 0; i < maxlen; i++) {
    str[i] = (char)getByte();
    if (str[i] == delim || err) {
      str[i] = '\0';
      return i;
    }
  }

  str[maxlen] = '\0';
  return maxlen;
}

#if BINIO_ENABLE_STRING
std::string binistream::readString(const char delim) {
  char buf[STRINGBUFSIZE + 1];
  std::string tempstr;
  unsigned long read;

  do {
    read = readString(buf, STRINGBUFSIZE, delim);
    tempstr.append(buf, read);
  } while (read == STRINGBUFSIZE);

  return tempstr;
}
#endif

binistream::Int binistream::peekInt(unsigned int size) {
  Int val = readInt(size);
  if (!err) seek(-(long)size, Add);
  return val;
}

binistream::Float binistream::peekFloat(FType ft) {
  Float val = readFloat(ft);

  if (!err) switch (ft) {
      case Single:
        seek(-4, Add);
        break;
      case Double:
        seek(-8, Add);
        break;
    }

  return val;
}

bool binistream::ateof() {
  Error olderr = err;  // Save current error state
  bool eof_then;

  peekInt(1);
  eof_then = eof();  // Get error state of next byte
  err = olderr;      // Restore original error state
  return eof_then;
}

void binistream::ignore(unsigned long amount) {
  unsigned long i;

  for (i = 0; i < amount; i++) getByte();
}

/***** binostream *****/

void binostream::writeInt(Int val, unsigned int size) {
  unsigned int i;

  // Check if 'size' doesn't exceed our system's biggest type.
  if (size > sizeof(Int)) {
    err |= Unsupported;
    throw runtime_error(
        "The size of the integer to be stored exceeds our system's biggest type");
    return;
  }

  for (i = 0; i < size; i++) {
    if (getFlag(BigEndian))
      putByte((unsigned char)(val >> ((size - i - 1) * 8)) & 0xff);
    else {
      putByte((unsigned char)val & 0xff);
      val >>= 8;
    }
  }
}

void binostream::writeFloat(Float f, FType ft) {
  unsigned int i = 0;
  unsigned int size = 0;
  Byte *out = nullptr;
  bool swap;

  auto outf = (float)f;
  auto outd = (double)f;

  // Hardware could be big or little endian, convert appropriately
  swap = getFlag(BigEndian) ^ (system_flags & BigEndian);

  switch (ft) {
    case Single:
      size = 4;
      break;  // 32 bits
    case Double:
      size = 8;
      break;  // 64 bits
  }

  if (!swap && ((size == sizeof(float)) || (size == sizeof(double)))) {
    if (size == 4) {
      putRaw((char *)&outf, size);
      return;
    }

    putRaw((char *)&outd, size);
    return;
  }
  // Determine appropriate size for given type and convert by hardware
  switch (ft) {
    case Single:
      out = (Byte *)&outf;
      break;  // 32 bits
    case Double:
      out = (Byte *)&outd;
      break;  // 64 bits
  }

  // Write the float byte by byte, converting endianess
  if (swap) out += size - 1;
  for (i = 0; i < size; i++) {
    putByte(*out);
    if (swap)
      out--;
    else
      out++;
  }

  return;  // We're done.
}

unsigned long binostream::writeString(const char *str, unsigned long amount) {
  unsigned int i;

  if (!amount) amount = strlen(str);

  for (i = 0; i < amount; i++) {
    putByte(str[i]);
    if (err) return i;
  }

  return amount;
}

#if BINIO_ENABLE_STRING
unsigned long binostream::writeString(const std::string &str) {
  return writeString(str.c_str());
}
#endif
