/* Copyright (C) 2002-2012 Oliver Lemke <olemke@core-dump.info>

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
  \file   xml_io.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/

#ifndef xml_io_h
#define xml_io_h

#include "absorption.h"
#include "mystring.h"

enum FileType {
  FILE_TYPE_ASCII = 0,
  FILE_TYPE_ZIPPED_ASCII = 1,
  FILE_TYPE_BINARY = 2
};

enum NumericType { NUMERIC_TYPE_FLOAT, NUMERIC_TYPE_DOUBLE };
enum EndianType { ENDIAN_TYPE_LITTLE, ENDIAN_TYPE_BIG };

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_xml(String& filename, const String& varname);

void filename_xml_with_index(String& filename,
                             const Index& file_index,
                             const String& varname,
                             const Index& digits = 0);

////////////////////////////////////////////////////////////////////////////
//   Generic IO routines for XML files
////////////////////////////////////////////////////////////////////////////

template <typename T>
void xml_read_from_file(const String& filename,
                        T& type,
                        const Verbosity& verbosity);

template <typename T>
void xml_write_to_file(const String& filename,
                       const T& type,
                       const FileType ftype,
                       const Index no_clobber,
                       const Verbosity& verbosity);

template <typename T>
void xml_read_from_file(const String&, T&, const Verbosity&);

template <typename T>
void xml_write_to_file(const String&, const T&, FileType, const Verbosity&);

#endif
