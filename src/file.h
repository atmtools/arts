/* Copyright (C) 2000-2008 Stefan Buehler <sbuehler@ltu.se>

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
/**
   \file  file.h

   This file contains basic functions to handle ASCII files.

   \author Patrick Eriksson
   \date 2000-10-28 
*/


#ifndef file_h
#define file_h

#include <fstream>
#include "matpackI.h"
#include "mystring.h"

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_ascii(
              String&  filename,
        const String&  varname );

void filename_bin(
              String&  filename,
        const String&  varname );



////////////////////////////////////////////////////////////////////////////
//   Functions to open and read ASCII files
////////////////////////////////////////////////////////////////////////////

void open_output_file(ofstream& file, const String& name);

void cleanup_output_file(ofstream& file, const String& name);

void open_input_file(ifstream& file, const String& name);

void read_text_from_stream(ArrayOfString& text, istream& is);

void read_text_from_file(ArrayOfString& text, const String& name);

String expand_path(const String& path);

void replace_all(String& s, const String& what, const String& with);

int check_newline(const String& s);

bool file_exists(const String& filename);

bool find_file(String& filename, const char* extension, const ArrayOfString& paths);

void get_dirname(String& dirname, const String& path);

#endif
