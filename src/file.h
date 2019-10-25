/* Copyright (C) 2000-2012 Stefan Buehler <sbuehler@ltu.se>

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
#include "messages.h"
#include "mystring.h"

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_ascii(String& filename, const String& varname);

void filename_bin(String& filename, const String& varname);

////////////////////////////////////////////////////////////////////////////
//   Functions to open and read ASCII files
////////////////////////////////////////////////////////////////////////////

void open_output_file(std::ofstream& file, const String& name);

void cleanup_output_file(std::ofstream& file, const String& name);

void open_input_file(std::ifstream& file, const String& name);

void read_text_from_stream(ArrayOfString& text, std::istream& is);

void read_text_from_file(ArrayOfString& text, const String& name);

void replace_all(String& s, const String& what, const String& with);

int check_newline(const String& s);

bool file_exists(const String& filename);

bool find_file(ArrayOfString& matches,
               const String& filename,
               const ArrayOfString& paths,
               const ArrayOfString& extensions = {""});

void find_xml_file(String& filename, const Verbosity& verbosity);

bool find_xml_file_existence(String& filename);

String expand_path(const String& path);

String add_basedir(const String& path);

void get_dirname(String& dirname, const String& path);

void list_directory(ArrayOfString& files, String dirname);

void make_filename_unique(String& filename, const String& extension = "");

////////////////////////////////////////////////////////////////////////////
//   IO manipulation classes for parsing nan and inf
////////////////////////////////////////////////////////////////////////////

// Not all compilers do support parsing of nan and inf.
// The code below is taken (and slightly modified) from the discussion at:
// http://stackoverflow.com/questions/11420263/is-it-possible-to-read-infinity-or-nan-values-using-input-streams

/** Input stream class for doubles that correctly handles nan and inf. */
class double_istream {
 public:
  double_istream(std::istream& i) : in(i) {}

  double_istream& parse_on_fail(double& x, bool neg);

  double_istream& operator>>(double& x) {
    bool neg = false;
    char c;
    if (!in.good()) return *this;
    while (isspace(c = (char)in.peek())) in.get();
    if (c == '-') {
      neg = true;
    }
    in >> x;
    if (!in.fail()) return *this;
    return parse_on_fail(x, neg);
  }

 private:
  std::istream& in;
};

/** Input manipulator class for doubles to enable nan and inf parsing. */
class double_imanip {
 public:
  const double_imanip& operator>>(double& x) const {
    double_istream(*in) >> x;
    return *this;
  }
  std::istream& operator>>(const double_imanip&) const { return *in; }

  friend const double_imanip& operator>>(std::istream& in,
                                         const double_imanip& dm);

 private:
  mutable std::istream* in;
};

const double_imanip& operator>>(std::istream& in, const double_imanip& dm);

#endif
