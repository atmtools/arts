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

/**
 * @file   messages.cc
 * @brief  Definitions having to do with the four output streams.
 *
 * See file messages.h for more explanations.
 *
 * @author Stefan Buehler
 * @date   2000-07-31
 */

#include "messages.h"
#include "array.h"
#include "arts.h"
#include "mystring.h"

// The global message verbosity settings:
Verbosity verbosity_at_launch;

/** The output path. For example for the report file. */
String out_path;

/** The basename for the report file and for all other output
    files. The files will have names like \<basename\>.rep (for the
    report file). */
String out_basename;

/** The report file. */
ofstream report_file;

ostream& operator<<(ostream& os, const Verbosity& v) {
  os << "Agenda Verbosity: " << v.get_agenda_verbosity() << "\n";
  os << "Screen Verbosity: " << v.get_screen_verbosity() << "\n";
  os << "File Verbosity  : " << v.get_file_verbosity() << "\n";

  return os;
}
