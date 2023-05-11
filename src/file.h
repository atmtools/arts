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

#include "double_imanip.h"
#include "mystring.h"

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_ascii(String& filename, const String& varname);

void filename_bin(String& filename, const std::string_view varname);

////////////////////////////////////////////////////////////////////////////
//   Functions to open and read ASCII files
////////////////////////////////////////////////////////////////////////////

void open_output_file(std::ofstream& file, const std::string_view name);

void cleanup_output_file(std::ofstream& file, const std::string_view name);

void open_input_file(std::ifstream& file, const std::string_view name);

[[nodiscard]] ArrayOfString read_text_from_stream(std::istream& is);

[[nodiscard]] ArrayOfString read_text_from_file(const std::string_view name);

void replace_all(String& s, const std::string_view what, const std::string_view with);

[[nodiscard]] int check_newline(const std::string_view s);

[[nodiscard]] bool file_exists(const std::string_view filename);

bool find_file(ArrayOfString& matches,
               const std::string_view filename,
               const ArrayOfString& paths,
               const ArrayOfString& extensions = {""});

void find_xml_file(String& filename);

[[nodiscard]] bool find_xml_file_existence(String& filename);

[[nodiscard]] String expand_path(String path);

[[nodiscard]] String add_basedir(const std::string_view path);

[[nodiscard]] String get_dirname(const std::string_view path);

[[nodiscard]] ArrayOfString list_directory(const std::string_view dirname);

[[nodiscard]] String make_filename_unique(const std::string_view filename, const String& extension = "");

#endif