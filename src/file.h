
/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

   This file contains basic functions to handle ASCII and binary (HDF)
   data files.

   \author Patrick Eriksson
   \date 2000-10-28 
*/


#ifndef file_h
#define file_h


////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "vecmat.h"

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_ascii(
              string&  filename,
        const string&  varname );

void filename_bin(
              string&  filename,
        const string&  varname );



////////////////////////////////////////////////////////////////////////////
//   Functions to open and read ASCII files
////////////////////////////////////////////////////////////////////////////

void open_output_file(ofstream& file, const string& name);

void open_input_file(ifstream& file, const string& name);

void read_text_from_stream(Array<string>& text, istream& is);

void read_text_from_file(Array<string>& text, const string& name);

void replace_all(string& s, const string& what, const string& with);



////////////////////////////////////////////////////////////////////////////
//   Matrix/Vector IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

void write_array_of_matrix_to_stream(ostream& os,
                                     const ArrayofMatrix& am);

void write_array_of_matrix_to_file(const string& filename,
                                   const ArrayofMatrix& am);

void read_array_of_matrix_from_stream(ArrayofMatrix& am,
                                      istream& is);

void read_array_of_matrix_from_file(ArrayofMatrix& am,
                                    const string& filename);



////////////////////////////////////////////////////////////////////////////
//   STRING IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

void write_array_of_string_to_stream(
              ostream&         os,
        const Arrayofstring&   as );

void write_array_of_string_to_file(
        const string&          filename,
        const Arrayofstring&   as );

void read_array_of_string_from_stream(
        Arrayofstring&   as,
        istream&         is );

void read_array_of_string_from_file(
           Arrayofstring&   as,
     const string&          filename );



////////////////////////////////////////////////////////////////////////////
//   Open and close binary
////////////////////////////////////////////////////////////////////////////

void binfile_open_out(
              int&      fid,
        const string&   filename );

void binfile_open_in(
              int&      fid,
        const string&   filename );

void binfile_close(
              int&      fid,
        const string&   filename );



////////////////////////////////////////////////////////////////////////////
//   Functions to read and write binary data for ARTS data types
////////////////////////////////////////////////////////////////////////////

void binfile_write_index(
        const string&   filename,
        const int&      fid,
        const size_t&   x,
        const string&   dataname );

void binfile_read_index(
              size_t&   x,
        const string&   filename,
        const int&      fid,
        const string&   dataname );

void binfile_write_numeric(
        const string&   filename,
        const int&      fid,
        const Numeric&  x,
        const string&   dataname );

void binfile_read_numeric(
              Numeric&  x,
        const string&   filename,
        const int&      fid,
        const string&   dataname );

void binfile_write_vector(
        const string&   filename,
        const int&      fid,
        const Vector&   x,
        const string&   dataname );

void binfile_read_vector(
              Vector&   x,
        const string&   filename,
        const int&      fid,
        const string&   dataname );

void binfile_write_matrix(
        const string&   filename,
        const int&      fid,
        const Matrix&   x,
        const string&   dataname );

void binfile_read_matrix(
              Matrix&   x,
        const string&   filename,
        const int&      fid,
        const string&   dataname );

void binfile_write_indexarray(
        const string&         filename,
        const int&            fid,
        const Arrayofsizet&   x,
        const string&         dataname );

void binfile_read_indexarray(
              Arrayofsizet&   x,
        const string&         filename,
        const int&            fid,
        const string&         dataname );

void binfile_write_vectorarray(
        const string&          filename,
        const int&             fid,
        const ArrayofVector&   x,
        const string&          dataname );

void binfile_read_vectorarray(
              ArrayofVector&   x,
        const string&          filename,
        const int&             fid,
	const string&          dataname );

void binfile_write_matrixarray(
        const string&          filename,
        const int&             fid,
        const ArrayofMatrix&   x,
        const string&          dataname );

void binfile_read_matrixarray(
              ArrayofMatrix&   x,
        const string&          filename,
        const int&             fid,
        const string&          dataname );

void binfile_write_string(
        const string&   filename,
        const int&      fid,
        const string&   s,
        const string&   dataname );

void binfile_read_string(
              string&   x,
        const string&   filename,
        const int&      fid,
        const string&   dataname );

void binfile_write_stringarray(
        const string&          filename,
        const int&             fid,
        const Arrayofstring&   x,
        const string&          dataname );

void binfile_read_stringarray(
              Arrayofstring&   x,
        const string&          filename,
        const int&             fid,
        const string&          dataname );


#endif
