
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

void open_input_file(ifstream& file, const String& name);

void read_text_from_stream(ArrayOfString& text, istream& is);

void read_text_from_file(ArrayOfString& text, const String& name);

void replace_all(String& s, const String& what, const String& with);



////////////////////////////////////////////////////////////////////////////
//   Matrix/Vector IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

void write_array_of_matrix_to_stream(ostream& os,
                                     const ArrayOfMatrix& am);

void write_array_of_matrix_to_file(const String& filename,
                                   const ArrayOfMatrix& am);

void read_array_of_matrix_from_stream(ArrayOfMatrix& am,
                                      istream& is);

void read_array_of_matrix_from_file(ArrayOfMatrix& am,
                                    const String& filename);



////////////////////////////////////////////////////////////////////////////
//   STRING IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

void write_array_of_String_to_stream(
              ostream&         os,
        const ArrayOfString&   as );

void write_array_of_String_to_file(
        const String&          filename,
        const ArrayOfString&   as );

void read_array_of_String_from_stream(
        ArrayOfString&   as,
        istream&         is );

void read_array_of_String_from_file(
           ArrayOfString&   as,
     const String&          filename );



////////////////////////////////////////////////////////////////////////////
//   Open and close binary
////////////////////////////////////////////////////////////////////////////

void binfile_open_out(
              int&      fid,
        const String&   filename );

void binfile_open_in(
              int&      fid,
        const String&   filename );

void binfile_close(
              int&      fid,
        const String&   filename );



////////////////////////////////////////////////////////////////////////////
//   Functions to read and write binary data for ARTS data types
////////////////////////////////////////////////////////////////////////////

void binfile_write_index(
        const String&   filename,
        const int&      fid,
        const Index&   x,
        const String&   dataname );

void binfile_read_index(
              Index&   x,
        const String&   filename,
        const int&      fid,
        const String&   dataname );

void binfile_write_numeric(
        const String&   filename,
        const int&      fid,
        const Numeric&  x,
        const String&   dataname );

void binfile_read_numeric(
              Numeric&  x,
        const String&   filename,
        const int&      fid,
        const String&   dataname );

void binfile_write_vector(
        const String&   filename,
        const int&      fid,
        const Vector&   x,
        const String&   dataname );

void binfile_read_vector(
              Vector&   x,
        const String&   filename,
        const int&      fid,
        const String&   dataname );

void binfile_write_matrix(
        const String&   filename,
        const int&      fid,
        const Matrix&   x,
        const String&   dataname );

void binfile_read_matrix(
              Matrix&   x,
        const String&   filename,
        const int&      fid,
        const String&   dataname );

void binfile_write_indexarray(
        const String&         filename,
        const int&            fid,
        const ArrayOfIndex&   x,
        const String&         dataname );

void binfile_read_indexarray(
              ArrayOfIndex&   x,
        const String&         filename,
        const int&            fid,
        const String&         dataname );

void binfile_write_vectorarray(
        const String&          filename,
        const int&             fid,
        const ArrayOfVector&   x,
        const String&          dataname );

void binfile_read_vectorarray(
              ArrayOfVector&   x,
        const String&          filename,
        const int&             fid,
	const String&          dataname );

void binfile_write_matrixarray(
        const String&          filename,
        const int&             fid,
        const ArrayOfMatrix&   x,
        const String&          dataname );

void binfile_read_matrixarray(
              ArrayOfMatrix&   x,
        const String&          filename,
        const int&             fid,
        const String&          dataname );

void binfile_write_String(
        const String&   filename,
        const int&      fid,
        const String&   s,
        const String&   dataname );

void binfile_read_String(
              String&   x,
        const String&   filename,
        const int&      fid,
        const String&   dataname );

void binfile_write_Stringarray(
        const String&          filename,
        const int&             fid,
        const ArrayOfString&   x,
        const String&          dataname );

void binfile_read_Stringarray(
              ArrayOfString&   x,
        const String&          filename,
        const int&             fid,
        const String&          dataname );


#endif
