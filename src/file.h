
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
   \file  file.cc

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
#include <hdf5.h>



////////////////////////////////////////////////////////////////////////////
//   Functions to open and read ASCII files
////////////////////////////////////////////////////////////////////////////

void open_output_file(ofstream& file, const string& name);

void open_input_file(ifstream& file, const string& name);

void read_text_from_stream(ARRAY<string>& text, istream& is);

void read_text_from_file(ARRAY<string>& text, const string& name);

void replace_all(string& s, const string& what, const string& with);



////////////////////////////////////////////////////////////////////////////
//   MATRIX/VECTOR IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

void write_array_of_matrix_to_stream(ostream& os,
                                     const ARRAYofMATRIX& am);

void write_array_of_matrix_to_file(const string& filename,
                                   const ARRAYofMATRIX& am);

void read_array_of_matrix_from_stream(ARRAYofMATRIX& am,
                                      istream& is);

void read_array_of_matrix_from_file(ARRAYofMATRIX& am,
                                    const string& filename);



////////////////////////////////////////////////////////////////////////////
//   Basic functions to handle HDF files
////////////////////////////////////////////////////////////////////////////

void binfile_open_out(
              hid_t&    fid,
        const string&   filename );

void binfile_open_in(
              hid_t&    fid,
        const string&   filename );

void binfile_close(
              hid_t&    fid,
        const string&   filename );



////////////////////////////////////////////////////////////////////////////
//   Functions to read and write binary data for ARTS data types
////////////////////////////////////////////////////////////////////////////

void binfile_write_size_t(
        const string&   filename,
        const hid_t&    fid,
        const size_t&   x,
        const string&   dataname );

void binfile_read_size_t(
              size_t&   x,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname );

void binfile_write_numeric(
        const string&   filename,
        const hid_t&    fid,
        const Numeric&  x,
        const string&   dataname );

void binfile_read_numeric(
              Numeric&  x,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname );

void binfile_write_vector(
        const string&   filename,
        const hid_t&    fid,
        const VECTOR&   x,
        const string&   dataname );

void binfile_read_vector(
              VECTOR&   x,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname );

void binfile_write_matrix(
        const string&   filename,
        const hid_t&    fid,
        const MATRIX&   x,
        const string&   dataname );

void binfile_read_matrix(
              MATRIX&   x,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname );


#endif
