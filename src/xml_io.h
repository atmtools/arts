/* Copyright (C) 2002 Oliver Lemke <olemke@uni-bremen.de>

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
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   Fri May 10 09:43:57 2002
  
  \brief This file contains basic functions to handle XML data files.
  
*/


#ifndef xml_io_h
#define xml_io_h

#include <fstream>
#include "matpackI.h"
#include "mystring.h"

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void
filename_xml (String&  filename,
              const String&  varname );


////////////////////////////////////////////////////////////////////////////
//   Functions to open and read XML files
////////////////////////////////////////////////////////////////////////////

void
open_output_xml (ofstream& file, const String& name);

void
open_input_xml (ifstream& file, const String& name);


////////////////////////////////////////////////////////////////////////////
//   Matrix/Vector IO routines for XML files
////////////////////////////////////////////////////////////////////////////

//void write_array_of_matrix_to_stream(ostream& os,
//                                     const ArrayOfMatrix& am);

void
xml_write_ArrayOfMatrix (const String&        filename,
                         const ArrayOfMatrix& am);

void
xml_read_ArrayOfMatrix (ArrayOfMatrix& am,
                        istream&       is);

void
xml_read_ArrayOfMatrix (ArrayOfMatrix& am,
                        const String&  filename);



////////////////////////////////////////////////////////////////////////////
//   STRING IO routines for XML files
////////////////////////////////////////////////////////////////////////////

void
xml_write_ArrayOfString (const String&        filename,
                         const ArrayOfString& as);

void
xml_read_ArrayOfString (ArrayOfString& as,
                        const String&  filename);


#endif
