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

/*!
  \file   m_xml.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   2002-06-18

  \brief  Workspace methods and template functions for supergeneric XML IO.

*/

#ifndef m_xml_h
#define m_xml_h

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#include "agenda_class.h"
#include "exceptions.h"
#include "workspace_ng.h"
#include "xml_io.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadXML(  // WS Generic Output:
    T& v,
    // WS Generic Output Names:
    const String& v_name,
    // WS Generic Input:
    const String& f,
    // WS Generic Input Names:
    const String& f_name _U_) {
  String filename = f;

  // Create default filename if empty
  filename_xml(filename, v_name);

  xml_read_from_file(filename, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadXML(Workspace& ws _U_,
             // WS Generic Output:
             T& v,
             // WS Generic Output Names:
             const String& v_name,
             // WS Generic Input:
             const String& f,
             // WS Generic Input Names:
             const String& f_name) {
  ReadXML(v, v_name, f, f_name);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadXMLIndexed(  // WS Generic Output:
    T& v,
    // WS Generic Output Names:
    const String& v_name,
    // WS Input:
    const Index& file_index,
    // WS Generic Input:
    const String& f,
    const Index& digits,
    // WS Generic Input Names:
    const String& f_name _U_,
    const String& digits_name _U_) {
  String filename = f;

  // Create default filename if empty
  filename_xml_with_index(filename, file_index, v_name, digits);

  xml_read_from_file(filename, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadXMLIndexed(Workspace& ws _U_,
                    // WS Generic Output:
                    T& v,
                    // WS Generic Output Names:
                    const String& v_name,
                    // WS Input:
                    const Index& file_index,
                    // WS Generic Input:
                    const String& f,
                    const Index& digits,
                    // WS Generic Input Names:
                    const String& f_name,
                    const String& digits_name) {
  ReadXMLIndexed(
      v, v_name, file_index, f, digits, f_name, digits_name);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteXML(  //WS Input:
    const String& file_format,
    // WS Generic Input:
    const T& v,
    const String& f,
    const Index& no_clobber,
    // WS Generic Input Names:
    const String& v_name,
    const String& f_name _U_,
    const String& no_clobber_name _U_)

{
  // If MPI is enabled make sure only master process performs the write.
#ifdef ENABLE_MPI
  int initialized;
  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init(nullptr, nullptr);
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return;
  }
#endif  // ENABLE_MPI

  String filename = f;

  // Create default filename if empty
  filename_xml(filename, v_name);

  const FileType ftype = string2filetype(file_format);

  String errmsg;

#pragma omp critical(WriteXML_critical_region)
  {
    try {
      xml_write_to_file(filename, v, ftype, no_clobber);
    } catch (const std::exception& e) {
      errmsg = e.what();
    }
  }

  ARTS_USER_ERROR_IF (errmsg.length(), errmsg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteXML(Workspace& ws _U_,
              //WS Input:
              const String& file_format,
              // WS Generic Input:
              const T& v,
              const String& f,
              const Index& no_clobber,
              // WS Generic Input Names:
              const String& v_name,
              const String& f_name,
              const String& no_clobber_name) {
  WriteXML(file_format,
           v,
           f,
           no_clobber,
           v_name,
           f_name,
           no_clobber_name);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteXMLIndexed(  //WS Input:
    const String& file_format,
    const Index& file_index,
    // WS Generic Input:
    const T& v,
    const String& f,
    const Index& digits,
    // WS Generic Input Names:
    const String& v_name,
    const String& f_name,
    const String& digits_name _U_) {
  String filename = f;

  // Create default filename if empty
  filename_xml_with_index(filename, file_index, v_name, digits);

  WriteXML(file_format, v, filename, 0, v_name, f_name, "");
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteXMLIndexed(Workspace& ws _U_,
                     //WS Input:
                     const String& file_format,
                     const Index& file_index,
                     // WS Generic Input:
                     const T& v,
                     const String& f,
                     const Index& digits,
                     // WS Generic Input Names:
                     const String& v_name,
                     const String& f_name,
                     const String& digits_name) {
  WriteXMLIndexed(file_format,
                  file_index,
                  v,
                  f,
                  digits,
                  v_name,
                  f_name,
                  digits_name);
}

#endif  // m_xml_h
