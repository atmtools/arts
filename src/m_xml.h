/* Copyright (C) 2002-2007 Oliver Lemke <olemke@core-dump.info>

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

#include "xml_io.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
ReadXML (// WS Generic Output:
         T&            v,
         // WS Generic Output Names:
         const String& v_name,
         // Control Parameters:
         const String& f)
{
  String filename = f;

  // Create default filename if empty
  filename_xml (filename, v_name);

  xml_read_from_file (filename, v);
}


/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
WriteXML (//WS Input:
          const String& file_format,
          // WS Generic Output:
          const T&            v,
          // WS Generic Output Names:
          const String& v_name,
          // Control Parameters:
          const String& f)
{
  String filename = f;

  // Create default filename if empty
  filename_xml (filename, v_name);

  xml_write_to_file
    (filename, v, file_format == "ascii" ? FILE_TYPE_ASCII : FILE_TYPE_BINARY);
}


/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
WriteXMLIndexed (//WS Input:
          const String& file_format,
          const Index&  file_index,
          // WS Generic Output:
          const T&            v,
          // WS Generic Output Names:
          const String& v_name,
          // Control Parameters:
          const String& f)
{
  String filename = f;

  // Create default filename if empty
  filename_xml_with_index( filename, file_index, v_name );

  WriteXML( file_format, v, v_name, filename );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
output_file_formatSetAscii (// WS Output:
                            String& file_format);


/* Workspace method: Doxygen documentation will be auto-generated */
void
output_file_formatSetBinary (// WS Output:
                             String& file_format);


#endif // m_xml_h

