/* Copyright (C) 2002
   Oliver Lemke <olemke@uni-bremen.de>

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
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   2002-06-18

  \brief  Workspace methods and template functions for supergeneric XML IO.

*/

#ifndef m_xml_h
#define m_xml_h

#include "xml_io.h"

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


template<typename T> void
WriteXML (//WS Input:
          const String &file_format,
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


void
output_file_formatSetAscii (// WS Output:
                            String &file_format);

void
output_file_formatSetBinary (// WS Output:
                             String &file_format);


#endif // m_xml_h

