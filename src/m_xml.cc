/*!
  \file   m_xml.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2013-03-04

  \brief  Non-template implementations of workspace methods for XML IO.

*/

#include "m_xml.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void output_file_formatSetAscii(  // WS Output:
    String& file_format,
    const Verbosity&) {
  file_format = "ascii";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void output_file_formatSetZippedAscii(  // WS Output:
    String& file_format,
    const Verbosity&) {
  file_format = "zascii";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void output_file_formatSetBinary(  // WS Output:
    String& file_format,
    const Verbosity&) {
  file_format = "binary";
}
