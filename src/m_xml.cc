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
  \file   m_xml.cc
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   Fri May 10 11:22:29 2002
  
  \brief  Workspace functions for reading and writing data from/to XML files.
*/


////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include "arts.h"
#include "matpackI.h"
#include "xml_io.h"
#include "messages.h"
#include "auto_md.h"
#include "make_array.h"



////////////////////////////////////////////////////////////////////////////
//   The functions
////////////////////////////////////////////////////////////////////////////


//**************************************************************************
//
//  2. XML file functions
//
//**************************************************************************

////////////////////////////////////////////////////////////////////////////
//   File workspace methods (sorted after workspace variable)
////////////////////////////////////////////////////////////////////////////

// The XML functions are created by Oliver Lemke


//=== Index ============================================================

// This function shall be modified to handle Index
void
IndexWriteXML (// WS Generic Output:
               const Index&  v,
               // WS Generic Output Names:
               const String& v_name,
               // Control Parameters:
               const String& f)
{
  String filename = f;

  // Create default filename if empty  
  filename_xml (filename, v_name);

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}


// This function shall be modified to handle Index
void
IndexReadXML (// WS Generic Output:
              Index&        v,
              // WS Generic Output Names:
              const String& v_name,
              // Control Parameters:
              const String& f )
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, v_name );

  xml_read_from_file (filename, v);
}



//=== NUMERIC ==========================================================

void
NumericWriteXML (// WS Generic Output:
                 const Numeric& v,
                 // WS Generic Output Names:
                 const String&  v_name,
                 // Control Parameters:
                 const String&  f)
{
  String filename = f;

  // Create default filename if empty  
  filename_xml( filename, v_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}


void
NumericReadXML (// WS Generic Output:
                Numeric&      v,
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



//=== Vector ==========================================================

void
VectorWriteXML (// WS Generic Output:
                const Vector& v,
                // WS Generic Output Names:
                const String& v_name,
                // Control Parameters:
                const String& f)
{
  String filename = f;

  // Create default filename if empty  
  filename_xml( filename, v_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}


void
VectorReadXML (// WS Generic Output:
               Vector&       v,
               // WS Generic Output Names:
               const String& v_name,
               // Control Parameters:
               const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, v_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



//=== Matrix ==========================================================

void
MatrixWriteXML (// WS Generic Output:
                const Matrix& m,
                // WS Generic Output Names:
                const String& m_name,
                // Control Parameters:
                const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, m_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



void
MatrixReadXML (// WS Generic Output:
               Matrix& m,
               // WS Generic Output Names:
               const String& m_name,
               // Control Parameters:
               const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, m_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



//=== ArrayOfIndex =====================================================

// This function shall be modified to handle ArrayOfIndex
void
ArrayOfIndexWriteXML (// WS Generic Output:
                      const ArrayOfIndex& v,
                      // WS Generic Output Names:
                      const String&       v_name,
                      // Control Parameters:
                      const String&       f)
{
  String filename = f;

  // Create default filename if empty  
  filename_xml( filename, v_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}


// This function shall be modified to handle Index
void
ArrayOfIndexReadXML (// WS Generic Output:
                     ArrayOfIndex& v,
                     // WS Generic Output Names:
                     const String& v_name,
                     // Control Parameters:
                     const String& f)
{
  // FIXME: This function is crap. Put the whole ASCII file stuff
  // should be changed in the future, so I leave it for now.
  
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, v_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



//=== ArrayOfVector ====================================================

void
ArrayOfVectorWriteXML (// WS Generic Output:
                       const ArrayOfVector& av,
                       // WS Generic Output Names:
                       const String&        av_name,
                       // Control Parameters:
                       const String&        f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, av_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



void
ArrayOfVectorReadXML (// WS Generic Output:
                      ArrayOfVector& av,
                      // WS Generic Output Names:
                      const String&  av_name,
                      // Control Parameters:
                      const String&  f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, av_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



//=== ArrayOfMatrix ====================================================

void
ArrayOfMatrixWriteXML (// WS Generic Output:
                       const ArrayOfMatrix& am,
                       // WS Generic Output Names:
                       const String&        am_name,
                       // Control Parameters:
                       const String&        f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, am_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



void
ArrayOfMatrixReadXML (// WS Generic Output:
                      ArrayOfMatrix& am,
                      // WS Generic Output Names:
                      const String&  am_name,
                      // Control Parameters:
                      const String&  f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, am_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



//=== String ===============================================================

void
StringWriteXML (// WS Generic Output:
                const String& s,
                // WS Generic Output Names:
                const String& s_name,
                // Control Parameters:
                const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, s_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}


void
StringReadXML (// WS Generic Output:
               String& s,
               // WS Generic Output Names:
               const String& s_name,
               // Control Parameters:
               const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, s_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



//=== ArrayOfString ====================================================

void
ArrayOfStringWriteXML (// WS Generic Input:
                       const ArrayOfString& as,
                       // WS Generic Input Names:
                       const String& as_name,
                       // Control Parameters:
                       const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, as_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



void
ArrayOfStringReadXML (// WS Generic Output:
                      ArrayOfString& as,
                      // WS Generic Output Names:
                      const String& as_name,
                      // Control Parameters:
                      const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_xml( filename, as_name );

  {
    ostringstream os;
    os << "Not yet implemented";
    throw runtime_error(os.str());
  }

}



