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
   \file   m_io.cc

   This file contains IO workspace methods. 

   \author Patrick Eriksson
   \date 2000-11-01 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "messages.h"
#include "file.h"
#include "vecmat.h"
#include "math_funcs.h"
#include "md.h"



////////////////////////////////////////////////////////////////////////////
//   Help  functions
////////////////////////////////////////////////////////////////////////////

// These functions gives, if filename is empty, the default file name for
// the different file types.

void filename_num_ascii(
              string&  filename,
        const string&  varname )
{
  if ( "" == filename )
  {
    extern const string basename;                       
    filename = basename+"."+varname+".am";
  }
}

void filename_string_ascii(
              string&  filename,
        const string&  varname )
{
  if ( "" == filename )
  {
    extern const string basename;                       
    filename = basename+"."+varname+".as";
  }
}

void filename_bin(
              string&  filename,
        const string&  varname )
{
  if ( "" == filename )
  {
    extern const string basename;                       
    filename = basename+"."+varname+".ab";
  }
}




////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////

//=== INDEX ============================================================

void IntSet(// WS Generic Output:
            int& x,
            // WS Generic Output Names:
            const string& x_name,
            // Control Parameters:
            const int& value)
{
  x = value;
  out3 << "  Setting " << x_name << " to " << value << ".\n";
}



//=== NUMERIC ==========================================================

void NumericSet(// WS Generic Output:
                Numeric& x,
                // WS Generic Output Names:
                const string& x_name,
                // Control Parameters:
                const Numeric& value)
{
  x = value;
  out3 << "  Setting " << x_name << " to " << value << ".\n";
}



void NumericWriteBinary(
        const Numeric&  v,
        const string&  var_name,
	const string&  f )
{
  hid_t  fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_numeric( filename, fid, v, "NUMERIC" );
  binfile_close( fid, filename );
}



void NumericReadBinary(
	      Numeric&  v,
        const string&  var_name,
	const string&  f )
{
  hid_t  fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_numeric( v, filename, fid, "NUMERIC" );
  binfile_close( fid, filename );
}



//=== VECTOR ==========================================================

void VectorSet(           VECTOR&  x, 
                    const string&  x_name,
                    const int& n,
                    const Numeric& value )
{
  x.newsize(n);
  x = value;
  out3 << "  Creating " << x_name << " as a constant vector\n"; 
  out3 << "         length: " << n << "\n";
  out3 << "          value: " << value << "\n";
}



void VectorLinSpace(      VECTOR&  x, 
                    const string&  x_name,
                    const Numeric& start,
                    const Numeric& stop,
                    const Numeric& step )
{
  x = linspace(start,stop,step);
  out3 << "  Creating " << x_name << " as linearly spaced vector\n";
  out3 << "         length: " << x.size() << "\n";
  out3 << "    first value: " << x(1) << "\n";
  if ( x.size() > 1 )
  {
    out3 << "      step size: " << x(2)-x(1) << "\n";
    out3 << "     last value: " << x(x.size()) << "\n";
  }
}



void VectorNLinSpace(     VECTOR&  x, 
                    const string&  x_name,
                    const Numeric& start,
                    const Numeric& stop,
                    const int& n )
{
  x = nlinspace(start,stop,n);
  out3 << "  Creating " << x_name << " as linearly spaced vector\n";
  out3 << "         length: " << n << "\n";
  out3 << "    first value: " << x(1) << "\n";
  if ( x.size() > 1 )
  {
    out3 << "      step size: " << x(2)-x(1) << "\n";
    out3 << "     last value: " << x(x.size()) << "\n";
  }
}



void VectorNLogSpace(     VECTOR&  x, 
                    const string&  x_name,
                    const Numeric& start,
                    const Numeric& stop,
                    const int& n )
{
  x = nlogspace(start,stop,n);
  out3 << "  Creating " << x_name << " as logarithmically spaced vector\n";
  out3 << "         length: " << n << "\n";
  out3 << "    first value: " << x(1) << "\n";
  if ( x.size() > 1 )
    out3 << "     last value: " << x(x.size()) << "\n";
}



void VectorCopy(
                      VECTOR&   y2,
                const string&   name_y2,
                const VECTOR&   y1,
                const string&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  y2 = y1;
}



void VectorWriteAscii(// WS Output:
		       const VECTOR& v,
		       // WS Variable Names:
		       const string& v_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;

  // Create default filename if empty  
  filename_num_ascii( filename, v_name );

  // Convert the vector to a matrix:
  MATRIX m;
  to_matrix(m,v);

  // Convert the matrix to an array of matrix:
  ARRAYofMATRIX am(1,m);

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void VectorReadAscii(// WS Generic Output:
			VECTOR& v,
			// WS Generic Output Names:
			const string& v_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_num_ascii( filename, v_name );

  // Read an array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to a matrix.
  if ( 1 != am.dim() )
    throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                        "but the dimension of the array is not 1.");
  MATRIX m(am[0]);

  // Convert the matrix to a vector:
  to_vector(v,m);
}



void VectorWriteBinary(
        const VECTOR&  v,
        const string&  var_name,
	const string&  f )
{
  hid_t  fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_vector( filename, fid, v, "VECTOR" );
  binfile_close( fid, filename );
}



void VectorReadBinary(
	      VECTOR&  v,
        const string&  var_name,
	const string&  f )
{
  hid_t  fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vector( v, filename, fid, "VECTOR" );
  binfile_close( fid, filename );
}



//=== MATRIX ==========================================================

void MatrixCopy(
                      MATRIX&   y2,
                const string&   name_y2,
                const MATRIX&   y1,
                const string&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  y2 = y1;
}



void MatrixWriteAscii(// WS Generic Input:
		       const MATRIX& m,
		       // WS Generic Input Names:
		       const string& m_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_num_ascii( filename, m_name );

  // Convert the matrix to an array of matrix:
  ARRAYofMATRIX am(1,m);

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void MatrixReadAscii(// WS Generic Output:
			MATRIX& m,
			// WS Generic Output Names:
			const string& m_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_num_ascii( filename, m_name );

  // Read the array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to a matrix.
  if ( 1 != am.dim() )
    throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                        "but the dimension of the array is not 1.");

  m = am[0];
}



void MatrixWriteBinary(
        const MATRIX&  v,
        const string&  var_name,
	const string&  f )
{
  hid_t  fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_matrix( filename, fid, v, "MATRIX" );
  binfile_close( fid, filename );
}



void MatrixReadBinary(
	      MATRIX&  v,
        const string&  var_name,
	const string&  f )
{
  hid_t  fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_matrix( v, filename, fid, "MATRIX" );
  binfile_close( fid, filename );
}



//=== ARRAYofINDEX =====================================================



//=== ARRAYofVECTOR ====================================================

void ArrayOfVectorWriteAscii(// WS Output:
			      const ARRAYofVECTOR& av,
			      // WS Variable Names:
			      const string& av_name,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_num_ascii( filename, av_name );

  // Convert the array of vector to an array of matrix:
  ARRAYofMATRIX am(av.dim());
  for (size_t i=0; i<av.dim(); ++i)
    {
      to_matrix(am[i],av[i]);
    }

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void ArrayOfVectorReadAscii(// WS Generic Output:
			       ARRAYofVECTOR& av,
			       // WS Generic Output Names:
			       const string& av_name,
			       // Control Parameters:
			       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_num_ascii( filename, av_name );

  // Read an array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to an array of vector.
  av.newsize(am.dim());
  for (size_t i=0; i<am.dim(); ++i)
    {
      to_vector(av[i],am[i]);
    }
}



//=== ARRAYofMATRIX ====================================================

void ArrayOfMatrixWriteAscii(// WS Generic Input:
			      const ARRAYofMATRIX& am,
			      // WS Generic Input Names:
			      const string& am_name,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_num_ascii( filename, am_name );

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void ArrayOfMatrixReadAscii(// WS Generic Output:
			       ARRAYofMATRIX& am,
			       // WS Generic Output Names:
			       const string& am_name,
			       // Control Parameters:
			       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_num_ascii( filename, am_name );

  // Read the array of matrix from the file:
  read_array_of_matrix_from_file(am,filename);
}



//=== MAYBESPARSE ====================================================

/** Generic function to read H matrices.

    \retval h        H matrix
    \param  h_name   which H matrix (h_total or h_data)
    \param  f        filename (if empty, the default name is used)

    \author Patrick Eriksson 
    \date   2000-10-06 
*/
void HmatrixReadAscii(// WS Generic Output:
			Hmatrix& h,
			// WS Generic Output Names:
			const string& h_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_num_ascii( filename, h_name );

  // Read the array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  h.issparse = 0;
  h.full     = am[0];
  // Remember to clear the other field
}


