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

/* Input / Output method functions. */

#include "arts.h"
#include "messages.h"
#include "file.h"
#include "vecmat.h"
#include "math_funcs.h"
#include "md.h"



//
//-------------------< Scalar initialization >--------------------------
//

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



//
//------------------------< Vector initialization >------------------------
//

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



//
//-----------< Array of Matrix and Array of Vector Write Methods >------------
//

void ArrayOfMatrixWriteToFile(// WS Generic Input:
			      const ARRAYofMATRIX& am,
			      // WS Generic Input Names:
			      const string& am_name,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+am_name+".am";
    }

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}

void ArrayOfVectorWriteToFile(// WS Output:
			      const ARRAYofVECTOR& av,
			      // WS Variable Names:
			      const string& av_name,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+av_name+".am";
    }

  // Convert the array of vector to an array of matrix:
  ARRAYofMATRIX am(av.dim());
  for (size_t i=0; i<av.dim(); ++i)
    {
      to_matrix(am[i],av[i]);
    }

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}


//
//--------------------< Matrix and Vector Write Methods >--------------------
//

void MatrixWriteToFile(// WS Generic Input:
		       const MATRIX& m,
		       // WS Generic Input Names:
		       const string& m_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+m_name+".am";
    }

  // Convert the matrix to an array of matrix:
  ARRAYofMATRIX am(1,m);

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}

void VectorWriteToFile(// WS Output:
		       const VECTOR& v,
		       // WS Variable Names:
		       const string& v_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+v_name+".am";
    }

  // Convert the vector to a matrix:
  MATRIX m;
  to_matrix(m,v);

  // Convert the matrix to an array of matrix:
  ARRAYofMATRIX am(1,m);

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



//
//-------------< Array of Matrix and Array of Vector Read Methods >--------
//

void ArrayOfMatrixReadFromFile(// WS Generic Output:
			       ARRAYofMATRIX& am,
			       // WS Generic Output Names:
			       const string& am_name,
			       // Control Parameters:
			       const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+am_name+".am";
    }

  // Read the array of matrix from the file:
  read_array_of_matrix_from_file(am,filename);
}

void ArrayOfVectorReadFromFile(// WS Generic Output:
			       ARRAYofVECTOR& av,
			       // WS Generic Output Names:
			       const string& av_name,
			       // Control Parameters:
			       const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+av_name+".am";
    }

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



//
//--------------------< Matrix and Vector Read Methods >--------------------
//

void MatrixReadFromFile(// WS Generic Output:
			MATRIX& m,
			// WS Generic Output Names:
			const string& m_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+m_name+".am";
    }

  // Read the array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to a matrix.
  if ( 1 != am.dim() )
    throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                        "but the dimension of the array is not 1.");

  m = am[0];
}

void VectorReadFromFile(// WS Generic Output:
			VECTOR& v,
			// WS Generic Output Names:
			const string& v_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+v_name+".am";
    }

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



//
//--------------------< Copy Methods >--------------------
//
void VectorCopy(
                      VECTOR&   y2,
                const string&   name_y2,
                const VECTOR&   y1,
                const string&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  y2 = y1;
}

void MatrixCopy(
                      MATRIX&   y2,
                const string&   name_y2,
                const MATRIX&   y1,
                const string&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  y2 = y1;
}



//
//--------------------< H matrix Methods >--------------------
//

/** Generic function to read H matrices.

    \retval h        H matrix
    \param  h_name   which H matrix (h_total or h_data)
    \param  f        filename (if empty, the default name is used)

    \author Patrick Eriksson 
    \date   2000-10-06 
*/
void HmatrixReadFromFile(// WS Generic Output:
			Hmatrix& h,
			// WS Generic Output Names:
			const string& h_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+"."+h_name+".am";
    }

  // Read the array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  h.issparse = 0;
  h.full     = am[0];
  // Remember to clear the other field
}


