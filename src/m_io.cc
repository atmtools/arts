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

   The functions are of two types:
     1. Reading and writing to/from files.
     2. Creation by workspace method keywords and generic input

   \author Patrick Eriksson
   \date 2000-11-01 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "messages.h"
#include "file.h"
#include "md.h"
#include "math_funcs.h"
#include "atm_funcs.h"          
#include <hdf.h>



////////////////////////////////////////////////////////////////////////////
//   General help functions
////////////////////////////////////////////////////////////////////////////

// setup_covmatrix ////////////////////////////////////////////////////////
/** 
   Core function to set-up covariance matrices from definition data.

   The function creates a covariance matrix from standard deviations and
   correlation lengths defined at some points (kp). Values at intermediate
   points are obtained by linear interpolation. The definition data must 
   cover all points of the retrieval/error grid.

   The correlation between different points is controled by the variables
   corrfun, cutoff and clength. The following correlation functions are
   defined:

     0: no correlation 
     1: linearly decreasing down to zero (tenth function)  
     2: exponential
     3: gaussian 

   The values of cutoff and clength are ignored when corrfun = 0.
 
   The clength is the distance to the point where the correlation has declined
   to exp(-1), the correlation length. The mean value of the correlation length
   at the two points of interest is used.

   The covariance is set to 0 for correlations below the cut-off value.

   \retval  s         covaraince matrix
   \param   kg        grid for the retrieval/error quantity
   \param   corrfun   correlation function (see above)
   \param   cutoff    cut-off for correlation coefficients (see above)
   \param   kp        abscissa for sdev and clength
   \param   sdev      standard deviation at the values of kp
   \param   clength   correlation length at the values of kp

   \author Patrick Eriksson
   \date   2000-12-01
*/
void setup_covmatrix(
                   MATRIX&    s,
             const VECTOR&    kg,
             const size_t&    corrfun,
             const Numeric&   cutoff,
             const VECTOR&    kp,
             const VECTOR&    sdev,
             const VECTOR&    clength )
{
  const size_t   n = kg.dim();
        size_t   row, col;
        VECTOR   sd, cl;
        Numeric  c;          // correlation

  if ( sdev.dim() != clength.dim() )
    throw runtime_error("The standard deviation and correlation length vectors must have the same length.");

  if ( (min(kg)<min(kp)) || (max(kg)>max(kp)) )
    throw runtime_error("The data defining the covariance do not cover all retrieval/error points.");

  if ( cutoff >= 1 )
    throw runtime_error("The correlation cut-off must be < 1.");

  // Interpolate to get standard deviation and correlation length at
  // each point of k_grid
  interp_lin( sd, kp, sdev, kg );
  interp_lin( cl, kp, clength, kg );

  // Resize s and fill with 0
  s.newsize(n,n);
  s = 0;  

  // Diagonal matrix
  if ( corrfun == 0 )
  {
    for ( row=1; row<=n; row++ )
      s(row,row) = sd(row)*sd(row); 
  }

  // Linearly decreasing (tenth function)
  else if ( corrfun == 1 )
  {
    for ( row=1; row<=n; row++ )
    {
      for ( col=row; col<=n; col++ )
      {
        c = 1 - 2*(1-exp(-1))*abs(kg(row)-kg(col))/(cl(row)+cl(col));
        if ( (c>0) && (c>cutoff) )
	{
          s(row,col) = c*sd(row)*sd(col); 
          s(col,row) = s(row,col);
        } 
      }
    }
  }

  // Exponential
  else if ( corrfun == 2 )
  {
    for ( row=1; row<=n; row++ )
    {
      for ( col=row; col<=n; col++ )
      {
        c = exp(-abs(kg(row)-kg(col))/((cl(row)+cl(col))/2));
        if ( c > cutoff )
	{
          s(row,col) = c*sd(row)*sd(col); 
          s(col,row) = s(row,col);
        } 
      }
    }
  }

  // Gaussian
  else if ( corrfun == 3 )
  {
    for ( row=1; row<=n; row++ )
    {
      for ( col=row; col<=n; col++ )
      {
        c = exp( -1.0 * pow( (kg(row)-kg(col))/((cl(row)+cl(col))/2), 2 ) );
        if ( c > cutoff )
	{
          s(row,col) = c*sd(row)*sd(col); 
          s(col,row) = s(row,col);
        } 
      }
    }
  }

  // Unknown correlation function
  else
  {
    ostringstream os;
    os << "Unknown correlation function flag (" << corrfun << ").";
    throw runtime_error(os.str());
  }
}





//**************************************************************************
//
//  1. File functions
//
//**************************************************************************

////////////////////////////////////////////////////////////////////////////
//   File help functions
////////////////////////////////////////////////////////////////////////////

// These functions gives, if filename is empty, the default file name for
// the different file types.

//// filename_ascii ////////////////////////////////////////////////////////
/**
   Gives the default file name for the ASCII formats.

   The default name is only used if the file name is empty.

   \retval   filename     file name
   \param    varname      variable name

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void filename_ascii(
              string&  filename,
        const string&  varname )
{
  if ( "" == filename )
  {
    extern const string basename;                       
    filename = basename+"."+varname+".aa";
  }
}



//// filename_bin ////////////////////////////////////////////////////////////
/**
   Gives the default file name for the binary format.

   The default name is only used if the file name is empty.

   \retval   filename     file name
   \param    varname      variable name

   \author Patrick Eriksson              
   \date   2000-11-01
*/
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
//   File workspace methods (sorted after workspace variable)
////////////////////////////////////////////////////////////////////////////

// The ASCII functions are created by Stefan Buehler around 00????
// The binary functions are created by Patrick Eriksson around 001102

//=== INDEX ============================================================

// This function shall be modified to handle INDEX
void IndexWriteAscii(
        const int&      v,
        const string&   v_name,
	const string&   f )
{
  string filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ARRAYofMATRIX and write to file
  ARRAYofMATRIX am(1);
  am(1).newsize(1,1);
  am(1) = (Numeric) v;
  write_array_of_matrix_to_file(filename,am);
}


// This function shall be modified to handle INDEX
void IndexReadAscii(
	      int&      v,
        const string&   v_name,
	const string&   f )
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ARRAYofMATRIX and move to Numeric
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.dim()!=1) || (am(1).dim(1)!=1) || (am(1).dim(2)!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains not a single value.";
    throw runtime_error(os.str());
  }

  Numeric a = am(1)(1,1);

  if ( (a-floor(a)) != 0 )
    throw runtime_error("The value in the file is not an integer.");
  if ( a < 0 )
    throw runtime_error("The value in the file is negative.");

  v = (int) a;
}


// This function shall be modified to handle INDEX
void IndexWriteBinary(
	const int&      v,
        const string&   var_name,
	const string&   f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_index( filename, fid, v, "INDEX" );
  binfile_close( fid, filename );
}


// This function shall be modified to handle INDEX
void IndexReadBinary(
	      int&      v,
        const string&   var_name,
	const string&   f )
{
  size_t vtemp;  //To be removed
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_index( vtemp, filename, fid, "INDEX" );
  v = (int) vtemp;
  binfile_close( fid, filename );
}



//=== NUMERIC ==========================================================

void NumericWriteAscii(
        const Numeric&  v,
        const string&   v_name,
	const string&   f )
{
  string filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ARRAYofMATRIX and write to file
  ARRAYofMATRIX am(1);
  am(1).newsize(1,1);
  am(1) = v;
  write_array_of_matrix_to_file(filename,am);
}


void NumericReadAscii(
	      Numeric&  v,
        const string&   v_name,
	const string&   f )
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ARRAYofMATRIX and move to Numeric
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.dim()!=1) || (am(1).dim(1)!=1) || (am(1).dim(2)!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains not a single numeric value";
    throw runtime_error(os.str());
  }

  v = am(1)(1,1);
}



void NumericWriteBinary(
        const Numeric&  v,
        const string&   var_name,
	const string&   f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_numeric( filename, fid, v, "NUMERIC" );
  binfile_close( fid, filename );
}



void NumericReadBinary(
	      Numeric&  v,
        const string&   var_name,
	const string&   f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_numeric( v, filename, fid, "NUMERIC" );
  binfile_close( fid, filename );
}



//=== VECTOR ==========================================================

void VectorWriteAscii(// WS Output:
		       const VECTOR& v,
		       // WS Variable Names:
		       const string& v_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

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
  filename_ascii( filename, v_name );

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
  int    fid;
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
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vector( v, filename, fid, "VECTOR" );
  binfile_close( fid, filename );
}



//=== MATRIX ==========================================================

void MatrixWriteAscii(// WS Generic Input:
		       const MATRIX& m,
		       // WS Generic Input Names:
		       const string& m_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, m_name );

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
  filename_ascii( filename, m_name );

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
  int    fid;
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
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_matrix( v, filename, fid, "MATRIX" );
  binfile_close( fid, filename );
}



//=== ARRAYofINDEX =====================================================

// This function shall be modified to handle ARRAYofINDEX
void ArrayOfIndexWriteAscii(
        const ARRAYofsizet&   v,
        const string&         v_name,
	const string&         f )
{
  string filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ARRAYofMATRIX and write to file
  const size_t  n = v.dim();
  ARRAYofMATRIX am(1);
  am(1).newsize(n,1);
  for ( size_t i=1; i<=n; i++ )
    am(1)(i,1) = (Numeric) v(i);
  write_array_of_matrix_to_file(filename,am);
}


// This function shall be modified to handle INDEX
void ArrayOfIndexReadAscii(
	      ARRAYofsizet&   v,
        const string&         v_name,
	const string&         f )
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ARRAYofMATRIX 
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.dim()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains more than one vector.";
    throw runtime_error(os.str());
  }

  // Move values to v using a vector
  VECTOR x;
  to_vector( x, am(1) );
  const size_t  n = x.dim();
  v.newsize(n);
  for ( size_t i=1; i<=n; i++ )
  {  
    if ( (x(i)-floor(x(i))) != 0 )
      throw runtime_error("A value in the file is not an integer.");
    if ( x(i) < 0 )
      throw runtime_error("A value in the file is negative.");
    v(i) = (int) x(i);
  }
}



void ArrayOfIndexWriteBinary(
        const ARRAYofsizet&  v,
        const string&        var_name,
	const string&        f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_indexarray( filename, fid, v, "INDEXARRAY" );
  binfile_close( fid, filename );
}



void ArrayOfIndexReadBinary(
	      ARRAYofsizet&  v,
        const string&        var_name,
	const string&        f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_indexarray( v, filename, fid, "INDEXARRAY" );
  binfile_close( fid, filename );
}



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
  filename_ascii( filename, av_name );

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
  filename_ascii( filename, av_name );

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



void ArrayOfVectorWriteBinary(
        const ARRAYofVECTOR&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_vectorarray( filename, fid, v, "VECTOR" );
  binfile_close( fid, filename );
}



void ArrayOfVectorReadBinary(
	      ARRAYofVECTOR&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vectorarray( v, filename, fid, "VECTOR" );
  binfile_close( fid, filename );
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
  filename_ascii( filename, am_name );

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
  filename_ascii( filename, am_name );

  // Read the array of matrix from the file:
  read_array_of_matrix_from_file(am,filename);
}



void ArrayOfMatrixWriteBinary(
        const ARRAYofMATRIX&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_matrixarray( filename, fid, v, "MATRIX" );
  binfile_close( fid, filename );
}



void ArrayOfMatrixReadBinary(
	      ARRAYofMATRIX&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_matrixarray( v, filename, fid, "MATRIX" );
  binfile_close( fid, filename );
}



//=== STRING ===============================================================

void StringWriteAscii( // WS Generic Input:
		       const string& s,
		       // WS Generic Input Names:
		       const string& s_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, s_name );

  // Convert the string to an array of string:
  ARRAYofstring as(1);
  as[0] = s;

  // Write the array of matrix to the file.
  write_array_of_string_to_file(filename,as);
}



void StringReadAscii(   // WS Generic Output:
			string& s,
			// WS Generic Output Names:
			const string& s_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, s_name );

  // Read the array of matrix from the file:
  ARRAYofstring as;
  read_array_of_string_from_file(as,filename);

  // Convert the array of string to a string.
  if ( 1 != as.dim() )
   throw runtime_error("You tried to convert an array of string to a string,\n"
                       "but the dimension of the array is not 1.");

  s = as[0];
}



void StringWriteBinary(
        const string&  v,
        const string&  var_name,
	const string&  f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_string( filename, fid, v, "STRING" );
  binfile_close( fid, filename );
}



void StringReadBinary(
	      string&  v,
        const string&  var_name,
	const string&  f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_string( v, filename, fid, "STRING" );
  binfile_close( fid, filename );
}



//=== ARRAYofstring ====================================================

void ArrayOfStringWriteAscii( // WS Generic Input:
			      const ARRAYofstring& as,
			      // WS Generic Input Names:
			      const string& as_name,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, as_name );

  // Write the array of matrix to the file.
  write_array_of_string_to_file(filename,as);
}



void ArrayOfStringReadAscii(   // WS Generic Output:
			       ARRAYofstring& as,
			       // WS Generic Output Names:
			       const string& as_name,
			       // Control Parameters:
			       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, as_name );

  // Read the array of matrix from the file:
  read_array_of_string_from_file(as,filename);
}



void ArrayOfStringWriteBinary(
        const ARRAYofstring&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_stringarray( filename, fid, v, "STRING" );
  binfile_close( fid, filename );
}



void ArrayOfStringReadBinary(
	      ARRAYofstring&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_stringarray( v, filename, fid, "STRING" );
  binfile_close( fid, filename );
}



//=== MAYBESPARSE ====================================================

void HmatrixReadAscii(// WS Generic Output:
			Hmatrix& h,
			// WS Generic Output Names:
			const string& h_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
 
  // Create default filename if empty  
  filename_ascii( filename, h_name );

  // Read the array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  h.issparse = 0;
  h.full     = am[0];
  // Remember to clear the other field
}



//=== LOS ================================================================

void LosWriteBinary(
        const LOS&      los,
        const string&   var_name,
	const string&   f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_vectorarray( filename, fid, los.p, "LOS.P" );
  binfile_write_vector( filename, fid, los.l_step, "LOS.L_STEP" );
  binfile_write_indexarray( filename, fid, los.ground, "LOS.GROUND" );
  binfile_write_indexarray( filename, fid, los.start, "LOS.START" );
  binfile_write_indexarray( filename, fid, los.stop, "LOS.STOP" );
  binfile_close( fid, filename );
}



void LosReadBinary(
	      LOS&      los,
        const string&   var_name,
	const string&   f )
{
  int   fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vectorarray( los.p, filename, fid, "LOS.P" );
  binfile_read_vector( los.l_step, filename, fid, "LOS.L_STEP" );
  binfile_read_indexarray( los.ground, filename, fid, "LOS.GROUND" );
  binfile_read_indexarray( los.start, filename, fid, "LOS.START" );
  binfile_read_indexarray( los.stop, filename, fid, "LOS.STOP" );
  binfile_close( fid, filename );
}





//**************************************************************************
//
//   2. Creation by workspace method keywords
//
//**************************************************************************

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



//=== VECTOR ==========================================================

void VectorSet(           VECTOR&  x, 
                    const string&  x_name,
                    const int&     n,
                    const Numeric& value )
{
  x.newsize(n);
  x = value;
  out2 << "  Creating " << x_name << " as a constant vector\n"; 
  out3 << "         length: " << n << "\n";
  out3 << "          value: " << value << "\n";
}



void VectorSet2(          VECTOR&  x, 
                    const string&  x_name,
                    const VECTOR&  z,
                    const string&  z_name,
                    const Numeric& value )
{
  const size_t  n = z.dim();
  x.newsize(n);
  x = value;
  out2 << "  Creating " << x_name << " as a constant vector\n"; 
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
  out2 << "  Creating " << x_name << " as linearly spaced vector\n";
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
  out2 << "  Creating " << x_name << " as linearly spaced vector\n";
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
                    const int&     n )
{
  x = nlogspace(start,stop,n);
  out2 << "  Creating " << x_name << " as logarithmically spaced vector\n";
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



void VectorCopyFromArrayOfVector(
                      VECTOR&          v,
                const string&          v_name,
                const ARRAYofVECTOR&   va,
                const string&          va_name,
                const int&             i )
{
  if ( i < 1 )
    throw runtime_error("The index must be > 0.");
  if ( size_t(i) > va.dim() )
  {
    ostringstream os;
    os << "The vector array has only " << va.dim() << "elements.";
    throw runtime_error(os.str());
  }

  out2 << "  Copies " << v_name << " from vector " << i << " in " << va_name
                                                                   << "\n";
  v = va(i);
}



void VectorPlanck(
                    VECTOR&   y,
              const string&   y_name,
              const VECTOR&   f,
              const string&   f_name,
              const Numeric&  t )
{
  if ( t > 0 )
  {
    planck( y, f, t );
    out2<<"  Setting " << y_name << " to blackbody radiation for "<<t<<" K.\n";
  }
  else
    throw runtime_error("The temperature must be > 0.");
}



void VectorRandUniform(
                    VECTOR&   y,
              const string&   y_name,
              const Numeric&  x_low,
              const Numeric&  x_high,
              const int&      n )
{
  out2<<"  Filling " << y_name << " with uniform random data.\n";
  rand_uniform( y, n, x_low, x_high );
}



void VectorRandNormal(
                    VECTOR&   y,
              const string&   y_name,
              const Numeric&  stddev,
              const int&      n )
{
  out2<<"  Filling " << y_name << " with Gaussian random data.\n";
  rand_normal( y, n, stddev );
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



//=== STRING ===============================================================

void StringSet(           string&  s, 
                    const string&  s_name,
                    const string&  s2 )
{
  s = s2;
  out3 << "  Setting " << s_name << " to " << s2 << "\n"; 
}



//=== COVARIANCE MATRIX =====================================================

void sDiagonal(
                   MATRIX&    s,
             const VECTOR&    k_grid,
             const Numeric&   stddev)
{
  const VECTOR   sdev(2,stddev);
  const VECTOR   clength(2,0.0);
        VECTOR   kp(2);

  kp(1) = k_grid(1);
  kp(2) = k_grid(k_grid.dim());
  setup_covmatrix( s, k_grid, 0, 0, kp, sdev, clength );
}



void sSimple(
                   MATRIX&    s,
             const VECTOR&    k_grid,
             const int&       corrfun,
             const Numeric&   cutoff,
             const Numeric&   stddev,
             const Numeric&   corrlength )
{
  const VECTOR   sdev(2,stddev);
  const VECTOR   clength(2,corrlength);
        VECTOR   kp(2);

  kp(1) = k_grid(1);
  kp(2) = k_grid(k_grid.dim());
  setup_covmatrix( s, k_grid, corrfun, cutoff, kp, sdev, clength );
}



void sFromFile(
                   MATRIX&    s,
             const VECTOR&    k_grid,
             const string&    filename )
{
  ARRAYofMATRIX   am;
         VECTOR   kp, sdev, clength;
         size_t   i, j, np;
 
  // Read the array of matrix from the file:
  read_array_of_matrix_from_file(am,filename);

  const size_t n = am.dim()-1;

  // Some checks of sizes
  if ( n < 1 )
    throw runtime_error("The file must contain > 1 matrix.");
  if ( am(1).dim(1) != 2 )
    throw runtime_error("The first matrix in the file must have 2 rows.");
  if ( am(1).dim(2) != n )
    throw runtime_error("The number of columns of the first matrix must equal the number of matrices - 1.");
  
  // Loop the different covariance matrices
  for ( i=1; i<=n; i++ )
  {
    // Check if the corrfun flag is an integer
    if ( (am(1)(1,i)-floor(am(1)(1,i))) != 0 )
      throw runtime_error("The first row of matrix 1 shall only contain integers..");

    // Move definition values to vectors
    np = am(i+1).dim(1);
    kp.newsize(np);
    sdev.newsize(np);
    clength.newsize(np);
    for ( j=1; j<=np; j++ )
    {
      kp(j)      = am(i+1)(j,1);
      sdev(j)    = am(i+1)(j,2);
      clength(j) = am(i+1)(j,3);
    }

    if ( i == 1 )
      setup_covmatrix( s, k_grid, size_t(am(1)(1,i)), am(1)(2,i), kp, 
                                                              sdev, clength );
      
    else
    {
      MATRIX stmp;
      setup_covmatrix( stmp, k_grid, size_t(am(1)(1,i)), am(1)(2,i), kp, 
                                                              sdev, clength );
      s = s + stmp;
    }
  }
}



void CovmatrixInit(
              MATRIX&   s,
        const string&   s_name)
{
  s.newsize(0,0);
}



void sxAppend(
              MATRIX&   sx,
        const MATRIX&   s )
{
  const size_t   ns  = s.dim(1); 
  const size_t   nsx = sx.dim(1); 
        MATRIX   stmp;
        size_t   row, col;

  if ( (ns!=s.dim(2)) || (nsx!=sx.dim(2)) )
    throw runtime_error("A covariance matrix must be square.");
    
  stmp = sx;
  sx.newsize(nsx+ns,nsx+ns);
  sx = 0;
  for ( row=1; row<=nsx; row++ )
  {
    for ( col=1; col<=nsx; col++ )
      sx(row,col) = stmp(row,col);
  }
  for ( row=1; row<=ns; row++ )
  {
    for ( col=1; col<=ns; col++ )
      sx(nsx+row,nsx+col) = s(row,col);
  }
}



void sbAppend(
              MATRIX&   sb,
        const MATRIX&   s )
{
  sxAppend( sb, s );
}
