/* Copyright (C) 2000, 2001 Patrick Eriksson <patrick@rss.chalmers.se>
                            Stefan Buehler <sbuehler@uni-bremen.de>

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
   \file   m_covmatrix.cc

   This file contains functions associated with covariance matrices. 

   Basic math operations are found elsewhere.

   \author Patrick Eriksson
   \date 2000-12-05 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "arts.h"
#include "file.h"
#include "math_funcs.h"
#include "auto_md.h"
#include "messages.h"
#include "make_vector.h"


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

   Input vectors are passed as VectorViews, so this function can be
   called with vector or matrix selections.

   \author Patrick Eriksson
   \date   2000-12-01
*/
void setup_covmatrix(
		     Matrix&            s,
		     ConstVectorView    kg,
		     const Index        corrfun,
		     const Numeric      cutoff,
		     ConstVectorView    kp,
		     ConstVectorView    sdev,
		     ConstVectorView    clength )
{
  const Index   n = kg.nelem();
        Index   row, col;
        Vector   sd(n), cl(n);
        Numeric  c;          // correlation

  if ( sdev.nelem() != clength.nelem() )
    throw runtime_error("The standard deviation and correlation length vectors must have the same length.");

  if ( (min(kg)<min(kp)) || (max(kg)>max(kp)) )
    throw runtime_error("The data defining the covariance do not cover all retrieval/error points.");

  if ( cutoff >= 1 )
    throw runtime_error("The correlation cut-off must be < 1.");

  if ( !corrfun )
    out2 << "  Creating a diagonal covariance matrix.\n";
  else
  {
    out2 << "  Creating a simple covariance matrix.\n";
    if ( corrfun == 1 )
      out2 << "  Tenth function correlation.\n";
    else if ( corrfun == 2 )
      out2 << "  Exponential correlation.\n";
    else if ( corrfun == 3 )
      out2 << "  Gaussian correlation.\n";
    out3 << "    Correlation cut-off : " << cutoff << "\n";
  }
  out3 << "    Size                : " << n << "\n";

  // Interpolate to get standard deviation and correlation length at
  // each point of k_grid
  //interp_lin( sd, kp, sdev, kg );
  //interp_lin( cl, kp, clength, kg );
  interp_lin_vector( sd, kp, sdev, kg );
  interp_lin_vector( cl, kp, clength, kg );

  // Resize s and fill with 0
  s.resize( n, n);
  s = 0;			// Matpack can set all elements like this.

  // Diagonal matrix
  if ( corrfun == 0 )
  {
    for ( row=0; row<n; row++ )
      s(row,row) = sd[row]*sd[row]; 
  }

  // Linearly decreasing (tenth function)
  else if ( corrfun == 1 )
  {
    for ( row=0; row<n; row++ )
    {
      for ( col=row; col<n; col++ )
      {
        c = 1 - 2*(1-exp(-1))*fabs(kg[row]-kg[col])/(cl[row]+cl[col]);
	// FIXME: Patrick, I changed abs in the above statement to
	// fabs. Please check.
        if ( (c>0) && (c>cutoff) )
	{
          s(row,col) = c*sd[row]*sd[col]; 
        } 
      }
    }
  }

  // Exponential
  else if ( corrfun == 2 )
  {
    for ( row=0; row<n; row++ )
    {
      for ( col=row; col<n; col++ )
      {
        c = exp(-fabs(kg[row]-kg[col])/((cl[row]+cl[col])/2));
	// FIXME: Patrick, I changed abs in the above statement to
	// fabs. Please check.
        if ( c > cutoff )
	{
          s(row,col) = c*sd[row]*sd[col]; 
        } 
      }
    }
  }

  // Gaussian
  else if ( corrfun == 3 )
  {
    for ( row=0; row<n; row++ )
    {
      for ( col=row; col<n; col++ )
      {
        c = exp( -1.0 * pow( (kg[row]-kg[col])/((cl[row]+cl[col])/2), 2 ) );
        if ( c > cutoff )
	{
          s(row,col) = c*sd[row]*sd[col]; 
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





////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sDiagonal(
	       Matrix&          s,
	       const Index&     n,
	       const Numeric&   stddev)
{
  // We can use the special Vector class MakeVector to create
  // two-element Vectors with particular contents:
  const MakeVector   sdev   ( stddev, stddev );
  const MakeVector   clength( 0.0,    0.0    );	 
  const MakeVector   kp     ( 0.0,    1.0    );		 
					 
  Vector kg(n);			// Kreate kg,
  kg = 0.5;			// and fill all elements with 0.5.

  setup_covmatrix( s, kg, 0, 0, kp, sdev, clength );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sDiagonalLengthFromVector(
			       Matrix&          s,
			       const Vector&    grid,
			       const String&    grid_name,
			       const Numeric&   stddev)
{
  // We can use the special Vector class MakeVector to create
  // two-element Vectors with particular contents:
  const MakeVector   sdev(    stddev,  stddev               );
  const MakeVector   clength( 0.0,     0.0                  );
  const MakeVector   kp(      grid[0], grid[grid.nelem()-1] );

  out2 << "  Creating a diagonal covariance matrix.\n";
  out3 << "    Standard deviation = " << stddev << "\n";

  setup_covmatrix( s, grid, 0, 0, kp, sdev, clength );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sDiagonalLengthFromVectors(
				Matrix&          s,
				const Vector&    grid1,
				const Vector&    grid2,
				const String&    grid1_name,
				const String&    grid2_name,
				const Numeric&   stddev)
{
  // We can use the special Vector class MakeVector to create
  // two-element Vectors with particular contents:
  const MakeVector   sdev(    stddev, stddev );
  const MakeVector   clength( 0.0,    0.0    );
  const MakeVector   kp(      0.0,    1.0)   ;

  Vector kg(grid1.nelem()*grid2.nelem()); // Kreate kg,                     
  kg = 0.5;				  // and fill all elements with 0.5.

  setup_covmatrix( s, kg, 0, 0, kp, sdev, clength );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sSimple(
	     Matrix&       s,
	     const Index&       n,
	     const Numeric&   stddev,
	     const Index&       corrfun,
	     const Numeric&   cutoff,
	     const Numeric&   corrlength )
{
  // We can use the special Vector class MakeVector to create
  // two-element Vectors with particular contents:
  const MakeVector   sdev(    stddev,     stddev     );
  const MakeVector   clength( corrlength, corrlength );
  const MakeVector   kp(      1,          n          );

  Vector ls;
  linspace(ls,1.0,n,1.0);

  setup_covmatrix( s, ls, corrfun, cutoff, kp, sdev, clength);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sSimpleLengthFromVector(
           Matrix&    s,
        const Vector&    grid,
        const String&    grid_name,
        const Numeric&   stddev,
        const Index&       corrfun,
        const Numeric&   cutoff,
        const Numeric&   corrlength )
{
  // We can use the special Vector class MakeVector to create
  // two-element Vectors with particular contents:
  const MakeVector   sdev(    stddev,     stddev               );
  const MakeVector   clength( corrlength, corrlength           );
  const MakeVector   kp(      grid[0],    grid[grid.nelem()-1] );

  setup_covmatrix( s, grid, corrfun, cutoff, kp, sdev, clength );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sSimpleLengthFromVectors(
			      Matrix&    s,
			      const Vector&    grid1,
			      const Vector&    grid2,
			      const String&    grid1_name,
			      const String&    grid2_name,
			      const Numeric&   stddev,
			      const Index&       corrfun,
			      const Numeric&   cutoff,
			      const Numeric&   corrlength )
{
  // Get matrix for one repitition of the vector
  Matrix   s0;
  sSimpleLengthFromVector( s0, grid1, grid1_name, stddev, corrfun, cutoff, 
                                                                 corrlength );
  // Create the total matrix
  const Index   n1 = grid1.nelem(); 
  const Index   n2 = grid2.nelem();
  Index   row, col, i, i0;
  
  s.resize( n1*n2, n1*n2 );
  s = 0;			// Matpack can set all elements like this.
  for ( i=0; i<n2; i++ )
  {
    i0 = (i-1)*n1;
    for ( row=0; row<n1; row++ )
    {
      for ( col=row; col<n1; col++ )
      {
        s(i0+row,i0+col) = s0(row,col);
      }
    }
  } 
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sFromFile(
	       Matrix&    s,
	       const Vector&    grid,
	       const String&    grid_name,
	       const String&    filename )
{
  ArrayOfMatrix am;
  Vector        kp, sdev, clength;
  Index         i, j, np;
 
  // Read the array of matrix from the file:
  read_array_of_matrix_from_file(am,filename);

  const Index n = am.nelem()-1;

  // Some checks of sizes
  if ( n < 1 )
    throw runtime_error("The file must contain > 1 matrix.");
  if ( am[0].nrows() != 2 )
    throw runtime_error("The first matrix in the file must have 2 rows.");
  if ( am[0].ncols() != n )
    throw runtime_error("The number of columns of the first matrix must equal the number of matrices - 1.");

  out2 << "  Summing " << n << " matrices.\n";
  
  // Loop the different covariance matrices
  for ( i=1; i<=n; i++ )
  {
    // Check if the corrfun flag is an integer
    if ( (am[0](0,i-1)-floor(am[0](0,i-1))) != 0 )
      throw runtime_error("The first row of matrix 1 shall only contain integers..");

    // Move definition values to vectors
    np = am[i].nrows();
    kp.resize(np);
    sdev.resize(np);
    clength.resize(np);
    for ( j=0; j<np; j++ )
    {
      kp[j]      = am[i](j,0);
      sdev[j]    = am[i](j,1);
      clength[j] = am[i](j,2);
    }

    if ( i == 1 )
      setup_covmatrix( s, grid, Index(am[0](0,i-1)), am[0](1,i-1), kp, 
                                                              sdev, clength );
      
    else
    {
      Matrix stmp;
      setup_covmatrix( stmp, grid, Index(am[0](0,i-1)), am[0](1,i-1), kp, 
                                                              sdev, clength );
      s += stmp;		// Matpack can add element-vise like this.
    }
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void CovmatrixInit(
		   Matrix&   s,
		   const String&   s_name)
{
  out2 << " Initializes " << s_name << "\n";
  s.resize(0,0);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sxAppend(
              Matrix&         sx,
	      const Matrix&   s )
{
  const Index   ns  = s.nrows(); 
  const Index   nsx = sx.nrows(); 
  Matrix   stmp(sx);		// Matpack can initilalize a Matrix
				// with another Matrix.  
  Index   row, col;

  if ( (ns!=s.ncols()) || (nsx!=sx.ncols()) )
    throw runtime_error("A covariance matrix must be square.");
    
  sx.resize( nsx+ns, nsx+ns );
  sx = 0;			// Matpack can set all elements like this.
  for ( row=0; row<nsx; row++ )
  {
    for ( col=0; col<nsx; col++ )
      sx(row,col) = stmp(row,col);
  }
  for ( row=0; row<ns; row++ )
  {
    for ( col=0; col<ns; col++ )
      sx(nsx+row,nsx+col) = s(row,col);
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-01
*/
void sbAppend(
              Matrix&         sb,
	      const Matrix&   s )
{
  sxAppend( sb, s );
}
