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
   \file   m_batch.cc

   This file contains for performing batch calculations. 

   \author Patrick Eriksson
   \date 2000-12-05 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "atm_funcs.h"
#include "absorption.h"
#include "math_funcs.h"
#include "messages.h"
#include "file.h"
#include "auto_md.h"



////////////////////////////////////////////////////////////////////////////
//   Help functions
////////////////////////////////////////////////////////////////////////////

//// filename_batch ////////////////////////////////////////////////////////
/**
   Composes name for batch files.

   \retval   filename     File name.
   \param    batchname    Default basename.
   \param    varname      Variable name.

   \author Patrick Eriksson              
   \date   2000-12-06
*/
void filename_batch(
              String&  filename,
        const String&  batchname,
        const String&  varname )
{
  if ( "" == filename )
    filename = batchname+"."+varname+".ab";
}



//// read_batchdata ////////////////////////////////////////////////////////
/**
   Reads batch data stored as a matrix and checks if the matrix has OK size.

   \retval   x          The read matrix.
   \param    batchname  Default basename.
   \param    filename   The filename specified in the control file.
   \param    varname    Variable name.
   \param    length     Expected length of data.
   \param    n          Number of spectra to calculate.

   \author Patrick Eriksson              
   \date   2000-12-06
*/
void read_batchdata( 
              Matrix&   x, 
        const String&   batchname, 
        const String&   filename, 
        const String&   varname,
        const Index&   length,
        const Index&   n )
{
  String fname = filename;
  filename_batch( fname, batchname, varname );
  MatrixReadBinary( x, "", fname );
  if ( x.nrows() != length )
  {
    ostringstream os;
    os << "The file " << fname << " contains data of length " << x.nrows() 
       << ", but a length of " << length << " is expected.\n";
      throw runtime_error(os.str());
  }
  if ( x.ncols() < n )
  {
    ostringstream os;
    os << "The file " << fname << " contains data for only " << x.ncols() 
       << " spectra when " << n << " spectra shall be calculated.\n";
      throw runtime_error(os.str());
  }    
  // Remove a possible excess of data
  if ( x.ncols() > n )
  {
    Matrix dummy(x(Range(joker),Range(0,n))); // Matpack can initialize a
					      // new Matrix from another
					      // Matrix. Used here in
					      // connection with
					      // Ranges to select a
					      // submatrix of the
					      // original Matrix.
    // Copy dummy back to x:
    x.resize( x.nrows(), n );
    x = dummy;			// Matpack can copy the contents of
				// matrices like this. The dimensions
				// must be the same! 
  }
}


// Part common for both versions of BatchdataGaussianTemperatureProfiles
void temperature_profiles(
              Matrix&   t,
              String&   fname,
        const Vector&   t_abs,
        const Matrix&   s,
        const String&   batchname,
        const Index&      n )
{
  fname = "";
  filename_batch( fname, batchname, "t_abs" );
  out2 << "  Creating " << n << " temperature profiles.\n";
  t.resize( t_abs.nelem(), n );

  // Check that s really is a square matrix:
  if ( s.ncols() != s.nrows() )
    throw runtime_error("The covariance matrix must be a square matrix."); 

  // Mean vector must be consistent with covariance matrix:
  if ( t_abs.nelem() != s.nrows() )
    throw runtime_error("The length of the mean vector and the size of the "
			"covariance matrix do not match.");  

  rand_data_gaussian( t, t_abs, s );
  MatrixWriteBinary( t, "", fname );
}



////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Carlos Jimenez
   \date   2001-26-01
*/
void BatchdataGaussianNoiseNoCorrelation(
     // WS Input
        const String&   batchname,
     // GInput
        const Vector&   grid1,
        const Vector&   grid2,
        const String&   grid1_name,
        const String&   grid2_name,
     // Control Parameters
        const Index&      n,
        const Numeric&  stddev)
{
  String fname = "";
  Matrix m(grid1.nelem()*grid2.nelem(),n);
  m = 0;			// Matpack can set all elements like this.

  filename_batch( fname,  batchname, "noise" );
  out2 << "  Creating " << n << " vectors with gaussian noise.\n";
  rand_matrix_gaussian( m, stddev);  
  MatrixWriteBinary(m,"",fname); 
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void BatchdataGaussianZeroMean(
      // WS Input:
        const String&   batchname,
        const Matrix&   s,
      // Control Parameters:
        const Index&      n,
        const String&   varname )
{
  const Index   l = s.nrows();
        String   fname = "";
        Matrix   zs;

  filename_batch( fname, batchname, varname );
  out2 << "  Creating " << n << " vectors with gaussian random data.\n";
  zs.resize(1,n);

  // Check that s really is a square matrix:
  if ( s.ncols() != s.nrows() )
    throw runtime_error("The covariance matrix must be a square matrix."); 

  // Create mean vector:
  Vector mean(l);
  mean = 0;

  rand_data_gaussian( zs, mean, s );
  MatrixWriteBinary( zs, "", fname );

}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void BatchdataGaussianTemperatureProfiles(
      // WS Input:
        const Vector&    p_abs,
        const Vector&    t_abs,
        const Vector&    z_abs,
        const Vector&    h2o_abs,
        const Matrix& s,
        const String&    batchname,
        const Numeric&   r_geoid,
        const Vector&    hse,
        const Index&       n )
{
  String fname;
  Matrix ts;
  temperature_profiles( ts, fname, t_abs, s, batchname, n );
  Matrix zs(t_abs.nelem(),n);
  Vector z, t;
  fname = "";
  filename_batch( fname, batchname, "z_abs" );
  out2 << "  Filling the file " << fname << "\n"
       << "  with vertical grids fulfilling hydrostatic eq.\n";

  t.resize( ts.nrows() );
  z.resize( z_abs.nelem() );
  for ( Index i=0; i<Index(n); i++ )
  {
    t = ts(Range(joker),i);	// Copy ith column of Matrix ts to
				// vector t.
    z = z_abs;			// Copy contents of Vector z_abs to
				// Vector z.

    hseCalc( z, p_abs, t, h2o_abs, r_geoid, hse );

    zs(Range(joker),i) = z;	// Copy contents of Vector z to ith
				// columns of Matrix zs.
  }

  MatrixWriteBinary( zs, "", fname );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void BatchdataGaussianTemperatureProfilesNoHydro(
        const Vector&     t_abs,
        const Matrix&  s,
        const String&     batchname,
        const Index&        n )
{
  String fname;
  Matrix t;
  temperature_profiles( t, fname, t_abs, s, batchname, n );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void BatchdataGaussianSpeciesProfiles(
      // WS Input:
        const TagGroups&       tgs,
        const Matrix&          vmrs,
        const Vector&          p_abs,
        const Vector&          t_abs,
        const Matrix&          s,
        const String&          batchname,
      // Control Parameters:
        const Index&           n,
        const ArrayOfString&   do_tags,
        const String&          unit )
{
  const Index   ntags = do_tags.nelem();    // Number of tags to do here 
  ArrayOfIndex  tagindex;                   // Index in tags for do_tags 
  String   	fname;
  Matrix   	x;
   
  // Check that s really is a square matrix:
  if ( s.ncols() != s.nrows() )
    throw runtime_error("The covariance matrix must be a square matrix."); 

  // Check if do_tags can be found in tag_groups and store indeces
  if ( ntags == 0 )
    throw runtime_error("No tags specified, no use in calling this function.");
  else
    get_tagindex_for_Strings( tagindex, tgs, do_tags );
  
  // Loop the tags
  for ( Index itag=0; itag<ntags; itag++ )
  {
    // Determine the name of the molecule for itag
    // The species lookup data:
    extern const Array<SpeciesRecord> species_data;
    String molname = species_data[tgs[tagindex[itag]][0].Species()].Name();

    // Create filename
    fname = "";
    filename_batch( fname, batchname, molname );

    out2 << "  Creating " << n << " profiles for " << molname << ".\n";

    // Make x the right size:
    x.resize( vmrs.ncols(), n );

    // Handle the different units
    if ( unit == "frac" )               // A priori fractions
    {
      Index   np = vmrs.ncols();
      Index   row, col;
      Numeric  a;

      // Create mean vector:
      Vector mean(np);
      mean = 1.0;

      // Check that mean vector is consisten with covariance matrix:
      if ( mean.nelem() != s.nrows() )
	throw runtime_error("The length of the mean vector and the size of the "
			    "covariance matrix do not match."); 

      rand_data_gaussian( x, mean, s );
      for ( row=0; row<np; row++ )
      {
        a = vmrs(tagindex[itag],row);
        for ( col=0; col<Index(n); col++ )
          x(row,col) *= a;
      }
    }

    else if ( unit == "vmr" )          // VMR
      {
	// Get handle on mean vector (for better readability):
	ConstVectorView mean = vmrs(tagindex[itag],Range(joker));
	// Range(joker) selects the entire row.    

	// Check that mean vector is consisten with covariance matrix:
	if ( mean.nelem() != s.nrows() )
	  throw runtime_error("The length of the mean vector and the size of the "
			      "covariance matrix do not match."); 

	rand_data_gaussian( x, mean, s );
      }
    else if ( unit == "nd" )           // Number density
    {
      Index   np = vmrs.ncols();
      Index   row, col;
      Numeric  a;

      // Create mean vector:
      Vector mean(np);
      mean = 0.0;

      // Check that mean vector is consisten with covariance matrix:
      if ( mean.nelem() != s.nrows() )
	throw runtime_error("The length of the mean vector and the size of the "
			    "covariance matrix do not match."); 

      rand_data_gaussian( x, mean, s );
      for ( row=0; row<np; row++ )
      {
        a = number_density ( p_abs[row], t_abs[row] );
        for ( col=0; col<Index(n); col++ )
          x(row,col) = vmrs(tagindex[itag],row) + x(row,col)/a;
      }
    }

    else
    {
      ostringstream os;
      os << "Selected unit not available: "<< unit <<"\n";
      throw runtime_error(os.str());
    }

    // Write to file
    MatrixWriteBinary( x, "", fname );
  }  
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void BatchdataGaussianOffSets(
      // WS Input:
        const String&  batchname,
      // WS Generic Input:
        const Vector&  z0,
      // WS Generic Input Names:
        const String&  z_name,
      // Control Parameters:
        const Index&     n,
        const Numeric& stddev)
{
  Vector   r(n);
  Matrix   a(1,n);

  String fname = "";
  filename_batch( fname, batchname, z_name );
  out2 << "  Creating a vector with " << n << " off-sets with gaussian PDF.\n";

  rand_gaussian( r, stddev );
  a(0,Range(joker)) = r;	// Copy Vector r to first row of
				// Matrix a.
  MatrixWriteBinary( a, "", fname );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void BatchdataUniformOffSets(
      // WS Input:
        const String&  batchname,
      // WS Generic Input:
        const Vector&  z0,
      // WS Generic Input Names:
        const String&  z_name,
      // Control Parameters:
        const Index&     n,
        const Numeric& low,
        const Numeric& high )
{
  Vector   r(n);
  Matrix   a(1,n);

  String fname = "";
  filename_batch( fname, batchname, z_name );
  out2 << "  Creating a vector with " << n << " off-sets with uniform PDF.\n";

  rand_uniform( r, low, high );
  a(0,Range(joker)) = r;	// Copy Vector r to first row of
				// Matrix a.
  MatrixWriteBinary( a, "", fname );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void BatchdataSinusoidalRippleNoCorrelations(
      // WS Input:
        const String&    batchname,
      // WS Generic Input:
        const Vector&    f,
        const Vector&    za,
      // WS Generic Input Names:
        const String&    f_name,
        const String&    za_name,
      // Control Parameters:
        const Index&       n,
        const Numeric&   period,
        const Numeric&   amplitude,
        const String&    pdf,
        const String&    varname )
{
  extern const Numeric PI;
  const Index   nza = za.nelem();
  const Index   nf  = f.nelem();
        Matrix   phase, amp;
        Matrix   rs(nf*nza,n);

  out2 << "  Creating " << n << " vectors sinusoidal baseline ripple.\n";

  // Make amp the right dize:
  amp.resize( nza, n );

  // Set amplitudes
  if ( pdf == "none" )
    amp = amplitude;		// Matpack can set all elements like this.
  else if ( pdf == "gaussian" )
    rand_matrix_gaussian( amp, amplitude );
  else if ( pdf == "uniform" )
    rand_matrix_uniform( amp, -amplitude, amplitude );
  else
  {
    ostringstream os;
    os << "The PDF " << pdf << " is not handled.\n";
      throw runtime_error(os.str());
  }    

  // Set phases:
  phase.resize( nza, n );
  rand_matrix_uniform( phase, 0, 2*PI );

  Index   i, iv, iza, if0;
  Vector   r(nf*nza);
  Numeric  a, p, b=2*PI/period;

  for ( i=0; i<Index(n); i++ )
  {
    for ( iza=0; iza<nza; iza++ )
    {
      if0 = iza*nf;
      a   = amp(iza,i);
      p   = phase(iza,i);
      for ( iv=0; iv<nf; iv++ )
        r[if0+iv] = a*sin(b*f[iv]+p);
    }
    //    put_in_col( rs, i+1, r );
    assert( r.nelem()==rs.nrows() );
    rs(Range(joker),i) = r;	// Copy Vector r to ith column of
				// Matrix rs.
  }

  // Write data 
  String fname = "";
  filename_batch( fname, batchname, varname );
  MatrixWriteBinary( rs, "", fname );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void ybatchAbsAndRte(
      // WS Output:
              Matrix&                     ybatch,
      // WS Input:
        const Vector&                     f_mono,
        const Vector&                     p_abs, 
        const Vector&                     t_abs,
	const Vector&                     n2_abs,
	const Vector&                     h2o_abs,
        const Matrix&                     vmrs,
        const ArrayOfArrayOfLineRecord&   lines_per_tag,
        const ArrayOfLineshapeSpec&       lineshape,
        const Vector&                     z_abs,
        const Numeric&                    z_plat,
        const Vector&                     za_pencil,
        const Numeric&                    l_step,
        const Index&                        refr,
        const Index&                        refr_lfac,
        const Vector&                     refr_index,
        const Numeric&                    z_ground,
        const Numeric&                    r_geoid,
	const Index&                        emission,
        const Vector&                     y_space,
        const Vector&                     e_ground,
        const Numeric&                    t_ground,
        const String&                     batchname,
        const TagGroups&                  tgs,
	const ArrayOfString&              cont_description_names,
	const ArrayOfVector& 		  cont_description_parameters,
      // Control Parameters:
        const Index&                        ncalc,
        const Index&                        do_t,
        const String&                     t_file,
        const Index&                        do_z,
        const String&                     z_file,
        const ArrayOfString&              do_tags,
        const ArrayOfString&              tag_files,
        const Index&                        do_f,
        const String&                     f_file,
        const Index&                        do_za,
        const String&                     za_file )
{
  const Index   np = p_abs.nelem();         // Number of pressure levels
  //  const Index   ntgs = tgs.nelem();         // Number of absorption tags
  const Index   ndo = do_tags.nelem();      // Number of tags to do here 
  ArrayOfIndex   tagindex;                  // Index in tgs for do_tags 
   
  if ( ndo != tag_files.nelem() )
  {
    ostringstream os;
    os << "There is " << ndo << " tga groups and only " << tag_files.nelem() 
       << " tag_files, even if they are empty they should match.\n";
      throw runtime_error(os.str());
  }

  // Check if do_tags can be found in tgs and store indeces
  if ( ndo > 0 )
    get_tagindex_for_Strings( tagindex, tgs, do_tags );

  out2 << "  Reading data from files.\n";

  // Read data from file(s) or use the workspace varaible
  
  // Temperature
  Vector t(t_abs);		// Matpack can initialize a
				// new Vector from another
				// Vector  
  Matrix Ts;
  if ( do_t )
    read_batchdata( Ts, batchname, t_file, "t_abs", np, ncalc );

  // Altitudes
  Vector z(z_abs);		// Matpack can initialize a
				// new Vector from another
				// Vector  

  Matrix Zs;
  if ( do_z )
    read_batchdata( Zs, batchname, z_file, "z_abs", np, ncalc );
  
  // Frequencies
  Vector f(f_mono);		// Matpack can initialize a
				// new Vector from another
				// Vector  

  Matrix f_oss;
  if ( do_f )
    read_batchdata( f_oss, batchname, f_file, "f_mono", 1, ncalc );
  
  // Zenith angles
  Vector za(za_pencil);		// Matpack can initialize a
				// new Vector from another
				// Vector  

  Matrix za_oss;
  if ( do_za )
    read_batchdata( za_oss, batchname, za_file, "za_pencil", 1, ncalc );
  
  // Species profiles
  Matrix vs(vmrs);		// Matpack can initalize a matrix from
				// another matrix.

  Index itag;
  ArrayOfMatrix VMRs(ndo);
  extern const Array<SpeciesRecord> species_data; // The species lookup data:
  for ( itag=0; itag<ndo; itag++ )
  {
    // Determine the name of the molecule for itag
    String molname = species_data[tgs[tagindex[itag]][0].Species()].Name();
    read_batchdata( VMRs[itag], batchname, tag_files[itag], molname, np, ncalc );
  }



  //--- Loop and calculate the spectra --------------------------------------
         Matrix   abs;
  ArrayOfMatrix   abs_per_tag;
            LOS   los;
  ArrayOfMatrix   trans, source;
         Vector   y, z_tan;

  out2 << "  Calculating spectra.\n";

  for ( Index i=0; i<Index(ncalc); i++ )
  {
    out2 << "  -------- Batch spectrum " << i << " --------\n";
    // Copy from file data for the "do" quantities
    if ( do_t )
      {
	assert( t.nelem()==Ts.nrows() );
	t = Ts(Range(joker),i);	// Copy ith column of Matrix Ts to
				// Vector t.
      }
    if ( do_z )
      {
	assert( z.nelem()==Zs.nrows() );
	z = Zs(Range(joker),i);	// copy ith column of Matrix Zs to
				// Vector z.
      }
    if ( do_f )
      {
	f = f_mono;		// Copy contents of f_mono to f.
	f += f_oss(0,i);	// Add f_oss(0,i) to all elements.
      }
    if ( do_za )
      {
	za = za_pencil;         // Copy contents of za_pencil to za.
	za += za_oss(0,i);      // Add za_oss(0,i) to all elements.
      }
    for ( itag=0; itag<ndo; itag++ )
      {
	assert( vs.ncols()==VMRs[itag].nrows() );
	vs(tagindex[itag],Range(joker)) = VMRs[itag](Range(joker),i); 
	// Set the row of vs for itag from ith column of VMRs[itag].
      }


    // Do the calculations
    //
    if ( (i==0) || do_t || ndo || do_f )
      absCalc( abs, abs_per_tag, tgs, f, p_abs, t, n2_abs, h2o_abs, vs, 
	       lines_per_tag, lineshape, 
	       cont_description_names, cont_description_parameters);

    if ( (i==0) || do_z || do_za )   
      losCalc( los, z_tan, z_plat, za, l_step, p_abs, z, refr, refr_lfac, 
                                               refr_index, z_ground, r_geoid );
    if ( (i==0) || do_t || do_z || do_f || do_za )   
      sourceCalc( source, emission, los, p_abs, t, f );
    transCalc( trans, los, p_abs, abs );
    yCalc ( y, emission, los, f, y_space, source, trans, e_ground, t_ground );
    
    // Move to ybatch
    if ( i == 0 )
      ybatch.resize( y.nelem(), ncalc );

    ybatch(Range(joker),i) = y;	// Copy to ith column of ybatch.
  }
  out2 << "  ------------------------------------\n";
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void ybatchLoadCalibration (
                    Matrix&   ybatch,
              const Vector&   i_cal1,
              const Vector&   i_cal2,
              const Vector&   y_cal1,
              const Vector&   y_cal2,
              const Vector&   za_sensor )
{
  Vector y(ybatch.nrows());
  for ( Index i=0; i<ybatch.ncols(); i++ )
  {
    y = ybatch(Range(joker),i);	// Copy ith column of ybatch to y.

    yLoadCalibration ( y, i_cal1, i_cal2, y_cal1, y_cal2, za_sensor );

    ybatch(Range(joker),i) = y;	// Copy to ith column of ybatch.
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void ybatchAdd (
                    Matrix&    ybatch,
              const String&    batchname,
              const String&    varname )
{
  const Index   l = ybatch.nrows();
  const Index   n = ybatch.ncols();
  Matrix x;
  read_batchdata( x, batchname, "", varname, l, n );
  ybatch += x;			// Matpack can add element-vise like this.
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void ybatchAddScaled (
                    Matrix&    ybatch,
              const String&    batchname,
              const String&    varname,
              const Numeric&   scalefac )
{
  const Index   l = ybatch.nrows();
  const Index   n = ybatch.ncols();
  Matrix x;
  read_batchdata( x, batchname, "", varname, l, n );

  // ybatch = ybatch + scalefac*x; 
  x *= scalefac;		
  ybatch += x;			// Matpack can perform these
				// element-vise operations.   
}
