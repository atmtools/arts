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

#if HAVE_CONFIG_H
#include "config.h"
#endif	

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


#ifdef HDF_SUPPORT

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



////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-12-06
*/
void ybatchCalc(
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
	      const ArrayOfString&              cont_description_models,
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
            Los   los;
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
	       cont_description_names, 
	       cont_description_models,
               cont_description_parameters);

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

#endif // HDF_SUPPORT

void ybatchFromRadiosonde(// WS Output:
                          Matrix& ybatch,
                          // WS Input:
			  const ArrayOfMatrix& radiosonde_data,
                          const Vector& f_mono,
                          const ArrayOfArrayOfLineRecord& lines_per_tg,
                          const ArrayOfLineshapeSpec& lineshape,
                          const Numeric& z_plat,
                          const Vector& za_pencil,
                          const Numeric& l_step,
                          const Index& refr,
			  const String& refr_model,
                          const Index& refr_lfac,
                          const Numeric& r_geoid,
                          const Index& emission,
                          const Vector& y_space,
                          const Vector& e_ground,
                          const TagGroups& tgs,
                          const ArrayOfString& cont_description_names,
			  const ArrayOfString& cont_description_models,
                          const ArrayOfVector& cont_description_parameters)
{
  // Initialize ybatch:
  ybatch.resize( f_mono.nelem()*za_pencil.nelem(), radiosonde_data.nelem() );
  ybatch = 0;

  // Loop over all radiosonde profiles:
  for ( Index i=0; i<radiosonde_data.nelem(); ++i )
    {
      
      out2 << "  Processing sonde profile " << i
	   << " of " << radiosonde_data.nelem() << "\n";

      const Matrix& rd = radiosonde_data[i];

      // Create p_abs, t_abs, z_abs:
      ConstVectorView p_abs = rd(Range(joker),0);
      ConstVectorView t_abs = rd(Range(joker),1);
      ConstVectorView z_abs = rd(Range(joker),2);

      // Create vmrs:
      Matrix vmrs(3,rd.nrows());
      vmrs(0,Range(joker)) = rd(Range(joker),3);     // H2O
      vmrs(1,Range(joker)) = 0.209;                  // O2 
      vmrs(2,Range(joker)) = 0.782;                  // N2

      // Set the physical H2O profile from the H2O profile in vmrs:
      Vector h2o_abs = vmrs(0,Range(joker));
      
      // Set the physical N2 profile from the N2 profile in vmrs:
      Vector n2_abs = vmrs(2,Range(joker));

      // Calculate absorption:
      Matrix        abs;
      ArrayOfMatrix abs_per_tg;
      // ... call the workspace method absCalc:
      absCalc(// Output:
             abs,
             abs_per_tg,
             // Input:		  
	     tgs,
             f_mono,
             p_abs,
             t_abs,
	     n2_abs,
	     h2o_abs,
             vmrs,
             lines_per_tg,
	     lineshape,
	     cont_description_names,
	     cont_description_models,
	     cont_description_parameters);

      // Set t_ground from lowest level of t_abs:
      Numeric t_ground = t_abs[0];

      // Set z_ground from lowest level of z_abs:
      Numeric z_ground = z_abs[0];

      // Calculate refractive index:
      Vector refr_index;
      refrCalc(// Output:
	       refr_index,
	       // Input:
	       p_abs,
	       t_abs,
	       h2o_abs,
	       refr,
	       refr_model);

      // Calculate the line of sight:
      Los los;
      Vector z_tan;
      losCalc(// Output:
	      los,
	      z_tan,
	      // Input:
	      z_plat,
	      za_pencil,
	      l_step,
	      p_abs,
	      z_abs,
	      refr,
	      refr_lfac,
	      refr_index,
	      z_ground,
	      r_geoid);

      // Calculate source:
      ArrayOfMatrix source;
      sourceCalc(// Output:
		 source,
		 // Input:
		 emission,
		 los,
		 p_abs,
		 t_abs,
		 f_mono);

      // Calculate transmittances:
      ArrayOfMatrix trans;
      transCalc(// Output:
		trans,
		//Input:
		los,
		p_abs,
		abs);

      // Calculate Spectrum:
      Vector y;
      yCalc(// Output:
	    y,
	    // Input:
	    emission,
	    los,
	    f_mono,
	    y_space,
	    source,
	    trans,
	    e_ground,
	    t_ground);

      // Assign y to this column of ybatch:
      ybatch(Range(joker),i) = y;
    }
}

void ybatchFromRadiosondeGlobal(// WS Output:
                          Matrix& ybatch,
                          // WS Input:
			  const ArrayOfMatrix& radiosonde_data,
                          const Vector& f_mono,
                          const ArrayOfArrayOfLineRecord& lines_per_tg,
                          const ArrayOfLineshapeSpec& lineshape,
                          const Numeric& z_plat,
                          const Vector& za_pencil,
                          const Numeric& l_step,
                          const Index& refr,
			  const String& refr_model,
                          const Index& refr_lfac,
                          const Numeric& r_geoid,
                          const Index& emission,
                          const Vector& y_space,
                          const Vector& e_ground,
                          const TagGroups& tgs,
                          const ArrayOfString& cont_description_names,
			  const ArrayOfString& cont_description_models,
                          const ArrayOfVector& cont_description_parameters)
{
  // Initialize ybatch:
  ybatch.resize( f_mono.nelem()*za_pencil.nelem(), radiosonde_data.nelem() );
  ybatch = 0;

  // Loop over all radiosonde profiles:
  for ( Index i=0; i<radiosonde_data.nelem(); ++i )
    {
      out2 << "  Processing sonde profile " << i
	   << " of " << radiosonde_data.nelem() << "\n";

      const Matrix& rd = radiosonde_data[i];

      // Create t_abs, z_abs:

      // For the interpolated profiles:
      
            
      ConstVectorView p_raw = rd(Range(joker),0);
      ConstVectorView z_raw = rd(Range(joker),1);
      ConstVectorView t_raw = rd(Range(joker),2);
      
      Vector p_abs;

      VectorNLogSpace (p_abs,"p_abs",p_raw[0],p_raw[p_raw.nelem()-1],100);
      Vector t_abs( p_abs.nelem() );
      Vector z_abs( p_abs.nelem() );
		       
      interpp (t_abs,p_raw,t_raw,p_abs);
      interpp (z_abs,p_raw,z_raw,p_abs);	
      
      // Create vmrs:
      
      ConstVectorView h20_raw = rd(Range(joker),4);
      
      Matrix vmrs(3,p_abs.nelem());
      interpp (vmrs(0,Range(joker)),p_raw,h20_raw,p_abs);
      //vmrs(0,Range(joker)) = rd(Range(joker),4);     // H2O
      vmrs(1,Range(joker)) = 0.209;                  // O2 
      vmrs(2,Range(joker)) = 0.782;                  // N2

      // Set the physical H2O profile from the H2O profile in vmrs:
      Vector h2o_abs = vmrs(0,Range(joker));
      
      // Set the physical N2 profile from the N2 profile in vmrs:
      Vector n2_abs = vmrs(2,Range(joker));

      // Calculate absorption:
      Matrix        abs;
      ArrayOfMatrix abs_per_tg;
      // ... call the workspace method absCalc:
      absCalc(// Output:
             abs,
             abs_per_tg,
             // Input:		  
	     tgs,
             f_mono,
             p_abs,
             t_abs,
	     n2_abs,
	     h2o_abs,
             vmrs,
             lines_per_tg,
	     lineshape,
	     cont_description_names,
	     cont_description_models,
	     cont_description_parameters);

      // Set t_ground from lowest level of t_abs:
      Numeric t_ground = t_raw[0];

      // Set z_ground from lowest level of z_abs:
      Numeric z_ground = z_raw[0];
      cout << z_raw[0] << endl;

      // Calculate refractive index:
      Vector refr_index;
      refrCalc(// Output:
	       refr_index,
	       // Input:
	       p_abs,
	       t_abs,
	       h2o_abs,
	       refr,
	       refr_model);


      // Calculate the line of sight:
      Los los;
      Vector z_tan;
      losCalc(// Output:
	      los,
	      z_tan,
	      // Input:
	      z_plat,
	      za_pencil,
	      l_step,
	      p_abs,
	      z_abs,
	      refr,
	      refr_lfac,
	      refr_index,
	      z_ground,
	      r_geoid);

      // Calculate source:
      ArrayOfMatrix source;
      sourceCalc(// Output:
		 source,
		 // Input:
		 emission,
		 los,
		 p_abs,
		 t_abs,
		 f_mono);

      // Calculate transmittances:
      ArrayOfMatrix trans;
      transCalc(// Output:
		trans,
		//Input:
		los,
		p_abs,
		abs);

      // Calculate Spectrum:
      Vector y;
      yCalc(// Output:
	    y,
	    // Input:
	    emission,
	    los,
	    f_mono,
	    y_space,
	    source,
	    trans,
	    e_ground,
	    t_ground);

      // Assign y to this column of ybatch:
      ybatch(Range(joker),i) = y;
    }
}

