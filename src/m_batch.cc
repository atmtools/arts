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
    x = dummy;                  // Matpack can copy the contents of
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
        const ArrayOfVector&              cont_description_parameters,
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
  Vector t(t_abs);              // Matpack can initialize a
                                // new Vector from another
                                // Vector  
  Matrix Ts;
  if ( do_t )
    read_batchdata( Ts, batchname, t_file, "t_abs", np, ncalc );

  // Altitudes
  Vector z(z_abs);              // Matpack can initialize a
                                // new Vector from another
                                // Vector  

  Matrix Zs;
  if ( do_z )
    read_batchdata( Zs, batchname, z_file, "z_abs", np, ncalc );
  
  // Frequencies
  Vector f(f_mono);             // Matpack can initialize a
                                // new Vector from another
                                // Vector  

  Matrix f_oss;
  if ( do_f )
    read_batchdata( f_oss, batchname, f_file, "f_mono", 1, ncalc );
  
  // Zenith angles
  Vector za(za_pencil);         // Matpack can initialize a
                                // new Vector from another
                                // Vector  

  Matrix za_oss;
  if ( do_za )
    read_batchdata( za_oss, batchname, za_file, "za_pencil", 1, ncalc );
  
  // Species profiles
  Matrix vs(vmrs);              // Matpack can initalize a matrix from
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
        t = Ts(Range(joker),i); // Copy ith column of Matrix Ts to
                                // Vector t.
      }
    if ( do_z )
      {
        assert( z.nelem()==Zs.nrows() );
        z = Zs(Range(joker),i); // copy ith column of Matrix Zs to
                                // Vector z.
      }
    if ( do_f )
      {
        f = f_mono;             // Copy contents of f_mono to f.
        f += f_oss(0,i);        // Add f_oss(0,i) to all elements.
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

    ybatch(Range(joker),i) = y; // Copy to ith column of ybatch.
  }
  out2 << "  ------------------------------------\n";
}

#endif // HDF_SUPPORT

void ybatchFromRadiosonde(// WS Output:
                          Matrix& ybatch,
                          ArrayOfMatrix& absbatch,
                          ArrayOfMatrix& jacbatch,
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
                          const ArrayOfVector& cont_description_parameters,
			  //Keyword
			  const Index& fine_abs_grid,
                          const Index& interpolation_in_rh,
                          const Index& za_batch,
                          const Index& e_ground_batch,
                          const Index& calc_abs,
                          const Index& calc_jac)
{
  // Check value of the keyword
  check_if_bool(fine_abs_grid, "Finegrid keyword" );
  check_if_bool(interpolation_in_rh, "Interpolation in RH keyword" );
  check_if_bool(za_batch, "za_batch keyword" );
  check_if_bool(e_ground_batch, "e_ground_batch keyword" );
  check_if_bool(calc_abs, "calc_abs keyword" );
  check_if_bool(calc_jac, "calc_jac keyword" );

  // this variable is keep the original za_pencil 
  Vector za_pencil_profile;
  Vector e_ground_profile;

  if (calc_abs) 
    {
      absbatch.resize( radiosonde_data.nelem() );
    }
  else
    {
      absbatch.resize( 1 );
      absbatch[0].resize( 1, 1 );
      absbatch[0] = -1.0;
    }
  
  if (calc_jac) 
    {
      jacbatch.resize( radiosonde_data.nelem() );
    }
  else
    {
      jacbatch.resize( 1 );
      jacbatch[0].resize( 1, 1 );
      jacbatch[0]  =  -1.0;
    }

  if (!za_batch)
    {
      // Initialize ybatch:
      ybatch.resize( f_mono.nelem()*za_pencil.nelem(), radiosonde_data.nelem() );
      ybatch = 0;
      za_pencil_profile.resize( za_pencil.nelem() );
      za_pencil_profile = za_pencil;
    }
  //The za_pencil vector now contains one za_pencil for each profile   
  else
    {
      //Check whether the number of matrices and the number of za_pencil
      //are the same
      if ( radiosonde_data.nelem() != za_pencil.nelem() )
        {
          ostringstream os;
          os << "The number of Radiosonde profiles is " << radiosonde_data.nelem() << "\n" 
             << "The number of zenith angles given is " << za_pencil.nelem() << "\n"
             << "But these two are expected to be the same.\n";
          throw runtime_error(os.str());
        }
      // Initialize ybatch: since there can be only one angle for each profile 
      // the size of ybatch does not depend on the angles.
      ybatch.resize( f_mono.nelem(), radiosonde_data.nelem() );
      ybatch = 0;
      za_pencil_profile.resize(1);

    }

  if (e_ground_batch)
    {
      if ( radiosonde_data.nelem() != e_ground.nelem() )
        {
          ostringstream os;
          os << "The number of Radiosonde profiles is " << radiosonde_data.nelem() << "\n" 
             << "The number of emissivities given is  " << e_ground.nelem() << "\n"
             << "But these two are expected to be the same, when e_ground_per_profile = 1 .\n";
          throw runtime_error(os.str());
        }
    }
  else
    {
      e_ground_profile.resize( f_mono.nelem() );
      e_ground_profile  =  e_ground;
    }
  
  // Loop over all radiosonde profiles:
  for ( Index i=0; i<radiosonde_data.nelem(); ++i )
    {
      const Matrix& rd = radiosonde_data[i];
      
      // When za_batch is set, it should be extracted from 
      // the za_pencil for each profile
      if (za_batch)
        {
          za_pencil_profile[0] = za_pencil[i];
        } 
    
      // When e_ground_batch is set, it should be extracted from 
      // the e_gound for each profile
      if (e_ground_batch)
        {
          e_ground_profile.resize( f_mono.nelem() );
          e_ground_profile  =  e_ground[i];
        }
      
      // Check whether the launch has reached upto 100 hpa
      Numeric min_p = rd(rd.nrows() - 1, 0);

      // The launch reached up to or above 100 hPa
      if (min_p <= 10000.0)
	{
	  Numeric t_ground, z_ground;
	  Vector  p_abs, t_abs, z_abs, h2o_abs, n2_abs;
	  Matrix  vmrs;
	  
	  // The absorption is calculated on a grid which is
	  // the same as the radiosonde levels
	  if (!fine_abs_grid)
	    {	  
	      // Create p_abs, t_abs, z_abs:
	      p_abs.resize(rd.nrows());
	      t_abs.resize(rd.nrows());
	      z_abs.resize(rd.nrows());
	      
	      p_abs = rd(Range(joker),0);
	      t_abs = rd(Range(joker),1);
	      z_abs = rd(Range(joker),2);
	    
              // FIXME: ARTS-RTTOV paper
              if (calc_abs)
                {
                  absbatch[i].resize( f_mono.nelem(), p_abs.nelem() );
                  absbatch[i] = 0;
                }
	      
	      // Create vmrs:
	      vmrs.resize(3, rd.nrows());
	      vmrs(0,Range(joker)) = rd(Range(joker),3);     // H2O
	      vmrs(1,Range(joker)) = 0.209;                  // O2 
	      vmrs(2,Range(joker)) = 0.782;                  // N2
	      
	      
	      // Set the physical H2O profile from the H2O profile in vmrs:
	      h2o_abs.resize(rd.nrows());
	      h2o_abs = vmrs(0,Range(joker));
	      
	      // Set the physical N2 profile from the N2 profile in vmrs:
	      n2_abs.resize(rd.nrows());
	      Vector n2_abs = vmrs(2,Range(joker));
	      
	      // Set t_ground from lowest level of t_abs:
	      t_ground = t_abs[0];
	      
	      // Set z_ground from lowest level of z_abs:
	      z_ground = z_abs[0];
	      
	    }
	  
	  // The absorption is calculated on a very fine grid
	  else
	    {
	      // Number of levels in the launch
	      Index n_rows = rd.nrows();
	      
	      // Create p_abs, t_abs, z_abs:
	      Vector p_raw(n_rows + 1);
	      Vector t_raw(n_rows + 1);
	      Vector z_raw(n_rows + 1);
	      
	      // Adding one more level at the bottom 
	      // Pressure is set as 1040 hpa
	      // Temperature is the lowermost value from the sonde
	      // Height is -99.0 m (a dummy value below ground)
	      p_raw[0] = 104000.0;
	      t_raw[0] = rd(0,1);
	      z_raw[0] = -99.0;
	      
	      p_raw[Range(1, n_rows)] = rd(Range(joker),0);
	      t_raw[Range(1, n_rows)] = rd(Range(joker),1);
	      z_raw[Range(1, n_rows)] = rd(Range(joker),2);
	      
	      // Creating p_abs 
	      VectorNLogSpace (p_abs,"p_abs", p_raw[0], 10000.0, 1001);
	      
	      t_abs.resize(p_abs.nelem());
	      z_abs.resize(p_abs.nelem());
	      
              // Interpolating the profiles on p_abs
	      interpp (t_abs, p_raw, t_raw, p_abs);
	      interpp (z_abs, p_raw, z_raw, p_abs);
              
              // Create vmrs:
	      Vector vmr_raw(n_rows + 1); 
	      
	      vmr_raw[0] = rd(0,3);
	      vmr_raw[Range(1, n_rows)] = rd(Range(joker),3);   // H2O	 
              
              vmrs.resize(3, p_abs.nelem());

              // Checking whether the interpolation should be done 
              // in vmr or RH
              if (interpolation_in_rh)
                {
                  // Calculates RH and interpolating RH on p_abs grid and 
                  // converting back to H2O VMR. This is done because the 
                  // interpolation in RH gives more realistic Tbs than
                  // interpolation VMR
                  Vector sat_pres_raw(n_rows + 1);
                  Vector sat_pres_abs(p_abs.nelem());
                  Vector rh_raw(n_rows + 1);
                  Vector rh_abs(p_abs.nelem());
                  
                  assert( sat_pres_raw.nelem() == t_raw.nelem() );
                  assert( sat_pres_abs.nelem() == t_abs.nelem() );                  
                  
                  e_eq_water(sat_pres_raw, t_raw);
                  e_eq_water(sat_pres_abs, t_abs);
                  
                  // Calculates RH for the raw profile
                  for ( Index j=0; j<rh_raw.nelem(); ++j )
                    {
                      rh_raw[j]  =  vmr_raw[j] / sat_pres_raw[j];
                    }
                  
                  // Interpolates raw RH on p_abs grid
                  interpp (rh_abs, p_raw, rh_raw, p_abs);
                  
                  // Converts RH back to VMR
                  for ( Index j=0; j<rh_abs.nelem(); ++j )
                    {
                      vmrs(0, j)  =  rh_abs[j] * sat_pres_abs[j];
                    }
                }
              else
                {
                  // Interpolating H2O VMR on p_abs grid
                  interpp (vmrs(0, Range(joker)), p_raw, vmr_raw, p_abs);
                }
	      
	      vmrs(1,Range(joker)) = 0.209;                     // O2 
	      vmrs(2,Range(joker)) = 0.782;                     // N2
	      
	      // Set the physical H2O profile from the H2O profile in vmrs:
	      h2o_abs.resize(p_abs.nelem());
	      h2o_abs = vmrs(0,Range(joker));
	      
	      // Set the physical N2 profile from the N2 profile in vmrs:
	      n2_abs.resize(p_abs.nelem());
	      n2_abs = vmrs(2,Range(joker));	  
	      
	      // Set t_ground from lowest level of launch:
	      t_ground = rd(0,1);
	      
	      // Set z_ground from lowest level of launch:
	      z_ground = rd(0,2);
	    }
	  
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

          // FIXME : ARTS-RTTOV paper
          if (calc_abs)
            {
              absbatch[i] = abs;
            }
      
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
                  za_pencil_profile,
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
		e_ground_profile,
		t_ground);
	  
          // Calculation of Jacobian
          if (calc_jac)
            {
              ArrayOfMatrix absloswfs;
              absloswfsCalc(// Output:
                            absloswfs,
                            // Input:
                            emission, 
                            los, 
                            source, 
                            trans, 
                            y, 
                            y_space, 
                            f_mono, 
                            e_ground_profile, 
                            t_ground);
              
              Vector k_grid = p_abs;
              
              TagGroups wfs_tgs( 1 );
              wfs_tgs = tgs[0];
                                          
              abs_per_tgReduce(// Output:
                               abs_per_tg,
                               // Input:
                               tgs, 
                               wfs_tgs);
              
              Matrix k;
              ArrayOfString k_names;
              Matrix k_aux;
              kSpecies(// Output:
                       k, 
                       k_names, 
                       k_aux,
                       // Input:
                       los, 
                       absloswfs, 
                       p_abs, 
                       t_abs, 
                       wfs_tgs, 
                       abs_per_tg, 
                       vmrs, 
                       k_grid,
                       // Keywords:
                       "frac");
              
              // Convert to RJ brightness temperatures
              for ( Index ik=0; ik<k.ncols(); ik++ )
                {
                  invrayjean( k(Range(joker),ik), f_mono, za_pencil_profile );
                }

              jacbatch[i] = k;
            }
          
          // Convert Radiance to Planck Tb
          yTB(y, f_mono, za_pencil_profile);
          
          // Assign y to this column of ybatch:
	  ybatch(Range(joker),i) = y;


	}
      // The launch did not reach upto or above 100 hPa
      else
	{
	  // Assign -1 to this column of ybatch:
	  ybatch(Range(joker),i) = -1.0;
	}
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
      const Matrix& rd = radiosonde_data[i];

      ConstVectorView p_raw = rd(Range(joker),0);
      ConstVectorView t_raw = rd(Range(joker),1);
      ConstVectorView z_raw = rd(Range(joker),2);

      Vector p_abs;
      
      VectorNLogSpace (p_abs,"p_abs",p_raw[0],p_raw[p_raw.nelem()-1],100);
      Vector t_abs( p_abs.nelem() );
      Vector z_abs( p_abs.nelem() );
                       
      interpp ( t_abs, p_raw, t_raw, p_abs );
      interpp ( z_abs, p_raw, z_raw, p_abs );        
      
      // Create vmrs:
      ConstVectorView h2o_raw = rd( Range(joker), 3 );
            
      Matrix vmrs( tgs.nelem(), p_abs.nelem() );
      interpp ( vmrs(0, Range( joker ) ), p_raw,h2o_raw, p_abs );
      vmrs( 1, Range( joker ) ) = 0.209;                  // O2 
      vmrs( 2, Range( joker ) ) = 0.782;                  // N2

      // Set the physical H2O profile from the H2O profile in vmrs:
      Vector h2o_abs = vmrs( 0, Range( joker ) );
      
      // Set the physical N2 profile from the N2 profile in vmrs:
      Vector n2_abs  = vmrs( 2, Range( joker ) );
      
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

