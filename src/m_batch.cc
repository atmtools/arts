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
#include "md.h"



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
              string&  filename,
        const string&  batchname,
        const string&  varname )
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
              MATRIX&   x, 
        const string&   batchname, 
        const string&   filename, 
        const string&   varname,
        const size_t&   length,
        const size_t&   n )
{
  string fname = filename;
  filename_batch( fname, batchname, varname );
  MatrixReadBinary( x, "", fname );
  if ( x.dim(1) != length )
  {
    ostringstream os;
    os << "The file " << fname << " contains data of length " << x.dim(1) 
       << ", but a length of " << length << " is expected.\n";
      throw runtime_error(os.str());
  }
  if ( x.dim(2) < n )
  {
    ostringstream os;
    os << "The file " << fname << " contains data for only " << x.dim(2) 
       << " spectra when " << n << " spectra shall be calculated.\n";
      throw runtime_error(os.str());
  }    
}


// Part common for both versions of BatchdataGaussianTemperatureProfiles
void temperature_profiles(
              MATRIX&   t,
              string&   fname,
        const VECTOR&   t_abs,
        const MATRIX&   s,
        const string&   batchname,
        const int&      n )
{
  fname = "";
  filename_batch( fname, batchname, "t_abs" );
  out2 << "  Filling the file " << fname << "\n"
       << "  with random temperature profiles.\n";
  out3 << "    Generating random data (" <<t_abs.dim()<<"x"<<n<<").\n";
  rand_data_gaussian( t, n, t_abs, s );
  MatrixWriteBinary( t, "", fname );
}



////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////

void BatchdataGaussianTemperatureProfiles(
      // WS Input:
        const VECTOR&    p_abs,
        const VECTOR&    t_abs,
        const VECTOR&    z_abs,
        const VECTOR&    h2o_abs,
        const MATRIX&    s,
        const string&    batchname,
      // Control Parameters:
        const int&       n,
        const Numeric&   g0,
        const Numeric&   pref,
        const Numeric&   zref,
        const int&       niter)
{
  string fname;
  MATRIX ts;
  temperature_profiles( ts, fname, t_abs, s, batchname, n );
  MATRIX zs(t_abs.dim(),n);
  VECTOR z, t;
  fname = "";
  filename_batch( fname, batchname, "z_abs" );
  out2 << "  Filling the file " << fname << "\n"
       << "  with vertical grids fulfilling hydrostatic eq.\n";
  for ( size_t i=1; i<=size_t(n); i++ )
  {
    col( t, i, ts );
    z = z_abs;
    z_absHydrostatic( z, p_abs, t, h2o_abs, g0, pref, zref, niter );
    put_in_col( zs, i, z );
  }
  MatrixWriteBinary( zs, "", fname );
}



void BatchdataGaussianTemperatureProfilesNoHydro(
        const VECTOR&   t_abs,
        const MATRIX&   s,
        const string&   batchname,
        const int&      n )
{
  string fname;
  MATRIX t;
  temperature_profiles( t, fname, t_abs, s, batchname, n );
}



void BatchdataGaussianSpeciesProfiles(
      // WS Input:
        const TagGroups&       tags,
        const ARRAYofVECTOR&   vmrs,
        const VECTOR&          p_abs,
        const VECTOR&          t_abs,
        const MATRIX&          s,
        const string&          batchname,
      // Control Parameters:
        const int&             n,
        const ARRAYofstring&   do_tags,
        const int&             unit )
{
  const size_t   ntags = do_tags.dim();    // Number of tags to do here 
  ARRAYofsizet   tagindex;                 // Index in tags for do_tags 
        string   fname;
        MATRIX   x;
   
  // Check if do_tags can be found in tag_groups and store indeces
  if ( ntags == 0 )
    throw runtime_error("No tags specified, no use in calling this function.");
  else
    get_tagindex_for_strings( tagindex, tags, do_tags );
  
  // Loop the tags
  for ( size_t itag=1; itag<=ntags; itag++ )
  {
    // Determine the name of the molecule for itag
    // The species lookup data:
    extern const ARRAY<SpeciesRecord> species_data;
    string molname = species_data[tags(tagindex(itag))(1).Species()].Name();

    // Create filename
    fname = "";
    filename_batch( fname, batchname, molname );
    out2 << "  Filling the file " << fname << "\n"
         << "  with random species profiles.\n";
    out3 << "    Generating random data (" <<p_abs.dim()<<"x"<<n<<").\n";

    // Handle the different units
    if ( unit == 1 )               // Relative
    {
      size_t   np = vmrs(tagindex(itag)).dim();
      size_t   row, col;
      Numeric  a;
      rand_data_gaussian( x, n, VECTOR(np,1.0), s );
      for ( row=1; row<=np; row++ )
      {
        a = vmrs(tagindex(itag))(row);
        for ( col=1; col<=size_t(n); col++ )
          x(row,col) *= a;
      }
    }

    else if ( unit == 2 )          // VMR
      rand_data_gaussian( x, n, vmrs(tagindex(itag)), s );

    else if ( unit == 3 )          // Number density
    {
      size_t   np = vmrs(tagindex(itag)).dim();
      size_t   row, col;
      Numeric  a;
      rand_data_gaussian( x, n, VECTOR(np,0.0), s );
      for ( row=1; row<=np; row++ )
      {
        a = number_density ( p_abs(row), t_abs(row) );
        for ( col=1; col<=size_t(n); col++ )
          x(row,col) = vmrs(tagindex(itag))(row) + x(row,col)/a;
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


void ybatchAbsAndRte(
      // WS Output:
              MATRIX&                     ybatch,
      // WS Input:
        const VECTOR&                     f_mono,
        const VECTOR&                     p_abs, 
        const VECTOR&                     t_abs,
        const ARRAYofVECTOR&              vmrs,
        const ARRAYofARRAYofLineRecord&   lines_per_tag,
        const ARRAYofsizet&               lineshape,
        const ARRAYofsizet&               lineshape_norm,
        const VECTOR&                     z_abs,
        const Numeric&                    z_plat,
        const VECTOR&                     za_pencil,
        const Numeric&                    l_step,
        const int&                        refr,
        const Numeric&                    l_step_refr,
        const VECTOR&                     refr_index,
        const Numeric&                    z_ground,
        const Numeric&                    r_geoid,
        const VECTOR&                     y_space,
        const VECTOR&                     e_ground,
        const Numeric&                    t_ground,
        const string&                     batchname,
        const TagGroups&                  tags,
      // Control Parameters:
        const int&                        ncalc,
        const int&                        do_t,
        const string&                     t_file,
        const int&                        do_z,
        const string&                     z_file,
        const ARRAYofstring&              do_tags,
        const ARRAYofstring&              tag_files,
        const int&                        do_f,
        const string&                     f_file,
        const int&                        do_za,
        const string&                     za_file )
{
  const size_t   np = p_abs.dim();         // Number of pressure levels
  const size_t   ntags = do_tags.dim();    // Number of tags to do here 
  ARRAYofsizet   tagindex;                 // Index in tags for do_tags 
   
  // Check if do_tags can be found in tag_groups and store indeces
  if ( ntags > 0 )
    get_tagindex_for_strings( tagindex, tags, do_tags );

  out2 << "  Reading data from files.\n";

  // Read data from file(s) or use the workspace varaible
  //
  // Temperature
  VECTOR t = t_abs;
  MATRIX Ts;
  if ( do_t )
    read_batchdata( Ts, batchname, t_file, "t_abs", np, ncalc );
  //
  // Altitudes
  VECTOR z = z_abs;
  MATRIX Zs;
  if ( do_z )
    read_batchdata( Zs, batchname, z_file, "z_abs", np, ncalc );
  //
  // Frequencies
  VECTOR f = f_mono;
  MATRIX Fs;
  if ( do_f )
    read_batchdata( Fs, batchname, f_file, "f_mono", f_mono.dim(), ncalc );
  //
  // Zenith angles
  VECTOR za = za_pencil;
  MATRIX ZAs;
  if ( do_za )
    read_batchdata( ZAs, batchname, za_file, "za_pencil", za_pencil.dim(), 
                                                                      ncalc );
  //
  // Species profiles
         size_t itag;
  ARRAYofVECTOR vs = vmrs;
  ARRAYofMATRIX VMRs(ntags);
  for ( itag=1; itag<=ntags; itag++ )
  {
    // The species lookup data:
    extern const ARRAY<SpeciesRecord> species_data;
    // Determine the name of the molecule for itag
    string molname = species_data[tags(tagindex(itag))(1).Species()].Name();
    read_batchdata( VMRs(itag), batchname, tag_files(itag), molname, np, 
                                                                      ncalc );
  }


  //--- Loop and calculate the spectra --------------------------------------
         MATRIX   abs;
  ARRAYofMATRIX   abs_per_tag;
            LOS   los;
  ARRAYofMATRIX   source;
  ARRAYofMATRIX   trans;
         VECTOR   y;

  out2 << "  Calculating spectra.\n";

  for ( size_t i=1; i<=size_t(ncalc); i++ )
  {
    out2 << "  -------- Batch spectrum " << i << " --------\n";
    // Copy from file data for the "do" quantities
    if ( do_t )
      col( t, i, Ts );
    if ( do_z )
      col( z, i, Zs );
    if ( do_f )
      col( f, i, Fs );
    if ( do_za )
      col( za, i, ZAs );
    for ( itag=1; itag<=ntags; itag++ )
      col( vs(tagindex(itag)), i, VMRs(itag) );

    // Do the calculations
    if ( (i==1) || do_t || ntags || do_f )
      absCalc( abs, abs_per_tag, f, p_abs, t, vs, lines_per_tag, lineshape, 
                                                              lineshape_norm);
    if ( (i==1) || do_z || do_za )   
      losCalc( los, z_plat, za, l_step, p_abs, z, refr, l_step_refr, 
                                              refr_index, z_ground, r_geoid );
    if (  (i==1) || do_t || do_f )   
      sourceCalc( source, los, p_abs, t, f );
    transCalc( trans, los, p_abs, abs );
    yRte ( y, los, f, y_space, source, trans, e_ground, t_ground );
    
    // Move to ybatch
    if ( i == 1 )
      ybatch.newsize( y.dim(), ncalc );
    put_in_col( ybatch, i, y );
  }
  out2 << "  ------------------------------------\n";
}
