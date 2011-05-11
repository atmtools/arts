/* Copyright (C) 2002-2008
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_atmosphere.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-16

  \brief  Workspace functions to set variables defining the atmosphere 
          (excluding the surface).

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <algorithm>
#include "agenda_class.h"
#include "arts.h"
#include "check_input.h"
#include "matpackIII.h"
#include "messages.h"
#include "special_interp.h"
#include "absorption.h"
#include "abs_species_tags.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "interpolation_poly.h"
#include "xml_io.h"

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_P_GRID;
extern const Index GFIELD4_LAT_GRID;
extern const Index GFIELD4_LON_GRID;

extern const Numeric DEG2RAD;

/*===========================================================================
 *=== Helper functions
 *===========================================================================*/

//! atm_fields_compactExpand
/*!
   Add a species to an *atm_fields_compact*. Does not add any content, but only
   resizes the data and adds a field to the *ArrayOfString* respresenting the
   species for this *GriddedField4*. This helper function is used by e.g.
   *atm_fields_compactAddSpecies* and *atm_fields_compactAddConstant*.
  

   \retval  af         The new atm_fields_compact 
   \retval  nf         The new number of fields
   \param   name       Name of new field

   \author Gerrit Holl
   \date   2011-05-04
*/

void atm_fields_compactExpand(GriddedField4& af,
                              Index& nf,
                              const String& name)
{
  // Number of fields already present:
  nf = af.get_string_grid(GFIELD4_FIELD_NAMES).nelem();

// Most of the functionality moved from atm_fields_compactAddConstant when
// atm_fields_compactAddSpecies was added, in order to share the code.
//
// Commented out by Gerrit 2011-05-04. I believe this check is not needed.
// Oliver agrees. We can still infer the dimensions even if there are zero
// fields; e.g., the data might have dimensions (0, 4, 3, 8).
// 
//  if (0==nf)
//    {
//      ostringstream os;
//      os << "The *atm_fields_compact* must already contain at least one field,\n"
//         << "so that we can infer the dimensions from that.";
//      throw runtime_error( os.str() );
//    }

  // Add name of new field to field name list:
  af.get_string_grid(GFIELD4_FIELD_NAMES).push_back(name);

  // Save original fields:
  const Tensor4 dummy = af.data;

  // Adjust size:
  af.resize( nf+1, dummy.npages(), dummy.nrows(), dummy.ncols() );
  nf++; // set to new number of fields

  // Copy back original field:
  af.data( Range(0,nf-1), Range(joker), Range(joker), Range(joker) ) = dummy;
}



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void atm_checkedCalc(
         Index&          atm_checked,
   const Index&          atmosphere_dim,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Tensor3&        z_field,
   const Tensor3&        t_field,
   const Tensor4&        vmr_field,
   const Matrix&         r_geoid,
   const Matrix&         z_surface )
{
  atm_checked = 1;

  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  // Consistency between dim, grids and atmospheric fields/surfaces
  //
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "z_field", z_field, atmosphere_dim, 
                                                  p_grid, lat_grid, lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, 
                                                  p_grid, lat_grid, lon_grid );
  // Ignore vmr_field if abs_species is empty
  if( abs_species.nelem() )
    chk_atm_field( "vmr_field", vmr_field, atmosphere_dim, abs_species.nelem(),
                                                  p_grid, lat_grid, lon_grid );
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_surface", z_surface, atmosphere_dim, lat_grid, lon_grid );

  // Check that z_field has strictly increasing pages.
  for( Index row=0; row<z_field.nrows(); row++ )
    {
      for( Index col=0; col<z_field.ncols(); col++ )
        {
          ostringstream os;
          os << "z_field (for latitude nr " << row << " and longitude nr " 
             << col << ")";
          chk_if_increasing( os.str(), z_field(joker,row,col) ); 
        }
    }

  // Check that there is no gap between the surface and lowest pressure 
  // level
  for( Index row=0; row<z_surface.nrows(); row++ )
    {
      for( Index col=0; col<z_surface.ncols(); col++ )
        {
          if( z_surface(row,col)<z_field(0,row,col) ||
                   z_surface(row,col)>=z_field(z_field.npages()-1,row,col) )
            {
              ostringstream os;
              os << "The surface altitude (*z_surface*) cannot be outside "
                 << "of the altitudes in *z_field*.";
              if( atmosphere_dim > 1 )
                os << "\nThis was found to be the case for:\n"
                   << "latitude " << lat_grid[row];
              if( atmosphere_dim > 2 )
                os << "\nlongitude " << lon_grid[col];
              throw runtime_error( os.str() );
            }
        }
    }
}



// Workspace method, doxygen header will be auto-generated.
// 2007-07-25 Stefan Buehler
void atm_fields_compactFromMatrix(// WS Output:
                                  GriddedField4& af, // atm_fields_compact
                                  // WS Input:
                                  const Index& atmosphere_dim,
                                  // WS Generic Input:
                                  const Matrix& im,
                                  // Control Parameters:
                                  const ArrayOfString& field_names )
{
  if (1!=atmosphere_dim)
    {
      ostringstream os; 
      os << "Atmospheric dimension must be one.";
      throw runtime_error( os.str() );
    }

  const Index np = im.nrows();   // Number of pressure levels.
  const Index nf = im.ncols()-1; // Total number of fields. 
  Index nf_1; // Number of required fields. 
                  // All fields called "ignore" are ignored.
  String fn_upper; // Temporary variable to hold upper case field_names.

  if (field_names.nelem()!=nf)
    {
      ostringstream os; 
      os << "Cannot extract fields from Matrix.\n"
         << "*field_names* must have one element less than there are\n"
         << "matrix columns.";
      throw runtime_error( os.str() );
    }


  // Remove additional fields from the field_names. All fields that
  // are flagged by 'ignore' in the field names, small or large letters,
  // are removed.
  nf_1 = 0;
  for(Index f = 0; f < field_names.nelem(); f++)
    {
      fn_upper = field_names[f];
      std::transform ( fn_upper.begin(),  fn_upper.end(), fn_upper.begin(), ::toupper);
      if(fn_upper != "IGNORE") nf_1++;
    }

  // Copy required field_names to a new variable called field_names_1
  ArrayOfString field_names_1(nf_1);
  for (Index f=0; f< nf_1; f++) field_names_1[f] = field_names[f];


  //  out3 << "Copying *" << im_name << "* to *atm_fields_compact*.\n";
  
  af.set_grid(GFIELD4_FIELD_NAMES, field_names_1);

  af.set_grid(GFIELD4_P_GRID, im(Range(joker),0));
  
  af.set_grid(GFIELD4_LAT_GRID, Vector());
  af.set_grid(GFIELD4_LON_GRID, Vector());
  
  af.resize(nf_1,np,1,1); // Resize it according to the required fields
  af.data(Range(joker),Range(joker),0,0) = transpose(im(Range(joker),Range(1,nf_1)));
}

// Workspace method, doxygen header will be auto-generated.
// 2011-01-24 Daniel Kreyling
void atm_fields_compactFromMatrixChevalAll(// WS Output:
                                  GriddedField4& af_all, // atm_fields_compact_all
				  GriddedField4& af_vmr, // atm_fields_compact
				  // WS Input:
                                  const Index& atmosphere_dim,
                                  // WS Generic Input:
                                  const Matrix& im,
                                  // Control Parameters:
                                  const ArrayOfString& field_names

 				    )
{
  // NOTE: follwoing section is HARD WIRED 
  // This method can only be applied to matrix data sets with a specific order of columns
  // See method documentation!
  
  if (1!=atmosphere_dim)
    {
      ostringstream os; 
      os << "Atmospheric dimension must be one.";
      throw runtime_error( os.str() );
    }

  const Index np = im.nrows();   // Number of pressure levels.
  const Index nf = im.ncols()-1; // Number of colums without pressure.
 
  
  Index nf_1, nf_2; // Number of required fields. 
                    // All fields called "ignore" are ignored.
  String fn_upper; // Temporary variable to hold upper case field_names.

  if (field_names.nelem()!= nf)
    {
      ostringstream os; 
      os << "Cannot extract fields from Matrix.\n"
         << "*field_names* must have one element less than there are\n"
         << "matrix columns.";
      throw runtime_error( os.str() );
    }


  // Remove additional fields from the field_names. All fields that
  // are flagged by 'ignore' in the field names, small or large letters,
  // are removed. The remaining number of fields is stored in *nf_1*.
  // Usage: total batch_atm_fields_compact_all
  nf_1 = 0;
  nf_2 = 0;
  ArrayOfIndex intarr;
  for(Index f = 0; f < field_names.nelem(); f++)
    {
      fn_upper = field_names[f];
      std::transform ( fn_upper.begin(),  fn_upper.end(), fn_upper.begin(), ::toupper);
      if(fn_upper != "IGNORE" ) nf_1++;
      
      // Remove all 'ignore' and massdensity fields in this step, so that only T, z and vmrs remain.
      // The Index of these fields in stroed in *intarr*, for later access.
      // The remaining number of fields is stored in *nf_2*.
      if ( fn_upper !="IGNORE" && field_names[f] != "LWC" && field_names[f] != "IWC" &&
	field_names[f] != "Rain" && field_names[f] != "Snow" ) 
      {
	nf_2++;
	intarr.push_back(f);
	
      }
    }
    
  //------- write batch_atm_fields_compact_all ----------------------------------------------------
  // including massdenity fields!
  
  // Copy required field_names to a new variable called field_names_1
  ArrayOfString field_names_1(nf_1); //f_names_2(2);
  for (Index f=0; f< nf_1; f++) field_names_1[f] = field_names[f];

  //  out3 << "Copying *" << im_name << "* to *atm_fields_compact*.\n";

  //cout<<nf<<"\t"<<nf_1<<"\t"<<field_names.nelem()<<endl;
   

  af_all.set_grid(GFIELD4_FIELD_NAMES, field_names_1);

  af_all.set_grid(GFIELD4_P_GRID, im(Range(joker),0));
  
  af_all.set_grid(GFIELD4_LAT_GRID, Vector());
  af_all.set_grid(GFIELD4_LON_GRID, Vector());
  
  af_all.resize(nf_1,np,1,1); // Resize it according to the required fields
  af_all.data(Range(joker),Range(joker),0,0) = transpose(im(Range(joker),Range(1,nf_1)));
  
 
  //------- write batch_atm_fields_compact -------------------------------------------------------------
  // excluding massdenity fields!
  
  ArrayOfString field_names_2(nf_2);
  for ( Index i=0; i<nf_2; i++ ) field_names_2[i] = field_names[intarr[i]] ;
  
  af_vmr.set_grid(GFIELD4_FIELD_NAMES, field_names_2);
  
  af_vmr.set_grid(GFIELD4_P_GRID, im(Range(joker),0));
  
  af_vmr.set_grid(GFIELD4_LAT_GRID, Vector());
  af_vmr.set_grid ( GFIELD4_LON_GRID, Vector() );

  af_vmr.resize ( nf_2,np,1,1 ); // Resize it according to the required fields
  for ( Index i=0; i<nf_2; i++ )
  {
    // write T, z and VMRs
    af_vmr.data ( Range(i,1) ,Range ( joker ),0,0 ) = transpose ( im ( Range ( joker ), Range(intarr[i]+1,1)) );
  }
}


// Workspace method, doxygen header is auto-generated.
// 2007-07-31 Stefan Buehler
// 2011-05-04 Adapted by Gerrit Holl
void atm_fields_compactAddConstant(// WS Output:
                                   GriddedField4& af,
                                   // Control Parameters:
                                   const String& name,
                                   const Numeric& value)
{
  Index nf; // Will hold new size

  // Add book
  atm_fields_compactExpand(af, nf, name);
  
  // Add the constant value:
  af.data( nf-1, Range(joker), Range(joker), Range(joker) ) = value;
}

// Workspace method, doxygen header is auto-generated
// 2011-05-02 Gerrit Holl
void atm_fields_compactAddSpecies(// WS Output:
                                GriddedField4& atm_fields_compact,
                                // WS Generic Input:
                                const String& name,
                                const GriddedField3& species)
{

  assert(atm_fields_compact.checksize());
  assert(species.checksize());

  ConstVectorView af_p_grid = atm_fields_compact.get_numeric_grid(GFIELD4_P_GRID);
  ConstVectorView af_lat_grid = atm_fields_compact.get_numeric_grid(GFIELD4_LAT_GRID);
  ConstVectorView af_lon_grid = atm_fields_compact.get_numeric_grid(GFIELD4_LON_GRID);
  ConstVectorView sp_p_grid = species.get_numeric_grid(GFIELD3_P_GRID);
  ConstVectorView sp_lat_grid = species.get_numeric_grid(GFIELD3_LAT_GRID);
  ConstVectorView sp_lon_grid = species.get_numeric_grid(GFIELD3_LON_GRID);

  Index new_n_fields; // To be set in next line
  atm_fields_compactExpand(atm_fields_compact, new_n_fields, name);


  // Interpolate species to atm_fields_compact

  // Common for all dim
  chk_interpolation_grids("species p_grid to atm_fields_compact p_grid",
                          sp_p_grid,
                          af_p_grid);
  ArrayOfGridPos p_gridpos(af_p_grid.nelem());
  // gridpos(p_gridpos, sp_p_grid, af_p_grid);
  p2gridpos(p_gridpos, sp_p_grid, af_p_grid);

  if (sp_lat_grid.nelem() > 1)
  {
      // Common for all dim>=2
      chk_interpolation_grids("species lat_grid to atm_fields_compact lat_grid",
                              sp_lat_grid,
                              af_lat_grid);
      ArrayOfGridPos lat_gridpos(af_lat_grid.nelem());
      gridpos(lat_gridpos, sp_lat_grid, af_lat_grid);

      if (sp_lon_grid.nelem() > 1)
      { // 3D-case
          chk_interpolation_grids("species lon_grid to atm_fields_compact lon_grid",
                                  sp_lon_grid,
                                  af_lon_grid);
          ArrayOfGridPos lon_gridpos(af_lon_grid.nelem());
          gridpos(lon_gridpos, sp_lon_grid, af_lon_grid);

          Tensor4 itw(p_gridpos.nelem(), lat_gridpos.nelem(), lon_gridpos.nelem(), 8);
          interpweights(itw, p_gridpos, lat_gridpos, lon_gridpos);

          Tensor3 newfield(af_p_grid.nelem(), af_lat_grid.nelem(), af_lon_grid.nelem());
          interp(newfield, itw, species.data, p_gridpos, lat_gridpos, lon_gridpos);

          atm_fields_compact.data(new_n_fields-1, joker, joker, joker) = newfield;
      } else { // 2D-case
          
          Tensor3 itw(p_gridpos.nelem(), lat_gridpos.nelem(), 4);
          interpweights(itw, p_gridpos, lat_gridpos);

          Matrix newfield(af_p_grid.nelem(), af_lat_grid.nelem());
          interp(newfield, itw, species.data(joker, joker, 0), p_gridpos, lat_gridpos);

          atm_fields_compact.data(new_n_fields-1, joker, joker, 0) = newfield;
      }
  } else { // 1D-case
      Matrix itw(p_gridpos.nelem(), 2);
      interpweights(itw, p_gridpos);

      Vector newfield(af_p_grid.nelem());
      interp(newfield, itw, species.data(joker, 0, 0), p_gridpos);

      atm_fields_compact.data(new_n_fields-1, joker, 0, 0) = newfield;
  }

}

// Workspace method, doxygen header is auto-generated
// 2011-05-09 Gerrit Holl
void batch_atm_fields_compactAddSpecies(// WS Output:
                                        ArrayOfGriddedField4& batch_atm_fields_compact,
                                        // WS Generic Input:
                                        const String& name,
                                        const GriddedField3& species)
{
    const Index nelem = batch_atm_fields_compact.nelem();

    // FIXME: can this loop be parallelised?
    for (Index i=0; i<nelem; i++)
    {
        atm_fields_compactAddSpecies(batch_atm_fields_compact[i], name, species);
    }
}

// Workspace method, doxygen header is auto-generated.
void batch_atm_fields_compactFromArrayOfMatrix(// WS Output:
                                               ArrayOfGriddedField4& batch_atm_fields_compact,
                                               // WS Input:
                                               const Index& atmosphere_dim,
                                               // WS Generic Input:
                                               const ArrayOfMatrix& am,
                                               // Control Parameters:
                                               const ArrayOfString& field_names,
                                               const ArrayOfString& extra_field_names,
                                               const Vector& extra_field_values)
{
  const Index amnelem = am.nelem();

  // We use the existing WSMs atm_fields_compactFromMatrix and
  // atm_fields_compactAddConstant to do most of the work.

  // Check that extra_field_names and extra_field_values have matching
  // dimensions:
  if (extra_field_names.nelem() != extra_field_values.nelem())
    {
      ostringstream os; 
      os << "The keyword arguments extra_field_names and\n"
         << "extra_field_values must have matching dimensions.";
      throw runtime_error( os.str() );
    }

  // Make output variable the proper size:
  batch_atm_fields_compact.resize(amnelem);

  // Loop the batch cases:
/*#pragma omp parallel for                                     \
  if(!arts_omp_in_parallel())                                \
  default(none)                                              \
  shared(am, atmosphere_dim, batch_atm_fields_compact,       \
         extra_field_names, extra_field_values, field_names) */
#pragma omp parallel for                                     \
  if(!arts_omp_in_parallel())
  for (Index i=0; i<amnelem; ++i)
    {
      // All the input variables are visible here, despite the
      // "default(none)". The reason is that they are return by
      // reference arguments of this function, which are shared
      // automatically. 

      // The try block here is necessary to correctly handle
      // exceptions inside the parallel region. 
      try
        {
          atm_fields_compactFromMatrix(batch_atm_fields_compact[i],
                                       atmosphere_dim,
                                       am[i],
                                       field_names);

          for (Index j=0; j<extra_field_names.nelem(); ++j)
            atm_fields_compactAddConstant(batch_atm_fields_compact[i],
                                          extra_field_names[j],
                                          extra_field_values[j]);
        }
      catch (runtime_error e)
        {
          exit_or_rethrow(e.what());
        }
    }    
}


// Workspace method, doxygen header is auto-generated.
// 2011-01-24 Daniel Kreyling
void batch_atm_fields_compactFromArrayOfMatrixChevalAll(// WS Output:
                                               ArrayOfGriddedField4& batch_atm_fields_compact,
					       ArrayOfGriddedField4& batch_atm_fields_compact_all,
                                               // WS Input:
                                               const Index& atmosphere_dim,
                                               // WS Generic Input:
                                               const ArrayOfMatrix& am,
                                               // Control Parameters:
                                               const ArrayOfString& field_names,
					       const ArrayOfString& extra_field_names,
                                               const Vector& extra_field_values)
{
  const Index amnelem = am.nelem();

  // We use the existing WSMs atm_fields_compactFromMatrix and
  // atm_fields_compactAddConstant to do most of the work.

  // Check that extra_field_names and extra_field_values have matching
  // dimensions:
  if (extra_field_names.nelem() != extra_field_values.nelem())
    {
      ostringstream os; 
      os << "The keyword arguments extra_field_names and\n"
         << "extra_field_values must have matching dimensions.";
      throw runtime_error( os.str() );
    }

  // Make output variable the proper size:
  batch_atm_fields_compact.resize(amnelem);
  batch_atm_fields_compact_all.resize(amnelem);
  // Loop the batch cases:
/*#pragma omp parallel for                                     \
  if(!arts_omp_in_parallel())                                \
  default(none)                                              \
  shared(am, atmosphere_dim, batch_atm_fields_compact,       \
         extra_field_names, extra_field_values, field_names) */
#pragma omp parallel for                                     \
  if(!arts_omp_in_parallel())
  for (Index i=0; i<amnelem; ++i)
    {
      // All the input variables are visible here, despite the
      // "default(none)". The reason is that they are return by
      // reference arguments of this function, which are shared
      // automatically. 

      // The try block here is necessary to correctly handle
      // exceptions inside the parallel region. 
      try
        {
          atm_fields_compactFromMatrixChevalAll( batch_atm_fields_compact_all[i],
					      batch_atm_fields_compact[i],
					      atmosphere_dim,
					      am[i],
					      field_names);


          for (Index j=0; j<extra_field_names.nelem(); ++j){
            atm_fields_compactAddConstant( batch_atm_fields_compact_all[i],
                                          extra_field_names[j],
                                          extra_field_values[j]);
	    
	    atm_fields_compactAddConstant(batch_atm_fields_compact[i],
                                          extra_field_names[j],
                                          extra_field_values[j]);}
        }
      catch (runtime_error e)
        {
          exit_or_rethrow(e.what());
        }
    }    
}

// Workspace method, doxygen header will be auto-generated.
// 2010-11-29 Daniel Kreyling
void AtmFieldsFromCompactChevalAll(// WS Output:
                          Vector& p_grid,
                          Vector& lat_grid,
                          Vector& lon_grid,
                          Tensor3& t_field,
                          Tensor3& z_field,
			  Tensor4& massdensity_field,
                          Tensor4& vmr_field,
		          // WS Input:
                          const ArrayOfArrayOfSpeciesTag& abs_species,
			  // const ArrayOfArrayOfSpeciesTag& 
                          const GriddedField4& atm_fields_compact_all,
                          const Index&  atmosphere_dim )
{
  // Make a handle on atm_fields_compact to save typing:
  const GriddedField4& c = atm_fields_compact_all;
  
  // Check if the grids in our data match atmosphere_dim
  // (throws an error if the dimensionality is not correct):
  chk_atm_grids(atmosphere_dim,
                c.get_numeric_grid(GFIELD4_P_GRID),
                c.get_numeric_grid(GFIELD4_LAT_GRID),
                c.get_numeric_grid(GFIELD4_LON_GRID));

  const Index nf   = c.get_grid_size(GFIELD4_FIELD_NAMES);
  const Index np   = c.get_grid_size(GFIELD4_P_GRID);
  const Index nlat = c.get_grid_size(GFIELD4_LAT_GRID);
  const Index nlon = c.get_grid_size(GFIELD4_LON_GRID);

  // Grids:
  p_grid = c.get_numeric_grid(GFIELD4_P_GRID);
  lat_grid = c.get_numeric_grid(GFIELD4_LAT_GRID);
  lon_grid = c.get_numeric_grid(GFIELD4_LON_GRID);

  // NOTE: follwoing section is HARD WIRED 
  // The order of the fields is:
  // T[K] z[m] LWC[kg/m^3] IWC[kg/m^3] Rain[kg/(m2*s)] Snow[kg/(m2*s)] VMR_1[1] ... VMR_n[1]

  // Number of VMR species:
  const Index ns = nf-6;
  
  // Check that there is at least one VMR species:
  if (ns<1)
    {
      ostringstream os;
      os << "There must be at least three fields in *atm_fields_compact*.\n"
         << "T, z, and at least one VMR.";
      throw runtime_error( os.str() );
    }

  // Check that first field is T:
  if (c.get_string_grid(GFIELD4_FIELD_NAMES)[0] != "T")
    {
      ostringstream os;
      os << "The first field must be \"T\", but it is:"
         << c.get_string_grid(GFIELD4_FIELD_NAMES)[0];
      throw runtime_error( os.str() );
    }

  // Check that second field is z:
  if (c.get_string_grid(GFIELD4_FIELD_NAMES)[1] != "z")
    {
      ostringstream os;
      os << "The second field must be \"z\"*, but it is:"
         << c.get_string_grid(GFIELD4_FIELD_NAMES)[1];
      throw runtime_error( os.str() );
    }

  // Check that third field is LWC:
  if (c.get_string_grid(GFIELD4_FIELD_NAMES)[2] != "LWC")
    {
      ostringstream os;
      os << "The third field must be \"LWC\"*, but it is:"
         << c.get_string_grid(GFIELD4_FIELD_NAMES)[2];
      throw runtime_error( os.str() );
    }

  // Check that fourth field is IWC:
  if (c.get_string_grid(GFIELD4_FIELD_NAMES)[3] != "IWC")
    {
      ostringstream os;
      os << "The fourth field must be \"IWC\"*, but it is:"
         << c.get_string_grid(GFIELD4_FIELD_NAMES)[3];
      throw runtime_error( os.str() );
    }

  // Check that fifth field is Rain:
  if (c.get_string_grid(GFIELD4_FIELD_NAMES)[4] != "Rain")
    {
      ostringstream os;
      os << "The fifth field must be \"Rain\"*, but it is:"
         << c.get_string_grid(GFIELD4_FIELD_NAMES)[4];
      throw runtime_error( os.str() );
    }

  // Check that sixth field is Snow:
  if (c.get_string_grid(GFIELD4_FIELD_NAMES)[5] != "Snow")
    {
      ostringstream os;
      os << "The sixth field must be \"IWC\"*, but it is:"
         << c.get_string_grid(GFIELD4_FIELD_NAMES)[5];
      throw runtime_error( os.str() );
    }
    
  // Check that the other fields are VMR fields and match abs_species:
  for (Index i=0; i<ns; ++i)
    {
      const String tf_species = c.get_string_grid(GFIELD4_FIELD_NAMES)[6+i];
      
      // Get name of species from abs_species:      
      extern const Array<SpeciesRecord> species_data;  // The species lookup data:
      const String as_species = species_data[abs_species[i][0].Species()].Name();

      // Species in field name and abs_species should be the same:
      if (tf_species != as_species)
        {
          ostringstream os;
          os << "Field name not valid: "
             << tf_species << "\n"
             << "Based on *abs_species*, the field name should be: "
             << as_species;
          throw runtime_error( os.str() );
        }
    }

  // Temperature field (first field):
  t_field.resize(np,nlat,nlon);
  t_field = c.data(0,Range(joker),Range(joker),Range(joker));

  // Altitude profile (second field):
  z_field.resize(np,nlat,nlon);
  z_field = c.data(1,Range(joker),Range(joker),Range(joker));

  //write all massdensity profile to one Tensor4
  massdensity_field.resize(4,np,nlat,nlon);
  massdensity_field = c.data(Range(2,4),Range(joker),Range(joker),Range(joker));

  // VMR profiles (remaining fields):
  vmr_field.resize(ns,np,nlat,nlon);
  vmr_field = c.data(Range(6,ns),Range(joker),Range(joker),Range(joker));
    
}


// Workspace method, doxygen header will be auto-generated.
// 2007-07-24 Stefan Buehler
void AtmFieldsFromCompact(// WS Output:
                          Vector& p_grid,
                          Vector& lat_grid,
                          Vector& lon_grid,
                          Tensor3& t_field,
                          Tensor3& z_field,
                          Tensor4& vmr_field,
                          // WS Input:
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const GriddedField4& atm_fields_compact,
                          const Index&  atmosphere_dim )
{
  // Make a handle on atm_fields_compact to save typing:
  const GriddedField4& c = atm_fields_compact;
  
  // Check if the grids in our data match atmosphere_dim
  // (throws an error if the dimensionality is not correct):
  chk_atm_grids(atmosphere_dim,
                c.get_numeric_grid(GFIELD4_P_GRID),
                c.get_numeric_grid(GFIELD4_LAT_GRID),
                c.get_numeric_grid(GFIELD4_LON_GRID));

  const Index nf   = c.get_grid_size(GFIELD4_FIELD_NAMES);
  const Index np   = c.get_grid_size(GFIELD4_P_GRID);
  const Index nlat = c.get_grid_size(GFIELD4_LAT_GRID);
  const Index nlon = c.get_grid_size(GFIELD4_LON_GRID);

  // Grids:
  p_grid = c.get_numeric_grid(GFIELD4_P_GRID);
  lat_grid = c.get_numeric_grid(GFIELD4_LAT_GRID);
  lon_grid = c.get_numeric_grid(GFIELD4_LON_GRID);

  // The order of the fields is:
  // T[K] z[m] VMR_1[1] ... VMR_n[1]

  // Number of VMR species:
  const Index ns = nf-2;
  
  // Check that there is at least one VMR species:
  if (ns<1)
    {
      ostringstream os;
      os << "There must be at least three fields in *atm_fields_compact*.\n"
         << "T, z, and at least one VMR.";
      throw runtime_error( os.str() );
    }

  // Check that first field is T:
  if (c.get_string_grid(GFIELD4_FIELD_NAMES)[0] != "T")
    {
      ostringstream os;
      os << "The first field must be \"T\", but it is:"
         << c.get_string_grid(GFIELD4_FIELD_NAMES)[0];
      throw runtime_error( os.str() );
    }

  // Check that second field is z:
  if (c.get_string_grid(GFIELD4_FIELD_NAMES)[1] != "z")
    {
      ostringstream os;
      os << "The second field must be \"z\"*, but it is:"
         << c.get_string_grid(GFIELD4_FIELD_NAMES)[1];
      throw runtime_error( os.str() );
    }

  // Check that the other fields are VMR fields and match abs_species:
  for (Index i=0; i<ns; ++i)
    {
      const String tf_species = c.get_string_grid(GFIELD4_FIELD_NAMES)[2+i];
      
      // Get name of species from abs_species:      
      extern const Array<SpeciesRecord> species_data;  // The species lookup data:
      const String as_species = species_data[abs_species[i][0].Species()].Name();

      // Species in field name and abs_species should be the same:
      if (tf_species != as_species)
        {
          ostringstream os;
          os << "Field name not valid: "
             << tf_species << "\n"
             << "Based on *abs_species*, the field name should be: "
             << as_species;
          throw runtime_error( os.str() );
        }
    }

  // Temperature field (first field):
  t_field.resize(np,nlat,nlon);
  t_field = c.data(0,Range(joker),Range(joker),Range(joker));

  // Altitude profile (second field):
  z_field.resize(np,nlat,nlon);
  z_field = c.data(1,Range(joker),Range(joker),Range(joker));

  // VMR profiles (remaining fields):
  vmr_field.resize(ns,np,nlat,nlon);
  vmr_field = c.data(Range(2,ns),Range(joker),Range(joker),Range(joker));
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmosphereSet1D(
        // WS Output:
              Index&    atmosphere_dim,
              Vector&   lat_grid,
              Vector&   lon_grid )
{
  out2 << "  Sets the atmospheric dimensionality to 1.\n";
  out3 << "    atmosphere_dim = 1\n";
  out3 << "    lat_grid is set to be an empty vector\n";
  out3 << "    lon_grid is set to be an empty vector\n";
  atmosphere_dim = 1;
  lat_grid.resize(0);
  lon_grid.resize(0);
 
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmosphereSet2D(
        // WS Output:
              Index&    atmosphere_dim,
              Vector&   lon_grid )
{
  out2 << "  Sets the atmospheric dimensionality to 2.\n";
  out3 << "    atmosphere_dim = 2\n";
  out3 << "    lon_grid is set to be an empty vector\n";
  atmosphere_dim = 2;
  lon_grid.resize(0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmosphereSet3D(
        // WS Output:
              Index&    atmosphere_dim )
{
  out2 << "  Sets the atmospheric dimensionality to 3.\n";
  out3 << "    atmosphere_dim = 3\n";
  atmosphere_dim = 3;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsCalc(//WS Output:
                   Tensor3& t_field,
                   Tensor3& z_field,
                   Tensor4& vmr_field,
                   //WS Input
                   const Vector&         p_grid,
                   const Vector&         lat_grid,
                   const Vector&         lon_grid,
                   const GriddedField3&        t_field_raw,
                   const GriddedField3&        z_field_raw,
                   const ArrayOfGriddedField3& vmr_field_raw,
                   const Index&          atmosphere_dim,
                   // WS Generic Input:
                   const Index& interp_order
                   )
{
  const ConstVectorView tfr_p_grid   = t_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView tfr_lat_grid = t_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView tfr_lon_grid = t_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
  const ConstVectorView zfr_p_grid   = z_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView zfr_lat_grid = z_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView zfr_lon_grid = z_field_raw.get_numeric_grid(GFIELD3_LON_GRID);

  out2 << "  Interpolation order: " << interp_order << "\n";
  
  // Basic checks of input variables
  //
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  

  // Note that we are using the special function p2gridpos_poly below for
  // all pressure interpolations. This does the usual ARTS pressure
  // interpolation in log(p). We don't have to take logs here
  // explicitly, since it is done by p2gridpos.


  //==========================================================================
  if ( atmosphere_dim == 1)
    {
      if( !( tfr_lat_grid.nelem() == 1 &&
             tfr_lon_grid.nelem() == 1 ))
        throw runtime_error(
                            "Temperature data (T_field) has wrong dimension "
                            "(2D or 3D).\n"
                            );

      if( !( zfr_lat_grid.nelem() == 1 &&
             zfr_lon_grid.nelem() == 1 ))
        throw runtime_error(
                            "Altitude data (z_field) has wrong dimension "
                            "(2D or 3D).\n"
                            );

      //Resize variables
      t_field.resize(p_grid.nelem(), 1, 1);
      z_field.resize(p_grid.nelem(), 1, 1);
      vmr_field.resize(vmr_field_raw.nelem(), p_grid.nelem(), 1, 1);

      // Gridpositions:
      ArrayOfGridPosPoly gp_p(p_grid.nelem());
  
      // Interpolate t_field:
      
      // Check that interpolation grids are ok (and throw a detailed
      // error message if not):
      chk_interpolation_grids("Raw temperature to p_grid, 1D case",
                              tfr_p_grid,
                              p_grid,
                              interp_order);
 
      // Calculate grid positions:
      p2gridpos_poly( gp_p, tfr_p_grid, p_grid, interp_order );

      // Interpolation weights:
      Matrix itw(p_grid.nelem(), interp_order+1);
      // (2 interpolation weights are required for 1D interpolation)
      interpweights( itw, gp_p);
  
      // Interpolate:
      interp( t_field(joker, 0, 0), itw, 
              t_field_raw.data(joker, 0, 0),  gp_p);

  
      // Interpolate z_field:
      
      // Check that interpolation grids are ok (and throw a detailed
      // error message if not):
      chk_interpolation_grids("Raw z to p_grid, 1D case",
                              zfr_p_grid,
                              p_grid,
                              interp_order);

      // Calculate grid positions:
      p2gridpos_poly( gp_p, zfr_p_grid, p_grid, interp_order );
     
      // Interpolation weights:
      interpweights( itw, gp_p );
      
      // Interpolate:
      interp( z_field(joker, 0, 0), itw,
              z_field_raw.data(joker, 0, 0), gp_p);
      
  
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          ostringstream os; 

          if( !( vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID).nelem() == 1 &&
                 vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LON_GRID).nelem() == 1 ))
            {
              os << "VMR data of the " << gas_i << "th species has "
                 << "wrong dimension (2D or 3D). \n";
              throw runtime_error( os.str() );
            }
          
          // Check that interpolation grids are ok (and throw a detailed
          // error message if not):
          os << "Raw VMR[" << gas_i << "] to p_grid, 1D case";
          chk_interpolation_grids(os.str(),
                                  vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID),
                                  p_grid,
                                  interp_order);

          // Calculate grid positions:
          p2gridpos_poly(gp_p, 
                         vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID), 
                         p_grid, 
                         interp_order);
          
          // Interpolation weights:
          interpweights( itw, gp_p);
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, 0, 0),
                  itw, vmr_field_raw[gas_i].data(joker, 0, 0), gp_p);
        }
      
    }

  //=========================================================================
  else if(atmosphere_dim == 2)
    {
      if( tfr_lat_grid.nelem() == 1 &&
          tfr_lon_grid.nelem() == 1 )
        throw runtime_error(
                            "Raw data has wrong dimension (1D). "
                            "You have to use \n"
                            "AtmFieldsCalcExpand1D instead of AtmFieldsCalc."
                            );
      
      //Resize variables
      t_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
      z_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
      vmr_field.resize(vmr_field_raw.nelem(), p_grid.nelem(), lat_grid.nelem(),
                       1);
      
      
      // Gridpositions:
      ArrayOfGridPosPoly gp_p(p_grid.nelem());
      ArrayOfGridPosPoly gp_lat(lat_grid.nelem());
      
      // Interpolate t_field:
      
      // Check that interpolation grids are ok (and throw a detailed
      // error message if not):
      chk_interpolation_grids("Raw temperature to p_grid, 2D case",
                              tfr_p_grid,
                              p_grid,
                              interp_order);
      chk_interpolation_grids("Raw temperature to lat_grid, 2D case",
                              tfr_lat_grid,
                              lat_grid,
                              interp_order);

      // Calculate grid positions:
      p2gridpos_poly( gp_p, tfr_p_grid, p_grid, interp_order );
      gridpos_poly( gp_lat, tfr_lat_grid, lat_grid, interp_order );
            
      // Interpolation weights:
      Tensor3 itw(p_grid.nelem(), lat_grid.nelem(), (interp_order+1)^2);
      // (4 interpolation weights are required for example for linear 2D interpolation)
      interpweights( itw, gp_p, gp_lat);
      
      // Interpolate:
      interp( t_field(joker, joker, 0 ), itw,
              t_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
      
      
      // Interpolate z_field:

      // Check that interpolation grids are ok (and throw a detailed
      // error message if not):
      chk_interpolation_grids("Raw z to p_grid, 2D case",
                              zfr_p_grid,
                              p_grid,
                              interp_order);
      chk_interpolation_grids("Raw z to lat_grid, 2D case",
                              zfr_lat_grid,
                              lat_grid,
                              interp_order);

      // Calculate grid positions:
      p2gridpos_poly( gp_p, zfr_p_grid, p_grid, interp_order );
      gridpos_poly( gp_lat, zfr_lat_grid, lat_grid, interp_order );
            
      // Interpolation weights:
      interpweights( itw, gp_p, gp_lat);
      
      // Interpolate:
      interp( z_field(joker, joker, 0), itw, 
              z_field_raw.data(joker, joker, 0), gp_p, gp_lat);
      
      
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          ostringstream os; 

          if( !( vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID).nelem() != 1 &&
                 vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LON_GRID).nelem() == 1 ))
            {
              os << "VMR data of the " << gas_i << "th species has "
                 << "wrong dimension (1D or 3D). \n";
              throw runtime_error( os.str() );
            }

          // Check that interpolation grids are ok (and throw a detailed
          // error message if not):
          os << "Raw VMR[" << gas_i << "] to p_grid, 2D case";
          chk_interpolation_grids(os.str(),
                                  vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID),
                                  p_grid,
                                  interp_order);
          os.str("");
          os << "Raw VMR[" << gas_i << "] to lat_grid, 2D case";
          chk_interpolation_grids(os.str(),
                                  vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID),
                                  lat_grid,
                                  interp_order);

          // Calculate grid positions:
          p2gridpos_poly(gp_p, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID), 
                         p_grid, interp_order);
          gridpos_poly(gp_lat, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID), 
                       lat_grid, interp_order);
                  
          // Interpolation weights:
          interpweights( itw, gp_p, gp_lat);
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, joker, 0),
                  itw, vmr_field_raw[gas_i].data(joker, joker, 0),
                  gp_p, gp_lat);
        }
    }

  //================================================================
  // atmosphere_dim = 3    
  else if(atmosphere_dim == 3)
    {
      if( tfr_lat_grid.nelem() == 1 &&
          tfr_lon_grid.nelem() == 1 )
        throw runtime_error(
                            "Raw data has wrong dimension. You have to use \n"
                            "AtmFieldsCalcExpand1D instead of AtmFieldsCalc."
                            );

      //Resize variables
      t_field.resize(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem());
      z_field.resize(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem());
      vmr_field.resize(vmr_field_raw.nelem(), p_grid.nelem(), lat_grid.nelem(),
                       lon_grid.nelem());
      
      
      // Gridpositions:
      ArrayOfGridPosPoly gp_p(p_grid.nelem());
      ArrayOfGridPosPoly gp_lat(lat_grid.nelem());
      ArrayOfGridPosPoly gp_lon(lon_grid.nelem());
      
      
      // Interpolate t_field:
      
      // Check that interpolation grids are ok (and throw a detailed
      // error message if not):
      chk_interpolation_grids("Raw temperature to p_grid, 3D case",
                              tfr_p_grid,
                              p_grid,
                              interp_order);
      chk_interpolation_grids("Raw temperature to lat_grid, 3D case",
                              tfr_lat_grid,
                              lat_grid,
                              interp_order);
      chk_interpolation_grids("Raw temperature to lon_grid, 3D case",
                              tfr_lon_grid,
                              lon_grid,
                              interp_order);

      // Calculate grid positions:
      p2gridpos_poly( gp_p, tfr_p_grid, p_grid, interp_order );
      gridpos_poly( gp_lat, tfr_lat_grid, lat_grid, interp_order );
      gridpos_poly( gp_lon, tfr_lon_grid, lon_grid, interp_order );
      
      // Interpolation weights:
      Tensor4 itw(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem(), (interp_order+1)^3);
      // (8 interpolation weights are required for example for linear 3D interpolation)
      interpweights( itw, gp_p, gp_lat, gp_lon );
      
      // Interpolate:
      interp( t_field, itw, t_field_raw.data,  gp_p, gp_lat, gp_lon);
      
      
      // Interpolate z_field:
      
      // Check that interpolation grids are ok (and throw a detailed
      // error message if not):
      chk_interpolation_grids("Raw z to p_grid, 3D case",
                              zfr_p_grid,
                              p_grid,
                              interp_order);
      chk_interpolation_grids("Raw z to lat_grid, 3D case",
                              zfr_lat_grid,
                              lat_grid,
                              interp_order);
      chk_interpolation_grids("Raw z to lon_grid, 3D case",
                              zfr_lon_grid,
                              lon_grid,
                              interp_order);

      // Calculate grid positions:
      p2gridpos_poly( gp_p, zfr_p_grid, p_grid, interp_order );
      gridpos_poly( gp_lat, zfr_lat_grid, lat_grid, interp_order );
      gridpos_poly( gp_lon, zfr_lon_grid, lon_grid, interp_order );
      
      // Interpolation weights:
      interpweights( itw, gp_p, gp_lat, gp_lon );
      
      // Interpolate:
      interp( z_field, itw, z_field_raw.data, gp_p, gp_lat, gp_lon);
      
      
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          ostringstream os; 

          if( !( vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID).nelem() != 1 &&
                 vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LON_GRID).nelem() != 1 ))
            {
              os << "VMR data of the " << gas_i << "th species has "
                 << "wrong dimension (1D or 2D). \n";
              throw runtime_error( os.str() );
            }

          // Check that interpolation grids are ok (and throw a detailed
          // error message if not):
          os << "Raw VMR[" << gas_i << "] to p_grid, 3D case";
          chk_interpolation_grids(os.str(),
                                  vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID),
                                  p_grid,
                                  interp_order);
          os.str("");
          os << "Raw VMR[" << gas_i << "] to lat_grid, 3D case";
          chk_interpolation_grids(os.str(),
                                  vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID),
                                  lat_grid,
                                  interp_order);
          os.str("");
          os << "Raw VMR[" << gas_i << "] to lon_grid, 3D case";
          chk_interpolation_grids(os.str(),
                                  vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LON_GRID),
                                  lon_grid,
                                  interp_order);

          // Calculate grid positions:
          p2gridpos_poly(gp_p, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID), 
                         p_grid, interp_order);
          gridpos_poly(gp_lat, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID), 
                       lat_grid, interp_order);
          gridpos_poly(gp_lon, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LON_GRID), 
                       lon_grid, interp_order);
          
          // Interpolation weights:
          interpweights( itw, gp_p, gp_lat, gp_lon );
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, joker, joker),
                  itw, vmr_field_raw[gas_i].data, gp_p, gp_lat, gp_lon);
        }
    }
  else
  {
    // We can never get here, since there was a runtime 
    // error check for atmosphere_dim at the beginning.
    assert(false);
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsCalcExpand1D(
            Tensor3&              t_field,
            Tensor3&              z_field,
            Tensor4&              vmr_field,
      const Vector&               p_grid,
      const Vector&               lat_grid,
      const Vector&               lon_grid,
      const GriddedField3&        t_field_raw,
      const GriddedField3&        z_field_raw,
      const ArrayOfGriddedField3& vmr_field_raw,
      const Index&                atmosphere_dim,
      const Index&                interp_order )
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  if( atmosphere_dim == 1 )
    { throw runtime_error( 
     "This function is intended for 2D and 3D. For 1D, use *AtmFieldsCalc*.");}

  // Make 1D interpolation using some dummy variables
  Vector    vempty(0);
  Tensor3   t_temp, z_temp;
  Tensor4   vmr_temp;
  AtmFieldsCalc( t_temp, z_temp, vmr_temp, p_grid, vempty, vempty, 
                 t_field_raw, z_field_raw, vmr_field_raw, 1, interp_order );

  // Move values from the temporary tensors to the return arguments
  const Index   np = p_grid.nelem();
  const Index   nlat = lat_grid.nelem();
        Index   nlon = lon_grid.nelem();
  if( atmosphere_dim == 2 )
    { nlon = 1; }
  const Index   nspecies = vmr_temp.nbooks();
  //
  assert( t_temp.npages() == np );
  //
  t_field.resize( np, nlat, nlon );
  z_field.resize( np, nlat, nlon );
  vmr_field.resize( nspecies, np, nlat, nlon );
  //
  for( Index ilon=0; ilon<nlon; ilon++ )
    {
      for( Index ilat=0; ilat<nlat; ilat++ )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              t_field(ip,ilat,ilon) = t_temp(ip,0,0);
              z_field(ip,ilat,ilon) = z_temp(ip,0,0);
              for( Index is=0; is<nspecies; is++ )
                { vmr_field(is,ip,ilat,ilon) = vmr_temp(is,ip,0,0); }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsExpand1D(
            Tensor3&              t_field,
            Tensor3&              z_field,
            Tensor4&              vmr_field,
      const Vector&               p_grid,
      const Vector&               lat_grid,
      const Vector&               lon_grid,
      const Index&                atmosphere_dim )
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Sizes
  const Index   np = p_grid.nelem();
  const Index   nlat = lat_grid.nelem();
  const Index   nlon = max( Index(1), lon_grid.nelem() );
  const Index   nspecies = vmr_field.nbooks();

  if( atmosphere_dim == 1 )
    { throw runtime_error( "No use in calling this method for 1D.");}
  chk_atm_field( "t_field", t_field, 1, p_grid, Vector(0), Vector(0) );
  chk_atm_field( "z_field", z_field, 1, p_grid, Vector(0), Vector(0) );
  if( nspecies )
    chk_atm_field( "vmr_field", vmr_field, 1, nspecies, p_grid, Vector(0), 
                                                                Vector(0) );

  // Temporary containers
  Tensor3 t_temp = t_field, z_temp = z_field;
  Tensor4 vmr_temp = vmr_field;

  // Resize and fill
  t_field.resize( np, nlat, nlon );
  z_field.resize( np, nlat, nlon );
  vmr_field.resize( nspecies, np, nlat, nlon );
  //
  for( Index ilon=0; ilon<nlon; ilon++ )
    {
      for( Index ilat=0; ilat<nlat; ilat++ )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              t_field(ip,ilat,ilon) = t_temp(ip,0,0);
              z_field(ip,ilat,ilon) = z_temp(ip,0,0);
              for( Index is=0; is<nspecies; is++ )
                { vmr_field(is,ip,ilat,ilon) = vmr_temp(is,ip,0,0); }
            }
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsRefinePgrid(// WS Output:
                          Vector& p_grid,
                          Tensor3& t_field,
                          Tensor3& z_field,
                          Tensor4& vmr_field,
                          // WS Input:
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const Index& atmosphere_dim,
                          // Control Parameters:
                          const Numeric& p_step)
{
  // Checks on input parameters:
  
  // We don't actually need lat_grid and lon_grid, but we have them as
  // input variables, so that we can use the standard functions to
  // check atmospheric fields and grids. A bit cheesy, but I don't
  // want to program all the checks explicitly.

  // Check grids:
  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);

  // Check T field:
  chk_atm_field("t_field", t_field, atmosphere_dim,
                p_grid, lat_grid, lon_grid);
 
  // Check z field:
  chk_atm_field("z_field", z_field, atmosphere_dim,
                p_grid, lat_grid, lon_grid);
 
  // Check VMR field (and abs_species):
  chk_atm_field("vmr_field", vmr_field, atmosphere_dim,
                vmr_field.nbooks(), p_grid, lat_grid, lon_grid);

  // Check the keyword arguments:
  if ( p_step <= 0  )
    {
      ostringstream os;
      os << "The keyword argument p_step must be >0.";
      throw runtime_error( os.str() );
    }

  // Ok, all input parameters seem to be reasonable.

  // We will need the log of the pressure grid:
  Vector log_p_grid(p_grid.nelem());
  transform(log_p_grid, log, p_grid);

  //  const Numeric epsilon = 0.01 * p_step; // This is the epsilon that
  //                                         // we use for comparing p grid spacings.

  // Construct abs_p
  // ---------------

  ArrayOfNumeric log_abs_p_a;  // We take log_abs_p_a as an array of
                             // Numeric, so that we can easily 
                             // build it up by appending new elements to the end. 

  // Check whether there are pressure levels that are further apart
  // (in log(p)) than p_step, and insert additional levels if
  // necessary:

  log_abs_p_a.push_back(log_p_grid[0]);

  for (Index i=1; i<log_p_grid.nelem(); ++i)
    {
      const Numeric dp =  log_p_grid[i-1] - log_p_grid[i]; // The grid is descending.

      const Numeric dp_by_p_step = dp/p_step;
      //          cout << "dp_by_p_step: " << dp_by_p_step << "\n";

      // How many times does p_step fit into dp?
      const Index n = (Index) ceil(dp_by_p_step); 
      // n is the number of intervals that we want to have in the
      // new grid. The number of additional points to insert is
      // n-1. But we have to insert the original point as well.
      //          cout << n << "\n";

      const Numeric ddp = dp/(Numeric)n;
      //          cout << "ddp: " << ddp << "\n";

      for (Index j=1; j<=n; ++j)
        log_abs_p_a.push_back(log_p_grid[i-1] - (Numeric)j*ddp);          
    }

  // Copy to a proper vector, we need this also later for
  // interpolation: 
  Vector log_abs_p(log_abs_p_a.nelem());
  for (Index i=0; i<log_abs_p_a.nelem(); ++i)
    log_abs_p[i] = log_abs_p_a[i];

  // Copy the new grid to abs_p, removing the log:
  Vector abs_p(log_abs_p.nelem());
  transform(abs_p, exp, log_abs_p);


  // We will also have to interpolate T and VMR profiles to the new
  // pressure grid. We interpolate in log(p), as usual in ARTS.

  // Grid positions:
  ArrayOfGridPos gp(log_abs_p.nelem());
  gridpos(gp, log_p_grid, log_abs_p);

  // Interpolation weights:
  Matrix itw(gp.nelem(),2);
  interpweights(itw,gp);

  // Extent of latitude and longitude grids:
  Index nlat = lat_grid.nelem();
  Index nlon = lon_grid.nelem();

  // This here is needed so that the method works for 1D and 2D
  // atmospheres as well, not just for 3D:
  if (0 == nlat) nlat = 1;
  if (0 == nlon) nlon = 1;  

  // Define variables for output fields:
  Tensor3 abs_t(log_abs_p.nelem(), nlat, nlon);
  Tensor3 abs_z(log_abs_p.nelem(), nlat, nlon);
  Tensor4 abs_vmr(vmr_field.nbooks(),
                  log_abs_p.nelem(), nlat, nlon);

  for (Index ilat=0; ilat<nlat; ++ilat)
    for (Index ilon=0; ilon<nlon; ++ilon)
      {
        interp(abs_t(joker, ilat, ilon),
               itw,
               t_field(joker, ilat, ilon),
               gp);

        interp(abs_z(joker, ilat, ilon),
               itw,
               z_field(joker, ilat, ilon),
               gp);

        for (Index ivmr=0; ivmr<vmr_field.nbooks(); ++ivmr)
          interp(abs_vmr(ivmr, joker, ilat, ilon),
                 itw,
                 vmr_field(ivmr, joker, ilat, ilon),
                 gp);
      }


  // Copy back the new fields to p_grid, t_field, z_field, vmr_field:
  p_grid    = abs_p;
  t_field   = abs_t;
  z_field   = abs_z; 
  vmr_field = abs_vmr;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AtmRawRead(//WS Output:
                GriddedField3&        t_field_raw,
                GriddedField3&        z_field_raw,
                ArrayOfGriddedField3& vmr_field_raw,
                //WS Input:
                const ArrayOfArrayOfSpeciesTag& abs_species,
                //Keyword:
                const String&   basename)
{
  // Read the temperature field:
  String file_name = basename + ".t.xml";
  xml_read_from_file( file_name, t_field_raw);
  
  out3 << "Temperature field read from file: " << file_name << "\n";  

  // Read geometrical altitude field:
  file_name = basename + ".z.xml";
  xml_read_from_file( file_name, z_field_raw);

  out3 << "Altitude field read from file: " << file_name << "\n";  


  // The species lookup data:

  extern const Array<SpeciesRecord> species_data;
  
  // We need to read one profile for each tag group.
  for ( Index i=0; i<abs_species.nelem(); i ++)
    {
      // Determine the name.
      file_name =
        basename + "." +
        species_data[abs_species[i][0].Species()].Name() + ".xml";
      
      // Add an element for this tag group to the vmr profiles:
      GriddedField3 vmr_field_data;
      vmr_field_raw.push_back(vmr_field_data);
      
      // Read the VMR:
      xml_read_from_file( file_name, vmr_field_raw[vmr_field_raw.nelem()-1]);
      
      // state the source of profile.
      out3 << "  " << species_data[abs_species[i][0].Species()].Name()
           << " profile read from file: " << file_name << "\n";
    }
}
  


/* Workspace method: Doxygen documentation will be auto-generated */
void InterpAtmFieldToRteGps(
                 Numeric&   outvalue,
           const Index&     atmosphere_dim,
           const GridPos&   rte_gp_p,
           const GridPos&   rte_gp_lat,
           const GridPos&   rte_gp_lon,
           const Tensor3&   field )
{
  // Interpolate
  outvalue = interp_atmfield_by_gp( atmosphere_dim, field, 
                                    rte_gp_p, rte_gp_lat, rte_gp_lon );

  out3 << "    Result = " << outvalue << "\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridFromAtmRaw(//WS Output
                      Vector& p_grid,
                      //WS Input
                      const GriddedField3& z_field_raw
                      )
{
  
  Index i=0; 
  while ( z_field_raw.data(i,0,0)< 0.0 ) i++;

  Vector p_grid_raw=z_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  p_grid=p_grid_raw[Range(i,p_grid_raw.nelem()-1)];
}




// A small help function
void z2g(
               Numeric& g,
         const Numeric& r,
         const Numeric& g0,
         const Numeric& z )
{
  g = g0 * pow( r/(r+z), 2 );
}

/* Workspace method: Doxygen documentation will be auto-generated */
void z_fieldFromHSE(
         Tensor3&        z_field,
   const Index&          atmosphere_dim,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Tensor3&        t_field,
   const Tensor4&        vmr_field,
   const Matrix&         r_geoid,
   const Index&          atm_checked,
   const Numeric&        p_hse,
   const Numeric&        z_hse_accuracy )
{
  // Some general variables
  const Index np   = p_grid.nelem();
  const Index nlat = t_field.nrows();
  const Index nlon = t_field.ncols();
  //
  const Index firstH2O = find_first_species_tg( abs_species,
                                      species_index_from_species_name("H2O") );

  // Input checks
  //
  if( !atm_checked )
    throw runtime_error( "The atmosphere must be flagged to have passed a "
                         "consistency check (atm_checked=1)." );
  //
  if( atmosphere_dim == 1  &&  lat_grid.nelem() != 1 )
    { throw runtime_error(
                "The method requires that, for 1D, *lat_grid* has length 1." );
    }
  if( min(lat_grid) < -90  ||  max(lat_grid) > 90 )
    { throw runtime_error(
                       "Values of *lat_grid* must be in the range [-90,90]." );
    }
  //
  if( firstH2O < 0 )
    throw runtime_error(
       "Water vapour is a requiered (must be a tag group in *abs_species*)." );
  //
  if( p_hse>p_grid[0]  ||  p_hse < p_grid[np-1] )
    {
      ostringstream os;
      os << "The value of *p_hse* must be inside the range of *p_grid*:"
         << "  p_hse  = " << p_hse << " Pa\n"
         << "  p_grid = << p_grid[np-1]" << " - " << p_grid[0] << " Pa\n";
      throw runtime_error( os.str() );
    }
  //
  if( z_hse_accuracy <= 0 )
    { throw runtime_error( "The value of *z_hse_accuracy* must be > 0." ); }


  // Determine interpolation weights for p_hse
  //
  ArrayOfGridPos gp(1);
  Matrix itw( 1, 2);
  p2gridpos( gp, p_grid, Vector(1,p_hse) );
  interpweights ( itw, gp );
  

  // The calculations
  //
  for( Index ilat=0; ilat<nlat; ilat++ )
    {
      // "Small g" at geoid altitude, g0:
      // Expression for g0 taken from Wikipedia page "Gravity of Earth", that
      // is stated to be: International Gravity Formula 1967, the 1967 Geodetic
      // Reference System Formula, Helmert's equation or Clairault's formula.
      const Numeric x = fabs( lat_grid[ilat] );
      const Numeric g0 = 9.780327 * ( 1 + 5.3024e-3*pow(sin(DEG2RAD*x),2) + 
                                      5.8e-6*pow(sin(2*DEG2RAD*x),2) );


      for( Index ilon=0; ilon<nlon; ilon++ )
        {
          // Determine altitude for p_hse
          Vector z_hse(1);
          interp( z_hse, itw, z_field(joker,ilat,ilon), gp );

          Numeric z_acc = 2 * z_hse_accuracy;

          while( z_acc > z_hse_accuracy )
            {
              z_acc = 0;
              Numeric g1=g0, g2=g0;

              for( Index ip=0; ip<np-1; ip++ )
                {
                  // Calculate average g
                  if( ip == 0 )
                    {
                      z2g( g2, r_geoid(ilat,ilon), g0, z_field(ip,ilat,ilon) );
                    }
                  g1 = g2;
                  z2g( g2, r_geoid(ilat,ilon), g0, z_field(ip+1,ilat,ilon) );
                  //
                  const Numeric g = ( g1 + g2 ) / 2.0;

                  // Weight mixing ratio for water assuming constant average
                  // molecular weight of the air
                  // 0.3108 = 18/28.96 * 0.5
                  const Numeric r = 0.3108 * ( vmr_field(firstH2O,ip,ilat,ilon)
                                         + vmr_field(firstH2O,ip+1,ilat,ilon) );
  
                  // Virtual temperature (no liquid water)
                  const Numeric tv = 0.5 * ( 1 + 0.61*r) * (
                              t_field(ip,ilat,ilon) + t_field(ip+1,ilat,ilon) );
  
                  // Change in vertical altitude from ip to ip+1 
                  const Numeric dz = 287.053 * (tv/g) * 
                                                 log( p_grid[ip]/p_grid[ip+1] );

                  // New altitude
                  Numeric znew = z_field(ip,ilat,ilon) + dz;
                  z_acc = max( z_acc, fabs(znew-z_field(ip+1,ilat,ilon)) );
                  z_field(ip+1,ilat,ilon) = znew;
                }

              // Adjust to z_hse
              Vector z_tmp(1);
              interp( z_tmp, itw, z_field(joker,ilat,ilon), gp );
              z_field(joker,ilat,ilon) -= z_tmp[0] - z_hse[0];
            }
        }
    }
}




