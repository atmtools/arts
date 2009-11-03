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
#include "xml_io.h"

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_P_GRID;
extern const Index GFIELD4_LAT_GRID;
extern const Index GFIELD4_LON_GRID;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


// Workspace method, doxygen header will be auto-generated.
// 2007-07-25 Stefan Buehler
void atm_fields_compactFromMatrix(// WS Output:
                                  GField4& af, // atm_fields_compact
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
      os << "Cannot copy Matrix.\n"
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
  af(Range(joker),Range(joker),0,0) = transpose(im(Range(joker),Range(1,nf_1)));
}



// Workspace method, doxygen header is auto-generated.
// 2007-07-31 Stefan Buehler
void atm_fields_compactAddConstant(// WS Output:
                                   GField4& af,
                                   // Control Parameters:
                                   const String& name,
                                   const Numeric& value)
{
  // Number of fields already present:
  const Index nf = af.get_string_grid(GFIELD4_FIELD_NAMES).nelem();

  if (0==nf)
    {
      ostringstream os;
      os << "The *atm_fields_compact* must already contain at least one field,\n"
         << "so that we can infer the dimensions from that.";
      throw runtime_error( os.str() );
    }

  // Add name of new field to field name list:
  af.get_string_grid(GFIELD4_FIELD_NAMES).push_back(name);

  // Save original fields:
  const Tensor4 dummy = af;

  // Adjust size:
  af.resize( nf+1, dummy.npages(), dummy.nrows(), dummy.ncols() );

  // Copy back original field:
  af( Range(0,nf), Range(joker), Range(joker), Range(joker) ) = dummy;
  
  // Add the constant value:
  af( nf, Range(joker), Range(joker), Range(joker) ) = value;
}


// Workspace method, doxygen header is auto-generated.
void batch_atm_fields_compactFromArrayOfMatrix(// WS Output:
                                               ArrayOfGField4& batch_atm_fields_compact,
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
                          const GField4& atm_fields_compact,
                          const Index&  atmosphere_dim )
{
  // Make a handle on atm_fields_compact to save typing:
  const GField4& c = atm_fields_compact;
  
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
  // T[K] z[m] VMR_1[1] ... VMR[2]

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
  t_field = c(0,Range(joker),Range(joker),Range(joker));

  // Altitude profile (second field):
  z_field.resize(np,nlat,nlon);
  z_field = c(1,Range(joker),Range(joker),Range(joker));

  // VMR profiles (remaining fields):
  vmr_field.resize(ns,np,nlat,nlon);
  vmr_field = c(Range(2,ns),Range(joker),Range(joker),Range(joker));
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
              Vector&   lon_grid,
              Numeric&  lat_1d,
              Numeric&  meridian_angle_1d )
{
  out2 << "  Sets the atmospheric dimensionality to 2.\n";
  out3 << "    atmosphere_dim = 2\n";
  out3 << "    lon_grid is set to be an empty vector\n";
  out3 << "    lat_1d = -999\n";
  out3 << "    meridian_angle_1d = -999\n";
  atmosphere_dim = 2;
  lon_grid.resize(0);
  lat_1d = -999;
  meridian_angle_1d = -999;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmosphereSet3D(
        // WS Output:
              Index&    atmosphere_dim,
              Numeric&  latitude_1d,
              Numeric&  meridian_angle_1d )
{
  out2 << "  Sets the atmospheric dimensionality to 3.\n";
  out3 << "    atmosphere_dim = 3\n";
  out3 << "    lat_1d = -999\n";
  out3 << "    meridian_angle_1d = -999\n";
  atmosphere_dim = 3;
  latitude_1d = -999;
  meridian_angle_1d = -999;
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
                   const GField3&        t_field_raw,
                   const GField3&        z_field_raw,
                   const ArrayOfGField3& vmr_field_raw,
                   const Index&          atmosphere_dim
                   )
{
  const ConstVectorView tfr_p_grid = t_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView tfr_lat_grid = t_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView tfr_lon_grid = t_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
  const ConstVectorView zfr_p_grid = z_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView zfr_lat_grid = z_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView zfr_lon_grid = z_field_raw.get_numeric_grid(GFIELD3_LON_GRID);

  // Basic checks of input variables
  //
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  

  // Note that we are using the special function p2gridpos below for
  // all pressure interpolations. This does the usual ARTS pressure
  // interpolation: Linear in log(p). We don't have to take logs here
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
      ArrayOfGridPos gp_p(p_grid.nelem());
  
      // Interpolate t_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, tfr_p_grid, p_grid );

      // Interpolation weights:
      Matrix itw(p_grid.nelem(), 2);
      // (2 interpolation weights are required for 1D interpolation)
      interpweights( itw, gp_p);
  
      // Interpolate:
      interp( t_field(joker, 0, 0), itw, 
              t_field_raw(joker, 0, 0),  gp_p);

  
      // Interpolate z_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, zfr_p_grid, p_grid );
     
      // Interpolation weights:
      interpweights( itw, gp_p );
      
      // Interpolate:
      interp( z_field(joker, 0, 0), itw,
              z_field_raw(joker, 0, 0), gp_p);
      
  
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          if( !( vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID).nelem() == 1 &&
                 vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LON_GRID).nelem() == 1 ))
            {
              ostringstream os; 
              os << "VMR data of the " << gas_i << "th species has "
                 << "wrong dimension (2D or 3D). \n";
              throw runtime_error( os.str() );
            }
          
          // Calculate grid positions:
          p2gridpos(gp_p, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID), p_grid);
          
          // Interpolation weights:
          interpweights( itw, gp_p);
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, 0, 0),
                  itw, vmr_field_raw[gas_i](joker, 0, 0), gp_p);
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
      ArrayOfGridPos gp_p(p_grid.nelem());
      ArrayOfGridPos gp_lat(lat_grid.nelem());
      
      // Interpolate t_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, tfr_p_grid, p_grid );
      gridpos( gp_lat, tfr_lat_grid, lat_grid );
            
      // Interpolation weights:
      Tensor3 itw(p_grid.nelem(), lat_grid.nelem(), 4);
      // (8 interpolation weights are required for 3D interpolation)
      interpweights( itw, gp_p, gp_lat);
      
      // Interpolate:
      interp( t_field(joker, joker, 0 ), itw,
              t_field_raw(joker, joker, 0),  gp_p, gp_lat);
      
      
      // Interpolate z_field:
      // Calculate grid positions:
      p2gridpos( gp_p, zfr_p_grid, p_grid );
      gridpos( gp_lat, zfr_lat_grid, lat_grid );
            
      // Interpolation weights:
      interpweights( itw, gp_p, gp_lat);
      
      // Interpolate:
      interp( z_field(joker, joker, 0), itw, 
              z_field_raw(joker, joker, 0), gp_p, gp_lat);
      
      
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID), p_grid);
          gridpos(gp_lat, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID), lat_grid);
                  
          // Interpolation weights:
          interpweights( itw, gp_p, gp_lat);
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, joker, 0),
                  itw, vmr_field_raw[gas_i](joker, joker, 0),
                  gp_p, gp_lat);
        }
    }

  //================================================================
  // atmosphere_dim = 3    
  else
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
      ArrayOfGridPos gp_p(p_grid.nelem());
      ArrayOfGridPos gp_lat(lat_grid.nelem());
      ArrayOfGridPos gp_lon(lon_grid.nelem());
      
      
      // Interpolate t_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, tfr_p_grid, p_grid );
      gridpos( gp_lat, tfr_lat_grid, lat_grid );
      gridpos( gp_lon, tfr_lon_grid, lon_grid );
      
      // Interpolation weights:
      Tensor4 itw(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem(), 8);
      // (8 interpolation weights are required for 3D interpolation)
      interpweights( itw, gp_p, gp_lat, gp_lon );
      
      // Interpolate:
      interp( t_field, itw, t_field_raw,  gp_p, gp_lat, gp_lon);
      
      
      // Interpolate z_field:
      
      // Calculate grid positions:
      p2gridpos( gp_p, zfr_p_grid, p_grid );
      gridpos( gp_lat, zfr_lat_grid, lat_grid );
      gridpos( gp_lon, zfr_lon_grid, lon_grid );
      
      // Interpolation weights:
      interpweights( itw, gp_p, gp_lat, gp_lon );
      
      // Interpolate:
      interp( z_field, itw, z_field_raw, gp_p, gp_lat, gp_lon);
      
      
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          // Calculate grid positions:
          p2gridpos(gp_p, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID), p_grid);
          gridpos(gp_lat, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID), lat_grid);
          gridpos(gp_lon, vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LON_GRID), lon_grid);
          
          // Interpolation weights:
          interpweights( itw, gp_p, gp_lat, gp_lon );
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, joker, joker),
                  itw, vmr_field_raw[gas_i], gp_p, gp_lat, gp_lon);
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsCalcExpand1D(
            Tensor3&        t_field,
            Tensor3&        z_field,
            Tensor4&        vmr_field,
      const Vector&         p_grid,
      const Vector&         lat_grid,
      const Vector&         lon_grid,
      const GField3&        t_field_raw,
      const GField3&        z_field_raw,
      const ArrayOfGField3& vmr_field_raw,
      const Index&          atmosphere_dim )
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
                                  t_field_raw, z_field_raw, vmr_field_raw, 1 );

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
                GField3&        t_field_raw,
                GField3&        z_field_raw,
                ArrayOfGField3& vmr_field_raw,
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
      GField3 vmr_field_data;
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
           const Vector&    p_grid,
           const Vector&    lat_grid,
           const Vector&    lon_grid,
           const GridPos&   rte_gp_p,
           const GridPos&   rte_gp_lat,
           const GridPos&   rte_gp_lon,
           const Tensor3&   field)
{
  // Interpolate
  outvalue = interp_atmfield_by_gp( atmosphere_dim, p_grid, lat_grid, 
                lon_grid, field, rte_gp_p, rte_gp_lat, rte_gp_lon );

  out3 << "    Result = " << outvalue << "\n";
}


