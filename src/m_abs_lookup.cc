/* Copyright (C) 2002-2007 Stefan Buehler <sbuehler@ltu.se>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   m_abs_lookup.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Wed Nov 20 18:04:20 2002
  
  \brief  Methods related to absorption, lookup table, etc.
*/

#include <algorithm> 
#include <map>
#include "auto_md.h"
#include "arts.h"
#include "messages.h"
#include "gas_abs_lookup.h"
#include "agenda_class.h"
#include "check_input.h"
#include "matpackV.h"
#include "physics_funcs.h"

//! Creates an empty gas absorption lookup table.
/*! 
  This is mainly there to help developers. For example, you can write
  the empty table to an XML file, to see the file format.
*/
  /* param x  Absorption lookup table.*/
void abs_lookupInit(GasAbsLookup& /* x */)
{
  // Nothing to do here.
  // That means, we rely on the default constructor.

  out2 << "  Created an empty gas absorption lookup table.\n";
}


//! Creates a gas absorption lookup table.
// FIXME Doc header!
// - works only with 1D atmospheric fields.
void abs_lookupCreate(// WS Output:
                      GasAbsLookup& gal,
                      Index& abs_lookup_is_adapted,
                      // WS Input:
                      const ArrayOfArrayOfSpeciesTag& abs_species,
                      const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                      const ArrayOfLineshapeSpec&     abs_lineshape,
                      const ArrayOfIndex&             abs_nls,
                      const Vector&                   f_grid,
                      const Vector&                   abs_p,
                      const Matrix&                   abs_vmrs,
                      const Vector&                   abs_t,
                      const Vector&                   abs_t_pert,
                      const Vector&                   abs_nls_pert, 
                      const Vector&                   abs_n2,            
                      const Vector&                   abs_h2o,           
                      const ArrayOfString&            abs_cont_names,    
                      const ArrayOfString&            abs_cont_models,   
                      const ArrayOfVector&            abs_cont_parameters )
{
  // We will be calling an absorption agenda one species at a
  // time. This is better than doing all simultaneously, because is
  // saves memory and allows for consistent treatment of nonlinear
  // species. But it means we need local copies of species, line list,
  // and line shapes for agenda communication.
  
  // 1. Output of absorption calculations:

  // Absorption coefficients:
  Matrix these_abs_coef;

  // "Cross sections" per tag group. These are not real
  // cross-sections, but the old definition, abs/VMR!
  ArrayOfMatrix abs_xsec_per_species;

  // 2. Determine various important sizes:
  const Index n_species = abs_species.nelem();   // Number of abs species
  const Index n_nls = abs_nls.nelem();           // Number of nonlinear species
  const Index n_f_grid = f_grid.nelem();         // Number of frequency grid points
  const Index n_p_grid = abs_p.nelem();         // Number of presure grid points
  const Index n_t_pert = abs_t_pert.nelem();     // Number of temp. perturbations
  const Index n_nls_pert = abs_nls_pert.nelem(); // Number of VMR pert. for NLS

  // 3. Input to absorption calculations:

  // Absorption vmrs and temperature:
  Matrix this_vmr(1,n_p_grid);
  Vector this_t(n_p_grid);

  // Species list, lines, and line shapes, all with only 1 element:
  ArrayOfArrayOfSpeciesTag this_species(1);
  ArrayOfArrayOfLineRecord these_lines(1);
  ArrayOfLineshapeSpec this_lineshape(1);

  // Local copy of nls_pert and t_pert:
  Vector these_nls_pert;        // Is resized every time when used
  Vector these_t_pert;          // Is resized later on

  // 3.a Set up a logical array for the nonlinear species
  ArrayOfIndex non_linear(n_species,0);
  for ( Index s=0; s<n_nls; ++s )
    {
      non_linear[abs_nls[s]] = 1;
    }

  // 4. Checks of input parameter correctness:

  // abs_species, f_grid, and p_grid should not be empty:
  if ( 0==n_species ||
       0==n_f_grid ||
       0==n_p_grid )
    {
      ostringstream os;
      os << "One of the required input variables is empty:\n"
         << "abs_species.nelem() = " << n_species << ",\n"
         << "f_grid.nelem() = " << n_f_grid << ",\n"
         << "abs_p.nelem() = " << n_p_grid << ".";
      throw runtime_error( os.str() );
    }

  // abs_nls must not be longer than abs_species, strictly increasing,
  // and point inside abs_species: 
  if (0!=n_nls)
    {
      chk_if_in_range( "abs_nls.nelem()",
                       n_nls, 
                       0, 
                       n_species );

      chk_if_increasing( "abs_nls",
                         abs_nls ); 

      chk_if_in_range( "min(abs_nls)",
                       min(abs_nls), 
                       0, 
                       n_species-1 );

      chk_if_in_range( "max(abs_nls)",
                       max(abs_nls), 
                       0, 
                       n_species-1 );
    }

  // VMR matrix must match species list and pressure levels:
  chk_size( "abs_vmrs",
            abs_vmrs,
            n_species,
            n_p_grid );

  // Temperature vector must match number of pressure levels:
  chk_size( "abs_t",
            abs_t,
            n_p_grid ); 

  // abs_nls_pert should only be not-empty if we have nonlinear species:
  if ( ( 0==n_nls && 0 != n_nls_pert ) ||
       ( 0!=n_nls && 0 == n_nls_pert ))
    {
      ostringstream os;
      os << "You have to set both abs_nls and abs_nls_pert, or none.";
      throw runtime_error( os.str() );
    }


  // 5. Set general lookup table properties:
  gal.species = abs_species;    // Species list
  gal.nonlinear_species = abs_nls;  // Nonlinear species (H2O)
  gal.f_grid = f_grid;           // Frequency grid
  gal.p_grid = abs_p;          // Pressure grid
  gal.vmrs_ref = abs_vmrs;
  gal.t_ref = abs_t;
  gal.t_pert = abs_t_pert;
  gal.nls_pert = abs_nls_pert;

  // 6. Create gal.xsec with the right dimensions:
  {
    Index a,b,c,d;

    if ( 0 == n_t_pert ) a = 1;
    else a = n_t_pert;

    b = n_species + n_nls * ( n_nls_pert - 1 );

    c = n_f_grid;

    d = n_p_grid;

    gal.xsec.resize( a, b, c, d );
  }


  // 6.a. Set up these_t_pert. This is done so that we can use the
  // same loop over the perturbations, independent of
  // whether we have temperature perturbations or not.
  if ( 0!=n_t_pert)
    {
      out2 << "  With temperature perturbations.\n";
      these_t_pert.resize(n_t_pert);
      these_t_pert = abs_t_pert;
    }
  else
    {
      out2 << "  No temperature perturbations.\n";
      these_t_pert.resize(1);
      these_t_pert = 0;
    }
  

  // 7. Now we have to fill gal.xsec with the right values!

  // Loop species:
  for ( Index i=0,spec=0; i<n_species; ++i )
    {
      // spec is the index for the second dimension of gal.xsec.
      
      // Prepare absorption agenda input for this species:
      out2 << "  Doing species " << i+1 << " of " << n_species << ": "
           << abs_species[i] << ".\n";

      // Get a dummy list of tag groups with only the current element:
      this_species[0].resize(abs_species[i].nelem());
      this_species[0] = abs_species[i];

      // List of lines:
      these_lines[0].resize(abs_lines_per_species[i].nelem());
      these_lines[0] = abs_lines_per_species[i];
      
      // List of lineshapes:
      this_lineshape[0] = abs_lineshape[i];

      // Set up these_nls_pert. This is done so that we can use the
      // same loop over the perturbations, independent of
      // whether we have nonlinear species or not.
      if ( non_linear[i] )
        {
          out2 << "  This is a species with non-linear treatment.\n";
          these_nls_pert.resize(n_nls_pert);
          these_nls_pert = abs_nls_pert;
        }
      else
        {
          these_nls_pert.resize(1);
          these_nls_pert = 1;
        }
      
      // Loop these_nls_pert:
      for ( Index s=0; s<these_nls_pert.nelem(); ++s,++spec )
        {
          // Remeber, spec is the index for the second dimension of gal.xsec
          
          if ( non_linear[i] )
            {
              out2 << "  Doing VMR variant " << s << " of " << n_nls_pert << ": "
                   << abs_nls_pert[s] << ".\n";
            }

          // VMR for this species:
          this_vmr(0,joker) = abs_vmrs(i,joker);
          this_vmr(0,joker) *= these_nls_pert[s]; // Add perturbation



          // Loop temperature perturbations:
          for ( Index j=0; j<these_t_pert.nelem(); ++j )
            {
              if ( 0!=n_t_pert )
                {
                  out2 << "  Doing temperature variant " << j+1
                       << " of " << n_t_pert << ": "
                       << these_t_pert[j] << ".\n";
                }
              
              // Create perturbed temperature profile:
              this_t = gal.t_ref;
              this_t += these_t_pert[j];
      
              // The sequence of function calls here is inspired from
              // abs_coefCalcSaveMemory. 

              abs_xsec_per_speciesInit( abs_xsec_per_species, this_species,
                                        f_grid, abs_p );

              abs_xsec_per_speciesAddLines( abs_xsec_per_species,
                                            this_species,
                                            f_grid,
                                            abs_p,
                                            this_t,
                                            abs_h2o,
                                            this_vmr,
                                            these_lines,
                                            this_lineshape );

              abs_xsec_per_speciesAddConts( abs_xsec_per_species,
                                            this_species,
                                            f_grid,
                                            abs_p,
                                            this_t,
                                            abs_n2,
                                            abs_h2o,
                                            this_vmr,
                                            abs_cont_names,
                                            abs_cont_parameters,
                                            abs_cont_models);

              // Store in the right place:
              // Loop through all altitudes
              for ( Index p=0; p<n_p_grid; ++p )
                {
                  // We have to multiply with VMR here, to get
                  // absorption coefficient from the arts-1-0
                  // "xsec" (which are not proper cross sections).
                  // Then divide by n*VMR to get the lookup table
                  // cross sections. The VMR cancels out, so we
                  // just divide by n (total number density).
                  
                  // IMPORTANT: There was a bug in my old Matlab
                  // function "create_lookup.m" to generate the lookup
                  // table. (The number density was always the
                  // reference one, and did not change for different
                  // temperatures.) Patricks Atmlab function
                  // "arts_abstable_from_arts1.m" did *not* have this bug.

                  // Calculate the number density for the given pressure and
                  // temperature: 
                  // n = n0*T0/p0 * p/T or n = p/kB/t, ideal gas law
                  const Numeric n = number_density( gal.p_grid[p],
                                                    this_t[p]   );

                  gal.xsec( j, spec, Range(joker), p )
                    = abs_xsec_per_species[0](Range(joker),p);

                  gal.xsec( j, spec, Range(joker), p ) /= n;
                }
            }
        }
    }

  // Set the abs_lookup_is_adapted flag. After all, the table fits the
  // current frequency grid and species selection.
  abs_lookup_is_adapted = 1;
}

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Stefan Buehler
   \date   2007-05-21
*/
void abs_nlsSet(// WS Output:
                ArrayOfIndex& abs_nls,
                // WS Input:
                const ArrayOfArrayOfSpeciesTag& abs_species,
                // Control Parameters:
                const ArrayOfString& names)
{
  Index n_nls = names.nelem();    // Number of non-linear species 

  if (0==n_nls)
    {
      out3 << "  No species with non-linear treatment.\n";
      abs_nls.resize(n_nls);
    }
  else
    {
      get_tagindex_for_Strings(abs_nls,
                               abs_species,
                               names);
    }

  // FIXME: Remove this.
//   // This loop is only executed if we do have non-linear species:
//   for ( Index i=0; i<n_nls; ++i )
//     {
//       ArrayOfSpeciesTag this_nls; // Tags for this non-linear species

//       out3 << " Non-linear treatment for " << names[i] << ".\n";
//       // Convert strings to species tags:
//       array_species_tag_from_string( this_nls, names[i] );  
      
//       // Check if this species is contained in abs_species:
//       get_tag_group_index_for_tag_group(abs_nls[i],
//                                         abs_species,
//                                         nls_species[]);
//     }

}

void abs_speciesAdd(// WS Output:
                    ArrayOfArrayOfSpeciesTag& abs_species,
                    // Control Parameters:
                    const ArrayOfString& names)
{
  // Size of initial array
  Index n_gs = abs_species.nelem();
  
  // Temporary ArrayOfSpeciesTag
  ArrayOfSpeciesTag temp;
    
  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for ( Index i=0; i<names.nelem(); ++i )
    {
      array_species_tag_from_string( temp, names[i] );  
      abs_species.push_back(temp);
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Added tag groups:";
  for ( Index i=n_gs; i<abs_species.nelem(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( Index s=0; s<abs_species[i].nelem(); ++s )
        {
          out3 << " " << abs_species[i][s].Name();
        }
    }
  out3 << '\n';
}



//! jacobianAddGas
/*!
   See the online help (arts -d FUNCTION_NAME)
   
   \author Patrick Eriksson
   \date   2006-08-29
*/
void abs_speciesAdd2(// WS Output:
                    ArrayOfArrayOfSpeciesTag& abs_species,
                    ArrayOfRetrievalQuantity& jq,
                    Agenda&                   jacobian_agenda,
                    // WS Input:
                    const Matrix&             jac,
                    const Index&              atmosphere_dim,
                    const Vector&             p_grid,
                    const Vector&             lat_grid,
                    const Vector&             lon_grid,
                    // WS Generic Input:
                    const Vector&             rq_p_grid,
                    const Vector&             rq_lat_grid,
                    const Vector&             rq_lon_grid,
                    // WS Generic Input Names:
                    const String&             rq_p_grid_name,
                    const String&             rq_lat_grid_name,
                    const String&             rq_lon_grid_name,
                    // Control Parameters:
                    const String&             species,
                    const String&             method,
                    const String&             mode,
                    const Numeric&            dx)
{
  // Add species to *abs_species*
  ArrayOfSpeciesTag tags;
  array_species_tag_from_string( tags, species );
  abs_species.push_back( tags );

  // Print list of added tag group to the most verbose output stream:
  out3 << "  Appended tag group:";
  out3 << "\n  " << abs_species.nelem()-1 << ":";
  for ( Index s=0; s<tags.nelem(); ++s )
  {
    out3 << " " << tags[s].Name();
  }
  out3 << '\n';

  // Do retrieval part
  jacobianAddAbsSpecies( jq, jacobian_agenda, jac, atmosphere_dim, 
                         p_grid, lat_grid, lon_grid, rq_p_grid, rq_lat_grid, 
                         rq_lon_grid, rq_p_grid_name, rq_lat_grid_name, 
                         rq_lon_grid_name, species, method, mode, dx);
}



void abs_speciesInit( ArrayOfArrayOfSpeciesTag& abs_species )
{
  abs_species.resize(0);
}



void abs_speciesSet(// WS Output:
                    ArrayOfArrayOfSpeciesTag& abs_species,
                    // Control Parameters:
                    const ArrayOfString& names)
{
  abs_species.resize(names.nelem());

  //cout << "Names: " << names << "\n";

  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for ( Index i=0; i<names.nelem(); ++i )
    {
      // This part has now been moved to array_species_tag_from_string.
      // Call this function.
      array_species_tag_from_string( abs_species[i], names[i] );  
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Defined tag groups:";
  for ( Index i=0; i<abs_species.nelem(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( Index s=0; s<abs_species[i].nelem(); ++s )
        {
          out3 << " " << abs_species[i][s].Name();
        }
    }
  out3 << '\n';
}

void abs_lookupAdapt( GasAbsLookup&                   abs_lookup,
                          Index&                          abs_lookup_is_adapted,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const Vector&                   f_grid)
{
  abs_lookup.Adapt( abs_species, f_grid );
  abs_lookup_is_adapted = 1;
}

void abs_scalar_gasExtractFromLookup( Matrix&             abs_scalar_gas,
                                      const GasAbsLookup& abs_lookup,
                                      const Index&        abs_lookup_is_adapted, 
                                      const Index&        f_index,
                                      const Numeric&      a_pressure,
                                      const Numeric&      a_temperature,
                                      const Vector&       a_vmr_list)
{
  // Check if the table has been adapted:
  if ( 1!=abs_lookup_is_adapted )
    throw runtime_error("Gas absorption lookup table must be adapted,\n"
                        "use method abs_lookupAdapt.");

  // The function we are going to call here is one of the few helper
  // functions that adjust the size of their output argument
  // automatically. 
  abs_lookup.Extract( abs_scalar_gas,
                          f_index,
                          a_pressure,
                          a_temperature,
                          a_vmr_list );
}

//! Calculate scalar gas absorption for all points in the atmosphere.
/*! 
  This is mainly for testing and plotting gas absorption. For RT
  calculations, gas absorption is calculated or extracted locally,
  therefore there is no need to calculate a global field. But this
  method is handy for easy plotting of absorption vs. pressure, for
  example.

  The calculation itself is performed by the
  *abs_scalar_gas_agenda*, which needs the input variables
  *a_pressure*, *a_temperature*, and *a_vmr_list*, and returns the
  output variable *abs_scalar_gas*.

  \param asg_field      Output: Scalar gas absorption field.

  \param asg            Agenda output: Local scalar gas absorption.
  \param a_pressure     Agenda input: Local pressure.
  \param a_temperature  Agenda input: Local temperature.
  \param a_vmr_list     Agenda input: Local list of VMR values.

  \param sga_agenda     Agenda to use to calculate local absorption.
  \param f_index        FIXME: Add documentation.
  \param f_grid         Frequency grid.
  \param atmosphere_dim Atmospheric dimensionality.
  \param p_grid         Global pressure grid.
  \param lat_grid       Global latitude grid.
  \param lon_grid       Global longitude grid.
  \param t_field        Global temperature field.
  \param vmr_field      Global VMR fields.

  \author Stefan Buehler
  \date   2002-12-20
*/
void abs_fieldCalc(// WS Output:
                   Tensor5& asg_field,
                   // WS Input:
                   const Agenda&  sga_agenda,
                   const Index&   f_index,
                   const Vector&  f_grid,
                   const Index&   atmosphere_dim,
                   const Vector&  p_grid,
                   const Vector&  lat_grid,
                   const Vector&  lon_grid,
                   const Tensor3& t_field,
                   const Tensor4& vmr_field )
{
  Matrix  asg;
  Numeric a_pressure;
  Numeric a_temperature;
  Vector a_vmr_list;
  // Get the number of species from the leading dimension of vmr_field:
  const Index n_species = vmr_field.nbooks();

  // Number of frequencies:
  const Index n_frequencies = f_grid.nelem();

  // Number of pressure levels:
  const Index n_pressures = p_grid.nelem();

  // Number of latitude grid points (must be at least one):
  const Index n_latitudes = max( Index(1), lat_grid.nelem() );

  // Number of longitude grid points (must be at least one):
  const Index n_longitudes = max( Index(1), lon_grid.nelem() );
  
  // Check grids:
  chk_atm_grids( atmosphere_dim,
                 p_grid,
                 lat_grid,
                 lon_grid );
  
  // Check if t_field is ok:
  chk_atm_field( "t_field",
                 t_field,
                 atmosphere_dim,
                 p_grid,
                 lat_grid,
                 lon_grid );

  // Check if vmr_field is ok.
  // (Actually, we are not checking the first dimension, since
  // n_species has been set from this.)
  chk_atm_field( "vmr_field",
                 vmr_field,
                 atmosphere_dim,
                 n_species,
                 p_grid,
                 lat_grid,
                 lon_grid );

  // We also set the start and extent for the frequency loop.
  Index f_extent;

  if ( f_index < 0 )
    {
      // This means we should extract for all frequencies.

      f_extent = n_frequencies;
    }
  else
    {
      // This means we should extract only for one frequency.

      // Make sure that f_index is inside f_grid:
      if ( f_index >= n_frequencies )
        {
          ostringstream os;
          os << "The frequency index f_index points to a frequency outside"
             << "the frequency grid. (f_index = " << f_index
             << ", n_frequencies = " << n_frequencies << ")";
          throw runtime_error( os.str() );
        }

      f_extent = 1;
    }

  // Resize output field.
  // The dimension in lat and lon must be at least one, even if these
  // grids are empty.
  out2 << "  Creating field with dimensions:\n"
       << "    " << n_species << "    gas species,\n"
       << "    " << f_extent << "     frequencies,\n"
       << "    " << n_pressures << "  pressures,\n"
       << "    " << n_latitudes << "  latitudes,\n"
       << "    " << n_longitudes << " longitudes.\n";

  asg_field.resize( n_species,
                    f_extent,
                    n_pressures,
                    n_latitudes,
                    n_longitudes );

  // Flag for first time agenda output:
  Index count = 0;

  // Now we have to loop all points in the atmosphere:
  for ( Index ipr=0; ipr<n_pressures; ++ipr )         // Pressure:  ipr
    {
      a_pressure = p_grid[ipr];

      out3 << "  p_grid[" << ipr << "] = " << a_pressure << "\n";

      for ( Index ila=0; ila<n_latitudes; ++ila )   // Latitude:  ila
        for ( Index ilo=0; ilo<n_longitudes; ++ilo ) // Longitude: ilo
          {
            a_temperature = t_field( ipr, ila, ilo );
            a_vmr_list    = vmr_field( Range(joker),
                                       ipr, ila, ilo );

            // Execute agenda to calculate local absorption.
            // Agenda input:  f_index, a_pressure, a_temperature, a_vmr_list
            // Agenda output: asg
            abs_scalar_gas_agendaExecute (asg, f_index, a_pressure,
                                                 a_temperature, a_vmr_list,
                                                 sga_agenda, (count != 0));

            // Verify, that the number of species in asg is
            // constistent with vmr_field:
            if ( n_species != asg.ncols() )
              {
                ostringstream os;
                os << "The number of gas species in vmr_field is "
                   << n_species << ",\n"
                   << "but the number of species returned by the agenda is "
                   << asg.ncols() << ".";
                throw runtime_error( os.str() );
              }

            // Verify, that the number of frequencies in asg is
            // constistent with f_extent:
            if ( f_extent != asg.nrows() )
              {
                ostringstream os;
                os << "The number of frequencies desired is "
                   << n_frequencies << ",\n"
                   << "but the number of frequencies returned by the agenda is "
                   << asg.nrows() << ".";
                throw runtime_error( os.str() );
              }

            // Store the result in output field.
            // We have to transpose asg, because the dimensions of the
            // two variables are:
            // asg_field: [ abs_species, f_grid, p_grid, lat_grid, lon_grid]
            // asg:       [ f_grid, abs_species ]
            asg_field( Range(joker),
                       Range(joker),
                       ipr, ila, ilo ) = transpose( asg );
            
            ++count;
          }
    }
}


//! f_gridFromGasAbsLookup
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2004-09-15
*/
void f_gridFromGasAbsLookup(
             Vector&         f_grid,
       const GasAbsLookup&   abs_lookup )
{
  abs_lookup.GetFgrid( f_grid );
}



//! p_gridFromGasAbsLookup
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2004-09-15
*/
void p_gridFromGasAbsLookup(
             Vector&         p_grid,
       const GasAbsLookup&   abs_lookup )
{
  abs_lookup.GetPgrid( p_grid );
}
