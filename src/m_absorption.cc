/*!
  \file   m_absorption.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Wed Nov 20 18:04:20 2002
  
  \brief  Methods related to absorption, lookup table, etc.
*/

#include <algorithm> 
#include "arts.h"
#include "messages.h"
#include "gas_abs_lookup.h"
#include "agenda_class.h"
#include "check_input.h"
#include "matpackV.h"

//! Creates an empty gas absorption lookup table.
/*! 
  This is mainly there to help developers. For example, you can write
  the empty table to an XML file, to see the file format.

  \param GasAbsLookup Absorption lookup table.
*/
void gas_abs_lookupInit(GasAbsLookup& x)
{
  // Nothing to do here.
  // That means, we rely on the default constructor.

  out2 << "  Created an empty gas absorption lookup table.\n";
}

void gas_speciesSet(// WS Output:
                    ArrayOfArrayOfSpeciesTag& gas_species,
                    // Control Parameters:
                    const ArrayOfString& names)
{
  gas_species.resize(names.nelem());

  //cout << "Names: " << names << "\n";

  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for ( Index i=0; i<names.nelem(); ++i )
    {
      // There can be a comma separated list of tag definitions, so we
      // need to break the String apart at the commas.
      ArrayOfString tag_def;

      bool go_on = true;
      String these_names = names[i];
      while (go_on)
        {
          //          Index n = find_first( these_names, ',' );
          Index n = these_names.find(',');
          if ( n == these_names.npos ) // Value npos indicates not found.
            {
              // There are no more commas.
              //              cout << "these_names: (" << these_names << ")\n";
              tag_def.push_back(these_names);
              go_on = false;
            }
          else
            {
              tag_def.push_back(these_names.substr(0,n));
              these_names.erase(0,n+1);
            }
        }

      // tag_def now holds the different tag Strings for this group.
      //      cout << "tag_def =\n" << tag_def << endl;

      // Set size to zero, in case the method has been called before.
      gas_species[i].resize(0);

      for ( Index s=0; s<tag_def.nelem(); ++s )
        {
          SpeciesTag this_tag(tag_def[s]);

          // Safety check: For s>0 check that the tags belong to the same species.
          if (s>0)
            if ( gas_species[i][0].Species() != this_tag.Species() )
              throw runtime_error("Tags in a tag group must belong to the same species.");

          gas_species[i].push_back(this_tag);
        }
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Defined tag groups:";
  for ( Index i=0; i<gas_species.nelem(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( Index s=0; s<gas_species[i].nelem(); ++s )
        {
          out3 << " " << gas_species[i][s].Name();
        }
    }
  out3 << '\n';
}

void gas_abs_lookupAdapt( GasAbsLookup&                   gas_abs_lookup,
                          Index&                          gas_abs_lookup_is_adapted,
                          const ArrayOfArrayOfSpeciesTag& gas_species,
                          const Vector&                   f_grid)
{
  gas_abs_lookup.Adapt( gas_species, f_grid );
  gas_abs_lookup_is_adapted = 1;
}

void abs_scalar_gasExtractFromLookup( Matrix&             abs_scalar_gas,
                                      const GasAbsLookup& gas_abs_lookup,
                                      const Index&        gas_abs_lookup_is_adapted, 
                                      const Index&        f_index,
                                      const Numeric&      a_pressure,
                                      const Numeric&      a_temperature,
                                      const Vector&       a_vmr_list)
{
  // Check if the table has been adapted:
  if ( 1!=gas_abs_lookup_is_adapted )
    throw runtime_error("Gas absorption lookup table must be adapted,\n"
                        "use method gas_abs_lookupAdapt.");

  // The function we are going to call here is one of the few helper
  // functions that adjust the size of their output argument
  // automatically. 
  gas_abs_lookup.Extract( abs_scalar_gas,
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
  *scalar_gas_absorption_agenda*, which needs the input variables
  *a_pressure*, *a_temperature*, and *a_vmr_list*, and returns the
  output variable *abs_scalar_gas*.

  \param asg_field      Output: Scalar gas absorption field.

  \param asg            Agenda output: Local scalar gas absorption.
  \param a_pressure     Agenda input: Local pressure.
  \param a_temperature  Agenda input: Local temperature.
  \param a_vmr_list     Agenda input: Local list of VMR values.

  \param sga_agenda     Agenda to use to calculate local absorption.
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
void abs_scalar_gas_fieldCalc(// WS Output:
                              Tensor5& asg_field,
                              Matrix&  asg,
                              Numeric& a_pressure,
                              Numeric& a_temperature,
                              Vector&  a_vmr_list,
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
            sga_agenda.execute(count);

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
            // asg_field: [ gas_species, f_grid, p_grid, lat_grid, lon_grid]
            // asg:       [ f_grid, gas_species ]
            asg_field( Range(joker),
                       Range(joker),
                       ipr, ila, ilo ) = transpose( asg );
            
            ++count;
          }
    }
}
