/*!
  \file   m_absorption.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Wed Nov 20 18:04:20 2002
  
  \brief  Methods related to absorption, lookup table, etc.
*/

#include "arts.h"
#include "messages.h"
#include "gas_abs_lookup.h"

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
