/*!
  \file   m_abs_lookup.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
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
                      // WS Input:
                      const Agenda& abs_coef_per_species_agenda,
//                      const Index& atmosphere_dim,
                      const ArrayOfArrayOfSpeciesTag& abs_species,
                      const Vector& f_grid,
                      const Vector& p_grid,
                      const Tensor4& vmr_field,
                      const Tensor3& t_field,
                      const Vector& t_pert,
                      const ArrayOfIndex& nls,
                      const Vector& nls_pert )
{
  // Absorption coefficients per species (output of absorption agenda):
  ArrayOfMatrix acps;

  // Absorption vmrs and temperature (input to absorption agenda):
  Matrix abs_vmrs;
  Vector abs_t;

  // Calling parameter for functions requireing lat/lon grids:
  const Vector empty_grid;      

  // FIXME: Which input parameters do we need?
  
  // FIXME: Checks of input parameter correctnes
  // Check Atmospheric fields:
  chk_atm_field("vmr_field",vmr_field,
                1,abs_species.nelem(),p_grid,
                empty_grid, empty_grid);

  chk_atm_field("t_field",t_field,
                1,p_grid,
                empty_grid, empty_grid);


  // FIXME: Check for simple case of no temperature variation?

  // FIXME: Perhaps go back to having t_pert, nls, nls_pert as WSVs,
  // so that we can pre-generate them by a specialized "to-match"
  // method? 
  // This would also allow more freedom to set up complicated
  // temperature perturbation vectors. (Explicit list could 
  // be a bit long and impractical.)

  // Set general lookup table properties:
  gal.species = abs_species;    // Species list
  gal.nonlinear_species = nls;  // Nonlinear species (H2O)
  gal.f_grid = f_grid;           // Frequency grid
  gal.p_grid = p_grid;          // Pressure grid
  gal.vmrs_ref = vmr_field(joker,joker,0,0);
  gal.t_ref = t_field(joker,0,0);
  gal.t_pert = t_pert;
  gal.nls_pert = nls_pert;

  // Now we have to fill gal.xsec with the right values!

  // Loop temperature perturbations:
  for ( Index i=0; i<gal.t_pert.nelem(); ++i )
    {
      // Create perturbed temperature profile:
      abs_t = gal.t_ref;
      abs_t += gal.t_pert[i];
      
      // FIXME: To do all species at once like this does not work for
      // nonlinear species. Three solutions:
      // a) Set up a longer species list with explicit nonlinear
      // species beforehand (plus perturbed VMR profiles). Not
      // difficult, but could blow memory?
      // b) Always loop over species, allowing consistent treatment of
      // nonlinear ones.
      // c) Do all regular species at once, then separate treatment of
      // nonlinear species. (Possibly also all nonlinear ones in one go.)
      
      // Call the agenda to calculate absorption coefficients:
      abs_coef_per_species_agendaExecute(
                                         // Output
                                         acps,
                                         // Input
                                         abs_vmrs,
                                         abs_t,
                                         // Wrapper Input
                                         abs_coef_per_species_agenda,
                                         false);
    }
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
