/*!
  \file   gas_abs_lookup.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Sep 19 17:25:07 2002
  
  \brief  Implementation of scalar gas absorption lookup table functions.
*/

#include <cmath>
#include "gas_abs_lookup.h"
#include "interpolation.h"
#include "make_vector.h"
#include "logic.h"

//! Adapt lookup table to current calculation.
/*! Verify that the lookup table contains all species used in the
  current calculation. Remove extra species which are not needed from
  the table and re-sort the table to be consistent with the actual
  list of species. */
void GasAbsLookup::Adapt()
{

}

//! Extract scalar gas absorption coefficients from the lookup table. 
/*! 
  This is the extraction routine intended for use inside the cloud
  box. It extracts absorption for only one frequency, and for only one
  atmospheric condition (one p/T combination).

  It carries out a simple interpolation in temperature and
  pressure. The interpolated value is then scaled by the ratio between
  actual VMR and reference VMR. In the case of nonlinear species the
  interpolation goes also over VMR.

  All input parameters (f_index, p, T) must be in the range coverd by
  the table.

  \retval sga A Matrix with scalar gas absorption coefficients
  [1/m]. Dimension: [n_frequencies, n_species]. Must already have the
  right dimension before function call!

  \param f_index The frequency index. If this is >=0, it means that
  absorption for this frequency will be extracted. (The leading
  dimension of sga must be 1.) If this is <0, it means that absorption
  for ALL frequencies is extracted. (The leading dimension fo sga must
  match the dimension of the tables frequency grid.)

  \param p The pressures [Pa].

  \param T The temperature [K].

  \param vmrs The VMRs [absolute number]. Dimension: [n_species]. 
*/
void GasAbsLookup::Extract( MatrixView sga,
                            Index      f_index,
                            Numeric    p,
                            Numeric    T,
                            VectorView vmrs )
{

  // Obtain some properties of the lookup table:
  
  // Number of gas species in the table:
  const Index n_species = species.nelem();

  // Number of nonlinear species:
  const Index n_nonlinear_species = nonlinear_species.nelem();

  // Number of frequencies in the table:
  const Index n_f_grid = f_grid.nelem();

  // Number of pressure grid points in the table:
  const Index n_p_grid = p_grid.nelem();

  // Number of temperature perturbations:
  const Index n_t_pert = t_pert.nelem();

  // Number of nonlinear species perturbations:
  const Index n_nls_pert = nls_pert.nelem();


  // First some checks on the lookup table itself:

  // Assert that log_p_grid has been initialized:
  assert( is_size( log_p_grid, n_p_grid) );

  // Check that the dimension of vmrs_ref is consistent with species and p_grid:
  assert( is_size( vmrs_ref, n_species, n_p_grid) );

  // Check dimension of t_ref:
  assert( is_size( t_ref, n_p_grid) );


  // Following are some checks on the input variables:

  // We also set the start and extent for the frequency loop.
  Index f_start, f_extent;

  if ( f_index < 0 )
    {
      // This means we should extract for all frequencies.

      // Dimension of output variable sga must match this:
      assert( is_size( sga, n_f_grid, n_species ) );

      f_start  = 0;
      f_extent = n_f_grid;
    }
  else
    {
      // This means we should extract only for one frequency.

      // Dimension of output variable sga must match this:
      assert( is_size( sga, 1, n_species ) );

      // Check that f_index is inside f_grid:
      assert( f_index < n_f_grid );

      f_start  = f_index;
      f_extent = 1;
    }

  // Assert that vmrs has the right dimension:
  assert( is_size( vmrs, n_species ) );

  // We need the log10 of the pressure:
  const Numeric log_p = log10( p );


  // Now we will start to do some business.

  // For sure, we need to store the vertical grid positions:
  // We need only one, because we want to extract only for a single
  // atmospheric condition.
  ArrayOfGridPos vgp(1);
  gridpos( vgp,
	   log_p_grid,
	   log_p );

  // Let's first treat the case with no nonlinear species.
  // This case is marked by n_nonlinear_species == 0.
  if ( 0 == n_nonlinear_species )
    {
      // In this case, n_nls_pert must also be zero:
      assert( is_size( nls_pert, 0 ) );
      
      // The simplest case is that the table contains neither temperature
      // nor VMR perturbations. This means it is not really a
      // lookup table, just an absorption profile. In this case we ignore
      // the temperature, and interpolate in pressure only.
      // This case is marked by t_pert being an empty vector, since we
      // already know that there are no non-linear species.
      if ( 0 == n_t_pert )
	{
	  // Verify, that abs has the right dimensions for this case:
          
	  assert( is_size( abs,
                           1,          // temperature perturbations
                           n_species,  // species  
                           n_f_grid,   // frequency grid points
                           n_p_grid    // pressure levels      
                           ) );
	  
	  // To store interpolation weights:
	  Matrix itw(1,2);
	  interpweights(itw,vgp);

          // Loop over frequency:
	  for ( Index s=f_start; s<f_extent; ++s )
	    {
              // Loop over species:
              for ( Index i=0; i<n_species; ++i )
                {
		  // Get the right view on abs. (Only a vector of
		  // pressure for this particular species and
		  // frequency):
		  ConstVectorView this_abs = abs( 0,
						  i,
						  s, 
						  Range(joker) );

		  // Get the right view on our result variable,
		  // sga. (Only a scalar):
		  Numeric& this_sga = sga(s,i);

		  // Do the interpolation:
		  interp( this_sga,
			  itw,
			  this_abs,
			  vgp );
                }
            }
	}
    }
}

