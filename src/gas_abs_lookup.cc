/*!
  \file   gas_abs_lookup.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Sep 19 17:25:07 2002
  
  \brief  Implementation of scalar gas absorption lookup table functions.
*/

#include "gas_abs_lookup.h"
#include "interpolation.h"
#include "make_vector.h"
#include <cmath>

//! Adapt lookup table to current calculation.
/*! Verify that the lookup table contains all species used in the
  current calculation. Remove extra species which are not needed from
  the table and re-sort the table to be consistent with the actual
  list of species. */
void GasAbsLookup::Adapt()
{

}

//! Extract scalar gas absorption coefficients from the lookup table. 
/*! This means a simple interpolation in pressure and
  temperature. The interpolated value is then scaled by the ratio
  between actual VMR and reference VMR. In the case of H2O the
  interpolation goes also over VMR.
    
  \retval sga A Tensor3 with scalar gas absorption
  coefficients [1/m]. Dimension: [Np, N_gas_tgs, N_f_grid]. Must already
  have the right dimension before function call! The first dimension
  corresponds to the number of points, for which absorption
  coefficients are wanted at the same time.

  \param p The vector of pressures for which absorption is wanted.
  \param T The vector of temperatures for which absorption is wanted.
    
*/
void GasAbsLookup::Extract( Tensor3 sga,
			    ConstVectorView& p,
			    ConstVectorView& T)
{
  // Number of points for which we want to extract absorption:
  const Index Np = p.nelem();

  // Number of gas species in the table:
  const Index Nspecies = gas_tgs.nelem();

  // Number of frequencies in the table:
  const Index Nf = f_grid.nelem();

  // The variable sga is the local scalar absorption coefficient in
  // [1/m]. It is stored for each gas species and as a function of
  // frequency. Hence the dimension is: [N_gas_tgs, N_f_grid].
  // We have to compute this here, as a function of the given input
  // variables pressure p and temperature T.
  // We assume that the user has already given sga the right
  // dimension. Check this:
  assert( sga.npages() == Np );
  assert( sga.nrows()  == Nspecies );
  assert( sga.ncols()  == Nf  );

  // Verify that vectors p and T have the same number of elements:
  assert( Np == T.nelem() );

  // We need the log10 of the pressure:
  Vector log_p( Np );
  transform( log_p, log10, p );

  // Assert that log_p_grid has been initialized:
  assert( p_grid.nelem() == log_p_grid.nelem() );

  // Assert that vmrs and t_ref have the right dimension:
  assert( Nspecies == vmrs.nrows() );
  assert( p_grid.nelem()  == vmrs.ncols() );
  assert( p_grid.nelem()  == t_ref.nelem() );

  // For sure, we need to store the vertical grid positions:
  ArrayOfGridPos vgp(Np);
  gridpos( vgp,
	   log_p_grid,
	   log_p );

  // Let's first treat the case with no water vapor perturbations.
  // This case is marked by h2o_pert being an empty vector.
  if ( 0 == h2o_pert.nelem() )
    {
      // The simplest case is that the table contains neither temperature
      // nor water vapor perturbations. This means it is not really a
      // lookup table, just an absorption profile. In this case we ignore
      // the temperature, and interpolate in pressure only.
      // This case is marked by t_pert being an empty vector, since we
      // already know that h2o_pert is also empty.
      if ( 0 == t_pert.nelem() )
	{
	  // Verify, that abs has the right dimensions for this case:
	  // 1. Dimension: Number of gas tags (no H2O perturbations)
	  // 2. Dimension: 1 (no T perturbations)
	  // 3. Dimension: Number of pressure levels
	  // 4. Dimension: Number of frequencies
	  assert( abs.nbooks() == Nspecies );
	  assert( abs.npages() == 1               );
	  assert( abs.nrows()  == p_grid.nelem()  );
	  assert( abs.ncols()  == Nf  );
	  
	  // To store interpolation weights:
	  Matrix itw(Np,2);
	  interpweights(itw,vgp);

	  // Loop over species:
	  for ( Index i=0; i<Nspecies; ++i )
	    {
	      // Loop over frequency:
	      for ( Index j=0; j<Nf; ++j )
		{
		  // Get the right view on abs. (Only a vector of
		  // pressure for this particular species and
		  // frequency):
		  ConstVectorView this_abs = abs( i,
						  0,
						  Range(joker),
						  j );

		  // Get the right view on our result variable,
		  // sga. (Only a Vector with dimension Np, because we
		  // selected species and frequency):
		  VectorView this_sga = sga( Range(joker),
					     i,
					     j );

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
