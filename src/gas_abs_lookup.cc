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

  \retval sga A Vector with scalar gas absorption coefficients
  [1/m]. Dimension: [n_species]. Must already have the right dimension
  before function call!

  \param f_index The frequency index.

  \param p The pressures [Pa].

  \param T The temperature [K].

  \param vmrs The VMRs [absolute number]. Dimension: [n_species]. 
*/
void GasAbsLookup::Extract( VectorView sga,
                            Index      f_index,
                            Numeric    p,
                            Numeric    T,
                            VectorView vmrs )
{
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

  // The variable sga is the local scalar absorption coefficient in
  // [1/m]. It is stored for each species.  We have to compute this
  // here, as a function of the given input variables pressure p and
  // temperature T.  We assume that the user has already given sga the
  // right dimension. Check this:
  assert( is_size( sga, n_species ) );

  // Check that f_index is inside f_grid:
  assert( 0 <= f_index );
  assert( f_index < n_f_grid );

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

  // Let's first treat the case with no water vapor perturbations.
  // This case is marked by h2o_pert being an empty vector.
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

	  // Loop over species:
	  for ( Index i=0; i<n_species; ++i )
	    {
		  // Get the right view on abs. (Only a vector of
		  // pressure for this particular species and
		  // frequency):
		  ConstVectorView this_abs = abs( 0,
						  i,
						  f_index, 
						  Range(joker) );

		  // Get the right view on our result variable,
		  // sga. (Only a scalar):
		  Numeric& this_sga = sga[i];

		  // Do the interpolation:
		  interp( this_sga,
			  itw,
			  this_abs,
			  vgp );
	    }
	}
    }
}

//! Extract scalar gas absorption coefficients from the lookup table. 
/*! This means a simple interpolation in pressure and
  temperature. The interpolated value is then scaled by the ratio
  between actual VMR and reference VMR. In the case of H2O the
  interpolation goes also over VMR.
    
  This function can do the extraction for several atmospheric
  conditions at the same time. Therefore, the first dimension of the
  output sga and the input vmr is n_p, the dimension of p and T.

  Because of the adaptability of Matpack, you can call this method
  also like this:

  GasAbsLookup::Extract( MatrixView sga,
                         Numeric    p,
			 Numeric    T)

  in which case only one atmospheric condition is done at a time. The
  dimension of sga is then [n_species, n_f_grid].

  \retval sga A Tensor3 with scalar gas absorption
  coefficients [1/m]. Dimension: [n_p, n_species, n_f_grid]. Must already
  have the right dimension before function call! The first dimension
  corresponds to the number of points, for which absorption
  coefficients are wanted at the same time.

  \param p The vector of pressures [Pa].
  \param T The vector of temperatures [K].
  \param vmrs The VMRs [absolute number]. This has to be a matrix of
  dimension [n_p,n_species]. 
*/
void GasAbsLookup::Extract( Tensor3View     sga,
			    ConstVectorView p,
			    ConstVectorView T)
{
  // Number of points for which we want to extract absorption:
  const Index n_p = p.nelem();

  // Number of gas species in the table:
  const Index n_species = species.nelem();

  // Number of nonlinear species:
  const Index n_nonlinear_species = nonlinear_species.nelem();

  // Number of frequencies in the table:
  const Index n_f_grid = f_grid.nelem();

  // Number of pressure grid points in the table:
  const Index n_p_grid = p_grid.nelem();

  // Number of H2O per

  // Note the important difference between n_p, the number of pressures
  // for which we want to extract absorption, and n_p_grid, the number
  // of pressures in the lookup table!

  // Following are some checks on the input variables:

  // The variable sga is the local scalar absorption coefficient in
  // [1/m]. It is stored for each atmospheric condition, each gas
  // species, and as a function of frequency. Hence the dimension is:
  // [n_species, n_f_grid].  We have to compute this here, as a
  // function of the given input variables pressure p and temperature
  // T.  We assume that the user has already given sga the right
  // dimension. Check this:
  assert( sga.npages() == n_p );
  assert( sga.nrows()  == n_species );
  assert( sga.ncols()  == n_f_grid  );

  // Verify that vectors p and T have the same number of elements:
  assert( n_p == T.nelem() );

  // Assert that vmrs_ref has the right dimension:
  assert( n_p        == vmrs_ref.nrows() );
  assert( n_species  == vmrs_ref.ncols() );

  // We need the log10 of the pressure:
  Vector log_p( n_p );
  transform( log_p, log10, p );

  // Now some checks on the lookup table:

  // Assert that log_p_grid has been initialized:
  assert( n_p_grid == log_p_grid.nelem() );

  // ...and that t_ref has the right dimension:
  assert( n_p_grid  == t_ref.nelem() );

  // Now we will start to do some business.

  // For sure, we need to store the vertical grid positions:
  ArrayOfGridPos vgp(n_p);
  gridpos( vgp,
	   log_p_grid,
	   log_p );

  // Let's first treat the case with no water vapor perturbations.
  // This case is marked by h2o_pert being an empty vector.
  if ( 0 == n_nonlinear_species )
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
	  assert( abs.nbooks() == n_species );
	  assert( abs.npages() == 1               );
	  assert( abs.nrows()  == n_p_grid  );
	  assert( abs.ncols()  == n_f_grid  );
	  
	  // To store interpolation weights:
	  Matrix itw(n_p,2);
	  interpweights(itw,vgp);

	  // Loop over species:
	  for ( Index i=0; i<n_species; ++i )
	    {
	      // Loop over frequency:
	      for ( Index j=0; j<n_f_grid; ++j )
		{
		  // Get the right view on abs. (Only a vector of
		  // pressure for this particular species and
		  // frequency):
		  ConstVectorView this_abs = abs( i,
						  0,
						  Range(joker),
						  j );

		  // Get the right view on our result variable,
		  // sga. (Only a Vector with dimension n_p, because we
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
