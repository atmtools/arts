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
    
  \retval sga A matrix with scalar gas absorption
  coefficients [1/m]. Dimension: [N_gas_tgs, N_f_grid]. Must already
  have the right dimension before function call!

  \param p The pressure for which absorption is wanted.
  \param T The temperature for which absorption is wanted.
    
*/
void GasAbsLookup::Extract( Matrix sga,
			    Numeric p,
			    Numeric T)
{
  // The variable sga is the local scalar absorption coefficient in
  // [1/m]. It is stored for each gas species and as a function of
  // frequency. Hence the dimension is: [N_gas_tgs, N_f_grid].
  // We have to compute this here, as a function of the given input
  // variables pressure p and temperature T.
  // We assume that the user has already given sga the right
  // dimension. Check this:
  assert( sga.nrows() == gas_tgs.nelem() );
  assert( sga.ncols() == f_grid.nelem()  );

  // Assert that log_p_grid has been initialized:
  assert( p_grid.nelem() == log_p_grid.nelem() );

  // Assert that vmrs and t_ref have the right dimension:
  assert( gas_tgs.nelem() == vmrs.nrows() );
  assert( p_grid.nelem()  == vmrs.ncols() );
  assert( p_grid.nelem()  == t_ref.nelem() );

  // For sure, we need to store the vertical grid position:
  ArrayOfGridPos vgp(1);
  gridpos( vgp,
	   log_p_grid,
	   MakeVector(log(p)));

  // Let's first treat the case with no water vapor perturbations:
  if ( 0 == h2o_pert.nelem() )
    {
      // The simplest case is that the table contains neither temperature
      // nor water vapor perturbations. This means it is not really a
      // lookup table, just an absorption profile. In this case we ignore
      // the temperature, and interpolate in pressure only.
      if ( 0 == t_pert.nelem() )
	{
	  // FIXME: Continue here...
	}

    }
}
