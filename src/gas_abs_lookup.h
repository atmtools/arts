/*!
  \file   gas_abs_lookup.h
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Sep 19 16:49:00 2002
  
  \brief  Declarations for the gas absorption lookup table.  
*/

#ifndef gas_abs_lookup_h
#define gas_abs_lookup_h

#include "matpackIV.h"
#include "absorption.h"

//! An absorption lookup table.
/*! This class holds an absorption lookup table, as well as all
    information that is necessary to use the table to extract
    absorption. Extraction routines are implemented as member functions. */
struct GasAbsLookup {
public:
  // Documentation is with the implementation!
  void Adapt();

  // Documentation is with the implementation!
  void Extract( Matrix sga,
		Numeric p,
		Numeric T);
private:

  //! This contains the tag groups for which the table is valid:
  TagGroups gas_tgs; 

  //! The frequency grid [Hz].
  /*! Must be sorted in ascending order. */
  Vector    f_grid;

  //! The pressure grid for the table [Pa].
  /*! Must be sorted in [FIXME: ascending or descending?] order. */
  Vector    p_grid;  

  //! The base 10 logarithm of the pressure grid.
  /*! This is not stored along with the table, but calculated when the
    table is initialized with Adapt. The reason to have this is that
    vertical interpolation should be linear in log(p). */
  Vector   log_p_grid;  

  //! The VMR profiles.
  /*! The VMRs for all species, associated with p_grid. Dimension:
    [N_gas_tgs, N_p_grid]. These VMRs are needed to scale the
    absorption coefficient to other VMRs. We are never working with
    "absorption cross-sections", always with real absorption coefficients,
    so we have to remember the associated VMR values. 

    Physical unit: Absolute value. */
  Matrix    vmrs;

  //! The reference temperature profile [K].
  /*! This is a temperature profile. The dimension must be the same as
    p_grid. */
  Vector    t_ref;

  //! The vector of temperature perturbations [K].
  /*! This can have any number of elements. Example:
    [-20,-10,0,10,20]. The actual temperatures for which absorption is
    stored are t_ref + t_pert for each level. The reference
    temperature itself should normally also be included, hence t_pert should
    always include 0. Must be sorted in ascending order!

    The vector t_pert may be an empty vector (nelem()=0), which marks
    the special case that no interpolation in temperature should be
    done. If t_pert is not empty, you will get an error message if you
    try to extract absorption for temperatures outside the range of
    t_pert. */
  Vector    t_pert;

  //! The vector of H2O VMR perturbations.
  /*! This can hold perturbations for the H2O profile, in analogy to
    t_pert. Fractional units are used! Example: [0,.5,1,10,100],
    meaning from H2O VMR 0 to 100 times the H2O profile given in
    vmrs. The reference value should normally be included, hence
    h2o_pert should always include the value 1.

    The vector h2o_pert may be an empty vector (nelem()=0), in which
    case H2O is treated like all other gases. As in the case of
    temperature, you cannot extract absorption for H2O VMR values
    outside h2o_pert, if h2o_pert is not empty.

    Per definition, h2o_pert refers to the first H2O tag in gas_tgs,
    should there be more than one. If there is no H2O tag in gas_tgs,
    h2o_pert must be an empty vector.*/
  Vector    h2o_pert;

  //! Absorption coefficients.
  /*! Physical unit: 1/m
    Dimension: [ N_gas_tgs (N_gas_tgs+N_h2o_pert if special treatment for H2O),
    N_t_pert (or 1 if N_t_pert=0),
    N_p_grid,
    N_f_grid] */
  Tensor4 abs;

};

#endif //  gas_abs_lookup_h
