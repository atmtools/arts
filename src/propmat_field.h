/**
  * @file   propmat_field.h
  * @author Richard Larsson
  * @date   2019-02-26
  * 
  * @brief Implements a propagation matrix field
  *
  * This file contains internal code to speed up propagation field
  * calculations by sacrificing memory and by interpolating from
  * potentially coarser resolution.
*/

#ifndef PROPAGATION_FIELD_HEADER
#define PROPAGATION_FIELD_HEADER

#include "atm.h"
#include "field.h"
#include "workspace_ng.h"

class Agenda;
class Workspace;

using FieldOfMuelmatVector = Field3D<MuelmatVector>;
using FieldOfPropmatVector = Field3D<PropmatVector>;
using FieldOfStokvecVector = Field3D<StokvecVector>;

/** Creates a field of propagation matrices, absorption vectors, and source vectors
 * 
 * @param[in] ws A workspace
 * @param[out] propmat_field A 3D field of propagation matrices
 * @param[out] absorption_field A 3D field of absorption vectors
 * @param[out] additional_source_field A 3D field of source vectors
 * @param[in] f_grid As WSV
 * @param[in] p_grid As WSV
 * @param[in] z_field As WSV
 * @param[in] t_field As WSV
 * @param[in] nlte_field As WSV
 * @param[in] vmr_field As WSV
 * @param[in] jacobian_quantities As WSV
 * @param[in] propmat_clearsky_agenda As WSA
 */
void field_of_propagation(Workspace& ws,
                          FieldOfPropmatVector& propmat_field,
                          FieldOfStokvecVector& absorption_field,
                          FieldOfStokvecVector& additional_source_field,
                          const Vector& f_grid,
                          const AtmField& atm_field,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const Agenda& propmat_clearsky_agenda);

/** Get a field of transmission matrices from the propagation matrix field
 * 
 * @param[in] propmat_field A 3D field of propagation matrices
 * @param[in] r The distance
 * @return FieldOfMuelmatVector 
 */
FieldOfMuelmatVector transmat_field_calc_from_propmat_field(
    const FieldOfPropmatVector& propmat_field, const Numeric& r = 1.0);

/** Computes the radiation and transmission from fields of atmospheric propagation
 * 
 * Only for 1D atmospheres for now.  Only works when the propagation matrix
 * is not polarizing
 * 
 * Computes The forward simulations by interpolating the fields of
 * radiative properties to the selected propagation path.
 * 
 * Not well-tested.
 * 
 * @param[in] ws A workspace
 * @param[out] lvl_rad Level by level radiation
 * @param[out] src_rad Level by level source function
 * @param[out] lyr_tra Layered transmission
 * @param[out] tot_tra Total transmission from layer to background
 * @param[in] propmat_field 3D field of propagation matrices
 * @param[in] absorption_field A 3D field of absorption vectors
 * @param[in] additional_source_field A 3D field of source vectors
 * @param[in] f_grid As WSV
 * @param[in] t_field As WSV
 * @param[in] nlte_field As WSV
 * @param[in] ppath As WSV
 * @param[in] iy_main_agenda As WSA
 * @param[in] iy_space_agenda As WSA
 * @param[in] iy_surface_agenda As WSA
 * @param[in] iy_cloudbox_agenda As WSA
 * @param[in] surface_props_data As WSV
 */
void emission_from_propmat_field(
    Workspace& ws,
    ArrayOfStokvecVector& lvl_rad,
    ArrayOfStokvecVector& src_rad,
    ArrayOfMuelmatVector& lyr_tra,
    ArrayOfMuelmatVector& tot_tra,
    const FieldOfPropmatVector& propmat_field,
    const FieldOfStokvecVector& absorption_field,
    const FieldOfStokvecVector& additional_source_field,
    const Vector& f_grid,
    const AtmField& atm_field,
    const Ppath& ppath,
    const Agenda& iy_main_agenda,
    const Agenda& iy_space_agenda,
    const Agenda& iy_surface_agenda,
    const Agenda& iy_cloudbox_agenda,
    const SurfaceField& surface_field);

#endif  // PROPAGATION_FIELD_HEADER
