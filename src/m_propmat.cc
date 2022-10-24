//! File contains ways to manipulate the propagation matrix

#include "agenda_class.h"
#include "auto_md.h"
#include "debug.h"
#include "energylevelmap.h"
#include "jacobian.h"
#include "propagationmatrix.h"
#include "species_tags.h"

void propmat_clearskyAddScaledSpecies(  // Workspace reference:
    Workspace& ws,
    // WS Output:
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    // WS Input:
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfSpeciesTag& select_abs_species,
    const Vector& f_grid,
    const Vector& rtp_mag,
    const Vector& rtp_los,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const EnergyLevelMap& rtp_nlte,
    const Vector& rtp_vmr,
    const Agenda& propmat_clearsky_agenda,
    // WS Generic Input:
    const ArrayOfSpeciesTag& target,
    const Numeric& scale,
    // Verbosity object:
    const Verbosity&) {
  ARTS_USER_ERROR_IF(jacobian_quantities.nelem(), "Cannot use with derivatives")

  if (select_abs_species not_eq target) {
    ARTS_USER_ERROR_IF(
        select_abs_species.nelem(),
        "Non-empty select_abs_species (lookup table calculations set select_abs_species)")

    PropagationMatrix pm;
    StokesVector sv;
    ArrayOfPropagationMatrix dpropmat_clearsky_dx;
    ArrayOfStokesVector dnlte_source_dx;

    propmat_clearsky_agendaExecute(ws,
                                   pm,
                                   sv,
                                   dpropmat_clearsky_dx,
                                   dnlte_source_dx,
                                   jacobian_quantities,
                                   target,
                                   f_grid,
                                   rtp_mag,
                                   rtp_los,
                                   rtp_pressure,
                                   rtp_temperature,
                                   rtp_nlte,
                                   rtp_vmr,
                                   propmat_clearsky_agenda);

    ARTS_USER_ERROR_IF(propmat_clearsky.Data().shape() not_eq pm.Data().shape(), "Mismatching sizes")
    ARTS_USER_ERROR_IF(nlte_source.Data().shape() not_eq sv.Data().shape(), "Mismatching sizes")
    
    propmat_clearsky += scale * pm;
    nlte_source += scale * sv;
  }
}