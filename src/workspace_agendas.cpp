#include "workspace_agendas.h"

std::unordered_map<std::string, WorkspaceAgendaInternalRecord>
internal_workspace_agendas() {
  std::unordered_map<std::string, WorkspaceAgendaInternalRecord> wsa_data;

  wsa_data["propmat_clearsky_agenda"] = {
      .desc = R"--(Calculate the absorption coefficient matrix.

This agenda calculates the absorption coefficient matrix for all
absorption species as a function of the given atmospheric state for
one point in the atmosphere. The result is returned in
*propmat_clearsky*.
)--",
      .output = {"propmat_clearsky",
                 "nlte_source",
                 "dpropmat_clearsky_dx",
                 "dnlte_source_dx"},
      .input = {"jacobian_targets",
                "select_abs_species",
                "f_grid",
                "path_point",
                "atm_point"}};

  wsa_data["spectral_radiance_background_space_agenda"] = {
      .desc = R"--(Spectral radiance as seen of space.

This agenda calculates the spectral radiance as seen of space. The
intent is to provide a background spectral radiance from space that
is input to the atmospheric radiative transfer calculations.

The input path point should be as if it is looking at space.

The output must be sized as:

- *spectral_radiance_background* : (*f_grid*)
- *spectral_radiance_background_jacobian* : (*jacobian_targets*, *f_grid*)
)--",
      .output = {"spectral_radiance_background",
                 "spectral_radiance_background_jacobian"},
      .input = {"f_grid", "jacobian_targets", "path_point"}};

  wsa_data["spectral_radiance_background_surface_agenda"] = {
      .desc = R"--(Spectral radiance as seen of the surface.

This agenda calculates the spectral radiance as seen of the surface.
The intent is to provide a background spectral radiance from the
surface that is input to the atmospheric radiative transfer
calculations.

The input path point should be as if it is looking at the surface.

The output must be sized as:

- *spectral_radiance_background* : (*f_grid*)
- *spectral_radiance_background_jacobian* : (*jacobian_targets*, *f_grid*)
)--",
      .output = {"spectral_radiance_background",
                 "spectral_radiance_background_jacobian"},
      .input = {"f_grid", "jacobian_targets", "path_point"}};

  return wsa_data;
}
