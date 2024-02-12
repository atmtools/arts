#include "workspace_agendas.h"

std::unordered_map<std::string, WorkspaceAgendaInternalRecord>
internal_workspace_agendas() {
  std::unordered_map<std::string, WorkspaceAgendaInternalRecord> wsa_data;

  wsa_data["propagation_matrix_agenda"] = {
      .desc = R"--(Compute the propagation matrix, the non-LTE source vector, and their derivatives
)--",
      .output = {"propagation_matrix",
                 "propagation_matrix_source_vector_nonlte",
                 "propagation_matrix_jacobian",
                 "propagation_matrix_source_vector_nonlte_jacobian"},
      .input = {"jacobian_targets",
                "propagation_matrix_select_species",
                "frequency_grid",
                "propagation_path_point",
                "atmospheric_point"},
  };

  wsa_data["propagation_path_observer_agenda"] = {
      .desc = R"--(Get the propagation path as it is obeserved.

The intent of this agenda is to provide a propagation path as seen from the observer
position and line of sight.
)--",
      .output = {"propagation_path"},
      .input = {"spectral_radiance_observer_position",
                "spectral_radiance_observer_line_of_sight"},
  };

  wsa_data["spectral_radiance_observer_agenda"] = {
      .desc =
          R"--(Spectral radiance as seen from the input position and environment

The intent of this agenda is to provide a spectral radiance as seen from the observer
position and line of sight.

The output must be sized as:

- *spectral_radiance* : (*frequency_grid*)
- *spectral_radiance_jacobian* : (*jacobian_targets*, *frequency_grid*)
)--",
      .output = {"spectral_radiance", "spectral_radiance_jacobian"},
      .input = {"frequency_grid",
                "jacobian_targets",
                "spectral_radiance_observer_position",
                "spectral_radiance_observer_line_of_sight",
                "atmospheric_field",
                "surface_field"},
  };

  wsa_data["spectral_radiance_space_agenda"] = {
      .desc = R"--(Spectral radiance as seen of space.

This agenda calculates the spectral radiance as seen of space.
One common use-case us to provide a background spectral radiance.

The input path point should be as if it is looking at space.

The output must be sized as:

- *spectral_radiance* : (*frequency_grid*)
- *spectral_radiance_jacobian* : (*jacobian_targets*, *frequency_grid*)
)--",
      .output = {"spectral_radiance", "spectral_radiance_jacobian"},
      .input = {"frequency_grid", "jacobian_targets", "propagation_path_point"},
  };

  wsa_data["spectral_radiance_surface_agenda"] = {
      .desc = R"--(Spectral radiance as seen of the surface.

This agenda calculates the spectral radiance as seen of the surface.
One common use-case us to provide a background spectral radiance.

The input path point should be as if it is looking at the surface.

The output must be sized as:

- *spectral_radiance* : (*frequency_grid*)
- *spectral_radiance_jacobian* : (*jacobian_targets*, *frequency_grid*)
)--",
      .output = {"spectral_radiance", "spectral_radiance_jacobian"},
      .input = {"frequency_grid",
                "jacobian_targets",
                "propagation_path_point",
                "surface_field"},
  };

  return wsa_data;
}
