#include "workspace_agendas.h"

#include <format>

using namespace std::literals;

std::unordered_map<std::string, WorkspaceAgendaInternalRecord>
internal_workspace_agendas_creator() {
  std::unordered_map<std::string, WorkspaceAgendaInternalRecord> wsa_data;

  wsa_data["propagation_matrix_agenda"] = {
      .desc =
          R"--(Compute the propagation matrix, the non-LTE source vector, and their derivatives.

The intent of this agenda is to be the workhorse for the propagation matrix
calculations that are happening deep in your ARTS method calls.

.. tip::
    Use *propagation_matrix_agendaAuto* after having defined
    your absorption data to create this agenda.  It covers most use-cases.
)--",
      .output       = {"propagation_matrix",
                       "propagation_matrix_source_vector_nonlte",
                       "propagation_matrix_jacobian",
                       "propagation_matrix_source_vector_nonlte_jacobian"},
      .input        = {"frequency_grid",
                       "frequency_grid_wind_shift_jacobian",
                       "jacobian_targets",
                       "select_species",
                       "ray_path_point",
                       "atmospheric_point"},
      .enum_options = {"Empty"},
      .output_constraints =
          {
              {"propagation_matrix.size() == frequency_grid.size()",
               "On output, *propagation_matrix* has the size of *frequency_grid*.",
               "propagation_matrix.size()",
               "frequency_grid.size()"},
              {"propagation_matrix_source_vector_nonlte.size() == frequency_grid.size()",
               "On output, *propagation_matrix_source_vector_nonlte* has the size of *frequency_grid*.",
               "propagation_matrix_source_vector_nonlte.size()",
               "frequency_grid.size()"},
              {"matpack::same_shape<2>({jacobian_targets.target_count(), frequency_grid.size()}, propagation_matrix_jacobian)",
               "On output, *propagation_matrix_jacobian* has the shape of the target-count of *jacobian_targets* times the size of *frequency_grid*.",
               "propagation_matrix_jacobian.shape()",
               "frequency_grid.size()",
               "jacobian_targets.target_count()"},
              {"matpack::same_shape<2>({jacobian_targets.target_count(), frequency_grid.size()}, propagation_matrix_source_vector_nonlte_jacobian)",
               "On output, *propagation_matrix_source_vector_nonlte_jacobian* has the shape of the target-count of *jacobian_targets* times the size of *frequency_grid*.",
               "propagation_matrix_source_vector_nonlte_jacobian.shape()",
               "frequency_grid.size()",
               "jacobian_targets.target_count()"},
          },
  };

  wsa_data["propagation_matrix_scattering_spectral_agenda"] = {
      .desc =
          R"--(Get the scattering propagation matrix, the scattering absorption vector, and the scattering spectral phase matrix.
)--",
      .output = {"propagation_matrix_scattering",
                 "absorption_vector_scattering",
                 "phase_matrix_scattering_spectral"},
      .input  = {"frequency_grid", "atmospheric_point", "legendre_degree"},
      .enum_options = {"FromSpeciesTRO"},
      .enum_default = "FromSpeciesTRO",
      .output_constraints =
          {
              {"propagation_matrix_scattering.size() == frequency_grid.size()",
               "On output, *propagation_matrix_scattering* has the size of *frequency_grid*.",
               "propagation_matrix_scattering.size()",
               "frequency_grid.size()"},
              {"absorption_vector_scattering.size() == frequency_grid.size()",
               "On output, *absorption_vector_scattering* has the size of *frequency_grid*.",
               "absorption_vector_scattering.size()",
               "frequency_grid.size()"},
              {"matpack::same_shape<2>({frequency_grid.size(), legendre_degree + 1}, phase_matrix_scattering_spectral)",
               "On output, *phase_matrix_scattering_spectral* has the shape of <*legendre_degree* + 1> times the size of *frequency_grid*.",
               "phase_matrix_scattering_spectral.shape()",
               "frequency_grid.size()",
               "legendre_degree"},
          },
  };

  wsa_data["propagation_matrix_scattering_agenda"] = {
      .desc =
          R"--(Compute the part of the propagation matrix that relates to scattering.
)--",
      .output       = {"propagation_matrix_scattering"},
      .input        = {"frequency_grid", "atmospheric_point"},
      .enum_options = {"AirSimple"},
      .enum_default = "AirSimple",
      .output_constraints =
          {
              {"propagation_matrix_scattering.size() == frequency_grid.size()",
               "On output, *propagation_matrix_scattering* has the size of *frequency_grid*.",
               "propagation_matrix_scattering.size()",
               "frequency_grid.size()"},
          },
  };

  wsa_data["ray_path_observer_agenda"] = {
      .desc   = R"--(Get the propagation path as it is obeserved.

The intent of this agenda is to provide a propagation path as seen from the observer
position and line of sight.

.. tip::
    The perhaps easiest way to set this agenda up is to use the *ray_path_observer_agendaSetGeometric* method.
)--",
      .output = {"ray_path"},
      .input  = {"spectral_radiance_observer_position",
                 "spectral_radiance_observer_line_of_sight"},
  };

  wsa_data["spectral_radiance_observer_agenda"] = {
      .desc =
          R"--(Spectral radiance as seen from the input position and environment.

The intent of this agenda is to provide the spectral radiance as seen from the observer
position and line of sight.

It also outputs the *ray_path* as seen from the observer position and line of sight.
This is useful in-case a call to the destructive *spectral_radianceApplyUnitFromSpectralRadiance*
is warranted
)--",
      .output = {"spectral_radiance", "spectral_radiance_jacobian", "ray_path"},
      .input  = {"frequency_grid",
                 "jacobian_targets",
                 "spectral_radiance_observer_position",
                 "spectral_radiance_observer_line_of_sight",
                 "atmospheric_field",
                 "surface_field",
                 "subsurface_field"},
      .enum_options = {"Emission", "EmissionNoSensor"},
      .enum_default = "Emission",
      .output_constraints =
          {
              {"spectral_radiance.size() == frequency_grid.size()",
               "On output, *spectral_radiance* has the size of *frequency_grid*.",
               "spectral_radiance.size()",
               "frequency_grid.size()"},
              {"matpack::same_shape<2>({jacobian_targets.x_size(), frequency_grid.size()}, spectral_radiance_jacobian)",
               "On output, *spectral_radiance_jacobian* has the shape of the expected *model_state_vector* (i.e., the x-size of *jacobian_targets*) times the size of *frequency_grid*.",
               "spectral_radiance_jacobian.shape()",
               "frequency_grid.size()",
               "jacobian_targets.x_size()"},
          },
  };

  wsa_data["spectral_radiance_space_agenda"] = {
      .desc         = R"--(Spectral radiance as seen of space.

This agenda calculates the spectral radiance as seen of space.
One common use-case is to provide a background spectral radiance.

The input path point should be as if it is looking at space.
)--",
      .output       = {"spectral_radiance", "spectral_radiance_jacobian"},
      .input        = {"frequency_grid", "jacobian_targets", "ray_path_point"},
      .enum_options = {"UniformCosmicBackground",
                       "SunOrCosmicBackground",
                       "Transmission"},
      .enum_default = "UniformCosmicBackground",
      .output_constraints = {
          {"spectral_radiance.size() == frequency_grid.size()",
           "On output, *spectral_radiance* has the size of *frequency_grid*.",
           "spectral_radiance.size()",
           "frequency_grid.size()"},
          {"matpack::same_shape<2>({jacobian_targets.x_size(), frequency_grid.size()}, spectral_radiance_jacobian)",
           "On output, *spectral_radiance_jacobian* has the shape of the expected *model_state_vector* (i.e., the x-size of *jacobian_targets*) times the size of *frequency_grid*.",
           "spectral_radiance_jacobian.shape()",
           "frequency_grid.size()",
           "jacobian_targets.x_size()"},
      }};

  wsa_data["spectral_radiance_surface_agenda"] = {
      .desc               = R"--(Spectral radiance as seen of the surface.

This agenda calculates the spectral radiance as seen of the surface.
One common use-case us to provide a background spectral radiance.

The input path point should be as if it is looking at the surface.

Subsurface calculations are also supported through this agenda.
)--",
      .output             = {"spectral_radiance", "spectral_radiance_jacobian"},
      .input              = {"frequency_grid",
                             "jacobian_targets",
                             "ray_path_point",
                             "surface_field",
                             "subsurface_field"},
      .enum_options       = {"Blackbody",
                             "Transmission",
                             "FlatScalarReflectance",
                             "FlatFresnelReflectance"},
      .enum_default       = "Blackbody",
      .output_constraints = {
          {"spectral_radiance.size() == frequency_grid.size()",
           "On output, *spectral_radiance* has the size of *frequency_grid*.",
           "spectral_radiance.size()",
           "frequency_grid.size()"},
          {"matpack::same_shape<2>({jacobian_targets.x_size(), frequency_grid.size()}, spectral_radiance_jacobian)",
           "On output, *spectral_radiance_jacobian* has the shape of the expected *model_state_vector* (i.e., the x-size of *jacobian_targets*) times the size of *frequency_grid*.",
           "spectral_radiance_jacobian.shape()",
           "frequency_grid.size()",
           "jacobian_targets.x_size()"},
      }};

  wsa_data["inversion_iterate_agenda"] = {
      .desc               = R"--(Work in progress ...

See *OEM*.

.. note::
    The output *measurement_jacobian* size may depend on the *do_jacobian* input.
)--",
      .output             = {"atmospheric_field",
                             "absorption_bands",
                             "measurement_sensor",
                             "surface_field",
                             "subsurface_field",
                             "measurement_vector_fitted",
                             "measurement_jacobian"},
      .input              = {"atmospheric_field",
                             "absorption_bands",
                             "measurement_sensor",
                             "surface_field",
                             "subsurface_field",
                             "jacobian_targets",
                             "model_state_vector",
                             "do_jacobian",
                             "inversion_iterate_agenda_counter"},
      .enum_options       = {"Full"},
      .enum_default       = "Full",
      .output_constraints = {
          {"(do_jacobian == 0 and measurement_jacobian.size() == 0) or (measurement_vector_fitted.size() == static_cast<Size>(measurement_jacobian.nrows()))",
           "On output, the measurement vector and Jacobian must match expected size.",
           "measurement_vector_fitted.size()",
           "measurement_jacobian.nrows()",
           "do_jacobian == 0"},
          {"(do_jacobian == 0 and measurement_jacobian.size() == 0) or (model_state_vector.size() == static_cast<Size>(measurement_jacobian.ncols()) and jacobian_targets.x_size() == model_state_vector.size())",
           "On output, the model state vector and Jacobian must match expected size.",
           "model_state_vector.size()",
           "measurement_jacobian.ncols()",
           "jacobian_targets.x_size()",
           "do_jacobian == 0"},
      }};

  wsa_data["measurement_inversion_agenda"] = {
      .desc =
          R"--(This is a helper *Agenda* intended for use within *inversion_iterate_agenda*.

It outputs the *measurement_vector_fitted* and *measurement_jacobian* for the
current iteration of the inversion. The *measurement_vector_fitted* is the
fitted measurement vector, i.e., the measurement vector that is expected to be
observed given the current *atmospheric_field*, *absorption_bands*, *measurement_sensor*,
and *surface_field*.  It does not take these as explicit input but via the Workspace
mechanism.  Within the *inversion_iterate_agenda*, these will be the local variables.

What is special about this Agenda is that it enforces that the *measurement_jacobian*
is empty on output if *do_jacobian* evaluates false.  Do not use this Agenda if you
do not mind having a non-empty *measurement_jacobian* on output even if *do_jacobian*
evaluates false.  Also do not use this Agenda if you wish to squeeze out performance,
it does a lot of unnecessary checks and operations that are not always needed.
)--",
      .output       = {"measurement_vector_fitted", "measurement_jacobian"},
      .input        = {"jacobian_targets", "do_jacobian"},
      .enum_options = {"Standard"},
      .enum_default = "Standard",
      .output_constraints =
          {
              {"do_jacobian != static_cast<Index>(measurement_jacobian.size() == 0)",
               "When *do_jacobian* evaluates as true, the *measurement_jacobian* must be non-empty.",
               "do_jacobian != 0",
               "measurement_jacobian.shape()"},
          },
  };

  wsa_data["disort_settings_agenda"] = {
      .desc   = R"--(An agenda for setting up Disort.

See *disort_settings_agendaSetup* for prepared agenda settings.

The only intent of this Agenda is to simplify the setup of Disort for different
scenarios.  The output of this Agenda is just that setting.
)--",
      .output = {"disort_settings"},
      .input  = {"frequency_grid",
                 "ray_path",
                 "disort_quadrature_dimension",
                 "disort_fourier_mode_dimension",
                 "disort_legendre_polynomial_dimension"},
  };

  // Add information about all automatically generated code
  for (auto& [name, record] : wsa_data) {
    record.desc += std::format(R"(
.. rubric:: Execution and customization
)");

    if (not record.enum_options.empty()) {
      record.desc += std::format(
          "See *{}Set* for builtin options that selects execution options.\n",
          name);
    }

    record.desc += std::format(R"(
You can execute *{0}* directly from the workspace by calling *{0}Execute*.

As all agendas in ARTS, it is also customizable via its operator helper class: *{0}Operator*.
See it, *{0}SetOperator*, and *{0}ExecuteOperator* for more details.

Also see the :func:`~pyarts.workspace.arts_agenda` property for how to fully define an agenda in python.
)",
                               name);
    if (not record.output_constraints.empty()) {
      record.desc +=
          std::format(R"(
.. rubric:: Constraint{0}

)",
                      record.output_constraints.size() > 1 ? "s"sv : ""sv);
      for (auto& c : record.output_constraints) {
        record.desc += std::format("#. {}\n", c.constraint);
      }
    }
    record.desc += "\n";
  }

  return wsa_data;
}

const std::unordered_map<std::string, WorkspaceAgendaInternalRecord>&
internal_workspace_agendas() {
  static const auto out = internal_workspace_agendas_creator();
  return out;
}

std::unordered_map<std::string, std::string> internal_workspace_agenda_names() {
  std::unordered_map<std::string, std::string> out;

  out["spectral_radiance_subsurface_agenda"] =
      "spectral_radiance_surface_agenda";

  return out;
}
