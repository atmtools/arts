#include "workspace_agendas.h"

#include <unique_unordered_map.h>

#include <format>

using namespace std::literals;

namespace {
std::unordered_map<std::string, WorkspaceAgendaInternalRecord> internal_workspace_agendas_creator() {
  UniqueMap<std::string, WorkspaceAgendaInternalRecord> wsa_data;

  wsa_data["spectral_propmat_agenda"] = {
      .desc =
          R"--(Computes the propagation matrix, the non-LTE source vector, and their derivatives.

The intent of this agenda is to be the workhorse for the propagation matrix
calculations that are happening deep in your ARTS method calls.

.. tip::
    Use *spectral_propmat_agendaAuto* after having defined
    your absorption data to create this agenda.  It covers most use-cases.
)--",
      .output       = {"spectral_propmat", "spectral_nlte_srcvec", "spectral_propmat_jac", "spectral_nlte_srcvec_jac"},
      .input        = {"freq_grid", "freq_wind_shift_jac", "jac_targets", "select_species", "ray_point", "atm_point"},
      .enum_options = {"Empty"},
  };

  wsa_data["spectral_propmat_and_atm_path_agenda"] = {
      .desc =
          R"--(Computes several path parameters along the path.

The main use of this agenda is to allow adapting the path points
based on spectral parameters.
)--",
      .output       = {"spectral_propmat_path",
                       "spectral_nlte_srcvec_path",
                       "spectral_propmat_jac_path",
                       "spectral_nlte_srcvec_jac_path",
                       "freq_grid_path",
                       "freq_wind_shift_jac_path",
                       "atm_path",
                       "ray_path"},
      .input        = {"ray_path", "jac_targets", "freq_grid", "atm_field", "surf_field"},
      .enum_options = {"Default", "AdaptiveHalfPath", "Profile2Path"},
      .enum_default = "Default",
  };

  wsa_data["single_propmat_agenda"] = {
      .desc =
          R"--(Computes the propagation matrix, the non-LTE source vector, the dispersion, and their derivatives.

The intent of this agenda is to be the workhorse for the propagation matrix
calculations that are happening deep in your ARTS method calls.  The methods
in question here only compute a single frequency point at a time.

If you do not need single-frequency-point calculations, consider using
*spectral_propmat_agenda* instead as it will likely be more efficient.

Convenience setters provide *single_dispersion* for microwave gases in Earth
or planetary atmospheres and for visible/near-infrared water or steam.
)--",
      .output = {"single_propmat",
                 "single_nlte_srcvec",
                 "single_dispersion",
                 "single_propmat_jac",
                 "single_nlte_srcvec_jac",
                 "single_dispersion_jac"},
      .input  = {"freq", "freq_wind_shift_jac", "jac_targets", "select_species", "ray_point", "atm_point"},
  };

  wsa_data["spectral_propmat_scat_spectral_agenda"] = {
      .desc =
          R"--(Gets the scattering propagation matrix, the scattering absorption vector, and the scattering spectral phase matrix.
)--",
      .output       = {"spectral_propmat_scat", "spectral_absvec_scat", "spectral_phamat_spectral"},
      .input        = {"freq_grid", "atm_point", "legendre_degree"},
      .enum_options = {"FromSpeciesTRO"},
      .enum_default = "FromSpeciesTRO",
  };

  wsa_data["spectral_propmat_scat_agenda"] = {
      .desc =
          R"--(Computes the part of the propagation matrix that relates to scattering.
)--",
      .output       = {"spectral_propmat_scat"},
      .input        = {"freq_grid", "atm_point"},
      .enum_options = {"AirSimple"},
      .enum_default = "AirSimple",
  };

  wsa_data["ray_path_observer_agenda"] = {
      .desc         = R"--(Gets the propagation path as it is observed.

The intent of this agenda is to provide a propagation path as seen from the observer
position and line of sight.

.. tip::
    The perhaps easiest way to set this agenda up is to use the *ray_path_observer_agendaSetGeometric* method.
)--",
      .output       = {"ray_path"},
      .input        = {"obs_pos", "obs_los"},
      .enum_options = {"GeometricDefault", "GeometricProfile"},
      .enum_default = "GeometricDefault",
  };

  wsa_data["ray_point_back_propagation_agenda"] = {
      .desc =
          R"--(Gets the next past point along a propagation path.

*ray_path* must have a point already.  This point is propagated backwards.

It is up to internal methods if they respect *single_dispersion* or not.

It is up to internal methods if they respect *max_stepsize* or not.

The ``RefractiveStepwise`` option consumes *single_dispersion*.  It can, for
example, be paired with *single_propmat_agendaSetGasMicrowavesEarth*,
*single_propmat_agendaSetGasMicrowavesGeneral*, or
*single_propmat_agendaSetWaterVisibleNIRHarvey98*.

A special exception may be made for a 1-size *ray_path* that is in space or at the surface,
where the next point may be the same point as the input.

The end of the path is reached when the last point in *ray_path* is
at *PathPositionType* ``space`` or ``surface``.
)--",
      .output       = {"ray_point"},
      .input        = {"ray_path", "single_dispersion", "single_propmat", "max_stepsize"},
      .enum_options = {"GeometricStepwise", "RefractiveStepwise"},
      .enum_default = "GeometricStepwise",
  };

  wsa_data["spectral_rad_observer_agenda"] = {
      .desc =
          R"--(Computes spectral radiance as seen from the input position and environment.

The intent of this agenda is to provide the spectral radiance as seen from the observer
position and line of sight.

It also outputs the *ray_path* as seen from the observer position and line of sight.
This is useful in-case a call to the destructive *spectral_radApplyUnitFromSpectralRadiance*
is warranted.
)--",
      .output       = {"spectral_rad", "spectral_rad_jac", "ray_path"},
      .input        = {"freq_grid", "jac_targets", "obs_pos", "obs_los", "atm_field", "surf_field", "subsurf_field"},
      .enum_options = {"Emission", "EmissionAdaptiveHalfsteps", "EmissionNoSensor", "MonteCarlo"},
      .enum_default = "Emission",
  };

  wsa_data["single_rad_space_agenda"] = {
      .desc =
          R"--(Gets spectral radiance as seen of space for a single frequency.

Otherwise same as *spectral_rad_space_agenda*.
)--",
      .output       = {"single_rad", "single_rad_jac"},
      .input        = {"freq", "jac_targets", "ray_point"},
      .enum_options = {"WrapGrid"},
      .enum_default = "WrapGrid",
  };

  wsa_data["spectral_rad_space_agenda"] = {
      .desc         = R"--(Gets spectral radiance as seen of space.

This agenda calculates the spectral radiance as seen of space.
One common use-case is to provide a background spectral radiance.

The input path point should be as if it is looking at space.
)--",
      .output       = {"spectral_rad", "spectral_rad_jac"},
      .input        = {"freq_grid", "jac_targets", "ray_point"},
      .enum_options = {"UniformCosmicBackground", "SunOrCosmicBackground", "Transmission"},
      .enum_default = "UniformCosmicBackground",
  };

  wsa_data["spectral_rad_surface_agenda"] = {
      .desc           = R"--(Computes spectral radiance as seen of the surface.

This agenda calculates the spectral radiance as seen of the surface.
One common use-case us to provide a background spectral radiance.

The input path point should be as if it is looking at the surface.

Subsurface calculations are also supported through this agenda,
but might require setting *spectral_rad_closed_surface_agenda*
as well.
)--",
      .output         = {"spectral_rad", "spectral_rad_jac"},
      .input          = {"freq_grid", "jac_targets", "ray_point", "surf_field", "subsurf_field"},
      .enum_options   = {"Blackbody", "Transmission", "SurfaceReflectance"},
      .enum_default   = "Blackbody",
      .named_operator = "SpectralRadianceSurfaceAgendaOperator"};

  wsa_data["spectral_rad_closed_surface_agenda"] = {
      .desc           = R"--(A closed surface agenda.

It behave exactly like *spectral_rad_surface_agenda*.  It exists
to allow chaining surface agendas.  The idea is that the main
*spectral_rad_surface_agenda* variable is the first interface
and can chain into another surface agenda - this one.

Thus this agenda must be "closed".  It cannot call another *spectral_rad_surface_agenda*,
whereas *spectral_rad_surface_agenda* can call this agenda.  Imagine a chain where
the *spectral_rad_surface_agenda* gets the reflectance from a land surface model
and calls the *spectral_rad_observer_agenda* to compute the downwelling radiation at the surface.
It can in turn call *spectral_rad_closed_surface_agenda* to get the upwelling radiation from the surface
that is being emitted.  That's the type of use case this agenda is made for and why it exists!
)--",
      .output         = wsa_data.at("spectral_rad_surface_agenda").output,
      .input          = wsa_data.at("spectral_rad_surface_agenda").input,
      .enum_options   = {"Blackbody"},
      .enum_default   = "Blackbody",
      .named_operator = wsa_data.at("spectral_rad_surface_agenda").named_operator};

  wsa_data["single_rad_surface_agenda"] = {
      .desc =
          R"--(Gets spectral radiance as seen of the surface for a single frequency.

Otherwise same as *spectral_rad_surface_agenda*.
)--",
      .output       = {"single_rad", "single_rad_jac"},
      .input        = {"freq", "jac_targets", "ray_point", "surf_field", "subsurf_field"},
      .enum_options = {"WrapGrid"},
      .enum_default = "WrapGrid",
  };

  wsa_data["inversion_iterate_agenda"] = {
      .desc         = R"--(Evaluate a retrieval state.  See *oemCalc*.

*model_state_targets* always contains the complete state mapping used by
*UpdateModelStates* and measurement-error values.  *jac_targets* contains
the derivative targets, or is empty for a value-only evaluation.
)--",
      .output       = {"atm_field",
                       "abs_bands",
                       "measurement_sensor",
                       "surf_field",
                       "subsurf_field",
                       "measurement_vec_fit",
                       "measurement_jac"},
      .input        = {"atm_field",
                       "abs_bands",
                       "measurement_sensor",
                       "surf_field",
                       "subsurf_field",
                       "model_state_targets",
                       "jac_targets",
                       "model_state_vec"},
      .enum_options = {"Full"},
      .enum_default = "Full",

      // Wraps *measurement_inversion_agenda*, so it inherits the empty
      // *measurement_jac* that a pass with empty *jac_targets* produces
      .output_constraints = false,
  };

  wsa_data["measurement_inversion_agenda"] = {
      .desc         = R"--(Simulate the fitted measurement for the current physical model.

Apply *UpdateModelStates* before this helper, as the predefined
*inversion_iterate_agenda* does.  *model_state_targets* supplies the full
mapping for measurement-error values.  *jac_targets* controls all derivatives;
when it is empty, *measurement_jac* is empty.  Both target sets are read-only.
)--",
      .output       = {"measurement_vec_fit", "measurement_jac"},
      .input        = {"model_state_targets", "jac_targets"},
      .enum_options = {"LowMemory", "HighPerformance"},
      .enum_default = "LowMemory",

      // Returning an empty *measurement_jac* when *jac_targets* is empty is the point
      // of this agenda, so its outputs have no shape to verify
      .output_constraints = false,
  };

  wsa_data["spectral_surf_refl_agenda"] = {
      .desc         = R"--(An agenda to compute the surface reflectance.
)--",
      .output       = {"spectral_surf_refl", "spectral_surf_refl_jac"},
      .input        = {"freq_grid", "surf_field", "ray_point", "jac_targets"},
      .enum_options = {"FlatScalar", "FlatRealFresnel", "Tessem", "Telsem"},
  };

  wsa_data["disort_settings_agenda"] = {.desc           = R"--(An agenda for setting up Disort.

See *disort_settings_agendaSetup* for prepared agenda settings.

The only intent of this Agenda is to simplify the setup of Disort for different
scenarios.  The output of this Agenda is just that setting.
)--",
                                        .output         = {"disort_settings"},
                                        .input          = {"freq_grid",
                                                           "ray_path",
                                                           "disort_quadrature_dimension",
                                                           "disort_fourier_mode_dimension",
                                                           "disort_legendre_polynomial_dimension"},
                                        .named_operator = "DisortSettingsAgendaOperator"};

  wsa_data["atm_disort_settings_agenda"] = {
      .desc =
          R"--(A specialization of *disort_settings_agenda* for atmospheric calculations.

.. seealso::
    *subsurf_disort_settings_agenda* for a similar agenda for subsurface calculations..
)--",
      .output         = wsa_data.at("disort_settings_agenda").output,
      .input          = wsa_data.at("disort_settings_agenda").input,
      .named_operator = wsa_data.at("disort_settings_agenda").named_operator};

  wsa_data["subsurf_disort_settings_agenda"] = {
      .desc =
          R"--(A specialization of *disort_settings_agenda* for subsurface calculations.

.. seealso::
    *atm_disort_settings_agenda* for a similar agenda for atmospheric calculations.
)--",
      .output         = wsa_data.at("disort_settings_agenda").output,
      .input          = wsa_data.at("disort_settings_agenda").input,
      .named_operator = wsa_data.at("disort_settings_agenda").named_operator};

  wsa_data["disort_settings_downwelling_wrapper_agenda"] = {
      .desc         = R"--(An wrapper agenda for calling *disort_settings_agenda*.

This agenda wraps the *disort_settings_agenda* to provide a simpler interface
for the common case of calculating downwelling radiation.  The idea is that a
call to *disort_settings_agenda* is made, and then a follow-up calculation of
the down-welling radiation is done to set the boundary condition at the top
of the tau-range covered by the ray path.

One use-case is to use this agenda to give downwelling atmospheric radiation
as a boundary condition to subsurface radiance calculation.
)--",
      .output       = {"disort_settings"},
      .input        = {"freq_grid",
                       "ray_path",
                       "atm_field",
                       "surf_field",
                       "subsurf_field",
                       "disort_quadrature_dimension",
                       "disort_fourier_mode_dimension",
                       "disort_legendre_polynomial_dimension",
                       "disort_settings_agenda"},
      .enum_options = {"Standard" /*, "Disort"*/},
      .enum_default = "Standard"};

  // Add information about all automatically generated code
  for (auto& [name, record] : wsa_data) {
    record.desc += std::format(R"(
.. rubric:: Execution and customization

)");

    if (not record.enum_options.empty()) {
      record.desc += std::format("See *{}Set* for builtin options that selects execution options.\n", name);
    }

    record.desc += std::format(R"(
You can execute *{0}* directly from the workspace by calling *{0}Execute*.

As all agendas in ARTS, it is also customizable via its operator helper class: *{1}*.
See it, *{0}SetOperator*, and *{0}ExecuteOperator* for more details.

Also see the :func:`~pyarts3.workspace.arts_agenda` decorator for how to fully define an agenda in python.
)",
                               name,
                               record.named_operator.empty() ? name + "Operator" : record.named_operator);

    // The output size constraints are documented where the agenda is documented.
    // They cannot be added here because reading them needs the workspace variables,
    // which are created after the agendas and from them.

    record.desc += "\n";
  }

  return std::move(wsa_data.map);
}
}  // namespace

const std::unordered_map<std::string, WorkspaceAgendaInternalRecord>& internal_workspace_agendas() {
  static const auto out = internal_workspace_agendas_creator();
  return out;
}
