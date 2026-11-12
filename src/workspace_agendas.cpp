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
      .enum_desc =
          {R"--(No absorption.  The propagation matrix, the non-LTE source vector, and their derivatives are all zero.)--"},
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
      .enum_desc =
          {R"--(The atmosphere of each point of *ray_path* is that of *atm_field*, and its frequencies are *freq_grid* shifted by the wind there.  *spectral_propmat_agenda* computes the propagation matrices of every point.  *ray_path* is kept as it is.)--",
           R"--(As ``Default``, but points are added to *ray_path* where a step is optically too thick, see *spectral_propmat_pathAddAdaptiveHalfPath*.)--",
           R"--(The path data are taken from a profile computed beforehand by *ProfileFromAltitude*, see *Profile2Path*.  *ray_path* is kept as it is.)--"},
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
      .enum_desc =
          {R"--(The scattering of *scat_species* as totally randomly oriented particles, see *spectral_propmat_scatAddSpectralScatteringSpeciesTRO*.)--"},
      .enum_default = "FromSpeciesTRO",
  };

  wsa_data["spectral_propmat_scat_agenda"] = {
      .desc =
          R"--(Computes the part of the propagation matrix that relates to scattering.
)--",
      .output       = {"spectral_propmat_scat"},
      .input        = {"freq_grid", "atm_point"},
      .enum_options = {"AirSimple"},
      .enum_desc    = {R"--(Simple scattering by air, see *spectral_propmat_scatAirSimple*.)--"},
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
      .enum_desc =
          {R"--(A geometric path, as *ray_path_observer_agendaSetGeometric* sets it with its default settings.)--",
           R"--(A geometric path with points only where it crosses *alt_grid*, for methods that work on altitude profiles, such as *Profile2Path*.  The points outside *alt_grid*, such as an observer above it, are removed, see *ray_pathRemoveOutsideAltitudeGrid*.  A limb path gets no tangent point, so the spacing of *alt_grid* limits its accuracy.)--"},
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
      .enum_desc =
          {R"--(The previous point along a straight line, see *ray_pointPastGeometric*.)--",
           R"--(The previous point along a line refracted by *single_dispersion*, see *ray_pointPastRefractive*.)--"},
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
      .enum_options = {"Emission", "EmissionAdaptiveHalfsteps", "EmissionProfile2Path", "MonteCarlo"},
      .enum_desc =
          {R"--(Clear-sky emission along the path of *ray_path_observer_agenda*, see *spectral_radClearskyEmission*.  *spectral_propmat_and_atm_path_agenda* gives the path data.)--",
           R"--(Clear-sky emission along the path of *ray_path_observer_agenda*, with points added where a step is optically too thick, see *spectral_propmat_pathAddAdaptiveHalfPath*.  The radiance is computed one frequency at a time, see *spectral_radSinglePathEmissionFrequencyLoop*.)--",
           R"--(Clear-sky emission along the path of *ray_path_observer_agenda*, with the path data taken from a profile computed beforehand by *ProfileFromAltitude*, see *Profile2Path*.  Points of the path between the altitudes of *alt_grid*, such as the tangent points of limb paths, are computed by *spectral_propmat_agenda*.)--",
           R"--(Passive Monte Carlo radiance for a pencil beam, see *spectral_radMonteCarlo*.)--"},
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
      .enum_desc    = {R"--(*spectral_rad_space_agenda* at the single frequency *freq*.)--"},
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
      .enum_desc =
          {R"--(The cosmic microwave background, the same in all directions, see *spectral_radUniformCosmicBackground*.)--",
           R"--(A sun if the ray hits one, otherwise the cosmic microwave background, see *spectral_radSunsOrCosmicBackground*.)--",
           R"--(A unit background without Jacobian, for transmission calculations, see *spectral_radDefaultTransmission*.)--"},
      .enum_default = "UniformCosmicBackground",
  };

  wsa_data["spectral_rad_surface_agenda"] = {
      .desc         = R"--(Computes spectral radiance as seen of the surface.

This agenda calculates the spectral radiance as seen of the surface.
One common use-case us to provide a background spectral radiance.

The input path point should be as if it is looking at the surface.

Subsurface calculations are also supported through this agenda,
but might require setting *spectral_rad_closed_surface_agenda*
as well.
)--",
      .output       = {"spectral_rad", "spectral_rad_jac"},
      .input        = {"freq_grid", "jac_targets", "ray_point", "surf_field", "subsurf_field"},
      .enum_options = {"Blackbody", "Transmission", "SurfaceReflectance"},
      .enum_desc =
          {R"--(Blackbody emission at the surface temperature, see *spectral_radSurfaceBlackbody*.)--",
           R"--(A unit background without Jacobian, for transmission calculations, see *spectral_radDefaultTransmission*.)--",
           R"--(Emission and reflection at the surface, see *spectral_radSurfaceReflectance*.)--"},
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
      .enum_desc      = {R"--(Blackbody emission at the surface temperature, see *spectral_radSurfaceBlackbody*.)--"},
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
      .enum_desc    = {R"--(*spectral_rad_surface_agenda* at the single frequency *freq*.)--"},
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
                       "measurement_sensors",
                       "surf_field",
                       "subsurf_field",
                       "measurement_vec_fit",
                       "measurement_jac"},
      .input        = {"atm_field",
                       "abs_bands",
                       "measurement_sensors",
                       "surf_field",
                       "subsurf_field",
                       "model_state_targets",
                       "jac_targets",
                       "model_state_vec"},
      .enum_options = {"Full"},
      .enum_desc =
          {R"--(Updates the model to *model_state_vec*, see *UpdateModelStates*, and simulates the measurement with *measurement_inversion_agenda*.)--"},
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
      .enum_options = {"LowMemory", "HighPerformance", "Indexed"},
      .enum_desc =
          {R"--(The measurement is simulated with the ``"Low Memory"`` kernel of *MeasurementVectorSumupKernel*, which keeps the least in memory.)--",
           R"--(The measurement is simulated with the ``"High Performance"`` kernel of *MeasurementVectorSumupKernel*, which keeps the spectral radiance of every geometry in memory to sum it faster.)--",
           R"--(The measurement is simulated with the ``"Indexed"`` kernel of *MeasurementVectorSumupKernel*, which lists the weights of each geometry before the agenda calls.)--"},
      .enum_default = "Indexed",

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
      .enum_desc =
          {R"--(A flat surface that reflects all polarizations alike, see *spectral_surf_reflFlatScalar*.)--",
           R"--(A flat surface with the Fresnel reflectance of a real refractive index, see *spectral_surf_reflFlatRealFresnel*.)--",
           R"--(Specular ocean reflectance from TESSEM2, see *spectral_surf_reflTessem*.)--",
           R"--(Specular land reflectance from the TELSEM2 atlas, see *spectral_surf_reflTelsem*.)--"},
  };

  wsa_data["disort_settings_agenda"] = {.desc           = R"--(An agenda for setting up Disort.

See *disort_settings_agendaSetup* for prepared agenda settings.

The only intent of this Agenda is to simplify the setup of Disort for different
scenarios.  The output of this Agenda is just that setting.
)--",
                                        .output         = {"disort_settings"},
                                        .input          = {"freq_grid",
                                                           "ray_path",
                                                           "quadrature_dimension",
                                                           "fourier_mode_dimension",
                                                           "legendre_polynomial_dimension"},
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
                       "quadrature_dimension",
                       "fourier_mode_dimension",
                       "legendre_polynomial_dimension",
                       "disort_settings_agenda"},
      .enum_options = {"Standard" /*, "Disort"*/},
      .enum_desc =
          {R"--(*disort_settings_agenda* with the downwelling boundary condition from *spectral_rad_observer_agenda*, see *disort_settingsDownwellingObserver*.)--"},
      .enum_default = "Standard"};

  wsa_data["vdisort_settings_agenda"] = {.desc           = R"--(An agenda for setting up VDISORT.

See *vdisort_settings_agendaSetup* for prepared agenda settings.

This is the polarized counterpart of *disort_settings_agenda*.  The only intent
of this Agenda is to simplify the setup of VDISORT for different scenarios.  The
output of this Agenda is just that setting.
)--",
                                         .output         = {"vdisort_settings"},
                                         .input          = {"freq_grid",
                                                            "ray_path",
                                                            "quadrature_dimension",
                                                            "fourier_mode_dimension",
                                                            "legendre_polynomial_dimension"},
                                         .named_operator = "VDisortSettingsAgendaOperator"};

  wsa_data["vdisort_settings_downwelling_wrapper_agenda"] = {
      .desc         = R"--(An wrapper agenda for calling *vdisort_settings_agenda*.

The polarized counterpart of *disort_settings_downwelling_wrapper_agenda*.
A call to *vdisort_settings_agenda* is made, and then a follow-up calculation
of the down-welling radiation is done to set the boundary condition at the top
of the tau-range covered by the ray path.
)--",
      .output       = {"vdisort_settings"},
      .input        = {"freq_grid",
                       "ray_path",
                       "atm_field",
                       "surf_field",
                       "subsurf_field",
                       "quadrature_dimension",
                       "fourier_mode_dimension",
                       "legendre_polynomial_dimension",
                       "vdisort_settings_agenda"},
      .enum_options = {"Standard"},
      .enum_desc =
          {R"--(*vdisort_settings_agenda* with the downwelling boundary condition from *spectral_rad_observer_agenda*, see *vdisort_settingsDownwellingObserver*.)--"},
      .enum_default = "Standard"};

  // Add information about all automatically generated code
  for (auto& [name, record] : wsa_data) {
    if (record.enum_desc.size() != record.enum_options.size())
      throw std::runtime_error(std::format("Agenda {} has {} enum_options but {} enum_desc, it needs one per option",
                                           name,
                                           record.enum_options.size(),
                                           record.enum_desc.size()));

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
