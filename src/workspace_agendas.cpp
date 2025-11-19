#include "workspace_agendas.h"

#include <format>

using namespace std::literals;

namespace {
std::unordered_map<std::string, WorkspaceAgendaInternalRecord>
internal_workspace_agendas_creator() {
  std::unordered_map<std::string, WorkspaceAgendaInternalRecord> wsa_data;

  wsa_data["spectral_propmat_agenda"] = {
      .desc =
          R"--(Computes the propagation matrix, the non-LTE source vector, and their derivatives.

The intent of this agenda is to be the workhorse for the propagation matrix
calculations that are happening deep in your ARTS method calls.

.. tip::
    Use *spectral_propmat_agendaAuto* after having defined
    your absorption data to create this agenda.  It covers most use-cases.
)--",
      .output       = {"spectral_propmat",
                       "spectral_nlte_srcvec",
                       "spectral_propmat_jac",
                       "spectral_nlte_srcvec_jac"},
      .input        = {"freq_grid",
                       "freq_wind_shift_jac",
                       "jac_targets",
                       "select_species",
                       "ray_point",
                       "atm_point"},
      .enum_options = {"Empty"},
      .output_constraints =
          {
              {"spectral_propmat.size() == freq_grid.size()",
               "On output, *spectral_propmat* has the size of *freq_grid*.",
               "spectral_propmat.size()",
               "freq_grid.size()"},
              {"spectral_nlte_srcvec.size() == freq_grid.size()",
               "On output, *spectral_nlte_srcvec* has the size of *freq_grid*.",
               "spectral_nlte_srcvec.size()",
               "freq_grid.size()"},
              {"matpack::same_shape<2>({jac_targets.target_count(), freq_grid.size()}, spectral_propmat_jac)",
               "On output, *spectral_propmat_jac* has the shape of the target-count of *jac_targets* times the size of *freq_grid*.",
               "spectral_propmat_jac.shape()",
               "freq_grid.size()",
               "jac_targets.target_count()"},
              {"matpack::same_shape<2>({jac_targets.target_count(), freq_grid.size()}, spectral_nlte_srcvec_jac)",
               "On output, *spectral_nlte_srcvec_jac* has the shape of the target-count of *jac_targets* times the size of *freq_grid*.",
               "spectral_nlte_srcvec_jac.shape()",
               "freq_grid.size()",
               "jac_targets.target_count()"},
          },
  };

  wsa_data["single_propmat_agenda"] = {
      .desc =
          R"--(Computes the propagation matrix, the non-LTE source vector, the dispersion, and their derivatives.

The intent of this agenda is to be the workhorse for the propagation matrix
calculations that are happening deep in your ARTS method calls.  The methods
in question here only compute a single frequency point at a time.

If you do not need single-frequency-point calculations, consider using
*spectral_propmat_agenda* instead as it will likely be more efficient.
)--",
      .output = {"single_propmat",
                 "single_nlte_srcvec",
                 "single_dispersion",
                 "single_propmat_jac",
                 "single_nlte_srcvec_jac",
                 "single_dispersion_jac"},
      .input  = {"freq",
                 "freq_wind_shift_jac",
                 "jac_targets",
                 "select_species",
                 "ray_point",
                 "atm_point"},
  };

  wsa_data["spectral_propmat_scat_spectral_agenda"] = {
      .desc =
          R"--(Gets the scattering propagation matrix, the scattering absorption vector, and the scattering spectral phase matrix.
)--",
      .output       = {"spectral_propmat_scat",
                       "spectral_absvec_scat",
                       "spectral_phamat_spectral"},
      .input        = {"freq_grid", "atm_point", "legendre_degree"},
      .enum_options = {"FromSpeciesTRO"},
      .enum_default = "FromSpeciesTRO",
      .output_constraints =
          {
              {"spectral_propmat_scat.size() == freq_grid.size()",
               "On output, *spectral_propmat_scat* has the size of *freq_grid*.",
               "spectral_propmat_scat.size()",
               "freq_grid.size()"},
              {"spectral_absvec_scat.size() == freq_grid.size()",
               "On output, *spectral_absvec_scat* has the size of *freq_grid*.",
               "spectral_absvec_scat.size()",
               "freq_grid.size()"},
              {"matpack::same_shape<2>({freq_grid.size(), legendre_degree + 1}, spectral_phamat_spectral)",
               "On output, *spectral_phamat_spectral* has the shape of <*legendre_degree* + 1> times the size of *freq_grid*.",
               "spectral_phamat_spectral.shape()",
               "freq_grid.size()",
               "legendre_degree"},
          },
  };

  wsa_data["spectral_propmat_scat_agenda"] = {
      .desc =
          R"--(Computes the part of the propagation matrix that relates to scattering.
)--",
      .output       = {"spectral_propmat_scat"},
      .input        = {"freq_grid", "atm_point"},
      .enum_options = {"AirSimple"},
      .enum_default = "AirSimple",
      .output_constraints =
          {
              {"spectral_propmat_scat.size() == freq_grid.size()",
               "On output, *spectral_propmat_scat* has the size of *freq_grid*.",
               "spectral_propmat_scat.size()",
               "freq_grid.size()"},
          },
  };

  wsa_data["ray_path_observer_agenda"] = {
      .desc   = R"--(Gets the propagation path as it is obeserved.

The intent of this agenda is to provide a propagation path as seen from the observer
position and line of sight.

.. tip::
    The perhaps easiest way to set this agenda up is to use the *ray_path_observer_agendaSetGeometric* method.
)--",
      .output = {"ray_path"},
      .input  = {"obs_pos", "obs_los"},
  };

  wsa_data["ray_point_back_propagation_agenda"] = {
      .desc =
          R"--(Gets the next past point along a propagation path.

*ray_path* must have a point already.  This point is propagated backwards.

It is up to internal methods if they respect *single_dispersion* or not.

It is up to internal methods if they respect *max_stepsize* or not.

A special exception may be made for a 1-size *ray_path* that is in space or at the surface,
where the next point may be the same point as the input.

The end of the path is reached when the last point in *ray_path* is
at *PathPositionType* ``space`` or ``surface``.
)--",
      .output       = {"ray_point"},
      .input        = {"ray_path",
                       "single_dispersion",
                       "single_propmat",
                       "max_stepsize"},
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
      .input        = {"freq_grid",
                       "jac_targets",
                       "obs_pos",
                       "obs_los",
                       "atm_field",
                       "surf_field",
                       "subsurf_field"},
      .enum_options = {"Emission",
                       "EmissionAdaptiveHalfsteps",
                       "EmissionNoSensor"},
      .enum_default = "Emission",
      .output_constraints =
          {
              {"spectral_rad.size() == freq_grid.size()",
               "On output, *spectral_rad* has the size of *freq_grid*.",
               "spectral_rad.size()",
               "freq_grid.size()"},
              {"matpack::same_shape<2>({jac_targets.x_size(), freq_grid.size()}, spectral_rad_jac)",
               "On output, *spectral_rad_jac* has the shape of the expected *model_state_vec* (i.e., the x-size of *jac_targets*) times the size of *freq_grid*.",
               "spectral_rad_jac.shape()",
               "freq_grid.size()",
               "jac_targets.x_size()"},
          },
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
      .desc               = R"--(Gets spectral radiance as seen of space.

This agenda calculates the spectral radiance as seen of space.
One common use-case is to provide a background spectral radiance.

The input path point should be as if it is looking at space.
)--",
      .output             = {"spectral_rad", "spectral_rad_jac"},
      .input              = {"freq_grid", "jac_targets", "ray_point"},
      .enum_options       = {"UniformCosmicBackground",
                             "SunOrCosmicBackground",
                             "Transmission"},
      .enum_default       = "UniformCosmicBackground",
      .output_constraints = {
          {"spectral_rad.size() == freq_grid.size()",
           "On output, *spectral_rad* has the size of *freq_grid*.",
           "spectral_rad.size()",
           "freq_grid.size()"},
          {"matpack::same_shape<2>({jac_targets.x_size(), freq_grid.size()}, spectral_rad_jac)",
           "On output, *spectral_rad_jac* has the shape of the expected *model_state_vec* (i.e., the x-size of *jac_targets*) times the size of *freq_grid*.",
           "spectral_rad_jac.shape()",
           "freq_grid.size()",
           "jac_targets.x_size()"},
      }};

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
      .input        = {"freq_grid",
                       "jac_targets",
                       "ray_point",
                       "surf_field",
                       "subsurf_field"},
      .enum_options = {"Blackbody", "Transmission", "SurfaceReflectance"},
      .enum_default = "Blackbody",
      .output_constraints = {
          {"spectral_rad.size() == freq_grid.size()",
           "On output, *spectral_rad* has the size of *freq_grid*.",
           "spectral_rad.size()",
           "freq_grid.size()"},
          {"matpack::same_shape<2>({jac_targets.x_size(), freq_grid.size()}, spectral_rad_jac)",
           "On output, *spectral_rad_jac* has the shape of the expected *model_state_vec* (i.e., the x-size of *jac_targets*) times the size of *freq_grid*.",
           "spectral_rad_jac.shape()",
           "freq_grid.size()",
           "jac_targets.x_size()"},
      }};

  wsa_data["single_rad_surface_agenda"] = {
      .desc =
          R"--(Gets spectral radiance as seen of the surface for a single frequency.

Otherwise same as *spectral_rad_surface_agenda*.
)--",
      .output = {"single_rad", "single_rad_jac"},
      .input =
          {"freq", "jac_targets", "ray_point", "surf_field", "subsurf_field"},
      .enum_options = {"WrapGrid"},
      .enum_default = "WrapGrid",
  };

  wsa_data["inversion_iterate_agenda"] = {
      .desc               = R"--(Work in progress ...

See *OEM*.

.. note::
    The output *measurement_jac* size may depend on the *do_jac* input.
)--",
      .output             = {"atm_field",
                             "abs_bands",
                             "measurement_sensor",
                             "surf_field",
                             "subsurf_field",
                             "measurement_vec_fit",
                             "measurement_jac"},
      .input              = {"atm_field",
                             "abs_bands",
                             "measurement_sensor",
                             "surf_field",
                             "subsurf_field",
                             "jac_targets",
                             "model_state_vec",
                             "do_jac",
                             "inversion_iterate_agenda_counter"},
      .enum_options       = {"Full"},
      .enum_default       = "Full",
      .output_constraints = {
          {"(do_jac == 0 and measurement_jac.size() == 0) or (measurement_vec_fit.size() == static_cast<Size>(measurement_jac.nrows()))",
           "On output, the measurement vector and Jacobian must match expected size.",
           "measurement_vec_fit.size()",
           "measurement_jac.nrows()",
           "do_jac == 0"},
          {"(do_jac == 0 and measurement_jac.size() == 0) or (model_state_vec.size() == static_cast<Size>(measurement_jac.ncols()) and jac_targets.x_size() == model_state_vec.size())",
           "On output, the model state vector and Jacobian must match expected size.",
           "model_state_vec.size()",
           "measurement_jac.ncols()",
           "jac_targets.x_size()",
           "do_jac == 0"},
      }};

  wsa_data["measurement_inversion_agenda"] = {
      .desc =
          R"--(This is a helper *Agenda* intended for use within *inversion_iterate_agenda*.

It outputs the *measurement_vec_fit* and *measurement_jac* for the
current iteration of the inversion. The *measurement_vec_fit* is the
fitted measurement vector, i.e., the measurement vector that is expected to be
observed given the current *atm_field*, *abs_bands*, *measurement_sensor*,
and *surf_field*.  It does not take these as explicit input but via the Workspace
mechanism.  Within the *inversion_iterate_agenda*, these will be the local variables.

What is special about this Agenda is that it enforces that the *measurement_jac*
is empty on output if *do_jac* evaluates false.  Do not use this Agenda if you
do not mind having a non-empty *measurement_jac* on output even if *do_jac*
evaluates false.  Also do not use this Agenda if you wish to squeeze out performance,
it does a lot of unnecessary checks and operations that are not always needed.
)--",
      .output       = {"measurement_vec_fit", "measurement_jac"},
      .input        = {"jac_targets", "do_jac"},
      .enum_options = {"Standard"},
      .enum_default = "Standard",
      .output_constraints =
          {
              {"do_jac != static_cast<Index>(measurement_jac.size() == 0)",
               "When *do_jac* evaluates as true, the *measurement_jac* must be non-empty.",
               "do_jac != 0",
               "measurement_jac.shape()"},
          },
  };

  wsa_data["spectral_surf_refl_agenda"] = {
      .desc         = R"--(An agenda to compute the surface reflectance.
)--",
      .output       = {"spectral_surf_refl", "spectral_surf_refl_jac"},
      .input        = {"freq_grid", "surf_field", "ray_point", "jac_targets"},
      .enum_options = {"FlatScalar", "FlatRealFresnel"},
      .output_constraints = {
          {"spectral_surf_refl.size() == freq_grid.size()",
           "*spectral_surf_refl* match *freq_grid* size",
           "spectral_surf_refl.size()",
           "freq_grid.size()"},
          {"matpack::same_shape<2>({jac_targets.target_count(), freq_grid.size()}, spectral_surf_refl_jac)",
           "*spectral_surf_refl_jac* match *jac_targets* target count and *freq_grid* size",
           "jac_targets.target_count()",
           "freq_grid.size()",
           "spectral_surf_refl_jac.shape()"}}};

  wsa_data["disort_settings_agenda"] = {
      .desc   = R"--(An agenda for setting up Disort.

See *disort_settings_agendaSetup* for prepared agenda settings.

The only intent of this Agenda is to simplify the setup of Disort for different
scenarios.  The output of this Agenda is just that setting.
)--",
      .output = {"disort_settings"},
      .input  = {"freq_grid",
                 "ray_path",
                 "disort_quadrature_dimension",
                 "disort_fourier_mode_dimension",
                 "disort_legendre_polynomial_dimension"}};

  wsa_data["disort_settings_downwelling_wrapper_agenda"] = {
      .desc   = R"--(An wrapper agenda for calling *disort_settings_agenda*.

This agenda wraps the *disort_settings_agenda* to provide a simpler interface
for the common case of calculating downwelling radiation.  The idea is that a
call to *disort_settings_agenda* is made, and then a follow-up calculation of
the down-welling radiation is done to set the boundary condition at the top
of the tau-range covered by the ray path.

One use-case is to use this agenda to give downwelling atmospheric radiation
as a boundary condition to subsurface radiance calculation.
)--",
      .output = {"disort_settings"},
      .input  = {"freq_grid",
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
      record.desc += std::format(
          "See *{}Set* for builtin options that selects execution options.\n",
          name);
    }

    record.desc += std::format(R"(
You can execute *{0}* directly from the workspace by calling *{0}Execute*.

As all agendas in ARTS, it is also customizable via its operator helper class: *{0}Operator*.
See it, *{0}SetOperator*, and *{0}ExecuteOperator* for more details.

Also see the :class:`~pyarts3.workspace.arts_agenda` property for how to fully define an agenda in python.
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
}  // namespace

const std::unordered_map<std::string, WorkspaceAgendaInternalRecord>&
internal_workspace_agendas() {
  static const auto out = internal_workspace_agendas_creator();
  return out;
}

std::unordered_map<std::string, std::string> internal_workspace_agenda_names() {
  std::unordered_map<std::string, std::string> out;

  out["spectral_rad_closed_surface_agenda"] = "spectral_rad_surface_agenda";

  return out;
}
