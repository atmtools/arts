#include "workspace_meta_methods.h"

#include <sorting.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "workspace_agendas.h"
#include "workspace_methods.h"
#include "workspace_variables.h"

namespace stdr = std::ranges;

namespace {
std::vector<WorkspaceMethodInternalMetaRecord> internal_meta_methods_creator() {
  std::vector<WorkspaceMethodInternalMetaRecord> wsm_meta;

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "measurement_sensorSimple",
      .desc    = R"(Creates a single simple Dirac-opening sensor

This means that the sensor has no bandwidth per channel, i.e., only one
frequency point is used to simulate the spectral radiance before being
averaged into a channel.
)",
      .author  = {"Richard Larsson"},
      .methods = {"measurement_sensorInit", "measurement_sensorAddSimple"},
      .out     = {"measurement_sensor"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "measurement_sensorSimpleGaussian",
      .desc    = R"(Creates a single simple Gaussian-opening sensor.

This means that the sensor has a Gaussian bandwidth per channel.
That is, multiple frequency points are used to simulate the spectral
radiance before being averaged into a channel.  The bandwidth of each
channel is the same.
)",
      .author  = {"Richard Larsson"},
      .methods = {"measurement_sensorInit",
                  "measurement_sensorAddSimpleGaussian"},
      .out     = {"measurement_sensor"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "measurement_sensorVectorGaussian",
      .desc    = R"(Creates a single simple Gaussian-opening sensor.

This means that the sensor has a Gaussian bandwidth per channel.
That is, multiple frequency points are used to simulate the spectral
radiance before being averaged into a channel.  The bandwidth of each
channel is independent.
)",
      .author  = {"Richard Larsson"},
      .methods = {"measurement_sensorInit",
                  "measurement_sensorAddVectorGaussian"},
      .out     = {"measurement_sensor"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name   = "disort_spectral_flux_fieldFromAgenda",
      .desc   = R"(Use Disort for clearsky calculations of spectral flux field.

The agenda is used to setup Disort, i.e., to compute the *disort_settings*
that governs how the solver is run.
)",
      .author = {"Richard Larsson"},
      .methods = {"disort_settings_agendaExecute",
                  "disort_spectral_flux_fieldCalc"},
      .out     = {"disort_spectral_flux_field"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "disort_spectral_flux_fieldProfile",
      .desc =
          R"(Extract a 1D path through the atmosphere and calculate spectral flux using Disort.

This wrapper helps setting up a downlooking ray path through the atmosphere to form the
basis for the agenda to setup the Disort calculations.
)",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pathGeometricDownlooking",
                  "disort_spectral_flux_fieldFromAgenda"},
      .out     = {"disort_spectral_flux_field", "ray_path"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "disort_spectral_rad_fieldFromAgenda",
      .desc =
          R"(Use Disort for clearsky calculations of spectral radiance field.

The agenda is used to setup Disort, i.e., to compute the *disort_settings*
that governs how the solver is run.
)",
      .author  = {"Richard Larsson"},
      .methods = {"disort_settings_agendaExecute",
                  "disort_spectral_rad_fieldCalc"},
      .out     = {"disort_spectral_rad_field", "disort_quadrature"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "disort_spectral_rad_fieldProfile",
      .desc =
          R"(Extract a 1D path through the atmosphere and calculate spectral radiance using Disort.

This wrapper helps setting up a downlooking ray path through the atmosphere to form the
basis for the agenda to setup the Disort calculations.
)",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pathGeometricDownlooking",
                  "disort_spectral_rad_fieldFromAgenda"},
      .out     = {"disort_spectral_rad_field", "disort_quadrature", "ray_path"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "disort_spectral_rad_fieldDepthProfile",
      .desc =
          R"(Sets a ray path from a point and depth profile and calculates spectral radiance using Disort.

This wrapper helps setting up a downlooking ray path to form the
basis for the agenda to setup the Disort calculations.
)",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pathFromPointAndDepth",
                  "disort_spectral_rad_fieldFromAgenda"},
      .out     = {"disort_spectral_rad_field", "disort_quadrature", "ray_path"},
  });

#ifdef ENABLE_CDISORT
  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name   = "disort_spectral_rad_fieldFromAgendaCdisort",
      .desc   = "Use the disort settings agenda to calculate spectral radiance",
      .author = {"Oliver Lemke"},
      .methods = {"disort_settings_agendaExecute",
                  "atm_pathFromPath",
                  "freq_grid_pathFromPath",
                  "disort_spectral_rad_fieldCalcCdisort"},
      .out     = {"disort_spectral_rad_field", "disort_quadrature"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "disort_spectral_rad_fieldProfileCdisort",
      .desc =
          "Extract a 1D path through the atmospheric field and calculate spectral radiance using Disort",
      .author  = {"Oliver Lemke"},
      .methods = {"ray_pathGeometricDownlooking",
                  "disort_spectral_rad_fieldFromAgendaCdisort"},
      .out     = {"disort_spectral_rad_field", "disort_quadrature", "ray_path"},
  });
#endif

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radApplyUnitFromSpectralRadiance",
      .desc    = R"(Helper method for calling *spectral_radApplyUnit*.

It is common that *ray_path* is defined but not *ray_point*.
This method simply is a convenience wrapper for that use case.
)",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pointForeground", "spectral_radApplyUnit"},
      .out     = {"spectral_rad", "spectral_rad_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radClearskyEmission",
      .desc    = "Computes clearsky emission of spectral radiances",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pointBackground",
                  "spectral_rad_bkgAgendasAtEndOfPath",
                  "atm_pathFromPath",
                  "freq_grid_pathFromPath",
                  "spectral_propmat_pathFromPath",
                  "spectral_tramat_pathFromPath",
                  "spectral_tramat_cumulative_pathFromPath",
                  "spectral_rad_srcvec_pathFromPropmat",
                  "spectral_tramat_bkgFromPathPropagationBack",
                  "spectral_radStepByStepEmission",
                  "spectral_rad_jacFromBackground",
                  "spectral_rad_jacAddPathPropagation"},
      .out     = {"spectral_rad", "spectral_rad_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radClearskyLinearInTauEmission",
      .desc    = "[WIP] Computes clearsky emission of spectral radiances using linear-in-tau approximation",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pointBackground",
                  "spectral_rad_bkgAgendasAtEndOfPath",
                  "atm_pathFromPath",
                  "freq_grid_pathFromPath",
                  "spectral_propmat_pathFromPath",
                  "spectral_tramat_pathLinearInTauFromPath",
                  "spectral_tramat_cumulative_pathFromPath",
                  "spectral_rad_srcvec_pathFromPropmat",
                  "spectral_tramat_bkgFromPathPropagationBack",
                  "spectral_radLinearEvolutionStepByStepEmission",
                  "spectral_rad_jacFromBackground",
                  "spectral_rad_jacAddPathPropagation"},
      .out     = {"spectral_rad", "spectral_rad_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radClearskyLinearInTauAndPropEmission",
      .desc    = "[WIP] Computes clearsky emission of spectral radiances using linear-in-tau and linear-in-propmat approximation",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pointBackground",
                  "spectral_rad_bkgAgendasAtEndOfPath",
                  "atm_pathFromPath",
                  "freq_grid_pathFromPath",
                  "spectral_propmat_pathFromPath",
                  "spectral_tramat_pathLinearInTauAndPropFromPath",
                  "spectral_tramat_cumulative_pathFromPath",
                  "spectral_rad_srcvec_pathFromPropmat",
                  "spectral_tramat_bkgFromPathPropagationBack",
                  "spectral_radLinearEvolutionStepByStepEmission",
                  "spectral_rad_jacFromBackground",
                  "spectral_rad_jacAddPathPropagation"},
      .out     = {"spectral_rad", "spectral_rad_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radClearskyEmissionParFreq",
      .desc    = "Computes clearsky emission of spectral radiances",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pointBackground",
                  "spectral_rad_bkgAgendasAtEndOfPath",
                  "atm_pathFromPath",
                  "freq_grid_pathFromPath",
                  "spectral_propmat_pathFromPath",
                  "spectral_radSetToBackground",
                  "spectral_radSinglePathEmissionFrequencyLoop"},
      .out     = {"spectral_rad", "spectral_rad_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "spectral_radClearskyRayleighScattering",
      .desc =
          "Computes clearsky emission of spectral radiances with solar Rayleigh scattering",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pointBackground",
                  "spectral_rad_bkgAgendasAtEndOfPath",
                  "atm_pathFromPath",
                  "freq_grid_pathFromPath",
                  "spectral_propmat_pathFromPath",
                  "spectral_propmat_scat_pathFromPath",
                  "spectral_propmat_pathAddScattering",
                  "spectral_tramat_pathFromPath",
                  "spectral_tramat_cumulative_pathFromPath",
                  "spectral_rad_srcvec_pathFromPropmat",
                  "spectral_rad_scat_pathSunsFirstOrderRayleigh",
                  "spectral_rad_srcvec_pathAddScattering",
                  "spectral_tramat_bkgFromPathPropagationBack",
                  "spectral_radStepByStepEmission",
                  "spectral_rad_jacFromBackground",
                  "spectral_rad_jacAddPathPropagation"},
      .out     = {"spectral_rad", "spectral_rad_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radClearskyTransmission",
      .desc    = "Computes clearsky transmission of spectral radiances",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pointBackground",
                  "spectral_rad_bkgAgendasAtEndOfPath",
                  "atm_pathFromPath",
                  "freq_grid_pathFromPath",
                  "spectral_propmat_pathFromPath",
                  "spectral_tramat_pathFromPath",
                  "spectral_tramat_cumulative_pathFromPath",
                  "spectral_tramat_bkgFromPathPropagationBack",
                  "spectral_radCumulativeTransmission",
                  "spectral_rad_jacFromBackground",
                  "spectral_rad_jacAddPathPropagation"},
      .out     = {"spectral_rad", "spectral_rad_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_flux_profilePseudo2D",
      .desc    = "Computes the spectral flux profile using pseudo-2D geometry",
      .author  = {"Richard Larsson"},
      .methods = {"zen_gridProfilePseudo2D",
                  "spectral_rad_fieldProfilePseudo2D",
                  "spectral_flux_profileFromSpectralRadianceField"},
      .out     = {"spectral_flux_profile"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radClearskyBackgroundTransmission",
      .desc    = "Computes clearsky transmission of spectral radiances",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pointBackground",
                  "atm_pathFromPath",
                  "freq_grid_pathFromPath",
                  "spectral_propmat_pathFromPath",
                  "spectral_tramat_pathFromPath",
                  "spectral_tramat_cumulative_pathFromPath",
                  "spectral_tramat_bkgFromPathPropagationBack",
                  "spectral_radCumulativeTransmission",
                  "spectral_rad_jacFromBackground",
                  "spectral_rad_jacAddPathPropagation"},
      .out     = {"spectral_rad", "spectral_rad_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name             = "atm_fieldRead",
      .desc             = "Read atmospheric data files from a directory",
      .author           = {"Richard Larsson"},
      .methods          = {"atm_fieldInit",
                           "atm_fieldAppendBaseData",
                           "atm_fieldAppendAuto"},
      .out              = {"atm_field"},
      .preset_gin       = {"replace_existing"},
      .preset_gin_value = {Index{0}},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "spectral_radSubsurfaceDisortEmission",
      .desc =
          "Get the spectral radiance from subsurface emission simulated using Disort",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pathFromPointAndDepth",
                  "disort_settings_downwelling_wrapper_agendaExecute",
                  "disort_spectral_rad_fieldCalc",
                  "spectral_radFromDisort"},
      .out     = {"spectral_rad",
                  "disort_settings",
                  "ray_path",
                  "disort_spectral_rad_field",
                  "disort_quadrature"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "atm_fieldFitNonLTE",
      .desc    = "Fits non-LTE atmospheric field values",
      .author  = {"Richard Larsson"},
      .methods = {"atm_profileFromGrid",
                  "atm_profileFitNonLTE",
                  "atm_fieldFromProfile"},
      .out     = {"atm_field"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "UpdateModelStates",
      .desc =
          "Update state of the model in preparation for a forward model run",
      .author  = {"Richard Larsson"},
      .methods = {"abs_bandsFromModelState",
                  "surf_fieldFromModelState",
                  "subsurf_fieldFromModelState",
                  "atm_fieldFromModelState",
                  "measurement_sensorFromModelState"},
      .out     = {"abs_bands",
                  "surf_field",
                  "subsurf_field",
                  "atm_field",
                  "measurement_sensor"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "model_state_vecFromData",
      .desc    = "Get *model_state_vec* from available data",
      .author  = {"Richard Larsson"},
      .methods = {"model_state_vecInit",
                  "model_state_vecFromAtmosphere",
                  "model_state_vecFromSurface",
                  "model_state_vecFromSubsurface",
                  "model_state_vecFromBands",
                  "model_state_vecFromSensor"},
      .out     = {"model_state_vec"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "measurement_jacTransformations",
      .desc =
          "Apply all transformations to the Jacobian related to states in *model_state_vecFromData*",
      .author  = {"Richard Larsson"},
      .methods = {"measurement_jacAtmosphereTransformation",
                  "measurement_jacSurfaceTransformation",
                  "measurement_jacSubsurfaceTransformation",
                  "measurement_jacBandTransformation",
                  "measurement_jacSensorTransformation"},
      .out     = {"measurement_jac"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "model_state_vec_aprioriFromData",
      .desc    = "Get *model_state_vec_apriori* from available data",
      .author  = {"Richard Larsson"},
      .methods = {"model_state_vecFromData",
                  "model_state_vec_aprioriFromState"},
      .out     = {"model_state_vec_apriori"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "abs_lookup_dataCalc",
      .desc    = R"(Get *abs_lookup_data* from available data.

This method will use the *atm_field* and *abs_bands* to
calculate the *abs_lookup_data*.  The atmospheric field is first
gridded using *atm_profileExtract*.
)",
      .author  = {"Richard Larsson"},
      .methods = {"atm_profileExtract",
                  "abs_lookup_dataInit",
                  "abs_lookup_dataPrecomputeAll"},
      .out     = {"abs_lookup_data"},
  });

  return wsm_meta;
}
}  // namespace

const std::vector<WorkspaceMethodInternalMetaRecord>& internal_meta_methods() {
  static const auto out = internal_meta_methods_creator();
  return out;
}

WorkspaceMethodInternalRecord WorkspaceMethodInternalMetaRecord::create(
    const std::unordered_map<std::string, WorkspaceMethodInternalRecord>& wsms)
    const try {
  if (wsms.contains(name)) {
    throw std::runtime_error("Meta-function " + name + " already exists");
  }

  WorkspaceMethodInternalRecord wsm{
      .desc   = desc,
      .author = author,
      .out    = out,
  };

  std::vector<std::string> first_out;
  std::vector<std::string> first_inout;

  wsm.desc += "\n\nWrapper calling Methods (in order):\n\n";
  for (const auto& m : methods) {
    wsm.desc += "  - *" + m + "*\n";

    const auto ptr = wsms.find(m);
    if (ptr == wsms.end()) {
      throw std::runtime_error(std::format(R"(Method "{}" not found)", m));
    }

    const auto& wm = ptr->second;

    if (wm.has_any() or wm.has_overloads()) {
      throw std::runtime_error(std::format(
          R"(Method "{}"  has overloads and does not work with meta-functions)",
          m));
    }

    wsm.author.insert(wsm.author.end(), wm.author.begin(), wm.author.end());
    wsm.gout.insert(wsm.gout.end(), wm.gout.begin(), wm.gout.end());
    wsm.gout_type.insert(
        wsm.gout_type.end(), wm.gout_type.begin(), wm.gout_type.end());
    wsm.gout_desc.insert(
        wsm.gout_desc.end(), wm.gout_desc.begin(), wm.gout_desc.end());
    wsm.gin.insert(wsm.gin.end(), wm.gin.begin(), wm.gin.end());
    wsm.gin_type.insert(
        wsm.gin_type.end(), wm.gin_type.begin(), wm.gin_type.end());
    wsm.gin_value.insert(
        wsm.gin_value.end(), wm.gin_value.begin(), wm.gin_value.end());
    wsm.gin_desc.insert(
        wsm.gin_desc.end(), wm.gin_desc.begin(), wm.gin_desc.end());
    wsm.pass_workspace = wsm.pass_workspace or wm.pass_workspace;

    stdr::copy_if(wm.in, std::back_inserter(wsm.in), [&](const std::string& i) {
      const auto cmp = Cmp::eq(i);
      return stdr::none_of(wsm.in, cmp) and stdr::none_of(first_out, cmp) and
             stdr::none_of(first_inout, cmp) and not stdr::any_of(wm.out, cmp);
    });

    stdr::copy_if(
        wm.out, std::back_inserter(first_out), [&](const std::string& o) {
          const auto cmp = Cmp::eq(o);
          return stdr::none_of(wsm.in, cmp) and
                 stdr::none_of(first_out, cmp) and
                 stdr::none_of(first_inout, cmp) and
                 not stdr::any_of(wm.in, cmp);
        });

    stdr::copy_if(
        wm.out, std::back_inserter(first_inout), [&](const std::string& o) {
          const auto cmp = Cmp::eq(o);
          return stdr::none_of(wsm.in, cmp) and
                 stdr::none_of(first_out, cmp) and
                 stdr::none_of(first_inout, cmp) and stdr::any_of(wm.in, cmp);
        });
  }

  wsm.desc += R"(
Equivalent (mostly) Python code:

.. code-block:: python3
  :linenos:

  ws = pyarts.Workspace()

  # ...

)";
  for (const auto& m : methods) {
    wsm.desc += "   ws." + m + "()\n";
  }
  wsm.desc += "\n";

  stdr::sort(wsm.author);
  wsm.author.erase(std::unique(wsm.author.begin(), wsm.author.end()),
                   wsm.author.end());

  if (not wsm.gout.empty()) {
    throw std::runtime_error(
        "Output generic variables not supported in meta-functions");
  }

  bubble_sort_by([&wsm](Size i, Size j) { return wsm.gin[i] < wsm.gin[j]; },
                 wsm.gin,
                 wsm.gin_type,
                 wsm.gin_value,
                 wsm.gin_desc);

  for (auto ptr = stdr::adjacent_find(wsm.gin); ptr != wsm.gin.end();
       ptr      = stdr::adjacent_find(wsm.gin)) {
    const auto idx0 = std::distance(wsm.gin.begin(), ptr);
    const auto idx1 = idx0 + 1;

    if (wsm.gin_type[idx0] != wsm.gin_type[idx1]) {
      throw std::runtime_error("Incompatible types for generic input " +
                               wsm.gin[idx0]);
    }

    if (wsm.gin_desc[idx0] != wsm.gin_desc[idx1]) {
      throw std::runtime_error("Incompatible descriptions for generic input " +
                               wsm.gin[idx0]);
    }

    //! Keep preferably nullopt value
    wsm.gin_value.erase(wsm.gin_value.begin() +
                        (wsm.gin_value[idx0] == std::nullopt ? idx0 : idx1));

    wsm.gin.erase(ptr);
    wsm.gin_type.erase(wsm.gin_type.begin() + idx0);
    wsm.gin_desc.erase(wsm.gin_desc.begin() + idx0);
  }

  ARTS_USER_ERROR_IF(preset_gin.size() != preset_gin_value.size(),
                     "preset_gin and preset_gin_value must have the same size")
  for (auto& preset : preset_gin) {
    for (auto ptr = stdr::find(wsm.gin, preset); ptr != wsm.gin.end();
         ptr      = stdr::adjacent_find(wsm.gin)) {
      const auto idx0 = std::distance(wsm.gin.begin(), ptr);

      wsm.gin.erase(ptr);
      wsm.gin_value.erase(wsm.gin_value.begin() + idx0);
      wsm.gin_type.erase(wsm.gin_type.begin() + idx0);
      wsm.gin_desc.erase(wsm.gin_desc.begin() + idx0);
    }
  }

  stdr::sort(first_out);
  stdr::sort(first_inout);
  stdr::sort(wsm.in);

  for (Size i = 0; i < out.size(); i++) {
    const auto& o = out[i];

    if (stdr::binary_search(first_inout, o)) {
      wsm.in.insert(wsm.in.begin() + i, o);
    }
  }

  for (auto& o : first_inout) {
    if (stdr::none_of(out, Cmp::eq(o))) {
      throw std::runtime_error("Input output variable " + o +
                               " is not declared as output");
    }
  }

  return wsm;
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error creating meta-function \"{}\":\n\n{}",
                  name,
                  std::string_view(e.what())));
}

std::string WorkspaceMethodInternalMetaRecord::call(
    const std::unordered_map<std::string, WorkspaceMethodInternalRecord>& wsms)
    const try {
  const auto& wsvs = internal_workspace_variables();
  const auto& wsas = internal_workspace_agendas();

  const auto ptr = wsms.find(name);
  if (ptr == wsms.end()) {
    throw std::runtime_error("Meta-function " + name + " not found");
  }
  const auto& wsm = ptr->second;
  std::vector<std::string> first_out;

  std::stringstream code;

  code << wsm.header(name, 0) << " try {\n  ARTS_TIME_REPORT\n\n";

  for (Size i = 0; i < preset_gin.size(); i++) {
    const auto t = preset_gin_value[i].type_name();
    std::println(code,
                 R"(  static const {0}& {1} =
    stdr::find_if(_wsmmeta, [](auto& _wsmmeta_var) {{
      return _wsmmeta_var.name=="{2}";
    }}) -> preset_gin_value[{3}].get_unsafe<{0}>();)",
                 t,
                 preset_gin[i],
                 name,
                 i);
  }

  for (const auto& m : methods) {
    const auto& wm = wsms.at(m);
    for (const auto& outvar : wm.out) {
      const auto cmp = Cmp::eq(outvar);
      if (stdr::any_of(first_out, cmp) or stdr::any_of(wsm.in, cmp) or
          stdr::any_of(out, cmp)) {
        continue;
      }
      first_out.push_back(outvar);

      auto wsvptr = wsvs.find(outvar);
      if (wsvptr != wsvs.end()) {
        code << "  " << wsvptr->second.type << " " << outvar << ";\n";
        continue;
      }

      auto wsaptr = wsas.find(outvar);
      if (wsaptr != wsas.end()) {
        code << "  Agenda " << outvar << ";\n";
        continue;
      }

      throw std::runtime_error("Output " + outvar +
                               " not found as Agenda or Variable");
    }

    code << "  " << wm.call(m) << '\n';
  }
  code << "} ARTS_METHOD_ERROR_CATCH";

  return code.str();
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error creating call for meta-function \"{}\":\n\n{}",
                  name,
                  std::string_view(e.what())));
}