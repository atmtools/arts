#include "workspace_meta_methods.h"

#include <sorting.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "workspace_agendas.h"
#include "workspace_methods.h"
#include "workspace_variables.h"

std::vector<WorkspaceMethodInternalMetaRecord> internal_meta_methods_creator() {
  std::vector<WorkspaceMethodInternalMetaRecord> wsm_meta;

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "measurement_sensorSimple",
      .desc    = "Wrapper for a single simple dirac-opening sensor",
      .author  = {"Richard Larsson"},
      .methods = {"measurement_sensorInit", "measurement_sensorAddSimple"},
      .out     = {"measurement_sensor"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "measurement_sensorSimpleGaussian",
      .desc    = "Wrapper for a single simple Gaussian-opening sensor",
      .author  = {"Richard Larsson"},
      .methods = {"measurement_sensorInit",
                  "measurement_sensorAddSimpleGaussian"},
      .out     = {"measurement_sensor"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "measurement_sensorVectorGaussian",
      .desc    = "Wrapper for a single simple Gaussian-opening sensor",
      .author  = {"Richard Larsson"},
      .methods = {"measurement_sensorInit",
                  "measurement_sensorAddVectorGaussian"},
      .out     = {"measurement_sensor"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "disort_spectral_flux_fieldFromAgenda",
      .desc    = "Use Disort for clearsky calculations of spectral flux field",
      .author  = {"Richard Larsson"},
      .methods = {"disort_settings_agendaExecute",
                  "disort_spectral_flux_fieldCalc"},
      .out     = {"disort_spectral_flux_field"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "disort_spectral_flux_fieldProfile",
      .desc =
          "Extract a 1D path through the atmospheric field and calculate spectral flux using Disort",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pathGeometricDownlooking",
                  "disort_spectral_flux_fieldFromAgenda"},
      .out     = {"disort_spectral_flux_field", "ray_path"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name   = "disort_spectral_radiance_fieldFromAgenda",
      .desc   = "Use the disort settings agenda to calculate spectral radiance",
      .author = {"Richard Larsson"},
      .methods = {"disort_settings_agendaExecute",
                  "disort_spectral_radiance_fieldCalc"},
      .out     = {"disort_spectral_radiance_field",
                  "disort_quadrature_angles",
                  "disort_quadrature_weights"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "disort_spectral_radiance_fieldProfile",
      .desc =
          "Extract a 1D path through the atmospheric field and calculate spectral radiance using Disort",
      .author  = {"Richard Larsson"},
      .methods = {"ray_pathGeometricDownlooking",
                  "disort_spectral_radiance_fieldFromAgenda"},
      .out     = {"disort_spectral_radiance_field",
                  "disort_quadrature_angles",
                  "disort_quadrature_weights",
                  "ray_path"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radianceApplyUnitFromSpectralRadiance",
      .desc    = R"(Apply unit changes to spectral radiance and its Jacobian

.. warning::
  This is a destructive method.  Any use of it means that it is undefined behavior
  to use *spectral_radiance* or *spectral_radiance_jacobian* in future methods.
)",
      .author  = {"Richard Larsson"},
      .methods = {"ray_path_pointForeground",
                  "spectral_radiance_jacobianApplyUnit",
                  "spectral_radianceApplyUnit"},
      .out     = {"spectral_radiance", "spectral_radiance_jacobian"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radianceClearskyEmission",
      .desc    = "Computes clearsky emission of spectral radiances",
      .author  = {"Richard Larsson"},
      .methods = {"ray_path_pointBackground",
                  "spectral_radiance_backgroundAgendasAtEndOfPath",
                  "ray_path_atmospheric_pointFromPath",
                  "ray_path_frequency_gridFromPath",
                  "ray_path_propagation_matrixFromPath",
                  "ray_path_transmission_matrixFromPath",
                  "ray_path_transmission_matrix_cumulativeFromPath",
                  "ray_path_spectral_radiance_sourceFromPropmat",
                  "transmission_matrix_backgroundFromPathPropagationBack",
                  "spectral_radianceStepByStepEmission",
                  "spectral_radiance_jacobianFromBackground",
                  "spectral_radiance_jacobianAddPathPropagation"},
      .out     = {"spectral_radiance", "spectral_radiance_jacobian"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radianceClearskyRayleighScattering",
      .desc    = "Computes clearsky emission of spectral radiances",
      .author  = {"Richard Larsson"},
      .methods = {"ray_path_pointBackground",
                  "spectral_radiance_backgroundAgendasAtEndOfPath",
                  "ray_path_atmospheric_pointFromPath",
                  "ray_path_frequency_gridFromPath",
                  "ray_path_propagation_matrixFromPath",
                  "ray_path_propagation_matrix_scatteringFromPath",
                  "ray_path_propagation_matrixAddScattering",
                  "ray_path_transmission_matrixFromPath",
                  "ray_path_transmission_matrix_cumulativeFromPath",
                  "ray_path_spectral_radiance_sourceFromPropmat",
                  "ray_path_spectral_radiance_scatteringSunsFirstOrderRayleigh",
                  "ray_path_spectral_radiance_sourceAddScattering",
                  "transmission_matrix_backgroundFromPathPropagationBack",
                  "spectral_radianceStepByStepEmission",
                  "spectral_radiance_jacobianFromBackground",
                  "spectral_radiance_jacobianAddPathPropagation"},
      .out     = {"spectral_radiance", "spectral_radiance_jacobian"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radianceClearskyTransmission",
      .desc    = "Computes clearsky transmission of spectral radiances",
      .author  = {"Richard Larsson"},
      .methods = {"ray_path_pointBackground",
                  "spectral_radiance_backgroundAgendasAtEndOfPath",
                  "ray_path_atmospheric_pointFromPath",
                  "ray_path_frequency_gridFromPath",
                  "ray_path_propagation_matrixFromPath",
                  "ray_path_transmission_matrixFromPath",
                  "ray_path_transmission_matrix_cumulativeFromPath",
                  "transmission_matrix_backgroundFromPathPropagationBack",
                  "spectral_radianceCumulativeTransmission",
                  "spectral_radiance_jacobianFromBackground",
                  "spectral_radiance_jacobianAddPathPropagation"},
      .out     = {"spectral_radiance", "spectral_radiance_jacobian"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_flux_profilePseudo2D",
      .desc    = "Computes the spectral flux profile using pseudo-2D geometry",
      .author  = {"Richard Larsson"},
      .methods = {"zenith_gridProfilePseudo2D",
                  "spectral_radiance_fieldProfilePseudo2D",
                  "spectral_flux_profileFromSpectralRadianceField"},
      .out     = {"spectral_flux_profile"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "spectral_radianceClearskyBackgroundTransmission",
      .desc    = "Computes clearsky transmission of spectral radiances",
      .author  = {"Richard Larsson"},
      .methods = {"ray_path_pointBackground",
                  "ray_path_atmospheric_pointFromPath",
                  "ray_path_frequency_gridFromPath",
                  "ray_path_propagation_matrixFromPath",
                  "ray_path_transmission_matrixFromPath",
                  "ray_path_transmission_matrix_cumulativeFromPath",
                  "transmission_matrix_backgroundFromPathPropagationBack",
                  "spectral_radianceCumulativeTransmission",
                  "spectral_radiance_jacobianFromBackground",
                  "spectral_radiance_jacobianAddPathPropagation"},
      .out     = {"spectral_radiance", "spectral_radiance_jacobian"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name             = "atmospheric_fieldRead",
      .desc             = "Reads absorption file from a directory",
      .author           = {"Richard Larsson"},
      .methods          = {"atmospheric_fieldInit",
                           "atmospheric_fieldAppendBaseData",
                           "atmospheric_fieldAppendAbsorptionData"},
      .out              = {"atmospheric_field"},
      .preset_gin       = {"replace_existing"},
      .preset_gin_value = {Index{0}},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "atmospheric_fieldFitNonLTE",
      .desc    = "Fits ",
      .author  = {"Richard Larsson"},
      .methods = {"atmospheric_profileFromGrid",
                  "atmospheric_profileFitNonLTE",
                  "atmospheric_fieldFromProfile"},
      .out     = {"atmospheric_field"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "UpdateModelStates",
      .desc =
          "Update state of the model in preparation for a forward model run",
      .author  = {"Richard Larsson"},
      .methods = {"absorption_bandsFromModelState",
                  "surface_fieldFromModelState",
                  "atmospheric_fieldFromModelState",
                  "measurement_sensorFromModelState"},
      .out     = {"absorption_bands",
                  "surface_field",
                  "atmospheric_field",
                  "measurement_sensor"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "model_state_vectorFromData",
      .desc    = "Get *model_state_vector* from available data",
      .author  = {"Richard Larsson"},
      .methods = {"model_state_vectorSize",
                  "model_state_vectorZero",
                  "model_state_vectorFromAtmosphere",
                  "model_state_vectorFromSurface",
                  "model_state_vectorFromBands",
                  "model_state_vectorFromSensor"},
      .out     = {"model_state_vector"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "model_state_vector_aprioriFromData",
      .desc    = "Get *model_state_vector_apriori* from available data",
      .author  = {"Richard Larsson"},
      .methods = {"model_state_vectorFromData",
                  "model_state_vector_aprioriFromState"},
      .out     = {"model_state_vector_apriori"},
  });

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name    = "absorption_lookup_tableCalc",
      .desc    = R"(Get *absorption_lookup_table* from available data.

This method will use the *atmospheric_field* and *absorption_bands* to
calculate the *absorption_lookup_table*.  The atmospheric field is first
gridded using *atmospheric_profileExtract*.
)",
      .author  = {"Richard Larsson"},
      .methods = {"atmospheric_profileExtract",
                  "absorption_lookup_tableInit",
                  "absorption_lookup_tablePrecomputeAll"},
      .out     = {"absorption_lookup_table"},
  });

  return wsm_meta;
}
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

    std::ranges::copy_if(
        wm.in, std::back_inserter(wsm.in), [&](const std::string& i) {
          const auto cmp = Cmp::eq(i);
          return std::ranges::none_of(wsm.in, cmp) and
                 std::ranges::none_of(first_out, cmp) and
                 std::ranges::none_of(first_inout, cmp) and
                 not std::ranges::any_of(wm.out, cmp);
        });

    std::ranges::copy_if(
        wm.out, std::back_inserter(first_out), [&](const std::string& o) {
          const auto cmp = Cmp::eq(o);
          return std::ranges::none_of(wsm.in, cmp) and
                 std::ranges::none_of(first_out, cmp) and
                 std::ranges::none_of(first_inout, cmp) and
                 not std::ranges::any_of(wm.in, cmp);
        });

    std::ranges::copy_if(
        wm.out, std::back_inserter(first_inout), [&](const std::string& o) {
          const auto cmp = Cmp::eq(o);
          return std::ranges::none_of(wsm.in, cmp) and
                 std::ranges::none_of(first_out, cmp) and
                 std::ranges::none_of(first_inout, cmp) and
                 std::ranges::any_of(wm.in, cmp);
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

  std::ranges::sort(wsm.author);
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

  for (auto ptr = std::ranges::adjacent_find(wsm.gin); ptr != wsm.gin.end();
       ptr      = std::ranges::adjacent_find(wsm.gin)) {
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
    for (auto ptr = std::ranges::find(wsm.gin, preset); ptr != wsm.gin.end();
         ptr      = std::ranges::adjacent_find(wsm.gin)) {
      const auto idx0 = std::distance(wsm.gin.begin(), ptr);

      wsm.gin.erase(ptr);
      wsm.gin_value.erase(wsm.gin_value.begin() + idx0);
      wsm.gin_type.erase(wsm.gin_type.begin() + idx0);
      wsm.gin_desc.erase(wsm.gin_desc.begin() + idx0);
    }
  }

  std::ranges::sort(first_out);
  std::ranges::sort(first_inout);
  std::ranges::sort(wsm.in);

  for (Size i = 0; i < out.size(); i++) {
    const auto& o = out[i];

    if (std::ranges::binary_search(first_inout, o)) {
      wsm.in.insert(wsm.in.begin() + i, o);
    }
  }

  for (auto& o : first_inout) {
    if (std::ranges::none_of(out, Cmp::eq(o))) {
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
    code << "  static const " << t << "& " << preset_gin[i]
         << " =\n    "
            "std::ranges::find_if(_wsmmeta, [](auto& _wsmmeta_var) {"
            "\n      return _wsmmeta_var.name==\""
         << name
         << "\";\n    "
            "}) -> preset_gin_value["
         << i << "].get_unsafe<" << t << ">();\n";
  }

  for (const auto& m : methods) {
    const auto& wm = wsms.at(m);
    for (const auto& outvar : wm.out) {
      const auto cmp = Cmp::eq(outvar);
      if (std::ranges::any_of(first_out, cmp) or
          std::ranges::any_of(wsm.in, cmp) or std::ranges::any_of(out, cmp)) {
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