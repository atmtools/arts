#include "workspace_meta_methods.h"

#include <sorting.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "workspace_agendas.h"
#include "workspace_methods.h"
#include "workspace_variables.h"

std::vector<WorkspaceMethodInternalMetaRecord> internal_meta_methods_creator() {
  std::vector<WorkspaceMethodInternalMetaRecord> out;

  out.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "spectral_radianceStandardEmission",
      .desc = "Computes standard emission of spectral radiances",
      .author = {"Richard Larsson"},
      .methods = {"spectral_radiance_backgroundAgendasAtEndOfPath",
                  "propagation_path_atmospheric_pointFromPath",
                  "propagation_path_frequency_gridFromPath",
                  "propagation_path_propagation_matrixFromPath",
                  "propagation_path_transmission_matrixFromPath",
                  "propagation_path_transmission_matrix_cumulativeForward",
                  "propagation_path_spectral_radiance_sourceFromPropmat",
                  "propagation_path_spectral_radianceCalcEmission",
                  "background_transmittanceFromPathPropagationBack",
                  "spectral_radianceFromPathPropagation",
                  "spectral_radiance_jacobianFromBackground",
                  "spectral_radiance_jacobianAddPathPropagation"},
      .out = {"spectral_radiance", "spectral_radiance_jacobian"},
  });

  out.push_back(WorkspaceMethodInternalMetaRecord{
      .name = "atmospheric_fieldReadCatalog",
      .desc = "Reads absorption file from a directory",
      .author = {"Richard Larsson"},
      .methods = {"atmospheric_fieldInit",
                  "atmospheric_fieldAppendBaseData",
                  "atmospheric_fieldAppendAbsorptionData"},
      .out = {"atmospheric_field"},
      .preset_gin = {"replace_existing"},
      .preset_gin_value = {Index{0}},
  });

  return out;
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
      .desc = desc,
      .author = author,
      .out = out,
  };

  std::vector<std::string> first_out;
  std::vector<std::string> first_inout;

  wsm.desc += "\n\nWrapper calling Methods (in order):\n\n";
  for (const auto& m : methods) {
    wsm.desc += "  - *" + m + "*\n";

    const auto& wm = wsms.at(m);

    if (wm.has_any() or wm.has_overloads()) {
      throw std::runtime_error(
          "Method " + m +
          " has overloads and does not work with meta-functions");
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
       ptr = std::ranges::adjacent_find(wsm.gin)) {
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
         ptr = std::ranges::adjacent_find(wsm.gin)) {
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
    if (not std::ranges::binary_search(first_out, o) and
        not std::ranges::binary_search(first_inout, o)) {
      throw std::runtime_error("Output " + o + " not found in any method");
    }

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
  throw std::runtime_error(var_string(
      "Error creating meta-function ", std::quoted(name), ":\n\n", e.what()));
}

std::string WorkspaceMethodInternalMetaRecord::call(
    const std::unordered_map<std::string, WorkspaceMethodInternalRecord>& wsms)
    const try {
  const static auto wsvs = internal_workspace_variables();
  const static auto wsas = internal_workspace_agendas();

  const auto& wsm = wsms.at(name);
  std::vector<std::string> first_out;

  std::stringstream code;

  code << wsm.header(name, 0) << " {\n";

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
        if (wsaptr->second.array)
          code << "  ArrayOfAgenda " << outvar << ";\n";
        else
          code << "  Agenda " << outvar << ";\n";
        continue;
      }

      throw std::runtime_error("Output " + outvar +
                               " not found as Agenda or Variable");
    }

    code << "  " << wm.call(m) << '\n';
  }
  code << '}';

  return code.str();
} catch (std::exception& e) {
  throw std::runtime_error(var_string("Error creating call for meta-function ",
                                      std::quoted(name),
                                      ":\n\n",
                                      e.what()));
}