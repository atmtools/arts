#include "workspace_method_extra_doc.h"

#include <configtypes.h>
#include <enums.h>
#include <enumsSpeciesEnum.h>
#include <isotopologues.h>
#include <workspace.h>

#include <ranges>

namespace {
std::string& docstr(std::unordered_map<std::string, std::string>& map, const std::string& key) {
  const auto& orig = workspace_methods();
  if (orig.find(key) == orig.end()) { throw std::runtime_error(std::format(R"(No method "{}" in ARTS methods)", key)); }
  return map[key];
}

std::string replace(std::string str, const std::string& from, const std::string& to) {
  size_t start_pos = str.find(from);
  if (start_pos == std::string::npos) return str;
  str.replace(start_pos, from.length(), to);
  return str;
}

/* The documentation of an agenda setup method of DISORT (disort_settings_agendaSetup) or VDISORT
   (vdisort_settings_agendaSetup): every combination of the five options that the method accepts, with the agenda
   it makes.  A combination that the method refuses is not listed. */
template <typename Setup> std::string get_agendaSetup_doc(const std::string& method_name, Setup&& setup) {
  auto& method = workspace_methods().at(method_name);

  Size        count = 0;
  std::string combinations{};
  for (auto layer_emission_setting : enumstrs::settings_agenda_setup_layer_emission_typeNames<>) {
    for (auto scattering_setting : enumstrs::settings_agenda_setup_scattering_typeNames<>) {
      for (auto space_setting : enumstrs::settings_agenda_setup_space_typeNames<>) {
        for (auto sun_setting : enumstrs::settings_agenda_setup_sun_typeNames<>) {
          for (auto surface_setting : enumstrs::settings_agenda_setup_surface_typeNames<>) {
            Agenda x{};
            Vector y{-99, 0, 99};
            try {
              setup(x,
                    std::string{layer_emission_setting},
                    std::string{scattering_setting},
                    std::string{space_setting},
                    std::string{sun_setting},
                    std::string{surface_setting},
                    y,
                    123);
            } catch (std::exception&) { continue; }

            count++;
            combinations += std::format(
                R"(
.. rubric:: ``{}({}="{}", {}="{}", {}="{}", {}="{}", {}="{}", {}=lambertian_reflection, min_optical_depth=min_optical_depth)``

{}
)",
                method_name,
                method.in[0].substr(1),
                layer_emission_setting,
                method.in[1].substr(1),
                scattering_setting,
                method.in[2].substr(1),
                space_setting,
                method.in[3].substr(1),
                sun_setting,
                method.in[4].substr(1),
                surface_setting,
                method.in[5].substr(1),
                replace(replace(x.sphinx_list(), "[-99, 0, 99]", "lambertian_reflection"), "123", "min_optical_depth"));
          }
        }
      }
    }
  }

  return std::format(
             R"(
There are {0} possible combinations for calling the ``{1}`` method.

Below, these are all listed with the generated agenda-call order for each combination in full.

Before that, a concise overview of what each option does is available by the types in this table:

.. list-table::
  :widths: auto
  :align: left
  :header-rows: 1

  * - Input variable
    - pyarts class
  * - ``{2}``
    - :class:`~pyarts3.arts.settings_agenda_setup_layer_emission_type`
  * - ``{3}``
    -  :class:`~pyarts3.arts.settings_agenda_setup_scattering_type`
  * - ``{4}``
    -  :class:`~pyarts3.arts.settings_agenda_setup_space_type`
  * - ``{5}``
    -  :class:`~pyarts3.arts.settings_agenda_setup_sun_type`
  * - ``{6}``
    -  :class:`~pyarts3.arts.settings_agenda_setup_surface_type`
)",
             count,
             method_name,
             method.in[0].substr(1),
             method.in[1].substr(1),
             method.in[2].substr(1),
             method.in[3].substr(1),
             method.in[4].substr(1)) +
         combinations;
}

std::string get_disort_settings_agendaSetup_doc() {
  return get_agendaSetup_doc("disort_settings_agendaSetup",
                             [](auto&&... args) { disort_settings_agendaSetup(args...); });
}

std::string get_vdisort_settings_agendaSetup_doc() {
  return get_agendaSetup_doc("vdisort_settings_agendaSetup",
                             [](auto&&... args) { vdisort_settings_agendaSetup(args...); });
}

std::string get_abs_speciesSet_doc() {
  const std::string method_name = "abs_speciesSet";

  std::string doc = R"(
Below follows a complete list of all single species tags that can be set using this method.

.. list-table::
   :widths: auto
   :align: left
   :header-rows: 1

   * - Isotopologue name
     - Type
     - Mass
     - :math:`g_i`)";

  for (auto& x : Species::Isotopologues | stdv::drop(1)) {
    doc += std::format(R"(
   * - ``{}``
     - {}
     - {}
     - {})",
                       x.FullName(),
                       x.is_predefined() ? "Predefined Model"sv
                       : x.is_normal()   ? "Normal Isotopologue"sv
                                         : "All Isotopologues"sv,
                       x.is_normal() ? std::format("{}", x.mass) : ""s,
                       x.is_normal() ? std::format("{}", x.gi) : ""s);
  }

  doc += "\n\nAnd in short: ";
  for (auto& x : Species::Isotopologues | stdv::drop(1)) { doc += std::format("``{}``, ", x.FullName()); }
  doc += '\n';

  return doc;
}

std::unordered_map<std::string, std::string> workspace_method_extra_doc_internal() {
  std::unordered_map<std::string, std::string> doc{};

  docstr(doc, "disort_settings_agendaSetup") = get_disort_settings_agendaSetup_doc();

  docstr(doc, "vdisort_settings_agendaSetup") = get_vdisort_settings_agendaSetup_doc();

  docstr(doc, "abs_speciesSet") = get_abs_speciesSet_doc();

  return doc;
}
}  // namespace

const std::unordered_map<std::string, std::string>& workspace_method_extra_doc() {
  static auto doc = workspace_method_extra_doc_internal();
  return doc;
}
