#include <enums.h>
#include <workspace.h>

std::string& docstr(std::unordered_map<std::string, std::string>& map,
                    const std::string& key) {
  const auto& orig = workspace_methods();
  if (orig.find(key) == orig.end()) {
    throw std::runtime_error(
        std::format(R"(No method "{}" in ARTS methods)", key));
  }
  return map[key];
}

std::string replace(std::string str,
                    const std::string& from,
                    const std::string& to) {
  size_t start_pos = str.find(from);
  if (start_pos == std::string::npos) return str;
  str.replace(start_pos, from.length(), to);
  return str;
}

std::string get_disort_settings_agendaSetup_doc() {
  const std::string method_name = "disort_settings_agendaSetup";

  auto& method = workspace_methods().at(method_name);

  std::string doc = std::format(
      R"(
There are {0} possible combinations for calling the ``{1}`` method.

Below, these are all listed with the generated agenda-call order for each combination in full.

Before that, a concise overview of what each option do is available by the types in this table:

.. list-table::
  :name: Setup options for ``{1}``
  :widths: auto
  :align: left
  :header-rows: 1

  * - Input variable
    - pyarts class
  * - ``{2}``
    - :class:`~pyarts.arts.disort_settings_agenda_setup_layer_emission_type`
  * - ``{3}``
    -  :class:`~pyarts.arts.disort_settings_agenda_setup_scattering_type`
  * - ``{4}``
    -  :class:`~pyarts.arts.disort_settings_agenda_setup_space_type` 
  * - ``{5}``
    -  :class:`~pyarts.arts.disort_settings_agenda_setup_sun_type`
  * - ``{6}``
    -  :class:`~pyarts.arts.disort_settings_agenda_setup_surface_type`
)",
      enumstrs::disort_settings_agenda_setup_layer_emission_typeNames<>.size() *
          enumstrs::disort_settings_agenda_setup_scattering_typeNames<>.size() *
          enumstrs::disort_settings_agenda_setup_space_typeNames<>.size() *
          enumstrs::disort_settings_agenda_setup_sun_typeNames<>.size() *
          enumstrs::disort_settings_agenda_setup_surface_typeNames<>.size(),
      method_name,
      
      method.in[0].substr(1),
      method.in[1].substr(1),
      method.in[2].substr(1),
      method.in[3].substr(1),
      method.in[4].substr(1));

  for (auto layer_emission_setting :
       enumstrs::disort_settings_agenda_setup_layer_emission_typeNames<>) {
    for (auto scattering_setting :
         enumstrs::disort_settings_agenda_setup_scattering_typeNames<>) {
      for (auto space_setting :
           enumstrs::disort_settings_agenda_setup_space_typeNames<>) {
        for (auto sun_setting :
             enumstrs::disort_settings_agenda_setup_sun_typeNames<>) {
          for (auto surface_setting :
               enumstrs::disort_settings_agenda_setup_surface_typeNames<>) {
            Agenda x{};
            Vector y{-99, 0, 99};
            std::string docstr{};
            try {
              disort_settings_agendaSetup(x,
                                          std::string{layer_emission_setting},
                                          std::string{scattering_setting},
                                          std::string{space_setting},
                                          std::string{sun_setting},
                                          std::string{surface_setting},
                                          y);
              docstr = replace(
                  x.sphinx_list(), "[-99, 0, 99]", "lambertian_reflection");
            } catch (std::exception& e) {
              docstr = std::format("Invalid combination: {}", e.what());
            }

            doc += std::format(
                R"(
------------------------------------------------------------

``{}({}="{}", {}="{}", {}="{}", {}="{}", {}="{}", {}=lambertian_reflection)``

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
                docstr);
          }
        }
      }
    }
  }

  return doc;
}

std::unordered_map<std::string, std::string>
workspace_method_extra_doc_internal() {
  std::unordered_map<std::string, std::string> doc{};

  docstr(doc, "disort_settings_agendaSetup") =
      get_disort_settings_agendaSetup_doc();

  return doc;
}

const std::unordered_map<std::string, std::string>&
workspace_method_extra_doc() {
  static auto doc = workspace_method_extra_doc_internal();
  return doc;
}
