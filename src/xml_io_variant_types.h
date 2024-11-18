#pragma once

#include <workspace.h>
#include <xml_io.h>

#include <variant>

template <WorkspaceGroup... Types>
void xml_write(std::ostream& os_xml,
               const std::variant<Types...>& var,
               bofstream* pbofs,
               const String& name) {
  constexpr std::array<std::string_view, sizeof...(Types)> valtype{
      WorkspaceGroupInfo<Types>::name...};

  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Variant");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", String{valtype[var.index()]});

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  std::visit(
      [&os_xml, pbofs](const auto& value) {
        xml_write_to_stream(os_xml, value, pbofs, "");
      },
      var);

  close_tag.set_name("/Variant");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

template <Size j = 0, WorkspaceGroup... Types>
std::variant<Types...> init(Size i) {
  if (i == j) return std::variant_alternative_t<j, std::variant<Types...>>{};

  if constexpr (j + 1 < sizeof...(Types)) {
    return init<j + 1, Types...>(i);
  } else {
    throw std::runtime_error("Cannot understand the variant type");
  }
}

template <WorkspaceGroup... Types>
void xml_read(std::istream& is_xml,
              std::variant<Types...>& var,
              bifstream* pbifs) {
  constexpr std::array<std::string_view, sizeof...(Types)> valtype{
      WorkspaceGroupInfo<Types>::name...};

  ArtsXMLTag tag;
  String type;

  tag.read_from_stream(is_xml);
  tag.check_name("Variant");
  tag.get_attribute_value("type", type);

  var = init<0, Types...>(std::distance(
      valtype.begin(), std::find(valtype.begin(), valtype.end(), type)));

  std::visit([&is_xml, pbifs](
                 auto& value) { xml_read_from_stream(is_xml, value, pbifs); },
             var);

  tag.read_from_stream(is_xml);
  tag.check_name("/Variant");
}
