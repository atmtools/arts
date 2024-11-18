#pragma once

#include <workspace.h>
#include <xml_io.h>

#include <unordered_map>

template <WorkspaceGroup Key, WorkspaceGroup Type>
void xml_write(std::ostream& os_xml,
               const std::unordered_map<Key, Type>& map,
               bofstream* pbofs,
               const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Map");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", String{WorkspaceGroupInfo<Type>::name});
  open_tag.add_attribute("key", String{WorkspaceGroupInfo<Key>::name});
  open_tag.add_attribute("nelem", static_cast<Index>(map.size()));

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (const auto& [key, value] : map) {
    xml_write_to_stream(os_xml, key, pbofs, "");
    xml_write_to_stream(os_xml, value, pbofs, "");
  }

  close_tag.set_name("/Map");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

template <WorkspaceGroup Key, WorkspaceGroup Type>
void xml_read(std::istream& is_xml,
              std::unordered_map<Key, Type>& map,
              bifstream* pbifs) {
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Map");
  tag.check_attribute("type", String{WorkspaceGroupInfo<Type>::name});
  tag.check_attribute("key", String{WorkspaceGroupInfo<Key>::name});

  tag.get_attribute_value("nelem", nelem);
  map.reserve(nelem);

  Index n;
  for (n = 0; n < nelem; n++) {
    Key key;
    xml_read_from_stream(is_xml, key, pbifs);
    xml_read_from_stream(is_xml, map[std::move(key)], pbifs);
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Map");
}
