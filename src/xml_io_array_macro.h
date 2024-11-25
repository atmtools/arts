#pragma once

#include "xml_io.h"

//! Both T and T{}[0] are ARTS groups exposed to the user if this is true
template <typename T>
concept array_of_group = WorkspaceGroup<std::remove_cvref_t<T>> and
                         WorkspaceGroup<std::remove_cvref_t<decltype(T{}[0])>>;

template <array_of_group T>
void xml_read(std::istream &is_xml, T &at, bifstream *pbifs) try {
  const static String subtype{
      WorkspaceGroupInfo<std::remove_cvref_t<decltype(T{}[0])>>::name};

  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", subtype);

  tag.get_attribute_value("nelem", nelem);
  at.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++) xml_read_from_stream(is_xml, at[n], pbifs);
  } catch (const std::runtime_error &e) {
    std::ostringstream os;
    os << "Error reading " << WorkspaceGroupInfo<std::remove_cvref_t<T>>::name
       << ": "
       << "\n Element: " << n << "\n"
       << e.what();
    throw std::runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
} catch (std::runtime_error &e) {
  throw std::runtime_error(
      var_string("Failed reading routine for ",
                 WorkspaceGroupInfo<std::remove_cvref_t<T>>::name,
                 "\nError reads:\n",
                 e.what()));
}

template <array_of_group T>
void xml_write(std::ostream &os_xml,
               const T &at,
               bofstream *pbofs,
               const String &name) try {
  const static String subtype{
      WorkspaceGroupInfo<std::remove_cvref_t<decltype(T{}[0])>>::name};

  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", subtype);
  open_tag.add_attribute("nelem", static_cast<Index>(at.size()));

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Size n = 0; n < at.size(); n++)
    xml_write_to_stream(os_xml, at[n], pbofs, std::format("{}", n));

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
} catch (std::runtime_error &e) {
  throw std::runtime_error(
      var_string("Failed saving routine for ",
                 WorkspaceGroupInfo<std::remove_cvref_t<T>>::name,
                 "\nError reads:\n",
                 e.what()));
}

//! Helper macro for when both Array<T> and T are ARTS groups
#define TMPL_XML_READ_WRITE_STREAM_ARRAY(T)                                  \
  void xml_read_from_stream(std::istream &is_xml, T &at, bifstream *pbifs) { \
    xml_read(is_xml, at, pbifs);                                             \
  }                                                                          \
                                                                             \
  void xml_write_to_stream(std::ostream &os_xml,                             \
                           const T &at,                                      \
                           bofstream *pbofs,                                 \
                           const String &name) {                             \
    xml_write(os_xml, at, pbofs, name);                                      \
  }
