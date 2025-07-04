#include "xml_io_stream_core.h"

#include <double_imanip.h>

#include "xml_io_base.h"

void xml_io_stream<std::complex<Numeric>>::put(
    std::span<const std::complex<Numeric>> v, bofstream* pbofs) {
  pbofs->putRaw(reinterpret_cast<const char*>(v.data()),
                v.size() * sizeof(std::complex<Numeric>));
}

void xml_io_stream<std::complex<Numeric>>::write(std::ostream& os,
                                                 const std::complex<Numeric>& n,
                                                 bofstream* pbofs,
                                                 std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(std::span{&n, 1}, pbofs);
  else
    std::println(os, "{} {}", n.real(), n.imag());

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Numeric>::put(std::span<const Numeric> v, bofstream* pbofs) {
  pbofs->putRaw(reinterpret_cast<const char*>(v.data()),
                v.size() * sizeof(Numeric));
}

void xml_io_stream<Numeric>::write(std::ostream& os,
                                   const Numeric& n,
                                   bofstream* pbofs,
                                   std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(std::span{&n, 1}, pbofs);
  else
    std::println(os, "{}", n);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Index>::put(std::span<const Index> v, bofstream* pbofs) {
  pbofs->putRaw(reinterpret_cast<const char*>(v.data()),
                v.size() * sizeof(Index));
}

void xml_io_stream<Index>::write(std::ostream& os,
                                 const Index& n,
                                 bofstream* pbofs,
                                 std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(std::span{&n, 1}, pbofs);
  else
    std::println(os, "{}", n);

  std::println(os, R"(</{0}>)", type_name);
}
void xml_io_stream<bool>::write(std::ostream& os,
                                const bool& n,
                                bofstream*,
                                std::string_view name) {
  std::println(os, R"(<{0} name="{1}"> {2} </{0}>)", type_name, name, Index{n});
}

void xml_io_stream<Size>::put(std::span<const Size> v, bofstream* pbofs) {
  pbofs->putRaw(reinterpret_cast<const char*>(v.data()),
                v.size() * sizeof(Size));
}

void xml_io_stream<Size>::write(std::ostream& os,
                                const Size& n,
                                bofstream* pbofs,
                                std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(std::span{&n, 1}, pbofs);
  else
    std::println(os, "{}", n);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<String>::write(std::ostream& os_xml,
                                  const String& str,
                                  bofstream* pbofs,
                                  std::string_view name) {
  XMLTag open_tag;
  XMLTag close_tag;

  open_tag.set_name("String");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);

  std::print(os_xml, R"("{}")", str);

  close_tag.set_name("/String");
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);
}

void xml_io_stream<std::complex<Numeric>>::get(
    std::span<std::complex<Numeric>> v, bifstream* pbifs) {
  pbifs->getRaw(reinterpret_cast<char*>(v.data()),
                v.size() * sizeof(std::complex<Numeric>));
}

void xml_io_stream<std::complex<Numeric>>::parse(
    std::span<std::complex<Numeric>> v, std::istream& is) {
  xml_io_stream<Numeric>::parse(
      std::span{reinterpret_cast<Numeric*>(v.data()), 2 * v.size()}, is);
}

void xml_io_stream<std::complex<Numeric>>::read(std::istream& is,
                                                std::complex<Numeric>& n,
                                                bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(std::span{&n, 1}, pbifs);
  } else {
    parse(std::span{&n, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Numeric>::get(std::span<Numeric> v, bifstream* pbifs) {
  pbifs->getRaw(reinterpret_cast<char*>(v.data()), v.size() * sizeof(Numeric));
}

void xml_io_stream<Numeric>::parse(std::span<Numeric> v, std::istream& is) {
  for (auto& x : v) is >> double_imanip() >> x;
}

void xml_io_stream<Numeric>::read(std::istream& is,
                                  Numeric& n,
                                  bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(std::span{&n, 1}, pbifs);
  } else {
    parse(std::span{&n, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Index>::get(std::span<Index> v, bifstream* pbifs) {
  pbifs->getRaw(reinterpret_cast<char*>(v.data()), v.size() * sizeof(Index));
}

void xml_io_stream<Index>::parse(std::span<Index> v, std::istream& is) {
  for (auto& x : v) is >> x;
}

void xml_io_stream<Index>::read(std::istream& is, Index& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(std::span{&n, 1}, pbifs);
  } else {
    parse(std::span{&n, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<bool>::read(std::istream& is, bool& n, bifstream*) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  is >> n;

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Size>::get(std::span<Size> v, bifstream* pbifs) {
  pbifs->getRaw(reinterpret_cast<char*>(v.data()), v.size() * sizeof(Size));
}

void xml_io_stream<Size>::parse(std::span<Size> v, std::istream& is) {
  for (auto& x : v) is >> x;
}

void xml_io_stream<Size>::read(std::istream& is, Size& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(std::span{&n, 1}, pbifs);
  } else {
    parse(std::span{&n, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<String>::read(std::istream& is_xml,
                                 String& str,
                                 bifstream* pbifs) {
  XMLTag tag;
  char dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("String");

  // Skip whitespaces
  bool string_starts_with_quotes = true;
  do {
    is_xml >> dummy;
    switch (dummy) {
      case ' ':
      case '\"':
      case '\n':
      case '\r':
      case '\t': break;
      default:   string_starts_with_quotes = false;
    }
  } while (is_xml.good() && dummy != '"' && string_starts_with_quotes);

  // Throw exception if first char after whitespaces is not a quote
  if (!string_starts_with_quotes) {
    xml_parse_error("String must begin with \"");
  }

  //catch case where string is empty. CPD 29/8/05
  dummy = (char)is_xml.peek();
  if (dummy == '"') {
    str = "";
  } else {
    std::stringbuf strbuf;

    is_xml.get(strbuf, '"');
    if (is_xml.fail()) {
      xml_parse_error("String must end with \"");
    }
    str = strbuf.str();
  }

  // Ignore quote
  is_xml >> dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("/String");
}

void xml_io_stream<Any>::write(std::ostream& os,
                               const Any&,
                               bofstream*,
                               std::string_view) {
  std::println(os, "<{0}> </{0}>", type_name);
}

void xml_io_stream<Any>::read(std::istream& is, Any&, bifstream*) {
  std::string x;
  is >> x >> x;
}
