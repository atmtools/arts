#include "xml_io_stream_core.h"

#include <double_imanip.h>

#include "xml_io_base.h"

//! NUMERIC

void xml_io_stream<Numeric>::put(std::span<const Numeric> v, bofstream* pbofs) {
  for (auto& elem : v) *pbofs << elem;
}

void xml_io_stream<Numeric>::write(std::ostream& os,
                                   const Numeric& n,
                                   bofstream* pbofs,
                                   std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  if (pbofs)
    xml_io_stream<Numeric>::put({&n, 1}, pbofs);
  else
    std::println(os, "{}", n);

  tag.write_to_end_stream(os);
}

void xml_io_stream<Numeric>::get(std::span<Numeric> v, bifstream* pbifs) {
  pbifs->readDoubleArray(v.data(), v.size());
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
    xml_io_stream<Numeric>::get({&n, 1}, pbifs);
  } else {
    xml_io_stream<Numeric>::parse({&n, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

//! COMPLEX

void xml_io_stream<std::complex<Numeric>>::put(
    std::span<const std::complex<Numeric>> v, bofstream* pbofs) {
  xml_io_stream<Numeric>::put(
      {reinterpret_cast<const Numeric*>(v.data()), 2 * v.size()}, pbofs);
}

void xml_io_stream<std::complex<Numeric>>::write(std::ostream& os,
                                                 const std::complex<Numeric>& n,
                                                 bofstream* pbofs,
                                                 std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  if (pbofs)
    xml_io_stream<std::complex<Numeric>>::put({&n, 1}, pbofs);
  else
    std::println(os, "{} {}", n.real(), n.imag());

  tag.write_to_end_stream(os);
}

void xml_io_stream<std::complex<Numeric>>::get(
    std::span<std::complex<Numeric>> v, bifstream* pbifs) {
  xml_io_stream<Numeric>::get(
      {reinterpret_cast<Numeric*>(v.data()), 2 * v.size()}, pbifs);
}

void xml_io_stream<std::complex<Numeric>>::parse(
    std::span<std::complex<Numeric>> v, std::istream& is) {
  xml_io_stream<Numeric>::parse(
      {reinterpret_cast<Numeric*>(v.data()), 2 * v.size()}, is);
}

void xml_io_stream<std::complex<Numeric>>::read(std::istream& is,
                                                std::complex<Numeric>& n,
                                                bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    xml_io_stream<std::complex<Numeric>>::get({&n, 1}, pbifs);
  } else {
    xml_io_stream<std::complex<Numeric>>::parse({&n, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

//! INDEX

void xml_io_stream<Index>::put(std::span<const Index> v, bofstream* pbofs) {
  for (auto& elem : v) *pbofs << elem;
}

void xml_io_stream<Index>::write(std::ostream& os,
                                 const Index& n,
                                 bofstream* pbofs,
                                 std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  if (pbofs)
    xml_io_stream<Index>::put({&n, 1}, pbofs);
  else
    std::println(os, "{}", n);

  tag.write_to_end_stream(os);
}

void xml_io_stream<Index>::get(std::span<Index> v, bifstream* pbifs) {
  for (auto& elem : v) *pbifs >> elem;
}

void xml_io_stream<Index>::parse(std::span<Index> v, std::istream& is) {
  for (auto& x : v) is >> x;
}

void xml_io_stream<Index>::read(std::istream& is, Index& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    xml_io_stream<Index>::get({&n, 1}, pbifs);
  } else {
    xml_io_stream<Index>::parse({&n, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

//! SIZE --- NOTE IT BELIEVES IT IS AN INDEX IN BINARY XML

void xml_io_stream<Size>::put(std::span<const Size> v, bofstream* pbofs) {
  //! FIXME: This is a workaround because we cannot write Size directly
  for (auto& elem : v) {
    if (elem > std::numeric_limits<Index>::max()) {
      throw std::runtime_error("Cannot store binary.  Size is larger than Index allows.");
    }

    *pbofs << static_cast<Index>(elem);
  }
}

void xml_io_stream<Size>::write(std::ostream& os,
                                const Size& n,
                                bofstream* pbofs,
                                std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  if (pbofs)
    xml_io_stream<Size>::put({&n, 1}, pbofs);
  else
    std::println(os, "{}", n);

  tag.write_to_end_stream(os);
}

void xml_io_stream<Size>::get(std::span<Size> v, bifstream* pbifs) {
  for (auto& elem : v) {
    Index workaround;
    *pbifs >> workaround;
    elem = static_cast<Size>(workaround);
  }
}

void xml_io_stream<Size>::parse(std::span<Size> v, std::istream& is) {
  for (auto& x : v) is >> x;
}

void xml_io_stream<Size>::read(std::istream& is, Size& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    xml_io_stream<Size>::get({&n, 1}, pbifs);
  } else {
    xml_io_stream<Size>::parse({&n, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

//! STRING

void xml_io_stream<String>::write(std::ostream& os_xml,
                                  const String& str,
                                  bofstream*,
                                  std::string_view name) {
  XMLTag open_tag;
  XMLTag close_tag;

  open_tag.name = "String";
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);

  std::print(os_xml, R"("{}")", str);

  close_tag.name = "/String";
  close_tag.write_to_stream(os_xml);
}

void xml_io_stream<String>::read(std::istream& is_xml,
                                 String& str,
                                 bifstream*) {
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

//! BOOL

void xml_io_stream<bool>::read(std::istream& is, bool& n, bifstream*) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  is >> n;

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

//! ANY

void xml_io_stream<bool>::write(std::ostream& os,
                                const bool& n,
                                bofstream*,
                                std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  std::println(os, "{}", Index{n});
  tag.write_to_end_stream(os);
}

void xml_io_stream<Any>::write(std::ostream& os,
                               const Any&,
                               bofstream*,
                               std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  tag.write_to_end_stream(os);
}

void xml_io_stream<Any>::read(std::istream& is, Any&, bifstream*) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name("Any");
  tag.read_from_stream(is);
  tag.check_end_name("Any");
}
