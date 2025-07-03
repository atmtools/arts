#include "xml_io_stream_core.h"

#include <double_imanip.h>

#include "xml_io_base.h"

void xml_io_stream<std::complex<Numeric>>::put(
    std::span<const std::complex<Numeric>> v, bofstream* pbofs) {
  assert(pbofs);
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
  assert(pbofs);
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
  assert(pbofs);
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
  assert(pbofs);
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

void xml_io_stream<String>::write(std::ostream& os,
                                  const String& n,
                                  bofstream* pbofs,
                                  std::string_view name) {
  if (pbofs) {
    std::println(os,
                 R"(<{0} name="{2}" nelem="{3}"></{0}>)",
                 type_name,
                 n,
                 name,
                 n.size());
    pbofs->putRaw(n.data(), n.size());
  } else {
    std::println(os,
                 R"(<{0} name="{2}" nelem="{3}">{1}</{0}>)",
                 type_name,
                 n,
                 name,
                 n.size());
  }
}

void xml_io_stream<std::complex<Numeric>>::get(
    std::span<std::complex<Numeric>> v, bifstream* pbifs) {
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v.data()),
                v.size() * sizeof(std::complex<Numeric>));
}

void xml_io_stream<std::complex<Numeric>>::parse(
    std::span<std::complex<Numeric>> v, std::istream& is) {
  assert(pbifs);
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
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v.data()), v.size() * sizeof(Numeric));
}

void xml_io_stream<Numeric>::parse(std::span<Numeric> v, std::istream& is) {
  assert(pbifs);
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
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v.data()), v.size() * sizeof(Index));
}

void xml_io_stream<Index>::parse(std::span<Index> v, std::istream& is) {
  assert(pbifs);
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
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v.data()), v.size() * sizeof(Size));
}

void xml_io_stream<Size>::parse(std::span<Size> v, std::istream& is) {
  assert(pbifs);
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

void xml_io_stream<String>::read(std::istream& is,
                                 String& n,
                                 bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);
  Size nelem = 0;
  tag.get_attribute_value("nelem", nelem);
  n.resize(nelem);

  if (pbifs) {
    pbifs->getRaw(n.data(), nelem);
  } else {
    is.read(n.data(), nelem);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
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
