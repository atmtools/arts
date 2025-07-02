#include "xml_io_stream_core.h"

#include "double_imanip.h"
#include "xml_io_base.h"

void xml_io_stream<std::complex<Numeric>>::put(
    const std::complex<Numeric>* const v, bofstream* pbofs, Size n) {
  assert(pbofs);
  pbofs->putRaw(reinterpret_cast<const char*>(v),
                n * sizeof(std::complex<Numeric>));
}

void xml_io_stream<std::complex<Numeric>>::write(std::ostream& os,
                                                 const std::complex<Numeric>& n,
                                                 bofstream* pbofs,
                                                 std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(&n, pbofs);
  else
    std::println(os, "{} {}", n.real(), n.imag());

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Numeric>::put(const Numeric* const v,
                                 bofstream* pbofs,
                                 Size n) {
  assert(pbofs);
  pbofs->putRaw(reinterpret_cast<const char*>(v), n * sizeof(Numeric));
}

void xml_io_stream<Numeric>::write(std::ostream& os,
                                   const Numeric& n,
                                   bofstream* pbofs,
                                   std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(&n, pbofs);
  else
    std::println(os, "{}", n);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Index>::put(const Index* const v, bofstream* pbofs, Size n) {
  assert(pbofs);
  pbofs->putRaw(reinterpret_cast<const char*>(v), n * sizeof(Index));
}

void xml_io_stream<Index>::write(std::ostream& os,
                                 const Index& n,
                                 bofstream* pbofs,
                                 std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(&n, pbofs);
  else
    std::println(os, "{}", n);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<bool>::put(const bool* const v, bofstream* pbofs, Size n) {
  assert(pbofs);
  pbofs->putRaw(reinterpret_cast<const char*>(v), n * sizeof(bool));
}

void xml_io_stream<bool>::write(std::ostream& os,
                                const bool& n,
                                bofstream* pbofs,
                                std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(&n, pbofs);
  else
    std::println(os, "{}", Index{n});

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Size>::put(const Size* const v, bofstream* pbofs, Size n) {
  assert(pbofs);
  pbofs->putRaw(reinterpret_cast<const char*>(v), n * sizeof(Size));
}

void xml_io_stream<Size>::write(std::ostream& os,
                                const Size& n,
                                bofstream* pbofs,
                                std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    put(&n, pbofs);
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

void xml_io_stream<std::complex<Numeric>>::get(std::complex<Numeric>* v,
                                               bifstream* pbifs,
                                               Size n) {
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v), n * sizeof(std::complex<Numeric>));
}

void xml_io_stream<std::complex<Numeric>>::read(std::istream& is,
                                                std::complex<Numeric>& n,
                                                bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(&n, pbifs);
  } else {
    Numeric r{}, i{};
    is >> double_imanip() >> r >> i;
    n = std::complex<Numeric>{r, i};
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Numeric>::get(Numeric* v, bifstream* pbifs, Size n) {
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v), n * sizeof(Numeric));
}

void xml_io_stream<Numeric>::read(std::istream& is,
                                  Numeric& n,
                                  bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(&n, pbifs);
  } else {
    is >> double_imanip() >> n;
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Index>::get(Index* v, bifstream* pbifs, Size n) {
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v), n * sizeof(Index));
}

void xml_io_stream<Index>::read(std::istream& is, Index& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(&n, pbifs);
  } else {
    is >> n;
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<bool>::get(bool* v, bifstream* pbifs, Size n) {
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v), n * sizeof(bool));
}

void xml_io_stream<bool>::read(std::istream& is, bool& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(&n, pbifs);
  } else {
    is >> n;
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Size>::get(Size* v, bifstream* pbifs, Size n) {
  assert(pbifs);
  pbifs->getRaw(reinterpret_cast<char*>(v), n * sizeof(Size));
}

void xml_io_stream<Size>::read(std::istream& is, Size& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(&n, pbifs);
  } else {
    is >> n;
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
