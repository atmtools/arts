/**
 * @file   rational.cc
 * @author Richard Larsson
 * @date   2012-10-31
 * 
 * @brief  Contains the rational class implmentations
 **/

#include "rational.h"

#include <xml_io_base.h>

#include <ostream>

#include "debug.h"
#include "mystring.h"

std::ostream& operator<<(std::ostream& os, const Rational& a) {
  Rational r = reduce_by_gcd(a);
  r.fixSign();

  if (r.denom == 1)
    os << r.numer;
  else
    os << r.numer << "/" << r.denom;
  return os;
}

std::istream& operator>>(std::istream& is, Rational& a) {
  String s;

  is >> s;
  a = Rational(s);

  ARTS_USER_ERROR_IF(a.isUndefined(), "Cannot read {} as rational", s)

  return is;
}

Rational::Rational(const std::string_view s) {
  auto len = s.length();

  if (len) {
    auto dot_pos   = s.find('.');
    auto slash_pos = s.find('/');
    if (len > dot_pos) {
      *this = numeric2rational(std::stod(std::string{s}), len - dot_pos - 1);
    } else if (len > slash_pos) {
      const String a{s.substr(0, slash_pos)};
      const String b{s.substr(slash_pos + 1, len)};
      try {
        *this = Rational(std::stoi(a), std::stoi(b));
      } catch (...) {
        ARTS_USER_ERROR(
            "Cannot interpret either '{}' or '{}' as an integer (or neither)",
            a,
            b);
      }
    } else {
      try {
        *this = Rational(std::stoi(std::string{s}));
      } catch (...) {
        ARTS_USER_ERROR("Cannot interpret '{}' as an integer", s);
      }
    }
  } else {
    *this = RATIONAL_UNDEFINED;
  }
}

void Rational::simplify_in_place() noexcept {
  Rational a = reduce_by_gcd(*this);
  numer      = a.numer;
  denom      = a.denom;
  fixSign();
}

bifstream& Rational::read(bifstream& bif) {
  bif >> numer >> denom;
  return bif;
}

/** Binary write for Rational */
bofstream& Rational::write(bofstream& bof) const {
  bof << numer << denom;
  return bof;
}

void xml_io_stream<Rational>::write(std::ostream& os_xml,
                                    const Rational& rational,
                                    bofstream* pbofs,
                                    std::string_view name) {
  std::println(os_xml, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs)
    *pbofs << rational;
  else
    os_xml << rational;

  std::println(os_xml, R"(</{0}>)", type_name);
}

void xml_io_stream<Rational>::read(std::istream& is_xml,
                                   Rational& rational,
                                   bifstream* pbifs) {
  XMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  if (pbifs) {
    *pbifs >> rational;
    if (pbifs->fail()) {
      xml_data_parse_error(tag, "");
    }
  } else {
    is_xml >> rational;
    if (is_xml.fail()) {
      xml_data_parse_error(tag, "");
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
}
