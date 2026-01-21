/**
 * @file   rational.cc
 * @author Richard Larsson
 * @date   2012-10-31
 * 
 * @brief  Contains the rational class implmentations
 **/

#include "rational.h"

#include <debug.h>
#include <mystring.h>
#include <xml_io_base.h>

#include <ostream>

namespace {
constexpr Rational fixSign(Rational r) noexcept {
  if (r.denom < 0) {
    r.numer = -r.numer;
    r.denom = -r.denom;
  }
  return r;
}

constexpr Rational numeric2rational(Numeric x, size_t maxdec = 4) noexcept {
  Index nom = 0, denom = 1;

  // Keep track of sign independently
  const bool signchange = x < 0;
  x                     = signchange ? -x : x;

  // Add numbers by keeping the floor
  size_t i = 0;
  do {
    const auto xi  = Index(x);
    nom           += xi;
    x              = 10 * (x - Numeric(xi));
    nom           *= 10;
    denom         *= 10;
    i++;
  } while (i <= maxdec);

  // Fix possible rounding error
  if (x >= 5) nom += 10;

  // Change sign or not
  return signchange ? Rational{-nom, denom} : Rational{nom, denom};
}
}  // namespace

Numeric sqrt(const Rational r) { return std::sqrt(r.toNumeric()); }

Numeric pow(const Rational base, Numeric exp) {
  return std::pow(base.toNumeric(), exp);
}

Numeric pow(Numeric base, const Rational exp) {
  return std::pow(base, exp.toNumeric());
}

Numeric pow(const Rational base, const Rational exp) {
  return pow(base, exp.toNumeric());
}

std::ostream& operator<<(std::ostream& os, const Rational& a) {
  Rational r = Rational::reduce_by_gcd(a);
  r          = fixSign(r);

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

  return is;
}

Rational::Rational(const std::string_view s) {
  auto len = s.length();

  ARTS_USER_ERROR_IF(0 == len, "Cannot create Rational from empty string");

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
}

void simplify_in_place(Rational& r) noexcept {
  r = fixSign(Rational::reduce_by_gcd(r));
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
