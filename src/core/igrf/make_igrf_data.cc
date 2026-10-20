#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <ostream>
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace {
//! The content of an IGRF coefficient file as distributed by IAGA
struct IgrfFile {
  //! The comment lines at the top of the file, without the leading '#'
  std::vector<std::string> comments;

  //! The model of each epoch, "IGRF" or "DGRF"
  std::vector<std::string> models;

  //! The epochs [decimal year]
  std::vector<double> epochs;

  //! The end of the secular variation, which starts at the last epoch [decimal year]
  double sv_end{};

  //! The largest degree
  int nmax{};

  //! The coefficients of each epoch, and the secular variation last, keyed by ('g' or 'h', n, m)
  std::map<std::tuple<char, int, int>, std::vector<double>> rows;
};

std::vector<std::string_view> split(std::string_view line) {
  std::vector<std::string_view> out;

  constexpr std::string_view ws = " \t\r";
  while (true) {
    const auto first = line.find_first_not_of(ws);
    if (first == std::string_view::npos) break;
    line              = line.substr(first);
    const auto length = line.find_first_of(ws);
    out.push_back(line.substr(0, length));
    if (length == std::string_view::npos) break;
    line = line.substr(length);
  }

  return out;
}

template <typename T> T parse(std::string_view s, std::size_t iline) {
  T x{};

  const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), x);
  if (ec != std::errc{} or ptr != s.data() + s.size()) {
    throw std::runtime_error(std::format("Line {}: cannot read a number from \"{}\"", iline, s));
  }

  if constexpr (std::floating_point<T>) {
    if (not std::isfinite(x)) throw std::runtime_error(std::format("Line {}: \"{}\" is not finite", iline, s));
  }

  return x;
}

//! Parses the secular variation interval, e.g., "2025-30", and returns its end, e.g., 2030
double parse_sv_end(std::string_view s, double last_epoch, std::size_t iline) {
  const auto dash = s.find('-');
  if (dash == std::string_view::npos) {
    throw std::runtime_error(std::format("Line {}: \"{}\" is not a secular variation interval", iline, s));
  }

  const int start = parse<int>(s.substr(0, dash), iline);
  int       end   = parse<int>(s.substr(dash + 1), iline);
  if (s.size() - dash - 1 == 2) end += start - start % 100 + (end < start % 100 ? 100 : 0);

  if (start != last_epoch or end <= start) {
    throw std::runtime_error(std::format(
        "Line {}: the secular variation interval \"{}\" does not start at the last epoch {}", iline, s, last_epoch));
  }

  return end;
}

IgrfFile read(const std::filesystem::path& path) {
  std::ifstream file(path);
  if (not file) throw std::runtime_error(std::format("Cannot open {}", path.string()));

  IgrfFile out;

  std::size_t ncols = 0;
  std::size_t iline = 0;
  std::string line;
  while (std::getline(file, line)) {
    ++iline;

    const auto tokens = split(line);
    if (tokens.empty()) continue;

    if (tokens.front().starts_with('#')) {
      if (ncols == 0) out.comments.emplace_back(std::string_view{line}.substr(line.find('#') + 1));
      continue;
    }

    if (tokens.front() == "c/s") {
      if (tokens.size() < 5 or tokens[1] != "deg" or tokens[2] != "ord" or tokens.back() != "SV") {
        throw std::runtime_error(std::format("Line {}: expected \"c/s deg ord <model>... SV\"", iline));
      }

      for (std::size_t i = 3; i + 1 < tokens.size(); i++) {
        if (tokens[i] != "IGRF" and tokens[i] != "DGRF") {
          throw std::runtime_error(std::format("Line {}: unknown model \"{}\"", iline, tokens[i]));
        }
        out.models.emplace_back(tokens[i]);
      }

      ncols = tokens.size();
      continue;
    }

    if (ncols == 0) throw std::runtime_error(std::format("Line {}: the \"c/s\" line must come first", iline));

    if (tokens.size() != ncols) {
      throw std::runtime_error(std::format("Line {}: expected {} columns, got {}", iline, ncols, tokens.size()));
    }

    if (tokens.front() == "g/h") {
      if (not out.epochs.empty()) throw std::runtime_error(std::format("Line {}: a second \"g/h\" line", iline));

      for (std::size_t i = 3; i + 1 < ncols; i++) {
        out.epochs.push_back(parse<double>(tokens[i], iline));
        if (out.epochs.size() > 1 and out.epochs.back() <= out.epochs[out.epochs.size() - 2]) {
          throw std::runtime_error(std::format("Line {}: the epochs must be increasing", iline));
        }
      }

      out.sv_end = parse_sv_end(tokens.back(), out.epochs.back(), iline);
      continue;
    }

    if (out.epochs.empty()) throw std::runtime_error(std::format("Line {}: the \"g/h\" line must come first", iline));

    if (tokens.front() != "g" and tokens.front() != "h") {
      throw std::runtime_error(std::format("Line {}: expected 'g' or 'h', got \"{}\"", iline, tokens.front()));
    }

    const char type = tokens.front().front();
    const int  n    = parse<int>(tokens[1], iline);
    const int  m    = parse<int>(tokens[2], iline);
    if (n < 1 or m < (type == 'g' ? 0 : 1) or m > n) {
      throw std::runtime_error(std::format("Line {}: invalid degree and order of {} ({}, {})", iline, type, n, m));
    }

    std::vector<double> values;
    for (std::size_t i = 3; i < ncols; i++) values.push_back(parse<double>(tokens[i], iline));

    if (not out.rows.try_emplace({type, n, m}, std::move(values)).second) {
      throw std::runtime_error(std::format("Line {}: {} ({}, {}) is given twice", iline, type, n, m));
    }

    out.nmax = std::max(out.nmax, n);
  }

  if (out.epochs.empty()) throw std::runtime_error(std::format("No coefficients in {}", path.string()));

  for (int n = 1; n <= out.nmax; n++) {
    for (int m = 0; m <= n; m++) {
      if (not out.rows.contains({'g', n, m})) throw std::runtime_error(std::format("g ({}, {}) is missing", n, m));
      if (m > 0 and not out.rows.contains({'h', n, m})) {
        throw std::runtime_error(std::format("h ({}, {}) is missing", n, m));
      }
    }
  }

  return out;
}

//! A literal that is a floating point number in C++, e.g., "-31543.0" and not "-31543"
std::string literal(double x) {
  std::string out = std::format("{}", x);
  if (out.find_first_of(".e") == std::string::npos) out += ".0";
  return out;
}

//! The coefficient of column i, with zeros for the n = 0 term and for h with m = 0
double value(const IgrfFile& data, char type, int n, int m, std::size_t i) {
  if (n == 0 or (type == 'h' and m == 0)) return 0.0;
  return data.rows.at({type, n, m})[i];
}

//! Writes g and h of column i as a pair of triangular matrices, one line per degree, followed by end
void write_pair(std::ostream&    os,
                const IgrfFile&  data,
                std::size_t      i,
                std::string_view title,
                std::string_view indent,
                std::string_view end) {
  constexpr std::size_t width = 10;

  std::println(os, "{}// {}", indent, title);
  for (const char type : {'g', 'h'}) {
    std::println(os, "{}{}Coefficients{{{{", indent, type == 'g' ? "{" : " ");
    for (int n = 0; n <= data.nmax; n++) {
      std::print(os, "{}   ", indent);
      for (int m = 0; m <= n; m++) std::print(os, " {:>{}},", literal(value(data, type, n, m, i)), width);
      std::println(os);
    }
    std::println(os, "{}  }}}}{}", indent, type == 'g' ? std::string{","} : std::format("}}{}", end));
  }
}

void write_comments(std::ostream& os, const IgrfFile& data, std::string_view filename) {
  std::println(os, "//! auto-generated file from {}", filename);
  for (auto& comment : data.comments) std::println(os, "//{}", comment);
}

void write_header(std::ostream& os, const IgrfFile& data, std::string_view filename) {
  std::string epochs;
  for (const double epoch : data.epochs) epochs += std::format("{}{}", epochs.empty() ? "" : ", ", literal(epoch));

  write_comments(os, data, filename);

  std::println(os,
               R"(
#pragma once

#include <configtypes.h>
#include <matpack.h>

#include <array>
#include <utility>

namespace IGRF::Data {{
//! The largest degree n
inline constexpr Size nmax = {};

//! The number of epochs
inline constexpr Size nepoch = {};

//! The epochs [decimal year]
inline constexpr std::array<Numeric, nepoch> epochs{{{}}};

//! The secular variation is valid from the last epoch to this time [decimal year]
inline constexpr Numeric sv_end = {};

//! The first and the last valid time of the coefficients [decimal year], where 2001.0 is the start of 2001
inline constexpr std::pair<Numeric, Numeric> valid_range{{{}, {}}};

//! The g or h coefficients, the (n, m) coefficient is at [n, m], zero for n = 0 and for h with m = 0
using Coefficients = matpack::triangular_cdata_t<Numeric, nmax + 1>;

//! The g and h coefficients [nT] of epoch i
[[nodiscard]] const std::pair<Coefficients, Coefficients>& coefficients(Size i);

//! The secular variation of g and h [nT/year], from the last epoch to sv_end
[[nodiscard]] const std::pair<Coefficients, Coefficients>& secular_variation();
}}  // namespace IGRF::Data)",
               data.nmax,
               data.epochs.size(),
               epochs,
               literal(data.sv_end),
               literal(data.epochs.front()),
               literal(data.sv_end));
}

void write_source(std::ostream& os, const IgrfFile& data, std::string_view filename) {
  write_comments(os, data, filename);

  std::println(os, R"(
#include "auto_igrf_data.h"

#include <cassert>

namespace IGRF::Data {{
namespace {{
constexpr std::array<std::pair<Coefficients, Coefficients>, nepoch> gh{{{{)");

  for (std::size_t i = 0; i < data.epochs.size(); i++) {
    write_pair(os, data, i, std::format("{} {}", literal(data.epochs[i]), data.models[i]), "    ", ",");
  }

  std::println(os, R"(}}}};

constexpr std::pair<Coefficients, Coefficients> dgh =)");

  write_pair(os,
             data,
             data.epochs.size(),
             std::format("secular variation {}-{}", literal(data.epochs.back()), literal(data.sv_end)),
             "  ",
             ";");

  std::println(os, R"(}}  // namespace

const std::pair<Coefficients, Coefficients>& coefficients(Size i) {{
  assert(i < nepoch);
  return gh[i];
}}

const std::pair<Coefficients, Coefficients>& secular_variation() {{ return dgh; }}
}}  // namespace IGRF::Data)");
}
}  // namespace

int main(int argc, char** argv) try {
  if (argc != 2) {
    std::println(stderr, "Usage: {} <igrf.txt>", argv[0]);
    return EXIT_FAILURE;
  }

  const std::filesystem::path path = argv[1];
  const IgrfFile              data = read(path);

  std::ofstream h("auto_igrf_data.h", std::ios::out);
  write_header(h, data, path.filename().string());

  std::ofstream cpp("auto_igrf_data.cpp", std::ios::out);
  write_source(cpp, data, path.filename().string());

  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "Error in main: {}", e.what());
  return EXIT_FAILURE;
}
