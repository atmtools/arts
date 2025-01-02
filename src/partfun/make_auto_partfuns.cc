#include <debug.h>
#include <enumsSpeciesEnum.h>
#include <matpack.h>
#include <spec/isotopologues.h>
#include <spec/species.h>
#include <xml_io_base.h>

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "xml_io_partfun.h"

namespace {
struct func_body {
  std::string data{};
  std::string main{};
  std::string deriv{};
  std::vector<std::string> includes{};
};

func_body make_cpp_string_interp(const Matrix& data) {
  assert(data.size() > 3);
  const auto Tv = data[joker, 0];
  const auto q  = data[joker, 1];

  for (Size i = 1; i < Tv.size(); i++) {
    if (Tv[i] <= Tv[i - 1])
      throw std::runtime_error("Temperature grid must be increasing");
  }

  return {.data = std::format(
              R"(
inline constexpr std::array<Numeric, {0}> coef{{ {1:,} }};
inline constexpr std::array<Numeric, {0}> grid{{ {2:,} }};
)",
              q.size(),
              q,
              Tv),
          .main     = std::format(R"(
  const Index i_low =
      std::distance(grid.cbegin(), std::lower_bound(grid.cbegin(), grid.cend(), T));
  const Size i = std::min<Size>(i_low - (i_low > 0), Size{{{}}});

  return coef[i] + (T - grid[i]) * (coef[i + 1] - coef[i]) / (grid[i + 1] - grid[i]);
)",
                              q.size() - 2),
          .deriv    = std::format(R"(
  const Index i_low =
      std::distance(grid.cbegin(), std::lower_bound(grid.cbegin(), grid.cend(), T));
  const Size i = std::min<Size>(i_low - (i_low > 0), Size{{{}}});

  return (coef[i + 1] - coef[i]) / (grid[i + 1] - grid[i]);
)",
                               q.size() - 2),
          .includes = {"<algorithm>", "<array>"}};
}

func_body make_cpp_string_coeff(const Matrix& data) {
  assert(data.size() > 1);

  const auto q = data[joker, 0];

  return {.data = std::format(
              R"(
inline constexpr std::array<Numeric, {0}> Q{{ {1:,} }};
)",
              q.size(),
              q),
          .main     = std::format(R"(
  Numeric result = Q[0];
  Numeric TN     = 1.0;

  for (int i = 1; i < {0}; i++) {{
    TN     *= T;
    result += TN * Q[i];
  }}

  return result;
)",
                              q.size()),
          .deriv    = std::format(R"(
  Numeric result = Q[1];
  Numeric TN     = 1.0;

  for (int i = 2; i < {0}; i++) {{
    TN     *= T;
    result += static_cast<Numeric>(i) * TN * Q[i];
  }}

  return result;
)",
                               q.size()),
          .includes = {"<array>"}};
}

func_body make_cpp_string_static_interp(const Matrix& data) {
  assert(data.size() > 3);
  const auto Tv = data[joker, 0];
  const auto q  = data[joker, 1];

  const Numeric r_dT = 1.0 / (Tv[1] - Tv[0]);
  for (Size i = 1; i < Tv.size(); i++) {
    if (Tv[i] <= Tv[i - 1])
      throw std::runtime_error("Temperature grid must be increasing");
    if (not is_same_within_epsilon(r_dT, 1.0 / (Tv[i] - Tv[i - 1])))
      throw std::runtime_error("Temperature grid must be equidistant");
  }

  return {.data = std::format(
              R"(
inline constexpr std::array<Numeric, {0}> Q{{ {1:,} }};
)",
              q.size(),
              q),
          .main     = std::format(R"(
  const Numeric Tx = (T - {1}) * {0};
  const Size iTx   = static_cast<Size>(Tx);
  const Size i     = (iTx > {2}) ? {2} : iTx;
  const Numeric To = Tx - static_cast<Numeric>(i);

  return Q[i] + To * (Q[i + 1] - Q[i]);
)",
                              r_dT,
                              Tv[0],
                              Tv.size() - 2),
          .deriv    = std::format(R"(
  const Numeric Tx = (T - {1}) * {0};
  const Size iTx   = static_cast<Size>(Tx);
  const Size i     = (iTx > {2}) ? {2} : iTx;

  return (Q[i + 1] - Q[i]) * {0};
)",
                               r_dT,
                               Tv[0],
                               Tv.size() - 2),
          .includes = {"<array>"}};
}

std::string make_call_signature(const std::string_view filename) {
  const Size N = filename.size();
  assert(N > 4);

  std::string out;

  out.reserve(N);
  for (Size i = 0; i < N - 4; i++) {
    if (filename[i] == '-') continue;
    if (filename[i] == '+') {
      out.push_back('p');
      out.push_back('l');
      out.push_back('u');
      out.push_back('s');
    } else {
      out.push_back(filename[i]);
    }
  }

  return out;
}

//! take XY+-X.xml and return XYZ.cpp (change 'xml' to 'cpp')
std::string make_cpp_filename(const std::string_view filename) {
  const Size N = filename.size();

  assert(N > 4);

  std::string out;

  out        = filename;
  out[N - 3] = 'c';
  out[N - 2] = 'p';
  out[N - 1] = 'p';

  return out;
}

struct head_data {
  SpeciesEnum spec;
  std::string isot;
  std::string call;
};

head_data make_header_data(const std::string_view filename) {
  const Size N = filename.size();

  assert(N > 4);

  SpeciesEnum spec;
  std::string isot;
  isot.reserve(N);
  for (Size i = 0; i < N - 4; i++) {
    if (filename[i] == '-') {
      spec = to<SpeciesEnum>(isot);
      isot.clear();
      continue;
    }
    isot.push_back(filename[i]);
  }

  return {.spec = spec, .isot = isot, .call = make_call_signature(filename)};
}

//! take XY+-X.xml and return XYplusZ (drop '-' and extension, replace '+' with "plus")
std::string make_cpp_string_impl(const std::string_view filename,
                                 const func_body& body) {
  std::string includes{"#include <configtypes.h>\n\n"};

  for (const auto& include : body.includes) {
    includes += std::format("#include {}\n", include);
  }

  return std::format(R"(//! auto-generated file

{4}

namespace {{{1}}}  // namespace

namespace PartitionFunctions {{
Numeric {0}(Numeric T) noexcept {{{2}}}

Numeric d{0}(Numeric T) noexcept {{{3}}}
}}  // namespace PartitionFunctions
)",
                     make_call_signature(filename),
                     body.data,
                     body.main,
                     body.deriv,
                     includes);
}

std::string make_cpp_string(const std::string_view filename,
                            const PartitionFunctionsData& data) {
  switch (data.type) {
    using enum PartitionFunctionsType;
    case Interp:
      return make_cpp_string_impl(filename, make_cpp_string_interp(data.data));
    case Coeff:
      return make_cpp_string_impl(filename, make_cpp_string_coeff(data.data));
    case StaticInterp:
      return make_cpp_string_impl(filename,
                                  make_cpp_string_static_interp(data.data));
  }

  return {};
}

std::string make_h_string(
    std::map<SpeciesEnum, std::vector<std::pair<std::string, std::string>>>&
        data) {
  std::string exists;
  std::string calls;
  std::string compute;
  exists.reserve(1025);
  calls.reserve(1025);
  compute.reserve(1025);

  for (auto& spec : enumtyps::SpeciesEnumTypes) {
    const auto& v                  = data[spec];
    const std::string_view specstr = toString<0>(spec);

    exists += std::format(
        "inline constexpr std::array<std::string_view, {0}> has{1}{{",
        v.size(),
        specstr);

    compute += std::format(
        R"(template<Derivatives deriv>
Numeric compute{0}(Numeric T [[maybe_unused]], const std::string_view isot) {{
  using enum Derivatives;
)",
        specstr);

    for (Size i = 0; i < v.size(); i++) {
      auto& [isot, call] = v[i];

      exists += std::format(R"( "{}"sv,
)",
                            isot);

      calls += std::format(
          R"(Numeric  {0}(Numeric T) noexcept;
Numeric d{0}(Numeric T) noexcept;
)",
          call);

      compute += std::format(
          R"(
  if (isot == std::get<{0}>(has{1})) {{
    if constexpr (deriv == Yes) return  {2}(T);
    if constexpr (deriv == No)  return d{2}(T);
  }}
)",
          i,
          specstr,
          call);
    }
    exists += "};\n\n";

    compute += std::format(
        R"(
  ARTS_USER_ERROR("No partition function for {}-{{}}", isot);
}}
)",
        spec);
  }

  return std::format(
      R"(//! auto-generated file

#pragma once

#include <configtypes.h>
#include <debug.h>

#include <array>
#include <string_view>

enum class Derivatives : bool {{
  Yes = true,
  No  = false
}};

namespace PartitionFunctions {{
{}

{}

{}
}}  // namespace PartitionFunctions
)",
      exists,
      calls,
      compute);
}

void make_files(const std::vector<std::string>& full_filenames) {
  std::map<SpeciesEnum, std::vector<std::pair<std::string, std::string>>> data;

  for (auto& full_filename : full_filenames) {
    auto filename = std::filesystem::path(full_filename).filename().string();

    const std::string cpp = make_cpp_string(
        filename, PartitionFunctions::data_read_file(full_filename));

    std::ofstream(make_cpp_filename(filename)) << cpp;

    const auto [s, i, c] = make_header_data(filename);
    data[s].emplace_back(i, c);
  }

  std::ofstream("auto_partfun.h") << make_h_string(data);
}
}  // namespace

int main(int argn, char** argv) try {
  if (argn < 2)
    throw std::runtime_error(
        "USAGE: PROG INPUT1.XML INPUT2.xml ... INPUTX.xml");
  const std::vector<std::string> xmlfiles{argv + 1, argv + argn};

  make_files(xmlfiles);

  // for (auto& fn : xmlfiles) {
  //   auto xmlfile = std::filesystem::canonical(fn);
  //   make_cc(xmlfile);
  // }

  // make_h(xmlfiles);

  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Failed:\n" << e.what() << '\n';
  return EXIT_FAILURE;
}