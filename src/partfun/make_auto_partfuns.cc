#include <debug.h>
#include <gridded_fields.h>
#include <isotopologues.h>
#include <species.h>
#include <template_partfun.h>
#include <xml_io_base.h>
#include <xml_io_partfun.h>

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

struct file_wrap {
  std::ofstream os;
  template <typename... Headers>
  file_wrap(const std::string& x, Headers&&... hs) : os(x) {
    if (std::filesystem::path(x).extension() == ".h") os << "#pragma once\n\n";
    ((os << "#include <" << std::forward<Headers>(hs) << ">\n"), ...);
    os << '\n';
    os << "namespace PartitionFunctions {\n";
  }
  ~file_wrap() { os << "}  //namespace PartitionFunctions \n"; }

  template <typename T>
  friend std::ostream& operator<<(file_wrap& lhs, T&& rhs) {
    return lhs.os << std::forward<T>(rhs);
  }
};

void print_data(const PartitionFunctionsData& data, auto& os) {
  constexpr int cutline = 10;
  const Index n = data.data.nrows();

  using enum PartitionFunctions::Type;
  switch (data.type) {
    case Interp:
      os << "static constexpr std::array<Numeric, " << n << "> data{";
      for (Index i = 0; i < n; i++) {
        if (i % cutline == 0) {
          os << '\n';
        }
        os << data.data(i, 1) << ',' << ' ';
      }
      os << "};\n\n";

      os << "static constexpr std::array<Numeric, " << n << "> grid{";
      for (Index i = 0; i < n; i++) {
        if (i % cutline == 0) {
          os << '\n';
        }
        os << data.data(i, 0) << ',' << ' ';
      }
      os << "};\n";
      break;
    case Coeff:
      os << "static constexpr std::array<Numeric, " << n << "> coeff{";
      for (Index i = 0; i < n; i++) {
        if (i % cutline == 0) {
          os << '\n';
        }
        os << data.data(i, 0) << ',' << ' ';
      }
      os << "};\n";
      break;
    case StaticInterp:
        os << "static constexpr std::array<Numeric, " << n << "> data{";
        for (Index i = 0; i < n; i++) {
          if (i % cutline == 0) {
            os << '\n';
          }
          os << data.data(i, 1) << ',' << ' ';
        }
        os << "};\n\n";

        os << "constexpr Numeric dT = " << data.data(0, 0) << ";\n";
        os << "constexpr Numeric T0 = " << data.data(1, 0) - data.data(0, 0) << ";\n";
      break;
    case FINAL:
      throw std::logic_error("invalid");
  }
}

void print_method(const PartitionFunctions::Type& type, auto& os) {
  using enum PartitionFunctions::Type;
  switch (type) {
    case Interp:
      os << "return linterp<derivative>(grid, data, T);\n";
      break;
    case Coeff:
      os << "return polynom<derivative>(coeff, T);\n";
      break;
    case StaticInterp:
      os << "return static_linterp<derivative, dT, T0>(data, T);\n";
      break;
    case FINAL:
      throw std::logic_error("invalid");
  }
}

std::string spec_from_xml(const std::filesystem::path& xmlfile) {
  return std::filesystem::path(xmlfile).filename().stem().c_str();
}

std::string func_name(std::string spec, bool deriv) {
  spec.replace(spec.find('-'), 1, "");
  if (auto plus = spec.find('+'); plus not_eq spec.npos)
    spec.replace(plus, 1, "plus");
  if (deriv) spec = "d" + spec;
  return spec;
}

void make_cc(const std::filesystem::path& xmlfile) {
  if (xmlfile.extension() not_eq ".xml")
    throw std::runtime_error("Not an xml file");

  // Read data
  const PartitionFunctionsData data = [&] {
    return PartitionFunctions::data_read_file(xmlfile);
  }();

  // Names of species and function calls
  const std::string spec{spec_from_xml(xmlfile)};
  file_wrap os{spec + ".cc", "template_partfun.h"};

  // Write static data
  print_data(data, os);

  // Write function call
  os << '\n'
     << "Numeric " << func_name(spec, false)
     << "(Numeric T) noexcept {\n"
        "  constexpr auto derivative = Derivatives::No;\n  ";
  print_method(data.type, os);
  os << '}' << '\n';

  // Write derivative of function call
  os << '\n'
     << "Numeric " << func_name(spec, true)
     << "(Numeric T) noexcept {\n"
        "  constexpr auto derivative = Derivatives::Yes;\n  ";
  print_method(data.type, os);
  os << '}' << '\n';
}

std::string species_name(Species::Species spec, const std::string& isot) {
  return std::string{Species::toShortName(spec)} + "-" + isot;
}

void make_h(const std::vector<std::string>& xmlfiles) {
  // Make a complete species list and isotopologues:
  const auto data = [&] {
    std::map<Species::Species, std::vector<std::string>> x;
    for (auto& xmlfile : xmlfiles) {
      const auto spec_name = spec_from_xml(xmlfile);
      const auto delim = spec_name.find('-');
      ARTS_USER_ERROR_IF(delim == spec_name.npos,
                         "Cannot find isotopologue split")
      const auto spec = spec_name.substr(0, delim);
      const auto isot = spec_name.substr(delim + 1);
      x[Species::fromShortName(spec)].push_back(isot);
    }

    for (auto& isot: Species::Isotopologues) {
      if (x.end() == x.find(isot.spec)) x[isot.spec];
    }

    return x;
  }();

  // Open file
  file_wrap os("auto_partfun.h",
               "array",
               "debug.h",
               "matpack.h",
               "string_view",
               "template_partfun.h");

  // Write if the data exist
  for (auto& spec : data) {
    os << "inline constexpr std::array<std::string_view, " << spec.second.size()
       << "> has" << spec.first << "{\n";
    for (auto& isot : spec.second) os << "  \"" << isot << "\",\n";
    os << "};\n\n";
  }

  // Declare existence of functions
  for (auto& spec : data) {
    for (auto& isot : spec.second) {
      os << "Numeric " << func_name(species_name(spec.first, isot), false)
         << "(Numeric T) noexcept;\n";
    }

    os << '\n';

    for (auto& isot : spec.second) {
      os << "Numeric " << func_name(species_name(spec.first, isot), true)
         << "(Numeric T) noexcept;\n";
    }
  os << '\n';
  }

  // Write the call operator
  for (auto& spec : data) {
    Index i = 0;

    os << "\ntemplate<Derivatives deriv>\nNumeric compute" << spec.first
       << "(Numeric T, const std::string_view isot) {\n";
    for (auto& isot : spec.second) {
      os << "  if (isot == std::get<" << i++ << ">(has" << spec.first
         << ")) {\n"
            "    if constexpr (deriv == Derivatives::No) return "
         << func_name(species_name(spec.first, isot), false)
         << "(T);\n"
            "    if constexpr (deriv == Derivatives::Yes) return "
         << func_name(species_name(spec.first, isot), true) << "(T);\n  }\n\n";
    }
    os << "  if (isot == \"*\") return T * std::numeric_limits<Numeric>::signaling_NaN();\n";
    os << R"--(  ARTS_USER_ERROR("Cannot find ", isot, " for species )--"
       << Species::toShortName(spec.first) << "\")\n";
    os << "}\n";
  }
}

int main(int argn, char** argv) try {
  if (argn < 2)
    throw std::runtime_error(
        "USAGE: PROG INPUT1.XML INPUT2.xml ... INPUTX.xml");
  const std::vector<std::string> xmlfiles{argv + 1, argv + argn};

  for (auto& fn : xmlfiles) {
    auto xmlfile = std::filesystem::canonical(fn);
    make_cc(xmlfile);
  }

  make_h(xmlfiles);

  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Failed:\n" << e.what() << '\n';
  return EXIT_FAILURE;
}