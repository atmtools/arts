#include <filesystem>
#include <fstream>
#include <map>
#include <vector>

#include <pugixml.hpp>

#include "debug.h"
#include "gridded_fields.h"
#include "species.h"
#include "template_partfun.h"

//! Contains the partition function data in a simple format
struct PartitionFunctionData {
  PartitionFunctions::Type typ;
  String name;
  Vector T;
  Vector Q;
};


std::pair<std::string_view, std::string_view> split_filename_species(const std::string_view isot) {
  auto x = std::min(isot.find('-'), isot.size());
  return {isot.substr(0, x),
          isot.substr(std::min(x+1, isot.size()), isot.size())};
}


std::pair<Species::Species, PartitionFunctionData> gen_fil(std::filesystem::path fil) {
  if (not std::filesystem::exists(fil)) {
    ARTS_USER_ERROR("file does not exist: ", fil)
  }
  
  PartitionFunctionData out;
  
  auto filname_no_ext = fil.filename().replace_extension().native();
  const auto [specname, isotname] = split_filename_species(filname_no_ext);
  const Species::Species spec = Species::fromShortName(specname);
  out.name = isotname;
  ARTS_USER_ERROR_IF(not good_enum(spec), "Got Species: ", spec, "\nFrom isotope: ", out.name, "\nIn filename of: ", fil)
  
  pugi::xml_document doc;
  doc.load_file(fil.c_str());
  pugi::xml_node root = doc.document_element();
  ARTS_USER_ERROR_IF(std::string_view(root.name()) not_eq "arts", "Expects \"arts\" as root of XML document, got: ", root.name())
  
  pugi::xml_node gf1 = root.first_child();
  ARTS_USER_ERROR_IF(std::string_view(gf1.name()) not_eq "GriddedField1", "Expects \"GriddedField1\" as first type of arts-xml, got: ", root.name())
  out.typ = PartitionFunctions::toTypeOrThrow(std::string_view(gf1.attribute("name").as_string()));
  
  pugi::xml_node grid = gf1.first_child();
  ARTS_USER_ERROR_IF(std::string_view(grid.name()) not_eq "Vector", "Expects \"Vector\" as first type of GriddedField1 document, got: ", grid.name())
  
  pugi::xml_node data = grid.next_sibling();
  ARTS_USER_ERROR_IF(std::string_view(data.name()) not_eq "Vector", "Expects \"Vector\" as second type of GriddedField1 document, got: ", data.name())
  ARTS_USER_ERROR_IF(std::string_view(data.attribute("name").as_string()) not_eq "Data", "Expects name=\"Data\", got: ", data.attribute("name").as_string())
  
  const Index ng = grid.attribute("nelem").as_llong();
  const Index nd = data.attribute("nelem").as_llong();
  switch (out.typ) {
    case PartitionFunctions::Type::Interp:
      ARTS_USER_ERROR_IF(nd not_eq ng or nd < 2, "Expects the two vectors to have the same length")
      ARTS_USER_ERROR_IF(std::string_view(grid.attribute("name").as_string()) not_eq "Temperature", "Expects name=\"Temperature\", got: ", grid.attribute("name").as_string())
      break;
    case PartitionFunctions::Type::Coeff:
      ARTS_USER_ERROR_IF(nd < 1, "Expects the some values")
      ARTS_USER_ERROR_IF(ng not_eq nd, "Expects the same values")
      ARTS_USER_ERROR_IF(std::string_view(grid.attribute("name").as_string()) not_eq "Coeff", "Expects name=\"Coeff\", got: ", grid.attribute("name").as_string())
      break;
    case PartitionFunctions::Type::FINAL: {/* leave last
      */ }
  }
  
  out.Q.resize(nd);
  std::istringstream is_data{data.text().as_string()};
  for (Index i=0; i<nd; i++) {
    is_data >> out.Q[i];
  }
  
  out.T.resize(ng);
  std::istringstream is_grid{grid.text().as_string()};
  for (Index i=0; i<ng; i++) {
    is_grid >> out.T[i];
  }
  
  return {spec, out};
}


std::map<Species::Species, std::vector<PartitionFunctionData>> gen_dir(std::filesystem::path dir) {
  if (not std::filesystem::is_directory(dir)) {
    ARTS_USER_ERROR("not a directory: ", dir)
  }
  
  // Make sure all species exists in the data-map
  std::map<Species::Species, std::vector<PartitionFunctionData>> out;
  for (auto spec: Species::enumtyps::SpeciesTypes) {
    out[spec] = {};
  }
  
  // Fill out the data map per file in the directory
  for(auto& p: std::filesystem::directory_iterator(dir)) {
    if (p.path().extension() == ".xml") {
      const auto [spec, data] = gen_fil(p.path());
      out.at(spec).push_back(data);
    }
  }
  for (auto spec: Species::enumtyps::SpeciesTypes) {
    std::sort(out[spec].begin(), out[spec].end(), [](auto& a, auto& b){return a.name.compare(b.name) < 0;});
  }
  
  return out;
}


void print_data(const PartitionFunctionData& data) {
  constexpr int cutline = 10;
  
  switch (data.typ) {
    case PartitionFunctions::Type::Interp:
      std::cout << "constexpr std::array<Numeric, " << data.T.nelem() << "> data{";
      for (Index i=0; i<data.T.nelem(); i++) {
        if (i % cutline == 0) std::puts("");
        std::cout << data.Q[i] <<','<<' ';
      }
      std::puts("};");
      
      std::cout << "constexpr std::array<Numeric, " << data.T.nelem() << "> grid{";
      for (Index i=0; i<data.T.nelem(); i++) {
        if (i % cutline == 0) std::puts("");
        std::cout << data.T[i] <<','<<' ';
      }
      std::puts("};");
      break;
    case PartitionFunctions::Type::Coeff:
      std::cout << "constexpr std::array<Numeric, " << data.T.nelem() << "> coeff{";
      for (Index i=0; i<data.T.nelem(); i++) {
        if (i % cutline == 0) std::puts("");
        std::cout << data.Q[i] <<','<<' ';
      }
      std::puts("};");
      break;
    case PartitionFunctions::Type::FINAL: {/* leave last
      */ }
  }
}


void print_method(const PartitionFunctionData& data) {
  switch (data.typ) {
    case PartitionFunctions::Type::Interp:
      std::cout << "return linterp<derivative>(grid, data, T);\n";
      break;
    case PartitionFunctions::Type::Coeff:
      std::cout << "return polynom<derivative>(coeff, T);\n";
      break;
    case PartitionFunctions::Type::FINAL: {/* leave last
      */ }
  }
}


void print_auto_partfun_h(const std::map<Species::Species, std::vector<PartitionFunctionData>>& data) {
  std::puts(R"AUTO_PARTFUN(//! Auto-generated partition function data
#ifndef auto_partfun_h
#define auto_partfun_h

#include "template_partfun.h"
#include "debug.h"

template <std::size_t N>
std::string column_string(const std::array<std::string_view, N>& list) {
  if constexpr (N) {
    std::ostringstream os;
    for (auto& x: list) os << x << '\n';
    return os.str();
  } else return "Nothing Was Compiled For This Species";
}

namespace PartitionFunctions {)AUTO_PARTFUN");
  
  std::cout << std::setprecision(15);
  for (auto& [spec, vec_data]: data) {
    
    std::cout << "constexpr std::array<std::string_view, " << vec_data.size() << "> has"<<spec<<"{\n";
    for (auto& pfdata: vec_data) {
      std::cout << "\"" << pfdata.name << "\",\n";
    }
    std::cout << "};\n";
    
    std::cout << '\n';
    std::cout << "template <Derivatives derivative>\nconstexpr Numeric compute"<<spec<<"(const Numeric T, const std::string_view name) {\n";
    
    bool first = true;
    for (auto& pfdata: vec_data) {
      if (not first) {
        std::cout << "else if";
      } else {
        std::cout << "if";
        first = false;
      }
      
      std::cout << " (\"" << pfdata.name << "\" == name) {\n";
      print_data(pfdata);
      print_method(pfdata);
      std::cout << "} ";
    }
    
    if (not first) {
      std::cout << "else if";
    } else {
      std::cout << "if";
    }
    
    std::cout << " (\"*\" == name) {\n"
                 "return T * std::numeric_limits<Numeric>::signaling_NaN();\n"
                 "} else ARTS_USER_ERROR(\"Cannot understand species: \", name, \" of "
              << spec << "\\n\"\n"
                 "                       \"Did you compile without supporting it?\\n\"\n"
                 "                       \"Supported species are:\\n\", column_string(has" << spec << "))"
                 "\n}\n\n";
  }
  
  std::puts(R"AUTO_PARTFUN(
}  // PartitionFunctions

#endif  // auto_partfun_h
)AUTO_PARTFUN");
}


int main(int argc, char ** argv) try {
  ARTS_USER_ERROR_IF(argc not_eq 2, "Bad call, must be used as ", argv[0], " DIR")
  
  //! All calculations are done on the directory data
  print_auto_partfun_h(gen_dir(argv[1]));
  
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
