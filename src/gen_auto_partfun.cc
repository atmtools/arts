#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <vector>

#include <pugixml.hpp>

#include "debug.h"
#include "gridded_fields.h"
#include "species.h"
#include "template_partfun.h"

#include "xml_io_partfun.h"

//! Contains the partition function data in a simple format
struct PrintData {
  String name;
  PartitionFunctions::Data data;
};


std::pair<std::string_view, std::string_view> split_filename_species(const std::string_view isot) {
  auto x = std::min(isot.find('-'), isot.size());
  return {isot.substr(0, x),
          isot.substr(std::min(x+1, isot.size()), isot.size())};
}


std::pair<Species::Species, PrintData> gen_fil(std::filesystem::path fil) {
  if (not std::filesystem::exists(fil)) {
    ARTS_USER_ERROR("file does not exist: ", fil)
  }
  
  PrintData out;
  
  auto filname_no_ext = fil.filename().replace_extension().native();
  const auto [specname, isotname] = split_filename_species(filname_no_ext);
  const Species::Species spec = Species::fromShortName(specname);
  out.name = isotname;
  ARTS_USER_ERROR_IF(not good_enum(spec), "Got Species: ", specname, "\nFrom isotope: ", out.name, "\nIn filename of: ", fil)
  
  out.data = PartitionFunctions::data_read_file(fil);
  
  return {spec, out};
}


std::map<Species::Species, std::vector<PrintData>> gen_dir(std::filesystem::path dir) {
  if (not std::filesystem::is_directory(dir)) {
    ARTS_USER_ERROR("not a directory: ", dir)
  }
  
  // Make sure all species exists in the data-map
  std::map<Species::Species, std::vector<PrintData>> out;
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


void print_data(const PrintData& data) {
  data.data.print_data();
}


void print_method(const PrintData& data) {
  data.data.print_method();
}


void print_auto_partfun_h(const std::map<Species::Species, std::vector<PrintData>>& data) {
  std::puts(R"AUTO_PARTFUN(//! Auto-generated partition function data
#ifndef auto_partfun_h
#define auto_partfun_h

#include <limits>

#include "template_partfun.h"
#include "debug.h"

template <std::size_t N>
std::string column_string(const std::array<std::string_view, N>& list) {
  if constexpr (static_cast<bool>(N)) {
    std::ostringstream os;
    for (auto& x: list) os << x << '\n';
    return os.str();
  } else return "Nothing Was Compiled For This Species";
}

namespace PartitionFunctions {)AUTO_PARTFUN");
  
  std::cout << std::setprecision(std::numeric_limits<Numeric>::digits10 + 1);
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
