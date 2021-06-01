#include "debug.h"
#include "xml_io_partfun.h"

#include <pugixml.hpp>

namespace PartitionFunctions {
void Data::print_data() const {
  constexpr int cutline = 10;
  const Index n = data.nrows();
  
  switch (type) {
    case PartitionFunctions::Type::Interp:
      std::cout << "constexpr std::array<Numeric, " << n << "> data{";
      for (Index i=0; i<n; i++) {
        if (i % cutline == 0) {
          std::cout << '\n';
        }
        std::cout << data(i, 0) << ',' << ' ';
      }
      std::cout << "};\n";
      
      std::cout << "constexpr std::array<Numeric, " << n << "> grid{";
      for (Index i=0; i<n; i++) {
        if (i % cutline == 0)  {
          std::cout << '\n';
        }
        std::cout << data(i, 1) << ',' << ' ';
      }
      std::cout << "};\n";
      break;
    case PartitionFunctions::Type::Coeff:
      std::cout << "constexpr std::array<Numeric, " << n << "> coeff{";
      for (Index i=0; i<n; i++) {
        if (i % cutline == 0)  {
          std::cout << '\n';
        }
        std::cout << data(i, 0) << ',' << ' ';
      }
      std::cout << "};\n";
      break;
    case PartitionFunctions::Type::FINAL: {/* leave last
      */ }
  }
}

void Data::print_method() const {
  switch (type) {
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
  
Data data_read_file(const std::filesystem::path& path) {
  Data out;
  
  pugi::xml_document doc;
  
  // Must start with ARTS tag
  doc.load_file(path.c_str());
  pugi::xml_node root = doc.document_element();
  ARTS_USER_ERROR_IF(std::string_view(root.name()) not_eq "arts", "Expects \"arts\" as root of XML document, got: ", root.name(), " in file: ", path)
  
  // Must continue with partition functions data, this must have a type understandable by ARTS
  pugi::xml_node pfd = root.first_child();
  ARTS_USER_ERROR_IF(std::string_view(pfd.name()) not_eq "PartitionFunctionsData", "Expects \"PartitionFunctionsData\" as first type of arts-xml, got: ", root.name(), " in file: ", path)
  out.type = toTypeOrThrow(std::string_view(pfd.attribute("type").as_string()));
  
  // FIXME: from here, this should use a proper matrix-reader
  
  // Must continue with partition functions data
  pugi::xml_node mat = pfd.first_child();
  ARTS_USER_ERROR_IF(std::string_view(mat.name()) not_eq "Matrix", "Expects \"Matrix\" as first type of arts-xml, got: ", root.name(), " in file: ", path)
  const Index nr = Index(mat.attribute("nrows").as_llong());
  const Index nc = Index(mat.attribute("ncols").as_llong());
  
  switch (out.type) {
    case PartitionFunctions::Type::Interp:
      ARTS_USER_ERROR_IF(nc not_eq 2, "Must have form: TEMP DATA", " in file: ", path)
      ARTS_USER_ERROR_IF(nr < 2, "Must have some data", " in file: ", path)
      break;
    case PartitionFunctions::Type::Coeff:
      ARTS_USER_ERROR_IF(nc not_eq 1, "Must have form: COEFF", " in file: ", path)
      ARTS_USER_ERROR_IF(nr < 1, "Must have some data", " in file: ", path)
      break;
    case PartitionFunctions::Type::FINAL: {/* leave last
      */ }
  }
  
  // Read all the data
  out.data.resize(nr, nc);
  std::istringstream is_data{mat.text().as_string()};
  for (Index i=0; i<nr; i++) {
    for (Index j=0; j<nc; j++) {
      is_data >> out.data(i, j);
    }
  }
  
  return out;
}
}
