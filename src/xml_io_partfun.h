#ifndef xml_io_partfun_h
#define xml_io_partfun_h

#include <filesystem>

#include "matpackI.h"
#include "template_partfun.h"

namespace PartitionFunctions {
struct Data {
  Type type;
  Matrix data;
  
  void print_data() const;
  
  void print_method() const;
};

Data data_read_file(const std::filesystem::path& path);
}

using PartitionFunctionsData = PartitionFunctions::Data;

#endif  // xml_io_partfun_h
