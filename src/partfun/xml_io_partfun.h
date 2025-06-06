#ifndef xml_io_partfun_h
#define xml_io_partfun_h

#include <filesystem>

#include "partition_function_data.h"
namespace PartitionFunctions {
Data data_read_file(const std::filesystem::path& path);
}  // namespace PartitionFunctions

#endif  // xml_io_partfun_h
