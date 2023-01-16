#ifndef xml_io_partfun_h
#define xml_io_partfun_h

#include <filesystem>

#include "matpack_data.h"
#include "template_partfun.h"

namespace PartitionFunctions {
Data data_read_file(const std::filesystem::path& path);
}

#endif  // xml_io_partfun_h
