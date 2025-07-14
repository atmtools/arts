#ifndef xml_io_partfun_h
#define xml_io_partfun_h

#include <xml.h>

#include <filesystem>

#include "partition_function_data.h"
namespace PartitionFunctions {
Data data_read_file(const std::filesystem::path& path);
}  // namespace PartitionFunctions

template <>
struct xml_io_stream<PartitionFunctionsData> {
  constexpr static std::string_view type_name   = "PartitionFunctionsData"sv;

  static void write(std::ostream&,
                    const PartitionFunctionsData&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv);

  static void read(std::istream&,
                   PartitionFunctionsData&,
                   bifstream* = nullptr);
};

#endif  // xml_io_partfun_h
