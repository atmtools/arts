#include <math_funcs.h>
#include <mystring.h>
#include <partfun.h>
#include <partition_function_data.h>
#include <workspace.h>
#include <xml_io.h>
#include <xml_io_partfun.h>

#include <filesystem>

void WriteBuiltinPartitionFunctionsXML(
    const String& fileformat, const String& dir, const Numeric& Tlow, const Numeric& Tupp, const Index& N) {
  if (Tupp <= Tlow) arts::user_error("Need a range [low, high], has [{}, {}]", Tlow, Tupp);
  if (N < 2) arts::user_error("Need a positive step counter 2 or larger, has: {}", N);

  const auto d = std::filesystem::path(dir.c_str());
  if (not std::filesystem::is_directory(d)) arts::user_error("dir: {} is not a directory", dir);

  const Vector T = [&] {
    Vector x;
    nlinspace(x, Tlow, Tupp, N);
    return x;
  }();
  const FileType ftype = to<FileType>(fileformat);

  for (auto& ir : Species::Isotopologues) {
    if (PartitionFunctions::has_partfun(ir)) {
      PartitionFunctionsData data{PartitionFunctionsType::Interp, Matrix(N, 2)};
      for (Index i = 0; i < N; i++) {
        data.data[i, 0] = T[i];
        data.data[i, 1] = PartitionFunctions::Q(T[i], ir);
      }

      xml_write_to_file((d / (ir.FullName() + ".xml")).string(), data, ftype, 0);
    }
  }
}
