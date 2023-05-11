#include <filesystem>

#include "math_funcs.h"
#include "matpack_concepts.h"
#include "mystring.h"
#include "partfun.h"

#include "xml_io.h"

void WriteBuiltinPartitionFunctionsXML(
  const String& fileformat,
  const String& dir,
  const Numeric& Tlow,
  const Numeric& Tupp,
  const Index& N) {
  ARTS_USER_ERROR_IF(Tupp <= Tlow, "Need a range [low, high], has [", Tlow, ", ", Tupp, "]")
  ARTS_USER_ERROR_IF(N < 2, "Need a positive step counter 2 or larger, has: ", N)
  
  const auto d = std::filesystem::path(dir.c_str());
  ARTS_USER_ERROR_IF(not std::filesystem::is_directory(d), "dir: ", dir, " is not a directory")
  
  const Vector T = [&]{Vector x; nlinspace(x, Tlow, Tupp, N); return x;}();
  const FileType ftype = string2filetype(fileformat);
  
  for (auto& ir: Species::Isotopologues) {
    if (PartitionFunctions::has_partfun(ir)) {
      
      PartitionFunctionsData data{PartitionFunctions::Type::Interp, Matrix(N, 2)};
      for (Index i=0; i<N; i++) {
        data.data(i, 0) = T[i];
        data.data(i, 1) = PartitionFunctions::Q(T[i], ir);
      }
      
      xml_write_to_file(String((d / (ir.FullName() + ".xml").c_str()).native()),
                        data, ftype, 0);
    }
  }
}
