#include <matpack.h>
#include <xml.h>

template <typename T>
void test_basic(std::string fn) try {
  const auto n = static_cast<T>(-3.14);
  xml_write_to_file_base(fn + fn + "test.xml", n, FileType::ascii);
  T n_read{};
  xml_read_from_file_base(fn + fn + "test.xml", n_read);
  std::println("{} --- {} vs {}", xml_io_stream<T>::type_name, n, n_read);

  if (n != n_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", n, FileType::binary);
  T n_read_bin{};
  xml_read_from_file_base(fn + "testbin.xml", n_read_bin);
  std::println("{} --- {} vs {}", xml_io_stream<T>::type_name, n, n_read_bin);

  if (n != n_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_string(std::string fn) try {
  const String s = "Hello, XML!";
  xml_write_to_file_base(fn + "test.xml", s, FileType::ascii);
  String s_read;
  xml_read_from_file_base(fn + "test.xml", s_read);
  std::println(
      R"({} --- "{}" vs "{}")", xml_io_stream<String>::type_name, s, s_read);

  if (s != s_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", s, FileType::binary);
  String s_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", s_read_bin);
  std::println(R"({} --- "{}" vs "{}")",
               xml_io_stream<String>::type_name,
               s,
               s_read_bin);
  if (s != s_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

template <typename T>
void test_basic_array(std::string fn) try {
  const Array<T> arr{
      static_cast<T>(1.2), static_cast<T>(2.4), static_cast<T>(-3.8)};
  xml_write_to_file_base(fn + "test.xml", arr, FileType::ascii);
  Array<T> arr_read;
  xml_read_from_file_base(fn + "test.xml", arr_read);
  std::println(R"({}<{}> --- {:B,} vs {:B,})",
               xml_io_stream<Array<T>>::type_name,
               xml_io_stream<T>::type_name,
               arr,
               arr_read);

  if (arr != arr_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", arr, FileType::binary);
  Array<T> arr_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", arr_read_bin);
  std::println(R"({}<{}> --- {:B,} vs {:B,})",
               xml_io_stream<Array<T>>::type_name,
               xml_io_stream<T>::type_name,
               arr,
               arr_read_bin);

  if (arr != arr_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_array_string(std::string fn) try {
  const Array<String> s = {"Hello", "my", "XML!"};
  xml_write_to_file_base(fn + "test.xml", s, FileType::ascii);
  Array<String> s_read;
  xml_read_from_file_base(fn + "test.xml", s_read);
  std::println(R"(Array<{}> --- "{}" vs "{}")",
               xml_io_stream<String>::type_name,
               s,
               s_read);

  if (s != s_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", s, FileType::binary);
  Array<String> s_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", s_read_bin);
  std::println(R"(Array<{}> --- "{}" vs "{}")",
               xml_io_stream<String>::type_name,
               s,
               s_read_bin);
  if (s != s_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_vector(std::string fn) try {
  const Vector x{1, 2, 3, 4, 5.5};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  Vector x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(
      R"({} --- {:B,} vs {:B,})", xml_io_stream<Vector>::type_name, x, x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  Vector x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<Vector>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_matrix(std::string fn) try {
  const Matrix x = Vector{1, 2, 3, 4, 5.5, 3.2}.reshape(3, 2);
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  Matrix x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(
      "{} ---\n{:B,}\nvs\n{:B,}", xml_io_stream<Matrix>::type_name, x, x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  Matrix x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println("{} ---\n{:B,}\nvs\n{:B,}",
               xml_io_stream<Matrix>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_tensor3(std::string fn) try {
  const Tensor3 x =
      Vector{1, 2, 3, 4, 5.5, 3.2, 1, 2, 3, 4, 5.5, 3.2}.reshape(3, 2, 2);
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  Tensor3 x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(
      "{} ---\n{:B,}\nvs\n{:B,}", xml_io_stream<Tensor3>::type_name, x, x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  Tensor3 x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println("{} ---\n{:B,}\nvs\n{:B,}",
               xml_io_stream<Tensor3>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_complex_vector(std::string fn) try {
  const ComplexVector x{1 + 2i, 2 + 5i, 3 + 5i, 4 + 12i, 5.5 - 32i};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  ComplexVector x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<ComplexVector>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  ComplexVector x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<ComplexVector>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_complex_matrix(std::string fn) try {
  const ComplexMatrix x =
      ComplexVector{1 + 2i, 2 + 5i, 3 + 5i, 4 + 12i, 5.5 - 32i, 32 - 90i}
          .reshape(3, 2);
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  ComplexMatrix x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println("{} ---\n{:B,}\nvs\n{:B,}",
               xml_io_stream<ComplexMatrix>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  ComplexMatrix x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println("{} ---\n{:B,}\nvs\n{:B,}",
               xml_io_stream<ComplexMatrix>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_complex_tensor3(std::string fn) try {
  const ComplexTensor3 x = ComplexVector{1 + 2i,
                                         2 + 5i,
                                         3 + 5i,
                                         4 + 12i,
                                         5.5 - 32i,
                                         32 - 90i,
                                         1 + 2i,
                                         2 + 5i,
                                         3 + 5i,
                                         4 + 12i,
                                         5.5 - 32i,
                                         32 - 90i}
                               .reshape(3, 2, 2);
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  ComplexTensor3 x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println("{} ---\n{:B,}\nvs\n{:B,}",
               xml_io_stream<ComplexTensor3>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  ComplexTensor3 x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println("{} ---\n{:B,}\nvs\n{:B,}",
               xml_io_stream<ComplexTensor3>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_vector2(std::string fn) try {
  const Vector2 x{1, 2};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  Vector2 x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(
      R"({} --- {:B,} vs {:B,})", xml_io_stream<Vector2>::type_name, x, x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  Vector2 x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<Vector2>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_matrix44(std::string fn) try {
  using Matrix44 = matpack::cdata_t<Numeric, 4, 4>;
  const Matrix44 x{1, 2, 3, 4, 5, 6, 7, 8, 9, 321, 312, 31, 23, 123, 12};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  Matrix44 x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<Matrix44>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  Matrix44 x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<Matrix44>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_array_matrix44(std::string fn) try {
  using Matrix44 = matpack::cdata_t<Numeric, 4, 4>;
  const Matrix44 v{1, 2, 3, 4, 5, 6, 7, 8, 9, 321, 312, 31, 23, 123, 12};
  Array<Matrix44> x  = {v, v, v};
  x[0][0, 0]        += 1;
  x[2][2, 0]        += 1;
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  Array<Matrix44> x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<Array<Matrix44>>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  Array<Matrix44> x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<Array<Matrix44>>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_complex_vector2(std::string fn) try {
  using ComplexVector2 = matpack::cdata_t<Complex, 2>;
  const ComplexVector2 x{1 - 3i, 2 + 8i};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  ComplexVector2 x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(
      R"({} --- {:B,} vs {:B,})", xml_io_stream<Vector2>::type_name, x, x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  ComplexVector2 x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<Vector2>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_complex_matrix44(std::string fn) try {
  using ComplexMatrix44 = matpack::cdata_t<Complex, 4, 4>;
  const ComplexMatrix44 x{
      1 - 6i, 2 + 8i, 3 - 99i, 4, 5, 6, 7, 8, 9, 321, 312, 31, 23, 123, 12};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  ComplexMatrix44 x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<ComplexMatrix44>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  ComplexMatrix44 x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<ComplexMatrix44>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_string_vector(std::string fn) try {
  const matpack::data_t<String, 1> x{
      "Hello"s, "World"s, "how"s, "are"s, "you"s, "doing"s};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  matpack::data_t<String, 1> x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<matpack::data_t<String, 1>>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  matpack::data_t<String, 1> x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<matpack::data_t<String, 1>>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}
ARTS_METHOD_ERROR_CATCH

void test_ascending_grid(std::string fn) {
  const AscendingGrid grid{0, 1, 2, 3, 4, 5};
  xml_write_to_file_base(fn + "test.xml", grid, FileType::ascii);
  AscendingGrid grid_read;
  xml_read_from_file_base(fn + "test.xml", grid_read);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<AscendingGrid>::type_name,
               grid,
               grid_read);
  if (grid != grid_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", grid, FileType::binary);
  AscendingGrid grid_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", grid_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<AscendingGrid>::type_name,
               grid,
               grid_read_bin);
  if (grid != grid_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}

void test_descending_grid(std::string fn) {
  const DescendingGrid grid{5, 4, 3, 2, 1};
  xml_write_to_file_base(fn + "test.xml", grid, FileType::ascii);
  DescendingGrid grid_read;
  xml_read_from_file_base(fn + "test.xml", grid_read);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<DescendingGrid>::type_name,
               grid,
               grid_read);
  if (grid != grid_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", grid, FileType::binary);
  DescendingGrid grid_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", grid_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<DescendingGrid>::type_name,
               grid,
               grid_read_bin);
  if (grid != grid_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}

template <typename T>
auto grid() {
  if constexpr (std::same_as<T, AscendingGrid>) {
    return AscendingGrid{0, 1, 2, 3, 4};
  } else if constexpr (std::same_as<T, DescendingGrid>) {
    return DescendingGrid{4, 3, 2, 1, 0};
  } else if constexpr (std::same_as<T, ArrayOfString>) {
    return ArrayOfString{"A", "B", "C", "D", "E"};
  } else {
    return Vector{0, 1, 2, 3, 4};
  }
}

template <class Grid>
void test_gridded_field(std::string fn) {
  const matpack::gridded_data_t<Numeric, Grid> x{
      .data_name  = "TestField",
      .data       = {1.0, 2.0, 3.0, 4.0, 5.0},
      .grid_names = {"Grid1"},
      .grids      = grid<Grid>()};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  matpack::gridded_data_t<Numeric, Grid> x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<matpack::gridded_data_t<Numeric, Grid>>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  matpack::gridded_data_t<Numeric, Grid> x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println(R"({} --- {:B,} vs {:B,})",
               xml_io_stream<matpack::gridded_data_t<Numeric, Grid>>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}

template <typename... Ts>
void test_variant(auto&& data, std::string fn) {
  const std::variant<Ts...> x{std::move(data)};
  xml_write_to_file_base(fn + "test.xml", x, FileType::ascii);
  std::variant<Ts...> x_read;
  xml_read_from_file_base(fn + "test.xml", x_read);
  std::println("{} ---\n{:B,}\nvs\n{:B,}",
               xml_io_stream<std::variant<Ts...>>::type_name,
               x,
               x_read);
  if (x != x_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", x, FileType::binary);
  std::variant<Ts...> x_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", x_read_bin);
  std::println("{} ---\n{:B,}\nvs\n{:B,}",
               xml_io_stream<std::variant<Ts...>>::type_name,
               x,
               x_read_bin);
  if (x != x_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}

void test_unordered_map(std::string fn) {
  const std::unordered_map<String, Numeric> map{
      {"A", 1.0}, {"B", 2.0}, {"C", 3.0}, {"D", 4.0}};
  xml_write_to_file_base(fn + "test.xml", map, FileType::ascii);
  std::unordered_map<String, Numeric> map_read;
  xml_read_from_file_base(fn + "test.xml", map_read);
  for (auto& [k, v] : map) {
    std::println("{} {} {}", k, map.at(k), map_read.at(k));
  }
  if (map != map_read) {
    throw std::runtime_error("Read from XML does not match original.");
  }

  xml_write_to_file_base(fn + "testbin.xml", map, FileType::binary);
  std::unordered_map<String, Numeric> map_read_bin;
  xml_read_from_file_base(fn + "testbin.xml", map_read_bin);
  for (auto& [k, v] : map) {
    std::println("{} {} {}", k, map.at(k), map_read_bin.at(k));
  }
  if (map != map_read_bin) {
    throw std::runtime_error("Read from XML does not match original.");
  }
}

int main() try {
  test_basic<Index>("index");
  test_basic<Size>("size");
  test_basic<Numeric>("numeric");
  test_string("string");
  test_basic_array<Index>("aindex");
  test_basic_array<Size>("asize");
  test_basic_array<Numeric>("anumeric");
  test_array_string("astring");
  test_vector("vector");
  test_matrix("matrix");
  test_tensor3("tensor3");
  test_complex_vector("cvector");
  test_complex_matrix("cmatrix");
  test_complex_tensor3("ctensor3");
  test_vector2("vector2");
  test_matrix44("matrix44");
  test_array_matrix44("amatrix44");
  test_complex_vector2("cvector2");
  test_complex_matrix44("cmatrix44");
  test_string_vector("svector");
  test_ascending_grid("ascend");
  test_descending_grid("descend");
  test_gridded_field<Vector>("gfvec");
  test_gridded_field<AscendingGrid>("gfasc");
  test_gridded_field<DescendingGrid>("gfdesc");
  test_gridded_field<ArrayOfString>("gfastr");
  test_variant<ArrayOfString, Vector, Matrix>(grid<Vector>(), "varvec");
  test_variant<ArrayOfString, Vector, Matrix>(grid<ArrayOfString>(), "varastr");
  test_variant<ArrayOfString, Vector, Matrix>(Vector{1, 2, 3, 4}.reshape(2, 2),
                                              "varmat");
  test_unordered_map("unmap");
} catch (const std::runtime_error& e) {
  std::println("Error:\n{}", e.what());
  return EXIT_FAILURE;
}

static_assert(std::is_trivially_copyable_v<Numeric>);
static_assert(std::is_trivially_copyable_v<Vector2>);
static_assert(std::is_trivially_copyable_v<Index>);
static_assert(std::is_trivially_copyable_v<Size>);
static_assert(std::is_trivially_copyable_v<Complex>);
static_assert(not std::is_trivially_copyable_v<String>);
static_assert(not std::is_trivially_copyable_v<std::vector<Numeric>>);
