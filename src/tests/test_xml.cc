#include <iostream>

#include <matpack.h>
#include <xml.h>

int main(int /*argc*/, char* /*argv*/[]) {
  Vector v1{1, 2, 3, 4, 5};
  String filename{"test.xml"};
  xml_write_to_file_base(filename, v1, FileType::ascii);
  std::cout << std::format("{}", v1) << '\n';

  Vector v2;
  xml_read_from_file_base(filename, v2);
  std::cout << std::format("{}", v2) << '\n';

  return (0);
}
