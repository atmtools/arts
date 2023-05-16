#include <iostream>
#include "matpack_data.h"
#include "xml_io_base.h"

int main(int /*argc*/, char * /*argv*/[]) {

  Vector v1{1, 2, 3, 4, 5};
  String filename{"test.xml"};
  xml_write_to_file_base(filename, v1, FILE_TYPE_ASCII, Verbosity());
  std::cout << v1 << '\n';

  Vector v2;
  xml_read_from_file_base(filename, v2, Verbosity());
  std::cout << v2 << '\n';

  return (0);
}
