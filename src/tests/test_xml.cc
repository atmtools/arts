#include <matpack.h>
#include <xml.h>

#include <iostream>

int main(int /*argc*/, char* /*argv*/[]) {
  Vector v1{1, 2, 3, 4, 5};
  String filename{"test.xml"};
  xml_write_to_file_base(filename, v1, FileType::ascii);
  arts::println("{}", v1);

  Vector v2;
  xml_read_from_file_base(filename, v2);
  arts::println("{}", v2);

  return (0);
}
