#include <cstdlib>
#include <iostream>

#include "matpack.h"
#include "xml_io.h"

int main(int /* argc */, char* /* argv */[]) {
  // Create binary file
  Tensor4 v(4, 4, 4, 4);

  for (Index i = 0; i < 4; i++)
    for (Index j = 0; j < 4; j++)
      for (Index k = 0; k < 4; k++)
        for (Index l = 0; l < 4; l++)
          v[i, j, k, l] = double(i * 4 * 4 * 4 + j * 4 * 4 + k * 4 + l);

  xml_write_to_file("outfile.xml", v, FileType::binary, 0);

  // Read binary file
  Tensor4 w;

  xml_read_from_file("outfile.xml", w);

  std::cout << std::format("{}", w) << '\n';

  return (EXIT_SUCCESS);
}
