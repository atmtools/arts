#include <cstdlib>
#include <iostream>
#include "bifstream.h"
#include "matpack_data.h"

//#define VERBOSE

typedef long PPHeader[65];

void readppheader(bifstream& is, PPHeader& pph);
void readppdata(bifstream& is, PPHeader& pph, Vector& v);

void readppheader(bifstream& is, PPHeader& pph) {
  for (int j = 0; j < 65 && is.good(); j++) {
    pph[j] = is.readInt(4);
#ifdef VERBOSE
    std::cout << j << " " << pph[j] << " " << is.error() << std::endl;
#endif
  }
}

void readppdata(bifstream& is, PPHeader& pph, Vector& v) {
  const int EXTRA_DATA = 3;
  v.resize(pph[15] + EXTRA_DATA);

  for (int j = 0; j < pph[15] + EXTRA_DATA && is.good(); j++) {
    v[j] = is.readFloat(binio::Single);
  }

  if (!is.good()) {
    std::cerr << "Error: " << is.error() << std::endl;
  }
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " PPFILE [MAXFIELDS]" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  long maxfields = -1;
  if (argc > 2) maxfields = strtol(argv[2], NULL, 10);

  bifstream is(argv[1]);
  if (is.good()) {
    PPHeader pph;
    Vector v;
    long field = 0;

    is.setFlag(binio::BigEndian, true);
    while (is.good() && field != maxfields) {
      binio::Error e;
      field++;

      // Check if more data is available in file
      is.peekInt(4);
      if (is.error() & 32) return (EXIT_FAILURE);

      readppheader(is, pph);
      if ((e = is.error())) {
        std::cerr << "Reading " << field << ". header failed with error " << e
             << std::endl;
        std::exit(EXIT_FAILURE);
      } else {
        std::cout << "Field # " << std::setw(5) << field << " -- STASH code " << std::setw(5)
             << pph[42] << " -- PP code " << std::setw(4) << pph[23] << " -- PP VCT "
             << std::setw(3) << pph[26] << std::endl;
      }

      readppdata(is, pph, v);
      if ((e = is.error())) {
        std::cerr << "Reading " << field << ". data failed with error " << e << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  } else {
    std::cerr << "Error reading from " << argv[1] << std::endl;
    return (EXIT_FAILURE);
  }

  return (EXIT_SUCCESS);
}
