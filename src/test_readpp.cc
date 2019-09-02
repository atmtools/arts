/* Copyright (C) 2004-2012 Oliver Lemke <olemke@core-dump.info>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

#include <cstdlib>
#include <iostream>
#include "bifstream.h"
#include "matpackI.h"

//#define VERBOSE

typedef long PPHeader[65];

void readppheader(bifstream& is, PPHeader& pph);
void readppdata(bifstream& is, PPHeader& pph, Vector& v);

void readppheader(bifstream& is, PPHeader& pph) {
  for (int j = 0; j < 65 && is.good(); j++) {
    pph[j] = is.readInt(4);
#ifdef VERBOSE
    cout << j << " " << pph[j] << " " << is.error() << endl;
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
    cerr << "Error: " << is.error() << endl;
  }
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " PPFILE [MAXFIELDS]" << endl;
    exit(EXIT_FAILURE);
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
        cerr << "Reading " << field << ". header failed with error " << e
             << endl;
        exit(EXIT_FAILURE);
      } else {
        cout << "Field # " << setw(5) << field << " -- STASH code " << setw(5)
             << pph[42] << " -- PP code " << setw(4) << pph[23] << " -- PP VCT "
             << setw(3) << pph[26] << endl;
      }

      readppdata(is, pph, v);
      if ((e = is.error())) {
        cerr << "Reading " << field << ". data failed with error " << e << endl;
        exit(EXIT_FAILURE);
      }
    }
  } else {
    cerr << "Error reading from " << argv[1] << endl;
    return (EXIT_FAILURE);
  }

  return (EXIT_SUCCESS);
}
