/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

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

#include "arts.h"
#include "matpack_data.h"
#include "xml_io.h"

int main(int /* argc */, char* /* argv */[]) {
  // Create binary file
  Tensor4 v(4, 4, 4, 4);

  for (Index i = 0; i < 4; i++)
    for (Index j = 0; j < 4; j++)
      for (Index k = 0; k < 4; k++)
        for (Index l = 0; l < 4; l++)
          v(i, j, k, l) = double(i * 4 * 4 * 4 + j * 4 * 4 + k * 4 + l);

  xml_write_to_file("outfile.xml", v, FILE_TYPE_BINARY, 0, Verbosity());

  // Read binary file
  Tensor4 w;

  xml_read_from_file("outfile.xml", w, Verbosity());

  cout << w << endl;

  return (EXIT_SUCCESS);
}
