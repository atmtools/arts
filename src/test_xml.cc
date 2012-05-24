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

#include <iostream>
#include "arts.h"
#include "matpackII.h"
#include "xml_io.h"
#include "exceptions.h"
#include "absorption.h"

extern Array<SpeciesRecord> species_data;

int
main (int /*argc*/, char * /*argv*/ [])
{
  define_species_data ();
  try
    {
      xml_write_to_file ("sdata1.xml", species_data, FILE_TYPE_ASCII, Verbosity());
      cout << "Wrote species_data: " << endl;

      species_data.clear ();

      xml_read_from_file ("sdata1.xml", species_data, Verbosity());
      cout << "Read species_data: " << endl;

      xml_write_to_file ("sdata2.xml", species_data, FILE_TYPE_ASCII, Verbosity());
      cout << "Wrote species_data: " << endl;
    }
  catch (runtime_error e)
    {
      cerr << e.what ();
    }

  return (0);
}
