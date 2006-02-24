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
      xml_write_to_file ("sdata1.xml", species_data);
      cout << "Wrote species_data: " << endl;

      species_data.clear ();

      xml_read_from_file ("sdata1.xml", species_data);
      cout << "Read species_data: " << endl;

      xml_write_to_file ("sdata2.xml", species_data);
      cout << "Wrote species_data: " << endl;
    }
  catch (runtime_error e)
    {
      cerr << e.what ();
    }

  return (0);
}
