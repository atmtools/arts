using namespace std;

#include <iostream>
#include "matpackII.h"
#include "xml_io.h"
#include "exceptions.h"
#include "absorption.h"

extern Array<SpeciesRecord> species_data;

int
main (int /* argc */, char * /* argv */ [])
{
  define_species_data ();
  try
    {
      xml_write_to_file ("sdata.xml", species_data);
      cout << "Wrote species_data: " << endl;

      //xml_read_from_file ("sparse.xml", b);
      //cout << "Read: " << b << endl;
    }
  catch (runtime_error e)
    {
      cerr << e.what ();
    }

  return (0);
}
