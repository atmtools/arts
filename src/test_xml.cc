using namespace std;

#include <iostream>
#include "matpackII.h"
#include "xml_io.h"
#include "exceptions.h"

int
main (int /* argc */, char * /* argv */ [])
{
  Sparse a (5, 3);
  Sparse b;

  try
    {
      a (1, 1) = 6.;
      a (2, 2) = 5.;

      xml_write_to_file ("sparse.xml", a);
      cout << "Wrote: " << a << endl;

      xml_read_from_file ("sparse.xml", b);
      cout << "Read: " << b << endl;
    }
  catch (runtime_error e)
    {
      cerr << e.what ();
    }

  return (0);
}
