using namespace std;

#include <iostream>
#include <string>
#include "xml_io.h"
#include "exceptions.h"

int
main (int /* argc */, char * /* argv */ [])
{
  String str;

  try
    {
      cout << "Testing StringWriteXML" << endl;
      str = "Hello World";
      xml_write_to_file ("string.xml", str);
      cout << "Wrote: " << str << endl;
      
      cout << endl << "Testing StringReadXML" << endl;
      xml_read_from_file ("string.xml", str);
      cout << "Read: " << str << endl;
    }
  catch (runtime_error e)
    {
      cerr << e.what ();
    }

  return (0);
}
