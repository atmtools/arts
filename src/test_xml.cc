using namespace std;

#include <iostream>
#include <string>
#include "xml_io.h"

int
main (int argc, char *argv[])
{
  String str;

  cout << "Testing StringWriteXML" << endl;
  str = "Hello World";
  xml_write_to_file ("string.xml", str);
  cout << "Wrote: " << str << endl;

  cout << endl << "Testing StringReadXML" << endl;
  xml_read_from_file ("string.xml", str);
  cout << "Read: " << str << endl;

  return (0);
}
