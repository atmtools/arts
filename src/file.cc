#include "arts.h"
#include "exceptions.h"
#include "file.h"

void open_output_file(ofstream& file, const string& name)
{
  // c_str explicitly converts to c string.
  file.open(name.c_str() );

  // See if the file is ok.
  if (!file)
    {
      std::ostrstream os;
      os << "Cannot open output file: " << name << '\n'
	 << "Maybe you don't have write access "
	 << "to the directory or the file?";
      throw IOError(os.str());
    }
}


void open_input_file(ifstream& file, const string& name)
{
  // c_str explicitly converts to c string.
  file.open(name.c_str() );

  // See if the file is ok.
  if (!file)
    {
      std::ostrstream os;
      os << "Cannot open input file: " << name;
      throw IOError(os.str());
    }
}


void read_text_from_stream(ARRAY<string>& text, istream& is)
{
  string linebuffer;

  // Read as long as `is' is good:
  while (is)
    {
      // Read line from file into linebuffer:
      getline(is,linebuffer);

      // Append to end of text:
      text.push_back(linebuffer);
    }
  
  // Check for error:
  if ( !is.eof() ) {
    std::ostrstream os;
    os << "Read Error. Last line read:\n" << linebuffer;
    throw IOError(os.str());
  }

}


void read_text_from_file(ARRAY<string>& text, const string& name)
{
  ifstream ifs;

  // Open input stream:
  open_input_file(ifs, name);
  // No need to check for error, because open_input_file throws an
  // IOError with an appropriate error message.

  // Read the text from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the 
  // filename.
  try
    {
      read_text_from_stream(text,ifs);
    }
  catch (runtime_error x)
    {
      std::ostrstream os;
      os << "Error reading file: " << name << '\n'
	 << x.what();
      throw IOError(os.str());
    }
}
