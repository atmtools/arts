#include "arts.h"
#include "messages.h"
#include "file.h"

void open_output_file(ofstream& file, const string& name)
{
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // FIXME: This does not yet work in  egcs-2.91.66, try again later.
  file.exceptions(ios::badbit |
		  ios::failbit);
  
  // c_str explicitly converts to c string.
  file.open(name.c_str() );

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly. (In that case we would not
  // get here if there really was a problem, because of the exception
  // thrown by open().)
  if (!file)
    {
      std::ostrstream os;
      os << "Cannot open output file: " << name << '\n'
	 << "Maybe you don't have write access "
	 << "to the directory or the file?";
      throw runtime_error(os.str());
    }
}


void open_input_file(ifstream& file, const string& name)
{
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  file.exceptions(ios::badbit |
		  ios::failbit);

  // c_str explicitly converts to c string.
  file.open(name.c_str() );

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly.
  if (!file)
    {
      std::ostrstream os;
      os << "Cannot open input file: " << name << '\n'
	 << "Maybe the file does not exist?";
      throw runtime_error(os.str());
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
  // FIXME: This should not be necessary anymore when stream
  // exceptions work properly.
  if ( !is.eof() ) {
    std::ostrstream os;
    os << "Read Error. Last line read:\n" << linebuffer;
    throw runtime_error(os.str());
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
      throw runtime_error(os.str());
    }
}

void replace_all(string& s, const string& what, const string& with)
{
  string::size_type j = s.find(what);
  while ( j != string::npos )
    {
      //		cout << "j = " << j << '\n';
      s.replace(j,1,with);
      j = s.find(what,j+with.size());
    }
}

