#include "arts.h"
#include "parameters.h"
#include "messages.h"
#include "exceptions.h"
#include "file.h"
#include "wsv.h"		// Automatically generated!
#include "workspace.h"
#include "methods.h"
#include "parser.h"
#include "md.h"

/*
   This is the main function of ARTS. (You never guessed that, did you?)
   The getopt_long function is used to parse the command line parameters.

   @return    0=ok, 1=error
   @param     argc Number of command line parameters 
   @param     argv Values of command line parameters
   @author    Stefan Buehler
   @version   1
*/
int main (int argc, char **argv)
{
  extern const Parameters parameters; // Global variable that holds
                                      // all command line parameters. 

  //--------------------< Get command line parameters >--------------------
  if ( get_parameters(argc, argv) )
    {
      // Print an error message and exit:
      cerr << "Try `arts --help' for help.\n";
      exit(1);
    }

  if (parameters.help)
    {
      // Just print a help message and then exit.
      cerr << '\n' << parameters.usage << "\n\n";
      cerr << parameters.helptext << "\n\n";
      return(0);
    }

  if (parameters.version)
    {
      extern const string subversion;
      // Just print version information and then exit.
      cerr << "This is " 
	   << PACKAGE << " " << VERSION << "." << subversion << '\n';
      return(0);
    }

  //--------------------< Open report file >--------------------
  // This one needs its own little try block, because we have to
  // write error messages to cerr directly since the report file
  // will not exist.
  try
    {
      extern ofstream report_file;	// Report file pointer

      //      cout << "rep = " << parameters.reportFile << '\n';
      open_output_file(report_file, parameters.reportFile);
    }
  catch (runtime_error x)
    {
      cerr << x.what() << '\n'
	   << "I have to be able to write to my report file.";
      exit(1);
    }

  // Now comes the global try block. Exceptions caught after this
  // one are general stuff like file opening errors.
  // FIXME: Maybe this is not really necessary, since methods
  // using files could always check for these errors? Think about
  // which way is easier.
  try
    {
      // Some global variables that we need:
      extern WorkSpace workspace;
      //      extern ARRAY<WsvRecord> wsv_data;
      extern ARRAY<MdRecord> md_data;

      {
	// Output program name and version number: 
	// The name (PACKAGE) and the major and minor version number
	// (VERSION) are set in configure.in. The configuration tools
	// place them in the file config.h, which is included in arts.h.
	// The subminor number is set in version.cc, which is linked with
	// arts.

	extern const string subversion;
  
	out1 << PACKAGE << " " << VERSION << "." << subversion << '\n';
      }

      // Initialize the wsv data:
      define_wsv_data();

      {
	// Quick test that the pointers stored in wsv_data can be
	// really used to access the workspace variables.
	// YES!! It works.
	// FIXME: Remove all this.
// 	WsvP   *wp = wsv_data[basename_].Pointer(); 
// 	string *s = *wp;
// 	*s = "test";
// 	cout << "workspace.basename = " << workspace.basename << '\n';
      }

      // Initialize the md data:
      define_md_data();

      {
	// Quick test:
// 	for (size_t i=0; i<md_data.size(); ++i)
// 	  {
// 	    cout << md_data[i].Name() << '\n';
// 	  }
      }

      // The list of methods to execute and their keyword data from
      // the control file. 
      ARRAY<MRecord> tasklist;

      // The text of the controlfile.
      SourceText text;
	
      // Read the control text from the control files:
      out3 << "\nReading control files:\n";
      for ( size_t i=0; i<parameters.controlfiles.size(); ++i )
	{
	  out3 << "- " << parameters.controlfiles[i] << '\n';
	  text.AppendFile(parameters.controlfiles[i]);
	}

      // Call the parser to parse the control text:
      parse_main(tasklist, text);

      // Execute the methods in tasklist.
      out3 << "\nExecuting methods:\n";
      for (size_t i=0; i<tasklist.size(); ++i)
	try
	  {
	    // The array holding the pointers to the getaway functions:
	    extern void (*getaways[])(WorkSpace&, const ARRAY<TokVal>&);

	    out1 << "- " << md_data[tasklist[i].Id()].Name() << '\n';
	    
	    getaways[tasklist[i].Id()]
	      ( workspace, tasklist[i].Values() );
	  }
	catch (runtime_error x)
	  {
	    out0 << "Error in method: " << md_data[tasklist[i].Id()].Name() << '\n'
		 << x.what() << '\n';
	    exit(1);
	  }
    }
  catch (runtime_error x)
    {
      out0 << x.what() << '\n';
      exit(1);
    }

  out1 << "Goodbye.\n";
  return(0);
}
