// #include "arts.h"
#include "arts.h"
#include <getopt.h>
#include "parameters.h"

/// Holds the command line parameters.
Parameters parameters;

bool get_parameters(int argc, char **argv)
{
  /*
    The header file getopt.h declares some external variables:

    extern char *optarg; (argument value for options that take an argument)
    extern int optind;   (index in ARGV of the next element to be scanned.)
    extern int opterr;   (we don´t use this)
    extern int optopt;   (set to an option caracter that was recoginized)
  */

  /*
    From the GNU documentation:

    Describe the long-named options requested by the application.
    The LONG_OPTIONS argument to getopt_long or getopt_long_only is a vector
    of `struct option' terminated by an element containing a name which is
    zero.

    The field `has_arg' is:
    no_argument          (or 0) if the option does not take an argument,
    required_argument    (or 1) if the option requires an argument,
    optional_argument    (or 2) if the option takes an optional argument.

    If the field `flag' is not NULL, it points to a variable that is set
    to the value given in the field `val' when the option is found, but
    left unchanged if the option is not found.

    To have a long-named option do something other than set an `int' to
    a compiled-in constant, such as set a value from `optarg', set the
    option's `flag' field to zero and its `val' field to a nonzero
    value (the equivalent single-letter option character, if there is
    one).  For long options that have a zero `flag' field, `getopt'
    returns the contents of the `val' field. 

    struct option
    {
       #if defined (__STDC__) && __STDC__
          const char *name;
       #else
          char *name;
       #endif
       int has_arg;
       int *flag;
       int val;
    };              
  */
  struct option longopts[] =
  {
    { "help",        no_argument, 	NULL, 'h' },
    { "version",     no_argument, 	NULL, 'v' },
    { "basename",    required_argument, NULL, 'b' },
    { NULL,          no_argument,       NULL, 0   }
  };

  parameters.usage =
    "Usage: arts [-hvb] [--help] [--version] [--basename <name>] file1.arts file2.arts ...";

  parameters.helptext =
    "The Atmospheric Radiative Transfer System.\n\n"
    "-h, --help          Print this message.\n"
    "-v, --version       Show version information.\n"
    "-b, --basename      Set the basename for the report\n"
    "                    file and for other output files.\n";

  // Set the short options automatically from the last columns of
  // longopts.
  //
  // Watch out! We also have to put in colons after options that take
  // an argument. Please do not use optional arguments (would be two
  // colons) since that works only with GNU getopt.
  //
  string shortopts;
  {
    int i=0;
    while (NULL != longopts[i].name )
      {
	// FIXME: Should there be a cast here? The real type of val is int.
	char c = longopts[i].val;
	shortopts += c;
	//	cout << "name = " << longopts[i].name << "\n";
	//	cout << "val  = " << longopts[i].val << "\n";

	// Check if we need to insert a colon
	if ( required_argument == longopts[i].has_arg )
	  {
	    shortopts += ':';
	  }
	++i;
      }
    shortopts += '\0';
  }
  //  cout << "shortopts: " << shortopts << '\n';

  int optc;

  while ( EOF != (optc = getopt_long (argc, argv, shortopts.c_str(),
				      longopts, (int *) 0) ) )
    {
      //      cout << "optc = " << optc << '\n';
      switch (optc)
	{
	case 'h':
	  parameters.help = true;
	  break;
	case 'v':
	  parameters.version = true;
	  break;
	case 'b':
	  //	  cout << "optarg = " << optarg << endl;
	  parameters.basename = optarg;
	  break;
	default:
	  // There were strange options.
	  return(1);
	  break;
	}
    }

  // Remaining things on the command line must be control file names. 
  // Get them one by one.
  if ( (optind==argc) &&
       !parameters.help &&
       !parameters.version)
    {
      cerr << "You must specify at least one control file name.\n";
      return(1);
    }
  else
    {
      while ( optind < argc )
	{
	  string dummy=argv[optind];
	  parameters.controlfiles.push_back(dummy);
	  optind++;
	}
    }
  //  cout << "alle:\n" << parameters.controlfiles << '\n';



  return false;
}
