#ifndef parameters_h
#define parameters_h

#include "vecmat.h"
#include "string.h"

/**
   Structure to hold all command line Parameters. This holds all the
   command line parameters, plut the usage message and the helptext
   message. The messages are in the same structure, because they need
   to be changed whenever the parameters are changed, so it is better
   to have them in the same place. 
   @author SAB
*/
struct Parameters {
  /** Default constructor. Care has to be taken to properly initialize
      all variables, e.g., bool options to false. */
  Parameters() :
    help(false),
    version(false),
    basename(""),
    controlfiles()
  { /* Nothing to be done here */ }
  /** Short message how to call the program. */
  string usage;
  /** Longer message explaining the options. */
  string helptext;
  /** Only display the help text. */
  bool help;			
  /** Display version information. */
  bool version;			
  /** If this is specified (with the -b --basename option), it is used
      as the base name for the report file and for other output
      files. */ 
  string basename;
  /** The filenames of the controlfiles. Can be only one or as many as
      you want. */
  ARRAY<string> controlfiles;
};


/**
   Get the command line parameters. They are stored in the global
   variable parameters which is a structure of type Parameters. If
   needed, this variable should be declared like this:

   extern const Parameters parameters

   @return    false=ok, true=error
   @param     argc Number of command line parameters 
   @param     argv Values of command line parameters
   @author    Stefan Buehler
   @version   1
 */
bool get_parameters(int argc, char **argv);


#endif // parameters_h
