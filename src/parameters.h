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
    controlfiles(),
    reporting(-1),
    methods(""),
    workspacevariables(""),
    describe(""),
    groups(false)
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
  /** This should be a two digit integer. The first digit specifies
      the output level for stdout (stderr for error messages), the
      second digit the output level for the report file. The levels
      can reach from 0 (show only error messages) to 3 (show
      everything). Example:

      03 = only errors to the screen, everything to the file. */
  int reporting;
  /** If this is given the argument `all', it simply prints a list of 
      all methods. If it is given the name of a variable, it
      prints all methods that produce this variable as output. */
  string methods;
  /** If this is given the argument `all', it simply prints a list of 
      all workspace variables. If it is given the name of a method,
      it prints all variables needed by that method. */
  string workspacevariables; 
  /** Print the description string of the given workspace variable or
      method. */
  string describe;
  /** Print a list of all workspace variable groups. */
  bool groups;
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
