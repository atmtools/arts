/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
   \file   main.cc

   This file contains the main function of ARTS, as well as functions
   to deal with command line parameters. It also contains the
   executor, which is the `engine' that executes workspace methods in
   a controlfile one by one, in order to carry out an ARTS
   calculations. 

   \author Stefan Buehler
   \date   2001-07-24
*/

#include "arts.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <map>
#include "parameters.h"
#include "messages.h"
#include "exceptions.h"
#include "file.h"
#include "auto_wsv.h"
#include "methods.h"
#include "parser.h"
#include "auto_md.h"
#include "absorption.h"
#include "wsv_aux.h"

// This must be here rather than in arts.h, because arts.h does not
// know about Array.
void define_wsv_pointers(Array<WsvP*>&    wsv_pointers,
                         WorkSpace&       workspace);


/** Print the error message and exit. */
void give_up(const String& message)
{
  out0 << message << '\n';
  exit(1);
}


/** The arts executor. This executes the methods specified in tasklist
    on the given workspace. It also checks for errors during the
    method execution and stops the program if an error has
    occured. 

    \param workspace Output. The workspace to act on.
    \param tasklist The list of methods to execute (including keyword data).

    \author Stefan Buehler  */
void executor(WorkSpace& workspace, const Array<MRecord>& tasklist)
{
  // The method description lookup table:
  extern const Array<MdRecord> md_data;

  // The workspace variable lookup table:
  extern const Array<WsvRecord> wsv_data;
  
  // The array holding the pointers to the getaway functions:
  extern const void (*getaways[])(WorkSpace&, const MRecord&);

  // We need a place to remember which workspace variables are
  // occupied and which aren't.
  // 
  // For some weird reason, Arrays of bool do not work, although all
  // other types seem to work fine. So in this case, I'll use
  // the stl vector directly. The other place where this is done is in
  // the function lines_per_tgCreateFromLines in m_abs.cc.
  // FIXME: Fix this when Array<bool> works.
  std::vector<bool> occupied(wsv_data.nelem(),false);

  out3 << "\nExecuting methods:\n";

//   for (Index i=0; i<tasklist.nelem(); ++i)
//     {
//       const MRecord&  mrr = tasklist[i];
//       cout << "id, input: " << mrr.Id() << ", ";
//       print_vector(mrr.Input());
//       cout << "id, output: " << mrr.Id() << ", ";
//       print_vector(mrr.Output());
//     }


  for (Index i=0; i<tasklist.nelem(); ++i)
    {
      // Runtime method data for this method:
      const MRecord&  mrr = tasklist[i];
      // Method data for this method:
      const MdRecord& mdd = md_data[mrr.Id()];
      
      try
	{
	  out1 << "- " << mdd.Name() << '\n';

	
	  { // Check if all specific input variables are occupied:
	    const ArrayOfIndex& v(mdd.Input());
	    for (Index s=0; s<v.nelem(); ++s)
	      if (!occupied[v[s]])
		give_up("Method "+mdd.Name()+" needs input variable: "+
			wsv_data[v[s]].Name());
	  }

	  { // Check if all generic input variables are occupied:
	    const ArrayOfIndex& v(mrr.Input());
	    //	    cout << "v.nelem(): " << v.nelem() << endl;
	    for (Index s=0; s<v.nelem(); ++s)
	      if (!occupied[v[s]])
		give_up("Generic Method "+mdd.Name()+" needs input variable: "+
			wsv_data[v[s]].Name());
	  }

	  // Call the getaway function:
	  getaways[mrr.Id()]
	    ( workspace, mrr );

	  { // Flag the specific output variables as occupied:
	    const ArrayOfIndex& v(mdd.Output());
	    for (Index s=0; s<v.nelem(); ++s) occupied[v[s]] = true;
	  }

	  { // Flag the generic output variables as occupied:
	    const ArrayOfIndex& v(mrr.Output());
	    for (Index s=0; s<v.nelem(); ++s) occupied[v[s]] = true;
	  }

	}
      catch (runtime_error x)
	{
	  out0 << "Run-time error in method: " << mdd.Name() << '\n'
	       << x.what() << '\n';
	  exit(1);
	}
      catch (logic_error x)
	{
	  out0 << "Logic error in method: " << mdd.Name() << '\n'
	       << x.what() << '\n';
	  exit(1);
	}
    }
}


/** Remind the user of --help and exit return value 1. */
void polite_goodby()
{
  cerr << "Try `arts --help' for help.\n";
  exit(1);
}

/** Set the reporting level, either the default or based on
    reporting. If reporting was specified, check if the values make
    sense. The value -1 for reporting means that it was (probably)
    not given on the command line, since this is the initialization
    value.
    \param r Reporting level from Command line.
    \author Stefan Buehler */
void set_reporting_level(Index r)
{
  // The global variable that holds the message levels for screen and file.
  extern Messages messages;

  if ( -1 == r )
    {
      // Reporting was not specified, set default. (Only the
      // important stuff to the screen, everything to the file.)
      messages.screen = 1;
      messages.file   = 3;
    }
  else
    {
      // Reporting was specified. Check consistency and set report
      // level accordingly. 
	
	// Separate the two digits by taking modulo 10:
      Index s = r / 10;
      Index f = r % 10;
      //	cout << "s=" << s << " f=" << f << '\n';

      if ( s<0 || s>3 || f<0 || f>3 )
	{
	  cerr << "Illegal value specified for --reporting (-r).\n"
	       << "The specified value is " << r << ", which would be\n"
	       << "interpreted as screen=" << s << ", file=" << f << ".\n"
	       << "Only values of 0-3 are allowed for screen and file.\n";
	  exit(1);
	}
      messages.screen = s;
      messages.file   = f;
    }
}


/** React to option `methods'. If given the argument `all', it
    should simply prints a list of all methods. If given the name of
    a variable, it should print all methods that produce this
    variable as output.

    \param methods All or name of a variable.
    \author Stefan Buehler */
void option_methods(const String& methods)
{
  // Make global data visible:
  extern const Array<MdRecord>  md_data;
  extern const Array<WsvRecord> wsv_data;
  //  extern const std::map<String, Index> MdMap;
  extern const std::map<String, Index> WsvMap;
  extern const ArrayOfString wsv_group_names;

  // This is used to count the number of matches to a query, so
  // that `none' can be output if necessary
  Index hitcount;

  // First check if the user gave the special name `all':

  if ( "all" == methods )
    {
      cout
	<< "\n*--------------------------------------------------------------*\n"
	<< "Complete list of ARTS workspace methods:\n"
	<< "----------------------------------------------------------------\n";
      for ( Index i=0; i<md_data.nelem(); ++i )
	{
	  cout << "- " << md_data[i].Name() << '\n';
	}
      cout
	<< "*--------------------------------------------------------------*\n\n";
      return;
    }

  // Ok, so the user has probably specified a workspace variable or
  // workspace variable group.

  // Check if the user gave the name of a specific variable.
  map<String, Index>::const_iterator mi =
    WsvMap.find(methods);
  if ( mi != WsvMap.end() )
    {
      // If we are here, then the given name matches a variable.
      Index wsv_key = mi->second;

      // List generic methods:
      hitcount = 0;
      cout 
	<< "\n*--------------------------------------------------------------*\n"
	<< "Generic methods that can generate " << wsv_data[wsv_key].Name() 
	<< ":\n"
	<< "----------------------------------------------------------------\n";
      for ( Index i=0; i<md_data.nelem(); ++i )
	{
	  // This if statement checks whether GOutput, the list
	  // of output variable types contains the group of the
	  // requested variable.
	  if ( count( md_data[i].GOutput().begin(),
		      md_data[i].GOutput().end(),
		      wsv_data[wsv_key].Group() ) )
	    {
	      cout << "- " << md_data[i].Name() << '\n';
	      ++hitcount;
	    }
	}
      if ( 0==hitcount )
	cout << "none\n";

      // List specific methods:
      hitcount = 0;
      cout 
	<< "\n----------------------------------------------------------------\n"
	<< "Specific methods that can generate " << wsv_data[wsv_key].Name() 
	<< ":\n"
	<< "----------------------------------------------------------------\n";
      for ( Index i=0; i<md_data.nelem(); ++i )
	{
	  // This if statement checks whether Output, the list
	  // of output variables contains the workspace
	  // variable key.
	  if ( count( md_data[i].Output().begin(),
		      md_data[i].Output().end(),
		      wsv_key ) ) 
	    {
	      cout << "- " << md_data[i].Name() << '\n';
	      ++hitcount;
	    }
	}
      if ( 0==hitcount )
	cout << "none\n";

      cout
	<< "*--------------------------------------------------------------*\n\n";

      return;
    }

  // Check if the user gave the name of a variable group.

  // We use the find algorithm from the STL to do this. It
  // returns an iterator, so to get the index we take the
  // difference to the begin() iterator.
  Index group_key =
    find( wsv_group_names.begin(),
	  wsv_group_names.end(),
	  methods ) - wsv_group_names.begin();

  // group_key == wsv_goup_names.nelem() indicates that a
  // group with this name was not found.
  if ( group_key != wsv_group_names.nelem() )
    {
      // List generic methods:
      hitcount = 0;
      cout 
	<< "\n*--------------------------------------------------------------*\n"
	<< "Generic methods that can generate variables of group " 
	<< wsv_group_names[group_key] << ":\n"
	<< "----------------------------------------------------------------\n";
      for ( Index i=0; i<md_data.nelem(); ++i )
	{
	  // This if statement checks whether GOutput, the list
	  // of output variable types contains the
	  // requested group.
	  if ( count( md_data[i].GOutput().begin(),
		      md_data[i].GOutput().end(),
		      group_key ) )
	    {
	      cout << "- " << md_data[i].Name() << '\n';
	      ++hitcount;
	    }
	}
      if ( 0==hitcount )
	cout << "none\n";

      cout
	<< "*--------------------------------------------------------------*\n\n";

      return;
    }

  // If we are here it means that what the user specified is neither
  // `all', nor a variable, nor a variable group.
  cerr << "The name " << methods << " matches neither `all',\n"
       << "nor the name of a workspace variable, nor the name\n"
       << "of a workspace variable group.\n";
  exit(1);
}

/** React to option `input'. Given the name of
    a variable, it should print all methods that need this
    variable as input.

    \param methods Name of a variable.
    \author Stefan Buehler
    \date   2001-07-24 */
void option_input(const String& input)
{
  // Make global data visible:
  extern const Array<MdRecord>  md_data;
  extern const Array<WsvRecord> wsv_data;
  //  extern const std::map<String, Index> MdMap;
  extern const std::map<String, Index> WsvMap;
  extern const ArrayOfString wsv_group_names;

  // This is used to count the number of matches to a query, so
  // that `none' can be output if necessary
  Index hitcount;

  // Ok, so the user has probably specified a workspace variable or
  // workspace variable group.

  // Check if the user gave the name of a specific variable.
  map<String, Index>::const_iterator mi =
    WsvMap.find(input);
  if ( mi != WsvMap.end() )
    {
      // If we are here, then the given name matches a variable.
      Index wsv_key = mi->second;

      // List generic methods:
      hitcount = 0;
      cout 
      << "\n*--------------------------------------------------------------*\n"
      << "Generic methods that can use " << wsv_data[wsv_key].Name() << ":\n"
      << "----------------------------------------------------------------\n";
      for ( Index i=0; i<md_data.nelem(); ++i )
	{
	  // This if statement checks whether GInput, the list
	  // of input variable types contains the group of the
	  // requested variable.
	  if ( count( md_data[i].GInput().begin(),
		      md_data[i].GInput().end(),
		      wsv_data[wsv_key].Group() ) )
	    {
	      cout << "- " << md_data[i].Name() << '\n';
	      ++hitcount;
	    }
	}
      if ( 0==hitcount )
	cout << "none\n";

      // List specific methods:
      hitcount = 0;
      cout 
      << "\n----------------------------------------------------------------\n"
      << "Specific methods that require " << wsv_data[wsv_key].Name() 
      << ":\n"
      << "----------------------------------------------------------------\n";
      for ( Index i=0; i<md_data.nelem(); ++i )
	{
	  // This if statement checks whether Output, the list
	  // of output variables contains the workspace
	  // variable key.
	  if ( count( md_data[i].Input().begin(),
		      md_data[i].Input().end(),
		      wsv_key ) ) 
	    {
	      cout << "- " << md_data[i].Name() << '\n';
	      ++hitcount;
	    }
	}
      if ( 0==hitcount )
	cout << "none\n";

      cout
	<< "*--------------------------------------------------------------*\n\n";

      return;
    }

  // Check if the user gave the name of a variable group.

  // We use the find algorithm from the STL to do this. It
  // returns an iterator, so to get the index we take the
  // difference to the begin() iterator.
  Index group_key =
    find( wsv_group_names.begin(),
	  wsv_group_names.end(),
	  input ) - wsv_group_names.begin();

  // group_key == wsv_goup_names.nelem() indicates that a
  // group with this name was not found.
  if ( group_key != wsv_group_names.nelem() )
    {
      // List generic methods:
      hitcount = 0;
      cout
      << "\n*--------------------------------------------------------------*\n"
      << "Generic methods that require a variable of group " 
      << wsv_group_names[group_key] << ":\n"
      << "----------------------------------------------------------------\n";
      for ( Index i=0; i<md_data.nelem(); ++i )
	{
	  // This if statement checks whether GOutput, the list
	  // of output variable types contains the
	  // requested group.
	  if ( count( md_data[i].GInput().begin(),
		      md_data[i].GInput().end(),
		      group_key ) )
	    {
	      cout << "- " << md_data[i].Name() << '\n';
	      ++hitcount;
	    }
	}
      if ( 0==hitcount )
	cout << "none\n";

      cout
	<< "*--------------------------------------------------------------*\n\n";

      return;
    }

  // If we are here it means that what the user specified is neither
  // a variable nor a variable group.
  cerr << "The name " << input << " matches neither the name of a\n"
       << "workspace variable, nor the name of a workspace variable group.\n";
  exit(1);
}


/** React to option `workspacevariables'. If given the argument `all',
    it should simply prints a list of all variables. If given the
    name of a method, it should print all variables that are needed
    by that method.

    \param  workspacevariables All or name of a method.
    \author Stefan Buehler */
void option_workspacevariables(const String& workspacevariables)
{
  // Make global data visible:
  extern const Array<MdRecord>  md_data;
  extern const Array<WsvRecord> wsv_data;
  extern const std::map<String, Index> MdMap;
  //  extern const std::map<String, Index> WsvMap;
  extern const ArrayOfString wsv_group_names;

  // This is used to count the number of matches to a query, so
  // that `none' can be output if necessary
  Index hitcount;

  // First check for `all':

  if ( "all" == workspacevariables )
    {
      cout
	<< "\n*--------------------------------------------------------------*\n"
	<< "Complete list of ARTS workspace variables:\n"
	<< "----------------------------------------------------------------\n";
      for ( Index i=0; i<wsv_data.nelem(); ++i )
	{
	  cout << "- " << wsv_data[i].Name() << '\n';
	}
      cout
	<< "*--------------------------------------------------------------*\n\n";
      return;
    }


  // Now check if the user gave the name of a method.
  map<String, Index>::const_iterator mi =
    MdMap.find(workspacevariables);
  if ( mi != MdMap.end() )
    {
      // If we are here, then the given name matches a method.
      // Assign the data record for this method to a local
      // variable for easier access:
      const MdRecord& mdr = md_data[mi->second];
      
      // List generic variables required by this method.
      hitcount = 0;
      cout
      << "\n*--------------------------------------------------------------*\n"
      << "Generic workspace variables required by " << mdr.Name()
      << " are of type:\n"
      << "----------------------------------------------------------------\n";
      for ( Index i=0; i<mdr.GInput().nelem(); ++i )
	{
	  cout << "- " << wsv_group_names[mdr.GInput()[i]] << "\n";
	  ++hitcount;
	}
      if ( 0==hitcount )
	cout << "none\n";

      // List specific variables required by this method.
      hitcount = 0;
      cout 
      << "\n----------------------------------------------------------------\n"
      << "Specific workspace variables required by " << mdr.Name() << ":\n"
      << "----------------------------------------------------------------\n";
      for ( Index i=0; i<mdr.Input().nelem(); ++i )
	{
	  cout << "- " << wsv_data[mdr.Input()[i]].Name() << '\n';
	  ++hitcount;
	}
      if ( 0==hitcount )
	cout << "none\n";

      cout
	<< "*--------------------------------------------------------------*\n\n";

      return;
    }

  // If we are here, then the user specified nothing that makes sense.
  cerr << "The name " << workspacevariables << " matches neither `all',\n" 
       << "nor the name of a workspace method.\n";
  exit(1);
}


/** React to option `describe'. This should print the description
    String of the given workspace variable or method.

    \param describe What to describe.
    \author Stefan Buehler */
void option_describe(const String& describe)
{
  // Make global data visible:
  extern const Array<MdRecord>  md_data;
  extern const Array<WsvRecord> wsv_data;
  extern const std::map<String, Index> MdMap;
  extern const std::map<String, Index> WsvMap;
  //  extern const ArrayOfString wsv_group_names;

  // Let's first assume it is a method that the user wants to have
  // described.

  // Find method id:
  map<String, Index>::const_iterator i =
    MdMap.find(describe);
  if ( i != MdMap.end() )
    {
      // If we are here, then the given name matches a method.
      cout << md_data[i->second] << '\n';
      return;
    }

  // Ok, let's now assume it is a variable that the user wants to have
  // described.

  // Find wsv id:
  i = WsvMap.find(describe);
  if ( i != WsvMap.end() )
    {
      // If we are here, then the given name matches a workspace
      // variable. 
      cout << wsv_data[i->second] << '\n';
      return;
    }

  // If we are here, then the given name does not match anything.
  cerr << "The name " << describe
       << " matches neither method nor variable.\n";
  exit(1);      
}


//! Checks the dimensions of stuff in generated headers.
/*!
  The header files wsv.h and md.h are generated automatically from the
  files workspace.cc and methods.cc. This function checks, if at least
  the dimensions agree between the .h files and the .cc files.

  FIXME: Update this, add automatic wsv groups.

  \author Stefan Buehler
  \date   2000-08-04
*/
void check_built_headers()
{
  // Make lookup data visible:
  //  extern const Array<MdRecord>  md_data;
  extern const ArrayOfString wsv_group_names;
  extern const Array<WsvRecord> wsv_data;

  // Checks:
  assert( N_WSV_GROUPS == wsv_group_names.nelem() );
  assert( N_WSV        == wsv_data.nelem()        );

}

/** This is the main function of ARTS. (You never guessed that, did you?)
    The getopt_long function is used to parse the command line parameters.
 
    \verbatim
    Overview:
    1. Get command line parameters.
    2. Evaluate the command line parameters. (This also checks if the
       parameters make sense, where necessary.) 
    \endverbatim

    \return    0=ok, 1=error
    \param     argc Number of command line parameters 
    \param     argv Values of command line parameters
    \author    Stefan Buehler */
int main (int argc, char **argv)
{
  extern const Parameters parameters; // Global variable that holds
                                      // all command line parameters. 

  //---------------< 1. Get command line parameters >---------------
  if ( get_parameters(argc, argv) )
    {
      // Print an error message and exit:
      polite_goodby();
    }

  //----------< 2. Evaluate the command line parameters >----------
  if (parameters.help)
    {
      // Just print a help message and then exit.
      cerr << '\n' << parameters.usage << "\n\n";
      cerr << parameters.helptext << "\n\n";
      return(0);
    }

  if (parameters.version)
    {
      extern const String full_name;
      // Just print version information and then exit.
      cerr << "This is " << full_name
#ifdef HDF_SUPPORT
           << " with HDF support."
#endif // HDF_SUPPORT
           << '\n'
           << "Compiled on " << OS_NAME << " " << OS_VERSION
           << " with "
           << ((sizeof (Numeric) == sizeof (double)) ? "double" : "float")
           << " precision." << endl;
      return(0);
    }



  // For the next couple of options we need to have the workspce and
  // method lookup data.


  // Initialize the md data:
  define_md_data();

  // Initialize the wsv group name array:
  define_wsv_group_names();

  // Initialize the wsv data:
  define_wsv_data();

  // Initialize wsv pointers:
  {
    // We need to make the workspace and the wsv_pointers visible for
    // a moment. This is because the pointers will be inititalized for
    // this particular workspace. There can be other (local)
    // workspaces, which have their own pointers associated with
    // them. 
    extern WorkSpace workspace;
    extern Array<WsvP*> wsv_pointers;
    define_wsv_pointers(wsv_pointers,workspace);
  }

  // Initialize MdMap
  define_md_map();

  // Initialize WsvMap
  define_wsv_map();


  // While we are at it, we can also initialize the molecular data and
  // the coefficients of the partition function that we need for the
  // absorption part, and check that the inputs are sorted the same way:
  define_species_data();

  // And also the species map:
  define_species_map();

  // And the lineshape lookup data:
  define_lineshape_data();
  define_lineshape_norm_data();

  // Make all these data visible:
  //  extern const Array<MdRecord>  md_data;
  //  extern const Array<WsvRecord> wsv_data;
  //  extern const std::map<String, Index> MdMap;
  //  extern const std::map<String, Index> WsvMap;
  extern const ArrayOfString wsv_group_names;

  // Now we are set to deal with the more interesting command line
  // switches. 


  // React to option `methods'. If given the argument `all', it
  // should simply prints a list of all methods. If given the name of
  // a variable, it should print all methods that produce this
  // variable as output.
  if ( "" != parameters.methods )
    {
      option_methods(parameters.methods);
      return(0);
    }

  // React to option `input'. Given the name of a variable (or group)
  // it should print all methods that need this variable (or group) as
  // input.
  if ( "" != parameters.input )
    {
      option_input(parameters.input);
      return(0);
    }
  
  // React to option `workspacevariables'. If given the argument `all',
  // it should simply prints a list of all variables. If given the
  // name of a method, it should print all variables that are needed
  // by that method.
  if ( "" != parameters.workspacevariables )
    {
      option_workspacevariables(parameters.workspacevariables);
      return(0);
    }

  // React to option `describe'. This should print the description
  // String of the given workspace variable or method.
  if ( "" != parameters.describe )
    {
      option_describe(parameters.describe);
      return(0);
    }

  
  // React to option `groups'. This should simply print a list of all
  // workspace variable groups.
  if ( parameters.groups )
    {
      cout
	<< "\n*--------------------------------------------------------------*\n"
	<< "Complete list of ARTS workspace variable groups:\n"
	<< "----------------------------------------------------------------\n";
      for ( Index i=0; i<wsv_group_names.nelem(); ++i )
	{
	  cout << "- " << wsv_group_names[i] << '\n';
	}
      cout
	<< "*--------------------------------------------------------------*\n\n";
      return(0);
    }


  // Ok, we are past all the special options. This means the user
  // wants to get serious and really do a calculation. Check if we
  // have at least one control file:
  if ( 0 == parameters.controlfiles.nelem() )
    {
      cerr << "You must specify at least one control file name.\n";
      polite_goodby();
    }

  // Set the basename according to the first control file, if not
  // explicitly specified.
  if ( "" == parameters.basename )
    {
      extern String out_basename;
      out_basename = parameters.controlfiles[0];
      // Find the last . in the name
      String::size_type p = out_basename.rfind('.');
      // Kill everything starting from the `.'
      out_basename.erase(p);
    }

  // Set the reporting level, either from reporting command line
  // option or default.  
  set_reporting_level(parameters.reporting);


  //--------------------< Open report file >--------------------
  // This one needs its own little try block, because we have to
  // write error messages to cerr directly since the report file
  // will not exist.
  try
    {
      extern const String out_basename;     // Basis for file name
      extern ofstream report_file;	// Report file pointer

      //      cout << "rep = " << out_basename+".rep" << '\n';
      open_output_file(report_file, out_basename+".rep");
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
      //      extern Array<WsvRecord> wsv_data;
      //      extern Array<MdRecord> md_data;

      {
	// Output program name and version number: 
	// The name (PACKAGE) and the major and minor version number
	// (VERSION) are set in configure.in. The configuration tools
	// place them in the file config.h, which is included in arts.h.
  
	extern const String full_name;

	out1 << full_name << '\n';
      }


      // The list of methods to execute and their keyword data from
      // the control file. 
      Array<MRecord> tasklist;

      // The text of the controlfile.
      SourceText text;
	
      // Read the control text from the control files:
      out3 << "\nReading control files:\n";
      for ( Index i=0; i<parameters.controlfiles.nelem(); ++i )
	{
	  out3 << "- " << parameters.controlfiles[i] << '\n';
	  text.AppendFile(parameters.controlfiles[i]);
	}

      // Call the parser to parse the control text:
      parse_main(tasklist, text);

      // Call the executor:
      executor(workspace, tasklist);

    }
  catch (runtime_error x)
    {
      out0 << x.what() << '\n';
      exit(1);
    }

  out1 << "Goodbye.\n";
  return(0);
}
