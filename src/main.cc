/* Copyright (C) 2000-2008 Stefan Buehler <sbuehler@ltu.se>

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

#include <algorithm>
#include <map>
#include "arts.h"
#include "parameters.h"
#include "messages.h"
#include "exceptions.h"
#include "file.h"
#include "methods.h"
#include "parser.h"
#include "auto_md.h"
#include "absorption.h"
#include "wsv_aux.h"
#include "agenda_record.h"
#include "mystring.h"
#include "workspace_ng.h"
#include "arts_omp.h"


/** Remind the user of --help and exit return value 1. */
void polite_goodby()
{
  cerr << "Try `arts --help' for help.\n";
  arts_exit ();
}

/**
   Set the reporting level.

   Set the reporting level, either the default or based on
   reporting. If reporting was specified, check if the values make
   sense. The value -1 for reporting means that it was (probably)
   not given on the command line, since this is the initialization
   value.   

   \param r Reporting level from Command line.
   \author Stefan Buehler 
*/
void set_reporting_level(Index r)
{
  extern Messages arts_messages;

  if ( -1 == r )
    {
      // Reporting was not specified, set default. (No output from
      // agendas (except main of course), only the important stuff to
      // the screen, nothing to the file.)
      arts_messages.va = 0;     // agendas
      arts_messages.vs = 1;     // screen
      arts_messages.vf = 0;     // file
    }
  else
    {
      // Reporting was specified. Check consistency and set report
      // level accordingly. 
        
      // Separate the three digits:
      arts_messages.va = r/100;
      arts_messages.vs = (r%100)/10;
      arts_messages.vf = r%10;

      if ( !arts_messages.valid() )
        {
          cerr << "Illegal value specified for --reporting (-r).\n"
               << "The specified value is " << r << ", which would be\n"
               << "interpreted as:\n"
               << "Verbosity for agendas:     " << arts_messages.va << "\n"
               << "Verbosity for screen:      " << arts_messages.vs << "\n"
               << "Verbosity for report file: " << arts_messages.vf << "\n"
               << "Only values of 0-3 are allowed for each verbosity.\n";
          arts_exit ();
        }
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
  Workspace workspace;
  workspace.initialize();
  // Make global data visible:
  extern const Array<MdRecord>  md_data_raw;
  extern const Parameters parameters;
  //  extern const map<String, Index> MdMap;
  extern const ArrayOfString wsv_group_names;

  // This is used to count the number of matches to a query, so
  // that `none' can be output if necessary
  Index hitcount;

  // First check if the user gave the special name `all':

  if ( "all" == methods )
    {
      if (!parameters.plain)
        {
          cout
            << "\n*-------------------------------------------------------------------*\n"
            << "Complete list of ARTS workspace methods:\n"
            << "---------------------------------------------------------------------\n";
        }
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
        {
          if (!parameters.plain) cout << "- ";
          cout << md_data_raw[i].Name() << "\n";
        }

      if (!parameters.plain)
        cout << "*-------------------------------------------------------------------*\n\n";

      return;
    }

  // Ok, so the user has probably specified a workspace variable or
  // workspace variable group.

  // Check if the user gave the name of a specific variable.
  map<String, Index>::const_iterator mi =
    Workspace::WsvMap.find(methods);
  if ( mi != Workspace::WsvMap.end() )
    {
      // If we are here, then the given name matches a variable.
      Index wsv_key = mi->second;

      // List generic methods:
      hitcount = 0;
      cout 
        << "\n*-------------------------------------------------------------------*\n"
        << "Generic and supergeneric methods that can generate " << Workspace::wsv_data[wsv_key].Name() 
        << ":\n"
        << "---------------------------------------------------------------------\n";
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
        {
          // Get handle on method record:
          const MdRecord& mdd = md_data_raw[i];

          // This if statement checks whether GOutType, the list
          // of output variable types contains the group of the
          // requested variable.
          // The else clause picks up methods with supergeneric input.
          if ( count( mdd.GOutType().begin(),
                      mdd.GOutType().end(),
                      Workspace::wsv_data[wsv_key].Group() ) )
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
          else if  ( count( mdd.GOutType().begin(),
                      mdd.GOutType().end(),
                      get_wsv_group_id("Any") ) )
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
        }
      if ( 0==hitcount )
        cout << "none\n";

      // List specific methods:
      hitcount = 0;
      cout 
        << "\n---------------------------------------------------------------------\n"
        << "Specific methods that can generate " << Workspace::wsv_data[wsv_key].Name() 
        << ":\n"
        << "---------------------------------------------------------------------\n";
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
        {
          // Get handle on method record:
          const MdRecord& mdd = md_data_raw[i];

          // This if statement checks whether Output, the list
          // of output variables contains the workspace
          // variable key.
          if ( count( mdd.Output().begin(),
                      mdd.Output().end(),
                      wsv_key ) ) 
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
        }
      if ( 0==hitcount )
        cout << "none\n";

      cout
        << "*-------------------------------------------------------------------*\n\n";

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
        << "\n*-------------------------------------------------------------------*\n"
        << "Generic and supergeneric methods that can generate variables of group " 
        << wsv_group_names[group_key] << ":\n"
        << "---------------------------------------------------------------------\n";
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
        {
          // Get handle on method record:
          const MdRecord& mdd = md_data_raw[i];

          // This if statement checks whether GOutType, the list
          // of output variable types contains the
          // requested group.
          // The else clause picks up methods with supergeneric input.
          if ( count( mdd.GOutType().begin(),
                      mdd.GOutType().end(),
                      group_key ) )
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
          else if  ( count( mdd.GOutType().begin(),
                      mdd.GOutType().end(),
                      get_wsv_group_id("Any") ) )
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
        }
      if ( 0==hitcount )
        cout << "none\n";

      cout
        << "*-------------------------------------------------------------------*\n\n";

      return;
    }

  // If we are here it means that what the user specified is neither
  // `all', nor a variable, nor a variable group.
  cerr << "The name " << methods << " matches neither `all',\n"
       << "nor the name of a workspace variable, nor the name\n"
       << "of a workspace variable group.\n";
  arts_exit ();
}

/** React to option `input'. Given the name of
    a variable, it should print all methods that need this
    variable as input.

    \param input Name of a variable.
    \author Stefan Buehler
    \date   2001-07-24 */
void option_input(const String& input)
{
  // Make global data visible:
  extern const Array<MdRecord>  md_data_raw;
  //  extern const map<String, Index> MdMap;
  extern const ArrayOfString wsv_group_names;

  // This is used to count the number of matches to a query, so
  // that `none' can be output if necessary
  Index hitcount;

  // Ok, so the user has probably specified a workspace variable or
  // workspace variable group.

  // Check if the user gave the name of a specific variable.
  map<String, Index>::const_iterator mi =
    Workspace::WsvMap.find(input);
  if ( mi != Workspace::WsvMap.end() )
    {
      // If we are here, then the given name matches a variable.
      Index wsv_key = mi->second;

      // List generic methods:
      hitcount = 0;
      cout 
      << "\n*-------------------------------------------------------------------*\n"
      << "Generic and supergeneric methods that can use " << Workspace::wsv_data[wsv_key].Name() << ":\n"
      << "---------------------------------------------------------------------\n";
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
        {
          // Get handle on method record:
          const MdRecord& mdd = md_data_raw[i];
          
          // This if statement checks whether GInType, the list
          // of input variable types contains the group of the
          // requested variable.
          // The else clause picks up methods with supergeneric input.
          if ( count( mdd.GInType().begin(),
                      mdd.GInType().end(),
                      Workspace::wsv_data[wsv_key].Group() ) )
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
          else if  ( count( mdd.GInType().begin(),
                      mdd.GInType().end(),
                      get_wsv_group_id("Any") ) )
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
        }
      if ( 0==hitcount )
        cout << "none\n";

      // List specific methods:
      hitcount = 0;
      cout 
      << "\n---------------------------------------------------------------------\n"
      << "Specific methods that require " << Workspace::wsv_data[wsv_key].Name() 
      << ":\n"
      << "---------------------------------------------------------------------\n";
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
        {
          // Get handle on method record:
          const MdRecord& mdd = md_data_raw[i];

          // This if statement checks whether Output, the list
          // of output variables contains the workspace
          // variable key.
          if ( count( mdd.Input().begin(),
                      mdd.Input().end(),
                      wsv_key ) ) 
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
        }
      if ( 0==hitcount )
        cout << "none\n";

      cout
        << "*-------------------------------------------------------------------*\n\n";

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
      << "\n*-------------------------------------------------------------------*\n"
      << "Generic and supergeneric methods that require a variable of group " 
      << wsv_group_names[group_key] << ":\n"
      << "---------------------------------------------------------------------\n";
      for ( Index i=0; i<md_data_raw.nelem(); ++i )
        {
          // Get handle on method record:
          const MdRecord& mdd = md_data_raw[i];

          // This if statement checks whether GOutType, the list
          // of output variable types contains the
          // requested group.
          // The else clause picks up methods with supergeneric input.
          if ( count( mdd.GInType().begin(),
                      mdd.GInType().end(),
                      group_key ) )
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
          else if  ( count( mdd.GInType().begin(),
                      mdd.GInType().end(),
                      get_wsv_group_id("Any") ) )
            {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }   }
      if ( 0==hitcount )
        cout << "none\n";

      cout
        << "*-------------------------------------------------------------------*\n\n";

      return;
    }

  // If we are here it means that what the user specified is neither
  // a variable nor a variable group.
  cerr << "The name " << input << " matches neither the name of a\n"
       << "workspace variable, nor the name of a workspace variable group.\n";
  arts_exit ();
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
  extern const map<String, Index> MdMap;
  extern const Parameters parameters;
  //  extern const map<String, Index> WsvMap;
  extern const ArrayOfString wsv_group_names;

  // This is used to count the number of matches to a query, so
  // that `none' can be output if necessary
  Index hitcount;

  // First check for `all':

  if ( "all" == workspacevariables )
    {
      if (!parameters.plain)
        {
          cout
            << "\n*-------------------------------------------------------------------*\n"
            << "Complete list of ARTS workspace variables:\n"
            << "---------------------------------------------------------------------\n";
        }

      for ( Index i=0; i<Workspace::wsv_data.nelem(); ++i )
        {
          if (!parameters.plain) cout << "- ";
          cout << Workspace::wsv_data[i].Name() << "\n";
        }

      if (!parameters.plain)
        cout << "*-------------------------------------------------------------------*\n\n";
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
      << "\n*-------------------------------------------------------------------*\n"
      << "Generic workspace variables required by " << mdr.Name()
      << " are of type:\n"
      << "---------------------------------------------------------------------\n";
      for ( Index i=0; i<mdr.GInType().nelem(); ++i )
        {
          cout << "- " << wsv_group_names[mdr.GInType()[i]] << "\n";
          ++hitcount;
        }
      if ( 0==hitcount )
        cout << "none\n";

      // List specific variables required by this method.
      hitcount = 0;
      cout 
      << "\n---------------------------------------------------------------------\n"
      << "Specific workspace variables required by " << mdr.Name() << ":\n"
      << "---------------------------------------------------------------------\n";
      for ( Index i=0; i<mdr.Input().nelem(); ++i )
        {
          cout << "- " << Workspace::wsv_data[mdr.Input()[i]].Name() << "\n";
          ++hitcount;
        }
      if ( 0==hitcount )
        cout << "none\n";

      cout
        << "*-------------------------------------------------------------------*\n\n";

      return;
    }

  // If we are here, then the user specified nothing that makes sense.
  cerr << "The name " << workspacevariables << " matches neither `all',\n" 
       << "nor the name of a workspace method.\n";
  arts_exit ();
}


/** React to option `describe'. This should print the description
    String of the given workspace variable or method.

    \param describe What to describe.
    \author Stefan Buehler */
void option_describe(const String& describe)
{
  // Make global data visible:
  extern const Array<MdRecord>  md_data_raw;
  extern const map<String, Index> MdRawMap;
  //  extern const ArrayOfString wsv_group_names;

  // Let's first assume it is a method that the user wants to have
  // described.

  // Find method id:
  map<String, Index>::const_iterator i =
    MdRawMap.find(describe);
  if ( i != MdRawMap.end() )
    {
      // If we are here, then the given name matches a method.
      cout << md_data_raw[i->second] << "\n";
      return;
    }

  // Ok, let's now assume it is a variable that the user wants to have
  // described.

  // Find wsv id:
  i = Workspace::WsvMap.find(describe);
  if ( i != Workspace::WsvMap.end() )
    {
      // If we are here, then the given name matches a workspace
      // variable.
      cout << Workspace::wsv_data[i->second] << "\n";
      return;     
    }

  // If we are here, then the given name does not match anything.
  cerr << "The name " << describe
       << " matches neither method nor variable.\n";
  arts_exit ();      
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
      cout << "\n" << parameters.usage << "\n\n";
      cout << parameters.helptext << "\n\n";
      arts_exit (EXIT_SUCCESS);
    }

  if (parameters.version)
    {
      extern const String full_name;
      // Just print version information and then exit.
      cout << full_name
        << " (compiled " << __DATE__ << " " << __TIME__
        << " on " << OS_NAME << " " << OS_VERSION << ")" << endl
        << "Compile flags: " << COMPILE_FLAGS << endl
        << "Features in this build: " << endl
        << "   Numeric precision:  "
        << ((sizeof (Numeric) == sizeof (double)) ? "double" : "float") << endl
        << "   OpenMP support:     "
#ifdef _OPENMP
        << "enabled" << endl
#else
        << "disabled" << endl
#endif
        << "   Zipped XML support: "
#ifdef ENABLE_ZLIB
        << "enabled" << endl
#else
        << "disabled" << endl
#endif
        << "   Disort algorithm:   "
#ifdef ENABLE_DISORT
        << "enabled" << endl
#else
        << "disabled" << endl
#endif
        << "";
      arts_exit (EXIT_SUCCESS);
    }

  if (parameters.numthreads)
    {
#ifdef _OPENMP
      omp_set_num_threads (parameters.numthreads);
#else
      cerr << "Ignoring commandline option --numthreads/-n.\n"
           << "This option only works with an OpenMP enabled ARTS build.\n";
#endif
    }


  // For the next couple of options we need to have the workspce and
  // method lookup data.

  // Initialize the wsv group name array:
  define_wsv_group_names();

  // Initialize the wsv data:
  Workspace::define_wsv_data();

  // Initialize WsvMap:
  Workspace::define_wsv_map();

  // Initialize the md data:
  define_md_data_raw();

  // Expand supergeneric methods:
  expand_md_data_raw_to_md_data();

  // Initialize MdMap:
  define_md_map();

  // Initialize MdRawMap (needed by parser and online docu).
  define_md_raw_map();

  // Initialize the agenda lookup data:
  define_agenda_data();

  // Initialize AgendaMap:
  define_agenda_map();

  // Check that agenda information in wsv_data and agenda_data is consistent: 
  assert( check_agenda_data() );

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
  //  extern const map<String, Index> MdMap;
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
      arts_exit (EXIT_SUCCESS);
    }

  // React to option `input'. Given the name of a variable (or group)
  // it should print all methods that need this variable (or group) as
  // input.
  if ( "" != parameters.input )
    {
      option_input(parameters.input);
      arts_exit (EXIT_SUCCESS);
    }
  
  // React to option `workspacevariables'. If given the argument `all',
  // it should simply prints a list of all variables. If given the
  // name of a method, it should print all variables that are needed
  // by that method.
  if ( "" != parameters.workspacevariables )
    {
      option_workspacevariables(parameters.workspacevariables);
      arts_exit (EXIT_SUCCESS);
    }

  // React to option `describe'. This should print the description
  // String of the given workspace variable or method.
  if ( "" != parameters.describe )
    {
      option_describe(parameters.describe);
      arts_exit (EXIT_SUCCESS);
    }

  
  // React to option `groups'. This should simply print a list of all
  // workspace variable groups.
  if ( parameters.groups )
    {
      if (!parameters.plain)
        {
          cout
            << "\n*-------------------------------------------------------------------*\n"
            << "Complete list of ARTS workspace variable groups:\n"
            << "---------------------------------------------------------------------\n";
        }

      for ( Index i=0; i<wsv_group_names.nelem(); ++i )
        {
          if (!parameters.plain) cout << "- ";
          cout << wsv_group_names[i] << "\n";
        }

      if (!parameters.plain)
        cout << "*-------------------------------------------------------------------*\n\n";
      arts_exit (EXIT_SUCCESS);
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
      String::size_type p = out_basename.rfind(".arts");

      if (String::npos==p)
        {
          // This is an error handler for the case that somebody gives
          // a supposed file name that does not contain the extension
          // ".arts"

          cerr << "The controlfile must have the extension .arts.\n";
          polite_goodby();
        }
      
      // Kill everything starting from the `.'
      out_basename.erase(p);
    }
  else
    {
      extern String out_basename;
      out_basename = parameters.basename;
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
      extern String out_basename;       // Basis for file name
      extern ofstream report_file;      // Report file pointer
      ostringstream report_file_ext;

      report_file_ext << ".rep";
      open_output_file(report_file, out_basename + report_file_ext.str ());
    }
  catch (runtime_error x)
    {
      cerr << x.what() << "\n"
           << "I have to be able to write to my report file.\n";
      arts_exit ();
    }
  catch (ios_base::failure x)
    {
      cerr << x.what() << "\n"
           << "I have to be able to write to my report file.\n"
           << "Make sure you have write permissions for the directory where\n"
           << "the report file is written.\n";
      arts_exit ();
    }

  // Now comes the global try block. Exceptions caught after this
  // one are general stuff like file opening errors.
  try
    {
      {
        // Output program name and version number: 
        // The name (PACKAGE) and the major and minor version number
        // (VERSION) are set in configure.in. The configuration tools
        // place them in the file config.h, which is included in arts.h.
  
        extern const String full_name;
        out1 << full_name << "\n";
      }

      // Output some OpenMP specific information on output level 2:
#ifdef _OPENMP
      out2 << "Running with OpenMP, "
           << "maximum number of threads = "
           << arts_omp_get_max_threads() << ".\n";
#else
      out2 << "Running without OpenMP.\n";        
#endif

      {
        // Output verbosity settings. This is not too interesting, it
        // goes only to out3.
        extern Messages arts_messages;        
        out3 << "Verbosity settings: Agendas:     " << arts_messages.va << "\n"
             << "                    Screen:      " << arts_messages.vs << "\n"
             << "                    Report file: " << arts_messages.vf << "\n";
      }

      out3 << "\nReading control files:\n";
      for ( Index i=0; i<parameters.controlfiles.nelem(); ++i )
        {
          out3 << "- " << parameters.controlfiles[i] << "\n";

          // The list of methods to execute and their keyword data from
          // the control file. 
          Agenda tasklist;

          Workspace workspace;

          // Call the parser to parse the control text:
          ArtsParser arts_parser(tasklist, parameters.controlfiles[i]);

          arts_parser.parse_tasklist();

          tasklist.set_name("Arts");

          tasklist.set_main_agenda();

          workspace.initialize ();

          // Execute main agenda:
          Arts(workspace, tasklist);
        }
    }
  catch (runtime_error x)
    {
      arts_exit_with_error_message(x.what());
    }

  out1 << "Goodbye.\n";
  arts_exit (EXIT_SUCCESS);
}
