/* Copyright (C) 2000-2012 Stefan Buehler <sbuehler@ltu.se>

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
#include <memory>

#ifdef ENABLE_DOCSERVER
#include <unistd.h>
#endif

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <map>

#include "absorption.h"
#include "agenda_record.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "docserver.h"
#include "exceptions.h"
#include "file.h"
#include "global_data.h"
#include "methods.h"
#include "mystring.h"
#include "parameters.h"
#include "parser.h"
#include "workspace_global_data.h"
#include "workspace.h"
#include "workspace_ng.h"
#include "wsv_aux.h"
#include "wsv_aux_operator.h"

/** Remind the user of --help and exit return value 1. */
void polite_goodby() {
  cerr << "Try `arts --help' for help.\n";
  arts_exit();
}

/**
   Set the reporting level.

   Set the global reporting level, either the default or based on
   reporting. If reporting was specified, check if the values make
   sense. The value -1 for reporting means that it was (probably)
   not given on the command line, since this is the initialization
   value.   

   \param[in] r Reporting level from Command line.
 
   \author Stefan Buehler 
*/

/** React to option `methods'. If given the argument `all', it
    should simply prints a list of all methods. If given the name of
    a variable, it should print all methods that produce this
    variable as output.

    \param methods All or name of a variable.
    \author Stefan Buehler */
void option_methods(Workspace& workspace, const String& methods) {
  // Make global data visible:
  using global_data::md_data_raw;
  extern const Parameters parameters;
  using global_data::wsv_groups;

  // This is used to count the number of matches to a query, so
  // that `none' can be output if necessary
  Index hitcount;

  // First check if the user gave the special name `all':

  if ("all" == methods) {
    if (!parameters.plain) {
      cout
          << "\n*-------------------------------------------------------------------*\n"
          << "Complete list of ARTS workspace methods:\n"
          << "---------------------------------------------------------------------\n";
    }
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      if (!parameters.plain) cout << "- ";
      cout << md_data_raw[i].Name() << "\n";
    }

    if (!parameters.plain)
      cout
          << "*-------------------------------------------------------------------*\n\n";

    return;
  }

  // Ok, so the user has probably specified a workspace variable or
  // workspace variable group.

  // Check if the user gave the name of a specific variable.
  auto mi = workspace.WsvMap_ptr->find(methods);
  if (mi != workspace.WsvMap_ptr->end()) {
    // If we are here, then the given name matches a variable.
    Index wsv_key = mi->second;

    // List generic methods:
    hitcount = 0;
    cout
        << "\n*-------------------------------------------------------------------*\n"
        << "Generic and supergeneric methods that can generate "
        << (*workspace.wsv_data_ptr)[wsv_key].Name() << ":\n"
        << "---------------------------------------------------------------------\n";
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether GOutType, the list
      // of output variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if (count(mdd.GOutType().begin(),
                mdd.GOutType().end(),
                (*workspace.wsv_data_ptr)[wsv_key].Group())) {
        cout << "- " << mdd.Name() << "\n";
        ++hitcount;
      } else if (count(mdd.GOutType().begin(),
                       mdd.GOutType().end(),
                       get_wsv_group_id("Any"))) {
        for (Index j = 0; j < mdd.GOutType().nelem(); j++) {
          if (mdd.GOutType()[j] == get_wsv_group_id("Any")) {
            if (mdd.GOutSpecType()[j].nelem()) {
              if (count(mdd.GOutSpecType()[j].begin(),
                        mdd.GOutSpecType()[j].end(),
                        (*workspace.wsv_data_ptr)[wsv_key].Group())) {
                cout << "- " << mdd.Name() << "\n";
                ++hitcount;
              }
            } else {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
          }
        }
      }
    }
    if (0 == hitcount) cout << "none\n";

    // List specific methods:
    hitcount = 0;
    cout
        << "\n---------------------------------------------------------------------\n"
        << "Specific methods that can generate "
        << (*workspace.wsv_data_ptr)[wsv_key].Name() << ":\n"
        << "---------------------------------------------------------------------\n";
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if (count(mdd.Out().begin(), mdd.Out().end(), wsv_key)) {
        cout << "- " << mdd.Name() << "\n";
        ++hitcount;
      }
    }
    if (0 == hitcount) cout << "none\n";

    cout
        << "*-------------------------------------------------------------------*\n\n";

    return;
  }

  // Check if the user gave the name of a variable group.

  // We use the find algorithm from the STL to do this. It
  // returns an iterator, so to get the index we take the
  // difference to the begin() iterator.
  Index group_key =
      find(wsv_groups.begin(), wsv_groups.end(), methods) -
      wsv_groups.begin();

  // group_key == wsv_goup_names.nelem() indicates that a
  // group with this name was not found.
  if (group_key != wsv_groups.nelem()) {
    // List generic methods:
    hitcount = 0;
    cout
        << "\n*-------------------------------------------------------------------*\n"
        << "Generic and supergeneric methods that can generate variables of group "
        << wsv_groups[group_key] << ":\n"
        << "---------------------------------------------------------------------\n";
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether GOutType, the list
      // of output variable types contains the
      // requested group.
      // The else clause picks up methods with supergeneric input.
      if (count(mdd.GOutType().begin(), mdd.GOutType().end(), group_key)) {
        cout << "- " << mdd.Name() << "\n";
        ++hitcount;
      } else if (count(mdd.GOutType().begin(),
                       mdd.GOutType().end(),
                       get_wsv_group_id("Any"))) {
        for (Index j = 0; j < mdd.GOutType().nelem(); j++) {
          if (mdd.GOutType()[j] == get_wsv_group_id("Any")) {
            if (mdd.GOutSpecType()[j].nelem()) {
              if (count(mdd.GOutSpecType()[j].begin(),
                        mdd.GOutSpecType()[j].end(),
                        group_key)) {
                cout << "- " << mdd.Name() << "\n";
                ++hitcount;
              }
            } else {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
          }
        }
      }
    }
    if (0 == hitcount) cout << "none\n";

    cout
        << "*-------------------------------------------------------------------*\n\n";

    return;
  }

  // If we are here it means that what the user specified is neither
  // `all', nor a variable, nor a variable group.
  cerr << "The name " << methods << " matches neither `all',\n"
       << "nor the name of a workspace variable, nor the name\n"
       << "of a workspace variable group.\n";
  arts_exit();
}

/** React to option `input'. Given the name of
    a variable, it should print all methods that need this
    variable as input.

    \param input Name of a variable.
    \author Stefan Buehler
    \date   2001-07-24 */
void option_input(const String& input) {
  // Make global data visible:
  using global_data::md_data_raw;
  using global_data::wsv_groups;

  // Ok, so the user has probably specified a workspace variable or
  // workspace variable group.

  // Check if the user gave the name of a specific variable.
  auto mi = global_data::WsvMap.find(input);
  if (mi != global_data::WsvMap.end()) {
    // This is used to count the number of matches to a query, so
    // that `none' can be output if necessary
    Index hitcount = 0;

    // If we are here, then the given name matches a variable.
    Index wsv_key = mi->second;

    // List generic methods:
    cout
        << "\n*-------------------------------------------------------------------*\n"
        << "Generic and supergeneric methods that can use "
        << global_data::wsv_data[wsv_key].Name() << ":\n"
        << "---------------------------------------------------------------------\n";
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether GInType, the list
      // of input variable types contains the group of the
      // requested variable.
      // The else clause picks up methods with supergeneric input.
      if (count(mdd.GInType().begin(),
                mdd.GInType().end(),
                global_data::wsv_data[wsv_key].Group())) {
        cout << "- " << mdd.Name() << "\n";
        ++hitcount;
      } else if (count(mdd.GInType().begin(),
                       mdd.GInType().end(),
                       get_wsv_group_id("Any"))) {
        for (Index j = 0; j < mdd.GInType().nelem(); j++) {
          if (mdd.GInType()[j] == get_wsv_group_id("Any")) {
            if (mdd.GInSpecType()[j].nelem()) {
              if (count(mdd.GInSpecType()[j].begin(),
                        mdd.GInSpecType()[j].end(),
                        global_data::wsv_data[wsv_key].Group())) {
                cout << "- " << mdd.Name() << "\n";
                ++hitcount;
              }
            } else {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
          }
        }
      }
    }
    if (0 == hitcount) cout << "none\n";

    // List specific methods:
    hitcount = 0;
    cout
        << "\n---------------------------------------------------------------------\n"
        << "Specific methods that require "
        << global_data::wsv_data[wsv_key].Name() << ":\n"
        << "---------------------------------------------------------------------\n";
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether Output, the list
      // of output variables contains the workspace
      // variable key.
      if (count(mdd.In().begin(), mdd.In().end(), wsv_key)) {
        cout << "- " << mdd.Name() << "\n";
        ++hitcount;
      }
    }
    if (0 == hitcount) cout << "none\n";

    cout
        << "*-------------------------------------------------------------------*\n\n";

    return;
  }

  // Check if the user gave the name of a variable group.

  // We use the find algorithm from the STL to do this. It
  // returns an iterator, so to get the index we take the
  // difference to the begin() iterator.
  Index group_key =
      find(wsv_groups.begin(), wsv_groups.end(), input) -
      wsv_groups.begin();

  // group_key == wsv_goup_names.nelem() indicates that a
  // group with this name was not found.
  if (group_key != wsv_groups.nelem()) {
    // This is used to count the number of matches to a query, so
    // that `none' can be output if necessary
    Index hitcount = 0;

    // List generic methods:
    cout
        << "\n*-------------------------------------------------------------------*\n"
        << "Generic and supergeneric methods that require a variable of group "
        << wsv_groups[group_key] << ":\n"
        << "---------------------------------------------------------------------\n";
    for (Index i = 0; i < md_data_raw.nelem(); ++i) {
      // Get handle on method record:
      const MdRecord& mdd = md_data_raw[i];

      // This if statement checks whether GOutType, the list
      // of output variable types contains the
      // requested group.
      // The else clause picks up methods with supergeneric input.
      if (count(mdd.GInType().begin(), mdd.GInType().end(), group_key)) {
        cout << "- " << mdd.Name() << "\n";
        ++hitcount;
      } else if (count(mdd.GInType().begin(),
                       mdd.GInType().end(),
                       get_wsv_group_id("Any"))) {
        for (Index j = 0; j < mdd.GInType().nelem(); j++) {
          if (mdd.GInType()[j] == get_wsv_group_id("Any")) {
            if (mdd.GInSpecType()[j].nelem()) {
              if (count(mdd.GInSpecType()[j].begin(),
                        mdd.GInSpecType()[j].end(),
                        group_key)) {
                cout << "- " << mdd.Name() << "\n";
                ++hitcount;
              }
            } else {
              cout << "- " << mdd.Name() << "\n";
              ++hitcount;
            }
          }
        }
      }
    }
    if (0 == hitcount) cout << "none\n";

    cout
        << "*-------------------------------------------------------------------*\n\n";

    return;
  }

  // If we are here it means that what the user specified is neither
  // a variable nor a variable group.
  cerr << "The name " << input << " matches neither the name of a\n"
       << "workspace variable, nor the name of a workspace variable group.\n";
  arts_exit();
}

/** React to option `workspacevariables'. If given the argument `all',
    it should simply prints a list of all variables. If given the
    name of a method, it should print all variables that are needed
    by that method.

    \param  workspacevariables All or name of a method.
    \author Stefan Buehler */
void option_workspacevariables(const String& workspacevariables) {
  // Make global data visible:
  using global_data::md_data;
  using global_data::MdMap;
  extern const Parameters parameters;
  using global_data::wsv_groups;

  // This is used to count the number of matches to a query, so
  // that `none' can be output if necessary
  Index hitcount;

  // First check for `all':

  if ("all" == workspacevariables) {
    if (!parameters.plain) {
      cout
          << "\n*-------------------------------------------------------------------*\n"
          << "Complete list of ARTS workspace variables:\n"
          << "---------------------------------------------------------------------\n";
    }

    for (Index i = 0; i < global_data::wsv_data.nelem(); ++i) {
      if (!parameters.plain) cout << "- ";
      cout << global_data::wsv_data[i].Name() << "\n";
    }

    if (!parameters.plain)
      cout
          << "*-------------------------------------------------------------------*\n\n";
    return;
  }

  // Now check if the user gave the name of a method.
  auto mi = MdMap.find(workspacevariables);
  if (mi != MdMap.end()) {
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
    for (Index i = 0; i < mdr.GInType().nelem(); ++i) {
      cout << "- " << wsv_groups[mdr.GInType()[i]] << "\n";
      ++hitcount;
    }
    if (0 == hitcount) cout << "none\n";

    // List specific variables required by this method.
    hitcount = 0;
    cout
        << "\n---------------------------------------------------------------------\n"
        << "Specific workspace variables required by " << mdr.Name() << ":\n"
        << "---------------------------------------------------------------------\n";
    for (Index i = 0; i < mdr.In().nelem(); ++i) {
      cout << "- " << global_data::wsv_data[mdr.In()[i]].Name() << "\n";
      ++hitcount;
    }
    if (0 == hitcount) cout << "none\n";

    cout
        << "*-------------------------------------------------------------------*\n\n";

    return;
  }

  // If we are here, then the user specified nothing that makes sense.
  cerr << "The name " << workspacevariables << " matches neither `all',\n"
       << "nor the name of a workspace method.\n";
  arts_exit();
}

/** React to option `describe'. This should print the description
    String of the given workspace variable or method.

    \param describe What to describe.
    \author Stefan Buehler */
void option_describe(const String& describe) {
  // Make global data visible:
  using global_data::md_data_raw;
  using global_data::MdRawMap;

  // Let's first assume it is a method that the user wants to have
  // described.

  // Find method id:
  auto i = MdRawMap.find(describe);
  if (i != MdRawMap.end()) {
    // If we are here, then the given name matches a method.
    cout << md_data_raw[i->second] << "\n";
    return;
  }

  // Ok, let's now assume it is a variable that the user wants to have
  // described.

  // Find wsv id:
  i = global_data::WsvMap.find(describe);
  if (i != global_data::WsvMap.end()) {
    // If we are here, then the given name matches a workspace
    // variable.
    cout << global_data::wsv_data[i->second] << "\n";
    return;
  }

  // If we are here, then the given name does not match anything.
  cerr << "The name " << describe << " matches neither method nor variable.\n";
  arts_exit();
}

template <typename TimePoint>
std::time_t to_time_t(TimePoint time_point) {
  auto tp = std::chrono::time_point_cast<std::chrono::system_clock::duration>(
      time_point - TimePoint::clock::now() + std::chrono::system_clock::now());
  return std::chrono::system_clock::to_time_t(tp);
}

/** This function returns the modification time of the arts executable
    as a string.

    \author Oliver Lemke
 */
String arts_mod_time(std::string_view filename) {
  ostringstream os;
  try {
    if (std::filesystem::is_regular_file(filename)) {
      auto modtime = to_time_t(std::filesystem::last_write_time(filename));
      os << " (compiled " << std::put_time(std::localtime(&modtime), "%F %T")
         << ")\n";
    }
  } catch (const std::exception&) {
    // Don't crash if we can't get the modification time.
  }
  return os.str();
}

double get_arts_runtime_in_sec(
    std::chrono::high_resolution_clock::time_point arts_realtime_start) {
  const auto arts_realtime_end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double, std::ratio<1>>(arts_realtime_end -
                                                      arts_realtime_start)
      .count();
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
int main(int argc, char** argv) {
  extern const Parameters parameters;  // Global variable that holds
                                       // all command line parameters.

  //---------------< -1. Time the arts run if possible >---------------
  const auto arts_realtime_start = std::chrono::high_resolution_clock::now();

  //---------------< 1. Get command line parameters >---------------
  if (get_parameters(argc, argv)) {
    // Print an error message and exit:
    polite_goodby();
  }

  //----------< 2. Evaluate the command line parameters >----------
  if (parameters.help) {
    // Just print a help message and then exit.
    cout << "\n" << parameters.usage << "\n\n";
    cout << parameters.helptext << "\n\n";
    arts_exit(EXIT_SUCCESS);
  }

  ostringstream osfeatures;
  {
    osfeatures << "Compiler: " << String(COMPILER) << endl;

    if (String(COMPILER) != "Xcode")
      osfeatures << "Compile flags: " << COMPILE_FLAGS << endl;

    osfeatures << "Features in this build: " << endl
               << "   Numeric precision:    "
               << ((sizeof(Numeric) == sizeof(double)) ? "double" : "float")
               << endl
               << "   OpenMP support:       "
#ifdef _OPENMP
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << "   Documentation server: "
#ifdef ENABLE_DOCSERVER
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << "   Zipped XML support:   "
#ifdef ENABLE_ZLIB
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << "   NetCDF support:       "
#ifdef ENABLE_NETCDF
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << "   Fortran support:      "
#ifdef FORTRAN_COMPILER
               << "enabled (" << FORTRAN_COMPILER << ")" << endl
#else
               << "disabled" << endl
#endif
               << "   Legacy Fortran Disort:"
#ifdef ENABLE_DISORT
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << "   RT4 support:          "
#ifdef ENABLE_RT4
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << "   FASTEM support:       "
#ifdef ENABLE_FASTEM
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << "   OEM support:          "
#ifdef OEM_SUPPORT
               << "enabled" << endl
               << "   MPI support for OEM:  "
#ifdef ENABLE_MPI
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
#else
               << "disabled" << endl
#endif
               << "   Refice support:       "
#ifdef ENABLE_REFICE
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << "   Tmatrix support:      "
#ifdef ENABLE_TMATRIX
#ifdef ENABLE_TMATRIX_QUAD
               << "enabled (quad-precision)" << endl
#else
               << "enabled (double-precision)" << endl
#endif
#else
               << "disabled" << endl
#endif
               << "   IPO/LTO support:      "
#ifdef IPO_SUPPORTED
               << "enabled" << endl
#else
               << "disabled" << endl
#endif
               << endl;

    osfeatures << "Include search paths: " << endl;
    for (auto& path : parameters.includepath) {
      osfeatures << "   " << path << endl;
    }

    osfeatures << "Data searchpaths: " << endl;
    for (auto& path : parameters.datapath) {
      osfeatures << "   " << path << endl;
    }
  }

  if (parameters.version) {
    cout << arts_get_version_string() << arts_mod_time(argv[0]) << endl;
    cout << osfeatures.str();
    arts_exit(EXIT_SUCCESS);
  }

  if (parameters.numthreads) {
#ifdef _OPENMP
    omp_set_num_threads((int)parameters.numthreads);
#else
    cerr << "Ignoring commandline option --numthreads/-n.\n"
         << "This option only works with an OpenMP enabled ARTS build.\n";
#endif
  }

  // For the next couple of options we need to have the workspce and
  // method lookup data.

  // Initialize the wsv group name array:
  define_wsv_groups();

  // Initialize the wsv data:
  define_wsv_data();

  // Initialize WsvMap:
  define_wsv_map();

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
  ARTS_ASSERT(check_agenda_data());

  // Initialize memory handler.
  global_data::workspace_memory_handler.initialize();

  auto workspace_shared = Workspace::create();

  Workspace& workspace=*workspace_shared;

  // Make all global data visible:
  using global_data::wsv_groups;

  // Now we are set to deal with the more interesting command line
  // switches.
#ifdef ENABLE_DOCSERVER
  if (parameters.check_docs) {
    // Check built-in docs and then exit
    const auto broken_links = Docserver::list_broken_description_links();
    const size_t nbroken = std::get<0>(broken_links);
    for (auto&& s : std::get<1>(broken_links)) {
      std::cout << s << std::endl;
    }
    std::cout << std::endl
              << nbroken << " broken link" << (nbroken == 1 ? "" :"s")
              << " found." << std::endl;
    arts_exit(nbroken ? EXIT_FAILURE : EXIT_SUCCESS);
  }
#endif

  // React to option `methods'. If given the argument `all', it
  // should simply prints a list of all methods. If given the name of
  // a variable, it should print all methods that produce this
  // variable as output.
  if ("" != parameters.methods) {
    option_methods(workspace, parameters.methods);
    arts_exit(EXIT_SUCCESS);
  }

  // React to option `input'. Given the name of a variable (or group)
  // it should print all methods that need this variable (or group) as
  // input.
  if ("" != parameters.input) {
    option_input(parameters.input);
    arts_exit(EXIT_SUCCESS);
  }

  // React to option `workspacevariables'. If given the argument `all',
  // it should simply prints a list of all variables. If given the
  // name of a method, it should print all variables that are needed
  // by that method.
  if ("" != parameters.workspacevariables) {
    option_workspacevariables(parameters.workspacevariables);
    arts_exit(EXIT_SUCCESS);
  }

  // React to option `describe'. This should print the description
  // String of the given workspace variable or method.
  if ("" != parameters.describe) {
    option_describe(parameters.describe);
    arts_exit(EXIT_SUCCESS);
  }

  // React to option `groups'. This should simply print a list of all
  // workspace variable groups.
  if (parameters.groups) {
    if (!parameters.plain) {
      cout
          << "\n*-------------------------------------------------------------------*\n"
          << "Complete list of ARTS workspace variable groups:\n"
          << "---------------------------------------------------------------------\n";
    }

    for (Index i = 0; i < wsv_groups.nelem(); ++i) {
      if (!parameters.plain) cout << "- ";
      cout << wsv_groups[i] << "\n";
    }

    if (!parameters.plain)
      cout
          << "*-------------------------------------------------------------------*\n\n";
    arts_exit(EXIT_SUCCESS);
  }

#ifdef ENABLE_DOCSERVER
  if (0 != parameters.docserver) {
    run_docserver(parameters.docserver, parameters.baseurl, parameters.daemon);
    arts_exit(0);
  }
#endif

  // Ok, we are past all the special options. This means the user
  // wants to get serious and really do a calculation. Check if we
  // have at least one control file:
  if (0 == parameters.controlfiles.nelem()) {
    cerr << "You must specify at least one control file name.\n";
    polite_goodby();
  }

  //--------------------< Open report file >--------------------
  // This one needs its own little try block, because we have to
  // write error messages to cerr directly since the report file
  // will not exist.
  try {
    ostringstream report_file_ext;

    report_file_ext << ".rep";
    //open_output_file(report_file, report_file_ext.str());
  } catch (const std::runtime_error& x) {
    cerr << x.what() << "\n"
         << "I have to be able to write to my report file.\n";
    arts_exit();
  }

  // Now comes the global try block. Exceptions caught after this
  // one are general stuff like file opening errors.
  try {
    // Output a short hello from each thread to out3:
#ifdef _OPENMP
#pragma omp parallel default(none)
    {
      ostringstream os;
      int tn = arts_omp_get_thread_num();
      os << "   Thread " << tn << ": ready.\n";
    }
#endif

    time_t rawtime;
    struct tm* timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    for (Index i = 0; i < parameters.controlfiles.nelem(); ++i) {
      try {

        // The list of methods to execute and their keyword data from
        // the control file.
        Agenda tasklist{workspace};

        // Call the parser to parse the control text:
        ArtsParser arts_parser(tasklist, parameters.controlfiles[i]);

        arts_parser.parse_tasklist();

        tasklist.set_name("Arts");

        tasklist.set_main_agenda();

        // Execute main agenda:
        Arts2(workspace, tasklist);
      } catch (const std::exception& x) {
        ostringstream os;
        os << "Run-time error in controlfile: " << parameters.controlfiles[i]
           << '\n'
           << x.what();
        throw runtime_error(os.str());
      }
    }
  } catch (const std::runtime_error& x) {
    arts_exit_with_error_message(x.what());
  }
  arts_exit(EXIT_SUCCESS);
}
