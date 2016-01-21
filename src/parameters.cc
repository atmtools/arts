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
   \file   parameters.cc
   
   This file contains the function get_parameters, which reads command
   line parameters. Standard GNU functions are used for this.

   \author Stefan Buehler
   \date   2001-07-24
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "arts.h"
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#else
#include <arts_getopt.h>
#endif
#include "parameters.h"
#include "file.h"

/// Holds the command line parameters.
Parameters parameters;


//! Parse path environment variable
/** 
 Parse a colon separated list of paths from the given environment variable
 into an ArrayOfString.
 
 \param[in]  envvar  Name of environment variable.
 \param[out] paths   ArrayOfString of paths.
 
 \author Oliver Lemke
 */
void parse_path_from_environment (String envvar, ArrayOfString& paths)
{
  char *envval = getenv (envvar.c_str());
  if (envval)
  {
    String pathstring(envval);
    
    // Skip delimiters at beginning.
    String::size_type lastPos = pathstring.find_first_not_of(":", 0);
    // Find first "non-delimiter".
    String::size_type pos     = pathstring.find_first_of(":", lastPos);
    
    while (String::npos != pos || String::npos != lastPos)
    {
      paths.push_back (pathstring.substr (lastPos, pos - lastPos));
      lastPos = pathstring.find_first_not_of(":", pos);
      pos = pathstring.find_first_of(":", lastPos);
    }
  }
}


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
    { "basename",           required_argument, NULL, 'b' },
    { "describe",           required_argument, NULL, 'd' },
    { "groups",             no_argument,       NULL, 'g' },
#ifdef ENABLE_GUI
    { "gui",                no_argument,       NULL, 'G' },
#endif
    { "help",               no_argument,       NULL, 'h' },
    { "includepath",        required_argument, NULL, 'I' },
    { "datapath",           required_argument, NULL, 'D' },
    { "input",              required_argument, NULL, 'i' },
    { "methods",            required_argument, NULL, 'm' },
    { "numthreads",         required_argument, NULL, 'n' },
    { "outdir",             required_argument, NULL, 'o' },
    { "plain",              no_argument,       NULL, 'p' },
    { "reporting",          required_argument, NULL, 'r' },
#ifdef ENABLE_DOCSERVER
    { "docserver",          optional_argument, NULL, 's' },
    { "docdaemon",          optional_argument, NULL, 'S' },
    { "baseurl",            required_argument, NULL, 'U' },
#endif
    { "workspacevariables", required_argument, NULL, 'w' },
    { "version",            no_argument,       NULL, 'v' },
    { NULL,                 no_argument,       NULL, 0   }
  };

  parameters.usage =
    "Usage: arts [-bBdghimnrsSvw]\n"
    "       [--basename <name>]\n"
    "       [--describe <method or variable>]\n"
    "       [--groups]\n"
#ifdef ENABLE_GUI
    "       [--gui]\n"
#endif
    "       [--help]\n"
    "       [--includepath <path>]\n"
    "       [--datapath <path>]\n"
    "       [--input <variable>]\n"
    "       [--methods all|<variable>]\n"
    "       [--numthreads <#>\n"
    "       [--outdir <name>]\n"
    "       [--plain]\n"
    "       [--reporting <xyz>]\n"
#ifdef ENABLE_DOCSERVER
    "       [--docserver[=<port>] --baseurl=BASEURL]\n"
    "       [--docdaemon[=<port>] --baseurl=BASEURL]\n"
#endif
    "       [--workspacevariables all|<method>]\n"
    "       file1.arts file2.arts ...";

  parameters.helptext =
    "The Atmospheric Radiative Transfer Simulator.\n\n"
    "-b, --basename      Set the basename for the report\n"
    "                    file and for other output files.\n"
    "-d, --describe      Print the description String of the given\n"
    "                    workspace variable or method.\n"
    "-g  --groups        List all workspace variable groups.\n"
#ifdef ENABLE_GUI
    "-G  --gui           Start with graphical user interface.\n"
#endif
    "-h, --help          Print this message.\n"
    "-i, --input         This is complementary to the --methods switch.\n"
    "                    It must be given the name of a variable (or group).\n"
    "                    Then it lists all methods that take this variable\n"
    "                    (or group) as input.\n"
    "-I  --includepath   Search path for include files. Can be given more\n"
    "                    than once to add several paths.\n"
    "                    Include paths can also be added by setting the\n"
    "                    environment variable ARTS_INCLUDE_PATH. Multiple\n"
    "                    paths have to be separated by colons.\n"
    "                    Paths specified on the commandline have precedence\n"
    "                    over the environment variable and will be searched\n"
    "                    first.\n"
    "-D  --datapath      Additional search path for data files. Directories\n"
    "                    specified here will be searched after the includepath.\n"
    "                    Data paths can also be added by setting the\n"
    "                    environment variable ARTS_DATA_PATH.\n"
    "-m, --methods       If this is given the argument 'all',\n"
    "                    it simply prints a list of all methods.\n"
    "                    If it is given the name of a variable\n"
    "                    (or variable group), it prints all\n"
    "                    methods that produce this\n"
    "                    variable (or group) as output.\n"
    "-n, --numthreads    If arts was compiled with OpenMP support this option\n"
    "                    can be used to set the maximum number of threads.\n"
    "                    By default OpenMP uses all processors/cores.\n"
    "-o, --outdir        Set the output directory for the report\n"
    "                    file and for other output files with relative paths.\n"
    "                    Default is the current directory.\n"
    "-p  --plain         Generate plain help output suitable for\n"
    "                    script processing.\n"
    "-r, --reporting     Three digit integer. Sets the reporting\n"
    "                    level for agenda calls (first digit),\n"
    "                    screen (second digit) and file (third \n"
    "                    digit). All reporting levels can reach from 0\n"
    "                    (only error messages) to 3 (everything).\n"
    "                    The agenda setting applies in addition to both\n"
    "                    screen and file output.\n"
    "                    Default is 010.\n"
#ifdef ENABLE_DOCSERVER
    "-s, --docserver     Start documentation server. Optionally, specify\n"
    "                    the port number the server should listen on,\n"
    "                    e.g. arts -s9999 or arts --docserver=9999.\n"
    "                    Default is 9000.\n"
    "-S, --docdaemon     Start documentation server in the background.\n"
    "-U, --baseurl       Base URL for the documentation server.\n"
#endif
    "-v, --version       Show version information.\n"
    "-w, --workspacevariables  If this is given the argument 'all',\n"
    "                    it simply prints a list of all variables.\n"
    "                    If it is given the name of a method, it\n"
    "                    prints all variables needed by this method.";

  // Set the short options automatically from the last columns of
  // longopts.
  //
  // Watch out! We also have to put in colons after options that take
  // an argument, and a double colon for options with optional
  // arguments. (Watch out, optional arguments only work with GNU
  // getopt. But since the GNU getopt source code is included with
  // ARTS, why not use this feature?)
  //
  String shortopts;
  {
    int i=0;
    while (NULL != longopts[i].name )
      {
        char c = (char)longopts[i].val;
        shortopts += c;
        //      cout << "name = " << longopts[i].name << "\n";
        //      cout << "val  = " << longopts[i].val << "\n";

        // Check if we need to insert a colon
        if ( required_argument == longopts[i].has_arg )
          {
            shortopts += ":";
          }

        // Or a double colon maybe?
        if ( optional_argument == longopts[i].has_arg )
          {
            shortopts += "::";
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
        case 'b':
          parameters.basename = optarg;
          break;
        case 'd':
          parameters.describe = optarg;
          break;
        case 'g':
          parameters.groups = true;
          break;
        case 'G':
          parameters.gui = true;
          break;
        case 'i':
          parameters.input = optarg;
          break;
        case 'I':
          parameters.includepath.push_back (optarg);
          break;
        case 'D':
          parameters.datapath.push_back (optarg);
          break;
        case 'm':
          parameters.methods = optarg;
          break;
        case 'n':
        {
          istringstream iss(optarg);
          iss >> std::dec >> parameters.numthreads;
          if ( iss.bad() || !iss.eof() )
          {
            cerr << "Argument to --numthreads (-n) must be an integer!\n";
            arts_exit ();
          }
          break;
        }
        case 'o':
          parameters.outdir = optarg;
          break;
        case 'p':
          parameters.plain = true;
          break;
        case 'r':
        {
          //      cout << "optarg = " << optarg << endl;
          istringstream iss(optarg);
          iss >> parameters.reporting;
          ws(iss);
          // This if statement should cover all cases: If there is
          // no integer at all, is becomes bad (first condition). If 
          // there is something else behind the integer, ws does not 
          // reach the end of is (second condition).
          if ( iss.bad() || !iss.eof() )
          {
            cerr << "Argument to --reporting (-r) must be an integer!\n";
            arts_exit ();
          }
          break;
        }
        case 's':
          {
            if (optarg)
            {
              istringstream iss(optarg);
              iss >> std::dec >> parameters.docserver;
              if ( iss.bad() || !iss.eof() )
              {
                cerr << "Argument to --docserver (-s) must be an integer!\n";
                arts_exit ();
              }
            }
            else
              parameters.docserver = -1;
            break;
          }
        case 'S':
          {
            if (optarg)
            {
              istringstream iss(optarg);
              iss >> std::dec >> parameters.docserver;
              if ( iss.bad() || !iss.eof() )
              {
                cerr << "Argument to --docdaemon (-S) must be an integer!\n";
                arts_exit ();
              }
            }
            else
              parameters.docserver = -1;

            parameters.daemon = true;
            break;
          }
        case 'U':
          parameters.baseurl = optarg;
          break;
        case 'v':
          parameters.version = true;
          break;
        case 'w':
          parameters.workspacevariables = optarg;
          break;
        default:
          // There were strange options.
          return(1);
          break;
        }
    }

  // Remaining things on the command line must be control file names. 
  // Get them one by one.
  while ( optind < argc )
    {
      String dummy=argv[optind];
      if (dummy.nelem())
          parameters.controlfiles.push_back(dummy);
      optind++;
    }

  // Look for include paths in the ARTS_INCLUDE_PATH environment variable and
  // append them to parameters.includepath


  parse_path_from_environment(String("ARTS_INCLUDE_PATH"), parameters.includepath);
  parse_path_from_environment(String("ARTS_DATA_PATH"), parameters.datapath);
  
#ifdef ARTS_DEFAULT_INCLUDE_DIR
  String arts_default_include_path (ARTS_DEFAULT_INCLUDE_DIR);
  if (arts_default_include_path != "" && !parameters.includepath.nelem())
  {
    // Skip delimiters at beginning.
    String::size_type lastPos = arts_default_include_path.find_first_not_of(":", 0);
    // Find first "non-delimiter".
    String::size_type pos = arts_default_include_path.find_first_of(":", lastPos);
    
    while (String::npos != pos || String::npos != lastPos)
    {
      parameters.includepath.push_back (arts_default_include_path.substr (lastPos, pos - lastPos));
      lastPos = arts_default_include_path.find_first_not_of(":", pos);
      pos = arts_default_include_path.find_first_of(":", lastPos);
    }
  }
#endif

  parameters.includepath.insert(parameters.includepath.begin(), ".");
  parameters.datapath.insert(parameters.datapath.begin(), ".");
  
  if (parameters.outdir.nelem()) 
    parameters.datapath.insert(parameters.datapath.begin(), parameters.outdir);

  if (parameters.controlfiles.nelem())
  {
    String cfdirname;
    get_dirname(cfdirname, parameters.controlfiles[0]);
    if (cfdirname.nelem()) parameters.includepath.push_back(cfdirname);
  }
  
  return false;
}
