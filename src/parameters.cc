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
   \file   parameters.cc
   
   This file contains the function get_parameters, which reads command
   line parameters. Standard GNU functions are used for this.

   \author Stefan Buehler
   \date   2001-07-24
*/

#include <iostream>
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
    { "help",               no_argument,       NULL, 'h' },
    { "version",            no_argument,       NULL, 'v' },
    { "basename",           required_argument, NULL, 'b' },
    { "reporting",          required_argument, NULL, 'r' },
    { "methods",            required_argument, NULL, 'm' },
    { "input",              required_argument, NULL, 'i' },
    { "workspacevariables", required_argument, NULL, 'w' },
    { "describe",           required_argument, NULL, 'd' },
    { "groups",             no_argument,       NULL, 'g' },
    { "plain",              no_argument,       NULL, 'p' },
    { NULL,                 no_argument,       NULL, 0   }
  };

  parameters.usage =
    "Usage: arts [-hvbrmiwdg] [--help] [--version] [--basename <name>]\n"
    "       [--reporting xy]\n"
    "       [--methods all|<variable>]\n"
    "       [--input <variable>]\n"
    "       [--workspacevariables all|<method>]\n"
    "       [--describe <method or variable>]\n"
    "       [--groups]\n"
    "       [--plain]\n"
    "       file1.arts file2.arts ...";

  parameters.helptext =
    "The Atmospheric Radiative Transfer System.\n\n"
    "-h, --help          Print this message.\n"
    "-v, --version       Show version information.\n"
    "-b, --basename      Set the basename for the report\n"
    "                    file and for other output files.\n"
    "-r, --reporting     Two digit integer. Sets the reporting\n"
    "                    level for screen (first digit) anf file\n"
    "                    (second digit). Levels can reach from 0\n"
    "                    (only error messages) to 3 (everything).\n"
    "-m, --methods       If this is given the argument `all',\n"
    "                    it simply prints a list of all methods.\n"
    "                    If it is given the name of a variable\n"
    "                    (or variable group), it prints all\n"
    "                    methods that produce this\n"
    "                    variable (or group) as output.\n"
    "-i, --input         This is complementary to the --methods switch.\n"
    "                    It must be given the name of a variable (or group).\n"
    "                    Then it lists all methods that take this variable\n"
    "                    (or group) as input.\n"
    "-w, --workspacevariables  If this is given the argument `all',\n"
    "                    it simply prints a list of all variables.\n"
    "                    If it is given the name of a method, it\n"
    "                    prints all variables needed by this method.\n"
    "-d, --describe      Print the description String of the given\n"
    "                    workspace variable or method.\n"
    "-g  --groups        List all workspace variable groups.\n"
    "-p  --plain         Generate plain help output suitable for\n"
    "                    script processing.";

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
        // FIXME: Should there be a cast here? The real type of val is int.
        char c = longopts[i].val;
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
        case 'v':
          parameters.version = true;
          break;
        case 'b':
          parameters.basename = optarg;
          break;
        case 'r':
          {
            //      cout << "optarg = " << optarg << endl;
            istringstream is(optarg);
            is >> parameters.reporting;
            ws(is);
            // This if statement should cover all cases: If there is
            // no integer at all, is becomes bad (first condition). If 
            // there is something else behind the integer, ws does not 
            // reach the end of is (second condition).
            if ( !is || !is.eof() )
              {
                cerr << "Argument to --reporting (-r) must be an integer!\n";
                arts_exit ();
              }
            break;
          }
        case 'm':
          parameters.methods = optarg;
          break;
        case 'i':
          parameters.input = optarg;
          break;
        case 'w':
          parameters.workspacevariables = optarg;
          break;
        case 'd':
          parameters.describe = optarg;
          break;
        case 'g':
          parameters.groups = true;
          break;
        case 'p':
          parameters.plain = true;
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
      parameters.controlfiles.push_back(dummy);
      optind++;
    }

  //  cout << "alle:\n" << parameters.controlfiles << '\n';



  return false;
}
