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
   \file   parameters.h

   This file contains header information for the dealing with command
   line parameters.

   \author Stefan Buehler
   \date   2001-07-24
*/

#ifndef parameters_h
#define parameters_h

#include "array.h"
#include "mystring.h"

/**
   Structure to hold all command line Parameters. This holds all the
   command line parameters, plut the usage message and the helptext
   message. The messages are in the same structure, because they need
   to be changed whenever the parameters are changed, so it is better
   to have them in the same place. 
   @author SAB
*/
class Parameters {
 public:
  /** Default constructor. Care has to be taken to properly initialize
      all variables, e.g., bool options to false. */
  Parameters()
      : usage(),
        helptext(),
        help(false),
        version(false),
        basename(""),
        outdir(""),
        controlfiles(),
        reporting(-1),
        methods(""),
        numthreads(0),
        includepath(),
        datapath(),
        input(""),
        workspacevariables(""),
        describe(""),
        groups(false),
        plain(false),
        docserver(0),
        baseurl(""),
        daemon(false),
        gui(false) { /* Nothing to be done here */
  }

  /** Short message how to call the program. */
  String usage;
  /** Longer message explaining the options. */
  String helptext;
  /** Only display the help text. */
  bool help;
  /** Display version information. */
  bool version;
  /** If this is specified (with the -b --basename option), it is used
      as the base name for the report file and for other output
      files. */
  String basename;
  /** If this is specified (with the -o --outdir option), it is used
   as the base directory for the report file and for other output
   files. If a full path is given for an output file it will not
   be affected by this. */
  String outdir;
  /** The filenames of the controlfiles. Can be only one or as many as
      you want. */
  ArrayOfString controlfiles;
  /** This should be a two digit integer. The first digit specifies
      the output level for stdout (stderr for error messages), the
      second digit the output level for the report file. The levels
      can reach from 0 (show only error messages) to 3 (show
      everything). Example:

      03 = only errors to the screen, everything to the file. */
  Index reporting;
  /** If this is given the argument `all', it simply prints a list of 
      all methods. If it is given the name of a variable (or group), it
      prints all methods that produce this variable (or group) as output. */
  String methods;
  /** The maximum number of threads to use. */
  Index numthreads;
  /** List of paths to search for include files. */
  ArrayOfString includepath;
  /** List of paths to search for data files. */
  ArrayOfString datapath;
  /** This is complementary to the methods switch. It must be given
      the name of a variable (or group). Then it lists all methods that take this
      variable (or group) as input. */
  String input;
  /** If this is given the argument `all', it simply prints a list of 
      all workspace variables. If it is given the name of a method,
      it prints all variables needed by that method. */
  String workspacevariables;
  /** Print the description String of the given workspace variable or
      method. */
  String describe;
  /** Print a list of all workspace variable groups. */
  bool groups;
  /** Generate plain help out suitable for script processing. */
  bool plain;
  /** Port to use for the docserver. */
  Index docserver;
  /** Baseurl for the docserver. */
  String baseurl;
  /** Flag to run the docserver in the background. */
  bool daemon;
  /** Flag to run with graphical user interface. */
  bool gui;
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

#endif  // parameters_h
