/**
   \file   parameters.cc
   
   This file contains the function get_parameters, which reads command
   line parameters. Standard GNU functions are used for this.

   \author Stefan Buehler
   \date   2001-07-24
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "file.h"
#include "parameters.h"

/// Holds the command line parameters.
extern Parameters parameters{};

//! Parse path environment variable
/** 
 Parse a colon separated list of paths from the given environment variable
 into an ArrayOfString.
 
 \param[in]  envvar  Name of environment variable.
 \param[out] paths   ArrayOfString of paths.
 
 \author Oliver Lemke
 */
void parse_path_from_environment(String envvar, ArrayOfString &paths) {
  char *envval = getenv(envvar.c_str());
  if (envval) {
    String pathstring(envval);

    // Skip delimiters at beginning.
    String::size_type lastPos = pathstring.find_first_not_of(":", 0);
    // Find first "non-delimiter".
    String::size_type pos = pathstring.find_first_of(":", lastPos);

    while (String::npos != pos || String::npos != lastPos) {
      paths.push_back(pathstring.substr(lastPos, pos - lastPos));
      lastPos = pathstring.find_first_not_of(":", pos);
      pos = pathstring.find_first_of(":", lastPos);
    }
  }
}


