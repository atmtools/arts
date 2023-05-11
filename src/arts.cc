/**
   \file   arts.cc

   This file contains global functions.

   \author Oliver Lemke
   \date   2003-05-07
*/

#include "arts.h"
#include <unistd.h>
#include <cstdlib>
#include <stdexcept>
#include "file.h"

/** This is the exit function of ARTS. Whenever arts has to be terminated
  at some point, call this function.

  You can call without any parameters, since the exit status then
  defaults to EXIT_FAILURE.  

  \param  status  Exit code. EXIT_FAILURE if omitted.
*/
void arts_exit(int status) {
  exit(status);
}

//! Print error message and exit.
/*!
  This function is intended for use in catch blocks.

  \param m Error message.
  \param m ArtsOut stream to use.

  \author Stefan Buehler
  \date   2008-05-09
*/
void arts_exit_with_error_message(const String &m) {
  std::cerr << m << "\n"
            << "Stopping ARTS execution.\n"
            << "Goodbye.\n";
  arts_exit(); // No argument means failure.
}
