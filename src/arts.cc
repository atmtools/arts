/* Copyright (C) 2003-2012 Stefan Buehler <sbuehler@ltu.se>

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
   \file   arts.cc

   This file contains global functions.

   \author Oliver Lemke
   \date   2003-05-07
*/

#include <cstdlib>
#include <stdexcept>
#include <unistd.h>
#include "arts.h"
#include "file.h"
#include "messages.h"


/** This is the exit function of ARTS. Whenever arts has to be terminated
  at some point, call this function.

  You can call without any parameters, since the exit status then
  defaults to EXIT_FAILURE.  

  \param  status  Exit code. EXIT_FAILURE if omitted.
*/
void arts_exit(int status)
{
  extern ofstream report_file;
  extern String out_basename;       // Basis for file name
  ostringstream report_file_ext;
  
  report_file_ext << ".rep";
  cleanup_output_file(report_file, add_basedir(out_basename + report_file_ext.str()));
  
  exit (status);
}

//! Print error message and exit.
/*!
  This function is intended for use in catch blocks.
  
  \param m Error message.
  \param m ArtsOut stream to use.

  \author Stefan Buehler
  \date   2008-05-09
*/
void arts_exit_with_error_message(const String& m, ArtsOut &out)
{
  ostringstream os;

  os << m << "\n";

  if (out.get_verbosity().get_screen_verbosity() < out.get_verbosity().get_file_verbosity())
  {
    extern String out_basename;
    char *buf = getcwd(nullptr, 0);
    os << "Detailed report file: " << buf << "/"
       << add_basedir(out_basename + ".rep") << "\n";
    free(buf);
  }

  os << "Stopping ARTS execution.\n"
     << "Goodbye.\n";
  out << os.str();
  
  arts_exit();              // No argument means failure.
}
