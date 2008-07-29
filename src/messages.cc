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

/*!
  \file   messages.cc
  \brief  Definitions having to do with the four output streams.

  See file messages.h for more explanations.

  \author Stefan Buehler
  \date   2000-07-31
*/

#include "arts.h"
#include "messages.h"
#include "mystring.h"
#include "array.h"

/** The output path. For example for the report file. */
String out_path;

/** The basename for the report file and for all other output
    files. The files will have names like \<basename\>.rep (for the
    report file). */
String out_basename;

/** The report file. */
ofstream report_file;

/** Verbosity levels.
    @see Messages 

    This has to be declared copyin for Open MP loops that involve agendas.
*/
Messages artsmessages;
#pragma omp threadprivate(artsmessages)


bool sufficient_priority_screen(Index priority)
{
  return artsmessages.screen >= priority;
}


bool sufficient_priority_file(Index priority)
{
  return artsmessages.file >= priority;
}

//--------------------< The different output streams >--------------------

/** Level 0 output stream. @see OutStream */
Out0 out0;
/** Level 1 output stream. @see OutStream */
Out1 out1;
/** Level 2 output stream. @see OutStream */
Out2 out2;
/** Level 3 output stream. @see OutStream */
Out3 out3;
