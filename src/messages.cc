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
#include "mystring.h"
#include "array.h"

/** A thread private flag that says whether we are in the main agenda,
    or not. If we are not, then agenda output is suppressed or
    enabled, according to the verbosity setting.

    I also tried to implement this functionality by making the flag an
    array over the number of threads. However, that approach does not
    work, since the thread number is only valid inside the same
    parallel region. Outside it will always be zero, even if we are
    not in the main thread.

    Threadprivate variables are not initialized, so we have to
    explicitly initialize this in main for all threads. */
bool in_main_agenda=true;
/** This is the actual thread index, a threadprivate variable. The OMP 
    function for the thread index works only in the same lexical scope where 
    the split has happened. (Not, for example, in functions called from there.) 
    Like in_main_agenda, this has to be explicitly set for all threads in main. */
int actual_thread_index;
#ifdef THREADPRIVATE_SUPPORTED
#pragma omp threadprivate(in_main_agenda)
#pragma omp threadprivate(actual_thread_index)
#endif

#include "messages.h"

// The global message verbosity settings:
Messages arts_messages;

/** The output path. For example for the report file. */
String out_path;

/** The basename for the report file and for all other output
    files. The files will have names like \<basename\>.rep (for the
    report file). */
String out_basename;

/** The report file. */
ofstream report_file;

/**
   Check if artsmessages contains valid message levels.

   \return True if ok. */
bool Messages::valid()
{
  if (va<0 || va>3) return false;
  if (vs<0 || vs>3) return false;
  if (vf<0 || vf>3) return false;

  return true;
}


//! Does the current message have sufficient priority for agenda?
/*! 
  \param priority Priority of current message.

  \return true if priority is sufficient, otherwise false. */
bool Messages::sufficient_priority_agenda(Index priority)
{
  return va >= priority;
}


//! Does the current message have sufficient priority for screen?
/*! 
  \param priority Priority of current message.

  \return true if priority is sufficient, otherwise false. */
bool Messages::sufficient_priority_screen(Index priority)
{
  return vs >= priority;
}


//! Does the current message have sufficient priority for file?
/*! 
  \param priority Priority of current message.

  \return true if priority is sufficient, otherwise false. */
bool Messages::sufficient_priority_file(Index priority)
{
  return vf >= priority;
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
