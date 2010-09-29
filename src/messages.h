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
  \file   messages.h
  \brief  Declarations having to do with the four output streams.

  ARTS uses four output streams: \a out0 to \a out3, where \a out0 has
  the highest priority, \a out3 the lowest. These are global
  variables. They are intended to be
  used as follows:

  \arg \a  out0: Error messages
  \arg \a  out1: Output of the `engine'
  \arg \a  out2: Important workspace method output
  \arg \a  out3: Unimportant workspace method output

  The classes associated with the four output stream variables have
  the same name, but start with a capital letter: Out0, Out1, Out2,
  and Out3.

  \author Stefan Buehler
  \date 2000-07-31 */

#ifndef messages_h
#define messages_h

#include <iostream>
#include <fstream>

#include "arts.h"
#include "array.h"
#include "arts_omp.h"

extern bool in_main_agenda;
extern int actual_thread_index;
#ifdef THREADPRIVATE_SUPPORTED
#pragma omp threadprivate(in_main_agenda)
#pragma omp threadprivate(actual_thread_index)
#endif

//! For global ARTS verbosity settings.
/*! 
  This class controls the ARTS verbosity.

  There are four different output streams: out0, out1, out2, and
  out3. They have different priority, out0 the highest, out3 the
  lowest.

  The verbosity level can be set separately for file and screen output. In
  both cases the level can range from 0 to 3, where 0 = no output
  (except error messages), 1 = only out1, 2 = out1+out2, and 3 = all
  output.
   
  For agenda output, the verbosity level can also be between 0 and
  3. The condition for agenda output is evaluated in addition to the
  screen or file output condition. In other words, output from a
  sub-agenda will only be shown on the screen, if its priority is high
  enough for both the agenda verbosity setting, and the screen
  verbosity setting.

  \author Stefan Buehler 
  \date   2008-07-29  */
class Messages {
public:
  bool valid();
  bool sufficient_priority_agenda(Index priority);
  bool sufficient_priority_screen(Index priority);
  bool sufficient_priority_file(Index priority);

/** Print a message to stream and report file. The message is printed
    only if the priority is higher than specified in messages. (Low
    number means high priority.)
  
    \param os       Stream to print to (cout or cerr).
    \param priority Priority of this message (0-3, 0=highest).
    \param t        The stuff to print (can be of any type).
    \author Stefan Buehler 
    \see Messages  */
  template<class T> 
  void Print(ostream& os, Index priority, const T& t)
  {
    extern ofstream report_file;
    
    // cout << "Printing object of type: " << typeid(t).name() << "\n";
    
    // If we are not in the main agenda, then the condition for agenda
    // output must be fulfilled in addition to the condition for
    // screen or file. 

    if (in_main_agenda || sufficient_priority_agenda(priority))
      {
        // We are marking the actual output operations as omp
        // critical, to somewhat reduce the mess when several threads
        // output simultaneously. 

        // This works well if the output operations themselves are
        // atomic, that is if a string is prepared beforehand and then
        // put to outx with a single << operation.
        
        if (sufficient_priority_screen(priority))
          {
#pragma omp critical
            {
              os << t;
            }
          }

        if (sufficient_priority_file(priority))
          {
#pragma omp critical
            {
              //    if (report_file)              // Check if report file is good
              report_file << t << flush;
              // The flush here is necessary to make the output really
              // appear in the report file. We are not producing a huge
              // amount of output to the report file, so I think the
              // performance penalty here is acceptable.
            }
          }
      }
  }
  
  //! Verbosity for agenda output. Can be 0-3.
  Index va;
  //! Verbosity for output to screen. Can be 0-3.
  Index vs;
  //! Verbosity for output to file. Can be 0-3.
  Index vf;
};




/** Highest priority output stream. This stream is only used for error 
    messages and can not be turned off. You should not use this
    directly. Instead, throw an error with an appropriate message.
    \see Out1 Out2 Out3 */
class Out0 {
};

/** Engine output stream. Used  by the engine to report which methods it is 
    calling. Do not use this inside your methods.
    \see Out0 Out2 Out3 */
class Out1 {
};

/** Medium priority output stream. Use this for normal informational
    messages. 
    \see Out0 Out1 Out3  */
class Out2 {
};

/** Lowest priority output stream. This should be used for stuff that
    you would not normally want to see. 
    \see Out0 Out1 Out2  */
class Out3 {
};


//--------------------< Output Operators >--------------------

/** Output operator for Out0. */
template<class T>
Out0& operator<<(Out0& os, const T& t)
{
  extern Messages arts_messages;
  arts_messages.Print(cerr,0,t);
  return os;
}

/** Output operator For Out1. */
template<class T>
Out1& operator<<(Out1& os, const T& t)
{
  //  cout << "Outing Object Of Type: " << Typeid(T).Name() << Endl;
  extern Messages arts_messages;
  arts_messages.Print(cout,1,t);
  return os;
}

/** Output operator For Out2. */
template<class T>
Out2& operator<<(Out2& os, const T& t)
{
  extern Messages arts_messages;
  arts_messages.Print(cout,2,t);
  return os;
}

/** Output operator For Out3. */
template<class T>
Out3& operator<<(Out3& os, const T& t)
{
  extern Messages arts_messages;
  arts_messages.Print(cout,3,t);
  return os;
}

//----------< Make output streams globally visible >----------

extern Out0 out0;
extern Out1 out1;
extern Out2 out2;
extern Out3 out3; 



#endif // messages_h
