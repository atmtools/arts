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

#include <fstream>
#include "arts.h"

/**
   The verbosity level for screen and file output. There are four
   different output streams: out0, out1, out2, and out3. They have different
   priority, out0 the highest, out3 the lowest.

   The verbosity level is stored in the workspace variable messages of
   this type and can be set separately for file and screen output. In
   both cases the level can range from 0 to 3, where 0 = no output
   (except error messages), 1 = only out1, 2 = out1+out2, 3 = all
   output.

   \author Stefan Buehler 
   \date   1999-06-30

   Introduced Out0, which prints to stderr rather than stdout.
   \author Stefan Buehler 
   \date   1999-11-03 */
class Messages {
public:
  /** Default constructor. Set both levels to -1. The default output
      level is set explicitly in main. */
  Messages() : screen(-1), file(-1) { /* nothing to do here */ }
  /** Verbosity of screen output. */
  Index screen;
  /** Verbosity of file output. */
  Index file;
};


/** Print a message to stream and report file. The message is printed
    only if the priority is higher 
    than specified in messages.
  
    \param os       Stream to print to (cout or cerr).
    \param priority Priority of this message (0-3, 0=highest).
    \param t        The stuff to print (can be of any type).
    \author Stefan Buehler 
    \see Messages  */
template<class T> 
void MessagePrint(ostream& os, Index priority, const T& t)
{
  extern Messages messages;
  extern ofstream report_file;

  // cout << "Printing object of type: " << typeid(t).name() << endl;

  if (messages.screen >= priority)
    os << t;

  if (messages.file >= priority)
    //    if (report_file)              // Check if report file is good
    {
      report_file << t;
      report_file.flush ();
    }
}

/* Parent class for all output streams. These control the level of detail
    for your output. There are four output streams altogether, out0 to
    out3. Out0 means highest priority, out3 lowest. Out0 and out1
    should not be used within methods, only by the ARTS engine.

    Nothing is done in this class, it just servers to group the four
    output stream classes together.

    \see Messages */
//class OutStream {
//};


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
  MessagePrint(cerr,0,t);
  return os;
}

/** Output operator For Out1. */
template<class T>
Out1& operator<<(Out1& os, const T& t)
{
  //  cout << "Outing Object Of Type: " << Typeid(T).Name() << Endl;
  MessagePrint(cout,1,t);
  return os;
}

/** Output operator For Out2. */
template<class T>
Out2& operator<<(Out2& os, const T& t)
{
  MessagePrint(cout,2,t);
  return os;
}

/** Output operator For Out3. */
template<class T>
Out3& operator<<(Out3& os, const T& t)
{
  MessagePrint(cout,3,t);
  return os;
}

//----------< Make output streams globally visible >----------

extern Out0 out0;
extern Out1 out1;
extern Out2 out2;
extern Out3 out3; 



#endif // messages_h
