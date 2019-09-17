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
 * @file   messages.h
 * @brief  Declarations having to do with the four output streams.
 *
 * ARTS uses four output streams: \a out0 to \a out3, where \a out0 has
 * the highest priority, \a out3 the lowest. These are global variables.
 * They are intended to be used as follows:
 *
 * @arg \a  out0: Error messages.
 * @arg \a  out1: Output of the 'engine'.
 * @arg \a  out2: Important workspace method output.
 * @arg \a  out3: Unimportant workspace method output.
 *
 * The classes associated with the four output stream variables have
 * the same name, but start with a capital letter: Out0, Out1, Out2,
 * and Out3.
 *
 * @author Stefan Buehler
 * @date 2000-07-31
 */

#ifndef messages_h
#define messages_h

#include <fstream>
#include <iostream>

#include "array.h"
#include "arts.h"
#include "arts_omp.h"

class Verbosity {
 public:
  Verbosity() : va(0), vs(0), vf(0), in_main_agenda(false) {}

  Verbosity(Index vagenda, Index vscreen, Index vfile)
      : va(vagenda), vs(vscreen), vf(vfile), in_main_agenda(false) {}

  /**
   * Check if artsmessages contains valid message levels.
   *
   * @return True if ok.
   */
  bool valid() const {
    return (va >= 0 && va <= 3) && (vs >= 0 && vs <= 3) && (vf >= 0 || vf <= 3);
  }

  Index get_agenda_verbosity() const { return va; }
  Index get_screen_verbosity() const { return vs; }
  Index get_file_verbosity() const { return vf; }
  bool is_main_agenda() const { return in_main_agenda; }

  void set_agenda_verbosity(Index v) { va = v; }
  void set_screen_verbosity(Index v) { vs = v; }
  void set_file_verbosity(Index v) { vf = v; }
  void set_main_agenda(bool main_agenda) { in_main_agenda = main_agenda; }

  friend ostream& operator<<(ostream& os, const Verbosity& v);

 private:
  /** Verbosity for agenda output. Can be 0-3.*/
  Index va;
  /** Verbosity for output to screen. Can be 0-3.*/
  Index vs;
  /** Verbosity for output to file. Can be 0-3.*/
  Index vf;
  bool in_main_agenda;
};

class ArtsOut {
 public:
  ArtsOut(const int p, const Verbosity& v) : verbosity(v), priority(p) {}

  int get_priority() const { return priority; }
  const Verbosity& get_verbosity() const { return verbosity; }

  /** Does the current message have sufficient priority for output?
   *
   * @return true if priority is sufficient, otherwise false.
   */
  bool sufficient_priority() const {
    return (sufficient_priority_agenda() &&
            (sufficient_priority_screen() || sufficient_priority_file()));
  }

  /** Does the current message have sufficient priority for agenda?
   *
   * @return true if priority is sufficient, otherwise false.
   */
  bool sufficient_priority_agenda() const {
    return (in_main_agenda() || verbosity.get_agenda_verbosity() >= priority);
  }

  /** Does the current message have sufficient priority for screen?
   *
   * @return true if priority is sufficient, otherwise false.
   */
  bool sufficient_priority_screen() const {
    return verbosity.get_screen_verbosity() >= priority;
  }

  /** Does the current message have sufficient priority for file?
   *
   * @return true if priority is sufficient, otherwise false.
   */
  bool sufficient_priority_file() const {
    return verbosity.get_file_verbosity() >= priority;
  }

  /** Are we in the main agenda?
   *
   * @return true if in main agenda, otherwise false.
   */
  bool in_main_agenda() const { return verbosity.is_main_agenda(); }

 private:
  const Verbosity& verbosity;
  int priority;
};

class ArtsOut0 : public ArtsOut {
 public:
  ArtsOut0(const Verbosity& v) : ArtsOut(0, v) {}
};

class ArtsOut1 : public ArtsOut {
 public:
  ArtsOut1(const Verbosity& v) : ArtsOut(1, v) {}
};

class ArtsOut2 : public ArtsOut {
 public:
  ArtsOut2(const Verbosity& v) : ArtsOut(2, v) {}
};

class ArtsOut3 : public ArtsOut {
 public:
  ArtsOut3(const Verbosity& v) : ArtsOut(3, v) {}
};

/** Output operator for ArtsOut. */
template <class T>
ArtsOut& operator<<(ArtsOut& aos, const T& t) {
  extern ofstream report_file;

  // cout << "Printing object of type: " << typeid(t).name() << "\n";

  // If we are not in the main agenda, then the condition for agenda
  // output must be fulfilled in addition to the condition for
  // screen or file.

  if (aos.sufficient_priority_agenda()) {
    // We are marking the actual output operations as omp
    // critical, to somewhat reduce the mess when several threads
    // output simultaneously.

    // This works well if the output operations themselves are
    // atomic, that is if a string is prepared beforehand and then
    // put to outx with a single << operation.

    if (aos.sufficient_priority_screen()) {
#pragma omp critical(ArtsOut_screen)
      {
        if (aos.get_priority() == 0)
          cerr << t << flush;
        else
          cout << t << flush;
      }
    }

    if (aos.sufficient_priority_file()) {
#pragma omp critical(ArtsOut_file)
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

  return aos;
}

#define CREATE_OUT0 ArtsOut0 out0(verbosity)
#define CREATE_OUT1 ArtsOut1 out1(verbosity)
#define CREATE_OUT2 ArtsOut2 out2(verbosity)
#define CREATE_OUT3 ArtsOut3 out3(verbosity)

#define CREATE_OUTS         \
  ArtsOut0 out0(verbosity); \
  ArtsOut1 out1(verbosity); \
  ArtsOut2 out2(verbosity); \
  ArtsOut3 out3(verbosity)

#endif  // messages_h
