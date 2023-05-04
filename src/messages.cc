/**
 * @file   messages.cc
 * @brief  Definitions having to do with the four output streams.
 *
 * See file messages.h for more explanations.
 *
 * @author Stefan Buehler
 * @date   2000-07-31
 */

#include "messages.h"
#include "array.h"
#include "arts.h"
#include "mystring.h"

/** The global message verbosity settings: */
Verbosity verbosity_at_launch;

/** The output path. For example for the report file. */
String out_path;

/** The basename for the report file and for all other output
    files. The files will have names like \<basename\>.rep (for the
    report file). */
String out_basename;

/** The report file. */
ofstream report_file;

ostream& operator<<(ostream& os, const Verbosity& v) {
  os << "Agenda Verbosity: " << v.get_agenda_verbosity() << "\n";
  os << "Screen Verbosity: " << v.get_screen_verbosity() << "\n";
  os << "File Verbosity  : " << v.get_file_verbosity() << "\n";

  return os;
}
