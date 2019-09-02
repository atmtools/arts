#include "timings.h"
#include <iostream>

using std::endl;
using std::left;
using std::setw;

//! Print timing results.
/*!

  Prints the timing results to the given output stream.

  \param os The output stream to write the timin results to.
  \param timer The Timings object.

  \return The output stream.
*/
std::ostream &operator<<(std::ostream &os, const Timings &timings) {
  os << endl << setw(40) << left << "TOTAL TIME:" << timings.total_time << endl;

  for (unsigned int i = 0; i < timings.times.size(); i++) {
    os << setw(40) << left << timings.names[i] << timings.times[i];
    os << endl;
  }

  os << endl << endl;

  return os;
}
