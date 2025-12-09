////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file  double_imanip.h

   Fast double input stream with support for parsing nan and inf.

   \author Oliver Lemke, Richard Larsson
   \date 2022-05-23
*/

#ifndef double_imanip_h
#define double_imanip_h

/** Input manipulator class for doubles to enable nan and inf parsing. */
#include <iosfwd>

class double_imanip {
 public:
  const double_imanip& operator>>(double& x) const;

  std::istream& operator>>(const double_imanip&) const;

  friend const double_imanip& operator>>(std::istream& in,
                                         const double_imanip& dm);

 private:
  mutable std::istream* in;
};

#endif
