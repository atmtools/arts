/*!
  \file   supergeneric.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Thu Jun 13 14:44:02 2002
  
  \brief  Declarations for supergeneric methods.

  The class Any can be used to mark supergenerio IO.
*/

#ifndef supergeneric_h
#define supergeneric_h

#include <ostream>

//! A placeholder for any type.
/*!  Used to mark supergeneric methods in file
  methods.cc.
*/
class Any {
  // Nothing to do here.
  friend std::ostream& operator<<(std::ostream& os, const Any&) {return os;}
};

#endif  // supergeneric_h
