/*!
  \file   m_gridded_fields.h
  \author Lukas Kluft  <lukas.kluft@gmail.com>
  \date   2017-02-01

  \brief  Implementation of GriddedField workspace methods.

  This file contains the implementation of the supergeneric methods.
*/

#ifndef m_gridded_fields_h
#define m_gridded_fields_h

#include "gridded_fields.h"
#include "mystring.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void GriddedFieldGetName(  // WS Generic Output:
    String& name,
    // WS Generic Input:
    const T& gf,
    const Verbosity&) {
  // Return the name of the given GriddedField.
  name = gf.get_name();
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ArrayOfGriddedFieldGetNames(  // WS Generic Output:
    ArrayOfString& names,
    // WS Generic Input:
    const Array<T>& aogf,
    const Verbosity&) {
  // Return the name of the given GriddedField.
  names.resize(aogf.nelem());
  for (Index i = 0; i < aogf.nelem(); i++) {
    names[i] = aogf[i].get_name();
  }
}

#endif  // m_gridded_fields_h
