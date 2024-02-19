#include "atm.h"
#include "debug.h"
#include "xml_io.h"

//! AtmFunctionalData

void xml_read_from_stream(std::istream&, AtmFunctionalData&, bifstream*) {
  ARTS_USER_ERROR("Cannot read AtmFunctionalData from XML.");
}

void xml_write_to_stream(std::ostream&,
                         const AtmFunctionalData&,
                         bofstream*,
                         const String&) {
  ARTS_USER_ERROR(
      "Cannot write AtmFunctionalData to XML, please convert it to valid type first.");
}