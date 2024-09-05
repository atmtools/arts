#include "callback.h"
#include "debug.h"
#include "operators.h"
#include "xml_io.h"

//=== CallbackOperator =========================================

void xml_read_from_stream(std::istream&,
                          CallbackOperator&,
                          bifstream* /* pbifs */) {
  ARTS_USER_ERROR("Reading operators from XML is not supported.");
}

void xml_write_to_stream(std::ostream&,
                         const CallbackOperator&,
                         bofstream* /* pbofs */,
                         const String& /* name */) {
  ARTS_USER_ERROR("Writing operators to XML is not supported.");
}

//=== NumericUnaryOperator =========================================

void xml_read_from_stream(std::istream&,
                          NumericUnaryOperator&,
                          bifstream* /* pbifs */) {
  ARTS_USER_ERROR("Reading operators from XML is not supported.");
}

void xml_write_to_stream(std::ostream&,
                         const NumericUnaryOperator&,
                         bofstream* /* pbofs */,
                         const String& /* name */) {
  ARTS_USER_ERROR("Writing operators to XML is not supported.");
}

//=== NumericBinaryOperator =========================================

void xml_read_from_stream(std::istream&,
                          NumericBinaryOperator&,
                          bifstream* /* pbifs */) {
  ARTS_USER_ERROR("Reading operators from XML is not supported.");
}

void xml_write_to_stream(std::ostream&,
                         const NumericBinaryOperator&,
                         bofstream* /* pbofs */,
                         const String& /* name */) {
  ARTS_USER_ERROR("Writing operators to XML is not supported.");
}

//=== NumericTernaryOperator =========================================

void xml_read_from_stream(std::istream&,
                          NumericTernaryOperator&,
                          bifstream* /* pbifs */) {
  ARTS_USER_ERROR("Reading operators from XML is not supported.");
}

void xml_write_to_stream(std::ostream&,
                         const NumericTernaryOperator&,
                         bofstream* /* pbofs */,
                         const String& /* name */) {
  ARTS_USER_ERROR("Writing operators to XML is not supported.");
}
