#ifndef wsv_aux_operator_h
#define wsv_aux_operator_h

#include "wsv_aux.h"

/** Output operator for WsvRecord.

    \author Stefan Buehler */
std::ostream& operator<<(std::ostream& os, const WsvRecord& wr);

#endif  // wsv_aux_operator_h
