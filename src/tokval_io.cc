#include "agenda_class.h"
#include "tokval.h"
#include "tokval_io.h"
#include "tokval_variant.h"

std::ostream& operator<<(std::ostream& os, const TokValPrinter& tvp) {
    return std::visit(
      [&os](auto&& val) -> std::ostream& { return os << *val; },
      *tokval_type(tvp.ref.data()));
}
