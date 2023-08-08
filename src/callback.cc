#include "callback.h"

CallbackFunction::CallbackFunction(std::function<void(const Workspace&)> x)
    : std::function<void(const Workspace&)>(std::move(x)) {}

std::ostream& operator<<(std::ostream& os, const CallbackFunction&) {
  return os << "Callback";
}
