#include "callback.h"

CallbackFunction::CallbackFunction()
    : std::function<void(Workspace&)>([](Workspace&) {}) {}

CallbackFunction::CallbackFunction(std::function<void(Workspace&)> x)
    : std::function<void(Workspace&)>(std::move(x)) {}

std::ostream& operator<<(std::ostream& os, const CallbackFunction&) {
  return os << "Callback";
}
