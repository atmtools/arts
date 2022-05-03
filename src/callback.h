#ifndef callback_h
#define callback_h

#include "workspace_ng.h"

#include <functional>

struct CallbackFunction : public std::function<void(Workspace&)> {
  CallbackFunction() : std::function<void(Workspace&)>([](Workspace&){}) {}
  CallbackFunction(std::function<void(Workspace&)> x) : std::function<void(Workspace&)>(x) {}
  friend std::ostream& operator<<(std::ostream& os, const CallbackFunction&) {return os << "Callback";}
};

#endif
