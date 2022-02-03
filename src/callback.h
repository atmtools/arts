#ifndef callback_h
#define callback_h

#include "workspace_ng.h"

#include <functional>

struct CallbackFunction : public std::function<void(void)> {
  CallbackFunction() : std::function<void(void)>([](){}) {}
  CallbackFunction(std::function<void(void)> x) : std::function<void(void)>(x) {}
  friend std::ostream& operator<<(std::ostream& os, const CallbackFunction&) {return os << "Callback";}
};

#endif
