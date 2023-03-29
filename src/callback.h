#ifndef callback_h
#define callback_h

class Workspace;

#include <functional>
#include <ostream>

struct CallbackFunction : public std::function<void(Workspace&)> {
  CallbackFunction();
  CallbackFunction(std::function<void(Workspace&)> x);
  friend std::ostream& operator<<(std::ostream& os, const CallbackFunction&);
};

#endif
