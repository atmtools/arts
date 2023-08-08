#ifndef callback_h
#define callback_h

class Workspace;

#include <functional>
#include <ostream>
#include <stdexcept>

struct CallbackFunction : public std::function<void(const Workspace&)> {
  CallbackFunction() : std::function<void(const Workspace&)>([](const Workspace&) {throw std::runtime_error("Not a function yet");}) {}
  CallbackFunction(std::function<void(const Workspace&)> x);
  friend std::ostream& operator<<(std::ostream& os, const CallbackFunction&);
};

#endif
