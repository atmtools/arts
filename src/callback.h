#pragma once

#include <functional>
#include <ostream>

class Workspace;

struct CallbackFunction : public std::function<void(Workspace&)> {
  CallbackFunction();
  CallbackFunction(std::function<void(Workspace&)> x);
  friend std::ostream& operator<<(std::ostream& os, const CallbackFunction&);
};

struct CallbackOperator {
  CallbackFunction callback;
  std::vector<std::string> inputs;
  std::vector<std::string> outputs;

  CallbackOperator() = default;

  CallbackOperator(CallbackFunction cb,
                   std::vector<std::string> i,
                   std::vector<std::string> o) : callback(std::move(cb)), inputs(std::move(i)), outputs(std::move(o)) {}

  void operator()(Workspace& ws) const;

  friend std::ostream& operator<<(std::ostream& os, const CallbackOperator&);
};
