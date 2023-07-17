#pragma once

#include <functional>
#include <memory>
#include <ostream>

class Workspace;

struct CallbackFunction : public std::function<void(Workspace&)> {
  CallbackFunction();
  CallbackFunction(std::function<void(Workspace&)> x);
  friend std::ostream& operator<<(std::ostream& os, const CallbackFunction&);
};

struct CallbackOperator {
  std::function<void(std::shared_ptr<Workspace>)> callback;
  std::vector<std::string> inputs;
  std::vector<std::string> outputs;

  void operator()(Workspace& ws) const;

  friend std::ostream& operator<<(std::ostream& os, const CallbackOperator&);
};
