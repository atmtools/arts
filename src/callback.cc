#include "callback.h"

#include "workspace_ng.h"

CallbackFunction::CallbackFunction()
    : std::function<void(Workspace&)>([](Workspace&) {}) {}

CallbackFunction::CallbackFunction(std::function<void(Workspace&)> x)
    : std::function<void(Workspace&)>(std::move(x)) {}

std::ostream& operator<<(std::ostream& os, const CallbackFunction&) {
  return os << "Callback";
}

std::ostream& operator<<(std::ostream& os, const CallbackOperator& op) {
  os << "CallbackOperator\nInputs: [";
  for (auto& n: op.inputs) {
    os << n << ", ";
  }
  os << "]\nOutputs: [";
  for (auto& n: op.outputs) {
    os << n << ", ";
  }
  return os << ']';
}

void CallbackOperator::operator()(Workspace& ws_in) const {
  auto ws = Workspace::create();

  for (auto& n : inputs) {
    ws->set_wsv(ws_in.copy_wsv(n));
  }

  callback(*ws);

  for (auto& n : outputs) {
    ws_in.set_wsv(ws->copy_wsv(n));
  }
}
