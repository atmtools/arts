#include "callback.h"
#include "workspace_class.h"

#include <exception>
#include <stdexcept>
#include <ostream>

#include <workspace.h>

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

void CallbackOperator::operator()(Workspace& ws_in) const try {
  ARTS_USER_ERROR_IF(not callback, "No callback function set for operator:\n{}", *this);
  
  Workspace ws(WorkspaceInitialization::Empty);

  for (auto& n : inputs) ws.set(n, ws_in.share(n));
  for (auto& n : outputs) {
    if (ws_in.contains(n)) {
      ws.set(n, ws_in.share(n));
    } else {
      ws.init(n);
      ws_in.set(n, ws.share(n));
    }
  }
  callback(ws);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("Error in callback operator:\n{}\n{}", *this, e.what()));
}
