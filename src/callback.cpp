#include "callback.h"
#include "compare.h"

#include <exception>
#include <stdexcept>

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
  ARTS_USER_ERROR_IF(not callback, "No callback function set for operator:\n", *this);
  
  auto ws = std::make_shared<Workspace>(WorkspaceInitialization::Empty);

  for (auto& n : inputs) {
    if (std::ranges::any_of(outputs, Cmp::eq(n))) {
      ws -> set(n, ws_in.share(n));
    } else {
      ws -> set(n, ws_in.copy(n));
    }
  }

  callback(ws);

  for (auto& n : outputs) {
    if (std::ranges::none_of(inputs, Cmp::eq(n))) ws_in.set(n, ws -> share(n));
  }
} catch (std::exception& e) {
  throw std::runtime_error(var_string("Error in callback operator:\n", *this, '\n', e.what()));
}
