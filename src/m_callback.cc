#include <exception>
#include <stdexcept>
#include "callback.h"
#include "messages.h"

void CallbackFunctionExecute(Workspace& ws,
                             const CallbackFunction& function,
                             const Verbosity&) {
  function(ws);
}

void CallbackOperatorExecute(Workspace& ws,
                             const CallbackOperator& op,
                             const Verbosity&) try {
  op(ws);
} catch (std::exception& e) {
  throw std::runtime_error(var_string(e.what(), "\n\n", op));
}
