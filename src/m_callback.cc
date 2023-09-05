#include "callback.h"
#include "debug.h"

void CallbackOperatorExecute(const Workspace& ws,
                             const CallbackOperator& op) try {
  //op(ws);
} catch (std::exception& e) {
  throw std::runtime_error(var_string(e.what(), "\n\n", op));
}
