#include "callback.h"

void CallbackFunctionExecute(const Workspace& ws,
                             const CallbackFunction& function) {
  function(ws);
}
