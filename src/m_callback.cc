#include "callback.h"
#include "messages.h"

void CallbackFunctionExecute(Workspace& ws,
                             const CallbackFunction& function,
                             const Verbosity&) {
  function(ws);
}
