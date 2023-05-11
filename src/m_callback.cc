#include "callback.h"

void CallbackFunctionExecute(Workspace& ws,
                             const CallbackFunction& function) {
  function(ws);
}
