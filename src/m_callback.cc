#include "callback.h"
#include "messages.h"

void CallbackFunctionExecute(const CallbackFunction& function,
                             const Verbosity&) {
  function();
}
