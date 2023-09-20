#include <string>

#include "enums.h"
#include "gui.h"

namespace ARTSGUI {
enum class ErrorStatus : char {OnHold, Exit, Continue};
[[nodiscard]] ErrorStatus error(const std::string& errmsg);
}  // namespace ARTSGUI
