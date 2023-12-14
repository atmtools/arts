#include <string>

namespace gui {
enum class ErrorStatus : char {OnHold, Exit, Continue};
[[nodiscard]] ErrorStatus error(const std::string& errmsg);
}  // namespace gui
