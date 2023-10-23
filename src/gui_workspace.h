#pragma once

#include <workspace.h>

namespace gui {
//! Experimental support for GUI workspace manipulation
bool change_item(const std::string& name,
                std::shared_ptr<Wsv>& wsv,
                const std::shared_ptr<Wsv>& wsv_old);
}  // namespace gui
