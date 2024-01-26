#pragma once

#include <string_view>

bool workspace_variables_keywords_match(const std::string_view var,
                                        const std::string_view exclude);