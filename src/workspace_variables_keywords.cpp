#include "workspace_variables_keywords.h"

#include <algorithm>
#include <ranges>
#include <string_view>
#include <vector>

#include "workspace_variable_shortnames.h"

bool workspace_variables_keywords_match(const std::string_view some_wsv,
                                        const std::string_view this_name) {
  static const std::vector<std::string> keywords = {
      std::from_range, workspace_variables_shortnames() | std::views::keys};

  return some_wsv != this_name and
         (some_wsv.find(this_name) != some_wsv.npos or
          std::ranges::any_of(keywords,
                              [a1 = some_wsv, a2 = this_name](auto& keyname) {
                                return a1.find(keyname) != a1.npos and
                                       a2.find(keyname) != a2.npos;
                              }));
}
