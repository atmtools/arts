#include <algorithm>
#include <array>
#include <ranges>
#include <string_view>

static constexpr std::array keywords{"jacobian",
                                     "spectral_radiance",
                                     "propagation_matrix",
                                     "source_vector",
                                     "nonlte",
                                     "propagation_path",
                                     "grid"};

bool workspace_variables_keywords_match(const std::string_view some_wsv,
                                        const std::string_view this_name) {
  return some_wsv != this_name and
         (some_wsv.find(this_name) != some_wsv.npos or
          std::ranges::any_of(keywords,
                              [a1 = some_wsv, a2 = this_name](auto& keyname) {
                                return a1.find(keyname) != a1.npos and
                                       a2.find(keyname) != a2.npos;
                              }));
}
