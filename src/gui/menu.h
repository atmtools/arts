#include <species_tags.h>

#include <limits>

#include "debug.h"
#include "enums.h"
#include "gui.h"
#include "imgui.h"
#include "jacobian.h"

namespace ARTSGUI::MainMenu {
ENUMCLASS(VMR, char, exact, percent, ppmv)

struct Options {
  VMR vmr{VMR::exact};
};

void fullscreen(Config &cfg, GLFWwindow *window);

void quitscreen(const Config &cfg, GLFWwindow *window);

bool exportdata(const Config &cfg, ImGui::FileBrowser &fileBrowser, const char * dialog = " Export Data ", bool shortcut=true);

[[nodiscard]] bool change_item(const char *);

[[nodiscard]] bool change_item(const char *, ArrayOfRetrievalQuantity&);

[[nodiscard]] bool change_item(const char *,
                               Vector &,
                               const ArrayOfString &keys);
[[nodiscard]] bool change_item(const char *,
                               Vector &,
                               const ArrayOfArrayOfSpeciesTag &,
                               Options &);
[[nodiscard]] bool change_item(
    const char *,
    Vector &,
    Numeric min = std::numeric_limits<Numeric>::lowest(),
    Numeric max = std::numeric_limits<Numeric>::max());

[[nodiscard]] bool change_item(
    const char *,
    const char *,
    Numeric &,
    Numeric min = std::numeric_limits<Numeric>::lowest(),
    Numeric max = std::numeric_limits<Numeric>::max());

[[nodiscard]] bool change_item(const char *,
                               ArrayOfSpeciesTag &,
                               const ArrayOfArrayOfSpeciesTag &);

template <class T, size_t N>
void select_option(T &current,
                   const std::array<T, N> &opt,
                   const std::array<std::string, N> &disp) {
  for (size_t i = 0; i < N; i++) {
    if (ImGui::Selectable(disp[i].c_str(),
                          current == opt[i],
                          ImGuiSelectableFlags_DontClosePopups)) {
      current = opt[i];
    }
  }
}

void select_option(Index&, const ArrayOfRetrievalQuantity& jac);
}  // namespace ARTSGUI::MainMenu
