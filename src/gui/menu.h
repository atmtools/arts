#include <species_tags.h>

#include <limits>

#include "atm.h"
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

[[nodiscard]] std::string change_item_name(const JacobianTarget&);
[[nodiscard]] bool change_item(const char *, ArrayOfRetrievalQuantity&, ArrayOfRetrievalQuantity&);

[[nodiscard]] bool change_item(const char *,
                               Vector &,Vector &,
                               const ArrayOfString &keys);
[[nodiscard]] bool change_item(const char *,
                               std::array<Numeric, 3> &,std::array<Numeric, 3> &,
                               const std::array<String, 3> &keys);
[[nodiscard]] bool change_item(const char *,
                               AtmPoint &,AtmPoint &,
                               const ArrayOfArrayOfSpeciesTag &,
                               Options &);
[[nodiscard]] bool change_item(
    const char *,
    Vector &,Vector &,
    Numeric min = std::numeric_limits<Numeric>::lowest(),
    Numeric max = std::numeric_limits<Numeric>::max());
[[nodiscard]] bool change_item(
    const char *,
    std::array<Numeric, 3> &,std::array<Numeric, 3> &,
    Numeric min = std::numeric_limits<Numeric>::lowest(),
    Numeric max = std::numeric_limits<Numeric>::max());

[[nodiscard]] bool change_item(
    const char *,
    const char *,
    Numeric &,Numeric &,
    Numeric min = std::numeric_limits<Numeric>::lowest(),
    Numeric max = std::numeric_limits<Numeric>::max());

[[nodiscard]] bool change_item(const char *,
                               ArrayOfSpeciesTag &,ArrayOfSpeciesTag &,
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
    ImGui::Separator();
  }
}

[[nodiscard]] std::string absunit(const Jacobian::Target& target);
void select_option(Index&, const ArrayOfRetrievalQuantity& jac);

void tooltip(const char*, const Config& config);

struct scoped_disable final {
  bool disabled;
  scoped_disable(bool x) : disabled(x) {
    if (disabled) {
      ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
      ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    }
  }
  ~scoped_disable() noexcept {
    if (disabled) {
      ImGui::PopItemFlag();
      ImGui::PopStyleVar();
    }
  }
};

void show_plot_controls();
void show_propmat_controls();
}  // namespace ARTSGUI::MainMenu
