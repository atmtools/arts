#include "gui_workspace.h"
#include "gui.h"

#include <memory>
#include <type_traits>

namespace gui {
template <WorkspaceGroup T>
struct GuiWsg {
  static constexpr bool enabled = false;
  static constexpr bool change_item(T&) {
    return false;
  }
};

template <>
struct GuiWsg<Index> {
  static constexpr bool enabled = true;
  static bool change_item(Index& n) {
    return ImGui::InputScalar("\tValue:\t", ImGuiDataType_S64, &n);
  }
};

template <>
struct GuiWsg<Numeric> {
  static constexpr bool enabled = true;
  static bool change_item(Numeric& n) {
    return ImGui::InputDouble("\tValue:\t", &n, 0, 0, "%g");
  }
};

bool change_item(const std::string& name,
                 std::shared_ptr<Wsv>& wsv,
                 const std::shared_ptr<Wsv>& wsv_old) {
  const auto enabler = [](auto& wsg) {
    using impl = GuiWsg<std::remove_cvref_t<decltype(*wsg)>>;
    return impl::enabled;
  };
  const auto changer = [](auto& wsg) {
    using impl = GuiWsg<std::remove_cvref_t<decltype(*wsg)>>;
    return impl::change_item(*wsg);
  };

  bool change = false;
  if (ImGui::BeginMenu(name.c_str(), std::visit(enabler, wsv->value))) {
    change = std::visit(changer, wsv->value);
    if (ImGui::MenuItem("Reset")) {
      wsv = std::make_shared<Wsv>(wsv_old->copy()); //! FIXME: Better to get original than to copy?
      change = true;
    }
    ImGui::Separator();
    ImGui::EndMenu();
  }
  return change;
}

bool change_item(Workspace& ws, const Workspace& ws_old) {
  bool change = false;
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Workspace")) {
      for (auto& wsv : ws) {
        change = change or change_item(wsv.first, wsv.second, ws_old.share(wsv.first));
      }
    ImGui::Separator();
    ImGui::EndMenu();
    }
  ImGui::EndMainMenuBar();
  }
  return change;
}
}  // namespace gui
