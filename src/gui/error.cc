#include "error.h"
#include "imgui.h"

namespace ARTSGUI {
ErrorStatus error(const std::string& errmsg) {
  ImGui::OpenPopup("Error");

  if (ImGui::BeginPopupModal("Error")) {
    ImGui::Text("\tThe following internal error message was caught:\t");
    ImGui::Text("\t----------------------------ERROR-----------------------------\t");
    ImGui::Text("%s", errmsg.c_str());
    ImGui::Text("\t----------------------------ERROR-----------------------------\t");
    ImGui::Text("\tNote that this error might be fatal.  If you believe you can\t");
    ImGui::Text("\tnot recover, press Exit.  Otherwise, OK should leave you in\t");
    ImGui::Text("\tfull control of the GUI to fix the error\t");
    if (ImGui::Button("\tOK\t", {-1.0f, 30.0f})) {
      ImGui::CloseCurrentPopup();
      return ErrorStatus::Continue;
    }
    if (ImGui::Button("\tExit\t", {-1.0f, 30.0f})) {
      ImGui::CloseCurrentPopup();
      return ErrorStatus::Exit;
    }
  }

  return ErrorStatus::OnHold;
}
}  // namespace ARTSGUI