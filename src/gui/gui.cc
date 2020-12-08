#include "gui.h"
#include <xml_io.h>

namespace ARTSGUI {
  void LayoutAndStyleSettings() {
    auto &style = ImGui::GetStyle();
    style.FramePadding = {0.0f, 0.0f};
    style.FrameRounding = 0.0f;
    style.FrameBorderSize = 0.0f;
    style.WindowPadding = {0.0f, 0.0f};
    style.WindowRounding = 0.0f;
    style.WindowBorderSize = 0.1f;
    style.TabRounding = 0.0f;
    style.TabBorderSize = 1.0f;
    style.AntiAliasedLinesUseTex = true;
    style.AntiAliasedFill = true;
    style.ChildRounding = 0.0f;
    
    style.Colors[ImGuiCol_Text] = ImVec4(0.92f, 0.92f, 0.92f, 1.00f);
    style.Colors[ImGuiCol_TextDisabled] = ImVec4(0.44f, 0.44f, 0.44f, 1.00f);
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.06f, 0.06f, 0.06f, 1.00f);
    style.Colors[ImGuiCol_ChildBg] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    style.Colors[ImGuiCol_PopupBg] = ImVec4(0.08f, 0.08f, 0.08f, 0.94f);
    style.Colors[ImGuiCol_Border] = ImVec4(0.51f, 0.36f, 0.15f, 1.00f);
    style.Colors[ImGuiCol_BorderShadow] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    style.Colors[ImGuiCol_FrameBg] = ImVec4(0.11f, 0.11f, 0.11f, 1.00f);
    style.Colors[ImGuiCol_FrameBgHovered] = ImVec4(0.51f, 0.36f, 0.15f, 1.00f);
    style.Colors[ImGuiCol_FrameBgActive] = ImVec4(0.78f, 0.55f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_TitleBg] = ImVec4(0.51f, 0.36f, 0.15f, 1.00f);
    style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.91f, 0.64f, 0.13f, 1.00f);
    style.Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.00f, 0.00f, 0.00f, 0.51f);
    style.Colors[ImGuiCol_MenuBarBg] = ImVec4(0.11f, 0.11f, 0.11f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarBg] = ImVec4(0.06f, 0.06f, 0.06f, 0.53f);
    style.Colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.21f, 0.21f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarGrabHovered] =
    ImVec4(0.47f, 0.47f, 0.47f, 1.00f);
    style.Colors[ImGuiCol_ScrollbarGrabActive] =
    ImVec4(0.81f, 0.83f, 0.81f, 1.00f);
    style.Colors[ImGuiCol_CheckMark] = ImVec4(0.78f, 0.55f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.91f, 0.64f, 0.13f, 1.00f);
    style.Colors[ImGuiCol_SliderGrabActive] = ImVec4(0.91f, 0.64f, 0.13f, 1.00f);
    style.Colors[ImGuiCol_Button] = ImVec4(0.51f, 0.36f, 0.15f, 1.00f);
    style.Colors[ImGuiCol_ButtonHovered] = ImVec4(0.91f, 0.64f, 0.13f, 1.00f);
    style.Colors[ImGuiCol_ButtonActive] = ImVec4(0.78f, 0.55f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_Header] = ImVec4(0.51f, 0.36f, 0.15f, 1.00f);
    style.Colors[ImGuiCol_HeaderHovered] = ImVec4(0.91f, 0.64f, 0.13f, 1.00f);
    style.Colors[ImGuiCol_HeaderActive] = ImVec4(0.93f, 0.65f, 0.14f, 1.00f);
    style.Colors[ImGuiCol_Separator] = ImVec4(0.21f, 0.21f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_SeparatorHovered] = ImVec4(0.91f, 0.64f, 0.13f, 1.00f);
    style.Colors[ImGuiCol_SeparatorActive] = ImVec4(0.78f, 0.55f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_ResizeGrip] = ImVec4(0.21f, 0.21f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_ResizeGripHovered] = ImVec4(0.91f, 0.64f, 0.13f, 1.00f);
    style.Colors[ImGuiCol_ResizeGripActive] = ImVec4(0.78f, 0.55f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_Tab] = ImVec4(0.51f, 0.36f, 0.15f, 1.00f);
    style.Colors[ImGuiCol_TabHovered] = ImVec4(0.91f, 0.64f, 0.13f, 1.00f);
    style.Colors[ImGuiCol_TabActive] = ImVec4(0.78f, 0.55f, 0.21f, 1.00f);
    style.Colors[ImGuiCol_TabUnfocused] = ImVec4(0.07f, 0.10f, 0.15f, 0.97f);
    style.Colors[ImGuiCol_TabUnfocusedActive] =
    ImVec4(0.14f, 0.26f, 0.42f, 1.00f);
    style.Colors[ImGuiCol_PlotLines] = ImVec4(0.61f, 0.61f, 0.61f, 1.00f);
    style.Colors[ImGuiCol_PlotLinesHovered] = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
    style.Colors[ImGuiCol_PlotHistogram] = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
    style.Colors[ImGuiCol_PlotHistogramHovered] =
    ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
    style.Colors[ImGuiCol_TextSelectedBg] = ImVec4(0.26f, 0.59f, 0.98f, 0.35f);
    style.Colors[ImGuiCol_DragDropTarget] = ImVec4(1.00f, 1.00f, 0.00f, 0.90f);
    style.Colors[ImGuiCol_NavHighlight] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
    style.Colors[ImGuiCol_NavWindowingHighlight] =
    ImVec4(1.00f, 1.00f, 1.00f, 0.70f);
    style.Colors[ImGuiCol_NavWindowingDimBg] = ImVec4(0.80f, 0.80f, 0.80f, 0.20f);
    style.Colors[ImGuiCol_ModalWindowDimBg] = ImVec4(0.80f, 0.80f, 0.80f, 0.35f);
    
    auto &plot_style = ImPlot::GetStyle();
    plot_style.UseLocalTime = true;
    plot_style.UseISO8601 = true;
    plot_style.Use24HourClock = true;
}

GLFWmonitor *get_current_monitor(GLFWwindow *window) {
  int nmonitors, i;
  int wx, wy, ww, wh;
  int mx, my;
  int bestoverlap;
  GLFWmonitor *bestmonitor;
  GLFWmonitor **monitors;
  const GLFWvidmode *mode;
  auto mini = [](int x, int y){ return x < y ? x : y; };
  auto maxi = [](int x, int y){ return x > y ? x : y; };
  
  bestoverlap = 0;
  bestmonitor = NULL;
  
  glfwGetWindowPos(window, &wx, &wy);
  glfwGetWindowSize(window, &ww, &wh);
  monitors = glfwGetMonitors(&nmonitors);
  
  for (i = 0; i < nmonitors; i++) {
    mode = glfwGetVideoMode(monitors[i]);
    glfwGetMonitorPos(monitors[i], &mx, &my);
    
    int mw = mode->width;
    int mh = mode->height;
    
    int overlap = maxi(0, mini(wx + ww, mx + mw) - maxi(wx, mx)) *
    maxi(0, mini(wy + wh, my + mh) - maxi(wy, my));
    
    if (bestoverlap < overlap) {
      bestoverlap = overlap;
      bestmonitor = monitors[i];
    }
  }
  
  return bestmonitor;
}

namespace MainMenu {
void fullscreen(Config &cfg, GLFWwindow *window) {
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem(" Fullscreen ", "F11")) {
        if (not cfg.fullscreen) {
          glfwGetWindowSize(window, &cfg.width, &cfg.height);
          glfwGetWindowPos(window, &cfg.xpos, &cfg.ypos);
          auto *monitor = get_current_monitor(window);
          const auto *mode = glfwGetVideoMode(monitor);
          glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height,
                               0);
        } else
          glfwSetWindowMonitor(window, NULL, cfg.xpos, cfg.ypos, cfg.width,
                               cfg.height, 0);

        cfg.fullscreen = not cfg.fullscreen;
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  if (ImGui::IsKeyPressed(GLFW_KEY_F11) or
      (cfg.fullscreen and ImGui::IsKeyPressed(GLFW_KEY_ESCAPE))) {
    if (not cfg.fullscreen) {
      glfwGetWindowSize(window, &cfg.width, &cfg.height);
      glfwGetWindowPos(window, &cfg.xpos, &cfg.ypos);
      glfwGetWindowPos(window, &cfg.xpos, &cfg.ypos);
      auto *monitor = get_current_monitor(window);
      const auto *mode = glfwGetVideoMode(monitor);
      glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, 0);
    } else
      glfwSetWindowMonitor(window, NULL, cfg.xpos, cfg.ypos, cfg.width,
                           cfg.height, 0);
    cfg.fullscreen = not cfg.fullscreen;
  }
}

void quitscreen(const Config &cfg, GLFWwindow *window) {
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem(" Quit ", "Ctrl+X"))
        glfwSetWindowShouldClose(window, 1);
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  if (cfg.io.KeyCtrl and ImGui::IsKeyPressed(GLFW_KEY_X)) {
    glfwSetWindowShouldClose(window, 1);
  }
}

void exportdata(const Config &cfg, ImGui::FileBrowser& fileBrowser) {
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem(" Export Data ", "Ctrl+S")) {
        fileBrowser.Open();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }
  
  if (cfg.io.KeyCtrl and ImGui::IsKeyPressed(GLFW_KEY_S)) {
    fileBrowser.Open();
  }
}
}  // MainMenu

namespace Windows {
bool full(GLFWwindow *window, const ImVec2 origpos, const char *name) {
  return sub(window, origpos, name);
}

ImVec2 CurrentPosition() { return ImGui::GetCursorPos(); }
  
void end() { ImGui::End(); }
}  // Windows

namespace Files {
  ImGui::FileBrowser xmlfile_chooser() {
    ImGui::FileBrowser fileBrowser(ImGuiFileBrowserFlags_EnterNewFilename | ImGuiFileBrowserFlags_CloseOnEsc);
    fileBrowser.SetTitle("Choose Save File");
    fileBrowser.SetPwd();
    fileBrowser.SetTypeFilters({".xml", ".xml.bin", ".xml.gz"});
    return fileBrowser;
  }
  
  void fix_file_ext_xml(Config& config) {
    if(config.save_path.extension() == ".xml") {
    } else {
      config.save_path  += ".xml";
    }
  }
  
  template <class X>
  void save_data_impl(Config& config, ImGui::FileBrowser& fileBrowser, const X& data) {
    fileBrowser.Display();
    if (fileBrowser.HasSelected()) {
      config.save_path = fileBrowser.GetSelected();
      
      // Fix filename
      fix_file_ext_xml(config);
      
      // Open popup box to not overwrite files without consent
      if (not config.new_save_path and std::filesystem::exists(config.save_path)) {
        ImGui::OpenPopup("Overwrite?");
      }
      
      // Still, we are here so we have a new file to save
      config.new_save_path = true;
      
      // Now, if we have opened a popup, that means we had the file so we need to be careful
      bool save = false;
      if (ImGui::BeginPopupModal("Overwrite?")) {
        // Tell the user about the problem
        ImGui::Text("\n\t %s exists\t \n\t Overwrite?\t ", config.save_path .c_str());
        
        // Press OK and we save the file
        if (ImGui::Button(" OK ", {80.0f, 30.0f})) {
          ImGui::CloseCurrentPopup();
          save = true;
        }
        ImGui::SameLine();
        
        // If you cancel, we clear the selection already and won't save
        if (ImGui::Button(" Cancel ", {80.0f, 30.0f})) {
          config.new_save_path = false;
          fileBrowser.ClearSelected();
          ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
      } else {
        save = true;
      }
      
      // Without popup, we save but otherwise we skip saving and wait
      if (save) {
        xml_write_to_file(config.save_path.native(), data, FILE_TYPE_ASCII, 0, {});
        config.new_save_path = false;
        fileBrowser.ClearSelected();
      }
    }
  }
  
  void save_data(Config& config, ImGui::FileBrowser& fileBrowser, const ArrayOfVector& data) {
    save_data_impl(config, fileBrowser, data);
  }
  
  void save_data(Config& config, ImGui::FileBrowser& fileBrowser, const Vector& data) {
    save_data_impl(config, fileBrowser, data);
  }
}  // Files
}  // ARTSGUI
