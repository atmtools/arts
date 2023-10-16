#include "gui.h"
#include <xml_io.h>
#include "imgui.h"
#include "implot.h"
#include "propagationmatrix.h"

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

    const ImVec4 arts_blue{2.f/255.f, 113.f/255.f, 187.f/255.f, 1.0f};
    const ImVec4 arts_black{0.f, 0.f, 0.f, 1.0f};
    const ImVec4 arts_grey{0.44f, 0.44f, 0.44f, 1.0f};
    const ImVec4 arts_white{1.f, 1.f, 1.f, 1.0f};
    
    style.Colors[ImGuiCol_Text] = arts_white;
    style.Colors[ImGuiCol_TextDisabled] = arts_grey;
    style.Colors[ImGuiCol_WindowBg] = arts_black;
    style.Colors[ImGuiCol_Border] = arts_blue;
    
    auto &plot_style = ImPlot::GetStyle();
    plot_style.UseLocalTime = true;
    plot_style.UseISO8601 = true;
    plot_style.Use24HourClock = true;

    // Line colors ???
    const ImVec4 arts_yellow{235.f/255.f, 172.f/255.f, 35.f/255.f, 1.0f};
    const ImVec4 arts_red{184.f/255.f, 0.f/255.f, 88.f/255.f, 1.0f};
    const ImVec4 arts_green{0.f/255.f, 110.f/255.f, 0.f/255.f, 1.0f};
    const ImVec4 arts_pink{255.f/255.f, 146.f/255.f, 135.f/255.f, 1.0f};
    const ImVec4 arts_brown{178.f/255.f, 69.f/255.f, 2.f/255.f, 1.0f};
    const ImVec4 arts_teal{0.f/255.f, 187.f/255.f, 173.f/255.f, 1.0f};
    static const std::array<ImVec4, 7> arts_line_colors{arts_blue,
                                                        arts_yellow,
                                                        arts_red,
                                                        arts_green,
                                                        arts_pink,
                                                        arts_brown,
                                                        arts_teal};
    ImPlot::SetColormap(arts_line_colors.data(), arts_line_colors.size());
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
  
  void save_data(Config& config, ImGui::FileBrowser& fileBrowser, const PropagationMatrix& data) {
    save_data_impl(config, fileBrowser, data);
  }
  
  void save_data(Config& config, ImGui::FileBrowser& fileBrowser, const ArrayOfPropagationMatrix& data) {
    save_data_impl(config, fileBrowser, data);
  }
}  // Files
}  // ARTSGUI
