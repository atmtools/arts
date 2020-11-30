#ifndef gui_h
#define gui_h

#include <imgui.h>

// IMGUI header must be first
#include <imfilebrowser.h>
#include <imgui_stdlib.h>
#include <implot.h>
#include <implot_internal.h>

#include <string>
#include <vector>

#include "gui_macros.h"

namespace ARTSGUI {
  /** A global config for all things ARTSGUI */
struct Config {
  /** The io variable of ImGui, has a lot of key presses */
  ImGuiIO &io;

  /** Should we go fullscreen? */
  bool fullscreen;

  /** Are we dealing with the error? */
  int active_errors;

  /** Window positions; used when fullscreen is toggled on and off to return
   * window to good position and size */
  int width, height, xpos, ypos;

  /** List of tabs */
  size_t tabspos;
  std::vector<std::string> tabs;

  /** User input */
  bool new_save_path;

  Config(bool fullscreen_on = false)
      : io(ImGui::GetIO()),
        fullscreen(fullscreen_on),
        active_errors(0),
        width(1280),
        height(720),
        xpos(50),
        ypos(50),
        tabspos(0),
        tabs(0),
        new_save_path(false) {
    (void)io;
    io.ConfigFlags |=
        ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
    io.IniFilename = NULL;  // Don't save any initialization file
  }
};  // Config

void LayoutAndStyleSettings();

GLFWmonitor *get_current_monitor(GLFWwindow *window);

namespace MainMenu {
void fullscreen(Config &cfg, GLFWwindow *window);
void quitscreen(const Config &cfg, GLFWwindow *window);
}  // namespace MainMenu

namespace Windows {
template <unsigned WIDTH = 1, unsigned HEIGHT = 1, unsigned WIDTH_POS = 0,
          unsigned HEIGHT_POS = 0, unsigned WIDTH_EXTENT = 1,
          unsigned HEIGHT_EXTENT = 1>
bool sub(GLFWwindow *window, const ImVec2 origpos, const char *name) {
  static_assert(WIDTH, "None size window not allowed");
  static_assert(HEIGHT, "None size window not allowed");
  static_assert(WIDTH_EXTENT, "None extent window not allowed");
  static_assert(HEIGHT_EXTENT, "None extent window not allowed");
  static_assert(WIDTH > (WIDTH_POS + WIDTH_EXTENT - 1),
                "More width than possible");
  static_assert(HEIGHT > (HEIGHT_POS + HEIGHT_EXTENT - 1),
                "More height than possible");
  constexpr float wscale = 1.0f / float(WIDTH);
  constexpr float hscale = 1.0f / float(HEIGHT);

  // Cursors and sizes
  int width = 0, height = 0;
  glfwGetWindowSize(window, &width, &height);
  ImVec2 size = {(float(width) - origpos.x - 5.0f) * wscale,
                 (float(height) - 2 * origpos.y) * hscale};
  ImVec2 pos = {2.5f + origpos.x + size.x * WIDTH_POS,
                1.5f * origpos.y + size.y * HEIGHT_POS};
  size.x *= WIDTH_EXTENT;
  size.y *= HEIGHT_EXTENT;

  // Set a simple window frame
  ImGui::SetNextWindowPos(pos, 1);
  ImGui::SetNextWindowSize(size, 1);

  return ImGui::Begin(name, NULL,
                      ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                          ImGuiWindowFlags_NoSavedSettings |
                          ImGuiWindowFlags_NoResize |
                          ImGuiWindowFlags_NoScrollbar);
}

bool full(GLFWwindow *window, const ImVec2 origpos, const char *name);
ImVec2 CurrentPosition();
void end();
}  // namespace Windows
}  // namespace GUI

#endif  // gui_h
