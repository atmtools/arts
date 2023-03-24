#include "menu.h"

#include <debug.h>
#include <imgui.h>
#include <logic.h>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

#include "atm.h"
#include "jacobian.h"
#include "math_funcs.h"

namespace ARTSGUI::MainMenu {
void fullscreen(Config& cfg, GLFWwindow* window) {
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem(" Fullscreen ", "F11")) {
        if (not cfg.fullscreen) {
          glfwGetWindowSize(window, &cfg.width, &cfg.height);
          glfwGetWindowPos(window, &cfg.xpos, &cfg.ypos);
          auto* monitor = get_current_monitor(window);
          const auto* mode = glfwGetVideoMode(monitor);
          glfwSetWindowMonitor(
              window, monitor, 0, 0, mode->width, mode->height, 0);
        } else
          glfwSetWindowMonitor(
              window, nullptr, cfg.xpos, cfg.ypos, cfg.width, cfg.height, 0);

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
      auto* monitor = get_current_monitor(window);
      const auto* mode = glfwGetVideoMode(monitor);
      glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, 0);
    } else
      glfwSetWindowMonitor(
          window, nullptr, cfg.xpos, cfg.ypos, cfg.width, cfg.height, 0);
    cfg.fullscreen = not cfg.fullscreen;
  }
}

void quitscreen(const Config& cfg, GLFWwindow* window) {
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

bool exportdata(const Config& cfg,
                ImGui::FileBrowser& fileBrowser,
                const char* dialog,
                bool shortcut) {
  bool pressed = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("File")) {
      if (ImGui::MenuItem(dialog, shortcut ? "Ctrl+S" : "")) {
        fileBrowser.Open();
        pressed = true;
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  if (shortcut) {
    if (cfg.io.KeyCtrl and ImGui::IsKeyPressed(GLFW_KEY_S)) {
      fileBrowser.Open();
    }
  }

  return pressed;
}

bool change_item(const char* name) {
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      if (ImGui::BeginMenu(name, false)) {
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return false;
}

[[nodiscard]] std::string change_item_name(const JacobianTarget& target) {
  switch (target.type) {
    case Jacobian::Type::Atm:
      switch (target.atm) {
        case Jacobian::Atm::Temperature:
          return "Temperature";
          break;
        case Jacobian::Atm::WindMagnitude:
          return "Wind Magnitude";
          break;
        default:
          ARTS_USER_ERROR("Not implemented")
      }

      break;
    case Jacobian::Type::Line:
      switch (target.line) {
        case Jacobian::Line::VMR:
          return var_string("VMR: ", target.qid.Isotopologue());
          break;
        default:
          ARTS_USER_ERROR("Not implemented")
      }
      break;
    default:
      ARTS_USER_ERROR("Not implemented")
  }

  return "Bad Value";
}

bool change_item(const char* name,
                 ArrayOfRetrievalQuantity& jac,
                 ArrayOfRetrievalQuantity& old) {
  bool did_something = false;
  Jacobian::Target target;
  std::string target_name;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      auto tpred = [](auto& j) { return j == Jacobian::Atm::Temperature; };
      auto fpred = [](auto& j) { return j == Jacobian::Atm::WindMagnitude; };

      if (ImGui::BeginMenu(name)) {
        target = Jacobian::Target(Jacobian::Atm::Temperature);
        target_name = '\t' + change_item_name(target) + '\t';
        if (bool has = std::any_of(jac.begin(), jac.end(), tpred);
            ImGui::Selectable(target_name.c_str(),
                              has,
                              ImGuiSelectableFlags_DontClosePopups)) {
          if (has) {
            std::remove_if(jac.begin(), jac.end(), tpred);
            jac.pop_back();
          } else {
            target.perturbation = 0.1;
            jac.emplace_back().Target() = target;
          }
          did_something = true;
        }
        ImGui::Separator();

        target = Jacobian::Target(Jacobian::Atm::WindMagnitude);
        target_name = '\t' + change_item_name(target) + '\t';
        if (bool has = std::any_of(jac.begin(), jac.end(), fpred);
            ImGui::Selectable(target_name.c_str(),
                              has,
                              ImGuiSelectableFlags_DontClosePopups)) {
          if (has) {
            std::remove_if(jac.begin(), jac.end(), fpred);
            jac.pop_back();
          } else {
            target.perturbation = 100;
            jac.emplace_back().Target() = target;
          }
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tRestore Value\t", {-1, 0})) {
          jac = old;
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tStore Value\t", {-1, 0})) old = jac;
        ImGui::Separator();
        ImGui::EndMenu();
      }

      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return did_something;
}

bool change_item(const char* name,
                 Vector& vec,
                 Vector& old,
                 const ArrayOfString& keys) {
  const Index n = vec.nelem();
  ARTS_ASSERT(n == keys.nelem())
  bool did_something = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      if (ImGui::BeginMenu(name)) {
        for (Index i = 0; i < n; i++) {
          ImGui::Text("\t");
          ImGui::SameLine();
          if (ImGui::InputDouble(keys[i].c_str(), &vec[i], 0, 0, "%g"))
            did_something = true;
          ImGui::Separator();
        }
        ImGui::Separator();
        if (ImGui::Button("\tRestore Value\t", {-1, 0})) {
          vec = old;
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tStore Value\t", {-1, 0})) old = vec;
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return did_something;
}



bool change_item(const char* name,
                 std::array<Numeric, 3>& vec,
                 std::array<Numeric, 3>& old,
                 const std::array<String, 3>& keys) {
  constexpr Index n = 3;
  bool did_something = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      if (ImGui::BeginMenu(name)) {
        for (Index i = 0; i < n; i++) {
          ImGui::Text("\t");
          ImGui::SameLine();
          if (ImGui::InputDouble(keys[i].c_str(), &vec[i], 0, 0, "%g"))
            did_something = true;
          ImGui::Separator();
        }
        ImGui::Separator();
        if (ImGui::Button("\tRestore Value\t", {-1, 0})) {
          vec = old;
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tStore Value\t", {-1, 0})) old = vec;
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return did_something;
}

bool change_item(const char* name,
                 std::array<Numeric, 3>& vec,
                 std::array<Numeric, 3>& old,
                 const ArrayOfString& keys) {
  const Index n = 3;
  ARTS_ASSERT(n == keys.nelem())
  bool did_something = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      if (ImGui::BeginMenu(name)) {
        for (Index i = 0; i < n; i++) {
          ImGui::Text("\t");
          ImGui::SameLine();
          if (ImGui::InputDouble(keys[i].c_str(), &vec[i], 0, 0, "%g"))
            did_something = true;
          ImGui::Separator();
        }
        ImGui::Separator();
        if (ImGui::Button("\tRestore Value\t", {-1, 0})) {
          vec = old;
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tStore Value\t", {-1, 0})) old = vec;
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return did_something;
}

bool change_item(const char* name,
                 AtmPoint& vec,
                 AtmPoint& old,
                 const ArrayOfArrayOfSpeciesTag& spec,
                 Options& menu) {
  ARTS_ASSERT(vec.nelem() == spec.nelem())
  bool did_something = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      if (ImGui::BeginMenu(name)) {
        if (ImGui::BeginMenu("\tSelect VMR type\t")) {
          for (auto& x : enumtyps::VMRTypes) {
            if (ImGui::Selectable(
                    (' ' + std::string{toString(x)} + ' ').c_str(),
                    x == menu.vmr,
                    ImGuiSelectableFlags_DontClosePopups))
              menu.vmr = x;
          }
          ImGui::Separator();
          ImGui::EndMenu();
        }
        ImGui::Separator();

        const std::string vmr_type{toString(menu.vmr)};
        constexpr Numeric max = 1.5;
        const Numeric scale = [](VMR& vmr) {
          switch (vmr) {
            case VMR::exact:
              return 1.0;
            case VMR::percent:
              return 100.0;
            case VMR::ppmv:
              return 1'000'000.0;
            case VMR::FINAL: { /* leave last */
            }
          }
          return 0.0;
        }(menu.vmr);
        ImGui::Text("\tVMR (range: [0, %g]; type: %s):\t",
                    scale * max,
                    vmr_type.c_str());

        for (Index i = 0; i < spec.nelem(); i++) {
          const std::string spec_name{var_string('\t', spec[i], '\t')};
          Numeric val = scale * vec[spec[i]];
          ImGui::Text("\t");
          ImGui::SameLine();
          if (ImGui::InputDouble(spec_name.c_str(), &val, 0, 0, "%g")) {
            val /= scale;
            vec[spec[i]] = std::clamp(val, 0.0, max);
            did_something = true;
          }
          ImGui::Separator();
        }

        ImGui::Separator();
        if (ImGui::Button("\tRestore Value\t", {-1, 0})) {
          vec = old;
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tStore Value\t", {-1, 0})) old = vec;
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return did_something;
}

bool change_item(
    const char* name, Vector& vec, Vector& old, Numeric min, Numeric max) {
  Index n = vec.nelem();
  ARTS_ASSERT(is_sorted(vec))
  ARTS_ASSERT(min < max)
  ARTS_ASSERT(n > 1)
  bool did_something = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      if (ImGui::BeginMenu(name)) {
        Numeric start = vec[0];
        Numeric stop = vec[n - 1];

        bool change = false;
        ImGui::Text("\t");
        ImGui::SameLine();
        if (ImGui::InputDouble("\tStart\t", &start, 0, 0, "%g")) {
          start = std::clamp(start, min, stop);
          change = true;
        }
        ImGui::Text("\t");
        ImGui::SameLine();
        if (ImGui::InputDouble("\tStop\t", &stop, 0, 0, "%g")) {
          stop = std::clamp(stop, start, max);
          change = true;
        }
        ImGui::Text("\t");
        ImGui::SameLine();
        if (ImGui::InputScalar("\tnelem\t", ImGuiDataType_S64, &n)) {
          n = std::clamp(n, Index{2}, std::numeric_limits<Index>::max());
          change = true;
        }

        if (change) {
          nlinspace(vec, start, stop, n);
          did_something = true;
        }

        ImGui::Separator();
        if (ImGui::Button("\tRestore Value\t", {-1, 0})) {
          vec = old;
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tStore Value\t", {-1, 0})) old = vec;
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return did_something;
}

bool change_item(const char* name,
                 const char* value_name,
                 Numeric& val,
                 Numeric& old,
                 Numeric min,
                 Numeric max) {
  bool did_something = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      if (ImGui::BeginMenu(name)) {
        Numeric v = val;
        ImGui::Text("\t");
        ImGui::SameLine();
        if (ImGui::InputDouble(value_name, &v, 0, 0, "%g")) {
          val = std::clamp(v, min, max);
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tRestore Value\t", {-1, 0})) {
          val = old;
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tStore Value\t", {-1, 0})) old = val;
        ImGui::Separator();
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return did_something;
}

bool change_item(const char* name,
                 ArrayOfSpeciesTag& out,
                 ArrayOfSpeciesTag& old,
                 const ArrayOfArrayOfSpeciesTag& keys) {
  bool did_something = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Value")) {
      if (ImGui::BeginMenu(name)) {
        if (ImGui::Selectable(" *All* ",
                              out.nelem() == 0,
                              ImGuiSelectableFlags_DontClosePopups)) {
          out.resize(0);
          did_something = true;
        }
        for (auto& key : keys) {
          ImGui::Separator();
          const std::string str{var_string(' ', key, ' ')};
          if (ImGui::Selectable(str.c_str(),
                                key == out,
                                ImGuiSelectableFlags_DontClosePopups)) {
            out = key;
            did_something = true;
          }
        }
        ImGui::Separator();
        if (ImGui::Button("\tRestore Value\t", {-1, 0})) {
          out = old;
          did_something = true;
        }
        ImGui::Separator();
        if (ImGui::Button("\tStore Value\t", {-1, 0})) old = out;
        ImGui::EndMenu();
      }
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  return did_something;
}

std::string absunit(const Jacobian::Target& target) {
  switch (target.type) {
    case Jacobian::Type::Atm:
      switch (target.atm) {
        case Jacobian::Atm::Temperature:
          return "Absorption Partial Derivative [1/m K]";
          break;
        case Jacobian::Atm::WindMagnitude:
          return "Absorption Partial Derivative [1/m (m/s)]";
          break;
        default:
          ARTS_USER_ERROR("Not implemented")
      }

      break;
    case Jacobian::Type::Line:
      switch (target.line) {
        case Jacobian::Line::VMR:
          return var_string(
              "Absorption Partial Derivative [1/m ", target.qid, ']');
          break;
        default:
          ARTS_USER_ERROR("Not implemented")
      }
      break;
    default:
      ARTS_USER_ERROR("Not implemented")
  }
}

void select_option(Index& ind, const ArrayOfRetrievalQuantity& jac) {
  if (ImGui::Selectable("\tMain Calculations\t",
                        ind == -1,
                        ImGuiSelectableFlags_DontClosePopups)) {
    ind = -1;
  }

  for (Index i = 0; i < jac.nelem(); i++) {
    const std::string opt{var_string('\t', "Derivative: ", change_item_name(jac[i].Target()), '\t')};
    if (ImGui::Selectable(
            opt.c_str(), ind == i, ImGuiSelectableFlags_DontClosePopups)) {
      ind = i;
    }
  }
}

void tooltip(const char* tip, const Config& config) {
  if (ImGui::IsItemHovered()) {
    if (ImGui::GetCurrentContext()->HoveredIdTimer > config.hover_time_limit) {
      ImGui::SetTooltip(" \n\t%s\t\n ", tip);
      ImGui::BeginTooltip();
      ImGui::EndTooltip();
    }
  }
}

void show_plot_controls() {
  static bool doit = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Help")) {
      if (ImGui::MenuItem("Show Plot Controls")) doit = true;
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  if (doit) {
    ImGui::OpenPopup("Plot Controls");
    doit = false;
  }

  if (ImGui::BeginPopupModal("Plot Controls")) {
    ImGui::Text(R"--(    The following controls are available to manipulate the plot panel                
    Note that gui has to be selected for any interactions to work

    Inside the plot panel:
        Left-click, hold, and drag:    Moves the current axes
        Double left-click:             Recenter the plot axis (might miss very small values)
        Right-click, hold, and drag:   Zoom in on the highlighted area
                     Shift modifier:   Selects the full y-axis
                       Alt modifier:   Selects the full x-axis
        Double right-click:            Open a menu that allows modifying the plot layout
        Spin mouse scroll-wheel:       Zoom in and out
        Press the color of the legend: Turn on-and-off the plot line (useful mostly when stokes_dim>1)    

    On the x-axis and on the y-axis:
        Same as above but for just the axis, bar the right-click, hold, and drag features
    
    The double right-click menu:
        X-Axis/Y-Axis:
            Set min value manually, and/or lock the minimum value by ticking the box
            Set max value manually, and/or lock the maximum value by ticking the box
            Invert:    The order of the axis is inverted (min at top) if ticked
            Log Scale: The axis is shown logarithmically if ticked
            Time:      Interpret the X-Axis as a unix time stamp (do not use)
        Other settings: Various interaction best understood by using them
)--");
    if (ImGui::Button("\tOK\t", {-1.0f, 30.0f})) {
      ImGui::CloseCurrentPopup();
    }
  }
}

void show_propmat_controls() {
  static bool doit = false;

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Help")) {
      if (ImGui::MenuItem("Show Propagation Matrix Agenda Controls")) doit = true;
      ImGui::Separator();
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  if (doit) {
    ImGui::OpenPopup("Propagation Matrix Agenda Controls");
    doit = false;
  }

  if (ImGui::BeginPopupModal("Propagation Matrix Agenda Controls")) {
    ImGui::Text(R"--(    The following controls are available to manipulate agenda                              
    Note that gui has to be selected for any interactions to work

    The Panel Control Panel:
        Run Agenda Button:                 This button runs the agenda with the current
                                           values
        Update Automatic Tickbox:          If ticked, then the agenda will run as soon
                                           as a value is updated    
        Automatic Frequency Grid Tickbox:  If ticked, then the agenda is run when the
                                           drawn x-axis is updated    
        Last Run Information:              Shows the time it took to execute the last
                                           agenda run
        Not Running / Running Information: The first shows if we are running, the latter    
                                           shows for how long the current calculations
                                           have been going on
    
    The Panel Status Panel:
        This shows the input to the last run agenda.  If a star appears next to a value,
        then this value is different from what the current settings entail (indicating
        you might want to run the agenda)
    
    Value Menu:
        Allows manually setting the propmat_clearsky_agenda inputs, as well as the
        transmission distance

    Display Menu:
        Propagation Matrix Scale:     Scale the Y-Axis of the propagation matrix
        Transmission Matrix Scale:    Scale the Y-Axis of the transmission matrix
        Frequency Scale:              Scale the X-Axis to the new frequency unit
        Select Plot:                  Select main calculations or a derivative to show
        Binning Count:                Bin data to make it smoother (number of binned
                                      points)
        Show As Transmission Tickbox: Shows the propagation matrix when unticked and the
        transmission matrix when ticked    
    Note that the rescaling does not change the numbers you are looking at in the plot
    panel but that you will have to refit the plotting window to match the new scale.
    Some scales might even make the numbers so small that the automatic rescaling of
    the plot panel with double left-click misses the values.  Also note that automatic
    frequency grids do not work very well with the inverse frequency scales.
)--");
    if (ImGui::Button("\tOK\t", {-1.0f, 30.0f})) {
      ImGui::CloseCurrentPopup();
    }
  }
}
}  // namespace ARTSGUI::MainMenu
