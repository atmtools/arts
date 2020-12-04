#include "plot.h"

namespace ARTSGUI {
void plot(const Vector& ydata) {
  // Get Graphics data
  InitializeGUI("Plot");
  
  // Our global states are stored in config
  Config config;
  
  // Our style
  LayoutAndStyleSettings();
  
  // Internal states
  std::string x = "X";
  std::string y = "Y";
  auto fileBrowser = ARTSGUI::Files::xmlfile_chooser();
  
  // Main loop
  BeginWhileLoopGUI;
  
  // Main menu bar
  MainMenu::fullscreen(config, window);
  MainMenu::exportdata(config, fileBrowser);
  MainMenu::quitscreen(config, window);
  
  // Full screen plot
  if (Windows::full(window, Windows::CurrentPosition(), "Plot Window")) {
    if (ImPlot::BeginPlot("Plot Frame", x.c_str(), y.c_str(), {-1, -1})) {
      ImPlot::PlotLine("Line", ydata.get_c_array(), int(ydata.nelem()));
      ImPlot::EndPlot();
    }
  }
  Windows::end();
  
  ARTSGUI::Files::save_data(config, fileBrowser, ydata);
  
  // End of main loop
  EndWhileLoopGUI;
  
  // Clean Graphics data
  CleanupGUI;
}


void plot(const Vector& xdata, const Vector& ydata) {
  // Get Graphics data
  InitializeGUI("Plot");
  
  // Our global states are stored in config
  Config config;
  
  // Our style
  LayoutAndStyleSettings();
  
  // Internal states
  std::string x = "X";
  std::string y = "Y";
  auto fileBrowser = ARTSGUI::Files::xmlfile_chooser();
  
  // Is this valid?
  const bool valid = xdata.nelem() == ydata.nelem();
  
  // Main loop
  BeginWhileLoopGUI;
  
  // Main menu bar
  MainMenu::fullscreen(config, window);
  MainMenu::exportdata(config, fileBrowser);
  MainMenu::quitscreen(config, window);
  
  // Full screen plot if valid or just a warning if invalid data
  if (Windows::full(window, Windows::CurrentPosition(), "Plot Window")) {
    if (valid) {
      if (ImPlot::BeginPlot("Plot Frame", x.c_str(), y.c_str(), {-1, -1})) {
        ImPlot::PlotLine("Line", xdata.get_c_array(), ydata.get_c_array(), int(ydata.nelem()));
        ImPlot::EndPlot();
      }
    } else {
      ImGui::Text("Invalid sizes, xdata is %ld elements and ydata is %ld elements", xdata.nelem(), ydata.nelem());
    }
  }
  Windows::end();
  
  ARTSGUI::Files::save_data(config, fileBrowser, ydata);
  
  // End of main loop
  EndWhileLoopGUI;
  
  // Clean Graphics data
  CleanupGUI;
}


bool same_lengths(const ArrayOfVector& xdata, const ArrayOfVector& ydata) {
  for (std::size_t i=0; i<xdata.size(); i++)
    if (xdata[i].nelem() not_eq ydata[i].nelem())
      return false;
  return true;
}


void plot(const ArrayOfVector& xdata, const ArrayOfVector& ydata) {
  if (xdata.size() == 1 and ydata.size() == 1) {
    plot(xdata[0], ydata[0]);
    return;
  }
  
  // Get Graphics data
  InitializeGUI("Plot");
  
  // Our global states are stored in config
  Config config;
  
  // Our style
  LayoutAndStyleSettings();
  
  // Internal states
  std::string x = "X";
  std::string y = "Y";
  std::vector<std::string> lines(xdata.size());
  for (std::size_t i=0; i<lines.size(); i++) {
    lines[i] = std::string("Line ") + std::to_string(i+1);
  }
  auto fileBrowser = ARTSGUI::Files::xmlfile_chooser();
  
  // Is this valid?
  const bool valid = xdata.size() == ydata.size() and same_lengths(xdata, ydata);
  
  // Main loop
  BeginWhileLoopGUI;
  
  // Main menu bar
  MainMenu::fullscreen(config, window);
  MainMenu::exportdata(config, fileBrowser);
  MainMenu::quitscreen(config, window);
  
  // Full screen plot if valid or just a warning if invalid data
  if (Windows::full(window, Windows::CurrentPosition(), "Plot Window")) {
    if (valid) {
      if (ImPlot::BeginPlot("Plot Frame", x.c_str(), y.c_str(), {-1, -1})) {
        for (std::size_t i=0; i<lines.size(); i++) {
          ImPlot::PlotLine(lines[i].c_str(), xdata[i].get_c_array(), ydata[i].get_c_array(), int(ydata[i].nelem()));
        }
        ImPlot::EndPlot();
      }
    } else {
      ImGui::Text("Invalid sizes, xdata is %ld elements and ydata is %ld elements", xdata.size(), ydata.size());
    }
  }
  Windows::end();
  
  ARTSGUI::Files::save_data(config, fileBrowser, ydata);
  
  // End of main loop
  EndWhileLoopGUI;
  
  // Clean Graphics data
  CleanupGUI;
}
}  // ARTSGUI
