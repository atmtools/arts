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
  
  // Main loop
  BeginWhileLoopGUI;
  
  // Main menu bar
  MainMenu::fullscreen(config, window);
  MainMenu::quitscreen(config, window);
  
  // Full screen plot
  if (Windows::full(window, Windows::CurrentPosition(), "Plot Window")) {
    if (ImPlot::BeginPlot("Plot Frame", x.c_str(), y.c_str(), {-1, -1})) {
      ImPlot::PlotLine("Line", ydata.get_c_array(), int(ydata.nelem()));
      ImPlot::EndPlot();
    }
  }
  Windows::end();
  
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
  
  // Is this valid?
  const bool valid = xdata.nelem() == ydata.nelem();
  
  // Main loop
  BeginWhileLoopGUI;
  
  // Main menu bar
  MainMenu::fullscreen(config, window);
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
  
  // End of main loop
  EndWhileLoopGUI;
  
  // Clean Graphics data
  CleanupGUI;
}
}
