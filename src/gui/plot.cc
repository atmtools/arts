#include <cinttypes>
#include <mutex>

#include "plot.h"
#include "menu.h"

namespace ARTSGUI {
// Defaults
std::string PlotConfig::Frame = "Plot";
std::string PlotConfig::X = "X";
std::string PlotConfig::Y = "Y";
std::string PlotConfig::Title = "Plot Frame";
std::vector<std::string> PlotConfig::Legend = {};
  
bool same_lengths(const ArrayOfVector& xdata, const ArrayOfVector& ydata) {
  for (std::size_t i=0; i<xdata.size(); i++)
    if (xdata[i].nelem() not_eq ydata[i].nelem())
      return false;
  return true;
}

std::mutex mtx;

void plot(const ArrayOfVector& xdata, const ArrayOfVector& ydata) {
  auto one_at_a_time = std::lock_guard(mtx);
  
  // Get Graphics data
  InitializeGUI(PlotConfig::Frame.c_str(), 1280, 720);
  
  // Our global states are stored in config
  Config config{};
  
  // Our style
  LayoutAndStyleSettings();
  
  // Internal states
  std::string x = PlotConfig::X;
  std::string y = PlotConfig::Y;
  std::vector<std::string> lines(xdata.size());
  if (lines.size() == PlotConfig::Legend.size()) {
    lines = PlotConfig::Legend;
  } else {
    if (lines.size() > 1) {
      for (std::size_t i=0; i<lines.size(); i++) {
        lines[i] = std::string("Line ") + std::to_string(i+1);
      }
    } else if (lines.size() == 1) {
      lines[0] = std::string("Line");
    }
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
      if (ImPlot::BeginPlot(PlotConfig::Title.c_str(), x.c_str(), y.c_str(), {-1, -1})) {
        for (std::size_t i=0; i<lines.size(); i++) {
          ImPlot::PlotLine(lines[i].c_str(), xdata[i].get_c_array(), ydata[i].get_c_array(), int(ydata[i].nelem()));
        }
        ImPlot::EndPlot();
      }
    } else {
      if (xdata.size() not_eq ydata.size())
        ImGui::Text("Invalid sizes, xdata is %ld elements and ydata is %ld elements", xdata.size(), ydata.size());
      else {
        for (Index i = 0; i < xdata.nelem(); i++) {
          ImGui::Text("xdata[%" PRId64 "] is %" PRId64 " elements and ydata[%" PRId64
                      "] is %" PRId64 " elements",
                      i,
                      xdata[i].size(),
                      i,
                      ydata[i].size());
        }
      }
    }
  }
  Windows::end();
  
  // Save the data to file?
  if (lines.size() == 1) {
    ARTSGUI::Files::save_data(config, fileBrowser, ydata[0]);
  } else {
    ARTSGUI::Files::save_data(config, fileBrowser, ydata);
  }
  
  // End of main loop
  EndWhileLoopGUI;
  
  // Clean Graphics data
  CleanupGUI;
}

void plot(const Vector& xdata, const ArrayOfVector& y) {
  const ArrayOfVector x(y.nelem(), xdata);
  plot(x, y);
}

void plot(const Vector& xdata, const Vector& ydata) {
  const ArrayOfVector x(1, xdata);
  const ArrayOfVector y(1, ydata);
  plot(x, y);
}

void plot(const Vector& ydata) {
  const Vector xdata(0.0, ydata.size(), 1.0);
  plot(xdata, ydata);
}
} // namespace ARTSGUI
