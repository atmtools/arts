#include "propmat.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "artstime.h"
#include "constants.h"
#include "debug.h"
#include "gui/gui.h"
#include "imgui.h"
#include "implot.h"
#include "jacobian.h"
#include "math_funcs.h"
#include "matpackI.h"
#include "menu.h"
#include "propagationmatrix.h"
#include "species_tags.h"

namespace ARTSGUI {
namespace PropmatClearsky {
constexpr Numeric xscale(Numeric x, Scaling scale) noexcept {
  switch (scale) {
    case Scaling::Hz:
      return x;
    case Scaling::GHz:
      return 1e-9 * x;
    case Scaling::THz:
      return 1e-12 * x;
    case Scaling::Angcm:
      return Conversion::freq2angcm(x);
    case Scaling::Kaycm:
      return Conversion::freq2kaycm(x);
    case Scaling::m:
      return Conversion::freq2wavelen(x);
    case Scaling::nm:
      return 1e9 * Conversion::freq2wavelen(x);
    case Scaling::Angfreq:
      return Conversion::freq2angfreq(x);
    case Scaling::FINAL: { /* leave last */
    }
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

constexpr std::pair<Numeric, Numeric> xunscale(Numeric min,
                                               Numeric max,
                                               Scaling scale) noexcept {
  switch (scale) {
    case Scaling::Hz:
      return {min, max};
    case Scaling::GHz:
      return {1e9 * min, 1e9 * max};
    case Scaling::THz:
      return {1e12 * min, 1e12 * max};
    case Scaling::Angcm:
      return {Conversion::angcm2freq(min), Conversion::angcm2freq(max)};
    case Scaling::Kaycm:
      return {Conversion::kaycm2freq(min), Conversion::kaycm2freq(max)};
    case Scaling::m:
      return {Conversion::wavelen2freq(max), Conversion::wavelen2freq(min)};
    case Scaling::nm:
      return {Conversion::wavelen2freq(1e-9 * max),
              Conversion::wavelen2freq(1e-9 * min)};
    case Scaling::Angfreq:
      return {Conversion::angfreq2freq(min), Conversion::angfreq2freq(max)};
    case Scaling::FINAL: { /* leave last */
    }
  }
  return {std::numeric_limits<Numeric>::quiet_NaN(),
          std::numeric_limits<Numeric>::quiet_NaN()};
}

constexpr std::string_view xunit(Scaling scale) noexcept {
  switch (scale) {
    case Scaling::Hz:
      return "Frequency [Hz]";
    case Scaling::GHz:
      return "Frequency [GHz]";
    case Scaling::THz:
      return "Frequency [THz]";
    case Scaling::Angcm:
      return "Angular Wavenumber [cm-1]";
    case Scaling::Kaycm:
      return "Kayser Wavenumber [cm-1]";
    case Scaling::m:
      return "Wavelength [m]";
    case Scaling::nm:
      return "Wavelength [nm]";
    case Scaling::Angfreq:
      return "Angular Frequency [Hz]";
    case Scaling::FINAL: { /* leave last */
    }
  }
  return "BadUnit";
}

std::array<std::string, PropmatClearsky::enumtyps::ScalingTypes.size()>
xscale_option() {
  std::array<std::string, PropmatClearsky::enumtyps::ScalingTypes.size()> out;
  std::array<char, 100> buf;
  for (size_t i = 0; i < out.size(); i++) {
    buf.fill('\0');
    std::sprintf(buf.data(),
                 "%g",
                 xscale(1e12, PropmatClearsky::enumtyps::ScalingTypes[i]));

    out[i] = var_string("Scale: ",
                        xunit(PropmatClearsky::enumtyps::ScalingTypes[i]),
                        "; 1e12 Hz -> ",
                        buf.data());
  }
  return out;
}

struct DataHolder {
  const PropagationMatrix& pm;
  const Vector& f_grid;
  Scaling scale;

  static ImPlotPoint Kjj(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->scale),
            data_ptr->pm.Kjj()[i]};
  }

  static ImPlotPoint K12(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->scale),
            data_ptr->pm.K12()[i]};
  }

  static ImPlotPoint K13(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->scale),
            data_ptr->pm.K13()[i]};
  }

  static ImPlotPoint K14(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->scale),
            data_ptr->pm.K14()[i]};
  }

  static ImPlotPoint K23(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->scale),
            data_ptr->pm.K23()[i]};
  }

  static ImPlotPoint K24(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->scale),
            data_ptr->pm.K24()[i]};
  }

  static ImPlotPoint K34(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->scale),
            data_ptr->pm.K34()[i]};
  }
};

ImPlotLimits draw_propmat(const ComputeValues& v, const DisplayOptions& opts) {
  ImPlotLimits out{};
  DataHolder main_data{(opts.jacobian_target < 0 or
                        opts.jacobian_target >= v.jacobian_quantities.nelem())
                           ? v.pm
                           : v.aopm[opts.jacobian_target],
                       v.f_grid,
                       opts.scale};

  std::string y_axis{"Absorption [1/m]"};

  if (ImPlot::BeginPlot("Propagation Matrix",
                        xunit(opts.scale).data(),
                        y_axis.c_str(),
                        {-1, -1})) {
    if (v.pm.StokesDimensions() > 0) {
      ImPlot::PlotLineG(
          "Kjj", DataHolder::Kjj, &main_data, int(v.f_grid.nelem()));
    }

    if (v.pm.StokesDimensions() > 1) {
      ImPlot::PlotLineG(
          "K12", DataHolder::K12, &main_data, int(v.f_grid.nelem()));
    }

    if (v.pm.StokesDimensions() > 2) {
      ImPlot::PlotLineG(
          "K13", DataHolder::K13, &main_data, int(v.f_grid.nelem()));
      ImPlot::PlotLineG(
          "K23", DataHolder::K23, &main_data, int(v.f_grid.nelem()));
    }

    if (v.pm.StokesDimensions() > 3) {
      ImPlot::PlotLineG(
          "K24", DataHolder::K24, &main_data, int(v.f_grid.nelem()));
      ImPlot::PlotLineG(
          "K34", DataHolder::K34, &main_data, int(v.f_grid.nelem()));
    }

    out = ImPlot::GetPlotLimits();
    ImPlot::EndPlot();
  }

  return out;
}

void start_run(Index pos, Control& ctrl, Time& start_time) {
  if (not ctrl.run.load()) {
    ctrl.pos.store(int(pos));
    ctrl.run.store(true);
    start_time = Time{};
  }
}
}  // namespace PropmatClearsky

void propmat(PropmatClearsky::ResultsArray& res,
             PropmatClearsky::Control& ctrl,
             ArrayOfRetrievalQuantity& jacobian_quantities,
             ArrayOfSpeciesTag& select_abs_species,
             Vector& f_grid,
             Vector& rtp_mag,
             Vector& rtp_los,
             Numeric& rtp_pressure,
             Numeric& rtp_temperature,
             EnergyLevelMap&,
             Vector& rtp_vmr,
             const ArrayOfArrayOfSpeciesTag&& abs_species) {
  // Get Graphics data
  InitializeGUI("Propagation Agenda Output", 1280, 720);

  // Our global states are stored in config
  Config config{};
  MainMenu::Options menu_opt{};
  PropmatClearsky::DisplayOptions disp_options{};

  // Our style
  LayoutAndStyleSettings();

  // Local configurations
  const auto xscale_display = PropmatClearsky::xscale_option();

  // Local buffers
  ImPlotLimits r{};
  Time start_time{};
  Time end_time{};
  std::size_t curpos = PropmatClearsky::n;
  auto fileBrowser = ARTSGUI::Files::xmlfile_chooser();

  // Main loop
  BeginWhileLoopGUI;

  // Main menu bar
  MainMenu::fullscreen(config, window);
  if (MainMenu::exportdata(config, fileBrowser, " Export Propagation Matrix ", false)) config.save_type=0;
  if (MainMenu::exportdata(config, fileBrowser, " Export Propagation Matrix Derivatives ", false)) config.save_type=1;
  MainMenu::quitscreen(config, window);

  // Allow changing the values
  bool updated = false;
  {
    std::lock_guard lock{ctrl.copy};

    updated =
        MainMenu::change_item("\tjacobian_quantities\t", jacobian_quantities) or
        updated;

    updated = MainMenu::change_item(
                  "\tselect_abs_species\t", select_abs_species, abs_species) or
              updated;

    updated = MainMenu::change_item("\tf_grid\t", f_grid, 1) or updated;

    updated =
        MainMenu::change_item(
            "\trtp_mag\t", rtp_mag, {"\tU [T]\t", "\tV [T]\t", "\tW [T]\t"}) or
        updated;

    updated =
        MainMenu::change_item("\trtp_los\t",
                              rtp_los,
                              {"\tZenith [deg]\t", "\tAzimuth [deg]\t"}) or
        updated;

    updated = MainMenu::change_item("\trtp_pressure\t",
                                    "\t[Pa]\t",
                                    rtp_pressure,
                                    std::numeric_limits<Numeric>::min()) or
              updated;

    updated = MainMenu::change_item("\trtp_temperature\t",
                                    "\t[K]\t",
                                    rtp_temperature,
                                    std::numeric_limits<Numeric>::min()) or
              updated;

    updated = MainMenu::change_item("\trtp_nlte\t") or updated;

    updated =
        MainMenu::change_item("\trtp_vmr\t", rtp_vmr, abs_species, menu_opt) or
        updated;
  }

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Panels")) {
      for (std::size_t i = 0; i < PropmatClearsky::n; i++) {
        std::string opt{var_string("\tSettings plot panel ", i, '\t')};
        if (ImGui::BeginMenu(opt.c_str())) {
          if (ImGui::Button("\tUpdate Manually\t")) {
            start_run(i, ctrl, start_time);
          }
          if (ImGui::IsItemHovered()) {
            if (ImGui::GetCurrentContext()->HoveredIdTimer >
                config.hover_time_limit) {
              ImGui::SetTooltip(" \n\tUpdates the calculations\t\n ");
              ImGui::BeginTooltip();
              ImGui::EndTooltip();
            }
          }
          ImGui::Separator();
          ImGui::Checkbox("\tUpdate Automatic\t", &res[i].auto_update);
          if (ImGui::IsItemHovered()) {
            if (ImGui::GetCurrentContext()->HoveredIdTimer >
                config.hover_time_limit) {
              ImGui::SetTooltip(
                  " \n\tUpdates the calculations as the atmospheric state is updated\t\n ");
              ImGui::BeginTooltip();
              ImGui::EndTooltip();
            }
          }
          ImGui::Separator();
          ImGui::Checkbox("\tAutomatic Frequency Grid\t", &res[i].auto_f_grid);
          if (ImGui::IsItemHovered()) {
            if (ImGui::GetCurrentContext()->HoveredIdTimer >
                config.hover_time_limit) {
              ImGui::SetTooltip(
                  " \n\tUpdates the calculations as you move the x-axis\t\n"
                  "\n"
                  "\tNote that this is very unstable with inverse frequency scales (wavelengths)\t\n ");
              ImGui::BeginTooltip();
              ImGui::EndTooltip();
            }
          }

          ImGui::EndMenu();
        }
        ImGui::Separator();
      }
      ImGui::EndMenu();
    }
    if (ImGui::IsItemHovered()) {
      if (ImGui::GetCurrentContext()->HoveredIdTimer >
          config.hover_time_limit) {
        ImGui::SetTooltip(" \n\tSettings related to individual panels\t\n ");
        ImGui::BeginTooltip();
        ImGui::EndTooltip();
      }
    }

    if (ImGui::BeginMenu("Display")) {
      if (ImGui::BeginMenu("\tX Scale\t")) {
        MainMenu::select_option(disp_options.scale,
                                PropmatClearsky::enumtyps::ScalingTypes,
                                xscale_display);
        ImGui::EndMenu();
      }
      ImGui::Separator();

      if (ImGui::BeginMenu("\tSelect Plot\t")) {
        MainMenu::select_option(disp_options.jacobian_target,
                                jacobian_quantities);
        ImGui::EndMenu();
      }
      ImGui::Separator();
      if (ImGui::IsItemHovered()) {
        if (ImGui::GetCurrentContext()->HoveredIdTimer >
            config.hover_time_limit) {
          ImGui::SetTooltip(" \n\tSelect main or derivative to display \t\n ");
          ImGui::BeginTooltip();
          ImGui::EndTooltip();
        }
      }

      ImGui::EndMenu();
    }
    if (ImGui::IsItemHovered()) {
      if (ImGui::GetCurrentContext()->HoveredIdTimer >
          config.hover_time_limit) {
        ImGui::SetTooltip(" \n\tScale of the X-axis\t\n ");
        ImGui::BeginTooltip();
        ImGui::EndTooltip();
      }
    }

    ImGui::EndMainMenuBar();
  }

  if (Windows::full(window, Windows::CurrentPosition(), "DrawingWindow")) {
    if (ImGui::BeginTabBar("TabBar")) {
      for (std::size_t i = 0; i < PropmatClearsky::n; i++) {
        std::string opt{var_string(' ', "Plot panel ", i, ' ')};
        if (ImGui::BeginTabItem(opt.c_str())) {
          curpos = i;

          // Run the update if input has be updated
          if (updated and res[i].auto_update) {
            start_run(i, ctrl, start_time);
          }

          if (res[i].ok.load()) {
            std::lock_guard lock{ctrl.copy};
            auto g = PropmatClearsky::draw_propmat(res[i].value, disp_options);

            // If automatic f_grid
            if (g.X.Min not_eq r.X.Min or r.X.Max not_eq g.X.Max) {
              r = g;
              if (res[i].auto_f_grid) {
                auto [min, max] = PropmatClearsky::xunscale(
                    g.X.Min, g.X.Max, disp_options.scale);

                if (std::isnormal(min) and std::isnormal(max)) {
                  min = std::clamp<Numeric>(
                      min, 1, std::numeric_limits<Numeric>::max());
                  max = std::clamp<Numeric>(
                      max, g.X.Min + 1, std::numeric_limits<Numeric>::max());
                  nlinspace(f_grid, min, max, f_grid.nelem());
                  start_run(i, ctrl, start_time);
                }
              }
            }
          }
          ImGui::EndTabItem();
        }
      }
      ImGui::EndTabBar();
    }
  }
  Windows::end();

  // Timer
  if (ctrl.run.load()) {
    end_time = Time{};
  }

  if (std::any_of(
          res.begin(), res.end(), [](auto& x) { return x.ok.load(); })) {
    if (ImGui::BeginMainMenuBar()) {
      std::string msg{
          var_string((ctrl.run.load() ? std::string{"Running: "}
                                      : std::string{"Last run took: "}),
                     end_time.Seconds() - start_time.Seconds(),
                     " s")};
      ImGui::Text("\t%s", msg.c_str());
      ImGui::EndMainMenuBar();
    }
  }

  // Save the data to file?
  if (curpos < res.size() and res[curpos].ok.load()) {
    if (config.save_type == 0)
      ARTSGUI::Files::save_data(config, fileBrowser, res[curpos].value.pm);
    else if (config.save_type == 1)
      ARTSGUI::Files::save_data(config, fileBrowser, res[curpos].value.aopm);
  }

  if (ctrl.exit.load()) glfwSetWindowShouldClose(window, 1);
  // End of main loop
  EndWhileLoopGUI;

  // Clean Graphics data
  CleanupGUI;

  // Set the exit parameter
  ctrl.exit.store(true);
}
}  // namespace ARTSGUI