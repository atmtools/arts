#include "propmat.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
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
#include "physics_funcs.h"
#include "propagationmatrix.h"
#include "species_tags.h"

namespace ARTSGUI {
namespace PropmatClearsky {
constexpr Numeric xscale(Numeric x, XScaling scale) noexcept {
  switch (scale) {
    case XScaling::Hz:
      return x;
    case XScaling::GHz:
      return 1e-9 * x;
    case XScaling::THz:
      return 1e-12 * x;
    case XScaling::Angcm:
      return Conversion::freq2angcm(x);
    case XScaling::Kaycm:
      return Conversion::freq2kaycm(x);
    case XScaling::m:
      return Conversion::freq2wavelen(x);
    case XScaling::nm:
      return 1e9 * Conversion::freq2wavelen(x);
    case XScaling::Angfreq:
      return Conversion::freq2angfreq(x);
    case XScaling::FINAL: { /* leave last */
    }
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

constexpr std::pair<Numeric, Numeric> xunscale(Numeric min,
                                               Numeric max,
                                               XScaling scale) noexcept {
  switch (scale) {
    case XScaling::Hz:
      return {min, max};
    case XScaling::GHz:
      return {1e9 * min, 1e9 * max};
    case XScaling::THz:
      return {1e12 * min, 1e12 * max};
    case XScaling::Angcm:
      return {Conversion::angcm2freq(min), Conversion::angcm2freq(max)};
    case XScaling::Kaycm:
      return {Conversion::kaycm2freq(min), Conversion::kaycm2freq(max)};
    case XScaling::m:
      return {Conversion::wavelen2freq(max), Conversion::wavelen2freq(min)};
    case XScaling::nm:
      return {Conversion::wavelen2freq(1e-9 * max),
              Conversion::wavelen2freq(1e-9 * min)};
    case XScaling::Angfreq:
      return {Conversion::angfreq2freq(min), Conversion::angfreq2freq(max)};
    case XScaling::FINAL: { /* leave last */
    }
  }
  return {std::numeric_limits<Numeric>::quiet_NaN(),
          std::numeric_limits<Numeric>::quiet_NaN()};
}

constexpr std::string_view xunit(XScaling scale) noexcept {
  switch (scale) {
    case XScaling::Hz:
      return "Frequency [Hz]";
    case XScaling::GHz:
      return "Frequency [GHz]";
    case XScaling::THz:
      return "Frequency [THz]";
    case XScaling::Angcm:
      return "Angular Wavenumber [cm-1]";
    case XScaling::Kaycm:
      return "Kayser Wavenumber [cm-1]";
    case XScaling::m:
      return "Wavelength [m]";
    case XScaling::nm:
      return "Wavelength [nm]";
    case XScaling::Angfreq:
      return "Angular Frequency [Hz]";
    case XScaling::FINAL: { /* leave last */
    }
  }
  return "BadUnit";
}

std::array<std::string, PropmatClearsky::enumtyps::XScalingTypes.size()>
xscale_option() {
  std::array<std::string, PropmatClearsky::enumtyps::XScalingTypes.size()> out;
  std::array<char, 100> buf;
  for (size_t i = 0; i < out.size(); i++) {
    buf.fill('\0');
    std::sprintf(buf.data(),
                 "%g",
                 xscale(1e12, PropmatClearsky::enumtyps::XScalingTypes[i]));

    out[i] = var_string("\tScale: ",
                        xunit(PropmatClearsky::enumtyps::XScalingTypes[i]),
                        "; 1e12 Hz -> ",
                        buf.data(),
                        '\t');
  }
  return out;
}

Numeric yscale(YScaling scale,
               Numeric T,
               Numeric P,
               const PropagationMatrix& pm) {
  switch (scale) {
    case YScaling::None:
      return 1.0;
    case YScaling::CrossSection:
      return number_density(P, T);
    case YScaling::Normalize:
      return max(pm.Data());
    case YScaling::FINAL: { /* leave last */
    }
  }
  return 1.0;
}

std::array<std::string, PropmatClearsky::enumtyps::YScalingTypes.size()>
yscale_option() {
  std::array<std::string, PropmatClearsky::enumtyps::YScalingTypes.size()> out;
  for (size_t i = 0; i < out.size(); i++) {
    out[i] =
        var_string('\t', PropmatClearsky::enumstrs::YScalingNames[i], '\t');
  }
  return out;
}

struct DataHolder {
  const PropagationMatrix& pm;
  const Vector& f_grid;
  XScaling xscale_fun;
  Numeric yscale_const;

  static ImPlotPoint Kjj(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->xscale_fun),
            data_ptr->yscale_const * data_ptr->pm.Kjj()[i]};
  }

  static ImPlotPoint K12(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->xscale_fun),
            data_ptr->yscale_const * data_ptr->pm.K12()[i]};
  }

  static ImPlotPoint K13(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->xscale_fun),
            data_ptr->yscale_const * data_ptr->pm.K13()[i]};
  }

  static ImPlotPoint K14(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->xscale_fun),
            data_ptr->yscale_const * data_ptr->pm.K14()[i]};
  }

  static ImPlotPoint K23(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->xscale_fun),
            data_ptr->yscale_const * data_ptr->pm.K23()[i]};
  }

  static ImPlotPoint K24(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->xscale_fun),
            data_ptr->yscale_const * data_ptr->pm.K24()[i]};
  }

  static ImPlotPoint K34(void* self, int i) {
    auto* data_ptr = reinterpret_cast<DataHolder*>(self);
    return {xscale(data_ptr->f_grid[i], data_ptr->xscale_fun),
            data_ptr->yscale_const * data_ptr->pm.K34()[i]};
  }
};

ImPlotLimits draw_propmat(const ComputeValues& v, const DisplayOptions& opts) {
  ImPlotLimits out{};

  const bool select_jac = opts.jacobian_target >= 0 and
                          opts.jacobian_target < v.jacobian_quantities.nelem();

  const PropagationMatrix& pm =
      select_jac ? v.aopm[opts.jacobian_target] : v.pm;

  DataHolder main_data{
      pm,
      v.f_grid,
      opts.xscale,
      (opts.inverse_yscale ? 1.0 / opts.yscale_const : opts.yscale_const) /
          yscale(opts.yscale, v.rtp_temperature, v.rtp_pressure, pm)};

  std::array<char, 100> yscale_str;
  yscale_str.fill('\0');
  sprintf(yscale_str.data(), " (Scale: 1:%g)", main_data.yscale_const);

  const std::string y_axis = var_string(
      select_jac ? MainMenu::absunit(
                       v.jacobian_quantities[opts.jacobian_target].Target())
                 : std::string{"Absorption [1/m]"},
      main_data.yscale_const == 1.0 ? "" : yscale_str.data());

  const std::string_view title = select_jac
                                     ? "Propagation Matrix Partial Derivative"
                                     : "Propagation Matrix";

  if (ImPlot::BeginPlot(
          title.data(), xunit(opts.xscale).data(), y_axis.c_str(), {-1, -1})) {
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
  const auto yscale_display = PropmatClearsky::yscale_option();

  // Local buffers
  ImPlotLimits r{};
  Time start_time{};
  Time end_time{};
  Numeric last_runtime = std::numeric_limits<Numeric>::quiet_NaN();
  std::size_t curpos = PropmatClearsky::n;
  auto fileBrowser = ARTSGUI::Files::xmlfile_chooser();

  // Old copies
  auto old_jacobian_quantities = jacobian_quantities;
  auto old_select_abs_species = select_abs_species;
  auto old_f_grid = f_grid;
  auto old_rtp_mag = rtp_mag;
  auto old_rtp_los = rtp_los;
  auto old_rtp_pressure = rtp_pressure;
  auto old_rtp_temperature = rtp_temperature;
  auto old_rtp_vmr = rtp_vmr;

  // Main loop
  BeginWhileLoopGUI;

  // Main menu bar
  MainMenu::fullscreen(config, window);
  if (MainMenu::exportdata(
          config, fileBrowser, " Export Propagation Matrix ", false))
    config.save_type = 0;
  if (MainMenu::exportdata(config,
                           fileBrowser,
                           " Export Propagation Matrix Derivatives ",
                           false))
    config.save_type = 1;
  MainMenu::quitscreen(config, window);

  // Allow changing the values
  bool updated = false;
  {
    std::lock_guard lock{ctrl.copy};

    updated = MainMenu::change_item("\tjacobian_quantities\t",
                                    jacobian_quantities,
                                    old_jacobian_quantities) or
              updated;

    updated = MainMenu::change_item("\tselect_abs_species\t",
                                    select_abs_species,
                                    old_select_abs_species,
                                    abs_species) or
              updated;

    updated =
        MainMenu::change_item("\tf_grid\t", f_grid, old_f_grid, 1) or updated;

    updated = MainMenu::change_item("\trtp_mag\t",
                                    rtp_mag,
                                    old_rtp_mag,
                                    {"\tU [T]\t", "\tV [T]\t", "\tW [T]\t"}) or
              updated;

    updated =
        MainMenu::change_item("\trtp_los\t",
                              rtp_los,
                              old_rtp_los,
                              {"\tZenith [deg]\t", "\tAzimuth [deg]\t"}) or
        updated;

    updated = MainMenu::change_item("\trtp_pressure\t",
                                    "\t[Pa]\t",
                                    rtp_pressure,
                                    old_rtp_pressure,
                                    std::numeric_limits<Numeric>::min()) or
              updated;

    updated = MainMenu::change_item("\trtp_temperature\t",
                                    "\t[K]\t",
                                    rtp_temperature,
                                    old_rtp_temperature,
                                    std::numeric_limits<Numeric>::min()) or
              updated;

    updated = MainMenu::change_item("\trtp_nlte\t") or updated;

    updated = MainMenu::change_item(
                  "\trtp_vmr\t", rtp_vmr, old_rtp_vmr, abs_species, menu_opt) or
              updated;
  }

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Display")) {
      if (ImGui::BeginMenu("\tY Scale\t")) {
        MainMenu::select_option(disp_options.yscale,
                                PropmatClearsky::enumtyps::YScalingTypes,
                                yscale_display);
        if (disp_options.yscale not_eq PropmatClearsky::YScaling::Normalize) {
          ImGui::InputDouble(
              "\tY Scale Constant\t", &disp_options.yscale_const, 0, 0, "%g");
          ImGui::Checkbox("\tInverse Y Scale Constant\t",
                          &disp_options.inverse_yscale);
          ImGui::Separator();
        } else {
          disp_options.inverse_yscale = false;
          disp_options.yscale_const = 1.0;
        }
        ImGui::EndMenu();
      }
      ImGui::Separator();

      if (ImGui::BeginMenu("\tX Scale\t")) {
        MainMenu::select_option(disp_options.xscale,
                                PropmatClearsky::enumtyps::XScalingTypes,
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

  if (Windows::sub<5, 1, 0, 0, 4, 1>(
          window, Windows::CurrentPosition(), "DrawingWindow")) {
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
                    g.X.Min, g.X.Max, disp_options.xscale);

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
  } else if (end_time not_eq start_time) {
    last_runtime = end_time.Seconds() - start_time.Seconds();
  }

  // Control panel
  if (Windows::sub<5, 5, 4, 0, 1, 1>(
          window, Windows::CurrentPosition(), "ControlWindow")) {
    ImGui::Separator();
    ImGui::Text("\tPanel Control\t");
    ImGui::Separator();
    ImGui::Separator();

    {
      MainMenu::disable_lock disabled{ctrl.run.load()};
      if (ImGui::Button("\tRun Agenda\t", {-1, 0})) {
        start_run(curpos, ctrl, start_time);
      }
      MainMenu::tooltip(
          "Updates the calculations using the current Agenda input", config);
    }

    ImGui::Separator();
    {
      MainMenu::disable_lock disabled{not res[curpos].ok.load()};
      ImGui::Text(" ");
      ImGui::SameLine();
      ImGui::Checkbox("\tUpdate Automatic\t", &res[curpos].auto_update);
      MainMenu::tooltip(
          "Updates the calculations as the Agenda input is updated", config);
      ImGui::Separator();
      ImGui::Text(" ");
      ImGui::SameLine();
      ImGui::Checkbox("\tAutomatic Frequency Grid\t", &res[curpos].auto_f_grid);
      MainMenu::tooltip("Updates the calculations as you move the x-axis",
                        config);
    }
    ImGui::Separator();

    ImGui::Text(" Last run: %g seconds", last_runtime);
    ImGui::Separator();
    if (ctrl.run.load()) {
      ImGui::Text(" Running: %g seconds",
                  end_time.Seconds() - start_time.Seconds());
    } else {
      ImGui::Text(" Not Running");
    }
    ImGui::Separator();
  }
  Windows::end();

  // Information about previous run
  if (Windows::sub<5, 5, 4, 1, 1, 4>(
          window, Windows::CurrentPosition(), "InformationWindow")) {
    ImGui::Separator();
    ImGui::Text("\tPanel Status\t");
    ImGui::Separator();
    ImGui::Separator();
    if (res[curpos].ok.load()) {
      std::lock_guard lock{ctrl.copy};
      auto& v = res[curpos].value;

      // Display current pressure
      ImGui::Text("\tPressure:    %g Pa%c\t",
                  v.rtp_pressure,
                  rtp_pressure == v.rtp_pressure ? ' ' : '*');
      ImGui::Separator();

      // Display current temperature
      ImGui::Text("\tTemperature: %g K%c\t",
                  v.rtp_temperature,
                  rtp_temperature == v.rtp_temperature ? ' ' : '*');
      ImGui::Separator();

      // Display current frequency grid
      ImGui::Text(
          "\tFrequency Grid:\n\t  Start: %g Hz%c\t\n\t  Stop: %g Hz%c\t\n\t  nelem: %ld%c\t",
          v.f_grid[0],
          v.f_grid[0] == f_grid[0] ? ' ' : '*',
          v.f_grid[v.f_grid.nelem() - 1],
          v.f_grid[v.f_grid.nelem() - 1] == f_grid[f_grid.nelem() - 1] ? ' '
                                                                       : '*',
          f_grid.nelem(),
          f_grid.nelem() == v.f_grid.nelem() ? ' ' : '*');
      ImGui::Separator();

      // Display current magnetic field
      ImGui::Text(
          "\tMagnetic Field:\n\t  U: %g T%c\t\n\t  V: %g T%c\t\n\t  W: %g T%c\t",
          v.rtp_mag[0],
          v.rtp_mag[0] == rtp_mag[0] ? ' ' : '*',
          v.rtp_mag[1],
          v.rtp_mag[1] == rtp_mag[1] ? ' ' : '*',
          v.rtp_mag[2],
          v.rtp_mag[2] == rtp_mag[2] ? ' ' : '*');
      ImGui::Separator();

      // Display current LOS
      ImGui::Text(
          "\tLine Of Sight:\n\t  Zenith: %g deg%c\t\n\t  Azimuth: %g deg%c\t",
          v.rtp_los[0],
          v.rtp_los[0] == rtp_los[0] ? ' ' : '*',
          v.rtp_los[1],
          v.rtp_los[1] == rtp_los[1] ? ' ' : '*');
      ImGui::Separator();

      // Display current VMR
      ImGui::Text("\tVMR\t");
      for (Index i = 0; i < abs_species.nelem(); i++) {
        const std::string spec{var_string(abs_species[i])};
        ImGui::Text("\t  %s:\t\n\t    %g%c",
                    spec.c_str(),
                    v.rtp_vmr[i],
                    v.rtp_vmr[i] == rtp_vmr[i] ? ' ' : '*');
      }
      ImGui::Separator();

      // Display current Jacobian
      bool jac_agree =
          jacobian_quantities.nelem() == v.jacobian_quantities.nelem();
      for (Index i = 0; i < jacobian_quantities.nelem() and jac_agree; i++) {
        jac_agree = jac_agree and
                    jacobian_quantities[i].Target().type ==
                        v.jacobian_quantities[i].Target().type and
                    jac_agree and
                    jacobian_quantities[i].Target().atm ==
                        v.jacobian_quantities[i].Target().atm and
                    jac_agree and
                    jacobian_quantities[i].Target().line ==
                        v.jacobian_quantities[i].Target().line and
                    jac_agree and
                    jacobian_quantities[i].Target().qid ==
                        v.jacobian_quantities[i].Target().qid;
      }
      ImGui::Text("\tPartial Derivatives:%c", jac_agree ? ' ' : '*');
      for (auto& jac : v.jacobian_quantities) {
        const auto str = MainMenu::change_item_name(jac.Target());
        ImGui::Text("\t  %s\t", str.c_str());
      }
      ImGui::Separator();

      // Display current select species
      auto spec_str = var_string(v.select_abs_species);
      ImGui::Text("\tSelect Species:%c\t\n\t  %s",
                  spec_str == var_string(select_abs_species) ? ' ' : '*',
                  spec_str.length() == 0 ? "All" : spec_str.c_str());
      ImGui::Separator();
    }
  }
  Windows::end();

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