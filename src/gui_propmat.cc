#include "gui_propmat.h"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "arts_conversions.h"
#include "artstime.h"
#include "debug.h"
#include "gui_error.h"
#include "gui.h"
#include "imgui.h"
#include "implot.h"
#include "jacobian.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "matpack_math.h"
#include "gui_menu.h"
#include "path_point.h"
#include "physics_funcs.h"
#include <rtepack.h>
#include "species_tags.h"

namespace gui {
Numeric no_inf(Numeric x) noexcept {
  return std::isfinite(x) ? x : std::numeric_limits<Numeric>::quiet_NaN();
}

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
    std::snprintf(buf.data(),
                  buf.size(),
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

Numeric yscale(PropmatScaling scale,
               Numeric T,
               Numeric P,
               const PropmatConstVectorView& pm) {
  switch (scale) {
    case PropmatScaling::None:
      return 1.0;
    case PropmatScaling::CrossSection:
      return number_density(P, T);
    case PropmatScaling::Normalize:
      return max(pm).A();
    case PropmatScaling::FINAL: { /* leave last */
    }
  }
  return 1.0;
}

std::array<std::string, PropmatClearsky::enumtyps::PropmatScalingTypes.size()>
propmat_scale_option() {
  std::array<std::string, PropmatClearsky::enumtyps::PropmatScalingTypes.size()>
      out;
  for (size_t i = 0; i < out.size(); i++) {
    out[i] = var_string(
        '\t', PropmatClearsky::enumstrs::PropmatScalingNames[i], '\t');
  }
  return out;
}

Numeric yscale(TramatScaling scale, bool inverse, Numeric x) {
  switch (scale) {
    case TramatScaling::None:
      if (inverse) return 1 - x;
      return x;
    case TramatScaling::dB:
      return (inverse ? -1 : 1) *
             (x > 0 ? 10 * std::log10(x)
                    : std::numeric_limits<Numeric>::quiet_NaN());
    case TramatScaling::FINAL: { /* leave last */
    }
  }
  return x;
}

std::string_view yscale_str(TramatScaling scale, bool inverse) {
  switch (scale) {
    case TramatScaling::None:
      if (inverse) return "One minus Transmission [-]";
      return "Transmission [-]";
    case TramatScaling::dB:
      if (inverse) return "Negative Transmission [dB]";
      return "Transmission [dB]";
    case TramatScaling::FINAL: { /* leave last */
    }
  }
  return "Bad Axis Data";
}

std::array<std::string, PropmatClearsky::enumtyps::TramatScalingTypes.size()>
tramat_scale_option() {
  std::array<std::string, PropmatClearsky::enumtyps::TramatScalingTypes.size()>
      out;
  for (size_t i = 0; i < out.size(); i++) {
    out[i] = var_string(
        '\t', PropmatClearsky::enumstrs::TramatScalingNames[i], '\t');
  }
  return out;
}

constexpr Numeric avg(Numeric a, Numeric b) noexcept {
  return a + 0.5 * (b - a);
}

template <size_t pos>
Numeric range_mean(const PropmatConstVectorView& vec, Index start, Index last) {
  static_assert(pos < 7 and pos >= 0);
  const Numeric scale = 1.0 / static_cast<Numeric>(last - start + 1);
  Numeric out = 0;
  for (Index i = start; i <= last; i++) out += scale * vec[i][pos];
  return out;
}

template <size_t stokes, size_t r, size_t c>
auto range_mean(const MuelmatConstVectorView& t, Index start, Index last) {
  const Numeric scale = 1.0 / static_cast<Numeric>(last - start + 1);
  Numeric out = 0;
  for (Index i = start; i <= last; i++)
    out += scale * t[i](r, c);
  return out;
}

struct TraMatDataHolder {
  const MuelmatConstVectorView& tm;
  const Vector& f_grid;
  XScaling xscale_fun;
  TramatScaling yscale_fun;
  bool inverse_y;
  int running_average;

  [[nodiscard]] int size() const noexcept {
    const Index size = f_grid.size();
    return static_cast<int>(size / running_average) +
           bool(size % running_average);
  }

  [[nodiscard]] std::pair<Index, Index> range(int i) const noexcept {
    const Index f_grid_size{f_grid.size()};
    const Index start{i * running_average};
    const Index last{std::min(start + running_average, f_grid_size) - 1};
    return {start, last};
  }

  template <size_t stokes, size_t r, size_t c>
  static ImPlotPoint T(void* self, int i) {
    auto* data_ptr = static_cast<TraMatDataHolder*>(self);

    const auto [start, last] = data_ptr->range(i);
    const Numeric x = avg(xscale(data_ptr->f_grid[start], data_ptr->xscale_fun),
                          xscale(data_ptr->f_grid[last], data_ptr->xscale_fun));

    return {
        x,
        no_inf(yscale(data_ptr->yscale_fun,
                      data_ptr->inverse_y,
                      range_mean<stokes, r, c>(data_ptr->tm, start, last)))};
  }
};

struct PropMatDataHolder {
  const PropmatConstVectorView& pm;
  const Vector& f_grid;
  XScaling xscale_fun;
  Numeric yscale_const;
  int running_average;

  [[nodiscard]] int size() const noexcept {
    const Index size = f_grid.size();
    return static_cast<int>(size / running_average) +
           bool(size % running_average);
  }

  [[nodiscard]] std::pair<Index, Index> range(int i) const noexcept {
    const Index f_grid_size{f_grid.size()};
    const Index start{i * running_average};
    const Index last{std::min(start + running_average, f_grid_size) - 1};
    return {start, last};
  }

  static ImPlotPoint A(void* self, int i) {
    auto* data_ptr = static_cast<PropMatDataHolder*>(self);

    const auto [start, last] = data_ptr->range(i);
    const Numeric x = avg(xscale(data_ptr->f_grid[start], data_ptr->xscale_fun),
                          xscale(data_ptr->f_grid[last], data_ptr->xscale_fun));
    const Numeric y = no_inf(data_ptr->yscale_const *
                             range_mean<0>(data_ptr->pm, start, last));
    return {x, y};
  }

  static ImPlotPoint B(void* self, int i) {
    auto* data_ptr = static_cast<PropMatDataHolder*>(self);

    const auto [start, last] = data_ptr->range(i);
    const Numeric x = avg(xscale(data_ptr->f_grid[start], data_ptr->xscale_fun),
                          xscale(data_ptr->f_grid[last], data_ptr->xscale_fun));
    const Numeric y = no_inf(data_ptr->yscale_const *
                             range_mean<1>(data_ptr->pm, start, last));
    return {x, y};
  }

  static ImPlotPoint C(void* self, int i) {
    auto* data_ptr = static_cast<PropMatDataHolder*>(self);

    const auto [start, last] = data_ptr->range(i);
    const Numeric x = avg(xscale(data_ptr->f_grid[start], data_ptr->xscale_fun),
                          xscale(data_ptr->f_grid[last], data_ptr->xscale_fun));
    const Numeric y = no_inf(data_ptr->yscale_const *
                             range_mean<2>(data_ptr->pm, start, last));
    return {x, y};
  }

  static ImPlotPoint D(void* self, int i) {
    auto* data_ptr = static_cast<PropMatDataHolder*>(self);

    const auto [start, last] = data_ptr->range(i);
    const Numeric x = avg(xscale(data_ptr->f_grid[start], data_ptr->xscale_fun),
                          xscale(data_ptr->f_grid[last], data_ptr->xscale_fun));
    const Numeric y = no_inf(data_ptr->yscale_const *
                             range_mean<3>(data_ptr->pm, start, last));
    return {x, y};
  }

  static ImPlotPoint U(void* self, int i) {
    auto* data_ptr = static_cast<PropMatDataHolder*>(self);

    const auto [start, last] = data_ptr->range(i);
    const Numeric x = avg(xscale(data_ptr->f_grid[start], data_ptr->xscale_fun),
                          xscale(data_ptr->f_grid[last], data_ptr->xscale_fun));
    const Numeric y = no_inf(data_ptr->yscale_const *
                             range_mean<4>(data_ptr->pm, start, last));
    return {x, y};
  }

  static ImPlotPoint V(void* self, int i) {
    auto* data_ptr = static_cast<PropMatDataHolder*>(self);

    const auto [start, last] = data_ptr->range(i);
    const Numeric x = avg(xscale(data_ptr->f_grid[start], data_ptr->xscale_fun),
                          xscale(data_ptr->f_grid[last], data_ptr->xscale_fun));
    const Numeric y = no_inf(data_ptr->yscale_const *
                             range_mean<5>(data_ptr->pm, start, last));
    return {x, y};
  }

  static ImPlotPoint W(void* self, int i) {
    auto* data_ptr = static_cast<PropMatDataHolder*>(self);

    const auto [start, last] = data_ptr->range(i);
    const Numeric x = avg(xscale(data_ptr->f_grid[start], data_ptr->xscale_fun),
                          xscale(data_ptr->f_grid[last], data_ptr->xscale_fun));
    const Numeric y = no_inf(data_ptr->yscale_const *
                             range_mean<6>(data_ptr->pm, start, last));
    return {x, y};
  }
};

ImPlotLimits draw_propmat(const ComputeValues& v, const DisplayOptions& opts) {
  ImPlotLimits out{};

  const PropmatConstVectorView pm = v.pm;

  PropMatDataHolder main_data{
      pm,
      v.f_grid,
      opts.xscale,
      (opts.inverse_propmat_scale ? 1.0 / opts.propmat_scale_const
                                  : opts.propmat_scale_const) /
          yscale(opts.propmat_scale, v.atm_point.temperature, v.atm_point.pressure, pm),
      opts.smooth_counter};

  std::array<char, 100> yscale_str;
  yscale_str.fill('\0');
  snprintf(yscale_str.data(),
           yscale_str.size(),
           " (Scale: 1:%g)",
           main_data.yscale_const);

  const std::string y_axis = var_string("Absorption [1/m]",
      main_data.yscale_const == 1.0 ? "" : yscale_str.data());

  const std::string_view title = "Propagation Matrix";

  if (ImPlot::BeginPlot(
          title.data(), xunit(opts.xscale).data(), y_axis.c_str(), {-1, -1})) {
    ImPlot::PlotLineG(
        "A", PropMatDataHolder::A, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "B", PropMatDataHolder::B, &main_data, int(v.f_grid.size()));
    ImPlot::PlotLineG(
        "C", PropMatDataHolder::C, &main_data, int(v.f_grid.size()));
    ImPlot::PlotLineG(
        "D", PropMatDataHolder::D, &main_data, int(v.f_grid.size()));
    ImPlot::PlotLineG(
        "U", PropMatDataHolder::U, &main_data, int(v.f_grid.size()));
    ImPlot::PlotLineG(
        "V", PropMatDataHolder::V, &main_data, int(v.f_grid.size()));
    ImPlot::PlotLineG(
        "W", PropMatDataHolder::W, &main_data, int(v.f_grid.size()));

    out = ImPlot::GetPlotLimits();
    ImPlot::EndPlot();
  }

  return out;
}

ImPlotLimits draw_tramat(const ComputeValues& v, const DisplayOptions& opts) {
  ImPlotLimits out{};

  const MuelmatConstVectorView tm = v.tm;

  TraMatDataHolder main_data{tm,
                             v.f_grid,
                             opts.xscale,
                             opts.tramat_scale,
                             opts.inverse_tramat_scale,
                             opts.smooth_counter};

  const std::string_view title = "Transmission Matrix";

  if (ImPlot::BeginPlot(
          title.data(),
          xunit(opts.xscale).data(),
          yscale_str(main_data.yscale_fun, main_data.inverse_y).data(),
          {-1, -1})) {
    ImPlot::PlotLineG(
        "T11", TraMatDataHolder::T<4, 0, 0>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T12", TraMatDataHolder::T<4, 0, 1>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T13", TraMatDataHolder::T<4, 0, 2>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T14", TraMatDataHolder::T<4, 0, 3>, &main_data, main_data.size());

    ImPlot::PlotLineG(
        "T21", TraMatDataHolder::T<4, 1, 0>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T22", TraMatDataHolder::T<4, 1, 1>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T23", TraMatDataHolder::T<4, 1, 2>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T24", TraMatDataHolder::T<4, 1, 3>, &main_data, main_data.size());

    ImPlot::PlotLineG(
        "T31", TraMatDataHolder::T<4, 2, 0>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T32", TraMatDataHolder::T<4, 2, 1>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T33", TraMatDataHolder::T<4, 2, 2>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T34", TraMatDataHolder::T<4, 2, 3>, &main_data, main_data.size());

    ImPlot::PlotLineG(
        "T41", TraMatDataHolder::T<4, 3, 0>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T42", TraMatDataHolder::T<4, 3, 1>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T43", TraMatDataHolder::T<4, 3, 2>, &main_data, main_data.size());
    ImPlot::PlotLineG(
        "T44", TraMatDataHolder::T<4, 3, 3>, &main_data, main_data.size());

    out = ImPlot::GetPlotLimits();
    ImPlot::EndPlot();
  }

  return out;
}

ImPlotLimits draw(const ComputeValues& v, const DisplayOptions& opts) {
  return opts.transmission ? draw_tramat(v, opts) : draw_propmat(v, opts);
}

void start_run(Index pos, Control& ctrl, Time& start_time, Time& end_time) {
  if (not ctrl.run.load()) {
    ctrl.pos.store(int(pos));
    ctrl.run.store(true);
    start_time = Time{};
    end_time = start_time;
  }
}
}  // namespace PropmatClearsky

void propmat(PropmatClearsky::ResultsArray& res,
             PropmatClearsky::Control& ctrl,
             ArrayOfSpeciesTag& select_abs_species,
             Vector& f_grid,
             PropagationPathPoint& path_point,
             AtmPoint& atm_point,
             Numeric& transmission_distance,
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
  const auto propmat_scale_display = PropmatClearsky::propmat_scale_option();
  const auto tramat_scale_display = PropmatClearsky::tramat_scale_option();

  // Local buffers
  ImPlotLimits limits{};
  Time start_time{};
  Time end_time{start_time};
  Numeric last_runtime = std::numeric_limits<Numeric>::quiet_NaN();
  std::size_t curpos = PropmatClearsky::n;
  auto fileBrowser = Files::xmlfile_chooser();

  // Old copies
  auto old_select_abs_species = select_abs_species;
  auto old_f_grid = f_grid;
  auto old_rtp_mag = atm_point.mag;
  auto old_path_point = path_point;
  auto old_rtp_pressure = atm_point.pressure;
  auto old_rtp_temperature = atm_point.temperature;
  auto old_rtp_vmr = atm_point;
  auto old_transmission_distance = transmission_distance;

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

    updated = MainMenu::change_item("\tselect_abs_species\t",
                                    select_abs_species,
                                    old_select_abs_species,
                                    abs_species) or
              updated;

    updated =
        MainMenu::change_item("\tf_grid\t", f_grid, old_f_grid, 1) or updated;

    updated = MainMenu::change_item("\trtp_mag\t",
                                    atm_point.mag,
                                    old_rtp_mag,
                                    {"\tU [T]\t", "\tV [T]\t", "\tW [T]\t"}) or
              updated;

    updated =
        MainMenu::change_item("\tpath_point\t",
                              path_point,
                              old_path_point) or
        updated;

    updated = MainMenu::change_item("\trtp_pressure\t",
                                    "\t[Pa]\t",
                                    atm_point.pressure,
                                    old_rtp_pressure,
                                    std::numeric_limits<Numeric>::min()) or
              updated;

    updated = MainMenu::change_item("\trtp_temperature\t",
                                    "\t[K]\t",
                                    atm_point.temperature,
                                    old_rtp_temperature,
                                    std::numeric_limits<Numeric>::min()) or
              updated;

    updated = MainMenu::change_item("\trtp_nlte\t") or updated;

    updated = MainMenu::change_item(
                  "\trtp_vmr\t", atm_point, old_rtp_vmr, abs_species, menu_opt) or
              updated;

    updated = MainMenu::change_item("\tTransmission Distance\t",
                                    "\t[m]\t",
                                    transmission_distance,
                                    old_transmission_distance,
                                    0) or
              updated;
  }

  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Display")) {
      // Scale tPropagation
      if (ImGui::BeginMenu("\tPropagation Matrix Scale\t")) {
        MainMenu::select_option(disp_options.propmat_scale,
                                PropmatClearsky::enumtyps::PropmatScalingTypes,
                                propmat_scale_display);
        if (disp_options.propmat_scale not_eq
            PropmatClearsky::PropmatScaling::Normalize) {
          ImGui::Text(" ");
          ImGui::SameLine();
          ImGui::InputDouble("\tScale Constant\t",
                             &disp_options.propmat_scale_const,
                             0,
                             0,
                             "%g");
          ImGui::Text(" ");
          ImGui::SameLine();
          ImGui::Checkbox("\tInverse Scale Constant\t",
                          &disp_options.inverse_propmat_scale);
          ImGui::Separator();
        } else {
          disp_options.inverse_propmat_scale = false;
          disp_options.propmat_scale_const = 1.0;
        }
        ImGui::EndMenu();
      }
      MainMenu::tooltip("Various ways to scale the Propagation Matrix", config);
      ImGui::Separator();

      // Scale Transmission
      if (ImGui::BeginMenu("\tTransmission Matrix Scale\t")) {
        MainMenu::select_option(disp_options.tramat_scale,
                                PropmatClearsky::enumtyps::TramatScalingTypes,
                                tramat_scale_display);
        ImGui::Text(" ");
        ImGui::SameLine();
        ImGui::Checkbox("\tInverse?\t", &disp_options.inverse_tramat_scale);
        ImGui::Separator();
        ImGui::EndMenu();
      }
      MainMenu::tooltip("Various ways to scale the Transmission Matrix",
                        config);
      ImGui::Separator();

      // Scale X
      if (ImGui::BeginMenu("\tFrequency Scale\t")) {
        MainMenu::select_option(disp_options.xscale,
                                PropmatClearsky::enumtyps::XScalingTypes,
                                xscale_display);
        ImGui::EndMenu();
      }
      MainMenu::tooltip("Various ways to scale the X-axis", config);
      ImGui::Separator();

      // Running average
      ImGui::Text(" ");
      ImGui::SameLine();
      ImGui::InputInt("\tBinning Count\t", &disp_options.smooth_counter);
      MainMenu::tooltip("Display the results by averaging neighboring grid",
                        config);
      disp_options.smooth_counter = std::clamp(
          disp_options.smooth_counter, 1, std::numeric_limits<int>::max());
      ImGui::Separator();

      // Toggle transmission
      ImGui::Text(" ");
      ImGui::SameLine();
      ImGui::Checkbox("Show As Transmission", &disp_options.transmission);
      ImGui::Separator();
      MainMenu::tooltip(
          "Displays the results as transmission instead of propagation",
          config);

      ImGui::EndMenu();
    }
    MainMenu::tooltip("Options for changing how values are displayed", config);

    ImGui::EndMainMenuBar();
  }
  MainMenu::show_plot_controls();
  MainMenu::show_propmat_controls();

  if (Windows::sub<5, 1, 0, 0, 4, 1>(
          window, Windows::CurrentPosition(), "DrawingWindow")) {
    if (ImGui::BeginTabBar("TabBar")) {
      for (std::size_t i = 0; i < PropmatClearsky::n; i++) {
        std::string opt{var_string(' ', "Plot panel ", i, ' ')};
        if (ImGui::BeginTabItem(opt.c_str())) {
          curpos = i;

          // Run the update if input has be updated
          if (updated and res[i].auto_update) {
            start_run(i, ctrl, start_time, end_time);
          }

          if (res[i].ok.load()) {
            std::lock_guard lock{ctrl.copy};
            auto new_limits = PropmatClearsky::draw(res[i].value, disp_options);

            // If automatic f_grid
            if (new_limits.X.Min not_eq limits.X.Min or
                new_limits.X.Max not_eq limits.X.Max) {
              if (not ctrl.run.load()) limits = new_limits;
              if (res[i].auto_f_grid) {
                auto [min, max] = PropmatClearsky::xunscale(
                    new_limits.X.Min, new_limits.X.Max, disp_options.xscale);

                if (std::isnormal(min) and std::isnormal(max)) {
                  min = std::clamp<Numeric>(
                      min, 1, std::numeric_limits<Numeric>::max());
                  max =
                      std::clamp<Numeric>(max,
                                          new_limits.X.Min + 1,
                                          std::numeric_limits<Numeric>::max());
                  nlinspace(f_grid, min, max, f_grid.size());
                  start_run(i, ctrl, start_time, end_time);
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
      MainMenu::scoped_disable disabled{ctrl.run.load()};
      if (ImGui::Button("\tRun Agenda\t", {-1, 0})) {
        start_run(curpos, ctrl, start_time, end_time);
      }
      MainMenu::tooltip(
          "Updates the calculations using the current Agenda input", config);
    }

    ImGui::Separator();
    {
      MainMenu::scoped_disable disabled{not res[curpos].ok.load()};
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
                  v.atm_point.pressure,
                  atm_point.pressure == v.atm_point.pressure ? ' ' : '*');
      ImGui::Separator();

      // Display current temperature
      ImGui::Text("\tTemperature: %g K%c\t",
                  v.atm_point.temperature,
                  atm_point.temperature == v.atm_point.temperature ? ' ' : '*');
      ImGui::Separator();

      // Display current frequency grid
      ImGui::Text(
          "\tFrequency Grid:\n\t  Start: %g Hz%c\t\n\t  Stop: %g Hz%c\t\n\t  size: %" PRId64
          "%c\t",
          v.f_grid[0],
          v.f_grid[0] == f_grid[0] ? ' ' : '*',
          v.f_grid[v.f_grid.size() - 1],
          v.f_grid[v.f_grid.size() - 1] == f_grid[f_grid.size() - 1] ? ' '
                                                                       : '*',
          f_grid.size(),
          f_grid.size() == v.f_grid.size() ? ' ' : '*');
      ImGui::Separator();

      // Display current magnetic field
      ImGui::Text(
          "\tMagnetic Field:\n\t  U: %g T%c\t\n\t  V: %g T%c\t\n\t  W: %g T%c\t",
          v.atm_point.mag[0],
          v.atm_point.mag[0] == atm_point.mag[0] ? ' ' : '*',
          v.atm_point.mag[1],
          v.atm_point.mag[1] == atm_point.mag[1] ? ' ' : '*',
          v.atm_point.mag[2],
          v.atm_point.mag[2] == atm_point.mag[2] ? ' ' : '*');
      ImGui::Separator();

      // Display current LOS
      ImGui::Text(
          "\tLine Of Sight:\n\t  Zenith: %g deg%c\t\n\t  Azimuth: %g deg%c\t",
          v.path_point.los[0],
          v.path_point.los[0] == path_point.los[0] ? ' ' : '*',
          v.path_point.los[1],
          v.path_point.los[1] == path_point.los[1] ? ' ' : '*');
      ImGui::Separator();

      // Display current VMR
      ImGui::Text("\tVMR\t");
      for (Size i = 0; i < abs_species.size(); i++) {
        const std::string spec{var_string(abs_species[i])};
        ImGui::Text("\t  %s:\t\n\t    %g%c",
                    spec.c_str(),
                    v.atm_point[abs_species[i].Species()],
                    v.atm_point[abs_species[i].Species()] == atm_point[abs_species[i].Species()] ? ' ' : '*');
      }
      ImGui::Separator();

      // Display current select species
      auto spec_str = var_string(v.select_abs_species);
      ImGui::Text("\tSelect Species:%c\t\n\t  %s",
                  spec_str == var_string(select_abs_species) ? ' ' : '*',
                  spec_str.length() == 0 ? "All" : spec_str.c_str());
      ImGui::Separator();

      // Display current transmissio distance
      ImGui::Text("\tTransmission Distance: %g m%c\t",
                  v.transmission_distance,
                  v.transmission_distance == transmission_distance ? ' ' : '*');
      ImGui::Separator();
    }
  }
  Windows::end();

  // Save the data to file?
  if (curpos < res.size() and res[curpos].ok.load()) {
    Files::save_data(config, fileBrowser, res[curpos].value.pm);
  }

  if (ctrl.error.load()) {
    res[curpos].auto_f_grid = false;
    res[curpos].auto_update = false;
    switch (error(ctrl.errmsg)) {
      case ErrorStatus::Exit:
        ctrl.exit.store(true);
        [[fallthrough]];
      case ErrorStatus::Continue:
        ctrl.error.store(false);
        [[fallthrough]];
      case ErrorStatus::OnHold: {
      }
    }
  }

  if (ctrl.exit.load()) glfwSetWindowShouldClose(window, 1);
  // End of main loop
  EndWhileLoopGUI;

  // Clean Graphics data
  CleanupGUI;

  // Set the exit parameter
  ctrl.exit.store(true);
}
}  // namespace gui
