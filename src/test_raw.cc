#include <autoarts.h>

#include "artstime.h"
#include "gui/plot.h"
#include "raw.h"

constexpr Index NumSpecChannel=4096;
struct XmlStruct {
    Time t;
    int chopperpos;
    double Cts1Temperature1;
    double Cts1Temperature2;
    double Cts2Temperature1;
    double Cts2Temperature2;
    double cold_temp;
    double hemt_temp;
    double hot_temp;
    double room_temp1;
    double room_temp2;
    std::array<float, NumSpecChannel> waspam3;
    std::array<float, NumSpecChannel> waspam7;
};

std::ostream& operator<<(std::ostream& os, const XmlStruct& x) {
  return  os << x.t << ' '
             << x.chopperpos << ' '
             << x.Cts1Temperature1 << ' '
             << x.Cts1Temperature2 << ' '
             << x.Cts2Temperature1 << ' '
             << x.Cts2Temperature2 << ' '
             << x.cold_temp << ' '
             << x.hemt_temp << ' '
             << x.hot_temp << ' '
             << x.room_temp1 << ' '
             << x.room_temp2 << ' '
             << '[' <<  x.waspam3.front() << ", ..., " << x.waspam3.back() << ']' << ' '
             << '[' <<  x.waspam7.front() << ", ..., " << x.waspam7.back() << ']';
}

int main() {
  std::fstream file = std::fstream(
    "/home/larsson/Work/kiruna/WASPAM.2021-01-19.0.xml.bin",
    std::fstream::in | std::fstream::binary);
  
  // Read raw data
  std::vector<XmlStruct> data;
  while (not file.eof()) {
    XmlStruct tmp;
    file.read((char *)&tmp, sizeof(Time) + sizeof(int));
    file.read((char *)&tmp.Cts1Temperature1,
              sizeof(XmlStruct) - sizeof(Time) - 2*sizeof(int));
    if (not file.eof()) data.emplace_back(tmp);
  }
  
  // Translation to ARTS style
  ArrayOfVector waspam3(data.size(), Vector(NumSpecChannel));
  ArrayOfVector waspam7(data.size(), Vector(NumSpecChannel));
  Vector tc(data.size()), th(data.size());
  ArrayOfTime time_grid(data.size());
  Index first_c=-1;
  for (std::size_t i=0; i<data.size(); i++) {
    for (Index j=0; j<NumSpecChannel; j++) waspam3[i][j] = Numeric(data[i].waspam3[j]);
    for (Index j=0; j<NumSpecChannel; j++) waspam7[i][j] = Numeric(data[i].waspam7[j]);
    tc[i] = data[i].cold_temp;
    th[i] = data[i].hot_temp;
    time_grid[i] = data[i].t;
    if (first_c < 0 and data[i].chopperpos == 0) first_c = i;
  }
  
  // ARTS f_grid
  Vector f_grid(0, NumSpecChannel, 40e6/NumSpecChannel);
  
  // CAHA calibrations of the data
  const ArrayOfVector calib_waspam3 = Raw::Calibration::caha(waspam3, tc, th, first_c);
  const ArrayOfVector calib_waspam7 = Raw::Calibration::caha(waspam7, tc, th, first_c);
  
  // Daily averages
  Vector fullwaspam3(NumSpecChannel), fullwaspam7(NumSpecChannel);
  Raw::Average::avg(fullwaspam3, calib_waspam3);
  Raw::Average::avg(fullwaspam7, calib_waspam7);
  
  // Tropospheric correction of daily averages
  const Vector fullwaspam3_c = fullwaspam3, fullwaspam7_c = fullwaspam7;
  const Numeric tau_waspam3 = Raw::Correction::naive_tropospheric_singletau_median(fullwaspam3, 273, 2.73);
  const Numeric tau_waspam7 = Raw::Correction::naive_tropospheric_singletau_median(fullwaspam7, 273, 2.73);
  Raw::Correction::naive_tropospheric(fullwaspam3, tau_waspam3, 273);
  Raw::Correction::naive_tropospheric(fullwaspam7, tau_waspam7, 273);
  
  // Apply a mask for values that seem bad
  const Vector fullwaspam3_c2 = fullwaspam3, fullwaspam7_c2 = fullwaspam7;
  const std::vector<bool> bad_waspam3 = Raw::Mask::out_of_bounds(fullwaspam3, nanmean(fullwaspam3) - 1, nanmean(fullwaspam3) + 1);
  const std::vector<bool> bad_waspam7 = Raw::Mask::out_of_bounds(fullwaspam7, nanmean(fullwaspam7) - 1, nanmean(fullwaspam7) + 1);
  Raw::Mask::mask(fullwaspam3, bad_waspam3);
  Raw::Mask::mask(fullwaspam7, bad_waspam7);
  
  // Apply a reduction to the data, resampling close to the line center
  const Vector fullwaspam3_c3 = fullwaspam3, fullwaspam7_c3 = fullwaspam7, f_grid_c = f_grid;
  const ArrayOfIndex red = Raw::Reduce::focus_doublescaling(f_grid, mean(f_grid) + 250e3, 20 * (f_grid[1] - f_grid[0]));
  fullwaspam3 = Raw::Reduce::nanfocus(fullwaspam3, red);
  fullwaspam7 = Raw::Reduce::nanfocus(fullwaspam7, red);
  f_grid = Raw::Reduce::nanfocus(f_grid, red);
  ARTSGUI::plot(f_grid_c, fullwaspam3_c3,
                f_grid_c, fullwaspam7_c3,
                f_grid, fullwaspam3,
                f_grid, fullwaspam7);
}
