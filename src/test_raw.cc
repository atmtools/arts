#include <autoarts.h>

#include "artstime.h"
#include "gui/plot.h"
#include "raw.h"
#include "test_utils.h"

#include <Faddeeva/Faddeeva.hh>

Numeric time_sign(const Time t) {
  return 200 + 5 * std::sin(Constant::two_pi*t.Seconds()/86400);
}

Numeric time_srcs(const Time t) {
  return 220 + 25 * std::sin(Constant::two_pi*t.Seconds()/86400);
}

struct MockData {
  Vector cold;
  Vector meas;
  Vector hot;
  Numeric tc;
  Numeric th;
};

MockData mock_data(const Vector& f, const Time& t, const Numeric f0, const Numeric g, const Numeric gd) {
  MockData out;
  out.cold = Vector(f.nelem(), 100);
  add_noise(out.cold, 20);
  out.meas = Vector(f.nelem());
  out.hot = Vector(f.nelem(), 500);
  add_noise(out.hot, 200);
  out.th = 300;
  out.tc = 3;
  
  const Numeric I = time_sign(t);
  const Numeric J = time_srcs(t);
  Vector mock_atm_meas(f.nelem());
  for (Index i=0; i<f.nelem(); i++) {
    const Numeric T = std::exp(- 5e-1 * Faddeeva::w(Complex(f[i] - f0, g) / gd).real());
    mock_atm_meas[i] = (I - J) * T + J;
  }
  add_noise(mock_atm_meas, 1);
  
  for (Index i=0; i<f.nelem(); i++) {
    out.meas[i] = out.cold[i] + (mock_atm_meas[i] - out.tc) * (out.hot[i] - out.cold[i]) / (out.th - out.tc);
  }
  
  return out;
}

int main() {
  const Time now = Time() - Time().seconds_into_day();
  const TimeStep DT(60);
  
  constexpr Index NumSpecChannel=10'000;
  constexpr Numeric F0 = 1.5e9;
  constexpr Numeric g = 1e7;
  constexpr Numeric gd = 1e6;
  
  const Vector f_grid_raw(1e9, NumSpecChannel, 1e9/(NumSpecChannel-1));
  
  // Mock some data
  ArrayOfTime times;
  ArrayOfVector rawdata;
  ArrayOfNumeric tc_arr, th_arr;
  for (Index i=0; i<1440; i+=4) {
    const auto mockraw = mock_data(f_grid_raw, now+i*DT, F0, g, gd);
    for (Index j=0; j<4; j++) {
      times.emplace_back(now + (i+j)*DT);
      tc_arr.emplace_back(mockraw.tc);
      th_arr.emplace_back(mockraw.th);
    }
    rawdata.emplace_back(mockraw.cold);
    rawdata.emplace_back(mockraw.meas);
    rawdata.emplace_back(mockraw.hot);
    rawdata.emplace_back(mockraw.meas);
  }
  const Vector th(th_arr);
  const Vector tc(tc_arr);

  
  // This is the ARTS interface to use the raw data 
  auto ws = ARTS::init();
  ARTS::Var::level0_cold_temperature(ws).value() = tc;
  ARTS::Var::level0_data(ws).value() = rawdata;
  ARTS::Var::level0_hot_temperature(ws).value() = th;
  ARTS::Var::level0_time(ws).value() = times;
  
  // Do CAHA
  ARTS::Method::ybatchCAHA(ws, 0);
  std::cout << "Plotting first, central and last CAHA calibrated data on averaged sensor_time\n";
  ARTSGUI::plot(f_grid_raw, ARTS::Var::ybatch(ws).value()[0],
                f_grid_raw, ARTS::Var::ybatch(ws).value()[ARTS::Var::ybatch(ws).value().nelem()/2],
                f_grid_raw, ARTS::Var::ybatch(ws).value()[ARTS::Var::ybatch(ws).value().nelem()-1]);
  
  // Do time averaging
  ARTS::Method::ybatchTimeAveraging(ws, String{"120 min"});
  std::cout << "Plotting CAHA calibrated data on averaged sensor_time\n";
  ARTSGUI::plot(f_grid_raw, ARTS::Var::ybatch(ws).value());
  
  // Do tropospheric correction
  ARTS::Method::ybatchTroposphericCorrectionNaiveMedianForward(ws, Vector{273});
  std::cout << "Plotting CAHA calibrated data on averaged sensor_time after tropospheric correction\n";
  ARTSGUI::plot(f_grid_raw, ARTS::Var::ybatch(ws).value());
  
  // Focus the data around F0
  ARTS::Method::VectorSet(ws, ARTS::Var::f_grid(ws), f_grid_raw);
  ARTS::Method::ybatchDoublingMeanFocus(ws, F0, 10*gd);
  std::cout << "Plotting CAHA calibrated data on averaged sensor_time after tropospheric correction with simple focus\n";
  ARTSGUI::plot(ARTS::Var::f_grid(ws).value(), ARTS::Var::ybatch(ws).value());
}
