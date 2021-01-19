#ifndef arts_gui_plot_h
#define arts_gui_plot_h

#include <type_traits>

#include <matpackI.h>
#include <complex.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "gui.h"
#pragma GCC diagnostic pop

namespace ARTSGUI {
void plot(const Vector& y);
void plot(const Vector& xdata, const Vector& ydata);

void plot(const ArrayOfVector& xdata, const ArrayOfVector& ydata);

template <typename ... Lines> void plot(Lines ... lines) {
  static_assert(sizeof...(Lines) % 2 == 0, "Must have even combinations of X1, Y1, X2, Y2, ... , XN, YN");
  const std::array<Vector, sizeof...(Lines)> l{Vector(lines)...};
  ArrayOfVector x, y;
  for (std::size_t i=0; i<sizeof...(Lines); i+=2) {
    x.emplace_back(std::move(l[i]));
    y.emplace_back(std::move(l[i+1]));
  }
  
  plot(x, y);
}
}

#endif  // arts_gui_plot_h
