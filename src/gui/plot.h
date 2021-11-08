#ifndef arts_gui_plot_h
#define arts_gui_plot_h

#include <type_traits>
#include <vector>

#include <matpackI.h>
#include <matpack_complex.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "gui.h"
#pragma GCC diagnostic pop

namespace ARTSGUI {

//! Static data for making the plots nicer to look at
struct PlotConfig {
  //! Frame title
  static std::string Frame;
  
  //! X label
  static std::string X;
  
  //! Y label
  static std::string Y;
  
  //! Plot title
  static std::string Title;
  
  //! Line legend (must be exactly the right length, or will default to standard legend)
  static std::vector<std::string> Legend;
};

/** Main plotting function
 * 
 * This function is the main plotting function and will
 * display all x and y.  Special case when y is 1-long
 * in legends and data-saver
 * 
 * @param[in] x All x-vectors
 * @param[in] y All y-vectors
 */
void plot(const ArrayOfVector& x, const ArrayOfVector& y);

/** Plots XY
 * 
 * Wraps plot(const ArrayOfVector& x, const ArrayOfVector& y)
 * with 1-long inputs for x
 * 
 * @param[in] x An x-vector
 * @param[in] y All y-vectors
 */
void plot(const Vector& x, const ArrayOfVector& y);

/** Plots XY
 * 
 * Wraps plot(const ArrayOfVector& x, const ArrayOfVector& y)
 * with 1-long inputs for x and y
 * 
 * @param[in] x An x-vector
 * @param[in] y An y-vector
 */
void plot(const Vector& x, const Vector& y);

/** Plots Y
 * 
 * Wraps void plot(const Vector& x, const Vector& y) with
 * the X-vector defined as Vector(0.0, y.size(), 1.0)
 * 
 * @param[in] y An y-vector
 */
void plot(const Vector& y);

/** Plots any numbers of X-Y pairs
 * 
 * First calls the converter to Vector for each xy-pair
 * 
 * Second moves the X and Y pair into ArrayOfVector before
 * calling plot(const ArrayOfVector& x, const ArrayOfVector& y)
 * 
 * @param[in] xy ... Any number of pairs of variables convertible to Vector
 */
template <typename ... Vectors> void plot(Vectors ... xy) {
  constexpr std::size_t N = sizeof...(Vectors);
  static_assert((N % 2 == 0) or N == 1,
                "Must have a single entry or an even "
                "combinations of X1, Y1, X2, Y2, ... , XN, YN");
  if constexpr (N == 1) {
    const ArrayOfVector y = {Vector(xy)...};
    plot(y[0]);
  } else {
    const std::array<Vector, N> l{Vector(xy)...};
    ArrayOfVector x, y;
    for (std::size_t i=0; i<N; i+=2) {
      x.emplace_back(std::move(l[i]));
      y.emplace_back(std::move(l[i+1]));
    }
    plot(x, y);
  }
}
} // namespace ARTSGUI

#endif  // arts_gui_plot_h
