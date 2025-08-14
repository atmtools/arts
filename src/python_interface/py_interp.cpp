#include <lagrange_interp.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>

#include <variant>

#include "nanobind/nanobind.h"

namespace Python {
using id = lagrange_interp::lag_t<-1, lagrange_interp::identity>;
using lc = lagrange_interp::lag_t<-1, lagrange_interp::loncross>;

using PythonLags = std::variant<id, lc>;

template <lagrange_interp::lagrange_type T>
py::class_<T>& interp_class_(py::class_<T>& cl) {
  cl.def(
        "__init__",
        [](T* l, const Numeric x, const Vector& xi, const Index polyorder) {
          new (l) T(xi, x, polyorder);
        },
        "val"_a,
        "vec"_a,
        "polyorder"_a = Index{1},
        "Construct a Lagrange interpolation object")
      .def_prop_ro(
          "order",
          [](const T& l) { return l.size() - 1; },
          "The order of interpolation")
      .def_rw("pos", &T::indx, "The interpolation positions")
      .def_rw("weights", &T::data, "The interpolation weights")
      .def("__getstate__",
           [](const T& self) { return std::make_tuple(self.indx, self.data); })
      .def("__setstate__",
           [](T* l,
              const std::tuple<std::vector<Index>, std::vector<Numeric>>&
                  state) {
             new (l) T{};

             l->indx = std::get<0>(state);
             l->data = std::get<1>(state);
           });

  return cl;
}

template <lagrange_interp::lagrange_type_list T, typename... Ts>
auto& interp_vec_class(py::class_<T, Ts...>& cl) {
  cl.def(
      "__init__",
      [](T* l,
         const Vector& xi,
         const Vector& xn,
         const Index polyorder,
         Numeric extrapol) {
        new (l)
            T(lagrange_interp::make_lags<typename T::value_type::transform_t>(
                xi, xn, polyorder, extrapol));
      },
      "grid"_a,
      "new_grid"_a,
      "polyorder"_a     = Index{1},
      "extrapolation"_a = 0.0,
      "Construct a vector of Lagrange interpolation objects");

  return cl;
}

void py_interp(py::module_& m) {
  auto interp = m.def_submodule("interp", "Interpolation methods");

  py::class_<id> id_(interp, "Lagrange");
  py::class_<lc> lc_(interp, "LagrangeCyclic");
  interp_class_<id>(id_).doc() = "Lagrange interpolation object";
  interp_class_<lc>(lc_).doc() =
      "Lagrange cyclic interpolation object for [-180, 180)";

  auto vid_ =
      py::bind_vector<std::vector<id>, py::rv_policy::reference_internal>(
          interp, "ArrayOfLagrange");
  auto vlc_ =
      py::bind_vector<std::vector<lc>, py::rv_policy::reference_internal>(
          interp, "ArrayOfLagrangeCyclic");
  interp_vec_class(vid_).doc() = "List of :class:`~pyarts.arts.Lagrange`";
  interp_vec_class(vlc_).doc() = "List of :class:`~pyarts.arts.LagrangeCyclic`";

  interp.def(
      "interp1d",
      [](const Vector& yi,
         const Vector& xi,
         const Vector& x,
         const Index polyorder,
         bool cyclic) {
        if (cyclic)
          return reinterp(yi,
                          lagrange_interp::make_lags<lagrange_interp::loncross>(
                              xi, x, polyorder, 0.0));
        return reinterp(yi,
                        lagrange_interp::make_lags<>(xi, x, polyorder, 0.0));
      },
      "yi"_a,
      "xi"_a,
      "x"_a,
      "polyorder"_a = 1,
      "cyclic"_a    = false,
      "Interpolate a 1D vector");
}
}  // namespace Python
