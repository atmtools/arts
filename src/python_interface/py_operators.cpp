#include <enumsSpectralRadianceUnitType.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <operators.h>

#include "henyey_greenstein.h"
#include "hpy_arts.h"
#include "hpy_numpy.h"
namespace Python {
namespace {
template <typename> using as_object = py::object;

template <typename R, typename... Args>
py::class_<CustomOperator<R, Args...>>& defaults_to_python(py::class_<CustomOperator<R, Args...>>& cop, auto... py_args)
    requires(sizeof...(Args) == sizeof...(py_args)) {
  using self_t    = CustomOperator<R, Args...>;
  using func_t    = self_t::func_t;
  using structs_t = xml_io_stream_functional<func_t>::structs_t;

  // The structs are kept as they are (not wrapped as Python callables), so they can be saved.
  // These must come before the callable, as the structs are themselves callable.
  if constexpr (not std::same_as<void, structs_t>) {
    [&cop]<typename... Ts>(std::type_identity<std::variant<Ts...>>) {
      (cop.def("__init__", [](self_t* op, const Ts& x) { new (op) self_t{.f = x}; }, "f"_a), ...);
      (py::implicitly_convertible<Ts, self_t>(), ...);
    }(std::type_identity<structs_t>{});
  }

  cop.def(
      "__init__",
      [](self_t* op, func_t f) {
        new (op) self_t([f = std::move(f)](Args... args) {
          py::gil_scoped_acquire gil{};
          return f(args...);
        });
      },
      "f"_a);
  py::implicitly_convertible<func_t, self_t>();

  if constexpr ((std::same_as<Args, Numeric> and ...)) {
    // Purely numeric operators also take numpy arrays, element by element
    cop.def(
        "__call__",
        [](const self_t& self, as_object<Args>... args) {
          return vectorize([op = &self](Args... x) { return (*op)(x...); }, args...);
        },
        py_args...);
  } else {
    cop.def("__call__", [](const self_t& self, Args... args) { return self(args...); }, py_args...);
  }

  generic_interface(cop);
  return cop;
}
}  // namespace

void py_operators(py::module_& m) {
  py::class_<WaterEquivalentPressureMK05> mk05(m, "WaterEquivalentPressureMK05");
  generic_interface(mk05);
  mk05.def(
          "__init__",
          [](WaterEquivalentPressureMK05* x, bool only_liquid) { new (x) WaterEquivalentPressureMK05{only_liquid}; },
          "only_liquid"_a)
      .def_rw("only_liquid",
              &WaterEquivalentPressureMK05::only_liquid,
              "Use the pressure over liquid water also below 0 C, otherwise over ice.\n\n.. :class:`bool`")
      .def(
          "__call__",
          [](const WaterEquivalentPressureMK05& self, py::object t) {
            return vectorize([&self](Numeric t) { return self(t); }, t);
          },
          "t"_a,
          "The water equivalent saturation pressure [Pa] at the temperature [K]");
  mk05.doc() = "Water equivalent saturation pressure [Pa] from temperature [K], by Murphy and Koop (2005)";

  py::class_<NumericUnaryOperator> nuop(m, "NumericUnaryOperator");
  defaults_to_python(nuop, "x"_a);

  py::class_<ExtSSACallback> esop(m, "ExtSSACallback");
  defaults_to_python(esop, "frequency"_a, "atm_point"_a);
  esop.doc() = "A callback to get [extinction, ssa] from a single scattering albedo and extinction field.";

  py::class_<NumericBinaryOperator> nbop(m, "NumericBinaryOperator");
  defaults_to_python(nbop, "x"_a, "y"_a);

  py::class_<NumericTernaryOperator> ntop(m, "NumericTernaryOperator");
  defaults_to_python(ntop, "x"_a, "y"_a, "z"_a);

  py::class_<NumericQuaternaryOperator> nqop(m, "NumericQuaternaryOperator");
  defaults_to_python(nqop, "w"_a, "x"_a, "y"_a, "z"_a);

  py::class_<ScatteringGeneralSpectralTROFunc> sgstro(m, "ScatteringGeneralSpectralTROFunc");
  defaults_to_python(sgstro, "atm_point"_a, "f_grid"_a, "degree"_a);
  sgstro.doc() = "A callback to get the bulk scattering properties of a species.";

  // Indirect operators
  py::class_<SpectralRadianceTransformOperator> srtop(m, "SpectralRadianceTransformOperator");
  srtop.def(py::init_implicit<const std::string_view&>())
      .def("__init__",
           [](SpectralRadianceTransformOperator* op, SpectralRadianceTransformOperator::Op::func_t f) {
             new (op) SpectralRadianceTransformOperator(SpectralRadianceTransformOperator::Op(
                 [f = std::move(f)](
                     StokvecVector& x, StokvecMatrix& y, const AscendingGrid& z, const PropagationPathPoint& c) {
                   py::gil_scoped_acquire gil{};
                   f(x, y, z, c);
                 }));
           })
      .def(
          "__call__",
          [](const SpectralRadianceTransformOperator& f,
             StokvecVector&                           x,
             StokvecMatrix&                           y,
             const AscendingGrid&                     z,
             const PropagationPathPoint&              c) { return f(x, y, z, c); },
          "spectral_rad"_a,
          "spectral_rad_jac"_a,
          "freq_grid"_a,
          "ray_point"_a);
  generic_interface(srtop);
  py::implicitly_convertible<SpectralRadianceTransformOperator::Op::func_t, SpectralRadianceTransformOperator>();
}
}  // namespace Python
