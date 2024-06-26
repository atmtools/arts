#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <obsel.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "sorted_grid.h"

namespace Python {
void py_sensor(py::module_& m) try {
  py::class_<SensorPosLos> splos(m, "SensorPosLos");
  workspace_group_interface(splos);
  common_ndarray(splos);
  splos
      .def(
          "__array__",
          [](SensorPosLos& x, py::object dtype, py::object copy) {
            std::array<size_t, 1> shape = {5};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<5>, py::c_contig>(
                    &x, 1, shape.data(), py::handle());
            return np.attr("asarray")(
                w, py::arg("dtype") = dtype, py::arg("copy") = copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")("copy"_a = false); },
          [](SensorPosLos& a, const SensorPosLos& b) { a = b; })
      .def(py::init<Vector3, Vector2>(), "From pos and los")
      .def_rw("pos", &SensorPosLos::pos, "Position")
      .def_rw("los", &SensorPosLos::los, "Line of sight");

  py::class_<SensorPosLosVector> vsplos(m, "SensorPosLosVector");
  workspace_group_interface(vsplos);
  common_ndarray(vsplos);
  vsplos
      .def(
          "__array__",
          [](SensorPosLosVector& x, py::object dtype, py::object copy) {
            std::array<size_t, 2> shape = {static_cast<size_t>(x.size()), 5};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<-1, 5>, py::c_contig>(
                    x.data_handle(), 2, shape.data(), py::handle());
            return np.attr("asarray")(
                w, py::arg("dtype") = dtype, py::arg("copy") = copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")("copy"_a = false); },
          [](SensorPosLosVector& x, Matrix& y) {
            if (y.ncols() != 5) {
              throw std::runtime_error("Bad shape");
            }
            x.resize(y.nrows());
            std::transform(y.begin(), y.end(), x.begin(), [](const auto& row) {
              return SensorPosLos{.pos = {row[0], row[1], row[2]},
                                  .los = {row[3], row[4]}};
            });
          },
          ":class:`~pyarts.arts.Matrix`")
      .def("__getstate__",
           [](const py::object& self) {
             return py::make_tuple(self.attr("value"));
           })
      .def("__setstate__", [](py::object& self, const py::tuple& state) {
        self.attr("value") = state[0];
      });

  py::class_<SensorObsel> so(m, "SensorObsel");
  workspace_group_interface(so);
  so.def(py::init<Vector,
                  AscendingGrid,
                  MuelmatVector,
                  SensorPosLosVector,
                  Stokvec>(),
         "From pos and los")
      .def_rw("f_grid_w", &SensorObsel::f_grid_w, "Frequency weights")
      .def_rw("f_grid", &SensorObsel::f_grid, "Frequency grid")
      .def_rw("poslos_w",
              &SensorObsel::poslos_grid_w,
              "Position and line of sight weights")
      .def_rw("poslos",
              &SensorObsel::poslos_grid,
              "Position and line of sight grid")
      .def_rw(
          "polarization", &SensorObsel::polarization, "Polarization sampling")
      .def(
          "set_frequency_gaussian",
          [](SensorObsel& s,
             const Numeric& f0,
             const Numeric& fwhm,
             const Index& Nfwhm,
             const Index& Nhwhm) {
            s.set_frequency_gaussian(f0, fwhm, Nfwhm, Nhwhm);
          },
          py::arg("f0"),
          py::arg("fwhm"),
          py::arg("Nfwhm") = Index{5},
          py::arg("Nhwhm") = Index{3},
          R"--(Gaussian frequency grid
  
Parameters
----------
f0 : float
    Center frequency
fwhm : float
    Full width at half maximum
Nfwhm : int
    Number of fwhm to include
Nhwhm : int
    Number of half width at half maximum to include
)--")
      .def(
          "set_frequency_lochain",
          [](SensorObsel& s,
             const Vector& f0s,
             const Numeric& width,
             const Index& N,
             const String& filter) {
            s.set_frequency_lochain(DescendingGrid{f0s}, width, N, filter);
          },
          py::arg("f0s"),
          py::arg("width"),
          py::arg("N"),
          py::arg("filter") = String{},
          R"--(Local oscillator style channel selection frequency grid
  
Parameters
----------
f0s : list of descending floats
    Effectively, the local oscillator frequencies for the channels
width : float
    Boxcar width of each sub-channel
N : int
    Number of points per sub-channel
filter : list of int, optional
    Selection of sub-channels to include - one shorter than f0s
    with each element being U for upper bandpass, L for lower bandpass,
    or anything else for full bandpass.  Default is full bandpass.
)--")
      .def("ok", &SensorObsel::ok, "Check if the obsel is valid")
      .def("cutoff_frequency_weights",
           &SensorObsel::cutoff_frequency_weights,
           py::arg("cutoff"),
           py::arg("relative") = true,
           "Cuts out parts of the frequency grid with low weights");

  auto a1 =
      py::bind_vector<ArrayOfSensorObsel, py::rv_policy::reference_internal>(
          m, "ArrayOfSensorObsel");
  workspace_group_interface(a1);
  vector_interface(a1);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize sensors\n", e.what()));
}
}  // namespace Python
