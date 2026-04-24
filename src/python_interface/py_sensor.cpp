#include <arts_omp.h>
#include <antenna_pattern.h>
#include <frequency_bandpass_filters.h>
#include <frequency_channel_selection.h>
#include <frequency_range_selection.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <obsel.h>
#include <python_interface.h>
#include <sensor_meta_info.h>

#include <algorithm>
#include <memory>
#include <stdexcept>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"

namespace Python {
void py_sensor(py::module_& m) try {
  py::class_<SensorPosLos> splos(m, "SensorPosLos");
  generic_interface(splos);
  common_ndarray(splos);
  splos
      .def(
          "__array__",
          [](SensorPosLos& x, py::object dtype, py::object copy)
              -> std::variant<
                  py::ndarray<py::numpy, Numeric, py::shape<5>, py::c_contig>,
                  py::object> {
            std::array<size_t, 1> shape = {5};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<5>, py::c_contig>(
                    &x, 1, shape.data(), py::cast(&x));

            if (not dtype.is_none()) {
              return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy);
            }

            return w.cast((copy.is_none() or not py::bool_(copy))
                              ? py::rv_policy::automatic_reference
                              : py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
          [](SensorPosLos& a, const SensorPosLos& b) { a = b; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`")
      .def(py::init<Vector3, Vector2>(), "From pos and los")
      .def_rw("pos", &SensorPosLos::pos, "Position\n\n.. :class:`Vector3`")
      .def_rw(
          "los", &SensorPosLos::los, "Line of sight\n\n.. :class:`Vector2`");

  py::class_<SensorPosLosVector> vsplos(m, "SensorPosLosVector");
  generic_interface(vsplos);
  common_ndarray(vsplos);
  vsplos.def(
      "sort",
      [](SensorPosLosVector& v, const SensorKeyType x, bool reverse) {
        using enum SensorKeyType;

        switch (x) {
          case alt:
            stdr::sort(v, {}, [](const SensorPosLos& a) { return a.alt(); });
            break;
          case lat:
            stdr::sort(v, {}, [](const SensorPosLos& a) { return a.lat(); });
            break;
          case lon:
            stdr::sort(v, {}, [](const SensorPosLos& a) { return a.lon(); });
            break;
          case zen:
            stdr::sort(v, {}, [](const SensorPosLos& a) { return a.zen(); });
            break;
          case azi:
            stdr::sort(v, {}, [](const SensorPosLos& a) { return a.azi(); });
            break;
          case freq: break;
        }

        if (reverse) {
          if (freq == x) throw std::runtime_error("Cannot reverse frequency");
          stdr::reverse(v);
        }
      },
      "key"_a,
      "reverse"_a = false,
      "Sorts the vector by the given key.");
  vsplos
      .def(
          "__array__",
          [](SensorPosLosVector& x,
             py::object dtype,
             py::object copy) -> std::variant<py::ndarray<py::numpy,
                                                          Numeric,
                                                          py::shape<-1, 5>,
                                                          py::c_contig>,
                                              py::object> {
            std::array<size_t, 2> shape = {static_cast<size_t>(x.size()), 5};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<-1, 5>, py::c_contig>(
                    x.data_handle(), 2, shape.data(), py::cast(&x));

            if (not dtype.is_none()) {
              return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy);
            }

            return w.cast((copy.is_none() or not py::bool_(copy))
                              ? py::rv_policy::automatic_reference
                              : py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
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
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");

  py::class_<sensor::SparseStokvecMatrix> ssm(m, "SparseStokvecMatrix");
  ssm.doc() = "A sparse matrix of Stokvec";
  generic_interface(ssm);
  ssm.def("__getitem__",
          [](const sensor::SparseStokvecMatrix& self, std::tuple<Size, Size> ij)
              -> Stokvec { return self[std::get<0>(ij), std::get<1>(ij)]; });
  ssm.def("__getitem__",
          [](const sensor::SparseStokvecMatrix& self,
             std::tuple<Size, Size, Size> ijk) -> Numeric {
            return self[std::get<0>(ijk), std::get<1>(ijk)][std::get<2>(ijk)];
          });
  ssm.def("__setitem__",
          [](sensor::SparseStokvecMatrix& self,
             std::tuple<Size, Size> ij,
             Stokvec v) { self[std::get<0>(ij), std::get<1>(ij)] = v; });
  ssm.def("__setitem__",
          [](sensor::SparseStokvecMatrix& self,
             std::tuple<Size, Size, Size> ijk,
             Numeric v) {
            self[std::get<0>(ijk), std::get<1>(ijk)][std::get<2>(ijk)] = v;
          });
  ssm.def_prop_ro(
      "dense",
      [](const sensor::SparseStokvecMatrix& self) {
        return StokvecMatrix{self};
      },
      "Get a dense copy.\n\n.. :class:`~pyarts3.arts.StokvecMatrix`");
  ssm.def_prop_ro(
      "value",
      [](const py::object& self) {
        return self.attr("dense").attr("__array__")("copy"_a = true);
      },
      "Allow using numpy on a copy of the class.\n\n.. :class:`~numpy.ndarray`");
  ssm.def(
      "__array__",
      [](const py::object& self,
         const py::object& dtype,
         const py::object& copy) {
        if (not copy.is_none() and not py::bool_(copy))
          throw std::invalid_argument("Must copy dense array");

        return self.attr("value").attr("__array__")("dtype"_a = dtype,
                                                    "copy"_a  = true);
      },
      "dtype"_a = py::none(),
      "copy"_a  = py::none(),
      "Allows :func:`~numpy.array` to be called with the object.");
  common_ndarray(ssm);
  ssm.def(
      "reduce",
      [](const sensor::SparseStokvecMatrix& self,
         Stokvec v,
         bool along_poslos,
         bool along_freq) -> std::variant<Matrix, Vector, Numeric> {
        if (along_poslos and along_freq) {
          Numeric sum = 0.0;
          for (const auto& elem : self.vector()) {
            sum += dot(elem.data, v);
          }
          return sum;
        }

        if (along_freq) {
          Vector out(self.nrows(), 0.0);
          for (const auto& elem : self.vector()) {
            out[elem.irow] += dot(elem.data, v);
          }
          return out;
        }

        if (along_poslos) {
          Vector out(self.ncols(), 0.0);
          for (const auto& elem : self.vector()) {
            out[elem.icol] += dot(elem.data, v);
          }
          return out;
        }

        Matrix out(self.nrows(), self.ncols(), 0.0);
        for (const auto& elem : self.vector()) {
          out[elem.irow, elem.icol] = dot(elem.data, v);
        }
        return out;
      },
      "pol"_a          = Stokvec{1., 0., 0., 0.},
      "along_poslos"_a = false,
      "along_freq"_a   = false,
      R"(Reduces the matrix by polarization.

The reduction happens either along the position/line-of-sight axis,
the frequency axis, or both.  The reduction is the equivalent of performing
``np.einsum("ijk,k->ij", weight_matrix, pol)``

where ``i`` is the position/line-of-sight index, ``j`` is the frequency index,
and ``k`` is the Stokes index if both ``along_poslos`` and ``along_freq``
are ``False``.  If both are ``True``, the result is a single number.
If only one of them is ``True``, the result is a vector the size
of the remaining axis.

Parameters
----------
pol : Stokvec
    The polarization to reduce with.  Default is unpolarized, i.e., ``[1, 0, 0, 0]``.
along_poslos : bool
    If ``True``, reduce along the position/line-of-sight axis.  Default is ``False``.
along_freq : bool
    If ``True``, reduce along the frequency axis.  Default is ``False``.

Returns
-------
Numeric, Vector, or Matrix
    The reduced value, depending on the parameters.
)");

  py::class_<SensorObsel> so(m, "SensorObsel");
  generic_interface(so);
  so.def(py::init<const AscendingGrid&,
                  const SensorPosLosVector&,
                  StokvecMatrix>())
      .def_prop_ro("f_grid",
                   &SensorObsel::f_grid,
                   "Frequency grid\n\n.. :class:`AscendingGrid`")
      .def_prop_ro(
          "weight_matrix",
          [](const SensorObsel& self) { return self.weight_matrix(); },
          "Weights matrix\n\n.. :class:`SparseStokvecMatrix`")
      .def_prop_ro(
          "poslos",
          &SensorObsel::poslos_grid,
          "Position and line of sight grid\n\n.. :class:`SensorPosLosVector`");

  auto a0 = py::bind_vector<Array<SensorPosLosVector>,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfSensorPosLosVector");
  a0.doc() = "Array of SensorPosLosVector";
  generic_interface(a0);
  vector_interface(a0);

  auto a1 =
      py::bind_vector<ArrayOfSensorObsel, py::rv_policy::reference_internal>(
          m, "ArrayOfSensorObsel");
  generic_interface(a1);
  vector_interface(a1);

  a1.def(
      "unique_freq_grids",
      [](const ArrayOfSensorObsel& x) {
        return std::vector{std::from_range,
                           unique_frequency_grids(x) |
                               stdv::transform([](auto& ptr) { return *ptr; })};
      },
      "List of the unique frequency grids");

  a1.def(
      "unique_poslos_grids",
      [](const ArrayOfSensorObsel& x) {
        return std::vector{std::from_range,
                           unique_poslos_grids(x) |
                               stdv::transform([](auto& ptr) { return *ptr; })};
      },
      "List of the unique poslos grids");

  a1.def(
      "collect_freq_grids",
      [](ArrayOfSensorObsel& x) {
        if (x.empty()) return;

        const auto tmp = [&x]() {
          std::set<std::shared_ptr<const AscendingGrid>> freq_set;
          for (const auto& obsel : x) freq_set.insert(obsel.f_grid_ptr());
          if (freq_set.size() == 1)
            return std::make_pair(true, *freq_set.begin());

          std::vector<Numeric> fs;
          fs.reserve(std::accumulate(
              freq_set.begin(),
              freq_set.end(),
              Size(0),
              [](Size a, const auto& b) { return a + b->size(); }));
          for (auto& freqs : freq_set) fs.append_range(*freqs);

          stdr::sort(fs);
          fs.erase(std::unique(fs.begin(), fs.end()), fs.end());
          return std::make_pair(
              false, std::make_shared<const AscendingGrid>(std::move(fs)));
        }();

        //! Clang bug workaround
        const auto& early    = tmp.first;
        const auto& freq_ptr = tmp.second;

        if (early) return;

#pragma omp parallel for if (not arts_omp_in_parallel())
        for (Size i = 0; i < x.size(); i++) {
          sensor::SparseStokvecMatrix w(x[i].poslos_grid().size(),
                                        freq_ptr->size());
          for (const auto& we : x[i].weight_matrix()) {
            const Index ifreq = stdr::distance(
                freq_ptr->begin(),
                stdr::lower_bound(*freq_ptr, x[i].f_grid()[we.icol]));
            w[we.irow, ifreq] = we.data;
          }

          x[i] = {freq_ptr, x[i].poslos_grid_ptr(), std::move(w)};
        }
      },
      R"(Collect all frequency grids in a single grid.

.. tip::
      This generally speeds up radiative transfer calculations.
      The internal weighting matrix is updated internally to
      give the same forward model result.  This matrix is stored
      sparsely, so memory usage should not be a big issue.
      Beware that accessing the weight matrix creates a dense matrix copy,
      which might be large if many frequency grids are combined.

.. note::
      If you need to consider sensor error corrections in frequency,
      such as frequency shifts, stretches or standing waves,
      this method might not have the desired effect.  All channels are considered
      as part of the same "physical" spectrometer, i.e., a standing wave in one channel
      will stretch to all the rest.
)");

  a1.def(
      "collect_poslos_grids",
      [](ArrayOfSensorObsel& x) {
        std::set<std::shared_ptr<const SensorPosLosVector>> poslos_set;
        for (const auto& obsel : x) poslos_set.insert(obsel.poslos_grid_ptr());

        const auto all_poslos_ptr = [&poslos_set]() {
          SensorPosLosVector all_poslos;
          for (const auto& pl : poslos_set) {
            for (auto& poslos : *pl) {
              if (stdr::find(all_poslos, poslos) == all_poslos.end()) {
                all_poslos.push_back(poslos);
              }
            }
          }
          return std::make_shared<const SensorPosLosVector>(
              std::move(all_poslos));
        }();

#pragma omp parallel for if (not arts_omp_in_parallel())
        for (Size i = 0; i < x.size(); i++) {
          sensor::SparseStokvecMatrix w(all_poslos_ptr->size(),
                                        x[i].f_grid().size());

          for (const auto& we : x[i].weight_matrix()) {
            const Index iposlos = stdr::distance(
                all_poslos_ptr->begin(),
                stdr::find(*all_poslos_ptr, x[i].poslos_grid()[we.irow]));
            w[iposlos, we.icol] = we.data;
          }

          x[i] = {x[i].f_grid_ptr(), all_poslos_ptr, std::move(w)};
        }
      },
      R"(Collect all poslos grids in a single grid.

.. tip::
      This generally speeds up radiative transfer calculations.
      The internal weighting matrix is updated internally to
      give the same forward model result.  This matrix is stored
      sparsely, so memory usage should not be a big issue.
      Beware that accessing the weight matrix creates a dense matrix copy,
      which might be large if many poslos grids are combined.

.. note::
      If you need to consider sensor error corrections in line-of-sight or position,
      this method might not have the desired effect.  All channels are considered
      as part of the same "physical" sensor, i.e., a correction in position or line-of-sight
      of one channel will affect all channels.
)");

  py::class_<SensorMetaInfo> smi(m, "SensorMetaInfo");
  generic_interface(smi);
  smi.def_rw(
         "data",
         &SensorMetaInfo::data,
         "The variant gridded field data\n\n.. :class:`~pyarts3.arts.GriddedField1` or :class:`~pyarts3.arts.CameraGriddedField`")
      .def_prop_ro(
          "count",
          &SensorMetaInfo::count,
          "Total number of scalar elements (product of all grid sizes)\n\n.. :class:`int`");

  auto a2 =
      py::bind_vector<ArrayOfSensorMetaInfo, py::rv_policy::reference_internal>(
          m, "ArrayOfSensorMetaInfo");
  a2.doc() = "Array of SensorMetaInfo";
  generic_interface(a2);
  vector_interface(a2);

  auto sch  = py::class_<sensor::Channel>(m, "SensorChannel");
  sch.doc() = "Base class for relative spectrometer channel responses.";
  generic_interface(sch);
  sch.def_prop_ro(
         "response",
         [](const sensor::Channel& self) -> const SortedGriddedField1& {
           return self.channel;
         },
         py::rv_policy::reference_internal,
         "Relative channel response as a gridded field."
         "\n\n.. :class:`~pyarts3.arts.SortedGriddedField1`")
      .def_prop_ro("freq_grid",
                   &sensor::Channel::freq_grid,
                   py::rv_policy::reference_internal,
                   "Relative frequency grid.\n\n.. :class:`AscendingGrid`")
      .def_prop_ro("weights",
                   &sensor::Channel::weights,
                   py::rv_policy::reference_internal,
                   "Channel weights on the relative frequency grid."
                   "\n\n.. :class:`Vector`")
      .def("is_always_relative",
           &sensor::Channel::is_always_relative,
           "Whether the channel grid is anchored at or below zero.");

  auto sbox =
      py::class_<sensor::BoxChannel, sensor::Channel>(m, "SensorBoxChannel");
  sbox.doc() =
      "A channel with uniform weights across a finite relative-frequency interval.";
  sbox.def(py::init<Numeric, Numeric, Size>(),
           "lower"_a,
           "upper"_a,
           "n"_a,
           "Construct a box channel on ``[lower, upper]`` with ``n`` points.")
      .def(
          py::init<Numeric, Size>(),
          "half_width"_a,
          "n"_a,
          "Construct a symmetric box channel on ``[-half_width, half_width]``.")
      .def(py::init<AscendingGrid>(),
           "freq_grid"_a,
           "Construct a box channel directly from a relative frequency grid.");
  generic_interface(sbox);

  auto sdirac = py::class_<sensor::DiracChannel, sensor::Channel>(
      m, "SensorDiracChannel");
  sdirac.doc() = "A single-frequency relative channel.";
  sdirac
      .def(py::init<Numeric>(),
           "frequency"_a,
           "Construct a Dirac channel at one relative frequency.")
      .def(py::init<>(), "Construct a Dirac channel at 0 relative frequency.");
  generic_interface(sdirac);

  auto sgauss = py::class_<sensor::GaussianChannel, sensor::Channel>(
      m, "SensorGaussianChannel");
  sgauss.doc() = "A Gaussian relative channel response.";
  sgauss
      .def(py::init<AscendingGrid, Numeric, Numeric>(),
           "freq_grid"_a,
           "center"_a,
           "std"_a,
           "Construct a Gaussian channel on a custom relative frequency grid.")
      .def(
          py::init<Numeric, Numeric, Size, Size>(),
          "center"_a,
          "std"_a,
          "n"_a,
          "m"_a,
          "Construct a Gaussian channel on ``center +/- m * std`` with ``n`` points.")
      .def(py::init<AscendingGrid, Numeric>(),
           "freq_grid"_a,
           "std"_a,
           "Construct a zero-centered Gaussian channel on a custom grid.")
      .def(py::init<Numeric, Size, Size>(),
           "std"_a,
           "n"_a,
           "m"_a,
           "Construct a zero-centered Gaussian channel on ``+/- m * std``.");
  generic_interface(sgauss);

        auto sap = py::class_<sensor::AntennaPattern>(m, "SensorAntennaPattern");
        sap.doc() =
         "A local 2D antenna pattern on zenith and azimuth offsets from boresight.";
        sap.def(py::init<>(),
          "Construct a pencil-beam antenna pattern.")
         .def_static("pencil",
               &sensor::AntennaPattern::pencil,
               "Construct a pencil-beam antenna pattern.")
         .def_static("gaussian",
               &sensor::AntennaPattern::gaussian,
               "zenith_std"_a,
               "azimuth_std"_a,
               "Construct a Gaussian antenna pattern from zenith and azimuth standard deviations.")
         .def_static("gaussian_fwhm",
               &sensor::AntennaPattern::gaussian_fwhm,
               "zenith_fwhm"_a,
               "azimuth_fwhm"_a,
               "Construct a Gaussian antenna pattern from zenith and azimuth FWHM values.")
         .def_static("lookup",
               [](const AscendingGrid& zenith_grid,
                  const AscendingGrid& azimuth_grid,
                  const Matrix& response_lookup) {
                 return sensor::AntennaPattern::lookup(
                     zenith_grid, azimuth_grid, response_lookup);
               },
               "zenith_grid"_a,
               "azimuth_grid"_a,
               "response_lookup"_a,
               "Construct an antenna pattern from a lookup table on local zenith and azimuth grids.")
         .def_static("lookup",
               [](const ZenGrid& zenith_grid,
                  const AziGrid& azimuth_grid,
                  const Matrix& response_lookup) {
                 return sensor::AntennaPattern::lookup(
                     zenith_grid, azimuth_grid, response_lookup);
               },
               "zenith_grid"_a,
               "azimuth_grid"_a,
               "response_lookup"_a,
               "Construct an antenna pattern from lookup tables defined on typed zenith and azimuth angle grids.")
         .def("set_pencil_beam",
           &sensor::AntennaPattern::set_pencil_beam,
           "Reset the antenna pattern to a pencil beam.")
         .def("set_gaussian",
           &sensor::AntennaPattern::set_gaussian,
           "zenith_std"_a,
           "azimuth_std"_a,
           "Set the antenna pattern to a Gaussian described by standard deviations.")
         .def("set_gaussian_fwhm",
           &sensor::AntennaPattern::set_gaussian_fwhm,
           "zenith_fwhm"_a,
           "azimuth_fwhm"_a,
           "Set the antenna pattern to a Gaussian described by FWHM values.")
         .def("set_lookup",
           [](sensor::AntennaPattern& self,
              const AscendingGrid& zenith_grid,
              const AscendingGrid& azimuth_grid,
              const Matrix& response_lookup) {
             self.set_lookup(zenith_grid, azimuth_grid, response_lookup);
           },
           "zenith_grid"_a,
           "azimuth_grid"_a,
           "response_lookup"_a,
           "Set the antenna pattern from a lookup table.")
         .def("set_lookup",
           [](sensor::AntennaPattern& self,
              const ZenGrid& zenith_grid,
              const AziGrid& azimuth_grid,
              const Matrix& response_lookup) {
             self.set_lookup(zenith_grid, azimuth_grid, response_lookup);
           },
           "zenith_grid"_a,
           "azimuth_grid"_a,
           "response_lookup"_a,
           "Set the antenna pattern from typed zenith and azimuth angle grids.")
         .def("__call__",
           [](const sensor::AntennaPattern& self,
              Numeric delta_zenith,
              Numeric delta_azimuth) {
             return self(delta_zenith, delta_azimuth);
           },
           "delta_zenith"_a,
           "delta_azimuth"_a,
           "Evaluate the antenna pattern at one local zenith and azimuth offset.")
         .def("response",
           [](const sensor::AntennaPattern& self,
              const AscendingGrid& zenith_grid,
              const AscendingGrid& azimuth_grid,
              bool normalize) {
             return normalize ? self.normalized_response(zenith_grid, azimuth_grid)
                  : self.response(zenith_grid, azimuth_grid);
           },
           "zenith_grid"_a,
           "azimuth_grid"_a,
           "normalize"_a = false,
           "Evaluate the antenna pattern on a 2D local angular grid.")
         .def("raw_sensor",
           [](const sensor::AntennaPattern& self,
              const AscendingGrid& dzen_grid,
              const AscendingGrid& dazi_grid,
              bool normalize) {
             return normalize ? self.normalized_raw_sensor(dzen_grid, dazi_grid)
                  : self.raw_sensor(dzen_grid, dazi_grid);
           },
           "dzen_grid"_a,
           "dazi_grid"_a,
           "normalize"_a = false,
           R"(Evaluate the antenna pattern on local angular offsets and return a raw-sensor gridded field.

      The returned field uses grid names ``dzen`` and ``dazi`` and can be passed to
      ``measurement_sensorAddRawSensor``.)")
         .def_ro("type",
              &sensor::AntennaPattern::type,
              "Antenna pattern family.")
         .def_ro("sigma_zenith",
              &sensor::AntennaPattern::sigma_zenith,
              "Zenith standard deviation used by Gaussian patterns.")
         .def_ro("sigma_azimuth",
              &sensor::AntennaPattern::sigma_azimuth,
              "Azimuth standard deviation used by Gaussian patterns.")
         .def_ro("lookup_zenith_grid",
              &sensor::AntennaPattern::lookup_zenith_grid,
              "Zenith grid used by lookup antenna patterns.")
         .def_ro("lookup_azimuth_grid",
              &sensor::AntennaPattern::lookup_azimuth_grid,
              "Azimuth grid used by lookup antenna patterns.")
         .def_ro("lookup_response",
              &sensor::AntennaPattern::lookup_response,
              "Response matrix used by lookup antenna patterns.");
        generic_interface(sap);

  auto shdfr = py::class_<sensor::HeterodyneFrequencyRange>(
      m, "SensorHeterodyneFrequencyRange");
  shdfr.def(py::init<>(),
            R"(Construct an empty staged heterodyne response.

  Stages can then be applied in sequence using :func:`lowpass`, :func:`highpass`,
  :func:`bandpass`, :func:`filter`, and :func:`mix`.)");
  shdfr.def(
      py::init<Numeric, Vector2>(),
      "lo"_a,
      "bandpass"_a,
      R"(Construct a heterodyne response from one ideal bandpass and one LO stage.

  This is shorthand for creating an empty object, applying :func:`bandpass`, and
  then applying :func:`mix`.)");
  shdfr.def(
      "__init__",
      [](sensor::HeterodyneFrequencyRange* v,
         const std::vector<Numeric>& clock_frequencies,
         const std::vector<Vector2>& bandpasses) {
        new (v) sensor::HeterodyneFrequencyRange(clock_frequencies, bandpasses);
      },
      "lo"_a,
      "bandpasses"_a,
      R"(Construct a heterodyne response from a sequence of ideal bandpass and LO stages.

  The sequence is applied as ``bandpass[0] -> lo[0] -> bandpass[1] -> lo[1] -> ...``.)");
  shdfr
      .def_ro("global_ranges",
              &sensor::HeterodyneFrequencyRange::global_ranges,
              "Global frequency range\n\n.. :class:`list[Vector2]`")
      .def_ro("local_ranges",
              &sensor::HeterodyneFrequencyRange::local_ranges,
              "Local frequency range\n\n.. :class:`list[Vector2]`")
      .def("lowpass",
           &sensor::HeterodyneFrequencyRange::apply_lowpass,
           "upper"_a,
           "Apply an ideal lowpass filter on the current local frequency axis.")
      .def(
          "highpass",
          &sensor::HeterodyneFrequencyRange::apply_highpass,
          "lower"_a,
          "Apply an ideal highpass filter on the current local frequency axis.")
      .def(
          "bandpass",
          [](sensor::HeterodyneFrequencyRange& self, const Vector2& bandpass) {
            self.apply_bandpass(bandpass);
          },
          "bandpass"_a,
          "Apply an ideal bandpass filter on the current local frequency axis.")
      .def(
          "filter",
          [](sensor::HeterodyneFrequencyRange& self,
             const SortedGriddedField1& bandpass_filter) {
            self.apply_bandpass(bandpass_filter);
          },
          "bandpass_filter"_a,
          R"(Apply a weighted bandpass filter on the current local frequency axis.

  The filter weights are interpreted on the filter's relative frequency grid and
  are zero outside that grid.)")
      .def("mix",
           &sensor::HeterodyneFrequencyRange::apply_mixer,
           "lo"_a,
           "Apply one heterodyne LO mixing stage.")
      .def(
          "local_response",
          [](const sensor::HeterodyneFrequencyRange& self,
             const Vector& f,
             Size path_index) { return self.local_response(f, path_index); },
          "f"_a,
          "path_index"_a = 0,
          "Evaluate one path response on the current local frequency axis.")
      .def(
          "global_response",
          [](const sensor::HeterodyneFrequencyRange& self,
             const Vector& f,
             Size path_index) { return self.global_response(f, path_index); },
          "f"_a,
          "path_index"_a = 0,
          "Evaluate one path response on the original real-frequency axis.")
      .def(
          "channel_response",
          [](const sensor::HeterodyneFrequencyRange& self,
             const sensor::Channel& channel) {
            return sensor::FrequencyRangeBandpassFilter(
                       self, std::vector<sensor::Channel>{channel})
                .filters.front();
          },
          "channel"_a,
          R"(Compute the real-frequency response for one spectrometer channel.

  The returned gridded field is aggregated across all active mixer paths.)")
      .def(
          "channel_responses",
          [](const sensor::HeterodyneFrequencyRange& self,
             const std::vector<sensor::Channel>& channels) {
            return sensor::FrequencyRangeBandpassFilter(self, channels).filters;
          },
          "channels"_a,
          R"(Compute the real-frequency response for multiple spectrometer channels.

  Each returned gridded field is aggregated across all active mixer paths for the
  matching input channel.)");
  shdfr.doc() = "A staged heterodyne mixer and filter response builder.";
  generic_interface(shdfr);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize sensors\n{}", e.what()));
}
}  // namespace Python
