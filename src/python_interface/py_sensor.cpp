#include <antenna_pattern.h>
#include <arts_omp.h>
#include <frequency_bandpass_filters.h>
#include <frequency_channel_selection.h>
#include <frequency_range_selection.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>
#include <sensor.h>
#include <sensor_builder.h>
#include <sensor_meta_info.h>

#include <algorithm>
#include <memory>
#include <stdexcept>

#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "xml_io_stream.h"

namespace Python {
namespace {
//! The layout and the sensor of a channel through the antenna pattern of the bound class c
template <typename T, typename... Ts> void antenna_calls(py::class_<T, Ts...>& c) {
  c.def_prop_ro("pixels", &T::pixels, "The layout of the pixels.\n\n.. :class:`~pyarts3.arts.sensors.PixelGrid`");
  c.def(
       "__call__",
       [](const T&               self,
          const sensor::Channel& channel,
          const Vector3&         pos,
          const Vector2&         bore_los,
          const Vector2&         ell) {
         if (channel.polarization.scan_rotation != 0.0) {
           throw std::runtime_error(
               "The channel's polarization turns with the scan, so it needs a scan angle: "
               "pass the response at the boresight, or use the sensor builder with a scan geometry");
         }
         return self(std::span{&channel, 1}, std::span{&channel.polarization.state, 1}, pos, bore_los, ell);
       },
       "channel"_a,
       "pos"_a,
       "bore_los"_a,
       "ell"_a,
       R"(The sensor of a channel through the antenna pattern, one element per pixel by row and then column.

The channel's polarization must not depend on the scan.)")
      .def(
          "__call__",
          [](const T&               self,
             const sensor::Channel& channel,
             const Stokvec&         response,
             const Vector3&         pos,
             const Vector2&         bore_los,
             const Vector2&         ell) { return self(std::span{&channel, 1}, std::span{&response, 1}, pos, bore_los, ell); },
          "channel"_a,
          "response"_a,
          "pos"_a,
          "bore_los"_a,
          "ell"_a,
          R"(The sensor of a channel through the antenna pattern, one element per pixel by row and then column.

``response`` is the channel's polarization response at the boresight, e.g. ``channel.polarization.at(scan_angle)``.)");
}
}  // namespace

void py_sensor(py::module_& m) try {
  py::class_<SensorPosLos> splos(m, "SensorPosLos");
  generic_interface(splos);
  common_ndarray(splos);
  splos
      .def(
          "__array__",
          [](SensorPosLos& x,
             py::object    dtype,
             py::object copy) -> std::variant<py::ndarray<py::numpy, Numeric, py::shape<5>, py::c_contig>, py::object> {
            std::array<size_t, 1> shape = {5};
            auto                  np    = py::module_::import_("numpy");
            auto w = py::ndarray<py::numpy, Numeric, py::shape<5>, py::c_contig>(&x, 1, shape.data(), py::cast(&x));

            if (not dtype.is_none()) { return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy); }

            if (copy.is_none() or not py::bool_(copy)) { return w.cast(py::rv_policy::automatic_reference); }
            return w.cast(py::rv_policy::copy);
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
      .def_rw("pos", &SensorPosLos::pos, "Position\n\n.. :class:`~pyarts3.arts.Vector3`")
      .def_rw("los", &SensorPosLos::los, "Line of sight\n\n.. :class:`~pyarts3.arts.Vector2`");

  py::class_<SensorPosLosVector> vsplos(m, "SensorPosLosVector");
  generic_interface(vsplos);
  common_ndarray(vsplos);
  vsplos.def(
      "sort",
      [](SensorPosLosVector& v, const SensorKeyType x, bool reverse) {
        using enum SensorKeyType;

        switch (x) {
          case alt:  stdr::sort(v, {}, [](const SensorPosLos& a) { return a.alt(); }); break;
          case lat:  stdr::sort(v, {}, [](const SensorPosLos& a) { return a.lat(); }); break;
          case lon:  stdr::sort(v, {}, [](const SensorPosLos& a) { return a.lon(); }); break;
          case zen:  stdr::sort(v, {}, [](const SensorPosLos& a) { return a.zen(); }); break;
          case azi:  stdr::sort(v, {}, [](const SensorPosLos& a) { return a.azi(); }); break;
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
          [](SensorPosLosVector& x, py::object dtype, py::object copy)
              -> std::variant<py::ndarray<py::numpy, Numeric, py::shape<-1, 5>, py::c_contig>, py::object> {
            std::array<size_t, 2> shape = {static_cast<size_t>(x.size()), 5};
            auto                  np    = py::module_::import_("numpy");
            auto                  w     = py::ndarray<py::numpy, Numeric, py::shape<-1, 5>, py::c_contig>(
                x.data_handle(), 2, shape.data(), py::cast(&x));

            if (not dtype.is_none()) { return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy); }

            if (copy.is_none() or not py::bool_(copy)) { return w.cast(py::rv_policy::automatic_reference); }
            return w.cast(py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
          [](SensorPosLosVector& x, Matrix& y) {
            if (y.ncols() != 5) { throw std::runtime_error("Bad shape"); }
            x.resize(y.nrows());
            std::transform(y.begin(), y.end(), x.begin(), [](const auto& row) {
              return SensorPosLos{.pos = {row[0], row[1], row[2]}, .los = {row[3], row[4]}};
            });
          },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");

  py::class_<sensor::SparseStokvecMatrix> ssm(m, "SparseStokvecMatrix");
  ssm.doc() = "A sparse matrix of Stokvec";
  generic_interface(ssm);
  ssm.def("__getitem__", [](const sensor::SparseStokvecMatrix& self, std::tuple<Size, Size> ij) -> Stokvec {
    return self[std::get<0>(ij), std::get<1>(ij)];
  });
  ssm.def("__getitem__", [](const sensor::SparseStokvecMatrix& self, std::tuple<Size, Size, Size> ijk) -> Numeric {
    return self[std::get<0>(ijk), std::get<1>(ijk)][std::get<2>(ijk)];
  });
  ssm.def("__setitem__", [](sensor::SparseStokvecMatrix& self, std::tuple<Size, Size> ij, Stokvec v) {
    self[std::get<0>(ij), std::get<1>(ij)] = v;
  });
  ssm.def("__setitem__", [](sensor::SparseStokvecMatrix& self, std::tuple<Size, Size, Size> ijk, Numeric v) {
    self[std::get<0>(ijk), std::get<1>(ijk)][std::get<2>(ijk)] = v;
  });
  ssm.def_prop_ro(
      "dense",
      [](const sensor::SparseStokvecMatrix& self) { return StokvecMatrix{self}; },
      "Get a dense copy.\n\n.. :class:`~pyarts3.arts.StokvecMatrix`");
  ssm.def_prop_ro(
      "value",
      [](const py::object& self) { return self.attr("dense").attr("__array__")("copy"_a = true); },
      "Allow using numpy on a copy of the class.\n\n.. :class:`~numpy.ndarray`");
  ssm.def(
      "__array__",
      [](const py::object& self, const py::object& dtype, const py::object& copy) {
        if (not copy.is_none() and not py::bool_(copy)) throw std::invalid_argument("Must copy dense array");

        return self.attr("value").attr("__array__")("dtype"_a = dtype, "copy"_a = true);
      },
      "dtype"_a = py::none(),
      "copy"_a  = py::none(),
      "Allows :func:`~numpy.array` to be called with the object.");
  common_ndarray(ssm);
  ssm.def(
      "reduce",
      [](const sensor::SparseStokvecMatrix& self, Stokvec v, bool along_poslos, bool along_freq)
          -> std::variant<Matrix, Vector, Numeric> {
        if (along_poslos and along_freq) {
          Numeric sum = 0.0;
          for (const auto& elem : self.vector()) { sum += dot(elem.data, v); }
          return sum;
        }

        if (along_freq) {
          Vector out(self.nrows(), 0.0);
          for (const auto& elem : self.vector()) { out[elem.irow] += dot(elem.data, v); }
          return out;
        }

        if (along_poslos) {
          Vector out(self.ncols(), 0.0);
          for (const auto& elem : self.vector()) { out[elem.icol] += dot(elem.data, v); }
          return out;
        }

        Matrix out(self.nrows(), self.ncols(), 0.0);
        for (const auto& elem : self.vector()) { out[elem.irow, elem.icol] = dot(elem.data, v); }
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

  auto assm = py::bind_vector<std::vector<sensor::SparseStokvecMatrix>, py::rv_policy::reference_internal>(
      m, "ArrayOfSparseStokvecMatrix");
  assm.doc() = "Array of SparseStokvecMatrix";
  generic_interface(assm);
  vector_interface(assm);

  py::class_<Sensor> so(m, "Sensor");
  generic_interface(so);
  so.def(
        "__init__",
        [](Sensor*                                         self,
           const AscendingGrid&                            freq_grid,
           const SensorPosLosVector&                       geometry,
           const std::vector<sensor::SparseStokvecMatrix>& elements,
           const SensorMetaInfo&                           meta) {
          new (self) Sensor{.freq_grid = freq_grid, .geometry = geometry, .elements = elements, .meta = meta};
          self->check();
        },
        "freq_grid"_a,
        "geometry"_a,
        "elements"_a,
        "meta"_a,
        "Construct a sensor from all its parts, which must agree, see :meth:`check`.  Its name is that of its meta.")
      .def_prop_rw(
          "name",
          [](const Sensor& self) { return self.name(); },
          [](Sensor& self, const String& name) { self.name() = name; },
          "The name of the sensor, the data name of its meta: a label for messages, not necessarily unique."
          "\n\n.. :class:`~pyarts3.arts.String`")
      .def_rw("freq_grid",
              &Sensor::freq_grid,
              "The frequencies at which the radiance is simulated.\n\n.. :class:`~pyarts3.arts.AscendingGrid`")
      .def_rw("geometry",
              &Sensor::geometry,
              "The positions and lines of sight at which the radiance is simulated, one agenda call each."
              "\n\n.. :class:`~pyarts3.arts.SensorPosLosVector`")
      .def_rw("elements",
              &Sensor::elements,
              "The observation elements, each a geometry by frequency sparse matrix of Stokes weights."
              "\n\n.. :class:`~pyarts3.arts.ArrayOfSparseStokvecMatrix`")
      .def_rw("meta",
              &Sensor::meta,
              "The structure of the elements, one value per element in their order, named as the sensor.\n\n"
              ".. :class:`~pyarts3.arts.SortedGriddedField1 | ~pyarts3.arts.SortedGriddedField2 | "
              "~pyarts3.arts.CameraGriddedField | ~pyarts3.arts.SortedGriddedField4`")
      .def("check",
           &Sensor::check,
           "Raise an error if the weights of an element do not have one row per geometry and one column per "
           "frequency, or if the meta does not describe one value per element.")
      .def("normalize",
           &Sensor::normalize,
           "pol"_a = Stokvec{1., 0., 0., 0.},
           R"(Scale the weights of each element by one common factor so that their summed I is the I of ``pol``.

When either has no I, the magnitudes of their summed Q, U, and V are made to agree instead.
The polarization of the weights is not changed: it is an error to give a polarized ``pol`` to
weights without polarization.  Set the polarization of a channel when building the sensor,
see :attr:`~pyarts3.arts.sensors.Channel.polarization`.)")
      .def(
          "reshape",
          [](const Sensor& self, const Vector& y) {
            if (y.size() != self.elements.size())
              throw std::invalid_argument(std::format(
                  "Sensor \"{}\" has {} elements, got {} values", self.name(), self.elements.size(), y.size()));

            SensorMetaInfo out = self.meta;
            std::visit(
                [&y](auto& gf) {
                  for (Size i = 0; i < y.size(); ++i) gf.data.elem_at(static_cast<Index>(i)) = y[i];
                },
                out);
            return out;
          },
          "y"_a,
          R"(The meta with its data set to ``y``, one value per element, e.g., this sensor's part of *measurement_vec*.

See :meth:`~pyarts3.arts.ArrayOfSensor.split` for the part of each sensor.)")
      .def("extract",
           &Sensor::extract,
           "spectral_rad"_a,
           R"(One value per element: its weights times the spectral radiance, summed.

``spectral_rad`` is the spectral radiance of each geometry (rows) at each frequency of
:attr:`freq_grid` (columns).)")
      .def(
          "extract",
          [](const Sensor& self, const StokvecVector& spectral_rad) {
            if (self.geometry.size() != 1)
              throw std::invalid_argument(
                  std::format("Sensor \"{}\" has {} geometries, so it needs the spectral radiance of each of them",
                              self.name(),
                              self.geometry.size()));
            return self.extract(StokvecMatrix{spectral_rad.view_as(1, spectral_rad.size())});
          },
          "spectral_rad"_a,
          R"(One value per element: its weights times the spectral radiance, summed.

``spectral_rad`` is the spectral radiance at each frequency of :attr:`freq_grid`, for a
sensor with one geometry.)")
      .def(
          "cost_report",
          [](const Sensor& self) { return cost_report({&self, 1}); },
          "A table of what simulating the sensor is expected to cost.");

  auto a0 =
      py::bind_vector<Array<SensorPosLosVector>, py::rv_policy::reference_internal>(m, "ArrayOfSensorPosLosVector");
  a0.doc() = "Array of SensorPosLosVector";
  generic_interface(a0);
  vector_interface(a0);

  auto a1 = py::bind_vector<ArrayOfSensor, py::rv_policy::reference_internal>(m, "ArrayOfSensor");
  generic_interface(a1);
  vector_interface(a1);

  a1.def(
      "split",
      [](const ArrayOfSensor& self, const Vector& y) {
        if (y.size() != measurement_size(self))
          throw std::invalid_argument(
              std::format("The sensors have {} elements, got {} values", measurement_size(self), y.size()));

        std::vector<Vector> out;
        out.reserve(self.size());
        Size start = 0;
        for (const auto& sensor : self) {
          out.emplace_back(y[Range(start, sensor.elements.size())]);
          start += sensor.elements.size();
        }
        return out;
      },
      "y"_a,
      R"(The part of ``y`` of each sensor, e.g., of *measurement_vec*.

The elements of the sensors, in order, are the elements of ``y``.  See
:meth:`~pyarts3.arts.Sensor.reshape` to give a part the structure of its sensor.)");

  a1.def(
      "report",
      [](const ArrayOfSensor& self) { return cost_report(self); },
      R"(A table of what simulating the sensors is expected to cost.

Each sensor makes one agenda call per geometry, over its whole frequency grid.
Per sensor, the table lists its elements, the shape of its meta, its agenda
calls, its frequencies, the spectral points those calls simulate, and its
non-zero weights, followed by the totals.  Geometries that more than one sensor
simulates are reported last, as each such sensor makes its own agenda call.)");

  auto sen = m.def_submodule("sensors", "Classes and functions to help create a sensor");

  auto spg  = py::class_<sensor::PixelGrid>(sen, "PixelGrid");
  spg.doc() = R"(The layout of an antenna's pixels.

An antenna makes one observation element per pixel for each channel.  Its pixels form a grid of rows
and columns, given by the angular offsets [deg] of their boresights from the antenna's boresight, and
its elements are ordered by row and then by column.  A kind with more than one pixel defines the
directions of its rows and columns.  An antenna that is a single beam has one pixel at no offset.)";
  generic_interface(spg);
  spg.def(
         "__init__",
         [](sensor::PixelGrid* self, const AscendingGrid& rows, const AscendingGrid& cols) {
           new (self) sensor::PixelGrid{.rows = rows, .cols = cols};
         },
         "rows"_a,
         "cols"_a,
         "Construct a pixel grid from the angular offsets of its rows and columns [deg].")
      .def_rw("rows",
              &sensor::PixelGrid::rows,
              "Angular offsets of the rows [deg].\n\n.. :class:`~pyarts3.arts.AscendingGrid`")
      .def_rw("cols",
              &sensor::PixelGrid::cols,
              "Angular offsets of the columns [deg].\n\n.. :class:`~pyarts3.arts.AscendingGrid`")
      .def("__len__", &sensor::PixelGrid::size, "The number of pixels.");

  py::class_<sensor::AntennaPatternField> sapf(sen, "AntennaPatternField");
  sapf.doc() = "A 2D gridded antenna response on local zenith and azimuth offsets.";
  gridded_data_interface(sapf);
  generic_interface(sapf);

  auto sgp  = py::class_<sensor::GriddedAntennaPattern>(sen, "GriddedAntennaPattern");
  sgp.doc() = "A gridded angular antenna response on local zenith and azimuth offsets.";
  sgp.def(py::init<>(), "Construct an empty gridded antenna pattern.")
      .def(
          "__init__",
          [](sensor::GriddedAntennaPattern* self, const sensor::AntennaPatternField& response) {
            new (self) sensor::GriddedAntennaPattern{};
            self->data = response;
          },
          "response"_a,
          "Construct a gridded antenna pattern from a local zenith/azimuth response field.")
      .def_rw("response",
              &sensor::GriddedAntennaPattern::data,
              py::rv_policy::reference_internal,
              "Local antenna response field.\n\n.. :class:`~pyarts3.arts.sensors.AntennaPatternField`");
  antenna_calls(sgp);
  generic_interface(sgp);

  auto spp  = py::class_<sensor::PencilBeamAntenna>(sen, "PencilBeamAntenna");
  spp.doc() = "A 1x1 antenna response centered on the bore line of sight.";
  spp.def(py::init<>(), "Construct a pencil-beam antenna that samples only the bore LOS.");
  antenna_calls(spp);
  generic_interface(spp);

  auto sga  = py::class_<sensor::GaussianAntenna>(sen, "GaussianAntenna");
  sga.doc() = "A Gaussian antenna response defined by angular offset from the bore LOS on a local zenith/azimuth grid.";
  sga.def(py::init<ZenGrid, Numeric, Size>(),
          "grid"_a,
          "std"_a,
          "azi_grid_size"_a = Size{2},
          "Construct a Gaussian antenna gain exp(-(offset / std)^2 / 2) of the angular offset [deg] from the bore LOS, "
          "sampled on the supplied zenith offsets and an azimuth grid of the provided size.");
  sga.def_rw("std",
             &sensor::GaussianAntenna::std,
             "Standard deviation of the gain [deg].\n\n.. :class:`~pyarts3.arts.Numeric`");
  sga.def_rw("zen_grid",
             &sensor::GaussianAntenna::zen_grid,
             "Local zenith grid of the antenna response.\n\n.. :class:`~pyarts3.arts.ZenGrid`");
  sga.def_rw("azi_grid",
             &sensor::GaussianAntenna::azi_grid,
             "Local azimuth grid of the antenna response.\n\n.. :class:`~pyarts3.arts.AziGrid`");
  antenna_calls(sga);
  generic_interface(sga);

  auto sgairy = py::class_<sensor::GaussianAiryAntenna>(sen, "GaussianAiryAntenna");
  sgairy.doc() =
      "A Gaussian antenna gain with the half-power width of the Airy pattern of a uniformly illuminated "
      "circular aperture, so its width scales with wavelength over aperture diameter.";
  sgairy.def(py::init<ZenGrid, Numeric, Size>(),
             "grid"_a,
             "aperture_diameter"_a,
             "azi_grid_size"_a = Size{2},
             "Construct a Gaussian Airy antenna on the supplied local grids."
             " The channel frequency grid must be strictly positive because the"
             " current builder path does not carry a separate reference frequency.");
  sgairy.def_rw("aperture_diameter",
                &sensor::GaussianAiryAntenna::aperture_diameter,
                "Aperture diameter [m].\n\n.. :class:`~pyarts3.arts.Numeric`");
  sgairy.def_rw("zen_grid",
                &sensor::GaussianAiryAntenna::zen_grid,
                "Local zenith grid of the antenna response.\n\n.. :class:`~pyarts3.arts.ZenGrid`");
  sgairy.def_rw("azi_grid",
                &sensor::GaussianAiryAntenna::azi_grid,
                "Local azimuth grid of the antenna response.\n\n.. :class:`~pyarts3.arts.AziGrid`");
  sgairy.def(
      "std",
      [](const sensor::GaussianAiryAntenna& self, py::object& f_) {
        return vectorize([&self](Numeric f) { return self.std(f); }, f_);
      },
      "Standard deviation of the antenna response at a given frequency.");
  antenna_calls(sgairy);
  generic_interface(sgairy);
  static_assert(arts_xml_ioable<sensor::GaussianAiryAntenna>);

  auto scam  = py::class_<sensor::CameraAntenna>(sen, "CameraAntenna");
  scam.doc() = R"(A camera: a grid of pixels in CCD coordinates, each with the same beam.

``grid`` holds the angular offsets [deg] of the pixels from the optical axis, the boresight, in CCD
coordinates: rows towards increasing zenith and, because the lens inverts the image, columns towards
decreasing azimuth.  A pixel at row offset ``r`` and column offset ``c`` looks along
``g + tan(r) e_zen - tan(c) e_azi``, where ``g`` is the boresight and ``e_zen`` and ``e_azi`` point
towards increasing zenith and azimuth there.  Each pixel's polarization response is the boresight
response carried along the great circle to the pixel, and ``beam`` is centred on each pixel.)";
  generic_interface(scam);
  scam.def(
          "__init__",
          [](sensor::CameraAntenna* self, const sensor::PixelGrid& grid, const sensor::BeamPattern& beam) {
            new (self) sensor::CameraAntenna{.grid = grid, .beam = beam};
          },
          "grid"_a,
          "beam"_a = sensor::BeamPattern{sensor::PencilBeamAntenna{}},
          "Construct a camera from its pixel grid and the beam of each pixel (defaults to pencil-beam).")
      .def_static("thin_lens",
                  &sensor::CameraAntenna::thin_lens,
                  "n_rows"_a,
                  "n_cols"_a,
                  "ccd_height"_a,
                  "ccd_width"_a,
                  "focal_length"_a,
                  "focus_distance"_a,
                  "beam"_a = sensor::BeamPattern{sensor::PencilBeamAntenna{}},
                  R"(A thin-lens camera.

``n_rows`` by ``n_cols`` pixels on a CCD of ``ccd_height`` by ``ccd_width`` [m], behind a lens of
``focal_length`` [m] focused at ``focus_distance`` [m].  The CCD sits at the image distance of the
thin-lens equation, and each pixel's offset is the angle of its centre seen from the lens.)")
      .def_rw("grid",
              &sensor::CameraAntenna::grid,
              "The pixels in CCD coordinates.\n\n.. :class:`~pyarts3.arts.sensors.PixelGrid`")
      .def_rw("beam",
              &sensor::CameraAntenna::beam,
              "The beam of each pixel.\n\n"
              ".. :class:`~pyarts3.arts.sensors.PencilBeamAntenna | ~pyarts3.arts.sensors.GaussianAntenna | "
              "~pyarts3.arts.sensors.GaussianAiryAntenna | ~pyarts3.arts.sensors.GriddedAntennaPattern`");
  antenna_calls(scam);

  auto spol  = py::class_<sensor::Polarization>(sen, "Polarization");
  spol.doc() = R"(The polarization response of a receiver channel.

``state`` is the Stokes vector of weights at the boresight, in the polarization basis of the
received radiation (that of its propagation direction, opposite to the line of sight), so that
``"Iv"`` receives vertical and ``"Ih"`` horizontal polarization.

``scan_rotation`` turns the linear part of the response about the ray by that many degrees per
degree of boresight nadir angle (180 minus the zenith angle of the line of sight).  It is 0 when
the response does not depend on the pointing and +1 or -1 for a cross-track scanner with a
rotating reflector, whose polarization is ``state`` at nadir.  The sign belongs to the instrument
and the side of the scan and only matters when the observed U is not zero.)";
  spol.def(py::init<>(), "Construct an unpolarized response.")
      .def(
          "__init__",
          [](sensor::Polarization* self, const Stokvec& state, Numeric scan_rotation) {
            new (self) sensor::Polarization{.state = state, .scan_rotation = scan_rotation};
          },
          "state"_a,
          "scan_rotation"_a = 0.0,
          "Construct a response from its boresight state and scan rotation.")
      .def(
          "__init__",
          [](sensor::Polarization* self, const String& state) {
            new (self) sensor::Polarization{.state = rtepack::to_stokvec(to<PolarizationChoice>(state))};
          },
          "state"_a,
          "Construct a response that does not depend on the pointing from a named polarization, e.g. ``\"Iv\"``.")
      .def_rw(
          "state", &sensor::Polarization::state, "The response at the boresight.\n\n.. :class:`~pyarts3.arts.Stokvec`")
      .def_rw("scan_rotation",
              &sensor::Polarization::scan_rotation,
              "Turn of the linear polarization per degree of boresight nadir angle.\n\n.. :class:`float`")
      .def("at",
           &sensor::Polarization::at,
           "scan_angle"_a,
           "The response at a signed scan angle [deg], positive to the right of the flight direction.");
  generic_interface(spol);
  py::implicitly_convertible<Stokvec, sensor::Polarization>();
  py::implicitly_convertible<String, sensor::Polarization>();

  auto sch  = py::class_<sensor::Channel>(sen, "Channel");
  sch.doc() = R"(A receiver channel: its frequency response and its polarization response.

The frequency response is a sum of pieces.  Each piece is linear between its frequencies and zero
outside them; a piece with a single frequency is a delta function of that weight.  The frequencies
of a piece are where the radiance is sampled, so their density is chosen with the channel.  The
measurement is the integral of the response times the radiance over the integral of the response,
see :meth:`quadrature`.)";
  generic_interface(sch);
  sch.def(
         "__init__",
         [](sensor::Channel* self, const SortedGriddedField1& piece, const sensor::Polarization& polarization) {
           new (self) sensor::Channel{.response = {piece}, .polarization = polarization};
         },
         "response"_a,
         "polarization"_a = sensor::Polarization{},
         "Construct a channel from one response piece and its polarization response.")
      .def(
          "__init__",
          [](sensor::Channel*                        self,
             const std::vector<SortedGriddedField1>& pieces,
             const sensor::Polarization&             polarization) {
            new (self) sensor::Channel{.response = pieces, .polarization = polarization};
          },
          "response"_a,
          "polarization"_a = sensor::Polarization{},
          "Construct a channel from its response pieces and its polarization response.");
  py::implicitly_convertible<SortedGriddedField1, sensor::Channel>();
  sch.def_rw("polarization",
             &sensor::Channel::polarization,
             "What the channel receives of the Stokes vector at the boresight."
             "\n\n.. :class:`~pyarts3.arts.sensors.Polarization`")
      .def_rw("response",
              &sensor::Channel::response,
              "The pieces of the frequency response."
              "\n\n.. :class:`list[~pyarts3.arts.SortedGriddedField1]`")
      .def("quadrature",
           &sensor::Channel::quadrature,
           R"(The integration weights of the response on the frequencies of its pieces, not normalized.

The trapezoid rule within each piece, delta weights as they are, summed where pieces share a frequency.)");

  auto sbox  = py::class_<sensor::BoxChannel, sensor::Channel>(sen, "BoxChannel");
  sbox.doc() = "A channel with a uniform response across a finite frequency interval, sampled at its points.";
  sbox.def(py::init<Numeric, Numeric, Size>(),
           "lower"_a,
           "upper"_a,
           "n"_a,
           "Construct a box channel on ``[lower, upper]`` with ``n`` points.")
      .def(py::init<Numeric, Size>(),
           "half_width"_a,
           "n"_a,
           "Construct a symmetric box channel on ``[-half_width, half_width]``.")
      .def(
          py::init<AscendingGrid>(), "freq_grid"_a, "Construct a box channel directly from a relative frequency grid.");
  generic_interface(sbox);

  auto sdirac  = py::class_<sensor::DiracChannel, sensor::Channel>(sen, "DiracChannel");
  sdirac.doc() = "A single-frequency relative channel.";
  sdirac.def(py::init<Numeric>(), "frequency"_a, "Construct a Dirac channel at one relative frequency.")
      .def(py::init<>(), "Construct a Dirac channel at 0 relative frequency.");
  generic_interface(sdirac);

  auto sgauss  = py::class_<sensor::GaussianChannel, sensor::Channel>(sen, "GaussianChannel");
  sgauss.doc() = "A Gaussian relative channel response.";
  sgauss
      .def(py::init<AscendingGrid, Numeric, Numeric>(),
           "freq_grid"_a,
           "center"_a,
           "std"_a,
           "Construct a Gaussian channel on a custom relative frequency grid.")
      .def(py::init<Numeric, Numeric, Size, Size>(),
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

  auto sspec  = py::class_<sensor::Spectrometer>(sen, "Spectrometer");
  sspec.doc() = "A spectrometer made from channels arranged on one shared relative-frequency grid.";
  sspec
      .def(py::init<const sensor::Channel&, const AscendingGrid&>(),
           "base_channel"_a,
           "freq_offsets"_a,
           "Construct a spectrometer by shifting one base channel across the supplied relative-frequency offsets.")
      .def(py::init_implicit<const AscendingGrid&>(),
           "freq_offsets"_a,
           "Construct a spectrometer with Dirac channels across the supplied relative-frequency offsets.")
      .def(py::init<std::vector<sensor::Channel>>(), "channels"_a, "Construct a spectrometer from explicit channels.")
      .def_prop_ro(
          "channels",
          [](const sensor::Spectrometer& self) -> const std::vector<sensor::Channel>& { return self.channels; },
          py::rv_policy::reference_internal,
          "Spectrometer channels.\n\n.. :class:`list[~pyarts3.arts.sensors.Channel]`")
      .def("__len__", [](const sensor::Spectrometer& self) { return self.channels.size(); })
      .def(
          "__getitem__",
          [](const sensor::Spectrometer& self, Size i) -> const sensor::Channel& { return self.channels.at(i); },
          py::rv_policy::reference_internal)
      .def(
          "__iter__",
          [](const sensor::Spectrometer& self) {
            return py::make_iterator(py::type<sensor::Spectrometer>(),
                                     "sensor-spectrometer-iterator",
                                     self.channels.begin(),
                                     self.channels.end());
          },
          py::rv_policy::reference_internal,
          "Iterate over the spectrometer channels.");
  generic_interface(sspec);

  auto ssb = py::class_<SensorBuilder>(m, "SensorBuilder");
  ssb.def(py::init<>(), "Construct an empty sensor builder.")
      .def(
          "__init__",
          [](sensor::Builder*                    self,
             const std::vector<sensor::Channel>& channels,
             const sensor::AntennaPattern&       antenna) {
            new (self) sensor::Builder{.channels = channels, .antenna = antenna};
          },
          "channels"_a,
          "antenna"_a = sensor::AntennaPattern{sensor::PencilBeamAntenna{}},
          "Construct a sensor builder from channels and one antenna pattern (defaults to pencil-beam).")
      .def(
          "__init__",
          [](sensor::Builder* self, const sensor::Spectrometer& spectrometer, const sensor::AntennaPattern& antenna) {
            new (self) sensor::Builder{sensor::Builder::from(spectrometer, antenna)};
          },
          "spectrometer"_a,
          "antenna"_a = sensor::AntennaPattern{sensor::PencilBeamAntenna{}},
          "Construct a sensor builder from a spectrometer and one antenna pattern (defaults to pencil-beam).")
      .def(
          "__init__",
          [](sensor::Builder*                        self,
             const sensor::AntennaPattern&           antenna,
             const sensor::Spectrometer&             spectrometer,
             const sensor::HeterodyneFrequencyRange& backend) {
            new (self) sensor::Builder{sensor::Builder::from(spectrometer, backend, antenna)};
          },
          "antenna"_a,
          "spectrometer"_a,
          "backend"_a,
          "Construct a sensor builder from an antenna, a spectrometer, and a backend frequency response.")
      .def_rw("channels",
              &sensor::Builder::channels,
              "Spectrometer channels.\n\n.. :class:`list[~pyarts3.arts.sensors.Channel]`")
      .def_rw("antenna",
              &sensor::Builder::antenna,
              "Angular antenna response and pixel layout.\n\n"
              ".. :class:`~pyarts3.arts.sensors.PencilBeamAntenna | ~pyarts3.arts.sensors.GaussianAntenna | "
              "~pyarts3.arts.sensors.GaussianAiryAntenna | ~pyarts3.arts.sensors.GriddedAntennaPattern | "
              "~pyarts3.arts.sensors.CameraAntenna`")
      .def(
          "__call__",
          [](const sensor::Builder&      self,
             const std::vector<Vector3>& pos,
             const std::vector<Vector2>& los,
             const Vector2&              ell,
             const String&               name) {
            if (pos.size() != los.size()) {
              throw std::invalid_argument(std::format(
                  "Sensor builder expects matching position and LOS counts. Got {} positions and {} LOS values.",
                  pos.size(),
                  los.size()));
            }

            return self(pos, los, ell, name);
          },
          "pos"_a,
          "los"_a,
          "ell"_a,
          "name"_a = String{},
          R"(Build a sensor from paired position and bore-LOS sequences.

Each ``pos[i]`` with ``los[i]`` is one placement of the instrument.  The elements
of the returned :class:`~pyarts3.arts.Sensor` are ordered by placement, then pixel
of the antenna, then channel, and its meta has these axes.  The sensor is named ``name``, or
after its number of channels if ``name`` is empty.)")
      .def(
          "__call__",
          [](const sensor::Builder& self, const Vector3 pos, const Vector2 los, const Vector2 ell, const String& name) {
            return self(std::span{&pos, 1}, std::span{&los, 1}, ell, name);
          },
          "pos"_a,
          "los"_a,
          "ell"_a,
          "name"_a = String{},
          R"(Build a sensor from one position and bore LOS.)")
      .def(
          "__call__",
          [](const sensor::Builder& self, const SensorPosLosVector& poslos, const Vector2& ell, const String& name) {
            std::vector<Vector3> pos(poslos.size());
            std::vector<Vector2> los(poslos.size());

            for (Size i = 0; i < poslos.size(); ++i) {
              pos[i] = poslos[i].pos;
              los[i] = poslos[i].los;
            }

            return self(pos, los, ell, name);
          },
          "poslos"_a,
          "ell"_a,
          "name"_a = String{},
          "Build a sensor from paired position/LOS samples, one placement each.");
  ssb.def(
      "cross_track",
      [](const sensor::Builder&      self,
         const std::vector<Vector3>& pos,
         const std::vector<Numeric>& heading,
         const Vector&               scan_angles,
         const Vector2&              ell,
         const String&               name) { return self.cross_track(pos, heading, scan_angles, ell, name); },
      "pos"_a,
      "heading"_a,
      "scan_angles"_a,
      "ell"_a,
      "name"_a = String{},
      R"(Build a sensor for cross-track scans.

Every ``pos[i]`` scans all ``scan_angles`` [deg] across its flight direction ``heading[i]`` [deg,
clockwise from north].  A scan angle is the nadir angle of the boresight, positive to the right of
the flight direction, so its bore line of sight is ``[180 - |s|, heading + 90]`` for ``s >= 0`` and
``[180 - |s|, heading - 90]`` for ``s < 0``.  Each channel's polarization is taken at ``s``, see
:meth:`~pyarts3.arts.sensors.Polarization.at`.

Each position and scan angle is one placement of the instrument: the elements of the returned
:class:`~pyarts3.arts.Sensor` are ordered by position, then scan angle, then channel.  The antenna
must have one pixel.  The sensor is named ``name``, or after its number of channels if
``name`` is empty.)");
  generic_interface(ssb);

  auto shdfr = py::class_<sensor::HeterodyneFrequencyRange>(sen, "HeterodyneFrequencyRange");
  shdfr.def(py::init<>(),
            R"(Construct an empty staged heterodyne response.

  Stages can then be applied in sequence using :func:`lowpass`, :func:`highpass`,
  :func:`bandpass`, :func:`filter`, and :func:`mix`.)");
  shdfr.def(py::init<Numeric, Vector2>(),
            "lo"_a,
            "bandpass"_a = Vector2{0.0, std::numeric_limits<Numeric>::infinity()},
            R"(Construct a heterodyne response from one ideal bandpass and one LO stage.

  This is shorthand for creating an empty object, applying :func:`bandpass`, and
  then applying :func:`mix`.)");
  shdfr.def(
      "__init__",
      [](sensor::HeterodyneFrequencyRange* v,
         const std::vector<Numeric>&       clock_frequencies,
         const std::vector<Vector2>&       bandpasses) {
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
      .def("highpass",
           &sensor::HeterodyneFrequencyRange::apply_highpass,
           "lower"_a,
           "Apply an ideal highpass filter on the current local frequency axis.")
      .def(
          "bandpass",
          [](sensor::HeterodyneFrequencyRange& self, const Vector2& bandpass) { self.apply_bandpass(bandpass); },
          "bandpass"_a,
          "Apply an ideal bandpass filter on the current local frequency axis.")
      .def(
          "filter",
          [](sensor::HeterodyneFrequencyRange& self, const SortedGriddedField1& bandpass_filter) {
            self.apply_bandpass(bandpass_filter);
          },
          "bandpass_filter"_a,
          R"(Apply a weighted bandpass filter on the current local frequency axis.

  The filter weights are interpreted on the filter's relative frequency grid and
  are zero outside that grid.)")
      .def("mix", &sensor::HeterodyneFrequencyRange::apply_mixer, "lo"_a, "Apply one heterodyne LO mixing stage.")
      .def(
          "local_response",
          [](const sensor::HeterodyneFrequencyRange& self, const Vector& f, Size path_index) {
            return self.local_response(f, path_index);
          },
          "f"_a,
          "path_index"_a = 0,
          "Evaluate one path response on the current local frequency axis.")
      .def(
          "global_response",
          [](const sensor::HeterodyneFrequencyRange& self, const Vector& f, Size path_index) {
            return self.global_response(f, path_index);
          },
          "f"_a,
          "path_index"_a = 0,
          "Evaluate one path response on the original real-frequency axis.")
      .def(
          "channel_response",
          [](const sensor::HeterodyneFrequencyRange& self, const sensor::Channel& channel) {
            return sensor::make_bandpass_channels(self, std::span<const sensor::Channel>{&channel, 1}).front();
          },
          "channel"_a,
          R"(The channel as seen at the real frequencies.

Each response piece of the channel becomes one piece per mixer path that passes it, so the
response is the sum over the sidebands.  The polarization is kept.)")
      .def(
          "channel_responses",
          [](const sensor::HeterodyneFrequencyRange& self, const std::vector<sensor::Channel>& channels) {
            return sensor::make_bandpass_channels(self,
                                                  std::span<const sensor::Channel>{channels.data(), channels.size()});
          },
          "channels"_a,
          "The channels as seen at the real frequencies, see :meth:`channel_response`.")
      .def(
          "channel_responses",
          [](const sensor::HeterodyneFrequencyRange& self, const sensor::Spectrometer& spectrometer) {
            return sensor::make_bandpass_channels(self, spectrometer);
          },
          "spectrometer"_a,
          "The spectrometer's channels as seen at the real frequencies, see :meth:`channel_response`.");
  shdfr.doc() = "A staged heterodyne mixer and filter response builder.";
  generic_interface(shdfr);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize sensors\n{}", e.what()));
}
}  // namespace Python
