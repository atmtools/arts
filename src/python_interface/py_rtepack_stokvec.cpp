#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "py_rtepack_helpers.h"

namespace Python {
void py_rtepack_stokvec(py::module_ &m) {
  py::class_<Stokvec> sv(m, "Stokvec");
  sv.def(py::init_implicit<Numeric>())
      .def("__init__", [](Stokvec *s, const PolarizationChoice p) { new (s) Stokvec{rtepack::to_stokvec(p)}; })
      .def("__init__",
           [](Stokvec *s, const String &p) { new (s) Stokvec{rtepack::to_stokvec(to<PolarizationChoice>(p))}; })
      .def(py::init_implicit<std::array<Numeric, 4>>())
      .def_static(
          "linpol",
          [](const Numeric angle) {
            return Stokvec{1.0, Conversion::cosd(2.0 * angle), Conversion::sind(2.0 * angle), 0.0};
          },
          "Returns [1.0, cos(2*angle), sin(2*angle), 0.0], the linear polarization vector for a given angle",
          "angle"_a)
      .def_static(
          "cirpol",
          [](const Numeric angle) { return Stokvec{1.0, 0.0, 0.0, Conversion::sind(angle)}; },
          "Returns [1.0, 0.0, 0.0, sin(angle)], the circular polarization vector for a given phase delay angle",
          "angle"_a)
      .def(
          "__array__",
          [](Stokvec   &v,
             py::object dtype,
             py::object copy) -> std::variant<py::ndarray<py::numpy, Numeric, py::shape<4>, py::c_contig>, py::object> {
            std::array<size_t, 1> shape = {4};
            auto                  np    = py::module_::import_("numpy");
            auto                  x     = py::ndarray<py::numpy, Numeric, py::shape<4>, py::c_contig>(
                v.data.data(), 1, shape.data(), py::cast(&v));

            if (not dtype.is_none()) { return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy); }

            if (copy.is_none() or not py::bool_(copy)) { return x.cast(py::rv_policy::automatic_reference); }
            return x.cast(py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object &x) { return x.attr("__array__")(); },
          [](Stokvec &x, Stokvec &y) { x = y; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");
  common_ndarray(sv);
  generic_interface(sv);
  py::implicitly_convertible<PolarizationChoice, Stokvec>();
  py::implicitly_convertible<String, Stokvec>();

  auto asv = py::bind_vector<std::vector<Stokvec>>(m, "ArrayOfStokvec");
  generic_interface(asv);
  vector_interface(asv);
  asv.doc() = "A list of :class:`~pyarts3.arts.Stokvec`";

  py::class_<StokvecVector> vsv(m, "StokvecVector");
  vsv.def(
      "__init__",
      [](StokvecVector *v, const std::vector<Numeric> &a) {
        new (v) StokvecVector(a.size());
        stdr::transform(a, v->begin(), [](const Numeric &x) { return Stokvec{x}; });
      },
      "a"_a);
  vsv.def(
      "__init__",
      [](StokvecVector *v, const std::vector<Stokvec> &a) {
        new (v) StokvecVector(a.size());
        stdr::transform(a, v->begin(), [](const Stokvec &x) { return x; });
      },
      "a"_a);
  rtepack_array<Stokvec, 1, 4>(vsv);
  generic_interface(vsv);
  py::implicitly_convertible<std::vector<Numeric>, StokvecVector>();
  py::implicitly_convertible<std::vector<Stokvec>, StokvecVector>();

  py::class_<StokvecMatrix> msv(m, "StokvecMatrix");
  rtepack_array<Stokvec, 2, 4>(msv);
  generic_interface(msv);

  py::class_<StokvecTensor3> t3sv(m, "StokvecTensor3");
  rtepack_array<Stokvec, 3, 4>(t3sv);
  generic_interface(t3sv);

  py::class_<StokvecTensor4> t4sv(m, "StokvecTensor4");
  rtepack_array<Stokvec, 4, 4>(t4sv);
  generic_interface(t4sv);

  py::class_<StokvecTensor5> t5sv(m, "StokvecTensor5");
  rtepack_array<Stokvec, 5, 4>(t5sv);
  generic_interface(t5sv);

  py::class_<StokvecTensor6> t6sv(m, "StokvecTensor6");
  rtepack_array<Stokvec, 6, 4>(t6sv);
  generic_interface(t6sv);

  py::class_<CompactPlanarStokvec> csv(m, "CompactPlanarStokvec");
  csv.def(py::init<Numeric, Numeric>(), "I"_a, "Q"_a = 0.0)
      .def(py::init_implicit<std::array<Numeric, 2>>())
      .def("expand", &CompactPlanarStokvec::expand, "The full :class:`~pyarts3.arts.Stokvec`")
      .def(
          "__array__",
          [](CompactPlanarStokvec &x,
             py::object            dtype,
             py::object copy) -> std::variant<py::ndarray<py::numpy, Numeric, py::shape<2>, py::c_contig>, py::object> {
            std::array<size_t, 1> shape = {2};
            auto                  np    = py::module_::import_("numpy");
            auto                  w     = py::ndarray<py::numpy, Numeric, py::shape<2>, py::c_contig>(
                x.data.data(), 1, shape.data(), py::cast(&x));

            if (not dtype.is_none()) { return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy); }

            if (copy.is_none() or not py::bool_(copy)) { return w.cast(py::rv_policy::automatic_reference); }
            return w.cast(py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object, [I, Q].")
      .def_prop_rw(
          "value",
          [](py::object &x) { return x.attr("__array__")(); },
          [](CompactPlanarStokvec &x, CompactPlanarStokvec &y) { x = y; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");
  for (const auto &[name, i] : std::initializer_list<std::pair<const char *, Size>>{{"I", 0}, {"Q", 1}}) {
    csv.def_prop_rw(
        name,
        [i](const CompactPlanarStokvec &x) { return x.data[i]; },
        [i](CompactPlanarStokvec &x, Numeric v) { x.data[i] = v; },
        "An element of the Stokes vector\n\n.. :class:`float`");
  }
  common_ndarray(csv);
  generic_interface(csv);
  csv.doc() =
      R"(The Stokes vector of a medium with a plane of symmetry through the propagation direction, in compact form.

In the basis of that plane it is :math:`[I, Q, 0, 0]`, and only :math:`[I, Q]`, the elements of the same
names of a :class:`~pyarts3.arts.Stokvec`, are stored.  For the absorption vector of azimuthally randomly
oriented particles they are :math:`a_1` and :math:`a_2`; for totally randomly oriented particles :math:`Q = 0`.
)";

  py::class_<CompactPlanarStokvecVector> vcsv(m, "CompactPlanarStokvecVector");
  rtepack_array<CompactPlanarStokvec, 1, 2>(vcsv);
  generic_interface(vcsv);
  vcsv.doc() = "A vector of :class:`~pyarts3.arts.CompactPlanarStokvec`; as an array [n, 2]";

  py::class_<CompactPlanarStokvecMatrix> mcsv(m, "CompactPlanarStokvecMatrix");
  rtepack_array<CompactPlanarStokvec, 2, 2>(mcsv);
  generic_interface(mcsv);
  mcsv.doc() = "A matrix of :class:`~pyarts3.arts.CompactPlanarStokvec`; as an array [n, m, 2]";

  py::class_<CompactPlanarStokvecTensor3> tcsv(m, "CompactPlanarStokvecTensor3");
  rtepack_array<CompactPlanarStokvec, 3, 2>(tcsv);
  generic_interface(tcsv);
  tcsv.doc() = "A 3-tensor of :class:`~pyarts3.arts.CompactPlanarStokvec`; as an array [n, m, k, 2]";
}
}  // namespace Python
