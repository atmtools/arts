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
}
}  // namespace Python
