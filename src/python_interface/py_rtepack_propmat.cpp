#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "py_rtepack_helpers.h"

namespace Python {
void py_rtepack_propmat(py::module_ &m) {
  py::class_<Propmat> pm(m, "Propmat");
  pm.def(py::init_implicit<Numeric>())
      .def(py::init_implicit<std::array<Numeric, 7>>())
      .def(
          "__array__",
          [](Propmat   &x,
             py::object dtype,
             py::object copy) -> std::variant<py::ndarray<py::numpy, Numeric, py::shape<7>, py::c_contig>, py::object> {
            std::array<size_t, 1> shape = {7};
            auto                  np    = py::module_::import_("numpy");
            auto                  w     = py::ndarray<py::numpy, Numeric, py::shape<7>, py::c_contig>(
                x.data.data(), 1, shape.data(), py::cast(&x));

            if (not dtype.is_none()) { return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy); }

            if (copy.is_none() or not py::bool_(copy)) { return w.cast(py::rv_policy::automatic_reference); }
            return w.cast(py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def(
          "as_matrix", [](Propmat &x) { return to_matrix(x); }, "Returns the Propmat as a matrix.")
      .def_prop_rw(
          "value",
          [](py::object &x) { return x.attr("__array__")(); },
          [](Propmat &x, Propmat &y) { x = y; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`")
      .def(
          "exp",
          [](const Propmat &k, Numeric r) { return exp(k, -r); },
          "r"_a = -1.0,
          "Returns the matrix exponential of the propagation matrix scaled by r.)")
      .def("inv", [](const Propmat &k) { return inv(k); }, "Returns the inverse of the propagation matrix.");

  common_ndarray(pm);
  generic_interface(pm);

  auto apm = py::bind_vector<std::vector<Propmat>>(m, "ArrayOfPropmat");
  generic_interface(apm);
  vector_interface(apm);
  apm.doc() = "A list of :class:`~pyarts3.arts.Propmat`";

  py::class_<PropmatVector> vpm(m, "PropmatVector");
  vpm.def(
      "__init__",
      [](PropmatVector *v, const std::vector<Numeric> &a) {
        new (v) PropmatVector(a.size());
        stdr::transform(a, v->begin(), [](const Numeric &x) { return Propmat{x}; });
      },
      "a"_a);
  vpm.def(
      "__init__",
      [](PropmatVector *v, const std::vector<Propmat> &a) {
        new (v) PropmatVector(a.size());
        stdr::transform(a, v->begin(), [](const Propmat &x) { return x; });
      },
      "a"_a);
  rtepack_array<Propmat, 1, 7>(vpm);
  generic_interface(vpm);
  py::implicitly_convertible<std::vector<Numeric>, PropmatVector>();
  py::implicitly_convertible<std::vector<Propmat>, PropmatVector>();

  py::class_<PropmatMatrix> mpm(m, "PropmatMatrix");
  rtepack_array<Propmat, 2, 7>(mpm);
  generic_interface(mpm);
}
}  // namespace Python
