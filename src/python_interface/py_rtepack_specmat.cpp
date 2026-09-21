#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "py_rtepack_helpers.h"

namespace Python {
void py_rtepack_specmat(py::module_ &m) {
  py::class_<Specmat> cmm(m, "Specmat");
  cmm.def(py::init_implicit<Complex>())
      .def(py::init_implicit<std::array<Complex, 16>>())
      .def(
          "__array__",
          [](Specmat &x, py::object dtype, py::object copy)
              -> std::variant<py::ndarray<py::numpy, Complex, py::shape<4, 4>, py::c_contig>, py::object> {
            std::array<size_t, 2> shape = {4, 4};
            auto                  np    = py::module_::import_("numpy");
            auto                  w     = py::ndarray<py::numpy, Complex, py::shape<4, 4>, py::c_contig>(
                x.data.data(), 2, shape.data(), py::cast(&x));

            if (not dtype.is_none()) { return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy); }

            if (copy.is_none() or not py::bool_(copy)) { return w.cast(py::rv_policy::automatic_reference); }
            return w.cast(py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object &x) { return x.attr("__array__")(); },
          [](Specmat &x, Specmat &y) { x = y; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");
  common_ndarray(cmm);
  generic_interface(cmm);

  auto asp  = py::bind_vector<std::vector<Specmat>, py::rv_policy::reference_internal>(m, "ArrayOfSpecmat");
  asp.doc() = "A list of :class:`~pyarts3.arts.Specmat`";
  vector_interface(asp);
  generic_interface(asp);

  py::class_<SpecmatVector> vcmm(m, "SpecmatVector");
  vcmm.def(
      "__init__",
      [](SpecmatVector *v, const std::vector<Numeric> &a) {
        new (v) SpecmatVector(a.size());
        stdr::transform(a, v->begin(), [](const Numeric &x) { return Specmat{x}; });
      },
      "a"_a);
  vcmm.def(
      "__init__",
      [](SpecmatVector *v, const std::vector<Specmat> &a) {
        new (v) SpecmatVector(a.size());
        stdr::transform(a, v->begin(), [](const Specmat &x) { return x; });
      },
      "a"_a);
  rtepack_array<Specmat, 1, 4, 4>(vcmm);
  vcmm.doc() = "A vector of :class:`~pyarts3.arts.Specmat`";
  generic_interface(vcmm);
  py::implicitly_convertible<std::vector<Numeric>, SpecmatVector>();
  py::implicitly_convertible<std::vector<Specmat>, SpecmatVector>();

  py::class_<SpecmatMatrix> mcmm(m, "SpecmatMatrix");
  rtepack_array<Specmat, 2, 4, 4>(mcmm);
  generic_interface(mcmm);

  py::class_<SpecmatTensor3> cmt3(m, "SpecmatTensor3");
  rtepack_array<Specmat, 3, 4, 4>(cmt3);
  cmt3.doc() = "A 3-tensor of :class:`~pyarts3.arts.Specmat`";
  generic_interface(cmt3);
}
}  // namespace Python
