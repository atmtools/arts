#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "py_rtepack_helpers.h"

namespace Python {
void py_rtepack_muelmat(py::module_ &m) {
  py::class_<Muelmat> mm(m, "Muelmat");
  mm.def(py::init_implicit<Numeric>())
      .def(py::init_implicit<std::array<Numeric, 16>>())
      .def("is_polarized", &Muelmat::is_polarized, "Check if the Mueller matrix represents a polarized state.")
      .def(
          "__array__",
          [](Muelmat &x, py::object dtype, py::object copy)
              -> std::variant<py::ndarray<py::numpy, Numeric, py::shape<4, 4>, py::c_contig>, py::object> {
            std::array<size_t, 2> shape = {4, 4};
            auto                  np    = py::module_::import_("numpy");
            auto                  w     = py::ndarray<py::numpy, Numeric, py::shape<4, 4>, py::c_contig>(
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
          [](Muelmat &x, Muelmat &y) { x = y; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");
  common_ndarray(mm);
  generic_interface(mm);

  auto amm = py::bind_vector<std::vector<Muelmat>>(m, "ArrayOfMuelmat");
  generic_interface(amm);
  vector_interface(amm);
  amm.doc() = "A list of :class:`~pyarts3.arts.Muelmat`";

  py::class_<MuelmatVector> vmm(m, "MuelmatVector");
  vmm.def(
      "is_polarized",
      [](const MuelmatVector &m) { return stdr::any_of(m, [](const Muelmat &mm) { return mm.is_polarized(); }); },
      "Check if the Mueller matrix represents a polarized state.");
  vmm.def(
      "__init__",
      [](MuelmatVector *v, const std::vector<Numeric> &a) {
        new (v) MuelmatVector(a.size());
        stdr::transform(a, v->begin(), [](const Numeric &x) { return Muelmat{x}; });
      },
      "a"_a);
  vmm.def(
      "__init__",
      [](MuelmatVector *v, const std::vector<Muelmat> &a) {
        new (v) MuelmatVector(a.size());
        stdr::transform(a, v->begin(), [](const Muelmat &x) { return x; });
      },
      "a"_a);
  rtepack_array<Muelmat, 1, 4, 4>(vmm);
  generic_interface(vmm);
  py::implicitly_convertible<std::vector<Numeric>, MuelmatVector>();
  py::implicitly_convertible<std::vector<Muelmat>, MuelmatVector>();

  py::class_<MuelmatMatrix> mmm(m, "MuelmatMatrix");
  rtepack_array<Muelmat, 2, 4, 4>(mmm);
  generic_interface(mmm);

  py::class_<MuelmatTensor3> mt3(m, "MuelmatTensor3");
  rtepack_array<Muelmat, 3, 4, 4>(mt3);
  generic_interface(mt3);

  py::class_<MuelmatTensor4> mt4(m, "MuelmatTensor4");
  mt4.doc() = "A 4-tensor of :class:`~pyarts3.arts.Muelmat`";
  rtepack_array<Muelmat, 4, 4, 4>(mt4);
  generic_interface(mt4);

  py::class_<MuelmatTensor5> mt5(m, "MuelmatTensor5");
  mt5.doc() = "A 5-tensor of :class:`~pyarts3.arts.Muelmat`";
  rtepack_array<Muelmat, 5, 4, 4>(mt5);
  generic_interface(mt5);
}
}  // namespace Python
