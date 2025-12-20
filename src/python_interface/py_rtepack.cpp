#include <debug.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>
#include <rtepack.h>

#include <algorithm>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"

namespace Python {
template <typename T, Index M, size_t... N>
void rtepack_array(py::class_<matpack::data_t<T, M>> &c) {
  using U = T::value_type;

  c.def(
      "__init__",
      [](matpack::data_t<T, M> *y, const matpack::data_t<U, M> &x) {
        new (y) matpack::data_t<T, M>(x.shape());
        std::transform(
            x.elem_begin(), x.elem_end(), y->elem_begin(), [](const U &z) {
              return T(z);
            });
      },
      "x"_a);
  py::implicitly_convertible<matpack::data_t<U, M>, matpack::data_t<T, M>>();
  c.def(
      "__init__",
      [](matpack::data_t<T, M> *y, const matpack::data_t<U, M + 1> &x) {
        const Size sz = T{}.size();

        if (x.size() != 0 and static_cast<Size>(x.shape().back()) != sz) {
          throw std::invalid_argument(std::format(
              "Last dimension of input must be equal to the size of the "
              "underlying rtepack type.  "
              "The shape is {:B,} and the size of the underlying rtepack type is {}",
              x.shape(),
              sz));
        }

        std::array<Index, M> shape{};
        for (Size i = 0; i < M; i++) shape[i] = x.extent(i);

        new (y) matpack::data_t<T, M>(shape);

        auto outview = y->view_as(y->size());
        auto inview  = x.view_as(y->size(), sz);
        for (Size i = 0; i < y->size(); i++) {
          for (Size j = 0; j < sz; j++) {
            outview[i][j] = inview[i, j];
          }
        }
      },
      "x"_a);
  py::implicitly_convertible<matpack::data_t<U, M + 1>,
                             matpack::data_t<T, M>>();

  using nd =
      py::ndarray<py::numpy, Numeric, py::ndim<M + sizeof...(N)>, py::c_contig>;
  c.def(
      "__array__",
      [](matpack::data_t<T, M> &v,
         py::object dtype,
         py::object copy) -> std::variant<nd, py::object> {
        constexpr auto n = M + sizeof...(N);
        std::array<size_t, n> shape{};
        stdr::copy(v.shape(), shape.begin());
        stdr::copy(std::array{N...}, shape.begin() + M);
        auto np = py::module_::import_("numpy");
        auto x  = nd(v.data_handle(), n, shape.data(), py::cast(&v));

        if (not dtype.is_none()) {
          return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
        }

        return x.cast((copy.is_none() or not py::bool_(copy))
                          ? py::rv_policy::automatic_reference
                          : py::rv_policy::copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none(),
      "Returns a :class:`~numpy.ndarray` of the object.  The last dimensions are fixed to the underlying rtepack type's size.");

  c.def_prop_rw(
      "value",
      [](py::object &x) { return x.attr("__array__")(); },
      [](matpack::data_t<T, M> &x, matpack::data_t<T, M> &y) { x = y; },
      "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");

  common_ndarray(c);
}

void py_rtepack(py::module_ &m) try {
  py::class_<Stokvec> sv(m, "Stokvec");
  sv.def(py::init_implicit<Numeric>())
      .def("__init__",
           [](Stokvec *s, const PolarizationChoice p) {
             new (s) Stokvec{rtepack::to_stokvec(p)};
           })
      .def("__init__",
           [](Stokvec *s, const String &p) {
             new (s) Stokvec{rtepack::to_stokvec(to<PolarizationChoice>(p))};
           })
      .def(py::init_implicit<std::array<Numeric, 4>>())
      .def_static(
          "linpol",
          [](const Numeric angle) {
            return Stokvec{1.0,
                           Conversion::cosd(2.0 * angle),
                           Conversion::sind(2.0 * angle),
                           0.0};
          },
          "Returns [1.0, cos(2*angle), sin(2*angle), 0.0], the linear polarization vector for a given angle",
          "angle"_a)
      .def_static(
          "cirpol",
          [](const Numeric angle) {
            return Stokvec{1.0, 0.0, 0.0, Conversion::sind(angle)};
          },
          "Returns [1.0, 0.0, 0.0, sin(angle)], the circular polarization vector for a given phase delay angle",
          "angle"_a)
      .def(
          "__array__",
          [](Stokvec &v, py::object dtype, py::object copy)
              -> std::variant<
                  py::ndarray<py::numpy, Numeric, py::shape<4>, py::c_contig>,
                  py::object> {
            std::array<size_t, 1> shape = {4};
            auto np                     = py::module_::import_("numpy");
            auto x =
                py::ndarray<py::numpy, Numeric, py::shape<4>, py::c_contig>(
                    v.data.data(), 1, shape.data(), py::cast(&v));

            if (not dtype.is_none()) {
              return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
            }

            return x.cast((copy.is_none() or not py::bool_(copy))
                              ? py::rv_policy::automatic_reference
                              : py::rv_policy::copy);
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
  vsv.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Stokvec>>());
  rtepack_array<Stokvec, 1, 4>(vsv);
  generic_interface(vsv);

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

  py::class_<Propmat> pm(m, "Propmat");
  pm.def(py::init_implicit<Numeric>())
      .def(py::init_implicit<std::array<Numeric, 7>>())
      .def(
          "__array__",
          [](Propmat &x, py::object dtype, py::object copy)
              -> std::variant<
                  py::ndarray<py::numpy, Numeric, py::shape<7>, py::c_contig>,
                  py::object> {
            std::array<size_t, 1> shape = {7};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<7>, py::c_contig>(
                    x.data.data(), 1, shape.data(), py::cast(&x));

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
      .def(
          "as_matrix",
          [](Propmat &x) { return to_matrix(x); },
          "Returns the Propmat as a matrix.")
      .def_prop_rw(
          "value",
          [](py::object &x) { return x.attr("__array__")(); },
          [](Propmat &x, Propmat &y) { x = y; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");
  common_ndarray(pm);
  generic_interface(pm);

  auto apm = py::bind_vector<std::vector<Propmat>>(m, "ArrayOfPropmat");
  generic_interface(apm);
  vector_interface(apm);
  apm.doc() = "A list of :class:`~pyarts3.arts.Propmat`";

  py::class_<PropmatVector> vpm(m, "PropmatVector");
  vpm.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Propmat>>());
  rtepack_array<Propmat, 1, 7>(vpm);
  generic_interface(vpm);

  py::class_<PropmatMatrix> mpm(m, "PropmatMatrix");
  rtepack_array<Propmat, 2, 7>(mpm);
  generic_interface(mpm);

  py::class_<Muelmat> mm(m, "Muelmat");
  mm.def(py::init_implicit<Numeric>())
      .def(py::init_implicit<std::array<Numeric, 16>>())
      .def("is_polarized",
           &Muelmat::is_polarized,
           "Check if the Mueller matrix represents a polarized state.")
      .def(
          "__array__",
          [](Muelmat &x,
             py::object dtype,
             py::object copy) -> std::variant<py::ndarray<py::numpy,
                                                          Numeric,
                                                          py::shape<4, 4>,
                                                          py::c_contig>,
                                              py::object> {
            std::array<size_t, 2> shape = {4, 4};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<4, 4>, py::c_contig>(
                    x.data.data(), 2, shape.data(), py::cast(&x));

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
      [](const MuelmatVector &m) {
        return stdr::any_of(
            m, [](const Muelmat &mm) { return mm.is_polarized(); });
      },
      "Check if the Mueller matrix represents a polarized state.");
  vmm.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Muelmat>>());
  rtepack_array<Muelmat, 1, 4, 4>(vmm);
  generic_interface(vmm);

  py::class_<MuelmatMatrix> mmm(m, "MuelmatMatrix");
  rtepack_array<Muelmat, 2, 4, 4>(mmm);
  generic_interface(mmm);

  py::class_<MuelmatTensor3> mt3(m, "MuelmatTensor3");
  rtepack_array<Muelmat, 3, 4, 4>(mt3);
  generic_interface(mt3);

  py::class_<Specmat> cmm(m, "Specmat");
  cmm.def(py::init_implicit<Complex>())
      .def(py::init_implicit<std::array<Complex, 16>>())
      .def(
          "__array__",
          [](Specmat &x,
             py::object dtype,
             py::object copy) -> std::variant<py::ndarray<py::numpy,
                                                          Complex,
                                                          py::shape<4, 4>,
                                                          py::c_contig>,
                                              py::object> {
            std::array<size_t, 2> shape = {4, 4};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Complex, py::shape<4, 4>, py::c_contig>(
                    x.data.data(), 2, shape.data(), py::cast(&x));

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
          [](py::object &x) { return x.attr("__array__")(); },
          [](Specmat &x, Specmat &y) { x = y; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");
  common_ndarray(cmm);
  generic_interface(cmm);

  auto asp =
      py::bind_vector<std::vector<Specmat>, py::rv_policy::reference_internal>(
          m, "ArrayOfSpecmat");
  asp.doc() = "A list of :class:`~pyarts3.arts.Specmat`";
  vector_interface(asp);
  generic_interface(asp);

  py::class_<SpecmatVector> vcmm(m, "SpecmatVector");
  vcmm.def(py::init_implicit<std::vector<Complex>>())
      .def(py::init_implicit<std::vector<Specmat>>());
  rtepack_array<Specmat, 1, 4, 4>(vcmm);
  vcmm.doc() = "A vector of :class:`~pyarts3.arts.Specmat`";
  generic_interface(vcmm);

  py::class_<SpecmatMatrix> mcmm(m, "SpecmatMatrix");
  rtepack_array<Specmat, 2, 4, 4>(mcmm);
  generic_interface(mcmm);

  py::class_<SpecmatTensor3> cmt3(m, "SpecmatTensor3");
  rtepack_array<Specmat, 3, 4, 4>(cmt3);
  cmt3.doc() = "A 3-tensor of :class:`~pyarts3.arts.Specmat`";
  generic_interface(cmt3);

  auto a1 =
      py::bind_vector<ArrayOfPropmatVector, py::rv_policy::reference_internal>(
          m, "ArrayOfPropmatVector");
  generic_interface(a1);
  vector_interface(a1);
  auto a2 = py::bind_vector<ArrayOfArrayOfPropmatVector,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfPropmatVector");
  generic_interface(a2);
  vector_interface(a2);
  auto a3 =
      py::bind_vector<ArrayOfPropmatMatrix, py::rv_policy::reference_internal>(
          m, "ArrayOfPropmatMatrix");
  generic_interface(a3);
  vector_interface(a3);
  auto a4 = py::bind_vector<ArrayOfArrayOfPropmatMatrix,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfPropmatMatrix");
  generic_interface(a4);
  vector_interface(a4);

  auto b1 =
      py::bind_vector<ArrayOfMuelmatVector, py::rv_policy::reference_internal>(
          m, "ArrayOfMuelmatVector");
  b1.def(
      "is_polarized",
      [](const ArrayOfMuelmatVector &m) {
        return stdr::any_of(m, [](const MuelmatVector &mm) {
          return stdr::any_of(
              mm, [](const Muelmat &m) { return m.is_polarized(); });
        });
      },
      "Check if the Mueller matrix represents a polarized state.");
  generic_interface(b1);
  vector_interface(b1);
  auto b2 = py::bind_vector<ArrayOfArrayOfMuelmatVector,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfMuelmatVector");
  generic_interface(b2);
  vector_interface(b2);
  auto b3 =
      py::bind_vector<ArrayOfMuelmatMatrix, py::rv_policy::reference_internal>(
          m, "ArrayOfMuelmatMatrix");
  generic_interface(b3);
  vector_interface(b3);
  auto b4 = py::bind_vector<ArrayOfArrayOfMuelmatMatrix,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfMuelmatMatrix");
  generic_interface(b4);
  vector_interface(b4);
  auto b5 =
      py::bind_vector<ArrayOfMuelmatTensor3, py::rv_policy::reference_internal>(
          m, "ArrayOfMuelmatTensor3");
  generic_interface(b5);
  vector_interface(b5);

  auto c1 =
      py::bind_vector<ArrayOfStokvecVector, py::rv_policy::reference_internal>(
          m, "ArrayOfStokvecVector");
  generic_interface(c1);
  vector_interface(c1);
  auto c2 = py::bind_vector<ArrayOfArrayOfStokvecVector,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfStokvecVector");
  generic_interface(c2);
  vector_interface(c2);
  auto c3 =
      py::bind_vector<ArrayOfStokvecMatrix, py::rv_policy::reference_internal>(
          m, "ArrayOfStokvecMatrix");
  generic_interface(c3);
  vector_interface(c3);
  auto c4 = py::bind_vector<ArrayOfArrayOfStokvecMatrix,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfStokvecMatrix");
  generic_interface(c4);
  vector_interface(c4);
  auto c5 =
      py::bind_vector<ArrayOfStokvecTensor3, py::rv_policy::reference_internal>(
          m, "ArrayOfStokvecTensor3");
  generic_interface(c5);
  vector_interface(c5);

  auto d3 =
      py::bind_vector<ArrayOfSpecmatMatrix, py::rv_policy::reference_internal>(
          m, "ArrayOfSpecmatMatrix");
  generic_interface(d3);
  vector_interface(d3);

  auto rtepack  = m.def_submodule("rtepack");
  rtepack.doc() = "Interface to some of the core RTE functionality";
  rtepack.def(
      "two_level_exp",
      [](const ArrayOfPropmatVector &K,
         const ArrayOfPropmatMatrix &dK,
         const Vector &r,
         const Tensor3 &dr) {
        ArrayOfMuelmatVector T;
        ArrayOfMuelmatTensor3 dT;

        rtepack::two_level_exp(T, dT, K, dK, r, dr);

        return std::pair{T, dT};
      },
      "K"_a,
      "dK"_a,
      "r"_a,
      "dr"_a,
      "Returns the two-level exponential of the input matrices");

  rtepack.def(
      "two_level_radiative_transfer",
      [](const ArrayOfMuelmatVector &Ts,
         const ArrayOfMuelmatTensor3 &dTs,
         const ArrayOfStokvecVector &Js,
         const ArrayOfStokvecMatrix &dJs,
         const StokvecVector &I0) {
        StokvecVector I;
        ArrayOfStokvecMatrix dI;

        const auto Pi = forward_cumulative_transmission(Ts);
        rtepack::two_level_linear_emission_step_by_step_full(
            I, dI, Ts, Pi, dTs, Js, dJs, I0);

        return std::pair{I, dI};
      },
      "Ts"_a,
      "dTs"_a,
      "Js"_a,
      "dJs"_a,
      "I0"_a,
      "Returns the two-level radiative transfer of the input matrices");

  auto rp = m.def_submodule("rtepack");
  auto tr = py::class_<rtepack::tran>(rp, "tran");
  generic_interface(tr);
  tr.doc() =
      "Class for computing the transmission Mueller matrix and its derivative";
  tr.def(py::init<Propmat, Propmat, Numeric>(), "k1"_a, "k2"_a, "r"_a)
      .def("__call__", &rtepack::tran::operator(), "Returns the Mueller matrix")
      .def("deriv",
           &rtepack::tran::deriv,
           "t"_a,
           "k1"_a,
           "k2"_a,
           "dk"_a,
           "r"_a,
           "dr"_a,
           "Returns the derivative of the Mueller matrix")
      .def("expm1",
           &rtepack::tran::expm1,
           "Returns the Mueller matrix minus the identity matrix");
  rp.def(
        "sqrt",
        &rtepack::sqrt,
        "K"_a,
        R"(Returns the square root of the Propagation matrix as a spectral matrix

.. math::
    \mathrm{K} = \mathrm{S} \mathrm{S},

The return value of this method is the spectral matrix :math:`\mathrm{S}`.

Parameters
----------
K : Specmat
    The propagation matrix

Returns
-------
Specmat
    The square root of the Propagation matrix as a spectral matrix.
)")
      .def(
          "logK",
          &rtepack::logK,
          "m"_a,
          R"(Returns the logarithm of the Mueller matrix as a Propagation matrix

This only works for a Mueller matrix that has been generated from a Propagation matrix using

.. math::
    \mathrm{M} = e^{-\mathrm{K} r},

where :math:`\mathrm{M}` is the Mueller matrix, :math:`\mathrm{K}` is the Propagation matrix, and :math:`r` is the distance.

The return value of this method is then not the original Propagation matrix, but rather the logarithm of the Mueller matrix
as if it were generated from a Propagation matrix using the same equation.  That is

.. math::
    -\mathrm{K}r = \log(\mathrm{M})

Parameters
----------
m : Muelmat
    The Mueller matrix to compute the logarithm of.

Returns
-------
Propmat
    The logarithm of the Mueller matrix as a Propagation matrix.
)");
} catch (std::exception &e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize rtepack\n{}", e.what()));
}
}  // namespace Python
