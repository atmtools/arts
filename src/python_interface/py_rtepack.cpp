#include <enums.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <python_interface.h>
#include <rtepack.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"

namespace Python {
template <typename T, Index M, size_t... N>
void rtepack_array(py::class_<matpack::matpack_data<T, M>> &c) {
  c.def(
      "__array__",
      [](matpack::matpack_data<T, M> &v, py::object dtype, py::object copy) {
        constexpr auto n = M + sizeof...(N);
        std::array<size_t, n> shape{};
        std::ranges::copy(v.shape(), shape.begin());
        std::ranges::copy(std::array{N...}, shape.begin() + M);
        auto np = py::module_::import_("numpy");
        auto x  = py::ndarray<py::numpy, Numeric, py::ndim<n>, py::c_contig>(
            v.data_handle(), n, shape.data(), py::handle());
        return np.attr("asarray")(
            x, py::arg("dtype") = dtype, py::arg("copy") = copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none());

  c.def_prop_rw(
      "value",
      [](matpack::matpack_data<T, M> &x) {
        py::object np = py::module_::import_("numpy");
        return np.attr("asarray")(x, py::arg("copy") = false);
      },
      [](matpack::matpack_data<T, M> &x, matpack::matpack_data<T, M> &y) {
        x = y;
      });

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
      .def_static(
          "linpol",
          [](const Numeric angle) {
            return Stokvec{1.0,
                           Conversion::cosd(2.0 * angle),
                           Conversion::sind(2.0 * angle),
                           0.0};
          },
          "Returns [1.0, cos(2*angle), sin(2*angle), 0.0], the linear polarization vector for a given angle",
          py::arg("angle"))
      .def_static(
          "cirpol",
          [](const Numeric angle) {
            return Stokvec{1.0, 0.0, 0.0, Conversion::sind(angle)};
          },
          "Returns [1.0, 0.0, 0.0, sin(angle)], the circular polarization vector for a given phase delay angle",
          py::arg("angle"))
      .def(
          "__array__",
          [](Stokvec &v, py::object dtype, py::object copy) {
            std::array<size_t, 1> shape = {4};
            auto np                     = py::module_::import_("numpy");
            auto x =
                py::ndarray<py::numpy, Numeric, py::shape<4>, py::c_contig>(
                    v.data.data(), 1, shape.data(), py::handle());
            return np.attr("asarray")(
                x, py::arg("dtype") = dtype, py::arg("copy") = copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](Stokvec &x) {
            py::object np = py::module_::import_("numpy");
            return np.attr("asarray")(x, py::arg("copy") = false);
          },
          [](Stokvec &x, Stokvec &y) { x = y; });
  common_ndarray(sv);
  workspace_group_interface(sv);
  py::implicitly_convertible<PolarizationChoice, Stokvec>();
  py::implicitly_convertible<String, Stokvec>();

  py::class_<StokvecVector> vsv(m, "StokvecVector");
  vsv.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Stokvec>>());
  rtepack_array<Stokvec, 1, 4>(vsv);
  workspace_group_interface(vsv);

  py::class_<StokvecMatrix> msv(m, "StokvecMatrix");
  rtepack_array<Stokvec, 2, 4>(msv);
  workspace_group_interface(msv);

  py::class_<StokvecTensor3> t3sv(m, "StokvecTensor3");
  rtepack_array<Stokvec, 3, 4>(t3sv);
  workspace_group_interface(t3sv);

  py::class_<StokvecTensor4> t4sv(m, "StokvecTensor4");
  rtepack_array<Stokvec, 4, 4>(t4sv);
  workspace_group_interface(t4sv);

  py::class_<StokvecTensor5> t5sv(m, "StokvecTensor5");
  rtepack_array<Stokvec, 5, 4>(t5sv);
  workspace_group_interface(t5sv);

  py::class_<StokvecTensor6> t6sv(m, "StokvecTensor6");
  rtepack_array<Stokvec, 6, 4>(t6sv);
  workspace_group_interface(t6sv);

  py::class_<Propmat> pm(m, "Propmat");
  pm.def(py::init_implicit<Numeric>())
      .def(py::init<Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric>())
      .def(
          "__array__",
          [](Propmat &x, py::object dtype, py::object copy) {
            std::array<size_t, 1> shape = {7};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<7>, py::c_contig>(
                    x.data.data(), 1, shape.data(), py::handle());
            return np.attr("asarray")(
                w, py::arg("dtype") = dtype, py::arg("copy") = copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](Propmat &x) {
            py::object np = py::module_::import_("numpy");
            return np.attr("asarray")(x, py::arg("copy") = false);
          },
          [](Propmat &x, Propmat &y) { x = y; });
  common_ndarray(pm);
  workspace_group_interface(pm);

  py::class_<PropmatVector> vpm(m, "PropmatVector");
  vpm.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Propmat>>());
  rtepack_array<Propmat, 1, 7>(vpm);
  workspace_group_interface(vpm);

  py::class_<PropmatMatrix> mpm(m, "PropmatMatrix");
  rtepack_array<Propmat, 2, 7>(mpm);
  workspace_group_interface(mpm);

  py::class_<Muelmat> mm(m, "Muelmat");
  mm.def(py::init_implicit<Numeric>())
      .def(py::init<Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric>())
      .def(
          "__array__",
          [](Muelmat &x, py::object dtype, py::object copy) {
            std::array<size_t, 2> shape = {4, 4};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<4, 4>, py::c_contig>(
                    x.data.data(), 2, shape.data(), py::handle());
            return np.attr("asarray")(
                w, py::arg("dtype") = dtype, py::arg("copy") = copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](Muelmat &x) {
            py::object np = py::module_::import_("numpy");
            return np.attr("asarray")(x, py::arg("copy") = false);
          },
          [](Muelmat &x, Muelmat &y) { x = y; });
  common_ndarray(mm);
  workspace_group_interface(mm);

  py::class_<MuelmatVector> vmm(m, "MuelmatVector");
  vmm.def(py::init_implicit<std::vector<Numeric>>())
      .def(py::init_implicit<std::vector<Muelmat>>());
  rtepack_array<Muelmat, 1, 4, 4>(vmm);
  workspace_group_interface(vmm);

  py::class_<MuelmatMatrix> mmm(m, "MuelmatMatrix");
  rtepack_array<Muelmat, 2, 4, 4>(mmm);
  workspace_group_interface(mmm);

  auto a1 = py::bind_vector<ArrayOfPropmatVector, py::rv_policy::reference_internal>(m, "ArrayOfPropmatVector");
  workspace_group_interface(a1);
  vector_interface(a1);
  auto a2 = py::bind_vector<ArrayOfArrayOfPropmatVector, py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfPropmatVector");
  workspace_group_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfPropmatMatrix, py::rv_policy::reference_internal>(m, "ArrayOfPropmatMatrix");
  workspace_group_interface(a3);
  vector_interface(a3);
  auto a4 = py::bind_vector<ArrayOfArrayOfPropmatMatrix, py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfPropmatMatrix");
  workspace_group_interface(a4);
  vector_interface(a4);

  auto b1 = py::bind_vector<ArrayOfMuelmatVector, py::rv_policy::reference_internal>(m, "ArrayOfMuelmatVector");
  workspace_group_interface(b1);
  vector_interface(b1);
  auto b2 = py::bind_vector<ArrayOfArrayOfMuelmatVector, py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfMuelmatVector");
  workspace_group_interface(b2);
  vector_interface(b2);
  auto b3 = py::bind_vector<ArrayOfMuelmatMatrix, py::rv_policy::reference_internal>(m, "ArrayOfMuelmatMatrix");
  workspace_group_interface(b3);
  vector_interface(b3);
  auto b4 = py::bind_vector<ArrayOfArrayOfMuelmatMatrix, py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfMuelmatMatrix");
  workspace_group_interface(b4);
  vector_interface(b4);

  auto c1 = py::bind_vector<ArrayOfStokvecVector, py::rv_policy::reference_internal>(m, "ArrayOfStokvecVector");
  workspace_group_interface(c1);
  vector_interface(c1);
  auto c2 = py::bind_vector<ArrayOfArrayOfStokvecVector, py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfStokvecVector");
  workspace_group_interface(c2);
  vector_interface(c2);
  auto c3 = py::bind_vector<ArrayOfStokvecMatrix, py::rv_policy::reference_internal>(m, "ArrayOfStokvecMatrix");
  workspace_group_interface(c3);
  vector_interface(c3);
  auto c4 = py::bind_vector<ArrayOfArrayOfStokvecMatrix, py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfStokvecMatrix");
  workspace_group_interface(c4);
  vector_interface(c4);

} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize rtepack\n", e.what()));
}
}  // namespace Python
