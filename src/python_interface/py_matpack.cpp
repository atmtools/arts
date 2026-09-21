#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/vector.h>

#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_vector.h"
#include "python_interface.h"

namespace Python {
// Defined in py_matpack_fixed.cpp and py_matpack_complex.cpp respectively.
// py_matpack.cpp used to bind all matpack types in one function; splitting
// the independent groups into their own TUs cuts the compile time and peak
// RAM of the single biggest offender without changing what gets bound (see
// the same rationale in py_griddedfield.cpp).
void py_matpack_fixed(py::module_& m);
void py_matpack_complex(py::module_& m);

void py_matpack(py::module_& m) try {
  py_matpack_fixed(m);

  py::class_<IndexVector> iv1(m, "IndexVector");
  iv1.doc() = "A vector of indices";
  matpack_interface(iv1);
  generic_interface(iv1);

  py::class_<Vector>  v1(m, "Vector");
  py::class_<Matrix>  v2(m, "Matrix");
  py::class_<Tensor3> v3(m, "Tensor3");
  py::class_<Tensor4> v4(m, "Tensor4");
  py::class_<Tensor5> v5(m, "Tensor5");
  py::class_<Tensor6> v6(m, "Tensor6");
  py::class_<Tensor7> v7(m, "Tensor7");

  matpack_interface(v1);
  matpack_interface(v2);
  matpack_interface(v3);
  matpack_interface(v4);
  matpack_interface(v5);
  matpack_interface(v6);
  matpack_interface(v7);
  generic_interface(v1);
  generic_interface(v2);
  generic_interface(v3);
  generic_interface(v4);
  generic_interface(v5);
  generic_interface(v6);
  generic_interface(v7);

  auto a1 = py::bind_vector<ArrayOfVector, py::rv_policy::reference_internal>(m, "ArrayOfVector");
  generic_interface(a1);
  vector_interface(a1);
  auto a2 = py::bind_vector<ArrayOfArrayOfVector, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfVector");
  generic_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfMatrix, py::rv_policy::reference_internal>(m, "ArrayOfMatrix");
  generic_interface(a3);
  vector_interface(a3);
  auto a4 = py::bind_vector<ArrayOfArrayOfMatrix, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfMatrix");
  generic_interface(a4);
  vector_interface(a4);
  auto a5 = py::bind_vector<ArrayOfTensor3, py::rv_policy::reference_internal>(m, "ArrayOfTensor3");
  generic_interface(a5);
  vector_interface(a5);
  auto a6 = py::bind_vector<ArrayOfArrayOfTensor3, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfTensor3");
  generic_interface(a6);
  vector_interface(a6);
  auto a7 = py::bind_vector<ArrayOfTensor4, py::rv_policy::reference_internal>(m, "ArrayOfTensor4");
  generic_interface(a7);
  vector_interface(a7);
  auto a8 = py::bind_vector<ArrayOfTensor5, py::rv_policy::reference_internal>(m, "ArrayOfTensor5");
  generic_interface(a8);
  vector_interface(a8);
  auto a9 = py::bind_vector<ArrayOfTensor6, py::rv_policy::reference_internal>(m, "ArrayOfTensor6");
  generic_interface(a9);
  vector_interface(a9);
  auto a10 = py::bind_vector<ArrayOfArrayOfTensor6, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfTensor6");
  generic_interface(a10);
  vector_interface(a10);
  auto a11 = py::bind_vector<ArrayOfTensor7, py::rv_policy::reference_internal>(m, "ArrayOfTensor7");
  generic_interface(a11);
  vector_interface(a11);
  auto a12 = py::bind_vector<ArrayOfVector2, py::rv_policy::reference_internal>(m, "ArrayOfVector2");
  generic_interface(a12);
  vector_interface(a12);
  auto a13 = py::bind_vector<ArrayOfVector3, py::rv_policy::reference_internal>(m, "ArrayOfVector3");
  generic_interface(a13);
  vector_interface(a13);

  py_matpack_complex(m);

  py::class_<AscendingGrid> g1(m, "AscendingGrid");
  matpack_grid_interface(g1);
  generic_interface(g1);
  g1.def(py::init_implicit<Vector>());
  v1.def(py::init_implicit<AscendingGrid>());

  py::class_<DescendingGrid> g2(m, "DescendingGrid");
  matpack_grid_interface(g2);
  generic_interface(g2);
  g2.def(py::init_implicit<Vector>());
  v1.def(py::init_implicit<DescendingGrid>());

  auto b1 = py::bind_vector<ArrayOfAscendingGrid, py::rv_policy::reference_internal>(m, "ArrayOfAscendingGrid");
  generic_interface(b1);
  vector_interface(b1);

  py::class_<LatGrid> gr1(m, "LatGrid");
  matpack_grid_interface(gr1);
  generic_interface(gr1);
  gr1.def(py::init_implicit<Vector>());
  v1.def(py::init_implicit<LatGrid>());

  py::class_<LonGrid> gr2(m, "LonGrid");
  matpack_grid_interface(gr2);
  generic_interface(gr2);
  gr2.def(py::init_implicit<Vector>());
  v1.def(py::init_implicit<LonGrid>());

  py::class_<ZenGrid> gr3(m, "ZenGrid");
  matpack_grid_interface(gr3);
  generic_interface(gr3);
  gr3.def(py::init_implicit<Vector>());
  v1.def(py::init_implicit<ZenGrid>());

  py::class_<AziGrid> gr4(m, "AziGrid");
  matpack_grid_interface(gr4);
  generic_interface(gr4);
  gr4.def(py::init_implicit<Vector>());
  v1.def(py::init_implicit<AziGrid>());
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize matpack\n{}", e.what()));
}
}  // namespace Python
