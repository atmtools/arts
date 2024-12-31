#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/variant.h>

#include <memory>
#include <utility>
#include <variant>

#include "covariance_matrix.h"
#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "nanobind/nanobind.h"
#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_sparse(py::module_& m) try {
  py::class_<Sparse> sp(m, "Sparse");
  workspace_group_interface(sp);
  sp.def(
        "__init__",
        [](Sparse* s, Eigen::SparseMatrix<Numeric, Eigen::RowMajor> es) {
          new (s) Sparse{};
          s->matrix.swap(es);
        },
        "From :class:`scipy.sparse.csr_matrix`")
      .def(
          "__getitem__",
          [](Sparse& x, std::tuple<Index, Index> ind) -> Numeric& {
            const auto [r, c] = ind;
            if (r < 0 or r >= x.nrows())
              throw std::out_of_range(std::format("row {}", r));
            if (c < 0 or c >= x.ncols())
              throw std::out_of_range(std::format("col {}", c));
            return x.rw(r, c);
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](Sparse& x, std::tuple<Index, Index> ind, Numeric y) {
             const auto [r, c] = ind;
             if (r < 0 or r >= x.nrows())
               throw std::out_of_range(std::format("row {}", r));
             if (c < 0 or c >= x.ncols())
               throw std::out_of_range(std::format("col {}", c));
             x.rw(r, c) = y;
           })
      .def(
          "tocsr",
          [](Sparse& sp) { return sp.matrix; },
          R"(Make a :class:`scipy.sparse.csr_matrix` from the sparse matrix

Returns
-------
arr : :class:`scipy.sparse.csr_matrix`
    A sparse array
)")
      .def("__getstate__",
           [](const Sparse& self) {
             return std::tuple<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>>{
                 self.matrix};
           })
      .def("__setstate__",
           [](Sparse* self,
              const std::tuple<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>>&
                  state) {
             new (self) Sparse{};
             self->matrix = std::get<0>(state);
           });
  py::implicitly_convertible<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>,
                             Sparse>();

  sp.def(
      "arcsin",
      [](py::object& self) {
        return self.attr("tocsr")().attr("arcsin")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "arcsinh",
      [](py::object& self) {
        return self.attr("tocsr")().attr("arcsinh")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "arctan",
      [](py::object& self) {
        return self.attr("tocsr")().attr("arctan")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "arctanh",
      [](py::object& self) {
        return self.attr("tocsr")().attr("arctanh")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "argmax",
      [](py::object& self) {
        return self.attr("tocsr")().attr("argmax")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "argmin",
      [](py::object& self) {
        return self.attr("tocsr")().attr("argmin")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "asformat",
      [](py::object& self) {
        return self.attr("tocsr")().attr("asformat")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "asfptype",
      [](py::object& self) {
        return self.attr("tocsr")().attr("asfptype")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "astype",
      [](py::object& self) {
        return self.attr("tocsr")().attr("astype")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "ceil",
      [](py::object& self) {
        return self.attr("tocsr")().attr("ceil")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "check_format",
      [](py::object& self) {
        return self.attr("tocsr")().attr("check_format")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "conj",
      [](py::object& self) {
        return self.attr("tocsr")().attr("conj")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "conjugate",
      [](py::object& self) {
        return self.attr("tocsr")().attr("conjugate")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "copy",
      [](py::object& self) {
        return self.attr("tocsr")().attr("copy")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "count_nonzero",
      [](py::object& self) {
        return self.attr("tocsr")().attr("count_nonzero")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "data",
      [](py::object& self) {
        return self.attr("tocsr")().attr("data")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "deg2rad",
      [](py::object& self) {
        return self.attr("tocsr")().attr("deg2rad")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "diagonal",
      [](py::object& self) {
        return self.attr("tocsr")().attr("diagonal")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "dot",
      [](py::object& self) {
        return self.attr("tocsr")().attr("dot")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "dtype",
      [](py::object& self) {
        return self.attr("tocsr")().attr("dtype")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "eliminate_zeros",
      [](py::object& self) {
        return self.attr("tocsr")().attr("eliminate_zeros")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "expm1",
      [](py::object& self) {
        return self.attr("tocsr")().attr("expm1")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "floor",
      [](py::object& self) {
        return self.attr("tocsr")().attr("floor")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "format",
      [](py::object& self) {
        return self.attr("tocsr")().attr("format")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "getH",
      [](py::object& self) {
        return self.attr("tocsr")().attr("getH")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "get_shape",
      [](py::object& self) {
        return self.attr("tocsr")().attr("get_shape")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "getcol",
      [](py::object& self) {
        return self.attr("tocsr")().attr("getcol")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "getformat",
      [](py::object& self) {
        return self.attr("tocsr")().attr("getformat")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "getmaxprint",
      [](py::object& self) {
        return self.attr("tocsr")().attr("getmaxprint")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "getnnz",
      [](py::object& self) {
        return self.attr("tocsr")().attr("getnnz")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "getrow",
      [](py::object& self) {
        return self.attr("tocsr")().attr("getrow")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "has_canonical_format",
      [](py::object& self) {
        return self.attr("tocsr")().attr("has_canonical_format")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "has_sorted_indices",
      [](py::object& self) {
        return self.attr("tocsr")().attr("has_sorted_indices")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "imag",
      [](py::object& self) {
        return self.attr("tocsr")().attr("imag")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "indices",
      [](py::object& self) {
        return self.attr("tocsr")().attr("indices")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "indptr",
      [](py::object& self) {
        return self.attr("tocsr")().attr("indptr")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "log1p",
      [](py::object& self) {
        return self.attr("tocsr")().attr("log1p")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "max",
      [](py::object& self) {
        return self.attr("tocsr")().attr("max")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "maximum",
      [](py::object& self) {
        return self.attr("tocsr")().attr("maximum")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "maxprint",
      [](py::object& self) {
        return self.attr("tocsr")().attr("maxprint")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "mean",
      [](py::object& self) {
        return self.attr("tocsr")().attr("mean")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "min",
      [](py::object& self) {
        return self.attr("tocsr")().attr("min")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "minimum",
      [](py::object& self) {
        return self.attr("tocsr")().attr("minimum")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "multiply",
      [](py::object& self) {
        return self.attr("tocsr")().attr("multiply")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "nanmax",
      [](py::object& self) {
        return self.attr("tocsr")().attr("nanmax")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "nanmin",
      [](py::object& self) {
        return self.attr("tocsr")().attr("nanmin")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "ndim",
      [](py::object& self) {
        return self.attr("tocsr")().attr("ndim")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "nnz",
      [](py::object& self) {
        return self.attr("tocsr")().attr("nnz")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "nonzero",
      [](py::object& self) {
        return self.attr("tocsr")().attr("nonzero")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "power",
      [](py::object& self) {
        return self.attr("tocsr")().attr("power")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "prune",
      [](py::object& self) {
        return self.attr("tocsr")().attr("prune")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "rad2deg",
      [](py::object& self) {
        return self.attr("tocsr")().attr("rad2deg")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "real",
      [](py::object& self) {
        return self.attr("tocsr")().attr("real")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "reshape",
      [](py::object& self) {
        return self.attr("tocsr")().attr("reshape")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "resize",
      [](py::object& self) {
        return self.attr("tocsr")().attr("resize")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "rint",
      [](py::object& self) {
        return self.attr("tocsr")().attr("rint")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "set_shape",
      [](py::object& self) {
        return self.attr("tocsr")().attr("set_shape")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "setdiag",
      [](py::object& self) {
        return self.attr("tocsr")().attr("setdiag")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "shape",
      [](py::object& self) {
        return self.attr("tocsr")().attr("shape")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "sign",
      [](py::object& self) {
        return self.attr("tocsr")().attr("sign")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "sin",
      [](py::object& self) {
        return self.attr("tocsr")().attr("sin")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "sinh",
      [](py::object& self) {
        return self.attr("tocsr")().attr("sinh")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "size",
      [](py::object& self) {
        return self.attr("tocsr")().attr("size")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "sort_indices",
      [](py::object& self) {
        return self.attr("tocsr")().attr("sort_indices")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "sorted_indices",
      [](py::object& self) {
        return self.attr("tocsr")().attr("sorted_indices")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "sqrt",
      [](py::object& self) {
        return self.attr("tocsr")().attr("sqrt")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "sum",
      [](py::object& self) {
        return self.attr("tocsr")().attr("sum")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "sum_duplicates",
      [](py::object& self) {
        return self.attr("tocsr")().attr("sum_duplicates")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "tan",
      [](py::object& self) {
        return self.attr("tocsr")().attr("tan")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "tanh",
      [](py::object& self) {
        return self.attr("tocsr")().attr("tanh")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "tobsr",
      [](py::object& self) {
        return self.attr("tocsr")().attr("tobsr")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "tocoo",
      [](py::object& self) {
        return self.attr("tocsr")().attr("tocoo")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "tocsc",
      [](py::object& self) {
        return self.attr("tocsr")().attr("tocsc")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "todense",
      [](py::object& self) {
        return self.attr("tocsr")().attr("todense")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "todia",
      [](py::object& self) {
        return self.attr("tocsr")().attr("todia")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "todok",
      [](py::object& self) {
        return self.attr("tocsr")().attr("todok")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "tolil",
      [](py::object& self) {
        return self.attr("tocsr")().attr("tolil")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "trace",
      [](py::object& self) {
        return self.attr("tocsr")().attr("trace")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "transpose",
      [](py::object& self) {
        return self.attr("tocsr")().attr("transpose")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");
  sp.def(
      "trunc",
      [](py::object& self) {
        return self.attr("tocsr")().attr("trunc")();
      },
      "See :class:`scipy.sparse.csr_matrix`, but for no input arguments.");

  auto a1 = py::bind_vector<ArrayOfSparse, py::rv_policy::reference_internal>(
      m, "ArrayOfSparse");
  workspace_group_interface(a1);
  vector_interface(a1);

  py::class_<Block>(m, "Block")
      .def(py::init<Range, Range, IndexPair, std::shared_ptr<Matrix>>(),
           "By value, dense")
      .def(py::init<Range, Range, IndexPair, std::shared_ptr<Sparse>>(),
           "By value, sparse")
      .PythonInterfaceCopyValue(Block)
      .def_prop_rw(
          "matrix",
          [](Block& x) -> std::variant<Matrix*, Sparse*> {
            if (x.is_dense()) return &x.get_dense();
            return &x.get_sparse();
          },
          [](Block& x, std::variant<Matrix*, Sparse*> y) {
            if (std::holds_alternative<Matrix*>(y)) {
              x.set_matrix(
                  std::make_shared<Matrix>(**std::get_if<Matrix*>(&y)));
            } else {
              x.set_matrix(
                  std::make_shared<Sparse>(**std::get_if<Sparse*>(&y)));
            }
          },
          ":class:`~pyarts.arts.Matrix` or :class:`~pyarts.arts.Sparse` The matrix held inside the instance")
      .def("__getstate__",
           [](const Block& self) {
             if (self.is_sparse())
               return std::
                   tuple<Range, Range, IndexPair, std::variant<Matrix, Sparse>>{
                       self.get_row_range(),
                       self.get_column_range(),
                       self.get_indices(),
                       self.get_sparse()};
             return std::
                 tuple<Range, Range, IndexPair, std::variant<Matrix, Sparse>>{
                     self.get_row_range(),
                     self.get_column_range(),
                     self.get_indices(),
                     self.get_dense()};
           })
      .def(
          "__setstate__",
          [](Block* self,
             const std::
                 tuple<Range, Range, IndexPair, std::variant<Matrix, Sparse>>&
                     state) {
            auto row_range    = std::get<0>(state);
            auto column_range = std::get<1>(state);
            auto indices      = std::get<2>(state);
            auto mat          = std::get<3>(state);

            if (std::holds_alternative<Sparse>(mat))
              new (self) Block(row_range,
                               column_range,
                               indices,
                               std::make_shared<Sparse>(std::get<Sparse>(mat)));
            else
              new (self) Block(row_range,
                               column_range,
                               indices,
                               std::make_shared<Matrix>(std::get<Matrix>(mat)));
          })
      .doc() = "A single block matrix";

  py::class_<BlockMatrix> bm(m, "BlockMatrix");
  bm.def(py::init_implicit<Matrix>());
  bm.def(py::init_implicit<Sparse>());
  bm.def(
      "__init__",
      [](BlockMatrix* s, Eigen::SparseMatrix<Numeric, Eigen::RowMajor> es) {
        new (s) BlockMatrix{};
        *s = py::cast<Sparse>(py::type<Sparse>()(es));
      },
      "From :class:`scipy.sparse.csr_matrix`");
  py::implicitly_convertible<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>,
                             BlockMatrix>();
  bm.def(
      "__init__",
      [](BlockMatrix* v,
         const py::ndarray<py::numpy, Numeric, py::ndim<2>, py::c_contig>& a) {
        new (v) BlockMatrix{};
        *v = py::cast<Matrix>(py::type<Matrix>()(a));
      },
      "a"_a);
  py::implicitly_convertible<
      py::ndarray<py::numpy, Numeric, py::ndim<2>, py::c_contig>,
      BlockMatrix>();
  bm.def_prop_rw(
      "matrix",
      [](BlockMatrix& bm) -> std::variant<Matrix, Sparse> {
        if (bm.not_null()) {
          if (bm.is_dense()) return bm.dense();
          return bm.sparse();
        }

        return Matrix{};
      },
      [](BlockMatrix& bm, const std::variant<Matrix, Sparse>& mat) {
        std::visit([&bm](auto& m) { bm = m; }, mat);
      },
      "The matrix of the block");
  bm.def(
      "__array__",
      [](py::object& v, py::object dtype, py::object copy) {
        return v.attr("matrix").attr("__array__")(dtype, copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none(),
      "Returns a :class:`~numpy.ndarray` of the object.");
  bm.def_prop_rw(
      "value",
      [](py::object& x) { return x.attr("__array__")("copy"_a = false); },
      [](BlockMatrix& a, const std::variant<Matrix, Sparse>& b) {
        std::visit([&a](auto& c) { a = c; }, b);
      },
      "A :class:`~numpy.ndarray` or :class:`scipy.sparse.csr_matrix` of the object.");
  common_ndarray(bm);
  workspace_group_interface(bm);

  py::class_<CovarianceMatrix> covm(m, "CovarianceMatrix");
  workspace_group_interface(covm);
  covm.def_prop_rw(
          "blocks",
          [](CovarianceMatrix& x) { return x.get_blocks(); },
          [](CovarianceMatrix& x, std::vector<Block> y) {
            x.get_blocks() = std::move(y);
          },
          ":class:`list` of :class:`~pyarts.arts.Block`")
      .def("__getstate__",
           [](CovarianceMatrix& self) {
             return std::tuple<std::vector<Block>, std::vector<Block>>(
                 self.get_blocks(), self.get_inverse_blocks());
           })
      .def("__setstate__",
           [](CovarianceMatrix* self,
              const std::tuple<std::vector<Block>, std::vector<Block>>& state) {
             new (self) CovarianceMatrix{};
             self->get_blocks()         = std::get<0>(state);
             self->get_inverse_blocks() = std::get<1>(state);
           });
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize sparse\n{}", e.what()));
}
}  // namespace Python
