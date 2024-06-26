#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/vector.h>
#include <workspace.h>

#include "python_interface_value_type.h"

NB_MAKE_OPAQUE(
    std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>);
NB_MAKE_OPAQUE(std::vector<lbl::line_shape::species_model>);
NB_MAKE_OPAQUE(std::vector<lbl::line>)
NB_MAKE_OPAQUE(Array<LagrangeInterpolation>);
NB_MAKE_OPAQUE(Array<AbsorptionSingleLine>);
NB_MAKE_OPAQUE(Array<Array<SpeciesTag>>);

namespace Python {
namespace py = nanobind;
using namespace py::literals;

template <typename T>
void vector_interface(py::class_<Array<T>> &c) {
  using Vec = Array<T>;

  c.def("__getstate__", [](const Vec &v) { return std::tuple<Vec>{v}; });

  c.def("__setstate__",
        [](Vec &v, const std::tuple<Vec> &x) { new (&v) Vec{std::get<0>(x)}; });
}

// THIS IS A DIRECT COPY OF THE FUNCTION FROM bind_vector.h
// SOME MODIFICATIONS ARE DONE, SPECIFICALLY TO GET THE INIT
// FUNCTION TO WORK WITH THE ARRAY OF STRING USING OUR VALUEHOLDER
// IF THIS FAILS, WE MUST PROBABLY REVERT TO THE ORIGINAL FUNCTION
// OR REWRITE THIS ENTIRERLY.
// IT ALSO WRAPS vector_interface FOR SIMPLICITY
template <typename Value>
void value_holder_vector_interface(py::class_<Array<Value>> &cl) {
  using namespace nanobind;
  using ValueRef = Value &;
  using Vector   = Array<Value>;

  cl.def(init<>(), "Default constructor")

      .def("__len__", [](const Vector &v) { return v.size(); })

      .def(
          "__bool__",
          [](const Vector &v) { return !v.empty(); },
          "Check whether the vector is nonempty")

      .def("__repr__",
           [](handle_t<Vector> h) {
             return steal<str>(detail::repr_list(h.ptr()));
           })

      .def(
          "__iter__",
          [](Vector &v) {
            return make_iterator(
                type<Vector>(), "Iterator", v.begin(), v.end());
          },
          keep_alive<0, 1>())

      .def(
          "__getitem__",
          [](Vector &v, Py_ssize_t i) -> ValueRef {
            return v[detail::wrap(i, v.size())];
          },
          rv_policy::automatic_reference)

      .def(
          "clear", [](Vector &v) { v.clear(); }, "Remove all items from list.");

  if constexpr (detail::is_copy_constructible_v<Value>) {
    cl.def(init<const Vector &>(), "Copy constructor");

    cl.def(
        "__init__",
        [](Vector *v, typed<iterable, ValueHolder<Value>> seq) {
          new (v) Vector();
          v->reserve(len_hint(seq));
          for (handle h : seq) v->push_back(cast<ValueHolder<Value>>(h));
        },
        "Construct from an iterable object");

    implicitly_convertible<iterable, Vector>();

    cl.def(
          "append",
          [](Vector &v, const Value &value) { v.push_back(value); },
          "Append `arg` to the end of the list.")

        .def(
            "insert",
            [](Vector &v, Py_ssize_t i, const Value &x) {
              if (i < 0) i += (Py_ssize_t)v.size();
              if (i < 0 || (size_t)i > v.size()) throw index_error();
              v.insert(v.begin() + i, x);
            },
            "Insert object `arg1` before index `arg0`.")

        .def(
            "pop",
            [](Vector &v, Py_ssize_t i) {
              size_t index = detail::wrap(i, v.size());
              Value result = std::move(v[index]);
              v.erase(v.begin() + index);
              return result;
            },
            arg("index") = -1,
            "Remove and return item at `index` (default last).")

        .def(
            "extend",
            [](Vector &v, const Vector &src) {
              v.insert(v.end(), src.begin(), src.end());
            },
            "Extend `self` by appending elements from `arg`.")

        .def("__setitem__",
             [](Vector &v, Py_ssize_t i, const Value &value) {
               v[detail::wrap(i, v.size())] = value;
             })

        .def("__delitem__",
             [](Vector &v, Py_ssize_t i) {
               v.erase(v.begin() + detail::wrap(i, v.size()));
             })

        .def("__getitem__",
             [](const Vector &v, const slice &slice) -> Vector * {
               auto [start, stop, step, length] = slice.compute(v.size());
               auto *seq                        = new Vector();
               seq->reserve(length);

               for (size_t i = 0; i < length; ++i) {
                 seq->push_back(v[start]);
                 start += step;
               }

               return seq;
             })

        .def("__setitem__",
             [](Vector &v, const slice &slice, const Vector &value) {
               auto [start, stop, step, length] = slice.compute(v.size());

               if (length != value.size())
                 throw index_error(
                     "The left and right hand side of the slice "
                     "assignment have mismatched sizes!");

               for (size_t i = 0; i < length; ++i) {
                 v[start]  = value[i];
                 start    += step;
               }
             })

        .def("__delitem__", [](Vector &v, const slice &slice) {
          auto [start, stop, step, length] = slice.compute(v.size());
          if (length == 0) return;

          stop = start + (length - 1) * step;
          if (start > stop) {
            std::swap(start, stop);
            step = -step;
          }

          if (step == 1) {
            v.erase(v.begin() + start, v.begin() + stop + 1);
          } else {
            for (size_t i = 0; i < length; ++i) {
              v.erase(v.begin() + stop);
              stop -= step;
            }
          }
        });
  }

  if constexpr (detail::is_equality_comparable_v<Value>) {
    cl.def(self == self)
        .def(self != self)

        .def("__contains__",
             [](const Vector &v, const Value &x) {
               return std::find(v.begin(), v.end(), x) != v.end();
             })

        .def("__contains__",  // fallback for incompatible types
             [](const Vector &, handle) { return false; })

        .def(
            "count",
            [](const Vector &v, const Value &x) {
              return std::count(v.begin(), v.end(), x);
            },
            "Return number of occurrences of `arg`.")

        .def(
            "remove",
            [](Vector &v, const Value &x) {
              auto p = std::find(v.begin(), v.end(), x);
              if (p != v.end())
                v.erase(p);
              else
                throw value_error();
            },
            "Remove first occurrence of `arg`.");
  }

  vector_interface(cl);
}
}  // namespace Python
