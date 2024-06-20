#pragma once

#include <nanobind/nanobind.h>

namespace Python {
namespace py = nanobind;
using namespace py::literals;

void common_ndarray(auto& c) {
  c.def(
      "__abs__",
      [](py::object& self) { return self.attr("value").attr("__abs__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__abs__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__add__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__add__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__add__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__and__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__and__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__and__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__contains__",
      [](py::object& self, py::object& key_) {
        return self.attr("value").attr("__contains__")(key_);
      },
      "key"_a,
      "Wrapper for :attr:`~numpy.ndarray.__contains__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__divmod__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__divmod__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__divmod__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__eq__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__eq__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__eq__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__float__",
      [](py::object& self) { return self.attr("value").attr("__float__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__float__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__floordiv__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__floordiv__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__floordiv__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__ge__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__ge__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__ge__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__getitem__",
      [](py::object& self, py::object& key_) {
        return self.attr("value").attr("__getitem__")(key_);
      },
      "key"_a,
      "Wrapper for :attr:`~numpy.ndarray.__getitem__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__gt__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__gt__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__gt__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__iadd__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__iadd__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__iadd__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__iand__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__iand__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__iand__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__ifloordiv__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__ifloordiv__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__ifloordiv__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__ilshift__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__ilshift__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__ilshift__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__imatmul__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__imatmul__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__imatmul__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__imod__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__imod__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__imod__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__imul__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__imul__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__imul__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__index__",
      [](py::object& self) { return self.attr("value").attr("__index__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__index__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__int__",
      [](py::object& self) { return self.attr("value").attr("__int__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__int__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__invert__",
      [](py::object& self) { return self.attr("value").attr("__invert__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__invert__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__ior__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__ior__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__ior__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__ipow__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__ipow__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__ipow__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__irshift__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__irshift__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__irshift__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__isub__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__isub__")(value_);
        return self;
      },
      py::rv_policy::reference_internal,
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__isub__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__iter__",
      [](py::object& self) { return self.attr("value").attr("__iter__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__iter__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__itruediv__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__itruediv__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__itruediv__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__ixor__",
      [](py::object& self, py::object& value_) {
        self.attr("value").attr("__ixor__")(value_);
        return self;
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__ixor__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__le__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__le__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__le__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__len__",
      [](py::object& self) { return self.attr("value").attr("__len__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__len__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__lshift__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__lshift__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__lshift__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__lt__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__lt__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__lt__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__matmul__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__matmul__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__matmul__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__mod__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__mod__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__mod__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__mul__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__mul__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__mul__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__ne__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__ne__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__ne__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__neg__",
      [](py::object& self) { return self.attr("value").attr("__neg__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__neg__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__or__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__or__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__or__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__pos__",
      [](py::object& self) { return self.attr("value").attr("__pos__")(); },
      "Wrapper for :attr:`~numpy.ndarray.__pos__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__pow__",
      [](py::object& self, py::object& value_, py::object& mod_) {
        return self.attr("value").attr("__pow__")(value_, mod_);
      },
      "value"_a,
      "mod"_a = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.__pow__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__radd__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__radd__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__radd__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rand__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rand__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rand__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rdivmod__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rdivmod__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rdivmod__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rfloordiv__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rfloordiv__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rfloordiv__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rlshift__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rlshift__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rlshift__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rmatmul__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rmatmul__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rmatmul__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rmod__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rmod__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rmod__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rmul__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rmul__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rmul__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__ror__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__ror__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__ror__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rpow__",
      [](py::object& self, py::object& value_, py::object& mod_) {
        return self.attr("value").attr("__rpow__")(value_, mod_);
      },
      "value"_a,
      "mod"_a = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.__rpow__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rrshift__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rrshift__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rrshift__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rshift__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rshift__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rshift__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rsub__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rsub__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rsub__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rtruediv__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rtruediv__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rtruediv__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__rxor__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__rxor__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__rxor__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__setitem__",
      [](py::object& self, py::object& key_, py::object& value_) {
        return self.attr("value").attr("__setitem__")(key_, value_);
      },
      "key"_a,
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__setitem__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__sub__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__sub__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__sub__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__truediv__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__truediv__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__truediv__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "__xor__",
      [](py::object& self, py::object& value_) {
        return self.attr("value").attr("__xor__")(value_);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.__xor__` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "T",
      [](py::object& self) { return py::object(self.attr("value").attr("T")); },
      "Wrapper for :attr:`~numpy.ndarray.T` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  // Renamed to behave as data_ as the return type is similar
  c.def_prop_ro(
      "_base",
      [](py::object& self) {
        return py::object(self.attr("value").attr("base"));
      },
      py::keep_alive<0, 1>{},
      "Wrapper for :attr:`~numpy.ndarray.base` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  // Renamed because we often have properties with the same name
  c.def_prop_ro(
      "_data",
      [](py::object& self) {
        return py::object(self.attr("value").attr("data"));
      },
      py::keep_alive<0, 1>{},
      "Wrapper for :attr:`~numpy.ndarray.data` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "dtype",
      [](py::object& self) {
        return py::object(self.attr("value").attr("dtype"));
      },
      "Wrapper for :attr:`~numpy.ndarray.dtype` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  // Renamed to behave as data_ as the return type is similar
  c.def_prop_ro(
      "_flat",
      [](py::object& self) {
        return py::object(self.attr("value").attr("flat"));
      },
      py::keep_alive<0, 1>{},
      "Wrapper for :attr:`~numpy.ndarray.flat` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "imag",
      [](py::object& self) {
        return py::object(self.attr("value").attr("imag"));
      },
      py::keep_alive<0, 1>{},
      "Wrapper for :attr:`~numpy.ndarray.imag` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "itemsize",
      [](py::object& self) {
        return py::object(self.attr("value").attr("itemsize"));
      },
      "Wrapper for :attr:`~numpy.ndarray.itemsize` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "nbytes",
      [](py::object& self) {
        return py::object(self.attr("value").attr("nbytes"));
      },
      "Wrapper for :attr:`~numpy.ndarray.nbytes` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "ndim",
      [](py::object& self) {
        return py::object(self.attr("value").attr("ndim"));
      },
      "Wrapper for :attr:`~numpy.ndarray.ndim` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "real",
      [](py::object& self) {
        return py::object(self.attr("value").attr("real"));
      },
      py::keep_alive<0, 1>{},
      "Wrapper for :attr:`~numpy.ndarray.real` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "shape",
      [](py::object& self) {
        return py::object(self.attr("value").attr("shape"));
      },
      "Wrapper for :attr:`~numpy.ndarray.shape` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "size",
      [](py::object& self) {
        return py::object(self.attr("value").attr("size"));
      },
      "Wrapper for :attr:`~numpy.ndarray.size` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def_prop_ro(
      "strides",
      [](py::object& self) {
        return py::object(self.attr("value").attr("strides"));
      },
      "Wrapper for :attr:`~numpy.ndarray.strides` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "flatten",
      [](py::object& self, char order) {
        return self.attr("value").attr("flatten")(order);
      },
      "order"_a = 'C',
      "Wrapper for :attr:`~numpy.ndarray.flatten` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "copy",
      [](py::object& self, char order) {
        return self.attr("value").attr("copy")(order);
      },
      "order"_a = 'C',
      "Wrapper for :attr:`~numpy.ndarray.copy` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "fill",
      [](py::object& self, py::object& value) {
        return self.attr("value").attr("fill")(value);
      },
      "value"_a,
      "Wrapper for :attr:`~numpy.ndarray.fill` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "tolist",
      [](py::object& self) { return self.attr("value").attr("tolist")(); },
      "Wrapper for :attr:`~numpy.ndarray.tolist` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "conj",
      [](py::object& self) { return self.attr("value").attr("conj")(); },
      "Wrapper for :attr:`~numpy.ndarray.conj` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "conjugate",
      [](py::object& self) { return self.attr("value").attr("conjugate")(); },
      "Wrapper for :attr:`~numpy.ndarray.conjugate` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "nonzero",
      [](py::object& self) { return self.attr("value").attr("nonzero")(); },
      "Wrapper for :attr:`~numpy.ndarray.nonzero` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "cumprod",
      [](py::object& self,
         py::object& axis,
         py::object& dtype,
         py::object& out) {
        return self.attr("value").attr("cumprod")(axis, dtype, out);
      },
      "axis"_a  = py::none(),
      "dtype"_a = py::none(),
      "out"_a   = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.cumprod` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "cumsum",
      [](py::object& self,
         py::object& axis,
         py::object& dtype,
         py::object& out) {
        return self.attr("value").attr("cumsum")(axis, dtype, out);
      },
      "axis"_a  = py::none(),
      "dtype"_a = py::none(),
      "out"_a   = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.cumsum` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "diagonal",
      [](py::object& self,
         py::object& offset,
         py::object& axis1,
         py::object& axis2) {
        return self.attr("value").attr("diagonal")(offset, axis1, axis2);
      },
      "offset"_a = 0,
      "axis1"_a  = 0,
      "axis2"_a  = 1,
      "Wrapper for :attr:`~numpy.ndarray.diagonal` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "max",
      [](py::object& self,
         py::object& axis,
         py::object& out,
         py::object& keepdims) {
        return self.attr("value").attr("max")(axis, out, keepdims);
      },
      "axis"_a     = py::none(),
      "out"_a      = py::none(),
      "keepdims"_a = false,
      "Wrapper for :attr:`~numpy.ndarray.max` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "mean",
      [](py::object& self,
         py::object& axis,
         py::object& dtype,
         py::object& out,
         py::object& keepdims) {
        return self.attr("value").attr("mean")(axis, dtype, out, keepdims);
      },
      "axis"_a     = py::none(),
      "dtype"_a    = py::none(),
      "out"_a      = py::none(),
      "keepdims"_a = false,
      "Wrapper for :attr:`~numpy.ndarray.mean` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "min",
      [](py::object& self,
         py::object& axis,
         py::object& out,
         py::object& keepdims) {
        return self.attr("value").attr("min")(axis, out, keepdims);
      },
      "axis"_a     = py::none(),
      "out"_a      = py::none(),
      "keepdims"_a = false,
      "Wrapper for :attr:`~numpy.ndarray.min` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "prod",
      [](py::object& self,
         py::object& axis,
         py::object& dtype,
         py::object& out,
         py::object& keepdims) {
        return self.attr("value").attr("prod")(axis, dtype, out, keepdims);
      },
      "axis"_a     = py::none(),
      "dtype"_a    = py::none(),
      "out"_a      = py::none(),
      "keepdims"_a = false,
      "Wrapper for :attr:`~numpy.ndarray.prod` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "ravel",
      [](py::object& self, py::list& order) {
        return self.attr("value").attr("ravel")(order);
      },
      "order"_a,
      "Wrapper for :attr:`~numpy.ndarray.ravel` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "repeat",
      [](py::object& self, py::object& repeats, py::object& axis) {
        return self.attr("value").attr("repeat")(repeats, axis);
      },
      "repeats"_a,
      "axis"_a = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.repeat` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "reshape",
      [](py::object& self, py::object& shape, py::object& order) {
        return self.attr("value").attr("reshape")(shape, order);
      },
      "shape"_a,
      "order"_a = 'C',
      "Wrapper for :attr:`~numpy.ndarray.reshape` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "round",
      [](py::object& self, py::object& decimals, py::object& out) {
        return self.attr("value").attr("round")(decimals, out);
      },
      "decimals"_a = 0,
      "out"_a      = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.round` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "squeeze",
      [](py::object& self, py::object& axis) {
        return self.attr("value").attr("squeeze")(axis);
      },
      "axis"_a = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.squeeze` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "std",
      [](py::object& self,
         py::object& axis,
         py::object& dtype,
         py::object& out,
         py::object& ddof,
         py::object& keepdims) {
        return self.attr("value").attr("std")(axis, dtype, out, ddof, keepdims);
      },
      "axis"_a     = py::none(),
      "dtype"_a    = py::none(),
      "out"_a      = py::none(),
      "ddof"_a     = 0,
      "keepdims"_a = false,
      "Wrapper for :attr:`~numpy.ndarray.std` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "sum",
      [](py::object& self,
         py::object& axis,
         py::object& dtype,
         py::object& out,
         py::object& keepdims) {
        return self.attr("value").attr("sum")(axis, dtype, out, keepdims);
      },
      "axis"_a     = py::none(),
      "dtype"_a    = py::none(),
      "out"_a      = py::none(),
      "keepdims"_a = false,
      "Wrapper for :attr:`~numpy.ndarray.sum` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "trace",
      [](py::object& self,
         py::object& offset,
         py::object& axis1,
         py::object& axis2,
         py::object& dtype,
         py::object& out) {
        return self.attr("value").attr("trace")(
            offset, axis1, axis2, dtype, out);
      },
      "offset"_a = 0,
      "axis1"_a  = 0,
      "axis2"_a  = 1,
      "dtype"_a  = py::none(),
      "out"_a    = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.trace` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "transpose",
      [](py::object& self, py::list& axes) {
        return self.attr("value").attr("transpose")(axes);
      },
      "axes"_a = py::none(),
      "Wrapper for :attr:`~numpy.ndarray.transpose` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");

  c.def(
      "var",
      [](py::object& self,
         py::object& axis,
         py::object& dtype,
         py::object& out,
         py::object& ddof,
         py::object& keepdims) {
        return self.attr("value").attr("var")(axis, dtype, out, ddof, keepdims);
      },
      "axis"_a     = py::none(),
      "dtype"_a    = py::none(),
      "out"_a      = py::none(),
      "ddof"_a     = 0,
      "keepdims"_a = false,
      "Wrapper for :attr:`~numpy.ndarray.var` using ARTS types."
      "\n\n"
      "Use the original for greater control and more functionality.");
}
}  // namespace Python
