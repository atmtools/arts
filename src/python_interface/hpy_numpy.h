#pragma once

#include "hpy_numpy_impl.h"

namespace Python {
template <typename FUNC, typename... in>
auto vectorize(FUNC f, in&&... args) {
  auto np = py::module_::import_("numpy");
  return np.attr("vectorize")(py::cpp_function(f))(std::forward<in>(args)...);
}

void common_ndarray(auto& c) {
  c.def("__abs__", &npabs, "Allows `abs(self)`");

  c.def("__add__", &npadd, "value"_a, "Allows `self + value`");

  c.def("__and__", &npand, "value"_a, "Allows `self and value`");

  c.def("__contains__", &npcontains, "key"_a, "Allows `key in self`");

  c.def("__divmod__", &npdivmod, "value"_a, "Allows `divmod(self, value)`");

  c.def("__eq__", &npeq, "value"_a, "Allows `self == value`");

  c.def("__float__", &npfloat, "Allows float(self)");

  c.def("__floordiv__", &npfloordiv, "value"_a, "Allows `self // value`");

  c.def("__ge__", &npge, "value"_a, "Allows `self >= value`");

  c.def("__getitem__", &npgetitem, "key"_a, "Allows `self[key]`");

  c.def("__gt__", &npgt, "value"_a, "Allows `self > value`");

  c.def("__iadd__", &npiadd, "value"_a, "Allows `self += value`");

  c.def("__iand__", &npiand, "value"_a, "Allows `self &= value`");

  c.def("__ifloordiv__", &npifloordiv, "value"_a, "Allows `self //= value`");

  c.def("__ilshift__", &npilshift, "value"_a, "Allows `self <<= value`");

  c.def("__imatmul__", &npimatmul, "value"_a, "Allows `self @= value`");

  c.def("__imod__", &npimod, "value"_a, "Allows `self %= value`");

  c.def("__imul__", &npimul, "value"_a, "Allows `self *= value`");

  c.def("__index__", &npindex, "Allows `int(self)`");

  c.def("__int__", &npint, "Allows `int(self)`");

  c.def("__invert__", &npinvert, "Allows `~self`");

  c.def("__ior__", &npior, "value"_a, "Allows `self |= value`");

  c.def("__ipow__", &npipow, "value"_a, "Allows `self **= value`");

  c.def("__irshift__", &npirshift, "value"_a, "Allows `self >>= value`");

  c.def("__isub__",
        &npisub,
        py::rv_policy::reference_internal,
        "value"_a,
        "Allows `self -= value`");

  c.def("__iter__", &npiter, "Allows `iter(self)`");

  c.def("__itruediv__", &npitruediv, "value"_a, "Allows `self /= value`");

  c.def("__ixor__", &npixor, "value"_a, "Allows `self ^= value`");

  c.def("__le__", &nple, "value"_a, "Allows `self <= value`");

  c.def("__len__", &nplen, "Allows `len(self)`");

  c.def("__lshift__", &nplshift, "value"_a, "Allows `self << value`");

  c.def("__lt__", &nplt, "value"_a, "Allows `self < value`");

  c.def("__matmul__", &npmatmul, "value"_a, "Allows `self @ value`");

  c.def("__mod__", &npmod, "value"_a, "Allows `self % value`");

  c.def("__mul__", &npmul, "value"_a, "Allows `self * value`");

  c.def("__ne__", &npne, "value"_a, "Allows `self != value`");

  c.def("__neg__", &npneg, "Allows `-self`");

  c.def("__or__", &npor, "value"_a, "Allows `self | value`");

  c.def("__pos__", &nppos, "Allows `+self`");

  c.def("__pow__",
        &nppow,
        "value"_a,
        "mod"_a = py::none(),
        "Allows `self ** value`");

  c.def("__radd__", &npradd, "value"_a, "Allows `value + self`");

  c.def("__rand__", &nprand, "value"_a, "Allows `value & self`");

  c.def("__rdivmod__", &nprdivmod, "value"_a, "Allows `value / self`");

  c.def("__rfloordiv__", &nprfloordiv, "value"_a, "Allows `value // self`");

  c.def("__rlshift__", &nprlshift, "value"_a, "Allows `value << self`");

  c.def("__rmatmul__", &nprmatmul, "value"_a, "Allows `value @ self`");

  c.def("__rmod__", &nprmod, "value"_a, "Allows `value % self`");

  c.def("__rmul__", &nprmul, "value"_a, "Allows `value * self`");

  c.def("__ror__", &npror, "value"_a, "Allows `value | self`");

  c.def("__rpow__",
        &nprpow,
        "value"_a,
        "mod"_a = py::none(),
        "Allows `value ** self`");

  c.def("__rrshift__", &nprrshift, "value"_a, "Allows `value >> self`");

  c.def("__rshift__", &nprshift, "value"_a, "Allows `self >> value`");

  c.def("__rsub__", &nprsub, "value"_a, "Allows `value - self`");

  c.def("__rtruediv__", &nprtruediv, "value"_a, "Allows `value / self`");

  c.def("__rxor__", &nprxor, "value"_a, "Allows `value ^ self`");

  c.def("__setitem__",
        &npsetitem,
        "key"_a,
        "value"_a,
        "Allows `self[key] = value`");

  c.def("__sub__", &npsub, "value"_a, "Allows `self - value`");

  c.def("__truediv__", &nptruediv, "value"_a, "Allows `self / value`");

  c.def("__xor__", &npxor, "value"_a, "Allows `self ^ value`");

  c.def_prop_ro("T",
                &npT,
                "Wrapper for :attr:`numpy.ndarray.T` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  // Renamed to behave as data_ as the return type is similar
  c.def_prop_ro("_base",
                &np_base,
                py::keep_alive<0, 1>{},
                "Wrapper for :attr:`numpy.ndarray.base` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  // Renamed because we often have properties with the same name
  c.def_prop_ro("_data",
                &np_data,
                py::keep_alive<0, 1>{},
                "Wrapper for :attr:`numpy.ndarray.data` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("dtype",
                &npdtype,
                "Wrapper for :attr:`numpy.ndarray.dtype` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  // Renamed to behave as data_ as the return type is similar
  c.def_prop_ro("_flat",
                &np_flat,
                py::keep_alive<0, 1>{},
                "Wrapper for :attr:`numpy.ndarray.flat` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("imag",
                &npimag,
                py::keep_alive<0, 1>{},
                "Wrapper for :attr:`numpy.ndarray.imag` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("itemsize",
                &npitemsize,
                "Wrapper for :attr:`numpy.ndarray.itemsize` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("nbytes",
                &npnbytes,
                "Wrapper for :attr:`numpy.ndarray.nbytes` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("ndim",
                &npndim,
                "Wrapper for :attr:`numpy.ndarray.ndim` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("real",
                &npreal,
                py::keep_alive<0, 1>{},
                "Wrapper for :attr:`numpy.ndarray.real` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("shape",
                &npshape,
                "Wrapper for :attr:`numpy.ndarray.shape` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("size",
                &npsize,
                "Wrapper for :attr:`numpy.ndarray.size` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def_prop_ro("strides",
                &npstrides,
                "Wrapper for :attr:`numpy.ndarray.strides` using ARTS types."
                "\n\n"
                "Use the original for greater control and more functionality.");

  c.def("flatten",
        &npflatten,
        "order"_a = 'C',
        "Wrapper for :attr:`numpy.ndarray.flatten` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("copy",
        &npcopy,
        "order"_a = 'C',
        "Wrapper for :attr:`numpy.ndarray.copy` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("fill",
        &npfill,
        "value"_a,
        "Wrapper for :attr:`numpy.ndarray.fill` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("tolist",
        &nptolist,
        "Wrapper for :attr:`numpy.ndarray.tolist` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("conj",
        &npconj,
        "Wrapper for :attr:`numpy.ndarray.conj` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("conjugate",
        &npconjugate,
        "Wrapper for :attr:`numpy.ndarray.conjugate` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("nonzero",
        &npnonzero,
        "Wrapper for :attr:`numpy.ndarray.nonzero` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("cumprod",
        &npcumprod,
        "axis"_a  = py::none(),
        "dtype"_a = py::none(),
        "out"_a   = py::none(),
        "Wrapper for :attr:`numpy.ndarray.cumprod` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("cumsum",
        &npcumsum,
        "axis"_a  = py::none(),
        "dtype"_a = py::none(),
        "out"_a   = py::none(),
        "Wrapper for :attr:`numpy.ndarray.cumsum` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("diagonal",
        &npdiagonal,
        "offset"_a = 0,
        "axis1"_a  = 0,
        "axis2"_a  = 1,
        "Wrapper for :attr:`numpy.ndarray.diagonal` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("max",
        &npmax,
        "axis"_a     = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "Wrapper for :attr:`numpy.ndarray.max` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("mean",
        &npmean,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "Wrapper for :attr:`numpy.ndarray.mean` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("min",
        &npmin,
        "axis"_a     = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "Wrapper for :attr:`numpy.ndarray.min` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("prod",
        &npprod,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "Wrapper for :attr:`numpy.ndarray.prod` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("ravel",
        &npravel,
        "order"_a,
        "Wrapper for :attr:`numpy.ndarray.ravel` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("repeat",
        &nprepeat,
        "repeats"_a,
        "axis"_a = py::none(),
        "Wrapper for :attr:`numpy.ndarray.repeat` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("reshape",
        &npreshape,
        "shape"_a,
        "order"_a = 'C',
        "Wrapper for :attr:`numpy.ndarray.reshape` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("round",
        &npround,
        "decimals"_a = 0,
        "out"_a      = py::none(),
        "Wrapper for :attr:`numpy.ndarray.round` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("squeeze",
        &npsqueeze,
        "axis"_a = py::none(),
        "Wrapper for :attr:`numpy.ndarray.squeeze` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("std",
        &npstd,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "ddof"_a     = 0,
        "keepdims"_a = false,
        "Wrapper for :attr:`numpy.ndarray.std` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("sum",
        &npsum,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "Wrapper for :attr:`numpy.ndarray.sum` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("trace",
        &nptrace,
        "offset"_a = 0,
        "axis1"_a  = 0,
        "axis2"_a  = 1,
        "dtype"_a  = py::none(),
        "out"_a    = py::none(),
        "Wrapper for :attr:`numpy.ndarray.trace` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("transpose",
        &nptranspose,
        "axes"_a = py::none(),
        "Wrapper for :attr:`numpy.ndarray.transpose` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");

  c.def("var",
        &npvar,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "ddof"_a     = 0,
        "keepdims"_a = false,
        "Wrapper for :attr:`numpy.ndarray.var` using ARTS types."
        "\n\n"
        "Use the original for greater control and more functionality.");
}
}  // namespace Python
