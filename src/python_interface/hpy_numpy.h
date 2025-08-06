#pragma once

#include "hpy_numpy_impl.h"

namespace Python {
template <typename FUNC, typename... in>
auto vectorize(FUNC f, in&&... args) {
  auto np = py::module_::import_("numpy");
  return np.attr("vectorize")(py::cpp_function(f))(std::forward<in>(args)...);
}

void common_ndarray(auto& c) {
  c.def("__abs__", &npabs, "See :attr:`numpy.ndarray.__abs__`");

  c.def("__add__", &npadd, "value"_a, "See :attr:`numpy.ndarray.__add__`");

  c.def("__and__", &npand, "value"_a, "See :attr:`numpy.ndarray.__and__`");

  c.def("__contains__",
        &npcontains,
        "key"_a,
        "See :attr:`numpy.ndarray.__contains__`");

  c.def("__divmod__",
        &npdivmod,
        "value"_a,
        "See :attr:`numpy.ndarray.__divmod__`");

  c.def("__eq__", &npeq, "value"_a, "See :attr:`numpy.ndarray.__eq__`");

  c.def("__float__", &npfloat, "See :attr:`numpy.ndarray.__float__`");

  c.def("__floordiv__",
        &npfloordiv,
        "value"_a,
        "See :attr:`numpy.ndarray.__floordiv__`");

  c.def("__ge__", &npge, "value"_a, "See :attr:`numpy.ndarray.__ge__`");

  c.def("__getitem__",
        &npgetitem,
        "key"_a,
        "See :attr:`numpy.ndarray.__getitem__`");

  c.def("__gt__", &npgt, "value"_a, "See :attr:`numpy.ndarray.__gt__`");

  c.def("__iadd__", &npiadd, "value"_a, "See :attr:`numpy.ndarray.__iadd__`");

  c.def("__iand__", &npiand, "value"_a, "See :attr:`numpy.ndarray.__iand__`");

  c.def("__ifloordiv__",
        &npifloordiv,
        "value"_a,
        "See :attr:`numpy.ndarray.__ifloordiv__`");

  c.def("__ilshift__",
        &npilshift,
        "value"_a,
        "See :attr:`numpy.ndarray.__ilshift__`");

  c.def("__imatmul__",
        &npimatmul,
        "value"_a,
        "See :attr:`numpy.ndarray.__imatmul__`");

  c.def("__imod__", &npimod, "value"_a, "See :attr:`numpy.ndarray.__imod__`");

  c.def("__imul__", &npimul, "value"_a, "See :attr:`numpy.ndarray.__imul__`");

  c.def("__index__", &npindex, "See :attr:`numpy.ndarray.__index__`");

  c.def("__int__", &npint, "See :attr:`numpy.ndarray.__int__`");

  c.def("__invert__", &npinvert, "See :attr:`numpy.ndarray.__invert__`");

  c.def("__ior__", &npior, "value"_a, "See :attr:`numpy.ndarray.__ior__`");

  c.def("__ipow__", &npipow, "value"_a, "See :attr:`numpy.ndarray.__ipow__`");

  c.def("__irshift__",
        &npirshift,
        "value"_a,
        "See :attr:`numpy.ndarray.__irshift__`");

  c.def("__isub__",
        &npisub,
        py::rv_policy::reference_internal,
        "value"_a,
        "See :attr:`numpy.ndarray.__isub__`");

  c.def("__iter__", &npiter, "See :attr:`numpy.ndarray.__iter__`");

  c.def("__itruediv__",
        &npitruediv,
        "value"_a,
        "See :attr:`numpy.ndarray.__itruediv__`");

  c.def("__ixor__", &npixor, "value"_a, "See :attr:`numpy.ndarray.__ixor__`");

  c.def("__le__", &nple, "value"_a, "See :attr:`numpy.ndarray.__le__`");

  c.def("__len__", &nplen, "See :attr:`numpy.ndarray.__len__`");

  c.def("__lshift__",
        &nplshift,
        "value"_a,
        "See :attr:`numpy.ndarray.__lshift__`");

  c.def("__lt__", &nplt, "value"_a, "See :attr:`numpy.ndarray.__lt__`");

  c.def("__matmul__",
        &npmatmul,
        "value"_a,
        "See :attr:`numpy.ndarray.__matmul__`");

  c.def("__mod__", &npmod, "value"_a, "See :attr:`numpy.ndarray.__mod__`");

  c.def("__mul__", &npmul, "value"_a, "See :attr:`numpy.ndarray.__mul__`");

  c.def("__ne__", &npne, "value"_a, "See :attr:`numpy.ndarray.__ne__`");

  c.def("__neg__", &npneg, "See :attr:`numpy.ndarray.__neg__`");

  c.def("__or__", &npor, "value"_a, "See :attr:`numpy.ndarray.__or__`");

  c.def("__pos__", &nppos, "See :attr:`numpy.ndarray.__pos__`");

  c.def("__pow__",
        &nppow,
        "value"_a,
        "mod"_a = py::none(),
        "See :attr:`numpy.ndarray.__pow__`");

  c.def("__radd__", &npradd, "value"_a, "See :attr:`numpy.ndarray.__radd__`");

  c.def("__rand__", &nprand, "value"_a, "See :attr:`numpy.ndarray.__rand__`");

  c.def("__rdivmod__",
        &nprdivmod,
        "value"_a,
        "See :attr:`numpy.ndarray.__rdivmod__`");

  c.def("__rfloordiv__",
        &nprfloordiv,
        "value"_a,
        "See :attr:`numpy.ndarray.__rfloordiv__`");

  c.def("__rlshift__",
        &nprlshift,
        "value"_a,
        "See :attr:`numpy.ndarray.__rlshift__`");

  c.def("__rmatmul__",
        &nprmatmul,
        "value"_a,
        "See :attr:`numpy.ndarray.__rmatmul__`");

  c.def("__rmod__", &nprmod, "value"_a, "See :attr:`numpy.ndarray.__rmod__`");

  c.def("__rmul__", &nprmul, "value"_a, "See :attr:`numpy.ndarray.__rmul__`");

  c.def("__ror__", &npror, "value"_a, "See :attr:`numpy.ndarray.__ror__`");

  c.def("__rpow__",
        &nprpow,
        "value"_a,
        "mod"_a = py::none(),
        "See :attr:`numpy.ndarray.__rpow__`");

  c.def("__rrshift__",
        &nprrshift,
        "value"_a,
        "See :attr:`numpy.ndarray.__rrshift__`");

  c.def("__rshift__",
        &nprshift,
        "value"_a,
        "See :attr:`numpy.ndarray.__rshift__`");

  c.def("__rsub__", &nprsub, "value"_a, "See :attr:`numpy.ndarray.__rsub__`");

  c.def("__rtruediv__",
        &nprtruediv,
        "value"_a,
        "See :attr:`numpy.ndarray.__rtruediv__`");

  c.def("__rxor__", &nprxor, "value"_a, "See :attr:`numpy.ndarray.__rxor__`");

  c.def("__setitem__",
        &npsetitem,
        "key"_a,
        "value"_a,
        "See :attr:`numpy.ndarray.__setitem__`");

  c.def("__sub__", &npsub, "value"_a, "See :attr:`numpy.ndarray.__sub__`");

  c.def("__truediv__",
        &nptruediv,
        "value"_a,
        "See :attr:`numpy.ndarray.__truediv__`");

  c.def("__xor__", &npxor, "value"_a, "See :attr:`numpy.ndarray.__xor__`");

  c.def_prop_ro("T", &npT, "See :attr:`numpy.ndarray.T`");

  // Renamed to behave as data_ as the return type is similar
  c.def_prop_ro("_base",
                &np_base,
                py::keep_alive<0, 1>{},
                "See :attr:`numpy.ndarray.base`");

  // Renamed because we often have properties with the same name
  c.def_prop_ro("_data",
                &np_data,
                py::keep_alive<0, 1>{},
                "See :attr:`numpy.ndarray.data`");

  c.def_prop_ro("dtype", &npdtype, "See :attr:`numpy.ndarray.dtype`");

  // Renamed to behave as data_ as the return type is similar
  c.def_prop_ro("_flat",
                &np_flat,
                py::keep_alive<0, 1>{},
                "See :attr:`numpy.ndarray.flat`");

  c.def_prop_ro("imag",
                &npimag,
                py::keep_alive<0, 1>{},
                "See :attr:`numpy.ndarray.imag`");

  c.def_prop_ro("itemsize", &npitemsize, "See :attr:`numpy.ndarray.itemsize`");

  c.def_prop_ro("nbytes", &npnbytes, "See :attr:`numpy.ndarray.nbytes`");

  c.def_prop_ro("ndim", &npndim, "See :attr:`numpy.ndarray.ndim`");

  c.def_prop_ro("real",
                &npreal,
                py::keep_alive<0, 1>{},
                "See :attr:`numpy.ndarray.real`");

  c.def_prop_ro("shape", &npshape, "See :attr:`numpy.ndarray.shape`");

  c.def_prop_ro("size", &npsize, "See :attr:`numpy.ndarray.size`");

  c.def_prop_ro("strides", &npstrides, "See :attr:`numpy.ndarray.strides`");

  c.def("flatten",
        &npflatten,
        "order"_a = 'C',
        "See :attr:`numpy.ndarray.flatten`");

  c.def("copy", &npcopy, "order"_a = 'C', "See :attr:`numpy.ndarray.copy`");

  c.def("fill", &npfill, "value"_a, "See :attr:`numpy.ndarray.fill`");

  c.def("tolist", &nptolist, "See :attr:`numpy.ndarray.tolist`");

  c.def("conj", &npconj, "See :attr:`numpy.ndarray.conj`");

  c.def("conjugate", &npconjugate, "See :attr:`numpy.ndarray.conjugate`");

  c.def("nonzero", &npnonzero, "See :attr:`numpy.ndarray.nonzero`");

  c.def("cumprod",
        &npcumprod,
        "axis"_a  = py::none(),
        "dtype"_a = py::none(),
        "out"_a   = py::none(),
        "See :attr:`numpy.ndarray.cumprod`");

  c.def("cumsum",
        &npcumsum,
        "axis"_a  = py::none(),
        "dtype"_a = py::none(),
        "out"_a   = py::none(),
        "See :attr:`numpy.ndarray.cumsum`");

  c.def("diagonal",
        &npdiagonal,
        "offset"_a = 0,
        "axis1"_a  = 0,
        "axis2"_a  = 1,
        "See :attr:`numpy.ndarray.diagonal`");

  c.def("max",
        &npmax,
        "axis"_a     = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "See :attr:`numpy.ndarray.max`");

  c.def("mean",
        &npmean,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "See :attr:`numpy.ndarray.mean`");

  c.def("min",
        &npmin,
        "axis"_a     = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "See :attr:`numpy.ndarray.min`");

  c.def("prod",
        &npprod,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "See :attr:`numpy.ndarray.prod`");

  c.def("ravel", &npravel, "order"_a = 'C', "See :attr:`numpy.ndarray.ravel`");

  c.def("repeat",
        &nprepeat,
        "repeats"_a,
        "axis"_a = py::none(),
        "See :attr:`numpy.ndarray.repeat`");

  c.def("reshape",
        &npreshape,
        "shape"_a,
        "order"_a = 'C',
        "See :attr:`numpy.ndarray.reshape`");

  c.def("round",
        &npround,
        "decimals"_a = 0,
        "out"_a      = py::none(),
        "See :attr:`numpy.ndarray.round`");

  c.def("squeeze",
        &npsqueeze,
        "axis"_a = py::none(),
        "See :attr:`numpy.ndarray.squeeze`");

  c.def("std",
        &npstd,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "ddof"_a     = 0,
        "keepdims"_a = false,
        "See :attr:`numpy.ndarray.std`");

  c.def("sum",
        &npsum,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "keepdims"_a = false,
        "See :attr:`numpy.ndarray.sum`");

  c.def("trace",
        &nptrace,
        "offset"_a = 0,
        "axis1"_a  = 0,
        "axis2"_a  = 1,
        "dtype"_a  = py::none(),
        "out"_a    = py::none(),
        "See :attr:`numpy.ndarray.trace`");

  c.def("transpose",
        &nptranspose,
        "axes"_a,
        "See :attr:`numpy.ndarray.transpose`");

  c.def("var",
        &npvar,
        "axis"_a     = py::none(),
        "dtype"_a    = py::none(),
        "out"_a      = py::none(),
        "ddof"_a     = 0,
        "keepdims"_a = false,
        "See :attr:`numpy.ndarray.var`");
}
}  // namespace Python
