#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

namespace Python {
namespace py = nanobind;
using namespace py::literals;

py::object npabs(py::object& self) {
  return self.attr("value").attr("__abs__")();
}

py::object npadd(py::object& self, py::object& value_) {
  return self.attr("value").attr("__add__")(value_);
}

py::object npand(py::object& self, py::object& value_) {
  return self.attr("value").attr("__and__")(value_);
}

py::object npcontains(py::object& self, py::object& key_) {
  return self.attr("value").attr("__contains__")(key_);
}

py::object npdivmod(py::object& self, py::object& value_) {
  return self.attr("value").attr("__divmod__")(value_);
}

py::object npeq(py::object& self, py::object& value_) {
  return self.attr("value").attr("__eq__")(value_);
}

py::object npfloat(py::object& self) {
  return self.attr("value").attr("__float__")();
}

py::object npfloordiv(py::object& self, py::object& value_) {
  return self.attr("value").attr("__floordiv__")(value_);
}

py::object npge(py::object& self, py::object& value_) {
  return self.attr("value").attr("__ge__")(value_);
}

py::object npgetitem(py::object& self, py::object& key_) {
  return self.attr("value").attr("__getitem__")(key_);
}

py::object npgt(py::object& self, py::object& value_) {
  return self.attr("value").attr("__gt__")(value_);
}

py::object npiadd(py::object& self, py::object& value_) {
  self.attr("value").attr("__iadd__")(value_);
  return self;
}

py::object npiand(py::object& self, py::object& value_) {
  self.attr("value").attr("__iand__")(value_);
  return self;
}

py::object npifloordiv(py::object& self, py::object& value_) {
  self.attr("value").attr("__ifloordiv__")(value_);
  return self;
}

py::object npilshift(py::object& self, py::object& value_) {
  self.attr("value").attr("__ilshift__")(value_);
  return self;
}

py::object npimatmul(py::object& self, py::object& value_) {
  self.attr("value").attr("__imatmul__")(value_);
  return self;
}

py::object npimod(py::object& self, py::object& value_) {
  self.attr("value").attr("__imod__")(value_);
  return self;
}

py::object npimul(py::object& self, py::object& value_) {
  self.attr("value").attr("__imul__")(value_);
  return self;
}

py::object npindex(py::object& self) {
  return self.attr("value").attr("__index__")();
}

py::object npint(py::object& self) {
  return self.attr("value").attr("__int__")();
}

py::object npinvert(py::object& self) {
  return self.attr("value").attr("__invert__")();
}

py::object npior(py::object& self, py::object& value_) {
  self.attr("value").attr("__ior__")(value_);
  return self;
}

py::object npipow(py::object& self, py::object& value_) {
  self.attr("value").attr("__ipow__")(value_);
  return self;
}

py::object npirshift(py::object& self, py::object& value_) {
  self.attr("value").attr("__irshift__")(value_);
  return self;
}

py::object npisub(py::object& self, py::object& value_) {
  self.attr("value").attr("__isub__")(value_);
  return self;
}

py::object npiter(py::object& self) {
  return self.attr("value").attr("__iter__")();
}

py::object npitruediv(py::object& self, py::object& value_) {
  self.attr("value").attr("__itruediv__")(value_);
  return self;
}

py::object npixor(py::object& self, py::object& value_) {
  self.attr("value").attr("__ixor__")(value_);
  return self;
}

py::object nple(py::object& self, py::object& value_) {
  return self.attr("value").attr("__le__")(value_);
}

py::object nplen(py::object& self) {
  return self.attr("value").attr("__len__")();
}

py::object nplshift(py::object& self, py::object& value_) {
  return self.attr("value").attr("__lshift__")(value_);
}
py::object nplt(py::object& self, py::object& value_) {
  return self.attr("value").attr("__lt__")(value_);
}

py::object npmatmul(py::object& self, py::object& value_) {
  return self.attr("value").attr("__matmul__")(value_);
}

py::object npmod(py::object& self, py::object& value_) {
  return self.attr("value").attr("__mod__")(value_);
}

py::object npmul(py::object& self, py::object& value_) {
  return self.attr("value").attr("__mul__")(value_);
}

py::object npne(py::object& self, py::object& value_) {
  return self.attr("value").attr("__ne__")(value_);
}

py::object npneg(py::object& self) {
  return self.attr("value").attr("__neg__")();
}

py::object npor(py::object& self, py::object& value_) {
  return self.attr("value").attr("__or__")(value_);
}

py::object nppos(py::object& self) {
  return self.attr("value").attr("__pos__")();
}

py::object nppow(py::object& self, py::object& value_, py::object& mod_) {
  return self.attr("value").attr("__pow__")(value_, mod_);
}

py::object npradd(py::object& self, py::object& value_) {
  return self.attr("value").attr("__radd__")(value_);
}

py::object nprand(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rand__")(value_);
}

py::object nprdivmod(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rdivmod__")(value_);
}

py::object nprfloordiv(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rfloordiv__")(value_);
}

py::object nprlshift(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rlshift__")(value_);
}

py::object nprmatmul(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rmatmul__")(value_);
}

py::object nprmod(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rmod__")(value_);
}

py::object nprmul(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rmul__")(value_);
}

py::object npror(py::object& self, py::object& value_) {
  return self.attr("value").attr("__ror__")(value_);
}

py::object nprpow(py::object& self, py::object& value_, py::object& mod_) {
  return self.attr("value").attr("__rpow__")(value_, mod_);
}

py::object nprrshift(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rrshift__")(value_);
}

py::object nprshift(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rshift__")(value_);
}

py::object nprsub(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rsub__")(value_);
}

py::object nprtruediv(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rtruediv__")(value_);
}

py::object nprxor(py::object& self, py::object& value_) {
  return self.attr("value").attr("__rxor__")(value_);
}

py::object npsetitem(py::object& self, py::object& key_, py::object& value_) {
  return self.attr("value").attr("__setitem__")(key_, value_);
}

py::object npsub(py::object& self, py::object& value_) {
  return self.attr("value").attr("__sub__")(value_);
}

py::object nptruediv(py::object& self, py::object& value_) {
  return self.attr("value").attr("__truediv__")(value_);
}

py::object npxor(py::object& self, py::object& value_) {
  return self.attr("value").attr("__xor__")(value_);
}

py::object npT(py::object& self) {
  return py::object(self.attr("value").attr("T"));
}

// Renamed to behave as data_ as the return type is similar
py::object np_base(py::object& self) {
  return py::object(self.attr("value").attr("base"));
}

// Renamed because we often have properties with the same name
py::object np_data(py::object& self) {
  return py::object(self.attr("value").attr("data"));
}

py::object npdtype(py::object& self) {
  return py::object(self.attr("value").attr("dtype"));
}

// Renamed to behave as data_ as the return type is similar
py::object np_flat(py::object& self) {
  return py::object(self.attr("value").attr("flat"));
}

py::object npimag(py::object& self) {
  return py::object(self.attr("value").attr("imag"));
}

py::object npitemsize(py::object& self) {
  return py::object(self.attr("value").attr("itemsize"));
}

py::object npnbytes(py::object& self) {
  return py::object(self.attr("value").attr("nbytes"));
}

py::object npndim(py::object& self) {
  return py::object(self.attr("value").attr("ndim"));
}

py::object npreal(py::object& self) {
  return py::object(self.attr("value").attr("real"));
}

py::object npshape(py::object& self) {
  return py::object(self.attr("value").attr("shape"));
}

py::object npsize(py::object& self) {
  return py::object(self.attr("value").attr("size"));
}

py::object npstrides(py::object& self) {
  return py::object(self.attr("value").attr("strides"));
}

py::object npflatten(py::object& self, char order) {
  return self.attr("value").attr("flatten")(order);
}

py::object npcopy(py::object& self, char order) {
  return self.attr("value").attr("copy")(order);
}

py::object npfill(py::object& self, py::object& value) {
  return self.attr("value").attr("fill")(value);
}

py::object nptolist(py::object& self) {
  return self.attr("value").attr("tolist")();
}

py::object npconj(py::object& self) {
  return self.attr("value").attr("conj")();
}

py::object npconjugate(py::object& self) {
  return self.attr("value").attr("conjugate")();
}

py::object npnonzero(py::object& self) {
  return self.attr("value").attr("nonzero")();
}

py::object npcumprod(py::object& self,
                     py::object& axis,
                     py::object& dtype,
                     py::object& out) {
  return self.attr("value").attr("cumprod")(axis, dtype, out);
}

py::object npcumsum(py::object& self,
                    py::object& axis,
                    py::object& dtype,
                    py::object& out) {
  return self.attr("value").attr("cumsum")(axis, dtype, out);
}

py::object npdiagonal(py::object& self,
                      py::object& offset,
                      py::object& axis1,
                      py::object& axis2) {
  return self.attr("value").attr("diagonal")(offset, axis1, axis2);
}

py::object npmax(py::object& self,
                 py::object& axis,
                 py::object& out,
                 py::object& keepdims) {
  return self.attr("value").attr("max")(axis, out, keepdims);
}

py::object npmean(py::object& self,
                  py::object& axis,
                  py::object& dtype,
                  py::object& out,
                  py::object& keepdims) {
  return self.attr("value").attr("mean")(axis, dtype, out, keepdims);
}

py::object npmin(py::object& self,
                 py::object& axis,
                 py::object& out,
                 py::object& keepdims) {
  return self.attr("value").attr("min")(axis, out, keepdims);
}

py::object npprod(py::object& self,
                  py::object& axis,
                  py::object& dtype,
                  py::object& out,
                  py::object& keepdims) {
  return self.attr("value").attr("prod")(axis, dtype, out, keepdims);
}

py::object npravel(py::object& self, py::list& order) {
  return self.attr("value").attr("ravel")(order);
}

py::object nprepeat(py::object& self, py::object& repeats, py::object& axis) {
  return self.attr("value").attr("repeat")(repeats, axis);
}

py::object npreshape(py::object& self, py::object& shape, py::object& order) {
  return self.attr("value").attr("reshape")(shape, order);
}

py::object npround(py::object& self, py::object& decimals, py::object& out) {
  return self.attr("value").attr("round")(decimals, out);
}

py::object npsqueeze(py::object& self, py::object& axis) {
  return self.attr("value").attr("squeeze")(axis);
}

py::object npstd(py::object& self,
                 py::object& axis,
                 py::object& dtype,
                 py::object& out,
                 py::object& ddof,
                 py::object& keepdims) {
  return self.attr("value").attr("std")(axis, dtype, out, ddof, keepdims);
}

py::object npsum(py::object& self,
                 py::object& axis,
                 py::object& dtype,
                 py::object& out,
                 py::object& keepdims) {
  return self.attr("value").attr("sum")(axis, dtype, out, keepdims);
}

py::object nptrace(py::object& self,
                   py::object& offset,
                   py::object& axis1,
                   py::object& axis2,
                   py::object& dtype,
                   py::object& out) {
  return self.attr("value").attr("trace")(offset, axis1, axis2, dtype, out);
}

py::object nptranspose(py::object& self, py::list& axes) {
  return self.attr("value").attr("transpose")(axes);
}

py::object npvar(py::object& self,
                 py::object& axis,
                 py::object& dtype,
                 py::object& out,
                 py::object& ddof,
                 py::object& keepdims) {
  return self.attr("value").attr("var")(axis, dtype, out, ddof, keepdims);
}
}  // namespace Python
