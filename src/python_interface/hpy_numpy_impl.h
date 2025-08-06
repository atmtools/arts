#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>

namespace Python {
namespace py = nanobind;
using namespace py::literals;
py::object npabs(py::object& self);
py::object npadd(py::object& self, py::object& value_);
py::object npand(py::object& self, py::object& value_);
py::object npcontains(py::object& self, py::object& key_);
py::object npdivmod(py::object& self, py::object& value_);
py::object npeq(py::object& self, py::object& value_);
py::object npfloat(py::object& self);
py::object npfloordiv(py::object& self, py::object& value_);
py::object npge(py::object& self, py::object& value_);
py::object npgetitem(py::object& self, py::object& key_);
py::object npgt(py::object& self, py::object& value_);
py::object npiadd(py::object& self, py::object& value_);
py::object npiand(py::object& self, py::object& value_);
py::object npifloordiv(py::object& self, py::object& value_);
py::object npilshift(py::object& self, py::object& value_);
py::object npimatmul(py::object& self, py::object& value_);
py::object npimod(py::object& self, py::object& value_);
py::object npimul(py::object& self, py::object& value_);
py::object npindex(py::object& self);
py::object npint(py::object& self);
py::object npinvert(py::object& self);
py::object npior(py::object& self, py::object& value_);
py::object npipow(py::object& self, py::object& value_);
py::object npirshift(py::object& self, py::object& value_);
py::object npisub(py::object& self, py::object& value_);
py::object npiter(py::object& self);
py::object npitruediv(py::object& self, py::object& value_);
py::object npixor(py::object& self, py::object& value_);
py::object nple(py::object& self, py::object& value_);
py::object nplen(py::object& self);
py::object nplshift(py::object& self, py::object& value_);
py::object nplt(py::object& self, py::object& value_);
py::object npmatmul(py::object& self, py::object& value_);
py::object npmod(py::object& self, py::object& value_);
py::object npmul(py::object& self, py::object& value_);
py::object npne(py::object& self, py::object& value_);
py::object npneg(py::object& self);
py::object npor(py::object& self, py::object& value_);
py::object nppos(py::object& self);
py::object nppow(py::object& self, py::object& value_, py::object& mod_);
py::object npradd(py::object& self, py::object& value_);
py::object nprand(py::object& self, py::object& value_);
py::object nprdivmod(py::object& self, py::object& value_);
py::object nprfloordiv(py::object& self, py::object& value_);
py::object nprlshift(py::object& self, py::object& value_);
py::object nprmatmul(py::object& self, py::object& value_);
py::object nprmod(py::object& self, py::object& value_);
py::object nprmul(py::object& self, py::object& value_);
py::object npror(py::object& self, py::object& value_);
py::object nprpow(py::object& self, py::object& value_, py::object& mod_);
py::object nprrshift(py::object& self, py::object& value_);
py::object nprshift(py::object& self, py::object& value_);
py::object nprsub(py::object& self, py::object& value_);
py::object nprtruediv(py::object& self, py::object& value_);
py::object nprxor(py::object& self, py::object& value_);
py::object npsetitem(py::object& self, py::object& key_, py::object& value_);
py::object npsub(py::object& self, py::object& value_);
py::object nptruediv(py::object& self, py::object& value_);
py::object npxor(py::object& self, py::object& value_);
py::object npT(py::object& self);
py::object np_base(py::object& self);
py::object np_data(py::object& self);
py::object npdtype(py::object& self);
py::object np_flat(py::object& self);
py::object npimag(py::object& self);
py::object npitemsize(py::object& self);
py::object npnbytes(py::object& self);
py::object npndim(py::object& self);
py::object npreal(py::object& self);
py::object npshape(py::object& self);
py::object npsize(py::object& self);
py::object npstrides(py::object& self);
py::object npflatten(py::object& self, char order);
py::object npcopy(py::object& self, char order);
py::object npfill(py::object& self, py::object& value);
py::object nptolist(py::object& self);
py::object npconj(py::object& self);
py::object npconjugate(py::object& self);
py::object npnonzero(py::object& self);
py::object npcumprod(py::object& self,
               py::object& axis,
               py::object& dtype,
               py::object& out);
py::object npcumsum(py::object& self,
              py::object& axis,
              py::object& dtype,
              py::object& out);
py::object npdiagonal(py::object& self,
                py::object& offset,
                py::object& axis1,
                py::object& axis2);
py::object npmax(py::object& self,
           py::object& axis,
           py::object& out,
           py::object& keepdims);
py::object npmean(py::object& self,
            py::object& axis,
            py::object& dtype,
            py::object& out,
            py::object& keepdims);
py::object npmin(py::object& self,
           py::object& axis,
           py::object& out,
           py::object& keepdims);
py::object npprod(py::object& self,
            py::object& axis,
            py::object& dtype,
            py::object& out,
            py::object& keepdims);
py::object npravel(py::object& self, py::object& order);
py::object nprepeat(py::object& self, py::object& repeats, py::object& axis);
py::object npreshape(py::object& self, py::object& shape, py::object& order);
py::object npround(py::object& self, py::object& decimals, py::object& out);
py::object npsqueeze(py::object& self, py::object& axis);
py::object npstd(py::object& self,
           py::object& axis,
           py::object& dtype,
           py::object& out,
           py::object& ddof,
           py::object& keepdims);
py::object npsum(py::object& self,
           py::object& axis,
           py::object& dtype,
           py::object& out,
           py::object& keepdims);
py::object nptrace(py::object& self,
             py::object& offset,
             py::object& axis1,
             py::object& axis2,
             py::object& dtype,
             py::object& out);
py::object nptranspose(py::object& self, py::args& axes);
py::object npvar(py::object& self,
           py::object& axis,
           py::object& dtype,
           py::object& out,
           py::object& ddof,
           py::object& keepdims);
}  // namespace Python
