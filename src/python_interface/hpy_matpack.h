#pragma once

#include <matpack.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/variant.h>

#include "configtypes.h"
#include "hpy_numpy.h"

namespace Python {
namespace py = nanobind;
using namespace py::literals;

template <typename mtype>
void matpack_common_interface(py::class_<mtype>& c) {
  c.def(py::init<>());

  c.def(
      "__init__",
      [](mtype* v, const py::list& l) {
        auto np = py::module_::import_("numpy");
        auto m  = py::type<mtype>()(np.attr("asarray")(l));

        new (v) mtype(py::cast<mtype>(m));
      },
      "l"_a);

  c.def_prop_rw(
      "value",
      [](py::object& x) { return x.attr("__array__")(); },
      [](mtype& a, const mtype& b) { a = b; },
      py::for_getter(py::rv_policy::reference_internal),
      "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");

  common_ndarray(c);

  py::implicitly_convertible<py::list, mtype>();
}

template <typename T, Size ndim>
void matpack_interface(py::class_<matpack::data_t<T, ndim>>& c) {
  using mtype = matpack::data_t<T, ndim>;
  using nd    = py::ndarray<py::numpy, T, py::ndim<ndim>, py::c_contig>;
  using const_nd =
      py::ndarray<py::numpy, const T, py::ndim<ndim>, py::c_contig>;

  c.def(
      "__init__",
      [](mtype* v, const const_nd& a) {
        std::array<Index, ndim> shape;
        for (Size i = 0; i < ndim; i++) {
          shape[i] = static_cast<Index>(a.shape(i));
        }

        new (v) mtype(shape);
        std::copy(a.data(), a.data() + a.size(), v->elem_begin());
      },
      "a"_a);

  c.def(
      "__init__",
      [](mtype* v, const nd& a) {
        std::array<Index, ndim> shape;
        for (Size i = 0; i < ndim; i++) {
          shape[i] = static_cast<Index>(a.shape(i));
        }

        new (v) mtype(shape);
        std::copy(a.data(), a.data() + a.size(), v->elem_begin());
      },
      "a"_a);

  c.def(
      "__array__",
      [](mtype& v,
         py::object dtype,
         py::object copy) -> std::variant<nd, py::object> {
        std::array<size_t, ndim> shape;
        std::ranges::copy(v.shape(), shape.begin());

        auto np = py::module_::import_("numpy");
        auto x  = nd(v.data_handle(), ndim, shape.data(), py::cast(&v));

        if (not dtype.is_none()) {
          return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
        }

        return x.cast((copy.is_none() or not py::bool_(copy))
                          ? py::rv_policy::automatic_reference
                          : py::rv_policy::copy);
      },
      "dtype"_a = py::none(),
      "copy"_a  = py::none(),
      "Allows :func:`~numpy.array` to be called with the object.");

  py::implicitly_convertible<nd, mtype>();
  py::implicitly_convertible<const_nd, mtype>();

  matpack_common_interface(c);
}

template <typename T, Size... ndim>
void matpack_constant_interface(py::class_<matpack::cdata_t<T, ndim...>>& c) {
  using mtype = matpack::cdata_t<T, ndim...>;
  using nd    = py::ndarray<py::numpy, T, py::shape<ndim...>, py::c_contig>;
  using const_nd =
      py::ndarray<py::numpy, const T, py::shape<ndim...>, py::c_contig>;

  c.def(
      "__init__",
      [](mtype* v, const nd& a) {
        new (v) mtype();
        std::copy(a.data(), a.data() + a.size(), v->elem_begin());
      },
      "a"_a);

  c.def(
      "__init__",
      [](mtype* v, const const_nd& a) {
        new (v) mtype();
        std::copy(a.data(), a.data() + a.size(), v->elem_begin());
      },
      "a"_a);

  c.def_rw("data", &mtype::data, "The data itself\n\n.. :class:`object`");

  c.def(
      "norm",
      [](const mtype& v) { return v / hypot(v); },
      "Return normalized version of self");

  c.def(
      "__array__",
      [](mtype& v,
         py::object dtype,
         py::object copy) -> std::variant<nd, py::object> {
        constexpr std::array<size_t, sizeof...(ndim)> shape{
            static_cast<size_t>(ndim)...};

        auto np = py::module_::import_("numpy");
        auto x = nd(v.data.data(), sizeof...(ndim), shape.data(), py::cast(&v));

        if (not dtype.is_none()) {
          return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
        }

        return x.cast((copy.is_none() or not py::bool_(copy))
                          ? py::rv_policy::automatic_reference
                          : py::rv_policy::copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none(),
      "Allows :func:`~numpy.array` to be called with the object.");

  py::implicitly_convertible<nd, mtype>();
  py::implicitly_convertible<const_nd, mtype>();

  matpack_common_interface(c);
}

template <class Compare>
void matpack_grid_interface(py::class_<matpack::grid_t<Compare>>& c) {
  using mtype = Vector;
  using nd = py::ndarray<py::numpy, const Numeric, py::ndim<1>, py::c_contig>;
  using mut_nd = py::ndarray<py::numpy, Numeric, py::ndim<1>, py::c_contig>;

  c.def(
      "__init__",
      [](matpack::grid_t<Compare>* v, const nd& a) {
        auto m = py::type<Vector>()(a);
        new (v) matpack::grid_t<Compare>(py::cast<Vector>(m));
      },
      "a"_a);

  c.def(
      "__init__",
      [](matpack::grid_t<Compare>* v, const mut_nd& a) {
        auto m = py::type<Vector>()(a);
        new (v) matpack::grid_t<Compare>(py::cast<Vector>(m));
      },
      "a"_a);

  c.def(
      "__init__",
      [](matpack::grid_t<Compare>* v, const mtype& a) {
        new (v) matpack::grid_t<Compare>(a);
      },
      "vec"_a);

  c.def(
      "__array__",
      [](matpack::grid_t<Compare>& v,
         py::object dtype,
         const py::object& copy) -> std::variant<nd, py::object> {
        std::array<size_t, 1> shape{static_cast<size_t>(v.size())};

        auto np = py::module_::import_("numpy");
        auto x  = nd(v.vec().data_handle(), 1, shape.data(), py::cast(&v));

        if (not dtype.is_none()) {
          return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
        }

        return x.cast((copy.is_none() or not py::bool_(copy))
                          ? py::rv_policy::automatic_reference
                          : py::rv_policy::copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none(),
      "Allows :func:`~numpy.array` to be called with the object.");

  py::implicitly_convertible<nd, matpack::grid_t<Compare>>();
  py::implicitly_convertible<mut_nd, matpack::grid_t<Compare>>();
  py::implicitly_convertible<mtype, matpack::grid_t<Compare>>();

  matpack_common_interface(c);
}

template <class Compare, Numeric l, Numeric u, bool il, bool iu>
void matpack_grid_interface(
    py::class_<matpack::ranged_grid_t<Compare, l, u, il, iu>>& c) {
  using mtype = Vector;
  using nd = py::ndarray<py::numpy, const Numeric, py::ndim<1>, py::c_contig>;
  using mut_nd = py::ndarray<py::numpy, Numeric, py::ndim<1>, py::c_contig>;
  using T      = matpack::ranged_grid_t<Compare, l, u, il, iu>;

  c.def(
      "__init__",
      [](T* v, const mut_nd& a) {
        auto m = py::type<Vector>()(a);
        new (v) T(py::cast<Vector>(m));
      },
      "a"_a);

  c.def(
      "__init__",
      [](T* v, const nd& a) {
        auto m = py::type<Vector>()(a);
        new (v) T(py::cast<Vector>(m));
      },
      "a"_a);

  c.def("__init__", [](T* v, const mtype& a) { new (v) T(a); }, "vec"_a);

  c.def(
      "__array__",
      [](T& v,
         py::object dtype,
         const py::object& copy) -> std::variant<nd, py::object> {
        std::array<size_t, 1> shape{static_cast<size_t>(v.size())};

        auto np = py::module_::import_("numpy");
        auto x  = nd(v.vec().data_handle(), 1, shape.data(), py::cast(&v));

        if (not dtype.is_none()) {
          return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
        }

        return x.cast((copy.is_none() or not py::bool_(copy))
                          ? py::rv_policy::automatic_reference
                          : py::rv_policy::copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none(),
      "Allows :func:`~numpy.array` to be called with the object.");

  py::implicitly_convertible<nd, T>();
  py::implicitly_convertible<mut_nd, T>();
  py::implicitly_convertible<mtype, T>();

  matpack_common_interface(c);
}

template <typename T, typename... Grids>
void gridded_data_interface(
    py::class_<matpack::gridded_data_t<T, Grids...>>& c) {
  constexpr Size dim = sizeof...(Grids);
  using mtype        = matpack::gridded_data_t<T, Grids...>;
  using dtype        = matpack::data_t<T, dim>;

  c.def(py::init<>());

  c.def(py::init<String,
                 dtype,
                 std::array<String, dim>,
                 typename mtype::grids_t>(),
        "name"_a = String{},
        "data"_a,
        "grid_names"_a = std::array<String, dim>{},
        "grids"_a);

  c.def_rw(
      "dataname", &mtype::data_name, "Name of the data\n\n.. :class:`String`");

  c.def_rw("data", &mtype::data, "The data itself\n\n.. :class:`object`");

  c.def_rw("gridnames",
           &mtype::grid_names,
           "The grid names\n\n.. :class:`list[str]`");

  c.def_rw("grids",
           &mtype::grids,
           "The grid of the data\n\n.. :class:`list[object]`");

  c.def_prop_ro_static(
      "dim",
      [](const py::object&) { return mtype::dim; },
      "Dimension of the field\n\n.. :class:`int`");

  c.def("ok", &mtype::ok, "Check the field");

  c.def(
      "__array__",
      [](py::object& gd, py::object dtype_, py::object copy) {
        return gd.attr("data").attr("__array__")("dtype"_a = dtype_,
                                                 "copy"_a  = copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none(),
      "Allows :func:`~numpy.array` to be called with the object.");

  c.def_prop_rw(
      "value",
      [](py::object& x) { return x.attr("__array__")("copy"_a = false); },
      [](mtype& a, const mtype& b) { a = b; },
      py::for_getter(py::rv_policy::reference_internal),
      "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`numpy.ndarray`");

  c.def(
      "to_dict",
      [](const py::object& gd) {
        py::dict out;

        out["dims"]   = py::list{};
        out["coords"] = py::dict{};
        for (Size i = 0; i < mtype::dim; ++i) {
          auto gridname = gd.attr("gridnames").attr("__getitem__")(i);
          auto grid     = gd.attr("grids").attr("__getitem__")(i);

          const auto n =
              gridname ? gridname : py::str(std::format("dim{}", i).c_str());
          out["dims"].attr("append")(n);
          out["coords"][n]         = py::dict{};
          out["coords"][n]["data"] = grid;
          out["coords"][n]["dims"] = n;
        }

        out["name"] = gd.attr("dataname");
        out["data"] = gd.attr("data");

        return out;
      },
      "Convert object to :class:`dict`.");

  c.def(
      "to_xarray",
      [](const py::object& gd) {
        py::module_ xarray = py::module_::import_("xarray");
        return xarray.attr("DataArray").attr("from_dict")(gd.attr("to_dict")());
      },
      "Create :class:`xarray.DataArray` from the object.");

  c.def_static(
      "from_dict",
      [](const py::dict& d) {
        if (not d.contains("coords"))
          throw std::invalid_argument(
              "cannot convert dict without the key 'coords''");

        if (not d.contains("dims"))
          throw std::invalid_argument(
              "cannot convert dict without the key 'dims''");

        if (not d.contains("data"))
          throw std::invalid_argument(
              "cannot convert dict without the key 'data''");

        auto gd = mtype{};

        gd.data = py::cast<dtype>(d["data"]);

        if (d.contains("name")) py::try_cast(d["name"], gd.data_name);

        gd.grid_names =
            py::cast<std::array<std::string, mtype::dim>>(d["dims"]);

        py::list coords{};
        for (auto& n : gd.grid_names) {
          coords.append(d["coords"][py::str(n.c_str())]["data"]);
        }

        gd.grids = py::cast<typename mtype::grids_t>(coords);

        return gd;
      },
      "Create object from :class:`dict`.");

  c.def_static(
      "from_xarray",
      [](const py::object& xarray) {
        return py::type<mtype>().attr("from_dict")(xarray.attr("to_dict")());
      },
      "Create object from :class:`xarray.DataArray`.");

  c.def(
      "__eq__",
      [](const mtype& a, const mtype& b) { return a == b; },
      py::is_operator(),
      "value"_a,
      "Allows `self == value`");

  c.def(
      "__ne__",
      [](const mtype& a, const mtype& b) { return a != b; },
      py::is_operator(),
      "value"_a,
      "Allows `self != value`");

  common_ndarray(c);

  c.def("__init__", [](mtype* v, const py::dict& d) {
    auto m = py::type<mtype>().attr("from_dict")(d);
    new (v) mtype(py::cast<mtype>(m));
  });

  py::implicitly_convertible<py::dict, mtype>();
}
}  // namespace Python
