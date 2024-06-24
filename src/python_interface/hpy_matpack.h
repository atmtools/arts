#pragma once

#include <matpack.h>

#include "hpy_numpy.h"

#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

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
        auto m =
            py::type<mtype>()(np.attr("array")(l, py::arg("copy") = false));

        new (v) mtype(py::cast<mtype>(m));
      },
      "l"_a);

  c.def("__str__", [](const mtype& v) { return var_string(v); });

  c.def("__repr__", [](const mtype& v) { return var_string(v); });

  c.def("__iter__",
        [](py::object& v) { return v.attr("value").attr("__iter__")(); });

  c.def_prop_rw(
      "value",
      [](py::object& v) {
        auto np = py::module_::import_("numpy");
        return np.attr("array")(v, py::arg("copy") = false);
      },
      [](mtype& a, const mtype& b) { a = b; });

  c.def("__getstate__",
        [](const py::object& v) { return std::tuple{v.attr("__array__")()}; });

  c.def("__setstate__", [](mtype& v, const std::tuple<py::object>& x) {
    auto m = py::type<mtype>()(std::get<0>(x));
    new (&v) mtype{py::cast<mtype>(m)};
  });

  common_ndarray(c);

  py::implicitly_convertible<py::list, mtype>();
}

template <typename T, Index ndim>
void matpack_interface(py::class_<matpack::matpack_data<T, ndim>>& c) {
  using mtype = matpack::matpack_data<T, ndim>;
  using nd    = py::ndarray<py::numpy, T, py::ndim<ndim>, py::c_contig>;

  c.def(
      "__init__",
      [](mtype* v, const nd& a) {
        std::array<Index, ndim> shape;
        for (Index i = 0; i < ndim; i++) {
          shape[i] = static_cast<Index>(a.shape(i));
        }

        new (v) mtype(shape);
        std::copy(a.data(), a.data() + a.size(), v->elem_begin());
      },
      "a"_a);

  c.def(
      "__array__",
      [](mtype& v) {
        std::array<size_t, ndim> shape;
        std::ranges::copy(v.shape(), shape.begin());

        return nd(v.data_handle(), ndim, shape.data(), py::handle());
      },
      py::rv_policy::reference_internal);

  py::implicitly_convertible<nd, mtype>();

  matpack_common_interface(c);
}

template <typename T, Index... ndim>
void matpack_constant_interface(
    py::class_<matpack::matpack_constant_data<T, ndim...>>& c) {
  using mtype = matpack::matpack_constant_data<T, ndim...>;
  using nd    = py::ndarray<py::numpy, T, py::shape<ndim...>, py::c_contig>;

  c.def(
      "__init__",
      [](mtype* v, const nd& a) {
        new (v) mtype();
        std::copy(a.data(), a.data() + a.size(), v->elem_begin());
      },
      "a"_a);

  c.def_rw("data", &mtype::data);

  c.def(
      "__array__",
      [](mtype& v) {
        constexpr std::array<size_t, sizeof...(ndim)> shape{
            static_cast<size_t>(ndim)...};

        return nd(v.data.data(), sizeof...(ndim), shape.data(), py::handle());
      },
      py::rv_policy::reference_internal);

  py::implicitly_convertible<nd, mtype>();

  matpack_common_interface(c);
}

template <class Compare>
void matpack_grid_interface(py::class_<matpack::grid<Compare>>& c) {
  using mtype = Vector;
  using nd = py::ndarray<py::numpy, const Numeric, py::ndim<1>, py::c_contig>;

  c.def(
      "__init__",
      [](matpack::grid<Compare>* v, const nd& a) {
        auto m = py::type<Vector>()(a);
        new (v) matpack::grid<Compare>(py::cast<Vector>(m));
      },
      "a"_a);

  c.def(
      "__init__",
      [](matpack::grid<Compare>* v, const mtype& a) {
        new (v) matpack::grid<Compare>(a);
      },
      "vec"_a);

  c.def(
      "__array__",
      [](matpack::grid<Compare>& v) {
        std::array<size_t, 1> shape{static_cast<size_t>(v.size())};

        return nd(v.data_handle(), 1, shape.data(), py::handle());
      },
      py::rv_policy::reference_internal);

  py::implicitly_convertible<nd, matpack::grid<Compare>>();
  py::implicitly_convertible<mtype, matpack::grid<Compare>>();

  matpack_common_interface(c);
}

template <typename T, typename... Grids>
void gridded_data_interface(py::class_<matpack::gridded_data<T, Grids...>>& c) {
  constexpr Index dim = sizeof...(Grids);
  using mtype         = matpack::gridded_data<T, Grids...>;
  using dtype         = matpack::matpack_data<T, dim>;

  c.def(py::init<>());

  c.def(py::init<String,
                 dtype,
                 std::array<String, dim>,
                 typename mtype::grids_t>(),
        "name"_a = String{},
        "data"_a,
        "grid_names"_a = std::array<String, dim>{},
        "grids"_a);

  c.def_rw("dataname", &mtype::data_name);

  c.def_rw("data", &mtype::data);

  c.def_rw("gridnames", &mtype::grid_names);

  c.def_rw("grids", &mtype::grids);

  c.def_prop_ro_static(
      "dim",
      [](const py::object&) { return mtype::dim; },
      "Dimension of the field");

  c.def("ok", &mtype::ok, "Check the field");

  c.def("__array__", [](mtype& gd) {
    auto np = py::module_::import_("numpy");
    return np.attr("array")(gd.data, py::arg("copy") = false);
  });

  c.def_prop_rw(
      "value",
      [](py::object& v) {
        auto np = py::module_::import_("numpy");
        return np.attr("array")(v, py::arg("copy") = false);
      },
      [](mtype& a, const mtype& b) { a = b; });

  c.def("__repr__", [](const mtype& m) { return var_string(m); });

  c.def("__str__", [](const mtype& m) { return var_string(m); });

  c.def("to_dict", [](const py::object& gd) {
    py::dict out;

    out["dims"]   = py::list{};
    out["coords"] = py::dict{};
    for (Size i = 0; i < mtype::dim; ++i) {
      auto gridname = gd.attr("gridnames").attr("__getitem__")(i);
      auto grid     = gd.attr("grids").attr("__getitem__")(i);

      const auto n =
          gridname ? gridname : py::str(var_string("dim", i).c_str());
      out["dims"].attr("append")(n);
      out["coords"][n]         = py::dict{};
      out["coords"][n]["data"] = grid;
      out["coords"][n]["dims"] = n;
    }

    out["name"] = gd.attr("dataname");
    out["data"] = gd.attr("data");

    return out;
  });

  c.def("to_xarray", [](const py::object& gd) {
    py::module_ xarray = py::module_::import_("xarray");
    return xarray.attr("DataArray").attr("from_dict")(gd.attr("to_dict")());
  });

  c.def_static("from_dict", [](const py::dict& d) {
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

    gd.data = py::cast<typename mtype::data_t>(d["data"]);

    if (d.contains("name")) {
      gd.data_name = py::cast<String>(d["name"]);
    }

    gd.grid_names = py::cast<std::array<std::string, mtype::dim>>(d["dims"]);

    py::list coords{};
    for (auto& n : gd.grid_names) {
      coords.append(d["coords"][py::str(n.c_str())]["data"]);
    }
    gd.grids = py::cast<typename mtype::grids_t>(coords);

    return gd;
  });

  c.def_static("from_xarray", [](const py::object& xarray) {
    return py::type<mtype>().attr("from_dict")(xarray.attr("to_dict")());
  });

  common_ndarray(c);

  c.def("__init__", [](mtype* v, const py::dict& d) {
    auto m = py::type<mtype>().attr("from_dict")(d);
    new (v) mtype(py::cast<mtype>(m));
  });

  py::implicitly_convertible<py::dict, mtype>();
}
}  // namespace Python
