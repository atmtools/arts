/*
    pybind11/eigen_tensor.h: Conversion of multi-dimensional numpy arrays to Eigen
    tensors.

    Copyright (c) 2020 Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/
#pragma once

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include <iostream>

#if defined(__INTEL_COMPILER)
// implicit conversion of a 64-bit integral type to a smaller integral type
// (potential portability problem)
#  pragma warning(disable: 1682)
#elif defined(__GNUG__) || defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#  ifdef __clang__
//   Eigen generates a bunch of implicit-copy-constructor-is-deprecated warnings with -Wdeprecated
//   under Clang, so disable that warning here:
#    pragma GCC diagnostic ignored "-Wdeprecated"
#  endif
#  if __GNUC__ >= 7
#    pragma GCC diagnostic ignored "-Wint-in-bool-context"
#  endif
#endif

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable: 4127) // warning C4127: Conditional expression is constant
#  pragma warning(disable: 4996) // warning C4996: std::unary_negate is deprecated in C++17
#endif

#include <Eigen/Core>
#include <Eigen/CXX11/Tensor>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

////////////////////////////////////////////////////////////////////////////////
// Helper types and aliases
////////////////////////////////////////////////////////////////////////////////

PYBIND11_NAMESPACE_BEGIN(detail)

//
// is_tensor
//

struct is_tensor_impl {
    template <typename Scalar, int NumIndices, int Options, typename IndexType>
        static std::true_type check(Eigen::Tensor<Scalar, NumIndices, Options, IndexType> *);
    static std::false_type check(...);
};

template <typename T>
using is_eigen_tensor = decltype(is_tensor_impl::check((intrinsic_t<T>*)nullptr));

//
// is_tensor_map
//

struct is_tensor_map_impl {
    template <typename PlainObjectType, int Options, template <typename> typename MakePointer>
        static std::true_type check(Eigen::TensorMap<PlainObjectType, Options, MakePointer> *);
    static std::false_type check(...);
};

template <typename T>
using is_eigen_tensor_map = decltype(is_tensor_map_impl::check((intrinsic_t<T>*)nullptr));

//
// is_tensor_ref
//

struct is_tensor_ref_impl {
    template <typename PlainObjectType>
        static std::true_type check(Eigen::TensorRef<PlainObjectType> *);
    static std::false_type check(...);
};

template <typename T>
using is_eigen_tensor_ref = decltype(is_tensor_ref_impl::check((intrinsic_t<T>*)nullptr));

////////////////////////////////////////////////////////////////////////////////
// Tensor conformable
////////////////////////////////////////////////////////////////////////////////

/* Captures numpy/eigen conformability status */

template <bool EigenRowMajor> struct EigenTensorConformable {
    /** Can numpy type converted to Eigen type? */
    bool conformable = false;
    /** Size of the tensor along each dimension. */
    std::vector<size_t> dimensions;
    /** The rank of the tensor. */
    size_t rank = 0;

    EigenTensorConformable(bool fits = false) : conformable{fits} {}
    EigenTensorConformable(std::vector<size_t> dimensions_)
        : conformable{true}, dimensions(dimensions_) {}

    operator bool() const { return conformable; }

};

////////////////////////////////////////////////////////////////////////////////
// Eigen Tensor traits
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t i = 0>
    struct tensor_dimensions {
        static constexpr char c[] = {'m' + i, '\0'};
        static constexpr auto text = _(c) + _(", ") + tensor_dimensions<N-1, i + 1>::text;
    };

template <size_t i>
struct tensor_dimensions<1, i> {
    static constexpr char c[] = {'m' + i, '\0'};
    static constexpr auto text = _(c);
};

template <size_t i>
struct tensor_dimensions<0, i> {
    static constexpr auto text = _("");
};

template<typename Tensor>
struct get_options {
    static constexpr int value = Tensor::Options;
};

template<typename Tensor>
struct get_options<Eigen::TensorRef<Tensor>> {
    static constexpr int value = Tensor::Options;
};

/** Properties of Eigen Tensor types.
 *
 * Helper struct to extract and combine information for Eigen tensor types.
 */
template <typename Type_> struct EigenTensorProps {
    /** The Eigen tensor type. */
    using Type = Type_;
    /** The Scalar used for its coefficients */
    using Scalar = typename Type::Scalar;
    /** Integer type used for indices */
    using Index = typename Type::Index;
    /** Corresponding type that manages its own data. */
    using PlainObjectType = Eigen::Tensor<Scalar, Type::NumIndices, get_options<Type>::value, typename Type::Index>;
    /** Corresponding map type. */
    using EigenMapType = Eigen::TensorMap<PlainObjectType, get_options<Type>::value>;

    /** The rank of the tensor */
    static constexpr Eigen::Index rank = Type::NumIndices;
    /** Is data in row-major storage format? */
    static constexpr bool row_major = static_cast<int>(Type::Layout) == Eigen::RowMajor;
    /** Is the tensor writable? */
    static constexpr bool writable = !std::is_const<typename Type::CoeffReturnType>::value;

    /** Are we trying to avoid copying ? */
    static constexpr bool avoid_copy = is_eigen_tensor_map<Type>() || is_eigen_tensor_ref<Type>();
    /** Do we require numpy array to be writable?
     *
     * This is the case when we are trying to avoid copying but still want mutable
     * access.
     */
    static constexpr bool needs_writable = writable && avoid_copy;

    using NumpyArrayType = array_t<Scalar, array::forcecast | (row_major ? array::c_style : array::f_style)>;


    // Takes an input array and determines whether we can make it fit into the Eigen type.  If
    // the array is a vector, we attempt to fit it into either an Eigen 1xN or Nx1 vector
    // (preferring the latter if it will fit in either, i.e. for a fully dynamic matrix type).
    static bool conformable(const array &arr) {
      if (!arr) { 
          return false;
      }

      const auto dims = arr.ndim();

      // If we are allowed to copy, all we care about is the tensor rank.
      if (dims != rank) {
          return false;
      }
      return true;
    }

    // Takes an input array and determines whether we can make it fit into the Eigen type.  If
    // the array is a vector, we attempt to fit it into either an Eigen 1xN or Nx1 vector
    // (preferring the latter if it will fit in either, i.e. for a fully dynamic matrix type).
    static bool needs_copy(const array &arr) {
      if (!isinstance<NumpyArrayType>(arr)) {
        return true;
      }

      std::array<long int, rank> expected_strides;
      size_t cs = sizeof(Scalar);
      expected_strides[0] = cs;
      for (long int i = 0; i < rank - 1; ++i) {
        cs *= arr.shape(row_major ? rank - i - 1 : i);
        expected_strides[i + 1] = cs;
      }

      for (long int i = 0; i < rank; ++i) {
          if (arr.strides(i) != expected_strides[row_major ? rank - i - 1 : i]) {
          return true;
        }
      }
      return false;
    }

    /** Descriptor for Python-side output. */
    static constexpr auto descriptor =
        _("numpy.ndarray[") + npy_format_descriptor<Scalar>::name +
        _("[")  + tensor_dimensions<rank>::text +
        _("]") +
        _<needs_writable>(", flags.writable", "") +
        _<avoid_copy && row_major>(", flags.c_contiguous", "") +
        _<avoid_copy && !row_major>(", flags.f_contiguous", "") +
        _("]");
};

template<size_t rank>
std::array<Eigen::Index, rank> get_dimensions(array a) {
    std::array<Eigen::Index, rank> result{};
    for (size_t i = 0; i < rank; ++i) {
        result[i] = a.shape(i);
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Converting Eigen to numpy
////////////////////////////////////////////////////////////////////////////////

/** Cast Eigen Tensor to numpy array.
 *
 * Returns numpy array object referencing the data of the given Eigen tensor.
 *
 * @param src The Eigen tensor
 * @param base Python object to tie the lifetime of the newly created array.
 * @param Value of the writable flag of the newly created array.
 * @return Numpy array object handle that references the data of the
 *         given Eigen tensor.
 */
template <typename props>
handle tensor_array_cast(typename props::Type &src,
                         handle base = handle(),
                         bool writable = props::writable) {
  typename props::NumpyArrayType a{src.dimensions(), src.data(), base};

  if (!writable) {
    array_proxy(a.ptr())->flags &= !detail::npy_api::NPY_ARRAY_WRITEABLE_;
  }

  return a.release();
}

// Takes a pointer to some dense, plain Eigen type, builds a capsule around it, then returns a numpy
// array that references the encapsulated data with a python-side reference to the capsule to tie
// its destruction to that of any dependent python objects.  Const-ness is determined by whether or
// not the Type of the pointer given is const.

/** Encapsule Eigen tensor.
 *
 * Creates a capsule around given Eigen  to manage its life time and creates
 * a numpy array pointing to the tensor's data. Lifetime of the capsule is tied
 * to the newly created numpy array.
 *
 * @param src Pointer to Eigen tensor to return to Python.
 * @return A numpy array pointing to the data of the tensor, which is tied to the
 *         capsule managing the tensors lifetime.
 */
template <typename props,
          typename Type>
handle tensor_encapsulate(Type *src) {
  capsule base(src, [](void *o) { delete static_cast<Type *>(o); });
  return tensor_array_cast<props>(*src, base);
}

////////////////////////////////////////////////////////////////////////////////
// Converting Eigen to numpy
////////////////////////////////////////////////////////////////////////////////

template<typename Type>
struct tensor_type_caster {
    using Scalar = typename Type::Scalar;
    using props = EigenTensorProps<Type>;
    using MapType = typename props::EigenMapType;

    bool load(handle src, bool convert) {
        // If we're in no-convert mode, only load if given an array of the correct type
        if (!convert && !isinstance<array_t<Scalar>>(src))
            return false;

        auto buf = array::ensure(src);
        if (!buf) return false;

        auto fits = props::conformable(buf);
        if (!fits) {
            return false;
        }

        // Allocate the new type, then build a numpy reference into it
        value = Type{};
        auto dims = get_dimensions<props::rank>(buf);
        value.resize(dims);


        auto ref = reinterpret_steal<typename props::NumpyArrayType>(tensor_array_cast<props>(value, none()));

        int result = detail::npy_api::get().PyArray_CopyInto_(ref.ptr(), buf.ptr());
        if (result < 0) { // Copy failed!
            PyErr_Clear();
            return false;
        }
        return true;
    }

private:

    // Cast implementation
    template <typename CType>
    static handle cast_impl(CType *src, return_value_policy policy, handle parent) {
        switch (policy) {
            case return_value_policy::take_ownership:
            case return_value_policy::automatic:
                return tensor_encapsulate<props>(src);
            case return_value_policy::move:
                return tensor_encapsulate<props>(new CType(std::move(*src)));
            case return_value_policy::copy:
                return tensor_array_cast<props>(*src);
            case return_value_policy::reference:
            case return_value_policy::automatic_reference:
                return tensor_array_cast<props>(*src, none());
            case return_value_policy::reference_internal:
                return tensor_array_cast<props>(*src, parent);
            default:
                throw cast_error("unhandled return_value_policy: should not happen!");
        };
    }

public:

    // Normal returned non-reference, non-const value:
    static handle cast(Type &&src, return_value_policy /* policy */, handle parent) {
        return cast_impl(&src, return_value_policy::move, parent);
    }
    // If you return a non-reference const, we mark the numpy array readonly:
    static handle cast(const Type &&src, return_value_policy /* policy */, handle parent) {
        return cast_impl(&src, return_value_policy::move, parent);
    }
    // lvalue reference return; default (automatic) becomes copy
    static handle cast(Type &src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic || policy == return_value_policy::automatic_reference)
            policy = return_value_policy::copy;
        return cast_impl(&src, policy, parent);
    }
    // const lvalue reference return; default (automatic) becomes copy
    static handle cast(const Type &src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic || policy == return_value_policy::automatic_reference)
            policy = return_value_policy::copy;
        return cast(&src, policy, parent);
    }
    // non-const pointer return
    static handle cast(Type *src, return_value_policy policy, handle parent) {
        return cast_impl(src, policy, parent);
    }
    // const pointer return
    static handle cast(const Type *src, return_value_policy policy, handle parent) {
        return cast_impl(src, policy, parent);
    }

    static constexpr auto name = props::descriptor;

    operator Type*() { return &value; }
    operator Type&() { return value; }
    operator Type&&() && { return std::move(value); }
    template <typename T> using cast_op_type = movable_cast_op_type<T>;

private:
    Type value;
};



template <typename Type>
struct type_caster<Type, enable_if_t<is_eigen_tensor<Type>::value>>
    : public tensor_type_caster<Type> {};

template <typename Type>
struct type_caster<Type, enable_if_t<is_eigen_tensor_map<Type>::value || is_eigen_tensor_ref<Type>::value>>
    : public tensor_type_caster<Type> {

    using Base = tensor_type_caster<Type>;
    using props = typename Base::props;
    using MapType = typename props::EigenMapType;
    using Array = typename props::NumpyArrayType;
    using Tensor = typename props::PlainObjectType;


public:

    bool load(handle src, bool convert) {

        //No need to carry on if array isn't even conform.
        auto buf = array::ensure(src);
        if (!buf) return false;
        auto fits = props::conformable(buf);
        if (!fits) {
           return false;
        }

        bool need_copy = props::needs_copy(buf);
        if (need_copy) {
          if (!convert || props::avoid_copy) {
              return false;
          }
          copy_or_ref = std::move(Array::ensure(src));
          if (!copy_or_ref) {
              return false;
          }
          loader_life_support::add_patient(copy_or_ref);
        } else {
          Array aref = reinterpret_borrow<Array>(src);
          if (!aref) {
              return false;
          }
          copy_or_ref = aref;
        }

        map.reset(new MapType(copy_or_ref.mutable_data(), get_dimensions<props::rank>(copy_or_ref)));
        value.reset(new Type(*map));
        return true;
    }

    operator Type*() { return value.get(); }
    operator Type&() { return *value; }
    template <typename _T> using cast_op_type = pybind11::detail::cast_op_type<_T>;

private:

    std::unique_ptr<MapType> map = nullptr;
    std::unique_ptr<Type> value = nullptr;
    Array copy_or_ref;

};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)

#if defined(__GNUG__) || defined(__clang__)
#  pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#  pragma warning(pop)
#endif
