/**
 * @file traits.h
 * @author Simon Pfreundschuh
 * @date 2016-03-07
 * @brief Template aliases for type traits.
 *
 */

#include <type_traits>
#include <functional>

#ifndef TRAITS_H
#define TRAITS_H

namespace invlib
{

// ----------------------------- //
// Remove std::reference_wrapper //
// ----------------------------- //

template<typename T1>
    T1 & remove_reference_wrapper(T1 &t)
{
    return t;
}

template<typename T1>
T1 & remove_reference_wrapper(std::reference_wrapper<T1> t)
{
    return t.get();
}

// ------------------- //
//  Custom Decay Type  //
// ------------------- //

/*
 *  Extension of std::decay that also removes std::reference_wrapper
 *  around a type.
 */
template<typename T1>
class DecayType
{
public:
    using type = typename std::decay<T1>::type;
};

template<typename T1>
class DecayType<std::reference_wrapper<T1>>
{
public:
    using type = typename std::decay<T1>::type;
};

template<typename T1>
struct DecayType<std::reference_wrapper<T1> &>
{
public:
    using type = typename std::decay<T1>::type;
};

template<typename T1>
struct DecayType<const std::reference_wrapper<T1> &>
{
public:
    using type = typename std::decay<T1>::type;
};

template<typename T1>
using decay = typename DecayType<T1>::type;

// ------------------------------ //
//  Remove std::reference_wrapper //
// ------------------------------ //

template<typename T1>
struct RemoveReferenceWrapperType
{
public:
    using type = T1;
};

template<typename T1>
struct RemoveReferenceWrapperType<std::reference_wrapper<T1>>
{
public:
    using type = T1;
};

template<typename T1>
struct RemoveReferenceWrapperType<const std::reference_wrapper<T1>&>
{
public:
    using type = T1;
};

template<typename T1>
struct RemoveReferenceWrapperType<std::reference_wrapper<T1>&>
{
public:
    using type = T1;
};

template<typename T1>
using RemoveReferenceWrapper = typename RemoveReferenceWrapperType<T1>::type;

// --------------- //
//  Custom Traits  //
// --------------- //

template<typename B1>
using enable_if = typename std::enable_if<B1::value>::type;

template<typename B1, typename B2>
using enable_if_either = typename std::enable_if<B1::value || B2::value>::type;

template<typename B1>
using disable_if = typename std::enable_if<!B1::value>::type;

template<typename T1, typename T2>
using is_same = typename std::is_same<T1, T2>;

template<typename T1>
struct Not
{
    static constexpr bool value = !T1::value;
};

template<typename T1, typename T2>
using is_base = typename std::is_base_of<T1, T2>;

template<typename T1, typename T2>
using is_constructible = typename std::is_constructible<T1, T2>;

template<typename T1>
using is_default_constructible = typename std::is_default_constructible<T1>;

template<typename T1, typename T2>
using is_copy_constructible = typename std::is_copy_constructible<T1>;

template<typename T1, typename T2>
using is_assignable = typename std::is_assignable<T1, T2>;

template<typename T1>
using is_copy_assignable = typename std::is_copy_assignable<T1>;

template<typename T1>
using is_move_assignable = typename std::is_move_assignable<T1>;

template<typename T1>
using return_type = typename std::result_of<T1>::type;

template<typename T1>
using CopyWrapper = typename std::conditional<std::is_lvalue_reference<T1>::value,
                                              std::reference_wrapper<const decay<T1>>,
                                              decay<T1>>::type;

// ----------------------------- //
//      Template equality        //
// ----------------------------- //

template
<
template <typename> class TT1,
template <typename> class TT2
>
struct is_same_template
{
    static constexpr bool value = false;
};

template
<
template <typename> class TT1
>
struct is_same_template<TT1, TT1>
{
    static constexpr bool value = true;
};

}     // namespace::invlib

#endif // TRAITS_H
