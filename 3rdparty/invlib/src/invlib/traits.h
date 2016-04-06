/**
 * @file traits.h
 * @author Simon Pfreundschuh
 * @date 2016-03-07
 * @brief Template aliases for type traits.
 *
 */

#include <type_traits>

namespace invlib
{

template<typename T1>
using decay = typename std::decay<T1>::type;

template<typename B1>
using enable_if = typename std::enable_if<B1::value>::type;

template<typename B1, typename B2>
using enable_if_either = typename std::enable_if<B1::value || B2::value>::type;

template<typename B1>
using disable_if = typename std::enable_if<!B1::value>::type;

template<typename T1, typename T2>
using is_same = typename std::is_same<T1, T2>;

template<typename T1, typename T2>
using is_base = typename std::is_base_of<T1, T2>;

template<typename T1, typename T2>
using is_constructible = typename std::is_constructible<T1, T2>;

template<typename T1, typename T2>
using is_assignable = typename std::is_assignable<T1, T2>;

template<typename T1>
using return_type = typename std::result_of<T1>::type;

}
