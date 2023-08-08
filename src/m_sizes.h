#pragma once

#include <workspace.h>

#include <exception>

#define NGET_GENERIC(what)                                                 \
  template <WorkspaceGroup T>                                              \
  void what##Get(Index& what, const T& x) {                                \
    if constexpr (matpack::has_##what<T>)                                  \
      what = x.what();                                                     \
    else                                                                   \
      throw std::runtime_error(var_string(                                 \
          "No " #what " for ", std::quoted(WorkspaceGroupInfo<T>::name))); \
  }

NGET_GENERIC(nelem)
NGET_GENERIC(ncols)
NGET_GENERIC(nrows)
NGET_GENERIC(npages)
NGET_GENERIC(nbooks)
NGET_GENERIC(nshelves)
NGET_GENERIC(nvitrines)
NGET_GENERIC(nlibraries)

#undef NGET_GENERIC

template <WorkspaceGroup T>
void IndexSetToLast(Index& last, const T& x) try {
  Index n;
  nelemGet(n, x);
  last = n - 1;
} catch (...) {
  throw std::runtime_error(var_string(
      "No clear last argument for ", std::quoted(WorkspaceGroupInfo<T>::name)));
}
