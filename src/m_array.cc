/**
 * @file   m_array.cc
 * @author Richard Larsson
 * @date   2020-04-24
 * 
 * @brief  Stuff related to manipulating arrays
 */


#include "artstime.h"
#include "matpack_data.h"
#include "matpack_arrays.h"
#include "sorting.h"


template <class T>
Array<T> TimeSortTemplate(const Array<T>& arr, const ArrayOfTime& time_stamps)
{
  // Size of problem
  const Size n=time_stamps.size();
  if (arr.size() not_eq n)
    throw std::runtime_error("Cannot sort, time array does not agree with sorting array size");
  
  // Sorted index
  ArrayOfIndex sortings(n);
  get_sorted_indexes(sortings, time_stamps);
  
  // Fill the data into a new array
  Array<T> out(n);
  for (Size i=0; i<n; i++)
    out[i] = arr[sortings[i]];
  
  return out;
}

#define TIME_SORT_MACRO(VAR) \
void time_stampsSort(VAR & out, const ArrayOfTime& time_stamps, const VAR & in) \
{out = TimeSortTemplate(in, time_stamps);}

TIME_SORT_MACRO(ArrayOfTime)
TIME_SORT_MACRO(ArrayOfVector)

#undef TIME_SORT_MACRO

template <class T>
Array<T> FlattenArrayTemplate(const Array<Array<T>>& in)
{
  // Size of problem
  Index n=0;
  for (auto& array: in)
    n += array.size();
  
  // Allocate output
  Array<T> out(n);
  
  // Assignment
  Index i=0;
  for (auto& array: in) {
    for (auto& val: array) {
      out[i] = val;
      i++;
    }
  }
  
  return out;
}

#define FLATTEN_MACRO(VAR) \
void Flatten(VAR & out, const Array< VAR > & in) \
{out = FlattenArrayTemplate(in);}

FLATTEN_MACRO(ArrayOfTime)
FLATTEN_MACRO(ArrayOfVector)

#undef FLATTEN_MACRO

void Flatten(Matrix& out, const ArrayOfVector& in)
{
  if (in.size() == 0) {
    out = Matrix(0, 0);
  } else {
    const Index n = in.size();
    const Index m = in[0].nelem();
    
    if (not std::all_of(in.cbegin(), in.cend(), [m](auto& v){return m == v.nelem();}))
      throw std::runtime_error("Can only flatten array of same length data");
    
    out = Matrix(n, m);
    for (Index i=0; i<n; i++)
      out(i, joker) = in[i];
  }
}

void Flatten(Tensor3& out, const ArrayOfMatrix& in)
{
  if (in.size() == 0) {
    out = Tensor3(0, 0, 0);
  } else {
    const Index n = in.size();
    const Index c = in[0].ncols();
    const Index r = in[0].nrows();
    
    if (not std::all_of(in.cbegin(), in.cend(), [c](auto& v){return c == v.ncols();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [r](auto& v){return r == v.nrows();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    }
      
    out = Tensor3(n, r, c);
    for (Index i=0; i<n; i++)
      out(i, joker, joker) = in[i];
  }
}

void Flatten(Tensor4& out, const ArrayOfTensor3& in)
{
  if (in.size() == 0) {
    out = Tensor4(0, 0, 0, 0);
  } else {
    const Index n = in.size();
    const Index c = in[0].ncols();
    const Index r = in[0].nrows();
    const Index p = in[0].npages();
    
    if (not std::all_of(in.cbegin(), in.cend(), [c](auto& v){return c == v.ncols();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [r](auto& v){return r == v.nrows();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [p](auto& v){return p == v.npages();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    }
        
    out = Tensor4(n, p, r, c);
    for (Index i=0; i<n; i++)
      out(i, joker, joker, joker) = in[i];
  }
}

void Flatten(Tensor5& out, const ArrayOfTensor4& in)
{
  if (in.size() == 0) {
    out = Tensor5(0, 0, 0, 0, 0);
  } else {
    const Index n = in.size();
    const Index c = in[0].ncols();
    const Index r = in[0].nrows();
    const Index p = in[0].npages();
    const Index b = in[0].nbooks();
    
    if (not std::all_of(in.cbegin(), in.cend(), [c](auto& v){return c == v.ncols();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [r](auto& v){return r == v.nrows();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [p](auto& v){return p == v.npages();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [b](auto& v){return b == v.nbooks();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    }
    
    out = Tensor5(n, b, p, r, c);
    for (Index i=0; i<n; i++)
      out(i, joker, joker, joker, joker) = in[i];
  }
}

void Flatten(Tensor6& out, const ArrayOfTensor5& in)
{
  if (in.size() == 0) {
    out = Tensor6(0, 0, 0, 0, 0, 0);
  } else {
    const Index n = in.size();
    const Index c = in[0].ncols();
    const Index r = in[0].nrows();
    const Index p = in[0].npages();
    const Index b = in[0].nbooks();
    const Index s = in[0].nshelves();
    
    if (not std::all_of(in.cbegin(), in.cend(), [c](auto& v){return c == v.ncols();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [r](auto& v){return r == v.nrows();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [p](auto& v){return p == v.npages();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [b](auto& v){return b == v.nbooks();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [s](auto& v){return s == v.nshelves();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    }
    
    out = Tensor6(n, s, b, p, r, c);
    for (Index i=0; i<n; i++)
      out(i, joker, joker, joker, joker, joker) = in[i];
  }
}

void Flatten(Tensor7& out, const ArrayOfTensor6& in)
{
  if (in.size() == 0) {
    out = Tensor7(0, 0, 0, 0, 0, 0, 0);
  } else {
    const Index n = in.size();
    const Index c = in[0].ncols();
    const Index r = in[0].nrows();
    const Index p = in[0].npages();
    const Index b = in[0].nbooks();
    const Index s = in[0].nshelves();
    const Index w = in[0].nvitrines();
    
    if (not std::all_of(in.cbegin(), in.cend(), [c](auto& v){return c == v.ncols();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [r](auto& v){return r == v.nrows();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [p](auto& v){return p == v.npages();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [b](auto& v){return b == v.nbooks();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [s](auto& v){return s == v.nshelves();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    } else if (not std::all_of(in.cbegin(), in.cend(), [w](auto& v){return w == v.nvitrines();})) {
      throw std::runtime_error("Can only flatten array of same size data");
    }
    
    out = Tensor7(n, w, s, b, p, r, c);
    for (Index i=0; i<n; i++)
      out(i, joker, joker, joker, joker, joker, joker) = in[i];
  }
}
