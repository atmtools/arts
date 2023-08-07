#include <workspace.h>

template <WorkspaceGroup T>
void ArraySelect(T& needles, const T& haystack, const ArrayOfIndex& needleind) {
  needles = T(needleind.nelem());
  for (auto& i : needleind) {
    ARTS_USER_ERROR_IF(
        i < 0 || i >= haystack.nelem(), "Index ", i, " out of range in Select.")
    needles[i] = haystack[i];
  }
}

template <WorkspaceGroup T>
void MatpackSelect(T& needles,
                   const T& haystack,
                   const ArrayOfIndex& needleind) {
  auto sh = haystack.shape();
  sh[0] = needleind.nelem();
  needles = T(sh);

  for (auto& i : needleind) {
    ARTS_USER_ERROR_IF(
        i < 0 || i >= haystack.shape()[0], "Index ", i, " out of range in Select.")
    needles[i] = haystack[i];
  }
}

#define array_select(arr) void Select(arr& x, const arr& y, const ArrayOfIndex& z) { \
  ArraySelect(x, y, z); \
}

#define matpack_select(arr) void Select(arr& x, const arr& y, const ArrayOfIndex& z) { \
  MatpackSelect(x, y, z); \
}

array_select(ArrayOfIndex)
array_select(ArrayOfQuantumIdentifier)
array_select(ArrayOfAbsorptionLines)
array_select(ArrayOfTime)
array_select(ArrayOfAgenda)
array_select(ArrayOfGriddedField1)
array_select(ArrayOfGriddedField2)
array_select(ArrayOfGriddedField3)
array_select(ArrayOfGriddedField4)
array_select(ArrayOfVector)
array_select(ArrayOfMatrix)
array_select(ArrayOfTensor3)
array_select(ArrayOfTensor4)
array_select(ArrayOfTensor5)
array_select(ArrayOfTensor6)
array_select(ArrayOfTensor7)
array_select(ArrayOfScatteringMetaData)
array_select(ArrayOfJacobianTarget)
array_select(ArrayOfPpath)
array_select(ArrayOfSparse)
array_select(ArrayOfSingleScatteringData)
array_select(ArrayOfTelsemAtlas)
array_select(ArrayOfString)

matpack_select(Vector)
matpack_select(Matrix)
matpack_select(Tensor3)
matpack_select(Tensor4)
matpack_select(Tensor5)
matpack_select(Tensor6)
matpack_select(Tensor7)
