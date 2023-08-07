#include <workspace.h>

template <typename AoT, WorkspaceGroup T>
void ArrayExtracter(T& e, const AoT& arr, const Index& index) {
  ARTS_USER_ERROR_IF(index >= arr.nelem() or index < 0,
                     "The index ",
                     index,
                     " is outside the range of the array.")

  e = arr[index];
}

template <typename AoT, WorkspaceGroup T>
void MatpackExtracter(T& e, const AoT& arr, const Index& index) {
  ARTS_USER_ERROR_IF(index >= arr.shape()[0] or index < 0,
                     "The index ",
                     index,
                     " is outside the range of the matpack type.")

  e = arr[index];
}

#define extracts(Arr, Val)                                   \
  void Extract(Val& e, const Arr& arr, const Index& index) { \
    ArrayExtracter(e, arr, index);                           \
  }

#define array_extracts(Val) extracts(ArrayOf##Val, Val)

array_extracts(Index)
array_extracts(Numeric)
array_extracts(QuantumIdentifier)
array_extracts(AbsorptionLines)
array_extracts(Time)
array_extracts(Agenda)
array_extracts(GriddedField1)
array_extracts(GriddedField2)
array_extracts(GriddedField3)
array_extracts(GriddedField4)
array_extracts(Vector)
array_extracts(Matrix)
array_extracts(Tensor3)
array_extracts(Tensor4)
array_extracts(Tensor5)
array_extracts(Tensor6)
array_extracts(Tensor7)
array_extracts(ScatteringMetaData)
array_extracts(JacobianTarget)
array_extracts(Ppath)
array_extracts(Sparse)
array_extracts(SingleScatteringData)
array_extracts(TelsemAtlas)
array_extracts(String)

array_extracts(ArrayOfSpeciesTag)
array_extracts(ArrayOfAbsorptionLines)
array_extracts(ArrayOfGriddedField1)
array_extracts(ArrayOfGriddedField2)
array_extracts(ArrayOfGriddedField3)
array_extracts(ArrayOfScatteringMetaData)
array_extracts(ArrayOfSingleScatteringData)
array_extracts(ArrayOfString)
array_extracts(ArrayOfIndex)

#define matpack_extract(Arr, Val) \
void Extract(Val& e, const Arr& arr, const Index& index) {\
  MatpackExtracter(e, arr, index);\
}

matpack_extract(Vector, Numeric)
matpack_extract(Matrix, Vector)
matpack_extract(Tensor3, Matrix)
matpack_extract(Tensor4, Tensor3)
matpack_extract(Tensor5, Tensor4)
matpack_extract(Tensor6, Tensor5)
matpack_extract(Tensor7, Tensor6)