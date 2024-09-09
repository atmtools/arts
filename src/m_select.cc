#include <workspace.h>

template <WorkspaceGroup T>
void ArraySelect(T& needles, const T& haystack, const ArrayOfIndex& needleind) {
  needles = T(needleind.size());
  for (auto& i : needleind) {
    ARTS_USER_ERROR_IF(i < 0 || static_cast<Size>(i) >= haystack.size(),
                       "Index ",
                       i,
                       " out of range in Select.")
    needles[i] = haystack[i];
  }
}

template <WorkspaceGroup T>
void MatpackSelect(T& needles,
                   const T& haystack,
                   const ArrayOfIndex& needleind) {
  auto sh = haystack.shape();
  sh[0]   = needleind.size();
  needles = T(sh);

  for (auto& i : needleind) {
    ARTS_USER_ERROR_IF(i < 0 || i >= haystack.shape()[0],
                       "Index {} out of range in Select.",
                       i)
    needles[i] = haystack[i];
  }
}

#define array_select(arr)                                    \
  void Select(arr& x, const arr& y, const ArrayOfIndex& z) { \
    ArraySelect(x, y, z);                                    \
  }

#define matpack_select(arr)                                  \
  void Select(arr& x, const arr& y, const ArrayOfIndex& z) { \
    MatpackSelect(x, y, z);                                  \
  }

array_select(ArrayOfIndex) array_select(ArrayOfQuantumIdentifier) array_select(
    ArrayOfAbsorptionLines) array_select(ArrayOfTime) array_select(ArrayOfAgenda)
    array_select(ArrayOfGriddedField1) array_select(ArrayOfGriddedField2) array_select(
        ArrayOfGriddedField3) array_select(ArrayOfGriddedField4) array_select(ArrayOfVector)
        array_select(ArrayOfMatrix) array_select(ArrayOfTensor3) array_select(
            ArrayOfTensor4) array_select(ArrayOfTensor5) array_select(ArrayOfTensor6)
            array_select(ArrayOfTensor7) array_select(ArrayOfSparse) array_select(
                ArrayOfScatteringMetaData) array_select(ArrayOfSingleScatteringData)
                array_select(ArrayOfTelsemAtlas) array_select(ArrayOfString) array_select(
                    ArrayOfSun) array_select(ArrayOfCIARecord) array_select(ArrayOfAtmPoint)
                    array_select(ArrayOfMuelmatVector) array_select(ArrayOfMuelmatMatrix) array_select(
                        ArrayOfPropmatVector) array_select(ArrayOfPropmatMatrix)
                        array_select(ArrayOfStokvecVector) array_select(ArrayOfStokvecMatrix) array_select(
                            ArrayOfXsecRecord) array_select(ArrayOfSpeciesTag)

                            array_select(ArrayOfArrayOfGriddedField1) array_select(
                                ArrayOfArrayOfGriddedField2) array_select(ArrayOfArrayOfGriddedField3)
                                array_select(ArrayOfArrayOfString) array_select(
                                    ArrayOfArrayOfScatteringMetaData) array_select(ArrayOfArrayOfSingleScatteringData)
                                    array_select(ArrayOfArrayOfIndex) array_select(
                                        ArrayOfArrayOfAbsorptionLines) array_select(ArrayOfArrayOfTime)
                                        array_select(ArrayOfArrayOfMuelmatVector) array_select(
                                            ArrayOfArrayOfMuelmatMatrix) array_select(ArrayOfArrayOfPropmatVector)
                                            array_select(ArrayOfArrayOfPropmatMatrix) array_select(
                                                ArrayOfArrayOfStokvecVector) array_select(ArrayOfArrayOfStokvecMatrix)
                                                array_select(ArrayOfArrayOfVector) array_select(
                                                    ArrayOfArrayOfMatrix) array_select(ArrayOfArrayOfTensor3)
                                                    array_select(ArrayOfArrayOfTensor6) array_select(
                                                        ArrayOfArrayOfSpeciesTag)

                                                        matpack_select(Vector) matpack_select(
                                                            Matrix) matpack_select(Tensor3)
                                                            matpack_select(Tensor4) matpack_select(
                                                                Tensor5) matpack_select(Tensor6)
                                                                matpack_select(
                                                                    Tensor7)

                                                                    void Select(  // WS Generic Output:
                                                                        Sparse&
                                                                            needles,
                                                                        // WS Generic Input:
                                                                        const Sparse&
                                                                            haystack,
                                                                        const ArrayOfIndex&
                                                                            needleind) {
  // We construct the output in this dummy variable, so that the
  // method also works properly if needles and haystack are the same
  // variable.
  Sparse dummy(needleind.size(), haystack.ncols());

  // If needleind only contains -1 as the only element, copy the whole thing
  if (needleind.size() == 1 && needleind[0] == -1) {
    needles = haystack;
    return;
  }

  for (Size i = 0; i < needleind.size(); i++) {
    if (haystack.nrows() <= needleind[i]) {
      std::ostringstream os;
      os << "The input matrix only has " << haystack.nrows()
         << " rows. But one of the needle indexes is " << needleind[i] << "."
         << std::endl;
      os << "The indexes must be between 0 and " << haystack.nrows() - 1;
      throw std::runtime_error(os.str());
    }

    if (needleind[i] < 0) {
      std::ostringstream os;
      os << "One of the needle indexes is " << needleind[i] << "." << std::endl;
      os << "The indexes must be between 0 and " << haystack.nrows() - 1;
      throw std::runtime_error(os.str());
    }

    // Copy this row of the sparse matrix.
    // This code is inefficient for Sparse, but I leave it like
    // this to be consistent with the other data types for which
    // Select is implemented.
    for (Index j = 0; j < haystack.ncols(); ++j) {
      Numeric value = haystack(needleind[i], j);
      if (0 != value) dummy.rw(i, j) = value;
    }
  }

  needles = dummy;
}

void Select(AscendingGrid& needles,
            const AscendingGrid& haystack,
            const ArrayOfIndex& needleind) {
  Vector n;
  Select(n, haystack, needleind);
  needles = std::move(n);
}
