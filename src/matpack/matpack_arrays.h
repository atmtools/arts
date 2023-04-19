#pragma once

#include "array.h"
#include "matpack_constexpr.h"
#include "matpack_data.h"

//! An array of vectors of Numeric
using ArrayOfVector = Array<Vector>;

//! An array of matrices of Numeric
using ArrayOfMatrix = Array<Matrix>;

//! An array of tensors of Numeric of rank 3
using ArrayOfTensor3 = Array<Tensor3>;

//! An array of tensors of Numeric of rank 4
using ArrayOfTensor4 = Array<Tensor4>;

//! An array of tensors of Numeric of rank 5
using ArrayOfTensor5 = Array<Tensor5>;

//! An array of tensors of Numeric of rank 6
using ArrayOfTensor6 = Array<Tensor6>;

//! An array of tensors of Numeric of rank 7
using ArrayOfTensor7 = Array<Tensor7>;

//! An array of an array of vectors of Numeric
using ArrayOfArrayOfVector = Array<ArrayOfVector>;

//! An array of an array of matrices of Numeric
using ArrayOfArrayOfMatrix = Array<ArrayOfMatrix>;

//! An array of an array of tensors of Numeric of rank 3
using ArrayOfArrayOfTensor3 = Array<ArrayOfTensor3>;

//! An array of an array of tensors of Numeric of rank 4
using ArrayOfArrayOfTensor4 = Array<ArrayOfTensor4>;

//! An array of an array of tensors of Numeric of rank 5
using ArrayOfArrayOfTensor5 = Array<ArrayOfTensor5>;

//! An array of an array of tensors of Numeric of rank 6
using ArrayOfArrayOfTensor6 = Array<ArrayOfTensor6>;

//! An array of an array of tensors of Numeric of rank 7
using ArrayOfArrayOfTensor7 = Array<ArrayOfTensor7>;

//! An array of vectors of Complex
using ArrayOfComplexVector = Array<ComplexVector>;

//! An array of matrices of Complex
using ArrayOfComplexMatrix = Array<ComplexMatrix>;

//! An array of Vectors of length 2 of Numeric
using ArrayOfVector2 = Array<Vector2>;

//! An array of Vectors of length 3 of Numeric
using ArrayOfVector3 = Array<Vector3>;
