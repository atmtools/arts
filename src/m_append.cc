#include <workspace.h>

void Append(  // WS Generic Output:
    ArrayOfSpeciesTag& out,
    // WS Generic Input:
    const ArrayOfSpeciesTag& in,
    const String& direction [[maybe_unused]]) {
  const ArrayOfSpeciesTag* in_pnt;
  ArrayOfSpeciesTag in_copy;

  if (&in == &out) {
    in_copy = in;
    in_pnt = &in_copy;
  } else
    in_pnt = &in;

  const ArrayOfSpeciesTag& in_ref = *in_pnt;

  // Reserve memory in advance to avoid reallocations:
  out.reserve(out.size() + in_ref.size());
  // Append in to end of out:
  for (Size i = 0; i < in_ref.size(); ++i) out.push_back(in_ref[i]);
}

/* Implementation for array types to append single element */
void Append(  // WS Generic Output:
    ArrayOfSpeciesTag& out,
    // WS Generic Input:
    const SpeciesTag& in,
    const String& direction [[maybe_unused]]) {
  // Append in to end of out:
  out.push_back(in);
}

/* Implementation for array types to append single element */
void Append(
    // WS Generic Output:
    ArrayOfAgenda& out,
    // WS Generic Input:
    const Agenda& in,
    const String& direction [[maybe_unused]]) {
  // Append in to end of out:
  auto& newag = out.emplace_back(in);
  newag.finalize();
}

/* Implementation for array types to append single element */
void Append(
    // WS Generic Output:
    ArrayOfAgenda& out,
    // WS Generic Input:
    const ArrayOfAgenda& in,
    const String& direction [[maybe_unused]]) {
  // Append in to end of out:
  for (const auto& it : in) {
    auto& newag = out.emplace_back(it);
    newag.finalize();
  }
}

/* Implementation for Vector */
void Append(  // WS Generic Output:
    Vector& out,
    // WS Generic Input:
    const Vector& in,
    const String& direction [[maybe_unused]]) {
  const Vector* in_pnt;
  Vector in_copy;

  if (&in == &out) {
    in_copy = in;
    in_pnt = &in_copy;
  } else
    in_pnt = &in;

  const Vector& in_ref = *in_pnt;

  // Get backup of out:
  Vector dummy = out;

  // Make out the right size:
  out.resize(dummy.nelem() + in_ref.nelem());

  // Copy dummy to first part of out:
  if (dummy.nelem()) out[Range(0, dummy.nelem())] = dummy;

  // Copy in to last part of out:
  if (in_ref.nelem()) out[Range(dummy.nelem(), in_ref.nelem())] = in_ref;
}

/* Implementation for Matrix */
void Append(  // WS Generic Output:
    Matrix& out,
    // WS Generic Input:
    const Matrix& in,
    const String& direction) {
  const Matrix* in_pnt;
  Matrix in_copy;

  if (&in == &out) {
    in_copy = in;
    in_pnt = &in_copy;
  } else
    in_pnt = &in;

  const Matrix& in_ref = *in_pnt;

  // Get backup of out:
  Matrix dummy = out;

  if (!out.nrows() || !out.ncols()) {
    out = in_ref;
  } else if (direction == "leading") {
    if (out.ncols() != in_ref.ncols())
      throw std::runtime_error(
          "Input and output matrix must have the same number of columns.");

    out.resize(dummy.nrows() + in_ref.nrows(), dummy.ncols());

    if (dummy.nrows() && dummy.ncols())
      out(Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
    if (dummy.nrows() && in_ref.nrows() && in_ref.ncols())
      out(Range(dummy.nrows(), in_ref.nrows()), Range(0, in_ref.ncols())) =
          in_ref;
  } else if (direction == "trailing") {
    if (out.nrows() != in_ref.nrows())
      throw std::runtime_error(
          "Input and output matrix must have the same number of rows.");

    out.resize(dummy.nrows(), dummy.ncols() + in_ref.ncols());

    if (dummy.nrows() && dummy.ncols())
      out(Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
    if (dummy.ncols() && in_ref.nrows() && in_ref.ncols())
      out(Range(0, in_ref.nrows()), Range(dummy.ncols(), in_ref.ncols())) =
          in_ref;
  } else
    throw std::runtime_error(R"(Dimension must be either "leading" or "trailing".)");
}

/* Implementation for Matrix/Vector */
void Append(  // WS Generic Output:
    Matrix& out,
    // WS Generic Input:
    const Vector& in,
    const String& direction) {
  // Get backup of out:
  Matrix dummy = out;

  if (direction == "leading") {
    if (!out.nrows() || !out.ncols()) {
      out = ExhaustiveMatrixView{in};
    } else {
      if (out.ncols() != in.nelem())
        throw std::runtime_error(
            "Number of elements in the input Vector has to match "
            "the number of columns in the output Matrix.");

      out.resize(dummy.nrows() + 1, dummy.ncols());
      out(Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
      out(Range(dummy.nrows(), 1), Range(0, in.nelem())) =
          transpose(ExhaustiveMatrixView{in});
    }
  } else if (direction == "trailing") {
    if (!out.nrows() || !out.ncols()) {
      out = transpose(ExhaustiveMatrixView{in});
    } else if (in.nelem()) {
      if (out.nrows() != in.nelem() && out.nrows() && out.ncols())
        throw std::runtime_error(
            "Number of elements in the input Vector has to match "
            "the number of rows in the output Matrix.");

      out.resize(dummy.nrows(), dummy.ncols() + 1);
      out(Range(0, dummy.nrows()), Range(0, dummy.ncols())) = dummy;
      out(Range(0, in.nelem()), Range(dummy.ncols(), 1)) =
          ExhaustiveMatrixView{in};
    }
  } else
    throw std::runtime_error(R"(Dimension must be either "leading" or "trailing".)");
}

/* Implementation for Vector/Numeric */
void Append(  // WS Generic Output:
    Vector& out,
    // WS Generic Input:
    const Numeric& in,
    const String& direction [[maybe_unused]]) {
  // Get backup of out:
  Vector dummy = out;

  // Make out the right size:
  out.resize(dummy.nelem() + 1);

  // Copy dummy to first part of out:
  if (dummy.nelem()) out[Range(0, dummy.nelem())] = dummy;

  // Copy in to last part of out:
  out[Range(dummy.nelem(), 1)] = in;
}

/* Implementation for Tensor3/Matrix */
void Append(  // WS Generic Output:
    Tensor3& out,
    // WS Generic Input:
    const Matrix& in,
    //            const String& direction,
    const String& direction [[maybe_unused]]) {
  // Get backup of out:
  Tensor3 dummy = out;

  if (!out.npages() || !out.nrows() || !out.ncols()) {
    out.resize(1, in.nrows(), in.ncols());
    out(0, joker, joker) = in;
  } else {
    if (out.nrows() != in.nrows() || out.ncols() != in.ncols())
      throw std::runtime_error(
          "Number of rows and columns in the input Matrix have to match\n"
          "the number of rows and columns in the output Tensor3.");

    out.resize(dummy.npages() + 1, dummy.nrows(), dummy.ncols());
    out(Range(0, dummy.npages()),
        Range(0, dummy.nrows()),
        Range(0, dummy.ncols())) = dummy;
    out(dummy.npages(), Range(0, dummy.nrows()), Range(0, dummy.ncols())) = in;
  }
}

/* Implementation for Tensor3 */
void Append(  // WS Generic Output:
    Tensor3& out,
    // WS Generic Input:
    const Tensor3& in,
    //            const String& direction,
    const String& direction [[maybe_unused]]) {
  const Tensor3* in_pnt;
  Tensor3 in_copy;

  if (&in == &out) {
    in_copy = in;
    in_pnt = &in_copy;
  } else
    in_pnt = &in;

  const Tensor3& in_ref = *in_pnt;

  // Get backup of out:
  Tensor3 dummy = out;

  if (out.nrows() != in_ref.nrows() || out.ncols() != in_ref.ncols())
    throw std::runtime_error(
        "Tensor3 append is performed in pages dimension.\n"
        "All other dimensions (rows, columns) must have identical\n"
        "sizes in In and Out Tensor.");

  out.resize(dummy.npages() + in_ref.npages(), dummy.nrows(), dummy.ncols());

  if (dummy.npages() && dummy.nrows() && dummy.ncols())
    out(Range(0, dummy.npages()),
        Range(0, dummy.nrows()),
        Range(0, dummy.ncols())) = dummy;
  if (dummy.npages() && in_ref.npages() && in_ref.nrows() && in_ref.ncols())
    out(Range(dummy.npages(), in_ref.npages()),
        Range(0, in_ref.nrows()),
        Range(0, in_ref.ncols())) = in_ref;
}

/* Implementation for Tensor4/Tensor3 */
void Append(  // WS Generic Output:
    Tensor4& out,
    // WS Generic Input:
    const Tensor3& in,
    //            const String& direction,
    const String& direction [[maybe_unused]]) {
  // Get backup of out:
  Tensor4 dummy = out;

  if (!out.nbooks() || !out.npages() || !out.nrows() || !out.ncols()) {
    out.resize(1, in.npages(), in.nrows(), in.ncols());
    out(0, joker, joker, joker) = in;
  } else {
    if (out.npages() != in.npages() || out.nrows() != in.nrows() ||
        out.ncols() != in.ncols())
      throw std::runtime_error(
          "Dimensions of input Tensor3 have to match corresponding\n"
          "dimensions in the output Tensor4.");

    out.resize(
        dummy.nbooks() + 1, dummy.npages(), dummy.nrows(), dummy.ncols());
    out(Range(0, dummy.nbooks()),
        Range(0, dummy.npages()),
        Range(0, dummy.nrows()),
        Range(0, dummy.ncols())) = dummy;
    out(dummy.nbooks(),
        Range(0, dummy.npages()),
        Range(0, dummy.nrows()),
        Range(0, dummy.ncols())) = in;
  }
}

/* Implementation for Tensor4 */
void Append(  // WS Generic Output:
    Tensor4& out,
    // WS Generic Input:
    const Tensor4& in,
    //            const String& direction,
    const String& direction [[maybe_unused]]) {
  const Tensor4* in_pnt;
  Tensor4 in_copy;

  if (&in == &out) {
    in_copy = in;
    in_pnt = &in_copy;
  } else
    in_pnt = &in;

  const Tensor4& in_ref = *in_pnt;

  // Get backup of out:
  Tensor4 dummy = out;

  if (out.npages() != in_ref.npages() || out.nrows() != in_ref.nrows() ||
      out.ncols() != in_ref.ncols())
    throw std::runtime_error(
        "Tensor4 append is performed in books dimension.\n"
        "All other dimensions (pages, rows, columns) must have identical\n"
        "sizes in In and Out Tensor.");

  out.resize(dummy.nbooks() + in_ref.nbooks(),
             dummy.npages(),
             dummy.nrows(),
             dummy.ncols());

  if (dummy.nbooks() && dummy.npages() && dummy.nrows() && dummy.ncols())
    out(Range(0, dummy.nbooks()),
        Range(0, dummy.npages()),
        Range(0, dummy.nrows()),
        Range(0, dummy.ncols())) = dummy;
  if (dummy.nbooks() && in_ref.nbooks() && in_ref.npages() && in_ref.nrows() &&
      in_ref.ncols())
    out(Range(dummy.nbooks(), in_ref.nbooks()),
        Range(0, in_ref.npages()),
        Range(0, in_ref.nrows()),
        Range(0, in_ref.ncols())) = in_ref;
}

/* Implementation for String */
void Append(  // WS Generic Output:
    String& out,
    // WS Generic Input:
    const String& in,
    const String& direction [[maybe_unused]]) {
  // String stream for easy string operations:
  std::ostringstream os;

  os << out << in;

  out = os.str();
}

#define self_append(Arr)              \
  void Append(Arr& x,                 \
              const Arr& y,           \
              const String&) {        \
    for (auto& z : y) x.push_back(z); \
  }

#define once_append(Val)       \
  void Append(Array<Val>& x,   \
              const Val& y,    \
              const String&) { \
    x.push_back(y);            \
  }

once_append(Index)
once_append(QuantumIdentifier)
once_append(AbsorptionLines)
once_append(Time)
once_append(GriddedField1)
once_append(GriddedField2)
once_append(GriddedField3)
once_append(GriddedField4)
once_append(Vector)
once_append(Matrix)
once_append(Tensor3)
once_append(Tensor4)
once_append(Tensor5)
once_append(Tensor6)
once_append(Tensor7)
once_append(JacobianTargets)
once_append(Sparse)
once_append(ScatteringMetaData)
once_append(SingleScatteringData)
once_append(TelsemAtlas)
once_append(String)
once_append(Sun)
once_append(CIARecord)
once_append(AtmPoint)
once_append(MuelmatVector)
once_append(MuelmatMatrix)
once_append(PropmatVector)
once_append(PropmatMatrix)
once_append(StokvecVector)
once_append(StokvecMatrix)
once_append(XsecRecord)

once_append(ArrayOfIndex)
once_append(ArrayOfQuantumIdentifier)
once_append(ArrayOfAbsorptionLines)
once_append(ArrayOfTime)
once_append(ArrayOfGriddedField1)
once_append(ArrayOfGriddedField2)
once_append(ArrayOfGriddedField3)
once_append(ArrayOfGriddedField4)
once_append(ArrayOfVector)
once_append(ArrayOfMatrix)
once_append(ArrayOfTensor3)
once_append(ArrayOfTensor4)
once_append(ArrayOfTensor5)
once_append(ArrayOfTensor6)
once_append(ArrayOfTensor7)
once_append(ArrayOfSparse)
once_append(ArrayOfScatteringMetaData)
once_append(ArrayOfSingleScatteringData)
once_append(ArrayOfTelsemAtlas)
once_append(ArrayOfString)
once_append(ArrayOfSun)
once_append(ArrayOfCIARecord)
once_append(ArrayOfAtmPoint)
once_append(ArrayOfMuelmatVector)
once_append(ArrayOfMuelmatMatrix)
once_append(ArrayOfPropmatVector)
once_append(ArrayOfPropmatMatrix)
once_append(ArrayOfStokvecVector)
once_append(ArrayOfStokvecMatrix)
once_append(ArrayOfXsecRecord)
once_append(ArrayOfSpeciesTag)

self_append(ArrayOfIndex)
self_append(ArrayOfQuantumIdentifier)
self_append(ArrayOfAbsorptionLines)
self_append(ArrayOfTime)
self_append(ArrayOfGriddedField1)
self_append(ArrayOfGriddedField2)
self_append(ArrayOfGriddedField3)
self_append(ArrayOfGriddedField4)
self_append(ArrayOfVector)
self_append(ArrayOfMatrix)
self_append(ArrayOfTensor3)
self_append(ArrayOfTensor4)
self_append(ArrayOfTensor5)
self_append(ArrayOfTensor6)
self_append(ArrayOfTensor7)
self_append(ArrayOfSparse)
self_append(ArrayOfScatteringMetaData)
self_append(ArrayOfSingleScatteringData)
self_append(ArrayOfTelsemAtlas)
self_append(ArrayOfString)
self_append(ArrayOfSun)
self_append(ArrayOfCIARecord)
self_append(ArrayOfAtmPoint)
self_append(ArrayOfMuelmatVector)
self_append(ArrayOfMuelmatMatrix)
self_append(ArrayOfPropmatVector)
self_append(ArrayOfPropmatMatrix)
self_append(ArrayOfStokvecVector)
self_append(ArrayOfStokvecMatrix)
self_append(ArrayOfXsecRecord)

self_append(ArrayOfArrayOfGriddedField1)
self_append(ArrayOfArrayOfGriddedField2)
self_append(ArrayOfArrayOfGriddedField3)
self_append(ArrayOfArrayOfString)
self_append(ArrayOfArrayOfScatteringMetaData)
self_append(ArrayOfArrayOfSingleScatteringData)
self_append(ArrayOfArrayOfIndex)
self_append(ArrayOfArrayOfAbsorptionLines)
self_append(ArrayOfArrayOfTime)
self_append(ArrayOfArrayOfMuelmatVector)
self_append(ArrayOfArrayOfMuelmatMatrix)
self_append(ArrayOfArrayOfPropmatVector)
self_append(ArrayOfArrayOfPropmatMatrix)
self_append(ArrayOfArrayOfStokvecVector)
self_append(ArrayOfArrayOfStokvecMatrix)
self_append(ArrayOfArrayOfVector)
self_append(ArrayOfArrayOfMatrix)
self_append(ArrayOfArrayOfTensor3)
self_append(ArrayOfArrayOfTensor6)
self_append(ArrayOfArrayOfSpeciesTag)