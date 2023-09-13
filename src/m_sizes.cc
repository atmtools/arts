#include <auto_wsg.h>

#include <exception>

#define NELEM_GET(T)                                          \
  void nelemGet(Index& size, const T& x) { size = x.size(); } \
  void IndexSetToLast(Index& last, const T& x) {              \
    nelemGet(last, x);                                        \
    last--;                                                   \
  }

NELEM_GET(ArrayOfAbsorptionLines);
NELEM_GET(ArrayOfAgenda);
NELEM_GET(ArrayOfArrayOfAbsorptionLines);
NELEM_GET(ArrayOfArrayOfGriddedField1);
NELEM_GET(ArrayOfArrayOfGriddedField2);
NELEM_GET(ArrayOfArrayOfGriddedField3);
NELEM_GET(ArrayOfArrayOfIndex);
NELEM_GET(ArrayOfArrayOfMatrix);
NELEM_GET(ArrayOfArrayOfMuelmatMatrix);
NELEM_GET(ArrayOfArrayOfMuelmatVector);
NELEM_GET(ArrayOfArrayOfPropmatMatrix);
NELEM_GET(ArrayOfArrayOfPropmatVector);
NELEM_GET(ArrayOfArrayOfScatteringMetaData);
NELEM_GET(ArrayOfArrayOfSingleScatteringData);
NELEM_GET(ArrayOfArrayOfSpeciesTag);
NELEM_GET(ArrayOfArrayOfStokvecMatrix);
NELEM_GET(ArrayOfArrayOfStokvecVector);
NELEM_GET(ArrayOfArrayOfString);
NELEM_GET(ArrayOfArrayOfTensor3);
NELEM_GET(ArrayOfArrayOfTensor6);
NELEM_GET(ArrayOfArrayOfTime);
NELEM_GET(ArrayOfArrayOfVector);
NELEM_GET(ArrayOfAtmPoint);
NELEM_GET(ArrayOfCIARecord);
NELEM_GET(ArrayOfGriddedField1);
NELEM_GET(ArrayOfGriddedField2);
NELEM_GET(ArrayOfGriddedField3);
NELEM_GET(ArrayOfGriddedField4);
NELEM_GET(ArrayOfIndex);
NELEM_GET(ArrayOfJacobianTarget);
NELEM_GET(ArrayOfMatrix);
NELEM_GET(ArrayOfMuelmatMatrix);
NELEM_GET(ArrayOfMuelmatVector);
NELEM_GET(ArrayOfPpath);
NELEM_GET(ArrayOfPropmatMatrix);
NELEM_GET(ArrayOfPropmatVector);
NELEM_GET(ArrayOfQuantumIdentifier);
NELEM_GET(ArrayOfRetrievalQuantity);
NELEM_GET(ArrayOfScatteringMetaData);
NELEM_GET(ArrayOfSingleScatteringData);
NELEM_GET(ArrayOfSparse);
NELEM_GET(ArrayOfSpeciesTag);
NELEM_GET(ArrayOfStokvecMatrix);
NELEM_GET(ArrayOfStokvecVector);
NELEM_GET(ArrayOfString);
NELEM_GET(ArrayOfSun);
NELEM_GET(ArrayOfTelsemAtlas);
NELEM_GET(ArrayOfTensor3);
NELEM_GET(ArrayOfTensor4);
NELEM_GET(ArrayOfTensor5);
NELEM_GET(ArrayOfTensor6);
NELEM_GET(ArrayOfTensor7);
NELEM_GET(ArrayOfTime);
NELEM_GET(ArrayOfVector);
NELEM_GET(ArrayOfXsecRecord);
NELEM_GET(Vector);

#define NCOLS_GET(T) \
  void ncolsGet(Index& ncols, const T& x) { ncols = x.ncols(); }

NCOLS_GET(Matrix)
NCOLS_GET(Sparse)
NCOLS_GET(Tensor3)
NCOLS_GET(Tensor4)
NCOLS_GET(Tensor5)
NCOLS_GET(Tensor6)
NCOLS_GET(Tensor7)

#define NROWS_GET(T) \
  void nrowsGet(Index& nrows, const T& x) { nrows = x.nrows(); }

NROWS_GET(Matrix)
NROWS_GET(Sparse)
NROWS_GET(Tensor3)
NROWS_GET(Tensor4)
NROWS_GET(Tensor5)
NROWS_GET(Tensor6)
NROWS_GET(Tensor7)

#define NPAGES_GET(T) \
  void npagesGet(Index& npages, const T& x) { npages = x.npages(); }

NPAGES_GET(Tensor3)
NPAGES_GET(Tensor4)
NPAGES_GET(Tensor5)
NPAGES_GET(Tensor6)
NPAGES_GET(Tensor7)

#define NBOOKS_GET(T) \
  void nbooksGet(Index& nbooks, const T& x) { nbooks = x.nbooks(); }

NBOOKS_GET(Tensor4)
NBOOKS_GET(Tensor5)
NBOOKS_GET(Tensor6)
NBOOKS_GET(Tensor7)

#define NSHELVES_GET(T) \
  void nshelvesGet(Index& nshelves, const T& x) { nshelves = x.nshelves(); }

NSHELVES_GET(Tensor5)
NSHELVES_GET(Tensor6)
NSHELVES_GET(Tensor7)

#define NVITRINES_GET(T) \
  void nvitrinesGet(Index& nvitrines, const T& x) { nvitrines = x.nvitrines(); }

NVITRINES_GET(Tensor6)
NVITRINES_GET(Tensor7)

#define NLIBRARIES_GET(T) \
  void nlibrariesGet(Index& nlibraries, const T& x) { nlibraries = x.nlibraries(); }

NLIBRARIES_GET(Tensor7)
