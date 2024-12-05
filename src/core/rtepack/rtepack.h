#pragma once

#include "array.h"
#include "rtepack_rtestep.h"
#include "rtepack_scattering.h"
#include "rtepack_source.h"
#include "rtepack_transmission.h"
#include "rtepack_spectral_matrix.h"

using Propmat = rtepack::propmat;
using PropmatVector = rtepack::propmat_vector;
using PropmatVectorView = rtepack::propmat_vector_view;
using PropmatConstVectorView = rtepack::propmat_vector_const_view;
using PropmatMatrix = rtepack::propmat_matrix;
using PropmatMatrixView = rtepack::propmat_matrix_view;
using PropmatConstMatrixView = rtepack::propmat_matrix_const_view;
using ArrayOfPropmatVector = Array<PropmatVector>;
using ArrayOfPropmatMatrix = Array<PropmatMatrix>;
using ArrayOfArrayOfPropmatVector = Array<ArrayOfPropmatVector>;
using ArrayOfArrayOfPropmatMatrix = Array<ArrayOfPropmatMatrix>;

using Stokvec = rtepack::stokvec;
using StokvecVector = rtepack::stokvec_vector;
using StokvecVectorView = rtepack::stokvec_vector_view;
using StokvecConstVectorView = rtepack::stokvec_vector_const_view;
using StokvecMatrix = rtepack::stokvec_matrix;
using StokvecMatrixView = rtepack::stokvec_matrix_view;
using StokvecConstMatrixView = rtepack::stokvec_matrix_const_view;
using StokvecTensor3 = rtepack::stokvec_tensor3;
using StokvecTensor3View = rtepack::stokvec_tensor3_view;
using StokvecConstTensor3View = rtepack::stokvec_tensor3_const_view;
using StokvecTensor4 = rtepack::stokvec_tensor4;
using StokvecTensor4View = rtepack::stokvec_tensor4_view;
using StokvecConstTensor4View = rtepack::stokvec_tensor4_const_view;
using StokvecTensor5 = rtepack::stokvec_tensor5;
using StokvecTensor5View = rtepack::stokvec_tensor5_view;
using StokvecConstTensor5View = rtepack::stokvec_tensor5_const_view;
using StokvecTensor6 = rtepack::stokvec_tensor6;
using StokvecTensor6View = rtepack::stokvec_tensor6_view;
using StokvecConstTensor6View = rtepack::stokvec_tensor6_const_view;
using StokvecGriddedField6 = matpack::gridded_data<Stokvec,
                                                   AscendingGrid,
                                                   AscendingGrid,
                                                   AscendingGrid,
                                                   AscendingGrid,
                                                   AscendingGrid,
                                                   AscendingGrid>;
using ArrayOfStokvecVector = Array<StokvecVector>;
using ArrayOfStokvecMatrix = Array<StokvecMatrix>;
using ArrayOfStokvecTensor3 = Array<StokvecTensor3>;
using ArrayOfArrayOfStokvecVector = Array<ArrayOfStokvecVector>;
using ArrayOfArrayOfStokvecMatrix = Array<ArrayOfStokvecMatrix>;

using Muelmat = rtepack::muelmat;
using MuelmatVector = rtepack::muelmat_vector;
using MuelmatVectorView = rtepack::muelmat_vector_view;
using MuelmatConstVectorView = rtepack::muelmat_vector_const_view;
using MuelmatMatrix = rtepack::muelmat_matrix;
using MuelmatTensor3 = rtepack::muelmat_tensor3;
using MuelmatMatrixView = rtepack::muelmat_matrix_view;
using MuelmatConstMatrixView = rtepack::muelmat_matrix_const_view;
using ArrayOfMuelmatVector = Array<MuelmatVector>;
using ArrayOfMuelmatMatrix = Array<MuelmatMatrix>;
using ArrayOfMuelmatTensor3 = Array<MuelmatTensor3>;
using ArrayOfArrayOfMuelmatVector = Array<ArrayOfMuelmatVector>;
using ArrayOfArrayOfMuelmatMatrix = Array<ArrayOfMuelmatMatrix>;

using Specmat = rtepack::specmat;
using SpecmatVector = rtepack::specmat_vector;
using SpecmatVectorView = rtepack::specmat_vector_view;
using SpecmatConstVectorView = rtepack::specmat_vector_const_view;
using SpecmatMatrix = rtepack::specmat_matrix;
using SpecmatTensor3 = rtepack::specmat_tensor3;
using SpecmatMatrixView = rtepack::specmat_matrix_view;
using SpecmatConstMatrixView = rtepack::specmat_matrix_const_view;
using ArrayOfSpecmatVector = Array<SpecmatVector>;
using ArrayOfSpecmatMatrix = Array<SpecmatMatrix>;
using ArrayOfSpecmatTensor3 = Array<SpecmatTensor3>;
using ArrayOfArrayOfSpecmatVector = Array<ArrayOfSpecmatVector>;
using ArrayOfArrayOfSpecmatMatrix = Array<ArrayOfSpecmatMatrix>;

namespace rtepack {
std::ostream& operator<<(std::ostream& os, const Array<propmat>& a);
std::ostream& operator<<(std::ostream& os, const Array<Array<propmat>>& a);
std::ostream& operator<<(std::ostream& os, const Array<muelmat>& a);
std::ostream& operator<<(std::ostream& os, const Array<Array<muelmat>>& a);
std::ostream& operator<<(std::ostream& os, const Array<stokvec>& a);
std::ostream& operator<<(std::ostream& os, const Array<Array<stokvec>>& a);
}  // namespace rtepack
