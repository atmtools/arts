#pragma once

#include "array.h"
#include "rtepack_rtestep.h"
#include "rtepack_source.h"
#include "rtepack_transmission.h"

using Propmat = rtepack::propmat;
using PropmatVector = rtepack::propmat_vector;
using PropmatVectorView = rtepack::propmat_vector_view;
using PropmatConstVectorView = rtepack::propmat_vector_const_view;
using PropmatMatrix = rtepack::propmat_matrix;
using PropmatMatrixView = rtepack::propmat_matrix_view;
using PropmatConstMatrixView = rtepack::propmat_matrix_const_view;
using ArrayOfPropmatVector = Array<PropmatVector>;
using ArrayOfPropmatMatrix = Array<PropmatMatrix>;

using Stokvec = rtepack::stokvec;
using StokvecVector = rtepack::stokvec_vector;
using StokvecVectorView = rtepack::stokvec_vector_view;
using StokvecConstVectorView = rtepack::stokvec_vector_const_view;
using StokvecMatrix = rtepack::stokvec_matrix;
using StokvecMatrixView = rtepack::stokvec_matrix_view;
using StokvecConstMatrixView = rtepack::stokvec_matrix_const_view;
using ArrayOfStokvecVector = Array<StokvecVector>;
using ArrayOfStokvecMatrix = Array<StokvecMatrix>;
using ArrayOfArrayOfStokvecVector = Array<ArrayOfStokvecVector>;
using ArrayOfArrayOfStokvecMatrix = Array<ArrayOfStokvecMatrix>;

using Muelmat = rtepack::muelmat;
using MuelmatVector = rtepack::muelmat_vector;
using MuelmatVectorView = rtepack::muelmat_vector_view;
using MuelmatConstVectorView = rtepack::muelmat_vector_const_view;
using MuelmatMatrix = rtepack::muelmat_matrix;
using MuelmatMatrixView = rtepack::muelmat_matrix_view;
using MuelmatConstMatrixView = rtepack::muelmat_matrix_const_view;
using ArrayOfMuelmatVector = Array<MuelmatVector>;
using ArrayOfMuelmatMatrix = Array<MuelmatMatrix>;
using ArrayOfArrayOfMuelmatVector = Array<ArrayOfMuelmatVector>;
using ArrayOfArrayOfMuelmatMatrix = Array<ArrayOfMuelmatMatrix>;
