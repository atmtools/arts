#pragma once

#include "array.h"
#include "rtepack_compact_planar_mueller_matrix.h"
#include "rtepack_compact_planar_propagation_matrix.h"
#include "rtepack_compact_planar_stokes_vector.h"
#include "rtepack_jones_matrix.h"
#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_polarization.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_radar.h"
#include "rtepack_rtestep.h"
#include "rtepack_scattering.h"
#include "rtepack_source.h"
#include "rtepack_spectral_matrix.h"
#include "rtepack_stokes_vector.h"
#include "rtepack_surface.h"
#include "rtepack_transmission.h"

using CompactPlanarMuelmat                 = rtepack::compact_planar_muelmat;
using CompactPlanarMuelmatVector           = rtepack::compact_planar_muelmat_vector;
using CompactPlanarMuelmatVectorView       = rtepack::compact_planar_muelmat_vector_view;
using CompactPlanarMuelmatConstVectorView  = rtepack::compact_planar_muelmat_vector_const_view;
using CompactPlanarMuelmatMatrix           = rtepack::compact_planar_muelmat_matrix;
using CompactPlanarMuelmatMatrixView       = rtepack::compact_planar_muelmat_matrix_view;
using CompactPlanarMuelmatConstMatrixView  = rtepack::compact_planar_muelmat_matrix_const_view;
using CompactPlanarMuelmatTensor3          = rtepack::compact_planar_muelmat_tensor3;
using CompactPlanarMuelmatTensor3View      = rtepack::compact_planar_muelmat_tensor3_view;
using CompactPlanarMuelmatConstTensor3View = rtepack::compact_planar_muelmat_tensor3_const_view;
using ArrayOfCompactPlanarMuelmatMatrix    = Array<CompactPlanarMuelmatMatrix>;

using CompactPlanarPropmat                 = rtepack::compact_planar_propmat;
using CompactPlanarPropmatVector           = rtepack::compact_planar_propmat_vector;
using CompactPlanarPropmatVectorView       = rtepack::compact_planar_propmat_vector_view;
using CompactPlanarPropmatConstVectorView  = rtepack::compact_planar_propmat_vector_const_view;
using CompactPlanarPropmatMatrix           = rtepack::compact_planar_propmat_matrix;
using CompactPlanarPropmatMatrixView       = rtepack::compact_planar_propmat_matrix_view;
using CompactPlanarPropmatConstMatrixView  = rtepack::compact_planar_propmat_matrix_const_view;
using CompactPlanarPropmatTensor3          = rtepack::compact_planar_propmat_tensor3;
using CompactPlanarPropmatTensor3View      = rtepack::compact_planar_propmat_tensor3_view;
using CompactPlanarPropmatConstTensor3View = rtepack::compact_planar_propmat_tensor3_const_view;

using CompactPlanarStokvec                 = rtepack::compact_planar_stokvec;
using CompactPlanarStokvecVector           = rtepack::compact_planar_stokvec_vector;
using CompactPlanarStokvecVectorView       = rtepack::compact_planar_stokvec_vector_view;
using CompactPlanarStokvecConstVectorView  = rtepack::compact_planar_stokvec_vector_const_view;
using CompactPlanarStokvecMatrix           = rtepack::compact_planar_stokvec_matrix;
using CompactPlanarStokvecMatrixView       = rtepack::compact_planar_stokvec_matrix_view;
using CompactPlanarStokvecConstMatrixView  = rtepack::compact_planar_stokvec_matrix_const_view;
using CompactPlanarStokvecTensor3          = rtepack::compact_planar_stokvec_tensor3;
using CompactPlanarStokvecTensor3View      = rtepack::compact_planar_stokvec_tensor3_view;
using CompactPlanarStokvecConstTensor3View = rtepack::compact_planar_stokvec_tensor3_const_view;

using Propmat                     = rtepack::propmat;
using PropmatVector               = rtepack::propmat_vector;
using PropmatVectorView           = rtepack::propmat_vector_view;
using PropmatConstVectorView      = rtepack::propmat_vector_const_view;
using PropmatMatrix               = rtepack::propmat_matrix;
using PropmatMatrixView           = rtepack::propmat_matrix_view;
using PropmatConstMatrixView      = rtepack::propmat_matrix_const_view;
using ArrayOfPropmatVector        = Array<PropmatVector>;
using ArrayOfPropmatMatrix        = Array<PropmatMatrix>;
using ArrayOfArrayOfPropmatVector = Array<ArrayOfPropmatVector>;
using ArrayOfArrayOfPropmatMatrix = Array<ArrayOfPropmatMatrix>;

using Stokvec                 = rtepack::stokvec;
using StokvecVector           = rtepack::stokvec_vector;
using StokvecVectorView       = rtepack::stokvec_vector_view;
using StokvecConstVectorView  = rtepack::stokvec_vector_const_view;
using StokvecMatrix           = rtepack::stokvec_matrix;
using StokvecMatrixView       = rtepack::stokvec_matrix_view;
using StokvecConstMatrixView  = rtepack::stokvec_matrix_const_view;
using StokvecTensor3          = rtepack::stokvec_tensor3;
using StokvecTensor3View      = rtepack::stokvec_tensor3_view;
using StokvecConstTensor3View = rtepack::stokvec_tensor3_const_view;
using StokvecTensor4          = rtepack::stokvec_tensor4;
using StokvecTensor4View      = rtepack::stokvec_tensor4_view;
using StokvecConstTensor4View = rtepack::stokvec_tensor4_const_view;
using StokvecTensor5          = rtepack::stokvec_tensor5;
using StokvecTensor5View      = rtepack::stokvec_tensor5_view;
using StokvecConstTensor5View = rtepack::stokvec_tensor5_const_view;
using StokvecTensor6          = rtepack::stokvec_tensor6;
using StokvecTensor6View      = rtepack::stokvec_tensor6_view;
using StokvecConstTensor6View = rtepack::stokvec_tensor6_const_view;

using ArrayOfStokvecVector        = Array<StokvecVector>;
using ArrayOfStokvecMatrix        = Array<StokvecMatrix>;
using ArrayOfStokvecTensor3       = Array<StokvecTensor3>;
using ArrayOfArrayOfStokvecVector = Array<ArrayOfStokvecVector>;
using ArrayOfArrayOfStokvecMatrix = Array<ArrayOfStokvecMatrix>;

using Jones = rtepack::jones;

using Muelmat                     = rtepack::muelmat;
using MuelmatVector               = rtepack::muelmat_vector;
using MuelmatVectorView           = rtepack::muelmat_vector_view;
using MuelmatConstVectorView      = rtepack::muelmat_vector_const_view;
using MuelmatMatrix               = rtepack::muelmat_matrix;
using MuelmatTensor3              = rtepack::muelmat_tensor3;
using MuelmatMatrixView           = rtepack::muelmat_matrix_view;
using MuelmatConstMatrixView      = rtepack::muelmat_matrix_const_view;
using MuelmatTensor4              = rtepack::muelmat_tensor4;
using MuelmatTensor4View          = rtepack::muelmat_tensor4_view;
using MuelmatConstTensor4View     = rtepack::muelmat_tensor4_const_view;
using MuelmatTensor5              = rtepack::muelmat_tensor5;
using MuelmatTensor5View          = rtepack::muelmat_tensor5_view;
using MuelmatConstTensor5View     = rtepack::muelmat_tensor5_const_view;
using MuelmatTensor6              = rtepack::muelmat_tensor6;
using MuelmatTensor6View          = rtepack::muelmat_tensor6_view;
using MuelmatConstTensor6View     = rtepack::muelmat_tensor6_const_view;
using ArrayOfMuelmatVector        = Array<MuelmatVector>;
using ArrayOfMuelmatMatrix        = Array<MuelmatMatrix>;
using ArrayOfMuelmatTensor3       = Array<MuelmatTensor3>;
using ArrayOfArrayOfMuelmatVector = Array<ArrayOfMuelmatVector>;
using ArrayOfArrayOfMuelmatMatrix = Array<ArrayOfMuelmatMatrix>;

using Specmat                     = rtepack::specmat;
using SpecmatVector               = rtepack::specmat_vector;
using SpecmatVectorView           = rtepack::specmat_vector_view;
using SpecmatConstVectorView      = rtepack::specmat_vector_const_view;
using SpecmatMatrix               = rtepack::specmat_matrix;
using SpecmatTensor3              = rtepack::specmat_tensor3;
using SpecmatMatrixView           = rtepack::specmat_matrix_view;
using SpecmatConstMatrixView      = rtepack::specmat_matrix_const_view;
using ArrayOfSpecmatVector        = Array<SpecmatVector>;
using ArrayOfSpecmatMatrix        = Array<SpecmatMatrix>;
using ArrayOfSpecmatTensor3       = Array<SpecmatTensor3>;
using ArrayOfArrayOfSpecmatVector = Array<ArrayOfSpecmatVector>;
using ArrayOfArrayOfSpecmatMatrix = Array<ArrayOfSpecmatMatrix>;

using StokvecSortedGriddedField1 = matpack::gridded_data_t<Stokvec, AscendingGrid>;
using StokvecSortedGriddedField2 = matpack::gridded_data_t<Stokvec, AscendingGrid, AscendingGrid>;
using StokvecSortedGriddedField3 = matpack::gridded_data_t<Stokvec, AscendingGrid, AscendingGrid, AscendingGrid>;
using StokvecSortedGriddedField4 =
    matpack::gridded_data_t<Stokvec, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid>;
using StokvecSortedGriddedField5 =
    matpack::gridded_data_t<Stokvec, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid>;
using StokvecSortedGriddedField6 = matpack::
    gridded_data_t<Stokvec, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid>;
using GriddedSpectralField6 =
    matpack::gridded_data_t<Stokvec, AscendingGrid, LatGrid, LonGrid, ZenGrid, AziGrid, AscendingGrid>;
