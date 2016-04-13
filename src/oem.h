/*!
  \file   oem.h
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Fri March 25 15:53:54 2016

  \brief Optimal estimation method for retrieval.
*/

#ifndef oem_h
#define oem_h

#include "invlib/map.h"
#include "invlib/algebra.h"
#include "invlib/optimization.h"
#include "invlib/algebra/solvers.h"
#include "invlib/algebra/precision_matrix.h"
#include "invlib/interfaces/arts_wrapper.h"
#include "invlib/profiling/timer.h"

using OEMVector       = invlib::Vector<invlib::Timer<ArtsVector>>;
using OEMMatrix       = invlib::Matrix<invlib::Timer<ArtsMatrix>>;
using OEMSparse       = invlib::Matrix<invlib::Timer<ArtsSparse>>;

/* using OEMVector       = invlib::Vector<ArtsVector>; */
/* using OEMMatrix       = invlib::Matrix<ArtsMatrix>; */
/* using OEMSparse       = invlib::Matrix<ArtsSparse>; */

using Identity        = invlib::MatrixIdentity<OEMMatrix>;
using Transform       = invlib::NormalizeDiagonal<OEMMatrix>;
using PrecisionMatrix = invlib::PrecisionMatrix<OEMMatrix>;
using PrecisionSparse = invlib::PrecisionMatrix<OEMSparse>;

// OEM types.

using invlib::Formulation;

// Standard Form.

template <typename ForwardModel>
using OEM_S_S = invlib::MAP<ForwardModel, OEMMatrix, Sparse,
    Sparse, Formulation::STANDARD>;

template <typename ForwardModel>
using OEM_PS_PS = invlib::MAP<ForwardModel, OEMMatrix, PrecisionSparse,
    PrecisionSparse, Formulation::STANDARD>;

template <typename ForwardModel>
using OEM_D_D = invlib::MAP<ForwardModel, OEMMatrix, OEMMatrix,
    OEMMatrix, Formulation::STANDARD>;

template <typename ForwardModel>
using OEM_PD_PD = invlib::MAP<ForwardModel, OEMMatrix, PrecisionMatrix,
    PrecisionMatrix, Formulation::STANDARD>;

// N-Form.

template <typename ForwardModel>
using OEM_NFORM_S_S = invlib::MAP<ForwardModel, OEMMatrix, Sparse,
    Sparse, Formulation::NFORM>;

template <typename ForwardModel>
using OEM_NFORM_PS_PS = invlib::MAP<ForwardModel, OEMMatrix, PrecisionSparse,
    PrecisionSparse, Formulation::NFORM>;

template <typename ForwardModel>
using OEM_NFORM_D_D = invlib::MAP<ForwardModel, OEMMatrix, OEMMatrix,
    OEMMatrix, Formulation::NFORM>;

template <typename ForwardModel>
using OEM_NFORM_PD_PD = invlib::MAP<ForwardModel, OEMMatrix, PrecisionMatrix,
    PrecisionMatrix, Formulation::NFORM>;

// M-Form.

template <typename ForwardModel>
using OEM_MFORM_S_S = invlib::MAP<ForwardModel, OEMMatrix, Sparse,
    Sparse, Formulation::MFORM>;

template <typename ForwardModel>
using OEM_MFORM_PS_PS = invlib::MAP<ForwardModel, OEMMatrix, PrecisionSparse,
    PrecisionSparse, Formulation::MFORM>;

template <typename ForwardModel>
using OEM_MFORM_D_D = invlib::MAP<ForwardModel, OEMMatrix, OEMMatrix,
    OEMMatrix, Formulation::MFORM>;

template <typename ForwardModel>
using OEM_MFORM_PD_PD = invlib::MAP<ForwardModel, OEMMatrix, PrecisionMatrix,
    PrecisionMatrix, Formulation::MFORM>;

// Solvers.
using Pre     = invlib::NormalizeDiagonal<OEMMatrix>;
using Std     = invlib::Standard;
using Std_Pre = invlib::PreconditionedSolver<Std, Pre>;
using CG      = invlib::ConjugateGradient;
using CG_Pre  = invlib::PreconditionedSolver<CG, Pre>;

// Optimization Methods.
using GN        = invlib::GaussNewton<Numeric>;
using GN_CG     = invlib::GaussNewton<Numeric, CG>;
using GN_Pre    = invlib::GaussNewton<Numeric, Std_Pre>;
using GN_CG_Pre = invlib::GaussNewton<Numeric, CG_Pre>;
using LM_D        = invlib::LevenbergMarquardt<Numeric, OEMMatrix>;
using LM_CG_D     = invlib::LevenbergMarquardt<Numeric, OEMMatrix, CG>;
using LM_Pre_D    = invlib::LevenbergMarquardt<Numeric, OEMMatrix, Std_Pre>;
using LM_CG_Pre_D = invlib::LevenbergMarquardt<Numeric, OEMMatrix, CG_Pre>;
using LM_Sparse_D = invlib::LevenbergMarquardt<Numeric, OEMMatrix>;
using LM_S        = invlib::LevenbergMarquardt<Numeric, OEMSparse>;
using LM_CG_S     = invlib::LevenbergMarquardt<Numeric, OEMSparse, CG>;
using LM_Pre_S    = invlib::LevenbergMarquardt<Numeric, OEMSparse, Std_Pre>;
using LM_I        = invlib::LevenbergMarquardt<Numeric, Identity>;
using LM_CG_I     = invlib::LevenbergMarquardt<Numeric, Identity, CG>;
using LM_Pre_I    = invlib::LevenbergMarquardt<Numeric, Identity, Std_Pre>;
using LM_CG_Pre_S = invlib::LevenbergMarquardt<Numeric, OEMSparse, CG_Pre>;
using LM_Sparse_S = invlib::LevenbergMarquardt<Numeric, OEMSparse>;

#endif // oem_h
