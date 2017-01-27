/*!
  \file   oem_mpi.h
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Thu October 06 16:50:46 2016

  \brief Optimal estimation method for MPI-distributed retrieval.
*/

#ifndef oem_mpi_h
#define oem_mpi_h

#include "invlib/map.h"
#include "invlib/algebra.h"
#include "invlib/optimization.h"
#include "invlib/algebra/solvers.h"
#include "invlib/algebra/precision_matrix.h"
#include "invlib/mpi/mpi_matrix.h"
#include "invlib/mpi/mpi_vector.h"
#include "invlib/mpi/log.h"
#include "invlib/interfaces/arts_wrapper.h"
#include "invlib/profiling/timer.h"

// OEM types.
using MPIMatrix       = invlib::Matrix<invlib::MPIMatrix<invlib::Timer<ArtsMatrix>>>;
using MPISparse       = invlib::Matrix<invlib::MPIMatrix<invlib::Timer<ArtsMatrixReference<const Sparse>>>>;
using MPIVector       = invlib::Vector<invlib::MPIVector<invlib::Timer<ArtsVector>>>;

using PrecisionMPI    = invlib::PrecisionMatrix<MPISparse>;

// Standard Form.
template <typename ForwardModel>
using OEM_PS_PS_MPI = invlib::MAP<ForwardModel, OEMMatrix, PrecisionMPI,
    PrecisionMPI, OEMVector, Formulation::STANDARD>;

// Optimization Methods.
using LM_CG_S_MPI = invlib::LevenbergMarquardt<Numeric, MPISparse, CG>;

#endif // oem_h
