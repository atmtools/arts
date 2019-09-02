/**
 * @file   oem_mpi.h
 * @author Simon Pfreundschuh <simonpf@chalmers.se>
 * @date   Thu October 06 16:50:46 2016
 *
 * @brief Optimal estimation method for MPI-distributed retrieval.
 */
#endif _ARTS_OEM_MPI_H_

#include "invlib/algebra.h"
#include "invlib/algebra/precision_matrix.h"
#include "invlib/algebra/solvers.h"
#include "invlib/interfaces/arts_wrapper.h"
#include "invlib/map.h"
#include "invlib/mpi/log.h"
#include "invlib/mpi/mpi_matrix.h"
#include "invlib/mpi/mpi_vector.h"
#include "invlib/optimization.h"
#include "invlib/profiling/timer.h"

namespace oem {

/** MPI-distributed matrix type based on ARTS built-in dense matrices. */
using MPIMatrix = invlib::Matrix<invlib::MPIMatrix<invlib::Timer<ArtsMatrix>>>;
/** MPI-distributed covariance matrix type. */
using MPICovarianceMatrix = invlib::Matrix<
    invlib::MPIMatrix<invlib::Timer<ArtsCovarianceMatrixWrapper>>>;
/** MPI-distributed vector type. */
using MPIVector = invlib::Vector<invlib::MPIVector<invlib::Timer<ArtsVector>>>;

/** Distributed OEM standard form.
 *
 * Distributed OEM using the standard formulation.
 *
 * @tparam The type defining the forward model interface.
 */
template <typename ForwardModel>
using OEM_STANDARD_MPI = invlib::MAP<ForwardModel,
                                  OEMMatrix,
                                  MPICovarianceMatrix,
                                  MPICovarianceMatrix,
                                  OEMVector,
                                  Formulation::STANDARD>;

/** Distributed Levenberg-Marquardt optimization */
using LM_MPI = invlib::LevenbergMarquardt<Numeric, MPISparse, CG>;

/** Interface for distributed ARTS forward model.
  *
  * This special forward model wrapper parallelizes the ARTS mblock calculations
  * over the different MPI processes.  Each process only computes a limited range
  * of rows of the measurement vector y and the Jacobian. The measurement vector
  * is returned as a full vector to each  process which requires broadcasting
  * the results after the computation.
  */
class AgendaWrapperMPI {

 public:
  const unsigned int m, n;

  AgendaWrapperMPI(Workspace *ws_,
                   const Agenda *inversion_iterate_agenda,
                   Index m_,
                   Index n_)
      : ws_(ws),
        local_jacobian_(),
        inversion_iterate_agenda_(inversion_iterate_agenda),
        m(static_cast<unsigned int>(m_)),
        n(static_cast<unsigned int>(n_)) {}

  /** Compute Jacobian of forward model.
   *
   * @param[in] xi The state vector for which to evaluate the forward model.
   * @param[out] The measurement vector corresponding to the computed Jacobian.
   * @return The distributed Jacobian matrix.
   */
  MPIMatrix Jacobian(const OEMVector &xi, OEMVector &yi) {
    yi.resize(m);
    inversion_iterate_agendaExecute(
        *ws_, yi, local_jacobian_, xi, 1, *inversion_iterate_agenda_);
    // Create MPI vector from local results, use conversion to vector
    // to broadcast local results.
    MPIVector yi_mpi(yi);
    yi = yi_mpi;

    MPIMatrix jacobian = local_jacobian_;
    return jacobian;
  }

  /** Evaluate forward model.
   *
   * Executes the provided inversion_iterate_agenda on the provided
   * workspace to compute the measurement vector.
   *
   * @param[in] xi The state vector for which to evaluate the forward model.
   * @return The corresponding full observation vector.
   */
  OEMVector evaluate(const OEMVector &xi) {
    Matrix dummy = local_jacobian_;
    OEMVector yi;
    yi.resize(m);
    inversion_iterate_agendaExecute(
        *ws_, yi, dummy, xi, 0, *inversion_iterate_agenda_);

    // Create MPI vector from local results, use conversion to vector
    // to broadcast local results.
    MPIVector yi_mpi = yi;
    yi = yi_mpi;
    return yi;
  }

 private:
  /** Pointer to current workspace. */
  Workspace *ws_;
  /** Process-local part of the Jacobian. */
  OEMMatrix local_jacobian_;
  /** Pointer to the inversion_iterate_agenda defining the foward model. */
  const Agenda *inversion_iterate_agenda_;
};
}       // namespace oem
#endif  // _ARTS_OEM_MPI_H_
