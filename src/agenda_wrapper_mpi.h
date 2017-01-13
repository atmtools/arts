/*!
  \file   agenda_wrappers_mpi.h
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Thu October 06 17:11:00 2016

  \brief Wrappers for the inversion_iterate_agenda in order to them as invlib
  forward models for MPI-distributed retrievals.
*/

#ifndef agenda_wrappers_mpi_h
#define agenda_wrappers_mpi_h

#include "oem.h"
#include "oem_mpi.h"
#include "mpi.h"

//! Wrapper class for forward model using a distributed Jacobian with MPI.
/*!
  Each process only computes a limited range of rows of the measurement vector
  y and the Jacobian. The measurement vector is returned as a full vector to each
  process which requires broadcasting the results after the computation.
 */
class AgendaWrapperMPI
{
    Workspace *ws;
    OEMMatrix local_jacobian;
    const Agenda *inversion_iterate_agenda;

public:

    const unsigned int m,n;
    AgendaWrapperMPI(Workspace *ws_,
                     const Agenda *inversion_iterate_agenda_,
                     Index m_, Index n_) :
        ws(ws_), inversion_iterate_agenda(inversion_iterate_agenda_), m(m_),
        n(n_), local_jacobian()
    {}

    MPIMatrix Jacobian(const OEMVector &xi, OEMVector &yi)
    {
        yi.resize(m);
        inversion_iterate_agendaExecute( *ws, yi, local_jacobian, xi,
                                         1, *inversion_iterate_agenda );
        // Create MPI vector from local results, use conversion to vector
        // to broadcast local results.
        MPIVector yi_mpi(yi);
        yi = yi_mpi;

        MPIMatrix jacobian = local_jacobian;
        return jacobian;
    }

    OEMVector evaluate(const OEMVector &xi)
    {
        Matrix dummy = local_jacobian;
        OEMVector yi; yi.resize(m);
        inversion_iterate_agendaExecute( *ws, yi, dummy, xi, 0,
                                         *inversion_iterate_agenda );

        // Create MPI vector from local results, use conversion to vector
        // to broadcast local results.
        MPIVector yi_mpi = yi;
        yi = yi_mpi;
        return yi;
    }
};

#endif // agenda_wrapper_mpi_h
