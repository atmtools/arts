/** \file mpi/utility.h
 *
 * \brief Utility functions for handling distributed matrices.
 *
 */

#ifndef MPI_MPI_UTILITY_H
#define MPI_MPI_UTILITY_H

#include "traits.h"

namespace invlib
{

// --------------- //
//   Reductions    //
// --------------- //

template<typename T1>
T1 mpi_sum(T1 &t)
{
    T1 result;
    MPI_Allreduce(&t, &result, 1, MPIDataType<T1>::name, MPI_SUM, MPI_COMM_WORLD);
    return result;
}

}

      // namespace invlib

#endif // MPI_MPI_UTILITY_H
