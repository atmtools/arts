/** \file mpi/traits.h
 *
 * \brief Contains type traits for a generic MPI implementation.
 *
 */

#ifndef MPI_TRAITS_H
#define MPI_TRAITS_H

namespace invlib
{

// -------------- //
//  Storage Types //
// -------------- //

/*! Passed as second template argument to the MPIMatrix or MPIVector classes,
 * the ConstRef template template triggers the MPIMatrix to only hold a
 * const reference to an already existing matrix or vector.
 */
template
<
typename T
>
struct ConstRef
{
    using type = const T &;
};

/*! Passed as second template argument to the MPIMatrix or MPIVector classes,
 * the LValue template template triggers the MPIMatrix to hold an MPIMatrix as
 * lvalue. This is required if a distributed matrix should be maniulated
 * or created from scratch and not from an already locally existing matrix or vector.
 */
template
<
typename T
>
struct LValue
{
    using type = T;
};

}      // namespace invlib

#endif // MPI_TRATIS_H
