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

// ------------------- //
//   MPI Data Types    //
// ------------------- //

template <typename T>
struct MPIDataType;

template <>
struct MPIDataType<double>
{
public:
    static constexpr MPI_Datatype name = MPI_DOUBLE;
};

template <>
struct MPIDataType<float>
{
public:
    static constexpr MPI_Datatype value = MPI_FLOAT;
};

// ---------------- //
//      MPI Type    //
// ---------------- //

// Forward declarations.

template
<
typename Base
>
class Vector;

template
<
typename LocalType,
template <typename> typename StorageTrait
>
class MPIVector;

template
<
typename Base
>
class Matrix;

template
<
typename LocalType,
template <typename> typename StorageTrait
>
class MPIMatrix;

// MPIType struct.

template
<
typename T1,
template <typename> typename StorageType
>
struct MPITypeStruct;

template
<
typename T1,
template <typename> typename StorageType
>
struct MPITypeStruct<Vector<T1>, StorageType>
{
public:
    using type = Vector<MPIVector<T1, StorageType>>;
};

template
<
typename T1,
template <typename> typename StorageType
>
struct MPITypeStruct<Matrix<T1>, StorageType>
{
public:
    using type = Matrix<MPIMatrix<T1, StorageType>>;
};

// Type alias.

template
<
typename T1,
template <typename> typename StorageType
>
using MPIType = typename MPITypeStruct<T1, StorageType>::type;

}      // namespace invlib


#endif // MPI_TRATIS_H
