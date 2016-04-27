/** \file mpi/mpi_vector.h
 *
 * \brief Contains the MPIVector class, a generic class for vectors distributed
 * row-wise over nodes.
 *
 */

#ifndef MPI_MPI_VECTOR_H
#define MPI_MPI_VECTOR_H

#include "mpi.h"
#include "invlib/mpi/traits.h"

namespace invlib
{

template
<
typename LocalType,
template <typename> typename StorageTrait = ConstRef
>
class MPIVector
{

public:

    /*! The basic scalar type. */
    using RealType   = typename LocalType::RealType;
    /*! The basic vector type  */
    using VectorType = typename LocalType::VectorType;
    /*! The local Matrix type.  */
    using MatrixType = LocalType;
    /*!
     * Result type of an algebraic expression with MPIMatrix as right hand
     * operator.
     */
    using ResultType = LocalType;
    /*! The type used to store the local vector. */
    using StorageType = typename StorageTrait<LocalType>::type;

    MPIVector();

    template
    <
    typename T,
    typename = enable_if<is_constructible<StorageType, T>>
    >
    MPIVector(T &&local_vector);

    void resize(unsigned int i);

    unsigned int rows() const;

    LocalType & get_local();

    const LocalType & get_local() const;

    RealType operator()(unsigned int i) const;
    RealType& operator()(unsigned int i);

    LocalType broadcast() const;

    operator LocalType() const;

private:

    static constexpr MPI_Datatype mpi_data_type = MPI_DOUBLE;

    void broadcast_local_rows(int proc_rows[]) const;
    void broadcast_local_block(double *vector,
                               const double *block) const;

    std::vector<unsigned int> row_indices;
    std::vector<unsigned int> row_ranges;

    int rank, nprocs;

    unsigned int m;
    unsigned int local_rows;

    StorageType local;
    RealType local_element;

};

#include "mpi_vector.cpp"

}      // namespace invlib

#endif // MPI_MPI_MATRIX_H
