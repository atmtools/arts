#ifndef INTERFACES_EIGEN_H
#define INTERFACES_EIGEN_H

#include <Eigen/Sparse>
#include <utility>

#include "invlib/traits.h"
#include "invlib/dense/vector_data.h"
#include "invlib/sparse/sparse_data.h"

namespace invlib {

using EigenSparseBase   = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using EigenVectorBase   = Eigen::VectorXd;

class EigenSparse;

// ----------------- //
//   Eigen Vector    //
// ----------------- //

class EigenVector : public EigenVectorBase
{

public:

    // -------------- //
    //  Type Aliases  //
    // -------------- //

    using BaseType   = EigenVectorBase;
    using RealType   = double;
    using VectorType = EigenVector;
    using MatrixType = EigenSparse;
    using ResultType = EigenVector;

    // -------------- //
    //  Constructors  //
    // -------------- //

    EigenVector()                    = default;
    EigenVector(const EigenVector &) = default;

    template
    <
    typename T,
    typename = disable_if<is_invlib_expression<T>>
    >
    EigenVector(T &&t)
        : EigenVectorBase(std::forward<T>(t))
    {
        // Nothing to do here.
    }

    // ----------------- //
    //     Conversion    //
    // ----------------- //

    EigenVector(const VectorData<double> & data)
    {
        BaseType::resize(data.rows());
        for (size_t i = 0; i < data.rows(); i++)
        {
            BaseType::operator()(i) = data.get_element_pointer()[i];
        }
    }

    operator VectorData<double>() const
    {
        auto elements  = std::shared_ptr<double *>(new (double *),
                                                   ArrayDeleter<double *>());
        *elements = new double[rows()];
        for (size_t i = 0; i < rows(); i++)
        {
            (*elements)[i] = this->operator()(i);
        }
        return VectorData<double>(rows(), elements);
    }

    // ------------------- //
    //     Manipulation    //
    // ------------------- //

    unsigned int rows() const
    {
        return this->EigenVectorBase::rows();
    }

    void resize(unsigned int n)
    {
        this->EigenVectorBase::resize((int) n);
    }

    EigenVector get_block(unsigned int start, unsigned int extent) const
    {
        return EigenVector(this->block((int) start, 0, (int) extent, 1));
    }

    RealType * data_pointer()
    {
        return this->data();
    }

    const RealType * data_pointer() const
    {
        return this->data();
    }

    // ---------------------- //
    //  Arithmetic Operations //
    // ---------------------- //

    void accumulate(const EigenVector& v)
    {
        *this += v;
    }

    void subtract(const EigenVector& v)
    {
        *this -= v;
    }

    void scale(RealType c)
    {
        *this *= c;
    }

    VectorType element_multiply(const VectorType &v) const
    {
        VectorType w = cwiseProduct(v);
        return w;
    }

    void element_invert()
    {
        cwiseInverse();
    }

    RealType norm() const
    {
        return this->EigenVectorBase::norm();
    }

    friend EigenSparse;
    friend double dot(const EigenVector &v, const EigenVector &w);

};

double dot(const EigenVector &v, const EigenVector &w)
{
    return v.dot(w);
}

// ----------------- //
//   Eigen Sparse    //
// ----------------- //

class EigenSparse : public EigenSparseBase
{

public:

    // -------------- //
    //  Type Aliases  //
    // -------------- //

    using BaseType   = EigenSparseBase;
    using RealType   = double;
    using VectorType = EigenVector;
    using MatrixType = EigenSparse;
    using ResultType = EigenSparse;

    // -------------- //
    //  Constructors  //
    // -------------- //

    template
    <
    typename T,
    typename = invlib::disable_if<is_invlib_expression<T>>
    >
    EigenSparse(T &&t)
        : EigenSparseBase(std::forward<T>(t))
    {
        // Nothing to do here.
    }

    template <typename Index>
    EigenSparse(const SparseData<double, Index, Representation::Coordinates> & data)
    {
        std::vector<Eigen::Triplet<double>> triplets{};

        size_t nnz = data.non_zeros();
        triplets.reserve(nnz);
        for (size_t i = 0; i < nnz; i++)
        {
            triplets.emplace_back(data.get_row_index_pointer()[i],
                                  data.get_column_index_pointer()[i],
                                  data.get_element_pointer()[i]);
        }
        BaseType::resize(data.rows(), data.cols());
        BaseType::setFromTriplets(triplets.begin(), triplets.end());
    }

    // ------------------- //
    //     Manipulation    //
    // ------------------- //

    unsigned int rows() const
    {
        return this->EigenSparseBase::rows();
    }

    unsigned int cols() const
    {
        return this->EigenSparseBase::cols();
    }

    EigenSparse get_block(unsigned int row_start,
                          unsigned int col_start,
                          unsigned int row_extent,
                          unsigned int col_extent) const
    {
        return this->block((int) row_start,
                           (int) col_start,
                           (int) row_extent,
                           (int) col_extent);
    }

    // ---------------------- //
    //  Arithmetic Operations //
    // ---------------------- //

    VectorType multiply(const VectorType &v) const
    {
        VectorType w = *this * v;
        return w;
    }

    VectorType transpose_multiply(const VectorType &v) const
    {
        VectorType w = this->transpose() * v;
        return w;
    }

    VectorType transpose_multiply_block(
        const VectorType &v,
        size_t block_start,
        size_t block_length) const
    {
        VectorType w = this->transpose() * v.block(block_start, 0, block_length, 1);
        return w;
    }

};

}      // namespace invlib
#endif // INTERFACES_EIGEN_H
