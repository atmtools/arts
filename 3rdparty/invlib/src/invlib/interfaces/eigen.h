#include <Eigen/Sparse>
#include <utility>

#include "invlib/traits.h"

using EigenSparseBase   = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using EigenVectorBase   = Eigen::VectorXd;

class EigenSparse;

// ----------------- //
//   Eigen Vector    //
// ----------------- //

class EigenVector : protected EigenVectorBase
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
    typename = invlib::enable_if<invlib::is_constructible<EigenVectorBase, T>>
    >
    EigenVector(T &&t)
        : EigenVectorBase(std::forward<T>(t))
    {
        // Nothing to do here.
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
        return this->block((int) start, 0, (int) extent, 1);
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

class EigenSparse : protected EigenSparseBase
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
    typename = invlib::enable_if<invlib::is_constructible<EigenVectorBase, T>>
    >
    EigenSparse(T &&t)
        : EigenSparseBase(std::forward<T>(t))
    {
        // Nothing to do here.
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

};
