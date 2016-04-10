#ifndef FORWARD_MODELS_LINEAR_H
#define FORWARD_MODELS_LINEAR_H

#include <utility.h>

namespace invlib
{

/**
 * \brief Linear forward model for testing.
 *
 * On construction the linear model by generating a random m times n matrix
 * K, that is used as the linear forward model. Implementes all necesary
 * for the computation of the MAP estimator using the MAP class.
 *
 * \tparam MatrixType The matrix type to be used for the linear model.
 */
template
<
typename MatrixType
>
class Linear
{
public:

    Linear(unsigned int n_, unsigned int m_)
        : n(n_), m(m_)
    {
        K = random<MatrixType>(m,n);
    }

    template<typename VectorType>
    VectorType evaluate(const VectorType &x)
    {
        VectorType w = K * x;
        return w;
    }

    template<typename VectorType>
    MatrixType Jacobian(const VectorType &x, VectorType &y)
    {
        MatrixType J = K;
        y = K * x;
        return J;
    }

    const unsigned int m,n;

private:

    MatrixType K;

};

}

#endif // FORWARD_MODELS_LINEAR_H
