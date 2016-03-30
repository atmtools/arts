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
 * \tparam Matrix The matrix type to be used for the linear model.
 */
template
<
typename Matrix
>
class Linear
{
public:

    Linear(unsigned int n_, unsigned int m_)
        : n(n_), m(m_)
    {
        K = random<Matrix>(m,n);
    }

    template<typename Vector>
    Vector evaluate(const Vector &x)
    {
        Vector w = K * x;
        return w;
    }

    template<typename Vector>
    Matrix Jacobian(const Vector &x)
    {
        Matrix J = K;
        return J;
    }

    const unsigned int m,n;

private:

    Matrix K;

};

}

#endif // FORWARD_MODELS_LINEAR_H
