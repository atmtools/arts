#ifndef FORWARD_MODELS_SPHERE_H
#define FORWARD_MODELS_SPHERE_H

#include <cmath>
#include <iostream>

template
<
typename Matrix
>
class Sphere
{
public:

    Sphere(unsigned int n_)
        : n(n_), m(1) {}

    template<typename VectorType>
    VectorType evaluate(const VectorType& x)
    {
        VectorType w; w.resize(m);
        w(0) = dot(x,x);
        return w;
    }

    template<typename VectorType>
    Matrix Jacobian(const VectorType& x, VectorType &y)
    {
        y.resize(m);
        y(0) = dot(x,x);

        Matrix J; J.resize(m,n);
        for (unsigned int i = 0; i < n; i++)
            J(0,i) = 2.0 * x(i);
        return J;
    }

    const unsigned int n, m;
};

#endif // FORWARD_MODELS_SPHERE_H
