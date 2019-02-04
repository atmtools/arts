#ifndef CONVERGENCE_CRITERIA_H
#define CONVERGENCE_CRITERIA_H

#include <memory>

/** file convergence_criteria.h
 * \brief Functors to computer different convergence criteria.
 *
 * This file contains different convergence criteria to test for
 * convergence for different iteration methods.
 *
 * The convergence functor must be callable with the following
 * arguments:
 *   - xi: The state vector of the current iteration.
 *   - yi: The simulated measurement vector corresponding to yi.
 *   -  y: The observed measurement vector.
 *   -  g: The gradient of the loss function.
 *   -  K: The forward model jacobian.
 *   - Sa: The a priori covariance matrix.
 *   - Se: The measurement covariance matrix.
 *
 * And return a bool indicating if convergence has been reached or not.
 * That is, it must implement operator() with the following signature:
 *
 * bool operator(Vector xi, Vector yi, Vector yk
 *               Matrix K, Matrix Sa, Matrix Se)
 *
 */
namespace invlib
{

template
<
    typename VectorType
>
struct Rodgers531
{
    using RealType = typename VectorType::RealType;

    template<typename JacobianType, typename SaType, typename SeType>
    auto operator()(const VectorType & xi,
                    const VectorType & /*yi*/,
                    const VectorType & /*y*/,
                    const VectorType & g,
                    const JacobianType & /*K*/,
                    const SaType & /*Sa*/,
                    const SeType & /*Se*/)
        -> typename VectorType::RealType
    {
        if (!x_im1_ptr) {
            x_im1_ptr = std::unique_ptr<VectorType>(new VectorType(xi));
            return std::numeric_limits<RealType>::max();
        }

        // Order here must be reversed compared to Rodgers because
        // g here is the gradient of the cost function.
        VectorType dx = *x_im1_ptr - xi;
        auto conv = abs(dot(dx, g) / xi.rows());
        return conv;
    }

    std::unique_ptr<VectorType> x_im1_ptr = nullptr;

};

template
<
    typename VectorType
>
struct Rodgers530
{
    using RealType = typename VectorType::RealType;

    template<typename JacobianType, typename SaType, typename SeType>
    auto operator()(const VectorType & xi,
                    const VectorType & yi,
                    const VectorType & y,
                    const VectorType & g,
                    const JacobianType & K,
                    const SaType & Sa,
                    const SeType & Se)
    -> typename VectorType::RealType
    {
        if (!x_im1_ptr) {
            x_im1_ptr = std::unique_ptr<VectorType>(new VectorType(xi));
            return std::numeric_limits<RealType>::max();
        }


        VectorType dx = xi - *x_im1_ptr;
        VectorType Hdx = (transp(K) * inv(Se) * K + inv(Sa)) * dx;
        auto conv = dot(dx, Hdx) / xi.rows();

        *x_im1_ptr = xi;

        return conv;
    }

    std::unique_ptr<VectorType> x_im1_ptr = nullptr;
};

template
<
typename VectorType
>
struct Rodgers533
{
    using RealType = typename VectorType::RealType;

    template<typename JacobianType, typename SaType, typename SeType>
        auto operator()(VectorType & xi,
                        VectorType & yi,
                        VectorType & y,
                        VectorType & g,
                        JacobianType & K,
                        SaType & Sa,
                        SeType & Se)
        -> typename VectorType::RealType
        {
            if (!y_im1_ptr) {
                y_im1_ptr = std::unique_ptr<VectorType>(new VectorType(yi));
                return std::numeric_limits<RealType>::max();
            }

            VectorType dy = yi - *y_im1_ptr;
            auto Sdy_inv  = inv(Se) * (K * Sa * transp(K) + Se) * inv(Se);
            auto conv     = dot(dy, Sdy_inv * dy) / y.rows();

            *y_im1_ptr = yi;

            return conv;
        }

    std::unique_ptr<VectorType> y_im1_ptr = nullptr;
};

}      // namespace invlib
#endif // CONVERGENCE_CRITERIA_H
