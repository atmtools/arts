#include "minimize.h"

#ifndef _MSC_VER
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

#include <unsupported/Eigen/NonLinearOptimization>

#ifndef _MSC_VER
#pragma GCC diagnostic pop
#endif

namespace Minimize {
//! Functor for minimizing X0 + X1 * X + X2 * X**2 + ... XN * X**N - Y
struct Polynom {
  // typedef needed by the functional-style optimizer
  using Scalar       = Numeric;
  using InputType    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using ValueType    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using JacobianType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using QRSolver     = Eigen::ColPivHouseholderQR<JacobianType>;

  // Size: The number of parameters to optimize towards
  const int m_inputs;

  // Size: The number of inputs required by the
  const int m_values;

  //! Temperature grid
  const StridedConstVectorView& X;

  //! Spectroscopic values at grid
  const StridedConstVectorView& Y;

  //! Size: The number of inputs required by the
  [[nodiscard]] int inputs() const { return m_inputs; }

  //! Size: The number of inputs required by the
  [[nodiscard]] int values() const { return m_values; }

  /*! The only constructor
    * 
    * @param[in] x The grid of the problem
    * @param[in] y The measured/simulated values
    * @param[in] order The max degree of the polynominal (e.g., 3 means x**3 is largest factor)
    */
  Polynom(const StridedConstVectorView& x,
          const StridedConstVectorView& y,
          const Index order)
      : m_inputs(int(order + 1)), m_values(int(x.size())), X(x), Y(y) {}

  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] f [values()]-sized Vector for model values
    * @return 0
    */
  int operator()(const InputType& p, ValueType& f) const {
    for (Index i = 0; i < m_values; i++) {
      f[i]      = p[0] - Y[i];
      Numeric x = 1;
      for (int j = 1; j < m_inputs; j++) {
        x    *= X[i];
        f[i] += p[j] * x;
      }
    }
    return 0;
  }

  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] J [values(), inputs()]-sized Matrix for model Jacobian
    * @return 0
    */
  int df(const InputType&, JacobianType& J) const {
    for (Index i = 0; i < m_values; i++) {
      J(i, 0)   = 1;
      Numeric x = 1;
      for (int j = 1; j < m_inputs; j++) {
        x       *= X[i];
        J(i, j)  = x;
      }
    }
    return 0;
  }

  //! start values helper function, operator()(...) must be not too bad
  InputType x0() const {
    InputType out(m_inputs);
    for (int j = 0; j < m_inputs; j++) {
      out[j] = 1;
    }
    return out;
  }
};

//! Functor for minimizing (X0 + X1 (T0 / T - 1)) * (T0 / T) ** X2 - Y
struct T4 {
  // typedef needed by the functional-style optimizer
  using Scalar       = Numeric;
  using InputType    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using ValueType    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using JacobianType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using QRSolver     = Eigen::ColPivHouseholderQR<JacobianType>;

  // Size: The number of parameters to optimize towards
  static constexpr int m_inputs = 3;

  // Size: The number of inputs required by the
  const int m_values;

  //! Temperature grid
  const StridedConstVectorView& T;

  //! Spectroscopic values at grid
  const StridedConstVectorView& Y;

  //! Reference temperature
  const Numeric T0;

  //! Reference start exponent
  const Numeric EXP0;

  //! Size: The number of inputs required by the
  static constexpr int inputs() { return m_inputs; }

  //! Size: The number of inputs required by the
  [[nodiscard]] int values() const { return m_values; }

  /*! The only constructor
    * 
    * @param[in] x The grid of the problem (Temperatures)
    * @param[in] y The measured/simulated values (Spectroscopic parameter)
    * @param[in] t0 The model's reference temperature
    * @param[in] exp0 Some exponent to start the function at
    */
  T4(const StridedConstVectorView& x,
     const StridedConstVectorView& y,
     const Numeric t0,
     const Numeric exp0)
      : m_values(int(x.size())), T(x), Y(y), T0(t0), EXP0(exp0) {}

  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] f [values()]-sized Vector for model values
    * @return 0
    */
  int operator()(const InputType& p, ValueType& f) const {
    for (Index i = 0; i < m_values; i++) {
      const Numeric G  = T0 / T[i];
      const Numeric GX = nonstd::pow(G, p[2]);
      f[i]             = (p[0] + p[1] * (G - 1)) * GX - Y[i];
    }
    return 0;
  }

  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] J [values(), inputs()]-sized Matrix for model Jacobian
    * @return 0
    */
  int df(const InputType& p, JacobianType& J) const {
    for (Index i = 0; i < m_values; i++) {
      const Numeric G  = T0 / T[i];
      const Numeric GX = nonstd::pow(G, p[2]);
      J(i, 0)          = GX;
      J(i, 1)          = (G - 1) * GX;
      J(i, 2)          = (p[0] + p[1] * (G - 1)) * GX * std::log(G);
    }
    return 0;
  }

  //! start values helper function, operator()(...) must be not too bad
  InputType x0() const {
    const Numeric mean_y = mean(Y);
    InputType out(m_inputs);
    out << mean_y, -0.01 * mean_y, mean_y < 0 ? -EXP0 : EXP0;
    return out;
  }
};

//! Functor for minimizing X0 * (T0 / T) ** X1 + X2 * (T0 / T) ** X3 - Y
struct DPL {
  // typedef needed by the functional-style optimizer
  using Scalar       = Numeric;
  using InputType    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using ValueType    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using JacobianType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using QRSolver     = Eigen::ColPivHouseholderQR<JacobianType>;

  // Size: The number of parameters to optimize towards
  static constexpr int m_inputs = 4;

  // Size: The number of inputs required by the
  const int m_values;

  //! Temperature grid
  const StridedConstVectorView& T;

  //! Spectroscopic values at grid
  const StridedConstVectorView& Y;

  //! Reference temperature
  const Numeric T0;

  //! Reference start exponent
  const Numeric EXP0;

  //! Size: The number of inputs required by the
  static constexpr int inputs() { return m_inputs; }

  //! Size: The number of inputs required by the
  [[nodiscard]] int values() const { return m_values; }

  /*! The only constructor
    * 
    * @param[in] x The grid of the problem (Temperatures)
    * @param[in] y The measured/simulated values (Spectroscopic parameter)
    * @param[in] t0 The model's reference temperature
    * @param[in] exp0 Some exponent to start the function at
    */
  DPL(const StridedConstVectorView& x,
      const StridedConstVectorView& y,
      const Numeric t0,
      const Numeric exp0)
      : m_values(int(x.size())), T(x), Y(y), T0(t0), EXP0(exp0) {}

  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] f [values()]-sized Vector for model values
    * @return 0
    */
  int operator()(const InputType& p, ValueType& f) const {
    for (Index i = 0; i < m_values; i++) {
      const Numeric G   = T0 / T[i];
      const Numeric GX1 = nonstd::pow(G, p[1]);
      const Numeric GX3 = nonstd::pow(G, p[3]);
      f[i]              = p[0] * GX1 + p[2] * GX3 - Y[i];
    }
    return 0;
  }

  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] J [values(), inputs()]-sized Matrix for model Jacobian
    * @return 0
    */
  int df(const InputType& p, JacobianType& J) const {
    for (Index i = 0; i < m_values; i++) {
      const Numeric G   = T0 / T[i];
      const Numeric lG  = std::log(G);
      const Numeric GX1 = nonstd::pow(G, p[1]);
      const Numeric GX3 = nonstd::pow(G, p[3]);
      J(i, 0)           = GX1;
      J(i, 1)           = p[0] * GX1 * lG;
      J(i, 2)           = GX3;
      J(i, 3)           = p[2] * GX3 * lG;
    }
    return 0;
  }

  //! start values helper function, operator()(...) must be not too bad
  InputType x0() const {
    const Numeric mean_y = mean(Y);
    InputType out(m_inputs);
    out << mean_y, mean_y < 0 ? -EXP0 : EXP0, -0.01 * mean_y,
        mean_y < 0 ? EXP0 : -EXP0;
    return out;
  }
};

namespace {
bool goodStatus(int status) {
  return status != Eigen::LevenbergMarquardtSpace::ImproperInputParameters;
}

/*! Fit a curve to data values
 * 
 * The Functor is required to have this signature:
 * 
 * struct Functor {
 * 
 * typedef Numeric Scalar;
 * 
 * typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> InputType;
 * 
 * typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ValueType;
 * 
 * typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
 * 
 * typedef Eigen::ColPivHouseholderQR<JacobianType> QRSolver;
 * 
 * int operator()(const InputType& p, ValueType& f) const;
 * 
 * int df(const InputType& p, JacobianType& J) const;
 * 
 * InputType x0() const;
 * 
 * int inputs() const;
 * 
 * int values() const;
 * 
 * };
 *
 * Note that this is expensive to compile.  It is recommended to use the
 * same approach as for polyfit below to avoid recompilation.
 * 
 * @param[in] fun A Functor fulfilling the above conditions
 * @return An optional Vector if successful, otherwise an empty optional
 */
template <class Functor>
std::optional<Vector> curve_fit(const Functor& fun) {
  typename Functor::InputType p = fun.x0();
  Eigen::LevenbergMarquardt lm(fun);
  const auto status = lm.minimize(p);
  if (not goodStatus(status)) return std::nullopt;
  return Vector{p};
}
}  // namespace

std::optional<Vector> polyfit(const StridedConstVectorView& X,
                              const StridedConstVectorView& Y,
                              const Index& order) {
  return curve_fit(Polynom(X, Y, order));
}
}  // namespace Minimize
