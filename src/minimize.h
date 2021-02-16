#ifndef minimize_h
#define minimize_h

#include <array>
#include <utility>

#include "matpackI.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include <unsupported/Eigen/NonLinearOptimization>
#pragma GCC diagnostic pop

namespace Minimize {
//! Functor for minimizing X0 + X1 * SIN(S * T) + X2 * COS(S * T) - Y
template <typename Scalar> struct Wave {
  // typedef needed by the functional-style optimizer
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> InputType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  typedef Eigen::ColPivHouseholderQR<JacobianType> QRSolver;
  
  // Size: The number of parameters to optimize towards
  static constexpr int m_inputs=3;
  
  // Size: The number of inputs required by the 
  const int m_values;
  
  //! Temperature grid
  const Scalar * const T;
  
  //! Spectroscopic values at grid
  const Scalar * const Y;
  
  //! Frequency scaling
  const Scalar S;
  
  //! Save Jacobian
  mutable bool first_jac;
  mutable JacobianType internal_jac;
  
  //! Size: The number of inputs required by the 
  static constexpr int inputs() { return m_inputs; }
  
  //! Size: The number of inputs required by the 
  int values() const { return m_values; }
  
  /*! The only constructor
    * 
    * @param[in] x The grid of the problem
    * @param[in] y The measured/simulated values
    */
  Wave(const std::vector<Scalar>& x,
       const std::vector<Scalar>& y,
       Scalar s) :
       m_values(int(x.size())), T(x.data()), Y(y.data()), S(s), first_jac(true), internal_jac(m_values, m_inputs) {}
  
  /*! The only constructor
    * 
    * @param[in] x The grid of the problem
    * @param[in] y The measured/simulated values
    */
  Wave(const Scalar * const x,
       const Scalar * const y,
       Scalar s,
       int N) :
       m_values(N), T(x), Y(y), S(s), first_jac(true), internal_jac(m_values, m_inputs) {}
       
  template <std::size_t N>
  Wave(const std::array<Scalar, N>& x,
       const std::array<Scalar, N>& y,
       Scalar s) :
       m_values(int(N)), T(x.data()), Y(y.data()), S(s), first_jac(true), internal_jac(m_values, m_inputs) {}
  
  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] f [values()]-sized Vector for model values
    * @return 0
    */
  int operator()(const InputType& p, ValueType& f) const {
    for (int i=0; i<m_values; i++) {
      const double Sin = p[1] * std::sin(S * T[i]);
      const double Cos = p[2] * std::cos(S * T[i]);
      f[i] = p[0] + Sin + Cos - Y[i];
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
    if (first_jac) {
      for (int i=0; i<m_values; i++) {
        internal_jac(i, 0) = 1;
        internal_jac(i, 1) = std::sin(S * T[i]);
        internal_jac(i, 2) = std::cos(S * T[i]);
      }
      first_jac = false;
    }
    
    J = internal_jac;
    return 0;
  }
  
  //! start values helper function, operator()(...) must be not too bad
  InputType x0() const {
    std::vector<Scalar> copyY{Y, Y + m_values};
    std::sort(copyY.begin(), copyY.end());
    InputType x(m_inputs);
    x[0] = copyY[copyY.size() / 2];
    x[1] = x[0] - copyY[copyY.size() / 4];
    x[2] = -x[1];
    return x;
  }
};

//! Functor for minimizing (X0 + X1 (T0 / T - 1)) * (T0 / T) ** X2 - Y
struct T4 {
  // typedef needed by the functional-style optimizer
  typedef Numeric Scalar;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> InputType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  typedef Eigen::ColPivHouseholderQR<JacobianType> QRSolver;
  
  // Size: The number of parameters to optimize towards
  static constexpr int m_inputs=3;
  
  // Size: The number of inputs required by the 
  const int m_values;
  
  //! Temperature grid
  const ConstVectorView& T;
  
  //! Spectroscopic values at grid
  const ConstVectorView& Y;
  
  //! Reference temperature
  const Numeric T0;
  
  //! Reference start exponent
  const Numeric EXP0;
  
  //! Size: The number of inputs required by the 
  static constexpr int inputs() { return m_inputs; }
  
  //! Size: The number of inputs required by the 
  int values() const { return m_values; }
  
  /*! The only constructor
    * 
    * @param[in] x The grid of the problem (Temperatures)
    * @param[in] y The measured/simulated values (Spectroscopic parameter)
    * @param[in] t0 The model's reference temperature
    * @param[in] exp0 Some exponent to start the function at
    */
  T4(const ConstVectorView& x,
     const ConstVectorView& y,
     const Numeric t0,
     const Numeric exp0) :
    m_values(int(x.nelem())), T(x), Y(y), T0(t0), EXP0(exp0) {}
  
  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] f [values()]-sized Vector for model values
    * @return 0
    */
  int operator()(const InputType& p, ValueType& f) const;
    
  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] J [values(), inputs()]-sized Matrix for model Jacobian
    * @return 0
    */
  int df(const InputType& p, JacobianType& J) const;
  
  //! start values helper function, operator()(...) must be not too bad
  InputType x0() const;
};


//! Functor for minimizing X0 * (T0 / T) ** X1 + X2 * (T0 / T) ** X3 - Y
struct DPL {
  // typedef needed by the functional-style optimizer
  typedef Numeric Scalar;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> InputType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  typedef Eigen::ColPivHouseholderQR<JacobianType> QRSolver;
  
  // Size: The number of parameters to optimize towards
  static constexpr int m_inputs=4;
  
  // Size: The number of inputs required by the 
  const int m_values;
  
  //! Temperature grid
  const ConstVectorView& T;
  
  //! Spectroscopic values at grid
  const ConstVectorView& Y;
  
  //! Reference temperature
  const Numeric T0;
  
  //! Reference start exponent
  const Numeric EXP0;
  
  //! Size: The number of inputs required by the 
  static constexpr int inputs() { return m_inputs; }
  
  //! Size: The number of inputs required by the 
  int values() const { return m_values; }
  
  /*! The only constructor
    * 
    * @param[in] x The grid of the problem (Temperatures)
    * @param[in] y The measured/simulated values (Spectroscopic parameter)
    * @param[in] t0 The model's reference temperature
    * @param[in] exp0 Some exponent to start the function at
    */
  DPL(const ConstVectorView& x,
      const ConstVectorView& y,
      const Numeric t0,
      const Numeric exp0) :
    m_values(int(x.nelem())), T(x), Y(y), T0(t0), EXP0(exp0) {}
  
  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] f [values()]-sized Vector for model values
    * @return 0
    */
  int operator()(const InputType& p, ValueType& f) const;
    
  /*!  Opeartor evaluating the function
    * 
    * @param[in] p inputs()-sized parameter list
    * @param[in,out] J [values(), inputs()]-sized Matrix for model Jacobian
    * @return 0
    */
  int df(const InputType& p, JacobianType& J) const;
  
  //! start values helper function, operator()(...) must be not too bad
  InputType x0() const;
};


/*! Returns wether or not the Eigen minimize call worked
 * 
 * @param[in] status Eigen::LevenbergMarquardtSpace::Status value from minimize(x) call 
 */
constexpr bool goodStatus(int status) {
  if (status == Eigen::LevenbergMarquardtSpace::RelativeErrorAndReductionTooSmall) {
    return true;
  } else if (status == Eigen::LevenbergMarquardtSpace::RelativeReductionTooSmall) {
    return true;
  } else if (status == Eigen::LevenbergMarquardtSpace::RelativeErrorTooSmall) {
    return true;
  } else {
    return false;
  }
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
 * Note that the functions can also be declared as static and/or constexpr,
 * but the constexpr support and the constant-sized matrices/vectors did
 * not work within the solver at the time of implementing this.
 * 
 * @param[in] fun A Functor fulfilling the above conditions
 * @return {true if successful curve fit else false, list of parameters that minimizes the functor}
 */
template <class Functor> 
std::pair<bool, typename Functor::InputType>
  curve_fit(const Functor& fun) {
  typename Functor::InputType p = fun.x0();
  Eigen::LevenbergMarquardt lm(fun);
  const auto status = lm.minimize(p);
  return {goodStatus(status), p};
}


/*! Fit a curve to data values
 * 
 * See curve_fit(Functor(X, Y, std::forward<Inputs>(ins)...)) for more details.
 * 
 * In addition to the Functor definition of that signature,
 * the Functor must also have defined the implied constructor:
 * 
 * struct Functor {
 * 
 *  Functor(const ConstVectorView& X, const ConstVectorView& Y, std::forward<Inputs>(ins)...);
 * 
 * };
 *
 * The Functor itself is initialized inside curve_fit(X, Y, ...)
 * with a call to its constructor with all Inputs simply forwarded
 * as received.  This initialization is to a constant so the
 * interface functions must be acceptable from const-marked
 * calls.
 * 
 * @param[in] X The X-Grid
 * @param[in] Y The computed/measured/simulated values at X
 * @param[in] ins Inputs to forward to the constructor of Functor
 * @return See curve_fit(Functor(X, Y, std::forward<Inputs>(ins)...))
 */
template <class Functor, typename ... Inputs> 
std::pair<bool, typename Functor::InputType>
  curve_fit(const ConstVectorView& X,
            const ConstVectorView& Y,
            Inputs&& ... ins) {
  return curve_fit(Functor(X, Y, std::forward<Inputs>(ins)...));
}
}

#endif  // minimize_wrap_h
