/** \file scattering/integration.h
 *
 * Quadratures and integration functions.
 *
 * @author Simon Pfreundschuh, 2020
 */
#include <map>
#include <numbers>
#include "fftw3.h"
#include "scattering/maths.h"

#pragma once

namespace scattering {

using std::numbers::pi_v;

namespace detail {


////////////////////////////////////////////////////////////////////////////////
/// Type trait storing the desired precision for different FP types.
////////////////////////////////////////////////////////////////////////////////

template <std::floating_point Scalar>
struct Precision {
  static constexpr Scalar value = 1e-16;
};

template <>
struct Precision<float> {
  static constexpr float value = static_cast<float>(1e-6);
};

template <>
struct Precision<long double> {
  static constexpr long double value = 1e-19;
};
}  // namespace detail


/// Enum class to represent different quadrature classes.
enum class QuadratureType {
  Trapezoidal,
  GaussLegendre,
  DoubleGauss,
  Lobatto,
  ClenshawCurtis,
  Fejer
};

////////////////////////////////////////////////////////////////////////////////
// Gauss-Legendre Quadrature
////////////////////////////////////////////////////////////////////////////////

// pxx :: export
// pxx :: instance(["double"])
/** Gauss-Legendre Quadarature
 *
 * This class implements a Gauss-Legendre quadrature for the integration
 * of functions over the interval [-1, 1].
 */
template <std::floating_point Scalar>
class GaussLegendreQuadrature {
 private:
  /** Find Gauss-Legendre nodes and weights.
   *
   * Uses the Newton root finding algorithm to find the roots of the
   * Legendre polynomial of degree n. Legendre functions are evaluated
   * using a recursion relation.
   */
  // pxx :: hide
  void calculate_nodes_and_weights() {
    const long int n = degree_;
    const long int n_half_nodes = (n + 1) / 2;
    const long int n_max_iter = 10;
    Scalar x, x_old, p_l, p_l_1, p_l_2, dp_dx;

    for (int i = 1; i <= n_half_nodes; ++i) {
      p_l = pi_v<Scalar>;
      p_l_1 = 2 * n;
      //
      // Initial guess.
      //
      x = -(1.0 - (n - 1) / (p_l_1 * p_l_1 * p_l_1)) *
          cos((p_l * (4 * i - 1)) / (4 * n + 2));

      //
      // Evaluate Legendre Polynomial and its derivative at node.
      //
      for (int j = 0; j < n_max_iter; ++j) {
        p_l = x;
        p_l_1 = 1.0;
        for (int l = 2; l <= n; ++l) {
          // Legendre recurrence relation
          p_l_2 = p_l_1;
          p_l_1 = p_l;
          p_l = ((2.0 * l - 1.0) * x * p_l_1 - (l - 1.0) * p_l_2) / l;
        }
        dp_dx = ((1.0 - x) * (1.0 + x)) / (n * (p_l_1 - x * p_l));
        x_old = x;

        //
        // Perform Newton step.
        //
        x -= p_l * dp_dx;
        auto dx = x - x_old;
        if (math::small(std::abs(dx * (x + x_old)))) {
          break;
        }
      }
      nodes_[i - 1] = x;
      weights_[i - 1] = 2.0 * dp_dx * dp_dx / ((1.0 - x) * (1.0 + x));
      nodes_[n - i] = -x;
      weights_[n - i] = weights_[i - 1];
    }
  }

 public:
  GaussLegendreQuadrature() {}
  GaussLegendreQuadrature(int degree)
      : degree_(degree), nodes_(degree), weights_(degree) {
    calculate_nodes_and_weights();
  }

  static constexpr QuadratureType type = QuadratureType::GaussLegendre;

  const math::Vector<Scalar>& get_nodes() const { return nodes_; }
  const math::Vector<Scalar>& get_weights() const { return weights_; }

  int degree_;
  math::Vector<Scalar> nodes_;
  math::Vector<Scalar> weights_;
};

////////////////////////////////////////////////////////////////////////////////
// Double Gauss Quadrature
////////////////////////////////////////////////////////////////////////////////

// pxx :: export
// pxx :: instance(["double"])
/** F
 *
 * This class implements a Gauss-Legendre for the integration of
 * functions of the interval [-1, 1].
 */
template <typename Scalar>
class DoubleGaussQuadrature {
 private:

 public:
  DoubleGaussQuadrature() {}
  DoubleGaussQuadrature(int degree)
      : degree_(degree), nodes_(degree), weights_(degree) {
      assert(degree % 2 == 0);
      auto gq = GaussLegendreQuadrature<Scalar>(degree / 2);
      auto nodes = gq.get_nodes();
      auto weights = gq.get_weights();

      for (int i = 0; i < degree / 2; ++i) {
          nodes_[i] = -0.5 + nodes[i] / 2.0;
          nodes_[degree / 2 + i] = 0.5 + nodes[i] / 2.0;
          weights_[i] = 0.5 * weights[i];
          weights_[degree / 2 + i] = 0.5 *  weights[i];
      }
  }

  static constexpr QuadratureType type = QuadratureType::DoubleGauss;

  const math::Vector<Scalar>& get_nodes() const { return nodes_; }
  const math::Vector<Scalar>& get_weights() const { return weights_; }

  int degree_;
  math::Vector<Scalar> nodes_;
  math::Vector<Scalar> weights_;
};

// pxx :: export
// pxx :: instance(["double"])
/** Nodes and weights of the Lobatto quadrature.
 *
 * This class calculates the weights and nodes of a Lobatto quadrature. The
 * Lobatto quadrature includes the endpoints [-1, 1] of the integration interval
 * and the internal nodes are derived from Legendre polynomials.
 *
 * More information can be found here: https://mathworld.wolfram.com/LobattoQuadrature.html
 print(weights, weights_ref)
 *
 */
template <std::floating_point Scalar>
class LobattoQuadrature {

 private:
  // pxx :: hide
  /// Calculate nodes and weights of the Lobatto quadrature.
  void calculate_nodes_and_weights() {
    const long int n = degree_;
    const long int n_half_nodes = (n + 1) / 2;
    const long int n_max_iter = 10;
    Scalar x, x_old, p_l, p_l_1, p_l_2, dp_dx, d2p_dx;

    const long int left = (n + 1) / 2 - 1;
    const long int right = n / 2;

    for (int i = 0; i < n_half_nodes - 1; ++i) {

      //
      // Initial guess.
      //
      Scalar d_i = ((n % 2) == 0) ? 0.5 : 0.0;
      x = sin(pi_v<Scalar> * (i + d_i) / (n - 0.5));

      //
      // Evaluate Legendre Polynomial and its derivative at node.
      //
      for (int j = 0; j < n_max_iter; ++j) {
        p_l = x;
        p_l_1 = 1.0;
        for (int l = 2; l < n; ++l) {
          // Legendre recurrence relation
          p_l_2 = p_l_1;
          p_l_1 = p_l;
          p_l = ((2.0 * l - 1.0) * x * p_l_1 - (l - 1.0) * p_l_2) / l;
        }
        dp_dx = (n - 1) * (x * p_l - p_l_1) / (x * x - 1.0);
        d2p_dx = (2.0 * x * dp_dx - (n - 1) * n * p_l) / (1.0 - x * x);

        x_old = x;

        //
        // Perform Newton step.
        //
        x -= dp_dx / d2p_dx;
        auto dx = x - x_old;
        if (math::small(std::abs(dx * (x + x_old)))) {
          break;
        }
      }
      nodes_[right + i] = x;
      weights_[right + i] = 2.0 / (n * (n - 1) * p_l * p_l);
      nodes_[left - i] = -x;
      weights_[left - i] = weights_[right + i];
    }
    nodes_[0] = -1.0;
    weights_[0] = 2.0 / (n * (n - 1));
    nodes_[n - 1] = - nodes_[0];
    weights_[n - 1] = weights_[0];
  }

 public:
  LobattoQuadrature() {}
  LobattoQuadrature(int degree)
      : degree_(degree), nodes_(degree), weights_(degree) {
    calculate_nodes_and_weights();
  }

  static constexpr QuadratureType type = QuadratureType::Lobatto;

  const math::Vector<Scalar>& get_nodes() const { return nodes_; }
  const math::Vector<Scalar>& get_weights() const { return weights_; }

  int degree_;
  math::Vector<Scalar> nodes_;
  math::Vector<Scalar> weights_;
};

////////////////////////////////////////////////////////////////////////////////
// Regular grid
////////////////////////////////////////////////////////////////////////////////

// pxx :: export
// pxx :: instance(["double"])
/** Clenshaw-Curtis quadrature
 *
 * This class implements the Clenshaw-Curtis quadrature rule, which uses
 * the points
 *  $x_i = arccos(2 * \pi * i / N)$ for $i = 0,\dots, N$.
 *  as integration nodes.
 */
template <std::floating_point Scalar>
class ClenshawCurtisQuadrature {
 private:
  // pxx :: hide
  void calculate_nodes_and_weights() {
    long int n = degree_ - 1;
    fftw_plan ifft;
    double* weights = reinterpret_cast<double*>(
        fftw_malloc(2 * (n / 2 + 1) * sizeof(double)));
    std::complex<double>* coeffs =
        reinterpret_cast<std::complex<double>*>(weights);

    ifft = fftw_plan_dft_c2r_1d(n,
                                reinterpret_cast<double(*)[2]>(coeffs),
                                weights,
                                FFTW_ESTIMATE);
    // Calculate DFT input.
    for (int i = 0; i < n / 2 + 1; ++i) {
      coeffs[i] = 2.0 / (1.0 - 4.0 * i * i);
    }
    fftw_execute_dft_c2r(ifft, reinterpret_cast<double(*)[2]>(coeffs), weights);

    weights[0] *= 0.5;
    for (int i = 0; i < n; ++i) {
      weights_[i] = weights[i] / n;
    }
    weights_[n] = weights[0];
    fftw_destroy_plan(ifft);
    fftw_free(weights);

    // Calculate nodes.
    for (long int i = 0; i <= n; i++) {
      nodes_[i] = -cos((pi_v<Scalar> * i) / n);
    }
  }

 public:
  ClenshawCurtisQuadrature() {}
  ClenshawCurtisQuadrature(int degree)
      : degree_(degree), nodes_(degree), weights_(degree) {
    calculate_nodes_and_weights();
  }

  static constexpr QuadratureType type = QuadratureType::ClenshawCurtis;

  const math::Vector<Scalar>& get_nodes() const { return nodes_; }
  const math::Vector<Scalar>& get_weights() const { return weights_; }

  int degree_;
  math::Vector<Scalar> nodes_;
  math::Vector<Scalar> weights_;
};

////////////////////////////////////////////////////////////////////////////////
// Fejer quadrature
////////////////////////////////////////////////////////////////////////////////

// pxx :: export
// pxx :: instance(["double"])
/** Fejer quadrature.
 *
 * This class implements the Fejer quadrature, which is similar to the
 * Clenshaw-Curtis quadrature but uses the points
 *  $x_i = arccos(2 * \pi * (i + 0.5) / N)$ for $i = 0,\dots, N-1$.
 * as integration nodes.
 */
template <std::floating_point Scalar>
class FejerQuadrature {
 private:

  // pxx :: hide
  void calculate_nodes_and_weights() {

      int n = degree_;
      fftw_plan ifft;
      double *weights = reinterpret_cast<double*>(fftw_malloc(2 * (n + 1) * sizeof(double)));
      std::complex<double> *coeffs = reinterpret_cast<std::complex<double>*>(weights);

      ifft = fftw_plan_dft_c2r_1d(n,
                                  reinterpret_cast<double (*)[2]>(coeffs),
                                  weights,
                                  FFTW_ESTIMATE);
      // Calculate DFT input.
      for (int i = 0; i < n / 2 + 1; ++i) {
          Scalar x = (pi_v<Scalar> * i) / static_cast<Scalar>(n);
          coeffs[i] = std::complex<double>(cos(x), sin(x));
          coeffs[i] *= 2.0 / (1.0 - 4.0 * i * i);
      }
      fftw_execute_dft_c2r(ifft,
                           reinterpret_cast<double (*)[2]>(coeffs),
                           weights);
      for (long int i = 0; i < n; ++i) {
          weights_[n - i - 1] = weights[i] / n;
      }

      fftw_destroy_plan(ifft);
      fftw_free(weights);

      // Calculate nodes.
      for (long int i = 0; i < n; i++) {
          nodes_[i] = -cos(pi_v<Scalar> * (static_cast<double>(i) + 0.5) / static_cast<double>(n));
      }
  }

 public:

  static constexpr QuadratureType type = QuadratureType::Fejer;

  FejerQuadrature() {}
  FejerQuadrature(int degree)
    : degree_(degree), nodes_(degree), weights_(degree) {
    calculate_nodes_and_weights();
  }

  const math::Vector<Scalar>& get_nodes() const { return nodes_; }
  const math::Vector<Scalar>& get_weights() const { return weights_; }


  int degree_;
  math::Vector<Scalar> nodes_;
  math::Vector<Scalar> weights_;
};

template <std::floating_point Scalar, template<typename> typename Quadrature>
class QuadratureProvider {
 public:
  QuadratureProvider() {}

  Quadrature<Scalar> get_quadrature(int degree) {
    auto found = quadratures_.find(degree);
    if (found != quadratures_.end()) {
      return found->second;
    } else {
      quadratures_.insert({degree, Quadrature<Scalar>(degree)});
      return quadratures_[degree];
    }
  }

 private:
  std::map<int, Quadrature<Scalar>> quadratures_;
};

/** Base class for latitude grids.
 *
 * This class defines the basic interface for latitude grids.
 * It is used to store the co-latitude values of the grid points and,
 * in addition to that, the integration weights of the corresponding
 * quadrature.
 */
template <std::floating_point Scalar>
class LatitudeGrid : public math::Vector<Scalar> {
 public:
  LatitudeGrid() : math::Vector<Scalar>() {}
  LatitudeGrid(const math::Vector<Scalar>& latitudes)
      : math::Vector<Scalar>(latitudes) {}

  virtual ~LatitudeGrid() {};

  /// The latitude grid points in radians.
  virtual const math::Vector<Scalar>& get_colatitudes() const = 0;
  /// The latitude grid points in radians.
  virtual const math::Vector<Scalar>& get_latitudes() const {return *this;}
  /// The integration weights.
  virtual const math::Vector<Scalar>& get_weights() const = 0;

  /// The type of quadrature.
  virtual QuadratureType get_type() = 0;
};

template <std::floating_point Scalar>
using LatitudeGridPtr = std::shared_ptr<LatitudeGrid<Scalar>>;

// pxx :: instance(["double"])
template <std::floating_point Scalar>
class IrregularLatitudeGrid : public LatitudeGrid<Scalar> {

 public:
    using math::Vector<Scalar>::operator[];
  IrregularLatitudeGrid() {}
  /** Create new latitude grid.
   * @param latitudes Vector containing the latitude grid points in radians.
   * @param weights The integration weight corresponding to each grid point.
   */
  IrregularLatitudeGrid(const math::Vector<Scalar>& latitudes)
      : LatitudeGrid<Scalar>(latitudes),
        weights_(latitudes.size()),
      colatitudes_(std::make_unique<math::Vector<Scalar>>(-Eigen::cos(latitudes.array()))),
        type_(QuadratureType::Trapezoidal) {
    weights_.setConstant(0.0);
    int n = math::Vector<Scalar>::size();
    for (int i = 0; i < n - 1; ++i) {
        auto dx = 0.5 * (this->operator[](i + 1) - this->operator[](i));
        weights_[i] += dx;
        weights_[i + 1] += dx;
    }
    weights_[0] += this->operator[](0) + 1.0;
    weights_[n - 1] += 1.0 - this->operator[](n - 1);
  }

  /// The latitude grid points in radians.
  const math::Vector<Scalar>& get_colatitudes() const { return *colatitudes_; }

  /// The latitude grid points in radians.
  const math::Vector<Scalar>& get_weights() const { return weights_; }

  /// The type of quadrature.
  QuadratureType get_type() { return QuadratureType::Trapezoidal; }

 protected:
  math::Vector<Scalar> weights_;
  std::unique_ptr <math::Vector<Scalar>> colatitudes_;
  QuadratureType type_;
};

// pxx :: export
// pxx :: instance("GaussLegendreLatitudeGrid", ["scattering::GaussLegendreQuadrature<double>", "double"])
// pxx :: instance("DoubleGaussLatitudeGrid", ["scattering::DoubleGaussQuadrature<double>", "double"])
// pxx :: instance("LobattoLatitudeGrid", ["scattering::LobattoQuadrature<double>", "double"])
template <typename Quadrature, std::floating_point Scalar>
class QuadratureLatitudeGrid : public LatitudeGrid<Scalar>
{
public:
  /** Create new quadrature latitude grid with given number of points.
   *
   * Creates a latitude grid using the nodes and weights of the given quadrature
   * class as grid points.
   *
   * @tparam Quadrature The quadrature class to use to calculate the nodes and
   * weights of the quadrature.
   * @param degree The number of points of the quadrature.
   */
  QuadratureLatitudeGrid() : LatitudeGrid<Scalar>() {}
  QuadratureLatitudeGrid(int n_points) : quadrature_(n_points) {
    math::Vector<Scalar>::operator=(
        Eigen::acos(quadrature_.get_nodes().array() * -1.0));
  }

 QuadratureLatitudeGrid(size_t n_points, size_t /*unused*/)
     : QuadratureLatitudeGrid(n_points) {}

 /// The co-latitude grid points in radians.
 const math::Vector<Scalar>& get_colatitudes() const {
   return quadrature_.get_nodes();
 }

  /// The integration weights.
  const math::Vector<Scalar>& get_weights() const { return quadrature_.get_weights(); }

  /// The type of quadrature.
  QuadratureType get_type() { return Quadrature::type; }

 protected:
  Quadrature quadrature_;

};

template <std::floating_point Scalar>
using GaussLegendreGrid = QuadratureLatitudeGrid<GaussLegendreQuadrature<Scalar>, Scalar>;

template <std::floating_point Scalar>
using DoubleGaussGrid = QuadratureLatitudeGrid<DoubleGaussQuadrature<Scalar>, Scalar>;

template <std::floating_point Scalar>
using LobattoGrid = QuadratureLatitudeGrid<LobattoQuadrature<Scalar>, Scalar>;


////////////////////////////////////////////////////////////////////////////////
// Integration functions
////////////////////////////////////////////////////////////////////////////////

static QuadratureProvider<double, FejerQuadrature> quadratures = QuadratureProvider<double, FejerQuadrature>();

template <std::floating_point Scalar>
Scalar integrate_latitudes(math::ConstVectorRef<Scalar> data,
                           const LatitudeGrid<Scalar>& grid) {
  auto weights = grid.get_weights();
  return weights.dot(data);
}

template <std::floating_point Scalar>
Scalar integrate_angles(math::ConstMatrixRef<Scalar> data,
                        math::ConstVectorRef<Scalar> longitudes,
                        const LatitudeGrid<Scalar>& latitude_grid) {
  Scalar result = 0.0;
  math::Index n = longitudes.size();

  Scalar latitude_integral_first = integrate_latitudes<Scalar>(data.row(0), latitude_grid);
  Scalar latitude_integral_left = latitude_integral_first;
  Scalar latitude_integral_right = latitude_integral_first;

  for (math::Index i = 0; i < n - 1; ++i) {
    latitude_integral_right =
        integrate_latitudes<Scalar>(data.row(i + 1), latitude_grid);
    Scalar dl = longitudes[i + 1] - longitudes[i];
    result += 0.5 * (latitude_integral_left + latitude_integral_right) * dl;
    latitude_integral_left = latitude_integral_right;
  }

  Scalar dl = 2.0 * pi_v<Scalar> + longitudes[0] - longitudes[n - 1];
  result += 0.5 * (latitude_integral_first + latitude_integral_right) * dl;

  return result;
}

}  // namespace scattering
