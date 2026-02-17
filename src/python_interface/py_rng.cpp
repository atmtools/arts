#include <arts_omp.h>
#include <nanobind/stl/function.h>
#include <rng.h>

#include "configtypes.h"
#include "hpy_arts.h"
#include "python_interface.h"

namespace Python {
void py_rng(py::module_& m) try {
  using RNG = RandomNumberGenerator<>;
  py::class_<RNG> rngs(m, "RandomNumberGenerator");

  rngs.def(py::init<Time>(),
           "Initialize the RNG with the time seed (at startup by default)",
           "time"_a = Time{})
      .def(py::init<std::mt19937_64::result_type>(),
           "Initialize the RNG with an integer seed",
           "seed"_a);

  rngs.def(
      "uniform_int_distribution",
      [](RNG& r,
         Index lower_bound,
         Index upper_bound) -> std::function<Index()> {
        return r.get<std::uniform_int_distribution>(lower_bound, upper_bound);
      },
      "Generate random numbers following a uniform distribution",
      "lower_bound"_a = 0,
      "upper_bound"_a = 10);

  rngs.def(
      "uniform_real_distribution",
      [](RNG& r,
         Numeric lower_bound,
         Numeric upper_bound) -> std::function<Numeric()> {
        return r.get<std::uniform_real_distribution>(lower_bound, upper_bound);
      },
      "Generate random numbers following a uniform distribution",
      "lower_bound"_a = 0.0,
      "upper_bound"_a = 1.0);

  rngs.def(
      "normal_distribution",
      [](RNG& r, Numeric mean, Numeric stddev) -> std::function<Numeric()> {
        return r.get<std::normal_distribution>(mean, stddev);
      },
      "Generate random numbers following a normal distribution",
      "mean"_a   = 0.0,
      "stddev"_a = 1.0);

  rngs.def(
      "exponential_distribution",
      [](RNG& r, Numeric lambda) -> std::function<Numeric()> {
        return r.get<std::exponential_distribution>(lambda);
      },
      "Generate random numbers following an exponential distribution",
      "lambda"_a = 1.0);

  rngs.def(
      "gamma_distribution",
      [](RNG& r, Numeric alpha, Numeric beta) -> std::function<Numeric()> {
        return r.get<std::gamma_distribution>(alpha, beta);
      },
      "Generate random numbers following a gamma distribution",
      "alpha"_a = 1.0,
      "beta"_a  = 1.0);

  rngs.def(
      "lognormal_distribution",
      [](RNG& r, Numeric mean, Numeric stddev) -> std::function<Numeric()> {
        return r.get<std::lognormal_distribution>(mean, stddev);
      },
      "Generate random numbers following a lognormal distribution",
      "mean"_a   = 0.0,
      "stddev"_a = 1.0);

  rngs.def(
      "chi_squared_distribution",
      [](RNG& r, Numeric n) -> std::function<Numeric()> {
        return r.get<std::chi_squared_distribution>(n);
      },
      "Generate random numbers following a chi-squared distribution",
      "n"_a = 1.0);

  rngs.def(
      "cauchy_distribution",
      [](RNG& r, Numeric location, Numeric scale) -> std::function<Numeric()> {
        return r.get<std::cauchy_distribution>(location, scale);
      },
      "Generate random numbers following a cauchy distribution",
      "location"_a = 0.0,
      "scale"_a    = 1.0);

  rngs.def(
      "weibull_distribution",
      [](RNG& r, Numeric a, Numeric b) -> std::function<Numeric()> {
        return r.get<std::weibull_distribution>(a, b);
      },
      "Generate random numbers following a weibull distribution",
      "a"_a = 1.0,
      "b"_a = 1.0);

  rngs.def(
      "fisher_f_distribution",
      [](RNG& r, Numeric d1, Numeric d2) -> std::function<Numeric()> {
        return r.get<std::fisher_f_distribution>(d1, d2);
      },
      "Generate random numbers following a fisher f distribution",
      "d1"_a = 1.0,
      "d2"_a = 1.0);

  rngs.def(
      "poisson_distribution",
      [](RNG& r, Numeric mean) -> std::function<Index()> {
        return r.get<std::poisson_distribution>(mean);
      },
      "Generate random numbers following a poisson distribution",
      "mean"_a = 1.0);

  rngs.def(
      "geometric_distribution",
      [](RNG& r, Numeric p) -> std::function<Index()> {
        return r.get<std::geometric_distribution>(p);
      },
      "Generate random numbers following a geometric distribution",
      "p"_a = 0.5);

  rngs.def(
      "binomial_distribution",
      [](RNG& r, Index n, Numeric p) -> std::function<Index()> {
        return r.get<std::binomial_distribution>(n, p);
      },
      "Generate random numbers following a binomial distribution",
      "n"_a = 1,
      "p"_a = 0.5);

  rngs.def(
      "student_t_distribution",
      [](RNG& r, Numeric n) -> std::function<Numeric()> {
        return r.get<std::student_t_distribution>(n);
      },
      "Generate random numbers following a student t distribution",
      "n"_a = 1.0);

  rngs.doc() =
      "Random number generator interface. Create an instance of this class and call the distribution methods to get a function that generates random numbers following the specified distribution.";

  generic_interface(rngs);

  auto rng = m.def_submodule("random", "Random Number Interface");

  rng.def(
      "uniform_int_distribution",
      [](Size N, Index lower_bound, Index upper_bound) {
        auto f = RandomNumberGenerator<>{}.get<std::uniform_int_distribution>(
            lower_bound, upper_bound);

        IndexVector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a uniform distribution

Wraps `std::uniform_int_distribution`_ from the C++ standard library.

.. _std::uniform_int_distribution: https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.uniform_int_distribution(10000, 0, 9), bins=10, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
lower_bound : int
    Lower bound of the uniform distribution. Default is 0.
upper_bound : int
    Upper bound of the uniform distribution. Default is 10.

Returns
-------
IndexVector
    The distribution
)",
      "N"_a           = Size{1},
      "lower_bound"_a = 0,
      "upper_bound"_a = 10);

  rng.def(
      "uniform_real_distribution",
      [](Size N, Numeric lower_bound, Numeric upper_bound) {
        auto f = RandomNumberGenerator<>{}.get<std::uniform_real_distribution>(
            lower_bound, upper_bound);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a uniform distribution

Wraps `std::uniform_real_distribution`_ from the C++ standard library.

.. _std::uniform_real_distribution: https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.uniform_real_distribution(10000, 0, 10), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
lower_bound : float
    Lower bound of the uniform distribution. Default is 0.0.
upper_bound : float
    Upper bound of the uniform distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a           = Size{1},
      "lower_bound"_a = 0.0,
      "upper_bound"_a = 1.0);

  rng.def(
      "normal_distribution",
      [](Size N, Numeric mean, Numeric stddev) {
        auto f = RandomNumberGenerator<>{}.get<std::normal_distribution>(
            mean, stddev);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a normal distribution

Wraps `std::normal_distribution`_ from the C++ standard library.

.. _std::normal_distribution: https://en.cppreference.com/w/cpp/numeric/random/normal_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.normal_distribution(10000, 0, 10), bins=51, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
mean : float
    Mean of the normal distribution. Default is 0.0.
stddev : float
    Standard deviation of the normal distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a      = Size{1},
      "mean"_a   = 0.0,
      "stddev"_a = 1.0);

  rng.def(
      "exponential_distribution",
      [](Size N, Numeric lambda) {
        auto f = RandomNumberGenerator<>{}.get<std::exponential_distribution>(
            lambda);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following an exponential distribution

Wraps `std::exponential_distribution`_ from the C++ standard library.

.. _std::exponential_distribution: https://en.cppreference.com/w/cpp/numeric/random/exponential_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.exponential_distribution(10000, 1.0), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
lambda : float
    Rate parameter of the exponential distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a      = Size{1},
      "lambda"_a = 1.0);

  rng.def(
      "gamma_distribution",
      [](Size N, Numeric alpha, Numeric beta) {
        auto f =
            RandomNumberGenerator<>{}.get<std::gamma_distribution>(alpha, beta);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a gamma distribution

Wraps `std::gamma_distribution`_ from the C++ standard library.

.. _std::gamma_distribution: https://en.cppreference.com/w/cpp/numeric/random/gamma_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.gamma_distribution(10000, 1.0, 1.0), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
alpha : float
    Shape parameter of the gamma distribution. Default is 1.0.
beta : float
    Scale parameter of the gamma distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a     = Size{1},
      "alpha"_a = 1.0,
      "beta"_a  = 1.0);

  rng.def(
      "lognormal_distribution",
      [](Size N, Numeric mean, Numeric stddev) {
        auto f = RandomNumberGenerator<>{}.get<std::lognormal_distribution>(
            mean, stddev);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a lognormal distribution

Wraps `std::lognormal_distribution`_ from the C++ standard library.

.. _std::lognormal_distribution: https://en.cppreference.com/w/cpp/numeric/random/lognormal_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.lognormal_distribution(10000, 0.0, 1.0), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
mean : float
    Mean of the underlying normal distribution. Default is 0.0.
stddev : float
    Standard deviation of the underlying normal distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a      = Size{1},
      "mean"_a   = 0.0,
      "stddev"_a = 1.0);

  rng.def(
      "chi_squared_distribution",
      [](Size N, Numeric n) {
        auto f =
            RandomNumberGenerator<>{}.get<std::chi_squared_distribution>(n);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a :math:`\chi^2` distribution

Wraps `std::chi_squared_distribution`_ from the C++ standard library.

.. _std::chi_squared_distribution: https://en.cppreference.com/w/cpp/numeric/random/chi_squared_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.chi_squared_distribution(10000, 1.0), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
n : float
    Degrees of freedom of the chi-squared distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a = Size{1},
      "n"_a = 1.0);

  rng.def(
      "cauchy_distribution",
      [](Size N, Numeric location, Numeric scale) {
        auto f = RandomNumberGenerator<>{}.get<std::cauchy_distribution>(
            location, scale);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a Cauchy/Lorentz distribution

Wraps `std::cauchy_distribution`_ from the C++ standard library.

.. _std::cauchy_distribution: https://en.cppreference.com/w/cpp/numeric/random/cauchy_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.cauchy_distribution(10000, 0.0, 1.0), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
location : float
    Location parameter of the cauchy distribution. Default is 0.0.
scale : float
    Scale parameter of the cauchy distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a        = Size{1},
      "location"_a = 0.0,
      "scale"_a    = 1.0);

  rng.def(
      "weibull_distribution",
      [](Size N, Numeric shape, Numeric scale) {
        auto f = RandomNumberGenerator<>{}.get<std::weibull_distribution>(
            shape, scale);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a Weibull distribution

Wraps `std::weibull_distribution`_ from the C++ standard library.

.. _std::weibull_distribution: https://en.cppreference.com/w/cpp/numeric/random/weibull_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.weibull_distribution(10000, 1.0, 1.0), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
shape : float
    Shape parameter of the weibull distribution. Default is 1.0.
scale : float
    Scale parameter of the weibull distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a     = Size{1},
      "shape"_a = 1.0,
      "scale"_a = 1.0);

  rng.def(
      "fisher_f_distribution",
      [](Size N, Numeric d1, Numeric d2) {
        auto f =
            RandomNumberGenerator<>{}.get<std::fisher_f_distribution>(d1, d2);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a Fisher F distribution

Wraps `std::fisher_f_distribution`_ from the C++ standard library.

.. _std::fisher_f_distribution: https://en.cppreference.com/w/cpp/numeric/random/fisher_f_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa
   import numpy as np

   plt.hist(np.log(pa.arts.random.fisher_f_distribution(10000, 1.0, 1.0)), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Log of Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
d1 : float
    First degrees of freedom of the fisher f distribution. Default is 1.0.
d2 : float
    Second degrees of freedom of the fisher f distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a  = Size{1},
      "d1"_a = 1.0,
      "d2"_a = 1.0);

  rng.def(
      "poisson_distribution",
      [](Size N, Numeric mean) {
        auto f = RandomNumberGenerator<>{}.get<std::poisson_distribution>(mean);

        IndexVector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a Poisson distribution

Wraps `std::poisson_distribution`_ from the C++ standard library.

.. _std::poisson_distribution: https://en.cppreference.com/w/cpp/numeric/random/poisson_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.poisson_distribution(10000, 250), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
mean : float
    Mean of the poisson distribution. Default is 1.0.

Returns
-------
IndexVector
    The distribution
)",
      "N"_a    = Size{1},
      "mean"_a = 1.0);

  rng.def(
      "geometric_distribution",
      [](Size N, Numeric p) {
        auto f = RandomNumberGenerator<>{}.get<std::geometric_distribution>(p);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a geometric distribution

Wraps `std::geometric_distribution`_ from the C++ standard library.

.. _std::geometric_distribution: https://en.cppreference.com/w/cpp/numeric/random/geometric_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.geometric_distribution(10000, 0.01), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
p : float
    Probability of success in the geometric distribution. Default is 0.5.

Returns
-------
Vector
    The distribution
)",
      "N"_a = Size{1},
      "p"_a = 0.5);

  rng.def(
      "binomial_distribution",
      [](Size N, Index n, Numeric p) {
        auto f =
            RandomNumberGenerator<>{}.get<std::binomial_distribution>(n, p);

        IndexVector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a binomial distribution

Wraps `std::binomial_distribution`_ from the C++ standard library.

.. _std::binomial_distribution: https://en.cppreference.com/w/cpp/numeric/random/binomial_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.binomial_distribution(10000, 10, 0.5), bins=11, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
n : int
    Number of trials in the binomial distribution. Default is 1.
p : float
    Probability of success in the binomial distribution. Default is 0.5.

Returns
-------
IndexVector
    The distribution
)",
      "N"_a = Size{1},
      "n"_a = 1,
      "p"_a = 0.5);

  rng.def(
      "student_t_distribution",
      [](Size N, Numeric n) {
        auto f = RandomNumberGenerator<>{}.get<std::student_t_distribution>(n);

        Vector x(N);

#pragma omp parallel for if (arts_omp_parallel(N))
        for (auto& a : x) a = f();

        return x;
      },
      R"(Generate random numbers following a Student :math:`t` distribution

Wraps `std::student_t_distribution`_ from the C++ standard library.

.. _std::student_t_distribution: https://en.cppreference.com/w/cpp/numeric/random/student_t_distribution

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import pyarts3 as pa

   plt.hist(pa.arts.random.student_t_distribution(10000, 1.0), bins=50, density=True)
   plt.title("A plotting example")
   plt.xlabel("Value")
   plt.ylabel("Density")
   plt.yscale("log")

Parameters
----------
N : int
    Number of random numbers to generate. Default is 1.
n : float
    Degrees of freedom of the student t distribution. Default is 1.0.

Returns
-------
Vector
    The distribution
)",
      "N"_a = Size{1},
      "n"_a = 1.0);

} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize RNG\n{}", e.what()));
}
}  // namespace Python
