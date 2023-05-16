#pragma once

#include "artstime.h"
#include "debug.h"
#include "matpack_concepts.h"

#include <algorithm>
#include <limits>
#include <mutex>
#include <ostream>
#include <random>
#include <vector>

/** A C++ standards dependent random number generator class
 *
 * This class once created provides the get<>() method to return a
 * haystack of random samples following a custom distribution model
 * 
 * @tparam Generator A C++ standards uniform random bit generator (defaults to mt19937_64)
 */
template <std::uniform_random_bit_generator Generator = std::mt19937_64>
class RandomNumberGenerator {
  using int_t = typename Generator::result_type;

  static constexpr int_t min_v = std::numeric_limits<int_t>::lowest();
  static constexpr int_t max_v = std::numeric_limits<int_t>::max();

  // A mutable generator that
  mutable Generator v;

  // The
  mutable std::mutex g;

  // Keep the seed so we can "copy" the RandomNumberGenerator
  int_t current_seed;

public:
  /** Construct a new Random Number Generator object
   *
   * @param t Time, defaults to the time when the object was initiated
   */
  RandomNumberGenerator(const Time &t = {}) { seed(t); }

  /** Construct a new Random Number Generator object
   *
   * @param s A generator natural integer seed
   */
  RandomNumberGenerator(const int_t &s) { seed(s); }

  /**
   * @brief Copy a Random Number Generator object
   *
   * @param rng Old object, the state is preserved
   */
  RandomNumberGenerator(const RandomNumberGenerator &rng) {
    std::lock_guard lock{rng.g};
    seed(rng.current_seed);
  }

  /** Returns a new Generator 
   *
   * The new generator is seeded by a uniform_int_distribution created
   * from the current generator
   */
  Generator new_generator() const {
    const std::lock_guard my_lock(this->g);
    return Generator{std::uniform_int_distribution<int_t>(min_v, max_v)(v)};
  }

  /** Returns a random number generator of some random distribution
   *
   * The return type is a function-like object so that calling it generates
   * a random number
   *
   * Note that the random distribution you choose must be constructible
   * and have its output type determined by the input.  One caveat that
   * is very important to be aware of is that if you input Index or int
   * types instead of Numeric when calling this function, the output is
   * not necessarily good, and it can often lead to undefined behavior
   *
   * Some examples of random distribution and call options to this
   * function are:
   *
   * Uniform distribution between [0.0, 1.0):
   * auto gen = rng.get<std::uniform_real_distribution>(0.0, 1.0);
   *
   * Uniform distribution between [-180., 180.):
   * auto gen = rng.get<std::uniform_real_distribution>(-180., 180.);
   *
   * Normal distribution around 0.0 with a standard deviation of 1.0:
   * auto gen = rng.get<std::normal_distribution>(0.0, 1.0);
   *
   * Logormal distribution about 0.0 with a standard deviation of 1.0:
   * auto gen = rng.get<std::lognormal_distribution>(0.0, 1.0);
   *
   * @tparam random_distribution A random number distribution
   * @param x The parameters to construct the random_distribution
   * @return auto A callable generator of random numbers
   */
  template <template <typename>
            class random_distribution = std::uniform_real_distribution,
            typename... Ts>
  auto get(Ts &&...x) const {
    return [v = new_generator(),
            draw = random_distribution(std::forward<Ts>(x)...)]() mutable {
      return draw(v);
    };
  }

  /** Seeds the random number engine with a new and unique seed
   *
   * If the seed is not used, it is safely added to the list of used seeds
   *
   * Otherwise the a selection process is used to determine a new seed
   *
   * WARNING: This function should generally only be called internally,
   * any manual seeding and reseeding messes with the randomness of the
   * distribution, and so it is very, very important you know what you
   * are doing if you call this method
   *
   * @param n The new seed
   * @return true if the seed was used as is, otherwise false
   */
  bool seed(int_t n) {
    static std::mutex seed_mtx;
    static std::vector<int_t> seeds{};

    std::lock_guard class_lock{g};
    std::lock_guard method_lock{seed_mtx};
    bool seed_as_is = true;

    constexpr auto add_one = [](int_t i) {
      return (i >= max_v - 1) ? min_v : i + 1;
    };

    if (std::find(seeds.begin(), seeds.end(), n) != seeds.end()) {
      seed_as_is = false;

      auto minmax = std::minmax_element(seeds.begin(), seeds.end());
      int_t min = *minmax.first;
      int_t max = *minmax.second;

      if (min not_eq min_v) {
        n = add_one(max);
      } else if (static_cast<std::size_t>(max - min) <= seeds.size()) {
        n = add_one(max);
        if (n == min)
          seeds.clear();
      } else {
        do {
          n = add_one(n);
        } while (std::find(seeds.begin(), seeds.end(), n) != seeds.end());
      }
    }
    seeds.push_back(n);

    force_seed(n);

    return seed_as_is;
  }

  /** As the main seed, but the number is extracted from the time
   *
   * The integer type is taken from the epoch counter modulos the
   * maximum possible value (limit ourselves to only positive seeds)
   *
   * @param t The time, defaults to current time
   * @return See main seed
   */
  bool seed(const Time &t = {}) {
    return seed(static_cast<int_t>(t.EpochTime().count() % max_v));
  }

  /** Force the seed to some integer, determining the randomness of the
   * distributions
   *
   * This is cannot be called directly from a parallel environment
   *
   * @param n The new seed
   */
  void force_seed(int_t n) {
    v.seed(n);
    current_seed = n;
  }

  /** output the current seed */
  friend std::ostream &operator<<(std::ostream &os,
                                  const RandomNumberGenerator &rng) {
    return os << rng.current_seed;
  }
};

/** Wraps the generation of a random number generator for a matpack type
 *
 * @tparam random_distribution As for RandomNumberGenerator::get
 * @tparam Generator As for RandomNumberGenerator
 * @param sz Size of output matpack data
 * @param x The parameters passed to get_par<>(x...)
 * @return matpack_data A new object of random numbers of size sz
 */
template <std::size_t N,
          template <typename>
          class random_distribution = std::uniform_real_distribution,
          std::uniform_random_bit_generator Generator = std::mt19937_64,
          typename... param_types, class RNG = RandomNumberGenerator<Generator>,
          class T = decltype(RNG{}.template get<random_distribution>(
              param_types{}...)())>
matpack::matpack_data<T, N> random_numbers(std::array<Index, N> sz,
                                           param_types &&...x) {
  matpack::matpack_data<T, N> out(sz);
  std::generate(out.elem_begin(), out.elem_end(),
                RNG{}.template get<random_distribution>(
                    std::forward<param_types>(x)...));
  return out;
}

/** Wraps more advanced random_numbers generator above for pure vectors
 *
 * See random_numbers above for more details
 */
template <template <typename>
          class random_distribution = std::uniform_real_distribution,
          std::uniform_random_bit_generator Generator = std::mt19937_64,
          typename... param_types>
auto random_numbers(Index n, param_types &&...x) {
  return random_numbers<1, random_distribution, Generator>(
      std::array<Index, 1>{n}, std::forward<param_types>(x)...);
}
