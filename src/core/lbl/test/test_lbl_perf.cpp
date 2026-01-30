#include <lbl.h>
#include <rng.h>
#include <time_report.h>

#include "configtypes.h"
#include "enumsLineShapeModelType.h"
#include "enumsLineShapeModelVariable.h"
#include "enumsQuantumNumberType.h"

Numeric lbl_temperature_t0(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t0");

  using namespace lbl::temperature;
  // constexpr Numeric T0 = 296;
  // constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::T0(x[i, 0]);

  return result;
}

Numeric lbl_temperature_t1(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t1");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::T1(x[i, 0], x[i, 1], T0, T);

  return result;
}

Numeric lbl_temperature_t2(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t2");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += model::T2(x[i, 0], x[i, 1], x[i, 2], T0, T);

  return result;
}

Numeric lbl_temperature_t3(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t3");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::T3(x[i, 0], x[i, 1], T0, T);
  return result;
}

Numeric lbl_temperature_t4(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t4");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += model::T4(x[i, 0], x[i, 1], x[i, 2], T0, T);
  return result;
}

Numeric lbl_temperature_t5(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t5");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::T5(x[i, 0], x[i, 1], T0, T);
  return result;
}

Numeric lbl_temperature_dpl(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_dpl");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += model::DPL(x[i, 0], x[i, 1], x[i, 2], x[i, 3], T0, T);
  return result;
}

Numeric lbl_temperature_aer(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_aer");

  using namespace lbl::temperature;
  // constexpr Numeric T0 = 296;
  constexpr Numeric T = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += model::AER(x[i, 0], x[i, 1], x[i, 2], x[i, 3], T);
  return result;
}

Numeric lbl_temperature_poly(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_poly");

  using namespace lbl::temperature;
  // constexpr Numeric T0 = 296;
  constexpr Numeric T = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::POLY(x[i], T);
  return result;
}

Numeric lbl_data_line_s(const std::vector<lbl::line>& lines) {
  ARTS_NAMED_TIME_REPORT("lbl_data_line_s");

  const Index n = lines.size();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += lines[i].s(250.0, 182.0);
  return result;
}

int main() {
  {
    constexpr Index M = 10'000'000;
    constexpr Index n = 4;
    const Matrix x    = random_numbers<2>({M, n}, 0.0, 1.0);

    std::println("{}", lbl_temperature_t0(x));
    std::println("{}", lbl_temperature_t1(x));
    std::println("{}", lbl_temperature_t2(x));
    std::println("{}", lbl_temperature_t3(x));
    std::println("{}", lbl_temperature_t4(x));
    std::println("{}", lbl_temperature_t5(x));
    std::println("{}", lbl_temperature_dpl(x));
    std::println("{}", lbl_temperature_aer(x));
    std::println("{}", lbl_temperature_poly(x));
  }

  {
    constexpr Index M = 1'000'000;
    std::vector<lbl::line> lines(M);

    constexpr Index n = 5;
    const Matrix x    = random_numbers<2>({M, n}, 0.0, 1.0);

    for (Index i = 0; i < M; i++) {
      auto& line = lines[i];

      line.a      = 4.479289583303983e-09 * (1 + nonstd::abs(x[i, 0]));
      line.f0     = 118750348044.712 + nonstd::abs(x[i, 1]) * 1e9;
      line.e0     = 1e-23 * (1 + nonstd::abs(x[i, 2]));
      line.gu     = 3.0;
      line.gl     = 1.0;
      line.z.gu() = 1.0011;
      line.z.gl() = 0.0;
      line.ls.T0  = 296.0;
      auto& ls    = line.ls.single_models["O2"_spec];
      ls.data[LineShapeModelVariable::G0] = lbl::temperature::data(
          LineShapeModelType::T1, Vector{x[i, Range(3, 2)]});
      line.qn[QuantumNumberType::J] = {.upper = Rational(1),
                                       .lower = Rational(0)};
    }

    std::println("{}", lbl_data_line_s(lines));
  }

  arts::print_report();
}