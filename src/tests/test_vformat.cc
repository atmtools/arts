#include <workspace.h>

#include "quantum.h"

QuantumIdentifier get_quantum_identifier() {
  QuantumIdentifier qid{};
  qid.isot = "H2O-161"_isot;
  qid.state.emplace(QuantumNumberType::J,
                    Quantum::UpperLower{.upper = Quantum::Value{1},
                                        .lower = Quantum::Value{1}});
  qid.state.emplace(QuantumNumberType::K,
                    Quantum::UpperLower{.upper = Quantum::Value{0},
                                        .lower = Quantum::Value{0}});
  return qid;
}

lbl::line get_lbl_line() {
  lbl::line l{};
  l.f0 = 1.0;
  l.a  = 2.0;
  l.e0 = 3.0;
  l.gu = 4.0;
  l.gl = 5.0;
  l.qn.emplace(QuantumNumberType::J,
               Quantum::UpperLower{.upper = Quantum::Value{1},
                                   .lower = Quantum::Value{1}});

  l.ls.one_by_one = false;
  l.ls.T0         = 300.0;
  auto& mod       = l.ls.single_models["H2O"_spec];

  mod.data[LineShapeModelVariable::G0] =
      lbl::temperature::data{LineShapeModelType::T0, Vector{6.0}};
  mod.data[LineShapeModelVariable::D0] =
      lbl::temperature::data{LineShapeModelType::T1, Vector{7.0, 8.0}};

  l.z.gu(6.0);
  l.z.gl(7.0);
  l.z.on = true;

  return l;
}

lbl::band_data get_lbl_band_data() {
  lbl::band_data b{};

  b.lines.emplace_back(get_lbl_line());

  return b;
}

AbsorptionBands get_absorption_bands() {
  return {{get_quantum_identifier(), get_lbl_band_data()}};
}

#define TEST_MACRO(x)                         \
  void test_##x() {                           \
    auto data = get_##x();                    \
                                              \
    std::println("NORMAL    {}", data);       \
    std::println("FMT       {:B,qN}", data);  \
    std::println("FMT-SHORT {:B,qNs}", data); \
  }

TEST_MACRO(lbl_line)
TEST_MACRO(lbl_band_data)
TEST_MACRO(quantum_identifier)
TEST_MACRO(absorption_bands)

#undef TEST_MACRO

#define TEST_MACRO(x)                                                       \
  std::println("-----------------------\n{}\n-----------------------", #x); \
  test_##x();

int main() {
  TEST_MACRO(lbl_line)
  TEST_MACRO(lbl_band_data)
  TEST_MACRO(quantum_identifier)
  TEST_MACRO(absorption_bands)
  std::println("-----------------------");
}
