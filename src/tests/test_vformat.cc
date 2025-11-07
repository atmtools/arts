#include <workspace.h>

#include "quantum.h"

QuantumIdentifier get_quantum_identifier() {
  QuantumIdentifier qid{};
  qid.isot = "H2O-161"_isot;
  qid.state.emplace(QuantumNumberType::J,
                    Quantum::UpperLower{.upper = 1, .lower = 1});
  qid.state.emplace(QuantumNumberType::K,
                    Quantum::UpperLower{.upper = 0, .lower = 0});
  return qid;
}

AbsorptionLine get_lbl_line() {
  using enum LineShapeModelVariable;
  using enum LineShapeModelType;
  using enum QuantumNumberType;

  AbsorptionLine line{};
  line.f0 = 1.0;
  line.a  = 2.0;
  line.e0 = 3.0;
  line.gu = 4.0;
  line.gl = 5.0;
  line.qn.emplace(J, Quantum::UpperLower{.upper = 1, .lower = 1});

  line.ls.one_by_one = false;
  line.ls.T0         = 300.0;
  auto& mod          = line.ls.single_models["H2O"_spec];
  mod.data[G0]       = {T0, Vector{6.0}};
  mod.data[D0]       = {T1, Vector{7.0, 8.0}};
  line.z.gu(6.0);
  line.z.gl(7.0);
  line.z.on = true;

  return line;
}

AbsorptionBand get_lbl_band_data() {
  AbsorptionBand b{};

  b.lines.emplace_back(get_lbl_line());

  return b;
}

AbsorptionBands get_abs_bands() {
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
TEST_MACRO(abs_bands)

#undef TEST_MACRO

#define TEST_MACRO(x)                                                       \
  std::println("-----------------------\n{}\n-----------------------", #x); \
  test_##x();

int main() {
  TEST_MACRO(lbl_line)
  TEST_MACRO(lbl_band_data)
  TEST_MACRO(quantum_identifier)
  TEST_MACRO(abs_bands)
  std::println("-----------------------");
}
