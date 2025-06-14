#include <workspace.h>

QuantumIdentifier get_quantum_identifier() {
  QuantumIdentifier qid{};
  qid.isotopologue_index = "H2O-161"_isot_index;
  qid.val.add(QuantumNumberValue{QuantumNumberType::J, 1, 1});
  qid.val.add(QuantumNumberValue{QuantumNumberType::K, 0, 0});
  return qid;
}

lbl::line get_lbl_line() {
  lbl::line l{};
  l.f0 = 1.0;
  l.a  = 2.0;
  l.e0 = 3.0;
  l.gu = 4.0;
  l.gl = 5.0;
  l.qn.val.add(QuantumNumberValue{QuantumNumberType::J, 1, 1});

  l.ls.one_by_one = false;
  l.ls.T0         = 300.0;
  auto& mod       = l.ls.single_models.emplace_back("H2O"_spec);

  mod.data.emplace_back(
      LineShapeModelVariable::G0,
      lbl::temperature::data{LineShapeModelType::T0, Vector{6.0}});
  mod.data.emplace_back(
      LineShapeModelVariable::D0,
      lbl::temperature::data{LineShapeModelType::T1, Vector{7.0, 8.0}});

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
