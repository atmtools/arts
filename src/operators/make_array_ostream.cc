#include <iostream>

#define MACROEXECUTE(Type)                                                      \
  std::cout                                                                     \
      << "std::ostream& operator<<(std::ostream&, const Array< "#Type" >&);\n";      \
  std::cerr                                                                     \
      << "static_assert(coutable<"#Type">, \"" #Type " NOT PRINTABLE TO STREAM\");\n"; \
  std::cerr                                                                     \
      << "std::ostream& operator<<(std::ostream& os, const Array<"#Type">& x) {\n";     \
  std::cerr << "  for (auto& a : x) os << a << \'\\n\';\n"                      \
               "  return os;\n}\n\n";

int main() {
  std::cout << "#pragma once\n\n#include <ostream>\n\n#include <auto_wsg.h>\n\n";
  std::cerr << "#include <array_ostream.h>\n\ntemplate <typename T>\n"
               "concept coutable = requires(std::ostream& os, const T& t) {\n"
               "  os << t;\n};\n\n";

  std::cout << "namespace std {\n";
  std::cerr << "namespace std {\n";
  MACROEXECUTE(Index)
  MACROEXECUTE(ArrayOfIndex)
  MACROEXECUTE(String)
  MACROEXECUTE(ArrayOfString)
  std::cout << "}  // namespace std\n";
  std::cerr << "}  // namespace std\n";
  
  MACROEXECUTE(CIARecord)
  MACROEXECUTE(XsecRecord)
  MACROEXECUTE(SingleScatteringData)
  MACROEXECUTE(ArrayOfSingleScatteringData)
  MACROEXECUTE(ScatteringMetaData)
  MACROEXECUTE(ArrayOfScatteringMetaData)
  MACROEXECUTE(Agenda)
  MACROEXECUTE(GriddedField1)
  MACROEXECUTE(GriddedField2)
  MACROEXECUTE(GriddedField3)
  MACROEXECUTE(GriddedField4)
  MACROEXECUTE(GriddedField5)
  MACROEXECUTE(ArrayOfGriddedField1)
  MACROEXECUTE(ArrayOfGriddedField2)
  MACROEXECUTE(ArrayOfGriddedField3)
  MACROEXECUTE(ArrayOfGriddedField4)
  MACROEXECUTE(ArrayOfGriddedField5)
  MACROEXECUTE(Time)
  MACROEXECUTE(ArrayOfTime)
  MACROEXECUTE(Ppath)
  MACROEXECUTE(RetrievalQuantity)
  MACROEXECUTE(Sparse)
  MACROEXECUTE(Sun)
  MACROEXECUTE(TelsemAtlas)

  MACROEXECUTE(ArrayOfSpeciesTag)

  std::cout << "namespace Quantum::Number {\n";
  std::cerr << "namespace Quantum::Number {\n";
  MACROEXECUTE(GlobalState)
  std::cout << "}  // namespace Quantum::Number\n";
  std::cerr << "}  // namespace Quantum::Number\n";
  
  std::cout << "namespace Atm {\n";
  std::cerr << "namespace Atm {\n";
  MACROEXECUTE(Point)
  std::cout << "}  // namespace Atm\n";
  std::cerr << "}  // namespace Atm\n";
  
  std::cout << "namespace Absorption {\n";
  std::cerr << "namespace Absorption {\n";
  MACROEXECUTE(AbsorptionLines)
  MACROEXECUTE(ArrayOfAbsorptionLines)
  std::cout << "}  // namespace Absorption\n";
  std::cerr << "}  // namespace Absorption\n";

  std::cout << "namespace matpack {\n";
  std::cerr << "namespace matpack {\n";
  MACROEXECUTE(Vector)
  MACROEXECUTE(Matrix)
  MACROEXECUTE(Tensor3)
  MACROEXECUTE(Tensor4)
  MACROEXECUTE(Tensor5)
  MACROEXECUTE(Tensor6)
  MACROEXECUTE(Tensor7)
  MACROEXECUTE(ArrayOfVector)
  MACROEXECUTE(ArrayOfMatrix)
  MACROEXECUTE(ArrayOfTensor3)
  MACROEXECUTE(ArrayOfTensor4)
  MACROEXECUTE(ArrayOfTensor5)
  MACROEXECUTE(ArrayOfTensor6)
  MACROEXECUTE(ArrayOfTensor7)
  MACROEXECUTE(MuelmatVector)
  MACROEXECUTE(MuelmatMatrix)
  MACROEXECUTE(PropmatVector)
  MACROEXECUTE(PropmatMatrix)
  MACROEXECUTE(StokvecVector)
  MACROEXECUTE(StokvecMatrix)
  MACROEXECUTE(ArrayOfMuelmatVector)
  MACROEXECUTE(ArrayOfMuelmatMatrix)
  MACROEXECUTE(ArrayOfPropmatVector)
  MACROEXECUTE(ArrayOfPropmatMatrix)
  MACROEXECUTE(ArrayOfStokvecVector)
  MACROEXECUTE(ArrayOfStokvecMatrix)
  std::cout << "}  // namespace matpack\n";
  std::cerr << "}  // namespace matpack\n";

  std::cout << "namespace rtepack {\n";
  std::cerr << "namespace rtepack {\n";
  MACROEXECUTE(Muelmat)
  MACROEXECUTE(Propmat)
  MACROEXECUTE(Stokvec)
  std::cout << "}  // namespace rtepack\n";
  std::cerr << "}  // namespace rtepack\n";
}
