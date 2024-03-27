#include <iostream>

#include "single_scattering_data.h"
#include "test_utils.h"
#include "xml_io.h"

bool test_single_scattering_data_from_legacy_tro() {
  SingleScatteringData legacy;
  std::string path =
      "/home/simonpf/src/arts/controlfiles/testdata/"
      "scatData/MieAtmlab_Liquid_72.9um.xml.gz";
  xml_read_from_file(path, legacy);

  using SingleScatteringData =
      scattering::SingleScatteringData<Numeric,
                                       scattering::Format::TRO,
                                       scattering::Representation::Gridded,
                                       4>;
  auto ssd = SingleScatteringData::from_legacy_tro(legacy);

  Numeric err = max_error<Tensor3View>(
      (*ssd.phase_matrix)(joker, 0, joker, joker),
      legacy.pha_mat_data(0, joker, joker, 0, 0, 0, joker));
  if (err > 0) {
    return false;
  }

  err = max_error<MatrixView>(ssd.extinction_matrix(0, joker, joker),
                              legacy.ext_mat_data(joker, 0, 0, 0, joker));
  if (err > 0) {
    return false;
  }

  err = max_error<MatrixView>(ssd.absorption_vector(joker, 0, joker),
                              legacy.abs_vec_data(0, joker, 0, 0, joker));
  if (err > 0) {
    return false;
  }
  return true;
}

int main() {
  bool passed = false;
  std::cout << "Test conversion from legacy format (TRO): ";
  passed = test_single_scattering_data_from_legacy_tro();
  if (passed) {
    std::cout << "PASSED." << std::endl;
  } else {
    std::cout << "FAILED." << std::endl;
    return 1;
  }
  return 0;
}
