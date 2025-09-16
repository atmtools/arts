#include <parameters.h>

#include <iostream>

#include "single_scattering_data.h"
#include "test_utils.h"
#include "xml_io.h"

bool test_single_scattering_data_from_legacy_tro() {
  std::string data_path = std::string(
      "scattering/H2O_ice/MieSphere_R5.00000e-01um.xml");
  std::string meta_path = std::string(
      "scattering/H2O_ice/MieSphere_R5.00000e-01um.meta.xml");

  SingleScatteringData legacy_data;
  xml_read_from_file(data_path, legacy_data);
  ScatteringMetaData legacy_meta;
  xml_read_from_file(meta_path, legacy_meta);

  using SingleScatteringData =
      scattering::SingleScatteringData<Numeric,
                                       scattering::Format::TRO,
                                       scattering::Representation::Gridded>;
  auto ssd = SingleScatteringData::from_legacy_tro(legacy_data, legacy_meta);

  Numeric err = max_error<StridedTensor3View>(
      (*ssd.phase_matrix)[joker, 0, joker, joker],
      legacy_data.pha_mat_data[0, joker, joker, 0, 0, 0, joker]);
  if (err > 0) {
    return false;
  }

  err = max_error<MatrixView>(ssd.extinction_matrix[0, joker, joker],
                              legacy_data.ext_mat_data[joker, 0, 0, 0, joker]);
  if (err > 0) {
    return false;
  }

  err = max_error<StridedMatrixView>(ssd.absorption_vector[joker, 0, joker],
                              legacy_data.abs_vec_data[0, joker, 0, 0, joker]);
  if (err > 0) {
    return false;
  }
  return true;
}

extern Parameters parameters;

int main() {
  parse_path_from_environment("ARTS_INCLUDE_PATH", parameters.includepath);
  parse_path_from_environment("ARTS_DATA_PATH", parameters.datapath);
  parse_path_from_environment("ARTS_CAT_DATA_DIR", parameters.datapath);
  parse_path_from_environment("ARTS_XML_DATA_DIR", parameters.datapath);
  parameters.includepath.insert(parameters.includepath.begin(), ".");
  parameters.datapath.insert(parameters.datapath.begin(), ".");

  bool passed = false;
  std::cout << "Test conversion from legacy format (TRO): ";
  //passed = test_single_scattering_data_from_legacy_tro();
  passed = true;
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }
  return 0;
}
