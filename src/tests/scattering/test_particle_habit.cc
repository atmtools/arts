#include <iostream>

#include "particle_habit.h"
#include "test_utils.h"
#include "xml_io.h"

bool test_particle_habit_from_legacy_tro() {

  std::string meta_path = std::string(TEST_DATA_PATH) + std::string("arts-xml-data/scattering/H2O_ice/ScatteringMetaFile_allH2Oice.xml");
  std::string data_path = std::string(TEST_DATA_PATH) + std::string("arts-xml-data/scattering/H2O_ice/SingleScatteringFile_allH2Oice.xml");

  ArrayOfSingleScatteringData legacy_data;
  xml_read_from_file(data_path, legacy_data);
  ArrayOfScatteringMetaData legacy_meta;
  xml_read_from_file(meta_path, legacy_meta);

  auto habit = scattering::ParticleHabit::from_legacy_tro(legacy_data, legacy_meta);
  using SSD = scattering::SingleScatteringData<Numeric, scattering::Format::TRO, scattering::Representation::Gridded>;

  for (Index ind = 0; ind < habit.size(); ++ind) {

    Numeric err = max_error<Tensor3View>((*std::get<SSD>(habit[ind]).phase_matrix)(joker, 0, joker, joker),
                                         legacy_data[ind].pha_mat_data(0, joker, joker, 0, 0, 0, joker));
    if (err > 0) {
      return false;
    }

    err = max_error<MatrixView>(std::get<SSD>(habit[ind]).extinction_matrix(0, joker, joker),
                                legacy_data[ind].ext_mat_data(joker, 0, 0, 0, joker));
    if (err > 0) {
      return false;
    }

    err = max_error<MatrixView>(std::get<SSD>(habit[ind]).absorption_vector(joker, 0, joker),
                                legacy_data[ind].abs_vec_data(0, joker, 0, 0, joker));
    if (err > 0) {
      return false;
    }
  }
  return true;
}

int main() {
  bool passed = false;
  std::cout << "Test conversion from legacy format (TRO): ";
  passed = test_particle_habit_from_legacy_tro();
  if (passed) {
    std::cout << "PASSED." << std::endl;
  } else {
    std::cout << "FAILED." << std::endl;
    return 1;
  }
  return 0;
}
