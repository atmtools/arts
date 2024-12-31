#include <iostream>

#include <test_utils.h>
#include <scattering/sht.h>
#include <scattering/integration.h>
#include <scattering/xml_io_scattering.h>


using namespace scattering;

bool test_serialize_sht() {

  auto sht = sht::SHT(5, 5, 32, 32);

  std::ostringstream ostrm;
  xml_write_to_stream(ostrm, sht, nullptr, "unused");

  std::istringstream istrm{ostrm.str()};
  auto sht_other = sht::SHT(6, 6, 48, 48);
  xml_read_from_stream(istrm, sht_other, nullptr);


  if (sht_other.get_l_max() != sht.get_l_max()) {
    return false;
  }
  if (sht_other.get_m_max() != sht.get_m_max()) {
    return false;
  }
  if (sht_other.get_n_azimuth_angles() != sht.get_n_azimuth_angles()) {
    return false;
  }
  if (sht_other.get_n_zenith_angles() != sht.get_n_zenith_angles()) {
    return false;
  }
  return true;
}

bool test_serialize_irregular_zenith_angle_grid() {

  scattering::GaussLegendreQuadrature quad(3);
  scattering::IrregularZenithAngleGrid grid = quad.get_weights();

  std::ostringstream ostrm;
  xml_write_to_stream(ostrm, grid, nullptr, "unused");

  std::istringstream istrm{ostrm.str()};
  scattering::IrregularZenithAngleGrid new_grid = quad.get_weights();
  xml_read_from_stream(istrm, new_grid, nullptr);

  if (grid.size() != new_grid.size()) return false;

  auto error = max_error(grid, new_grid);
  if (error > 1e-6) return false;

  return true;
}

template <typename Quadrature>
bool test_serialize_quadrature_grid() {

  Quadrature grid(10);

  std::ostringstream ostrm;
  xml_write_to_stream(ostrm, grid, nullptr, "unused");

  std::istringstream istrm{ostrm.str()};
  Quadrature new_grid;
  xml_read_from_stream(istrm, new_grid, nullptr);

  if (grid.size() != new_grid.size()) return false;

  auto error = max_error(grid, new_grid);
  if (error > 1e-6) return false;

  return true;
}


bool test_serialize_zenith_angle_grid() {

  scattering::ZenithAngleGrid grid = scattering::GaussLegendreGrid(10);

  std::ostringstream ostrm;
  xml_write_to_stream(ostrm, grid, nullptr, "unused");
  std::istringstream istrm{ostrm.str()};
  scattering::ZenithAngleGrid new_grid{};
  xml_read_from_stream(istrm, new_grid, nullptr);

  if (grid_size(grid) != grid_size(new_grid)) return false;
  auto error = max_error(grid_vector(grid), grid_vector(new_grid));
  if (error > 1e-6) return false;

  return true;
}


int main() {

  bool passed = false;

#ifndef ARTS_NO_SHTNS
  std::cout << "Testing SHT serialization: ";
  passed = test_serialize_sht();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }
#endif

  std::cout << "Testing IrregularZenithAngleGridSerialization: ";
  passed = test_serialize_irregular_zenith_angle_grid();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing GaussLegendreGrid: ";
  passed = test_serialize_quadrature_grid<scattering::GaussLegendreGrid>();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing DoubleGaussGrid: ";
  passed = test_serialize_quadrature_grid<scattering::DoubleGaussGrid>();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing LobattoGrid: ";
  passed = test_serialize_quadrature_grid<scattering::LobattoGrid>();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing FejerGrid: ";
  passed = test_serialize_quadrature_grid<scattering::FejerGrid>();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing generic zenith angle grid: ";
  passed = test_serialize_zenith_angle_grid();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  return 0;
}
