#include "extinction_matrix.h"
#include "test_utils.h"

#include <iostream>

bool test_extinction_matrix_tro() {
  using ExtinctionMatrix =
      scattering::ExtinctionMatrixData<Numeric,
                                       scattering::Format::TRO,
                                       scattering::Representation::Gridded>;

  auto t_grid = std::make_shared<Vector>(Vector{210.0, 240.0, 270.0});
  auto f_grid = std::make_shared<Vector>(Vector{1e9, 10e9, 100e9});

  ExtinctionMatrix extinction_matrix{t_grid, f_grid};
  extinction_matrix = random_tensor<Tensor3>(extinction_matrix.shape());

  auto t_grid_new = std::make_shared<Vector>(Vector{210.0});
  auto f_grid_new = std::make_shared<Vector>(Vector{1e9});

  scattering::ScatteringDataGrids grids(t_grid_new, f_grid_new, nullptr);
  auto weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, nullptr, grids);

  auto extinction_matrix_interp = extinction_matrix.regrid(grids, weights);
  Numeric err =
      std::abs(extinction_matrix[0, 0, 0] - extinction_matrix_interp[0, 0, 0]);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);

  fill_along_axis<0>(
      reinterpret_cast<matpack::data_t<Numeric, 3> &>(extinction_matrix));
  (*t_grid_new)[0] = 1.2345;
  (*f_grid_new)[0] = 1.2346;
  grids = scattering::ScatteringDataGrids{t_grid_new, f_grid_new, nullptr};
  weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, nullptr, grids);
  extinction_matrix_interp = extinction_matrix.regrid(grids, weights);
  err = std::abs(extinction_matrix_interp[0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<1>(
      reinterpret_cast<matpack::data_t<Numeric, 3> &>(extinction_matrix));
  extinction_matrix_interp = extinction_matrix.regrid(grids, weights);
  err = std::abs(extinction_matrix_interp[0, 0, 0] - 1.2346);
  if (err > 1e-10) {
    return false;
  }

  return true;
}

bool test_extinction_matrix_aro() {
  using ExtinctionMatrix =
      scattering::ExtinctionMatrixData<Numeric,
                                       scattering::Format::ARO,
                                       scattering::Representation::Gridded>;

  auto t_grid = std::make_shared<Vector>(Vector{210.0, 240.0, 270.0});
  auto f_grid = std::make_shared<Vector>(Vector{1e9, 10e9, 100e9});
  auto za_inc_grid = std::make_shared<Vector>(matpack::uniform_grid(0.0, 10.0, 19));

  ExtinctionMatrix extinction_matrix{t_grid, f_grid, za_inc_grid};
  extinction_matrix = random_tensor<Tensor4>(extinction_matrix.shape());

  auto t_grid_new = std::make_shared<Vector>(Vector{210.0});
  auto f_grid_new = std::make_shared<Vector>(Vector{1e9});
  auto za_inc_grid_new = std::make_shared<Vector>(Vector{0.0});

  scattering::ScatteringDataGrids grids(
      t_grid_new, f_grid_new, za_inc_grid_new, nullptr, nullptr);
  auto weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, za_inc_grid, nullptr, nullptr, grids);

  auto extinction_matrix_interp = extinction_matrix.regrid(grids, weights);
  Numeric err = std::abs(extinction_matrix[0, 0, 0, 0] -
                         extinction_matrix_interp[0, 0, 0, 0]);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);
  fill_along_axis<0>(*za_inc_grid);

  fill_along_axis<0>(
      reinterpret_cast<matpack::data_t<Numeric, 4> &>(extinction_matrix));
  (*t_grid_new)[0] = 1.2345;
  (*f_grid_new)[0] = 1.2346;
  (*za_inc_grid_new)[0] = 1.2347;
  grids = scattering::ScatteringDataGrids{
      t_grid_new, f_grid_new, za_inc_grid_new, nullptr, nullptr};
  weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, za_inc_grid, nullptr, nullptr, grids);
  extinction_matrix_interp = extinction_matrix.regrid(grids, weights);
  err = std::abs(extinction_matrix_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<1>(
      reinterpret_cast<matpack::data_t<Numeric, 4> &>(extinction_matrix));
  extinction_matrix_interp = extinction_matrix.regrid(grids, weights);
  err = std::abs(extinction_matrix_interp[0, 0, 0, 0] - 1.2346);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<2>(
      reinterpret_cast<matpack::data_t<Numeric, 4> &>(extinction_matrix));
  extinction_matrix_interp = extinction_matrix.regrid(grids, weights);
  err = std::abs(extinction_matrix_interp[0, 0, 0, 0] - 1.2347);
  if (err > 1e-10) {
    return false;
  }

  return true;
}

int main() {
  bool passed = false;

  std::cout << "Testing extinction matrix (TRO): ";
  passed = test_extinction_matrix_tro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  passed = false;
  std::cout << "Testing extinction matrix (ARO): ";
  passed = test_extinction_matrix_aro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  return 0;
}
