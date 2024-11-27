#include "absorption_vector.h"
#include "test_utils.h"

#include <iostream>

bool test_absorption_vector_tro() {
  using AbsorptionVector =
      scattering::AbsorptionVectorData<Numeric,
                                       scattering::Format::TRO,
                                       scattering::Representation::Gridded>;

  auto t_grid = std::make_shared<Vector>(Vector{210.0, 240.0, 270.0});
  auto f_grid = std::make_shared<Vector>(Vector{1e9, 10e9, 100e9});

  AbsorptionVector absorption_vector{t_grid, f_grid};
  absorption_vector = random_tensor<Tensor3>(absorption_vector.shape());

  auto t_grid_new = std::make_shared<Vector>(Vector{210.0});
  auto f_grid_new = std::make_shared<Vector>(Vector{1e9});

  scattering::ScatteringDataGrids grids(t_grid_new, f_grid_new, nullptr);
  auto weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, nullptr, grids);

  auto absorption_vector_interp = absorption_vector.regrid(grids, weights);
  Numeric err =
      std::abs(absorption_vector(0, 0, 0) - absorption_vector_interp(0, 0, 0));
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);

  fill_along_axis<0>(
      reinterpret_cast<matpack::matpack_data<Numeric, 3> &>(absorption_vector));
  (*t_grid_new)[0] = 1.2345;
  (*f_grid_new)[0] = 1.2346;
  grids = scattering::ScatteringDataGrids{t_grid_new, f_grid_new, nullptr};
  weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, nullptr, grids);
  absorption_vector_interp = absorption_vector.regrid(grids, weights);
  err = std::abs(absorption_vector_interp(0, 0, 0) - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<1>(
      reinterpret_cast<matpack::matpack_data<Numeric, 3> &>(absorption_vector));
  absorption_vector_interp = absorption_vector.regrid(grids, weights);
  err = std::abs(absorption_vector_interp(0, 0, 0) - 1.2346);
  if (err > 1e-10) {
    return false;
  }

  return true;
}

bool test_absorption_vector_aro() {
  using AbsorptionVector =
      scattering::AbsorptionVectorData<Numeric,
                                       scattering::Format::ARO,
                                       scattering::Representation::Gridded>;

  auto t_grid = std::make_shared<Vector>(Vector{210.0, 240.0, 270.0});
  auto f_grid = std::make_shared<Vector>(Vector{1e9, 10e9, 100e9});
  auto za_inc_grid = std::make_shared<Vector>(uniform_grid(0.0, 10.0, 19));

  AbsorptionVector absorption_vector{t_grid, f_grid, za_inc_grid};
  absorption_vector = random_tensor<Tensor4>(absorption_vector.shape());

  auto t_grid_new = std::make_shared<Vector>(Vector{210.0});
  auto f_grid_new = std::make_shared<Vector>(Vector{1e9});
  auto za_inc_grid_new = std::make_shared<Vector>(Vector{0.0});

  scattering::ScatteringDataGrids grids(
      t_grid_new, f_grid_new, za_inc_grid_new, nullptr, nullptr);
  auto weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, za_inc_grid, nullptr, nullptr, grids);

  auto absorption_vector_interp = absorption_vector.regrid(grids, weights);
  Numeric err = std::abs(absorption_vector(0, 0, 0, 0) -
                         absorption_vector_interp(0, 0, 0, 0));
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);
  fill_along_axis<0>(*za_inc_grid);

  fill_along_axis<0>(
      reinterpret_cast<matpack::matpack_data<Numeric, 4> &>(absorption_vector));
  (*t_grid_new)[0] = 1.2345;
  (*f_grid_new)[0] = 1.2346;
  (*za_inc_grid_new)[0] = 1.2347;
  grids = scattering::ScatteringDataGrids{
      t_grid_new, f_grid_new, za_inc_grid_new, nullptr, nullptr};
  weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, za_inc_grid, nullptr, nullptr, grids);
  absorption_vector_interp = absorption_vector.regrid(grids, weights);
  err = std::abs(absorption_vector_interp(0, 0, 0, 0) - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<1>(
      reinterpret_cast<matpack::matpack_data<Numeric, 4> &>(absorption_vector));
  absorption_vector_interp = absorption_vector.regrid(grids, weights);
  err = std::abs(absorption_vector_interp(0, 0, 0, 0) - 1.2346);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<2>(
      reinterpret_cast<matpack::matpack_data<Numeric, 4> &>(absorption_vector));
  absorption_vector_interp = absorption_vector.regrid(grids, weights);
  err = std::abs(absorption_vector_interp(0, 0, 0, 0) - 1.2347);
  if (err > 1e-10) {
    return false;
  }

  return true;
}

int main() {
  bool passed = false;
  std::cout << "Testing absorption vector (TRO): ";
  passed = test_absorption_vector_tro();
  if (passed) {
    std::cout << "PASSED." << std::endl;
  } else {
    std::cout << "FAILED." << std::endl;
    return 1;
  }

  std::cout << "Testing absorption vector (ARO): ";
  passed = test_absorption_vector_aro();
  if (passed) {
    std::cout << "PASSED." << std::endl;
  } else {
    std::cout << "FAILED." << std::endl;
    return 1;
  }

  return 0;
}
