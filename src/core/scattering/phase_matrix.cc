#include "phase_matrix.h"

namespace scattering {
Matrix expand_phase_matrix(const StridedConstVectorView &compact) {
  Matrix mat{4, 4};
  mat[0, 0] = detail::f11(compact);
  mat[0, 1] = detail::f12(compact);
  mat[1, 0] = detail::f12(compact);
  mat[1, 1] = detail::f22(compact);
  mat[2, 2] = detail::f33(compact);
  mat[2, 3] = detail::f34(compact);
  mat[3, 2] = detail::f34(compact);
  mat[3, 3] = detail::f33(compact);
  return mat;
}

ComplexMatrix expand_phase_matrix(
    const StridedConstComplexVectorView &compact) {
  ComplexMatrix mat{4, 4};
  mat[0, 0] = detail::f11(compact);
  mat[0, 1] = detail::f12(compact);
  mat[1, 0] = detail::f12(compact);
  mat[1, 1] = detail::f22(compact);
  mat[2, 2] = detail::f33(compact);
  mat[2, 3] = detail::f34(compact);
  mat[3, 2] = detail::f34(compact);
  mat[3, 3] = detail::f33(compact);
  return mat;
}

RegridWeights calc_regrid_weights(
    std::shared_ptr<const Vector> t_grid,
    std::shared_ptr<const Vector> f_grid,
    std::shared_ptr<const Vector> aa_inc_grid,
    std::shared_ptr<const Vector> za_inc_grid,
    std::shared_ptr<const Vector> aa_scat_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
    ScatteringDataGrids new_grids) {
  RegridWeights res{};

  if (!t_grid) {
    ARTS_USER_ERROR(
        "The old t_grid must be provided for calculating regridding weights.");
  }
  if (!new_grids.t_grid) {
    ARTS_USER_ERROR(
        "The new t_grid must be provided for calculating regridding weights.");
  }
  if (!f_grid) {
    ARTS_USER_ERROR(
        "The old f_grid must be provided for calculating regridding weights.");
  }
  if (!new_grids.f_grid) {
    ARTS_USER_ERROR(
        "The new f_grid must be provided for calculating regridding weights.");
  }

  res.t_grid_weights = ArrayOfGridPos(new_grids.t_grid->size());
  gridpos(res.t_grid_weights, *t_grid, *new_grids.t_grid, 1e99);
  res.f_grid_weights = ArrayOfGridPos(new_grids.f_grid->size());
  gridpos(res.f_grid_weights, *f_grid, *new_grids.f_grid, 1e99);

  if ((aa_inc_grid) && (new_grids.aa_inc_grid)) {
    res.aa_inc_grid_weights = ArrayOfGridPos(new_grids.aa_inc_grid->size());
    gridpos(
        res.aa_inc_grid_weights, *aa_inc_grid, *new_grids.aa_inc_grid, 1e99);
  }
  if ((za_inc_grid) && (new_grids.za_inc_grid)) {
    res.za_inc_grid_weights = ArrayOfGridPos(new_grids.za_inc_grid->size());
    gridpos(
        res.za_inc_grid_weights, *za_inc_grid, *new_grids.za_inc_grid, 1e99);
  }
  if ((aa_scat_grid) && (new_grids.aa_scat_grid)) {
    res.aa_scat_grid_weights = ArrayOfGridPos(new_grids.aa_scat_grid->size());
    gridpos(
        res.aa_scat_grid_weights, *aa_scat_grid, *new_grids.aa_scat_grid, 1e99);
  }
  if ((za_scat_grid) && (new_grids.za_scat_grid)) {
    res.za_scat_grid_weights = ArrayOfGridPos(
        std::visit([](const auto &grd) { return grd.angles.size(); },
                   *new_grids.za_scat_grid));
    gridpos(res.za_scat_grid_weights,
            std::visit(
                [](const auto &grd) { return static_cast<Vector>(grd.angles); },
                *za_scat_grid),
            std::visit(
                [](const auto &grd) { return static_cast<Vector>(grd.angles); },
                *new_grids.za_scat_grid),
            1e99);
  }
  return res;
}

std::ostream &operator<<(std::ostream &out, Format format) {
  switch (format) {
    case Format::TRO:     out << "TRO"; break;
    case Format::ARO:     out << "ARO"; break;
    case Format::General: out << "General"; break;
  }
  return out;
}

std::ostream &operator<<(std::ostream &out, Representation repr) {
  switch (repr) {
    case Representation::Gridded:        out << "gridded"; break;
    case Representation::Spectral:       out << "spectral"; break;
    case Representation::DoublySpectral: out << "doubly-spectral"; break;
  }
  return out;
}

}  // namespace scattering
