#include "phase_matrix.h"

#include <ranges>

namespace scattering {
ScatteringDataGrids::ScatteringDataGrids(std::shared_ptr<const Vector> t_grid_, std::shared_ptr<const Vector> f_grid_)
    : t_grid(std::move(t_grid_)),
      f_grid(std::move(f_grid_)),
      aa_inc_grid(nullptr),
      za_inc_grid(nullptr),
      aa_scat_grid(nullptr),
      za_scat_grid(nullptr) {}

ScatteringDataGrids::ScatteringDataGrids(std::shared_ptr<const Vector>          t_grid_,
                                         std::shared_ptr<const Vector>          f_grid_,
                                         std::shared_ptr<const ZenithAngleGrid> za_scat_grid_)
    : t_grid(std::move(t_grid_)),
      f_grid(std::move(f_grid_)),
      aa_inc_grid(nullptr),
      za_inc_grid(nullptr),
      aa_scat_grid(nullptr),
      za_scat_grid(std::move(za_scat_grid_)) {}

ScatteringDataGrids::ScatteringDataGrids(std::shared_ptr<const Vector>          t_grid_,
                                         std::shared_ptr<const Vector>          f_grid_,
                                         std::shared_ptr<const Vector>          za_inc_grid_,
                                         std::shared_ptr<const Vector>          delta_aa_grid_,
                                         std::shared_ptr<const ZenithAngleGrid> za_scat_grid_)
    : t_grid(std::move(t_grid_)),
      f_grid(std::move(f_grid_)),
      aa_inc_grid(nullptr),
      za_inc_grid(std::move(za_inc_grid_)),
      aa_scat_grid(std::move(delta_aa_grid_)),
      za_scat_grid(std::move(za_scat_grid_)) {}

Matrix expand_phase_matrix(const StridedConstVectorView &compact) {
  return Matrix{rtepack::compact_planar_muelmat{compact}.expand().view()};
}

ComplexMatrix expand_phase_matrix(const StridedConstComplexVectorView &compact) {
  // The real and the imaginary parts expand alike
  const auto re = rtepack::compact_planar_muelmat{compact | std::views::transform([](const Complex &x) { return x.real(); })}.expand();
  const auto im = rtepack::compact_planar_muelmat{compact | std::views::transform([](const Complex &x) { return x.imag(); })}.expand();
  ComplexMatrix mat(4, 4);
  for (Index i = 0; i < 4; ++i)
    for (Index j = 0; j < 4; ++j) mat[i, j] = Complex{re[i, j], im[i, j]};
  return mat;
}

RegridWeights calc_regrid_weights(std::shared_ptr<const Vector>          t_grid,
                                  std::shared_ptr<const Vector>          f_grid,
                                  std::shared_ptr<const Vector>          aa_inc_grid,
                                  std::shared_ptr<const Vector>          za_inc_grid,
                                  std::shared_ptr<const Vector>          aa_scat_grid,
                                  std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                  ScatteringDataGrids                    new_grids) {
  RegridWeights res{};

  if (!t_grid) { ARTS_USER_ERROR("The old t_grid must be provided for calculating regridding weights."); }
  if (!new_grids.t_grid) { ARTS_USER_ERROR("The new t_grid must be provided for calculating regridding weights."); }
  if (!f_grid) { ARTS_USER_ERROR("The old f_grid must be provided for calculating regridding weights."); }
  if (!new_grids.f_grid) { ARTS_USER_ERROR("The new f_grid must be provided for calculating regridding weights."); }

  const auto positions = [](ArrayOfGridPos &out, const Vector &old_grid, const Vector &new_grid) {
    if (old_grid.size() == 1) {
      stdr::fill(out, GridPos{.idx = 0, .fd = {0.0, 1.0}});
    } else {
      gridpos(out, old_grid, new_grid, 1e99);
    }
  };

  res.t_grid_weights = ArrayOfGridPos(new_grids.t_grid->size());
  positions(res.t_grid_weights, *t_grid, *new_grids.t_grid);
  res.f_grid_weights = ArrayOfGridPos(new_grids.f_grid->size());
  positions(res.f_grid_weights, *f_grid, *new_grids.f_grid);

  if ((aa_inc_grid) && (new_grids.aa_inc_grid)) {
    res.aa_inc_grid_weights = ArrayOfGridPos(new_grids.aa_inc_grid->size());
    positions(res.aa_inc_grid_weights, *aa_inc_grid, *new_grids.aa_inc_grid);
  }
  if ((za_inc_grid) && (new_grids.za_inc_grid)) {
    res.za_inc_grid_weights = ArrayOfGridPos(new_grids.za_inc_grid->size());
    positions(res.za_inc_grid_weights, *za_inc_grid, *new_grids.za_inc_grid);
  }
  if ((aa_scat_grid) && (new_grids.aa_scat_grid)) {
    res.aa_scat_grid_weights = ArrayOfGridPos(new_grids.aa_scat_grid->size());
    positions(res.aa_scat_grid_weights, *aa_scat_grid, *new_grids.aa_scat_grid);
  }
  if ((za_scat_grid) && (new_grids.za_scat_grid)) {
    res.za_scat_grid_weights =
        ArrayOfGridPos(std::visit([](const auto &grd) { return grd.angles.size(); }, *new_grids.za_scat_grid));
    positions(res.za_scat_grid_weights,
              std::visit([](const auto &grd) { return grd.angles.vec(); }, *za_scat_grid),
              std::visit([](const auto &grd) { return grd.angles.vec(); }, *new_grids.za_scat_grid));
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
