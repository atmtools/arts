#include <algorithm>

#include "covariance_matrix.h"
#include "jacobian.h"

////////////////////////////////////////////////////////////////////////////////
// Forward declarations of WSVs defined in m_jacobian.cc
////////////////////////////////////////////////////////////////////////////////

void jacobianAddAbsSpecies(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Index&, const Vector&,
    const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
    const String&, const String&, const String&, const Index&, const Numeric&, const Verbosity&
    );
void jacobianClose(
    Workspace&, Index&, ArrayOfArrayOfIndex&, Agenda&, const ArrayOfRetrievalQuantity&,
    const Matrix&, const Sparse&, const Verbosity&
    );
void jacobianInit(
    ArrayOfRetrievalQuantity&, ArrayOfArrayOfIndex&, Agenda&, const Verbosity&
    );
void jacobianAddConstantVMRAbsSpecies(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const String&, const String&, const Index&,
    const Numeric&, const Verbosity&
    );
void jacobianAddBeamFlux(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Index&, const Vector&,
    const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
    const Verbosity&
    );
void jacobianAddFreqShift(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Vector&, const Matrix&,
    const Vector&, const Index&, const Numeric&, const Verbosity&
    );
void jacobianAddFreqStretch(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Vector&, const Matrix&,
    const Vector&, const Index&, const Numeric&, const Verbosity&
    );
void jacobianAddCatalogParameter(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const QuantumIdentifier&,
    const String&, const Verbosity&
    );
void jacobianAddCatalogParameters(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const ArrayOfQuantumIdentifier&,
    const ArrayOfString&, const Verbosity&
    );
void jacobianAddMagField(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Index&, const Vector&,
    const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
    const String&, const Numeric&, const Verbosity&
    );
void jacobianAddPointingZa(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Matrix&, const Vector&,
    const Index&, const String&, const Numeric&, const Verbosity&
    );
void jacobianAddPolyfit(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const ArrayOfIndex&, const Matrix&,
    const Matrix&, const Index&, const Index&, const Index&, const Index&, const Verbosity&
    );
void jacobianAddScatSpecies(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Index&, const Vector&,
    const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
    const String&, const String&, const Verbosity&
    );
void jacobianAddSinefit(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const ArrayOfIndex&, const Matrix&,
    const Matrix&, const Vector&, const Index&, const Index&, const Index&, const Verbosity&
    );
void jacobianAddSpecialSpecies(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Index&, const Vector&,
    const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
    const String&, const Verbosity&);
void jacobianAddWind(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Index&, const Vector&, const Vector&,
    const Vector&, const Vector&, const Vector&, const Vector&, const String&, const Numeric&,
    const Verbosity&
    );
void jacobianAddTemperature(
    Workspace&, ArrayOfRetrievalQuantity&, Agenda&, const Index&, const Vector&,
    const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
    const String&, const String&, const Numeric&, const Verbosity&
    );

////////////////////////////////////////////////////////////////////////////////
// Helper Functions
////////////////////////////////////////////////////////////////////////////////

void check_and_add_block(CovarianceMatrix &covmat,
                         const ArrayOfRetrievalQuantity &jacobian_quantities,
                         const Index grid_dimensions,
                         const Sparse &covmat_block,
                         const Sparse &covmat_inv_block)
{

    Index i = jacobian_quantities.size() - 1;
    Index start = covmat.nrows();
    Index extent = covmat_block.nrows();
    Range range(start, extent);
    Index n_gps = 1;
    for (Index j = 0; j < grid_dimensions; ++j) {
        n_gps *= jacobian_quantities.back().Grids()[j].nelem();
    }

    if ((n_gps == extent) && (n_gps == covmat_block.ncols())) {
        std::shared_ptr<Sparse> mat = make_shared<Sparse>(covmat_block);
        covmat.add_correlation(Block(range, range, std::make_pair(i,i), mat));
    } else {
        std::cout << n_gps << " " << extent << " " << covmat_block.ncols() << std::endl;
        throw runtime_error("Matrix in covmat_block is inconsistent with the retrieval grids.");
    }
    if (!covmat_inv_block.empty()) {
        if ((n_gps == covmat_inv_block.nrows()) && (n_gps == covmat_inv_block.ncols())) {
            std::shared_ptr<Sparse> mat = make_shared<Sparse>(covmat_inv_block);
            covmat.add_correlation_inverse(Block(range, range, std::make_pair(i,i), mat));
        } else {
            throw runtime_error("Matrix in covmat_inv_block is inconsistent with the retrieval"
                                "grids.");
        }
    }
}

void add_scalar_variance(CovarianceMatrix &covmat,
                         ArrayOfRetrievalQuantity &jacobian_quantities,
                         Numeric var)
{
    Index i = jacobian_quantities.size() - 1;
    Index start = covmat.nrows();
    Index extent = 1;
    Range range(start, extent);
    std::shared_ptr<Matrix> mat = make_shared<Matrix>(1,1);
    mat->operator()(0,0) = var;
    covmat.add_correlation(Block(range, range, std::make_pair(i,i), mat));
    mat = make_shared<Matrix>(1,1);
    mat->operator()(0,0) = 1.0 / var;
    covmat.add_correlation_inverse(Block(range, range, std::make_pair(i,i), mat));
}

////////////////////////////////////////////////////////////////////////////////
// Creation of Correlation Blocks
////////////////////////////////////////////////////////////////////////////////

void covmat_blockSetDiagonal(Sparse &covmat_block,
                             Sparse &covmat_inv_block,
                             const Vector &vars,
                             const Verbosity &)
{
    if (vars.empty()) {
        throw runtime_error("Cannot pass empty vector of variances to covmat_blockSetDiagonal");
    }
    Index n = vars.nelem();
    covmat_block = Sparse(n, n);
    covmat_inv_block = Sparse(n,n);

    ArrayOfIndex indices(n);
    Vector       elements(n), elements_inv(n);
    for (Index i = 0; i < n; ++i) {
        indices[i]      = i;
        elements[i]     = vars[i];
        elements_inv[i] = 1.0 / vars[i];
    }

    covmat_block.insert_elements(n, indices, indices, elements);
    covmat_inv_block.insert_elements(n, indices, indices, elements_inv);
}

void covmat_seSet(CovarianceMatrix &covmat_se,
                  const Sparse &covmat_block,
                  const Sparse &covmat_inv_block,
                  const Verbosity &)
{

    Index n =     covmat_block.nrows();

    if (n != covmat_block.ncols()) {
        throw runtime_error("*covmat_block* is not a square matrix.");
    }

    Index i      = 0;
    Index start  = 0;
    Index extent = n;
    Range range(start, extent);
    std::shared_ptr<Sparse> mat = std::make_shared<Sparse>(covmat_block);
    covmat_se.add_correlation(Block(range, range, std::make_pair(i,i), mat));


    if (!covmat_inv_block.empty()) {
        Index n_inv = covmat_inv_block.nrows();
        if (n != n_inv) {
            throw runtime_error("Sizes of *covmat_block* and *covmat_inv_block* don't agree.");
        }

        if (n != covmat_inv_block.ncols()) {
            throw runtime_error("*covmat_inv_block* is not a square matrix.");
        }
        mat = std::make_shared<Sparse>(covmat_inv_block);
        covmat_se.add_correlation_inverse(Block(range, range, std::make_pair(i,i), mat));
    }
}

void Covmat1DMarkov(Sparse& covmat_block,
                    Sparse& covmat_inv_block,
                    const Vector& grid,
                    const Vector& sigma,
                    const Numeric& lc,
                    const Numeric& co)
{
    Index n = grid.nelem();
    ArrayOfIndex row_indices, column_indices;
    Vector       elements, elements_inv;

    Numeric dz = abs(grid[1] - grid[0]);

    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < n; ++j) {
            Numeric e  = sigma[i] * sigma[j] * exp(-dz/lc);
            if (e > co) {
                row_indices.push_back(i);
                column_indices.push_back(j);
            }
        }
    }

    elements = Vector(row_indices.size());
    for (size_t i = 0; i < row_indices.size(); ++i) {
        Numeric e  = sigma[row_indices[i]] * sigma[column_indices[i]] * exp(-dz/lc);
        elements[i] = e;
    }

    covmat_block = Sparse(n,n);
    covmat_block.insert_elements(row_indices.size(), row_indices, column_indices, elements);

    row_indices    = ArrayOfIndex{};
    column_indices = ArrayOfIndex{};
    elements       = Vector(3 * n - 2);

    Numeric alpha = exp(-dz / lc);
    Numeric c1    = - alpha / (1.0 - alpha * alpha);
    Numeric c2    = 1.0 /  (1.0 - alpha * alpha);

    for (Index i = 0; i < n; ++i) {

        // Lower Tri-diagonal
        if (i > 0) {
            column_indices.push_back(i-1);
            row_indices.push_back(i);
            elements[i*3 - 1]  = c1;
            elements[i*3 - 1] /= (sigma[i] * sigma[i-1]);
        }

        // Diagonal Elements
        column_indices.push_back(i);
        row_indices.push_back(i);
        if (i == 0 || i == n-1) {
            elements[i*3] = c2;
        } else {
            elements[i*3] = c2 * (1.0 + alpha * alpha);
        }
        elements[i*3] /= (sigma[i] * sigma[i]);

        // Upper Tri-diagonal
        if (i < n-1) {
            column_indices.push_back(i);
            row_indices.push_back(i+1);
            elements[i*3 + 1] = c1;
            elements[i*3 + 1] /= (sigma[i] * sigma[i+1]);
        }
    }
    covmat_inv_block = Sparse(n,n);
    covmat_inv_block.insert_elements(row_indices.size(), row_indices, column_indices, elements);
}

void Covmat1D(Sparse& covmat_block,
              Sparse& covmat_inv_block,
              const Vector& grid_1,
              const Vector& grid_2,
              const Vector& sigma_1,
              const Vector& sigma_2,
              const Vector& cls_1,
              const Vector& cls_2,
              const Numeric& co,
              const String& fname)
{
    Index m = grid_1.nelem();
    assert(sigma_1.nelem() == m);
    assert(cls_1.nelem() == m);

    Index n = grid_1.nelem();
    assert(sigma_1.nelem() == n);
    assert(cls_1.nelem() == n);

    assert((n == m) || (n==0));

    ConstVectorView grid_view_1(grid_1), sigma_view_1(sigma_1), cls_view_1(cls_1);
    ConstVectorView grid_view_2(grid_2), sigma_view_2(sigma_2), cls_view_2(cls_2);

    if (n==0) {
        grid_view_2  = grid_view_1;
        sigma_view_2 = sigma_view_1;
        cls_view_2   = sigma_view_1;
    }

    // Correlation Functions
    auto f_lin = [&](Index i, Index j) {
        Numeric d  = abs(grid_view_1[i] - grid_view_2[j]);
        Numeric cl = 0.5 * (cls_view_1[i] + cls_view_2[j]);
        return 1.0 - (1.0 - exp(-1.0)) * (d / cl);
    };

    auto f_exp = [&](Index i, Index j) {
        Numeric d  = abs(grid_view_1[i] - grid_view_2[j]);
        Numeric cl = 0.5 * (cls_view_1[i] + cls_view_2[j]);
        return exp(-d/cl);
    };

    auto f_gau = [&](Index i, Index j) {
        Numeric d  = abs(grid_view_1[i] - grid_view_2[j]);
        Numeric cl = 0.5 * (cls_view_1[i] + cls_view_2[j]);
        return exp(-d/cl);
    };

    ArrayOfIndex row_indices;
    ArrayOfIndex column_indices;

    row_indices.reserve(n * m);
    column_indices.reserve(n * m);

    std::function<Numeric(Index, Index)> f;

    if (fname == "exp") {
        f = f_exp;
    } else if (fname == "lin") {
        f = f_lin;
    } else if (fname == "gau") {
        f = f_gau;
    }

    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < n; ++j) {
            Numeric e = sigma_view_1[i] * sigma_view_2[j] * f(i,j);
            if (e > co) {
                row_indices.push_back(i);
                column_indices.push_back(i);
            }
        }
    }

    Vector elements(row_indices.size());
    for (size_t i = 0; i < row_indices.size(); ++i) {
        Index ii = row_indices[i];
        Index jj = column_indices[i];
        elements[i] = sigma_view_1[ii] * sigma_view_2[jj] * f(ii,jj);
    }

    covmat_block = Sparse(m,n);
    covmat_block.insert_elements(column_indices.size(), row_indices, column_indices, elements);
    covmat_inv_block = Sparse(0,0);
}

// void covmat_seAddDiagonalBlock(CovarianceMatrix &covmat_se,
//                                const Sparse &covmat_block,
//                                const Sparse &covmat_inv_block)
// {

//     Index n =     covmat_block.nrows();
//     Index n_inv = covmat_inv_block.nrows();

//     if (n != n_inv) {
//         throw runtime_error("Sizes of *covmat_block* and *covmat_inv_block* don't agree.");
//     }

//     if (n != covmat_block.ncols()) {
//         throw runtime_error("*covmat_block* is not a square matrix.");
//     }

//     if (n != covmat_inv_block.ncols()) {
//         throw runtime_error("*covmat_inv_block* is not a square matrix.");
//     }

//     Index i      = covmat.ndiagblocks();
//     Index start  = covmat_se.nrows();
//     Index extent = n;
//     Range range(start, extent);
//     std::shared_ptr<Sparse> mat = std::make_shared<Sparse>(covmat_block);
//     covmat_se.add_block(Block(range, range, std::make_pair(i,i), mat));

//     mat = std::make_shared<Sparse>(covmat_inv_block);
//     covmat_se.add_inverse_block(Block(range, range, std::make_pair(i,i), mat));
// }

////////////////////////////////////////////////////////////////////////////////
// retrievalAdd Functions
////////////////////////////////////////////////////////////////////////////////

void retrievalAddAbsSpecies(
        Workspace&                  ws,
        CovarianceMatrix &          covmat_sx,
        ArrayOfRetrievalQuantity&   jacobian_quantities,
        Agenda&                     jacobian_agenda,
  const Index&                      atmosphere_dim,
  const Sparse&                     covmat_block,
  const Sparse&                     covmat_inv_block,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Vector&                     rq_p_grid,
  const Vector&                     rq_lat_grid,
  const Vector&                     rq_lon_grid,
  const String&                     species,
  const String&                     method,
  const String&                     mode,
  const Index&                      for_species_tag,
  const Numeric&                    dx,
  const Verbosity&                  verbosity )
{
    jacobianAddAbsSpecies(ws, jacobian_quantities, jacobian_agenda, atmosphere_dim, p_grid,
                          lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid, species,
                          method, mode, for_species_tag, dx, verbosity);
    std::cout << "jq size:" << jacobian_quantities.size() << std::endl;
    check_and_add_block(covmat_sx, jacobian_quantities, atmosphere_dim, covmat_block, covmat_inv_block);
}

void retrievalAddConstantVMRAbsSpecies(
    Workspace&                  ws,
    CovarianceMatrix&           covmat_sx,
    ArrayOfRetrievalQuantity&   jacobian_quantities,
    Agenda&                     jacobian_agenda,
    const String&               species,
    const String&               mode,
    const Index&                for_species_tag,
    const Numeric&              dx,
    const Numeric&              var,
    const Verbosity&            verbosity)
{
    jacobianAddConstantVMRAbsSpecies(ws, jacobian_quantities, jacobian_agenda, species,
                                     mode, for_species_tag, dx, verbosity);
    add_scalar_variance(covmat_sx, jacobian_quantities, var);
}

void retrievalAddFreqShift(Workspace& ws,
                           CovarianceMatrix& covmat_sx,
                           ArrayOfRetrievalQuantity& jacobian_quantities,
                           Agenda& jacobian_agenda,
                           const Sparse& covmat_block,
                           const Sparse& covmat_inv_block,
                           const Vector& f_grid,
                           const Matrix& sensor_pos,
                           const Vector& sensor_time,
                           const Index& poly_order,
                           const Numeric& df,
                           const Verbosity& verbosity)
{
    jacobianAddFreqShift(ws, jacobian_quantities, jacobian_agenda, f_grid, sensor_pos,
                         sensor_time, poly_order, df, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, 1, covmat_block, covmat_inv_block);
}

void retrievalAddFreqStretch(Workspace& ws,
                             CovarianceMatrix& covmat_sx,
                             ArrayOfRetrievalQuantity& jacobian_quantities,
                             Agenda& jacobian_agenda,
                             const Vector& f_grid,
                             const Matrix& sensor_pos,
                             const Vector& sensor_time,
                             const Sparse& covmat_block,
                             const Sparse& covmat_inv_block,
                             const Index& poly_order,
                             const Numeric& df,
                             const Verbosity& verbosity)
{
    jacobianAddFreqStretch(ws, jacobian_quantities, jacobian_agenda, f_grid, sensor_pos,
                         sensor_time, poly_order, df, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, 1, covmat_block, covmat_inv_block);
}

void retrievalAddBeamFlux(
    Workspace&                  ws,
    CovarianceMatrix&           covmat_sx,
    ArrayOfRetrievalQuantity&   jacobian_quantities,
    Agenda&                     jacobian_agenda,
    const Index&                atmosphere_dim,
    const Sparse&               covmat_block,
    const Sparse&               covmat_inv_block,
    const Vector&               p_grid,
    const Vector&               lat_grid,
    const Vector&               lon_grid,
    const Vector&               rq_p_grid,
    const Vector&               rq_lat_grid,
    const Vector&               rq_lon_grid,
    const Verbosity&            verbosity )
{
    jacobianAddBeamFlux(ws, jacobian_quantities, jacobian_agenda, atmosphere_dim, p_grid,
                        lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, atmosphere_dim, covmat_block, covmat_inv_block);
}


void retrievalAddCatalogParameter(
    Workspace&                  ws,
    CovarianceMatrix&           covmat_sx,
    ArrayOfRetrievalQuantity&   jacobian_quantities,
    Agenda&                     jacobian_agenda,
    const QuantumIdentifier&    catalog_identity,
    const String&               catalog_parameter,
    const Numeric&              var,
    const Verbosity&            verbosity )
{
    jacobianAddCatalogParameter(ws, jacobian_quantities, jacobian_agenda, catalog_identity,
                                catalog_parameter, verbosity);
    add_scalar_variance(covmat_sx, jacobian_quantities, var);
}

void retrievalAddCatalogParameters(Workspace& ws,
                                   CovarianceMatrix& covmat_sx,
                                   ArrayOfRetrievalQuantity& jacobian_quantities,
                                   Agenda& jacobian_agenda,
                                   const Sparse& covmat_block,
                                   const Sparse& covmat_inv_block,
                                   const ArrayOfQuantumIdentifier& catalog_identities,
                                   const ArrayOfString& catalog_parameters,
                                   const Verbosity& verbosity )
{
    jacobianAddCatalogParameters(ws, jacobian_quantities, jacobian_agenda,
                                 catalog_identities, catalog_parameters, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, 0, covmat_block, covmat_inv_block);
}

void retrievalAddMagField(Workspace& ws,
                          CovarianceMatrix& covmat_sx,
                          ArrayOfRetrievalQuantity& jacobian_quantities,
                          Agenda& jacobian_agenda,
                          const Index& atmosphere_dim,
                          const Sparse& covmat_block,
                          const Sparse& covmat_inv_block,
                          const Vector& p_grid,
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const Vector& rq_p_grid,
                          const Vector& rq_lat_grid,
                          const Vector& rq_lon_grid,
                          const String& component,
                          const Numeric& dB,
                          const Verbosity& verbosity )
{
    jacobianAddMagField(ws, jacobian_quantities, jacobian_agenda, atmosphere_dim, p_grid,
                        lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid, component, dB,
                        verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, atmosphere_dim, covmat_block, covmat_inv_block);
}

void retrievalAddPointingZa(Workspace& ws,
                            CovarianceMatrix& covmat_sx,
                            ArrayOfRetrievalQuantity& jacobian_quantities,
                            Agenda& jacobian_agenda,
                            const Sparse& covmat_block,
                            const Sparse& covmat_inv_block,
                            const Matrix& sensor_pos,
                            const Vector& sensor_time,
                            const Index& poly_order,
                            const String& calcmode,
                            const Numeric& dza,
                            const Verbosity& verbosity)
{
    jacobianAddPointingZa(ws, jacobian_quantities, jacobian_agenda, sensor_pos,
                          sensor_time, poly_order, calcmode, dza, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, 1, covmat_block, covmat_inv_block);
}

void retrievalAddPolyfit(Workspace& ws,
                        CovarianceMatrix& covmat_sx,
                        ArrayOfRetrievalQuantity& jacobian_quantities,
                        Agenda& jacobian_agenda,
                        const Sparse& covmat_block,
                        const Sparse& covmat_inv_block,
                        const ArrayOfIndex& sensor_response_pol_grid,
                        const Matrix& sensor_response_dlos_grid,
                        const Matrix& sensor_pos,
                        const Index& poly_order,
                        const Index& no_pol_variation,
                        const Index& no_los_variation,
                        const Index& no_mblock_variation,
                        const Verbosity& verbosity)
{
    jacobianAddPolyfit(ws, jacobian_quantities, jacobian_agenda, sensor_response_pol_grid,
                       sensor_response_dlos_grid, sensor_pos, poly_order, no_pol_variation,
                       no_los_variation, no_mblock_variation, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, 4, covmat_block, covmat_inv_block);
}

void retrievalAddScatSpecies(Workspace& ws,
                             CovarianceMatrix& covmat_sx,
                             ArrayOfRetrievalQuantity& jacobian_quantities,
                             Agenda& jacobian_agenda,
                             const Index& atmosphere_dim,
                             const Sparse& covmat_block,
                             const Sparse& covmat_inv_block,
                             const Vector& p_grid,
                             const Vector& lat_grid,
                             const Vector& lon_grid,
                             const Vector& rq_p_grid,
                             const Vector& rq_lat_grid,
                             const Vector& rq_lon_grid,
                             const String& species,
                             const String& quantity,
                             const Verbosity& verbosity)
{
    jacobianAddScatSpecies(ws, jacobian_quantities, jacobian_agenda, atmosphere_dim, p_grid,
                           lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid, species,
                           quantity, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, atmosphere_dim, covmat_block, covmat_inv_block);
}

void retrievalAddSinefit(Workspace& ws,
                         CovarianceMatrix& covmat_sx,
                         ArrayOfRetrievalQuantity& jacobian_quantities,
                         Agenda& jacobian_agenda,
                         const Sparse& covmat_block,
                         const Sparse& covmat_inv_block,
                         const ArrayOfIndex& sensor_response_pol_grid,
                         const Matrix& sensor_response_dlos_grid,
                         const Matrix& sensor_pos,
                         const Vector& period_lengths,
                         const Index& no_pol_variation,
                         const Index& no_los_variation,
                         const Index& no_mblock_variation,
                         const Verbosity& verbosity)
{
    jacobianAddSinefit(ws, jacobian_quantities, jacobian_agenda, sensor_response_pol_grid,
                         sensor_response_dlos_grid, sensor_pos, period_lengths,
                       no_pol_variation, no_los_variation, no_mblock_variation, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, 4, covmat_block, covmat_inv_block);
}

void retrievalAddSpecialSpecies(Workspace& ws,
                                CovarianceMatrix& covmat_sx,
                                ArrayOfRetrievalQuantity& jacobian_quantities,
                                Agenda& jacobian_agenda,
                                const Index& atmosphere_dim,
                                const Sparse& covmat_block,
                                const Sparse& covmat_inv_block,
                                const Vector& p_grid,
                                const Vector& lat_grid,
                                const Vector& lon_grid,
                                const Vector& rq_p_grid,
                                const Vector& rq_lat_grid,
                                const Vector& rq_lon_grid,
                                const String& species,
                                const Verbosity& verbosity)
{
    jacobianAddSpecialSpecies(ws, jacobian_quantities, jacobian_agenda, atmosphere_dim, p_grid,
                    lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid, species,
                    verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, atmosphere_dim, covmat_block, covmat_inv_block);
}

void retrievalAddWind(Workspace& ws,
                      CovarianceMatrix& covmat_sx,
                      ArrayOfRetrievalQuantity& jacobian_quantities,
                      Agenda& jacobian_agenda,
                      const Index& atmosphere_dim,
                      const Sparse& covmat_block,
                      const Sparse& covmat_inv_block,
                      const Vector& p_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const Vector& rq_p_grid,
                      const Vector& rq_lat_grid,
                      const Vector& rq_lon_grid,
                      const String& component,
                      const Numeric& dfrequency,
                      const Verbosity& verbosity)
{
    jacobianAddWind(ws, jacobian_quantities, jacobian_agenda, atmosphere_dim, p_grid,
                    lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid, component,
                    dfrequency, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, atmosphere_dim, covmat_block,
                        covmat_inv_block);
}

void retrievalAddTemperature(Workspace& ws,
                             CovarianceMatrix& covmat_sx,
                             ArrayOfRetrievalQuantity& jacobian_quantities,
                             Agenda& jacobian_agenda,
                             const Index& atmosphere_dim,
                             const Sparse& covmat_block,
                             const Sparse& covmat_inv_block,
                             const Vector& p_grid,
                             const Vector& lat_grid,
                             const Vector& lon_grid,
                             const Vector& rq_p_grid,
                             const Vector& rq_lat_grid,
                             const Vector& rq_lon_grid,
                             const String& hse,
                             const String& method,
                             const Numeric& dx,
                             const Verbosity& verbosity)
{
    jacobianAddTemperature(ws, jacobian_quantities, jacobian_agenda, atmosphere_dim, p_grid,
                           lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid, hse,
                           method, dx, verbosity);
    check_and_add_block(covmat_sx, jacobian_quantities, atmosphere_dim, covmat_block,
                        covmat_inv_block);
}

void retrievalDefClose(Workspace &ws,
                       Index &jacobian_do,
                       ArrayOfArrayOfIndex &jacobian_indices,
                       Agenda &jacobian_agenda,
                       Index  &retrieval_checked,
                       const CovarianceMatrix &/*covmat_se*/,
                       const CovarianceMatrix &covmat_sx,
                       const ArrayOfRetrievalQuantity &jacobian_quantities,
                       const Matrix& sensor_pos,
                       const Sparse& sensor_response,
                       const Verbosity& verbosity)
{
    jacobianClose(ws, jacobian_do, jacobian_indices, jacobian_agenda, jacobian_quantities,
                  sensor_pos, sensor_response, verbosity);

    if (!covmat_sx.has_diagonal_blocks(jacobian_indices)) {
        throw runtime_error("*covmat_sx* does not contain a diagonal block for each retrieval"
                            " quantity in the Jacobian.");
    }
    if (!covmat_sx.is_consistent(jacobian_indices)) {
        throw runtime_error("The blocks in *covmat_sx* are not consistent with the retrieval"
                            " quantities in the Jacobian.");
    }
    retrieval_checked = true;
}

void retrievalDefInit(CovarianceMatrix &covmat_se,
                      CovarianceMatrix &covmat_sx,
                      ArrayOfRetrievalQuantity &jacobian_quantities,
                      ArrayOfArrayOfIndex &jacobian_indices,
                      Agenda &jacobian_agenda,
                      const Verbosity &verbosity)
{
    jacobianInit(jacobian_quantities, jacobian_indices, jacobian_agenda, verbosity);
    covmat_sx = CovarianceMatrix();
    covmat_se = CovarianceMatrix();
}
