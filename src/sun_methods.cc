#include "sun_methods.h"

using Constant::pi;

void get_scattered_sunsource(const Workspace& ws,
                              StokvecVector& scattered_sunlight,
                              const Vector& f_grid,
                              const AtmPoint& atm_point,
                              const Matrix& transmitted_sunlight,
                              const Vector& gas_scattering_los_in,
                              const Vector& gas_scattering_los_out,
                              const Agenda& gas_scattering_agenda) {
  PropmatVector K_sca;
  MuelmatVector gas_scattering_mat;
  Vector sca_fct_dummy;

  // calculate gas scattering properties
  gas_scattering_agendaExecute(ws,
                               K_sca,
                               gas_scattering_mat,
                               sca_fct_dummy,
                               f_grid,
                               atm_point,
                               gas_scattering_los_in,
                               gas_scattering_los_out,
                               0,
                               gas_scattering_agenda);

  //some basic quantities
  Index ns = transmitted_sunlight.ncols();
  Index nf = f_grid.size();

  Matrix mat_temp(1, ns,0.);
  // Calculate the scattered radiation
  for (Index i_f = 0; i_f < nf; i_f++) {
    Stokvec scattered_sunlight_temp{transmitted_sunlight[i_f]};
    scattered_sunlight[i_f] = K_sca[i_f].A() / (4*pi) * gas_scattering_mat[i_f] * scattered_sunlight_temp;
  }

  //TODO: Include jacobian mechanism
}
