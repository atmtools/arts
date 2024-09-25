/**
 * @file   zeeman.cc
 * @author Richard Larsson <larsson (at) mps.mpg.de>
 * @date   2012-08-14
 * 
 * @brief Public methods of ARTS to compute Zeeman effects
 * 
 * Several methods to change and alter and in other way set up
 * Zeeman effect calculations are implemented in this file
 */

#include "auto_md.h"
#include "debug.h"
#include "global_data.h"
#include "matpack_data.h"
#include "messages.h"
#include "ppath_struct.h"
#include "propagationmatrix.h"
#include "zeeman.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeeman(
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_source_dx,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const EnergyLevelMap& rtp_nlte,
    const Vector& rtp_vmr,
    const Vector& rtp_mag,
    const Vector& ppath_los,
    const Index& atmosphere_dim,
    const Index& nlte_do,
    const Index& lbl_checked,
    const Index& manual_zeeman_tag,
    const Numeric& manual_zeeman_magnetic_field_strength,
    const Numeric& manual_zeeman_theta,
    const Numeric& manual_zeeman_eta,
    const Verbosity&) {
  if (abs_lines_per_species.nelem() == 0) return;

  ARTS_USER_ERROR_IF((atmosphere_dim not_eq 3) and (not manual_zeeman_tag),
                     "Only for 3D *atmosphere_dim* or a manual magnetic field")

  ARTS_USER_ERROR_IF((ppath_los.nelem() not_eq 2) and (not manual_zeeman_tag),
                     "Only for 2D *ppath_los* or a manual magnetic field");

  ARTS_USER_ERROR_IF(not lbl_checked,
                     "Please set lbl_checked true to use this function")

  // Change to LOS by radiation
  Vector rtp_los;
  if (not manual_zeeman_tag) mirror_los(rtp_los, ppath_los, atmosphere_dim);

  // Main computations
  zeeman_on_the_fly(propmat_clearsky,
                    nlte_source,
                    dpropmat_clearsky_dx,
                    dnlte_source_dx,
                    abs_species,
                    select_abs_species,
                    jacobian_quantities,
                    abs_lines_per_species,
                    isotopologue_ratios,
                    f_grid,
                    rtp_vmr,
                    rtp_nlte,
                    rtp_mag,
                    rtp_los,
                    rtp_pressure,
                    rtp_temperature,
                    nlte_do,
                    manual_zeeman_tag,
                    manual_zeeman_magnetic_field_strength,
                    manual_zeeman_theta,
                    manual_zeeman_eta);
}

void abs_linesZeemanCoefficients(ArrayOfAbsorptionLines& abs_lines,
                                 const ArrayOfQuantumIdentifier& qid,
                                 const Vector& gs,
                                 const Verbosity&) {
  ARTS_USER_ERROR_IF(qid.nelem() not_eq gs.nelem(),
                     "Inputs not matching in size");
  for (Index i = 0; i < qid.nelem(); i++) {
    const QuantumIdentifier& id = qid[i];
    const Numeric g = gs[i];

    for (AbsorptionLines& band : abs_lines) {
      if (id.isotopologue_index not_eq band.quantumidentity.isotopologue_index)
        continue;

      for (auto& line : band.lines) {
        auto test = id.part_of(band.quantumidentity, line.localquanta);

        if (test.upp) line.zeeman.gu(g);
        if (test.low) line.zeeman.gl(g);
      }
    }
  }
}

void abs_lines_per_speciesZeemanCoefficients(
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfQuantumIdentifier& qid,
    const Vector& gs,
    const Verbosity& verbosity) {
  for (auto& abs_lines : abs_lines_per_species) {
    for (Index i = 0; i < qid.nelem(); i++) {
      abs_linesZeemanCoefficients(abs_lines, qid, gs, verbosity);
    }
  }
}

void ppvar_magFromPath(Matrix& ppvar_mag,
                       const Tensor3& mag_u_field,
                       const Tensor3& mag_v_field,
                       const Tensor3& mag_w_field,
                       const Ppath& ppath,
                       const Verbosity&) {
  const Index np = ppath.np;
  const Index atmosphere_dim = ppath.dim;

  ARTS_USER_ERROR_IF(mag_u_field.shape() != mag_v_field.shape() or
                         mag_u_field.shape() != mag_w_field.shape() or
                         mag_u_field.size() == 0,
                     "Magnetic field sizes not correct");
  ARTS_USER_ERROR_IF(atmosphere_dim != 3, "Only for 3D atmospheres");

  Matrix itw_field;
  interp_atmfield_gp2itw(
      itw_field, atmosphere_dim, ppath.gp_p, ppath.gp_lat, ppath.gp_lon);

  // Magnetic field:
  ppvar_mag.resize(3, np);
  ppvar_mag = 0;
  //
  interp_atmfield_by_itw(ppvar_mag(0, joker),
                         atmosphere_dim,
                         mag_u_field,
                         ppath.gp_p,
                         ppath.gp_lat,
                         ppath.gp_lon,
                         itw_field);

  interp_atmfield_by_itw(ppvar_mag(1, joker),
                         atmosphere_dim,
                         mag_v_field,
                         ppath.gp_p,
                         ppath.gp_lat,
                         ppath.gp_lon,
                         itw_field);

  interp_atmfield_by_itw(ppvar_mag(2, joker),
                         atmosphere_dim,
                         mag_w_field,
                         ppath.gp_p,
                         ppath.gp_lat,
                         ppath.gp_lon,
                         itw_field);
}

void zeeman_magnetic_fieldCalc(Workspace& ws,
                               ArrayOfMatrix& zeeman_magnetic_field,
                               ArrayOfPpath& zeeman_magnetic_field_path,
                               const Agenda& ppath_agenda,
                               const Numeric& ppath_lmax,
                               const Numeric& ppath_lraytrace,
                               const Vector& f_grid,
                               const Index& cloudbox_on,
                               const Index& ppath_inside_cloudbox_do,
                               const Matrix& sensor_pos,
                               const Matrix& sensor_los,
                               const Tensor3& mag_u_field,
                               const Tensor3& mag_v_field,
                               const Tensor3& mag_w_field,
                               const Verbosity& verbosity) {
  using Conversion::rad2deg, Conversion::deg2rad;

  const Index npaths = sensor_pos.nrows();
  ARTS_USER_ERROR_IF(npaths != sensor_los.nrows(),
                     "Different path lengths in *sensor_pos* and *sensor_los*");
  ARTS_USER_ERROR_IF(sensor_pos.ncols() != 3, "Inputs not matching in size");
  ARTS_USER_ERROR_IF(sensor_los.ncols() != 2, "Inputs not matching in size");
  ARTS_USER_ERROR_IF(mag_u_field.shape() != mag_v_field.shape() or
                         mag_u_field.shape() != mag_w_field.shape() or
                         mag_u_field.size() == 0,
                     "Magnetic field sizes not correct");

  zeeman_magnetic_field.resize(npaths);
  zeeman_magnetic_field_path.resize(npaths);

  for (Index ipath = 0; ipath < npaths; ipath++) {
    Ppath& ppath = zeeman_magnetic_field_path[ipath];
    const Vector rte_pos{sensor_pos[ipath]};
    const Vector rte_los{sensor_los[ipath]};

    ppath_agendaExecute(ws,
                        ppath,
                        ppath_lmax,
                        ppath_lraytrace,
                        rte_pos,
                        rte_los,
                        {},
                        cloudbox_on,
                        ppath_inside_cloudbox_do,
                        f_grid,
                        ppath_agenda);

    ARTS_USER_ERROR_IF(ppath.dim != 3, "Only for 3D atmospheres");

    Matrix ppvar_mag;
    ppvar_magFromPath(
        ppvar_mag, mag_u_field, mag_v_field, mag_w_field, ppath, verbosity);

    zeeman_magnetic_field[ipath].resize(ppath.np, 12);
    for (Index ip = 0; ip < ppath.np; ip++) {
      Vector rtp_los;
      mirror_los(rtp_los, ppath.los[ip], 3);

      const auto [H,
                  theta,
                  eta,
                  dH_du,
                  dH_dv,
                  dH_dw,
                  dtheta_du,
                  dtheta_dv,
                  dtheta_dw,
                  deta_du,
                  deta_dv,
                  deta_dw] = Zeeman::FromGrids(ppvar_mag(0, ip),
                                               ppvar_mag(1, ip),
                                               ppvar_mag(2, ip),
                                               deg2rad(rtp_los[0]),
                                               deg2rad(rtp_los[1]));
      zeeman_magnetic_field[ipath](ip, 0) = H;
      zeeman_magnetic_field[ipath](ip, 1) = rad2deg(theta);
      zeeman_magnetic_field[ipath](ip, 2) = rad2deg(eta);
      zeeman_magnetic_field[ipath](ip, 3) = dH_du;
      zeeman_magnetic_field[ipath](ip, 4) = dH_dv;
      zeeman_magnetic_field[ipath](ip, 5) = dH_dw;
      zeeman_magnetic_field[ipath](ip, 6) = rad2deg(dtheta_du);
      zeeman_magnetic_field[ipath](ip, 7) = rad2deg(dtheta_dv);
      zeeman_magnetic_field[ipath](ip, 8) = rad2deg(dtheta_dw);
      zeeman_magnetic_field[ipath](ip, 9) = rad2deg(deta_du);
      zeeman_magnetic_field[ipath](ip, 10) = rad2deg(deta_dv);
      zeeman_magnetic_field[ipath](ip, 11) = rad2deg(deta_dw);
    }
  }
}
