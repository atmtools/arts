/* Copyright (C) 2020
 * Richard Larsson <larsson@mps.mpg.de>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/**
 * @file m_linemixing.cc
 * @author Richard Larsson
 * @date 2020-06-23
 * 
 * @brief User interface for dealing with pure line mixing calculations
 * 
 * Note: If defined using parameterized form, the normal line-functions
 * approach is faster and more appropriate.  These functions should first
 * compute the relaxation, not simply use the relaxation
 */

#include "global_data.h"
#include "linemixing_hitran.h"
#include "propagationmatrix.h"


void abs_hitran_relmat_dataReadHitranRelmatDataAndLines(HitranRelaxationMatrixData& abs_hitran_relmat_data,
                                                        ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                        const ArrayOfArrayOfSpeciesTag& abs_species,
                                                        const String& basedir,
                                                        const Numeric& linemixinglimit,
                                                        const Numeric& fmin,
                                                        const Numeric& fmax,
                                                        const Numeric& stot,
                                                        const String& mode,
                                                        const Verbosity&)
{
  using global_data::species_data;
  
  lm_hitran_2017::ModeOfLineMixing intmode;
  if (mode == "VP") intmode = lm_hitran_2017::ModeOfLineMixing::VP;
  else if (mode == "VP_Y") intmode = lm_hitran_2017::ModeOfLineMixing::VP_Y;
  else if (mode == "SDVP") intmode = lm_hitran_2017::ModeOfLineMixing::SDVP;
  else if (mode == "SDVP_Y") intmode = lm_hitran_2017::ModeOfLineMixing::SDVP_Y;
  else if (mode == "FullW") intmode = lm_hitran_2017::ModeOfLineMixing::FullW;
  else throw std::runtime_error("Bad mode, see method instruction for valid arguments");
  
  if (abs_species.nelem() not_eq abs_lines_per_species.nelem())
    throw std::runtime_error("Bad size of input species+lines");
  
  ArrayOfAbsorptionLines lines;
  lm_hitran_2017::read(abs_hitran_relmat_data, lines, basedir, linemixinglimit, fmin, fmax, stot, intmode);
  std::for_each(lines.begin(), lines.end(), [](auto& band){band.sort_by_frequency();});  // Sort so largest frequency is last
  ArrayOfIndex used(lines.nelem(), false);
  
  bool emptied=false;
  for (Index i=0; i<abs_species.nelem(); i++) {
    
    for (Index j=0; j<abs_species[i].nelem(); j++) {
      if (abs_species[i][j].Species() not_eq SpeciesTag("CO2").Species()) 
        continue;
      
      if (not emptied) {
        abs_lines_per_species[i].resize(0);
        emptied = true;
      }
      
      for (Index k=0; k<lines.nelem(); k++) {
        if (used[k]) continue;
        
        const Numeric lf{abs_species[i][j].Lf() > 0 ? abs_species[i][j].Lf() : -std::numeric_limits<Numeric>::max()};
        const Numeric uf{abs_species[i][j].Uf() > 0 ? abs_species[i][j].Uf() : std::numeric_limits<Numeric>::max()};
        
        // Select lines with correct Isotopologue and one line center within the range
        if ((abs_species[i][j].Isotopologue() == lines[k].Isotopologue() or 
             abs_species[i][j].Isotopologue() == species_data[SpeciesTag("CO2").Species()].Isotopologue().nelem()) and 
             (lines[k].AllLines().front().F0() <= uf and lines[k].AllLines().back().F0() >= lf)) {
          used[k] = true;  // The lines should not be copied into other places
          abs_lines_per_species[i].push_back(lines[k]);
        }
      }
    }
  }
}

void propmat_clearskyAddHitranLineMixingLines(ArrayOfPropagationMatrix& propmat_clearsky,
                                              const HitranRelaxationMatrixData& abs_hitran_relmat_data,
                                              const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                              const Vector& f_grid,
                                              const ArrayOfArrayOfSpeciesTag& abs_species,
                                              const ArrayOfRetrievalQuantity& jacobian_quantities,
                                              const SpeciesAuxData& partition_functions,
                                              const Numeric& rtp_pressure,
                                              const Numeric& rtp_temperature,
                                              const Vector& rtp_vmr,
                                              const Verbosity&)
{
  if (jacobian_quantities.nelem())
    throw std::runtime_error("Cannot support any Jacobian at this time");
  if (abs_species.nelem() not_eq abs_lines_per_species.nelem())
    throw std::runtime_error("Bad size of input species+lines");
  if (abs_species.nelem() not_eq rtp_vmr.nelem())
    throw std::runtime_error("Bad size of input species+vmrs");
  
  // vmrs should be [air, water, co2]
  Vector vmrs(3, 0);
  for (Index i=0; i<abs_species.nelem(); i++) {
    auto& specs = abs_species[i];
    for (auto& spec: specs) {
      if (SpeciesTag("H2O").Species() == spec.Species()) {
        vmrs[1] = rtp_vmr[i];
      }
      else if (SpeciesTag("CO2").Species() == spec.Species()) {
        vmrs[2] = rtp_vmr[i];
      }
    }
  }
  vmrs[0] = 1 - vmrs[1] - vmrs[2];
    
  for (Index i=0; i<abs_species.nelem(); i++) {
    if (abs_lines_per_species[i].nelem() and 
      (abs_lines_per_species[i].front().Population() == Absorption::PopulationType::ByHITRANFullRelmat or
       abs_lines_per_species[i].front().Population() == Absorption::PopulationType::ByHITRANRosenkranzRelmat))
      propmat_clearsky[i].Kjj() += lm_hitran_2017::compute(abs_hitran_relmat_data, abs_lines_per_species[i], rtp_pressure, rtp_temperature, vmrs, f_grid, partition_functions);
  }
}
