/* Copyright (C) 2012
   Richard Larsson <ric.larsson@gmail.com>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

#include "auto_md.h"
#include "zeeman.h"
#include "global_data.h"


/* Workspace method: Doxygen documentation will be auto-generated */
void zeeman_linerecord_precalcCreateFromLines(ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                                              const ArrayOfArrayOfSpeciesTag& abs_species,
                                              const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                              const Index& wigner_initialized,
                                              const Verbosity& verbosity)
{
  if(not wigner_initialized)
    throw std::runtime_error("Must initialize wigner calculations to compute Zeeman effect");
  
  if (abs_species.nelem() != abs_lines_per_species.nelem())
    throw std::runtime_error("Dimension of *abs_species* and *abs_lines_per_species* don't match.");
  
  zeeman_linerecord_precalc.resize(0);
  zeeman_linerecord_precalc.reserve(24); //will always be multiple of three, default is high
  
  // creating the ArrayOfArrayOfLineRecord
  create_Zeeman_linerecordarrays(zeeman_linerecord_precalc, abs_species, abs_lines_per_species, false, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void zeeman_linerecord_precalcCreateWithZeroSplitting(ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                                                      const ArrayOfArrayOfSpeciesTag& abs_species,
                                                      const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                                      const Index& wigner_initialized,
                                                      const Verbosity& verbosity)
{
  if(not wigner_initialized)
    throw std::runtime_error("Must initialize wigner calculations to compute Zeeman effect");
  
  if (abs_species.nelem() != abs_lines_per_species.nelem())
    throw std::runtime_error("Dimension of *abs_species* and *abs_lines_per_species* don't match.");
  
  zeeman_linerecord_precalc.resize(0);
  zeeman_linerecord_precalc.reserve(24); //will always be multiple of three, default is high
  
  // creating the ArrayOfArrayOfLineRecord
  create_Zeeman_linerecordarrays(zeeman_linerecord_precalc, abs_species, abs_lines_per_species, true, verbosity);
}


void zeeman_linerecord_precalcModifyFromData(ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                                             const ArrayOfQuantumIdentifier& keys,
                                             const Vector& data,
                                             const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  if(keys.nelem() not_eq data.nelem()) throw std::runtime_error("Mismatching data and identifier vector");
  
  for(ArrayOfLineRecord& lines : zeeman_linerecord_precalc) {
    Index i=0, j=0;
    for(LineRecord& line: lines) {
      Index upper=-1, lower=-1;
      for(Index k=0; k<keys.nelem(); k++) {
        const QuantumIdentifier& qid = keys[k];
        if(qid < line.QuantumIdentity().LowerQuantumId())
          lower = k;
        else if(qid < line.QuantumIdentity().UpperQuantumId())
          upper = k;
      }
      
      if(lower not_eq -1)
        line.ZeemanEffect().LowerG() = data[lower];
      if(upper not_eq -1)
        line.ZeemanEffect().UpperG() = data[upper];
      
      if(lower not_eq -1 or  upper not_eq -1) ++i;
      if(lower not_eq -1 and upper not_eq -1) ++j;
    }
    out2 << "Modified " << i <<"/"<<lines.nelem() << " lines of which " 
                        << j <<"/"<<lines.nelem() << " were fully modified.\n";
  }
}


void zeeman_linerecord_precalcPrintMissing(const ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                                           const ArrayOfQuantumIdentifier& keys,
                                           const Verbosity& verbosity)
{
  CREATE_OUT0;
  
  ArrayOfArrayOfIndex c(zeeman_linerecord_precalc.nelem());
  
  for(Index i=0; i<c.nelem(); i++) {
    auto& lines = zeeman_linerecord_precalc[i];
    for(Index j=0; j<lines.nelem(); j++) {
      auto& line = lines[j];
      bool found=false;
      for(auto& key: keys) {
        if(key.In(line.QuantumIdentity().LowerQuantumId()))
          found = true;
        if(key.In(line.QuantumIdentity().UpperQuantumId()))
          found = true;
        if(found)
          break;
      }
      
      if(not found)
        c[i].push_back(j);
    }
  }
  
  for(Index i=0; i<c.nelem(); i++) {
    for(auto& x: c[i])
      out0 << "Line is missing in keys: " << zeeman_linerecord_precalc[i][x] << "\n";
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeeman(ArrayOfPropagationMatrix& propmat_clearsky,
                               ArrayOfStokesVector& nlte_source,
                               ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                               ArrayOfStokesVector& dnlte_dx_source,
                               ArrayOfStokesVector& nlte_dsource_dx,
                               const ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                               const Vector& f_grid,
                               const ArrayOfArrayOfSpeciesTag& abs_species,
                               const ArrayOfRetrievalQuantity& jacobian_quantities,
                               const SpeciesAuxData& isotopologue_ratios,
                               const SpeciesAuxData& partition_functions,
                               const Numeric& rtp_pressure,
                               const Numeric& rtp_temperature,
                               const Vector& rtp_nlte,
                               const Vector& rtp_vmr,
                               const Vector& rtp_mag,
                               const Vector& ppath_los,
                               const Index& atmosphere_dim,
                               const Index& manual_zeeman_tag,
                               const Numeric& manual_zeeman_magnetic_field_strength,
                               const Numeric& manual_zeeman_theta,
                               const Numeric& manual_zeeman_eta,
                               const Verbosity&)
{
  // Check that correct isotopologue ratios are defined for the species
  // we want to calculate
  checkIsotopologueRatios(abs_species, isotopologue_ratios);
  
  const Index nzeeman = zeeman_linerecord_precalc.nelem();
  
  bool do_src = !nlte_source.empty();
  {// Begin TEST(s)
    if (abs_species.nelem() == 0)
        throw std::runtime_error("No Zeeman species have been defined.");
    if( propmat_clearsky[0].StokesDimensions()  != 4 )
        throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
    if( propmat_clearsky[0].NumberOfFrequencies() != f_grid.nelem() )
        throw std::runtime_error("Frequency dimension of *propmat_clearsky* not equal to length of *f_grid*.");
    if( propmat_clearsky.nelem() != abs_species.nelem() )
        throw std::runtime_error("Species dimension of *propmat_clearsky* not equal to length of *abs_species*.");
    if( rtp_mag.nelem() != 3 )
        throw std::runtime_error("*rtp_mag* must have length 3.");
    if( atmosphere_dim != 3 )
        throw std::runtime_error("*atmosphere_dim* must be 3.  Zeeman Effect is only implemented for 3D geometry.");
    if( ppath_los.nelem() != 2 )
        throw std::runtime_error("*ppath_los* is not set correctly.");
    if( zeeman_linerecord_precalc.nelem() % 3 != 0 )
        throw std::runtime_error("Length of *zeeman_linerecord_precalc* must be multiple of 3 for polarization states.  It is not.");
    if(do_src) {
        if(nlte_source.nelem() != abs_species.nelem())
        throw std::runtime_error("Species dimension of *nlte_source* not equal to length of *abs_species*.");
        if(nlte_source[0].NumberOfFrequencies() != f_grid.nelem())
        throw std::runtime_error("Frequency dimension of *nlte_source* not equal to length of *f_grid*.");
        if(nlte_source[0].StokesDimensions() != 4)
        throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
    }
  }// End   TEST(s)
  
  if(nzeeman==0)
      return;
  
  Vector rtp_los;
  mirror_los(rtp_los, ppath_los, atmosphere_dim);
  
  // NEW method
  zeeman_on_the_fly(propmat_clearsky, nlte_source, dpropmat_clearsky_dx, dnlte_dx_source, nlte_dsource_dx,
                    abs_species, jacobian_quantities, zeeman_linerecord_precalc, isotopologue_ratios,
                    partition_functions, f_grid, rtp_vmr, rtp_nlte, rtp_mag, rtp_los, rtp_pressure,
                    rtp_temperature, manual_zeeman_tag, manual_zeeman_magnetic_field_strength,
                    manual_zeeman_theta, manual_zeeman_eta);
}
