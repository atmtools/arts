/* Copyright (C) 2000-2012
   Stefan Buehler  <sbuehler@ltu.se>
   Axel von Engeln <engeln@uni-bremen.de>

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

/**
  \file   absorption.cc

  Physical absorption routines. 

  The absorption workspace methods are
  in file m_abs.cc

  This is the file from arts-1-0, back-ported to arts-1-1.

  \author Stefan Buehler and Axel von Engeln
*/

#include "arts.h"
#include "auto_md.h"
#include <map>
#include <cfloat>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "file.h"
#include "absorption.h"
#include "math_funcs.h"
#include "messages.h"
#include "logic.h"
#include "interpolation_poly.h"
#include "linemixingdata.h"
#include "linescaling.h"

#include "global_data.h"
#include "linefunctions.h"


/** Mapping of species auxiliary type names to SpeciesAuxData::AuxType enum */
static const char *SpeciesAuxTypeNames[] = {
    "NONE",
    "ISORATIO",   //Built-in type
    "ISOQUANTUM",
    "PART_TFIELD",
    "PART_COEFF", //Built-in type
    "PART_COEFF_VIBROT"
};


// member fct of isotopologuerecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopologueRecord::CalculatePartitionFctAtTempFromCoeff_OldUnused( Numeric
                                                    temperature ) const
{
  Numeric result = 0.;
  Numeric exponent = 1.;

  Vector::const_iterator it;

//      cout << "T: " << temperature << "\n";
  for (it=mqcoeff.begin(); it != mqcoeff.end(); ++it)
    {
      result += *it * exponent;
      exponent *= temperature;
//      cout << "it: " << it << "\n";
//      cout << "res: " << result << ", exp: " << exponent << "\n";
    }
  return result;
}


// member fct of isotopologuerecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopologueRecord::CalculatePartitionFctAtTempFromData_OldUnused( Numeric
temperature ) const
{
    GridPosPoly gp;
    gridpos_poly( gp, mqcoeffgrid, temperature, mqcoeffinterporder);
    Vector itw;
    interpweights(itw, gp );
    return interp(  itw, mqcoeff, gp );
}


void SpeciesAuxData::InitFromSpeciesData()
{
    using global_data::species_data;

    mparams.resize(species_data.nelem());
    mparam_type.resize(species_data.nelem());

    for (Index isp = 0; isp < species_data.nelem(); isp++)
    {
        const Index niso = species_data[isp].Isotopologue().nelem();
        mparams[isp].resize(niso);
        mparam_type[isp].resize(niso);
        for (Index iso = 0; iso < niso; iso++)
        {
            mparams[isp][iso].resize(0);
            mparam_type[isp][iso] = SpeciesAuxData::AT_NONE;
        }
    }
}


void SpeciesAuxData::setParam(const Index species,
                              const Index isotopologue,
                              const AuxType auxtype,
                              const ArrayOfGriddedField1& auxdata)
{
    mparam_type[species][isotopologue] = auxtype;
    mparams[species][isotopologue] = auxdata;
}


void SpeciesAuxData::setParam(const String& artstag,
                              const String& auxtype,
                              const ArrayOfGriddedField1& auxdata)
{
    // Global species lookup data:
    using global_data::species_data;

    // We need a species index sorted by Arts identifier. Keep this in a
    // static variable, so that we have to do this only once.  The ARTS
    // species index is ArtsMap[<Arts String>].
    static map<String, SpecIsoMap> ArtsMap;

    // Remember if this stuff has already been initialized:
    static bool hinit = false;

    if ( !hinit )
    {
        for ( Index i=0; i<species_data.nelem(); ++i )
        {
            const SpeciesRecord& sr = species_data[i];
            for ( Index j=0; j<sr.Isotopologue().nelem(); ++j)
            {
                SpecIsoMap indicies(i,j);
                String buf = sr.Name()+"-"+sr.Isotopologue()[j].Name();
                ArtsMap[buf] = indicies;
            }
        }
        hinit = true;
    }

    Index species;
    Index isotopologue;

    // ok, now for the cool index map:
    // is this arts identifier valid?
    const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artstag);
    if ( i == ArtsMap.end() )
    {
        ostringstream os;
        os << "ARTS Tag: " << artstag << " is unknown.";
        throw runtime_error(os.str());
    }

    SpecIsoMap id = i->second;

    // Set mspecies:
    species = id.Speciesindex();

    // Set misotopologue:
    isotopologue = id.Isotopologueindex();

    Index this_auxtype = 0;

    while (this_auxtype < AT_FINAL_ENTRY && auxtype != SpeciesAuxTypeNames[this_auxtype])
        this_auxtype++;

    if (this_auxtype != AT_FINAL_ENTRY)
    {
        setParam(species, isotopologue, (AuxType)this_auxtype, auxdata);
    }
    else
    {
        ostringstream os;
        os << "Unknown SpeciesAuxData type: " << auxtype;
        std::runtime_error(os.str());
    }
}


const ArrayOfGriddedField1& SpeciesAuxData::getParam(const Index species,
                                                     const Index isotopologue) const
{
    return mparams[species][isotopologue];
}


String SpeciesAuxData::getTypeString(const Index species, const Index isotopologue) const
{
    assert(mparam_type[species][isotopologue] < AT_FINAL_ENTRY);
    return SpeciesAuxTypeNames[mparam_type[species][isotopologue]];
}


bool SpeciesAuxData::ReadFromStream(String& artsid, istream& is, Index nparams, const Verbosity& verbosity)
{
    CREATE_OUT3;

    // Global species lookup data:
    using global_data::species_data;

    // We need a species index sorted by Arts identifier. Keep this in a
    // static variable, so that we have to do this only once.  The ARTS
    // species index is ArtsMap[<Arts String>].
    static map<String, SpecIsoMap> ArtsMap;

    // Remember if this stuff has already been initialized:
    static bool hinit = false;

    if ( !hinit )
    {

        out3 << "  ARTS index table:\n";

        for ( Index i=0; i<species_data.nelem(); ++i )
        {
            const SpeciesRecord& sr = species_data[i];


            for ( Index j=0; j<sr.Isotopologue().nelem(); ++j)
            {

                SpecIsoMap indicies(i,j);
                String buf = sr.Name()+"-"+sr.Isotopologue()[j].Name();

                ArtsMap[buf] = indicies;

                // Print the generated data structures (for debugging):
                // The explicit conversion of Name to a c-String is
                // necessary, because setw does not work correctly for
                // stl Strings.
                const Index& i1 = ArtsMap[buf].Speciesindex();
                const Index& i2 = ArtsMap[buf].Isotopologueindex();

                out3 << "  Arts Identifier = " << buf << "   Species = "
                << setw(10) << setiosflags(ios::left)
                << species_data[i1].Name().c_str()
                << "iso = "
                << species_data[i1].Isotopologue()[i2].Name().c_str()
                << "\n";

            }
        }
        hinit = true;
    }


    // This always contains the rest of the line to parse. At the
    // beginning the entire line. Line gets shorter and shorter as we
    // continue to extract stuff from the beginning.
    String line;

    // Look for more comments?
    bool comment = true;

    while (comment)
    {
        // Return true if eof is reached:
        if (is.eof()) return true;

        // Throw runtime_error if stream is bad:
        if (!is) throw runtime_error ("Stream bad.");

        // Read line from file into linebuffer:
        getline(is,line);

        // It is possible that we were exactly at the end of the file before
        // calling getline. In that case the previous eof() was still false
        // because eof() evaluates only to true if one tries to read after the
        // end of the file. The following check catches this.
        if (line.nelem() == 0 && is.eof()) return true;

        // @ as first character marks catalogue entry
        char c;
        extract(c,line,1);

        // check for empty line
        if (c == '@')
        {
            comment = false;
        }
    }


    // read the arts identifier String
    istringstream icecream(line);

    icecream >> artsid;

    if (artsid.length() != 0)
    {
        Index mspecies;
        Index misotopologue;

        // ok, now for the cool index map:
        // is this arts identifier valid?
        const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artsid);
        if ( i == ArtsMap.end() )
        {
            ostringstream os;
            os << "ARTS Tag: " << artsid << " is unknown.";
            throw runtime_error(os.str());
        }

        SpecIsoMap id = i->second;


        // Set mspecies:
        mspecies = id.Speciesindex();

        // Set misotopologue:
        misotopologue = id.Isotopologueindex();

        ArrayOfGriddedField1 ratios;
        ratios.resize(1);
        // Extract accuracies:
        try
        {
            Numeric p = NAN;
            std::vector<Numeric> aux;
            for (Index ip = 0; ip < nparams; ip++)
            {
                icecream >> double_imanip() >> p;
                aux.push_back(p);
            }

            Vector grid;
            if (aux.size() > 1)
                nlinspace(grid, 1, (Numeric)aux.size(), aux.size());
            else
                grid = Vector(1, .1);

            ratios[0].set_grid(0, grid);
            ratios[0].data = aux;
            mparams[mspecies][misotopologue] = ratios;
        }
        catch (runtime_error)
        {
            throw runtime_error("Error reading SpeciesAuxData.");
        }
    }
    
    // That's it!
    return false;
}


void checkIsotopologueRatios(const ArrayOfArrayOfSpeciesTag& abs_species,
                             const SpeciesAuxData& isoratios)
{
    using global_data::species_data;
    
    // Check total number of species:
    if (species_data.nelem() != isoratios.nspecies())
      {
        ostringstream os;
        os << "Number of species in SpeciesAuxData (" << isoratios.nspecies()
        << "does not fit builtin species data (" << species_data.nelem() << ").";
        throw runtime_error(os.str());
      }
    
    // For the selected species, we check all isotopes by looping over the
    // species data. (Trying to check only the isotopes actually used gets
    // quite complicated, actually, so we do the simple thing here.)
    
    // Loop over the absorption species:
    for (Index i=0; i<abs_species.nelem(); i++)
      {
        // sp is the index of this species in the internal lookup table
        const Index sp = abs_species[i][0].Species();
        
        // Get handle on species data for this species:
        const SpeciesRecord& this_sd = species_data[sp];
     
        // Check number of isotopologues:
        if (this_sd.Isotopologue().nelem() != isoratios.nisotopologues(sp))
          {
            ostringstream os;
            os << "Incorrect number of isotopologues in isotopologue data.\n"
            << "Species: " << this_sd.Name() << ".\n"
            << "Number of isotopes in SpeciesAuxData ("
            << isoratios.nisotopologues(sp) << ") "
            << "does not fit builtin species data (" << this_sd.Isotopologue().nelem() << ").";
            throw runtime_error(os.str());
          }
        
        for (Index iso = 0; iso < this_sd.Isotopologue().nelem(); ++iso)
          {
            // For "real" species (not representing continau) the isotopologue
            // ratio must not be NAN or below zero.
            if (!this_sd.Isotopologue()[iso].isContinuum()) {
                if (isnan(isoratios.getParam(sp, iso)[0].data[0]) ||
                    isoratios.getParam(sp, iso)[0].data[0] < 0.) {
                    
                    ostringstream os;
                    os << "Invalid isotopologue ratio.\n"
                    << "Species: " << this_sd.Name() << "-"
                    << this_sd.Isotopologue()[iso].Name() << "\n"
                    << "Ratio:   " << isoratios.getParam(sp, iso)[0].data[0];
                    throw runtime_error(os.str());
                }
            }
          }
      }
}


void fillSpeciesAuxDataWithIsotopologueRatiosFromSpeciesData(SpeciesAuxData& sad)
{
    using global_data::species_data;

    sad.InitFromSpeciesData();

    Vector grid(1, 1.);
    ArrayOfGriddedField1 ratios;
    ratios.resize(1);
    ratios[0].set_name("IsoRatios");
    ratios[0].set_grid_name(0, "Index");
    ratios[0].set_grid(0, grid);
    ratios[0].resize(1);

    for (Index isp = 0; isp < species_data.nelem(); isp++)
    {
        for (Index iiso = 0; iiso < species_data[isp].Isotopologue().nelem(); iiso++)
        {
            ratios[0].data[0] = species_data[isp].Isotopologue()[iiso].Abundance();
            sad.setParam(isp, iiso,
                         SpeciesAuxData::AT_ISOTOPOLOGUE_RATIO, ratios);
        }
    }
}


void fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(SpeciesAuxData& sad)
{
    using global_data::species_data;

    sad.InitFromSpeciesData();

    ArrayOfGriddedField1 pfuncs;
    pfuncs.resize(2);
    pfuncs[0].set_name("PartitionFunction");
    pfuncs[0].set_grid_name(0, "Coeff");
    pfuncs[1].set_grid_name(0, "Temperature");

    ArrayOfString tgrid;
    tgrid.resize(2);
    tgrid[0] = "Tlower";
    tgrid[1] = "Tupper";

    for (Index isp = 0; isp < species_data.nelem(); isp++)
    {
        for (Index iiso = 0; iiso < species_data[isp].Isotopologue().nelem(); iiso++)
        {
            Vector grid;
            const Vector& coeffs = species_data[isp].Isotopologue()[iiso].GetCoeff();

            assert(coeffs.nelem() >= 2);

            nlinspace(grid, 0, (Numeric)coeffs.nelem()-1., coeffs.nelem());
            pfuncs[0].set_grid(0, grid);
            pfuncs[0].data = coeffs;

            const Vector& temp_range = species_data[isp].Isotopologue()[iiso].GetCoeffGrid();

            // Temperature data should either contain two Ts, lower and upper value of
            // the valid range or be empty
            assert(temp_range.nelem() == 0 || temp_range.nelem() == 2);

            if (temp_range.nelem() == 2)
            {
                pfuncs[1].set_grid(0, tgrid);
                pfuncs[1].data = temp_range;
            }
            else
            {
                pfuncs[1].resize(0);
                pfuncs[1].set_grid(0, ArrayOfString());
            }

            sad.setParam(isp, iiso,
                         SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF, pfuncs);
        }
    }
}


ostream& operator<< (ostream& os, const SpeciesRecord& sr)
{
  for ( Index j=0; j<sr.Isotopologue().nelem(); ++j)
    {
      os << sr.Name() << "-" << sr.Isotopologue()[j].Name() << "\n";
    }
  return os;
}
 
      
ostream& operator<< (ostream& os, const SpeciesAuxData& sad)
{
    using global_data::species_data;
    for (Index sp = 0; sp < sad.nspecies(); sp++)
    {
        for (Index iso = 0; iso < sad.nisotopologues(sp); iso++)
        {
            os << species_name_from_species_index(sp) << "-"
            << global_data::species_data[sp].Isotopologue()[iso].Name();
            os << " " << sad.getTypeString(sp, iso) << std::endl;
            for (Index ip = 0; ip < sad.getParam(sp, iso).nelem(); ip++)
                os << "AuxData " << ip << " " << sad.getParam(sp, iso) << std::endl;
        }
    }

    return os;
}



/** Find the location of all broadening species in abs_species. 
 
 Set to -1 if
 not found. The length of array broad_spec_locations is the number of allowed
 broadening species (in ARTSCAT-4 N2, O2, H2O, CO2, H2, H2). The value
 means:
 <p> -1 = not in abs_species
 <p> -2 = in abs_species, but should be ignored because it is identical to Self
 <p> N  = species is number N in abs_species
 
 The catalogue contains also broadening parameters for "Self", but we ignore this 
 here, since we know the position of species "Self" anyway.
 
 \retval broad_spec_locations See above.
 \param abs_species List of absorption species
 \param this_species Index of the current species in abs_species

 \author Stefan Buehler
 \date   2012-09-04
 
**/
void find_broad_spec_locations(ArrayOfIndex& broad_spec_locations,
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const Index this_species)
{
  // Make sure this_species points to somewhere inside abs_species:
  assert(this_species>=0);
  assert(this_species<abs_species.nelem());
  
  // Number of broadening species:
  const Index nbs = LineRecord::NBroadSpec();
  
  // Resize output array:
  broad_spec_locations.resize(nbs);
  
//  // Set broadening species names, using the enums defined in absorption.h.
//  // This is hardwired here and quite primitive, but should do the job.
//  ArrayOfString broad_spec_names(nbs);
//  broad_spec_names[LineRecord::SPEC_POS_N2]  = "N2";
//  broad_spec_names[LineRecord::SPEC_POS_O2]  = "O2";
//  broad_spec_names[LineRecord::SPEC_POS_H2O] = "H2O";
//  broad_spec_names[LineRecord::SPEC_POS_CO2] = "CO2";
//  broad_spec_names[LineRecord::SPEC_POS_H2]  = "H2";
//  broad_spec_names[LineRecord::SPEC_POS_He]  = "He";
  
  // Loop over all broadening species and see if we can find them in abs_species.
  for (Index i=0; i<nbs; ++i) {
    // Find associated internal species index (we do the lookup by index, not by name).
    const Index isi = LineRecord::BroadSpecSpecIndex(i);
    
    // First check if this broadening species is the same as this_species
    if ( isi == abs_species[this_species][0].Species() )
      broad_spec_locations[i] = -2;
    else
    {
      // Find position of broadening species isi in abs_species. The called
      // function returns -1 if not found, which is already the correct
      // treatment for this case that we also want here.
      broad_spec_locations[i] = find_first_species_tg(abs_species,isi);
    }
  }
}



/** Calculate line width and pressure shift for artscat4.
 
 \retval gamma Line width [Hz].
 \retval deltaf Pressure shift [Hz].
 \param p Pressure [Pa].
 \param  t Temperature [K].
 \param  vmrs Vector of VMRs for different species [dimensionless].
 \param  this_species Index of current species in vmrs.
 \param  broad_spec_locations Has length of number of allowed broadening species
 (6 in artscat-4). Gives for each species the position
 in vmrs, or negative if it should be ignored. See
 function find_broad_spec_locations for details.
 \param  T0  The temperature at which broadening parameters are determined.
 \param  Sgam  The self pressure broadening.
 \param  Nself  The self broadening temperature exponent.
 \param  Gamma_foreign  The foregin pressure broadenings.
 \param  N_foreign  The foregin temperature exponents.
 \param  Delta_foreign  The pressure dependent frequency shift.
 \param  verbosity Verbosity flag.
 
 \author Stefan Buehler
 \date   2012-09-05
 
 Removed LineRecord dependency.
 
 \author Richard Larsson
 \date   2014-10-29
 */
void calc_gamma_and_deltaf_artscat4(Numeric& gamma,
                                    Numeric& deltaf,
                                    const Numeric p,
                                    const Numeric t,
                                    ConstVectorView vmrs,
                                    const Index this_species,
                                    const ArrayOfIndex& broad_spec_locations,
                                    const Numeric& T0,
                                    const Numeric& Sgam,
                                    const Numeric& Nself,
                                    const Vector&  Gamma_foreign,
                                    const Vector&  N_foreign,
                                    const Vector&  Delta_foreign,
                                    const Verbosity& verbosity)
{
    CREATE_OUT2;
    
    // Number of broadening species:
    static const Index nbs = LineRecord::NBroadSpec();
    assert(nbs==broad_spec_locations.nelem());
    
    // Theta is reference temperature divided by local temperature. Used in
    // several place by the broadening and shift formula.
    const Numeric theta = T0 / t;
    
    // Split total pressure in self and foreign part:
    const Numeric p_self    = vmrs[this_species] * p;
    const Numeric p_foreign = p-p_self;
    
    // Calculate sum of VMRs of all available foreign broadening species (we need this
    // for normalization). The species "Self" will not be included in the sum!
    Numeric broad_spec_vmr_sum = 0;
    
    // Gamma is the line width. We first initialize gamma with the self width
    gamma =  Sgam * pow(theta, Nself) * p_self;
    
    // and treat foreign width separately:
    Numeric gamma_foreign = 0;
    
    // There is no self shift parameter (or rather, we do not have it), so
    // we do not need separate treatment of self and foreign for the shift:
    deltaf = 0;
    
    // Add up foreign broadening species, where available:
    for (Index i=0; i<nbs; ++i) {
        if ( broad_spec_locations[i] < -1 ) {
            // -2 means that this broadening species is identical to Self.
            // Throw runtime errors if the parameters are not identical.
            if (Gamma_foreign[i]!=Sgam ||
                N_foreign[i]!=Nself)
            {
                std::ostringstream os;
                os << "Inconsistency in LineRecord, self broadening and line "
                << "broadening for " << LineRecord::BroadSpecName(i) << "\n"
                << "should be identical.\n";
                throw std::runtime_error(os.str());
            }
        } else if ( broad_spec_locations[i] >= 0 ) {
            
            // Add to VMR sum:
            broad_spec_vmr_sum += vmrs[broad_spec_locations[i]];
            
            // foreign broadening:
            gamma_foreign +=  Gamma_foreign[i] * pow(theta, N_foreign[i])
            * vmrs[broad_spec_locations[i]];
            
            // Pressure shift:
            // The T dependence is connected to that of the corresponding
            // broadening parameter by:
            // n_shift = .25 + 1.5 * n_gamma
            deltaf += Delta_foreign[i]
                      * pow( theta , (Numeric).25 + (Numeric)1.5*N_foreign[i] )
                      * vmrs[broad_spec_locations[i]];
        }
    }

    // Check that sum of self and all foreign VMRs is not too far from 1:
    if ( abs(vmrs[this_species]+broad_spec_vmr_sum-1) > 0.1
         && out2.sufficient_priority() )
      {
        std::ostringstream os;
        os << "Warning: The total VMR of all your defined broadening\n"
             << "species (including \"self\") is "
             << vmrs[this_species]+broad_spec_vmr_sum
             << ", more than 10% " << "different from 1.\n";
        out2 << os.str();
      }
    
    // Normalize foreign gamma and deltaf with the foreign VMR sum (but only if
    // we have any foreign broadening species):
    if (broad_spec_vmr_sum != 0.)
      {
        gamma_foreign /= broad_spec_vmr_sum;
        deltaf        /= broad_spec_vmr_sum;
      }
    else if (p_self > 0.)
    // If there are no foreign broadening species present, the best assumption
    // we can make is to use gamma_self in place of gamma_foreign. for deltaf
    // there is no equivalent solution, as we don't have a Delta_self and don't
    // know which other Delta we should apply (in this case delta_f gets 0,
    // which should be okayish):
      {
        gamma_foreign = gamma/p_self;
      }
    // It can happen that broad_spec_vmr_sum==0 AND p_self==0 (e.g., when p_grid
    // exceeds the given atmosphere and zero-padding is applied). In this case,
    // both gamma_foreign and deltaf are 0 and we leave it like that.

    // Multiply by pressure. For the width we take only the foreign pressure.
    // This is consistent with that we have scaled with the sum of all foreign
    // broadening VMRs. In this way we make sure that the total foreign broadening
    // scales with the total foreign pressure.
    gamma_foreign  *= p_foreign;
    // For the shift we simply take the total pressure, since there is no self part.
    deltaf *= p;

    // For the width, add foreign parts:
    gamma += gamma_foreign;
    
    // That's it, we're done.
}


/** Calculate line absorption cross sections for one tag group. All
    lines in the line list must belong to the same species. This must
    be ensured by abs_lines_per_speciesCreateFromLines, so it is only verified
    with assert. Also, the input vectors abs_p, and abs_t must all
    have the same dimension.

    This is mainly a copy of abs_species which is removed now, with
    the difference that the vmrs are removed from the absorption
    coefficient calculation. (the vmr is still used for the self
    broadening)

    Continua are not handled by this function, you have to call
    xsec_continuum_tag for those.

    \retval xsec_attenuation   Cross section of one tag group. This is now the
                               true absorption cross section in units of m^2.
    \retval xsec_phase         Cross section of one tag group. This is now the
                               true dispersion cross section in units of m^2.
    \param f_grid       Frequency grid.
    \param abs_p        Pressure grid.
    \param abs_t        Temperatures associated with abs_p.
    \param all_vmrs     Gas volume mixing ratios [nspecies, np].
    \param abs_species  Species tags for all species.
    \param this_species Index of the current species in abs_species.
    \param abs_lines    The spectroscopic line list.
    \param ind_ls       Index to lineshape function.
    \param ind_lsn      Index to lineshape norm.
    \param cutoff       Lineshape cutoff.
    \param isotopologue_ratios  Isotopologue ratios.

    \author Stefan Buehler and Axel von Engeln
    \date   2001-01-11 

    Changed from pseudo cross sections to true cross sections

    \author Stefan Buehler 
    \date   2007-08-08 

    Adapted to new Perrin line parameters, treating broadening by different 
    gases explicitly
 
    \author Stefan Buehler
    \date   2012-09-03
    
    Adapted to LineRecord independency in central parts of function.
    
    \author Richard Larsson
    \date   2014-10-29

*/
void xsec_species( MatrixView               xsec_attenuation,
                   MatrixView               xsec_source,
                   MatrixView               xsec_phase,
                   ConstVectorView          f_grid,
                   ConstVectorView          abs_p,
                   ConstVectorView          abs_t,
                   ConstMatrixView          abs_t_nlte,
                   ConstMatrixView          all_vmrs,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const Index              this_species,
                   const ArrayOfLineRecord& abs_lines,
                   const Index              ind_ls,
                   const Index              ind_lsn,
                   const Numeric            cutoff,
                   const SpeciesAuxData&    isotopologue_ratios,
                   const SpeciesAuxData&    partition_functions,
                   const Verbosity&         verbosity )
{
    
    extern const Numeric BOLTZMAN_CONST;
    extern const Numeric AVOGADROS_NUMB;
    extern const Numeric SPEED_OF_LIGHT;
    static const Numeric doppler_const = sqrt(2.0 * BOLTZMAN_CONST *
                                              AVOGADROS_NUMB) / SPEED_OF_LIGHT;    
    
    // dimension of f_grid, abs_lines
    const Index nf = f_grid.nelem();
    const Index nl = abs_lines.nelem();
    
    // number of pressure levels:
    const Index np = abs_p.nelem();
    
    // Define the vector for the line shape function and the
    // normalization factor of the lineshape here, so that we don't need
    // so many free store allocations.  the last element is used to
    // calculate the value at the cutoff frequency
    Vector ls_attenuation(nf+1);
    Vector ls_phase_dummy;
    Vector fac(nf+1);
    
    const bool cut = (cutoff != -1) ? true : false;
    
    const bool calc_phase = 0;
    
    // Check that the frequency grid is sorted in the case of lineshape
    // with cutoff. Duplicate frequency values are allowed.
    if (cut)
    {
        if ( ! is_sorted( f_grid ) )
        {
            std::ostringstream os;
            os << "If you use a lineshape function with cutoff, your\n"
               << "frequency grid *f_grid* must be sorted.\n"
               << "(Duplicate values are allowed.)";
            throw std::runtime_error(os.str());
        }
    }
    
    // Check that all temperatures are non-negative
    bool negative = false;
    
    for (Index i = 0; !negative && i < abs_t.nelem (); i++)
    {
        if (abs_t[i] < 0.)
            negative = true;
    }
    
    if (negative)
    {
        std::ostringstream os;
        os << "abs_t contains at least one negative temperature value.\n"
           << "This is not allowed.";
        throw std::runtime_error(os.str());
    }
    
    // We need a local copy of f_grid which is 1 element longer, because
    // we append a possible cutoff to it.
    // The initialization of this has to be inside the line loop!
    Vector f_local( nf + 1 );
    f_local[Range(0,nf)] = f_grid;
    
    // Voigt generally needs a different frequency grid. If we allocate
    // that in the outer loop, instead of in voigt, we don't have the
    // free store allocation at each lineshape call. Calculation is
    // still done in the voigt routine itself, this is just an auxillary
    // parameter, passed to lineshape. For selected lineshapes (e.g.,
    // Rosenkranz) it is used additionally to pass parameters needed in
    // the lineshape (e.g., overlap, ...). Consequently we have to
    // assure that aux has a dimension not less then the number of
    // parameters passed.
    Index ii = (nf+1 < 10) ? 10 : nf+1;
    Vector aux(ii);
    
    // Check that abs_p, abs_t, and abs_vmrs have consistent
    // dimensions. This could be a user error, so we throw a
    // runtime_error.
    
    if ( abs_t.nelem() != np )
    {
        std::ostringstream os;
        os << "Variable abs_t must have the same dimension as abs_p.\n"
           << "abs_t.nelem() = " << abs_t.nelem() << '\n'
           << "abs_p.nelem() = " << np;
        throw std::runtime_error(os.str());
    }
    
    // all_vmrs should have dimensions [nspecies, np]:
    
    if ( all_vmrs.ncols() != np )
    {
        std::ostringstream os;
        os << "Number of columns of all_vmrs must match abs_p.\n"
           << "all_vmrs.ncols() = " << all_vmrs.ncols() << '\n'
           << "abs_p.nelem() = " << np;
        throw std::runtime_error(os.str());
    }
    
    const Index nspecies = abs_species.nelem();
    
    if ( all_vmrs.nrows() != nspecies)
    {
        std::ostringstream os;
        os << "Number of rows of all_vmrs must match abs_species.\n"
           << "all_vmrs.nrows() = " << all_vmrs.nrows() << '\n'
           << "abs_species.nelem() = " << nspecies;
        throw std::runtime_error(os.str());
    }
    
    // Check that the dimension of xsec is indeed [f_grid.nelem(),
    // abs_p.nelem()]:
    if ( xsec_attenuation.nrows() != nf || xsec_attenuation.ncols() != np )
    {
        std::ostringstream os;
        os << "Variable xsec_attenuation must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
           << "[xsec_attenuation.nrows(),xsec_attenuation.ncols()] = [" << xsec_attenuation.nrows()
           << ", " << xsec_attenuation.ncols() << "]\n"
           << "f_grid.nelem() = " << nf << '\n'
           << "abs_p.nelem() = " << np;
        throw std::runtime_error(os.str());
    }
    
    // Check that the dimension of xsec source is [f_grid.nelem(), abs_p.nelem()]:
    bool calc_src;
    if ( xsec_source.nrows() != nf || xsec_source.ncols() != np )
    {
      if( xsec_source.nrows() != 0 || xsec_source.ncols() != 0 )
      {
        std::ostringstream os;
        os << "Variable xsec_source must have dimensions [f_grid.nelem(),abs_p.nelem()] or [0,0].\n"
           << "[xsec_source.nrows(),xsec_source.ncols()] = [" << xsec_source.nrows()
           << ", " << xsec_source.ncols() << "]\n"
           << "f_grid.nelem() = " << nf << '\n'
           << "abs_p.nelem() = " << np;
        throw std::runtime_error(os.str());
      }
      else 
      {
	calc_src = false;
      }
    }
    else
    {
      calc_src = true;
    }
    
    if ( xsec_phase.nrows() != 0 || xsec_phase.ncols() != 0 )
    {
        std::ostringstream os;
        os << "Variable xsec_phase must have dimensions [0,0] here.\n"
           << "[xsec_phase.nrows(),xsec_phase.ncols()] = [" << xsec_phase.nrows()
           << ", " << xsec_phase.ncols() << "]\n"
           << "f_grid.nelem() = " << nf << '\n'
           << "abs_p.nelem() = " << np;
        throw std::runtime_error(os.str());
    }
    
    // Find the location of all broadening species in abs_species. Set to -1 if
    // not found. The length of array broad_spec_locations is the number of allowed
    // broadening species (in ARTSCAT-4 Self, N2, O2, H2O, CO2, H2, He). The value
    // means:
    // -1 = not in abs_species
    // -2 = in abs_species, but should be ignored because it is identical to Self
    // N  = species is number N in abs_species
    ArrayOfIndex broad_spec_locations;
    find_broad_spec_locations(broad_spec_locations,
                              abs_species,
                              this_species);
    const Index h2o_index = find_first_species_tg( abs_species,
                                                   species_index_from_species_name("H2O") );
    
    String fail_msg;
    bool failed = false;
    
    // Loop all pressures:
    if (np)
#pragma omp parallel for                    \
if (!arts_omp_in_parallel()               \
&& np >= arts_omp_get_max_threads())  \
firstprivate(ls_attenuation, fac, f_local, aux)
        for ( Index i=0; i<np; ++i )
        {
            if (failed) continue;
            
            // holder when things can be empty
            Vector empty_vector;
            
            // Store input profile variables, this is perhaps slightly faster.
            const Numeric p_i       = abs_p[i];
            const Numeric t_i       = abs_t[i];
            ConstVectorView t_nlte_i  = calc_src ? abs_t_nlte(joker,i):empty_vector;
            const Numeric vmr_i     = all_vmrs(this_species,i);
            
            //out3 << "  p = " << p_i << " Pa\n";
            
            // Calculate total number density from pressure and temperature.
            // n = n0*T0/p0 * p/T or n = p/kB/t, ideal gas law
            //      const Numeric n = p_i / BOLTZMAN_CONST / t_i;
            // This is not needed anymore, since we now calculate true cross
            // sections, which do not contain the n.
            
            // For the pressure broadening, we also need the partial pressure:
            const Numeric p_partial = p_i * vmr_i;
            
            // Get handle on xsec for this pressure level i.
            // Watch out! This is output, we have to be careful not to
            // introduce race conditions when writing to it.
            VectorView xsec_i_attenuation = xsec_attenuation(Range(joker),i);
            VectorView xsec_i_source = calc_src ? xsec_source(Range(joker),i) : empty_vector;
            
            
            //       if (omp_in_parallel())
            //         cout << "omp_in_parallel: true\n";
            //       else
            //         cout << "omp_in_parallel: false\n";
            
            
            // Prepare a variable that can be used by the individual LBL
            // threads to add up absorption:
            Index n_lbl_threads;
            if (arts_omp_in_parallel())
            {
                // If we already are running parallel, then the LBL loop
                // will not be parallelized.
                n_lbl_threads = 1;
            }
            else
            {
                n_lbl_threads = arts_omp_get_max_threads();
            }
            Matrix xsec_accum_attenuation(n_lbl_threads, xsec_i_attenuation.nelem(), 0);
            Matrix xsec_accum_source(n_lbl_threads, xsec_i_attenuation.nelem());
            if(calc_src)
              xsec_accum_source=0.0;
            
            // Simple caching of partition function to avoid recalculating things.
            Numeric qt_cache=-1, qref_cache=-1;
            Index iso_cache = - 1;
            Numeric line_t_cache = - 1;
            
            ConstVectorView vmrs = all_vmrs(joker,i);
            
            // Loop all lines:
            if (nl)
            {
#pragma omp parallel for                   \
if (!arts_omp_in_parallel()               \
&& nl >= arts_omp_get_max_threads())  \
firstprivate(ls_attenuation, fac, f_local, aux, qt_cache, qref_cache, iso_cache, line_t_cache)
                for ( Index l=0; l< nl; ++l )
                {
                    // Skip remaining iterations if an error occurred
                    if (failed) continue;
                    
                    //           if (omp_in_parallel())
                    //             cout << "LBL: omp_in_parallel: true\n";
                    //           else
                    //             cout << "LBL: omp_in_parallel: false\n";
                    
                    
                    // The try block here is necessary to correctly handle
                    // exceptions inside the parallel region.
                    try
                    {
                        const LineRecord& l_l = abs_lines[l];
                        
                        // Prepare pressure broadening parameters
                        Numeric gamma_0,gamma_2,eta,df_0,df_2,f_VC;
                        l_l.PressureBroadening().GetPressureBroadeningParams(gamma_0,gamma_2,eta,df_0,df_2,f_VC,
                                                                             l_l.Ti0()/t_i,p_i,
                                                                             p_partial,this_species,h2o_index,
                                                                             broad_spec_locations,
                                                                             vmrs,verbosity);
                        
                        // Check the chache is the tempearture of the line and the isotope is the same to avoid recalculating the partition sum
                        if(iso_cache!=l_l.Isotopologue() || line_t_cache != l_l.Ti0())
                        {
                          iso_cache = l_l.Isotopologue();
                          line_t_cache = l_l.Ti0();
                          qref_cache=-1;// no need to reset qt since it is done internally.
                        }
                        
                        // Prepare line strength scaling
                        Numeric partition_ratio, K1, K2, abs_nlte_ratio, src_nlte_ratio;
                        GetLineScalingData(qt_cache,
                                           qref_cache,
                                           partition_ratio, 
                                           K1, 
                                           K2,
                                           abs_nlte_ratio, 
                                           src_nlte_ratio, 
                                           partition_functions.getParamType(l_l.Species(), l_l.Isotopologue()),
                                           partition_functions.getParam(l_l.Species(), l_l.Isotopologue()),
                                           t_i,
                                           l_l.Ti0(),
                                           l_l.F(),
                                           l_l.Elow(),
                                           calc_src,
                                           l_l.Evlow(),
                                           l_l.Evupp(),
                                           l_l.NLTELowerIndex(),
                                           l_l.NLTEUpperIndex(),
                                           t_nlte_i);
                       
                        // Dopple broadening
                        const Numeric sigma = l_l.F() * doppler_const *sqrt( t_i / l_l.IsotopologueData().Mass());
                        
                        const bool calc_partials=false;
                        Vector da_dF, dp_dF, da_dP, dp_dP;
                        
                        Range this_f_range(0,0);
                        
                        // Calculate line cross section
                        xsec_single_line(// OUTPUT
                                         xsec_accum_attenuation(arts_omp_get_thread_num(),joker),
                                         xsec_accum_source(arts_omp_get_thread_num(),joker),
                                         empty_vector,
                                         // HELPER:
                                         ls_attenuation,
                                         ls_phase_dummy,
                                         da_dF,
                                         dp_dF,
                                         da_dP,
                                         dp_dP,
                                         this_f_range,
                                         fac,
                                         aux,
                                         // FREQUENCY
                                         f_local,
                                         f_grid,
                                         nf,
                                         cutoff,
                                         l_l.F(),
                                         // LINE STRENGTH
                                         l_l.I0(),
                                         partition_ratio,
                                         K1*K2,
                                         abs_nlte_ratio,
                                         src_nlte_ratio,
                                         isotopologue_ratios.getParam(l_l.Species(),l_l.Isotopologue())[0].data[0],
                                         // ATMOSPHERIC TEMPERATURE
                                         t_i,
                                         // LINE SHAPE
                                         ind_ls,
                                         ind_lsn,
                                         // LINE BROADENING
                                         gamma_0,
                                         gamma_2,
                                         eta,
                                         df_0,
                                         df_2,
                                         sigma,
                                         f_VC,
                                         // LINE MIXING
                                         0,
                                         0,
                                         0,
                                         // FEATURE FLAGS
                                         cut,
                                         calc_phase,
                                         calc_partials,
                                         calc_src);
                        
                    } // end of try block
                    catch (runtime_error e)
                    {
                        #pragma omp critical (xsec_species_fail)
                        { fail_msg = e.what(); failed = true; }
                    }
                    
                } // end of parallel LBL loop
            }
                
            // Bail out if an error occurred in the LBL loop
            if (failed) continue;
            
            // Now we just have to add up all the rows of xsec_accum:
            for (Index j=0; j<xsec_accum_attenuation.nrows(); ++j)
            {
                xsec_i_attenuation += xsec_accum_attenuation(j, Range(joker));
                if(calc_src)
                    {
                        xsec_i_source += xsec_accum_source(j, Range(joker));
                    }
            }
    } // end of parallel pressure loop
    
    if (failed) throw std::runtime_error("Run-time error in function: xsec_species\n" + fail_msg);
}


/** 

   Calculate line absorption cross sections for one line at one layer
   and accumulates to the total phase and attenuation change.
   
   No dependency on LineRecord to increase speed for wrapper applications.
   
   \retval xsec_accum_attenuation   Cross section of one tag group. This is now the
                                    true absorption cross section in units of m^2.
                                    It has inputs of all previously calculated lines.
   \retval xsec_accum_source        Deviation of line cross section from LTE cross-section
   \retval xsec_accum_phase         Cross section of one tag group. This is now the
                                    true dispersion cross section in units of m^2.
                                    It has inputs of all previously calculated lines.
   \retval attenuation              Input only to increase speed.  Holds attenuation internally.
   \retval phase                    Input only to increase speed.  Holds phase internally.
   \retval fac                      Input only to increase speed.  Holds lineshape factor internally.
   \retval aux                      Input only to increase speed.  Holds f_grid factor internally.
   \retval f_local                  Input only to increase speed.  Holds f_grid internally.
   \param f_grid,                   Frequency grid
   \param nf,                       Number of frequencies to calculate
   \param cutoff,                   Lineshape cutoff.
   \param F0,                       Line center
   \param intensity,                Line intensity
   \param part_fct_ratio            Ratio of partition sums to atmospheric temperature
   \param boltzmann_ratio           Ratio of Boltzmann statistics to atmospheric temperature
   \param abs_nlte_ratio            Ratio of absorption intensity to LTE
   \param src_nlte_ratio            Ratio of emission intesity to LTE
   \param Isotopologue_Ratio,       Ratio of the isotopologue in the atmosphere
   \param temperature,              Atmospheric temperature
   \param ind_ls,                   Index to lineshape function.
   \param ind_lsn,                  Index to lineshape norm.
   \param gamma_0,                  Line pressure broadening
   \param gamma_2,                  Speed-dependent line pressure broadening
   \param eta,                      Correlation of pressure broadening parameters
   \param df_0,                     Line pressure shift
   \param df_2,                     Speed-dependent line pressure shift
   \param sigma,                    Doppler broadening
   \param f_VC,                     Collisional frequency limit
   \param LM_DF,                    Line mixing frequency shift
   \param LM_Y,                     Line mixing dispersion dependency
   \param LM_G,                     Line mixing added attenuation
   \param cut,                      Is cutoff applied?
   \param calc_phase,               Is dispersion calculated?
   \param calc_src,                 Is source calculated?
 
   \author Stefan Buehler and Axel von Engeln
   \date   2001-01-11 
   
   Changed from pseudo cross sections to true cross sections
   
   \author Stefan Buehler 
   \date   2007-08-08 
   
   Adapted to new Perrin line parameters, treating broadening by different 
   gases explicitly
   
   \author Stefan Buehler
   \date   2012-09-03
   
   Removed LineRecord dependency.
   
   \author Richard Larsson
   \date   2014-10-29
   
*/
void xsec_single_line(// Output:
                      VectorView xsec_accum_attenuation, 
                      VectorView xsec_accum_source, 
                      VectorView xsec_accum_phase, 
                      // Helper variables
                      Vector& attenuation, 
                      Vector& phase,
                      Vector& dFa_dx,
                      Vector& dFb_dx,
                      Vector& dFa_dy,
                      Vector& dFb_dy,
                      Range&  this_f_range,
                      Vector& fac, 
                      Vector& aux, 
                      // Frequency grid:
                      Vector& f_local, 
                      const Vector& f_grid, 
                      const Index nf, 
                      const Numeric cutoff,
                      Numeric F0, 
                      // Line strength:
                      Numeric intensity, 
                      const Numeric part_fct_ratio,  
                      const Numeric boltzmann_ratio,
                      const Numeric abs_nlte_ratio,
                      const Numeric src_nlte_ratio,
                      const Numeric Isotopologue_Ratio,
                      // Atmospheric state
                      const Numeric temperature, 
                      // Line shape:
                      const Index ind_ls, 
                      const Index ind_lsn,  
                      // Line broadening:
                      const Numeric gamma_0,
                      const Numeric gamma_2,
                      const Numeric eta,
                      const Numeric df_0,
                      const Numeric df_2,
                      const Numeric sigma,
                      const Numeric f_VC,
                      // Line mixing
                      const Numeric LM_DF,
                      const Numeric LM_Y, 
                      const Numeric LM_G,
                      // Feature flags
                      const bool calc_cut, 
                      const bool calc_phase,
                      const bool calc_partials,
                      const bool calc_src)
{
    
    // Copy f_grid to the beginning of f_local. There is one
    // element left at the end of f_local.
    // THIS HAS TO BE INSIDE THE LINE LOOP, BECAUSE THE CUTOFF
    // FREQUENCY IS ALWAYS PUT IN A DIFFERENT PLACE!
    
                                              
    // This will hold the actual number of frequencies to add to
    // xsec later on:
    Index nfl = nf;
    
    // This will hold the actual number of frequencies for the
    // call to the lineshape functions later on:
    Index nfls = nf;
        
    // intensity at temperature
    // (calculate the line intensity according to the standard
    // expression given in HITRAN)
    intensity *= part_fct_ratio * boltzmann_ratio;
    
    // Apply pressure shift:
    //
    // 2017-10-15: Stuart Fox had reported the bug that for mirror
    // lines (negative F0) the shift was in the wrong direction. The
    // if statement below makes sure that the shift goes in the other
    // direction for negative frequency mirror lines, so that they
    // really end up at the negative of the frequency of the original
    // line.
    if (F0 >= 0)
      F0 += ( df_0 + LM_DF );
    else
      F0 -= ( df_0 + LM_DF );
    // FIXME:  This is probably not compatible with HTP line shape
    
    // Indices pointing at begin/end frequencies of f_grid or at
    // the elements that have to be calculated in case of cutoff
    Index i_f_min = 0;
    Index i_f_max = nf-1;
    
    // cutoff ?
    if ( calc_cut )
    {
        // Check whether we have elements in ls that can be
        // ignored at lower frequencies of f_grid.
        //
        // Loop through all frequencies, finding min value and
        // set all values to zero on that way.
        while ( i_f_min < nf && (F0 - cutoff) > f_grid[i_f_min] )
        {
            //              ls[i_f_min] = 0;
            ++i_f_min;
        }
        
        
        // Check whether we have elements in ls that can be
        // ignored at higher frequencies of f_grid.
        //
        // Loop through all frequencies, finding max value and
        // set all values to zero on that way.
        while ( i_f_max >= 0 && (F0 + cutoff) < f_grid[i_f_max] )
        {
            //              ls[i_f_max] = 0;
            --i_f_max;
        }
        
        // Append the cutoff frequency to f_local:
        ++i_f_max;
        f_local[i_f_max] = F0 + cutoff;
        
        // Number of frequencies to calculate:
        nfls = i_f_max - i_f_min + 1; // Add one because indices
                                      // are pointing to first and
                                      // last valid element. This
                                      // is for the lineshape
                                      // calls.
        nfl = nfls -1;              // This is for xsec.
    }
    else
    {
        // Nothing to do here. Note that nfl and nfls are both still set to nf.
    }
    
    //          cout << "nf, nfl, nfls = " << nf << ", " << nfl << ", " << nfls << ".\n";
    
    // Maybe there are no frequencies left to compute?  Note that
    // the number that counts here is nfl, since only these are
    // the `real' frequencies, for which xsec is changed. nfls
    // will always be at least one, because it contains the cutoff.
    if ( nfl > 0 )
    {
        this_f_range = Range(i_f_min,nfl);
        Range that_f_range(i_f_min,nfls);
        
        // Calculate the line shape:
        global_data::lineshape_data[ind_ls].Function()(attenuation, phase,
                                                       dFa_dx, dFb_dx, dFa_dy, dFb_dy, //partial derivatives
                                                       aux, F0, gamma_0, gamma_2, eta, 0.0, // FIXME: H-T ls...
                                                       df_2, sigma, f_VC, f_local[that_f_range],
                                                       calc_phase, calc_partials);
        
        // Calculate the chosen normalization factor:
        global_data::lineshape_norm_data[ind_lsn].Function()(fac, F0,
                                                             f_local[that_f_range],
                                                             temperature);
        
        // Reset f_local for next loop through now that cutoff is calculated
        if( i_f_max < nf )
            f_local[i_f_max] = f_grid[i_f_max];
        
        // Attenuation output
        VectorView this_xsec_attenuation      = xsec_accum_attenuation[this_f_range];
        
        // Optional outputs
        Vector empty_vector;
        VectorView this_xsec_source           = calc_src ? xsec_accum_source[this_f_range] : empty_vector;
        VectorView this_xsec_phase            = calc_phase ? xsec_accum_phase[this_f_range] : empty_vector;
        
        // Store cutoff values
        const Numeric 
        cutoff_attenuation = calc_cut?attenuation[nfls-1]:0.0, 
        cutoff_phase       = calc_phase?calc_cut?phase[nfls-1]:0.0:0.0,
        cutoff_dFa_dx      = calc_partials?calc_cut?dFa_dx[nfls-1]:0.0:0.0, 
        cutoff_dFb_dx      = calc_phase?calc_partials?calc_cut?dFb_dx[nfls-1]:0.0:0.0:0.0, 
        cutoff_dFa_dy      = calc_partials?calc_cut?dFa_dy[nfls-1]:0.0:0.0, 
        cutoff_dFb_dy      = calc_phase?calc_partials?calc_cut?dFb_dy[nfls-1]:0.0:0.0:0.0;
        
        // If one of these is non-zero, then there are line mixing calculations required
        const bool calc_LM = LM_G!=0||LM_Y!=0;
        
        // Line strength is going to be stored in attenuation to ensure the right numeric value.
        const Numeric str_line = Isotopologue_Ratio*intensity*(calc_src?abs_nlte_ratio:1.0);
        
        // To make this readable, we loop over the frequencies
        for(Index jj =0; jj<this_xsec_attenuation.nelem();jj++)
        {
            const Numeric str_scale = fac[jj]*str_line;
            if(jj<=nfl)
            {   
                // This ensures that output attenuation is the physically correct attenuation
                if(calc_cut)
                    attenuation[jj]-=cutoff_attenuation;
                attenuation[jj]*=str_scale;
                
                // This ensures that output phase is the physically correct phase
                if(calc_phase||calc_LM)
                {
                    if(calc_cut)
                        phase[jj]-=cutoff_phase;
                    phase[jj]*=str_scale;
                }
                
                if(calc_partials)
                {
                    if(calc_cut)
                    {
                        dFa_dx[jj]-=cutoff_dFa_dx;
                        dFa_dy[jj]-=cutoff_dFa_dy;
                        if(calc_phase)
                        {
                            dFb_dx[jj]-=cutoff_dFb_dx;
                            dFb_dy[jj]-=cutoff_dFb_dy;
                        }
                    }
                    dFa_dx[jj]*=str_scale;
                    dFa_dy[jj]*=str_scale;
                    if(calc_phase)
                    {
                        dFb_dx[jj]*=str_scale;
                        dFb_dy[jj]*=str_scale;
                    }
                }
                
                // Attenuation is by tradition added to total attenuation from here
                if(calc_LM&&calc_src)
                {
                    const Numeric tmp =(1.+LM_G)*attenuation[jj]+LM_Y*phase[jj];
                    this_xsec_source[jj] += ( src_nlte_ratio/abs_nlte_ratio - 1.0 ) * tmp;
                    this_xsec_attenuation[jj] += tmp;
                }
                else if(calc_src)
                {
                    this_xsec_source[jj] += ( src_nlte_ratio/abs_nlte_ratio - 1.0 ) * attenuation[jj];
                    this_xsec_attenuation[jj] += attenuation[jj];
                }
                else if(calc_LM)
                    this_xsec_attenuation[jj] += (1.+LM_G)*attenuation[jj]+LM_Y*phase[jj];
                else
                    this_xsec_attenuation[jj] += attenuation[jj];
                
                // Phase follows the attenuation practice
                if(calc_phase&&calc_LM)
                    this_xsec_phase[jj] += (1.+LM_G)*phase[jj] - LM_Y*attenuation[jj];
                else if(calc_phase)
                    this_xsec_phase[jj] += phase[jj];
                
                // Partials (Assuming F(v,v0) = C*F'(v,v0), where F' is from our lineshape functions, then the below is 
                // only the C*dF'/dt part.  The dC/dt*F' part remains to be calculated.  The CF is output as 
                // attenuation and phase variables [NOTE: attenuation and phase are without line mixing corrections!])
                if(calc_partials&&calc_LM)
                {
                    const Numeric orig_dFa_dx=dFa_dx[jj];
                    const Numeric orig_dFa_dy=dFa_dy[jj];
                    
                    dFa_dx[jj] = (1.+LM_G)*orig_dFa_dx + LM_Y*dFb_dx[jj];
                    dFb_dx[jj] = (1.+LM_G)*dFb_dx[jj]  - LM_Y*orig_dFa_dx;
                    dFa_dy[jj] = (1.+LM_G)*orig_dFa_dy + LM_Y*dFb_dy[jj];
                    dFb_dy[jj] = (1.+LM_G)*dFb_dy[jj]  - LM_Y*orig_dFa_dy;
                }
            }
        }  
    }
}



/** A little helper function to convert energy from units of
    wavenumber (cm^-1) to Joule (J). 

    This is used when reading HITRAN or JPL catalogue files, which
    have the lower state energy in cm^-1.

    \return Energy in J.
    \param  e Energy in cm^-1.

    \author Stefan Buehler
    \date   2001-06-26 */
Numeric wavenumber_to_joule(Numeric e)
{
  // Planck constant [Js]
  extern const Numeric PLANCK_CONST;

  // Speed of light [m/s]
  extern const Numeric SPEED_OF_LIGHT;

  // Constant to convert lower state energy from cm^-1 to J
  const Numeric lower_energy_const = PLANCK_CONST * SPEED_OF_LIGHT * 1E2;

  return e*lower_energy_const; 
}


//======================================================================
//             Functions related to species
//======================================================================

//! Return species index for given species name.
/*! 
  This is useful in connection with other functions that need a species
  index.

  \see find_first_species_tg.

  \param name Species name.

  \return Species index, -1 means not found.

  \author Stefan Buehler
  \date   2003-01-13
*/
Index species_index_from_species_name( String name )
{
  using global_data::SpeciesMap;

  // For the return value:
  Index mspecies;

  // Trim whitespace
  name.trim();

  //  cout << "name / def = " << name << " / " << def << endl;

  // Look for species name in species map:
  map<String, Index>::const_iterator mi = SpeciesMap.find(name);
  if ( mi != SpeciesMap.end() )
    {
      // Ok, we've found the species. Set mspecies.
      mspecies = mi->second;
    }
  else
    {
      // The species does not exist!
      mspecies = -1;
    }

  return mspecies;
}


//! Return species name for given species index.
/*!
 This is useful in connection with other functions that use a species
 index.
 
 Does an assertion that the index really corresponds to a species.
 
 \param spec_ind Species index.
 
 \return Species name
 
 \author Stefan Buehler
 \date   2013-01-04
 */
String species_name_from_species_index( const Index spec_ind )
{
    // Species lookup data:
    using global_data::species_data;

    // Assert that spec_ind is inside species data. (This is an assertion,
    // because species indices should never be user input, but set by the
    // program automatically, based on species names.)
    assert( spec_ind>=0 );
    assert( spec_ind<species_data.nelem() );
    
    // A reference to the relevant record of the species data:
    const  SpeciesRecord& spr = species_data[spec_ind];
    
    return spr.Name();
}


//======================================================================
//        Functions to convert the accuracy index to ARTS units
//======================================================================

// ********* for HITRAN database *************
//convert HITRAN index for line position accuracy to ARTS
//units (Hz).

void convHitranIERF(     
                    Numeric&     mdf,
              const Index&       df 
                    )
{
  switch ( df )
    {
    case 0:
      { 
        mdf = -1;
        break; 
      }
    case 1:
      {
        mdf = 30*1E9;
        break; 
      }
    case 2:
      {
        mdf = 3*1E9;
        break; 
      }
    case 3:
      {
        mdf = 300*1E6;
        break; 
      }
    case 4:
      {
        mdf = 30*1E6;
        break; 
      }
    case 5:
      {
        mdf = 3*1E6;
        break; 
      }
    case 6:
      {
        mdf = 0.3*1E6;
        break; 
      }
    }
}
                    
//convert HITRAN index for intensity and halfwidth accuracy to ARTS
//units (relative difference).
void convHitranIERSH(     
                    Numeric&     mdh,
              const Index&       dh 
                    )
{
  switch ( dh )
    {
    case 0:
      { 
        mdh = -1;
        break; 
      }
    case 1:
      {
        mdh = -1;
        break; 
      }
    case 2:
      {
        mdh = -1;
        break; 
      }
    case 3:
      {
        mdh = 30;
        break; 
      }
    case 4:
      {
        mdh = 20;
        break; 
      }
    case 5:
      {
        mdh = 10;
        break; 
      }
    case 6:
      {
        mdh =5;
        break; 
      }
    case 7:
      {
        mdh =2;
        break; 
      }
    case 8:
      {
        mdh =1;
        break; 
      }
    }
  mdh=mdh/100;              
}

// ********* for MYTRAN database ************* 
//convert MYTRAN index for intensity and halfwidth accuracy to ARTS
//units (relative difference).
void convMytranIER(     
                    Numeric&     mdh,
              const Index  &      dh 
                    )
{
  switch ( dh )
    {
    case 0:
      { 
        mdh = 200;
        break; 
      }
    case 1:
      {
        mdh = 100;
        break; 
      }
    case 2:
      {
        mdh = 50;
        break; 
      }
    case 3:
      {
        mdh = 30;
        break; 
      }
    case 4:
      {
        mdh = 20;
        break; 
      }
    case 5:
      {
        mdh = 10;
        break; 
      }
    case 6:
      {
        mdh =5;
        break; 
      }
    case 7:
      {
        mdh = 2;
        break; 
      }
    case 8:
      {
        mdh = 1;
        break; 
      }
    case 9:
      {
        mdh = 0.5;
        break; 
      }
    }
  mdh=mdh/100;              
}

ostream& operator<< (ostream &os, const LineshapeSpec& lsspec)
{
    os << "LineshapeSpec Index: " << lsspec.Ind_ls() << ", NormIndex: " << lsspec.Ind_lsn()
    << ", Cutoff: " << lsspec.Cutoff() << endl;

    return os;
}


/**
 *  
 *  This will work as a wrapper for linemixing when abs_species contain relevant data.
 *  The funciton will only pass on arguments to xsec_species if there is no linemixing.
 *  
 *  \retval xsec_attenuation    Cross section of one tag group. This is now the
 *                              true attenuation cross section in units of m^2.
 *  \retval xsec_source         Cross section of one tag group. This is now the
 *                              true source cross section in units of m^2.
 *  \retval xsec_phase          Cross section of one tag group. This is now the
 *                              true phase cross section in units of m^2.
 *  \param f_grid               Frequency grid.
 *  \param abs_p                Pressure grid.
 *  \param abs_t                Temperatures associated with abs_p.
 *  \param all_vmrs             Gas volume mixing ratios [nspecies, np].
 *  \param abs_species          Species tags for all species.
 *  \param this_species         Index of the current species in abs_species.
 *  \param abs_lines            The spectroscopic line list.
 *  \param ind_ls               Index to lineshape function.
 *  \param ind_lsn              Index to lineshape norm.
 *  \param cutoff               Lineshape cutoff.
 *  \param isotopologue_ratios  Isotopologue ratios.
 * 
 *  \author Richard Larsson
 *  \date   2013-04-24
 * 
 */
void xsec_species_line_mixing_wrapper(  MatrixView               xsec_attenuation,
                                        MatrixView               xsec_source,
                                        MatrixView               xsec_phase,
                                        ArrayOfMatrix&           partial_xsec_attenuation,
                                        ArrayOfMatrix&           partial_xsec_source,
                                        ArrayOfMatrix&           partial_xsec_phase,
                                        const PropmatPartialsData& flag_partials,
                                        ConstVectorView          f_grid,
                                        ConstVectorView          abs_p,
                                        ConstVectorView          abs_t,
                                        ConstMatrixView          abs_t_nlte,
                                        ConstMatrixView          all_vmrs,
                                        const ArrayOfArrayOfSpeciesTag& abs_species,
                                        const Index              this_species,
                                        const ArrayOfLineRecord& abs_lines,
                                        const Vector&            Z_DF,
                                        const Numeric            H_magntitude_Zeeman,
                                        const Index              ind_ls,
                                        const Index              ind_lsn,
                                        const Numeric            lm_p_lim,
                                        const Numeric            cutoff,
                                        const SpeciesAuxData&    isotopologue_ratios,
                                        const SpeciesAuxData&    partition_functions,
                                        const Verbosity&         verbosity )
{
    
    // Optional paths through the code...
    const bool cut = (cutoff != -1) ? true : false;
    const bool calc_partials = flag_partials.nelem();
    const bool calc_partials_phase = partial_xsec_phase.nelem()==partial_xsec_attenuation.nelem();
    const bool do_zeeman = !Z_DF.empty();
    const bool do_lm     = abs_species[this_species][0].LineMixing() != SpeciesTag::LINE_MIXING_OFF;
    const bool no_ls_partials = flag_partials.supportsLBLwithoutPhase();
    
    const bool we_need_phase = do_zeeman || do_lm || calc_partials_phase;
    const bool we_need_partials =calc_partials&&!no_ls_partials;
    
    // Must have the phase for two of the options
    using global_data::lineshape_data;
    if (!lineshape_data[ind_ls].Phase() && we_need_phase )
    {
        std::ostringstream os;
        os <<  "This is an error message. You are using " << lineshape_data[ind_ls].Name() <<
        ".\n"<<"This line shape does not include phase in its calculations and\nis therefore invalid for " <<
        "line mixing, Zeeman effect, and certain partial derivatives.\nYou are using one of these or you should not see this error.\n";
        throw std::runtime_error(os.str());
    }
    
    if(we_need_partials&&!lineshape_data[ind_ls].Partials())
    {
        std::ostringstream os;
        os <<  "This is an error message. You are using " << lineshape_data[ind_ls].Name() <<".\n"
           <<  "Your selected *jacobian_quantities* requires that the line shape returns partial\n"
           <<  "derivatives.";
           
        throw std::runtime_error(os.str());
    }
    
    extern const Numeric BOLTZMAN_CONST;
    extern const Numeric AVOGADROS_NUMB;
    extern const Numeric SPEED_OF_LIGHT;
    static const Numeric doppler_const = sqrt(2.0 * BOLTZMAN_CONST *
                                              AVOGADROS_NUMB) / SPEED_OF_LIGHT;
    
    bool precalc_zeeman;
    if(Z_DF.nelem()==0)
        precalc_zeeman=false;
    else if(Z_DF.nelem()==abs_lines.nelem())
        precalc_zeeman=true;
    else
    {
        std::ostringstream os;
        os <<  "This is an error message. You have sent a strange Zeeman frequency shift vector to"<<
        "the xsec_species_line_mixing_wrapper function.  This causes the failure.\n";
        throw std::runtime_error(os.str());
    }
    
    // Check that the frequency grid is sorted in the case of lineshape
    // with cutoff. Duplicate frequency values are allowed.
    if (cut)
    {
        if ( ! is_sorted( f_grid ) )
        {
            std::ostringstream os;
            os << "If you use a lineshape function with cutoff, your\n"
            << "frequency grid *f_grid* must be sorted.\n"
            << "(Duplicate values are allowed.)";
            throw std::runtime_error(os.str());
        }
    }
    
    // Check that all temperatures are non-negative
    bool negative = false;
    
    for (Index i = 0; !negative && i < abs_t.nelem (); i++)
    {
        if (abs_t[i] < 0.)
            negative = true;
    }
    
    if (negative)
    {
        std::ostringstream os;
        os << "abs_t contains at least one negative temperature value.\n"
        << "This is not allowed.";
        throw std::runtime_error(os.str());
    }
    
    if ( abs_t.nelem() != abs_p.nelem() )
    {
        std::ostringstream os;
        os << "Variable abs_t must have the same dimension as abs_p.\n"
        << "abs_t.nelem() = " << abs_t.nelem() << '\n'
        << "abs_p.nelem() = " << abs_p.nelem();
        throw std::runtime_error(os.str());
    }
    
    // all_vmrs should have dimensions [nspecies, np]:
    
    if ( all_vmrs.ncols() != abs_p.nelem() )
    {
        std::ostringstream os;
        os << "Number of columns of all_vmrs must match abs_p.\n"
        << "all_vmrs.ncols() = " << all_vmrs.ncols() << '\n'
        << "abs_p.nelem() = " << abs_p.nelem();
        throw std::runtime_error(os.str());
    }
    
    // Check that the dimension of xsec is indeed [f_grid.nelem(),
    // abs_p.nelem()]:
    if ( xsec_attenuation.nrows() != f_grid.nelem() || xsec_attenuation.ncols() != abs_p.nelem() )
    {
        std::ostringstream os;
        os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
        << "[xsec_attenuation.nrows(),xsec_attenuation.ncols()] = [" << xsec_attenuation.nrows()
        << ", " << xsec_attenuation.ncols() << "]\n"
        << "f_grid.nelem() = " << f_grid.nelem() << '\n'
        << "abs_p.nelem() = " << abs_p.nelem();
        throw std::runtime_error(os.str());
    }
    
    // Check that the dimension of xsec source is [f_grid.nelem(), abs_p.nelem()]:
    bool calc_src;
    if ( xsec_source.nrows() != f_grid.nelem() || xsec_source.ncols() != abs_p.nelem() )
    {
      if( xsec_source.nrows() != 0 || xsec_source.ncols() != 0 )
      {
        std::ostringstream os;
        os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()] or [0,0].\n"
           << "[xsec_source.nrows(),xsec_source.ncols()] = [" << xsec_source.nrows()
           << ", " << xsec_source.ncols() << "]\n"
           << "f_grid.nelem() = " << f_grid.nelem() << '\n'
           << "abs_p.nelem() = " << abs_p.nelem();
        throw std::runtime_error(os.str());
      }
      else 
      {
	calc_src = 0;
      }
    }
    else
    {
      calc_src = 1;
    }
    
    if ( xsec_phase.nrows() != f_grid.nelem() || xsec_phase.ncols() != abs_p.nelem() )
    {
        std::ostringstream os;
        os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
        << "[xsec_phase.nrows(),xsec_phase.ncols()] = [" << xsec_phase.nrows()
        << ", " << xsec_phase.ncols() << "]\n"
        << "f_grid.nelem() = " << f_grid.nelem() << '\n'
        << "abs_p.nelem() = " << abs_p.nelem();
        throw std::runtime_error(os.str());
    }
    
    //Helper var 
    Vector attenuation(f_grid.nelem()+1),
    phase(f_grid.nelem()+1),
    aux(f_grid.nelem()+1),
    fac(f_grid.nelem()+1),
    f_local(f_grid.nelem()+1);
    f_local[Range(0,f_grid.nelem())]=f_grid;
    
    // Set this to zero in the case of VoigtKuntz6 or other lines hapes that work in special circumstances
    if (!we_need_phase&&we_need_partials)
        phase = 0.0;
    
    // Broadening species
    ArrayOfIndex broad_spec_locations;
    find_broad_spec_locations(broad_spec_locations,
                              abs_species,
                              this_species);
    // Water index
    const Index h2o_index = find_first_species_tg( abs_species,
                           species_index_from_species_name("H2O") );
    
#pragma omp parallel for        \
if (!arts_omp_in_parallel())    \
firstprivate(attenuation, phase, fac, f_local, aux)
    for(Index jj=0; jj<abs_p.nelem();jj++)
    {
        Vector empty_vector;
        const Numeric& t=abs_t[jj];
        const Numeric& p=abs_p[jj];
        ConstVectorView t_nlte = calc_src?abs_t_nlte(joker, jj):empty_vector;
        ConstVectorView vmrs = all_vmrs(joker,jj);
        
        
        // Setup for calculating the partial derivatives
        Vector dFa_dx, dFb_dx, dFa_dy, dFb_dy;
        if(calc_partials&&!no_ls_partials) //Only size them if partials are wanted
        {
            dFa_dx.resize(f_grid.nelem()+1);
            dFb_dx.resize(f_grid.nelem()+1);
            dFa_dy.resize(f_grid.nelem()+1); 
            dFb_dy.resize(f_grid.nelem()+1);
        }
        
        // Simple caching of partition function to avoid recalculating things.
        Numeric qt_cache=-1, qref_cache=-1;
        Index iso_cache=-1;
        Numeric line_t_cache=-1,  dqt_dT_cache=-1.;
        
        const Numeric p_partial = p * vmrs[this_species];
        
        for(Index ii=0; ii<abs_lines.nelem();ii++)
        {
            
            // These lines should be ignored by user request
            if( LineMixingData::LM_BYBAND == abs_lines[ii].LineMixing().Type() )
            {
                continue;
            }
            
            if(calc_src)
            {
              if(abs_lines[ii].GetLinePopulationType() not_eq LinePopulationType::ByLTE and
                 abs_lines[ii].GetLinePopulationType() not_eq LinePopulationType::ByVibrationalTemperatures)
              {
                throw std::runtime_error("Bad data seen in xsec_species_line_mixing_wrapper.  Please use more modern function.");
              }
            }
            
            // Pressure broadening parameters
            // Prepare pressure broadening parameters
            Numeric gamma_0,gamma_2,eta,df_0,df_2,f_VC;
            abs_lines[ii].PressureBroadening().GetPressureBroadeningParams(gamma_0,gamma_2,eta,df_0,df_2,f_VC,
                                                                abs_lines[ii].Ti0()/t,p,
                                                                p_partial,this_species,h2o_index,
                                                                broad_spec_locations,
                                                                vmrs,verbosity);
            
            // Line mixing parameters
            Numeric Y=0,DV=0,G=0; // Set to zero since case with 0 exist
            abs_lines[ii].LineMixing().GetLineMixingParams( Y,  G,  DV,  t, p, lm_p_lim, 1);
            
            // Check the cache is the temperature of the line and the isotope is the same to avoid recalculating the partition sum
            if(iso_cache!=abs_lines[ii].Isotopologue() || line_t_cache != abs_lines[ii].Ti0())
            {
              iso_cache = abs_lines[ii].Isotopologue();
              line_t_cache = abs_lines[ii].Ti0();
              
              qt_cache  =-1;
              qref_cache=-1;
              dqt_dT_cache=-1.;
            }
            
            // Prepare line strength scaling
            Numeric partition_ratio, K1, K2, abs_nlte_ratio, src_nlte_ratio;
            GetLineScalingData( qt_cache,
                                qref_cache,
                                partition_ratio, 
                                K1, 
                                K2,
                                abs_nlte_ratio, 
                                src_nlte_ratio, 
                                partition_functions.getParamType(abs_lines[ii].Species(), abs_lines[ii].Isotopologue()),
                                partition_functions.getParam(abs_lines[ii].Species(), abs_lines[ii].Isotopologue()),
                                t,
                                abs_lines[ii].Ti0(),
                                abs_lines[ii].F(),
                                abs_lines[ii].Elow(),
                                calc_src,
                                abs_lines[ii].Evlow(),
                                abs_lines[ii].Evupp(),
                                abs_lines[ii].NLTELowerIndex(),
                                abs_lines[ii].NLTEUpperIndex(),
                                t_nlte);
            
            // Doppler broadening
            const Numeric sigma = abs_lines[ii].F() * doppler_const *sqrt( t / abs_lines[ii].IsotopologueData().Mass());
            
            Range this_f_range(0,0);
            
            // Time to calculate the line cross section
            // Still an ugly case with non-resonant line near 0 frequency, since this is actually a band...
            if( LineMixingData::LM_LBLRTM_O2NonResonant != abs_lines[ii].LineMixing().Type() )
            {
                xsec_single_line(   // OUTPUT   
                                    xsec_attenuation(joker,jj), calc_src?xsec_source(joker,jj):
                                    empty_vector,xsec_phase(joker,jj), 
                                    // HELPER
                                    attenuation, phase, dFa_dx, dFb_dx, dFa_dy, dFb_dy, this_f_range, fac, aux, 
                                    // FREQUENCY
                                    f_local,  f_grid,  f_grid.nelem(),  cutoff,
                                    abs_lines[ii].F()+(precalc_zeeman?(Z_DF[ii]*H_magntitude_Zeeman):0), // Since vector is 0-length if no Zeeman pre-calculations
                                    // LINE STRENGTH
                                    abs_lines[ii].I0(), partition_ratio, K1*K2, abs_nlte_ratio, src_nlte_ratio,
                                    isotopologue_ratios.getParam(abs_lines[ii].Species(),
                                                                 abs_lines[ii].Isotopologue())[0].data[0],
                                    // ATMOSPHERIC TEMPERATURE
                                    t,
                                    // LINE SHAPE
                                    ind_ls, ind_lsn,
                                    // LINE BROADENING
                                    gamma_0, gamma_2, eta, df_0, df_2, sigma, f_VC,
                                    // LINE MIXING
                                    DV, Y,  G, 
                                    // FEATURE FLAGS
                                    cutoff>0, we_need_phase, calc_partials&&!no_ls_partials, calc_src);
            }
            else
            {//FIXME:  Is all of this really working?  Should I not have a VVH_that is half of VVH?
                //The below follows from Debye lineshape approximating half of VVH lineshape
                ArrayOfLineshapeSpec tmp;
                abs_lineshapeDefine( tmp, "Faddeeva_Algorithm_916", "VVH", 
                                     -1, verbosity );
                 
                xsec_single_line(   // OUTPUT   
                                    xsec_attenuation(joker,jj), calc_src?xsec_source(joker,jj):empty_vector, xsec_phase(joker,jj), 
                                    // HELPER
                                    attenuation, phase, 
                                    dFa_dx, dFb_dx, dFa_dy, dFb_dy, 
                                    this_f_range, fac, aux, 
                                    // FREQUENCY
                                    f_local, f_grid, f_grid.nelem(), cutoff, abs_lines[ii].F()+(precalc_zeeman?(Z_DF[ii]*H_magntitude_Zeeman):0), // Since vector is 0-length if no Zeeman pre-calculations
                                    // LINE STRENGTH
                                    abs_lines[ii].I0(), partition_ratio, K1*K2, abs_nlte_ratio, src_nlte_ratio,
                                    isotopologue_ratios.getParam(abs_lines[ii].Species(),
                                                                 abs_lines[ii].Isotopologue())[0].data[0],
                                    // ATMOSPHERIC TEMPERATURE
                                    t,
                                    // LINE SHAPE
                                    tmp[0].Ind_ls(), tmp[0].Ind_lsn(),
                                    // LINE BROADENING
                                    gamma_0, gamma_2, eta, df_0, df_2, sigma, f_VC,
                                    // LINE MIXING
                                    DV, Y, G, 
                                    // FEATURE FLAGS
                                    cutoff>0, we_need_phase, calc_partials&&!no_ls_partials, calc_src);
            }
            
            
            if(calc_partials)
            {
                // These needs to be calculated and returned when Temperature is in list
                Numeric dgamma_0_dT, ddf_0_dT,
                dgamma_2_dT, ddf_2_dT, deta_dT, df_VC_dT,
                dY_dT=0.0,dG_dT=0.0,dDV_dT=0.0,
                dQ_dT, dK2_dT, dabs_nlte_ratio_dT=0.0,
                atm_tv_low, atm_tv_upp;
                if(flag_partials.do_temperature())
                {
                    abs_lines[ii].LineMixing().GetLineMixingParams_dT(dY_dT, dG_dT, dDV_dT, t, flag_partials.Temperature_Perturbation(),
                                                                      p, lm_p_lim, 1);
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dT(dgamma_0_dT,dgamma_2_dT,
                                                                                      deta_dT, ddf_0_dT, 
                                                                                      ddf_2_dT, df_VC_dT, t, 
                                                                                      abs_lines[ii].Ti0(),p,
                                                                                      p_partial,this_species,h2o_index,
                                                                                      broad_spec_locations,
                                                                                      vmrs,verbosity);
                    GetLineScalingData_dT(dqt_dT_cache,
                                          dK2_dT,
                                          dQ_dT, 
                                          dabs_nlte_ratio_dT,
                                          atm_tv_low, 
                                          atm_tv_upp,
                                          qt_cache,
                                          abs_nlte_ratio,  
                                          partition_functions.getParamType(abs_lines[ii].Species(), abs_lines[ii].Isotopologue()),
                                          partition_functions.getParam(abs_lines[ii].Species(), abs_lines[ii].Isotopologue()),
                                          t,
                                          abs_lines[ii].Ti0(),
                                          flag_partials.Temperature_Perturbation(),
                                          abs_lines[ii].F(),
                                          calc_src,
                                          abs_lines[ii].Evlow(),
                                          abs_lines[ii].Evupp(),
                                          abs_lines[ii].NLTELowerIndex(),
                                          abs_lines[ii].NLTEUpperIndex(),
                                          t_nlte);
                }
                
                // These needs to be calculated when pressure broadening partial derivatives are needed
                // Note that this gives plenty of wasted calculations for all lines that are not specifically
                // requesting their individual partial derivatives...
                Numeric gamma_dSelf=0.0, gamma_dForeign=0.0, gamma_dWater=0.0,
                        psf_dSelf=0.0, psf_dForeign=0.0, psf_dWater=0.0, 
                        gamma_dSelfExponent=0.0, gamma_dForeignExponent=0.0, gamma_dWaterExponent=0.0,
                        psf_dSelfExponent=0.0, psf_dForeignExponent=0.0, psf_dWaterExponent=0.0;
                if(flag_partials.PressureBroadeningTerm(0)) // Self broadening gamma
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dSelfGamma(
                        gamma_dSelf,abs_lines[ii].Ti0()/t,p_partial);
                if(flag_partials.PressureBroadeningTerm(1)) // Foreign broadening gamma
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dForeignGamma(
                        gamma_dForeign,abs_lines[ii].Ti0()/t,p,p_partial,this_species,h2o_index,vmrs);
                if(flag_partials.PressureBroadeningTerm(2)) // Water broadening gamma
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dWaterGamma(
                        gamma_dWater,abs_lines[ii].Ti0()/t,p,this_species,h2o_index,vmrs,verbosity);
                if(flag_partials.PressureBroadeningTerm(3)) // Self broadening exponent
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dSelfExponent(
                        gamma_dSelfExponent, psf_dSelfExponent, abs_lines[ii].Ti0()/t,p_partial);
                if(flag_partials.PressureBroadeningTerm(4)) // Foreign broadening exponent
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dForeignExponent(
                        gamma_dForeignExponent,psf_dForeignExponent,abs_lines[ii].Ti0()/t,p,p_partial,this_species,h2o_index,vmrs);
                if(flag_partials.PressureBroadeningTerm(5)) // Water broadening exponent
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dWaterExponent(
                        gamma_dWaterExponent,psf_dWaterExponent,abs_lines[ii].Ti0()/t,p,this_species,h2o_index,vmrs,verbosity);
                if(flag_partials.PressureBroadeningTerm(6)) // Self broadening psf
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dSelfPsf(
                        psf_dSelf,abs_lines[ii].Ti0()/t,p_partial);
                if(flag_partials.PressureBroadeningTerm(7)) // Foreign broadening psf
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dForeignPsf(
                        psf_dForeign,abs_lines[ii].Ti0()/t,p,p_partial,this_species,h2o_index,vmrs);
                if(flag_partials.PressureBroadeningTerm(8)) // Water broadening psf
                    abs_lines[ii].PressureBroadening().GetPressureBroadeningParams_dWaterPsf(
                        psf_dWater,abs_lines[ii].Ti0()/t,p,this_species,h2o_index,vmrs,verbosity);
                
                Numeric dY0=0., dY1=0., dYexp=0., dG0=0., dG1=0., dGexp=0., dDV0=0., dDV1=0., dDVexp=0.;
                if(do_lm)
                {
                    if(flag_partials.ZerothTermLM())
                        abs_lines[ii].LineMixing().GetLineMixingParams_dZerothOrder(dY0, dG0, dDV0, t, p, lm_p_lim);
                    if(flag_partials.FirstTermLM())
                        abs_lines[ii].LineMixing().GetLineMixingParams_dFirstOrder(dY1, dG1, dDV1, t, p, lm_p_lim);
                    if(flag_partials.ExponentLM())
                        abs_lines[ii].LineMixing().GetLineMixingParams_dExponent(dYexp, dGexp, dDVexp, t, p, lm_p_lim);
                }
                    
                // Gather all new partial derivative calculations in this function
                partial_derivatives_lineshape_dependency(partial_xsec_attenuation,
                                                         partial_xsec_phase, 
                                                         partial_xsec_source,
                                                         flag_partials, 
                                                         attenuation,
                                                         phase,
                                                         fac,
                                                         dFa_dx, 
                                                         dFb_dx, 
                                                         dFa_dy, 
                                                         dFb_dy,
                                                         f_grid,
                                                         this_f_range,
                                                         //Temperature
                                                         t,
                                                         sigma,
                                                         K2,
                                                         dK2_dT,
                                                         abs_nlte_ratio,//K3
                                                         dabs_nlte_ratio_dT,
                                                         src_nlte_ratio,//K4
                                                         // Line parameters
                                                         abs_lines[ii].F(),
                                                         abs_lines[ii].I0(),
                                                         abs_lines[ii].Ti0(),
                                                         abs_lines[ii].Elow(),
                                                         abs_lines[ii].Evlow(),
                                                         abs_lines[ii].Evupp(),
                                                         abs_lines[ii].NLTELowerIndex() > -1?
                                                         t_nlte[abs_lines[ii].NLTELowerIndex()]:-1.0,
                                                         abs_lines[ii].NLTEUpperIndex() > -1?
                                                         t_nlte[abs_lines[ii].NLTEUpperIndex()]:-1.0,
                                                         Y,
                                                         dY_dT,
                                                         dY0,
                                                         dY1,
                                                         dYexp,
                                                         G,
                                                         dG_dT,
                                                         dG0,
                                                         dG1,
                                                         dGexp,
                                                         DV,
                                                         dDV_dT,
                                                         dDV0,
                                                         dDV1,
                                                         dDVexp,
                                                         abs_lines[ii].QuantumNumbers(),
                                                         abs_lines[ii].Species(),
                                                         abs_lines[ii].Isotopologue(),
                                                         // LINE SHAPE
                                                         ind_ls,
                                                         ind_lsn,
                                                         df_0,//FIXME: For H-T, this part is difficult
                                                         ddf_0_dT,
                                                         gamma_0,
                                                         dgamma_0_dT,
                                                         gamma_dSelf,
                                                         gamma_dForeign,
                                                         gamma_dWater,
                                                         psf_dSelf,
                                                         psf_dForeign,
                                                         psf_dWater,
                                                         gamma_dSelfExponent,
                                                         gamma_dForeignExponent,
                                                         gamma_dWaterExponent,
                                                         psf_dSelfExponent,
                                                         psf_dForeignExponent,
                                                         psf_dWaterExponent,
                                                         // Partition data parameters
                                                         dQ_dT,
                                                         // Magnetic variables
                                                         precalc_zeeman?Z_DF[ii]:0,
                                                         H_magntitude_Zeeman,
                                                         precalc_zeeman,
                                                         // Programming
                                                         jj, 
                                                         calc_partials_phase,
                                                         calc_src);
            }
        }
    }
}


/*! cross-section replacement computer 
 *  
 * This will work as the interface for all line-by-line computations 
 * lacking special demands
 *  
 *  \retval xsec                Cross section of one tag group. This is now the
 *                              true attenuation cross section in units of m^2.
 *  \retval source              Cross section of one tag group. This is now the
 *                              true source cross section in units of m^2.
 *  \retval phase               Cross section of one tag group. This is now the
 *                              true phase cross section in units of m^2.
 *  \retval dxsec               Partial derivatives of xsec.
 *  \retval dsource             Partial derivatives of source.
 *  \retval dphase              Partial derivatives of phase.
 * 
 *  \param flag_partials        Partial derivatives flags.
 *  \param f_grid               Frequency grid.
 *  \param abs_p                Pressure grid.
 *  \param abs_t                Temperatures associated with abs_p.
 *  \param abs_t_nlte           Non-lte temperatures for various energy levels.
 *  \param all_vmrs             Gas volume mixing ratios [nspecies, np].
 *  \param abs_species          Species tags for all species.
 *  \param this_species         Index of the current species in abs_species.
 *  \param abs_lines            The spectroscopic line list.
 *  \param Z_DF                 The Zeeman line center shift over the magnitude of the magnetic field.
 *  \param H_magntitude_Zeeman  The magnitude of the magnetic field required by Zeeman effect.
 *  \param lm_p_lim             Line mixing pressure limit
 *  \param isotopologue_ratios  Isotopologue ratios.
 *  \param partition_functions  Partition functions.
 *  \param verbosity            Verbosity level.
 * 
 *  \author Richard Larsson
 *  \date   2013-04-24
 * 
 */
void xsec_species2(MatrixView xsec,
                   MatrixView source,
                   MatrixView phase,
                   ArrayOfMatrix& dxsec_dx,
                   ArrayOfMatrix& dsource_dx,
                   ArrayOfMatrix& dphase_dx,
                   const PropmatPartialsData& flag_partials,
                   ConstVectorView f_grid,
                   ConstVectorView abs_p,
                   ConstVectorView abs_t,
                   ConstMatrixView abs_t_nlte,
                   ConstMatrixView all_vmrs,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const Index this_species,
                   const ArrayOfLineRecord& abs_lines,
                   const Vector& Z_DF,
                   const Numeric H_magntitude_Zeeman,
                   const Numeric lm_p_lim,
                   const SpeciesAuxData& isotopologue_ratios,
                   const SpeciesAuxData& partition_functions,
                   const Index& binary_speedup,
                   const Verbosity& verbosity)
{
  // Size of problem
  const Index np = abs_p.nelem();                      // number of pressure levels
  const Index nf = f_grid.nelem();                     // number of Dirac frequencies
  DEBUG_ONLY(const Index ns = abs_species.nelem();)    // number of species vmrs
  const Index nl = abs_lines.nelem();                  // number of lines in the catalog
  const Index nz = Z_DF.nelem();                       // number of lines affected by Zeeman effect
  const Index nj = flag_partials.nelem();              // number of partial derivatives
  const Index nt = abs_t_nlte.nrows();                 // number of energy levels in NLTE
  
  // Type of problem
  const bool do_phase = nz;                  // phase return is requested only by Zeeman effect calculations
  const bool do_nonlte = nt;                 // source return is requested only if there are energy levels in NLTE
  DEBUG_ONLY(const bool do_jacobians = nj;)  // Jacobian return is requested only if there are partial derivatives
  
  // Test if the size of the problem is 0
  if(not np or not nf or not nl)
    return;
  
  // Standard variables that must always have correct size upon calling this function
  assert(this_species < ns);
  assert(abs_t.nelem() == np);
  assert(all_vmrs.ncols() == np and all_vmrs.nrows() == ns);
  assert(xsec.nrows() == nf and xsec.ncols() == np);
  assert(dxsec_dx.nelem() == nj);
  
  // Non-LTE stuff must follow these rules if applied
  assert((source.nrows() == nf and source.ncols() == np) or not do_nonlte);
  assert(dsource_dx.nelem() == nj or not do_nonlte);
  assert(np == abs_t_nlte.ncols() or not do_nonlte);
  
  // Phase-related only works for single layer
  assert((np == 1 and nz == nl) or not do_phase);
  assert((phase.nrows() == nf and phase.ncols() == 1) or not do_phase);
  assert(dphase_dx.nelem() == nj or not do_phase);
  assert(H_magntitude_Zeeman >= 0); // When not using it it should be 0
  
  // Jacobian must be the correct size throughout
  DEBUG_ONLY(for(auto& m : dxsec_dx) {assert((m.nrows() == nf and m.ncols() == np) or not do_jacobians);})
  DEBUG_ONLY(for(auto& m : dphase_dx) {assert((m.nrows() == nf and m.ncols() == np) or not do_jacobians or not do_phase);})
  DEBUG_ONLY(for(auto& m : dsource_dx) {assert((m.nrows() == nf and m.ncols() == np) or not do_jacobians or not do_nonlte);})
  
  // Algorithmic stuff to test that variables are as expected
  //assert(not std::any_of(abs_t.begin(), abs_t.end(), [](Numeric n){return n <= 0.;}));
  //assert(not std::any_of(abs_p.begin(), abs_p.end(), [](Numeric n){return n <= 0.;}));
  
  // Water index in VMRS is constant for all levels and lines
  const Index h2o_index = find_first_species_tg(abs_species, species_index_from_species_name("H2O"));
  
  // Broadening species locations in VMRS are constant for all levels and lines
  ArrayOfIndex broad_spec_locations;
  find_broad_spec_locations(broad_spec_locations, abs_species, this_species);
  
  // Results vectors are initialized first and then copied to the threads later
  ComplexVector F(nf), N(do_nonlte?nf:0);
  ArrayOfComplexVector dF(nj), dN(do_nonlte?nj:0);
  for(auto& aocv : dF) aocv.resize(nf);
  for(auto& aocv : dN) aocv.resize(nf);
  Range this_xsec_range(joker);

  for(Index ip = 0; ip < np; ip++)
  {
    // Constants for this level
    const Numeric& temperature = abs_t[ip];
    const Numeric& pressure = abs_p[ip];
    const Numeric partial_pressure = pressure * all_vmrs(this_species, ip);
    
    // Quasi-constants for this level, defined here to speed up later computations
    Index this_iso = -1; // line isotopologue number
    Numeric t0 = -1; // line temperature
    Numeric dc, ddc_dT, qt, qt0, dqt_dT; // Doppler and partition functions
    
    for(Index il = 0; il < nl; il++)
    {
      const LineRecord& line = abs_lines[il];
      
      // Partition function depends on isotopologue and line temperatures.
      // Both are commonly constant in a single catalog.  They are, however,
      // allowed to change so we must check that they do not
      if(line.Isotopologue() not_eq this_iso or line.Ti0() not_eq t0)
      {
        // set new line temperature
        t0 = line.Ti0();
        
        // set new partition function
        partition_function(qt0, qt,
          t0, temperature,
          partition_functions.getParamType(line.Species(), line.Isotopologue()),
          partition_functions.getParam(line.Species(), line.Isotopologue()));
        
        // if needed, set the partition function derivative
        if(flag_partials.do_temperature())
          dpartition_function_dT(dqt_dT, qt,
            temperature, flag_partials.Temperature_Perturbation(),
            partition_functions.getParamType(line.Species(), line.Isotopologue()),
            partition_functions.getParam(line.Species(), line.Isotopologue()));
      }
      
      // Same rule as for partition function applies to the Doppler broadening.
      // It is however more simple since it only depends on the mass of the species
      // of the line, so only check for this
      if(line.Isotopologue() not_eq this_iso)
      {
        // set new isotopologue number
        this_iso = line.Isotopologue();
        
        // set the frequency-independent part of the line shape
        dc = Linefunctions::DopplerConstant(temperature, line.IsotopologueData().Mass());
        
        // if needed, set the partial derivative of the frequency-independent part of the line shape
        if(flag_partials.do_temperature())
          ddc_dT = Linefunctions::dDopplerConstant_dT(temperature, line.IsotopologueData().Mass());
      }
      
      // we now compute the line shape
      Linefunctions::set_cross_section_for_single_line(F, dF, N, dN, this_xsec_range,
        flag_partials, line, f_grid, all_vmrs(joker, ip), nt?abs_t_nlte(joker, ip):Vector(0),
        pressure, temperature, dc, partial_pressure, 
        isotopologue_ratios.getParam(line.Species(), this_iso)[0].data[0],
        H_magntitude_Zeeman, ddc_dT, lm_p_lim, do_phase?Z_DF[il]:0.0, qt, dqt_dT, qt0,
        broad_spec_locations, this_species, h2o_index, verbosity, 
        false, binary_speedup);
      
      // range-based arguments that need be made to work for both complex and numeric
      const Index extent = (this_xsec_range.get_extent()<0)     ?
                (nf-this_xsec_range.get_start())                :
                this_xsec_range.get_extent();
      const Range this_out_range(this_xsec_range.get_start(), extent);
        
      for(Index i = 0; i < extent; i++)
      {
        xsec(this_out_range, ip)[i] += F[this_xsec_range][i].real();
        
        if(do_phase)
          phase(this_out_range, ip)[i] += F[this_xsec_range][i].imag();
        
        if(do_nonlte)
          source(this_out_range, ip)[i] += N[this_xsec_range][i].real();
        
        for(Index j = 0; j < nj; j++)
        {
          dxsec_dx[j](this_out_range, ip)[i] += dF[j][this_xsec_range][i].real();
          
          if(do_phase)
            dphase_dx[j](this_out_range, ip)[i] += dF[j][this_xsec_range][i].imag();
          
          if(do_nonlte)
            dsource_dx[j](this_out_range, ip)[i] += dN[j][this_xsec_range][i].real();
        }
      }
    }
  }
}
