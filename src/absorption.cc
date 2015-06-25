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
#include "make_vector.h"

#include "global_data.h"


/** Mapping of species auxiliary type names to SpeciesAuxData::AuxType enum */
static const char *SpeciesAuxTypeNames[] = {
    "NONE",
    "ISORATIO",   //Built-in type
    "ISOQUANTUM",
    "PART_TFIELD",
    "PART_COEFF", //Built-in type
    "PART_COEFF_VIBROT"
};


/** The map associated with species_data. */
std::map<String, Index> SpeciesMap;


// member fct of isotopologuerecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopologueRecord::CalculatePartitionFctAtTempFromCoeff( Numeric
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
Numeric IsotopologueRecord::CalculatePartitionFctAtTempFromData( Numeric
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
            mparam_type[isp][iso] = SpeciesAuxData::PT_NONE;
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

    while (this_auxtype < PT_FINAL_ENTRY && auxtype != SpeciesAuxTypeNames[this_auxtype])
        this_auxtype++;

    if (this_auxtype != PT_FINAL_ENTRY)
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
    assert(mparam_type[species][isotopologue] < PT_FINAL_ENTRY);
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
        for (Index iiso = 0; iiso < species_data[isp].Isotopologue().nelem(); iiso++)
        {
            ratios[0].data[0] = species_data[isp].Isotopologue()[iiso].Abundance();
            sad.setParam(isp, iiso,
                         SpeciesAuxData::PT_ISOTOPOLOGUE_RATIO, ratios);
        }
}


/*! Define the species data map.

    \author Stefan Buehler */
void define_species_map()
{
  using global_data::species_data;

  for ( Index i=0 ; i<species_data.nelem() ; ++i)
    {
      SpeciesMap[species_data[i].Name()] = i;
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
                   const Verbosity&         verbosity )
{
    // Make lineshape and species lookup data visible:
    using global_data::lineshape_data;
    
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
    Vector ls_phase(nf+1);
    Vector fac(nf+1);
    
    const bool cut = (cutoff != -1) ? true : false;
    
    const bool calc_phase = lineshape_data[ind_ls].Phase();
    
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
        os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
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
        os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()] or [0,0].\n"
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
    
    if ( xsec_phase.nrows() != nf || xsec_phase.ncols() != np )
    {
        std::ostringstream os;
        os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
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
firstprivate(ls_attenuation, ls_phase, fac, f_local, aux)
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
            VectorView xsec_i_phase = xsec_phase(Range(joker),i);
            
            
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
	    Matrix xsec_accum_source(n_lbl_threads, xsec_i_attenuation.nelem(), 0);
            Matrix xsec_accum_phase(n_lbl_threads, xsec_i_phase.nelem(), 0);
            
            ConstVectorView vmrs = all_vmrs(joker,i);
            
            // Loop all lines:
            if (nl)
#pragma omp parallel for                   \
if (!arts_omp_in_parallel()               \
&& nl >= arts_omp_get_max_threads())  \
firstprivate(ls_attenuation, ls_phase, fac, f_local, aux)
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
                        Numeric gamma=0, deltaf=0, partition_ratio,boltzmann_ratio, abs_nlte_ratio, src_nlte_ratio;
                        l_l.PressureBroadening().GetPressureBroadeningParams(gamma,deltaf,
                                                                             l_l.Ti0()/t_i,p_i,
                                                                             p_partial,this_species,h2o_index,
                                                                             broad_spec_locations,
                                                                             vmrs,verbosity);
                        
                        l_l.GetLineScalingData(partition_ratio, boltzmann_ratio, abs_nlte_ratio, src_nlte_ratio, 
                                               t_i, calc_src, t_nlte_i);
                        
                        xsec_single_line(xsec_accum_attenuation(arts_omp_get_thread_num(),joker),
					 xsec_accum_source(arts_omp_get_thread_num(),joker),
                                         xsec_accum_phase(arts_omp_get_thread_num(),joker),
                                         ls_attenuation,
                                         ls_phase,
                                         fac,
                                         f_local,
                                         aux,
                                         isotopologue_ratios,
                                         f_grid,
                                         l_l.F(),
                                         l_l.I0(),
                                         partition_ratio,
                                         boltzmann_ratio,
                                         abs_nlte_ratio,
                                         src_nlte_ratio,
                                         l_l.IsotopologueData().Mass(),
                                         t_i,
                                         gamma,
                                         deltaf,
                                         cutoff,
                                         0,
                                         0,
                                         0,
                                         nf,
                                         ind_ls,
                                         ind_lsn,
                                         l_l.Species(),
                                         l_l.Isotopologue(),
                                         cut,
                                         calc_phase,
                                         calc_src);
                        
                    } // end of try block
                    catch (runtime_error e)
                    {
                        #pragma omp critical (xsec_species_fail)
                        { fail_msg = e.what(); failed = true; }
                    }
                    
                } // end of parallel LBL loop
                
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
                
                if (calc_phase)
                    for (Index j=0; j<xsec_accum_phase.nrows(); ++j)
                    {
                        xsec_i_phase += xsec_accum_phase(j, Range(joker));
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
   \retval xsec_accum_phase         Cross section of one tag group. This is now the
                                    true dispersion cross section in units of m^2.
                                    It has inputs of all previously calculated lines.
   \retval attenuation              Input only to increase speed.  Holds attenuation internally.
   \retval phase                    Input only to increase speed.  Holds phase internally.
   \retval fac                      Input only to increase speed.  Holds lineshape factor internally.
   \retval f_local                  Input only to increase speed.  Holds f_grid internally.
   \retval aux                      Input only to increase speed.  Holds f_grid factor internally.
   \param isotopologue_ratios       Isotopologue ratios.const ArrayOfIndex& broad_spec_locations,
   \param f_grid,                   Frequency grid
   \param vmrs,                     Gas volume mixing ratios [nspecies].
   \param Gamma_foreign,            Foreign line pressure broadening
   \param N_foreign,                Foregin line temperature exponent
   \param Delta_foreign,            Foreign line pressure shift
   \param F0,                       Line center
   \param intensity,                Line intensity
   \param part_fct_ratio,           Partition function ratio
   \param Isotopologue_Mass,        Mass of isotopologue
   \param e_lower,                  Line lower energy state
   \param T0,                       Line reference temperature
   \param Sgam,                     Line self broadening
   \param Nself,                    Line self temperature exponent
   \param Tgam,                     Line pressure broadening reference temperature
   \param Agam,                     Line air broadening
   \param Nair,                     Line air temperature exponent
   \param Psf,                      Line air pressure shift
   \param temperature,              Atmospheric temperature
   \param pressure,                 Atmospheric pressure
   \param p_partial,                Molecular pressure of line
   \param cutoff,                   Lineshape cutoff.
   \param LM_DF,                    Line mixing frequency shift
   \param LM_Y,                     Line mixing dispersion dependency
   \param LM_G,                     Line mixing added attenuation
   \param this_species,             Index of the current species in abs_species.
   \param nf,                       Number of frequencies to calculate
   \param ind_ls,                   Index to lineshape function.
   \param ind_lsn,                  Index to lineshape norm.
   \param LineRecord_Version,       Line data version (for pressure broadening scheme)
   \param LineRecord_Species,       Line data species number
   \param LineRecord_Isotopologue,  Line data isotopologue identifier
   \param quadratic_lineshape,      Is line shape quadratic?
   \param cut,                      Is cutoff applied?
   \param calc_phase,               Is dispersion calculated?
 
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
void xsec_single_line(VectorView xsec_accum_attenuation, 
		      VectorView xsec_accum_source, 
                      VectorView xsec_accum_phase, 
                      Vector& attenuation, 
                      Vector& phase,
                      Vector& fac, 
                      Vector& f_local, 
                      Vector& aux, 
                      const SpeciesAuxData& isotopologue_ratios,
                      const Vector& f_grid, 
                      Numeric F0, 
                      Numeric intensity, 
                      const Numeric part_fct_ratio,  
                      const Numeric boltzmann_ratio,
                      const Numeric abs_nlte_ratio,
                      const Numeric src_nlte_ratio,
                      const Numeric Isotopologue_Mass,
                      const Numeric temperature, 
                      const Numeric gamma,
                      const Numeric deltaf,
                      const Numeric cutoff,
                      const Numeric LM_DF,
                      const Numeric LM_Y, 
                      const Numeric LM_G, 
                      const Index nf, 
                      const Index ind_ls, 
                      const Index ind_lsn,  
                      const Index LineRecord_Species, 
                      const Index LineRecord_Isotopologue, 
                      const bool calc_cut, 
                      const bool calc_phase,
                      const bool calc_src)
{//asdasd;
    
    extern const Numeric BOLTZMAN_CONST;
    extern const Numeric AVOGADROS_NUMB;
    extern const Numeric SPEED_OF_LIGHT;
    static const Numeric doppler_const = sqrt(2.0 * BOLTZMAN_CONST *
                                              AVOGADROS_NUMB) / SPEED_OF_LIGHT;
    
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
    F0 += deltaf + LM_DF;
    
    // 3. Doppler broadening without the sqrt(ln(2)) factor, which
    // seems to be redundant.
    const Numeric sigma = F0 * doppler_const *
    sqrt( temperature / Isotopologue_Mass);
    
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
        
        // Calculate the line shape:
        global_data::lineshape_data[ind_ls].Function()(attenuation, phase,
                                                       aux, F0, gamma, sigma,
                                                       f_local[Range(i_f_min,nfls)]);
        
        // Calculate the chosen normalization factor:
        global_data::lineshape_norm_data[ind_lsn].Function()(fac, F0,
                                                             f_local[Range(i_f_min,nfls)],
                                                             temperature);
        
        // Reset f_local for next loop through now that cutoff is calculated
        if( i_f_max < nf )
            f_local[i_f_max] = f_grid[i_f_max];
        
        // Get a handle on the range of xsec that we want to change.
        // We use nfl here, which could be one less than nfls.
        VectorView this_xsec_attenuation      = xsec_accum_attenuation[Range(i_f_min,nfl)];
        Vector empty_vector;
	    VectorView this_xsec_source           = calc_src ? xsec_accum_source[Range(i_f_min,nfl)] : empty_vector;
        VectorView this_xsec_phase            = xsec_accum_phase[Range(i_f_min,nfl)];
        
        // Get handles on the range of ls and fac that we need.
        VectorView this_ls_attenuation  = attenuation[Range(0,nfl)];
        VectorView this_ls_phase  = phase[Range(0,nfl)];
        VectorView this_fac = fac[Range(0,nfl)];
        
        // Line Mixing
        if(LM_G!=0)
        {
            // Apply line mixing to both cutoff and other values
            this_ls_attenuation*=1+LM_G;
            if(calc_cut)
                attenuation[nfls-1]*=1+LM_G;
        }
        if(LM_Y!=0)
        {
            // Apply line mixing to both cutoff and other values
            Vector tmp = this_ls_phase;
            tmp *= LM_Y;
            this_ls_attenuation+=tmp;
            if(calc_cut)
                attenuation[nfls-1]+=phase[nfls-1]*LM_Y;
        }
        
        // cutoff ?
        if ( calc_cut )
        {
            // Subtract baseline for cutoff frequency
            // The index nfls-1 should be exactly the index pointing
            // to the value at the cutoff frequency.
            // Subtract baseline from xsec.
            this_ls_attenuation -= attenuation[nfls-1]; // cutoff is constant
        }
        
        // Add line to xsec.
        {
            // To make the loop a bit faster, precompute all constant
            // factors. These are:
            // 1. Total number density of the air. --> Not
            //    anymore, we now to real cross-sections
            // 2. Line intensity.
            // 3. Isotopologue ratio.
            //
            // The isotopologue ratio must be applied here, since we are
            // summing up lines belonging to different isotopologues.
            
            //                const Numeric factors = n * intensity * l_l.IsotopologueData().Abundance();
            //                    const Numeric factors = intensity * l_l.IsotopologueData().Abundance();
            const Numeric factors = isotopologue_ratios.getParam(
                                           LineRecord_Species,
                                           LineRecord_Isotopologue)[0].data[0];
            
            // We have to do:
            // xsec(j,i) += factors * ls[j] * fac[j];
            //
            // We use ls as a dummy to compute the product, then add it
            // to this_xsec.
            
            this_ls_attenuation *= this_fac;
	    if(calc_src)
	    {
	      this_ls_attenuation *= factors * intensity * abs_nlte_ratio;
	      this_xsec_attenuation += this_ls_attenuation;
	      this_ls_attenuation *= ( src_nlte_ratio/abs_nlte_ratio - 1.0 );//Only retain the additional source term for NLTE calculations
	      this_xsec_source += this_ls_attenuation; // note that the multiplication above solves 
	      
	      if (calc_phase)
              {
                  this_ls_phase *= this_fac;
                  this_ls_phase *= factors * intensity * abs_nlte_ratio; // Is this right?  
                  this_xsec_phase += this_ls_phase;
              }
	    }
	    else
	    {
	      this_ls_attenuation *= factors * intensity;
	      this_xsec_attenuation += this_ls_attenuation;
              
              if (calc_phase)
              {
                  this_ls_phase *= this_fac;
                  this_ls_phase *= factors*intensity;
                  this_xsec_phase += this_ls_phase;
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
                                        ConstVectorView          f_grid,
                                        ConstVectorView          abs_p,
                                        ConstVectorView          abs_t,
                                        ConstMatrixView          abs_t_nlte,
                                        ConstMatrixView          all_vmrs,
                                        const ArrayOfArrayOfSpeciesTag& abs_species,
                                        const Index              this_species,
                                        const ArrayOfLineRecord& abs_lines,
                                        const Vector&            Z_DF,
                                        const Index              ind_ls,
                                        const Index              ind_lsn,
                                        const Numeric            lm_p_lim,
                                        const Numeric            cutoff,
                                        const SpeciesAuxData&    isotopologue_ratios,
                                        const Verbosity&         verbosity )
{
    // Must have the phase
    using global_data::lineshape_data;
    if (! lineshape_data[ind_ls].Phase())
    {
        std::ostringstream os;
        os <<  "This is an error message. You are using " << lineshape_data[ind_ls].Name() <<
        ".\n"<<"This line shape does not include phase in its calculations and\nis therefore invalid for " <<
        "line mixing or Zeeman.\n\n";
        throw std::runtime_error(os.str());
    }
    
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
    
    const bool cut = (cutoff != -1) ? true : false;
    
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
    ArrayOfIndex tmp_none_mixed_lines_index;
    Vector attenuation(f_grid.nelem()+1),
    phase(f_grid.nelem()+1),
    aux(f_grid.nelem()+1),
    fac(f_grid.nelem()+1),
    f_local(f_grid.nelem()+1);
    f_local[Range(0,f_grid.nelem())]=f_grid;
    
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
        
        const Numeric p_partial = p * vmrs[this_species];
        
        for(Index ii=0; ii<abs_lines.nelem();ii++)
        {
            // Pressure broadening parameters
            Numeric gamma=0, deltaf_pressure=0, partition_ratio,boltzmann_ratio, abs_nlte_ratio, src_nlte_ratio; // Set to zero since case with 0 exist
            abs_lines[ii].PressureBroadening().GetPressureBroadeningParams(gamma,deltaf_pressure,
                                                                            abs_lines[ii].Ti0()/t,p,
                                                                            p_partial,this_species,h2o_index,
                                                                            broad_spec_locations,
                                                                            vmrs,verbosity);
            
            // Line mixing parameters
            Numeric Y=0,DV=0,G=0; // Set to zero since case with 0 exist
            abs_lines[ii].LineMixing().GetLineMixingParams( Y,  G,  DV,  t, p, lm_p_lim, 1);
            
            abs_lines[ii].GetLineScalingData(partition_ratio,boltzmann_ratio, abs_nlte_ratio, src_nlte_ratio, 
                                             t, calc_src, t_nlte);
            
            // Still an ugly case with non-resonant line near 0 frequency, since this is actually a band...
            if( LineMixingData::LM_LBLRTM_O2NonResonant != abs_lines[ii].LineMixing().Type() )
            {
                xsec_single_line(   xsec_attenuation(joker,jj), 
				    calc_src?xsec_source(joker,jj):empty_vector,
                                    xsec_phase(joker,jj), 
                                    attenuation, 
                                    phase,
                                    fac, 
                                    f_local, 
                                    aux, 
                                    isotopologue_ratios,
                                    f_grid, 
                                    abs_lines[ii].F()+(precalc_zeeman?Z_DF[ii]:0), // Since vector is 0-length if no Zeeman pre-calculations
                                    abs_lines[ii].I0(), 
                                    partition_ratio, 
                                    boltzmann_ratio,
                                    abs_nlte_ratio,
                                    src_nlte_ratio,
                                    abs_lines[ii].IsotopologueData().Mass(),
                                    t,
                                    gamma,
                                    deltaf_pressure,
                                    cutoff,
                                    DV,
                                    Y, 
                                    G, 
                                    f_grid.nelem(), 
                                    ind_ls, 
                                    ind_lsn, 
                                    abs_lines[ii].Species(), 
                                    abs_lines[ii].Isotopologue(), 
                                    cutoff!=-1,
                                    1,
				    calc_src);
            }
            else
            {
                //The below follows from Debye lineshape approximating half of VVH lineshape
                ArrayOfLineshapeSpec tmp;
                abs_lineshapeDefine( tmp, "Faddeeva_Algorithm_916", "VVH", 
                                     -1, verbosity );
                
                xsec_single_line(   xsec_attenuation(joker,jj), 
				    calc_src?xsec_source(joker,jj):Vector(0),
                                    xsec_phase(joker,jj), 
                                    attenuation, 
                                    phase,
                                    fac, 
                                    f_local, 
                                    aux, 
                                    isotopologue_ratios,
                                    f_grid, 
                                    abs_lines[ii].F()+(precalc_zeeman?Z_DF[ii]:0), // Since vector is 0-length if no Zeeman pre-calculations 
                                    abs_lines[ii].I0(), 
                                    partition_ratio,
                                    boltzmann_ratio,
                                    abs_nlte_ratio,
                                    src_nlte_ratio,
                                    abs_lines[ii].IsotopologueData().Mass(),
                                    t,
                                    gamma,
                                    deltaf_pressure,
                                    cutoff,
                                    DV,
                                    Y, 
                                    G, // Non-resonant lineshape is half VVH --- this scales
                                    f_grid.nelem(), 
                                    tmp[0].Ind_ls(), 
                                    tmp[0].Ind_lsn(), 
                                    abs_lines[ii].Species(), 
                                    abs_lines[ii].Isotopologue(), 
                                    cutoff!=-1,
                                    1,
				    calc_src);
            }
            
        }
    }
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
 
 \retval xsec   Cross section of one tag group. This is now the
 true absorption cross section in units of m^2.
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
 
 */
void xsec_species_old_unused( MatrixView               xsec_attenuation,
                              MatrixView               xsec_phase,
                              ConstVectorView          f_grid,
                              ConstVectorView          abs_p,
                              ConstVectorView          abs_t,
                              ConstMatrixView          all_vmrs,
                              const ArrayOfArrayOfSpeciesTag& abs_species,
                              const Index              this_species,
                              const ArrayOfLineRecord& abs_lines,
                              const Index              ind_ls,
                              const Index              ind_lsn,
                              const Numeric            cutoff,
                              const SpeciesAuxData&    isotopologue_ratios,
                              const Verbosity&         verbosity )
{
    // Make lineshape and species lookup data visible:
    using global_data::lineshape_data;
    using global_data::lineshape_norm_data;
    
    // speed of light constant
    extern const Numeric SPEED_OF_LIGHT;
    
    // Boltzmann constant
    extern const Numeric BOLTZMAN_CONST;
    
    // Avogadros constant
    extern const Numeric AVOGADROS_NUMB;
    
    // Planck constant
    extern const Numeric PLANCK_CONST;
    
    // sqrt(ln(2))
    // extern const Numeric SQRT_NAT_LOG_2;
    
    // Constant within the Doppler Broadening calculation:
    const Numeric doppler_const = sqrt( 2.0 * BOLTZMAN_CONST *
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
    Vector ls_phase(nf+1);
    Vector fac(nf+1);
    
    const bool cut = (cutoff != -1) ? true : false;
    
    const bool calc_phase = lineshape_data[ind_ls].Phase();
    
    // Check that the frequency grid is sorted in the case of lineshape
    // with cutoff. Duplicate frequency values are allowed.
    if (cut)
    {
        if ( ! is_sorted( f_grid ) )
        {
            ostringstream os;
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
        ostringstream os;
        os << "abs_t contains at least one negative temperature value.\n"
        << "This is not allowed.";
        throw std::runtime_error(os.str());
    }
    
    // We need a local copy of f_grid which is 1 element longer, because
    // we append a possible cutoff to it.
    // The initialization of this has to be inside the line loop!
    Vector f_local( nf + 1 );
    
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
        ostringstream os;
        os << "Variable abs_t must have the same dimension as abs_p.\n"
        << "abs_t.nelem() = " << abs_t.nelem() << '\n'
        << "abs_p.nelem() = " << np;
        throw std::runtime_error(os.str());
    }
    
    // all_vmrs should have dimensions [nspecies, np]:
    
    if ( all_vmrs.ncols() != np )
    {
        ostringstream os;
        os << "Number of columns of all_vmrs must match abs_p.\n"
        << "all_vmrs.ncols() = " << all_vmrs.ncols() << '\n'
        << "abs_p.nelem() = " << np;
        throw std::runtime_error(os.str());
    }
    
    const Index nspecies = abs_species.nelem();
    
    if ( all_vmrs.nrows() != nspecies)
    {
        ostringstream os;
        os << "Number of rows of all_vmrs must match abs_species.\n"
        << "all_vmrs.nrows() = " << all_vmrs.nrows() << '\n'
        << "abs_species.nelem() = " << nspecies;
        throw std::runtime_error(os.str());
    }
    
    // With abs_h2o it is different. We do not really need this in most
    // cases, only the Rosenkranz lineshape for oxygen uses it. There is
    // a global (scalar) default value of -1, that we are expanding to a
    // vector here if we find it. The Rosenkranz lineshape does a check
    // to make sure that the value is actually set, and not the default
    // value. 
    //  Vector abs_h2o(np);
    //  if ( abs_h2o_orig.nelem() == np )
    //    {
    //      abs_h2o = abs_h2o_orig;
    //    }
    //  else
    //    {
    //      if ( ( 1   == abs_h2o_orig.nelem()) && 
    //           ( -.99 > abs_h2o_orig[0]) )
    //        {
    //          // We have found the global default value. Expand this to a
    //          // vector with the right length, by copying -1 to all
    //          // elements of abs_h2o.
    //          abs_h2o = -1;         
    //        }
    //      else
    //        {
    //          ostringstream os;
    //          os << "Variable abs_h2o must have default value -1 or the\n"
    //             << "same dimension as abs_p.\n"
    //             << "abs_h2o.nelem() = " << abs_h2o.nelem() << '\n'
    //             << "abs_p.nelem() = " << np;
    //          throw runtime_error(os.str());
    //        }
    //    }
    
    // Check that the dimension of xsec is indeed [f_grid.nelem(),
    // abs_p.nelem()]:
    if ( xsec_attenuation.nrows() != nf || xsec_attenuation.ncols() != np )
    {
        ostringstream os;
        os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
        << "[xsec_attenuation.nrows(),xsec_attenuation.ncols()] = [" << xsec_attenuation.nrows()
        << ", " << xsec_attenuation.ncols() << "]\n"
        << "f_grid.nelem() = " << nf << '\n'
        << "abs_p.nelem() = " << np;
        throw std::runtime_error(os.str());
    }
    if ( xsec_phase.nrows() != nf || xsec_phase.ncols() != np )
    {
        ostringstream os;
        os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
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
    
    String fail_msg;
    bool failed = false;
    
    // Loop all pressures:
    if (np)
        #pragma omp parallel for                    \
        if (!arts_omp_in_parallel()               \
            && np >= arts_omp_get_max_threads())  \
            firstprivate(ls_attenuation, ls_phase, fac, f_local, aux)
            for ( Index i=0; i<np; ++i )
            {
                if (failed) continue;
                
                // Store input profile variables, this is perhaps slightly faster.
                const Numeric p_i       = abs_p[i];
                const Numeric t_i       = abs_t[i];
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
                VectorView xsec_i_phase = xsec_phase(Range(joker),i);
                
                
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
                Matrix xsec_accum_phase(n_lbl_threads, xsec_i_phase.nelem(), 0);
                
                
                // Loop all lines:
                if (nl)
                    #pragma omp parallel for                   \
                    if (!arts_omp_in_parallel()               \
                        && nl >= arts_omp_get_max_threads())  \
                        firstprivate(ls_attenuation, ls_phase, fac, f_local, aux)
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
                                // Copy f_grid to the beginning of f_local. There is one
                                // element left at the end of f_local.  
                                // THIS HAS TO BE INSIDE THE LINE LOOP, BECAUSE THE CUTOFF
                                // FREQUENCY IS ALWAYS PUT IN A DIFFERENT PLACE!
                                f_local[Range(0,nf)] = f_grid;
                                
                                // This will hold the actual number of frequencies to add to
                                // xsec later on:
                                Index nfl = nf;
                                
                                // This will hold the actual number of frequencies for the
                                // call to the lineshape functions later on:
                                Index nfls = nf;      
                                
                                // abs_lines[l] is used several times, this construct should be
                                // faster (Oliver Lemke)
                                const LineRecord& l_l = abs_lines[l];  // which line are we dealing with
                                
                                // Make sure that catalogue version is either 3 or 4 (no other
                                // versions are supported yet):
                                if ( 3!=l_l.Version() && 4!=l_l.Version() )
                                {
                                    ostringstream os;
                                    os << "Unknown spectral line catalogue version (artscat-"
                                    << l_l.Version() << ").\n"
                                    << "Allowed are artscat-3 and artscat-4.";
                                    throw std::runtime_error(os.str());
                                }
                                
                                // Center frequency in vacuum:
                                Numeric F0 = l_l.F();
                                
                                // Intensity is already in the right units (Hz*m^2). It also
                                // includes already the isotopologue ratio. Needs only to be
                                // coverted to the actual temperature and multiplied by total
                                // number density and lineshape.
                                Numeric intensity = l_l.I0();
                                
                                // Lower state energy is already in the right unit (Joule).
                                Numeric e_lower = l_l.Elow();
                                
                                // Upper state energy:
                                Numeric e_upper = e_lower + F0 * PLANCK_CONST;
                                
                                // Get the ratio of the partition function.
                                // This will throw a runtime error if no data exists.
                                // Important: This function needs both the reference
                                // temperature and the actual temperature, because the
                                // reference temperature can be different for each line,
                                // even of the same species.
                                Numeric part_fct_ratio =
                                l_l.IsotopologueData().CalculatePartitionFctRatio( l_l.Ti0(),
                                                                                   t_i );
                                
                                // Boltzmann factors
                                Numeric nom = exp(- e_lower / ( BOLTZMAN_CONST * t_i ) ) - 
                                exp(- e_upper / ( BOLTZMAN_CONST * t_i ) );
                                
                                Numeric denom = exp(- e_lower / ( BOLTZMAN_CONST * l_l.Ti0() ) ) - 
                                exp(- e_upper / ( BOLTZMAN_CONST * l_l.Ti0() ) );
                                
                                
                                // intensity at temperature
                                // (calculate the line intensity according to the standard 
                                // expression given in HITRAN)
                                intensity *= part_fct_ratio * nom / denom;
                                
                                if (lineshape_norm_data[ind_lsn].Name() == "quadratic")
                                {
                                    // in case of the quadratic normalization factor use the 
                                    // so called 'microwave approximation' of the line intensity 
                                    // given by 
                                    // P. W. Rosenkranz, Chapter 2, Eq.(2.16), in M. A. Janssen, 
                                    // Atmospheric Remote Sensing by Microwave Radiometry, 
                                    // John Wiley & Sons, Inc., 1993
                                    Numeric mafac = (PLANCK_CONST * F0) / (2.000e0 * BOLTZMAN_CONST
                                    * t_i);
                                    intensity     = intensity * mafac / sinh(mafac);
                                }
                                
                                // 2. Calculate the pressure broadened line width and the pressure
                                // shifted center frequency.
                                //
                                // Here there is a difference betweeen catalogue version 4
                                // (from ESA planetary study) and earlier catalogue versions.
                                Numeric gamma;     // The line width.
                                Numeric deltaf=0;    // Pressure shift.
                                if (l_l.Version() == 4)
                                {
                                    calc_gamma_and_deltaf_artscat4_old_unused(gamma,
                                                                              deltaf,
                                                                              p_i,
                                                                              t_i,
                                                                              all_vmrs(joker,i),
                                                                              this_species,
                                                                              broad_spec_locations,
                                                                              l_l,
                                                                              verbosity);
                                }
                                else if (l_l.Version() == 3)
                                {
                                    //                  Numeric Tgam;
                                    //                  Numeric Nair;
                                    //                  Numeric Agam;
                                    //                  Numeric Sgam;
                                    //                  Numeric Nself;
                                    //                  Numeric Psf;
                                    
                                    //              if (l_l.Version() == 3)
                                    //              {
                                    const Numeric Tgam = l_l.Ti0(); // FIXME: This is breaking old catalogs??
                                    const Numeric Agam = l_l.Agam();
                                    const Numeric Nair = l_l.Nair();
                                    const Numeric Sgam = l_l.Sgam();
                                    const Numeric Nself = l_l.Nself();
                                    const Numeric Psf = l_l.Psf();
                                    //              }
                                    //              else
                                    //              {
                                    //                static bool warn = false;
                                    //                if (!warn)
                                    //                {
                                    //                  CREATE_OUT0;
                                    //                  warn = true;
                                    //                  out0 << "  WARNING: Using artscat version 4 for calculations is currently\n"
                                    //                       << "           just a hack and results are most likely wrong!!!\n";
                                    //                }
                                    //                Tgam = l_l.Tgam();
                                    //                // Use hardcoded mixing ratios for air
                                    //                Agam = l_l.Gamma_N2() * 0.79 + l_l.Gamma_O2() * 0.21;
                                    //                Nair = l_l.Gam_N_N2() * 0.79 + l_l.Gam_N_O2() * 0.21;
                                    //                Sgam = l_l.Gamma_self();
                                    //                Nself = l_l.Gam_N_self();
                                    //                Psf = l_l.Delta_N2() * 0.79 + l_l.Delta_O2() * 0.21;
                                    //              }
                                    
                                    // Get pressure broadened line width:
                                    // (Agam is in Hz/Pa, abs_p is in Pa, gamma is in Hz)
                                    const Numeric theta = Tgam / t_i;
                                    const Numeric theta_Nair = pow(theta, Nair);
                                    
                                    gamma =   Agam * theta_Nair  * (p_i - p_partial)
                                    + Sgam * pow(theta, Nself) * p_partial;
                                    
                                    
                                    // Pressure shift:
                                    // The T dependence is connected to that of agam by:
                                    // n_shift = .25 + 1.5 * n_agam
                                    // Theta has been initialized above.
                                    deltaf = Psf * p_i *
                                    std::pow( theta , (Numeric).25 + (Numeric)1.5*Nair );
                                    
                                }
                                else
                                {
                                    // There is a runtime error check for allowed artscat versions
                                    // further up.
                                    assert(false);
                                }
                                
                                // Apply pressure shift:
                                F0 += deltaf;
                                
                                // 3. Doppler broadening without the sqrt(ln(2)) factor, which
                                // seems to be redundant.
                                const Numeric sigma = F0 * doppler_const *
                                sqrt( t_i / l_l.IsotopologueData().Mass());
                                
                                //              cout << l_l.IsotopologueData().Name() << " " << l_l.F() << " "
                                //                   << Nair << " " << Agam << " " << Sgam << " " << Nself << endl;
                                
                                
                                
                                // Indices pointing at begin/end frequencies of f_grid or at
                                // the elements that have to be calculated in case of cutoff
                                Index i_f_min = 0;            
                                Index i_f_max = nf-1;         
                                
                                // cutoff ?
                                if ( cut )
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
                                    //               cout << ls << endl
                                    //                    << "aux / F0 / gamma / sigma" << aux << "/" << F0 << "/" << gamma << "/" << sigma << endl
                                    //                    << f_local[Range(i_f_min,nfls)] << endl
                                    //                    << nfls << endl;
                                    
                                    // Calculate the line shape:
                                    lineshape_data[ind_ls].Function()(ls_attenuation,ls_phase,
                                                                      aux,F0,gamma,sigma,
                                                                      f_local[Range(i_f_min,nfls)]);
                                    
                                    // Calculate the chosen normalization factor:
                                    lineshape_norm_data[ind_lsn].Function()(fac,F0,
                                                                            f_local[Range(i_f_min,nfls)],
                                                                            t_i);
                                    
                                    // Get a handle on the range of xsec that we want to change.
                                    // We use nfl here, which could be one less than nfls.
                                    VectorView this_xsec_attenuation      = xsec_accum_attenuation(arts_omp_get_thread_num(), Range(i_f_min,nfl));
                                    VectorView this_xsec_phase            = xsec_accum_phase(arts_omp_get_thread_num(), Range(i_f_min,nfl));
                                    
                                    // Get handles on the range of ls and fac that we need.
                                    VectorView this_ls_attenuation  = ls_attenuation[Range(0,nfl)];
                                    VectorView this_ls_phase  = ls_phase[Range(0,nfl)];
                                    VectorView this_fac = fac[Range(0,nfl)];
                                    
                                    // cutoff ?
                                    if ( cut )
                                    {
                                        // Subtract baseline for cutoff frequency
                                        // The index nfls-1 should be exactly the index pointing
                                        // to the value at the cutoff frequency.
                                        // Subtract baseline from xsec. 
                                        // this_xsec -= base;
                                        this_ls_attenuation -= ls_attenuation[nfls-1];
                                        //if (calc_phase) this_ls_phase -= ls_phase[nfls-1];  PHASE is not compatible with cutoff
                                    }
                                    
                                    // Add line to xsec. 
                                    {
                                        // To make the loop a bit faster, precompute all constant
                                        // factors. These are:
                                        // 1. Total number density of the air. --> Not
                                        //    anymore, we now to real cross-sections
                                        // 2. Line intensity.
                                        // 3. Isotopologue ratio.
                                        //
                                        // The isotopologue ratio must be applied here, since we are
                                        // summing up lines belonging to different isotopologues.
                                        
                                        //                const Numeric factors = n * intensity * l_l.IsotopologueData().Abundance();
                                        //                    const Numeric factors = intensity * l_l.IsotopologueData().Abundance();
                                        const Numeric factors = intensity
                                        * isotopologue_ratios.getParam(l_l.Species(), l_l.Isotopologue())[0].data[0];
                                        
                                        // We have to do:
                                        // xsec(j,i) += factors * ls[j] * fac[j];
                                        //
                                        // We use ls as a dummy to compute the product, then add it
                                        // to this_xsec.
                                        
                                        this_ls_attenuation *= this_fac;
                                        this_ls_attenuation *= factors;
                                        this_xsec_attenuation += this_ls_attenuation;
                                        
                                        if (calc_phase)
                                        {
                                            this_ls_phase *= this_fac;
                                            this_ls_phase *= factors;
                                            this_xsec_phase += this_ls_phase;
                                        }
                                    }
                                }
                                
                            } // end of try block
                            catch (runtime_error e)
                            {
                                #pragma omp critical (xsec_species_fail)
                                { fail_msg = e.what(); failed = true; }
                            }
                            
                        } // end of parallel LBL loop
                        
                        // Bail out if an error occurred in the LBL loop
                        if (failed) continue;
                        
                        // Now we just have to add up all the rows of xsec_accum:
                        for (Index j=0; j<xsec_accum_attenuation.nrows(); ++j)
                        {
                            xsec_i_attenuation += xsec_accum_attenuation(j, Range(joker));
                        }
                        
                        if (calc_phase)
                            for (Index j=0; j<xsec_accum_phase.nrows(); ++j)
                            {
                                xsec_i_phase += xsec_accum_phase(j, Range(joker));
                            }
            } // end of parallel pressure loop
            
            if (failed) throw std::runtime_error("Run-time error in function: xsec_species\n" + fail_msg);
            
}


/** Calculate line width and pressure shift for artscat4.
 * 
 *   \retval gamma Line width [Hz].
 *   \retval deltaf Pressure shift [Hz].
 *   \param  p Pressure [Pa].
 *   \param  t Temperature [K].
 *   \param  vmrs Vector of VMRs for different species [dimensionless].
 *   \param  this_species Index of current species in vmrs.
 *   \param  broad_spec_locations Has length of number of allowed broadening species
 *                                (6 in artscat-4). Gives for each species the position
 *                                in vmrs, or negative if it should be ignored. See 
 *                function find_broad_spec_locations for details.
 *   \param  l_l Spectral line data record (a single line).
 *   \param  verbosity Verbosity flag.
 * 
 *   \author Stefan Buehler
 *   \date   2012-09-05
 */
void calc_gamma_and_deltaf_artscat4_old_unused(Numeric& gamma,
                                               Numeric& deltaf,
                                               const Numeric p,
                                               const Numeric t,
                                               ConstVectorView vmrs,
                                               const Index this_species,
                                               const ArrayOfIndex& broad_spec_locations,
                                               const LineRecord& l_l,
                                               const Verbosity& verbosity)
{
    CREATE_OUT2;
    
    // Number of broadening species:
    const Index nbs = LineRecord::NBroadSpec();
    assert(nbs==broad_spec_locations.nelem());
    
    // Theta is reference temperature divided by local temperature. Used in
    // several place by the broadening and shift formula.
    const Numeric theta = l_l.Ti0() / t;
    
    // Split total pressure in self and foreign part:
    const Numeric p_self    = vmrs[this_species] * p;
    const Numeric p_foreign = p-p_self;
    
    // Calculate sum of VMRs of all available foreign broadening species (we need this
    // for normalization). The species "Self" will not be included in the sum!
    Numeric broad_spec_vmr_sum = 0;
    
    // Gamma is the line width. We first initialize gamma with the self width
    gamma =  l_l.Sgam() * pow(theta, l_l.Nself()) * p_self;
    
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
            if (l_l.Gamma_foreign(i)!=l_l.Sgam() ||
                l_l.N_foreign(i)!=l_l.Nself())
            {
                ostringstream os;
                os << "Inconsistency in LineRecord, self broadening and line "
                << "broadening for " << LineRecord::BroadSpecName(i) << "\n"
                << "should be identical.\n"
                << "LineRecord:\n"
                << l_l;
                throw std::runtime_error(os.str());
            }
        } else if ( broad_spec_locations[i] >= 0 ) {
            
            // Add to VMR sum:
            broad_spec_vmr_sum += vmrs[broad_spec_locations[i]];
            
            // foreign broadening:
            gamma_foreign +=  l_l.Gamma_foreign(i) * pow(theta, l_l.N_foreign(i))
            * vmrs[broad_spec_locations[i]];
            
            // Pressure shift:
            // The T dependence is connected to that of the corresponding
            // broadening parameter by:
            // n_shift = .25 + 1.5 * n_gamma
            deltaf += l_l.Delta_foreign(i)
            * pow( theta , (Numeric).25 + (Numeric)1.5*l_l.N_foreign(i) )
            * vmrs[broad_spec_locations[i]];
        }
    }
    
    // Check that sum of self and all foreign VMRs is not too far from 1:
    if ( abs(vmrs[this_species]+broad_spec_vmr_sum-1) > 0.1
        && out2.sufficient_priority() )
    {
        ostringstream os;
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