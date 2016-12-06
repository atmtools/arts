/* Copyright 2013, The ARTS Developers.
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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "arts.h"
#include "absorption.h"
#include "file.h"
#include "linemixingrecord.h"
#include "complex.h"
#include "lin_alg.h"
#include "global_data.h"
#include "Faddeeva.hh"
#include "linescaling.h"
#include "jacobian.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void line_mixing_dataInit(// WS Output:
                          ArrayOfArrayOfLineMixingRecord& line_mixing_data,
                          ArrayOfArrayOfIndex& line_mixing_data_lut,
                          // WS Input:
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          // Verbosity object:
                          const Verbosity&)
{
    line_mixing_data.resize(abs_species.nelem());
    line_mixing_data_lut.resize(abs_species.nelem());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void line_mixing_dataMatch(// WS Output:
                           //ArrayOfArrayOfLineMixingRecord& line_mixing_data,
                           //ArrayOfArrayOfIndex& line_mixing_data_lut,
                           // WS Output/Input:
                           ArrayOfArrayOfLineRecord& abs_lines_per_species,
                           // WS Input:
                           const ArrayOfArrayOfSpeciesTag& abs_species,
                           const String& species_tag,
                           const String& line_mixing_tag,
                           const ArrayOfLineMixingRecord& line_mixing_records,
                           const Verbosity& verbosity)
{
    CREATE_OUT2;
    CREATE_OUT3;
/*
    if (abs_species.nelem() != line_mixing_data.nelem())
        throw std::runtime_error( "*line_mixing_data* doesn't match *abs_species*.\n"
                             "Make sure to call line_mixing_dataInit first." );
    if (abs_species.nelem() != line_mixing_data_lut.nelem())
        throw std::runtime_error( "*line_mixing_data_lut* doesn't match *abs_species*.\n"
                            "Make sure to call line_mixing_dataInit first." );*/
    
    LineMixingData line_mixing_data_holder;
    line_mixing_data_holder.StorageTag2SetType(line_mixing_tag);

    // Find index of species_tag in abs_species
    SpeciesTag this_species( species_tag );
    Index species_index = -1;
    for (Index i = 0; species_index == -1 && i < abs_species.nelem(); i++)
        // We only look at the first SpeciesTag in each group because
        // line mixing tags can not be combined with other tags
        if (this_species == abs_species[i][0])
            species_index = i;

    if (species_index == -1)
    {
        ostringstream os;
        os << "Can't find tag \"" << species_tag << "\" in *abs_species*.";
        throw std::runtime_error(os.str());
    }
ArrayOfArrayOfLineMixingRecord line_mixing_data(abs_species.nelem());//FIXME: Added for convenience to match old use case
    line_mixing_data[species_index] = line_mixing_records;

    ArrayOfIndex matches;

    // Now we use the quantum numbers to match the line mixing
    // data to lines in abs_lines_per_species
    Index nmatches = 0;
    for (Index i = 0; i < line_mixing_data[species_index].nelem(); i++)
    {
        const LineMixingRecord& this_lmr = line_mixing_data[species_index][i];
        
        find_matching_lines(matches,
                            abs_lines_per_species[species_index],
                            this_lmr.Species(),
                            this_lmr.Isotopologue(),
                            this_lmr.Quantum());
        
        line_mixing_data_holder.SetDataFromVectorWithKnownType(line_mixing_data[species_index][i].Data());
        
        if (!matches.nelem())
        {
            out3 << "  Found no matching lines for\n" << this_lmr.Quantum() << "\n";
        }
        else if (matches.nelem() == 1)
        {
            out3 << "  Found matching line for\n" << this_lmr.Quantum() << "\n";
            abs_lines_per_species[species_index][matches[0]].SetLineMixingData(line_mixing_data_holder);
            nmatches++;
        }
        else
        {
            ostringstream os;
            os << "  Found multiple lines for\n" << this_lmr.Quantum() << std::endl
            << "  Matching lines are: " << std::endl;
            for (Index m = 0; m < matches.nelem(); m++)
                os << "  " << abs_lines_per_species[species_index][matches[m]] << std::endl
                << "  " << abs_lines_per_species[species_index][matches[m]].QuantumNumbers()
                << std::endl;
            throw std::runtime_error(os.str());
        }
    }

    out2 << "  Matched " << nmatches << " lines out of " << line_mixing_data[species_index].nelem()
         << "\n";
    out2 << "  abs_lines_per_species contains " << abs_lines_per_species[species_index].nelem()
         << " lines for " << species_tag << ".\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfLineMixingRecordReadAscii(// Generic Output:
                                      ArrayOfLineMixingRecord& line_mixing_records,
                                      // Generic Input:
                                      const String& filename,
                                      const Verbosity& verbosity)
{
    ifstream ifs;
    open_input_file(ifs, filename);

    line_mixing_records.resize(0);

    // Read the line mixing file
    Index linenr = 0;
    String line;
    istringstream is;
    while (!ifs.eof())
    {
        getline(ifs, line);
        linenr++;

        line.trim();
        // Skip empty lines and comments
        if (!line.nelem()
            || (line.nelem() && line[0] == '#'))
            continue;

        is.clear();
        is.str(line);
        try {
            String species_string;
            SpeciesTag line_species;

            is >> species_string;
            line_species = SpeciesTag(species_string);

            LineMixingRecord lmr(line_species.Species(), line_species.Isotopologue());

            Rational r;
            is >> r; lmr.Quantum().SetLower(QN_v1, r);
            lmr.Quantum().SetUpper(QN_v1, r);
            is >> r; lmr.Quantum().SetUpper(QN_N,  r);
            is >> r; lmr.Quantum().SetLower(QN_N,  r);
            is >> r; lmr.Quantum().SetUpper(QN_J,  r);
            is >> r; lmr.Quantum().SetLower(QN_J,  r);

            vector<Numeric> temp_mixing_data;
            String s;
            char *c;
            while (is.good() && is)
            {
                is >> s;
                s.trim();
                if (s.nelem())
                {
                    temp_mixing_data.push_back(strtod(s.c_str(), &c));
                    if (c != s.c_str() + s.nelem())
                        throw std::runtime_error(line);
                }
            }

            lmr.Data() = temp_mixing_data;
            line_mixing_records.push_back(lmr);
        } catch (runtime_error e) {
            ostringstream os;

            os << "Error parsing line mixing file in line " << linenr << std::endl;
            os << e.what();
            throw std::runtime_error(os.str());
        }
    }

    CREATE_OUT2;
    out2 << "  Read " << line_mixing_records.nelem() << " lines from " << filename << ".\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_bandFromband_identifiers( ArrayOfArrayOfLineRecord&       abs_lines_per_band,
                                             ArrayOfArrayOfSpeciesTag&       abs_species_per_band,
                                             ArrayOfArrayOfLineRecord&       abs_lines_per_species,
                                             const ArrayOfArrayOfSpeciesTag& abs_species,
                                             const ArrayOfQuantumIdentifier& band_identifiers,
                                             const Verbosity&                verbosity)
{
  CREATE_OUT3;
  out3<<"Sets line mixing tag for provided bands.\n" <<
  "\tNB. Requires \"*-LM-*\" tag in abs_species.";
  
  abs_lines_per_band.resize(band_identifiers.nelem());
  abs_species_per_band.resize(band_identifiers.nelem());
  
  // This is a preallocated finding of a band
  LineMixingData lmd_byband;
  lmd_byband.SetByBand(); 
  
  #pragma omp parallel for        \
  if (!arts_omp_in_parallel())    
    for (Index qi = 0; qi < band_identifiers.nelem(); qi++)
    {
      const QuantumIdentifier& band_id = band_identifiers[qi];
      
      // Two variables that are used inside the loop
      ArrayOfIndex matches;
      ArrayOfQuantumMatchInfo match_info;
      
      for (Index s = 0; s < abs_lines_per_species.nelem(); s++)
      {
        
        // Skip this species if qi is not part of the species represented by this abs_lines
        if(abs_species[s][0].Species() != band_id.Species() || abs_species[s][0].LineMixing() == SpeciesTag::LINE_MIXING_OFF)
          continue;
        
        ArrayOfLineRecord& species_lines = abs_lines_per_species[s];
        
        // Copy this
        abs_species_per_band[qi] = abs_species[s];
        
        // Run internal mathcing routine
        match_lines_by_quantum_identifier(matches, match_info, band_id, species_lines);
        
        // Use info about mathced lines to tag the relevant parameter
        for (Index i = 0; i < matches.nelem(); i++)
        {
          QuantumMatchInfo& qm = match_info[i];
          
          LineRecord& lr = species_lines[matches[i]];
          
          // If any of the levels match partially or fully set the right quantum number
          if(qm.Upper()==QMI_NONE||qm.Lower()==QMI_NONE)
          {
            continue;
          }
          else // we will accept this match if both levels are at least partially matched
          {
            abs_lines_per_band[qi].push_back(lr);
            lr.SetLineMixingData(lmd_byband);
          }
        }
      }
    }
}


void calculate_xsec_from_W( VectorView  xsec,
                            ConstMatrixView Wmat,
                            ConstVectorView f_grid,
                            ConstVectorView f0,
                            ConstVectorView d0,
                            ConstVectorView rhoT,
                            const Numeric&  T,
                            const Numeric&  P,
                            const Numeric&  isotopologue_mass,
                            const Index&    n//lines
)
{
    // Physical constants
    extern const Numeric BOLTZMAN_CONST;
    extern const Numeric AVOGADROS_NUMB;
    extern const Numeric SPEED_OF_LIGHT;;
    
    // internal constant
    const Index nf  = f_grid.nelem();
    const Index nf0 = f0.nelem();
    //const Numeric f_mean = mean(f0);
    const Numeric doppler_const = sqrt( 2.0 * BOLTZMAN_CONST * AVOGADROS_NUMB * T / isotopologue_mass ) / SPEED_OF_LIGHT;
    
    // Setuo for the matrix
    ComplexMatrix W(n,n);
    for(Index i=0;i<n;i++)
    {
        for(Index j=0;j<n;j++)
        {
            if(j==i)
                W(i,j)=Complex(0,-P*Wmat(i,j)-f0[i]);
            else
                W(i,j)=Complex(0,-P*Wmat(i,j));
        }
    }
    
    // Setup so that W above is equivalent to D*diag(z_eigs)*invD, where diag
    // is a diagonal matrix with values of eigs in the diagonal
    ComplexVector z_eigs(nf0);
    ComplexMatrix D(nf0,nf0),invD(nf0,nf0);
    diagonalize(D,z_eigs,W);
    inv(invD,D);
    // Question:  If this fails and produce baloney, is this caught later on?  Switch to zgeevx_?
    
    // Equivalent line strength
    ComplexVector equivS0(nf0,0);
    
    // Doppler broadening
    Vector sigma(nf0);
    
    // Set starts and equivs
    for(Index if1=0;if1<nf0;if1++)
    {
        sigma[if1]= f0[if1] * doppler_const ;
        z_eigs[if1]/=sigma[if1];
        
        for(Index if2=0;if2<nf0;if2++)
        {
            equivS0[if2] += rhoT[if1]*d0[if1]*d0[if2]*D(if2,if1)*invD(if1,if2);
        }
    }
    
    // Set xsec (need to normalize?)
    for(Index if0=0;if0<nf0;if0++)
        for(Index iv=0;iv<nf;iv++)
            xsec+=((equivS0[if0]/sigma[if0])*Faddeeva::w(z_eigs[if0]+f_grid[iv]/sigma[if0])).real();
}

#ifdef ENABLE_RELMAT
extern "C"
{
    // This is the interface between the Fortran code that calculates W and ARTS
    extern void arts_relmat_interface(
        long   *nlines,
        double *fmin,
        double *fmax,
        long   *M,
        long   *I,
        double *v,
        double *S,
        double *gamma_air,
        double *E_double_prime,
        double *n_air,
        long   *upper,
        long   *lower,
        long   *g_prime,
        long   *g_double_prime,
        double *temperature,
        double *pressure,
        double *partition_function_t,
        double *partition_function_t0,
        double *isotopologue_mass,
        long   *number_of_perturbers,
        long   *molecule_code_perturber,
        long   *iso_code_perturber,
        double *vmr,
        //outputs
        double *W,
        double *d0,
        double *rhoT
    );
}


void abs_xsec_per_speciesAddLineMixedBands( // WS Output:
                                            ArrayOfMatrix&                   abs_xsec_per_species,
                                            ArrayOfArrayOfMatrix&            /*dabs_xsec_per_species_dx*/,
                                            // WS Input:                     
                                            const ArrayOfArrayOfLineRecord&  abs_lines_per_band,
                                            const ArrayOfArrayOfSpeciesTag&  abs_species_per_band,
                                            const ArrayOfArrayOfSpeciesTag&  abs_species,
                                            const SpeciesAuxData&            partition_functions,
                                            const ArrayOfRetrievalQuantity&  jacobian_quantities,
                                            const Vector&                    f_grid,
                                            const Vector&                    abs_p,
                                            const Vector&                    abs_t,
                                            const Verbosity&)
{
    /*throw std::runtime_error("\nabs_xsec_per_speciesAddLineMixedBands of src/m_linemixing.cc\n"
    "is still a work in progress and does not work as of yet.\n"
    "Please remove this runtime_error to proceed with debugging the reasons.\n");*/
    
    using global_data::species_data;
    
    // Physical constants
    extern const Numeric PLANCK_CONST;
    extern const Numeric SPEED_OF_LIGHT;
    extern const Numeric ATM2PA;
    
    // HITRAN to ARTS constants
    const Numeric w2Hz               = SPEED_OF_LIGHT *1E2;
    const Numeric lower_energy_const = PLANCK_CONST * w2Hz;
    const Numeric I0_hi2arts         = 1E-2 * SPEED_OF_LIGHT;
    const Numeric gamma_hi2arts      = w2Hz / ATM2PA;
    
    // size of atmosphere and input/output
    const Index nps         = abs_p.nelem();
    const Index nts         = abs_t.nelem();
    const Index nf          = f_grid.nelem();
    const Index nspecies    = abs_species.nelem();
    
    // These should be identical
    if(nps!=nts)
        throw std::runtime_error("Different lengths of atmospheric input than expected.");
    
    if(nspecies==0 || nspecies!=abs_xsec_per_species.nelem())
        throw std::runtime_error("Absorption species and abs_xsec_per_species are not from same source.");
    else
        for(Index i=0; i<nspecies; i++)
            if(abs_xsec_per_species[i].ncols()!=nts)
                throw std::runtime_error("Unexpected size of xsec matrix not matching abs_t and abs_p length.");
            else if(abs_xsec_per_species[i].nrows()!=nf)
                throw std::runtime_error("Unexpected size of xsec matrix not matching f_grid length.");
    
    // FIXME: the partial derivations are necessary...
    if( 0!=jacobian_quantities.nelem() )
      throw std::runtime_error("Presently no support for partial derivation.");
    
    const Index nbands = abs_lines_per_band.nelem();
    
    if(nbands!=abs_species_per_band.nelem())
        throw std::runtime_error("Error in definition of the bands.  Mismatching length of *_per_band arrays.");
    
    // Make constant input not constant and convert to wavenumber
    Vector v(nf); 
    v=f_grid; 
    v/=w2Hz;
    
    // Setting up thermal bath:  only in simplistic air for now
    // This means: 21% O2 and 79% N2
    //
    long    number_of_perturbers = 2;
    long   *molecule_code_perturber = new long[number_of_perturbers];
    long   *iso_code_perturber = new long[number_of_perturbers];
    Vector  vmr(number_of_perturbers);
    bool    done_o2=false, done_n2=false;
    for(Index ispecies=0;ispecies<species_data.nelem();ispecies++)
    {
        const SpeciesRecord& sr = species_data[ispecies];
        const IsotopologueRecord& ir = sr.Isotopologue()[0];
        const String& name = sr.Name();
        
        if(name=="O2"&&!done_o2)
        {
            vmr[0] = 0.21;
            const Index hitran_tag = ir.HitranTag();
            iso_code_perturber[0] = (long) (hitran_tag%10);
            molecule_code_perturber[0] = (((long)hitran_tag)-iso_code_perturber[0])/10;
            done_o2=true;
        }
        else if(name=="N2"&&!done_n2)
        {
            vmr[1] = 0.79;
            const Index hitran_tag = ir.HitranTag();
            iso_code_perturber[1] = (long) (hitran_tag%10);
            molecule_code_perturber[1] = (((long)hitran_tag)-iso_code_perturber[1])/10;
            done_n2=true;
        }
        
        if(done_n2&&done_o2)
            break;
    }
    
    // FIXME:  Can this loop be parallelized?
    for(Index iband=0;iband<nbands;iband++)
    {
        // band pointer
        const ArrayOfLineRecord& this_band = abs_lines_per_band[iband];
        
        long nlines = (long) this_band.nelem();
        
        // Worth doing anything?
        if(nlines==0) { continue; }
        
        // Send in frequency range
        Numeric fmin, fmax;
        
        // To store the xsec matrix we need to know where
        Index this_species=-1;
        
        // Allocation of band specific inputs (types: to be used as fortran input)
        long   *M              = new long[nlines];
        long   *I              = new long[nlines];
        long   *upper          = new long[4*nlines];
        long   *lower          = new long[4*nlines];
        long   *g_prime        = new long[nlines];
        long   *g_double_prime = new long[nlines];
        Vector v0(nlines);
        Vector S(nlines);
        Vector gamma_air(nlines);
        Vector E_double_prime(nlines);
        Vector n_air(nlines);
        Numeric mass, abundance;
        
        for( long iline=0; iline<nlines; iline++ )
        {
            // Line data
            const LineRecord& this_line = this_band[iline];
            const IsotopologueRecord& this_iso = this_line.IsotopologueData();
            
            // For first line do something special
            if( iline==0 )
            {
                // Isotopologue values
                mass  = this_iso.Mass();
                abundance = this_iso.Abundance();
                String iso_name = this_iso.Name();
                int isona;
                //isona << iso_name;
                extract(isona,iso_name,iso_name.nelem());
                
                const ArrayOfSpeciesTag& band_tags=abs_species_per_band[iband];
                
                // Finds the first species in abs_species_per_band that matches abs_species
                bool this_one=false;
                for(Index ispecies=0; ispecies<nspecies; ispecies++)
                {
                    const ArrayOfSpeciesTag& species_tags=abs_species[ispecies];
                    const Index nbandtags = band_tags.nelem();
                    const Index nspeciestags = species_tags.nelem();
                    
                    // Test if there is
                    if(nbandtags!=nspeciestags) { break; }
                    for(Index itags=0; itags<nspeciestags; itags++)
                    {
                        if(band_tags[itags]==species_tags[itags])
                        {
                            this_one = true;
                            break;
                        }
                    }
                    
                    if(this_one)
                    {
                        this_species = ispecies;
                        break;
                    }
                    
                }
                
                if(!this_one)
                    throw std::runtime_error("abs_species and abs_species_per_band disagrees"
                    " on absorption species");
                
            }
            else 
            {
                if( mass!=this_iso.Mass() )
                {
                    throw std::runtime_error("There are lines of different Isotopologues"
                    " in abs_lines_per_band,");
                }
            }
            
            // Hitran tags
            long hitran_tag       = (long) this_iso.HitranTag();
            
            // Line information converted to relmat format --- i.e. to HITRAN format
            I[iline]              = hitran_tag%10;
            M[iline]              = (hitran_tag-I[iline])/10;
            v0[iline]             = this_line.F()/w2Hz*abundance; // WARNING:  Necessity of abundance means vmr should also be scaled?
            S[iline]              = this_line.I0()/I0_hi2arts;
            gamma_air[iline]      = this_line.Agam()/gamma_hi2arts;
            n_air[iline]          = this_line.Nair();
            E_double_prime[iline] = this_line.Elow()/lower_energy_const;
            g_prime[iline]        = (long) this_line.G_upper(); // NB:  Numeric to long... why?
            g_double_prime[iline] = (long) this_line.G_lower(); // NB:  Numeric to long... why?
            
            // Quantum numbers converted to relmat format, again Numeric/Rational to long... why?
            Rational a;
            
            // l2 is for molecules like CO2
            a = this_line.QuantumNumbers().Lower()[QN_l2];
            a.Simplify();
            if(a.isUndefined())
                lower[0+4*iline] = -1;
            else  if(a.Denom()==1)
                lower[0+4*iline] = (long) a.toIndex();
            else 
                throw std::runtime_error("Half quantum numbers not supported in l2.");
            a = this_line.QuantumNumbers().Upper()[QN_l2];
            a.Simplify();
            if(a.isUndefined())
                upper[0+4*iline] = -1;
            else  if(a.Denom()==1)
                upper[0+4*iline] = (long) a.toIndex();
            else 
                throw std::runtime_error("Half quantum numbers not supported in l2.");
            
            // J is universally important for linear molecules
            a = this_line.QuantumNumbers().Lower()[QN_J];
            a.Simplify();
            if(a.isUndefined())
                lower[1+4*iline] = -1;
            else  if(a.Denom()==1)
                lower[1+4*iline] = (long) a.toIndex();
            else 
                throw std::runtime_error("Half quantum numbers not supported in J.");
            a = this_line.QuantumNumbers().Upper()[QN_J];
            a.Simplify();
            if(a.isUndefined())
                upper[1+4*iline] = -1;
            else  if(a.Denom()==1)
                upper[1+4*iline] = (long) a.toIndex();
            else 
                throw std::runtime_error("Half quantum numbers not supported in J.");
            
            // N is important for molecules with magnetic dipoles
            a = this_line.QuantumNumbers().Lower()[QN_N];
            a.Simplify();
            if(a.isUndefined())
                lower[2+4*iline] = -1;
            else  if(a.Denom()==1)
                lower[2+4*iline] = (long) a.toIndex();
            else 
                throw std::runtime_error("Half quantum numbers not supported in N.");
            a = this_line.QuantumNumbers().Upper()[QN_N];
            a.Simplify();
            if(a.isUndefined())
                upper[2+4*iline] = -1;
            else  if(a.Denom()==1)
                upper[2+4*iline] = (long) a.toIndex();
            else 
                throw std::runtime_error("Half quantum numbers not supported in N.");
            
            // S is important for molecules with magnetic dipoles
            a = this_line.QuantumNumbers().Lower()[QN_S];
            a.Simplify();
            if(a.isUndefined())
                lower[3+4*iline] = -1;
            else  if(a.Denom()==1)
                lower[3+4*iline] = (long) a.toIndex();
            else 
                throw std::runtime_error("Half quantum numbers not supported in S.");
            a = this_line.QuantumNumbers().Upper()[QN_S];
            a.Simplify();
            if(a.isUndefined())
                upper[3+4*iline] = -1;
            else  if(a.Denom()==1)
                upper[3+4*iline] = (long) a.toIndex();
            else 
                throw std::runtime_error("Half quantum numbers not supported in S.");
            
            // Set fmax and fmin.  Why do I need this again?
            if(iline==0)
            {
                fmin=v0[0];
                fmax=v0[0];
            }
            else 
            {
                if(fmin>v0[iline])
                    fmin=v0[iline];
                if(fmax<v0[iline])
                    fmax=v0[iline];
            }
            
        }
        
        // FIXME:  Or can this loop be parallelized?
        for(Index ip=0;ip<nps;ip++)
        {
            // Information on the lines will be here after relmat is done
            Matrix W(nlines,nlines);
            Vector d0(nlines);
            Vector rhoT(nlines);
            
            // Find the partition function
            Numeric QT0;
            Numeric QT;
            
            // Get partition function information
            partition_function( QT0,
                                QT,
                                abs_lines_per_band[iband][0].Ti0(),
                                abs_t[ip],
                                partition_functions.getParamType(abs_lines_per_band[iband][0].Species(), abs_lines_per_band[iband][0].Isotopologue()),
                                partition_functions.getParam(abs_lines_per_band[iband][0].Species(), abs_lines_per_band[iband][0].Isotopologue()),
                                false);
            
            // Cannot be constants for Fortran's sake
            Numeric t;
            t = abs_t[ip];
            Numeric p;
            p = abs_p[ip];
            
            std::cout<<"Starting the arts_relmat_interface!\n";
            
            // Calling Teresa's code
            arts_relmat_interface(
                &nlines, &fmin, &fmax,
                M, I, v0.get_c_array(), S.get_c_array(),
                gamma_air.get_c_array(),E_double_prime.get_c_array(),n_air.get_c_array(),
                upper, lower,
                g_prime, g_double_prime,
                &t, &p, &QT, &QT0, &mass,
                &number_of_perturbers, molecule_code_perturber, 
                iso_code_perturber, vmr.get_c_array(),
                W.get_c_array(), d0.get_c_array(), rhoT.get_c_array() );
            
            std::cout<<"Succesful run of arts_relmat_interface!\n";
            
            // Using Rodrigues method
            calculate_xsec_from_W( abs_xsec_per_species[this_species](joker, ip),
                                   W,
                                   v0,
                                   v,
                                   d0,
                                   rhoT,
                                   t,
                                   p,
                                   mass,
                                   nlines);
            
        }
        
        delete[] M;
        delete[] I;
        delete[] g_prime;
        delete[] g_double_prime;
        delete[] upper;
        delete[] lower;
    }
    
    delete[] iso_code_perturber;
    delete[] molecule_code_perturber;
    
}  

#else
void abs_xsec_per_speciesAddLineMixedBands( // WS Output:
                                            ArrayOfMatrix&                   /* abs_xsec_per_species, */,
                                            ArrayOfArrayOfMatrix&            /*dabs_xsec_per_species_dx*/,
                                            // WS Input:                     
                                            const ArrayOfArrayOfLineRecord&  /* abs_lines_per_band */,
                                            const ArrayOfArrayOfSpeciesTag&  /* abs_species_per_band */,
                                            const ArrayOfArrayOfSpeciesTag&  /* abs_species */,
                                            const SpeciesAuxData&            /* partition_functions */,
                                            const ArrayOfRetrievalQuantity&  /* jacobian_quantities */,
                                            const Vector&                    /* f_grid */,
                                            const Vector&                    /* abs_p */,
                                            const Vector&                    /* abs_t */,
                                            const Verbosity&)
{
    throw std::runtime_error("This version of ARTS was compiled without external line mixing support.");
}

#endif //ENABLE_RELMAT
