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
#include <map>
#include <cfloat>
#include <algorithm>
#include <cmath>
#include "file.h"
#include "absorption.h"
#include "math_funcs.h"
#include "messages.h"
#include "logic.h"


/** The map associated with species_data. */
std::map<String, Index> SpeciesMap;


template<class T>
void extract(T& x, String& line, Index  n);


// member fct of isotopologuerecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopologueRecord::CalculatePartitionFctAtTemp( Numeric
                                                    temperature ) const
{
  Numeric result = 0.;
  Numeric exponent = 1.;

  ArrayOfNumeric::const_iterator it;

//      cout << "T: " << temperature << "\n";
  for (it=mqcoeff.begin(); it != mqcoeff.end(); it++)
    {
      result += *it * exponent;
      exponent *= temperature;
//      cout << "it: " << it << "\n";
//      cout << "res: " << result << ", exp: " << exponent << "\n";
    }
  return result;
}


void SpeciesAuxData::initParams(Index nparams)
{
    extern const Array<SpeciesRecord> species_data;

    mparams.resize(species_data.nelem());

    for (Index isp = 0; isp < species_data.nelem(); isp++)
    {
        mparams[isp].resize(species_data[isp].Isotopologue().nelem(), nparams);
        mparams[isp] = NAN;
    }
}


bool SpeciesAuxData::ReadFromStream(String& artsid, istream& is, Index nparams, const Verbosity& verbosity)
{
    CREATE_OUT3;

    // Global species lookup data:
    extern const Array<SpeciesRecord> species_data;

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

        Matrix& params = mparams[mspecies];
        // Extract accuracies:
        try
        {
            Numeric p;
            for (Index ip = 0; ip < nparams; ip++)
            {
                icecream >> double_imanip() >> p;
                params(misotopologue, ip) = p;
            }
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
    extern const Array<SpeciesRecord> species_data;
    
    // Check total number of species:
    if (species_data.nelem() != isoratios.getParams().nelem())
      {
        ostringstream os;
        os << "Number of species in SpeciesAuxData (" << isoratios.getParams().nelem()
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
        if (this_sd.Isotopologue().nelem() != isoratios.getParams()[sp].nrows())
          {
            ostringstream os;
            os << "Incorrect number of isotopologues in isotopologue data.\n"
            << "Species: " << this_sd.Name() << ".\n"
            << "Number of isotopes in SpeciesAuxData ("
            << isoratios.getParams()[i].nrows() << ") "
            << "does not fit builtin species data (" << this_sd.Isotopologue().nelem() << ").";
            throw runtime_error(os.str());
          }
        
        for (Index iso = 0; iso < this_sd.Isotopologue().nelem(); ++iso)
          {
            // For "real" species (not representing continau) the isotopologue
            // ratio must not be NAN or below zero.
            if (!this_sd.Isotopologue()[iso].isContinuum()) {
                if (isnan(isoratios.getParam(sp, iso, 0)) ||
                    isoratios.getParam(sp, iso, 0) < 0.) {
                    
                    ostringstream os;
                    os << "Invalid isotopologue ratio.\n"
                    << "Species: " << this_sd.Name() << "-"
                    << this_sd.Isotopologue()[iso].Name() << "\n"
                    << "Ratio:   " << isoratios.getParam(sp, iso, 0);
                    throw runtime_error(os.str());
                }
            }
          }
      }
}


void fillSpeciesAuxDataWithIsotopologueRatiosFromSpeciesData(SpeciesAuxData& sad)
{
    extern const Array<SpeciesRecord> species_data;

    sad.initParams(1);

    for (Index isp = 0; isp < species_data.nelem(); isp++)
        for (Index iiso = 0; iiso < species_data[isp].Isotopologue().nelem(); iiso++)
        {
            sad.setParam(isp, iiso, 0, species_data[isp].Isotopologue()[iiso].Abundance());
        }
}


/*! Define the species data map.

    \author Stefan Buehler */
void define_species_map()
{
  extern const Array<SpeciesRecord> species_data;

  for ( Index i=0 ; i<species_data.nelem() ; ++i)
    {
      SpeciesMap[species_data[i].Name()] = i;
    }
}


ostream& operator<< (ostream& os, const LineRecord& lr)
{
  // Determine the precision, depending on whether Numeric is double
  // or float:  
  int precision;
#ifdef USE_FLOAT
  precision = FLT_DIG;
#else
#ifdef USE_DOUBLE
  precision = DBL_DIG;
#else
#error Numeric must be double or float
#endif
#endif

  switch (lr.Version()) {
    case 3:
      os << "@"
      << " " << lr.Name  ()
      << " "
      << setprecision(precision)
      <<        lr.F     ()
      << " " << lr.Psf   ()
      << " " << lr.I0    ()
      << " " << lr.Ti0   ()
      << " " << lr.Elow  ()
      << " " << lr.Agam  ()
      << " " << lr.Sgam  ()
      << " " << lr.Nair  ()
      << " " << lr.Nself ()
      << " " << lr.Tgam  ()
      << " " << lr.Naux  ()
      << " " << lr.dF    ()
      << " " << lr.dI0   ()
      << " " << lr.dAgam ()
      << " " << lr.dSgam ()
      << " " << lr.dNair ()
      << " " << lr.dNself()
      << " " << lr.dPsf  ();

      // Added new lines for the spectroscopic parameters accuracies.
      for ( Index i=0; i<lr.Naux(); ++i )
        os << " " << lr.Aux()[i];
      
      break;
      
    case 4: {
      ostringstream ls;

      ls << "@"
         << " " << lr.Name  ()
         << " "
         << setprecision(precision)
         <<        lr.F()
         << " " << lr.I0()
         << " " << lr.Ti0()
         << " " << lr.Elow()
         << " " << lr.A()
         << " " << lr.G_upper()
         << " " << lr.G_lower()
         << " " << lr.Sgam();
      
      for (Index s=0; s<6; ++s)
        ls << " " << lr.Gamma_foreign(s);
//      << " " << lr.Gamma_foreign(LineRecord::SPEC_POS_N2)
//      << " " << lr.Gamma_foreign(LineRecord::SPEC_POS_O2)
//      << " " << lr.Gamma_foreign(LineRecord::SPEC_POS_H2O)
//      << " " << lr.Gamma_foreign(LineRecord::SPEC_POS_CO2)
//      << " " << lr.Gamma_foreign(LineRecord::SPEC_POS_H2)
//      << " " << lr.Gamma_foreign(LineRecord::SPEC_POS_He)

      ls << " " << lr.Nself();
      for (Index s=0; s<6; ++s)
        ls << " " << lr.N_foreign(s);
//      << " " << lr.N_foreign(LineRecord::SPEC_POS_N2)
//      << " " << lr.N_foreign(LineRecord::SPEC_POS_O2)
//      << " " << lr.N_foreign(LineRecord::SPEC_POS_H2O)
//      << " " << lr.N_foreign(LineRecord::SPEC_POS_CO2)
//      << " " << lr.N_foreign(LineRecord::SPEC_POS_H2)
//      << " " << lr.N_foreign(LineRecord::SPEC_POS_He)

      for (Index s=0; s<6; ++s)
        ls << " " << lr.Delta_foreign(s);
//      << " " << lr.Delta_foreign(LineRecord::SPEC_POS_N2)
//      << " " << lr.Delta_foreign(LineRecord::SPEC_POS_O2)
//      << " " << lr.Delta_foreign(LineRecord::SPEC_POS_H2O)
//      << " " << lr.Delta_foreign(LineRecord::SPEC_POS_CO2)
//      << " " << lr.Delta_foreign(LineRecord::SPEC_POS_H2)
//      << " " << lr.Delta_foreign(LineRecord::SPEC_POS_He)
      ls << " " << lr.Upper_GQuanta()
         << " " << lr.Lower_GQuanta()
         << " " << lr.Upper_LQuanta()
         << " " << lr.Lower_LQuanta();

      os << ls.str();
      
      break;
    }
      
    default:
      os << "Unknown ARTSCAT version: " << lr.Version();
      break;
  }
  
  return os;
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
  os << sad.getParams();
  return os;
}


/** Extract something from a catalogue line. This is just a small helper
    function to safe some typing. 

    \retval x    What was extracted from the beginning of the line.
    \retval line What was extracted is also cut away from line.
    \param n     The width of the stuff to extract.

    \author Stefan Buehler */
template<class T>
void extract(T&      x,
             String& line,
             Index  n)
{
  // Initialize output to zero! This is important, because otherwise
  // the output variable could `remember' old values.
  x = T(0);
  
  // This will contain the short subString with the item to extract.
  // Make it a String stream, for easy parsing,
  // extracting subString of width n from line:
  istringstream item( line.substr(0,n) );

//  cout << "line = '" << line << "'\n";
//   cout << "line.substr(0,n) = " << line.substr(0,n) << endl;
//   cout << "item = " << item.str() << endl;

  // Shorten line by n:
  line.erase(0,n);
//  cout << "line = " << line << endl;

  // Convert with the aid of String stream item:
  item >> x;
}

/** The full name of the species and isotopologue. E.g., `O3-666'. The
  name is found by looking up the information in species_data,
  using the species and isotopologue index. */
String LineRecord::Name() const {
  // The species lookup data:
  extern const Array<SpeciesRecord> species_data;
  const SpeciesRecord& sr = species_data[mspecies];
  return sr.Name() + "-" + sr.Isotopologue()[misotopologue].Name();
}


/** The matching SpeciesRecord from species_data. To get at the
    species data of a LineRecord lr, you can use:
    <ul>
    <li>species_data[lr.Species()]</li>
    <li>lr.SpeciesData()</li>
    </ul>
    The only advantages of the latter are that the notation is
    slightly nicer and that you don't have to declare the external
    variable species_data. */
const SpeciesRecord& LineRecord::SpeciesData() const {
  // The species lookup data:
  extern const Array<SpeciesRecord> species_data;
  return species_data[mspecies];
}


/** The matching IsotopologueRecord from species_data. The IsotopologueRecord
    is a subset of the SpeciesRecord. To get at the isotopologue data of
    a LineRecord lr, you can use:
    <ul>
    <li>species_data[lr.Species()].Isotopologue()[lr.Isotopologue()]</li>
    <li>lr.SpeciesData().Isotopologue()[lr.Isotopologue()]</li>
    <li>lr.IsotopologueData()</li>
    </ul>
    The last option is clearly the shortest, and has the advantage
    that you don't have to declare the external variable
    species_data. */
const IsotopologueRecord& LineRecord::IsotopologueData() const {
  // The species lookup data:
  extern const Array<SpeciesRecord> species_data;
  return species_data[mspecies].Isotopologue()[misotopologue];
}


void LineRecord::ARTSCAT4FromARTSCAT3() {

    // Skip this line if it is already ARTSCAT-4
    if ( this->Version() == 4 ) return;

    // Check that this line really is ARTSCAT-3
    if ( this->Version() != 3 )
    {
        ostringstream os;
        os << "This line is not ARTSCAT-3, it is ARTSCAT-" << this->Version();
        throw runtime_error(os.str());
    }

    // Set version to 4:
    mversion = 4;

    const Index nbs = NBroadSpec();

    // Resize foreign parameter arrays:
    mgamma_foreign.resize(nbs);
    mn_foreign.resize(nbs);
    mdelta_foreign.resize(nbs);

    // Loop over broadening species:
    for (Index i=0; i<nbs; ++i) {

        // Find out if this broadening species is identical to the line species:
        if (this->Species() == BroadSpecSpecIndex(i)) {
            // We have to copy the self parameters here.
            mgamma_foreign[i] = msgam;
            mn_foreign[i] =     mnself;
            mdelta_foreign[i] = 0;
        } else {
            // We have to copy the foreign parameters here.
            mgamma_foreign[i] = magam;
            mn_foreign[i] =     mnair;
            mdelta_foreign[i] = mpsf;
        }
    }

    // Erase the ARTSCAT-3 foreign parameteres:
    ARTSCAT4UnusedToNaN();
}


// This function reads line data in the Hitran 1986-2001 format. For the Hitran
// 2004 data format use ReadFromHitran2004Stream.
//
// 2005/03/29: A check for the right data record length (100 characters) has
//             been added by Hermann Berg (h.berg@utoronto.ca)
//
bool LineRecord::ReadFromHitranStream(istream& is, const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Global species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotopologue lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotopologue.
  const Index missing = species_data.nelem() + 100;

  // We need a species index sorted by HITRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<HITRAN tag>]. 
  //
  // Allow for up to 100 species in HITRAN in the future.
  static Array< Index >        hspec(100);

  // This is  an array of arrays for each hitran tag. It contains the
  // ARTS indices of the HITRAN isotopologues. 
  static Array< ArrayOfIndex > hiso(100);

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  // Remember, about which missing species we have already issued a
  // warning: 
  static ArrayOfIndex warned_missing;

  if ( !hinit )
    {
      // Initialize hspec.
      // The value of missing means that we don't have this species.
      hspec = missing;  // Matpack can set all elements like this.

      for ( Index i=0; i<species_data.nelem(); ++i )
        {
          const SpeciesRecord& sr = species_data[i];

          // We have to be careful and check for the case that all
          // HITRAN isotopologue tags are -1 (this species is missing in HITRAN).

          if ( 0 < sr.Isotopologue()[0].HitranTag() )
            {
              // The HITRAN tags are stored as species plus isotopologue tags
              // (MO and ISO)
              // in the Isotopologue() part of the species record.
              // We can extract the MO part from any of the isotopologue tags,
              // so we use the first one. We do this by taking an integer
              // division by 10.
          
              Index mo = sr.Isotopologue()[0].HitranTag() / 10;
              //          cout << "mo = " << mo << endl;
              hspec[mo] = i; 
          
              // Get a nicer to handle array of HITRAN iso tags:
              Index n_iso = sr.Isotopologue().nelem();
              ArrayOfIndex iso_tags;
              iso_tags.resize(n_iso);
              for ( Index j=0; j<n_iso; ++j )
                {
                  iso_tags[j] = sr.Isotopologue()[j].HitranTag();
                }

              // Reserve elements for the isotopologue tags. How much do we
              // need? This depends on the largest HITRAN tag that we know
              // about!
              // Also initialize the tags to missing.
              //          cout << "iso_tags = " << iso_tags << endl;
              //          cout << "static_cast<Index>(max(iso_tags))%10 + 1 = "
              //               << static_cast<Index>(max(iso_tags))%10 + 1 << endl;
              hiso[mo].resize( max(iso_tags)%10 + 1 );
              hiso[mo] = missing; // Matpack can set all elements like this.


              // Set the isotopologue tags:
              for ( Index j=0; j<n_iso; ++j )
                {
                  if ( 0 < iso_tags[j] )                                  // ignore -1 elements
                    {
                      // To get the iso tags from HitranTag() we also have to take
                      // modulo 10 to get rid of mo.
                      hiso[mo][iso_tags[j] % 10] = j;
                    }
                }
            }
        }


      // Print the generated data structures (for debugging):
      out3 << "  HITRAN index table:\n";
      for ( Index i=0; i<hspec.nelem(); ++i )
        {
          if ( missing != hspec[i] )
            {
              // The explicit conversion of Name to a c-String is
              // necessary, because setw does not work correctly for
              // stl Strings.
              out3 << "  mo = " << i << "   Species = "
                   << setw(10) << setiosflags(ios::left)
                   << species_data[hspec[i]].Name().c_str()
                   << "iso = ";
              for ( Index j=1; j<hiso[i].nelem(); ++j )
                {
                  if ( missing==hiso[i][j] )
                    out3 << " " << "m";
                  else
                    out3 << " " << species_data[hspec[i]].Isotopologue()[hiso[i][j]].Name();
                }
              out3 << "\n";
            }
        }

      hinit = true;
    }


  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

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

      // If the catalogue is in dos encoding, throw away the
      // additional carriage return
      if (line[line.nelem () - 1] == 13)
        {
          line.erase (line.nelem () - 1, 1);
        }

      // Because of the fixed FORTRAN format, we need to break up the line
      // explicitly in apropriate pieces. Not elegant, but works!

      // Extract molecule number:
      mo = 0;
      // Initialization of mo is important, because mo stays the same
      // if line is empty.
      extract(mo,line,2);
      //      cout << "mo = " << mo << endl;
  
      // If mo == 0 this is just a comment line:
      if ( 0 != mo )
        {
          // See if we know this species. Exit with an error if the species is unknown. 
          if ( missing != hspec[mo] )
            {
              comment = false;

              // Check if data record has the right number of characters for the
              // in Hitran 1986-2001 format
              Index nChar = line.nelem() + 2; // number of characters in data record;
              if ( nChar != 100 )
                {
                  ostringstream os;
                  os << "Invalid HITRAN 1986-2001 line data record with " << nChar <<
                        " characters (expected: 100)." << endl << line << " n: " << line.nelem ();
                  throw runtime_error(os.str());
                }

            }
          else
            {
              // See if this is already in warned_missing, use
              // std::count for that:
              if ( 0 == std::count(warned_missing.begin(),
                                   warned_missing.end(),
                                   mo) )
                {
                  CREATE_OUT0;
                  out0 << "Error: HITRAN mo = " << mo << " is not "
                       << "known to ARTS.\n";
                  warned_missing.push_back(mo);
                }
            }
        }
    }

  // Ok, we seem to have a valid species here.

  // Set mspecies from my cool index table:
  mspecies = hspec[mo];

  // Extract isotopologue:
  Index iso;                              
  extract(iso,line,1);
  //  cout << "iso = " << iso << endl;


  // Set misotopologue from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases. 
  misotopologue = missing;
  if ( iso < hiso[mo].nelem() )
    if ( missing != hiso[mo][iso] )
      misotopologue = hiso[mo][iso];

  // Issue error message if misotopologue is still missing:
  if (missing == misotopologue)
    {
      ostringstream os;
      os << "Species: " << species_data[mspecies].Name()
         << ", isotopologue iso = " << iso
         << " is unknown.";
      throw runtime_error(os.str());
    }

  
  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 100.;

    // Extract HITRAN postion:
    extract(v,line,12);

    // ARTS position in Hz:
    mf = v * w2Hz;
//    cout << "mf = " << mf << endl;
  }

  // Intensity.
  {
    extern const Numeric SPEED_OF_LIGHT; // in [m/s]

    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function). 
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c).
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4).
    // 3. Take out the isotopologue ratio.

    const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s,line,10);
    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    mi0 = s * hi2arts;
    // Take out isotopologue ratio:
    mi0 /= species_data[mspecies].Isotopologue()[misotopologue].Abundance();  
  }  
  
  // Skip transition probability:
  {
    Numeric r;
    extract(r,line,10);
  }
  

  // Air broadening parameters.
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // External constant from constants.cc: Converts atm to
    // Pa. Multiply value in atm by this number to get value in Pa. 
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN AGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    magam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    msgam = gam * hi2arts;

    // If zero, set to agam:
    if (0==msgam)
      msgam = magam;

    //    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }


  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(melow,line,10);

    // Convert to Joule:
    melow = wavenumber_to_joule(melow);
  }

  
  // Temperature coefficient of broadening parameters.
  {
    // This is dimensionless, we can also extract directly.
    extract(mnair,line,4);

    // Set self broadening temperature coefficient to the same value:
    mnself = mnair;
//    cout << "mnair = " << mnair << endl;
  }


  // Pressure shift.
  {
    // HITRAN value in cm^-1 / atm. So the conversion goes exactly as
    // for the broadening parameters.
    Numeric d;
    // External constant from constants.cc: Converts atm to
    // Pa. Multiply value in atm by this number to get value in Pa. 
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN value:
    extract(d,line,8);

    // ARTS value in Hz/Pa
    mpsf = d * hi2arts;
  }
  // Set the accuracies using the definition of HITRAN 
  // indices. If some are missing, they are set to -1.

  //Skip upper state global quanta index
  {
    Index eu;
    extract(eu,line,3);
  }

 //Skip lower state global quanta index
  {
    Index el;
    extract(el,line,3);
  }

  //Skip upper state local quanta 
  {
    Index eul;
    extract(eul,line,9);
  }

  //Skip lower state local quanta 
  {
    Index ell;
    extract(ell,line,9);
  }

  // Accuracy index for frequency reference
  {
  Index df;
  // Extract HITRAN value:
  extract(df,line,1);
  // Convert it to ARTS units (Hz)
  convHitranIERF(mdf,df);
  }

  // Accuracy index for intensity reference
  {
  Index di0;
  // Extract HITRAN value:
    extract(di0,line,1);
    convHitranIERSH(mdi0,di0);
  }

  // Accuracy index for halfwidth reference
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam,line,1);
    //Convert to ARTS units (%)
    convHitranIERSH(mdagam,dgam);
    // convHitranIERSH(mdsgam,dgam);
    // convHitranIERSH(mdnair,dgam);
    // convHitranIERSH(mdnself,dgam);
  }

  // Accuracy for pressure shift
  // This is missing in HITRAN catalogue and it is set to -1.
    mdpsf =-1;

  // These were all the parameters that we can extract from
  // HITRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for HITRAN)
  mti0 = 296.0;

  // Reference temperature for AGAM and SGAM in K.
  // (This is also fix for HITRAN)
  mtgam = 296.0;

  // That's it!
  return false;
}


// This function reads line data in the Hitran 2004 format. For the Hitran
// 1986-2004 data format use ReadFromHitranStream().
//
// 2005/03/29: This function was added based on ReadFromHitranStream(). There
//             is a check for the right data record length (160 characters).
//             Hermann Berg (h.berg@utoronto.ca)
//
bool LineRecord::ReadFromHitran2004Stream(istream& is, const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Global species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotopologue lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotopologue.
  const Index missing = species_data.nelem() + 100;

  // We need a species index sorted by HITRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<HITRAN tag>].
  //
  // Allow for up to 100 species in HITRAN in the future.
  static Array< Index >        hspec(100);

  // This is  an array of arrays for each hitran tag. It contains the
  // ARTS indices of the HITRAN isotopologues.
  static Array< ArrayOfIndex > hiso(100);

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  // Remember, about which missing species we have already issued a
  // warning:
  static ArrayOfIndex warned_missing;

  if ( !hinit )
    {
      // Initialize hspec.
      // The value of missing means that we don't have this species.
      hspec = missing;  // Matpack can set all elements like this.

      for ( Index i=0; i<species_data.nelem(); ++i )
        {
          const SpeciesRecord& sr = species_data[i];

          // We have to be careful and check for the case that all
          // HITRAN isotopologue tags are -1 (this species is missing in HITRAN).

          if ( 0 < sr.Isotopologue()[0].HitranTag() )
            {
              // The HITRAN tags are stored as species plus isotopologue tags
              // (MO and ISO)
              // in the Isotopologue() part of the species record.
              // We can extract the MO part from any of the isotopologue tags,
              // so we use the first one. We do this by taking an integer
              // division by 10.

              Index mo = sr.Isotopologue()[0].HitranTag() / 10;
              //          cout << "mo = " << mo << endl;
              hspec[mo] = i;

              // Get a nicer to handle array of HITRAN iso tags:
              Index n_iso = sr.Isotopologue().nelem();
              ArrayOfIndex iso_tags;
              iso_tags.resize(n_iso);
              for ( Index j=0; j<n_iso; ++j )
                {
                  iso_tags[j] = sr.Isotopologue()[j].HitranTag();
                }

              // Reserve elements for the isotopologue tags. How much do we
              // need? This depends on the largest HITRAN tag that we know
              // about!
              // Also initialize the tags to missing.
              //          cout << "iso_tags = " << iso_tags << endl;
              //          cout << "static_cast<Index>(max(iso_tags))%10 + 1 = "
              //               << static_cast<Index>(max(iso_tags))%10 + 1 << endl;
              hiso[mo].resize( max(iso_tags)%10 + 1 );
              hiso[mo] = missing; // Matpack can set all elements like this.


              // Set the isotopologue tags:
              for ( Index j=0; j<n_iso; ++j )
                {
                  if ( 0 < iso_tags[j] )                                  // ignore -1 elements
                    {
                      // To get the iso tags from HitranTag() we also have to take
                      // modulo 10 to get rid of mo.
                      hiso[mo][iso_tags[j] % 10] = j;
                    }
                }
            }
        }


      // Print the generated data structures (for debugging):
      out3 << "  HITRAN index table:\n";
      for ( Index i=0; i<hspec.nelem(); ++i )
        {
          if ( missing != hspec[i] )
            {
              // The explicit conversion of Name to a c-String is
              // necessary, because setw does not work correctly for
              // stl Strings.
              out3 << "  mo = " << i << "   Species = "
                   << setw(10) << setiosflags(ios::left)
                   << species_data[hspec[i]].Name().c_str()
                   << "iso = ";
              for ( Index j=1; j<hiso[i].nelem(); ++j )
                {
                  if ( missing==hiso[i][j] )
                    out3 << " " << "m";
                  else
                    out3 << " " << species_data[hspec[i]].Isotopologue()[hiso[i][j]].Name();
                }
              out3 << "\n";
            }
        }

      hinit = true;
    }


  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

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

      // If the catalogue is in dos encoding, throw away the
      // additional carriage return
      if (line[line.nelem () - 1] == 13)
        {
          line.erase (line.nelem () - 1, 1);
        }

      // Because of the fixed FORTRAN format, we need to break up the line
      // explicitly in apropriate pieces. Not elegant, but works!

      // Extract molecule number:
      mo = 0;
      // Initialization of mo is important, because mo stays the same
      // if line is empty.
      extract(mo,line,2);
      //      cout << "mo = " << mo << endl;

      // If mo == 0 this is just a comment line:
      if ( 0 != mo )
        {
          // See if we know this species. 
          if ( missing != hspec[mo] )
            {
              comment = false;
              
              // Check if data record has the right number of characters for the
              // in Hitran 2004 format
              Index nChar = line.nelem() + 2; // number of characters in data record;
              if ( nChar != 160 )
                {
                  ostringstream os;
                  os << "Invalid HITRAN 2004 line data record with " << nChar <<
                        " characters (expected: 160).";
                  throw runtime_error(os.str());
                }
                
            }
          else
            {
              // See if this is already in warned_missing, use
              // std::count for that:
              if ( 0 == std::count(warned_missing.begin(),
                                   warned_missing.end(),
                                   mo) )
                {
                  CREATE_OUT1;
                  out1 << "Warning: HITRAN molecule number mo = " << mo << " is not "
                       << "known to ARTS.\n";
                  warned_missing.push_back(mo);
               }
            }
        }
    }

  // Ok, we seem to have a valid species here.

  // Set mspecies from my cool index table:
  mspecies = hspec[mo];

  // Extract isotopologue:
  Index iso;
  extract(iso,line,1);
  //  cout << "iso = " << iso << endl;


  // Set misotopologue from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases.
  misotopologue = missing;
  if ( iso < hiso[mo].nelem() )
    if ( missing != hiso[mo][iso] )
      misotopologue = hiso[mo][iso];

  // Issue error message if misotopologue is still missing:
  if (missing == misotopologue)
    {
      ostringstream os;
      os << "Species: " << species_data[mspecies].Name()
         << ", isotopologue iso = " << iso
         << " is unknown.";
      throw runtime_error(os.str());
    }


  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 100.;

    // Extract HITRAN postion:
    extract(v,line,12);

    // ARTS position in Hz:
    mf = v * w2Hz;
//    cout << "mf = " << mf << endl;
  }

  // Intensity.
  {
    extern const Numeric SPEED_OF_LIGHT; // in [m/s]

    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function).
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c).
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4).
    // 3. Take out the isotopologue ratio.

    const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s,line,10);
    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    mi0 = s * hi2arts;
    // Take out isotopologue ratio:
    mi0 /= species_data[mspecies].Isotopologue()[misotopologue].Abundance();
  }

  // Skip Einstein coefficient
  {
    Numeric r;
    extract(r,line,10);
  }


  // Air broadening parameters.
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // External constant from constants.cc: Converts atm to
    // Pa. Multiply value in atm by this number to get value in Pa.
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN AGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    magam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    msgam = gam * hi2arts;

    // If zero, set to agam:
    if (0==msgam)
      msgam = magam;

    //    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }


  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(melow,line,10);

    // Convert to Joule:
    melow = wavenumber_to_joule(melow);
  }


  // Temperature coefficient of broadening parameters.
  {
    // This is dimensionless, we can also extract directly.
    extract(mnair,line,4);

    // Set self broadening temperature coefficient to the same value:
    mnself = mnair;
//    cout << "mnair = " << mnair << endl;
  }


  // Pressure shift.
  {
    // HITRAN value in cm^-1 / atm. So the conversion goes exactly as
    // for the broadening parameters.
    Numeric d;
    // External constant from constants.cc: Converts atm to
    // Pa. Multiply value in atm by this number to get value in Pa.
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN value:
    extract(d,line,8);

    // ARTS value in Hz/Pa
    mpsf = d * hi2arts;
  }
  // Set the accuracies using the definition of HITRAN
  // indices. If some are missing, they are set to -1.

  // Upper state global quanta
  {
    mupper_gquanta = line.substr(0,15);
    Index eu;
    extract(eu,line,15);
  }

  // Lower state global quanta
  {
    mlower_gquanta = line.substr(0,15);
    Index el;
    extract(el,line,15);
  }

  // Upper state local quanta
  {
    mupper_lquanta = line.substr(0,15);
    Index eul;
    extract(eul,line,15);
  }

  // Lower state local quanta
  {
    mlower_lquanta = line.substr(0,15);
    Index ell;
    extract(ell,line,15);
  }

  // Assign the local quantum numbers.
  if (mlower_lquanta.nelem() == 15)
  {
      Index DN = 0, DJ = 0;
      if(species_data[mspecies].Name()=="O2")//O2 FIXME: Different mspecies versus HITRAN?
      {
        mlower_n = atoi(mlower_lquanta.substr(2,3).c_str());
        mlower_j = atoi(mlower_lquanta.substr(6,3).c_str());
        DJ =  -  mlower_lquanta.compare(5,1,"Q");
        DN =  -  mlower_lquanta.compare(1,1,"Q");
        mupper_n = mlower_n - DN;
        mupper_j = mlower_j - DJ;
      }
      else if(species_data[mspecies].Name()=="NO")//NO FIXME: Different mspecies versus HITRAN?
      {
        DJ = - mlower_lquanta.compare(3,1,"Q");
        
        int nom=-1, denom=0;
        if (mlower_lquanta.substr(4,3) == "   ")
        {
            nom   = atoi(mlower_lquanta.substr(4,3).c_str());
            denom = atoi(mlower_lquanta.substr(8,1).c_str());
        }
        else //out2 << "Hitran read error in mlower_lquanta J for species " <<  mspecies << " in line " << mf << "\n";
        
        if( denom == 5 ) mlower_j = Rational(nom*2+1,2);
        else mlower_j = nom;
        mupper_j = mlower_j - DJ;
        
        nom=-1;denom=0;
        if (mlower_gquanta.substr(4,3) == "   ")
        {
            nom   = atoi(mlower_lquanta.substr(8,1).c_str());
            denom = atoi(mlower_lquanta.substr(10,1).c_str());
        }
        else //out2 << "Hitran read error in mlower_gquanta i for species " <<  mspecies << " in line " << mf << "\n";
        
        // Note that OMEGA is nom/denom
        
        if( nom==3&&denom==2 ) mlower_n = mlower_j+Rational(1,2);
        else if( nom==1&&denom==2 ) mlower_n = mlower_j-Rational(1,2);
        else mlower_n = -1;
        
        nom=-1;denom=0;
        if (mupper_gquanta.substr(4,3) == "   ")
        {
            nom   = atoi(mupper_lquanta.substr(8,1).c_str());
            denom = atoi(mupper_lquanta.substr(10,1).c_str());
        }
        else //out2 << "Hitran read error in mlower_gquanta i for species " <<  mspecies << " in line " << mf << "\n";
        
        // Note that OMEGA is nom/denom
        
        if( nom==3&&denom==2 ) mupper_n = mupper_j+Rational(1,2);
        else if( nom==1&&denom==2 ) mupper_n = mupper_j-Rational(1,2);
        else mupper_n = -1;
        
      }
      else if(species_data[mspecies].Name()=="OH")//OH FIXME: Different mspecies versus HITRAN?
      {
        //pass for now
      }
      else if(species_data[mspecies].Name()=="ClO")//ClO FIXME: Different mspecies versus HITRAN?
      {
        //pass for now
      }
  }

  // Accuracy index for frequency
  {
    Index df;
    // Extract HITRAN value:
    extract(df,line,1);
    // Convert it to ARTS units (Hz)
    convHitranIERF(mdf,df);
  }

  // Accuracy index for intensity
  {
    Index di0;
    // Extract HITRAN value:
    extract(di0,line,1);
    //Convert to ARTS units (%)
    convHitranIERSH(mdi0,di0);
  }

  // Accuracy index for air-broadened halfwidth
  {
    Index dagam;
    // Extract HITRAN value:
    extract(dagam,line,1);
    //Convert to ARTS units (%)
    convHitranIERSH(mdagam,dagam);
  }

  // Accuracy index for self-broadened half-width
  {
    Index dsgam;
    // Extract HITRAN value:
    extract(dsgam,line,1);
    //Convert to ARTS units (%)
    convHitranIERSH(mdsgam,dsgam);
  }
  
  // Accuracy index for temperature-dependence exponent for agam
  {
    Index dnair;
    // Extract HITRAN value:
    extract(dnair,line,1);
    //Convert to ARTS units (%)
    convHitranIERSH(mdnair,dnair);
  }

  // Accuracy index for temperature-dependence exponent for sgam
  // This is missing in HITRAN catalogue and is set to -1.
    mdnself =-1;

  // Accuracy index for pressure shift
  {
    Index dpsf;
    // Extract HITRAN value (given in cm-1):
    extract(dpsf,line,1);
    // Convert it to ARTS units (Hz)
    convHitranIERF(mdpsf,dpsf);
    // ARTS wants this error in %
    mdpsf = mdpsf / mf;
  }

  // These were all the parameters that we can extract from
  // HITRAN 2004. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for HITRAN 2004)
  mti0 = 296.0;

  // Reference temperature for AGAM and SGAM in K.
  // (This is also fix for HITRAN 2004)
  mtgam = 296.0;

  // That's it!
  return false;
}


bool LineRecord::ReadFromMytran2Stream(istream& is, const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Global species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotopologue lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotopologue.
  const Index missing = species_data.nelem() + 100;

  // We need a species index sorted by MYTRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<MYTRAN tag>]. The value of
  // missing means that we don't have this species.
  //
  // Allow for up to 100 species in MYTRAN in the future.
  static Array< Index >        hspec(100,missing);      

  // This is  an array of arrays for each mytran tag. It contains the
  // ARTS indices of the MYTRAN isotopologues. 
  static Array< ArrayOfIndex > hiso(100);

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  // Remember, about which missing species we have already issued a
  // warning: 
  static ArrayOfIndex warned_missing;

  if ( !hinit )
    {
      for ( Index i=0; i<species_data.nelem(); ++i )
        {
          const SpeciesRecord& sr = species_data[i];

          // We have to be careful and check for the case that all
          // MYTRAN isotopologue tags are -1 (this species is missing in MYTRAN).

          if ( 0 < sr.Isotopologue()[0].MytranTag() )
            {
              // The MYTRAN tags are stored as species plus isotopologue tags
              // (MO and ISO)
              // in the Isotopologue() part of the species record.
              // We can extract the MO part from any of the isotopologue tags,
              // so we use the first one. We do this by taking an integer
              // division by 10.
          
              Index mo = sr.Isotopologue()[0].MytranTag() / 10;
              //          cout << "mo = " << mo << endl;
              hspec[mo] = i; 
          
              // Get a nicer to handle array of MYTRAN iso tags:
              Index n_iso = sr.Isotopologue().nelem();
              ArrayOfIndex iso_tags;
              iso_tags.resize(n_iso);
              for ( Index j=0; j<n_iso; ++j )
                {
                  iso_tags[j] = sr.Isotopologue()[j].MytranTag();
                }

              // Reserve elements for the isotopologue tags. How much do we
              // need? This depends on the largest MYTRAN tag that we know
              // about!
              // Also initialize the tags to missing.
              //          cout << "iso_tags = " << iso_tags << endl;
              //          cout << "static_cast<Index>(max(iso_tags))%10 + 1 = "
              //               << static_cast<Index>(max(iso_tags))%10 + 1 << endl;
              hiso[mo].resize( max(iso_tags)%10 + 1 );
              hiso[mo] = missing; // Matpack can set all elements like this.

              // Set the isotopologue tags:
              for ( Index j=0; j<n_iso; ++j )
                {
                  if ( 0 < iso_tags[j] )                                  // ignore -1 elements
                    {
                      // To get the iso tags from MytranTag() we also have to take
                      // modulo 10 to get rid of mo.
                      //                  cout << "iso_tags[j] % 10 = " << iso_tags[j] % 10 << endl;
                      hiso[mo][iso_tags[j] % 10] = j;
                    }
                }
            }
        }
      
//      cout << "hiso = " << hiso << endl << "***********" << endl;


      // Print the generated data structures (for debugging):
      out3 << "  MYTRAN index table:\n";
      for ( Index i=0; i<hspec.nelem(); ++i )
        {
          if ( missing != hspec[i] )
            {
              // The explicit conversion of Name to a c-String is
              // necessary, because setw does not work correctly for
              // stl Strings.
              out3 << "  mo = " << i << "   Species = "
                   << setw(10) << setiosflags(ios::left)
                   << species_data[hspec[i]].Name().c_str()
                   << "iso = ";
              for ( Index j=1; j<hiso[i].nelem(); ++j )
                {
                  if ( missing==hiso[i][j] )
                    out3 << " " << "m";
                  else
                    out3 << " " << species_data[hspec[i]].Isotopologue()[hiso[i][j]].Name();
                }
              out3 << "\n";
            }
        }

      hinit = true;
    }


  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

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

      // Because of the fixed FORTRAN format, we need to break up the line
      // explicitly in apropriate pieces. Not elegant, but works!

      // Extract molecule number:
      mo = 0;
      // Initialization of mo is important, because mo stays the same
      // if line is empty.
      extract(mo,line,2);
      //           cout << "mo = " << mo << endl;
  
      // If mo == 0 this is just a comment line:
      if ( 0 != mo )
        {
          // See if we know this species. We will give an error if a
          // species is not known. 
          if ( missing != hspec[mo] )       comment = false ;
          else
            {
              // See if this is already in warned_missing, use
              // std::count for that:
              if ( 0 == std::count(warned_missing.begin(),
                                   warned_missing.end(),
                                   mo) )
                {
                  CREATE_OUT0;
                  out0 << "Error: MYTRAN mo = " << mo << " is not "
                       << "known to ARTS.\n";
                  warned_missing.push_back(mo);
                }
            }
        }
    }

  // Ok, we seem to have a valid species here.

  // Set mspecies from my cool index table:
  mspecies = hspec[mo];

  // Extract isotopologue:
  Index iso;                              
  extract(iso,line,1);
  //  cout << "iso = " << iso << endl;


  // Set misotopologue from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases. 
  misotopologue = missing;
  if ( iso < hiso[mo].nelem() )
    if ( missing != hiso[mo][iso] )
      misotopologue = hiso[mo][iso];

  // Issue error message if misotopologue is still missing:
  if (missing == misotopologue)
    {
      ostringstream os;
      os << "Species: " << species_data[mspecies].Name()
         << ", isotopologue iso = " << iso
         << " is unknown.";
      throw runtime_error(os.str());
    }

  
  // Position.
  {
    // MYTRAN position in MHz:
    Numeric v;

    // Extract MYTRAN postion:
    extract(v,line,13);

    // ARTS position in Hz:
    mf = v * 1E6;
    //    cout << "mf = " << mf << endl;
  }

  // Accuracy for line position
  {
    // Extract MYTRAN postion accuracy:
    Numeric df;
    extract(df,line,8);
    //  ARTS accuracy of line position  in Hz:
    mdf = df * 1E6;
  } 
  
  // Intensity.
  {
    extern const Numeric SPEED_OF_LIGHT; // in [m/s]

    // MYTRAN2 intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // (just like HITRAN, only isotopologue ratio is not included)
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function). 
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c)
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4)

    const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s,line,10);

    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    mi0 = s * hi2arts;
  }

  
  // Air broadening parameters.
  {
    // MYTRAN parameter is in MHz/Torr at reference temperature
    // All parameters are HWHM
    Numeric gam;
    // External constant from constants.cc: Converts torr to
    // Pa. Multiply value in torr by this number to get value in Pa. 
    extern const Numeric TORR2PA;

    // Extract HITRAN AGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    magam = gam * 1E6 / TORR2PA;

    // Extract MYTRAN SGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    msgam = gam * 1E6 / TORR2PA;

//    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }


  // Lower state energy.
  {
    // MYTRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(melow,line,10);

    // Convert to Joule:
    melow = wavenumber_to_joule(melow);
  }

  
  // Temperature coefficient of broadening parameters.
  {
    // This is dimensionless, we can also extract directly.
    extract(mnair,line,4);

    // Extract the self broadening parameter:
    extract(mnself,line,4);
//    cout << "mnair = " << mnair << endl;
  }

  
  // Reference temperature for broadening parameter in K:
  {
    // correct units, extract directly
    extract(mtgam,line,7);
  }


  // Pressure shift.
  {
    // MYTRAN value in MHz/Torr
    Numeric d;
    // External constant from constants.cc: Converts torr to
    // Pa. Multiply value in torr by this number to get value in Pa. 
    extern const Numeric TORR2PA;

    // Extract MYTRAN value:
    extract(d,line,9);

    // ARTS value in Hz/Pa
    mpsf = d * 1E6 / TORR2PA;
  }
  // Set the accuracies using the definition of MYTRAN accuracy
  // indices. If some are missing, they are set to -1.

  //Skip upper state global quanta index
  {
    Index eu;
    extract(eu,line,3);
  }

 //Skip lower state global quanta index
  {
    Index el;
    extract(el,line,3);
  }

  //Skip upper state local quanta 
  {
    Index eul;
    extract(eul,line,9);
  }

  //Skip lower state local quanta 
  {
    Index ell;
    extract(ell,line,9);
  }
 // Accuracy index for intensity
  {
  Index di0;
  // Extract MYTRAN value:
  extract(di0,line,1);
  //convert to ARTS units (%)
  convMytranIER(mdi0,di0);
  }

  // Accuracy index for AGAM
  {
  Index dgam;
  // Extract MYTRAN value:
  extract(dgam,line,1);
  //convert to ARTS units (%)
  convMytranIER(mdagam,dgam);
  }

  // Accuracy index for NAIR 
  {
  Index dnair;
  // Extract MYTRAN value:
  extract(dnair,line,1); 
  //convert to ARTS units (%);
  convMytranIER(mdnair,dnair);
  }


  // These were all the parameters that we can extract from
  // MYTRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for MYTRAN2)
  mti0 = 296.0;

  // It is important that you intialize here all the new parameters that
  // you added to the line record. (This applies to all the reading
  // functions, also for ARTS, JPL, and HITRAN format.) Parameters
  // should be either set from the catalogue, or set to -1.)

  // That's it!
  return false;
}


bool LineRecord::ReadFromJplStream(istream& is, const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Global species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // We need a species index sorted by JPL tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is JplMap[<JPL tag>]. We need Index in this map,
  // because the tag array within the species data is an Index array.
  static map<Index, SpecIsoMap> JplMap;

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  if ( !hinit )
    {

      out3 << "  JPL index table:\n";

      for ( Index i=0; i<species_data.nelem(); ++i )
        {
          const SpeciesRecord& sr = species_data[i];


          for ( Index j=0; j<sr.Isotopologue().nelem(); ++j)
            {
              
              for (Index k=0; k<sr.Isotopologue()[j].JplTags().nelem(); ++k)
                {

                  SpecIsoMap indicies(i,j);

                  JplMap[sr.Isotopologue()[j].JplTags()[k]] = indicies;

                  // Print the generated data structures (for debugging):
                  // The explicit conversion of Name to a c-String is
                  // necessary, because setw does not work correctly for
                  // stl Strings.
                  const Index& i1 = JplMap[sr.Isotopologue()[j].JplTags()[k]].Speciesindex();
                  const Index& i2 = JplMap[sr.Isotopologue()[j].JplTags()[k]].Isotopologueindex();
                                         
                  out3 << "  JPL TAG = " << sr.Isotopologue()[j].JplTags()[k] << "   Species = "
                       << setw(10) << setiosflags(ios::left)
                       << species_data[i1].Name().c_str()
                       << "iso = " 
                       << species_data[i1].Isotopologue()[i2].Name().c_str()
                       << "\n";
                }

            }
        }
      hinit = true;
    }


  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
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

      // Because of the fixed FORTRAN format, we need to break up the line
      // explicitly in apropriate pieces. Not elegant, but works!

      // Extract center frequency:
      // Initialization of v is important, because v stays the same
      // if line is empty.
      // JPL position in MHz:
      Numeric v = 0.0;
      
      // Extract JPL position:
      extract(v,line,13);
        
      // check for empty line
      if (v != 0.0)
        {
          // ARTS position in Hz:
          mf = v * 1E6;
          
          comment = false;
        }
    }

  // Accuracy for line position 
  {
    Numeric df;
    extract(df,line,8);
    //convert to ARTS units (Hz)
    mdf = df * 1E6; 
  } 
  
  // Intensity.
  {
    // JPL has log (10) of intensity in nm2 MHz at 300 Kelvin.
    //
    // We need to do the following:
    // 1. take 10^intensity 
    // 2. convert to cm-1/(molecule * cm-2): devide by c * 1e10
    // 3. Convert frequency from wavenumber to Hz (factor 1e2 * c)
    // 4. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4)

    Numeric s;

    // Extract JPL intensity:
    extract(s,line,8);

    // remove log
    s = pow((Numeric)10.,s);

    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    mi0 = s / 1E12;
  }

  // Degrees of freedom
  {
    Index dr;

    // Extract degrees of freedom
    extract(dr,line,2);
  }

  // Lower state energy.
  {
    // JPL parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(melow,line,10);

    // Convert to Joule:
    melow = wavenumber_to_joule(melow);
  }

  // Upper state degeneracy
  {
    Index gup;

    // Extract upper state degeneracy
    extract(gup,line,3);
  }

  // Tag number
  Index tag;
  {
    // Extract Tag number
    extract(tag,line,7);

    // make sure tag is not negative (damned jpl cat):
    tag = tag > 0 ? tag : -tag;
  }

  // ok, now for the cool index map:

  // is this tag valid?
  const map<Index, SpecIsoMap>::const_iterator i = JplMap.find(tag);
  if ( i == JplMap.end() )
    {
      ostringstream os;
      os << "JPL Tag: " << tag << " is unknown.";
      throw runtime_error(os.str());
    }

  SpecIsoMap id = i->second;


  // Set mspecies:
  mspecies = id.Speciesindex();

  // Set misotopologue:
  misotopologue = id.Isotopologueindex();

  // Air broadening parameters: unknown to jpl, use old iup forward
  // model default values, which is mostly set to 0.0025 GHz/hPa, even
  // though for some lines the pressure broadening is given explicitly
  // in the program code. The explicitly given values are ignored and
  // only the default value is set. Self broadening was in general not
  // considered in the old forward model.
  {
    // ARTS parameter in Hz/Pa:
    magam = 2.5E4;

    // ARTS parameter in Hz/Pa:
    msgam = magam;
  }


  // Temperature coefficient of broadening parameters. Was set to 0.75
  // in old forward model, even though for some lines the parameter is
  // given explicitly in the program code. The explicitly given values
  // are ignored and only the default value is set. Self broadening
  // not considered.
  {
    mnair = 0.75;
    mnself = 0.0;
  }

  
  // Reference temperature for broadening parameter in K, was
  // generally set to 300 K in old forward model, with the exceptions
  // as already mentioned above:
  {
    mtgam = 300.0;
  }


  // Pressure shift: not given in JPL, set to 0
  {
    mpsf = 0.0;
  }


  // These were all the parameters that we can extract from
  // JPL. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for JPL)
  mti0 = 300.0;

  // That's it!
  return false;
}


bool LineRecord::ReadFromArtscat3Stream(istream& is, const Verbosity& verbosity)
{
  CREATE_OUT3;
 
  // Global species lookup data:
  extern const Array<SpeciesRecord> species_data;
  
  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>]. 
  static map<String, SpecIsoMap> ArtsMap;
  
  // Remember if this stuff has already been initialized:
  static bool hinit = false;
  
  mversion = 3;
  
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
  
  String artsid;
  icecream >> artsid;
  
  if (artsid.length() != 0)
  {
    
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
    
    
    // Extract center frequency:
    icecream >> mf;
    
    
    // Extract pressure shift:
    icecream >> mpsf;
    
    // Extract intensity:
    icecream >> mi0;
    
    
    // Extract reference temperature for Intensity in K:
    icecream >> mti0;
    
    
    // Extract lower state energy:
    icecream >> melow;
    
    
    // Extract air broadening parameters:
    icecream >> magam;
    icecream >> msgam;
    
    // Extract temperature coefficient of broadening parameters:
    icecream >> mnair;
    icecream >> mnself;
    
    
    // Extract reference temperature for broadening parameter in K:
    icecream >> mtgam;
    
    // Extract the aux parameters:
    Index naux;
    icecream >> naux;
    
    // resize the aux array and read it
    maux.resize(naux);
    
    for (Index j = 0; j<naux; j++)
    {
      icecream >> maux[j];
      //cout << "maux" << j << " = " << maux[j] << "\n";
    }
    
    // Extract accuracies:
    try
    {
      icecream >> mdf;
      icecream >> mdi0;
      icecream >> mdagam;
      icecream >> mdsgam;
      icecream >> mdnair;
      icecream >> mdnself;
      icecream >> mdpsf;
    }
    catch (runtime_error)
    {
      // Nothing to do here, the accuracies are optional, so we
      // just set them to -1 and continue reading the next line of
      // the catalogue
      mdf      = -1;
      mdi0     = -1;
      mdagam   = -1;
      mdsgam   = -1;
      mdnair   = -1;
      mdnself  = -1;
      mdpsf    = -1;            
    }
  }
  
  // That's it!
  return false;
}

bool LineRecord::ReadFromArtscat4Stream(istream& is, const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Global species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>]. 
  static map<String, SpecIsoMap> ArtsMap;

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  mversion = 4;
  
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

  String artsid;
  icecream >> artsid;

  if (artsid.length() != 0)
    {

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


      // Extract center frequency:
      icecream >> mf;


      // Extract intensity:
      icecream >> mi0;

      // Extract reference temperature for Intensity in K:
      icecream >> mti0;
      
      // Extract lower state energy:
      icecream >> melow;

      // Extract Einstein A-coefficient:
      icecream >> ma;
      
      // Extract upper state stat. weight:
      icecream >> mgupper;
      
      // Extract lower state stat. weight:
      icecream >> mglower;

      // Extract broadening parameters:
      icecream >> msgam;

      mgamma_foreign.resize(6);
      for (Index s=0; s<6; ++s)
          icecream >> mgamma_foreign[s];
//      icecream >> mgamma_n2;
//      icecream >> mgamma_o2;
//      icecream >> mgamma_h2o;
//      icecream >> mgamma_co2;
//      icecream >> mgamma_h2;
//      icecream >> mgamma_he;
      
      // Extract GAM temp. exponents:
      icecream >> mnself;
    
      mn_foreign.resize(6);
      for (Index s=0; s<6; ++s)
          icecream >> mn_foreign[s];
//      icecream >> mn_n2;
//      icecream >> mn_o2;
//      icecream >> mn_h2o;
//      icecream >> mn_co2;
//      icecream >> mn_h2;
//      icecream >> mn_he;

      // Extract F pressure shifts:
      mdelta_foreign.resize(6);
      for (Index s=0; s<6; ++s)
          icecream >> mdelta_foreign[s];
//      icecream >> mdelta_n2;
//      icecream >> mdelta_o2;
//      icecream >> mdelta_h2o;
//      icecream >> mdelta_co2;
//      icecream >> mdelta_h2;
//      icecream >> mdelta_he;

      // FIXME: OLE: Added this if to catch crash for species like CO, PH3
      // where the line in the catalog is too short. Better would be to
      // only read the n and j for Zeeman species, but we don't have that
      // information here
      if (line.nelem() >= 288+12+1+12)
      {
          String qstr1 = line.substr(288,      12);
          String qstr2 = line.substr(288+12+1, 12);
          ArrayOfIndex q(6);
          for (Index qi=0; qi<3; qi++)
              extract(q[qi], qstr1, 4);
          for (Index qi=3; qi<6; qi++)
              extract(q[qi], qstr2, 4);

          if (q[0] > 0) mupper_n = q[0];
          if (q[1] > 0) mlower_n = q[1];
          if (q[3] > 0) mupper_j = q[3];
          if (q[4] > 0) mlower_j = q[4];
      }

      LineRecord::ARTSCAT4UnusedToNaN();
    }

  // That's it!
  return false;
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
    \param  l_l Spectral line data record (a single line).
    \param  verbosity Verbosity flag.

    \author Stefan Buehler
    \date   2012-09-05
 */
void calc_gamma_and_deltaf_artscat4(Numeric& gamma,
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
                throw runtime_error(os.str());
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

//  {
//    CREATE_OUT3;
//    ostringstream os;
//    os << "  Broad_spec_vmr_sum: " << broad_spec_vmr_sum << "\n";
//    out3 << os.str();
//  }
    
    // Check that sum of self and all foreign VMRs is not too far from 1:
    if ( abs(vmrs[this_species]+broad_spec_vmr_sum-1) > 0.1 )
      {
        ostringstream os;
//        os << "Error: The total VMR of all your defined broadening\n"
        os << "Warning: The total VMR of all your defined broadening\n"
             << "species (including \"self\") is " << broad_spec_vmr_sum
             << ", more than 10% " << "different from 1.\n";
//        throw runtime_error(os.str());
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
void xsec_species( MatrixView               xsec_attenuation,
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
  extern const Array<LineshapeRecord> lineshape_data;
  extern const Array<LineshapeNormRecord> lineshape_norm_data;

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
          throw runtime_error(os.str());
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
      throw runtime_error(os.str());
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
      throw runtime_error(os.str());
    }

  // all_vmrs should have dimensions [nspecies, np]:
  
  if ( all_vmrs.ncols() != np )
    {
      ostringstream os;
      os << "Number of columns of all_vmrs must match abs_p.\n"
         << "all_vmrs.ncols() = " << all_vmrs.ncols() << '\n'
         << "abs_p.nelem() = " << np;
      throw runtime_error(os.str());
    }

  const Index nspecies = abs_species.nelem();
  
  if ( all_vmrs.nrows() != nspecies)
  {
    ostringstream os;
    os << "Number of rows of all_vmrs must match abs_species.\n"
    << "all_vmrs.nrows() = " << all_vmrs.nrows() << '\n'
    << "abs_species.nelem() = " << nspecies;
    throw runtime_error(os.str());
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
      throw runtime_error(os.str());
    }
   if ( xsec_phase.nrows() != nf || xsec_phase.ncols() != np )
    {
      ostringstream os;
      os << "Variable xsec must have dimensions [f_grid.nelem(),abs_p.nelem()].\n"
         << "[xsec_phase.nrows(),xsec_phase.ncols()] = [" << xsec_phase.nrows()
         << ", " << xsec_phase.ncols() << "]\n"
         << "f_grid.nelem() = " << nf << '\n'
         << "abs_p.nelem() = " << np;
      throw runtime_error(os.str());
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
#pragma omp parallel for                               \
  if(!arts_omp_in_parallel()                           \
     && (np >= arts_omp_get_max_threads() || np > nl)) \
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
#pragma omp parallel for                                            \
  if(!arts_omp_in_parallel())                                       \
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

              // The baseline to substract for cutoff frequency
              Numeric base_attenuation=0.0;
              Numeric base_phase=0.0;

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
                throw runtime_error(os.str());
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
              Numeric deltaf;    // Pressure shift.
              if (l_l.Version() == 4)
                {
                  calc_gamma_and_deltaf_artscat4(gamma,
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
                  const Numeric Tgam = l_l.Tgam();
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
                      // The index nfls-1 should be exactly the index pointing
                      // to the value at the cutoff frequency.
                      base_attenuation = ls_attenuation[nfls-1];
                      base_phase = ls_phase[nfls-1];

                      // Subtract baseline from xsec. 
                      // this_xsec -= base;
                      this_ls_attenuation -= base_attenuation;
                      this_ls_phase       -= base_phase;
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
                          * isotopologue_ratios.getParam(l_l.Species(), l_l.Isotopologue(), 0);

                    // We have to do:
                    // xsec(j,i) += factors * ls[j] * fac[j];
                    //
                    // We use ls as a dummy to compute the product, then add it
                    // to this_xsec.

                    this_ls_attenuation *= this_fac;
                    this_ls_attenuation *= factors;
                    this_ls_phase *= this_fac;
                    this_ls_phase *= factors;
                    
                    this_xsec_attenuation += this_ls_attenuation;
                    this_xsec_phase += this_ls_phase;

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
      for (Index j=0; j<xsec_accum_phase.nrows(); ++j)
        {
            xsec_i_phase += xsec_accum_phase(j, Range(joker));
        }
    } // end of parallel pressure loop

  if (failed) throw runtime_error("Run-time error in function: xsec_species\n" + fail_msg);

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

