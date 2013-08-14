/* Copyright (C) 2000-2013
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

/** \file

  LineRecord implementation.

  \author Stefan Buehler and Axel von Engeln
*/

#include <cfloat>
#include "linerecord.h"
#include "absorption.h"

#include "global_data.h"


String LineRecord::VersionString() const
{
    ostringstream os;
    os << "ARTSCAT-" << mversion;
    return os.str();
}


String LineRecord::Name() const {
  // The species lookup data:
  using global_data::species_data;
  const SpeciesRecord& sr = species_data[mspecies];
  return sr.Name() + "-" + sr.Isotopologue()[misotopologue].Name();
}


const SpeciesRecord& LineRecord::SpeciesData() const {
  // The species lookup data:
  using global_data::species_data;
  return species_data[mspecies];
}


const IsotopologueRecord& LineRecord::IsotopologueData() const {
  // The species lookup data:
  using global_data::species_data;
  return species_data[mspecies].Isotopologue()[misotopologue];
}


Index LineRecord::BroadSpecSpecIndex(const Index i)  {
    // No need for asserts of i here, since the default clause in
    // BroadSpecName catches everything.
    return species_index_from_species_name(BroadSpecName(i));
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


bool LineRecord::ReadFromHitran2001Stream(istream& is, const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Global species lookup data:
  using global_data::species_data;

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


bool LineRecord::ReadFromHitran2004Stream(istream& is, const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Global species lookup data:
  using global_data::species_data;

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

          if ( sr.Isotopologue().nelem() && 0 < sr.Isotopologue()[0].HitranTag() )
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
              if ( (nChar == 161 && line[158] != ' ') || nChar > 161 )
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
      if(species_data[mspecies].Name() == "O2")
      {
        mlower_n = atoi(mlower_lquanta.substr(2, 3).c_str());
        mlower_j = atoi(mlower_lquanta.substr(6, 3).c_str());
        DJ =  -  mlower_lquanta.compare(5, 1, "Q");
        DN =  -  mlower_lquanta.compare(1, 1, "Q");
        mupper_n = mlower_n - DN;
        mupper_j = mlower_j - DJ;

        mquantum_numbers.SetLower(QN_N, mlower_n);
        mquantum_numbers.SetLower(QN_J, mlower_j);
        mquantum_numbers.SetUpper(QN_N, mlower_n - DN);
        mquantum_numbers.SetUpper(QN_J, mlower_j - DJ);
        mquantum_numbers.SetLower(QN_v1, atoi(mlower_gquanta.substr(13, 2).c_str()));
        mquantum_numbers.SetUpper(QN_v1, atoi(mupper_gquanta.substr(13, 2).c_str()));

        // Parse lower local quanta F
        String qnf = mlower_lquanta.substr(9, 5);
        qnf.trim();
        if (qnf.nelem())
        {
          ArrayOfString as;
          qnf.split(as, ".");
          if (as.nelem() == 2)
          {
            Index nom;
            char* endptr;

            nom = strtol(as[0].c_str(), &endptr, 10);
            if (endptr != as[0].c_str()+as[0].nelem())
                throw std::runtime_error("Error parsing quantum number F");

            if (as[1] == "5")
                mquantum_numbers.SetLower(QN_F, Rational(nom * 2 + 1, 2));
            else if (as[1] == "0")
                mquantum_numbers.SetLower(QN_F, nom);
            else
                throw std::runtime_error("Error parsing quantum number F");
          }
        }
      }
      else if(species_data[mspecies].Name() == "NO2" || species_data[mspecies].Name() == "HO2")
      {
          mlower_n = atoi(mlower_lquanta.substr(0, 3).c_str());
          mupper_n = atoi(mupper_lquanta.substr(0, 3).c_str());
          mquantum_numbers.SetUpper(QN_N, mupper_n);
          mquantum_numbers.SetLower(QN_N, mlower_n);
          
          if (mupper_lquanta[14] == '+')
            mquantum_numbers.SetUpper(QN_J,mupper_n+Rational(1,2));
          else if (mupper_lquanta[14] == '-')
            mquantum_numbers.SetUpper(QN_J,mupper_n-Rational(1,2));
          else
          { 
            // The J will be undefined and we fail at another stage.
          }
          
          if (mlower_lquanta[14] == '+')
              mquantum_numbers.SetLower(QN_J,mlower_n+Rational(1,2));
          else if (mlower_lquanta[14] == '-')
              mquantum_numbers.SetLower(QN_J,mlower_n-Rational(1,2));
          else
          { 
            // The J will be undefined and we fail at another stage.
          }
          
          mquantum_numbers.SetLower(QN_K1, atoi(mlower_lquanta.substr(3, 3).c_str()));
          mquantum_numbers.SetUpper(QN_K1, atoi(mupper_lquanta.substr(3, 3).c_str()));
          mquantum_numbers.SetLower(QN_K2, atoi(mlower_lquanta.substr(6, 3).c_str()));
          mquantum_numbers.SetUpper(QN_K2, atoi(mupper_lquanta.substr(6, 3).c_str()));
          
      }
      else if(species_data[mspecies].Name()=="NO"||species_data[mspecies].Name()=="ClO")
      {
        String qnf;
        
        // Parse upper local quanta J
        qnf = mlower_lquanta.substr(4, 5);
        qnf.trim();
        if (qnf.nelem())
        {
            ArrayOfString as;
            qnf.split(as, ".");
            if (as.nelem() == 2)
            {
                Index nom;
                char* endptr;
                
                nom = strtol(as[0].c_str(), &endptr, 10) + mlower_lquanta.compare(3,1,"Q");
                if (endptr != as[0].c_str()+as[0].nelem())
                    throw std::runtime_error("Error parsing quantum number J");
                
                if (as[1] == "5")
                    mquantum_numbers.SetUpper(QN_J, Rational(nom * 2 + 1 , 2));
                else if (as[1] == "0")
                    mquantum_numbers.SetUpper(QN_J, nom);
                else
                    throw std::runtime_error("Error parsing quantum number J");
            }
        }
        
        // Parse lower local quanta J
        qnf = mlower_lquanta.substr(4, 5);
        qnf.trim();
        if (qnf.nelem())
        {
            ArrayOfString as;
            qnf.split(as, ".");
            if (as.nelem() == 2)
            {
                Index nom;
                char* endptr;
                
                nom = strtol(as[0].c_str(), &endptr, 10);
                if (endptr != as[0].c_str()+as[0].nelem())
                    throw std::runtime_error("Error parsing quantum number J");
                
                if (as[1] == "5")
                    mquantum_numbers.SetLower(QN_J, Rational(nom * 2 + 1, 2));
                else if (as[1] == "0")
                    mquantum_numbers.SetLower(QN_J, nom);
                else
                    throw std::runtime_error("Error parsing quantum number J");
            }
        }
        
        // Parse lower local quanta F
        qnf = mlower_lquanta.substr(10, 5);
        qnf.trim();
        if (qnf.nelem())
        {
            ArrayOfString as;
            qnf.split(as, ".");
            if (as.nelem() == 2)
            {
                Index nom;
                char* endptr;
                
                nom = strtol(as[0].c_str(), &endptr, 10);
                if (endptr != as[0].c_str()+as[0].nelem())
                    throw std::runtime_error("Error parsing quantum number F");
                
                if (as[1] == "5")
                    mquantum_numbers.SetLower(QN_F, Rational(nom * 2 + 1, 2));
                else if (as[1] == "0")
                    mquantum_numbers.SetLower(QN_F, nom);
                else
                    throw std::runtime_error("Error parsing quantum number F");
            }
        }
        
        // Parse upper local quanta F
        qnf = mupper_lquanta.substr(10, 5);
        qnf.trim();
        if (qnf.nelem())
        {
            ArrayOfString as;
            qnf.split(as, ".");
            if (as.nelem() == 2)
            {
                Index nom;
                char* endptr;
                
                nom = strtol(as[0].c_str(), &endptr, 10);
                if (endptr != as[0].c_str()+as[0].nelem())
                    throw std::runtime_error("Error parsing quantum number F");
                
                if (as[1] == "5")
                    mquantum_numbers.SetUpper(QN_F, Rational(nom * 2 + 1, 2));
                else if (as[1] == "0")
                    mquantum_numbers.SetUpper(QN_F, nom);
                else
                    throw std::runtime_error("Error parsing quantum number F");
            }
        }  
        
        // Parse upper global quanta Omega
        qnf=mupper_gquanta.substr(8,3);
        qnf.trim();
        if (qnf.nelem())
        {
            ArrayOfString as;
            qnf.split(as, "/");
            if (as.nelem() == 2)
            {
                Index nom,denom;
                char* endptr;
                
                nom = strtol(as[0].c_str(), &endptr, 10);
                if (endptr != as[0].c_str()+as[0].nelem())
                    throw std::runtime_error("Error parsing nominator quantum number Omega");
                denom = strtol(as[1].c_str(), &endptr, 10);
                if (endptr != as[1].c_str()+as[1].nelem())
                    throw std::runtime_error("Error parsing denominator quantum number Omega");
                mquantum_numbers.SetUpper(QN_Omega, Rational(nom,denom));
            }
        }  
        
        // Parse lower global quanta Omega
        qnf=mlower_gquanta.substr(8,3);
        qnf.trim();
        if (qnf.nelem())
        {
            ArrayOfString as;
            qnf.split(as, "/");
            if (as.nelem() == 2)
            {
                Index nom,denom;
                char* endptr;
                
                nom = strtol(as[0].c_str(), &endptr, 10);
                if (endptr != as[0].c_str()+as[0].nelem())
                    throw std::runtime_error("Error parsing nominator quantum number Omega");
                denom = strtol(as[1].c_str(), &endptr, 10);
                if (endptr != as[1].c_str()+as[1].nelem())
                    throw std::runtime_error("Error parsing denominator quantum number Omega");
                mquantum_numbers.SetLower(QN_Omega, Rational(nom,denom));
            }
        }   
      }
      else if(species_data[mspecies].Name()=="OH")
      {
          // FIXME: What is wrong with HITRAN2004 format rules here? Seems like not only 2A1 on Br but also on Sym. 
          String qnf;
          
          // Parse upper local quanta J
          qnf = mlower_lquanta.substr(3, 5);
          qnf.trim();
          if (qnf.nelem())
          {
              ArrayOfString as;
              qnf.split(as, ".");
              if (as.nelem() == 2)
              {
                  Index nom;
                  char* endptr;
                  
                  nom = strtol(as[0].c_str(), &endptr, 10) + mlower_lquanta.compare(2,1,"Q");
                  if (endptr != as[0].c_str()+as[0].nelem())
                      throw std::runtime_error("Error parsing quantum number J");
                  
                  if (as[1] == "5")
                      mquantum_numbers.SetUpper(QN_J, Rational(nom * 2 + 1 , 2));
                  else if (as[1] == "0")
                      mquantum_numbers.SetUpper(QN_J, nom);
                  else
                      throw std::runtime_error("Error parsing quantum number J");
              }
          }
          
          // Parse lower local quanta J
          qnf = mlower_lquanta.substr(3, 5);
          qnf.trim();
          if (qnf.nelem())
          {
              ArrayOfString as;
              qnf.split(as, ".");
              if (as.nelem() == 2)
              {
                  Index nom;
                  char* endptr;
                  
                  nom = strtol(as[0].c_str(), &endptr, 10);
                  if (endptr != as[0].c_str()+as[0].nelem())
                      throw std::runtime_error("Error parsing quantum number J");
                  
                  if (as[1] == "5")
                      mquantum_numbers.SetLower(QN_J, Rational(nom * 2 + 1, 2));
                  else if (as[1] == "0")
                      mquantum_numbers.SetLower(QN_J, nom);
                  else
                      throw std::runtime_error("Error parsing quantum number J");
              }
          }
          
          // Parse lower local quanta F
          qnf = mlower_lquanta.substr(10, 5);
          qnf.trim();
          if (qnf.nelem())
          {
              ArrayOfString as;
              qnf.split(as, ".");
              if (as.nelem() == 2)
              {
                  Index nom;
                  char* endptr;
                  
                  nom = strtol(as[0].c_str(), &endptr, 10);
                  if (endptr != as[0].c_str()+as[0].nelem())
                      throw std::runtime_error("Error parsing quantum number F");
                  
                  if (as[1] == "5")
                      mquantum_numbers.SetLower(QN_F, Rational(nom * 2 + 1, 2));
                  else if (as[1] == "0")
                      mquantum_numbers.SetLower(QN_F, nom);
                  else
                      throw std::runtime_error("Error parsing quantum number F");
              }
          }
          
          // Parse upper local quanta F
          qnf = mupper_lquanta.substr(10, 5);
          qnf.trim();
          if (qnf.nelem())
          {
              ArrayOfString as;
              qnf.split(as, ".");
              if (as.nelem() == 2)
              {
                  Index nom;
                  char* endptr;
                  
                  nom = strtol(as[0].c_str(), &endptr, 10);
                  if (endptr != as[0].c_str()+as[0].nelem())
                      throw std::runtime_error("Error parsing quantum number F");
                  
                  if (as[1] == "5")
                      mquantum_numbers.SetUpper(QN_F, Rational(nom * 2 + 1, 2));
                  else if (as[1] == "0")
                      mquantum_numbers.SetUpper(QN_F, nom);
                  else
                      throw std::runtime_error("Error parsing quantum number F");
              }
          }  
          
          // Parse upper global quanta Omega
          qnf=mupper_gquanta.substr(8,3);
          qnf.trim();
          if (qnf.nelem())
          {
              ArrayOfString as;
              qnf.split(as, "/");
              if (as.nelem() == 2)
              {
                  Index nom,denom;
                  char* endptr;
                  
                  nom = strtol(as[0].c_str(), &endptr, 10);
                  if (endptr != as[0].c_str()+as[0].nelem())
                      throw std::runtime_error("Error parsing nominator quantum number Omega");
                  denom = strtol(as[1].c_str(), &endptr, 10);
                  if (endptr != as[1].c_str()+as[1].nelem())
                      throw std::runtime_error("Error parsing denominator quantum number Omega");
                  mquantum_numbers.SetUpper(QN_Omega, Rational(nom,denom));
              }
          }  
          
          // Parse lower global quanta Omega
          qnf=mlower_gquanta.substr(8,3);
          qnf.trim();
          if (qnf.nelem())
          {
              ArrayOfString as;
              qnf.split(as, "/");
              if (as.nelem() == 2)
              {
                  Index nom,denom;
                  char* endptr;
                  
                  nom = strtol(as[0].c_str(), &endptr, 10);
                  if (endptr != as[0].c_str()+as[0].nelem())
                      throw std::runtime_error("Error parsing nominator quantum number Omega");
                  denom = strtol(as[1].c_str(), &endptr, 10);
                  if (endptr != as[1].c_str()+as[1].nelem())
                      throw std::runtime_error("Error parsing denominator quantum number Omega");
                  mquantum_numbers.SetLower(QN_Omega, Rational(nom,denom));
              }
          }
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
  using global_data::species_data;

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
  using global_data::species_data;

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
  using global_data::species_data;
  
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
  using global_data::species_data;

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

      // Remaining entries are the quantum numbers
      getline(icecream, mquantum_numbers_str);
      mquantum_numbers_str.trim();

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

        // Do not write quantas from Hitran into ARTSCAT-4
        // because they're not compatible with our format
        // Only quantum numbers in Agnes' format are valid
//      ls << " " << lr.Upper_GQuanta()
//         << " " << lr.Lower_GQuanta()
//         << " " << lr.Upper_LQuanta()
//         << " " << lr.Lower_LQuanta();

      if (lr.QuantumNumbersString().nelem())
          ls << " " << lr.QuantumNumbersString();
      os << ls.str();
      
      break;
    }
      
    default:
      os << "Unknown ARTSCAT version: " << lr.Version();
      break;
  }
  
  return os;
}


//======================================================================
//         Functions for searches inside the line catalog
//======================================================================


bool find_matching_lines(ArrayOfIndex& matches,
                         const ArrayOfLineRecord& abs_lines,
                         const Index species,
                         const Index isotopologue,
                         const QuantumNumberRecord qr,
                         const LineMatchingCriteria match_criteria)
{
    bool ret = true;
    matches.resize(0);
    matches.reserve(100);

    for (Index l = 0; l < abs_lines.nelem(); l++)
    {
        const LineRecord& this_line = abs_lines[l];

        if ((species == -1 || this_line.Species() == species)
            && (isotopologue == -1 || this_line.Isotopologue() == isotopologue)
            && qr.Lower().Compare(this_line.QuantumNumbers().Lower())
            && qr.Upper().Compare(this_line.QuantumNumbers().Upper()))
        {
            matches.push_back(l);

            if (match_criteria == LINE_MATCH_FIRST)
                break;
            if (match_criteria == LINE_MATCH_UNIQUE && matches.nelem() > 1)
            {
                ret = false;
                break;
            }
        }
    }

    if (!matches.nelem()) ret = false;

    return ret;
}
