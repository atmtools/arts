/* Copyright (C) 2000, 2001 Stefan Buehler  <sbuehler@uni-bremen.de>
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

  \author Stefan Buehler
*/

#include "arts.h"
#include <map>
#include <cfloat>
#include <algorithm>
#include <math.h>
#include "absorption.h"
#include "math_funcs.h"
#include "auto_md.h"
#include "messages.h"


/** The map associated with species_data. */
std::map<String, Index> SpeciesMap;


// member fct of isotoperecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopeRecord::CalculatePartitionFctAtTemp( Numeric
                                                    temperature ) const
{
  Numeric result = 0.;
  Numeric exponent = 1.;

  ArrayOfNumeric::const_iterator it;

  for (it=mqcoeff.begin(); it != mqcoeff.end(); it++)
    {
      result += *it * exponent;
      exponent *= temperature;
    }
  return result;
}

void define_species_map()
{
  extern const Array<SpeciesRecord> species_data;
  extern std::map<String, Index> SpeciesMap;

  for ( Index i=0 ; i<species_data.nelem() ; ++i)
    {
      SpeciesMap[species_data[i].Name()] = i;
    }
}


ostream& operator << (ostream& os, const LineRecord& lr)
{
  // Determine the precision, depending on whether Numeric is double
  // or float:  
  Index precision;
  switch (sizeof(Numeric)) {
  case sizeof(float)  : precision = FLT_DIG; break;
  case sizeof(double) : precision = DBL_DIG; break;
  default: out0 << "Numeric must be double or float\n"; exit(1);
  }

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


bool LineRecord::ReadFromHitranStream(istream& is)
{
  // Global species lookup data:
  extern Array<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotope lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotope.
  const Index missing = species_data.nelem() + 100;

  // We need a species index sorted by HITRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<HITRAN tag>]. 
  //
  // Allow for up to 100 species in HITRAN in the future.
  static Array< Index >        hspec(100);

  // This is  an array of arrays for each hitran tag. It contains the
  // ARTS indices of the HITRAN isotopes. 
  static Array< ArrayOfIndex > hiso(100);

  // Remeber if this stuff has already been initialized:
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
          // HITRAN isotope tags are -1 (this species is missing in HITRAN).

          if ( 0 < sr.Isotope()[0].HitranTag() )
            {
              // The HITRAN tags are stored as species plus isotope tags
              // (MO and ISO)
              // in the Isotope() part of the species record.
              // We can extract the MO part from any of the isotope tags,
              // so we use the first one. We do this by taking an integer
              // division by 10.
          
              Index mo = sr.Isotope()[0].HitranTag() / 10;
              //          cout << "mo = " << mo << endl;
              hspec[mo] = i; 
          
              // Get a nicer to handle array of HITRAN iso tags:
              Index n_iso = sr.Isotope().nelem();
              ArrayOfIndex iso_tags;
              iso_tags.resize(n_iso);
              for ( Index j=0; j<n_iso; ++j )
                {
                  iso_tags[j] = sr.Isotope()[j].HitranTag();
                }

              // Reserve elements for the isotope tags. How much do we
              // need? This depends on the largest HITRAN tag that we know
              // about!
              // Also initialize the tags to missing.
              //          cout << "iso_tags = " << iso_tags << endl;
              //          cout << "static_cast<Index>(max(iso_tags))%10 + 1 = "
              //               << static_cast<Index>(max(iso_tags))%10 + 1 << endl;
              hiso[mo].resize( max(iso_tags)%10 + 1 );
              hiso[mo] = missing; // Matpack can set all elements like this.


              // Set the isotope tags:
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
                    out3 << " " << species_data[hspec[i]].Isotope()[hiso[i][j]].Name();
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
          if ( missing != hspec[mo] )       comment = false ;
          else
            {
              // See if this is already in warned_missing, use
              // std::count for that:
              if ( 0 == std::count(warned_missing.begin(),
                                   warned_missing.end(),
                                   mo) )
                {
                  out0 << "Error: HITRAN mo = " << mo << " is not "
                       << "known to ARTS.\n";
                  warned_missing.push_back(mo);
                  exit(1);
                  // SAB 08.08.2000 If you want to make the program
                  // continue anyway, just comment out the exit
                  // line.
                }
            }
        }
    }

  // Ok, we seem to have a valid species here.

  // Set mspecies from my cool index table:
  mspecies = hspec[mo];

  // Extract isotope:
  Index iso;                              
  extract(iso,line,1);
  //  cout << "iso = " << iso << endl;


  // Set misotope from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases. 
  misotope = missing;
  if ( iso < hiso[mo].nelem() )
    if ( missing != hiso[mo][iso] )
      misotope = hiso[mo][iso];

  // Issue error message if misotope is still missing:
  if (missing == misotope)
    {
      ostringstream os;
      os << "Species: " << species_data[mspecies].Name()
         << ", isotope iso = " << iso
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
    // 3. Take out the isotopic ratio.

    const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s,line,10);

    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    mi0 = s * hi2arts;

    // Take out isotopic ratio:
    mi0 /= species_data[mspecies].Isotope()[misotope].Abundance();    
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


bool LineRecord::ReadFromMytran2Stream(istream& is)
{
  // Global species lookup data:
  extern Array<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotope lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotope.
  const Index missing = species_data.nelem() + 100;

  // We need a species index sorted by MYTRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<MYTRAN tag>]. The value of
  // missing means that we don't have this species.
  //
  // Allow for up to 100 species in MYTRAN in the future.
  static Array< Index >        hspec(100,missing);      

  // This is  an array of arrays for each mytran tag. It contains the
  // ARTS indices of the MYTRAN isotopes. 
  static Array< ArrayOfIndex > hiso(100);

  // Remeber if this stuff has already been initialized:
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
          // MYTRAN isotope tags are -1 (this species is missing in MYTRAN).

          if ( 0 < sr.Isotope()[0].MytranTag() )
            {
              // The MYTRAN tags are stored as species plus isotope tags
              // (MO and ISO)
              // in the Isotope() part of the species record.
              // We can extract the MO part from any of the isotope tags,
              // so we use the first one. We do this by taking an integer
              // division by 10.
          
              Index mo = sr.Isotope()[0].MytranTag() / 10;
              //          cout << "mo = " << mo << endl;
              hspec[mo] = i; 
          
              // Get a nicer to handle array of MYTRAN iso tags:
              Index n_iso = sr.Isotope().nelem();
              ArrayOfIndex iso_tags;
              iso_tags.resize(n_iso);
              for ( Index j=0; j<n_iso; ++j )
                {
                  iso_tags[j] = sr.Isotope()[j].MytranTag();
                }

              // Reserve elements for the isotope tags. How much do we
              // need? This depends on the largest MYTRAN tag that we know
              // about!
              // Also initialize the tags to missing.
              //          cout << "iso_tags = " << iso_tags << endl;
              //          cout << "static_cast<Index>(max(iso_tags))%10 + 1 = "
              //               << static_cast<Index>(max(iso_tags))%10 + 1 << endl;
              hiso[mo].resize( max(iso_tags)%10 + 1 );
              hiso[mo] = missing; // Matpack can set all elements like this.

              // Set the isotope tags:
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
                    out3 << " " << species_data[hspec[i]].Isotope()[hiso[i][j]].Name();
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
      if (is.eof()) {
        //      cout << "Eof" << endl;
        return true;
      }

      // Throw runtime_error if stream is bad:
      if (!is) throw runtime_error ("Stream bad.");

      // Read line from file into linebuffer:
      getline(is,line);
      //      cout << line << endl;

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
                  out0 << "Error: MYTRAN mo = " << mo << " is not "
                       << "known to ARTS.\n";
                  warned_missing.push_back(mo);
                  exit(1);
                  // SAB 08.08.2000 If you want to make the program
                  // continue anyway, just comment out the exit
                  // line.
                }
            }
        }
    }

  // Ok, we seem to have a valid species here.

  // Set mspecies from my cool index table:
  mspecies = hspec[mo];

  // Extract isotope:
  Index iso;                              
  extract(iso,line,1);
  //  cout << "iso = " << iso << endl;


  // Set misotope from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases. 
  misotope = missing;
  if ( iso < hiso[mo].nelem() )
    if ( missing != hiso[mo][iso] )
      misotope = hiso[mo][iso];

  // Issue error message if misotope is still missing:
  if (missing == misotope)
    {
      ostringstream os;
      os << "Species: " << species_data[mspecies].Name()
         << ", isotope iso = " << iso
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
    // (just like HITRAN, only isotopic ratio is not included)
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
    extract(d,line,8);

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


bool LineRecord::ReadFromJplStream(istream& is)
{
  // Global species lookup data:
  extern Array<SpeciesRecord> species_data;

  // We need a species index sorted by JPL tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is JplMap[<JPL tag>]. We need Index in this map,
  // because the tag array within the species data is an Index array.
  static map<Index, SpecIsoMap> JplMap;

  // Remeber if this stuff has already been initialized:
  static bool hinit = false;

  if ( !hinit )
    {

      out3 << "  JPL index table:\n";

      for ( Index i=0; i<species_data.nelem(); ++i )
        {
          const SpeciesRecord& sr = species_data[i];


          for ( Index j=0; j<sr.Isotope().nelem(); ++j)
            {
              
              for (Index k=0; k<sr.Isotope()[j].JplTags().nelem(); ++k)
                {

                  SpecIsoMap indicies(i,j);

                  JplMap[sr.Isotope()[j].JplTags()[k]] = indicies;

                  // Print the generated data structures (for debugging):
                  // The explicit conversion of Name to a c-String is
                  // necessary, because setw does not work correctly for
                  // stl Strings.
                  const Index& i1 = JplMap[sr.Isotope()[j].JplTags()[k]].Speciesindex();
                  const Index& i2 = JplMap[sr.Isotope()[j].JplTags()[k]].Isotopeindex();
                                         
                  out3 << "  JPL TAG = " << sr.Isotope()[j].JplTags()[k] << "   Species = "
                       << setw(10) << setiosflags(ios::left)
                       << species_data[i1].Name().c_str()
                       << "iso = " 
                       << species_data[i1].Isotope()[i2].Name().c_str()
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
      if (is.eof()) {
        //      cout << "Eof" << endl;
        return true;
      }

      // Throw runtime_error if stream is bad:
      if (!is) throw runtime_error ("Stream bad.");

      // Read line from file into linebuffer:
      getline(is,line);

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
    s = pow(10,s);

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
    tag = abs(tag);
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

  // Set misotope:
  misotope = id.Isotopeindex();

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


bool LineRecord::ReadFromArtsStream(istream& is)
{
  // Global species lookup data:
  extern Array<SpeciesRecord> species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>]. 
  static map<String, SpecIsoMap> ArtsMap;

  // Remeber if this stuff has already been initialized:
  static bool hinit = false;

  if ( !hinit )
    {

      out3 << "  ARTS index table:\n";

      for ( Index i=0; i<species_data.nelem(); ++i )
        {
          const SpeciesRecord& sr = species_data[i];


          for ( Index j=0; j<sr.Isotope().nelem(); ++j)
            {
              
              SpecIsoMap indicies(i,j);
              String buf = sr.Name()+"-"+sr.Isotope()[j].Name();
              
              ArtsMap[buf] = indicies;
              
              // Print the generated data structures (for debugging):
              // The explicit conversion of Name to a c-String is
              // necessary, because setw does not work correctly for
              // stl Strings.
              const Index& i1 = ArtsMap[buf].Speciesindex();
              const Index& i2 = ArtsMap[buf].Isotopeindex();
                                         
              out3 << "  Arts Identifier = " << buf << "   Species = "
                   << setw(10) << setiosflags(ios::left)
                   << species_data[i1].Name().c_str()
                   << "iso = " 
                   << species_data[i1].Isotope()[i2].Name().c_str()
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
      if (is.eof()) {
        //      cout << "Eof" << endl;
        return true;
      }

      // Throw runtime_error if stream is bad:
      if (!is) throw runtime_error ("Stream bad.");

      // Read line from file into linebuffer:
      getline(is,line);

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

      // Set misotope:
      misotope = id.Isotopeindex();


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

      for (Index i = 0; i<naux; i++)
        {
          icecream >> maux[i];
          //cout << "maux" << i << " = " << maux[i] << "\n";
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
      catch (runtime_error x)
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



OneTag::OneTag(String def) 
{
  // Species lookup data:
  extern const Array<SpeciesRecord> species_data;
  // The species map. This is used to find the species id.
  extern std::map<String, Index> SpeciesMap;
  // Name of species and isotope (aux variables):
  String name, isoname;
  // Aux index:
  Index n;


  // Set frequency limits to default values (no limits):
  mlf = -1;
  muf = -1;

  // We cannot set a default value for the isotope, because the
  // default should be `ALL' and the value for `ALL' depends on the
  // species. 
    

  // Extract the species name:
  n    = def.find('-');    // find the '-'
  if (n != def.npos )
    {
      name = def.substr(0,n);      // Extract before '-'
      def.erase(0,n+1);             // Remove from def
    }
  else
    {
      // n==def.npos means that def does not contain a '-'. In that
      // case we assume that it contains just a species name and
      // nothing else
      name = def;
      def  = "";
    }


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
      ostringstream os;
      os << "Species " << name << " is not a valid species.";
      throw runtime_error(os.str());
    }

  // Get a reference to the relevant Species Record:
  const SpeciesRecord& spr = species_data[mspecies];

  if ( 0 == def.nelem() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants all isotopes and no frequency limits.
      // Frequency defaults are already set. Set isotope defaults:
      misotope = spr.Isotope().nelem();
      // This means all isotopes.

      return;
    }
    
  // Extract the isotope name:
  n    = def.find('-');    // find the '-'
  if (n != def.npos )
    {
      isoname = def.substr(0,n);          // Extract before '-'
      def.erase(0,n+1);             // Remove from def
    }
  else
    {
      // n==def.npos means that def does not contain a '-'. In that
      // case we assume that it contains just the isotope name and
      // nothing else.
      isoname = def;
      def  = "";
    }
    
  // Check for joker:
  if ( "*" == isoname )
    {
      // The user wants all isotopes. Set this accordingly:
      misotope = spr.Isotope().nelem();
    }
  else
    {
      // Make an array containing the isotope names:
      ArrayOfString ins;
      for ( Index i=0; i<spr.Isotope().nelem(); ++i )
        ins.push_back( spr.Isotope()[i].Name() );

      // Use the find algorithm from the STL to find the isotope ID. It
      // returns an iterator, so to get the index we take the
      // difference to the begin() iterator.
      misotope = find( ins.begin(),
                       ins.end(),
                       isoname ) - ins.begin();

      // Check if we found a matching isotope:
      if ( misotope >= ins.nelem() ) 
        {
          ostringstream os;
          os << "Isotope " << isoname << " is not a valid isotope for "
             << "species " << name << ".\n"
             << "Valid isotopes are:";
          for ( Index i=0; i<ins.nelem(); ++i )
            os << " " << ins[i];
          throw runtime_error(os.str());
        }
    }

  if ( 0 == def.nelem() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants no frequency limits.  Frequency defaults are
      // already set, so we can return directly.

      return;
    }


  // Look for the two frequency limits:
  
  // Extract first frequency
  n    = def.find('-');    // find the '-'
  if (n != def.npos )
    {
      // Frequency as a String:
      String fname;
      fname = def.substr(0,n);              // Extract before '-'
      def.erase(0,n+1);               // Remove from def

      // Check for joker:
      if ( "*" == fname )
        {
          // The default for mlf is already -1, meaning `ALL'.
          // So there is nothing to do here.
        }
      else
        {
          // Convert to Numeric:
          istringstream is(fname);
          is >> mlf;
        }
    }
  else
    {
      // n==def.npos means that def does not contain a '-'. In this
      // case that is not allowed!
      throw runtime_error("You must either speciefy both frequency limits\n"
                          "(at least with jokers), or none.");
    }


  // Now there should only be the upper frequency left in def.
  // Check for joker:
  if ( "*" == def )
    {
      // The default for muf is already -1, meaning `ALL'.
      // So there is nothing to do here.
    }
  else
    {
      // Convert to Numeric:
      istringstream is(def);
      is >> muf;
    }
}


String OneTag::Name() const 
{
  // Species lookup data:
  extern const Array<SpeciesRecord> species_data;
  // A reference to the relevant record of the species data:
  const  SpeciesRecord& spr = species_data[mspecies];
  // For return value:
  ostringstream os;

  // First the species name:
  os << spr.Name() << "-";

  // Now the isotope. Can be a single isotope or ALL.
  if ( misotope == spr.Isotope().nelem() )
    {
      // One larger than allowed means all isotopes!
      os << "*-";
    }
  else
    {
      os << spr.Isotope()[misotope].Name() << "-";
    }

  // Now the frequency limits, if there are any. For this we first
  // need to determine the floating point precision.

  // Determine the precision, depending on whether Numeric is double
  // or float:  
  Index precision;
  switch (sizeof(Numeric)) {
  case sizeof(float)  : precision = FLT_DIG; break;
  case sizeof(double) : precision = DBL_DIG; break;
  default: out0 << "Numeric must be double or float\n"; exit(1);
  }

  if ( 0 > mlf )
    {
      // mlf < 0 means no lower limit.
      os << "*-";
    }
  else
    {
      os << setprecision(precision);
      os << mlf << "-";
    }

  if ( 0 > muf )
    {
      // muf < 0 means no upper limit.
      os << "*";
    }
  else
    {
      os << setprecision(precision);
      os << muf;
    }

  return os.str();
}

ostream& operator << (ostream& os, const OneTag& ot)
{
  return os << ot.Name();
}



/** 
   Returns the index among some tag groups for an array of tag Strings. 
   
   \begin{verbatim}
   For example, if tags1 correspond to the definition
     ["O3","H2O-161,H2O-162"]
   and the tag Strings are
     ["H2O-161,H2O-162","O3"]
   the tags1_index becomes
     [2,1]
   \end{verbatim}

   @exception runtime_error  Some String is not a valid tag item.
   @exception runtime_error  Not all Strings are not found among the tags.

   \retval tags1_index     Index in tags1 for tags2_Strings
   \param  tags1           The tags to search in.
   \param  tags2_Strings   The tag Strings for which indeces shall be found.

   \author Patrick Eriksson 
   \date 2000-12-06
*/
void get_tagindex_for_Strings( 
                              ArrayOfIndex&   tags1_index, 
                              const TagGroups&      tags1, 
                              const ArrayOfString&  tags2_Strings )
{
  const Index   n1 = tags1.nelem();
  const Index   n2 = tags2_Strings.nelem();
     TagGroups   tags2;                // Internal tag names for tag_Strings
        Index   i1, i2, nj, j, found, ok;

  tags1_index.resize(n2);
  //  cout << "tags2_Strings: " << tags2_Strings << "\n";
  tgsDefine( tags2, tags2_Strings );

  for ( i2=0; i2<n2; i2++ )
  {
    found = 0;
    for ( i1=0; (i1<n1) && !found; i1++ )
    {
      nj = tags2[i2].nelem(); 
      if ( nj  == tags1[i1].nelem() )
      {
        ok = 1;
        for ( j=0; j<nj; j++ )
        {
          if ( tags2[i2][j].Name() != tags1[i1][j].Name() )
            ok = 0;
        }
        if ( ok )
        {
          found = 1;
          tags1_index[i2] = i1;
        }
      }
    } 
    if ( !found )
    {
      ostringstream os;
      os << "The tag String \"" << tags2_Strings[i2] << 
            "\" does not match any of the given tags.\n";
      throw runtime_error(os.str());
    }
  }
}


/** Returns the index of the tag group tg2 within the array of tag
    groups tgs1. Slightly modified copy of get_tagindex_for_Strings.
   
    @exception runtime_error  Could not find tg2 in tgs1.
   
    \retval tgs1_index     Index in tgs1 for tg2
    \param  tgs1           The tags groups to search in.
    \param  tg2            The tag group for which the index shall be found.
    
    \author Patrick Eriksson, Axel von Engeln, and Stefan Buehler
    \date 2001-01-31 */
void get_tag_group_index_for_tag_group( Index&               tgs1_index, 
                                        const TagGroups&      tgs1, 
                                        const Array<OneTag>&  tg2 )
{
  bool found = false;

  for ( Index i=0;
        i<tgs1.nelem() && !found;
        ++i )
    {
      // Is at least the size correct?
      if ( tg2.nelem() == tgs1[i].nelem() )
        {
          bool ok = true;

          for ( Index j=0; j<tg2.nelem(); ++j )
            {
              if ( tg2[j].Name() != tgs1[i][j].Name() )
                ok = false;
            }

          if ( ok )
            {
              found = true;
              tgs1_index = i;
            }
        }
    }

  if ( !found )
    {
      ostringstream os;
      os << "The tag String \"" << tg2 << 
        "\" does not match any of the given tags.\n";
      throw runtime_error(os.str());
    }
}


/**
   Print the name of a tag group. 

   A tag group consists of several elementary OneTags. This functions
   returns a String with the name of the entire tag group. This is
   nice for informational output messages, for example in the
   absorption routines.

   \param  tg  The tag group in question.
   \return The full name of the tag group, as it could occur in the controlfile.

   \author Stefan Buehler
   \date   2001-03-13
*/
String get_tag_group_name( const Array<OneTag>& tg )
{
  String name;
  Index i;
  
  for ( i=0; i<tg.nelem()-1; ++i )
    {
      name += tg[i].Name() + ", ";
    }
  name += tg[i].Name();

  return name;
}


/** A helper function that writes lines in a line list to a stream. 

    \retval os      The stream to write to.
    \param  lines   The line list to write.

    \date 2000-06-12 
    \author Stefan Buehler 
*/
void write_lines_to_stream(ostream& os,
                           const ArrayOfLineRecord& lines)
{
  // We need this dummy line record, so that we can get the catalogue
  // version tag from dummy.Version, even if the line list is empty.
  LineRecord dummy;

  // Output version String first, followed by number of lines:
  os << dummy.Version() << " " << lines.nelem() << "\n";

  // Now output the line data:
  for ( Index i=0; i<lines.nelem(); ++i )
    {
      //      out3 << lines[i] << "\n";
      os << lines[i] << "\n";
    }
}


//! Checks if a vector is sorted in ascending order.
/*!
  Duplicated values are allowed.

  We need this to make sure that the frequency grid is sorted in the
  case of lineshape with cutoff. Duplicate values are allowed.

  THIS FUNCTION IS TAKEN FROM ARTS-1-1.

  \param   x   A vector.
  \return      True if sorted.
*/
bool is_sorted( ConstVectorView   x )
{
  if( x.nelem() > 1 )
    {
      for( Index i=1; i<x.nelem(); i++ )
        {
          if( x[i] < x[i-1] )
            return false;
        }
    }
  return true;
}


/** Calculate line absorption cross sections for one tag group. All
    lines in the line list must belong to the same species. This must
    be ensured by lines_per_tgCreateFromLines, so it is only verified
    with assert. Also, the input vectors p_abs, and t_abs must all
    have the same dimension.

    This is mainly a copy of abs_species which is removed now, with
    the difference that the vmrs are removed from the absorption
    coefficient calculation. (the vmr is still used for the self
    broadening)

    Continua are not handled by this function, you have to call
    xsec_continuum_tag for those.

    \retval xsec   Cross section of one tag group.
    \param f_mono  Frequency grid.
    \param p_abs   Pressure grid.
    \param t_abs   Temperatures associated with p_abs.
    \param h2o_abs Total volume mixing ratio of water vapor.
    \param vmr     Volume mixing ratio of the calculated species.
    \param lines   The spectroscopic line list.
    \param ind_ls  Lineshape specifications.

    \author Stefan Buehler and Axel von Engeln
    \date   2001-01-11 */
void xsec_species( MatrixView               xsec,
                   ConstVectorView          f_mono,
                   ConstVectorView          p_abs,
                   ConstVectorView          t_abs,
                   ConstVectorView          h2o_abs,
                   ConstVectorView          vmr,
                   const ArrayOfLineRecord& lines,
                   const Index              ind_ls,
                   const Index              ind_lsn,
                   const Numeric            cutoff)
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
  static const Numeric doppler_const = sqrt( 2.0 * BOLTZMAN_CONST *
                                             AVOGADROS_NUMB) / SPEED_OF_LIGHT; 

  // dimension of f_mono, lines
  Index nf = f_mono.nelem();
  Index nl = lines.nelem();

  // Define the vector for the line shape function and the
  // normalization factor of the lineshape here, so that we don't need
  // so many free store allocations.  the last element is used to
  // calculate the value at the cutoff frequency
  Vector ls(nf+1);
  Vector fac(nf+1);

  bool cut = (cutoff != -1) ? true : false;

  // Check that the frequency grid is sorted in the case of lineshape
  // with cutoff. Duplicate frequency values are allowed.
  if (cut)
    {
      if ( ! is_sorted( f_mono ) )
        {
          ostringstream os;
          os << "If you use a lineshape function with cutoff, your\n"
             << "frequency grid *f_mono* must be sorted.\n"
             << "(Duplicate values are allowed.)";
          throw runtime_error(os.str());
        }
    }
  
  // Check that all temperatures are non-negative
  bool negative = false;
  
  for (Index i = 0; !negative && i < t_abs.nelem (); i++)
    {
      if (t_abs[i] < 0.)
        negative = true;
    }
  
  if (negative)
    {
      ostringstream os;
      os << "t_abs contains at least one negative temperature value.\n"
         << "This is not allowed.";
      throw runtime_error(os.str());
    }
  
  // We need a local copy of f_mono which is 1 element longer, because
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

  // Check that p_abs, t_abs, and h2o_abs all have the same
  // dimension. This could be a user error, so we throw a
  // runtime_error.  FIXME: why do we do this for each tag? wouldn't a
  // check in abscalc be sufficient?

  if ( t_abs.nelem() != p_abs.nelem() )
    {
      ostringstream os;
      os << "Variable t_abs must have the same dimension as p_abs.\n"
         << "t_abs.nelem() = " << t_abs.nelem() << '\n'
         << "p_abs.nelem() = " << p_abs.nelem();
      throw runtime_error(os.str());
    }

  if ( vmr.nelem() != p_abs.nelem() )
    {
      ostringstream os;
      os << "Variable vmr must have the same dimension as p_abs.\n"
         << "vmr.nelem() = " << vmr.nelem() << '\n'
         << "p_abs.nelem() = " << p_abs.nelem();
      throw runtime_error(os.str());
    }

  if ( h2o_abs.nelem() != p_abs.nelem() )
    {
      ostringstream os;
      os << "Variable h2o_abs must have the same dimension as p_abs.\n"
         << "h2o_abs.nelem() = " << h2o_abs.nelem() << '\n'
         << "p_abs.nelem() = " << p_abs.nelem();
      throw runtime_error(os.str());
    }


  // Check that the dimension of xsec is indeed [f_mono.nelem(),
  // p_abs.nelem()]:
  if ( xsec.nrows() != nf || xsec.ncols() != p_abs.nelem() )
    {
      ostringstream os;
      os << "Variable xsec must have dimensions [f_mono.nelem(),p_abs.nelem()].\n"
         << "[xsec.nrows(),xsec.ncols()] = [" << xsec.nrows()
         << ", " << xsec.ncols() << "]\n"
         << "f_mono.nelem() = " << nf << '\n'
         << "p_abs.nelem() = " << p_abs.nelem();
      throw runtime_error(os.str());
    }


  // Loop all pressures:
  for ( Index i=0; i<p_abs.nelem(); ++i )
    {

      // store variables p_abs[i] and t_abs[i],
      // this is slightly faster
      Numeric p_i = p_abs[i];
      Numeric t_i = t_abs[i];

    
      //out3 << "  p = " << p_i << " Pa\n";

      // Calculate total number density from pressure and temperature.
      // n = n0*T0/p0 * p/T or n = p/kB/t, ideal gas law
      Numeric n = p_i / BOLTZMAN_CONST / t_i;

      // For the pressure broadening, we also need the partial pressure:
      const Numeric p_partial = p_i * vmr[i];


      // Loop all lines:
      for ( Index l=0; l< nl; ++l )
        {
          // Copy f_mono to the beginning of f_local. There is one
          // element left at the end of f_local.  
          // THIS HAS TO BE INSIDE THE LINE LOOP, BECAUSE THE CUTOFF
          // FREQUENCY IS ALWAYS PUT IN A DIFFERENT PLACE!
          f_local[Range(0,nf)] = f_mono;

          // This will hold the actual number of frequencies to add to
          // xsec later on:
          Index nfl = nf;

          // This will hold the actual number of frequencies for the
          // call to the lineshape functions later on:
          Index nfls = nf;      

          // The baseline to substract for cutoff frequency
          Numeric base=0.0;

          // lines[l] is used several times, this construct should be
          // faster (Oliver Lemke)
          LineRecord l_l = lines[l];  // which line are we dealing with
          // Center frequency in vacuum:
          Numeric F0 = l_l.F();

          // Intensity is already in the right units (Hz*m^2). It also
          // includes already the isotopic ratio. Needs only to be
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
            l_l.IsotopeData().CalculatePartitionFctRatio( l_l.Ti0(),
                                                          t_i );

          // Boltzmann factors
          Numeric nom = exp(- e_lower / ( BOLTZMAN_CONST * t_i ) ) - 
            exp(- e_upper / ( BOLTZMAN_CONST * t_i ) );

          Numeric denom = exp(- e_lower / ( BOLTZMAN_CONST * l_l.Ti0() ) ) - 
            exp(- e_upper / ( BOLTZMAN_CONST * l_l.Ti0() ) );


          // intensity at temperature
          intensity *= part_fct_ratio * nom / denom;


          // 2. Get pressure broadened line width:
          // (Agam is in Hz/Pa, p_abs is in Pa, gamma is in Hz)
          const Numeric theta = l_l.Tgam() / t_i;
          const Numeric theta_Nair = pow(theta, l_l.Nair());

          Numeric gamma
            = l_l.Agam() * theta_Nair  * (p_i - p_partial)
            + l_l.Sgam() * pow(theta, l_l.Nself()) * p_partial;

          // 3. Doppler broadening without the sqrt(ln(2)) factor, which
          // seems to be redundant FIXME: verify .
          Numeric sigma = F0 * doppler_const * 
            sqrt( t_i / l_l.IsotopeData().Mass());

          // 3.a. Put in pressure shift.
          // The T dependence is connected to that of agam by:
          // n_shift = .25 + 1.5 * n_agam
          // Theta has been initialized above.
          F0 += l_l.Psf() * p_i * 
              std::pow( theta , (Numeric).25 + (Numeric)1.5*l_l.Nair() );

          // 4. the rosenkranz lineshape for oxygen calculates the
          // pressure broadening, overlap, ... differently. Therefore
          // all required parameters are passed in the aux array, the
          // given order must agree with how they are accessed in the
          // used lineshape (currently only the Rosenkranz routines are
          // using this method of passing parameters). Parameters are
          // always passed, because the Rosenkranz lineshape function
          // can be used without overlap correction, which is then just
          // set to 0. I know, this is not very nice, but I suggested to
          // put the oxygen absorption into a different workspace
          // method, but nobody listened to me, Schnief.
          aux[0] = theta;
          aux[1] = theta_Nair;
          aux[2] = p_i;
          aux[3] = vmr[i];
          aux[4] = h2o_abs[i];
          aux[5] = l_l.Agam();
          aux[6] = l_l.Nair();
          // is overlap available, otherwise pass zero
          if (l_l.Naux() > 1)
            {
              aux[7] = l_l.Aux()[0];
              aux[8] = l_l.Aux()[1];
              //            cout << "aux7, aux8: " << aux[7] << " " << aux[8] << "\n";
            }
          else
            {
              aux[7] = 0.0;
              aux[8] = 0.0;
            }


          // Indices pointing at begin/end frequencies of f_mono or at
          // the elements that have to be calculated in case of cutoff
          Index i_f_min = 0;            
          Index i_f_max = nf-1;         

          // cutoff ?
          if ( cut )
            {
              // Check whether we have elements in ls that can be
              // ignored at lower frequencies of f_mono.
              //
              // Loop through all frequencies, finding min value and
              // set all values to zero on that way.
              while ( i_f_min < nf && (F0 - cutoff) > f_mono[i_f_min] )
                {
                  //              ls[i_f_min] = 0;
                  ++i_f_min;
                }
              

              // Check whether we have elements in ls that can be
              // ignored at higher frequencies of f_mono.
              //
              // Loop through all frequencies, finding max value and
              // set all values to zero on that way.
              while ( i_f_max >= 0 && (F0 + cutoff) < f_mono[i_f_max] )
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

          //      cout << "nf, nfl, nfls = " << nf << ", " << nfl << ", " << nfls << ".\n";
        
          // Maybe there are no frequencies left to compute?  Note that
          // the number that counts here is nfl, since only these are
          // the `real' frequencies, for which xsec is changed. nfls
          // will always be at least one, because it contains the cutoff.
          if ( nfl > 0 )
            {
              // Calculate the line shape:
              lineshape_data[ind_ls].Function()(ls,
                                                aux,F0,gamma,sigma,
                                                f_local[Range(i_f_min,nfls)],
                                                nfls);

              // Calculate the chosen normalization factor:
              lineshape_norm_data[ind_lsn].Function()(fac,F0,
                                                      f_local[Range(i_f_min,nfls)],
                                                      nfls);

              // Get a handle on the range of xsec that we want to change.
              // We use nfl here, which could be one less than nfls.
              VectorView this_xsec = xsec(Range(i_f_min,nfl),i);

	      // Get handles on the range of ls and fac that we need.
	      VectorView this_ls  = ls[Range(0,nfl)];
	      VectorView this_fac = fac[Range(0,nfl)];

              // cutoff ?
              if ( cut )
                {
                  // The index nfls-1 should be exactly the index pointing
                  // to the value at the cutoff frequency.
                  base = ls[nfls-1];

                  // Subtract baseline from xsec. 
                  // this_xsec -= base;
                  this_ls -= base;
                }

              // Add line to xsec. 
              {
                // To make the loop a bit faster, precompute all constant
                // factors. These are:
                // 1. Total number density of the air.
                // 2. Line intensity.
                // 3. Isotopic ratio.
                //
                // The isotopic ratio must be applied here, since we are
                // summing up lines belonging to different isotopes.

                const Numeric factors = n * intensity * l_l.IsotopeData().Abundance();

                // We have to do:
                // xsec(j,i) += factors * ls[j1] * fac[j1];
                //
                // We use ls as a dummy to compute the product, then add it
                // to this_xsec.

                // FIXME: Maybe try if this is faster with a good
                // old-fashioned simple loop.

                this_ls *= this_fac;
                this_ls *= factors;
                this_xsec += this_ls;
              }
            }
        }
    }
}




//======================================================================
//             Functions related to refraction
//======================================================================


//// refr_index_BoudourisDryAir ///////////////////////////////////////////////
/**
   Calculates the refractive index for dry air at microwave frequncies 
   following Boudouris 1963.

   The expression is also found in Chapter 5 of the Janssen book.

   The atmosphere is assumed to have no water vapour.

   \retval   refr_index  refractive index
   \param    p_abs       absorption pressure grid
   \param    t_abs       temperatures at p_abs

   \author Patrick Eriksson
   \date   2001-02-16
*/
void refr_index_BoudourisDryAir (
                                Vector&   refr_index,
                                ConstVectorView   p_abs,
                                ConstVectorView   t_abs )
{
  const Index   n = p_abs.nelem();
  refr_index.resize( n );

  assert ( n == t_abs.nelem() );

  // N = 77.593e-2 * p / t ppm
  for ( Index i=0; i<n; i++ )
    refr_index[i] = 1.0 + 77.593e-8 * p_abs[i] / t_abs[i];
}



//// refr_index_Boudouris ///////////////////////////////////////////////
/**
   Calculates the refractive index at microwave frequncies 
   following Boudouris 1963.

   The expression is also found in Chapter 5 of the Janssen book.

   \retval   refr_index  refractive index
   \param    p_abs       absorption pressure grid
   \param    t_abs       temperatures at p_abs
   \param    h2o_abs     H2O vmr at p_abs

   \author Patrick Eriksson
   \date   2001-02-16
*/
void refr_index_Boudouris (
                          Vector&   refr_index,
                          ConstVectorView   p_abs,
                          ConstVectorView   t_abs,
                          ConstVectorView   h2o_abs )
{
  const Index   n = p_abs.nelem();
  refr_index.resize( n );

  assert ( n == t_abs.nelem() );
  assert ( n == h2o_abs.nelem() );

  Numeric   e;     // Partial pressure of water in Pa
  Numeric   p;     // Partial pressure of the dry air: p = p_tot - e 

  for ( Index i=0; i<n; i++ )
  {
    e = p_abs[i] * h2o_abs[i];
    p = p_abs[i] - e;

    refr_index[i] = 1.0 + 77.593e-8 * p / t_abs[i] + 
                          72e-8 * e / t_abs[i] +
                          3.754e-3 * e / (t_abs[i]*t_abs[i]);
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
  static const Numeric lower_energy_const = PLANCK_CONST * SPEED_OF_LIGHT * 1E2;

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

