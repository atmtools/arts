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

#include <algorithm>
#include "arts.h"
#include <map>
#include <cfloat>
#include <cmath>
#include "absorption.h"
#include "math_funcs.h"
#include "auto_md.h"
#include "messages.h"


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

ostream& operator << (ostream& os, const LineRecord& lr)
{
  // Determine the precision, depending on whether Numeric is double
  // or float:  
  Index precision;
  switch (sizeof(Numeric)) {
  case sizeof(float)  : precision = FLT_DIG; break;
  case sizeof(double) : precision = DBL_DIG; break;
  default: out0 << "Numeric must be double or float\n"; arts_exit ();
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
     << " " << lr.Naux  ();

  for ( Index i=0; i<lr.Naux(); ++i )
    os << " " << lr.Aux()[i];

  return os;
}


/** Extract something from a catalogue line. This is just a small helper
    function to safe some typing. 

    \param x Output:    What was extracted from the beginning of the line.
    \param line Output: What was extracted is also cut away from line.
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
              // count for that:
              if ( 0 == count(warned_missing.begin(),
                                   warned_missing.end(),
                                   mo) )
                {
                  out0 << "Error: HITRAN mo = " << mo << " is not "
                       << "known to ARTS.\n";
                  warned_missing.push_back(mo);
                  arts_exit ();
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


  // These were all the parameters that we can extract from
  // HITRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for HITRAN)
  mti0 = 296.0;

  // calculate the partition fct at the reference temperature. This is
  // done for all lines read from the catalogue, even the ones never
  // used in calculation. FIXME: are there better ways to implement this?
  species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 );

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
              // count for that:
              if ( 0 == count(warned_missing.begin(),
                                   warned_missing.end(),
                                   mo) )
                {
                  out0 << "Error: MYTRAN mo = " << mo << " is not "
                       << "known to ARTS.\n";
                  warned_missing.push_back(mo);
                  arts_exit ();
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

  // skip transition frequency error
  {
    Numeric r;
    extract(r,line,8);
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


  // These were all the parameters that we can extract from
  // MYTRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for MYTRAN2)
  mti0 = 296.0;

  // calculate the partition fct at the reference temperature. This is
  // done for all lines read from the catalogue, even the ones never
  // used in calculation. FIXME: are there better ways to implement this?
  species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 );



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


  // skip transition frequency error
  {
    Numeric r;
    extract(r,line,8);
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

  // calculate the partition fct at the reference temperature. This is
  // done for all lines read from the catalogue, even the ones never
  // used in calculation. FIXME: are there better ways to implement this?
  species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 );



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


      // calculate the partition fct at the reference temperature. This is
      // done for all lines read from the catalogue, even the ones never
      // used in calculation. FIXME: are there better ways to implement this?
      species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 );

    }

  // That's it!
  return false;
}



/** 
   Returns the index among some tag groups for an array of tag Strings. 
   
   <pre>
   For example, if tags1 correspond to the definition
     ["O3","H2O-161,H2O-162"]
   and the tag Strings are
     ["H2O-161,H2O-162","O3"]
   the tags1_index becomes
     [2,1]
   </pre>

   @exception runtime_error  Some String is not a valid tag item.
   @exception runtime_error  Not all Strings are not found among the tags.

   \param tags1_index Output:     Index in tags1 for tags2_Strings
   \param  tags1           The tags to search in.
   \param  tags2_Strings   The tag Strings for which indeces shall be found.

   \author Patrick Eriksson 
   \date 2000-12-06
*/
void get_tagindex_for_Strings( 
                              ArrayOfIndex&   tags1_index, 
                              const ArrayOfArrayOfSpeciesTag&      tags1, 
                              const ArrayOfString&  tags2_Strings )
{
  const Index   n1 = tags1.nelem();
  const Index   n2 = tags2_Strings.nelem();
     ArrayOfArrayOfSpeciesTag   tags2;                // Internal tag names for tag_Strings
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
   
    FIXME: SAB 2002-11-29 I think we will not need this for ARTS-1-1,
    instead one can use the more general function chk_contains.

    @exception runtime_error  Could not find tg2 in tgs1.
   
    \param tgs1_index Output:     Index in tgs1 for tg2
    \param  tgs1           The tags groups to search in.
    \param  tg2            The tag group for which the index shall be found.
    
    \author Patrick Eriksson, Axel von Engeln, and Stefan Buehler
    \date 2001-01-31 */
void get_tag_group_index_for_tag_group( Index&               tgs1_index, 
                                        const ArrayOfArrayOfSpeciesTag&      tgs1, 
                                        const Array<SpeciesTag>&  tg2 )
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


/** A helper function that writes lines in a line list to a stream. 

    \param os Output:      The stream to write to.
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

    \param xsec Output:   Cross section of one tag group.
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

  // We need a local copy of f_mono which is 1 element longer, because
  // we append a possible cutoff to it.
  bool cut = (cutoff != -1) ? true : false;
  Index nfl = nf;                       // This will hold the actual number of
                                        // frequencies to add to xsec
                                        // later on.
  Index nfls = nf;                      // This will hold the actual
                                        // nubmer of frequencies for
                                        // the call to the lineshape
                                        // functions later on.

  Vector f_local( nfl );
  f_local[Range(0,nf)] = f_mono; // Copy f_mono to the beginning of
                                 // f_local. There could be one
                                 // element left at the end of f_local.

  // The baseline to substract for cutoff frequency
  Numeric base=0.0;

  // Voigt generally needs a different frequency grid. If we allocate
  // that in the outer loop, insted of in voigt, we don't have the
  // free store allocation at each lineshape call. Calculation is
  // still done in the voigt routine itself, this is just an auxillary
  // parameter, passed to lineshape. For selected lineshapes (e.g.,
  // Rosenkranz) it is used additionally to pass parameters needed in
  // the lineshape (e.g., overlap, ...). Consequently we have to
  // assure that aux has a dimension not less then the number of
  // parameters passed.
  Index ii = (nf < 10) ? 10 : nf;
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
          Numeric part_fct_ratio =
            l_l.IsotopeData().CalculatePartitionFctRatio( t_i );

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
              pow( theta , .25 + 1.5*l_l.Nair() );         

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
              // Check whether we have elements in ls that have to be
              // set to zero at lower frequencies of f_mono.
              //
              // Loop through all frequencies, finding min value and
              // set all values to zero on that way.
              while ( i_f_min < nf && (F0 - cutoff) > f_mono[i_f_min] )
                {
                  ls[i_f_min] = 0;
                  ++i_f_min;
                }
              

              // Check whether we have elements in ls that have to be
              // set to zero at higher frequencies of f_mono.
              //
              // Loop through all frequencies, finding max value and
              // set all values to zero on that way.
              while ( i_f_max >= 0 && (F0 + cutoff) < f_mono[i_f_max] )
                {
                  ls[i_f_max] = 0;
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

        
          // Maybe there are no frequencies left to compute?  Note that
          // the number that counts here is nfl, since only these are
          // the `real' frequencies, for which xsec is changed. nfls
          // will always be at least one, because it contains the cutoff.
          if ( nfl > 0 )
            {
              // Calculate the line shape:
              lineshape_data[ind_ls].Function()(ls,aux,F0,gamma,sigma,
                                                f_local[Range(i_f_min,nfl)],
                                                nfls);

              // Calculate the chosen normalization factor:
              lineshape_norm_data[ind_lsn].Function()(fac,F0,
                                                      f_local[Range(i_f_min,nfl)],
                                                      nfls);

              // Get a handle on the range of xsec that we want to change.
              // We use nfl here, which could be one less than nfls.
              VectorView this_xsec = xsec(Range(i_f_min,nfl),i);

              // cutoff ?
              if ( cut )
                {
                  // The index nfls-1 should be exactly the index pointing
                  // to the value at the cutoff frequency.
                  base = ls[nfls-1] * fac[nfls-1];

                  // Subtract baseline from xsec. 
                  this_xsec -= base;
                }

              // Get handles on the range of ls and fac that we need.
              VectorView this_ls  = ls[Range(0,nfl)];
              VectorView this_fac = fac[Range(0,nfl)];

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

   \param   refr_index Output:  refractive index
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

   \param   refr_index Output:  refractive index
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
