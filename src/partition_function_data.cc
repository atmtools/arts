/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>,
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

/*!
  \file   partition_function_data.cc
  \brief  Coefficients for 3rd order polynomial of partition function
  in temperature, plus safety check. 

  This file contains the definition of this function and a safety check of
  the input. You have to add the coefficients of new species here, along
  with the new entry into the species_data (file species_data.cc)
  record, if you want to extend the capability of ARTS. These entries
  have to be in the same order, which is assured by check_q_data.

  \author Axel von Engeln
  \date 2000-08-21 */

#include "arts.h"
#include "make_array.h"
#include "absorption.h"



/*! \name Some #defines for better readability */
//@{ 
#define Qcoeff	     make_array<Numeric>		
//@} 


/**
  Define partition function coefficients lookup data.

  <h1>General Remarks</h1>
  <!------------------------>

  This function contains the polynomial coefficients of the partition
  function for each isotope of each species. The sorting of the array
  has to match the species_data entries, a safety check routine
  assures this order.

  The polynomial coefficients are mainly taken from the tips program,
  which comes along with the HITRAN database. tips separates the
  partition function into 3 temperature ranges, each range has its own
  coefficients. Ranges are: 70-500 K, 500-1500 K, 1500-3005 K. Since
  atmospheric temperatures are generally between 150 and 300 K, only
  the first range has been implemented.

  The coefficients for species which are not covered in HITRAN are
  calculated from JPL values. The JPL catalogue has a different way to
  calculate the partition function, it gives the partition function at
  specific temperatures: 300, 225, 150, 75, 37.5, 18.75, 9.375 K and
  an interpolation scheme is given for values inbetween. The partition
  functions are proportional to temperature T^1.5 for non-linear
  molecules (degrees of freedom: 3) and proportional to T for linear
  molecules (degrees of freedom 2). Nothing is known so far about the
  way how to treat atoms.

  The polynomial coeffients for JPL catalogue values are calculated
  with an idl program that comes along with the arts distribution
  (partition_function.pro). This calculates the partition function for
  all temperatures between 70 and 500 K from the given JPL values, by
  using the recommended JPL interpolation scheme. JPL does not
  consider vibrational energy levels in the partition functions and a
  correction based on a literature research is performed within the
  idl program (for energy levels refer to that idl file). The
  literature research yielded only main isotope vibrational levels,
  since nothing better is available they were used for the isotopes as
  well. The idl program fits a polynomial to the data, and outputs the
  coefficients. The polynomial is fitted between temperatures of 70 to
  500 K. No big difference is found in the ratio of the partition
  functions Q(300K)/Q(T) (the relevant quantity for the calculation of
  the JPL intensities to other temperatures) for other fit ranges,
  e.g., 9 - 300 K (range where the JPL partition functions are given),
  except for HCOOH, where the 0 to 300 K fit is better and was
  therefore use for this species (see below).

  This idl program allows additionally the comparison of partition
  functions of the 2 catalogues, when the species is defined in both
  catalogues. Partition function differ sometimes by more than 10 %,
  important nevertheless is the ratio of the partition function at the
  temperature where the catalogue line intensity is given, and the
  desired temperature, e.g., JPL: Q(300K)/Q/T. This ratio was found to
  vary as well sometimes by more than 10 % for specific species, as
  given in the table below, which list the difference found for the
  partition function and the one found in the ratio (Q(300K)/Q/T)):

\verbatim
  NAME     JPL TAG    HITRAN    Difference: Q / Ratio [%]   
  CH4      17003       63     > 10.0 /  > 10.0       
  NH3      18002      112     > 10.0 /  > 10.0       
  HNO3     63001      121     no coeff defined in tips
  H2O2     34004      251     > 10.0 /  > 10.0       
  COF2     66001      291     > 10.0 /  > 10.0       
  HCOOH    46005      321     > 10.0 /  > 10.0       
  O        16001      341     no coeff defined in tips
  ClONO2   97002      351     no coeff defined in tips
  ClONO2   99001      352     no coeff defined in tips
\endverbatim

  Certain species in the HITRAN database were marked in the tips
  program with a partition function equals -1 for all temperatures,
  meaning a polynomial fit is not applicable. These species are:

  HNO3, SF6, O, ClONO2, C2H6

  HNO3 and ClONO2 are covered in the JPL database, and a polynomial
  fit to the partition function calculated with the recommended scheme
  was found to be accurate within 4 percents, therefore the JPL
  partition function was used for these 2 molecules. 

  The other molecules (SF6, C2H6) and atom (O) were marked as not
  appropriate for a polynomial fit following the tips convention, e.g.,
  first coefficient is -1, all others are zero.

  HCOOH is the only species that shows errors of up to 7 % for the
  polynomial fit to the JPL partition functions, where the partition
  functions were fitted for the range 0 to 300 K (different to all
  other species). The fit with range 70 to 500 K showed errors larger
  than 10 %. It is to some extend questionable whether one should use
  a polynomial fit for this species, nevertheless it was done since
  the agreement between HITRAN and JPL ratios for the first isotope of
  HCOOH was found to be larger than 10 %, thus even larger than the
  fit error. It should be noted however, that the polynomial fit to
  the JPL HCOOH partition functions works very well for the first
  isotope, the other isotopes show a different partition function for
  some reason that I do not understand.


  <h1>Source for entries</h1>
  <!------------------------>

  <dl>
  <dt> Molecule Name:
  <dd> Arts convention

  <dt> Isotope Name:
  <dd> Arts convention

  <dt> Coefficients: 
  <dd> Generally taken from the tips program, version September 23,
  1997. Partition functions for species included only in JPL were
  calculated according to the JPL recommended scheme, and a 3rd order
  polynomial fit was performed (with the idl program mentioned above).
  Some molecules are included in HITRAN, but no coefficients are
  given, because a polynomial fit to the partition function is not
  applicable according to the tips program. Here, if applicable, the
  JPL partition function were used, and fitted with a polynomial.  The
  difference between the polynomial fit and the values calculated
  with the JPL scheme were below 2 % for the interesting atmospheric
  temperature range (except for HCOOH, see discussion above).

  A number code in the header gives the source of the polynomials:

  <table>
  <tr>
  <td> 1: <td> tips program, version September 23, 1997

  <tr>
  <td >2: <td> polynomial fit to JPL data

  </table>

  <dt> Quality:
  <dd> For each source of the coefficients, a different quality code
  is given, both qualities are given for the temperature range from
  150 to 300 K:

  <table>
  <tr>
  <td> 1: <td> for species covered in both catalogues, the maximum
  difference in the ratio between the two catalogues is given in
  percent (range 150 to 300 K).

  <tr>
  <td >2: <td> the difference in ratio between the polynomial fit to
  the JPL data and the original data (interpolation with the
  recommended scheme) (range 150 to 300 K).

  </table>
  </dl>

  \author Axel von Engeln  
  \date 2000-08-21 */

void spec(ARRAY<SpeciesRecord>::iterator& is,
	  ARRAY<IsotopeRecord>::iterator& ii,
	  string name);

void iso(ARRAY<IsotopeRecord>::iterator& ii,
	 string name,
	 const ARRAY<Numeric>& coeff);


void define_partition_species_data()
{
  extern ARRAY<SpeciesRecord> species_data;

  ARRAY<SpeciesRecord>::iterator it_species = species_data.begin();
  ARRAY<IsotopeRecord>::iterator it_isotope;


  // H2O
  // Coeff:       1      1      1      1      2      2
  // Quality:    0.28   0.28   0.35   0.46   0.32   0.34
  spec(it_species, it_isotope, "H2O");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"161",	Qcoeff(	-4.4405E+00	,2.7678E-01	,1.2536E-03	,-4.8938E-07) );
  iso(it_isotope,	"181",	Qcoeff(	-4.3624E+00	,2.7647E-01	,1.2802E-03	,-5.2046E-07) );
  iso(it_isotope,	"171",	Qcoeff(	-2.5767E+01	,1.6458E+00	,7.6905E-03	,-3.1668E-06) );
  iso(it_isotope,	"162",	Qcoeff(	-2.3916E+01	,1.3793E+00	,6.1246E-03	,-2.1530E-06) );
  iso(it_isotope,	"182",	Qcoeff(	-5.1056E+00	,2.4408E-01	,1.0230E-03	,-4.2596E-07) );
  iso(it_isotope,	"262",	Qcoeff(	-3.6689E+01	,1.7119E+00	,7.2123E-03	,-3.0200E-06) );



  // CO2
  // Coeff:       1      1      1      1      1      1      1      1
  // Quality:    ----   ----   4.40   4.32   ----   ----   ----   ----
  spec(it_species, it_isotope, "CO2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"626",	Qcoeff(	-1.3617E+00	,9.4899E-01	,-6.9259E-04	,2.5974E-06) );
  iso(it_isotope,	"636",	Qcoeff(	-2.0631E+00	,1.8873E+00	,-1.3669E-03	,5.4032E-06) );
  iso(it_isotope,	"628",	Qcoeff(	-2.9175E+00	,2.0114E+00	,-1.4786E-03	,5.5941E-06) );
  iso(it_isotope,	"627",	Qcoeff(	-1.6558E+01	,1.1733E+01	,-8.5844E-03	,3.2379E-05) );
  iso(it_isotope,	"638",	Qcoeff(	-4.4685E+00	,4.0330E+00	,-2.9590E-03	,1.1770E-05) );
  iso(it_isotope,	"637",	Qcoeff(	-2.6263E+01	,2.3350E+01	,-1.7032E-02	,6.7532E-05) );
  iso(it_isotope,	"828",	Qcoeff(	-1.4811E+00	,1.0667E+00	,-7.8758E-04	,3.0133E-06) );
  iso(it_isotope,	"728",	Qcoeff(	-1.7600E+01	,1.2445E+01	,-9.1837E-03	,3.4915E-05) );



  // O3
  // Coeff:       1      1      1      1      1
  // Quality:    3.86   2.64   2.49   1.25   1.28
  spec(it_species, it_isotope, "O3");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"666",	Qcoeff(	-1.6443E+02	,6.9047E+00	,1.0396E-02	,2.6669E-05) );
  iso(it_isotope,	"668",	Qcoeff(	-3.5222E+02	,1.4796E+01	,2.1475E-02	,5.9891E-05) );
  iso(it_isotope,	"686",	Qcoeff(	-1.7466E+02	,7.2912E+00	,1.0093E-02	,2.9991E-05) );
  iso(it_isotope,	"667",	Qcoeff(	-2.0540E+03	,8.5998E+01	,1.2667E-01	,3.3026E-04) );
  iso(it_isotope,	"676",	Qcoeff(	-1.0148E+03	,4.2494E+01	,6.2586E-02	,1.6319E-04) );



  // N2O
  // Coeff:       1      1      1      1      1
  // Quality:    0.89   1.26   0.95   0.82   ----
  spec(it_species, it_isotope, "N2O");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"446",	Qcoeff(	2.4892E+01	,1.4979E+01	,-7.6213E-03	,4.6310E-05) );
  iso(it_isotope,	"456",	Qcoeff(	3.6318E+01	,9.5497E+00	,-2.3943E-03	,2.6842E-05) );
  iso(it_isotope,	"546",	Qcoeff(	2.4241E+01	,1.0179E+01	,-4.3002E-03	,3.0425E-05) );
  iso(it_isotope,	"448",	Qcoeff(	6.7708E+01	,1.4878E+01	,-1.0730E-03	,3.4254E-05) );
  iso(it_isotope,	"447",	Qcoeff(	5.0069E+02	,8.4526E+01	,8.3494E-03	,1.7154E-04) );



  // CO
  // Coeff:       1      1      1      1      1      1
  // Quality:    0.03   0.02   0.02   0.02   ----   ----
  spec(it_species, it_isotope, "CO");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"26",	Qcoeff(	2.7758E-01	,3.6290E-01	,-7.4669E-06	,1.4896E-08) );
  iso(it_isotope,	"36",	Qcoeff(	5.3142E-01	,7.5953E-01	,-1.7810E-05	,3.5160E-08) );
  iso(it_isotope,	"28",	Qcoeff(	2.6593E-01	,3.8126E-01	,-9.2083E-06	,1.8086E-08) );
  iso(it_isotope,	"27",	Qcoeff(	1.6376E+00	,2.2343E+00	,-4.9025E-05	,9.7389E-08) );
  iso(it_isotope,	"38",	Qcoeff(	5.1216E-01	,7.9978E-01	,-2.1784E-05	,4.2749E-08) );
  iso(it_isotope,	"37",	Qcoeff(	3.2731E+00	,4.6577E+00	,-6.9833E-05	,1.8853E-07) );



  // CH4
  // Coeff:       1      1      1
  // Quality:    ----   ----  21.43
  spec(it_species, it_isotope, "CH4");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"211",	Qcoeff(	-2.6479E+01	,1.1557E+00	,2.6831E-03	,1.5117E-06) );
  iso(it_isotope,	"311",	Qcoeff(	-5.2956E+01	,2.3113E+00	,5.3659E-03	,3.0232E-06) );
  iso(it_isotope,	"212",	Qcoeff(	-2.1577E+02	,9.3318E+00	,2.1779E-02	,1.2183E-05) );



  // O2
  // Coeff:       1      1      1
  // Quality:    0.11   0.62   0.65
  spec(it_species, it_isotope, "O2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"66",	Qcoeff(	3.5923E-01	,7.3534E-01	,-6.4870E-05	,1.3073E-07) );
  iso(it_isotope,	"68",	Qcoeff(	-4.0039E+00	,1.5595E+00	,-1.5357E-04	,3.0969E-07) );
  iso(it_isotope,	"67",	Qcoeff(	-2.3325E+01	,9.0981E+00	,-8.4435E-04	,1.7062E-06) );



  // NO
  // Coeff:       1      1      1
  // Quality:    1.56   ----   ----
  spec(it_species, it_isotope, "NO");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"46",	Qcoeff(	-2.5296E+01	,2.6349E+00	,5.8517E-03	,-5.2020E-06) );
  iso(it_isotope,	"56",	Qcoeff(	-1.4990E+01	,1.8240E+00	,4.0261E-03	,-3.5648E-06) );
  iso(it_isotope,	"48",	Qcoeff(	-2.6853E+01	,2.7816E+00	,6.1493E-03	,-5.4410E-06) );



  // SO2
  // Coeff:       1      1      2      2
  // Quality:    1.26   1.23   1.32   1.32
  spec(it_species, it_isotope, "SO2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"626",	Qcoeff(	-2.4056E+02	,1.1101E+01	,2.2164E-02	,5.2334E-05) );
  iso(it_isotope,	"646",	Qcoeff(	-2.4167E+02	,1.1151E+01	,2.2270E-02	,5.2550E-05) );
  iso(it_isotope,	"636",	Qcoeff(	5.8740E+02	,1.2472E+01	,2.9113E-01	,-1.6236E-04) );
  iso(it_isotope,	"628",	Qcoeff(	3.1299E+02	,6.6372E+00	,1.5485E-01	,-8.6343E-05) );



  // NO2
  // Coeff:       1
  // Quality:    0.98
  spec(it_species, it_isotope, "NO2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"646",	Qcoeff(	-5.3042E+02	,2.4216E+01	,6.6856E-02	,4.3823E-05) );



  // NH3
  // Coeff:       1      1      2
  // Quality:    1.09  21.12   0.78
  spec(it_species, it_isotope, "NH3");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"4111",	Qcoeff(	-6.2293E+01	,3.0915E+00	,9.4575E-03	,1.8416E-06) );
  iso(it_isotope,	"5111",	Qcoeff(	-4.2130E+01	,2.0569E+00	,6.3387E-03	,1.2127E-06) );
  iso(it_isotope,	"4112",	Qcoeff(	-6.5172E+01	,4.3080E+00	,3.3865E-02	,-1.6513E-05) );



  // HNO3
  // Coeff:       2
  // Quality:    4.00
  spec(it_species, it_isotope, "HNO3");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"146",	Qcoeff(	4.7202E+03	,-5.5810E+01	,6.6216E-01	,-4.2602E-04) );



  // OH
  // Coeff:       1      1      1
  // Quality:    0.85   0.84   1.04
  spec(it_species, it_isotope, "OH");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"61",	Qcoeff(	8.7390E+00	,1.5977E-01	,3.8291E-04	,-3.5669E-07) );
  iso(it_isotope,	"81",	Qcoeff(	8.6770E+00	,1.6175E-01	,3.8223E-04	,-3.5466E-07) );
  iso(it_isotope,	"62",	Qcoeff(	1.0239E+01	,4.3783E-01	,1.0477E-03	,-9.4570E-07) );



  // HF
  // Coeff:       1      2
  // Quality:    0.02   0.04
  spec(it_species, it_isotope, "HF");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"19",	Qcoeff(	1.5486E+00	,1.3350E-01	,5.9154E-06	,-4.6889E-09) );
  iso(it_isotope,	"29",	Qcoeff(	3.4435E-01	,6.3996E-02	,4.2331E-07	,-3.0344E-10) );



  // HCl
  // Coeff:       1      1      2      2
  // Quality:    2.32   2.30   0.03   0.03
  spec(it_species, it_isotope, "HCl");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"15",	Qcoeff(	2.8627E+00	,5.3122E-01	,6.7464E-06	,-1.6730E-09) );
  iso(it_isotope,	"17",	Qcoeff(	2.8617E+00	,5.3203E-01	,6.6553E-06	,-1.5168E-09) );
  iso(it_isotope,	"25",	Qcoeff(	1.2492E+00	,5.1693E-01	,-1.0216E-06	,2.4388E-09) );
  iso(it_isotope,	"27",	Qcoeff(	1.2054E+00	,5.1906E-01	,-3.1578E-06	,4.8085E-09) );



  // HBr
  // Coeff:       1      1
  // Quality:    1.89   1.89
  spec(it_species, it_isotope, "HBr");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"19",	Qcoeff(	2.7963E+00	,6.6532E-01	,3.4255E-06	,5.2274E-09) );
  iso(it_isotope,	"11",	Qcoeff(	2.7953E+00	,6.6554E-01	,3.2931E-06	,5.4823E-09) );



  // HI
  // Coeff:       1
  // Quality:    ----
  spec(it_species, it_isotope, "HI");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"17",	Qcoeff(	4.0170E+00	,1.3003E+00	,-1.1409E-05	,4.0026E-08) );



  // ClO
  // Coeff:       1      1
  // Quality:   19.45  19.46
  spec(it_species, it_isotope, "ClO");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"56",	Qcoeff(	9.0968E+01	,7.0918E+00	,1.1639E-02	,3.0145E-06) );
  iso(it_isotope,	"76",	Qcoeff(	9.2598E+01	,7.2085E+00	,1.1848E-02	,3.1305E-06) );



  // OCS
  // Coeff:       1      1      1      1
  // Quality:    1.49   1.66   2.43   2.15
  spec(it_species, it_isotope, "OCS");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"622",	Qcoeff(	-9.3697E-01	,3.6090E+00	,-3.4552E-03	,1.7462E-05) );
  iso(it_isotope,	"624",	Qcoeff(	-1.1536E+00	,3.7028E+00	,-3.5582E-03	,1.7922E-05) );
  iso(it_isotope,	"632",	Qcoeff(	-6.1015E-01	,7.2200E+00	,-7.0044E-03	,3.6708E-05) );
  iso(it_isotope,	"822",	Qcoeff(	-2.1569E-01	,3.8332E+00	,-3.6783E-03	,1.9177E-05) );



  // H2CO
  // Coeff:       1      1      1      2      2
  // Quality:    2.80   3.18   0.74   0.18   0.21
  spec(it_species, it_isotope, "H2CO");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"1126",	Qcoeff(	-1.1760E+02	,4.6885E+00	,1.5088E-02	,3.5367E-06) );
  iso(it_isotope,	"1136",	Qcoeff(	-2.4126E+02	,9.6134E+00	,3.0938E-02	,7.2579E-06) );
  iso(it_isotope,	"1128",	Qcoeff(	-1.1999E+02	,5.2912E+00	,1.4686E-02	,4.3505E-06) );
  iso(it_isotope,	"1226",	Qcoeff(	-6.3717E+01	,3.0566E+00	,1.3879E-02	,-5.9420E-06) );
  iso(it_isotope,	"2226",	Qcoeff(	-3.4485E+02	,1.7375E+01	,8.4988E-02	,-3.7365E-05) );



  // HOCl
  // Coeff:       1      1
  // Quality:    0.95   0.95
  spec(it_species, it_isotope, "HOCl");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"165",	Qcoeff(	-7.3640E+02	,3.4149E+01	,9.3554E-02	,6.7409E-05) );
  iso(it_isotope,	"167",	Qcoeff(	-7.4923E+02	,3.4747E+01	,9.5251E-02	,6.8523E-05) );



  // N2
  // Coeff:       1
  // Quality:    ----
  spec(it_species, it_isotope, "N2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"44",	Qcoeff(	1.3684E+00	,1.5756E+00	,-1.8511E-05	,3.8960E-08) );



  // HCN
  // Coeff:       1      1      1      2
  // Quality:    0.51   0.38   0.45   1.72
  spec(it_species, it_isotope, "HCN");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"124",	Qcoeff(	-1.3992E+00	,2.9619E+00	,-1.7464E-03	,6.5937E-06) );
  iso(it_isotope,	"134",	Qcoeff(	-2.5869E+00	,6.0744E+00	,-3.5719E-03	,1.3654E-05) );
  iso(it_isotope,	"125",	Qcoeff(	-1.1408E+00	,2.0353E+00	,-1.2159E-03	,4.6375E-06) );
  iso(it_isotope,	"224",	Qcoeff(	1.4949E+01	,1.4567E+00	,9.0970E-04	,8.0655E-07) );



  // CH3Cl
  // Coeff:       1      1
  // Quality:    1.65   1.71
  spec(it_species, it_isotope, "CH3Cl");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"215",	Qcoeff(	-9.1416E+02	,3.4081E+01	,7.5461E-03	,1.7933E-04) );
  iso(it_isotope,	"217",	Qcoeff(	-9.2868E+02	,3.4621E+01	,7.6674E-03	,1.8217E-04) );



  // H2O2
  // Coeff:       1
  // Quality:   14.46
  spec(it_species, it_isotope, "H2O2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"1661",	Qcoeff(	-3.6499E+02	,1.3712E+01	,3.8658E-02	,2.3052E-05) );



  // C2H2
  // Coeff:       1      1
  // Quality:    ----   ----
  spec(it_species, it_isotope, "C2H2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"1221",	Qcoeff(	-8.3088E+00	,1.4484E+00	,-2.5946E-03	,8.4612E-06) );
  iso(it_isotope,	"1231",	Qcoeff(	-6.6736E+01	,1.1592E+01	,-2.0779E-02	,6.7719E-05) );



  // C2H6
  // Coeff:       1
  // Quality:    ----
  spec(it_species, it_isotope, "C2H6");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"1221",	Qcoeff(	-1.0000E+00	,0.0000E+00	,0.0000E+00	,0.0000E+00) );



  // PH3
  // Coeff:       1
  // Quality:    2.06
  spec(it_species, it_isotope, "PH3");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"1111",	Qcoeff(	-1.5068E+02	,6.4718E+00	,1.2588E-02	,1.4759E-05) );



  // COF2
  // Coeff:       1
  // Quality:   13.92
  spec(it_species, it_isotope, "COF2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"269",	Qcoeff(	-5.4180E+03	,1.8868E+02	,-3.3139E-01	,1.8650E-03) );



  // SF6
  // Coeff:       1
  // Quality:    ----
  spec(it_species, it_isotope, "SF6");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"29",	Qcoeff(	-1.0000E+00	,0.0000E+00	,0.0000E+00	,0.0000E+00) );



  // H2S
  // Coeff:       1      1      1      2
  // Quality:    0.29   ----   ----   0.32
  spec(it_species, it_isotope, "H2S");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"121",	Qcoeff(	-1.5521E+01	,8.3130E-01	,3.3656E-03	,-8.5691E-07) );
  iso(it_isotope,	"141",	Qcoeff(	-1.5561E+01	,8.3337E-01	,3.3744E-03	,-8.5937E-07) );
  iso(it_isotope,	"131",	Qcoeff(	-6.2170E+01	,3.3295E+00	,1.3480E-02	,-3.4323E-06) );
  iso(it_isotope,	"122",	Qcoeff(	-1.5731E+01	,7.1204E-01	,2.9745E-03	,-1.2403E-06) );



  // HCOOH
  // Coeff:       1      2      2      2
  // Quality:   12.54   6.54   0.55   0.57
  spec(it_species, it_isotope, "HCOOH");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"1261",	Qcoeff(	-2.9550E+03	,1.0349E+02	,-1.3146E-01	,8.7787E-04) );
  iso(it_isotope,	"1361",	Qcoeff(	2.0082E+02	,-4.9620E+00	,2.3524E-01	,-4.0531E-04) );
  iso(it_isotope,	"2261",	Qcoeff(	-2.5614E+01	,8.4704E+00	,1.2417E-01	,-1.1758E-04) );
  iso(it_isotope,	"1262",	Qcoeff(	-5.7615E+01	,9.5919E+00	,1.0930E-01	,-1.0134E-04) );



  // HO2
  // Coeff:       1
  // Quality:    1.10
  spec(it_species, it_isotope, "HO2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"166",	Qcoeff(	-1.5684E+02	,7.4450E+00	,2.6011E-02	,-9.2704E-07) );



  // O
  // Coeff:       1
  // Quality:    ----
  spec(it_species, it_isotope, "O");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"6",	Qcoeff(	-1.0000E+00	,0.0000E+00	,0.0000E+00	,0.0000E+00) );



  // ClONO2
  // Coeff:       2      2
  // Quality:    1.87   1.85
  spec(it_species, it_isotope, "ClONO2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"5646",	Qcoeff(	-4.9350E+03	,2.7375E+02	,7.2257E+00	,-4.2583E-03) );
  iso(it_isotope,	"7646",	Qcoeff(	-4.9721E+03	,2.7978E+02	,7.4127E+00	,-4.3703E-03) );



  // NO+
  // Coeff:       1
  // Quality:    0.01
  spec(it_species, it_isotope, "NO+");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"46",	Qcoeff(	9.1798E-01	,1.0416E+00	,-1.1614E-05	,2.4499E-08) );



  // OClO
  // Coeff:       2      2
  // Quality:    1.09   1.13
  spec(it_species, it_isotope, "OClO");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"656",	Qcoeff(	1.0809E+02	,4.5813E+01	,6.1875E-01	,-3.4280E-04) );
  iso(it_isotope,	"676",	Qcoeff(	1.0495E+03	,3.1385E+01	,6.6403E-01	,-3.7240E-04) );



  // BrO
  // Coeff:       2      2
  // Quality:    0.09   0.09
  spec(it_species, it_isotope, "BrO");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"96",	Qcoeff(	-8.8945E+00	,1.3384E+01	,-1.6496E-03	,1.4828E-06) );
  iso(it_isotope,	"16",	Qcoeff(	-8.4100E+00	,1.3433E+01	,-1.6322E-03	,1.4620E-06) );



  // H2SO4
  // Coeff:       2
  // Quality:    0.39
  spec(it_species, it_isotope, "H2SO4");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"126",	Qcoeff(	-6.6207E+03	,2.6856E+02	,1.0187E+00	,-4.0717E-04) );



  // Cl2O2
  // Coeff:       2      2
  // Quality:    0.29   0.30
  spec(it_species, it_isotope, "Cl2O2");
  //			Name		c0		c1		c2		c3
  //			|		|		|		|		|
  iso(it_isotope,	"565",	Qcoeff(	-1.6546E+04	,7.0845E+02	,2.9249E+00	,-1.2097E-03) );
  iso(it_isotope,	"765",	Qcoeff(	-1.6893E+04	,7.2651E+02	,3.0185E+00	,-1.2520E-03) );


}


void spec(ARRAY<SpeciesRecord>::iterator& is,
	  ARRAY<IsotopeRecord>::iterator& ii,
	  string name)
{
  
#ifndef NDEBUG
  {
    assert( name == is->Name() );
  }
#endif

  ii = is->Isotope().begin();

  is++;

}


void iso(ARRAY<IsotopeRecord>::iterator& ii,
	 string name,
	 const ARRAY<Numeric>& coeff)
{
#ifndef NDEBUG
  {
    assert( name == ii->Name() );
  }
#endif

  ii->SetPartitionFctCoeff(coeff);

  ii++;

}

