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


/*! The lookup information for all the different species. */
ARRAY<QRecord> q_data;


/*! \name Some #defines for better readability */
//@{ 
#define NAME(x)      x				
#define ISOTOPE      make_array<QIsotopeRecord>	
#define REC	     QIsotopeRecord		
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
  using the recommended JPL interpolation scheme. The program then
  fits a polynomial to the data, and outputs the coefficients. The
  polynomial is fitted between temperatures of 70 to 500 K. No big
  difference is found in the ratio of the partition functions
  Q(300K)/Q(T) (the relevant quantity for the calculation of
  the JPL intensities to other temperatures) for other fit ranges,
  e.g., 9 - 300 K (range where the JPL partition functions are given),
  except for HCOOH, where the  0 to 300 K fit is better and was
  therefore use for this species (see below).

  This idl program allows additionally the comparison of partition
  functions of the 2 catalogues, when the species is defined in both
  catalogues. Partition function differ often by more than 10 %,
  important nevertheless is the ratio of the partition function at the
  temperature where the catalogue line intensity is given, and the
  desired temperature, e.g., JPL: Q(300K)/Q/T. This ratio was found
  to vary as well by more than 10 % for specific species, as given in
  the table below, which list the difference found for the partition
  function and the one found in the ratio (Q(300K)/Q/T)):

\verbatim
  NAME     JPL TAG    HITRAN    Difference: Q / Ratio [%]   
  CO2      46013       23     > 10.0 /  10.0         
  CO2      45012       24     > 10.0 /  10.0         
  N2O      44004       41     > 10.0 /  > 10.0       
  N2O      45007       42     > 10.0 /  > 10.0       
  N2O      45008       43     > 10.0 /  > 10.0       
  N2O      46007       44     > 10.0 /  > 10.0       
  CH4      17003       63     > 10.0 /  > 10.0       
  SO2      64002       91     > 10.0 / 10.0          
  SO2      66002       92     > 10.0 / 10.0          
  NH3      18002      112     > 10.0 /  > 10.0       
  HNO3     63001      121     no coeff defined in tips
  OCS      60001      191     > 10.0 /  > 10.0       
  OCS      62001      192     > 10.0 /  > 10.0       
  OCS      61001      193     > 10.0 /  > 10.0       
  OCS      62002      194     > 10.0 /  > 10.0       
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
  was found to be accurate within 2 percents, therefore the JPL
  partition function was used for these 2 molecules. The other
  molecules (SF6, C2H6) and atom (O) were marked as not appropriate for
  a polynomial fit following the tips convention, e.g., first
  coefficient is -1, all others are zero.

  HCOOH is the only species that shows errors of up to 7 % for the
  polynomial fit to the JPL partition functions, where the partition
  functions were fitted for the range 0 to 300 K (different to all
  other species). The fit with range 70 to 500 K showed errors larger
  than 10 %. It is to some extend questionable whether one should use
  a polynomial fit for this species, nevertheless it was done since
  the agreement between HITRAN and JPL ratios for the first isotope of
  HCOOH was found to be larger than 10 %, thus even larger than the
  fit error.


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
  \date   2000-08-21 
*/
void define_q_data()
{
  extern ARRAY<QRecord> q_data;

  // Initialize to zero, just in case:
  q_data.clear();

  /* Here's an empty template record entry:

  // species
  // Coeff:
  // Quality:
  q_data.push_back
    ( QRecord
      ( NAME(""),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"",	Qcoeff(			,		,		,) ),
	 REC(	"",	Qcoeff(			,		,		,) )
	 ) ) );

  */


  // H2O
  // Coeff:       1      1      1      1      2      2
  // Quality:    0.33   0.33   0.39   0.50   0.33   0.35
  q_data.push_back
    ( QRecord
      ( NAME("H2O"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"161",	Qcoeff(	-4.4405E+00	,2.7678E-01	,1.2536E-03	,-4.8938E-07) ),
	 REC(	"181",	Qcoeff(	-4.3624E+00	,2.7647E-01	,1.2802E-03	,-5.2046E-07) ),
	 REC(	"171",	Qcoeff(	-2.5767E+01	,1.6458E+00	,7.6905E-03	,-3.1668E-06) ),
	 REC(	"162",	Qcoeff(	-2.3916E+01	,1.3793E+00	,6.1246E-03	,-2.1530E-06) ),
	 REC(	"182",	Qcoeff(	-5.1419E+00	,2.4479E-01	,1.0199E-03	,-4.2468E-07) ),
	 REC(	"262",	Qcoeff(	-3.7362E+01	,1.7225E+00	,7.1688E-03	,-2.9865E-06) )
	 ) ) );



  // CO2
  // Coeff:       1      1      1      1      1      1      1      1
  // Quality:    ----   ----   8.76   8.67   ----   ----   ----   ----
  q_data.push_back
    ( QRecord
      ( NAME("CO2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"626",	Qcoeff(	-1.3617E+00	,9.4899E-01	,-6.9259E-04	,2.5974E-06) ),
	 REC(	"636",	Qcoeff(	-2.0631E+00	,1.8873E+00	,-1.3669E-03	,5.4032E-06) ),
	 REC(	"628",	Qcoeff(	-2.9175E+00	,2.0114E+00	,-1.4786E-03	,5.5941E-06) ),
	 REC(	"627",	Qcoeff(	-1.6558E+01	,1.1733E+01	,-8.5844E-03	,3.2379E-05) ),
	 REC(	"638",	Qcoeff(	-4.4685E+00	,4.0330E+00	,-2.9590E-03	,1.1770E-05) ),
	 REC(	"637",	Qcoeff(	-2.6263E+01	,2.3350E+01	,-1.7032E-02	,6.7532E-05) ),
	 REC(	"828",	Qcoeff(	-1.4811E+00	,1.0667E+00	,-7.8758E-04	,3.0133E-06) ),
	 REC(	"728",	Qcoeff(	-1.7600E+01	,1.2445E+01	,-9.1837E-03	,3.4915E-05) )
	 ) ) );



  // O3
  // Coeff:       1      1      1      1      1
  // Quality:    1.30   5.82   5.88   5.25   5.29
  q_data.push_back
    ( QRecord
      ( NAME("O3"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"666",	Qcoeff(	-1.6443E+02	,6.9047E+00	,1.0396E-02	,2.6669E-05) ),
	 REC(	"668",	Qcoeff(	-3.5222E+02	,1.4796E+01	,2.1475E-02	,5.9891E-05) ),
	 REC(	"686",	Qcoeff(	-1.7466E+02	,7.2912E+00	,1.0093E-02	,2.9991E-05) ),
	 REC(	"667",	Qcoeff(	-2.0540E+03	,8.5998E+01	,1.2667E-01	,3.3026E-04) ),
	 REC(	"676",	Qcoeff(	-1.0148E+03	,4.2494E+01	,6.2586E-02	,1.6319E-04) )
	 ) ) );



  // N2O
  // Coeff:       1      1      1      1      1
  // Quality:   12.71  13.33  12.83  12.11   ----
  q_data.push_back
    ( QRecord
      ( NAME("N2O"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"446",	Qcoeff(	2.4892E+01	,1.4979E+01	,-7.6213E-03	,4.6310E-05) ),
	 REC(	"456",	Qcoeff(	3.6318E+01	,9.5497E+00	,-2.3943E-03	,2.6842E-05) ),
	 REC(	"546",	Qcoeff(	2.4241E+01	,1.0179E+01	,-4.3002E-03	,3.0425E-05) ),
	 REC(	"448",	Qcoeff(	6.7708E+01	,1.4878E+01	,-1.0730E-03	,3.4254E-05) ),
	 REC(	"447",	Qcoeff(	5.0069E+02	,8.4526E+01	,8.3494E-03	,1.7154E-04) )
	 ) ) );



  // CO
  // Coeff:       1      1      1      1      1      1
  // Quality:    0.03   0.01   0.01   0.02   ----   ----
  q_data.push_back
    ( QRecord
      ( NAME("CO"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"26",	Qcoeff(	2.7758E-01	,3.6290E-01	,-7.4669E-06	,1.4896E-08) ),
	 REC(	"36",	Qcoeff(	5.3142E-01	,7.5953E-01	,-1.7810E-05	,3.5160E-08) ),
	 REC(	"28",	Qcoeff(	2.6593E-01	,3.8126E-01	,-9.2083E-06	,1.8086E-08) ),
	 REC(	"27",	Qcoeff(	1.6376E+00	,2.2343E+00	,-4.9025E-05	,9.7389E-08) ),
	 REC(	"38",	Qcoeff(	5.1216E-01	,7.9978E-01	,-2.1784E-05	,4.2749E-08) ),
	 REC(	"37",	Qcoeff(	3.2731E+00	,4.6577E+00	,-6.9833E-05	,1.8853E-07) )
	 ) ) );



  // CH4
  // Coeff:       1      1      1
  // Quality:    ----   ----  22.18
  q_data.push_back
    ( QRecord
      ( NAME("CH4"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"211",	Qcoeff(	-2.6479E+01	,1.1557E+00	,2.6831E-03	,1.5117E-06) ),
	 REC(	"311",	Qcoeff(	-5.2956E+01	,2.3113E+00	,5.3659E-03	,3.0232E-06) ),
	 REC(	"212",	Qcoeff(	-2.1577E+02	,9.3318E+00	,2.1779E-02	,1.2183E-05) )
	 ) ) );



  // O2
  // Coeff:       1      1      1
  // Quality:    0.06   0.68   0.71
  q_data.push_back
    ( QRecord
      ( NAME("O2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"66",	Qcoeff(	3.5923E-01	,7.3534E-01	,-6.4870E-05	,1.3073E-07) ),
	 REC(	"68",	Qcoeff(	-4.0039E+00	,1.5595E+00	,-1.5357E-04	,3.0969E-07) ),
	 REC(	"67",	Qcoeff(	-2.3325E+01	,9.0981E+00	,-8.4435E-04	,1.7062E-06) )
	 ) ) );



  // NO
  // Coeff:       1      1      1
  // Quality:    1.56   ----   ----
  q_data.push_back
    ( QRecord
      ( NAME("NO"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"46",	Qcoeff(	-2.5296E+01	,2.6349E+00	,5.8517E-03	,-5.2020E-06) ),
	 REC(	"56",	Qcoeff(	-1.4990E+01	,1.8240E+00	,4.0261E-03	,-3.5648E-06) ),
	 REC(	"48",	Qcoeff(	-2.6853E+01	,2.7816E+00	,6.1493E-03	,-5.4410E-06) )
	 ) ) );



  // SO2
  // Coeff:       1      1      2      2
  // Quality:    9.28   9.25   0.35   0.35
  q_data.push_back
    ( QRecord
      ( NAME("SO2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"626",	Qcoeff(	-2.4056E+02	,1.1101E+01	,2.2164E-02	,5.2334E-05) ),
	 REC(	"646",	Qcoeff(	-2.4167E+02	,1.1151E+01	,2.2270E-02	,5.2550E-05) ),
	 REC(	"636",	Qcoeff(	-8.9492E+02	,3.8997E+01	,1.6619E-01	,-7.0199E-05) ),
	 REC(	"628",	Qcoeff(	-4.8125E+02	,2.0825E+01	,8.8091E-02	,-3.6976E-05) )
	 ) ) );



  // NO2
  // Coeff:       1
  // Quality:    3.46
  q_data.push_back
    ( QRecord
      ( NAME("NO2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"646",	Qcoeff(	-5.3042E+02	,2.4216E+01	,6.6856E-02	,4.3823E-05) )
	 ) ) );



  // NH3
  // Coeff:       1      1      2
  // Quality:    2.07  22.22   0.52
  q_data.push_back
    ( QRecord
      ( NAME("NH3"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"4111",	Qcoeff(	-6.2293E+01	,3.0915E+00	,9.4575E-03	,1.8416E-06) ),
	 REC(	"5111",	Qcoeff(	-4.2130E+01	,2.0569E+00	,6.3387E-03	,1.2127E-06) ),
	 REC(	"4112",	Qcoeff(	-1.0529E+02	,4.9832E+00	,3.0955E-02	,-1.4409E-05) )
	 ) ) );



  // HNO3
  // Coeff:       2
  // Quality:    0.38
  q_data.push_back
    ( QRecord
      ( NAME("HNO3"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"146",	Qcoeff(	-1.0936E+03	,4.6423E+01	,1.9143E-01	,-7.9439E-05) )
	 ) ) );



  // OH
  // Coeff:       1      1      1
  // Quality:    0.85   0.84   1.04
  q_data.push_back
    ( QRecord
      ( NAME("OH"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"61",	Qcoeff(	8.7390E+00	,1.5977E-01	,3.8291E-04	,-3.5669E-07) ),
	 REC(	"81",	Qcoeff(	8.6770E+00	,1.6175E-01	,3.8223E-04	,-3.5466E-07) ),
	 REC(	"62",	Qcoeff(	1.0239E+01	,4.3783E-01	,1.0477E-03	,-9.4570E-07) )
	 ) ) );



  // HF
  // Coeff:       1      2
  // Quality:    0.02   0.04
  q_data.push_back
    ( QRecord
      ( NAME("HF"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"19",	Qcoeff(	1.5486E+00	,1.3350E-01	,5.9154E-06	,-4.6889E-09) ),
	 REC(	"29",	Qcoeff(	3.4435E-01	,6.3996E-02	,4.2331E-07	,-3.0344E-10) )
	 ) ) );



  // HCl
  // Coeff:       1      1      2      2
  // Quality:    2.32   2.30   0.03   0.03
  q_data.push_back
    ( QRecord
      ( NAME("HCl"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"15",	Qcoeff(	2.8627E+00	,5.3122E-01	,6.7464E-06	,-1.6730E-09) ),
	 REC(	"17",	Qcoeff(	2.8617E+00	,5.3203E-01	,6.6553E-06	,-1.5168E-09) ),
	 REC(	"25",	Qcoeff(	1.2492E+00	,5.1693E-01	,-1.0216E-06	,2.4388E-09) ),
	 REC(	"27",	Qcoeff(	1.2017E+00	,5.1911E-01	,-3.3197E-06	,4.9817E-09) )
	 ) ) );



  // HBr
  // Coeff:       1      1
  // Quality:    1.89   1.89
  q_data.push_back
    ( QRecord
      ( NAME("HBr"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"19",	Qcoeff(	2.7963E+00	,6.6532E-01	,3.4255E-06	,5.2274E-09) ),
	 REC(	"11",	Qcoeff(	2.7953E+00	,6.6554E-01	,3.2931E-06	,5.4823E-09) )
	 ) ) );



  // HI
  // Coeff:       1
  // Quality:    ----
  q_data.push_back
    ( QRecord
      ( NAME("HI"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"17",	Qcoeff(	4.0170E+00	,1.3003E+00	,-1.1409E-05	,4.0026E-08) )
	 ) ) );



  // ClO
  // Coeff:       1      1
  // Quality:    0.83   0.81
  q_data.push_back
    ( QRecord
      ( NAME("ClO"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"56",	Qcoeff(	9.0968E+01	,7.0918E+00	,1.1639E-02	,3.0145E-06) ),
	 REC(	"76",	Qcoeff(	9.2598E+01	,7.2085E+00	,1.1848E-02	,3.1305E-06) )
	 ) ) );



  // OCS
  // Coeff:       1      1      1      1
  // Quality:   19.14  19.14  20.48  20.07
  q_data.push_back
    ( QRecord
      ( NAME("OCS"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"622",	Qcoeff(	-9.3697E-01	,3.6090E+00	,-3.4552E-03	,1.7462E-05) ),
	 REC(	"624",	Qcoeff(	-1.1536E+00	,3.7028E+00	,-3.5582E-03	,1.7922E-05) ),
	 REC(	"632",	Qcoeff(	-6.1015E-01	,7.2200E+00	,-7.0044E-03	,3.6708E-05) ),
	 REC(	"822",	Qcoeff(	-2.1569E-01	,3.8332E+00	,-3.6783E-03	,1.9177E-05) )
	 ) ) );



  // H2CO
  // Coeff:       1      1      1      2      2
  // Quality:    3.55   3.94   1.39   0.34   0.34
  q_data.push_back
    ( QRecord
      ( NAME("H2CO"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"1126",	Qcoeff(	-1.1760E+02	,4.6885E+00	,1.5088E-02	,3.5367E-06) ),
	 REC(	"1136",	Qcoeff(	-2.4126E+02	,9.6134E+00	,3.0938E-02	,7.2579E-06) ),
	 REC(	"1128",	Qcoeff(	-1.1999E+02	,5.2912E+00	,1.4686E-02	,4.3505E-06) ),
	 REC(	"1226",	Qcoeff(	-7.8182E+01	,3.2960E+00	,1.2870E-02	,-5.2158E-06) ),
	 REC(	"2226",	Qcoeff(	-4.2971E+02	,1.8783E+01	,7.9045E-02	,-3.3104E-05) )
	 ) ) );



  // HOCl
  // Coeff:       1      1
  // Quality:    3.89   3.89
  q_data.push_back
    ( QRecord
      ( NAME("HOCl"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"165",	Qcoeff(	-7.3640E+02	,3.4149E+01	,9.3554E-02	,6.7409E-05) ),
	 REC(	"167",	Qcoeff(	-7.4923E+02	,3.4747E+01	,9.5251E-02	,6.8523E-05) )
	 ) ) );



  // N2
  // Coeff:       1
  // Quality:    ----
  q_data.push_back
    ( QRecord
      ( NAME("N2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"44",	Qcoeff(	1.3684E+00	,1.5756E+00	,-1.8511E-05	,3.8960E-08) )
	 ) ) );



  // HCN
  // Coeff:       1      1      1      2
  // Quality:    6.84   7.05   7.13   1.42
  q_data.push_back
    ( QRecord
      ( NAME("HCN"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"124",	Qcoeff(	-1.3992E+00	,2.9619E+00	,-1.7464E-03	,6.5937E-06) ),
	 REC(	"134",	Qcoeff(	-2.5869E+00	,6.0744E+00	,-3.5719E-03	,1.3654E-05) ),
	 REC(	"125",	Qcoeff(	-1.1408E+00	,2.0353E+00	,-1.2159E-03	,4.6375E-06) ),
	 REC(	"224",	Qcoeff(	-1.0261E+01	,1.9009E+00	,-1.1334E-03	,2.3164E-06) )
	 ) ) );



  // CH3Cl
  // Coeff:       1      1
  // Quality:    5.86   5.92
  q_data.push_back
    ( QRecord
      ( NAME("CH3Cl"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"215",	Qcoeff(	-9.1416E+02	,3.4081E+01	,7.5461E-03	,1.7933E-04) ),
	 REC(	"217",	Qcoeff(	-9.2868E+02	,3.4621E+01	,7.6674E-03	,1.8217E-04) )
	 ) ) );



  // H2O2
  // Coeff:       1
  // Quality:   14.46
  q_data.push_back
    ( QRecord
      ( NAME("H2O2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"1661",	Qcoeff(	-3.6499E+02	,1.3712E+01	,3.8658E-02	,2.3052E-05) )
	 ) ) );



  // C2H2
  // Coeff:       1      1
  // Quality:    ----   ----
  q_data.push_back
    ( QRecord
      ( NAME("C2H2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"1221",	Qcoeff(	-8.3088E+00	,1.4484E+00	,-2.5946E-03	,8.4612E-06) ),
	 REC(	"1231",	Qcoeff(	-6.6736E+01	,1.1592E+01	,-2.0779E-02	,6.7719E-05) )
	 ) ) );



  // C2H6
  // Coeff:       1
  // Quality:    ----
  q_data.push_back
    ( QRecord
      ( NAME("C2H6"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"1221",	Qcoeff(	-1.0000E+00	,0.0000E+00	,0.0000E+00	,0.0000E+00) )
	 ) ) );



  // PH3
  // Coeff:       1
  // Quality:    3.70
  q_data.push_back
    ( QRecord
      ( NAME("PH3"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"1111",	Qcoeff(	-1.5068E+02	,6.4718E+00	,1.2588E-02	,1.4759E-05) )
	 ) ) );



  // COF2
  // Coeff:       1
  // Quality:   16.72
  q_data.push_back
    ( QRecord
      ( NAME("COF2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"269",	Qcoeff(	-5.4180E+03	,1.8868E+02	,-3.3139E-01	,1.8650E-03) )
	 ) ) );



  // SF6
  // Coeff:       1
  // Quality:    ----
  q_data.push_back
    ( QRecord
      ( NAME("SF6"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"29",	Qcoeff(	-1.0000E+00	,0.0000E+00	,0.0000E+00	,0.0000E+00) )
	 ) ) );



  // H2S
  // Coeff:       1      1      1      2
  // Quality:    0.60   ----   ----   0.40
  q_data.push_back
    ( QRecord
      ( NAME("H2S"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"121",	Qcoeff(	-1.5521E+01	,8.3130E-01	,3.3656E-03	,-8.5691E-07) ),
	 REC(	"141",	Qcoeff(	-1.5561E+01	,8.3337E-01	,3.3744E-03	,-8.5937E-07) ),
	 REC(	"131",	Qcoeff(	-6.2170E+01	,3.3295E+00	,1.3480E-02	,-3.4323E-06) ),
	 REC(	"122",	Qcoeff(	-1.7216E+01	,7.3676E-01	,2.8697E-03	,-1.1650E-06) )
	 ) ) );



  // HCOOH
  // Coeff:       1      2      2      2
  // Quality:   12.54   6.53   0.55   0.57
  q_data.push_back
    ( QRecord
      ( NAME("HCOOH"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"1261",	Qcoeff(	-2.9550E+03	,1.0349E+02	,-1.3146E-01	,8.7787E-04) ),
	 REC(	"1361",	Qcoeff(	2.0081E+02	,-4.9619E+00	,2.3524E-01	,-4.0531E-04) ),
	 REC(	"2261",	Qcoeff(	-2.5614E+01	,8.4704E+00	,1.2417E-01	,-1.1758E-04) ),
	 REC(	"1262",	Qcoeff(	-5.7615E+01	,9.5919E+00	,1.0930E-01	,-1.0134E-04) )
	 ) ) );



  // HO2
  // Coeff:       1
  // Quality:    1.10
  q_data.push_back
    ( QRecord
      ( NAME("HO2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"166",	Qcoeff(	-1.5684E+02	,7.4450E+00	,2.6011E-02	,-9.2704E-07) )
	 ) ) );



  // O
  // Coeff:       1
  // Quality:   ----
  q_data.push_back
    ( QRecord
      ( NAME("O"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"6",	Qcoeff(	-1.0000E+00	,0.0000E+00	,0.0000E+00	,0.0000E+00) )
	 ) ) );



  // ClONO2
  // Coeff:       2      2
  // Quality:    1.87   1.85
  q_data.push_back
    ( QRecord
      ( NAME("ClONO2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"5646",	Qcoeff(	-4.9350E+03	,2.7375E+02	,7.2257E+00	,-4.2583E-03) ),
	 REC(	"7646",	Qcoeff(	-4.9554E+03	,2.7955E+02	,7.4136E+00	,-4.3714E-03) )
	 ) ) );



  // NO+
  // Coeff:       1
  // Quality:    0.01
  q_data.push_back
    ( QRecord
      ( NAME("NO+"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"46",	Qcoeff(	9.1798E-01	,1.0416E+00	,-1.1614E-05	,2.4499E-08) )
	 ) ) );



  // OClO
  // Coeff:       2      2
  // Quality:    1.55   0.75
  q_data.push_back
    ( QRecord
      ( NAME("OClO"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"656",	Qcoeff(	-4.4583E+03	,1.2893E+02	,2.1655E-01	,-4.4496E-05) ),
	 REC(	"676",	Qcoeff(	-3.6048E+03	,1.1569E+02	,2.5805E-01	,-7.0900E-05) )
	 ) ) );



  // BrO
  // Coeff:       2      2
  // Quality:    0.09   0.09
  q_data.push_back
    ( QRecord
      ( NAME("BrO"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"96",	Qcoeff(	-8.8945E+00	,1.3384E+01	,-1.6496E-03	,1.4828E-06) ),
	 REC(	"16",	Qcoeff(	-8.4100E+00	,1.3433E+01	,-1.6322E-03	,1.4620E-06) )
	 ) ) );



  // H2SO4
  // Coeff:       2
  // Quality:    0.39
  q_data.push_back
    ( QRecord
      ( NAME("H2SO4"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"126",	Qcoeff(	-6.6207E+03	,2.6856E+02	,1.0187E+00	,-4.0717E-04) )
	 ) ) );



  // Cl2O2
  // Coeff:       2      2
  // Quality:    0.29   0.30
  q_data.push_back
    ( QRecord
      ( NAME("Cl2O2"),
	ISOTOPE
	(//	Name		c0		c1		c2		c3
	 //			|		|		|		|
	 REC(	"565",	Qcoeff(	-1.6546E+04	,7.0845E+02	,2.9249E+00	,-1.2097E-03) ),
	 REC(	"765",	Qcoeff(	-1.6863E+04	,7.2610E+02	,3.0201E+00	,-1.2538E-03) )
	 ) ) );


}


void check_q_data()
{
  /* routine assures that the species_data and the q_data are sorted
     the same way, for all molecules and isotopes, safety check. */
  
  extern ARRAY<QRecord> q_data;
  extern ARRAY<SpeciesRecord> species_data;

  for (size_t i=0; i<q_data.size(); ++i)
    {
      const QRecord& qd = q_data[i];
      const SpeciesRecord& sd = species_data[i];


      // check whether the species names are sorted the same way in
      // species_data and q_data
      assert( qd.Name() == sd.Name() );
      
      for (size_t j=0;j<qd.Isotope().size(); ++j)
	{

	  // check whether all isotopes of a species are sorted the
	  // same way in species_data and q_data
	  assert( qd.Isotope()[j].Name() == sd.Isotope()[j].Name() );
	}
    }
}

