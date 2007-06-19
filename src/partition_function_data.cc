/* Copyright (C) 2000-2007
   Stefan Buehler <sbuehler@ltu.se>,
   Axel von Engeln <engeln@uni-bremen.de>,
   Carmen Verdes <cverdes@uni-bremen.de>

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

  This is the file from arts-1-0, back-ported to arts-1-1.

  \author Axel von Engeln,  C. Verdes,
  \date 2000-08-21 
   */

#include "arts.h"
#include "make_array.h"
#include "absorption.h"

/*
  Default values for a linear molecule:

  c0        c1       c2       c3
  0.0       1.0      0.0      0.0

  and for non linear molecules:

  c0        c1       c2           c3
  -76.2339  6.95356  1.30510e-07  -1.17068e-11

  Here we assumed that either

  linear: Q = T
  non-linear: Q = T^1.5

  These numbers are from Axel von Engeln, May 2005.
*/

/*! \name Some #defines for better readability */
//@{ 
#define Qcoeff       MakeArray<Numeric>         
//@} 


/**
  Define partition function coefficients lookup data.

  <h1>General Remarks</h1>
  <!------------------------>

  This function contains the polynomial coefficients of the partition
  function for each isotope of each species. The sorting of the array
  has to match the species_data entries, a safety check routine
  assures this order.

  The information on the partition function are
  usually calculated using the latest version pf TIPS program, provided by B. Gamache.
  he latest TIPS code, TIPS_2003, can be obtained from the web page 
  (http://faculty.uml.edu/robert_gamache/) in the Software and Data section. 
  These include the data for all
  molecular species found in the 2000 HITRAN database, and, in addition, data for 13 other
  isotopomers/isotopologues of ozone and carbon dioxide.
  The calculations address the corrections suggested by Goldman et
  al. The calculations of the partition functions consider the temperature
  range 70--3000K, in steps of 25K to be applicable to a variety of remote sensing needs.
  The partition functions for any temperatures other than the ones at
  which the data is provided, are calculated by a Lagrange 4-point
  interpolation.  
  

  The coefficients for species which are not covered in HITRAN are
  calculated from JPL values. The JPL catalogue has a different way to
  calculate the partition function, it gives the partition function at
  specific temperatures: 300, 225, 150, 75, 37.5, 18.75, 9.375 K and
  an interpolation scheme is given for values inbetween. The partition
  functions are proportional to temperature T^1.5 for non-linear
  molecules (degrees of freedom: 3) and proportional to T for linear
  molecules (degrees of freedom 2). Nothing is known so far about the
  way how to treat atoms.

  The polynomial coefficients for JPL catalogue values are calculated
  with a matlab program. This calculates the partition function for
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

  In order to see the level of agreement of the to partition functions datasets, HITRAN 
  and JPL, a intercomparison has been made. 
  It was found that, with some exceptions, the two datasets are in a good
  agreement, if a vibrational correction is applied to the JPL partition
  functions. Only for 16 molecules (out of 66 considered molecules), the deviations in
  the two datasets are larger than 1\%. The largest discrepancies appear, for
  ClONO$_2$ (23\% at 150~K). Some deviations, however much smaller (around 2\%, or less),
  are found, e.g., for HCl, HOBr, and HBr.  

  <h1>Source for entries</h1>
  <!------------------------>

  <dl>
  <dt> Molecule Name:
  <dd> Arts convention

  <dt> Isotope Name:
  <dd> Arts convention

  <dt> Coefficients: 
  <dd> Generally taken from the new program provided by B. Gamache, version 2003. 
  Partition functions for species included only in JPL, except for BrO, were
  calculated according to the JPL recommended scheme, and a 3rd order
  polynomial fit was performed (with the idl program mentioned above).
  The partition function data for BrO is provided by Agnes Perrin, Orsay, France.
 
  \author Carmen Verdes  
  \date 2003-11-12 */

void spec(Array<SpeciesRecord>::iterator& is,
          Array<IsotopeRecord>::iterator& ii,
          String name);

void iso(Array<IsotopeRecord>::iterator& ii,
         String name,
         const ArrayOfNumeric& coeff);


void define_partition_species_data()
{
  extern Array<SpeciesRecord> species_data;

  Array<SpeciesRecord>::iterator it_species = species_data.begin();
  Array<IsotopeRecord>::iterator it_isotope;


  // H2O
  // Coeff:       1      1      1      1      1      1      2
  //
  // There are dummy entries here for the continuum tags. Of course,
  // continua need no partition functions, but the entries must be
  // present, so that the consistency check between species_data and
  // partition_function_data is successful. Thomas, you will also need
  // to add more dummy entries here in case you add more continuum tags.
  //
  // Excetp isotope 262,  data for all other isotopes is obtained by using the new code provided by Bob Gamasche (falged with HITRAN). The missing data is taken from JPL. 
  spec(it_species, it_isotope, "H2O");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"161" ,	Qcoeff(-6.065594e+00,  2.907027e-01,  1.246245e-03,  -5.606119e-07 )); //HITRAN 
iso(it_isotope,	"181" ,	Qcoeff(-7.220624e+00,  2.945347e-01,  1.250362e-03,  -5.554638e-07 )); //HITRAN 
iso(it_isotope,	"171" ,	Qcoeff(-4.668105e+01,  1.819186e+00,  7.137470e-03,  -2.670352e-06 )); //HITRAN 
iso(it_isotope,	"162" ,	Qcoeff(-4.084466e+01,  1.484533e+00,  5.953330e-03,  -2.359695e-06 )); //HITRAN 
iso(it_isotope,	"182" ,	Qcoeff(-3.529770e+01,  1.503267e+00,  6.020059e-03,  -2.389284e-06 )); //HITRAN 
iso(it_isotope,	"172" ,	Qcoeff(-2.098457e+02,  8.959286e+00,  3.593721e-02,  -1.428880e-05 )); //HITRAN 
 iso(it_isotope,	"262" ,	Qcoeff(-3.572493e+01,  1.652500e+00,  7.633309e-03,  -3.770940e-06 )); // JPL
  iso(it_isotope,       "SelfContStandardType",      Qcoeff( 0  ,0              ,0              ,0          ) );
  iso(it_isotope,       "ForeignContStandardType",   Qcoeff( 0  ,0              ,0              ,0          ) );
  iso(it_isotope,       "ForeignContMaTippingType",  Qcoeff( 0  ,0              ,0              ,0          ) ); 
  iso(it_isotope,       "ContMPM93",               Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "SelfContCKDMT100",        Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "ForeignContCKDMT100",     Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "SelfContCKD222",          Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "ForeignContCKD222",       Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "SelfContCKD242",          Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "ForeignContCKD242",       Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "SelfContCKD24",           Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "ForeignContCKD24",        Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "ForeignContATM01",        Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "CP98", Qcoeff( 0       ,0              ,0              ,0          ) );
  iso(it_isotope,       "MPM87",Qcoeff( 0       ,0              ,0              ,0          ) );
  iso(it_isotope,       "MPM89",Qcoeff( 0       ,0              ,0              ,0          ) );
  iso(it_isotope,       "MPM93",Qcoeff( 0       ,0              ,0              ,0          ) );
  iso(it_isotope,       "PWR98",Qcoeff( 0       ,0              ,0              ,0          ) );


  // CO2
  // Coeff:       1      1      1      1      1      1      1      1
  spec(it_species, it_isotope, "CO2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"626" ,	Qcoeff(-1.720718e+00,  9.669217e-01,  -8.277298e-04,  2.891070e-06 )); //HITRAN 
iso(it_isotope,	"636" ,	Qcoeff(-1.850250e+00,  1.912107e+00,  -1.599677e-03,  5.955462e-06 )); //HITRAN 
iso(it_isotope,	"628" ,	Qcoeff(-2.989446e+00,  2.041095e+00,  -1.732748e-03,  6.174831e-06 )); //HITRAN 
iso(it_isotope,	"627" ,	Qcoeff(-2.256240e+01,  1.197737e+01,  -1.036863e-02,  3.618820e-05 )); //HITRAN 
iso(it_isotope,	"638" ,	Qcoeff(-1.882095e+00,  4.025528e+00,  -3.299047e-03,  1.266725e-05 )); //HITRAN 
iso(it_isotope,	"637" ,	Qcoeff(-1.788894e+01,  2.358158e+01,  -1.957976e-02,  7.389405e-05 )); //HITRAN 
iso(it_isotope,	"828" ,	Qcoeff(-1.818540e+00,  1.086818e+00,  -9.427241e-04,  3.352295e-06 )); //HITRAN 
iso(it_isotope,	"728" ,	Qcoeff(-2.087791e+01,  1.266486e+01,  -1.091287e-02,  3.874034e-05 )); //HITRAN 
iso(it_isotope,       "CKD241",           Qcoeff( 0        ,0              ,0              ,0          ) );
iso(it_isotope,       "CKDMT100",         Qcoeff( 0        ,0              ,0              ,0          ) );
iso(it_isotope,       "SelfContPWR93",    Qcoeff( 0        ,0              ,0              ,0          ) );
iso(it_isotope,       "ForeignContPWR93", Qcoeff( 0        ,0              ,0              ,0          ) );



  // O3
  // Coeff:       1      1      1      1      1
  spec(it_species, it_isotope, "O3");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"666" ,	Qcoeff(-2.773214e+02,  8.175293e+00,  6.892651e-03,  2.842028e-05 )); //HITRAN 
iso(it_isotope,	"668" ,	Qcoeff(-5.978029e+02,  1.759117e+01,  1.353516e-02,  6.440030e-05 )); //HITRAN 
iso(it_isotope,	"686" ,	Qcoeff(-3.005190e+02,  8.726453e+00,  5.976672e-03,  3.241643e-05 )); //HITRAN 
iso(it_isotope,	"667" ,	Qcoeff(-3.454638e+03,  1.018144e+02,  8.249751e-02,  3.631247e-04 )); //HITRAN 
iso(it_isotope,	"676" ,	Qcoeff(-1.735693e+03,  5.072998e+01,  3.877763e-02,  1.821985e-04 )); //HITRAN 



  // N2O
  // Coeff:       1      1      1      1      1
  spec(it_species, it_isotope, "N2O");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"446" ,	Qcoeff(3.478254e+01,  1.530195e+01,  -1.120080e-02,  5.472145e-05 )); //HITRAN 
iso(it_isotope,	"456" ,	Qcoeff(3.479618e+01,  1.002537e+01,  -6.789834e-03,  3.681093e-05 )); //HITRAN 
iso(it_isotope,	"546" ,	Qcoeff(2.435117e+01,  1.055152e+01,  -7.756090e-03,  3.819981e-05 )); //HITRAN 
iso(it_isotope,	"448" ,	Qcoeff(4.066999e+01,  1.615921e+01,  -1.180945e-02,  5.883212e-05 )); //HITRAN 
iso(it_isotope,	"447" ,	Qcoeff(2.050163e+02,  9.473303e+01,  -7.029656e-02,  3.426216e-04 )); //HITRAN 


  // CO
  // Coeff:       1      1      1      1      1      1
  spec(it_species, it_isotope, "CO");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"26" ,	Qcoeff(3.243148e-01,  3.601229e-01,  1.538205e-06,  2.385704e-09 )); //HITRAN 
iso(it_isotope,	"36" ,	Qcoeff(4.632274e-01,  7.560062e-01,  -8.390593e-06,  2.229242e-08 )); //HITRAN 
iso(it_isotope,	"28" ,	Qcoeff(2.874382e-01,  3.786605e-01,  -5.551926e-07,  5.629838e-09 )); //HITRAN 
iso(it_isotope,	"27" ,	Qcoeff(1.697400e+00,  2.220079e+00,  -4.074631e-06,  3.291954e-08 )); //HITRAN 
iso(it_isotope,	"38" ,	Qcoeff(6.558005e-01,  7.928532e-01,  4.443750e-06,  3.520833e-09 )); //HITRAN 
iso(it_isotope,	"37" ,	Qcoeff(3.990599e+00,  4.641927e+00,  2.855732e-05,  1.499385e-08 )); //HITRAN 

  // CH4
  // Coeff:       1      1      1
  spec(it_species, it_isotope, "CH4");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"211" ,	Qcoeff(-3.640461e+01,  1.202398e+00,  3.005684e-03,  2.911372e-07 )); //HITRAN 
iso(it_isotope,	"311" ,	Qcoeff(-7.385939e+01,  2.419567e+00,  5.941999e-03,  6.864449e-07 )); //HITRAN 
iso(it_isotope,	"212" ,	Qcoeff(-3.003903e+02,  9.769371e+00,  2.411804e-02,  2.704667e-06 )); //HITRAN 



  // O2
  // Coeff:       1      1      1
  spec(it_species, it_isotope, "O2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"66" ,	Qcoeff(4.016432e-01,  7.315888e-01,  -3.313678e-05,  6.642877e-08 )); //HITRAN 
iso(it_isotope,	"68" ,	Qcoeff(-3.922253e+00,  1.551651e+00,  -8.580045e-05,  1.716056e-07 )); //HITRAN 
iso(it_isotope,	"67" ,	Qcoeff(-2.757545e+01,  9.118689e+00,  -7.483006e-04,  1.332269e-06 )); //HITRAN 
  iso(it_isotope,       "CIAfunCKDMT100", Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "v0v0CKDMT100",   Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "v1v0CKDMT100",   Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "SelfContStandardType",   Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "SelfContMPM93",  Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "SelfContPWR93",  Qcoeff( 0    ,0              ,0              ,0          ) );
  iso(it_isotope,       "PWR98",          Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "PWR93",          Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "PWR88",          Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "MPM93",          Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "MPM92",          Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "MPM89",          Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "MPM87",          Qcoeff( 0     ,0              ,0              ,0          ) );
  iso(it_isotope,       "MPM85",          Qcoeff( 0     ,0              ,0              ,0          ) );

  // NO
  // Coeff:       1      1      1
  spec(it_species, it_isotope, "NO");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"46" ,	Qcoeff(-5.824308e+01,  3.025484e+00,  4.976571e-03,  -5.060093e-06 )); //HITRAN 
iso(it_isotope,	"56" ,	Qcoeff(-4.036081e+01,  2.091668e+00,  3.435242e-03,  -3.490987e-06 )); //HITRAN 
iso(it_isotope,	"48" ,	Qcoeff(-6.255837e+01,  3.205744e+00,  5.176248e-03,  -5.223151e-06 )); //HITRAN 

  // SO2
  // Coeff:       1      1      2      2
  spec(it_species, it_isotope, "SO2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"626" ,	Qcoeff(-3.406710e+02,  1.214516e+01,  1.995262e-02,  5.157669e-05 )); //HITRAN 
iso(it_isotope,	"646" ,	Qcoeff(-3.389056e+02,  1.215747e+01,  2.023113e-02,  5.153272e-05 )); //HITRAN 
iso(it_isotope,       "636",  Qcoeff( 5.8740E+02      ,1.2472E+01     ,2.9113E-01     ,-1.6236E-04) ); //JPl
iso(it_isotope,       "628",  Qcoeff( 3.1299E+02      ,6.6372E+00     ,1.5485E-01     ,-8.6343E-05) ); //JPl



  // NO2
  // Coeff:       1
  spec(it_species, it_isotope, "NO2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"646" ,	Qcoeff(-8.761726e+02,  2.829842e+01,  5.398242e-02,  5.194329e-05 ));  //HITRAN 



  // NH3
  // Coeff:       1      1      2
  spec(it_species, it_isotope, "NH3");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"4111" ,	Qcoeff(-9.698124e+01,  3.402711e+00,  8.958578e-03,  1.157044e-06 )); //HITRAN 
iso(it_isotope,	"5111" ,	Qcoeff(-6.520038e+01,  2.279068e+00,  5.958356e-03,  8.170489e-07 )); //HITRAN 
 iso(it_isotope,	"4112" ,	Qcoeff(9.278991e+00,  4.053839e+00,  3.148529e-02,  -8.153125e-06 )); //JPL


  // HNO3
  // Coeff:       1
  spec(it_species, it_isotope, "HNO3");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"146" ,	Qcoeff(-3.402033e+04,  7.965238e+02,  -2.403160e+00,  8.593868e-03 )); //HITRAN 


  // OH
  // Coeff:       1      1      1
  spec(it_species, it_isotope, "OH");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"61" ,	Qcoeff(6.198722e+00,  1.870893e-01,  3.099551e-04,  -3.229806e-07 )); //HITRAN 
iso(it_isotope,	"81" ,	Qcoeff(6.173190e+00,  1.884492e-01,  3.126020e-04,  -3.263942e-07 )); //HITRAN 
iso(it_isotope,	"62" ,	Qcoeff(4.103720e+00,  5.095633e-01,  8.899807e-04,  -9.103002e-07 )); //HITRAN 


  // HF
  // Coeff:       1      2
  spec(it_species, it_isotope, "HF");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"19" ,	Qcoeff(1.472238e+00,  1.343685e-01,  3.150221e-06,  -2.120225e-09 )); //HITRAN 
 iso(it_isotope,	"29" ,	Qcoeff(3.375585e-01,  6.403473e-02,  3.134983e-07,  -3.970786e-11 )); //JPL


  // HCl
  // Coeff:       1      1      2      2
  spec(it_species, it_isotope, "HCl");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"15" ,	Qcoeff(2.729314e+00,  5.328097e-01,  8.234868e-07,  5.619026e-09 )); //HITRAN 
iso(it_isotope,	"17" ,	Qcoeff(2.719350e+00,  5.335676e-01,  2.054102e-06,  2.061213e-09 )); //HITRAN 
 iso(it_isotope,	"25" ,	Qcoeff(1.355208e+00,  5.155418e-01,  3.328246e-06,  1.718278e-12 )); //JPL
 iso(it_isotope,	"27" ,	Qcoeff(1.359929e+00,  5.170804e-01,  3.358101e-06,  -1.087936e-11 )); //JPL



  // HBr
  // Coeff:       1      1
  spec(it_species, it_isotope, "HBr");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"19" ,	Qcoeff(2.936148e+00,  6.629899e-01,  1.604872e-05,  -1.593934e-08 )); //HITRAN 
iso(it_isotope,	"11" ,	Qcoeff(2.875136e+00,  6.637710e-01,  1.449833e-05,  -1.498201e-08 )); //HITRAN 


  // HI
  // Coeff:       1
  spec(it_species, it_isotope, "HI");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"17" ,	Qcoeff(4.226561e+00,  1.295818e+00,  1.611346e-05,  -7.882228e-09 )); //HITRAN 


  // ClO
  // Coeff:       1      1
  spec(it_species, it_isotope, "ClO");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"56" ,	Qcoeff(1.290486e+02,  6.369550e+00,  1.441861e-02,  -1.211120e-07 )); //HITRAN 
iso(it_isotope,	"76" ,	Qcoeff(1.306461e+02,  6.492672e+00,  1.457301e-02,  1.142879e-07 )); //HITRAN


  // OCS
  // Coeff:       1      1      1      1      1
  spec(it_species, it_isotope, "OCS");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"622" ,	Qcoeff(1.199103e+01,  3.484349e+00,  -3.172632e-03,  1.757090e-05 )); //HITRAN 
iso(it_isotope,	"624" ,	Qcoeff(1.055761e+01,  3.598837e+00,  -3.406838e-03,  1.836238e-05 )); //HITRAN 
iso(it_isotope,	"632" ,	Qcoeff(3.246621e+01,  6.852374e+00,  -5.819381e-03,  3.599002e-05 )); //HITRAN 
iso(it_isotope,	"623" ,	Qcoeff(4.848356e+01,  1.411918e+01,  -1.292079e-02,  7.151233e-05 )); //HITRAN 
iso(it_isotope,	"822" ,	Qcoeff(1.444298e+01,  3.686311e+00,  -3.307686e-03,  1.920205e-05 )); //HITRAN 


  // H2CO
  // Coeff:       1      1      1      2      2
 spec(it_species, it_isotope, "H2CO");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"1126" ,	Qcoeff(-1.734031e+02,  5.682345e+00,  1.504875e-02,  7.509330e-07 )); //HITRAN 
iso(it_isotope,	"1136" ,	Qcoeff(-3.529337e+02,  1.160844e+01,  3.109193e-02,  1.153082e-06 )); //HITRAN 
iso(it_isotope,	"1128" ,	Qcoeff(-1.781662e+02,  5.905635e+00,  1.604851e-02,  3.936717e-07 )); //HITRAN 
 iso(it_isotope,	"1226" ,	Qcoeff(-5.332528e+01,  2.914098e+00,  1.444437e-02,  -6.565213e-06 )); //Jpl
iso(it_isotope,	"2226" ,	Qcoeff(-2.847116e+02,  1.672849e+01,  8.661739e-02,  -3.736935e-05 )); //Jpl



  // HOCl
  // Coeff:       1      1
  spec(it_species, it_isotope, "HOCl");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"165" ,	Qcoeff(-1.219223e+03,  3.989396e+01,  7.529869e-02,  8.046020e-05 )); //HITRAN 
iso(it_isotope,	"167" ,	Qcoeff(-1.215084e+03,  4.025848e+01,  7.807742e-02,  7.992701e-05 )); //HITRAN 


  // N2
  // Coeff:       1
  spec(it_species, it_isotope, "N2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"44" ,	Qcoeff(1.704255e+00,  1.562748e+00,  2.437406e-05,  -1.677703e-08 )); //HITRAN 
  iso(it_isotope,       "SelfContMPM93",        Qcoeff( 0       ,0              ,0              ,0   ) );
  iso(it_isotope,       "SelfContPWR93",        Qcoeff( 0       ,0              ,0              ,0   ) );
  iso(it_isotope,       "SelfContStandardType", Qcoeff( 0       ,0              ,0              ,0   ) );
  iso(it_isotope,       "SelfContBorysow",      Qcoeff( 0       ,0              ,0              ,0   ) );
  iso(it_isotope,       "CIArotCKDMT100",       Qcoeff( 0       ,0              ,0              ,0   ) );
  iso(it_isotope,       "CIAfunCKDMT100",       Qcoeff( 0       ,0              ,0              ,0   ) );
  iso(it_isotope,       "DryContATM01",         Qcoeff( 0       ,0              ,0              ,0   ) );


  // HCN
  // Coeff:       1      1      1      2
  spec(it_species, it_isotope, "HCN");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"124" ,	Qcoeff(-5.935227e+00,  3.077616e+00,  -2.476330e-03,  7.991253e-06 )); //HITRAN 
iso(it_isotope,	"134" ,	Qcoeff(-1.010578e+01,  6.290094e+00,  -4.988065e-03,  1.641309e-05 )); //HITRAN 
iso(it_isotope,	"125" ,	Qcoeff(-3.253498e+00,  2.118761e+00,  -1.680616e-03,  5.582555e-06 )); //HITRAN 
 iso(it_isotope,	"224" ,	Qcoeff(6.626957e+01,  8.670873e-01,  3.148000e-03,  -2.052228e-06 )); //JPL
 


  // CH3Cl
  // Coeff:       1      1
  spec(it_species, it_isotope, "CH3Cl");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"215" ,	Qcoeff(-1.140936e+04,  3.073757e+02,  4.383730e-02,  1.249421e-03 )); //HITRAN 
iso(it_isotope,	"217" ,	Qcoeff(-1.159736e+04,  3.123035e+02,  4.438509e-02,  1.269305e-03 )); //HITRAN 



  // H2O2
  // Coeff:       1
  spec(it_species, it_isotope, "H2O2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"1661" ,	Qcoeff(-3.865211e+02,  1.286868e+01,  3.910416e-02,  1.145394e-04 )); //HITRAN 
 



  // C2H2
  // Coeff:       1      1
  spec(it_species, it_isotope, "C2H2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"1221" ,	Qcoeff(-8.684002e+00,  1.453883e+00,  -2.597724e-03,  8.482153e-06 ));//HITRAN 
iso(it_isotope,	"1231" ,	Qcoeff(-3.468599e+01,  5.815575e+00,  -1.039390e-02,  3.393990e-05 ));//HITRAN 
 



  // C2H6
  // Coeff:       1
  spec(it_species, it_isotope, "C2H6");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"1221" ,	Qcoeff(-9.118157e+03,  2.088364e+02,  -4.404385e-01,  2.188428e-03 ));//HITRAN 



  // PH3
  // Coeff:       2
  spec(it_species, it_isotope, "PH3");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"1111" ,	Qcoeff(-2.426718e+02,  7.338849e+00,  1.131445e-02,  1.262244e-05 ));// JPL 



  // COF2
  // Coeff:       1
  spec(it_species, it_isotope, "COF2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"269" ,	Qcoeff(-8.322642e+03,  2.144407e+02,  -3.498616e-01,  1.755888e-03 ));//HITRAN 



  // SF6
  // Coeff:       1
  spec(it_species, it_isotope, "SF6");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"29" ,	Qcoeff(-1.668432e+06,  2.850128e+04,  -1.561230e+02,  3.288986e-01 ));//HITRAN 



  // H2S
  // Coeff:       1      1      1      2
  spec(it_species, it_isotope, "H2S");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |

iso(it_isotope,	"121" ,	Qcoeff(-2.308888e+01,  9.052647e-01,  3.237531e-03,  -9.823621e-07 )); //HITRAN 
iso(it_isotope,	"141" ,	Qcoeff(-2.333981e+01,  9.102537e-01,  3.233485e-03,  -9.665574e-07 )); //HITRAN 
iso(it_isotope,	"131" ,	Qcoeff(-9.329309e+01,  3.636877e+00,  1.291822e-02,  -3.864368e-06 )); //HITRAN 
 iso(it_isotope,	"122" ,	Qcoeff(-1.512671e+01,  6.851018e-01,  3.158080e-03,  -1.563931e-06 )); //JPL


  // HCOOH
  // Coeff:       1      2      2      2
  spec(it_species, it_isotope, "HCOOH");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"1261" ,	Qcoeff(-4.370811e+03,  1.141311e+02,  -1.217474e-01,  7.859656e-04 ));//HITRAN 
iso(it_isotope,	"1361" ,	Qcoeff(-4.910213e+03,  5.115094e+01,  3.433096e-02,  -1.340898e-04 ));// JPL 
iso(it_isotope,	"2261" ,	Qcoeff(3.823001e+02,  5.455419e+00,  1.108040e-01,  -5.086754e-06 ));// JPL 
iso(it_isotope,	"1262" ,	Qcoeff(8.193393e+02,  2.222546e+00,  1.070970e-01,  3.255965e-06 ));// JPL 




  // HO2
  // Coeff:       1
  spec(it_species, it_isotope, "HO2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"166" ,	Qcoeff(-2.341264e+02,  8.164256e+00,  2.506193e-02,  -3.012599e-06 )); //HITRAN 



  // O
  // Coeff:       1
  spec(it_species, it_isotope, "O");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
  iso(it_isotope,       "6",    Qcoeff( -1.0000E+00     ,0.0000E+00     ,0.0000E+00     ,0.0000E+00) ); //HITRAN 



  // ClONO2
  // Coeff:       1      1
  spec(it_species, it_isotope, "ClONO2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"5646" ,	Qcoeff(-2.052890e+06,  3.638094e+04,  -1.995279e+02,  5.224687e-01 ));//HITRAN 
iso(it_isotope,	"7646" ,	Qcoeff(-2.104484e+06,  3.729925e+04,  -2.045781e+02,  5.357327e-01 ));//HITRAN 


  // NO+
  // Coeff:       1
  spec(it_species, it_isotope, "NO+");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"46" ,	Qcoeff(1.125969e+00,  1.047028e+00,  1.174546e-05,  -1.519278e-08 ));//HITRAN 



  // OClO
  // Coeff:       2      2
  spec(it_species, it_isotope, "OClO");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"656" ,	Qcoeff(-1.617389e+03,  6.991068e+01,  5.003075e-01,  -1.442758e-04 ));// JPL 
iso(it_isotope,	"676" ,	Qcoeff(7.964396e+02,  4.768587e+01,  5.283347e-01,  -8.232128e-05 ));// JPL  

// The partition function for BrO is provided by Agned Perin, Orsay France!!!
  // BrO
  // Coeff:       3      3
  spec(it_species, it_isotope, "BrO");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"96" ,	Qcoeff(-4.084622e+01,  1.427999e+01,  -1.011647e-02,  2.783630e-05 )); //HITRAN 
iso(it_isotope,	"16" ,	Qcoeff(-4.118468e+01,  1.434034e+01,  -1.016302e-02,  2.795965e-05 )); //HITRAN 

  // H2SO4
  // Coeff:       2
  spec(it_species, it_isotope, "H2SO4");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"126" ,	Qcoeff(-5.913199e+03,  2.485770e+02,  1.140269e+00,  -5.679165e-04 ));// JPL 



  // Cl2O2
  // Coeff:       2      2
  spec(it_species, it_isotope, "Cl2O2");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"565" ,	Qcoeff(6.215326e+05,  -7.121447e+03,  2.784834e+01,  2.147458e-02 ));// JPL
iso(it_isotope,	"765" ,	Qcoeff(6.399192e+05,  -7.332314e+03,  2.866224e+01,  2.210953e-02 ));// JPL


  // HOBr
  // Coeff:       1      1
  spec(it_species, it_isotope, "HOBr");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"169" ,	Qcoeff(-1.665575e+03,  5.687767e+01,  9.982304e-02,  1.705212e-04 )); //HITRAN 
iso(it_isotope,	"161" ,	Qcoeff(-1.631140e+03,  5.625451e+01,  1.012339e-01,  1.676169e-04 )); //HITRAN 



  // C2H4
  // Coeff:       1      1
  spec(it_species, it_isotope, "C2H4");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"221" ,	Qcoeff(-1.379496e+03,  3.408740e+01,  -2.321387e-02,  1.682474e-04 )); //HITRAN 
iso(it_isotope,	"231" ,	Qcoeff(-5.653328e+03,  1.396050e+02,  -9.531910e-02,  6.891171e-04 )); //HITRAN 


  // CH3CN
  // Coeff:       2      2      2      2      2
  spec(it_species, it_isotope, "CH3CN");
  //                    Name                    c0              c1              c2              c3
  //                    |                       |               |               |               |
iso(it_isotope,	"211124" ,	Qcoeff(1.706820e+03,  1.093287e+00,  4.255850e-01,  3.367172e-05 ));// JPL 
iso(it_isotope,	"311124" ,	Qcoeff(-1.172596e+03,  4.973615e+01,  2.285735e-01,  -1.135942e-04 ));// JPL  
iso(it_isotope,	"211134" ,	Qcoeff(-1.139329e+03,  4.832504e+01,  2.220882e-01,  -1.103713e-04 ));// JPL
iso(it_isotope,	"211125" ,	Qcoeff(-3.861117e+02,  1.654635e+01,  7.638250e-02,  -3.776153e-05 ));// JPL   
iso(it_isotope,	"211224" ,	Qcoeff(-3.483734e+02,  1.464417e+01,  6.717486e-02,  -3.345710e-05 ));// JPL  


  // HNC
  // Coeff:       2      2      2      2
  spec(it_species, it_isotope, "HNC");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
iso(it_isotope,	"142" ,	Qcoeff(3.333499e-01,  4.595243e-01,  1.502307e-06,  2.413631e-13 ));// JPL 
iso(it_isotope,	"143" ,	Qcoeff(7.032558e-02,  4.813463e-01,  -6.651923e-06,  -1.954964e-11 ));// JPL  
iso(it_isotope,	"152" ,	Qcoeff(1.589700e-01,  4.708129e-01,  -3.866420e-06,  4.006830e-11 ));// JPL  
iso(it_isotope,	"242" ,	Qcoeff(2.865158e-01,  5.465990e-01,  2.142689e-07,  1.733211e-11 ));// JPL  

  // particle
  // Coeff:       1      1
  spec(it_species, it_isotope, "liquidcloud");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
  iso(it_isotope,       "MPM93", Qcoeff( 0      ,0              ,0              ,0  ) );

  // particle
  // Coeff:       1      1
  spec(it_species, it_isotope, "icecloud");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
  iso(it_isotope,       "MPM93",Qcoeff( 0              ,0               ,0              ,0  ) );

  // particle
  // Coeff:       1      1
  spec(it_species, it_isotope, "rain");
  //                    Name            c0              c1              c2              c3
  //                    |               |               |               |               |
  iso(it_isotope,       "MPM93",Qcoeff( 0              ,0               ,0              ,0  ) );

}


void spec(Array<SpeciesRecord>::iterator& is,
          Array<IsotopeRecord>::iterator& ii,
          String name)
{
  
  assert( name == is->Name() );

  if (name ==  is->Name() ){}
  ii = is->Isotope().begin();

  is++;

}


void iso(Array<IsotopeRecord>::iterator& ii,
         String name,
         const ArrayOfNumeric& coeff)
{
  assert( name == ii->Name() );

  if (name ==  ii->Name() ){}
  ii->SetPartitionFctCoeff(coeff);

  ii++;

}
