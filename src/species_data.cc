/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>,
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
  \file   species_data.cc
  \brief  Implementation of function define_species_data().

  This file contains the definition of this function and nothing
  else. You can add new records here if you want to extend the
  capability of ARTS.

  \author Stefan Buehler, Axel von Engeln
  \date 2000-08-10 */

#include "arts.h"
#include "make_array.h"
#include "absorption.h"


/*! The lookup information for all the different species. */
Array<SpeciesRecord> species_data;


/*! \name Some #defines for better readability */
//@{ 
#define NAME(x) x                               
#define DEGFR(x) x                              
#define ISOTOPES     MakeArray<IsotopeRecord>   
#define REC          IsotopeRecord              
#define TAGS         MakeArray<Index>           
//@} 


/*!
  \brief Define species lookup data.

  \author Stefan Buehler, Axel von Engeln

  <h1>Source for entries</h1>
  <!------------------------>

  <dl>
  <dt> Molecule Name:
  <dd> Generally Hitran convention, unless otherwise indicated.

  <dt> Degrees of freedom: 
  <dd> Generally taken from the forward_4_96 release,
       file glob_def.c, unless otherwise indicated.

  <dt> Isotopes: 
  <dd> Generally Hitran convention, e.g., for the main H2O
       isotope: 161, thus the last digit of the sum of the
       neutrons and protons of the individual atoms is
       taken. Nevertheless, molecules like H2CO have only 3
       digits in Hitran, only one for the two H atoms. This leads
       to problems with the jpl catalogue, if this includes DHCO
       and DDCO. For these cases, a different scheme is taken as
       indicated.

  <dt> Isotopic Ratio: 
  <dd> Generally Hitran convention, unless otherwise indicated.
       A number code in the header gives the source of the isotopic ratios 
       for each isotope of each species:
       
       <table>
       <tr>
       <td> 1: <td> hitran 96 isotopic ratio, taken from cd: file
                    software/generic/tables_96.txt or HITRAN 2000 edition
                    (default, even if jpl and hitran isotopic ratios are 
                    available).
       <tr>
       <td> 2: <td> jpl isotopic ratio, taken from the documentation coming
                    along with the catalogue. latest catalogue version
                    extracted 27.07.00, can be found at
                    /pool/lookup/jpl/cat7_00/doc/d<tag_nr>.cat 
       <tr>
       <td> 3: <td> jpl isotopic ratio is multiplied with the maximum isotopic ratio
                    of this species found in hitran. is only performed when
                    isotopic ratios of 1 were found in the jpl catalogue.
       </table>

  <dt> Mass: 
//  <dd> Rounded atomic weight. (SAB: From looking at the table in forward_4_96,
//       glob_def.c, the difference between the actual mass and the mass
//       simply estimated from the Atom number is only 0.001. This seems
//       not worth the trouble. Anyway, the field for the mass is
//       there. Should anybody feel like adding the true numbers, just go
//       ahead.)
  <dd> mass as given in file MOLPARAM.TXT of HITRAN2000. 
       If tag is not present in HITRAN, the rounded value is stated. 
       Note that the mass is given in units of [g/mol].

  <dt >MY-tags:
  <dd> Extracted from file glob_def.c.

  <dt> HI-tags:
  <dd> From hitran 96, software/generic/tables_96.txt.

  <dt> JPL-tags:
  <dd> Collected in file 
       arts/aux/abundancies/tag_species.jpl, taken from
       last issue of jpl catalogue (7/00).
  </dl>

  Some more information can be found at /pool/lookup/jpl/cat7_00/abundancies,
  where the idl script that reads/converts the isotopic ratios is located.

  \author Axel von Engeln  
  \date   2000-08-08 
*/


// prototyping
void define_basic_species_data();
extern void define_partition_species_data();

void define_species_data()
{
  define_basic_species_data();
  define_partition_species_data();
}


void define_basic_species_data()
{
  extern Array<SpeciesRecord> species_data;

  // Initialize to zero, just in case:
  species_data.resize(0);

  /* Here's an empty template record entry:

  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
        DEGFR(3),
        ISOTOPES
        (//   Name,     Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC( ""        ,               ,       ,       ,       ,TAGS() ),
         REC( ""        ,               ,       ,       ,       ,TAGS() )
         ) ) );

  For good looks, keep the commas on the marks!

  */


  // H2O
  // Isotopic Ratio: 1 1 1 1 1 1 3 
  //
  // Some tags relate to empirical continuum correction terms. Thomas,
  // this would be the place to add additional continuum tags, should
  // this be necessary. Continuum tags must come after the other tags
  // and have -1 for all data entries, like in my example below. Not
  // even the isotipic ratio is used for continuum tags.
  //
  // The isotopic ratio of -1 is used to identify continuum tags in
  // the absorption routines!
  // 
  // You also have to change the entry in the file
  // partition_function_data.cc consistently! 
  //
  // 2001-05-30 Stefan Buehler: 
  // - Added isotope 172 (HITRAN 2000 tag 16)
  // - Gave HITRAN tag 15 to isotope 182, which was already in ARTS,
  //   but previously not in HITRAN. Changed isotopic ratio of this
  //   one to the HITRAN value, which is 6.23e-7, instead of the
  //   previous value of 6.11e-7.
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
        DEGFR(3),
        ISOTOPES
        (//   Name,             Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //                     |               |       |       |       |
         REC( "161"             ,0.99731702     ,18.010565    ,11     ,11     ,TAGS(18003, 18005) ),
         REC( "181"             ,0.00199983     ,20.014811    ,12     ,12     ,TAGS(20003) ),
         REC( "171"             ,0.00037200     ,19.014780    ,13     ,13     ,TAGS(19003) ),
         REC( "162"             ,0.000310693     ,19.016740    ,14     ,14     ,TAGS(19002) ),
         REC( "182"             ,6.23003E-07    ,21.020985    ,-1     ,15     ,TAGS(21001) ),
         REC( "172"             ,1.15853E-07    ,20.020956    ,-1     ,16     ,TAGS()      ),
         REC( "262"             ,2.2430204E-08  ,20.000000    ,-1     ,-1     ,TAGS(20001) ),
         REC( "SelfContStandardType"    ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "ForeignContStandardType" ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "ForeignContMaTippingType",-1.    ,-1.    ,-1     ,-1     ,TAGS()      ), 
         REC( "ContMPM93"               ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "SelfContCKDMT100"        ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "ForeignContCKDMT100"     ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "SelfContCKD24"           ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "ForeignContCKD24"        ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "ForeignContATM01"        ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "CP98"                    ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "MPM87"                   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "MPM89"                   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "MPM93"                   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      ),
         REC( "PWR98"                   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      )
         ) ) );

  // CO2 
  // (missing mainly in JPL, latest version (7/00) includes some isotopes)
  // Degrees of freedom from Schanda:`Physical Fundamentals of Remote Sensing'
  // Isotopic Ratios: 1 1 1 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CO2"),
        DEGFR(2),
        ISOTOPES
        (//   Name,     Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC( "626"     ,0.984204       ,43.989830    ,21     ,21     ,TAGS() ),
         REC( "636"     ,1.10574E-02    ,44.993185    ,22     ,22     ,TAGS() ),
         REC( "628"     ,3.94707E-03    ,45.994076    ,23     ,23     ,TAGS(46013) ),
         REC( "627"     ,7.33989E-04   ,44.994045    ,24     ,24     ,TAGS(45012) ),
         REC( "638"     ,4.43446E-05   ,46.997431    ,25     ,25     ,TAGS() ),
         REC( "637"     ,8.24623E-06   ,45.997400    ,26     ,26     ,TAGS() ),
         REC( "828"     ,3.95734E-06   ,47.998322    ,27     ,27     ,TAGS() ),
         REC( "728"     ,1.47180E-06   ,46.998291    ,28     ,28     ,TAGS() ),
	 REC( "CKDMT100"         ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "SelfContPWR93"    ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "ForeignContPWR93" ,-1.    ,-1.    ,-1     ,-1     ,TAGS())
         ) ) );
  

  // O3
  // Isotopic Ratios: 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O3"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("666"      ,.992901E+00     ,47.984745    ,31     ,31     ,TAGS(48004, 48005, 48006, 48007, 48008)),
         REC("668"      ,3.98194E-03     ,49.988991    ,32     ,32     ,TAGS(50004, 50006)),
         REC("686"      ,1.99097E-03     ,49.988991    ,33     ,33     ,TAGS(50003, 50005)),
         REC("667"      ,7.40475E-04     ,48.988960    ,34     ,34     ,TAGS(49002)),
         REC("676"      ,3.70237E-04     ,48.988960    ,35     ,35     ,TAGS(49001))
         ) ) );


  // N2O
  // Isotopic Ratios: 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("N2O"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("446"      ,.990333E+00    ,44.001062    ,41     ,41     ,TAGS(44004, 44009, 44012)),
         REC("456"      ,3.64093E-03    ,44.998096    ,42     ,42     ,TAGS(45007)),
         REC("546"      ,3.64093E-03    ,44.998096    ,43     ,43     ,TAGS(45008)),
         REC("448"      ,1.98582E-03    ,46.005308    ,44     ,44     ,TAGS(46007)),
         REC("447"      ,3.69280E-04    ,45.005278    ,-1     ,45     ,TAGS() )
         ) ) );

  // CO
  // Isotopic Ratios: 1 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CO"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("26"       , .986544E+00      ,27.994915    ,51     ,51     ,TAGS(28001)),
         REC("36"       , 1.10836E-02      ,28.998270    ,52     ,52     ,TAGS(29001)),
         REC("28"       , 1.97822E-03      ,29.999161    ,53     ,53     ,TAGS(30001)),
         REC("27"       , 3.67867E-04      ,28.999130    ,-1     ,54     ,TAGS(29006)),
         REC("38"       , 2.22250E-05      ,31.002516    ,-1     ,55     ,TAGS() ),
         REC("37"       , 4.13292E-06      ,30.002485    ,-1     ,56     ,TAGS() )
         ) ) );

  // CH4
  // Degrees of freedom: jpl catalogue
  // Isotopic Ratios: 1 1 1
  // Note: - jpl isotopic ratio for tag 17003: 0.00014996848
  //       - CH4 is in official mytran list (6), but does not 
  //         seem to be included for calculation, as given
  //         by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CH4"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("211"      ,.988274E+00     ,16.031300    ,-1     ,61     ,TAGS()),
         REC("311"      ,1.11031E-02     ,17.034655    ,-1     ,62     ,TAGS()),
         REC("212"      ,6.15751E-04     ,17.037475    ,-1     ,63     ,TAGS(17003))
         ) ) );

  // O2
  // Isotopic Ratios: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O2"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("66"       , .995262E+00     ,31.989830    ,71     ,71     ,TAGS(32001, 32002)),
         REC("68"       , 3.99141E-03     ,33.994076    ,72     ,72     ,TAGS(34001)),
         REC("67"       , 7.42235E-04     ,32.994045    ,73     ,73     ,TAGS(33002)),
         REC( "CIAfunCKDMT100"  ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "v0v0CKDMT100"    ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "v1v0CKDMT100"    ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "SelfContStandardType"    ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "SelfContMPM93"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "SelfContPWR93"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "PWR98"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "PWR93"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "PWR88"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "MPM93"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "MPM92"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "MPM89"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "MPM87"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "MPM85"   ,-1.    ,-1.    ,-1     ,-1     ,TAGS())
         ) ) );

  // NO
  // Isotopic Ratios: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("46"       , .993974E+00    ,29.997989    ,81     ,81     ,TAGS(30008)),
         REC("56"       , 3.65431E-03    ,30.995023    ,-1     ,82     ,TAGS() ),
         REC("48"       , 1.99312E-03    ,32.002234    ,-1     ,83     ,TAGS() )
         ) ) );

  // SO2
  // Isotopic Ratios: 1 1 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("SO2"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("626"      ,.945678E+00    ,63.961901    ,91     ,91     ,TAGS(64002, 64005)),
         REC("646"      ,4.19503E-02    ,65.957695    ,92     ,92     ,TAGS(66002)),
         REC("636"      ,0.0074989421   ,65.00        ,93     ,-1     ,TAGS(65001)),
         REC("628"      ,0.0020417379   ,66.00        ,94     ,-1     ,TAGS(66004))  
         ) ) );

  // NO2
  // Isotopic Ratios: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO2"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("646"      ,.991616E+00      ,45.992904    ,101    ,101    ,TAGS(46006))
         ) ) );

  // NH3
  // Isotopic Ratios: 1 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NH3"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("4111"     ,.9958715E+00    ,17.026549    ,111    ,111    ,TAGS(17002, 17004)),
         REC("5111"     ,3.66129E-03    ,18.023583    ,112    ,112    ,TAGS(18002)),
         REC("4112"     ,0.00044792294  ,18.00        ,-1     ,-1     ,TAGS(18004))
         ) ) );

  // HNO3
  // Isotopic Ratios: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HNO3"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("146"      ,0.989110       ,62.995644    ,121    ,121    ,TAGS(63001, 63002, 63003, 63004, 63005, 63006))
         ) ) );

  // OH
  // Isotopic Ratios: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OH"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("61"       ,.997473E+00    ,17.002740    ,131    ,131    ,TAGS(17001)),
         REC("81"       ,2.00014E-03    ,19.006986    ,132    ,132    ,TAGS(19001)),
         REC("62"       ,1.55371E-04    ,18.008915    ,133    ,133    ,TAGS(18001))
         ) ) );

  // HF
  // Isotopic Ratios: 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HF"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("19"       ,0.99984425     ,20.006229    ,141    ,141    ,TAGS(20002)),
         REC("29"       ,0.00014994513  ,21.00        ,-1     ,-1     ,TAGS(21002))
         ) ) );

  // HCl
  // Isotopic Ratios: 1 1 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HCl"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("15"       ,0.757587       ,35.976678    ,151    ,151    ,TAGS(36001)),
         REC("17"       ,0.242257       ,37.973729    ,152    ,152    ,TAGS(38001)),
         REC("25"       ,0.00011324004  ,37.00        ,-1     ,-1     ,TAGS(37001)),
         REC("27"       ,3.6728230E-05  ,39.00        ,-1     ,-1     ,TAGS(39004))
         ) ) );

  // HBr
  // Isotopic Ratios: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HBr"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("19"       ,0.506781        ,79.926160    ,161    ,161    ,TAGS(80001)),
         REC("11"       ,0.493063        ,81.924115    ,162    ,162    ,TAGS(82001))
         ) ) );

  // HI
  // Degrees of freedom: guessed, since it seems to be linear
  // Isotopic Ratios: 1
  // Note: HI is in official mytran list (17), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c

  species_data.push_back
    ( SpeciesRecord
      ( NAME("HI"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("17"       ,0.99984425     ,127.912297   ,-1     ,171    ,TAGS( ))
         ) ) );

  // ClO
  // Isotopic Ratios: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("ClO"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("56"       ,.755908E+00       ,50.963768    ,181    ,181    ,TAGS(51002, 51003)),
         REC("76"       ,.241720E+00       ,52.960819    ,182    ,182    ,TAGS(53002, 53006))
         ) ) );

  // OCS
  // Isotopic Ratios: 1 1 1 1 1
  // Note: OCS-623 is new in Hitran 2000, and it is not the at least
  // abundant isotope. So what actually happend with the following
  // one, is that numbers changed in the hitran edition? I had a
  // look at the two editions, and they actually changed them, stupid
  // idiots. So what used to be 194 (hitran 96) is now 195 (hitran
  // 2000). This messes up the whole concept of reading a catalogue,
  // because the species depends now on the edition of the
  // catalogue.
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OCS"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("622"      ,.937395E+00     ,59.966986    ,191    ,191    ,TAGS(60001)),
         REC("624"      ,4.15828E-02     ,61.962780    ,192    ,192    ,TAGS(62001)),
         REC("632"      ,1.05315E-02     ,60.970341    ,193    ,193    ,TAGS(61001)),
	 REC("623"      ,7.39908E-03     ,60.966371    ,194    ,194    ,TAGS( )),
	 REC("822"      ,1.87967E-03     ,61.971231    ,195    ,195    ,TAGS(62002))
         ) ) );

  // H2CO
  // Isotopic Ratios: 1 1 1 3 3
  // Note: the isotope names differ from hitran convention, since the jpl catalogue has 
  //       isotopes HHCO, HDCO, DDCO.
  //       hitran convention  --  new convention  -- jpl species
  //              126                 1126             H2CO
  //              136                 1136             H2C-13-O
  //              128                 1128             H2CO-18
  //               -                  1226             HDCO
  //               -                  2226             D2CO
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2CO"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("1126"     ,.986237E+00    ,30.010565    ,201    ,201    ,TAGS(30004)),
         REC("1136"     ,1.10802E-02    ,31.013920    ,202    ,202    ,TAGS(31002)),
         REC("1128"     ,1.97761E-03    ,32.014811    ,203    ,203    ,TAGS(32004)),
         REC("1226"     ,0.00029578940  ,31.00        ,-1     ,-1     ,TAGS(31003)),
         REC("2226"     ,2.2181076E-08  ,32.00        ,-1     ,-1     ,TAGS(32006))
         ) ) );

  // HOCl
  // Isotopic Ratios: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HOCl"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("165"      ,.755790E+00       ,51.971593    ,211    ,211    ,TAGS(52006)),
         REC("167"      ,.241683E+00       ,53.968644    ,212    ,212    ,TAGS(54005))
         ) ) );

  // N2
  // Degrees of freedom: guessed, since it seems to be linear
  // Isotopic Ratios: 1
  // Note: N2 is in official mytran list (22), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("N2"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("44"       ,0.9926874      ,28.006147    ,-1     ,221    ,TAGS( )),
         REC( "SelfContMPM93"           ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "SelfContPWR93"           ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "SelfContStandardType"    ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "SelfContBorysow"         ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "CIArotCKDMT100"          ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "CIAfunCKDMT100"          ,-1.    ,-1.    ,-1     ,-1     ,TAGS()),
         REC( "DryContATM01"            ,-1.    ,-1.    ,-1     ,-1     ,TAGS())
         ) ) );

  // HCN
  // Isotopic Ratios: 1 1 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HCN"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("124"      ,.985114E+00      ,27.010899    ,231    ,231    ,TAGS(27001, 27003)),
         REC("134"      ,1.10676E-02      ,28.014254    ,232    ,232    ,TAGS(28002)),
         REC("125"      ,3.62174E-03      ,28.007933    ,233    ,233    ,TAGS(28003)),
         REC("224"      ,0.00014773545  ,28.00        ,-1     ,-1     ,TAGS(28004))
         ) ) );

  // CH3Cl
  // Isotopic Ratios: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CH3Cl"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("215"      ,.748937E+00        ,49.992328    ,241    ,241    ,TAGS(50007)),
         REC("217"      ,.239491E+00        ,51.989379    ,242    ,242    ,TAGS(52009))
         ) ) );

  // H2O2
  // Isotopic Ratios: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O2"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("1661"     ,.994952E+00       ,34.005480    ,251    ,251    ,TAGS(34004))
         ) ) );

  // C2H2
  // Degrees of freedom: guessed, since it seems to be non linear
  // Isotopic Ratios: 1 1
  // Note: C2H2 is in official mytran list (26), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("C2H2"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("1221"     ,.977599E+00       ,26.015650    ,-1     ,261    ,TAGS( )),
         REC("1231"     ,2.19663E-02       ,27.019005    ,-1     ,262    ,TAGS( ))
         ) ) );

  // C2H6
  // Degrees of freedom: guessed, since it seems to be non linear
  // Isotopic Ratios: 1
  // Note: C2H6 is in official mytran list (27), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("C2H6"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("1221"     ,.976990E+00        ,30.046950    ,-1     ,271    ,TAGS( ))
         ) ) );

  // PH3
  // Isotopic Ratios: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("PH3"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("1111"     ,0.99953283     ,33.997238    ,281    ,281    ,TAGS(34003))
         ) ) );

  // COF2
  // Isotopic Ratios: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("COF2"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("269"      ,.986544E+00        ,65.991722    ,291    ,291    ,TAGS(66001))
         ) ) );

  // SF6
  // Degrees of freedom: guessed, since it seems to be non linear
  // Isotopic Ratios: 1
  // Note: SF6 is in official mytran list (30), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("SF6"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("29"       ,0.95018        ,145.962492   ,-1     ,301    ,TAGS( ))
         ) ) );

  // H2S
  // Isotopic Ratios: 1 1 1 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2S"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("121"      ,.949884E+00    ,33.987721    ,311    ,311    ,TAGS(34002)),
         REC("141"      ,4.21369E-02    ,35.983515    ,-1     ,312    ,TAGS( )),
         REC("131"      ,7.49766E-03    ,34.987105    ,-1     ,313    ,TAGS( )),
         REC("122"      ,0.00029991625  ,35.00        ,-1     ,-1     ,TAGS(35001))
         ) ) );

  // HCOOH
  // Isotopic Ratios: 1 3 3 3
  // Note: the isotope names differ from hitran convention, since the jpl catalogue has 
  //       isotopes HCOOH, HC-13-OOH, DCOOH, HCOOD
  //       hitran convention  --  new convention  -- jpl species
  //              126                 1261             HCOOH
  //               -                  1361             HC-13-OOH
  //               -                  2261             DCOOH
  //               -                  1262             HCOOD
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HCOOH"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("1261"     ,0.983898E+00   ,46.005480    ,321    ,321    ,TAGS(46005)),
         REC("1361"     ,0.010913149    ,47.00        ,-1     ,-1     ,TAGS(47002)),
         REC("2261"     ,0.00014755369  ,47.00        ,-1     ,-1     ,TAGS(47003)),
         REC("1262"     ,0.00014755369  ,47.00        ,-1     ,-1     ,TAGS(47004))
         ) ) );

  // HO2
  // Isotopic Ratios: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HO2"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("166"      ,0.995107       ,32.997655    ,331    ,331    ,TAGS(33001))
         ) ) );

  // O
  // Isotopic Ratios: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O"),
        DEGFR(0),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("6"        ,0.997628       ,15.994915    ,341    ,341    ,TAGS(16001))
         ) ) );

  // ClONO2
  // Isotopic Ratios: 1 1
  // Note: ClONO2 in hitran is identical to ClNO3 in jpl (according to Johannes Orphal)
  species_data.push_back
    ( SpeciesRecord
      ( NAME("ClONO2"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("5646"     ,.749570E+00        ,96.956672    ,351    ,351    ,TAGS(97002)),
         REC("7646"     ,.239694E+00        ,98.953723    ,352    ,352    ,TAGS(99001))
         ) ) );

  // NO+
  // Degrees of freedom: guessed, since it seems to be linear
  // Isotopic Ratios: 1
  // Note: NO+ is in official mytran list (36), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO+"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("46"       ,0.993974       ,29.997989    ,-1     ,361    ,TAGS(30011))
         ) ) );

  // OClO
  // Isotopic Ratios: 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OClO"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("656"      ,0.75509223     ,67.00    ,431    ,-1     ,TAGS(67001)),
         REC("676"      ,0.24490632     ,69.00    ,432    ,-1     ,TAGS(69001))
         ) ) );

  // BrO
  // Isotopic Ratios: 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("BrO"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("96"       ,0.50582466     ,95.00    ,401    ,-1     ,TAGS(95001)),
         REC("16"       ,0.49431069     ,97.00    ,402    ,-1     ,TAGS(97001))
         ) ) );

  // H2SO4
  // Isotopic Ratios: 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2SO4"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("126"      ,0.95060479     ,98.00    ,481    ,-1     ,TAGS(98001))
         ) ) );

  // Cl2O2
  // Isotopic Ratios: 2 2
  // Note: refered to as Cl2O2 in mytran catalogue, in jpl cat: ClOOCl
  species_data.push_back
    ( SpeciesRecord
      ( NAME("Cl2O2"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("565"      ,0.57016427     ,102.00   ,491    ,-1     ,TAGS(102001)),
         REC("765"      ,0.36982818     ,104.00   ,492    ,-1     ,TAGS(104001))
         ) ) );

  // HOBr
  // Isotopic Ratios: 1 1
  // Note: latest addtion to Hitran 2000, DEGFR guessed
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HOBr"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("169"      ,.505579E+00    ,95.921076    ,371    ,371    ,TAGS(96001)),
         REC("161"      ,.491894E+00    ,97.919027    ,372    ,372    ,TAGS(98002))
         ) ) );

  // C2H4
  // Isotopic Ratios: 1 1
  // Note: latest addtion to Hitran 2000, DEGFR guessed
  species_data.push_back
    ( SpeciesRecord
      ( NAME("C2H4"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("221"      ,.977294E+00    ,28.031300    ,-1     ,381    ,TAGS( )),
         REC("231"      ,.219595E-01    ,29.034655    ,-1     ,382    ,TAGS( ))
         ) ) );

  // CH3CN
  // Isotopic Ratios: 2 2 2 2 2
  // Note: Isotopic ratio of 1.0 was found in JPL catalogue for main, 
  // isotope, the value given here is determinded by subtracting the other
  // isotopic ratios found in JPL from 1.0
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CH3CN"),
        DEGFR(3),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("211124"   ,0.97366840     ,41.00    ,-1     ,-1     ,TAGS(41001)),
         REC("311124"   ,0.011091748    ,42.00    ,-1     ,-1     ,TAGS(42006)),
         REC("211134"   ,0.011091748    ,42.00    ,-1     ,-1     ,TAGS(42007)),
         REC("211125"   ,0.0036982817   ,42.00    ,-1     ,-1     ,TAGS(42001)),
         REC("211224"   ,0.00044977985  ,42.00    ,-1     ,-1     ,TAGS(42008))
         ) ) );

  // HNC
  // Isotopic Ratios: 2 2 2 2
  // Note: Isotopic ratio of 1.0 was found in JPL catalogue for main, 
  // isotope, the value given here is determinded by subtracting the other
  // isotopic ratios found in JPL from 1.0
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HNC"),
        DEGFR(2),
        ISOTOPES
        (//  Name,      Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC("142"      ,0.98505998     ,27.00    ,-1     ,-1     ,TAGS(27002)),
         REC("143"      ,0.011091748    ,28.00    ,-1     ,-1     ,TAGS(28005)),
         REC("152"      ,0.0036982817   ,28.00    ,-1     ,-1     ,TAGS(28006)),
         REC("242"      ,0.00014996849  ,28.00    ,-1     ,-1     ,TAGS(28007))
         ) ) );


  // You also have to change the entry in the file
  // partition_function_data.cc consistently! 
  species_data.push_back
    ( SpeciesRecord
      ( NAME("liquidcloud"),
        DEGFR(0),
        ISOTOPES
        (//   Name,             Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //                     |               |       |       |       |
         REC( "MPM93"   ,-1.            ,-1.    ,-1     ,-1     ,TAGS())
         ) ) );

  // You also have to change the entry in the file
  // partition_function_data.cc consistently! 
  species_data.push_back
    ( SpeciesRecord
      ( NAME("icecloud"),
        DEGFR(0),
        ISOTOPES
        (//   Name,             Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //                     |               |       |       |       |
         REC( "MPM93"           ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      )
         ) ) );

  // You also have to change the entry in the file
  // partition_function_data.cc consistently! 
  species_data.push_back
    ( SpeciesRecord
      ( NAME("rain"),
        DEGFR(0),
        ISOTOPES
        (//   Name,             Isotopic Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //                     |               |       |       |       |
         REC( "MPM93"           ,-1.    ,-1.    ,-1     ,-1     ,TAGS()      )
         ) ) );

}
