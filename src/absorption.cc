/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

// Physical absorption routines. The absorption workspace methods are
// in file m_abs.cc 

#include "arts.h"
#include "make_array.h"
#include "absorption.h"
#include "math_funcs.h"
#include "messages.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DEGFR(x) x
#define ISOTOPES     make_array<IsotopeRecord>
#define REC	     IsotopeRecord
#define TAGS	     make_array<int>



void define_species_data()
{
  extern ARRAY<SpeciesRecord> species_data;

  // Initialize to zero, just in case:
  species_data.clear();

  /* Here's an empty template record entry:

  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( ""	,		,	,	,	,TAGS() ),
	 REC( ""	,		,	,	,	,TAGS() )
	 ) ) );

  For good looks, keep the commas on the marks!



  Source for entries:
  -------------------
  -------------------


  Molecule Name: Generally Hitran convention, unless otherwise indicate

  Degrees of freedom: Generally taken from the forward_4_96 release,
                      file glob_def.c, unless otherwise indicate

  Isotopes: Generally Hitran convention, e.g., for the main H2O
            isotope: 161, thus the last digit of the sum of the
            neutrons and protons of the individual atoms is
            taken. Nevertheless, molecules like H2CO have only 3
            digits in Hitran, only one for the two H atoms. This leads
            to problems with the jpl catalogue, if this includes DHCO
            and DDCO. For these cases, a different scheme is taken as
            indicated.

  Abundance: Generally Hitran convention, unless otherwise indicate.
             A number code in the header gives the source of the abundancies 
	     for each isotope of each species:

             1: hitran 96 abundance, taken from cd: file
                software/generic/tables_96.txt (default, even if jpl and
                hitran abundancies are available) 

             2: jpl abundance, taken from the documentation coming
                along with the catalogue. latest catalogue version
                extracted 27.07.00, can be found at
                /pool/lookup/jpl/cat7_00/doc/d<tag_nr>.cat 

             3: jpl abundance is multiplied with the maximum abundance
                of this species found in hitran. is only performed when
                abundancies of 1 were found in the jpl catalogue.

  Mass: Rounded atomic weight

  MY-tags: Extracted from file glob_def.c

  HI-tags: From hitran 96, software/generic/tables_96.txt

  JPL-tags: Collected in file 
            /pool/lookup/jpl/cat7_00/abundancies/tag_species.jpl, taken from
	    last issue of jpl catalogue (7/00)


  Some more information can be found at /pool/lookup/jpl/cat7_00/abundancies,
  where the idl script that reads/converts the abundancies is located.

  08.08.00 AvE

  */


  // H2O
  // Abundancies: 1 1 1 1 3 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( "161"	,0.99731702	,18.	,11	,11	,TAGS(18003, 18005) ),
	 REC( "181"	,0.00199983	,20.	,12	,12	,TAGS(20003) ),
	 REC( "171"	,0.00037200	,19.	,13	,13	,TAGS(19003) ),
	 REC( "162"	,0.00031069	,19.	,14	,14	,TAGS(19002) ),
	 REC( "182"	,6.1070746E-07	,21.	,-1	,-1	,TAGS(21001) ),
	 REC( "262"	,2.2430204E-08	,20.	,-1	,-1	,TAGS(20001) )
	 ) ) );

  // CO2 
  // (missing mainly in JPL, latest version (7/00) includes some isotopes)
  // Degrees of freedom from Schanda:`Physical Fundamentals of Remote Sensing'
  // Abundancies: 1 1 1 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CO2"),
	DEGFR(2),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( "626"	,0.98420	,44.	,21	,21	,TAGS() ),
	 REC( "636"	,0.01106	,45.	,22	,22	,TAGS() ),
	 REC( "628"	,0.0039471	,46.	,23	,23	,TAGS(46013) ),
	 REC( "627"	,0.000734	,45.	,24	,24	,TAGS(45012) ),
	 REC( "638"	,0.00004434	,47.	,25	,25	,TAGS() ),
	 REC( "637"	,0.00000825	,46.	,26	,26	,TAGS() ),
	 REC( "828"	,0.0000039573	,48.	,27	,27	,TAGS() ),
	 REC( "728"	,0.00000147	,47.	,28	,28	,TAGS() )
	 ) ) );
  

  // O3
  // Abundancies: 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("666"	,0.992901	,48.	,31	,31	,TAGS(48004, 48005, 48006, 48007, 48008)),
	 REC("668"	,0.00398194	,50.	,32	,32	,TAGS(50004, 50006)),
	 REC("686"	,0.00199097	,50.	,33	,33	,TAGS(50003, 50005)),
	 REC("667"	,0.000740	,49.	,34	,34	,TAGS(49002)),
	 REC("676"	,0.000370	,49.	,35	,35	,TAGS(49001))
	 ) ) );


  // N2O
  // Abundancies: 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("N2O"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("446"	,0.990333	,44.	,41	,41	,TAGS(44004, 44009, 44012)),
	 REC("456"	,0.0036409	,45.	,42	,42	,TAGS(45007)),
	 REC("546"	,0.0036409	,45.	,43	,43	,TAGS(45008)),
	 REC("448"	,0.00198582	,46.	,44	,44	,TAGS(46007)),
	 REC("447"	,0.000369	,45.	,-1	,45	,TAGS() )
	 ) ) );

  // CO
  // Abundancies: 1 1 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CO"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("26"	,0.98654	,28.	,51	,51	,TAGS(28001)),
	 REC("36"	,0.01108	,29.	,52	,52	,TAGS(29001)),
	 REC("28"	,0.0019782	,30.	,53	,53	,TAGS(30001)),
	 REC("27"	,0.000368	,29.	,-1	,54	,TAGS(29006)),
	 REC("38"	,0.00002222	,31.	,-1	,55	,TAGS() ),
	 REC("37"	,0.00000413	,30.	,-1	,56	,TAGS() )
	 ) ) );

  // CH4
  // Degrees of freedom: jpl catalogue
  // Abundancies: 1 1 1
  // Note: - jpl abundancies for tag 17003: 0.00014996848
  //       - CH4 is in official mytran list (6), but does not 
  //         seem to be included for calculation, as given
  //         by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CH4"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("211"	,0.98827	,16.	,-1	,61	,TAGS()),
	 REC("311"	,0.01110	,17.	,-1	,62	,TAGS()),
	 REC("212"	,0.00061575	,17.	,-1	,63	,TAGS(17003))
	 ) ) );

  // O2
  // Abundancies: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O2"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("66"	,0.995262	,32.	,71	,71	,TAGS(32001, 32002, 32005)),
	 REC("68"	,0.00399141	,34.	,72	,72	,TAGS(34001)),
	 REC("67"	,0.000742	,33.	,73	,73	,TAGS(33002))
	 ) ) );

  // NO
  // Abundancies: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("46"	,0.993974	,30.	,81	,81	,TAGS(30008)),
	 REC("56"	,0.0036543	,31.	,-1	,82	,TAGS() ),
	 REC("48"	,0.00199312	,32.	,-1	,83	,TAGS() )
	 ) ) );

  // SO2
  // Abundancies: 1 1 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("SO2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("626"	,0.94568	,64.	,91	,91	,TAGS(64002, 64005)),
	 REC("646"	,0.04195	,66.	,-1	,92	,TAGS(66002)),
	 REC("636"	,0.0074989421	,65.	,-1	,-1	,TAGS(65001)),
	 REC("628"	,0.0020417379	,66.	,-1	,-1	,TAGS(66004))
	 ) ) );

  // NO2
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("646"	,0.991616	,46.	,101	,101	,TAGS(46006))
	 ) ) );

  // NH3
  // Abundancies: 1 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NH3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("4111"	,0.9958715	,17.	,111	,111	,TAGS(17002, 17004)),
	 REC("5111"	,0.0036613	,18.	,111	,111	,TAGS(18002)),
	 REC("4112"	,0.00044792294	,18.	,-1	,-1	,TAGS(18004))
	 ) ) );

  // HNO3
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HNO3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("146"	,0.989110	,63.	,121	,121	,TAGS(63001, 63002, 63003, 63004, 63005, 63006))
	 ) ) );

  // OH
  // Abundancies: 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OH"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("61"	,0.997473	,17.	,131	,131	,TAGS(17001)),
	 REC("81"	,0.00200014	,19.	,132	,132	,TAGS(19001)),
	 REC("62"	,0.00015537	,18.	,133	,133	,TAGS(18001))
	 ) ) );

  // HF
  // Abundancies: 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HF"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("19"	,0.99984425	,20.	,141	,141	,TAGS(20002)),
	 REC("29"	,0.00014994513	,21.	,-1	,-1	,TAGS(21002))
	 ) ) );

  // HCl
  // Abundancies: 1 1 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HCl"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("15"	,0.757587	,36.	,151	,151	,TAGS(36001)),
	 REC("17"	,0.242257	,38.	,152	,152	,TAGS(38001)),
	 REC("25"	,0.00011324004	,37.	,-1	,-1	,TAGS(37001)),
	 REC("27"	,3.6728230E-05	,39.	,-1	,-1	,TAGS(39004))
	 ) ) );

  // HBr
  // Abundancies: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HBr"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("19"	,0.50678	,80.	,161	,161	,TAGS(80001)),
	 REC("11"	,0.49306	,82.	,162	,162	,TAGS(82001))
	 ) ) );

  // HI
  // Degrees of freedom: guessed, since it seems to be linear
  // Abundancies: 1
  // Note: HI is in official mytran list (17), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c

  species_data.push_back
    ( SpeciesRecord
      ( NAME("HI"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("17"	,0.99984425	,128.	,-1	,171	,TAGS( ))
	 ) ) );

  // ClO
  // Abundancies: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("ClO"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("56"	,0.75591	,51.	,181	,181	,TAGS(51002, 51003)),
	 REC("76"	,0.24172	,53.	,182	,182	,TAGS(53002, 53006))
	 ) ) );

  // OCS
  // Abundancies: 1 1 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OCS"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("622"	,0.93739	,60.	,191	,191	,TAGS(60001)),
	 REC("624"	,0.04158	,62.	,192	,192	,TAGS(62001)),
	 REC("632"	,0.01053	,61.	,193	,193	,TAGS(61001)),
	 REC("822"	,0.0018797	,62.	,194	,194	,TAGS(62002))
	 ) ) );

  // H2CO
  // Abundancies: 1 1 1 3 3
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
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1126"	,0.98624	,30.	,201	,201	,TAGS(30004)),
	 REC("1136"	,0.01108	,31.	,202	,202	,TAGS(31002)),
	 REC("1128"	,0.0019776	,32.	,203	,203	,TAGS(32004)),
	 REC("1226"	,0.00029578940	,31.	,-1	,-1	,TAGS(31003)),
	 REC("2226"	,2.2181076E-08	,32.	,-1	,-1	,TAGS(32006))
	 ) ) );

  // HOCl
  // Abundancies: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HOCl"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("165"	,0.75579	,52.	,211	,211	,TAGS(52006)),
	 REC("167"	,0.24168	,54.	,212	,212	,TAGS(54005))
	 ) ) );

  // N2
  // Degrees of freedom: guessed, since it seems to be linear
  // Abundancies: 1
  // Note: N2 is in official mytran list (22), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("N2"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("44"	,0.9926874	,28.	,-1	,221	,TAGS( ))
	 ) ) );

  // HCN
  // Abundancies: 1 1 1 3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HCN"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("124"	,0.98511	,27.	,231	,231	,TAGS(27001, 27003)),
	 REC("134"	,0.01107	,28.	,232	,232	,TAGS(28002)),
	 REC("125"	,0.0036217	,28.	,233	,233	,TAGS(28003)),
	 REC("224"	,0.00014773545	,28.	,-1	,-1	,TAGS(28004))
	 ) ) );

  // CH3Cl
  // Abundancies: 1 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CH3Cl"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("215"	,0.74894	,50.	,241	,241	,TAGS(50007)),
	 REC("217"	,0.23949	,52.	,242	,242	,TAGS(52009))
	 ) ) );

  // H2O2
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1661"	,0.994952	,34.	,251	,251	,TAGS(34004))
	 ) ) );

  // C2H2
  // Degrees of freedom: guessed, since it seems to be non linear
  // Abundancies: 1 1
  // Note: C2H2 is in official mytran list (26), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("C2H2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1221"	,0.97760	,26.	,-1	,261	,TAGS( )),
	 REC("1231"	,0.02197	,27.	,-1	,262	,TAGS( ))
	 ) ) );

  // C2H6
  // Degrees of freedom: guessed, since it seems to be non linear
  // Abundancies: 1
  // Note: C2H6 is in official mytran list (27), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("C2H6"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1221"	,0.97699	,30.	,-1	,271	,TAGS( ))
	 ) ) );

  // PH3
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("PH3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1111"	,0.99953283	,34.	,281	,281	,TAGS(34003))
	 ) ) );

  // COF2
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("COF2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("269"	,0.98654	,66.	,291	,291	,TAGS(66001))
	 ) ) );

  // SF6
  // Degrees of freedom: guessed, since it seems to be non linear
  // Abundancies: 1
  // Note: SF6 is in official mytran list (30), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("SF6"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("29"	,0.95018	,146.	,-1	,301	,TAGS( ))
	 ) ) );

  // H2S
  // Abundancies: 1 1 1 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2S"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("121"	,0.94988	,34.	,311	,311	,TAGS(34002)),
	 REC("141"	,0.04214	,36.	,-1	,312	,TAGS( )),
	 REC("131"	,0.007498	,35.	,-1	,313	,TAGS( )),
	 REC("122"	,0.00029991625	,35.	,-1	,-1	,TAGS(35001))
	 ) ) );

  // HCOOH
  // Abundancies: 1 3 3 3
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
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("1261"	,0.983898	,46.	,321	,321	,TAGS(46005)),
	 REC("1361"	,0.010913149	,47.	,-1	,-1	,TAGS(47002)),
	 REC("2261"	,0.00014755369	,47.	,-1	,-1	,TAGS(47003)),
	 REC("1262"	,0.00014755369	,47.	,-1	,-1	,TAGS(47004))
	 ) ) );

  // HO2
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("HO2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("166"	,0.995107	,33.	,331	,331	,TAGS(33001))
	 ) ) );

  // O
  // Abundancies: 1
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O"),
	DEGFR(0),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("6"	,0.997628	,16.	,341	,341	,TAGS(16001))
	 ) ) );

  // ClONO2
  // Abundancies: 1 1
  // Note: ClONO2 in hitran is identical to ClNO3 in jpl (according to Johannes Orphal)
  species_data.push_back
    ( SpeciesRecord
      ( NAME("ClONO2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("5646"	,0.74957	,97.	,351	,351	,TAGS(97002)),
	 REC("7646"	,0.23970	,99.	,352	,352	,TAGS(99001))
	 ) ) );

  // NO+
  // Degrees of freedom: guessed, since it seems to be linear
  // Abundancies: 1
  // Note: NO+ is in official mytran list (36), but does not 
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back
    ( SpeciesRecord
      ( NAME("NO+"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("46"	,0.993974	,30.	,-1	,361	,TAGS(30011))
	 ) ) );

  // OClO
  // Abundancies: 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("OClO"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("656"	,0.75509223	,67.	,431	,-1	,TAGS(67001)),
	 REC("676"	,0.24490632	,69.	,432	,-1	,TAGS(69001))
	 ) ) );

  // BrO
  // Abundancies: 2 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("BrO"),
	DEGFR(2),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("96"	,0.50582466	,95.	,461	,-1	,TAGS(95001)),
	 REC("16"	,0.49431069	,97.	,462	,-1	,TAGS(97001))
	 ) ) );

  // H2SO4
  // Abundancies: 2
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2SO4"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("126"	,0.95060479	,98.	,481	,-1	,TAGS(98001))
	 ) ) );

  // Cl2O2
  // Abundancies: 2 2
  // Note: refered to as Cl2O2 in mytran catalogue, in jpl cat: ClOOCl
  species_data.push_back
    ( SpeciesRecord
      ( NAME("Cl2O2"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("565"	,0.57016427	,102.	,491	,-1	,TAGS(102001)),
	 REC("765"	,0.36982818	,104.	,492	,-1	,TAGS(104001))
	 ) ) );


}


void define_species_map()
{
  extern const ARRAY<SpeciesRecord> species_data;
  extern std::map<string, size_t> SpeciesMap;

  for ( size_t i=0 ; i<species_data.size() ; ++i)
    {
      SpeciesMap[species_data[i].Name()] = i;
    }
}


ostream& operator << (ostream& os, const LineRecord& lr)
{
  // Determine the precision, depending on whether Numeric is double
  // or float:  
  int precision;
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
     << " " << lr.Naux  ();

  for ( size_t i=0; i<lr.Naux(); ++i )
    os << " " << lr.Aux()[i];

  return os;
}


/** Extract something from a HITRAN line. This is just a small helper
    function to safe some typing. 

    @param x Output.    What was extracted from the beginning of the line.
    @param line Output. What was extracted is also cut away from line.
    @param n            The width of the stuff to extract.

    @author Stefan Buehler */
template<class T>
void extract(T&      x,
	     string& line,
	     size_t  n)
{
  // This will contain the short substring with the item to extract.
  // Make it a string stream, for easy parsing,
  // extracting substring of width n from line:
  istringstream item( line.substr(0,n) );

//  cout << "line = '" << line << "'\n";
//   cout << "line.substr(0,n) = " << line.substr(0,n) << endl;
//   cout << "item = " << item.str() << endl;

  // Shorten line by n:
  line.erase(0,n);
//  cout << "line = " << line << endl;

  // Convert with the aid of string stream item:
  item >> x;
}


bool LineRecord::ReadFromHitranStream(istream& is)
{
  // Global species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotope lists. This is based on the assumption that there are a
  // lot more species than isotopes. Anyway, we also add 100, because
  // at the moment we have only 3 species and more than 3 isotopes.
  const size_t missing = species_data.size() + 100;

  // We need a species index sorted by HITRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<HITRAN tag>]. The value of
  // missing means that we don't have this species.
  //
  // Allow for up to 100 species in HITRAN in the future.
  static ARRAY< size_t >        hspec(100,missing);	

  // This is  an array of arrays for each hitran tag. It contains the
  // ARTS indices of the HITRAN isotopes. 
  static ARRAY< ARRAY<size_t> > hiso(100);

  // Remeber if this stuff has already been initialized:
  static bool hinit = false;

  // Remember, about which missing species we have already issued a
  // warning: 
  static ARRAY<size_t> warned_missing;

  if ( !hinit )
    {
      for ( size_t i=0; i<species_data.size(); ++i )
	{
	  const SpeciesRecord& sr = species_data[i];

	  // The HITRAN tags are stored as species plus isotope tags
	  // (MO and ISO)
	  // in the Isotope() part of the species record.
	  // We can extract the MO part from any of the isotope tags,
	  // so we use the first one. We do this by taking an integer
	  // division by 10.
	  size_t mo = sr.Isotope()[0].HitranTag() / 10;
//	  cout << "mo = " << mo << endl;
	  hspec[mo] = i; 
	  
	  // Get a nicer to handle array of HITRAN iso tags:
	  size_t n_iso = sr.Isotope().size();
	  ARRAY<int> iso_tags;
	  iso_tags.resize(n_iso);
	  for ( size_t j=0; j<n_iso; ++j )
	    {
	      iso_tags[j] = sr.Isotope()[j].HitranTag();
	    }

	  // Reserve elements for the isotope tags. How much do we
	  // need? This depends on the largest HITRAN tag that we know
	  // about!
	  // Also initialize the tags to missing.
// 	  cout << "iso_tags = " << iso_tags << endl;
// 	  cout << "static_cast<size_t>(max(iso_tags))%10 + 1 = "
// 	       << static_cast<size_t>(max(iso_tags))%10 + 1 << endl;
	  hiso[mo].resize( static_cast<size_t>(max(iso_tags))%10 + 1,
			   missing );
//	  cout << "hiso[mo].size() = " << hiso[mo].size() << endl;

	  // Set the isotope tags:
	  for ( size_t j=0; j<n_iso; ++j )
	    {
	      if ( 0 < iso_tags[j] )				  // ignore -1 elements
		{
		  // To get the iso tags from HitranTag() we also have to take
		  // modulo 10 to get rid of mo.
//		  cout << "iso_tags[j] % 10 = " << iso_tags[j] % 10 << endl;
		  hiso[mo][iso_tags[j] % 10] = j;
		}
	    }
	}
      
//      cout << "hiso = " << hiso << endl << "***********" << endl;


      // Print the generated data structures (for debugging):
      out3 << "  HITRAN index table:\n";
      for ( size_t i=0; i<hspec.size(); ++i )
	{
	  if ( missing != hspec[i] )
	    {
	      // The explicit conversion of Name to a c-string is
	      // necessary, because setw does not work correctly for
	      // stl strings.
	      out3 << "  mo = " << i << "   Species = "
		   << setw(10) << setiosflags(ios::left)
		   << species_data[hspec[i]].Name().c_str()
		   << "iso = ";
	      for ( size_t j=1; j<hiso[i].size(); ++j )
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
  string line;

  // The first item is the molecule number:
  size_t mo;

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
	  // See if we know this species. At the moment, because I
	  // only implemented 3 species, we will only give a warning
	  // if a species is unknown, and not an error. FIXME: This should
	  // probably be changed in the future.
	  if ( missing != hspec[mo] )	    comment = false ;
	  else
	    {
	      // See if this is already in warned_missing, use
	      // std::count for that:
	      if ( 0 == std::count(warned_missing.begin(),
				   warned_missing.end(),
				   mo) )
		{
		  out0 << "Warning: HITRAN mo = " << mo << " is not "
		       << "known to ARTS.\n";
		  warned_missing.push_back(mo);
		}
	    }
	}
    }

  // Ok, we seem to have a valid species here.

  // Set mspecies from my cool index table:
  mspecies = hspec[mo];

  // Extract isotope:
  size_t iso;				  
  extract(iso,line,1);
  //  cout << "iso = " << iso << endl;


  // Set misotope from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases. 
  misotope = missing;
  if ( iso < hiso[mo].size() )
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
    // HITRAN intensity in cm-1/(molec * cm-2) at 296 Kelvin.
    // Includes isotpic ratio!
    Numeric s;
    // Conversion from HITRAN intensities to JPL intensities 
    // (nm^2 MHz). This factor is given in the JPL documentation. 
    const Numeric hi2jpl = 2.99792458e18;
    // Because we want SI units m and Hz instead of nm and MHz, we
    // need to shift the decimals a bit.
    // More importantly, we have to take out the isotope ratio, which
    // is included in the HITRAN intensity.
    const Numeric hi2arts = hi2jpl * 1e-12 / IsotopeData().Abundance();

    // Extract HITRAN intensity:
    extract(s,line,10);

    // ARTS intensity in m^2/Hz:
    mi0 = s * hi2arts;
//    cout << "mi0 = " << mi0 << endl;
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

//    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }


  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // The same unit as in ARTS, we can extract directly!
    extract(melow,line,10);
//    cout << "mf, melow = " << mf << ", " << setprecision(20) << melow << endl;
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
    // HITRAN value in cm^-1
    // (I guess they really mean cm^-1 / atm. Anyway, that's what I'll
    // assume. So the conversion goes exactly as for the broadening
    // parameters. 
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

  // Reference temperature for AGAM and SGAM in K.
  // (This is also fix for HITRAN)
  mtgam = 296.0;


  // That's it!
  return false;
}

OneTag::OneTag(string def) 
{
  // Species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;
  // The species map. This is used to find the species id.
  extern std::map<string, size_t> SpeciesMap;
  // Name of species and isotope (aux variables):
  string name, isoname;
  // Aux index:
  size_t n;


  // Set frequency limits to default values (no limits):
  mlf = -1;
  muf = -1;

  // We cannot set a default value for the isotope, because the
  // default should be `ALL' and the value for `ALL' depends on the
  // species. 
    

  // Extract the species name:
  n    = def.find('-');    // find the '-'
  if (n < def.size() )
    {
      name = def.substr(0,n);      // Extract before '-'
      def.erase(0,n+1);		    // Remove from def
    }
  else
    {
      // n==def.size means that def does not contain a '-'. In that
      // case we assume that it contains just a species name and
      // nothing else
      name = def;
      def  = "";
    }


//  cout << "name / def = " << name << " / " << def << endl;

  // Look for species name in species map:
  map<string, size_t>::const_iterator mi = SpeciesMap.find(name);
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

  if ( 0 == def.size() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants all isotopes and no frequency limits.
      // Frequency defaults are already set. Set isotope defaults:
      misotope = spr.Isotope().size();
      // This means all isotopes.

      return;
    }
    
  // Extract the isotope name:
  n    = def.find('-');    // find the '-'
  if (n < def.size() )
    {
      isoname = def.substr(0,n);          // Extract before '-'
      def.erase(0,n+1);		    // Remove from def
    }
  else
    {
      // n==def.size means that def does not contain a '-'. In that
      // case we assume that it contains just the isotope name and
      // nothing else.
      isoname = def;
      def  = "";
    }
    
  // Check for joker:
  if ( "*" == isoname )
    {
      // The user wants all isotopes. Set this accordingly:
      misotope = spr.Isotope().size();
    }
  else
    {
      // Make an array containing the isotope names:
      ARRAY<string> ins;
      for ( size_t i=0; i<spr.Isotope().size(); ++i )
	ins.push_back( spr.Isotope()[i].Name() );

      // Use the find algorithm from the STL to find the isotope ID. It
      // returns an iterator, so to get the index we take the
      // difference to the begin() iterator.
      misotope = find( ins.begin(),
		       ins.end(),
		       isoname ) - ins.begin();

      // Check if we found a matching isotope:
      if ( misotope >= ins.size() ) 
	{
	  ostringstream os;
	  os << "Isotope " << isoname << " is not a valid isotope for "
	     << "species " << name << ".\n"
	     << "Valid isotopes are:";
	  for ( size_t i=0; i<ins.size(); ++i )
	    os << " " << ins[i];
	  throw runtime_error(os.str());
	}
    }

  if ( 0 == def.size() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants no frequency limits.  Frequency defaults are
      // already set, so we can return directly.

      return;
    }


  // Look for the two frequency limits:
  
  // Extract first frequency
  n    = def.find('-');    // find the '-'
  if (n < def.size() )
    {
      // Frequency as a string:
      string fname;
      fname = def.substr(0,n);              // Extract before '-'
      def.erase(0,n+1);		      // Remove from def

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
      // n==def.size means that def does not contain a '-'. In this
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


string OneTag::Name() const 
{
  // Species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;
  // A reference to the relevant record of the species data:
  const  SpeciesRecord& spr = species_data[mspecies];
  // For return value:
  ostringstream os;

  // First the species name:
  os << spr.Name() << "-";

  // Now the isotope. Can be a single isotope or ALL.
  if ( misotope == spr.Isotope().size() )
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
  int precision;
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


void write_lines_to_stream(ostream& os,
			   const ARRAYofLineRecord& lines)
{
  for ( size_t i=0; i<lines.size(); ++i )
    {
      os << lines[i] << '\n';
    }
}

/** The Lorentz line shape. This is a quick and dirty implementation.

    @param F     Output: The shape function.
    @param f0    Line center frequency.
    @param f_abs The frequency grid.
    @param gamma The pressure broadening parameter.

    @author Stefan Buehler 16.06.2000 */
void line_shape_lorentz(VECTOR&       ls,
			Numeric	      f0,
			Numeric       gamma,
			const VECTOR& f_abs)
{
  // FIXME: Maybe try if call by reference is faster for f0 and gamma?

  // 1/PI:
  extern const Numeric PI;
  static const Numeric invPI = 1. / PI;

  assert( ls.dim() == f_abs.dim() );

  for ( size_t i=0; i<f_abs.dim(); ++i )
    {
      ls[i] = invPI * gamma / ( pow( f_abs[i]-f0, 2) + pow(gamma,2) );
    }
}

void abs_species( MATRIX&                  abs,
		  const VECTOR&  	   f_abs,
		  const VECTOR&  	   p_abs,
		  const VECTOR&  	   t_abs,           
		  const VECTOR&            vmr,
		  const ARRAYofLineRecord& lines )
{
  // Define the vector for the line shape function here, so that we
  // don't need so many free store allocations.
  VECTOR ls(f_abs.dim());

  // Check that p_abs, t_abs, and vmr all have the same
  // dimension. This could be a user error, so we throw a
  // runtime_error. 

  if ( t_abs.dim() != p_abs.dim() )
    {
      ostringstream os;
      os << "Variable t_abs must have the same dimension as p_abs.\n"
	 << "t_abs.dim() = " << t_abs.dim() << '\n'
	 << "p_abs.dim() = " << p_abs.dim();
      throw runtime_error(os.str());
    }

  if ( vmr.dim() != p_abs.dim() )
    {
      ostringstream os;
      os << "Variable vmr must have the same dimension as p_abs.\n"
	 << "vmr.dim() = " << vmr.dim() << '\n'
	 << "p_abs.dim() = " << p_abs.dim();
      throw runtime_error(os.str());
    }

  // Check that the dimension of abs is indeed [f_abs.dim(),
  // p_abs.dim()]:
  if ( abs.dim(1) != f_abs.dim() || abs.dim(2) != p_abs.dim() )
    {
      ostringstream os;
      os << "Variable abs must have dimensions [f_abs.dim(),p_abs.dim()].\n"
	 << "[abs.dim(1),abs.dim(2)] = [" << abs.dim(1)
	 << ", " << abs.dim(2) << "]\n"
	 << "f_abs.dim() = " << f_abs.dim() << '\n'
	 << "p_abs.dim() = " << p_abs.dim();
      throw runtime_error(os.str());
    }

  // Loop all pressures:
  for ( size_t i=0; i<p_abs.dim(); ++i )
  {
    out3 << "  p = " << p_abs[i] << " Pa\n";
    // Loop all lines:
    for ( size_t l=0; l<lines.dim(); ++l )
      {
	// 1. Convert intensity to units of [ m^-1 * Hz / Pa ]
	Numeric intensity = lines[l].I0();
	// FIXME: Move the constants here to constants.cc?
	  {
	    // I have no time to find out about the correct intensity
	    // conversion. So, as a quick hack, we will convert to
	    // MYTRAN intensiy (HITRAN, except isotopic ratio is taken
	    // out, and then use the old conversion from abs_my.c. See
	    // LineRecord::ReadFromHitranStream for the below two
	    // conversions. 

	    // Conversion from HITRAN intensities to JPL intensities 
	    // (nm^2 MHz). This factor is given in the JPL documentation. 
	    const Numeric hi2jpl = 2.99792458e18;
	    // Because we want SI units m and Hz instead of nm and MHz, we
	    // need to shift the decimals a bit.
	    // More importantly, we have to take out the isotope ratio, which
	    // is included in the HITRAN intensity.
	    const Numeric hi2arts = hi2jpl * 1e-12;

	    intensity = intensity / hi2arts;

	    // intensity is now MYTRAN intensity.

	    // From abs_my.c:
	    const Numeric T_0_C = 273.15;  	       /* temp. of 0 Celsius in [K]  */
	    const Numeric p_0   = 101300.25; 	       /* standard p in [Pa]        */
	    const Numeric n_0   = 2.686763E19;         /* Loschmidt constant [cm^-3] */
	    const Numeric c     = 29.9792458;          /*  	     [ GHz / cm^-1 ] */

	    intensity = intensity * n_0 * T_0_C / (p_0*t_abs[i]);
	    /* n = n0*T0/p0 * p/T, ideal gas law
	       ==> [intensity] = [cm^2 / Pa] */

	    intensity = intensity * c; 			/* [cm^-1 * GHz / Pa] */
	    intensity = intensity * 1E5; 		/* [km^-1 * GHz / Pa] */
	    intensity = intensity * 1e6; 		/* [m^-1 * Hz / Pa] */
	  }

	// FIXME: Are the units of the intensity correct?
	
	// FIXME: We should calculate the temperature dependence of the
	// intensities here.

	// Scale with isotopic ratio. FIXME: This is inefficient. The
	// scaled intensities could be stored in lines.
	intensity = intensity * lines[l].IsotopeData().Abundance();

	// 2. Get pressure broadened line width:
	// (Agam is in Hz/Pa, p_abs is in Pa, f_abs is in Hz,
	// gamma is in Hz/Pa)
	// FIXME: This is inefficient, the scaling to Hz could be
	// stored in lines.

	const Numeric theta = lines[l].Tgam() / t_abs[i];

	const Numeric p_partial = p_abs[i] * vmr[i];

	Numeric gamma
	  = lines[l].Agam() * pow(theta, lines[l].Nair())  * (p_abs[i] - p_partial)
	  + lines[l].Sgam() * pow(theta, lines[l].Nself()) * p_partial;

	// Ignore Doppler broadening for now. FIXME: Add this.

	// Get line shape:
	line_shape_lorentz(ls,
			   lines[l].F(),
			   gamma,
			   f_abs);

	// Add line to abs:
	for ( size_t j=0; j<f_abs.dim(); ++j )
	  {
	    abs[j][i] = abs[j][i]
	      + p_partial * intensity * ls[j];
	  }
      }
  }
}
