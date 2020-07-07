/* Copyright (C) 2000-2012
   Stefan Buehler <sbuehler@ltu.se>,
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

  This is the file from arts-1-0, back-ported to arts-1-1.
  NOTE that adding new data here (new species or new isotopologues) requires
  equivalent additions in partition_function_data.cc. Also, equivalent data
  should be added to the isotopratio.xml for other planets residing in the
  planets/ folders of arts-xml-data should be done.

  \author Stefan Buehler, Axel von Engeln
  \date 2000-08-10 
*/

#include "absorption.h"
#include "arts.h"

/*! The lookup information for all the different species. */
namespace global_data {
Array<SpeciesRecord> species_data;
std::map<String, Index> SpeciesMap;
}  // namespace global_data

using global_data::species_data;

/*! \name Some #defines for better readability */
//@{
#define NAME(x) x
#define DEGFR(x) x
#define ISOTOPOLOGUES(...) \
  { __VA_ARGS__ }
#define REC IsotopologueRecord
#define TAGS(...) \
  { __VA_ARGS__ }
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

  <dt> Isotopologues: 
  <dd> Generally Hitran convention, e.g., for the main H2O
       isotopologue: 161, thus the last digit of the sum of the
       neutrons and protons of the individual atoms is
       taken. Nevertheless, molecules like H2CO have only 3
       digits in Hitran, only one for the two H atoms. This leads
       to problems with the jpl catalogue, if this includes DHCO
       and DDCO. For these cases, a different scheme is taken as
       indicated.

  <dt> Isotopologue Ratio: 
  <dd> Built-in isotoplogue ratios for Earth. Generally Hitran convention, unless
       otherwise indicated. A number code in the header gives the source of the
       isotopologue ratios for each isotopologue of each species:
       
       <table>
       <tr>
       <td> 1: <td> HITRAN isotopologue ratio (varying editions; usually the
                    recent one at time of addition)
                    (DEFAULT, even if jpl and hitran isotopologue ratios are 
                    available).
       <tr>
       <td> 2: <td> jpl isotopologue ratio, taken from the documentation coming
                    along with the catalogue (see d\<tag_nr\>.cat at
                    /storage3/data/pool/lookup2/jpl/cat7_00/doc/ or
                    http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/, where
                    log10(isotop.ratio) is documented as Isotope Corr.
       <tr>
       <td> 3: <td> jpl isotopologue ratio is multiplied with the maximum isotopologue
                    ratio of this species found in hitran, or (if species does
                    not exist in hitran) scaled with the sum of the relative
                    abundances of all the jpl isotopologues of this species. is
                    only performed when jpl main isotopologue ratio is given as 1.
       <tr>
       <td> 4: <td> isotopologue ratio calculated from (atomic) isotopic
                    abundances. Origin of isotopic ratio documented in individual
                    species headers.
  	 	 <tr> 
	 	   <td> 0: <td> Wild (but educated) guess (e.g. for inert gases).
       </table>

      On run-time, isotopologue ratios can be initialized from those built-in
      values or can be read in from file (replacing those given here).
      Pre-prepared xml files with isotopologue ratios of selected planets are available
      from arts-xml-data package.

  <dt> Mass: 
  <!-- <dd> Rounded atomic weight. (SAB: From looking at the table in forward_4_96,
       glob_def.c, the difference between the actual mass and the mass
       simply estimated from the Atom number is only 0.001. This seems
       not worth the trouble. Anyway, the field for the mass is
       there. Should anybody feel like adding the true numbers, just go
       ahead.) -->
  <dd> mass as given in file MOLPARAM.TXT of HITRAN2000. 
       If tag is not present in HITRAN, the rounded value is stated. 
       Note that the mass is given in units of [g/mol].

  <dt >MY-tags:
  <dd> Extracted from file glob_def.c.

  <dt> HI-tags:
  <dd> According to HITRAN's molparams.txt (order of isotopologues determines
       isotopologue number) found in
       ftp://cfa-ftp.harvard.edu/pub/HITRAN-XXXX/Global_Data/, where XXXX stands
       for the respective HITRAN version.

  <dt> JPL-tags:
  <dd> From JPL species list
       documented on http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html.
  </dl>

  Some more information can be found at svn:arts/branches/ARTS-1-0/misc/abundances/,
  where the idl script that reads/converts the isotopologue ratios is located.

  \author Axel von Engeln, Jana Mendrok  
  \date   2000-08-08 
*/

// prototyping
void define_basic_species_data();
extern void define_partition_species_data();

void define_species_data() {
  define_basic_species_data();
  define_partition_species_data();
}

void define_basic_species_data() {
  // Initialize to zero, just in case:
  species_data.resize(0);

  /* Here's an empty template record entry:

  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
        DEGFR(3),
        ISOTOPOLOGUES
        (//   Name,     Isotopologue Ratio, Mass,   MY-tag, HI-tag, JPL-tag
         //             |               |       |       |       |
         REC( ""        ,               ,       ,       ,       ,TAGS() ),
         REC( ""        ,               ,       ,       ,       ,TAGS() )
         ) ) );

  For good looks, keep the commas on the marks!

  */

  // H2O
  // Isotopologue Ratio: 1 1 1 1 1 1 3
  //
  // Some tags relate to empirical continuum correction terms.
  // This would be the place to add additional continuum tags, should
  // this be necessary. Continuum tags must come after the other tags
  // and have NAN for all data entries, like in my example below. Not
  // even the isotopic ratio is used for continuum tags.
  //
  // The isotopologue ratio of NAN is used to identify continuum tags in
  // the absorption routines!
  //
  // You also have to change the entry in the file
  // partition_function_data.cc consistently!
  //
  // 2001-05-30 Stefan Buehler:
  // - Added isotopologue 172 (HITRAN 2000 tag 16)
  // - Gave HITRAN tag 15 to isotopologue 182, which was already in ARTS,
  //   but previously not in HITRAN. Changed isotopologue ratio of this
  //   one to the HITRAN value, which is 6.23e-7, instead of the
  //   previous value of 6.11e-7.
  species_data.push_back(SpeciesRecord(
      NAME("H2O"),
      DEGFR(3),
      ISOTOPOLOGUES(  //   Name,             Isotop. Ratio,  Mass,          MY-tag, HI-tag, JPL-tag
          //                     |               |              |       |       |
          REC("161", .997317E+00, 18.010565, 11, 11, TAGS(18003, 18005)),
          REC("181", 1.99983E-03, 20.014811, 12, 12, TAGS(20003)),
          REC("171", 3.71884E-04, 19.014780, 13, 13, TAGS(19003)),
          REC("162", 3.10693E-04, 19.016740, 14, 14, TAGS(19002)),
          REC("182", 6.23003E-07, 21.020985, -1, 15, TAGS(21001)),
          REC("172", 1.15853E-07, 20.020956, -1, 16, TAGS()),
          REC("262", 2.41970E-08, 20.000000, -1, 17, TAGS(20001)),
          REC("SelfContStandardType", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContStandardType", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContMaTippingType", NAN, NAN, -1, -1, TAGS()),
          REC("ContMPM93", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContCKDMT100", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContCKDMT100", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContCKDMT252", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContCKDMT252", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContCKDMT320", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContCKDMT320", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContCKD222", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContCKD222", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContCKD242", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContCKD242", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContCKD24", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContCKD24", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContATM01", NAN, NAN, -1, -1, TAGS()),
          REC("CP98", NAN, NAN, -1, -1, TAGS()),
          REC("MPM87", NAN, NAN, -1, -1, TAGS()),
          REC("MPM89", NAN, NAN, -1, -1, TAGS()),
          REC("MPM93", NAN, NAN, -1, -1, TAGS()),
          REC("PWR98", NAN, NAN, -1, -1, TAGS()))));

  // CO2
  // (missing mainly in JPL, latest version (7/00) includes some isotopologues)
  // Degrees of freedom from Schanda:`Physical Fundamentals of Remote Sensing'
  // Isotopologue Ratios: 1 1 1 1 1 1 1 1 1 1 1
  // Note (JM): CO2-727 is new in Hitran 2012, and it is not the at least
  // abundant isotopologue. So hitran changed the isotopologue assignments of
  // the less abundant isotopologues: 838, which was 29 before now has hitran
  // tag 20!
  // Again (as for OCS from 2000 edition), this messes up the whole concept of
  // reading a catalogue, because the species depends now on the edition of the
  // catalogue.
  // ATTENTION: This version of species_data.cc works with HITRAN 2008 (and
  // earlier) editions. For use with 2012 (and later), outcomment the 2 lines
  // below "version for 2008 and earlier" and uncomment the corresponding 2
  // lines below "version for 2012 and later".
  // Also, isotopologue 837 is already in molparams.txt, too, as isotopologue
  // #11. It's unclear what will be done when this will get some transitions
  // assigned.
  species_data.push_back(SpeciesRecord(
      NAME("CO2"),
      DEGFR(2),
      ISOTOPOLOGUES(  //   Name,     Isotop. Ratio,  Mass,         MY-tag, HI-tag, JPL-tag
          //             |               |             |       |       |
          REC("626", .984204E+00, 43.989830, 21, 21, TAGS()),
          REC("636", 1.10574E-02, 44.993185, 22, 22, TAGS()),
          REC("628", 3.94707E-03, 45.994076, 23, 23, TAGS(46013)),
          REC("627", 7.33989E-04, 44.994045, 24, 24, TAGS(45012)),
          REC("638", 4.43446E-05, 46.997431, 25, 25, TAGS()),
          REC("637", 8.24623E-06, 45.997400, 26, 26, TAGS()),
          REC("828", 3.95734E-06, 47.998322, 27, 27, TAGS()),
          REC("728", 1.47180E-06, 46.998291, 28, 28, TAGS()),
#ifdef USE_HITRAN2008
          // version for 2008 and earlier
          REC("727", 1.36847E-07, 45.998262, -1, -1, TAGS()),
          REC("838", 4.44600E-08, 49.001675, -1, 29, TAGS()),
#else
          // version for 2012 and later
          REC("727", 1.36847E-07, 45.998262, -1, 29, TAGS()),
          REC("838", 4.44600E-08, 49.001675, -1, 20, TAGS()),
#endif
          REC("837", 1.65354E-08, 48.001646, -1, -1, TAGS()),
          REC("CKD241", NAN, NAN, -1, -1, TAGS()),
          REC("CKDMT100", NAN, NAN, -1, -1, TAGS()),
          REC("CKDMT252", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContPWR93", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContPWR93", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContHo66", NAN, NAN, -1, -1, TAGS()),
          REC("ForeignContHo66", NAN, NAN, -1, -1, TAGS()))));

  // O3
  // Isotopologue Ratios: 1 1 1 1 1
  species_data.push_back(SpeciesRecord(
      NAME("O3"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,    Mass,         MY-tag, HI-tag, JPL-tag
          //             |                 |             |       |       |
          REC("666",
              .992901E+00,
              47.984745,
              31,
              31,
              TAGS(48004, 48005, 48006, 48007, 48008)),
          REC("668", 3.98194E-03, 49.988991, 32, 32, TAGS(50004, 50006)),
          REC("686", 1.99097E-03, 49.988991, 33, 33, TAGS(50003, 50005)),
          REC("667", 7.40475E-04, 48.988960, 34, 34, TAGS(49002)),
          REC("676", 3.70237E-04, 48.988960, 35, 35, TAGS(49001)))));

  // N2O
  // Isotopologue Ratios: 1 1 1 1 1
  species_data.push_back(SpeciesRecord(
      NAME("N2O"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("446", .990333E+00, 44.001062, 41, 41, TAGS(44004, 44009, 44012)),
          REC("456", 3.64093E-03, 44.998096, 42, 42, TAGS(45007)),
          REC("546", 3.64093E-03, 44.998096, 43, 43, TAGS(45008)),
          REC("448", 1.98582E-03, 46.005308, 44, 44, TAGS(46007)),
          REC("447", 3.69280E-04, 45.005278, -1, 45, TAGS()))));

  // CO
  // Isotopologue Ratios: 1 1 1 1 1 1
  species_data.push_back(SpeciesRecord(
      NAME("CO"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,      Mass,         MY-tag, HI-tag, JPL-tag
          //             |                   |             |       |       |
          REC("26", .986544E+00, 27.994915, 51, 51, TAGS(28001)),
          REC("36", 1.10836E-02, 28.998270, 52, 52, TAGS(29001)),
          REC("28", 1.97822E-03, 29.999161, 53, 53, TAGS(30001)),
          REC("27", 3.67867E-04, 28.999130, -1, 54, TAGS(29006)),
          REC("38", 2.22250E-05, 31.002516, -1, 55, TAGS()),
          REC("37", 4.13292E-06, 30.002485, -1, 56, TAGS()))));

  // CH4
  // Degrees of freedom: jpl catalogue
  // Isotopologue Ratios: 1 1 1 1
  // Note: - jpl isotopologue ratio for tag 17003: 0.00014996848
  //       - CH4 is in official mytran list (6), but does not
  //         seem to be included for calculation, as given
  //         by table tag_table in file glob_def.c
  species_data.push_back(SpeciesRecord(
      NAME("CH4"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,       Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //              |                |             |       |       |
          REC("211", .988274E+00, 16.031300, -1, 61, TAGS()),
          // JM: the following version with MADE-UP JPL-TAG(!) is for SWI intercomparison ONLY!
          //         REC("211"      ,.988274E+00     ,16.031300    ,-1     ,61     ,TAGS(16002)),
          REC("311", 1.11031E-02, 17.034655, -1, 62, TAGS()),
          REC("212", 6.15751E-04, 17.037475, -1, 63, TAGS(17003)),
          REC("312", 6.91785E-06, 18.040830, -1, 64, TAGS()))));

  // O2
  // Isotopologue Ratios: 1 1 1
  species_data.push_back(SpeciesRecord(
      NAME("O2"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("66", .995262E+00, 31.989830, 71, 71, TAGS(32001, 32002)),
          REC("68", 3.99141E-03, 33.994076, 72, 72, TAGS(34001)),
          REC("67", 7.42235E-04, 32.994045, 73, 73, TAGS(33002)),
          REC("CIAfunCKDMT100", NAN, NAN, -1, -1, TAGS()),
          REC("v0v0CKDMT100", NAN, NAN, -1, -1, TAGS()),
          REC("v1v0CKDMT100", NAN, NAN, -1, -1, TAGS()),
          REC("visCKDMT252", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContStandardType", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContMPM93", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContPWR93", NAN, NAN, -1, -1, TAGS()),
          REC("PWR98", NAN, NAN, -1, -1, TAGS()),
          REC("PWR93", NAN, NAN, -1, -1, TAGS()),
          REC("PWR88", NAN, NAN, -1, -1, TAGS()),
          REC("MPM93", NAN, NAN, -1, -1, TAGS()),
          REC("TRE05", NAN, NAN, -1, -1, TAGS()),
          REC("MPM92", NAN, NAN, -1, -1, TAGS()),
          REC("MPM89", NAN, NAN, -1, -1, TAGS()),
          REC("MPM87", NAN, NAN, -1, -1, TAGS()),
          REC("MPM85", NAN, NAN, -1, -1, TAGS()),
          REC("MPM2020", NAN, NAN, -1, -1, TAGS()))));

  // NO
  // Isotopologue Ratios: 1 1 1
  species_data.push_back(SpeciesRecord(
      NAME("NO"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,    Mass,         MY-tag, HI-tag, JPL-tag
          //             |                 |             |       |       |
          REC("46", .993974E+00, 29.997989, 81, 81, TAGS(30008)),
          REC("56", 3.65431E-03, 30.995023, -1, 82, TAGS()),
          REC("48", 1.99312E-03, 32.002234, -1, 83, TAGS()))));

  // SO2
  // Isotopologue Ratios: 1 1 2 2
  species_data.push_back(SpeciesRecord(
      NAME("SO2"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("626", .945678E+00, 63.961901, 91, 91, TAGS(64002, 64005)),
          REC("646", 4.19503E-02, 65.957695, 92, 92, TAGS(66002)),
          REC("636", 0.0074989421, 65.00, 93, -1, TAGS(65001)),
          REC("628", 0.0020417379, 66.00, 94, -1, TAGS(66004)))));

  // NO2
  // Isotopologue Ratios: 1
  species_data.push_back(SpeciesRecord(
      NAME("NO2"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,     Mass,         MY-tag, HI-tag, JPL-tag
          //             |                  |             |       |       |
          REC("646", .991616E+00, 45.992904, 101, 101, TAGS(46006)))));

  // NH3
  // Isotopologue Ratios: 1 1 3
  species_data.push_back(SpeciesRecord(
      NAME("NH3"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("4111", .995872E+00, 17.026549, 111, 111, TAGS(17002, 17004)),
          REC("5111", 3.66129E-03, 18.023583, 112, 112, TAGS(18002)),
          REC("4112", 0.00044792294, 18.00, -1, -1, TAGS(18004)))));

  // HNO3
  // Isotopologue Ratios: 1 1
  species_data.push_back(SpeciesRecord(
      NAME("HNO3"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("146",
              .989110E+00,
              62.995644,
              121,
              121,
              TAGS(63001, 63002, 63003, 63004, 63005, 63006)),
          REC("156", 3.63600E-03, 63.992680, -1, 122, TAGS()))));

  // OH
  // Isotopologue Ratios: 1 1 1
  species_data.push_back(SpeciesRecord(
      NAME("OH"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("61", .997473E+00, 17.002740, 131, 131, TAGS(17001)),
          REC("81", 2.00014E-03, 19.006986, 132, 132, TAGS(19001)),
          REC("62", 1.55371E-04, 18.008915, 133, 133, TAGS(18001)))));

  // HF
  // Isotopologue Ratios: 1 1
  species_data.push_back(SpeciesRecord(
      NAME("HF"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("19", .999844E+00, 20.006229, 141, 141, TAGS(20002)),
          REC("29", 1.55741E-04, 21.012404, -1, 142, TAGS(21002)))));

  // HCl
  // Isotopologue Ratios: 1 1 1 1
  species_data.push_back(SpeciesRecord(
      NAME("HCl"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("15", .757587E+00, 35.976678, 151, 151, TAGS(36001)),
          REC("17", .242257E+00, 37.973729, 152, 152, TAGS(38001)),
          REC("25", 1.18005E-04, 36.982853, -1, 153, TAGS(37001)),
          REC("27", 3.77350E-05, 38.979904, -1, 154, TAGS(39004)))));

  // HBr
  // Isotopologue Ratios: 1 1 1 1
  species_data.push_back(SpeciesRecord(
      NAME("HBr"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,    Mass,         MY-tag, HI-tag, JPL-tag
          //             |                 |             |       |       |
          REC("19", .506781E+00, 79.926160, 161, 161, TAGS(80001)),
          REC("11", .493063E+00, 81.924115, 162, 162, TAGS(82001)),
          REC("29", 7.89384E-05, 80.932336, -1, 163, TAGS()),
          REC("21", 7.68016E-05, 82.930289, -1, 164, TAGS()))));

  // HI
  // Degrees of freedom: guessed, since it seems to be linear
  // Isotopologue Ratios: 1 1
  // Note: HI is in official mytran list (17), but does not
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c

  species_data.push_back(SpeciesRecord(
      NAME("HI"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("17", .999844E+00, 127.912297, -1, 171, TAGS()),
          REC("27", 1.55741E-04, 128.918472, -1, 172, TAGS()))));

  // ClO
  // Isotopologue Ratios: 1 1
  species_data.push_back(SpeciesRecord(
      NAME("ClO"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,      Mass,         MY-tag, HI-tag, JPL-tag
          //             |                   |             |       |       |
          REC("56", .755908E+00, 50.963768, 181, 181, TAGS(51002, 51003)),
          REC("76", .241720E+00, 52.960819, 182, 182, TAGS(53002, 53006)))));

  // OCS
  // Isotopologue Ratios: 1 1 1 1 1
  // Note (AvE): OCS-623 is new in Hitran 2000, and it is not the at least
  // abundant isotopologue. So what actually happend with the following
  // one, is that numbers changed in the hitran edition? I had a
  // look at the two editions, and they actually changed them, stupid
  // idiots. So what used to be 194 (hitran 96) is now 195 (hitran
  // 2000). This messes up the whole concept of reading a catalogue,
  // because the species depends now on the edition of the
  // catalogue.
  species_data.push_back(SpeciesRecord(
      NAME("OCS"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,    Mass,         MY-tag, HI-tag, JPL-tag
          //             |                 |             |       |       |
          REC("622", .937395E+00, 59.966986, 191, 191, TAGS(60001)),
          REC("624", 4.15828E-02, 61.962780, 192, 192, TAGS(62001)),
          REC("632", 1.05315E-02, 60.970341, 193, 193, TAGS(61001)),
          REC("623", 7.39908E-03, 60.966371, 194, 194, TAGS()),
          REC("822", 1.87967E-03, 61.971231, 195, 195, TAGS(62002)))));

  // H2CO
  // Isotopologue Ratios: 1 1 1 3 3
  // Note: the isotopologue names differ from hitran convention, since the jpl catalogue has
  //       isotopologues HHCO, HDCO, DDCO.
  //       hitran convention  --  new convention  -- jpl species
  //              126                 1126             H2CO
  //              136                 1136             H2C-13-O
  //              128                 1128             H2CO-18
  //               -                  1226             HDCO
  //               -                  2226             D2CO
  species_data.push_back(SpeciesRecord(
      NAME("H2CO"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("1126", .986237E+00, 30.010565, 201, 201, TAGS(30004)),
          REC("1136", 1.10802E-02, 31.013920, 202, 202, TAGS(31002)),
          REC("1128", 1.97761E-03, 32.014811, 203, 203, TAGS(32004)),
          REC("1226", 0.00029578940, 31.00, -1, -1, TAGS(31003)),
          REC("2226", 2.2181076E-08, 32.00, -1, -1, TAGS(32006)))));

  // HOCl
  // Isotopologue Ratios: 1 1
  species_data.push_back(SpeciesRecord(
      NAME("HOCl"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,      Mass,         MY-tag, HI-tag, JPL-tag
          //             |                   |             |       |       |
          REC("165", .755790E+00, 51.971593, 211, 211, TAGS(52006)),
          REC("167", .241683E+00, 53.968644, 212, 212, TAGS(54005)))));

  // N2
  // Degrees of freedom: guessed, since it seems to be linear
  // Isotopologue Ratios: 1 1
  // Note: N2 is in official mytran list (22), but does not
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back(SpeciesRecord(
      NAME("N2"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,   Mass,         MY-tag, HI-tag, JPL-tag
          //             |                |             |       |       |
          REC("44", .992687E+00, 28.006148, -1, 221, TAGS()),
          REC("45", 7.47809E-03, 29.003182, -1, 222, TAGS()),
          REC("SelfContMPM93", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContPWR93", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContStandardType", NAN, NAN, -1, -1, TAGS()),
          REC("SelfContBorysow", NAN, NAN, -1, -1, TAGS()),
          REC("CIArotCKDMT100", NAN, NAN, -1, -1, TAGS()),
          REC("CIAfunCKDMT100", NAN, NAN, -1, -1, TAGS()),
          REC("CIArotCKDMT252", NAN, NAN, -1, -1, TAGS()),
          REC("CIAfunCKDMT252", NAN, NAN, -1, -1, TAGS()),
          REC("DryContATM01", NAN, NAN, -1, -1, TAGS()))));

  // HCN
  // Isotopologue Ratios: 1 1 1 3
  species_data.push_back(SpeciesRecord(
      NAME("HCN"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,     Mass,         MY-tag, HI-tag, JPL-tag
          //             |                  |             |       |       |
          REC("124", .985114E+00, 27.010899, 231, 231, TAGS(27001, 27003)),
          REC("134", 1.10676E-02, 28.014254, 232, 232, TAGS(28002)),
          REC("125", 3.62174E-03, 28.007933, 233, 233, TAGS(28003)),
          REC("224", 0.00014773545, 28.00, -1, -1, TAGS(28004)))));

  // CH3Cl
  // Isotopologue Ratios: 1 1
  species_data.push_back(SpeciesRecord(
      NAME("CH3Cl"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,     Mass,         MY-tag, HI-tag, JPL-tag
          //             |                  |             |       |       |
          REC("215", .748937E+00, 49.992328, 241, 241, TAGS(50007)),
          REC("217", .239491E+00, 51.989379, 242, 242, TAGS(52009)))));

  // H2O2
  // Isotopologue Ratios: 1
  species_data.push_back(SpeciesRecord(
      NAME("H2O2"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,      Mass,         MY-tag, HI-tag, JPL-tag
          //             |                   |             |       |       |
          REC("1661", .994952E+00, 34.005480, 251, 251, TAGS(34004)))));

  // C2H2
  // Degrees of freedom: guessed, since it seems to be non linear
  // Isotopologue Ratios: 1 1 1
  // Note: C2H2 is in official mytran list (26), but does not
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back(SpeciesRecord(
      NAME("C2H2"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,      Mass,         MY-tag, HI-tag, JPL-tag
          //             |                   |             |       |       |
          REC("1221", .977599E+00, 26.015650, -1, 261, TAGS()),
          REC("1231", 2.19663E-02, 27.019005, -1, 262, TAGS()),
          REC("1222", 3.04550E-04, 27.021825, -1, 263, TAGS()))));

  // C2H6
  // Degrees of freedom: guessed, since it seems to be non linear
  // Isotopologue Ratios: 1 1
  // Note: C2H6 is in official mytran list (27), but does not
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back(SpeciesRecord(
      NAME("C2H6"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("1221", .976990E+00, 30.046950, -1, 271, TAGS()),
          REC("1231", 2.19526E-02, 31.050305, -1, 272, TAGS()))));

  // PH3
  // Isotopologue Ratios: 1
  species_data.push_back(SpeciesRecord(
      NAME("PH3"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("1111", .999533E+00, 33.997238, 281, 281, TAGS(34003)))));

  // COF2
  // Isotopologue Ratios: 1 1
  species_data.push_back(SpeciesRecord(
      NAME("COF2"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("269", .986544E+00, 65.991722, 291, 291, TAGS(66001)),
          REC("369", 1.10834E-02, 66.995083, -1, 292, TAGS()))));

  // SF6
  // Degrees of freedom: guessed, since it seems to be non linear
  // Isotopologue Ratios: 1
  // Note: SF6 is in official mytran list (30), but does not
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back(SpeciesRecord(
      NAME("SF6"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("29", .950180E+00, 145.962492, -1, 301, TAGS()))));

  // H2S
  // Isotopologue Ratios: 1 1 1 2
  species_data.push_back(SpeciesRecord(
      NAME("H2S"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("121", .949884E+00, 33.987721, 311, 311, TAGS(34002)),
          // A. Perrin found, HITRAN isotopologue assignments are wrong (inconsistent
          // between molparams.txt and the line parameter files *.par. This dates back to
          // at least HITRAN1996, so it seems safe to just swap around the HITRAN tags
          // here.
          //         REC("141"      ,4.21369E-02        ,35.983515    ,-1     ,312    ,TAGS( )),
          //         REC("131"      ,7.49766E-03        ,34.987105    ,-1     ,313    ,TAGS( )),
          REC("141", 4.21369E-02, 35.983515, -1, 313, TAGS()),
          REC("131", 7.49766E-03, 34.987105, -1, 312, TAGS()),
          REC("122", 0.00029991625, 35.00, -1, -1, TAGS(35001)))));

  // HCOOH
  // Isotopologue Ratios: 1 3 3 3
  // Note: the isotopologue names differ from hitran convention, since the jpl catalogue has
  //       isotopologues HCOOH, HC-13-OOH, DCOOH, HCOOD
  //       hitran convention  --  new convention  -- jpl species
  //              126                 1261             HCOOH
  //               -                  1361             HC-13-OOH
  //               -                  2261             DCOOH
  //               -                  1262             HCOOD
  species_data.push_back(SpeciesRecord(
      NAME("HCOOH"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("1261", .983898E+00, 46.005480, 321, 321, TAGS(46005)),
          REC("1361", 0.010913149, 47.00, -1, -1, TAGS(47002)),
          REC("2261", 0.00014755369, 47.00, -1, -1, TAGS(47003)),
          REC("1262", 0.00014755369, 47.00, -1, -1, TAGS(47004)))));

  // HO2
  // Isotopologue Ratios: 1
  species_data.push_back(SpeciesRecord(
      NAME("HO2"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("166", .995107E+00, 32.997655, 331, 331, TAGS(33001)))));

  // O
  // Isotopologue Ratios: 1
  species_data.push_back(SpeciesRecord(
      NAME("O"),
      DEGFR(0),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("6", .997628E+00, 15.994915, 341, 341, TAGS(16001)))));

  // ClONO2
  // Isotopologue Ratios: 1 1
  // Note: ClONO2 in hitran is identical to ClNO3 in jpl (according to Johannes Orphal)
  species_data.push_back(SpeciesRecord(
      NAME("ClONO2"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("5646", .749570E+00, 96.956672, 351, 351, TAGS(97002)),
          REC("7646", .239694E+00, 98.953723, 352, 352, TAGS(99001)))));

  // NO+
  // Degrees of freedom: guessed, since it seems to be linear
  // Isotopologue Ratios: 1
  // Note: NO+ is in official mytran list (36), but does not
  //       seem to be included for calculation, as given
  //       by table tag_table in file glob_def.c
  species_data.push_back(SpeciesRecord(
      NAME("NO+"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("46", .993974E+00, 29.997989, -1, 361, TAGS(30011)))));

  // OClO
  // Isotopologue Ratios: 2 2
  species_data.push_back(SpeciesRecord(
      NAME("OClO"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,     MY-tag, HI-tag, JPL-tag
          //             |                    |         |       |       |
          REC("656", 0.75509223, 67.00, 431, -1, TAGS(67001)),
          REC("676", 0.24490632, 69.00, 432, -1, TAGS(69001)))));

  // BrO
  // Isotopologue Ratios: 2 2
  species_data.push_back(SpeciesRecord(
      NAME("BrO"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,     MY-tag, HI-tag, JPL-tag
          //             |                    |         |       |       |
          REC("96", 0.50582466, 95.00, 401, -1, TAGS(95001)),
          REC("16", 0.49431069, 97.00, 402, -1, TAGS(97001)))));

  // H2SO4
  // Isotopologue Ratios: 2
  species_data.push_back(SpeciesRecord(
      NAME("H2SO4"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,     MY-tag, HI-tag, JPL-tag
          //             |                    |         |       |       |
          REC("126", 0.95060479, 98.00, 481, -1, TAGS(98001)))));

  // Cl2O2
  // Isotopologue Ratios: 2 2
  // Note: refered to as Cl2O2 in mytran catalogue, in jpl cat: ClOOCl
  species_data.push_back(SpeciesRecord(
      NAME("Cl2O2"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,     MY-tag, HI-tag, JPL-tag
          //             |                    |         |       |       |
          REC("565", 0.57016427, 102.00, 491, -1, TAGS(102001)),
          REC("765", 0.36982818, 104.00, 492, -1, TAGS(104001)))));

  // HOBr
  // Isotopologue Ratios: 1 1
  // Note: latest addition to Hitran 2000, DEGFR guessed
  species_data.push_back(SpeciesRecord(
      NAME("HOBr"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("169", .505579E+00, 95.921076, 371, 371, TAGS(96001)),
          REC("161", .491894E+00, 97.919027, 372, 372, TAGS(98002)))));

  // C2H4
  // Isotopologue Ratios: 1 1
  // Note: latest addition to Hitran 2000, DEGFR guessed
  species_data.push_back(SpeciesRecord(
      NAME("C2H4"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("221", .977294E+00, 28.031300, -1, 381, TAGS()),
          REC("231", .219595E-01, 29.034655, -1, 382, TAGS()))));

  // CH3OH
  // Isotopologue Ratios: 1
  // Note: added with Hitran 2008, DEGFR from JPL.
  species_data.push_back(SpeciesRecord(
      NAME("CH3OH"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("2161", .985930E+00, 32.026215, -1, 391, TAGS(32003))
          //         REC("2261"     ,1.0?               ,33.00        ,-1     ,-1     ,TAGS(33004)) //in JPL this is seems to be taken as separate species (IsoCorr=0)
          )));

  // CH3Br
  // Isotopologue Ratios: 1 1
  // Note: added with Hitran 2008, DEGFR guessed.
  species_data.push_back(SpeciesRecord(
      NAME("CH3Br"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("219", .500995E+00, 93.941811, -1, 401, TAGS()),
          REC("211", .487433E+00, 95.939764, -1, 402, TAGS()))));

  // CH3CN
  // Isotopologue Ratios: 1 3 3 3 3
  // Note: the isotopologue names have been set before CH3CH appeared in hitran,
  // hence naming is different from there
  //       hitran convention  --  our convention  -- jpl species
  //             2124                 212124             CH3CN
  //             3124                 311124             C-13-H3CN (?)
  //             2134                 211134             CH3C-13-N (?)
  //               -                  211125             CH3CN-15
  //               -                  211224             CH2DCN
  //             3134                 311134                -      (?)
  species_data.push_back(SpeciesRecord(
      NAME("CH3CN"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("211124", .973866E+00, 41.026549, -1, 411, TAGS(41001, 41010)),
          REC("311124", .102683e-01, 42.00, -1, -1, TAGS(42006)),
          REC("211134", .102683e-01, 42.00, -1, -1, TAGS(42007)),
          REC("211125", .347136e-02, 42.00, -1, -1, TAGS(42001)),
          REC("211224", .441185e-03, 42.00, -1, -1, TAGS(42008)))));

  // CF4
  // Isotopologue Ratios: 1
  // Note: added with Hitran 2008, DEGFR guessed.
  species_data.push_back(SpeciesRecord(
      NAME("CF4"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("29", .988890E+00, 87.993616, -1, 421, TAGS()))));

  // HC3N
  // Isotopologue Ratios: 1 3 3 3 3 3
  // Note: added with Hitran 2008, where lines are available, but no entry on
  // species parameters in molparams.txt yet. Name and Mass are temporary
  // guesses. DEGFR from JPL.
  // We took HCCCN from JPL, which has also HCCNC and HNCCC, but we assume those
  // to be different species.
  // (temporary) name convention:
  //       hitran convention  --  our convention  -- jpl species
  //             1224                 12224             HCCCN
  //               -                  12234             HCCC-13-N
  //               -                  12324             HCC-13-CN
  //               -                  13224             HC-13-CCN
  //               -                  12225             HCCCN-15
  //               -                  22224             DCCCN
  species_data.push_back(SpeciesRecord(
      NAME("HC3N"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("12224", .963346E+00, 51.010899, -1, 441, TAGS(51001)),
          REC("12234", .106852e-01, 52.00, -1, -1, TAGS(52001)),
          REC("12324", .106852e-01, 52.00, -1, -1, TAGS(52002)),
          REC("13224", .106852e-01, 52.00, -1, -1, TAGS(52003)),
          REC("12225", .356272e-02, 52.00, -1, -1, TAGS(52004)),
          REC("22224", .144472e-03, 52.00, -1, -1, TAGS(52005)))));

  // CS
  // Isotopologue Ratios: 1 1 1 1
  // Note: added with Hitran 2008, DEGFR from JPL.
  species_data.push_back(SpeciesRecord(
      NAME("CS"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("22", .939624E+00, 43.971036, -1, 461, TAGS(44001)),
          REC("24", .416817E-01, 45.966787, -1, 462, TAGS(46001)),
          REC("32", .105565E-01, 44.974368, -1, 463, TAGS(45001)),
          REC("23", .741668E-02, 44.970399, -1, 464, TAGS()))));

  // HNC
  // Isotopologue Ratios: 3 3 3 3
  species_data.push_back(SpeciesRecord(
      NAME("HNC"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("142", .985280e+00, 27.00, -1, -1, TAGS(27002, 27003)),
          REC("143", .109285e-01, 28.00, -1, -1, TAGS(28005)),
          REC("152", .364384e-02, 28.00, -1, -1, TAGS(28006)),
          REC("242", .147761e-03, 28.00, -1, -1, TAGS(28007)))));

  // SO
  // Isotopologue Ratios: 2 2 2
  // Note: Name and Mass guessed. DEGFR from JPL.
  species_data.push_back(SpeciesRecord(
      NAME("SO"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("26", .950605e+00, 48.00, -1, -1, TAGS(48001, 48002)),
          REC("46", .420727e-01, 50.00, -1, -1, TAGS(50001)),
          REC("28", .194089e-02, 50.00, -1, -1, TAGS(50002)))));

  // C3H8
  // Isotopologue Ratios: 4
  // Note: Name and Mass guessed. DEGFR and atomic isotopic ratios from JPL.
  species_data.push_back(SpeciesRecord(
      NAME("C3H8"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("21", 9.66290e-01, 44.00, -1, -1, TAGS(44013)))));

  // H2
  // Isotopologue Ratios: 1 1 4
  // Note: H2 is spectroscopically inert, but now needed as broadening/continuum
  // species for planetary atmospheres. Hence we need it to be defined as
  // absorption species.
  // Name and Mass guessed. DEGFR and atomic isotopic ratios from JPL.
  species_data.push_back(SpeciesRecord(
      NAME("H2"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("11", .999688E+00, 2.0156500, -1, 451, TAGS()),
          REC("12", 3.11432E-04, 3.0218250, -1, 452, TAGS(3001))
          //         REC("22"       ,.224838e-07        ,4.00         ,-1     ,-1     ,TAGS())
          )));

  // H
  // Isotopologue Ratios: 0
  // Note: H is spectroscopically inert, but now needed as broadening/continuum
  // species for planetary atmospheres. Hence we need it to be defined as
  // absorption species.
  // All parameters guessed.
  species_data.push_back(SpeciesRecord(
      NAME("H"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("1", 1.00, 1.00, -1, -1, TAGS()))));

  // He
  // Isotopologue Ratios: 0
  // Note: He is spectroscopically inert, but now needed as broadening/continuum
  // species for planetary atmospheres. Hence we need it to be defined as
  // absorption species.
  // All parameters guessed.
  species_data.push_back(SpeciesRecord(
      NAME("He"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("4", 1.00, 4.00, -1, -1, TAGS()))));

  // Ar
  // Isotopologue Ratios: 0
  // Note: Ar is spectroscopically inert, but now needed as broadening/continuum
  // species for planetary atmospheres. Hence it needs to be defined as
  // absorption species.
  // All parameters guessed.
  species_data.push_back(SpeciesRecord(
      NAME("Ar"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("8", 1.00, 18.00, -1, -1, TAGS()))));

  // C4H2
  // Isotopologue Ratios: 1
  // Note: added with HITRAN2012, Mass and DEGFR guessed.
  species_data.push_back(SpeciesRecord(
      NAME("C4H2"),
      DEGFR(2),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("2211", .955998E+00, 50.015650, -1, 431, TAGS()))));

  // SO3
  // Isotopologue Ratios: 1
  // Note: added with HITRAN2012, Mass and DEGFR guessed.
  species_data.push_back(SpeciesRecord(
      NAME("SO3"),
      DEGFR(3),
      ISOTOPOLOGUES(  //  Name,      Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("26", .943400E+00, 79.956820, -1, 471, TAGS()))));

  species_data.push_back(SpeciesRecord(
      NAME("liquidcloud"),
      DEGFR(0),
      ISOTOPOLOGUES(  //   Name,     Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("MPM93", NAN, NAN, -1, -1, TAGS()),
          REC("ELL07", NAN, NAN, -1, -1, TAGS()))));

  species_data.push_back(SpeciesRecord(
      NAME("icecloud"),
      DEGFR(0),
      ISOTOPOLOGUES(  //   Name,     Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("MPM93", NAN, NAN, -1, -1, TAGS()))));

  species_data.push_back(SpeciesRecord(
      NAME("rain"),
      DEGFR(0),
      ISOTOPOLOGUES(  //   Name,     Isotop. Ratio,       Mass,         MY-tag, HI-tag, JPL-tag
          //             |                    |             |       |       |
          REC("MPM93", NAN, NAN, -1, -1, TAGS()))));

  species_data.push_back(
      SpeciesRecord(NAME("free_electrons"), DEGFR(0), ISOTOPOLOGUES()));

  species_data.push_back(
      SpeciesRecord(NAME("particles"), DEGFR(0), ISOTOPOLOGUES()));

  // hitran cross section species
  for (const auto& s : {
           "C2F6",     "C3F8",      "C4F10",     "C5F12",     "C6F14",
           "C8F18",    "cC4F8",     "CCl4",      "CFC11",     "CFC113",
           "CFC114",   "CFC115",    "CFC12",     "CH2Cl2",    "CH3CCl3",
           "CHCl3",    "Halon1211", "Halon1301", "Halon2402", "HCFC141b",
           "HCFC142b", "HCFC22",    "HFC125",    "HFC134a",   "HFC143a",
           "HFC152a",  "HFC227ea",  "HFC23",     "HFC245fa",  "HFC32",
           "NF3",      "SO2F2",     "HFC4310mee"
       }) {
    species_data.push_back(SpeciesRecord(NAME(s), DEGFR(0), ISOTOPOLOGUES()));
  }

  // You also have to change the entry in the file
  // partition_function_data.cc consistently!
}

/*! Define the species data map.

 \author Stefan Buehler */
void define_species_map() {
  using global_data::species_data;
  using global_data::SpeciesMap;

  for (Index i = 0; i < species_data.nelem(); ++i) {
    SpeciesMap[species_data[i].Name()] = i;
  }
}
