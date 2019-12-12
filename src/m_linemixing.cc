/* Copyright 2013, The ARTS Developers.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <Eigen/Eigenvalues>
#include <Faddeeva/Faddeeva.hh>
#include "arts.h"
#include "auto_md.h"
#include "file.h"
#include "global_data.h"
#include "lin_alg.h"
#include "linefunctions.h"
#include "linemixing.h"
#include "linescaling.h"
#include "physics_funcs.h"
#include "wigner_functions.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void SetBandIdentifiersAuto(ArrayOfQuantumIdentifier& band_identifiers,
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const Verbosity& verbosity) {
  CREATE_OUT3;

  out3
      << "Will set band_identifiers to as many bands as is included by default in this setup.\n"
         "To extend these bands, please edit this function. Currently, supported species and \n"
         "bands are only for CO2 and O2-66, though not all bands are guaranteed to be available.\n";

  band_identifiers.resize(0);

  static const SpeciesTag CO2_626 = SpeciesTag("CO2-626");
  static const SpeciesTag CO2_627 = SpeciesTag("CO2-627");
  static const SpeciesTag CO2_628 = SpeciesTag("CO2-628");
  static const SpeciesTag CO2_636 = SpeciesTag("CO2-636");
  static const SpeciesTag CO2_638 = SpeciesTag("CO2-638");
  static const SpeciesTag CO2_828 = SpeciesTag("CO2-828");
  static const SpeciesTag O2_66 = SpeciesTag("O2-66");
  bool o2_66_done, co2_626_done, co2_638_done, co2_636_done, co2_627_done,
      co2_628_done, co2_828_done;
  o2_66_done = co2_626_done = co2_638_done = co2_636_done = co2_627_done =
      co2_628_done = co2_828_done = false;

  //for( const ArrayOfSpeciesTag& spec : abs_species)
  for (Index ispec = 0; ispec < abs_species.nelem(); ispec++) {
    const ArrayOfSpeciesTag& spec = abs_species[ispec];
    //for( const SpeciesTag& tag : spec )
    for (Index iabsorber = 0; iabsorber < spec.nelem(); iabsorber++) {
      const SpeciesTag& tag = spec[iabsorber];
      if (CO2_626.Species() == tag.Species()) {
        if (CO2_626.Isotopologue() == tag.Isotopologue() and
            not co2_626_done) {
          co2_626_done = true;
          QuantumIdentifier qi;
          Index nbands = 0;

          //Frequency is around 471 cm-1
          qi.SetFromStringForCO2Band("20003", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 479 cm-1
          qi.SetFromStringForCO2Band("13302", "12201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 544 cm-1
          qi.SetFromStringForCO2Band("11102", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 568 cm-1
          qi.SetFromStringForCO2Band("13302", "04401", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 594 cm-1
          qi.SetFromStringForCO2Band("20002", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 597 cm-1
          qi.SetFromStringForCO2Band("11102", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 608 cm-1
          qi.SetFromStringForCO2Band("10012", "01111", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 615 cm-1
          qi.SetFromStringForCO2Band("20003", "11102", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 618 cm-1
          qi.SetFromStringForCO2Band("10002", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 647 cm-1
          qi.SetFromStringForCO2Band("11102", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 654 cm-1
          qi.SetFromStringForCO2Band("01111", "00011", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 655 cm-1
          qi.SetFromStringForCO2Band("02211", "01111", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 667 cm-1
          qi.SetFromStringForCO2Band("01101", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 667 cm-1
          qi.SetFromStringForCO2Band("02201", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 668 cm-1
          qi.SetFromStringForCO2Band("03301", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 668 cm-1
          qi.SetFromStringForCO2Band("04401", "03301", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 681 cm-1
          qi.SetFromStringForCO2Band("13301", "12201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 683 cm-1
          qi.SetFromStringForCO2Band("12201", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 688 cm-1
          qi.SetFromStringForCO2Band("11101", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 703 cm-1
          qi.SetFromStringForCO2Band("21101", "20001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 710 cm-1
          qi.SetFromStringForCO2Band("10011", "01111", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 720 cm-1
          qi.SetFromStringForCO2Band("20001", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 720 cm-1
          qi.SetFromStringForCO2Band("10001", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 738 cm-1
          qi.SetFromStringForCO2Band("20002", "11102", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 739 cm-1
          qi.SetFromStringForCO2Band("21101", "12201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 741 cm-1
          qi.SetFromStringForCO2Band("11101", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 757 cm-1
          qi.SetFromStringForCO2Band("12201", "03301", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 770 cm-1
          qi.SetFromStringForCO2Band("13301", "04401", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 791 cm-1
          qi.SetFromStringForCO2Band("11101", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 828 cm-1
          qi.SetFromStringForCO2Band("12201", "11102", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 829 cm-1
          qi.SetFromStringForCO2Band("21101", "20002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 864 cm-1
          qi.SetFromStringForCO2Band("20001", "11102", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 898 cm-1
          qi.SetFromStringForCO2Band("02211", "12201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 917 cm-1
          qi.SetFromStringForCO2Band("10011", "20001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 927 cm-1
          qi.SetFromStringForCO2Band("01111", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 941 cm-1
          qi.SetFromStringForCO2Band("10012", "20002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 952 cm-1
          qi.SetFromStringForCO2Band("21101", "20003", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 960 cm-1
          qi.SetFromStringForCO2Band("00011", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1043 cm-1
          qi.SetFromStringForCO2Band("10011", "20002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1063 cm-1
          qi.SetFromStringForCO2Band("00011", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1064 cm-1
          qi.SetFromStringForCO2Band("10012", "20003", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1071 cm-1
          qi.SetFromStringForCO2Band("01111", "11102", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1880 cm-1
          qi.SetFromStringForCO2Band("20003", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1905 cm-1
          qi.SetFromStringForCO2Band("13302", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1932 cm-1
          qi.SetFromStringForCO2Band("11102", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2003 cm-1
          qi.SetFromStringForCO2Band("20002", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2076 cm-1
          qi.SetFromStringForCO2Band("11101", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2093 cm-1
          qi.SetFromStringForCO2Band("12201", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2107 cm-1
          qi.SetFromStringForCO2Band("13301", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2112 cm-1
          qi.SetFromStringForCO2Band("21101", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2129 cm-1
          qi.SetFromStringForCO2Band("20001", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2165 cm-1
          qi.SetFromStringForCO2Band("21101", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2170 cm-1
          qi.SetFromStringForCO2Band("11112", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2180 cm-1
          qi.SetFromStringForCO2Band("20012", "20001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2182 cm-1
          qi.SetFromStringForCO2Band("20013", "20002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2224 cm-1
          qi.SetFromStringForCO2Band("10012", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2288 cm-1
          qi.SetFromStringForCO2Band("13311", "13301", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2290 cm-1
          qi.SetFromStringForCO2Band("13312", "13302", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2299 cm-1
          qi.SetFromStringForCO2Band("04411", "04401", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2299 cm-1
          qi.SetFromStringForCO2Band("02221", "02211", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2301 cm-1
          qi.SetFromStringForCO2Band("12211", "12201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2301 cm-1
          qi.SetFromStringForCO2Band("10021", "10011", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2302 cm-1
          qi.SetFromStringForCO2Band("10022", "10012", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2302 cm-1
          qi.SetFromStringForCO2Band("20011", "20001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2305 cm-1
          qi.SetFromStringForCO2Band("20013", "20003", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2306 cm-1
          qi.SetFromStringForCO2Band("20012", "20002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2311 cm-1
          qi.SetFromStringForCO2Band("03311", "03301", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2311 cm-1
          qi.SetFromStringForCO2Band("01121", "01111", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2313 cm-1
          qi.SetFromStringForCO2Band("11111", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2315 cm-1
          qi.SetFromStringForCO2Band("11112", "11102", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2324 cm-1
          qi.SetFromStringForCO2Band("02211", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2324 cm-1
          qi.SetFromStringForCO2Band("00021", "00011", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2326 cm-1
          qi.SetFromStringForCO2Band("10011", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2327 cm-1
          qi.SetFromStringForCO2Band("10012", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2336 cm-1
          qi.SetFromStringForCO2Band("01111", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2349 cm-1
          qi.SetFromStringForCO2Band("00011", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2428 cm-1
          qi.SetFromStringForCO2Band("20011", "20002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2429 cm-1
          qi.SetFromStringForCO2Band("10011", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2429 cm-1
          qi.SetFromStringForCO2Band("20012", "20003", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2458 cm-1
          qi.SetFromStringForCO2Band("11111", "11102", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3465 cm-1
          qi.SetFromStringForCO2Band("20013", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3500 cm-1
          qi.SetFromStringForCO2Band("21101", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3527 cm-1
          qi.SetFromStringForCO2Band("30014", "20003", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3527 cm-1
          qi.SetFromStringForCO2Band("13312", "03301", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3550 cm-1
          qi.SetFromStringForCO2Band("30012", "20001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3556 cm-1
          qi.SetFromStringForCO2Band("30013", "20002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3566 cm-1
          qi.SetFromStringForCO2Band("10022", "00011", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3568 cm-1
          qi.SetFromStringForCO2Band("20013", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3580 cm-1
          qi.SetFromStringForCO2Band("11112", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3589 cm-1
          qi.SetFromStringForCO2Band("20012", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3612 cm-1
          qi.SetFromStringForCO2Band("10012", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3667 cm-1
          qi.SetFromStringForCO2Band("10021", "00011", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3676 cm-1
          qi.SetFromStringForCO2Band("30012", "20002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3679 cm-1
          qi.SetFromStringForCO2Band("30013", "20003", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3692 cm-1
          qi.SetFromStringForCO2Band("20012", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3705 cm-1
          qi.SetFromStringForCO2Band("30011", "20001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3711 cm-1
          qi.SetFromStringForCO2Band("20011", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3714 cm-1
          qi.SetFromStringForCO2Band("10011", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3723 cm-1
          qi.SetFromStringForCO2Band("11111", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3726 cm-1
          qi.SetFromStringForCO2Band("12211", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3727 cm-1
          qi.SetFromStringForCO2Band("13311", "03301", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3799 cm-1
          qi.SetFromStringForCO2Band("30012", "20003", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3814 cm-1
          qi.SetFromStringForCO2Band("20011", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4005 cm-1
          qi.SetFromStringForCO2Band("00021", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4687 cm-1
          qi.SetFromStringForCO2Band("30014", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4786 cm-1
          qi.SetFromStringForCO2Band("31113", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4790 cm-1
          qi.SetFromStringForCO2Band("30014", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4839 cm-1
          qi.SetFromStringForCO2Band("30013", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4853 cm-1
          qi.SetFromStringForCO2Band("20013", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4931 cm-1
          qi.SetFromStringForCO2Band("31113", "11102", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4942 cm-1
          qi.SetFromStringForCO2Band("30013", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4946 cm-1
          qi.SetFromStringForCO2Band("31112", "11101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4959 cm-1
          qi.SetFromStringForCO2Band("30012", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4977 cm-1
          qi.SetFromStringForCO2Band("20012", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5062 cm-1
          qi.SetFromStringForCO2Band("30012", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5099 cm-1
          qi.SetFromStringForCO2Band("20011", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5114 cm-1
          qi.SetFromStringForCO2Band("30011", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5217 cm-1
          qi.SetFromStringForCO2Band("30011", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5247 cm-1
          qi.SetFromStringForCO2Band("10022", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5291 cm-1
          qi.SetFromStringForCO2Band("02221", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5315 cm-1
          qi.SetFromStringForCO2Band("01121", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5349 cm-1
          qi.SetFromStringForCO2Band("10012", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5584 cm-1
          qi.SetFromStringForCO2Band("00031", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5687 cm-1
          qi.SetFromStringForCO2Band("00031", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6075 cm-1
          qi.SetFromStringForCO2Band("30014", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6196 cm-1
          qi.SetFromStringForCO2Band("31113", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6227 cm-1
          qi.SetFromStringForCO2Band("30013", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6347 cm-1
          qi.SetFromStringForCO2Band("30012", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6356 cm-1
          qi.SetFromStringForCO2Band("31112", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6503 cm-1
          qi.SetFromStringForCO2Band("30011", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6860 cm-1
          qi.SetFromStringForCO2Band("03331", "03301", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6897 cm-1
          qi.SetFromStringForCO2Band("02231", "02201", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6905 cm-1
          qi.SetFromStringForCO2Band("10031", "10001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6907 cm-1
          qi.SetFromStringForCO2Band("10032", "10002", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6935 cm-1
          qi.SetFromStringForCO2Band("01131", "01101", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6972 cm-1
          qi.SetFromStringForCO2Band("00031", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 8192 cm-1
          qi.SetFromStringForCO2Band("10032", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 8293 cm-1
          qi.SetFromStringForCO2Band("10031", "00001", "626");
          band_identifiers.push_back(qi);
          nbands++;

          out3 << "Set " << nbands << " for CO2-626";
        } else if (CO2_636.Isotopologue() == tag.Isotopologue() and
                    not co2_636_done) {
          co2_636_done = true;
          QuantumIdentifier qi;
          Index nbands = 0;

          //Frequency is around 617 cm-1
          qi.SetFromStringForCO2Band("10002", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 636 cm-1
          qi.SetFromStringForCO2Band("01111", "00011", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 648 cm-1
          qi.SetFromStringForCO2Band("01101", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 648 cm-1
          qi.SetFromStringForCO2Band("02201", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 649 cm-1
          qi.SetFromStringForCO2Band("03301", "02201", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 667 cm-1
          qi.SetFromStringForCO2Band("11101", "10001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 721 cm-1
          qi.SetFromStringForCO2Band("10001", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 739 cm-1
          qi.SetFromStringForCO2Band("11101", "02201", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 771 cm-1
          qi.SetFromStringForCO2Band("11101", "10002", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 883 cm-1
          qi.SetFromStringForCO2Band("01111", "11101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 913 cm-1
          qi.SetFromStringForCO2Band("00011", "10001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1017 cm-1
          qi.SetFromStringForCO2Band("00011", "10002", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2037 cm-1
          qi.SetFromStringForCO2Band("11101", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2157 cm-1
          qi.SetFromStringForCO2Band("10012", "10001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2248 cm-1
          qi.SetFromStringForCO2Band("01121", "01111", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2250 cm-1
          qi.SetFromStringForCO2Band("11111", "11101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2260 cm-1
          qi.SetFromStringForCO2Band("02211", "02201", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2260 cm-1
          qi.SetFromStringForCO2Band("00021", "00011", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2261 cm-1
          qi.SetFromStringForCO2Band("10012", "10002", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2262 cm-1
          qi.SetFromStringForCO2Band("10011", "10001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2271 cm-1
          qi.SetFromStringForCO2Band("01111", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2283 cm-1
          qi.SetFromStringForCO2Band("00011", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2367 cm-1
          qi.SetFromStringForCO2Band("10011", "10002", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2498 cm-1
          qi.SetFromStringForCO2Band("11112", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3450 cm-1
          qi.SetFromStringForCO2Band("13312", "03301", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3460 cm-1
          qi.SetFromStringForCO2Band("21113", "11102", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3473 cm-1
          qi.SetFromStringForCO2Band("12212", "02201", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3482 cm-1
          qi.SetFromStringForCO2Band("20013", "10002", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3482 cm-1
          qi.SetFromStringForCO2Band("21112", "11101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3497 cm-1
          qi.SetFromStringForCO2Band("30001", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3517 cm-1
          qi.SetFromStringForCO2Band("20012", "10001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3527 cm-1
          qi.SetFromStringForCO2Band("10012", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3621 cm-1
          qi.SetFromStringForCO2Band("20011", "10001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3621 cm-1
          qi.SetFromStringForCO2Band("20012", "10002", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3623 cm-1
          qi.SetFromStringForCO2Band("21112", "11102", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3625 cm-1
          qi.SetFromStringForCO2Band("21111", "11101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3632 cm-1
          qi.SetFromStringForCO2Band("10011", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3639 cm-1
          qi.SetFromStringForCO2Band("11111", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3641 cm-1
          qi.SetFromStringForCO2Band("12211", "02201", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3641 cm-1
          qi.SetFromStringForCO2Band("13311", "03301", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3725 cm-1
          qi.SetFromStringForCO2Band("20011", "10002", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4748 cm-1
          qi.SetFromStringForCO2Band("20013", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4887 cm-1
          qi.SetFromStringForCO2Band("20012", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 4991 cm-1
          qi.SetFromStringForCO2Band("20011", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5013 cm-1
          qi.SetFromStringForCO2Band("21111", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 5168 cm-1
          qi.SetFromStringForCO2Band("01121", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6745 cm-1
          qi.SetFromStringForCO2Band("01131", "01101", "636");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 6780 cm-1
          qi.SetFromStringForCO2Band("00031", "00001", "636");
          band_identifiers.push_back(qi);
          nbands++;

          out3 << "Set " << nbands << " for CO2-636";
        } else if (CO2_628.Isotopologue() == tag.Isotopologue() and
                    not co2_628_done) {
          co2_628_done = true;
          QuantumIdentifier qi;
          Index nbands = 0;

          //Frequency is around 597 cm-1
          qi.SetFromStringForCO2Band("10002", "01101", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 662 cm-1
          qi.SetFromStringForCO2Band("01101", "00001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 662 cm-1
          qi.SetFromStringForCO2Band("02201", "01101", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 703 cm-1
          qi.SetFromStringForCO2Band("10001", "01101", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 966 cm-1
          qi.SetFromStringForCO2Band("00011", "10001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1072 cm-1
          qi.SetFromStringForCO2Band("00011", "10002", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1259 cm-1
          qi.SetFromStringForCO2Band("10002", "00001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1365 cm-1
          qi.SetFromStringForCO2Band("10001", "00001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2205 cm-1
          qi.SetFromStringForCO2Band("10012", "10001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2307 cm-1
          qi.SetFromStringForCO2Band("02211", "02201", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2309 cm-1
          qi.SetFromStringForCO2Band("10011", "10001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2311 cm-1
          qi.SetFromStringForCO2Band("10012", "10002", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2319 cm-1
          qi.SetFromStringForCO2Band("01111", "01101", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2332 cm-1
          qi.SetFromStringForCO2Band("00011", "00001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2415 cm-1
          qi.SetFromStringForCO2Band("10011", "10002", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2614 cm-1
          qi.SetFromStringForCO2Band("20002", "00001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2618 cm-1
          qi.SetFromStringForCO2Band("21102", "01101", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2757 cm-1
          qi.SetFromStringForCO2Band("20001", "00001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3571 cm-1
          qi.SetFromStringForCO2Band("10012", "00001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3675 cm-1
          qi.SetFromStringForCO2Band("10011", "00001", "628");
          band_identifiers.push_back(qi);
          nbands++;

          out3 << "Set " << nbands << " for CO2-628";
        } else if (CO2_828.Isotopologue() == tag.Isotopologue() and
                    not co2_828_done) {
          co2_828_done = true;
          QuantumIdentifier qi;
          Index nbands = 0;

          //Frequency is around 2301 cm-1
          qi.SetFromStringForCO2Band("01111", "01101", "828");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2314 cm-1
          qi.SetFromStringForCO2Band("00011", "00001", "828");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3525 cm-1
          qi.SetFromStringForCO2Band("10012", "00001", "828");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3638 cm-1
          qi.SetFromStringForCO2Band("10011", "00001", "828");
          band_identifiers.push_back(qi);
          nbands++;

          out3 << "Set " << nbands << " for CO2-828";
        } else if (CO2_627.Isotopologue() == tag.Isotopologue() and
                    not co2_627_done) {
          co2_627_done = true;
          QuantumIdentifier qi;
          Index nbands = 0;

          //Frequency is around 664 cm-1
          qi.SetFromStringForCO2Band("01101", "00001", "627");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 711 cm-1
          qi.SetFromStringForCO2Band("10001", "01101", "627");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 963 cm-1
          qi.SetFromStringForCO2Band("00011", "10001", "627");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 1376 cm-1
          qi.SetFromStringForCO2Band("10001", "00001", "627");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2340 cm-1
          qi.SetFromStringForCO2Band("00011", "00001", "627");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 2641 cm-1
          qi.SetFromStringForCO2Band("20002", "00001", "627");
          band_identifiers.push_back(qi);
          nbands++;

          out3 << "Set " << nbands << " for CO2-627";
        } else if (CO2_638.Isotopologue() == tag.Isotopologue() and
                    not co2_638_done) {
          co2_638_done = true;
          QuantumIdentifier qi;
          Index nbands = 0;

          //Frequency is around 2265 cm-1
          qi.SetFromStringForCO2Band("00011", "00001", "638");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3490 cm-1
          qi.SetFromStringForCO2Band("10012", "00001", "638");
          band_identifiers.push_back(qi);
          nbands++;

          //Frequency is around 3587 cm-1
          qi.SetFromStringForCO2Band("10011", "00001", "638");
          band_identifiers.push_back(qi);
          nbands++;

          out3 << "Set " << nbands << " for CO2-638";
        }
      } else if (O2_66.Species() == tag.Species() and
                  O2_66.Isotopologue() == tag.Isotopologue() and
                  not o2_66_done) {
        o2_66_done = true;
        QuantumIdentifier qi;
        // The main band in the 60 GHz to 1.5 THz range.  LM mostly at 60 GHz, though the rest falls into this band.
        qi.SetFromString("O2-66 TR UP v1 0/1 LO v1 0/1");
        band_identifiers.push_back(qi);
        // The secondary band in the 60 GHz to 1.5 THz range.  LM mostly at 60 GHz, though the rest falls into this band.
        qi.SetFromString("O2-66 TR UP v1 1/1 LO v1 1/1");
        band_identifiers.push_back(qi);
      } else {
        throw std::runtime_error(
            "Unsupported species (edit m_linemixing.cc to fix or manually enter band_identifiers).\n");
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SetBandIdentifiersFromLines(ArrayOfQuantumIdentifier& band_identifiers,
                                 const ArrayOfAbsorptionLines& abs_lines,
                                 const QuantumIdentifier& band_quantums,
                                 const Index& global,
                                 const Verbosity&) {
  if (not band_quantums.IsEnergyLevelType())
    throw std::runtime_error("band_quantums must be energy type");

  auto& qns = band_quantums.EnergyLevelQuantumNumbers();

  if (global) {
    band_identifiers.resize(abs_lines.nelem());
    for (Index k=0; k<abs_lines.nelem(); k++) {
      band_identifiers[k] = abs_lines[k].QuantumIdentity();
    }
  } else {
    band_identifiers.resize(0);
    for (auto& band: abs_lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        QuantumIdentifier qid = band.QuantumIdentityOfLine(k);  // Copy
        auto& upper = qid.UpperQuantumNumbers();
        auto& lower = qid.LowerQuantumNumbers();

        for (Index i = 0; i < Index(QuantumNumberType::FINAL_ENTRY); i++) {
          if (qns[i].isUndefined()) {
            upper.Set(i, RATIONAL_UNDEFINED);  // Make undefined unwanted numbers
            lower.Set(i, RATIONAL_UNDEFINED);  // Make undefined unwanted numbers
          }
        }

        bool copy = false;
        for (auto& bid : band_identifiers) 
          if (bid.In(qid))
            copy = true;
        if (not copy)
          band_identifiers.push_back(qid);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetRelamtLineMixingToMatches(
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfQuantumIdentifier& band_identifiers,
    const String& relaxation_type,
    const Verbosity&) {
  const auto type = Absorption::string2populationtype(relaxation_type);
  
  if (not Absorption::relaxationtype_relmat(type)) {
    std::ostringstream os;
    os << "The input type is: " << relaxation_type << "\n";
    os << "\tThis is not a relaxation matrix type\n";
    throw std::runtime_error(os.str());
  }
  
  for (Index qi = 0; qi < band_identifiers.nelem(); qi++) {
    const QuantumIdentifier& id = band_identifiers[qi];

    for (auto& lines: abs_lines_per_species) {
      for (auto& band: lines) {
        for (Index k=0; k<band.NumLines() and (band.Population() not_eq type); k++) {
          if (Absorption::id_in_line(band, id, k)) {
            band.Population(type);
          }
        }
      }
    }
  }
}

//static constexpr Index hartman_tran_type = 0;
//static constexpr Index linear_type = 1;


// Ignore function arguments if compiled without RELMAT support
#ifdef ENABLE_RELMAT
#define _UU_
#else
#define _UU_ _U_
#endif

/* Workspace method: Doxygen documentation will be auto-generated */
// void abs_xsec_per_speciesAddLineMixedBands(  // WS Output:
//     ArrayOfMatrix& abs_xsec_per_species _UU_,
//     ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx _UU_,
//     ArrayOfArrayOfMatrix& relmat_per_band _UU_,
//     // WS Input:
//     const ArrayOfArrayOfLineRecord& abs_lines_per_band _UU_,
//     const ArrayOfArrayOfSpeciesTag& abs_species_per_band _UU_,
//     const ArrayOfQuantumIdentifier& band_identifiers _UU_,
//     const ArrayOfArrayOfSpeciesTag& abs_species _UU_,
//     const SpeciesAuxData& isotopologue_ratios _UU_,
//     const SpeciesAuxData& partition_functions _UU_,
//     const ArrayOfRetrievalQuantity& jacobian_quantities _UU_,
//     const Vector& f_grid _UU_,
//     const Vector& abs_p _UU_,
//     const Vector& abs_t _UU_,
//     const Numeric& lm_p_lim _UU_,
//     const ArrayOfIndex& relmat_type_per_band _UU_,
//     const Index& wigner_initialized _UU_,
//     const Numeric& pressure_rule_limit _UU_,
//     const Index& write_relmat_per_band _UU_,
//     const Index& error_handling _UU_,
//     const Index& order_of_linemixing _UU_,
//     const Index& use_adiabatic_factor _UU_,
//     const Verbosity& verbosity _UU_)
// #ifdef ENABLE_RELMAT
// {
//   CREATE_OUT3;
//   using global_data::species_data;
//   using global_data::SpeciesMap;
// 
//   // Physical constants
//   extern const Numeric SPEED_OF_LIGHT;
//   extern const Numeric DOPPLER_CONST;
//   extern const Numeric ATM2PA;
// 
//   // HITRAN to ARTS constants
//   const Numeric doppler_const = DOPPLER_CONST;
//   static const Numeric w2Hz = SPEED_OF_LIGHT * 1E2;
//   static const Numeric lower_energy_const = wavenumber_to_joule(1.0);
//   static const Numeric I0_hi2arts = 1E-2 * SPEED_OF_LIGHT;
//   static const Numeric gamma_hi2arts = w2Hz / ATM2PA;
// 
//   // size of atmosphere and input/output
//   const Index nps = abs_p.nelem();
//   const Index nts = abs_t.nelem();
//   const Index nf = f_grid.nelem();
//   const Index nspecies = abs_species.nelem();
// 
//   // Relmat constants
//   const Numeric relmat_T0 = 296.0;
// 
//   // Jacobian constants
//   const ArrayOfIndex jacobian_quantities_position =
//       equivalent_propmattype_indexes(jacobian_quantities);
// 
//   // Test that wigner is wigner is initialized
//   if (not wigner_initialized)
//     throw std::runtime_error("Run InitWigner6() before this function...");
// 
//   // These should be identical
//   if (nps not_eq nts)
//     throw std::runtime_error(
//         "Different lengths of atmospheric input than expected.");
// 
//   if (nspecies == 0 or nspecies not_eq abs_xsec_per_species.nelem())
//     throw std::runtime_error(
//         "Absorption species and abs_xsec_per_species are not from same source.");
//   else {
//     for (Index i = 0; i < nspecies; i++) {
//       if (abs_xsec_per_species[i].ncols() not_eq nts)
//         throw std::runtime_error(
//             "Unexpected size of xsec matrix not matching abs_t and abs_p length.");
//       else if (abs_xsec_per_species[i].nrows() not_eq nf)
//         throw std::runtime_error(
//             "Unexpected size of xsec matrix not matching f_grid length.");
//       if (relmat_type_per_band.nelem() not_eq abs_lines_per_band.nelem())
//         throw std::runtime_error(
//             "Mismatching relmat_type_per_band and abs_lines_per_band.\n");
//     }
//   }
// 
//   if (abs_lines_per_band.nelem() not_eq band_identifiers.nelem())
//     throw std::runtime_error("Mismatch between band_identifiers and bands\n");
// 
//   supports_relaxation_matrix(
//       jacobian_quantities);  // It should support temperature and wind, and maybe even the line mixing parameters as an error vector
// 
//   const Index nbands = abs_lines_per_band.nelem();
// 
//   if (write_relmat_per_band) {
//     relmat_per_band.resize(nps);
//     for (Index ip = 0; ip < nps; ip++) relmat_per_band[ip].resize(nbands);
//   }
// 
//   if (nbands not_eq abs_species_per_band.nelem())
//     throw std::runtime_error(
//         "Error in definition of the bands.  Mismatching length of *_per_band arrays.");
// 
//   // Setting up thermal bath:  only in simplistic air for now
//   // This means: 21% O2 and 79% N2
//   long number_of_perturbers = 2;
//   ArrayOfIndex molecule_code_perturber(number_of_perturbers);
//   ArrayOfIndex iso_code_perturber(number_of_perturbers);
//   Vector perturber_mass(number_of_perturbers);
//   Vector vmr(number_of_perturbers);
// 
//   // Setup air as background gas for now...
//   {
//     vmr[0] = 0.2095;
//     const SpeciesRecord& O2 = species_data[SpeciesMap.find("O2")->second];
//     const IsotopologueRecord& O2_66 = O2.Isotopologue()[0];
//     const Index O2_66_hitran_tag = O2_66.HitranTag();
//     iso_code_perturber[0] = O2_66_hitran_tag % 10;
//     molecule_code_perturber[0] = (O2_66_hitran_tag) / 10;
//     perturber_mass[0] = O2_66.Mass();
// 
//     vmr[1] = 1.0 - 0.2095;
//     const SpeciesRecord& N2 = species_data[SpeciesMap.find("N2")->second];
//     const IsotopologueRecord& N2_44 = N2.Isotopologue()[0];
//     const Index N2_44_hitran_tag = N2_44.HitranTag();
//     iso_code_perturber[1] = N2_44_hitran_tag % 10;
//     molecule_code_perturber[1] = (N2_44_hitran_tag) / 10;
//     perturber_mass[1] = N2_44.Mass();
//   }
// 
//   // Flags and rules
//   double tolerance_in_rule_nr2 = pressure_rule_limit;
//   bool bool_use_adiabatic_factor = use_adiabatic_factor;
// 
//   for (Index iband = 0; iband < nbands; iband++) {
//     // band pointer
//     const ArrayOfLineRecord& this_band = abs_lines_per_band[iband];
// 
//     long nlines = (long)this_band.nelem();
// 
//     // Worth doing anything?
//     if (nlines <= 0) {
//       continue;
//     }
// 
//     // Send in frequency range
//     Numeric fmin, fmax;
// 
//     // To store the xsec matrix we need to know where
//     Index this_species = -1;
// 
//     // Allocation of band specific inputs (types: to be used as fortran input)
//     long M = (this_band[0].IsotopologueData().HitranTag()) / 10;
//     long I = (this_band[0].IsotopologueData().HitranTag()) % 10;
//     ArrayOfIndex upper(4 * nlines);
//     ArrayOfIndex lower(4 * nlines);
//     ArrayOfIndex g_prime(nlines);
//     ArrayOfIndex g_double_prime(nlines);
//     Vector v0(nlines);
//     Vector S(nlines);
//     Vector gamma_air(nlines);
//     Vector delta_air(nlines);
//     Vector E_double_prime(nlines);
//     Vector n_air(nlines);
//     Numeric mass, iso_ratio;
// 
//     for (long iline = 0; iline < nlines; iline++) {
//       // Line data
//       const LineRecord& this_line = this_band[iline];
//       const IsotopologueRecord& this_iso = this_line.IsotopologueData();
// 
//       // For first line do something special with mass and name and such
//       if (iline == 0) {
//         // Isotopologue values
//         mass = this_iso.Mass();
//         String iso_name = this_iso.Name();
//         int isona;
//         extract(isona, iso_name, iso_name.nelem());
// 
//         const ArrayOfSpeciesTag& band_tags = abs_species_per_band[iband];
// 
//         // Finds the first species in abs_species_per_band that matches abs_species
//         bool this_one = false;
//         for (Index ispecies = 0; ispecies < nspecies; ispecies++) {
//           const ArrayOfSpeciesTag& species_tags = abs_species[ispecies];
//           const Index nbandtags = band_tags.nelem();
//           const Index nspeciestags = species_tags.nelem();
// 
//           // Test if there is
//           if (nbandtags not_eq nspeciestags) {
//             break;
//           }
//           for (Index itags = 0; itags < nspeciestags; itags++) {
//             if (band_tags[itags] == species_tags[itags]) {
//               this_one = true;
//               break;
//             }
//           }
// 
//           if (this_one) {
//             this_species = ispecies;
//             break;
//           }
//         }
// 
//         if (!this_one)
//           throw std::runtime_error(
//               "abs_species and abs_species_per_band disagrees"
//               " on absorption species");
// 
//         iso_ratio =
//             isotopologue_ratios
//                 .getParam(abs_lines_per_band[iband][iline].Species(),
//                           abs_lines_per_band[iband][iline].Isotopologue())[0]
//                 .data[0];
// 
//       } else {
//         if (mass not_eq this_iso.Mass()) {
//           throw std::runtime_error(
//               "There are lines of different Isotopologues"
//               " in abs_lines_per_band,");
//         }
//       }
// 
//       // Pressure broadening at relmat temperatures
//       if (not this_line.LineShapeModelHasAir()) {
//         std::ostringstream os;
//         os << "Line does not not have air broadening but this function only uses air broadening.\n";
//         os << "Its values are " << this_line.GetLineShapeModel();
//         throw std::runtime_error(os.str());
//       }
// 
//       // Ensure that temperatures are sufficiently close
//       if (1e-4 < abs(this_line.Ti0() - relmat_T0)) {
//         std::ostringstream os;
//         os << "Line is not of same standard temperature as relmat is expecting.\n";
//         os << "Expecting: " << relmat_T0 << " K.  Getting: " << this_line.Ti0()
//            << " K.";
//         throw std::runtime_error(os.str());
//       }
// 
//       gamma_air[iline] =
//           this_line.PressureBroadeningAirBroadeningAgam() / gamma_hi2arts;
//       delta_air[iline] = this_line.PressureBroadeningAirBroadeningPsf();
//       n_air[iline] = this_line.Nair();
// 
//       // Line information converted to relmat format --- i.e. to HITRAN format
//       v0[iline] = this_line.F() / w2Hz;
//       S[iline] = this_line.I0() / I0_hi2arts;
//       E_double_prime[iline] = this_line.Elow() / lower_energy_const;
//       g_prime[iline] =
//           (long)this_line.G_upper();  // NB:  Numeric to long... why?
//       g_double_prime[iline] =
//           (long)this_line.G_lower();  // NB:  Numeric to long... why?
// 
//       // Quantum numbers converted to relmat format, again Numeric/Rational to long... why?
//       Rational a;
// 
//       // l2 is for molecules like CO2
//       {
//         a = this_line.LowerQuantumNumbers()[QuantumNumberType::l2];
//         a.Simplify();
//         if (a.isUndefined())
//           lower[0 + 4 * iline] = -1;
//         else if (a.Denom() == 1)
//           lower[0 + 4 * iline] = (long)a.toIndex();
//         else
//           throw std::runtime_error("Half quantum numbers not supported in l2.");
//         a = this_line.UpperQuantumNumbers()[QuantumNumberType::l2];
//         a.Simplify();
//         if (a.isUndefined())
//           upper[0 + 4 * iline] = -1;
//         else if (a.Denom() == 1)
//           upper[0 + 4 * iline] = (long)a.toIndex();
//         else
//           throw std::runtime_error("Half quantum numbers not supported in l2.");
//       }
// 
//       // J is universally important for linear molecules
//       {
//         a = this_line.LowerQuantumNumbers()[QuantumNumberType::J];
//         a.Simplify();
//         if (a.isUndefined())
//           lower[1 + 4 * iline] = -1;
//         else if (a.Denom() == 1)
//           lower[1 + 4 * iline] = (long)a.toIndex();
//         else
//           throw std::runtime_error("Half quantum numbers not supported in J.");
//         a = this_line.UpperQuantumNumbers()[QuantumNumberType::J];
//         a.Simplify();
//         if (a.isUndefined())
//           upper[1 + 4 * iline] = -1;
//         else if (a.Denom() == 1)
//           upper[1 + 4 * iline] = (long)a.toIndex();
//         else
//           throw std::runtime_error("Half quantum numbers not supported in J.");
//       }
// 
//       // N is important for molecules with magnetic dipoles
//       {
//         a = this_line.LowerQuantumNumbers()[QuantumNumberType::N];
//         a.Simplify();
//         if (a.isUndefined())
//           lower[2 + 4 * iline] = -1;
//         else if (a.Denom() == 1)
//           lower[2 + 4 * iline] = (long)a.toIndex();
//         else
//           throw std::runtime_error("Half quantum numbers not supported in N.");
//         a = this_line.UpperQuantumNumbers()[QuantumNumberType::N];
//         a.Simplify();
//         if (a.isUndefined())
//           upper[2 + 4 * iline] = -1;
//         else if (a.Denom() == 1)
//           upper[2 + 4 * iline] = (long)a.toIndex();
//         else
//           throw std::runtime_error("Half quantum numbers not supported in N.");
//       }
// 
//       // S is important for molecules with magnetic dipoles
//       {
//         a = this_line.LowerQuantumNumbers()[QuantumNumberType::S];
//         a.Simplify();
//         if (a.isUndefined())
//           lower[3 + 4 * iline] = -1;
//         else if (a.Denom() == 1)
//           lower[3 + 4 * iline] = (long)a.toIndex();
//         else
//           throw std::runtime_error("Half quantum numbers not supported in S.");
//         a = this_line.UpperQuantumNumbers()[QuantumNumberType::S];
//         a.Simplify();
//         if (a.isUndefined())
//           upper[3 + 4 * iline] = -1;
//         else if (a.Denom() == 1)
//           upper[3 + 4 * iline] = (long)a.toIndex();
//         else
//           throw std::runtime_error("Half quantum numbers not supported in S.");
//       }
// 
//       // Set fmax and fmin.  Why do I need this again?
//       if (iline == 0) {
//         fmin = v0[0] - 1.0;
//         fmax = v0[0] + 1.0;
//       } else {
//         if (fmin > v0[iline]) fmin = v0[iline] - 1.0;
//         if (fmax < v0[iline]) fmax = v0[iline] + 1.0;
//       }
//     }
// 
//     Vector f0 = v0;
//     f0 *= w2Hz;
// 
//     for (Index ip = 0; ip < nps; ip++) {
//       // Temporary table
//       wig_temp_init(2 * int(wigner_initialized));
// 
//       // Information on the lines will be here after relmat is done
//       Matrix W(nlines, nlines);
//       Vector dipole(nlines);
//       Vector rhoT(nlines);
//       Vector Y(nlines), G(nlines), DV(nlines);
// 
//       Matrix W_dt;
//       Vector Y_dt, G_dt, DV_dt;
//       Vector dipole_dt;
//       Vector rhoT_dt;
//       if (do_temperature_jacobian(jacobian_quantities)) {
//         W_dt.resize(nlines, nlines);
//         dipole_dt.resize(nlines);
//         rhoT_dt.resize(nlines);
//         Y_dt.resize(nlines);
//         G_dt.resize(nlines);
//         DV_dt.resize(nlines);
//       }
// 
//       Vector psf(nlines), psf_dt;
//       if (do_temperature_jacobian(jacobian_quantities)) psf_dt.resize(nlines);
// 
//       if (abs_p[ip] > lm_p_lim or order_of_linemixing < 0) {
//         // Get partition function information
//         Numeric QT0;
//         Numeric QT;
//         partition_function(QT0,
//                            QT,
//                            abs_lines_per_band[iband][0].Ti0(),
//                            abs_t[ip],
//                            partition_functions.getParamType(
//                                abs_lines_per_band[iband][0].Species(),
//                                abs_lines_per_band[iband][0].Isotopologue()),
//                            partition_functions.getParam(
//                                abs_lines_per_band[iband][0].Species(),
//                                abs_lines_per_band[iband][0].Isotopologue()));
// 
//         // Cannot be constants for Fortran's sake
//         Numeric t;
//         t = abs_t[ip];
//         Numeric p;
//         p = abs_p[ip] / ATM2PA;  // HITRAN pressure unit is in atmospheres
// 
//         Index error_handling_type = error_handling;
// 
//         // If lm_p_lim < some-limit, the ordered approach instead of the full approach is used.
//         // This makes the computations work at low pressures where we need Voigt line shape
//         Index order_of_linemixing_type;
//         if (order_of_linemixing < 0 and not(abs_p[ip] > lm_p_lim))
//           order_of_linemixing_type = -order_of_linemixing;
//         else
//           order_of_linemixing_type = order_of_linemixing;
// 
//         // Calling Teresa's code
//         if (relmat_type_per_band[iband] == hartman_tran_type) {
//           arts_relmat_interface__hartmann_and_niro_type(
//               &nlines,
//               &fmin,
//               &fmax,
//               &M,
//               &I,
//               v0.get_c_array(),
//               S.get_c_array(),
//               gamma_air.get_c_array(),
//               E_double_prime.get_c_array(),
//               n_air.get_c_array(),
//               upper.data(),
//               lower.data(),
//               g_prime.data(),
//               g_double_prime.data(),
//               &t,
//               &p,
//               &QT,
//               &QT0,
//               &mass,
//               &number_of_perturbers,
//               molecule_code_perturber.data(),
//               iso_code_perturber.data(),
//               perturber_mass.get_c_array(),
//               vmr.get_c_array(),
//               &error_handling_type,
//               &order_of_linemixing_type,
//               &tolerance_in_rule_nr2,
//               &bool_use_adiabatic_factor,
//               W.get_c_array(),
//               dipole.get_c_array(),
//               rhoT.get_c_array(),
//               Y.get_c_array(),
//               G.get_c_array(),
//               DV.get_c_array());
//         }
// 
//         else if (relmat_type_per_band[iband] == linear_type) {
//           arts_relmat_interface__linear_type(&nlines,
//                                              &fmin,
//                                              &fmax,
//                                              &M,
//                                              &I,
//                                              v0.get_c_array(),
//                                              S.get_c_array(),
//                                              gamma_air.get_c_array(),
//                                              E_double_prime.get_c_array(),
//                                              n_air.get_c_array(),
//                                              upper.data(),
//                                              lower.data(),
//                                              g_prime.data(),
//                                              g_double_prime.data(),
//                                              &t,
//                                              &p,
//                                              &QT,
//                                              &QT0,
//                                              &mass,
//                                              &number_of_perturbers,
//                                              molecule_code_perturber.data(),
//                                              iso_code_perturber.data(),
//                                              perturber_mass.get_c_array(),
//                                              vmr.get_c_array(),
//                                              &error_handling_type,
//                                              &order_of_linemixing_type,
//                                              &tolerance_in_rule_nr2,
//                                              &bool_use_adiabatic_factor,
//                                              W.get_c_array(),
//                                              dipole.get_c_array(),
//                                              rhoT.get_c_array(),
//                                              Y.get_c_array(),
//                                              G.get_c_array(),
//                                              DV.get_c_array());
//         } else {
// #pragma omp critical
//           throw std::runtime_error(
//               "Unsupported relaxation matrix type encountered.\n");
//         }
// 
//         if (error_handling_type == 1) {
//           std::ostringstream os;
//           os << "Fatal error encountered in relmat calculations.  Check your input for sanity.\n"
//              << "\tTo check what relmat is doing, activate the debug flag of this code."
//              << "\tIdentity: " << band_identifiers[iband] << std::endl;
// #pragma omp critical
//           throw std::runtime_error(os.str());
//         } else if (error_handling_type == 2) {
//           std::ostringstream os;
//           os << "Band: " << band_identifiers[iband] << std::endl
//              << "Did not pass rule 1: you need more lines in this band. "
//              << "LBL without Line Mixing is performed\n"
//              << "Pressure: " << p << " atm. Temperature: " << t << " K";
// #pragma omp critical
//           out3 << os.str() << "\n";
//         } else if (error_handling_type == 3) {
//           std::ostringstream os;
//           os << "Band: " << band_identifiers[iband] << std::endl
//              << "Did not pass rule 2: the pressure check failed. "
//              << "LBL without Line Mixing is performed\n"
//              << "Pressure: " << p << " atm. Temperature: " << t << " K";
// 
// #pragma omp critical
//           out3 << os.str() << "\n";
//         } else if (error_handling_type == 4) {
//           std::ostringstream os;
//           os << "Band: " << band_identifiers[iband] << std::endl
//              << "Did not pass the sum rule: the band cannot be renormalized. "
//              << "LBL without Line Mixing is performed\n"
//              << "Pressure: " << p << " atm. Temperature: " << t << " K";
// 
// #pragma omp critical
//           out3 << os.str() << "\n";
//         }
// 
//         // Convert to SI-units
//         W *= w2Hz * 0.5;
//         dipole /= 100.0;  // sqrt(I0_hi2arts / w2Hz) = 1/100;
//         DV *= w2Hz;
// 
//         // The temperature derivatives are for now only possible to do with perturbations
//         if (do_temperature_jacobian(jacobian_quantities)) {
//           // Do not write debug information in this section... do not crash
//           Index e_tmp = -1;
// 
//           // Perturbed temperature
//           Numeric t_dt = t + temperature_perturbation(jacobian_quantities);
// 
//           Numeric QT_dt;
// 
//           // Perturbed partition functions
//           partition_function(QT0,
//                              QT_dt,
//                              abs_lines_per_band[iband][0].Ti0(),
//                              t_dt,
//                              partition_functions.getParamType(
//                                  abs_lines_per_band[iband][0].Species(),
//                                  abs_lines_per_band[iband][0].Isotopologue()),
//                              partition_functions.getParam(
//                                  abs_lines_per_band[iband][0].Species(),
//                                  abs_lines_per_band[iband][0].Isotopologue()));
// 
//           if (relmat_type_per_band[iband] == hartman_tran_type) {
//             arts_relmat_interface__hartmann_and_niro_type(
//                 &nlines,
//                 &fmin,
//                 &fmax,
//                 &M,
//                 &I,
//                 v0.get_c_array(),
//                 S.get_c_array(),
//                 gamma_air.get_c_array(),
//                 E_double_prime.get_c_array(),
//                 n_air.get_c_array(),
//                 upper.data(),
//                 lower.data(),
//                 g_prime.data(),
//                 g_double_prime.data(),
//                 &t_dt,
//                 &p,
//                 &QT_dt,
//                 &QT0,
//                 &mass,
//                 &number_of_perturbers,
//                 molecule_code_perturber.data(),
//                 iso_code_perturber.data(),
//                 perturber_mass.get_c_array(),
//                 vmr.get_c_array(),
//                 &e_tmp,
//                 &order_of_linemixing_type,
//                 &tolerance_in_rule_nr2,
//                 &bool_use_adiabatic_factor,
//                 W_dt.get_c_array(),
//                 dipole_dt.get_c_array(),
//                 rhoT_dt.get_c_array(),
//                 Y_dt.get_c_array(),
//                 G_dt.get_c_array(),
//                 DV_dt.get_c_array());
//           } else if (relmat_type_per_band[iband] == linear_type) {
//             arts_relmat_interface__linear_type(&nlines,
//                                                &fmin,
//                                                &fmax,
//                                                &M,
//                                                &I,
//                                                v0.get_c_array(),
//                                                S.get_c_array(),
//                                                gamma_air.get_c_array(),
//                                                E_double_prime.get_c_array(),
//                                                n_air.get_c_array(),
//                                                upper.data(),
//                                                lower.data(),
//                                                g_prime.data(),
//                                                g_double_prime.data(),
//                                                &t_dt,
//                                                &p,
//                                                &QT_dt,
//                                                &QT0,
//                                                &mass,
//                                                &number_of_perturbers,
//                                                molecule_code_perturber.data(),
//                                                iso_code_perturber.data(),
//                                                perturber_mass.get_c_array(),
//                                                vmr.get_c_array(),
//                                                &e_tmp,
//                                                &order_of_linemixing_type,
//                                                &tolerance_in_rule_nr2,
//                                                &bool_use_adiabatic_factor,
//                                                W_dt.get_c_array(),
//                                                dipole_dt.get_c_array(),
//                                                rhoT_dt.get_c_array(),
//                                                Y_dt.get_c_array(),
//                                                G_dt.get_c_array(),
//                                                DV_dt.get_c_array());
//           }
// 
// #pragma omp critical
//           if (e_tmp == 1) {
//             std::ostringstream os;
//             os << "Fatal error encountered in relmat calculations.  Check your input for sanity.\n"
//                << "\tTo check what relmat is doing, activate the debug flag of this code.\n"
//                << "\tYou passed normal calculations but failed in the calculations of the partial derivatives...\n"
//                << "\tIdentity (first line): " << M << " " << I << " "
//                << upper[0] << " " << upper[1] << " " << upper[2] << " "
//                << upper[3] << " " << lower[0] << " " << lower[1] << " "
//                << lower[2] << " " << lower[3];
//             throw std::runtime_error(os.str());
//           }
// 
//           // Convert to SI-units
//           W_dt *= w2Hz * 0.5;
//           dipole_dt /= 100.0;
//           DV_dt *= w2Hz;
//         }
// 
//         // Use the provided pressure shift  NOTE: this might be a bad idea
//         if (do_temperature_jacobian(jacobian_quantities)) {
//           const Numeric t_dt =
//               abs_t[ip] + temperature_perturbation(jacobian_quantities);
//           for (Index ii = 0; ii < nlines; ii++) {
//             psf[ii] = delta_air[ii] * abs_p[ip] *
//                       pow((relmat_T0 / abs_t[ip]),
//                           (Numeric)0.25 + (Numeric)1.5 * n_air[ii]);
//             psf_dt[ii] = delta_air[ii] * abs_p[ip] *
//                          pow((relmat_T0 / t_dt),
//                              (Numeric)0.25 + (Numeric)1.5 * n_air[ii]);
//           }
// 
//         } else {
//           for (Index ii = 0; ii < nlines; ii++) {
//             psf[ii] = delta_air[ii] * abs_p[ip] *
//                       pow((relmat_T0 / abs_t[ip]),
//                           (Numeric)0.25 + (Numeric)1.5 * n_air[ii]);
//           }
//         }
// 
//         out3 << "Adding to band " << iband + 1 << "/" << nbands << " with "
//              << nlines << " lines at T-P level " << ip + 1 << "/" << nps
//              << ". It "
//              << ((error_handling_type > 0) ? "fails some" : "passes all")
//              << " LM test.\n";
// 
//         // Using Rodrigues etal method
//         if (order_of_linemixing_type < 0) {
//           calculate_xsec_from_full_relmat(abs_xsec_per_species,
//                                           dabs_xsec_per_species_dx,
//                                           this_band,
//                                           jacobian_quantities,
//                                           jacobian_quantities_position,
//                                           W,
//                                           W_dt,
//                                           f0,
//                                           f_grid,
//                                           dipole,
//                                           rhoT,
//                                           rhoT_dt,
//                                           psf,
//                                           psf_dt,
//                                           abs_t[ip],
//                                           iso_ratio,
//                                           this_species,
//                                           ip,
//                                           nlines);
// 
//           if (write_relmat_per_band not_eq 0) relmat_per_band[ip][iband] = W;
//         } else {
//           ConstVectorView pressure_broadening = W.diagonal();
//           ConstVectorView pressure_broadening_dt =
//               do_temperature_jacobian(jacobian_quantities) ? W_dt.diagonal()
//                                                            : Vector(0);
// 
//           calculate_xsec_from_relmat_coefficients(abs_xsec_per_species,
//                                                   dabs_xsec_per_species_dx,
//                                                   jacobian_quantities,
//                                                   jacobian_quantities_position,
//                                                   pressure_broadening,
//                                                   pressure_broadening_dt,
//                                                   f0,
//                                                   f_grid,
//                                                   dipole,
//                                                   rhoT,
//                                                   rhoT_dt,
//                                                   psf,
//                                                   psf_dt,
//                                                   Y,
//                                                   Y_dt,
//                                                   G,
//                                                   G_dt,
//                                                   DV,
//                                                   DV_dt,
//                                                   abs_t[ip],
//                                                   mass,
//                                                   iso_ratio,
//                                                   this_species,
//                                                   ip,
//                                                   nlines);
// 
//           if (write_relmat_per_band not_eq 0) {
//             if (order_of_linemixing == 0) { /* do nothing here */
//             } else if (order_of_linemixing_type == 1) {
//               relmat_per_band[ip][iband].resize(1, nlines);
//               relmat_per_band[ip][iband](0, joker) = Y;
// 
//             } else if (order_of_linemixing == 2) {
//               relmat_per_band[ip][iband].resize(3, nlines);
//               relmat_per_band[ip][iband](0, joker) = Y;
//               relmat_per_band[ip][iband](1, joker) = G;
//               relmat_per_band[ip][iband](2, joker) = DV;
//             }
// 
//             else
//               relmat_per_band[ip][iband] = W;
//           }
//         }
//       } else {
//         Numeric QT = -1, QT0 = -1, part_ratio;
//         Eigen::VectorXcd F(nf), N(0);
//         Eigen::
//             Matrix<Complex, Eigen::Dynamic, Linefunctions::ExpectedDataSize()>
//                 data(nf, Linefunctions::ExpectedDataSize());
//         Eigen::MatrixXcd dF(nf, jacobian_quantities_position.nelem()), dN(0, 0);
//         const Numeric GD_div_F0 = doppler_const * sqrt(abs_t[ip] / mass);
// 
//         for (long iline = 0; iline < nlines; iline++) {
//           Numeric K1, K2, tmp1;
//           static const Vector v_tmp(0);
// 
//           GetLineScalingData(QT,
//                              QT0,
//                              part_ratio,
//                              K1,
//                              K2,
//                              tmp1,
//                              tmp1,
//                              partition_functions.getParamType(
//                                  abs_lines_per_band[iband][0].Species(),
//                                  abs_lines_per_band[iband][0].Isotopologue()),
//                              partition_functions.getParam(
//                                  abs_lines_per_band[iband][0].Species(),
//                                  abs_lines_per_band[iband][0].Isotopologue()),
//                              abs_t[ip],
//                              abs_lines_per_band[iband][iline].Ti0(),
//                              abs_lines_per_band[iband][iline].F(),
//                              abs_lines_per_band[iband][iline].Elow(),
//                              false,
//                              -1.0,
//                              -1.0,
//                              -1,
//                              -1,
//                              v_tmp);
// 
//           //abs_lines_per_band[iband][iline].SetAirPressureBroadening(W(iline, iline), psf[iline], abs_t[ip], abs_p[ip], 0.0);
//           // FIXME:  Update this entire section... the below is a temporary workaround
//           auto X = abs_lines_per_band[iband][iline].GetShapeParams(
//               abs_t[ip], abs_p[ip], {400e-6}, {{SpeciesTag("CO2")}});
//           W(iline, iline) = X.G0;
//           psf[iline] = X.D0;
// 
//           // TODO: Add derivatives here
// 
//           Linefunctions::set_voigt(F,
//                                    dF,
//                                    data,
//                                    MapToEigen(f_grid),
//                                    0.0,
//                                    0.0,
//                                    abs_lines_per_band[iband][iline].F(),
//                                    GD_div_F0,
//                                    X);  // Derivatives need to be added...
// 
//           Linefunctions::apply_linestrength_scaling_by_lte(
//               F,
//               dF,
//               N,
//               dN,
//               abs_lines_per_band[iband][iline],
//               abs_t[ip],
//               iso_ratio,
//               QT,
//               QT0);
// 
//           for (Index ii = 0; ii < nf; ii++) {
//             const Numeric& y = F[ii].real();
// #pragma omp atomic
//             abs_xsec_per_species[this_species](ii, ip) += y;
// 
//             for (Index jj = 0; jj < jacobian_quantities_position.nelem();
//                  jj++) {
//               const Numeric& dy_dx = dF(jj, ii).real();
// #pragma omp atomic
//               dabs_xsec_per_species_dx[this_species][jj](ii, ip) += dy_dx;
//             }
//           }
//         }
//       }
// 
//       // Remove the temporary table
//       wig_temp_free();
//     }
//   }
// }
// #else
// {
//   throw std::runtime_error(
//       "This version of ARTS was compiled without external line mixing support.");
// }
// #endif  //ENABLE_RELMAT
// #undef _UU_
// 
// 
/* Workspace method: Doxygen documentation will be auto-generated */
// void SetLineMixingCoefficinetsFromRelmat(  // WS Input And Output:
//     ArrayOfArrayOfLineRecord& abs_lines_per_band,
//     ArrayOfArrayOfMatrix& relmat_per_band,
//     // WS Input:
//     const ArrayOfArrayOfSpeciesTag& abs_species_per_band,
//     const ArrayOfQuantumIdentifier& band_identifiers,
//     const ArrayOfArrayOfSpeciesTag& abs_species,
//     const SpeciesAuxData& isotopologue_ratios,
//     const SpeciesAuxData& partition_functions,
//     const Numeric& rtp_pressure,
//     const Vector& abs_t,
//     const ArrayOfIndex& relmat_type_per_band,
//     const Index& wigner_initialized,
//     const Numeric& pressure_rule_limit,
//     const Index& error_handling,
//     const Index& order_of_linemixing,
//     const Index& use_adiabatic_factor,
//     const Verbosity& verbosity) {
//   const Index nband = abs_lines_per_band.nelem();
//   const Index nlevl = abs_t.nelem();
//   const Index ndata = (order_of_linemixing == 2) ? 3 : 1;
// 
//   // Other numbers are accepted by abs_xsec_per_speciesAddLineMixedBands
//   if (order_of_linemixing not_eq 1 and order_of_linemixing not_eq 2)
//     throw std::runtime_error(
//         "Need first or second order linemixing in order_of_linemixing");
// 
//   const static Vector f_grid(1, 1.0);
//   const static ArrayOfRetrievalQuantity jacobian_quantities(0);
//   const Vector abs_p(nlevl, rtp_pressure);
//   ArrayOfMatrix _tmp1, _tmp2;
//   ArrayOfArrayOfMatrix _tmp3, _tmp4;
//   ArrayOfIndex abs_species_active(abs_species.nelem());
//   for (Index i = 0; i < abs_species.nelem(); i++) abs_species_active[i] = i;
// 
//   // Initialize dummy variables
//   abs_xsec_per_speciesInit(_tmp1,
//                            _tmp2,
//                            _tmp3,
//                            _tmp4,
//                            abs_species,
//                            jacobian_quantities,
//                            abs_species_active,
//                            f_grid,
//                            abs_p,
//                            1,
//                            0,
//                            verbosity);
// 
//   // Main variable --- its size will be [pressure][band] after next function call
//   // Largest write:  relmat_per_band[ip][iband](0, joker) = Y;
//   // Largest write:  relmat_per_band[ip][iband](1, joker) = G;
//   // Largest write:  relmat_per_band[ip][iband](2, joker) = DV;
// 
//   abs_xsec_per_speciesAddLineMixedBands(_tmp1,
//                                         _tmp3,
//                                         relmat_per_band,
//                                         abs_lines_per_band,
//                                         abs_species_per_band,
//                                         band_identifiers,
//                                         abs_species,
//                                         isotopologue_ratios,
//                                         partition_functions,
//                                         jacobian_quantities,
//                                         f_grid,
//                                         abs_p,
//                                         abs_t,
//                                         0.0,
//                                         relmat_type_per_band,
//                                         wigner_initialized,
//                                         pressure_rule_limit,
//                                         1,
//                                         error_handling,
//                                         order_of_linemixing,
//                                         use_adiabatic_factor,
//                                         verbosity);
// 
//   for (Index iband = 0; iband < nband; iband++) {
//     // Compute vectors (copying data to make life easier)
//     ArrayOfVector data(ndata, Vector(nlevl, 0));
//     Vector delta(nlevl);
// 
//     const Index nline = abs_lines_per_band[iband].nelem();
// 
//     for (Index iline = 0; iline < nline; iline++) {
//       for (Index ilevl = 0; ilevl < nlevl; ilevl++)
//         for (Index idata = 0; idata < ndata; idata++)
//           data[idata][ilevl] = relmat_per_band[ilevl][iband](idata, iline);
// 
//       // data vector for holding the results before setting the line
//       Vector dvec(10, 0);
// 
//       /* dvec structure:
//        * 0 :  y0
//        * 1 :  y1
//        * 2 :  g0
//        * 3 :  g1
//        * 4 :  d0
//        * 5 :  d1
//        * 6 :  T0
//        * 7 :  n
//        * 8 :  2*n
//        * 9 :  2*n
//        * 
//        * Note that while all n and T0 are known from the line catalog, for legacy purpose this methods needs to set them
//        * These legacy reasons are bad to have here because they might create the confusion that you can change these and
//        * still be OK with the calculations...  regardless though, the method using these numbers is worse than direct 
//        * line mixing so some errors should be expected
//        * 
//        * Inteded computations:
//        * 
//        * Fit data to 
//        * 
//        * X  =  (F0 + F1 (T0/T - 1)) ((T0 / T)^n  * P)^k,
//        * 
//        * where F is any of y, g, d above, and T is temperature, and P is the pressure.  k is 1 for y but 2 for the 
//        * others.  n is the air pressure broadening parameter  (FIXME:  generalize to other planets)
//        */
// 
//       const Numeric T0 = abs_lines_per_band[iband][iline].Ti0();
//       const Numeric n = abs_lines_per_band[iband][iline]
//                             .PressureBroadeningAirBroadeningNair();
// 
//       dvec[6] = T0;
//       dvec[7] = n;
//       dvec[8] = 2 * n;
//       dvec[9] = 2 * n;
// 
//       // Fitting parameters for computations solving A(T, P) dC = data(P, T)-f(C, P, T); C += dC;
//       Vector C(2), dC(2);
//       Matrix A(nlevl, 2);
// 
//       const static Index max_loop_count = 10;
//       const static Numeric res_limit = 1e-30;
//       Index pos = 0;
//       bool squared = false;
// 
//       // Do similar fitting for all data vectors
//       for (auto& v : data) {
//         // Take a value from the center of the data and use it as a starting point for the fitting
//         C[0] = v[nlevl / 2] / rtp_pressure / (squared ? rtp_pressure : 1.0);
//         C[1] = v[nlevl / 2] / rtp_pressure / (squared ? rtp_pressure : 1.0);
// 
//         Numeric res;
//         Index loop_count = 0;
//         do {
//           for (Index ilevl = 0; ilevl < nlevl; ilevl++) {
//             const Numeric theta = T0 / abs_t[ilevl];
//             const Numeric theta_n = pow(theta, n);
//             const Numeric TP = theta_n * rtp_pressure *
//                                (squared ? theta_n * rtp_pressure : 1.0);
//             A(ilevl, 0) = TP;
//             A(ilevl, 1) = (theta - 1.0) * TP;
//             const Numeric f = C[0] * TP + C[1] * (theta - 1.0) * TP;
//             delta[ilevl] = v[ilevl] - f;
//           }
//           res = lsf(dC, A, delta);
//           C += dC;
//           loop_count++;
//         } while (res > res_limit and loop_count < max_loop_count);
// 
//         dvec[pos] = C[0];
//         pos++;
//         dvec[pos] = C[1];
//         pos++;
// 
//         // Only the first loop is not squared
//         squared = true;
//       }
// 
//       // Overwrite the linemixing data NOTE: Conversion to modern LM-format... y0 y1 ny, g0...dv1 ndv
//       abs_lines_per_band[iband][iline].SetLineMixing2SecondOrderData(
//           Vector({dvec[0],
//                   dvec[1],
//                   dvec[7],
//                   dvec[2],
//                   dvec[3],
//                   dvec[8],
//                   dvec[4],
//                   dvec[5],
//                   dvec[9]}));
//     }
//   }
// }

/* Workspace method: Doxygen documentation will be auto-generated */
void PrintSelfLineMixingStatus(
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Numeric& T,
    const Verbosity& verbosity) {
  CREATE_OUT0;
  
  constexpr QuantumNumberType J = QuantumNumberType::J;
  
  for (auto& lines: abs_lines_per_species) {
    for (auto& band: lines) {
      out0 << "BAND: " << band.QuantumIdentity() << "\n"
           << "JF\tJI\tY(" << T << ")\tG(" << T << ")\tDV(" << T
           << ")\tF0\tE0\n";
      for (Index k=0; k<band.NumLines(); k++) {
        const Vector vmrs(band.BroadeningSpecies().nelem(), 1);
        out0 << band.LowerQuantumNumber(k, J) << '\t'
             << band.UpperQuantumNumber(k, J) << '\t';
        auto x = band.ShapeParameters(k, T, 101325, vmrs);
        out0 << x.Y << '\t' << x.G << '\t' << x.DV << '\t' << band.F0(k) << '\t'
             << band.E0(k) << '\n';
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void relmat_per_bandInAir(ArrayOfArrayOfMatrix& relmat_per_band,
                          const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                          const SpeciesAuxData& partition_functions,
                          const Index& wigner_initialized,
                          const Vector& temperatures,
                          const Verbosity&) try {
  const Index n = nelem(abs_lines_per_species);
  const Index tsize = temperatures.nelem();

  // Ensure increasing temperatures
  for (auto i = 1; i < tsize; i++)
    if (temperatures[i] <= temperatures[i - 1])
      throw "Must have strictly increasing GIN temperatures";

  // Size of relaxation matrix
  relmat_per_band.resize(n);
  for (auto& r : relmat_per_band) r.resize(tsize);

  Index ib=0;
  for (Index i=0; i<abs_lines_per_species.nelem(); i++) {
    for (Index j=0; j<abs_lines_per_species[i].nelem(); j++) {
      for (Index k=0; k<tsize; k++) {
        if (Absorption::relaxationtype_relmat(abs_lines_per_species[i][j].Population()))
          relmatInAir(relmat_per_band[ib][k],
                      abs_lines_per_species[i][j],
                      partition_functions,
                      wigner_initialized,
                      temperatures[k]);
        ib++;
      }
    }
  }
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *relmat_per_bandInAir*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *relmat_per_bandInAir*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLineMixingFromRelmat(
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfArrayOfMatrix& relmat_per_band,
    const SpeciesAuxData& partition_functions,
    const Vector& temperatures,
    const String& linemixing_type,
    const Index& do_g,
    const Index& do_dv,
    const Verbosity&) try {
  const Index lsize = nelem(abs_lines_per_species);
  auto tsize = temperatures.nelem();

  if (relmat_per_band.nelem() not_eq nelem(abs_lines_per_species))
    throw "Mismatch between number of bands in abs_lines_per_band and relmat_per_band";

  for (auto j = 0; j < lsize; j++) {
    if (relmat_per_band[j].nelem() not_eq temperatures.nelem())
      throw "Mismatch between number of temperatures in relmat_per_band and GIN temperatures";

    for (auto i = 0; i < tsize; i++) {
      if (relmat_per_band[j][i].nrows() not_eq relmat_per_band[j][i].ncols())
        throw "Non-square relaxation matrix; this invalidates relmat_per_band";
    }
  }

  Index iline=0;
  for (auto& lines: abs_lines_per_species) {
    for (auto& band: lines) {
      const Index n = band.NumLines();
      if (n < 1) continue;
      if (not Absorption::relaxationtype_relmat(band.Population())) continue;
      
      Matrix Y(tsize, n, 0);
      Matrix G(tsize, n, 0);
      Matrix DV(tsize, n, 0);
      Vector d0 = dipole_vector(band, partition_functions);
      
      for (auto j = 0; j < tsize; j++) {
        Y(j, joker) = rosenkranz_first_order(band, relmat_per_band[iline][j], d0);
        G(j, joker) =
        rosenkranz_scaling_second_order(band, relmat_per_band[iline][j], d0);
        DV(j, joker) =
        rosenkranz_shifting_second_order(band, relmat_per_band[iline][j]);
      }
      
      if (linemixing_type == "AER") {
        if (tsize not_eq 4 or temperatures[0] not_eq 200 or
          temperatures[1] not_eq 250 or temperatures[2] not_eq 296 or
          temperatures[3] not_eq 340)
          throw "Bad GIN temperatures.  Must be 4-long vector of [200, 250, 296, 340] for AER type interpolation.";
        
        Vector data(12, 0);
        data[0] = 200;
        data[1] = 250;
        data[2] = 296;
        data[3] = 340;
        for (auto k = 0; k < n; k++) {
          data[Range(4, 4)] = Y(joker, k);
          if (do_g)
            data[Range(8, 4)] = G(joker, k);
          band.Line(k).SetLineMixing2AER(data);
        }
      } else if (linemixing_type == "LM2") {
        Vector data(9, 0);
        
        for (auto k = 0; k < n; k++) {
          const Numeric exp = band.Line(k).LineShape().Data().back().G0().X1;
          
          auto X = compute_2nd_order_lm_coeff(
            Y(joker, k), temperatures, exp, band.T0());
          data[0] = X.y0;
          data[1] = X.y1;
          data[2] = exp;
          
          if (do_g) {
            X = compute_2nd_order_lm_coeff(
              G(joker, k), temperatures, 2 * exp, band.T0());
            data[3] = X.y0;
            data[4] = X.y1;
            data[5] = 2 * exp;
          }
          
          if (do_dv) {
            X = compute_2nd_order_lm_coeff(
              DV(joker, k), temperatures, 2 * exp, band.T0());
            data[6] = X.y0;
            data[7] = X.y1;
            data[8] = 2 * exp;
          }
          band.Line(k).SetLineMixing2SecondOrderData(data);
        }
      } else {
        throw "Cannot interpret type of line mixing";
      }
      
      iline++;
    }
  }
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *abs_lines_per_bandSetLineMixingFromRelmat*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *abs_lines_per_bandSetLineMixingFromRelmat*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddLineMixedLines(  // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    // WS Input:
    const Vector& f_grid,
    const Vector& abs_p,
    const Vector& abs_t,
    const ArrayOfArrayOfMatrix& relmat_per_band,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const SpeciesAuxData& isotopologue_ratios,
    const SpeciesAuxData& partition_functions,
    const Verbosity&) try {
  const Index nbands = nelem(abs_lines_per_species);
  const Index nf = f_grid.nelem();
  const Index np = abs_t.nelem();

  if (relmat_per_band.nelem() not_eq nbands)
    throw "Bad sizes of relaxation matrix.  Must be flat.";
  
  ComplexVector F(nf);
  for (Index ip = 0; ip < np; ip++) {
    Index iband = 0;
    for (Index iline=0; iline<abs_lines_per_species.nelem(); iline++) {
      for (auto& band: abs_lines_per_species[iline]) {
        if (not Absorption::relaxationtype_relmat(band.Population())) continue;
        
        const Matrix& W = relmat_per_band[iband][ip];
        const Index N = band.NumLines();
        Eigen::MatrixXcd M(N, N);
        for (auto i1 = 0; i1 < N; i1++) {
          for (auto i2 = 0; i2 < N; i2++) {
            if (i1 not_eq i2) {
              M(i1, i2) = Complex(0, abs_p[ip]) * W(i1, i2);
            } else {
              M(i1, i2) = band.F0(i1) + Complex(0, abs_p[ip]) * W(i1, i2);
            }
          }
        }
        
        const Vector population =
        population_density_vector(band, partition_functions, abs_t[ip]);
        const Vector dipole = dipole_vector(band, partition_functions);
        const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> decM(M, true);
        const auto& D = decM.eigenvalues();
        const ComplexVector B =
        equivalent_linestrengths(population, dipole, decM);
        
        F = 0;
        for (Index il = 0; il < N; il++) {
          const Numeric gammaD =
          Linefunctions::DopplerConstant(abs_t[ip], band.SpeciesMass()) *
          D[il].real();
          for (auto iv = 0; iv < nf; iv++) {
            const Complex z = (f_grid[iv] - conj(D[il])) / gammaD;
            const Complex w = Faddeeva::w(z);
            const Complex zm = (f_grid[iv] + D[il]) / gammaD;
            const Complex wm = Faddeeva::w(zm);
            
            F[iv] +=
            (Constant::inv_sqrt_pi / gammaD) * (w * conj(B[il]) + wm * B[il]);
          }
        }
        
        for (Index iv = 0; iv < nf; iv++) {
          abs_xsec_per_species[iline](iv, ip) +=
          isotopologue_ratios.getIsotopologueRatio(band.QuantumIdentity()) * F[iv].real();
        }
        
        iband++;
      }
    }
  }
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *abs_xsec_per_speciesAddLineMixedLines*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *abs_xsec_per_speciesAddLineMixedLines*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddLineMixedLinesInAir(  // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    // WS Input:
    const Vector& f_grid,
    const Vector& abs_p,
    const Vector& abs_t,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const SpeciesAuxData& isotopologue_ratios,
    const SpeciesAuxData& partition_functions,
    const Index& wigner_initialized,
    const Index& minimum_line_count,
    const Verbosity& verbosity) try {
  CREATE_OUT2;
  const auto nf = f_grid.nelem();
  const auto np = abs_t.nelem();

  ComplexVector F(nf);
  const ArrayOfArrayOfSpeciesTag pseudo_spec({ArrayOfSpeciesTag(1, SpeciesTag("O2-66")),
                                              ArrayOfSpeciesTag(1, SpeciesTag("N2-44"))});
  const Vector pseudo_vmrs({0.21, 0.79});
  
  for (Index ispec=0; ispec<abs_lines_per_species.nelem(); ispec++) {
    Matrix& sum_xsec = abs_xsec_per_species[ispec];
    
    for (const AbsorptionLines& band: abs_lines_per_species[ispec]) {
      if (not Absorption::relaxationtype_relmat(band.Population())) continue;
      
      const Index N=band.NumLines();
      Vector vmrs = LineShape::vmrs(pseudo_vmrs, pseudo_spec, band.QuantumIdentity(), band.BroadeningSpecies(), band.Self(), band.Bath(), band.LineShapeType());
      
      for (auto ip = 0; ip < np; ip++) {
        Eigen::MatrixXcd M(N, N);
        if (minimum_line_count > N) {
          Matrix W;
          relmatInAir(W, band,
                      partition_functions,
                      wigner_initialized,
                      abs_t[ip]);
          for (auto i1 = 0; i1 < N; i1++) {
            auto x = band.ShapeParameters(i1, abs_t[ip], abs_p[ip], vmrs);
            for (auto i2 = 0; i2 < N; i2++) {
              if (i1 not_eq i2)
                M(i1, i2) = +Complex(0, abs_p[ip]) * W(i1, i2);
              else
                M(i1, i2) =
                    band.F0(i1) + x.D0 + Complex(0, abs_p[ip]) * W(i1, i2);
            }
          }
        } else {
          M.setZero();
          for (Index il = 0; il < N; il++) {
            auto x = band.ShapeParameters(il, abs_t[ip], abs_p[ip], vmrs);
            M(il, il) = Complex(0, x.G0);
          }
        }
        
        const Vector population =
        population_density_vector(band, partition_functions, abs_t[ip]);
        const Vector dipole = dipole_vector(band, partition_functions);
        const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> decM(M, true);
        const auto& D = decM.eigenvalues();
        const ComplexVector B =
        equivalent_linestrengths(population, dipole, decM);
        
        F = 0;
        for (Index il = 0; il < N; il++) {
          const Numeric gammaD =
          Linefunctions::DopplerConstant(abs_t[ip], band.SpeciesMass()) *
          D[il].real();
          for (auto iv = 0; iv < nf; iv++) {
            const Complex z = (f_grid[iv] - conj(D[il])) / gammaD;
            const Complex w = Faddeeva::w(z);
            const Complex zm = (f_grid[iv] + D[il]) / gammaD;
            const Complex wm = Faddeeva::w(zm);
            
            F[iv] += (Constant::inv_sqrt_pi / gammaD) * (w * conj(B[il]) + wm * B[il]);
          }
        }
        
        for (Index iv = 0; iv < nf; iv++) {
          sum_xsec(iv, ip) +=
          isotopologue_ratios.getIsotopologueRatio(band.QuantumIdentity()) * F[iv].real();
        }
      }
    }
  }
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *abs_xsec_per_speciesAddLineMixedLinesInAir*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *abs_xsec_per_speciesAddLineMixedLinesInAir*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}
