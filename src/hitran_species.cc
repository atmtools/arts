#include <map>

#include "isotopologues.h"
#include "hitran_species.h"

#include "species_tags.h"
#include "species_tags.h"

namespace Hitran {
using HitranMap = std::map<Index, std::map<char, std::pair<Index, Numeric>>>;

/** In 2012 the order if isotopologues were changed in HITRAN
 * 
 * This version takes that into account.  Note that the other species
 * are left as they were in the latest_molparam_map map at the time
 * of creating this file (2021-03-10)
 */
const HitranMap pre2012co2change_molparam_map{
    {1, {  // H2O
      {'1', {Species::find_species_index("H2O", "161"), 9.97317E-01}},
      {'2', {Species::find_species_index("H2O", "181"), 1.99983E-03}},
      {'3', {Species::find_species_index("H2O", "171"), 3.71884E-04}},
      {'4', {Species::find_species_index("H2O", "162"), 3.10693E-04}},
      {'5', {Species::find_species_index("H2O", "182"), 6.23003E-07}},
      {'6', {Species::find_species_index("H2O", "172"), 1.15853E-07}},
      {'7', {Species::find_species_index("H2O", "262"), 2.41970E-08}},
    }},
    {2, {  // CO2
      {'1', {Species::find_species_index("CO2", "626"), 9.84204E-01}},
      {'2', {Species::find_species_index("CO2", "636"), 1.10574E-02}},
      {'3', {Species::find_species_index("CO2", "628"), 3.94707E-03}},
      {'4', {Species::find_species_index("CO2", "627"), 7.33989E-04}},
      {'5', {Species::find_species_index("CO2", "638"), 4.43446E-05}},
      {'6', {Species::find_species_index("CO2", "637"), 8.24623E-06}},
      {'7', {Species::find_species_index("CO2", "828"), 3.95734E-06}},
      {'8', {Species::find_species_index("CO2", "827"), 1.47180E-06}},
      {'9', {Species::find_species_index("CO2", "838"), 4.44600E-08}},  // This is different for HITRAN2008 and earlier cf original map
    }},
    {3, {  // O3
      {'1', {Species::find_species_index("O3", "666"), 9.92901E-01}},
      {'2', {Species::find_species_index("O3", "668"), 3.98194E-03}},
      {'3', {Species::find_species_index("O3", "686"), 1.99097E-03}},
      {'4', {Species::find_species_index("O3", "667"), 7.40475E-04}},
      {'5', {Species::find_species_index("O3", "676"), 3.70237E-04}},
    }},
    {4, {  // N2O
      {'1', {Species::find_species_index("N2O", "446"), 9.90333E-01}},
      {'2', {Species::find_species_index("N2O", "456"), 3.64093E-03}},
      {'3', {Species::find_species_index("N2O", "546"), 3.64093E-03}},
      {'4', {Species::find_species_index("N2O", "448"), 1.98582E-03}},
      {'5', {Species::find_species_index("N2O", "447"), 3.69280E-04}},
    }},
    {5, {  // CO
      {'1', {Species::find_species_index("CO", "26"), 9.86544E-01}},
      {'2', {Species::find_species_index("CO", "36"), 1.10836E-02}},
      {'3', {Species::find_species_index("CO", "28"), 1.97822E-03}},
      {'4', {Species::find_species_index("CO", "27"), 3.67867E-04}},
      {'5', {Species::find_species_index("CO", "38"), 2.22250E-05}},
      {'6', {Species::find_species_index("CO", "37"), 4.13292E-06}},
    }},
    {6, {  // CH4
      {'1', {Species::find_species_index("CH4", "211"), 9.88274E-01}},
      {'2', {Species::find_species_index("CH4", "311"), 1.11031E-02}},
      {'3', {Species::find_species_index("CH4", "212"), 6.15751E-04}},
      {'4', {Species::find_species_index("CH4", "312"), 6.91785E-06}},
    }},
    {7, {  // O2
      {'1', {Species::find_species_index("O2", "66"), 9.95262E-01}},
      {'2', {Species::find_species_index("O2", "68"), 3.99141E-03}},
      {'3', {Species::find_species_index("O2", "67"), 7.42235E-04}},
    }},
    {8, {  // NO
      {'1', {Species::find_species_index("NO", "46"), 9.93974E-01}},
      {'2', {Species::find_species_index("NO", "56"), 3.65431E-03}},
      {'3', {Species::find_species_index("NO", "48"), 1.99312E-03}},
    }},
    {9, {  // SO2
      {'1', {Species::find_species_index("SO2", "626"), 9.45678E-01}},
      {'2', {Species::find_species_index("SO2", "646"), 4.19503E-02}},
    }},
    {10, {  // NO2
      {'1', {Species::find_species_index("NO2", "646"), 9.91616E-01}},
      {'2', {Species::find_species_index("NO2", "656"), 3.64564E-03}},
    }},
    {11, {  // NH3
      {'1', {Species::find_species_index("NH3", "4111"), 9.95872E-01}},
      {'2', {Species::find_species_index("NH3", "5111"), 3.66129E-03}},
    }},
    {12, {  // HNO3
      {'1', {Species::find_species_index("HNO3", "146"), 9.89110E-01}},
      {'2', {Species::find_species_index("HNO3", "156"), 3.63600E-03}},
    }},
    {13, {  // OH
      {'1', {Species::find_species_index("OH", "61"), 9.97473E-01}},
      {'2', {Species::find_species_index("OH", "81"), 2.00014E-03}},
      {'3', {Species::find_species_index("OH", "62"), 1.55371E-04}},
    }},
    {14, {  // HF
      {'1', {Species::find_species_index("HF", "19"), 9.99844E-01}},
      {'2', {Species::find_species_index("HF", "29"), 1.55741E-04}},
    }},
    {15, {  // HCl
      {'1', {Species::find_species_index("HCl", "15"), 7.57587E-01}},
      {'2', {Species::find_species_index("HCl", "17"), 2.42257E-01}},
      {'3', {Species::find_species_index("HCl", "25"), 1.18005E-04}},
      {'4', {Species::find_species_index("HCl", "27"), 3.77350E-05}},
    }},
    {16, {  // HBr
      {'1', {Species::find_species_index("HBr", "19"), 5.06781E-01}},
      {'2', {Species::find_species_index("HBr", "11"), 4.93063E-01}},
      {'3', {Species::find_species_index("HBr", "29"), 7.89384E-05}},
      {'4', {Species::find_species_index("HBr", "21"), 7.68016E-05}},
    }},
    {17, {  // HI
      {'1', {Species::find_species_index("HI", "17"), 9.99844E-01}},
      {'2', {Species::find_species_index("HI", "27"), 1.55741E-04}},
    }},
    {18, {  // ClO
      {'1', {Species::find_species_index("ClO", "56"), 7.55908E-01}},
      {'2', {Species::find_species_index("ClO", "76"), 2.41720E-01}},
    }},
    {19, {  // OCS
      {'1', {Species::find_species_index("OCS", "622"), 9.37395E-01}},
      {'2', {Species::find_species_index("OCS", "624"), 4.15828E-02}},
      {'3', {Species::find_species_index("OCS", "632"), 1.05315E-02}},
      {'4', {Species::find_species_index("OCS", "623"), 7.39908E-03}},
      {'5', {Species::find_species_index("OCS", "822"), 1.87967E-03}},
      {'6', {Species::find_species_index("OCS", "634"), 4.67508E-04}},
    }},
    {20, {  // H2CO
      {'1', {Species::find_species_index("H2CO", "126"), 9.86237E-01}},
      {'2', {Species::find_species_index("H2CO", "136"), 1.10802E-02}},
      {'3', {Species::find_species_index("H2CO", "128"), 1.97761E-03}},
    }},
    {21, {  // HOCl
      {'1', {Species::find_species_index("HOCl", "165"), 7.55790E-01}},
      {'2', {Species::find_species_index("HOCl", "167"), 2.41683E-01}},
    }},
    {22, {  // N2
      {'1', {Species::find_species_index("N2", "44"), 9.92687E-01}},
      {'2', {Species::find_species_index("N2", "45"), 7.47809E-03}},
    }},
    {23, {  // HCN
      {'1', {Species::find_species_index("HCN", "124"), 9.85114E-01}},
      {'2', {Species::find_species_index("HCN", "134"), 1.10676E-02}},
      {'3', {Species::find_species_index("HCN", "125"), 3.62174E-03}},
    }},
    {24, {  // CH3Cl
      {'1', {Species::find_species_index("CH3Cl", "215"), 7.48937E-01}},
      {'2', {Species::find_species_index("CH3Cl", "217"), 2.39491E-01}},
    }},
    {25, {  // H2O2
      {'1', {Species::find_species_index("H2O2", "1661"), 9.94952E-01}},
    }},
    {26, {  // C2H2
      {'1', {Species::find_species_index("C2H2", "1221"), 9.77599E-01}},
      {'2', {Species::find_species_index("C2H2", "1231"), 2.19663E-02}},
      {'3', {Species::find_species_index("C2H2", "1222"), 3.04550E-04}},
    }},
    {27, {  // C2H6
      {'1', {Species::find_species_index("C2H6", "1221"), 9.76990E-01}},
      {'2', {Species::find_species_index("C2H6", "1231"), 2.19526E-02}},
    }},
    {28, {  // PH3
      {'1', {Species::find_species_index("PH3", "1111"), 9.99533E-01}},
    }},
    {29, {  // COF2
      {'1', {Species::find_species_index("COF2", "269"), 9.86544E-01}},
      {'2', {Species::find_species_index("COF2", "369"), 1.10834E-02}},
    }},
    {30, {  // SF6
      {'1', {Species::find_species_index("SF6", "29"), .950180E+00}},
    }},
    {31, {  // H2S
      {'1', {Species::find_species_index("H2S", "121"), 9.49884E-01}},
      {'2', {Species::find_species_index("H2S", "141"), 4.21369E-02}},
      {'3', {Species::find_species_index("H2S", "131"), 7.49766E-03}},
    }},
    {32, {  // HCOOH
      {'1', {Species::find_species_index("HCOOH", "126"), 9.83898E-01}},
    }},
    {33, {  // HO2
      {'1', {Species::find_species_index("HO2", "166"), 9.95107E-01}},
    }},
    {34, {  // O
      {'1', {Species::find_species_index("O", "6"), 9.97628E-01}},
    }},
    {35, {  // ClONO2
      {'1', {Species::find_species_index("ClONO2", "5646"), .749570E+00}},
      {'2', {Species::find_species_index("ClONO2", "7646"), .239694E+00}},
    }},
    {36, {  // NO+
      {'1', {Species::find_species_index("NO+", "46"), 9.93974E-01}},
    }},
    {37, {  // HOBr
      {'1', {Species::find_species_index("HOBr", "169"), 5.05579E-01}},
      {'2', {Species::find_species_index("HOBr", "161"), 4.91894E-01}},
    }},
    {38, {  // C2H4
      {'1', {Species::find_species_index("C2H4", "221"), 9.77294E-01}},
      {'2', {Species::find_species_index("C2H4", "231"), 2.19595E-02}},
    }},
    {39, {  // CH3OH
      {'1', {Species::find_species_index("CH3OH", "2161"), 9.85930E-01}},
    }},
    {40, {  // CH3Br
      {'1', {Species::find_species_index("CH3Br", "219"), 5.00995E-01}},
      {'2', {Species::find_species_index("CH3Br", "211"), 4.87433E-01}},
    }},
    {41, {  // CH3CN
      {'1', {Species::find_species_index("CH3CN", "2124"), 9.73866E-01}},
    }},
    {42, {  // CF4
      {'1', {Species::find_species_index("CF4", "29"), 9.88890E-01}},
    }},
    {43, {  // C4H2
      {'1', {Species::find_species_index("C4H2", "2211"), 9.55998E-01}},
    }},
    {44, {  // HC3N
      {'1', {Species::find_species_index("HC3N", "12224"), 9.63346E-01}},
    }},
    {45, {  // H2
      {'1', {Species::find_species_index("H2", "11"), 9.99688E-01}},
      {'2', {Species::find_species_index("H2", "12"), 3.11432E-04}},
    }},
    {46, {  // CS
      {'1', {Species::find_species_index("CS", "22"), 9.39624E-01}},
      {'2', {Species::find_species_index("CS", "24"), 4.16817E-02}},
      {'3', {Species::find_species_index("CS", "32"), 1.05565E-02}},
      {'4', {Species::find_species_index("CS", "23"), 7.41668E-03}},
    }},
    {47, {  // SO3
      {'1', {Species::find_species_index("SO3", "26"), 9.43400E-01}},
    }},
    {48, {  // C2N2
      {'1', {Species::find_species_index("C2N2", "4224"), 9.70752E-01}},
    }},
    {49, {  // COCl2
      {'1', {Species::find_species_index("COCl2", "2655"), 5.66392E-01}},
      {'2', {Species::find_species_index("COCl2", "2657"), 3.62235E-01}},
    }},
    {53, {  // CS2
      {'1', {Species::find_species_index("CS2", "222"), 8.92811E-01}},
      {'2', {Species::find_species_index("CS2", "224"), 7.92600E-02}},
      {'3', {Species::find_species_index("CS2", "223"), 1.40940E-02}},
      {'4', {Species::find_species_index("CS2", "232"), 1.03100E-02}},
    }},
  };

/** The latest version of the HITRAN online molparam.txt file as a map
 * 
 * Note that ARTS does not use AFGL notation as HITRAN so several species has
 * to be changed manually when recreating this map.  Please keep the comments
 * of this change around so that it is easy to do this recreation.  Thanks!
 * 
 * To recreate this map, run the pyarts.hitran.gen_latest_molparam_map.  If more
 * species names mismatch between ARTS and HITRAN, do please add a comment to indicate
 * this when you update this variable
 */
const HitranMap latest_molparam_map{
  {1, {  // H2O
    {'1', {Species::find_species_index("H2O", "161"), 9.97317E-01}},
    {'2', {Species::find_species_index("H2O", "181"), 1.99983E-03}},
    {'3', {Species::find_species_index("H2O", "171"), 3.71884E-04}},
    {'4', {Species::find_species_index("H2O", "162"), 3.10693E-04}},
    {'5', {Species::find_species_index("H2O", "182"), 6.23003E-07}},
    {'6', {Species::find_species_index("H2O", "172"), 1.15853E-07}},
    {'7', {Species::find_species_index("H2O", "262"), 2.41970E-08}},
  }},
  {2, {  // CO2
    {'1', {Species::find_species_index("CO2", "626"), 9.84204E-01}},
    {'2', {Species::find_species_index("CO2", "636"), 1.10574E-02}},
    {'3', {Species::find_species_index("CO2", "628"), 3.94707E-03}},
    {'4', {Species::find_species_index("CO2", "627"), 7.33989E-04}},
    {'5', {Species::find_species_index("CO2", "638"), 4.43446E-05}},
    {'6', {Species::find_species_index("CO2", "637"), 8.24623E-06}},
    {'7', {Species::find_species_index("CO2", "828"), 3.95734E-06}},
    {'8', {Species::find_species_index("CO2", "827"), 1.47180E-06}},
    {'9', {Species::find_species_index("CO2", "727"), 1.36847E-07}},
    {'0', {Species::find_species_index("CO2", "838"), 4.44600E-08}},
    {'A', {Species::find_species_index("CO2", "837"), 1.65354E-08}},
    {'B', {Species::find_species_index("CO2", "737"), 1.53750E-09}},
  }},
  {3, {  // O3
    {'1', {Species::find_species_index("O3", "666"), 9.92901E-01}},
    {'2', {Species::find_species_index("O3", "668"), 3.98194E-03}},
    {'3', {Species::find_species_index("O3", "686"), 1.99097E-03}},
    {'4', {Species::find_species_index("O3", "667"), 7.40475E-04}},
    {'5', {Species::find_species_index("O3", "676"), 3.70237E-04}},
  }},
  {4, {  // N2O
    {'1', {Species::find_species_index("N2O", "446"), 9.90333E-01}},
    {'2', {Species::find_species_index("N2O", "456"), 3.64093E-03}},
    {'3', {Species::find_species_index("N2O", "546"), 3.64093E-03}},
    {'4', {Species::find_species_index("N2O", "448"), 1.98582E-03}},
    {'5', {Species::find_species_index("N2O", "447"), 3.69280E-04}},
  }},
  {5, {  // CO
    {'1', {Species::find_species_index("CO", "26"), 9.86544E-01}},
    {'2', {Species::find_species_index("CO", "36"), 1.10836E-02}},
    {'3', {Species::find_species_index("CO", "28"), 1.97822E-03}},
    {'4', {Species::find_species_index("CO", "27"), 3.67867E-04}},
    {'5', {Species::find_species_index("CO", "38"), 2.22250E-05}},
    {'6', {Species::find_species_index("CO", "37"), 4.13292E-06}},
  }},
  {6, {  // CH4
    {'1', {Species::find_species_index("CH4", "211"), 9.88274E-01}},
    {'2', {Species::find_species_index("CH4", "311"), 1.11031E-02}},
    {'3', {Species::find_species_index("CH4", "212"), 6.15751E-04}},
    {'4', {Species::find_species_index("CH4", "312"), 6.91785E-06}},
  }},
  {7, {  // O2
    {'1', {Species::find_species_index("O2", "66"), 9.95262E-01}},
    {'2', {Species::find_species_index("O2", "68"), 3.99141E-03}},
    {'3', {Species::find_species_index("O2", "67"), 7.42235E-04}},
  }},
  {8, {  // NO
    {'1', {Species::find_species_index("NO", "46"), 9.93974E-01}},
    {'2', {Species::find_species_index("NO", "56"), 3.65431E-03}},
    {'3', {Species::find_species_index("NO", "48"), 1.99312E-03}},
  }},
  {9, {  // SO2
    {'1', {Species::find_species_index("SO2", "626"), 9.45678E-01}},
    {'2', {Species::find_species_index("SO2", "646"), 4.19503E-02}},
  }},
  {10, {  // NO2
    {'1', {Species::find_species_index("NO2", "646"), 9.91616E-01}},
    {'2', {Species::find_species_index("NO2", "656"), 3.64564E-03}},
  }},
  {11, {  // NH3
    {'1', {Species::find_species_index("NH3", "4111"), 9.95872E-01}},
    {'2', {Species::find_species_index("NH3", "5111"), 3.66129E-03}},
  }},
  {12, {  // HNO3
    {'1', {Species::find_species_index("HNO3", "146"), 9.89110E-01}},
    {'2', {Species::find_species_index("HNO3", "156"), 3.63600E-03}},
  }},
  {13, {  // OH
    {'1', {Species::find_species_index("OH", "61"), 9.97473E-01}},
    {'2', {Species::find_species_index("OH", "81"), 2.00014E-03}},
    {'3', {Species::find_species_index("OH", "62"), 1.55371E-04}},
  }},
  {14, {  // HF
    {'1', {Species::find_species_index("HF", "19"), 9.99844E-01}},
    {'2', {Species::find_species_index("HF", "29"), 1.55741E-04}},
  }},
  {15, {  // HCl
    {'1', {Species::find_species_index("HCl", "15"), 7.57587E-01}},
    {'2', {Species::find_species_index("HCl", "17"), 2.42257E-01}},
    {'3', {Species::find_species_index("HCl", "25"), 1.18005E-04}},
    {'4', {Species::find_species_index("HCl", "27"), 3.77350E-05}},
  }},
  {16, {  // HBr
    {'1', {Species::find_species_index("HBr", "19"), 5.06781E-01}},
    {'2', {Species::find_species_index("HBr", "11"), 4.93063E-01}},
    {'3', {Species::find_species_index("HBr", "29"), 7.89384E-05}},
    {'4', {Species::find_species_index("HBr", "21"), 7.68016E-05}},
  }},
  {17, {  // HI
    {'1', {Species::find_species_index("HI", "17"), 9.99844E-01}},
    {'2', {Species::find_species_index("HI", "27"), 1.55741E-04}},
  }},
  {18, {  // ClO
    {'1', {Species::find_species_index("ClO", "56"), 7.55908E-01}},
    {'2', {Species::find_species_index("ClO", "76"), 2.41720E-01}},
  }},
  {19, {  // OCS
    {'1', {Species::find_species_index("OCS", "622"), 9.37395E-01}},
    {'2', {Species::find_species_index("OCS", "624"), 4.15828E-02}},
    {'3', {Species::find_species_index("OCS", "632"), 1.05315E-02}},
    {'4', {Species::find_species_index("OCS", "623"), 7.39908E-03}},
    {'5', {Species::find_species_index("OCS", "822"), 1.87967E-03}},
    {'6', {Species::find_species_index("OCS", "634"), 4.67508E-04}},
  }},
  {20, {  // H2CO
    {'1', {Species::find_species_index("H2CO", "126"), 9.86237E-01}},
    {'2', {Species::find_species_index("H2CO", "136"), 1.10802E-02}},
    {'3', {Species::find_species_index("H2CO", "128"), 1.97761E-03}},
  }},
  {21, {  // HOCl
    {'1', {Species::find_species_index("HOCl", "165"), 7.55790E-01}},
    {'2', {Species::find_species_index("HOCl", "167"), 2.41683E-01}},
  }},
  {22, {  // N2
    {'1', {Species::find_species_index("N2", "44"), 9.92687E-01}},
    {'2', {Species::find_species_index("N2", "45"), 7.47809E-03}},
  }},
  {23, {  // HCN
    {'1', {Species::find_species_index("HCN", "124"), 9.85114E-01}},
    {'2', {Species::find_species_index("HCN", "134"), 1.10676E-02}},
    {'3', {Species::find_species_index("HCN", "125"), 3.62174E-03}},
  }},
  {24, {  // CH3Cl
    {'1', {Species::find_species_index("CH3Cl", "215"), 7.48937E-01}},
    {'2', {Species::find_species_index("CH3Cl", "217"), 2.39491E-01}},
  }},
  {25, {  // H2O2
    {'1', {Species::find_species_index("H2O2", "1661"), 9.94952E-01}},
  }},
  {26, {  // C2H2
    {'1', {Species::find_species_index("C2H2", "1221"), 9.77599E-01}},
    {'2', {Species::find_species_index("C2H2", "1231"), 2.19663E-02}},
    {'3', {Species::find_species_index("C2H2", "1222"), 3.04550E-04}},
  }},
  {27, {  // C2H6
    {'1', {Species::find_species_index("C2H6", "1221"), 9.76990E-01}},
    {'2', {Species::find_species_index("C2H6", "1231"), 2.19526E-02}},
  }},
  {28, {  // PH3
    {'1', {Species::find_species_index("PH3", "1111"), 9.99533E-01}},
  }},
  {29, {  // COF2
    {'1', {Species::find_species_index("COF2", "269"), 9.86544E-01}},
    {'2', {Species::find_species_index("COF2", "369"), 1.10834E-02}},
  }},
  {30, {  // SF6
    {'1', {Species::find_species_index("SF6", "29"), .950180E+00}},
  }},
  {31, {  // H2S
    {'1', {Species::find_species_index("H2S", "121"), 9.49884E-01}},
    {'2', {Species::find_species_index("H2S", "141"), 4.21369E-02}},
    {'3', {Species::find_species_index("H2S", "131"), 7.49766E-03}},
  }},
  {32, {  // HCOOH
    {'1', {Species::find_species_index("HCOOH", "126"), 9.83898E-01}},
  }},
  {33, {  // HO2
    {'1', {Species::find_species_index("HO2", "166"), 9.95107E-01}},
  }},
  {34, {  // O
    {'1', {Species::find_species_index("O", "6"), 9.97628E-01}},
  }},
  {35, {  // ClONO2
    {'1', {Species::find_species_index("ClONO2", "5646"), .749570E+00}},
    {'2', {Species::find_species_index("ClONO2", "7646"), .239694E+00}},
  }},
  {36, {  // NO+
    {'1', {Species::find_species_index("NO+", "46"), 9.93974E-01}},
  }},
  {37, {  // HOBr
    {'1', {Species::find_species_index("HOBr", "169"), 5.05579E-01}},
    {'2', {Species::find_species_index("HOBr", "161"), 4.91894E-01}},
  }},
  {38, {  // C2H4
    {'1', {Species::find_species_index("C2H4", "221"), 9.77294E-01}},
    {'2', {Species::find_species_index("C2H4", "231"), 2.19595E-02}},
  }},
  {39, {  // CH3OH
    {'1', {Species::find_species_index("CH3OH", "2161"), 9.85930E-01}},
  }},
  {40, {  // CH3Br
    {'1', {Species::find_species_index("CH3Br", "219"), 5.00995E-01}},
    {'2', {Species::find_species_index("CH3Br", "211"), 4.87433E-01}},
  }},
  {41, {  // CH3CN
    {'1', {Species::find_species_index("CH3CN", "2124"), 9.73866E-01}},
  }},
  {42, {  // CF4
    {'1', {Species::find_species_index("CF4", "29"), 9.88890E-01}},
  }},
  {43, {  // C4H2
    {'1', {Species::find_species_index("C4H2", "2211"), 9.55998E-01}},
  }},
  {44, {  // HC3N
    {'1', {Species::find_species_index("HC3N", "12224"), 9.63346E-01}},
  }},
  {45, {  // H2
    {'1', {Species::find_species_index("H2", "11"), 9.99688E-01}},
    {'2', {Species::find_species_index("H2", "12"), 3.11432E-04}},
  }},
  {46, {  // CS
    {'1', {Species::find_species_index("CS", "22"), 9.39624E-01}},
    {'2', {Species::find_species_index("CS", "24"), 4.16817E-02}},
    {'3', {Species::find_species_index("CS", "32"), 1.05565E-02}},
    {'4', {Species::find_species_index("CS", "23"), 7.41668E-03}},
  }},
  {47, {  // SO3
    {'1', {Species::find_species_index("SO3", "26"), 9.43400E-01}},
  }},
  {48, {  // C2N2
    {'1', {Species::find_species_index("C2N2", "4224"), 9.70752E-01}},
  }},
  {49, {  // COCl2
    {'1', {Species::find_species_index("COCl2", "2655"), 5.66392E-01}},
    {'2', {Species::find_species_index("COCl2", "2657"), 3.62235E-01}},
  }},
  {53, {  // CS2
    {'1', {Species::find_species_index("CS2", "222"), 8.92811E-01}},
    {'2', {Species::find_species_index("CS2", "224"), 7.92600E-02}},
    {'3', {Species::find_species_index("CS2", "223"), 1.40940E-02}},
    {'4', {Species::find_species_index("CS2", "232"), 1.03100E-02}},
  }},
};

using OurHitranMap = std::map<Index, std::map<char, Species::IsotopeRecord>>;

/** Turns the string-map required at compile time into a species-map
 * to be used as a static runtime map
 */
OurHitranMap to_species_map(const HitranMap& string_map) {
  OurHitranMap species_map;
  for (auto& specs: string_map) {
    ARTS_ASSERT(specs.second.find('1') not_eq specs.second.cend(),
                "Must have species '1' in map")
    for (auto& isot: specs.second) {
      ARTS_ASSERT(isot.second.first >= 0, "Undefined species in ARTS found in HITRAN data")
      species_map[specs.first][isot.first] = Species::Isotopologues[isot.second.first];
    }
  }
  return species_map;
}

/** Selects the map and returns it.  Uses templates to avoid weird warnings */
template <Type t>
OurHitranMap select_hitran_map() {
  static_assert(good_enum(t), "Bad enum encountered, something is amiss");
  
  /*
  // Debug code to uncomment to check that the output is as expected
  auto x = to_species_map(latest_molparam_map);
  for (auto y: x) {
    std::cout << y.first << '\n';
    for (auto z: y.second) {
      std:: cout << z.first << ' ' << z.second << '\n';
    }
  }
  */
  
  if constexpr (t == Type::Newest) {
    return to_species_map(latest_molparam_map);
  } else if constexpr (t == Type::Pre2012CO2Change) {
    return to_species_map(pre2012co2change_molparam_map);
  } else {
    std::terminate();
  }
}

/** Returns the species if possible or throws an error if it cannot be found
 * 
 * @param[in] molnum The hitran molecular number
 * @param[in] isonum The hitran character representing an isotopologue
 * @return An ARTS identifier of the species' as a transition
 */
template <Type t> QuantumIdentifier from_mol_iso(Index molnum, char isonum) {
  static_assert(good_enum(t), "Bad enum encountered, something is amiss. enum is: ");
  
  // Map of all HITRAN species
  static const OurHitranMap hitmap = select_hitran_map<t>();
  
  // Search the map with pointers so that we can throw manually if something is bad
  if (auto species_list = hitmap.find(molnum); species_list not_eq hitmap.cend()) {
    if (auto species_info = species_list -> second.find(isonum); species_info not_eq species_list -> second.cend()) {
      return QuantumIdentifier(species_info -> second, Quantum::IdentifierType::Transition);
    } else {
      ARTS_USER_ERROR("Species ", molnum, " has no isotopologue ", isonum,
                      " in ARTS' HITRAN implementation.\n",
                      "(Species is ", Species::toShortName(hitmap.at(molnum).at('1').spec), ")\n"
                      "If you are using a new version of HITRAN that has added the isotopologue, please consider\n"
                      "contacting the ARTS developers so we can append the species to our list and make this work.\n"
                      "Note that you are calling this templated function as the ", t, " template")
    }
  } else {
    ARTS_USER_ERROR("Species ", molnum, " does not exist in ARTS' HITRAN implementation\n"
                    "If you are using a new version of HITRAN that has added the species, please consider\n"
                    "contacting the ARTS developers so we can append the species to our list and make this work.\n"
                    "Note that you are calling this templated function as the ", t, " template");
  }
}

QuantumIdentifier id_from_lookup(Index mol, char isochar, Type type) {
  switch (type) {
    case Type::Pre2012CO2Change:
      return from_mol_iso<Type::Pre2012CO2Change>(mol, isochar);
    case Type::Newest:
      return from_mol_iso<Type::Newest>(mol, isochar);
    case Type::FINAL: {/* leave last */}
  }
  return {};
}

Numeric ratio_from_lookup(Index mol, char isochar, Type type) {
  switch (type) {
    case Type::Pre2012CO2Change:
      return pre2012co2change_molparam_map.at(mol).at(isochar).second;
    case Type::Newest:
      return latest_molparam_map.at(mol).at(isochar).second;
    case Type::FINAL: {/* leave last */}
  }
  return {};
}

SpeciesIsotopologueRatios isotopologue_ratios_impl(const HitranMap& data) {
  SpeciesIsotopologueRatios out;
  for (auto& x: data) {
    for (auto& y: x.second) {
      out.data[y.second.first] = y.second.second;
    }
  }
  return out;
}

SpeciesIsotopologueRatios isotopologue_ratios(Type type) {
  switch (type) {
    case Type::Pre2012CO2Change:
      return isotopologue_ratios_impl(pre2012co2change_molparam_map);
    case Type::Newest:
      return isotopologue_ratios_impl(latest_molparam_map);
    case Type::FINAL: {/* leave last */}
  }
  return {};
}
}  // namespace Hitran
