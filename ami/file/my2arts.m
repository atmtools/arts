function atag = my2arts(mytag)
% convert Mytran tag to ARTS tag
% fuction necessary for read/wrire_mytran.m
% Format atag = my2arts(mytag)
%
% Input    mytag: MYTRAN tag (a number)
% Output   atag:  ARTS tag (a string)

switch mytag 
 
 case 11 
  atag = 'H2O-161';
 case 12 
  atag = 'H2O-181';
 case 13 
  atag = 'H2O-171';
 case 14 
  atag = 'H2O-162';
 case 21 
  atag = 'CO2-626';
 case 22 
  atag = 'CO2-636';
 case 23
  atag = 'CO2-628';
 case 24 
  atag = 'CO2-627';
 case 25 
  atag = 'CO2-638';
 case 26 
  atag = 'CO2-637';
 case 27 
  atag = 'CO2-828';
 case 28 
  atag = 'CO2-728';
 case 31 
  atag = 'O3-666';
 case 32 
  atag = 'O3-668';
 case 33 
  atag = 'O3-686';
 case 34 
  atag = 'O3-667';
 case 35 
  atag = 'O3-676';
 case 41 
  atag = 'N2O-446';
 case 42 
  atag = 'N2O-456';
 case 43 
  atag = 'N2O-546';
 case 44 
  atag = 'N2O-448';
 case 51 
  atag = 'CO-26';
 case 52 
  atag = 'CO-36';
 case 53 
  atag = 'CO-28';
 case 71 
  atag = 'O2-66';
 case 72 
  atag = 'O2-68';
 case 73 
  atag = 'O2-67';
 case 81 
  atag = 'NO-46';
 case 91 
  atag = 'SO2-626';
 case 101 
  atag = 'NO2-646';
 case 111 
  atag = 'NH3-4111';
 case 112 
  atag = 'NH3-5111';
 case 121 
  atag = 'HNO3-146';
 case 131 
  atag = 'OH-61';
 case 132 
  atag = 'OH-81';
 case 133 
  atag = 'OH-62';
 case 141 
  atag = 'HF-19';
 case 151 
  atag = 'HCl-15';
 case 152 
  atag = 'HCl-17';
 case 161 
  atag = 'HBr-19';
 case 162 
  atag = 'HBr-11';
 case 181 
  atag = 'ClO-56';
 case 182 
  atag = 'ClO-76';
 case 191 
  atag = 'OCS-622';
 case 192 
  atag = 'OCS-624';
 case 193 
  atag = 'OCS-632';
 case 194 
  atag = 'OCS-822';
 case 201 
  atag = 'H2CO-1126';
 case 202 
  atag = 'H2CO-1136';
 case 203 
  atag = 'H2CO-1128';
 case 211 
  atag = 'HOCl-165';
 case 212 
  atag = 'HOCl-167';
 case 231 
  atag = 'HCN-124';
 case 232 
  atag = 'HCN-134';
 case 233 
  atag = 'HCN-125';
 case 241 
  atag = 'CH3Cl-215';
 case 242 
  atag = 'CH3Cl-217';
 case 251 
  atag = 'H2O2-1661';
 case 281 
  atag = 'PH3-1111';
 case 291 
  atag = 'COF2-269';
 case 311 
  atag = 'H2S-121';
 case 321 
  atag = 'HCOOH-1261';
 case 331 
  atag = 'HO2-166';
 case 341 
  atag = 'O-6';
 case 351 
  atag = 'ClONO2-5646';
 case 352 
  atag = 'ClONO2-7646';
 case 431 
  atag = 'OClO-656';
 case 432 
  atag = 'OClO-676';
 case 401 
  atag = 'BrO-96';
 case 402 
  atag = 'BrO-16';
 case 481 
  atag = 'H2SO4-126';
 case 491 
  atag = 'Cl2O2-565';
 case 492 
  atag = 'Cl2O2-765';     
 otherwise
  display(['tag ', mytag, ' is not a valid tag']);
  break
end

