;
; $Id: partition_function.pro,v 1.6 2001/01/10 15:29:08 axel Exp $
;
pro HAK, dummy, mesg=mesg
; NAME:
;       HAK
; PURPOSE:
;       HAK is a procedure that performs "Hit any key to continue".  It waits
;       for keyboard input, clears the type-ahead buffer and allows the
;       application to continue.


; Check for the message keyword.
;
if keyword_set(mesg) then begin
   sm = size(mesg)
   if sm(1) ne 7 then begin          ; Print default string
      print, 'Hit any key to continue...'
   endif else begin                  ; Print user-defined string
      print, mesg
   endelse
endif
;
; Wait for keyboard input before continuing (returning).
;
dumb = get_kbrd (1)
empty
repeat dumb = get_kbrd (0) until dumb eq ''
empty
;
return
end

FUNCTION calculate_partition_functions,jpl_part,vib
;; little helper function for the correction of the jpl partition
;; functions  
cm2t	= 1.43875 ; Conversion factor for energy levels from cm-1 to K

;; go through all partition function values and vib states

print, 'Correcting JPL partition functions.'

for i=0, n_elements(jpl_part[0,*])-1 do begin
    for j=0,n_elements(vib)-1 do begin
        a=-vib[j]*cm2t/jpl_part[0,i]
        if a lt -80 then jpl_part[1,i]=jpl_part[1,i] else $
          jpl_part[1,i]=jpl_part[1,i]/(1-exp(a))
    endfor
endfor
return, jpl_part
end


PRO correct_jpl_partition_functions,jpl_part,species
;; corrects for vibrational partition fuctions not included in the
;; jpl catalogue. Data are either from:
;; Patrick Eriksson 
;; Miriam von König
;; or taken from the JANAF Thermochemical Tables Third Edition,
;; M.W. Chase et. al., Journal of Physical and Chemical Reference
;; Data, Volume 14, 1985.

case species of 
    18003 : begin ; H2O main, JANAF
        vib=[1594.7,3651.1,3755.9] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    20003 : begin ; H2O isotope, taken H2O main values
        vib=[1594.7,3651.1,3755.9] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    19003 : begin ; H2O isotope, taken H2O main values
        vib=[1594.7,3651.1,3755.9] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    19002 : begin ; H2O isotope, taken H2O main values
        vib=[1594.7,3651.1,3755.9] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    21001 : begin ; H2O isotope, taken H2O main values
        vib=[1594.7,3651.1,3755.9] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    20001 : begin ; H2O isotope, taken H2O main values
        vib=[1594.7,3651.1,3755.9] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    46013 : begin ; CO2 isotope, taken CO2 main values from JANAF
        vib=[667.3, 1384.9, 2349.3] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    45012 : begin ; CO2 isotope, taken CO2 main values from JANAF
        vib=[667.3, 1384.9, 2349.3] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    48004 : begin ; O3 main, JANAF
        vib=[705, 1043, 1110]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    50004 : begin ; O3-asym18, Patrick Eriksson
        vib=[693.0] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    50003 : begin ; O3-sym18 Patrick Eriksson
        vib=[678.0] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    49002 : begin ; O3 isotope, taken O3 main from JANAF
        vib=[705, 1043, 1110] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    49001 : begin ; O3 isotope, taken O3 main from JANAF
        vib=[705, 1043, 1110] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    44004 : begin ; N2O main, JANAF
        vib=[588.8,588.8,1284.9,2223.8] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    45007 : begin ; N2O isotope, taken main N2O values
        vib=[588.8,588.8,1284.9,2223.8] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    45008 : begin ; N2O isotope, taken main N2O values
        vib=[588.8,588.8,1284.9,2223.8] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    46007 : begin ; N2O isotope, taken main N2O values
        vib=[588.8,588.8,1284.9,2223.8] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    28001 : begin ; CO main, Janssen - Rosenkranz
        vib=[2143.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    29001 : begin ; CO isotope, taken main
        vib=[2143.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    30001 : begin ; CO isotope, taken main
        vib=[2143.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    29006 : begin ; CO isotope, taken main
        vib=[2143.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    17003 : begin ; CH4 isotope, taken main CH4 JANAF values
        vib=[1306,1306,1306,1534,1534,2917,3019] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    32001 : begin ; O2 main, Patrick Eriksson
        vib=[1556.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    34001 : begin ; O2 isotope, taken O2 main values
        vib=[1556.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    33002 : begin ; O2 isotope, taken O2 main values
        vib=[1556.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

;;    30008 : begin ; NO main, no info found in JANAF
;;        vib=[] ; [cm^-1]
;;        jpl_part=calculate_partition_functions(jpl_part,vib)
;;    end

    64002 : begin ; SO2 main, JANAF
        vib=[517.7,1151.4,1361.8] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    66002 : begin ; SO2 isotope, taken SO2 main
        vib=[517.7,1151.4,1361.8] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    65001 : begin ; SO2 isotope, taken SO2 main
        vib=[517.7,1151.4,1361.8] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    66004 : begin ; SO2 isotope, taken SO2 main
        vib=[517.7,1151.4,1361.8] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    46006 : begin ; NO2 main, JANAF
        vib=[756.8,1357.8,1665.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    17002 : begin ; NH3 main, JANAF
        vib=[950, 1629, 1629, 3335, 3414, 3414] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    18002 : begin ; NH3 isotope, taken main values
        vib=[950, 1629, 1629, 3335, 3414, 3414 ] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    18004 : begin ; NH3 isotope, taken main values
        vib=[950, 1629, 1629, 3335, 3414, 3414 ] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

;    17002 : begin ; NH3  Patrick Eriksson
;        vib=[948.6] ; [cm^-1]
;        jpl_part=calculate_partition_functions(jpl_part,vib)
;    end

    63001 : begin ; HNO3 main, JANAF
        vib=[465,583,680,765,886,1320,1335,1710,3560] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
;    63001 : begin ; HNO3 main, Patrick Eriksson
;        vib=[458.2,579.0,646.8,763.2,878.6,912,1035,1103,1158,1218,$
;	      1226,1294,1303,1326,1335,1341,1368,1409,1710,3550] ; [cm^-1]
;        jpl_part=calculate_partition_functions(jpl_part,vib)
;    end

    51002 : begin ; ClO main, Janssen-Rosenkranz
        vib=[842.4] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    53002 : begin ; ClO isotope, taken main
        vib=[842.4] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end


    60001 : begin ; OCS main, JANAF
        vib=[524,524,859,2064] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    62001 : begin ; OCS isotope, taken main
        vib=[524,524,859,2064] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    61001 : begin ; OCS isotope, taken main
        vib=[524,524,859,2064] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    62002 : begin ; OCS isotope, taken main
        vib=[524,524,859,2064] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    30004 : begin ; H2CO main, JANAF
        vib=[1163.5,1247.4,1500.6,1746.1,2766.4,2843.4] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    31002 : begin ; H2CO isotope, taken main
        vib=[1163.5,1247.4,1500.6,1746.1,2766.4,2843.4] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    32004 : begin ; H2CO isotope, taken main
        vib=[1163.5,1247.4,1500.6,1746.1,2766.4,2843.4] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    31003 : begin ; H2CO isotope, taken main
        vib=[1163.5,1247.4,1500.6,1746.1,2766.4,2843.4] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    32006 : begin ; H2CO isotope, taken main
        vib=[1163.5,1247.4,1500.6,1746.1,2766.4,2843.4] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    52006 : begin ; HOCl main, JANAF
        vib=[725,1239.4,3609.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    54005 : begin ; HOCl isotope, taken main
        vib=[725,1239.4,3609.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    27001 : begin ; HCN main, JANAF
        vib=[713.5,713.5,2096.3,3311.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    28002 : begin ; HCN isotope, taken main
        vib=[713.5,713.5,2096.3,3311.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    28003 : begin ; HCN isotope, taken main
        vib=[713.5,713.5,2096.3,3311.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    28004 : begin ; HCN isotope, taken main
        vib=[713.5,713.5,2096.3,3311.5] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    50007 : begin ; CH3Cl main, JANAF
        vib=[732,1017,1017,1355,1455,1455,2968,3054,3054] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    52009 : begin ; CH3Cl isotope, taken main
        vib=[732,1017,1017,1355,1455,1455,2968,3054,3054] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    34003 : begin ; PH3 main, JANAF
        vib=[992,1122,1122,2323,2328,2328] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    66001 : begin ; COF6 main, JANAF
        vib=[381,381,1074,1978] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    34002 : begin ; H2S main, JANAF
        vib=[1183,2615,2627] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    35001 : begin ; H2S isotope, taken main
        vib=[1183,2615,2627] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    67001 : begin ; OClO main, JANAF
        vib=[447.4,945.3,1109] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end
    69001 : begin ; OClO isotope, taken main
        vib=[447.4,945.3,1109] ; [cm^-1]
        jpl_part=calculate_partition_functions(jpl_part,vib)
    end

    else : begin
        jpl_part=jpl_part
        print,'No partition function correction defined.'
    end
endcase

end



PRO make_plot,part_fct,data_def,t_arr,ind,title,leg_str
;; plot routine

window,1
plot,t_arr[ind],part_fct[0,ind],title=title, $
  xtitle='Temperature [K]', ytitle='Partition Function [1]',xstyle=1
for i=1,data_def-1 do oplot,t_arr[ind],part_fct[i,ind],linestyle=i
        
;; make a simple legend
arr1=indgen(data_def) - indgen(data_def)
arr2=indgen(data_def)
legend,leg_str,$
  arr1,arr2,arr1,$
  !x.crange(1)-(!x.crange(1)-!x.crange(0))/1.2, $
  !y.crange(1)-!y.crange(1)/7, (abs(!y.crange(1)) + abs(!y.crange(0)))/15
end




PRO gis, species_arr, str, i, j, k, l, array
;; find the first free location in species_arr and store input there
ind=where(species_arr[*].name eq '')
ind=ind[0]


species_arr[ind].name        = str
species_arr[ind].arts_tag    = i
species_arr[ind].hit_mol     = j
species_arr[ind].hit_tag     = k
species_arr[ind].dgf         = l

if n_elements(array) ne 0 then $
  species_arr[ind].jpl_tag[0:n_elements(array)-1]  = array
end


PRO initialize_array,species_arr
;; some general data, refering jpl and hitran tags with 
;; each other, extracted from the arts program

;; structure of the data
species_data = { species, $
                 name    : '',        $ ;  name of species within ARTS
                 arts_tag: 0,         $ ;  arts isotope convention
                 hit_mol : 0,         $ ;  hitran molecule number code, 
                                        ;      e.g. major isotope of H2O:161
                                        ;      for molecules not present in 
                                        ;      hitran -1
                 hit_tag : 0,         $ ;  hitran tag number
                                        ;      e.g. major isotope of H2O:11
                                        ;           -1: not present 
                 dgf     : 0,         $ ;  degrees of freedom
                 jpl_tag : lonarr(10) } ;  jpl tag numbers of this hitran isotope

;; replicate the structure, currently there are not more than 200
;; species/isotopes. will be reduced to actual size later.
species_arr  = replicate({species},200)


;; put the data into structure use gis procedure, not too elegant, I
;; must admitt. data is extracted from arts program
gis, species_arr, 'H2O',   161,   161,    11, 3,  long([18003, 18005])
gis, species_arr, 'H2O',   181,   181,    12, 3,  long([20003])
gis, species_arr, 'H2O',   171,   171,    13, 3,  long([19003])
gis, species_arr, 'H2O',   162,   162,    14, 3,  long([19002])
gis, species_arr, 'H2O',   182,    -1,    -1, 3,  long([21001])
gis, species_arr, 'H2O',   262,    -1,    -1, 3,  long([20001])
gis, species_arr, 'CO2',   626,   626,    21, 2
gis, species_arr, 'CO2',   636,   636,    22, 2
gis, species_arr, 'CO2',   628,   628,    23, 2,  long([46013])
gis, species_arr, 'CO2',   627,   627,    24, 2,  long([45012])
gis, species_arr, 'CO2',   638,   638,    25, 2
gis, species_arr, 'CO2',   637,   637,    26, 2
gis, species_arr, 'CO2',   828,   828,    27, 2
gis, species_arr, 'CO2',   728,   728,    28, 2
gis, species_arr, 'O3',    666,   666,    31, 3,  long([48004, 48005, 48006, 48007, 48008])
gis, species_arr, 'O3',    668,   668,    32, 3,  long([50004, 50006])
gis, species_arr, 'O3',    686,   686,    33, 3,  long([50003, 50005])
gis, species_arr, 'O3',    667,   667,    34, 3,  long([49002])
gis, species_arr, 'O3',    676,   676,    35, 3,  long([49001])
gis, species_arr, 'N2O',   446,   446,    41, 2,  long([44004, 44009, 44012])
gis, species_arr, 'N2O',   456,   456,    42, 2,  long([45007])
gis, species_arr, 'N2O',   546,   546,    43, 2,  long([45008])
gis, species_arr, 'N2O',   448,   448,    44, 2,  long([46007])
gis, species_arr, 'N2O',   447,   447,    45, 2
gis, species_arr, 'CO',     26,   26,     51, 2,  long([28001])
gis, species_arr, 'CO',     36,   36,     52, 2,  long([29001])
gis, species_arr, 'CO',     28,   28,     53, 2,  long([30001])
gis, species_arr, 'CO',     27,   27,     54, 2,  long([29006])
gis, species_arr, 'CO',     38,   38,     55, 2
gis, species_arr, 'CO',     37,   37,     56, 2
gis, species_arr, 'CH4',   211,   211,    61, 3
gis, species_arr, 'CH4',   311,   311,    62, 3
gis, species_arr, 'CH4',   212,   212,    63, 3,  long([17003])
gis, species_arr, 'O2',     66,   66,     71, 2,  long([32001, 32002])
gis, species_arr, 'O2',     68,   68,     72, 2,  long([34001])
gis, species_arr, 'O2',     67,   67,     73, 2,  long([33002])
gis, species_arr, 'NO',     46,   46,     81, 2,  long([30008])
gis, species_arr, 'NO',     56,   56,     82, 2
gis, species_arr, 'NO',     48,   48,     83, 2
gis, species_arr, 'SO2',   626,   626,    91, 3,  long([64002, 64005])
gis, species_arr, 'SO2',   646,   646,    92, 3,  long([66002])
gis, species_arr, 'SO2',   636,    -1,    -1, 3,  long([65001])
gis, species_arr, 'SO2',   628,    -1,    -1, 3,  long([66004])
gis, species_arr, 'NO2',   646,   646,   101, 3,  long([46006])
gis, species_arr, 'NH3',  4111,   4111,  111, 3,  long([17002, 17004])
gis, species_arr, 'NH3',  5111,   5111,  112, 3,  long([18002])
gis, species_arr, 'NH3',  4112,     -1,   -1, 3,  long([18004])
gis, species_arr, 'HNO3',  146,   146,   121, 3,  long([63001, 63002, 63003, 63004, 63005, 63006])
gis, species_arr, 'OH',     61,   61,    131, 2,  long([17001])
gis, species_arr, 'OH',     81,   81,    132, 2,  long([19001])
gis, species_arr, 'OH',     62,   62,    133, 2,  long([18001])
gis, species_arr, 'HF',     19,   19,    141, 2,  long([20002])
gis, species_arr, 'HF',     29,   -1,     -1, 2,  long([21002])
gis, species_arr, 'HCl',    15,   15,    151, 2,  long([36001])
gis, species_arr, 'HCl',    17,   17,    152, 2,  long([38001])
gis, species_arr, 'HCl',    25,   -1,     -1, 2,  long([37001])
gis, species_arr, 'HCl',    27,   -1,     -1, 2,  long([39004])
gis, species_arr, 'HBr',    19,   19,    161, 2,  long([80001])
gis, species_arr, 'HBr',    11,   11,    162, 2,  long([82001])
gis, species_arr, 'HI',     17,   17,    171, 2
gis, species_arr, 'ClO',    56,   56,    181, 2,  long([51002, 51003])
gis, species_arr, 'ClO',    76,   76,    182, 2,  long([53002, 53006])
gis, species_arr, 'OCS',   622,   622,   191, 2,  long([60001])
gis, species_arr, 'OCS',   624,   624,   192, 2,  long([62001])
gis, species_arr, 'OCS',   632,   632,   193, 2,  long([61001])
gis, species_arr, 'OCS',   822,   822,   194, 2,  long([62002])
gis, species_arr, 'H2CO', 1126,   126,   201, 3,  long([30004])
gis, species_arr, 'H2CO', 1136,   136,   202, 3,  long([31002])
gis, species_arr, 'H2CO', 1128,   128,   203, 3,  long([32004])
gis, species_arr, 'H2CO', 1226,    -1,    -1, 3,  long([31003])
gis, species_arr, 'H2CO', 2226,    -1,    -1, 3,  long([32006])
gis, species_arr, 'HOCl',  165,   165,   211, 3,  long([52006])
gis, species_arr, 'HOCl',  167,   167,   212, 3,  long([54005])
gis, species_arr, 'N2',     44,   44,    221, 2
gis, species_arr, 'HCN',   124,   124,   231, 3,  long([27001, 27003])
gis, species_arr, 'HCN',   134,   134,   232, 3,  long([28002])
gis, species_arr, 'HCN',   125,   125,   233, 3,  long([28003])
gis, species_arr, 'HCN',   224,    -1,    -1, 3,  long([28004])
gis, species_arr, 'CH3Cl', 215,   215,   241, 3,  long([50007])
gis, species_arr, 'CH3Cl', 217,   217,   242, 3,  long([52009])
gis, species_arr, 'H2O2', 1661,   1661,  251, 3,  long([34004])
gis, species_arr, 'C2H2', 1221,   1221,  261, 3
gis, species_arr, 'C2H2', 1231,   1231,  262, 3
gis, species_arr, 'C2H6', 1221,   1221,  271, 3
gis, species_arr, 'PH3',  1111,   1111,  281, 3,  long([34003])
gis, species_arr, 'COF2',  269,   269,   291, 3,  long([66001])
gis, species_arr, 'SF6',    29,   29,    301, 3
gis, species_arr, 'H2S',   121,   121,   311, 3,  long([34002])
gis, species_arr, 'H2S',   141,   141,   312, 3
gis, species_arr, 'H2S',   131,   131,   313, 3
gis, species_arr, 'H2S',   122,    -1,    -1, 3,  long([35001])
gis, species_arr, 'HCOOH', 1261,  126,   321, 3,  long([46005])
gis, species_arr, 'HCOOH', 1361,   -1,    -1, 3,  long([47002])
gis, species_arr, 'HCOOH', 2261,   -1,    -1, 3,  long([47003])
gis, species_arr, 'HCOOH', 1262,   -1,    -1, 3,  long([47004])
gis, species_arr, 'HO2',   166,   166,   331, 3,  long([33001])
gis, species_arr, 'O',        6,   6,    341, 0,  long([16001])
gis, species_arr, 'ClONO2',5646, 5646,   351, 3,  long([97002])
gis, species_arr, 'ClONO2',7646, 7646,   352, 3,  long([99001])
gis, species_arr, 'NO+',     46,   46,   361, 2,  long([30011])
gis, species_arr, 'OClO'   ,656,    -1,   -1, 3,  long([67001])
gis, species_arr, 'OClO',   676,    -1,   -1, 3,  long([69001])
gis, species_arr, 'BrO',     96,    -1,   -1, 2,  long([95001])
gis, species_arr, 'BrO',     16,    -1,   -1, 2,  long([97001])
gis, species_arr, 'H2SO4',  126,    -1,   -1, 3,  long([98001])
gis, species_arr, 'Cl2O2',  565,    -1,   -1, 3,  long([102001])
gis, species_arr, 'Cl2O2',  765,    -1,   -1, 3,  long([104001])

;; resize array
ind=where(species_arr[*].name ne '',count)
species_arr = species_arr[0:count-1]
end



FUNCTION strip,str
;; removes all the preceding, following, and multiple blanks inbetween
;; of the string str
return, strtrim(strcompress(str),2)
end



;; MAIN ROUTINE
;----------------------------------------------------------------------------------------
;----------------------------------------------------------------------------------------

PRO partition_function,jpl_tag=jpl_tag,hit_tag=hit_tag,jpl_all=jpl_all,$
                       make_arts_entry=make_arts_entry,$
                       no_vib_correction=no_vib_correction,no_windows=no_windows

;; this procedure compares jpl and hitran partition functions (Q). 
;;
;; Calculation Schemes:
;;
;; JPL:
;; Qs are calculated according to the scheme given in the jpl
;; catalogue, e.g., the Qs tabulated are interpolated by assuming
;; either :
;;         linear molecule      : Q proportional to temp T
;;         non linear molecules : Q proportional to T^1.5
;;         atoms                : unknown
;;
;; HITRAN:
;; Qs are calculated as specified in the TIPS program which comes
;; along with the hitran catalogue. TIPS calculates Q as a third
;; order polynomical in T, coefficients are given in hitran cat
;;
;;
;; INPUT:
;;    - jpl_tag          : string of jpl tag, e.g., jpl_tag='32001'
;;    - hit_tag          : string of hit tag, e.g., hit_tag='O2-71'
;;    - jpl_all          : does the calculation for all jpl molecules
;;    - make_arts_entry  : creates a file with the partition
;;                         coefficients, ... which can be copied and
;;                         might has to be slightly modified for arts
;;    - no_vib_correction: flag, use no correction of the partition
;;                         functions. Corrects for not considered
;;                         vibrational energies in the JPL partition
;;                         funcitons. Only some of the species are
;;                         corrected, the vibrational energies are
;;                         provided by Patrick Eriksson, Miriam von
;;                         König, and a book research.
;;    - no_windows       : only printed output
;;
;; OUTPUT:
;;    - plots with Qs of both catalogues, for hitran the State
;;      independant degeneracy factor is either included or excluded 
;;    - Difference between Qs of both catalogues
;;    - Difference of ratio of Qs in both catalogues. this is the
;;      important quantaty, since the intensity at a certain
;;      temperarture is calculated with the ratio of Q. the ratio is
;;      calculated as: Q(300K)/Q(T) for both catalogues
;;    - 3rd order polynomial fit to jpl data, fit is only performed
;;      for temperatures between 70 and 500K, one of the temperature
;;      ranges given in the TIPS program.
;;    - partition function file for arts
;;
;;
;; HISTORY:
;;    15.08.00 AvE: created
;;



if n_elements(jpl_tag) eq 0 and n_elements(hit_tag) eq 0 $
  and n_elements(jpl_all) eq 0 and n_elements(make_arts_entry) eq 0 then begin
    print,''
    print,'Error: No JPL or HITRAN tag number given.'
    print,'Examples for Program Call: '
    print,'   For first isotope of oxygen:'
    print,"        partition_function,jpl_tag='32001'"
    print,"        partition_function,hit_tag='O2-71'"
    print,'   For all JPL molecules defined in ARTS:'
    print,'        partition_function,/jpl_all'
    print,'   No vibraional correction of JPL partition functions,'
    print,'        only available for certain species.'
    print,'        partition_function,.....,/no_vib_correction'
    print,'   For arts entry file:'
    print,'        partition_function,/make_arts_entry'
    print,'   For only text output:'
    print,'        partition_function,...,/no_windows'
    print,''
    return
endif

if keyword_set(no_windows) then no_windows=1 else no_windows=0

;; initialize array that refers jpl tags to hitran tags. this is
;; extracted from the arts program
initialize_array,species_arr

;; generate temperature array 0 to 500 K, T^2, T^3, this is the
;; important area of temperatures
t_arr=indgen(501)*1.0
t_arr2=t_arr*t_arr
t_arr3=t_arr*t_arr*t_arr

;; what is the index to 300 K, since this is used later in the ratio
;; calculation of partition functions
ind300=where(t_arr eq 300)

;; restrict fit to the TIPS area of 70 to 500 K, this is the first
;; temperature range given in TIPS, the other ones (> 500 K) are in
;; general not relevant to atmospheric temperatures
ind=where(t_arr ge 70)

;; jpl partition functions are only given between 9 and 300 K
ind_le_300=where(t_arr le 300)

;; index to all elements between 150 and 300 K, e.g., the atmospheric
;; range of interest
ind_150_300=where((t_arr ge 150) and (t_arr le 300))

;; which species should be plotted ?

;; keyword jpl_all set, check all jpl species
if keyword_set(jpl_all) then begin
    loop_index = where(species_arr[*].jpl_tag[0] ne 0)
endif

;; keyword make_arts_enty set, all species
if keyword_set(make_arts_entry) then begin
    loop_index = where(species_arr[*].name ne '')

    ;; these molecules do not exist in hitran, but a polynominal is
    ;; fittable to the jpl results, therefore these jpl coefficients
    ;; are  used in the output file 
    ;; identified by name and arts isotope
    use_jpl=['HNO3-146','ClONO2-5646','ClONO2-7646']

    ;; for these molecules, the polynomial fit to the JPL partition
    ;; functions for the temperature range 70 to 500 K does not work
    ;; that well, therefore the 0 to 300 K range is used instead
    use_jpl_300=['47002','47003','47004']

    ;; make stucture to hold arts entries
    arts_arr = { arts, $
                 name     : '',         $ ;  name of species within ARTS
                 sr_coeff : intarr(10), $ ;  source of coefficients, not more than 10 isotopes
                 quality  : strarr(10), $ ;  quality of fit
                 iso      : strarr(10)}   ;  isotope records

    ;; replicate the structure, currently there are not more than 200
    ;; species/isotopes. will be reduced to actual size later.
    arts_arr  = replicate({arts},200)

endif


;; keyword jpl_tag set, plot just one species
if keyword_set(jpl_tag) then begin
    loop_index = where(species_arr[*].jpl_tag[*] eq jpl_tag,count)
    if count eq 0 then begin
        print,'Error: JPL tag not found: '+jpl_tag
        stop
    endif else begin
        if count gt 1 then begin
            print,'Error: More than one JPL tag found.'
            stop
        endif
    endelse
    loop_index=floor(loop_index/n_elements(species_arr[0].jpl_tag)*1.0)
endif

;; keyword hit_tag set, plot just one species
if keyword_set(hit_tag) then begin
    loop_index = where(species_arr[*].name+'-'+$
                       strip(string(species_arr[*].hit_tag,format='(I3)')) eq hit_tag,count)
    if count eq 0 then begin
        print,'Error: HITRAN tag not found: '+hit_tag
        stop
    endif else begin
        if count gt 1 then begin
            print,'Error: More than one HITRAN tag found.'
            stop
        endif
    endelse
endif

;; open a standard file for arts entries
if keyword_set(make_arts_entry) then begin
    print,'Generating arts_partition_fct.txt file with arts input'
    last_mol='dummy'
    ind1=-1
    openw,10,'arts_partition_fct.txt'
endif


;; now cycle through all species identified
for ii=0,n_elements(loop_index)-1 do begin

    ;; initializations
    jpl=0                       ; tag not present
    hitran=0                    ; tag not present
    hit_ok = 0                  ; partition fct are ok
    jpl_ok = 0                  ; partition fct are ok

    ;; initialize arrays that hold the polynomical coefficients
    hit_coeff=dblarr(4)
    jpl_coeff=dblarr(4)
    jpl_coeff_full=dblarr(4)

    print,'-----------------------------------------------------------------'
        
    ;; find the corresponding tag of the other catalogue
    if species_arr[loop_index[ii]].jpl_tag[0] ne 0 then begin
        jpl = 1
        if keyword_set(jpl_all) or keyword_set(hit_tag) or $
          keyword_set(make_arts_entry) then $
          jpl_tag=species_arr[loop_index[ii]].jpl_tag[0]
    endif

    if species_arr[loop_index[ii]].hit_tag ne -1 then begin
        hitran = 1
    endif


    ;; output molecule
    print,''
    print,'Molecule: ',species_arr[loop_index[ii]].name
    if jpl then print,'         JPL tag   : ',strip(jpl_tag)
    if hitran then $
      print,'         HITRAN tag: ',strip(species_arr[loop_index[ii]].name)+'-'+$
      strip(species_arr[loop_index[ii]].hit_tag)
    print,'         ARTS tag : ',strip(species_arr[loop_index[ii]].name)+'-'+$
      strip(species_arr[loop_index[ii]].arts_tag)
    print,''


    if jpl then begin
        ;; search for the jpl_tag in the catdir.cat file
        str=''
        openr,1,'catdir.cat'
        ;; tag number is given at position 0 to 6
        while (strip(strmid(str,0,6)) ne strip(long(jpl_tag))) and not eof(1) do begin
            readf,1,str
        endwhile
        close,1
        
        ;; extract the partition function values out of str, first define all
        ;; variables
        tagnr=0 & jplname='' & nline=0 & qlog300=0.0 & qlog225=0.0 & qlog150=0.0 & qlog75=0.0
        qlog37=0.0 & qlog18=0.0 & qlog9=0.0 & version=0
        ;; then read them
        reads,str,tagnr,jplname,nline,qlog300,qlog225,qlog150,qlog75,qlog37,qlog18,$
          qlog9,version,format='(I6,A15,I5,7(F7.4),I2)'
            
        ;; store the partition functions in array, first one is temp, second
        ;; the partition function
        par_arr_jpl=[ [  9.375  , 10.0^qlog9  ],$
                      [  18.75  , 10.0^qlog18 ],$
                      [  37.5   , 10.0^qlog37 ],$
                      [ 75.0    , 10.0^qlog75 ],$
                      [150.0    , 10.0^qlog150],$
                      [225.0    , 10.0^qlog225],$
                      [300.0    , 10.0^qlog300  ] ]

        ;; do we have to correct the jpl part fct?
        if keyword_set(no_vib_correction) then $
          correct_jpl = 0 $ 
        else correct_jpl = species_arr[loop_index[ii]].jpl_tag[0]
        correct_jpl_partition_functions,par_arr_jpl,$
          correct_jpl
    endif
    
    
    if hitran then begin
        ;; now find the hitran partition function coefficients in file
        ;; tips_coeff.txt, which list the block data of fortran code of the
        ;; tips program in slightly modified way, plus the State
        ;; independent degeneracy factors 

        ;; this is the format of the species name in that file
        spe=species_arr[loop_index[ii]].name+' -- '+$
          strip(string(species_arr[loop_index[ii]].hit_mol,format='(I4)'))
        str=''
        openr,1,'tips_coeff.txt'
        while (strip(str) ne spe) and (not eof(1)) do begin
            readf,1,str
        endwhile

        if eof(1) then begin
            print,'Error: Species not found in file tips_coeff.txt.'
            close,1
            stop
        endif else begin
            ;; now read the tips coefficients and the state
            ;; independant degeneracy factor sidf
            readf,1,c0,c1,c2,c3,sidf
            hit_coeff[0]=c0 & hit_coeff[1]=c1
            hit_coeff[2]=c2 & hit_coeff[3]=c3
            close,1
        endelse
    endif


    if jpl then begin
        ;; check the degrees of freedom and set exponent for interpolation
        ;; routine, exponent comes form jpl catalogue, unfortunately,
        ;; we don't know how to handle atoms so far
        case species_arr[loop_index[ii]].dgf of
            2: begin
                nn=1.0
                jpl_ok = 1
            end
            3: begin
                nn=1.5
                jpl_ok = 1
            end
            else: begin
                print,'Warning: Degrees of freedom interpolation scheme not specified.'
                nn=10.0 ; set to dummy
                jpl_ok = 0 ; ; partition functions of jpl are not OK
            end
        endcase

        ;; interpolation scheme not specified in JPL, follow hitran convention
        if not jpl_ok then begin
            jpl_coeff[0] = -1.
            jpl_coeff[1] =  0.
            jpl_coeff[2] =  0.
            jpl_coeff[3] =  0.

            jpl_coeff_full[0] = -1.
            jpl_coeff_full[1] =  0.
            jpl_coeff_full[2] =  0.
            jpl_coeff_full[3] =  0.
        endif else begin ;; partition functions defined

            ;; calculate jpl partition function at all temperatures up to 500 K
            ;; using interpolation scheme given in jpl catalogue
            Q_temp_jpl=interpol(par_arr_jpl[1,*],par_arr_jpl[0,*]^nn,t_arr^nn)

            ;; now fit a polynomical to jpl partition function, use only
            ;; temperature range 70 to 500 K
            jpl_coeff = poly_fit(t_arr[ind],Q_temp_jpl[ind],3)

            ;; full range fit of available JPL data (0-300 K)
            jpl_coeff_full = poly_fit(t_arr[ind_le_300],Q_temp_jpl[ind_le_300],3)

            ;; calculate the jpl partition functions with polynomical
            ;; coeff with fit range up to 500 K 
            Q_temp_jpl_poly = jpl_coeff[0] + $
              jpl_coeff[1] * t_arr + $
              jpl_coeff[2] * t_arr2 + $
              jpl_coeff[3] * t_arr3
        
            ;; calculate the jpl partition functions with polynomical
            ;; coeff with fit range up to 300 K 
            Q_temp_jpl_poly_full = jpl_coeff_full[0] + $
              jpl_coeff_full[1] * t_arr + $
              jpl_coeff_full[2] * t_arr2 + $
              jpl_coeff_full[3] * t_arr3
        endelse
        
    endif
    
    if hitran then begin
        ;; calculate the hitran partition function at all temp up
        ;; to 500 K, with and without the state independant
        ;; degeneracy factor considered
        Q_temp_hit=(hit_coeff[0] + $
                    hit_coeff[1] * t_arr + $
                    hit_coeff[2] * t_arr2 + $
                    hit_coeff[3] * t_arr3)/sidf
        Q_temp_hit_in_sidf = Q_temp_hit * sidf

        ;; polynomical fit not applicable, as specified in TIPS
        ;; program
        if Q_temp_hit[ind300[0]] lt 0 then begin
            print,''
            print,'Warning: TIPS-Polynomial not applicable to this species.'
            print,''
            hit_ok = 0
        endif else hit_ok = 1
    endif
    

    if not keyword_set(make_arts_entry) then begin
    
        ;; store all partition functions in this array, we generally
        ;; do not plot more than 3 curves
        part_fct = dblarr(3,n_elements(t_arr))
    
        ;; make a plot of the data
        if jpl and jpl_ok and hitran and hit_ok then begin
            ;; all data are defined and partition fct are OK
            part_fct[0,*] = Q_temp_jpl
            part_fct[1,*] = Q_temp_hit
            part_fct[2,*] = Q_temp_hit_in_sidf


            data_def=3

            title='HITRAN/JPL Partition Functions'

            leg_str=['JPL Name: '+strip(jplname)+$
                     ',    JPL Tag: '+strip(jpl_tag),$
                     'HITRAN Name: '+species_arr[loop_index[ii]].name+$
                     ',    HITRAN Tag: '+$
                     strip(string(species_arr[loop_index[ii]].hit_tag))+$
                     ', SIDF : '+string(sidf,format='(I3)'),$
                     'HITRAN Name: '+species_arr[loop_index[ii]].name+$
                 ',    HITRAN Tag: '+$
                     strip(string(species_arr[loop_index[ii]].hit_tag))+', SIDF in' ]
        endif 
        
        if jpl and jpl_ok and not (hitran and hit_ok) then begin
            ;; jpl is OK, not hitran
            part_fct[0,*] = Q_temp_jpl
            
            data_def=1

            title='JPL Partition Function'
            leg_str=['JPL Name: '+strip(jplname)+',    JPL Tag: '+strip(jpl_tag)]
        endif 
        
        if not jpl and hitran and hit_ok then begin
            ;; hitran is OK, jpl not
            part_fct[0,*] = Q_temp_hit

            data_def=1

            title='HITRAN Partition Function'
            leg_str=['HITRAN Species: '+$
                     strip(string(species_arr[loop_index[ii]].name)) +'-'+$
                     strip(string(species_arr[loop_index[ii]].hit_tag))]
        endif 

        ;; call plot routine
        if not no_windows then begin

            make_plot,part_fct,data_def,t_arr,ind,title,leg_str
    
        
            if hitran and hit_ok and jpl and jpl_ok then begin
                !P.multi=[0,1,2,0]
            
                ;; make a plot of the differences
                window,2
                plot,t_arr[ind],(Q_temp_jpl[ind] - Q_temp_hit[ind])/Q_temp_jpl[ind] $
                  * 100.0,$
                  title='Difference of Partition Functions, '+'HITRAN Tag: '+$
                  species_arr[loop_index[ii]].name+$
                  '-'+strip(string(species_arr[loop_index[ii]].hit_tag))+$
                  ', JPL Tag: '+strip(jpl_tag),$
                  xtitle='Temperature [K]', ytitle='Difference [%]',$
                  yrange=[-10,10],xstyle=1
                oplot,t_arr[ind],(Q_temp_jpl[ind] - Q_temp_jpl_poly[ind])/$
                  Q_temp_jpl[ind] * 100.0,linestyle=1
                oplot,t_arr[ind],(Q_temp_jpl_poly[ind] - Q_temp_hit[ind])/$
                  Q_temp_jpl_poly[ind] * 100.0,linestyle=2
                oplot,t_arr[ind],(Q_temp_jpl[ind] - Q_temp_hit_in_sidf[ind])/$
                  Q_temp_jpl[ind] * 100.0,linestyle=3
                oplot,t_arr[ind],(Q_temp_jpl[ind] - Q_temp_jpl_poly_full[ind])/$
                  Q_temp_jpl[ind] * 100.0,linestyle=4
            
                ;; make a simple legend
                legend,['(JPL - HITRAN)/JPL','(JPL - JPL!dpoly!n)/JPL',$
                        '(JPL!dpoly!n-HITRAN)/JPL!dpoly!n','(JPL - HITRAN!dSIDF in!n)/JPL',$
                        '(JPL - JPL!dpoly full fit!n)/JPL'],$
                  [0,0,0,0,0],[0,1,2,3,4],[0,0,0,0,0],$
                  !x.crange(1)-(!x.crange(1)-!x.crange(0))/1.2, $
                  !y.crange(1)-!y.crange(1)/7, (abs(!y.crange(1)) + abs(!y.crange(0)))/15
                
                ;; compare ratio of Q at 300 K to Q at other values
                plot,t_arr,( (Q_temp_jpl[ind300[0]]/Q_temp_jpl) - $
                             (Q_temp_hit[ind300[0]]/Q_temp_hit)) /$
                  (Q_temp_jpl[ind300[0]]/Q_temp_jpl) * 100.0 ,$
                  title='Difference in Ratios: (Q!d1!n(300 K)/Q!d1!n(T)'+$
                  '-Q!d2!n(300 K)/Q!d2!n(T))/Q!d1!n(300 K)/Q!d1!n(T)',$
                  xtitle='Temperature [K]', ytitle='Difference [%]',yrange=[-10,10],xstyle=1
                oplot,t_arr,( (Q_temp_jpl[ind300[0]]/Q_temp_jpl) - $
                              (Q_temp_jpl_poly[ind300[0]]/Q_temp_jpl_poly)) /$
                  (Q_temp_jpl[ind300[0]]/Q_temp_jpl) * 100.0,linestyle=1
                oplot,t_arr,( (Q_temp_jpl[ind300[0]]/Q_temp_jpl) - $
                              (Q_temp_jpl_poly_full[ind300[0]]/Q_temp_jpl_poly_full)) /$
                  (Q_temp_jpl[ind300[0]]/Q_temp_jpl) * 100.0,linestyle=2
                
                ;; make a simple legend
                legend,['1: JPL, 2: HITRAN','1: JPL, 2: JPL poly 70-500K',$
                        '1: JPL, 2: JPL poly 0-300K'],$
                  [0,0,0],[0,1,2],[0,0,0],$
                  !x.crange(1)-(!x.crange(1)-!x.crange(0))/1.2, $
                  !y.crange(1)-!y.crange(1)/7, (abs(!y.crange(1)) + abs(!y.crange(0)))/15
                
                
                !P.multi=0
            endif else begin
                if not hitran and jpl and jpl_ok then begin
                    
                    ;; make a plot of the differences between the poly fit
                    ;; to jpl and the jpl recommended part fct calculation
                    window,2
                    ;; compare ratio of Q at 300 K to Q at other values
                    plot,t_arr,( (Q_temp_jpl[ind300[0]]/Q_temp_jpl) - $
                                 (Q_temp_jpl_poly[ind300[0]]/Q_temp_jpl_poly)) /$
                      (Q_temp_jpl[ind300[0]]/Q_temp_jpl) * 100.0 ,$
                      title='Difference in Ratios: (Q!d1!n(300 K)/Q!d1!n(T)'+$
                      '-Q!d2!n(300 K)/Q!d2!n(T))/Q!d1!n(300 K)/Q!d1!n(T)',$
                      xtitle='Temperature [K]', ytitle='Difference [%]',yrange=[-10,10],xstyle=1
                    oplot,t_arr,( (Q_temp_jpl[ind300[0]]/Q_temp_jpl) - $
                                  (Q_temp_jpl_poly_full[ind300[0]]/Q_temp_jpl_poly_full)) /$
                      (Q_temp_jpl[ind300[0]]/Q_temp_jpl) * 100.0,linestyle=1
                    
                    ;; make a simple legend
                    legend,['1: JPL, 2: JPL poly 70-500K',$
                            '1: JPL, 2: JPL poly 0-300K'],$
                      [0,0],[0,1],[0,0],$
                      !x.crange(1)-(!x.crange(1)-!x.crange(0))/1.2, $
                      !y.crange(1)-!y.crange(1)/7, (abs(!y.crange(1)) + abs(!y.crange(0)))/15
                    
                    
                endif else            window, 2 ; just empty to avoid confusion
            endelse
        endif

        
    
        ;; output some stuff on the screen
        print,''
        print,''
        print,'Found HITRAN Polynomial Coefficients of Partition Function'
        print,'      and fitted JPL values:'
        print,'               HITRAN                 JPL (70-500K)       JPL (0-300K)'
        for i=0,3 do print,'c'+string(i,format='(I1)'),':  ',+$
          string(hit_coeff[i],format='(E20.4)'),'   ',$
          string(jpl_coeff[i],format='(E20.4)')+$
          string(jpl_coeff_full[i],format='(E20.4)')
        print,''
        print,''

        ;; this output is only useful if polynomial is applicable
        if jpl and jpl_ok then begin
            ;; give the maximum difference found between the atmospheric
            ;; temperature of interest (150-300 K) found between the original
            ;; JPL partition functions and the polynomial fit
            a=Q_temp_jpl[ind300[0]]/Q_temp_jpl[ind_150_300]
            b=Q_temp_jpl_poly[ind300[0]]/Q_temp_jpl_poly[ind_150_300]
            c=Q_temp_jpl_poly_full[ind300[0]]/Q_temp_jpl_poly_full[ind_150_300]

            diff1 =  max( abs( (a-b)/a * 100.0 ) )
            diff2 =  max( abs( (a-c)/a * 100.0 ) )


            print,''
            print,''
            print,'Maximum difference in ratio of Q(300K)/Q(T) between original JPL '
            print,'and JPL polynominal fit for temperature range 150 to 300 K:'
            print,'            '+string(diff1,format='(F5.2)') +' % for 70 to 500 K fit'
            print,'            '+string(diff2,format='(F5.2)') +' % for  0 to 300 K fit'
            print,''
            if hitran and hit_ok then begin

                d=Q_temp_hit[ind300[0]]/Q_temp_hit[ind_150_300]
                diff3 =  max( abs( (a-d)/a * 100.0 ) )
                print,'Maximum difference in ratio of Q(300K)/Q(T) between original JPL '
                print,'and HITRAN for temperature range 150 to 300 K:'
                print,'            '+string(diff3,format='(F5.2)') +' %'
            endif
            print,''
        endif

        print,'-----------------------------------------------------------------'
    
    
        if keyword_set(jpl_all) then begin
            hak,/mesg
        endif

    endif else begin

        ;; make the entry for the arts program
        ;;------------------------------------


        ;; check wether we still have the same molecule
        if last_mol eq species_arr[loop_index[ii]].name then begin
            same=1
            ind2=ind2+1
        endif else begin
            same=0
            ind1=ind1+1
            ind2=0
        endelse
        last_mol = species_arr[loop_index[ii]].name

        ;; put data into arts_arr
        if not same then begin
            arts_arr[ind1].name=species_arr[loop_index[ii]].name
        endif 
    
        ;; check whether this species is not applicable in hitran and
        ;; jpl should be used
        str=strip(species_arr[loop_index[ii]].name)+'-'+$
          strip(string(species_arr[loop_index[ii]].arts_tag))
        dummy=where(str eq use_jpl,count)
        if count ne 0 then hitran=0

        ;; check whether the 0 to 300 K fit range should be used
        ;; instead of the 0 to 500 K range
        str=strip(species_arr[loop_index[ii]].jpl_tag[0])
        dummy=where(str eq use_jpl_300,count)
        if count ne 0 then use_300=1 else use_300=0

        if hitran and hit_ok and jpl and jpl_ok then begin
            sr=1
            ;; quality
            a=Q_temp_jpl[ind300[0]]/Q_temp_jpl[ind_150_300]
            b=Q_temp_hit[ind300[0]]/Q_temp_hit[ind_150_300]
            c=(a-b)/a * 100.0
            qual = string(max( abs(c) ),format='(F5.2)')
            ;; entry string
            arts_ent = '"'+$
              strip(string(species_arr[loop_index[ii]].arts_tag))+'",	Qcoeff(	'+ $
              strip(string(hit_coeff[0],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[1],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[2],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[3],format='(E11.4)'))+') );'
        endif 

        if hitran and hit_ok and not jpl then begin
            sr=1
            ;; quality not defined
            qual = '----'
            ;; entry string
            arts_ent = '"'+$
              strip(string(species_arr[loop_index[ii]].arts_tag))+'",	Qcoeff(	'+ $
              strip(string(hit_coeff[0],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[1],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[2],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[3],format='(E11.4)'))+') );'
        endif

        if jpl and jpl_ok and not hitran then begin
            sr=2

            ;; quality
            if use_300 then begin
                a=Q_temp_jpl[ind300[0]]/Q_temp_jpl[ind_150_300]
                b=Q_temp_jpl_poly_full[ind300[0]]/Q_temp_jpl_poly_full[ind_150_300]
                c=(a-b)/a * 100.0
                entry=jpl_coeff_full
            endif else begin
                a=Q_temp_jpl[ind300[0]]/Q_temp_jpl[ind_150_300]
                b=Q_temp_jpl_poly[ind300[0]]/Q_temp_jpl_poly[ind_150_300]
                c=(a-b)/a * 100.0
                entry=jpl_coeff
            endelse

            qual = string(max( abs(c) ),format='(F5.2)')

            ;; entry string
            arts_ent = '"'+$
              strip(string(species_arr[loop_index[ii]].arts_tag))+'",	Qcoeff(	'+ $
              strip(string(entry[0],format='(E11.4)'))+'	,'+$
              strip(string(entry[1],format='(E11.4)'))+'	,'+$
              strip(string(entry[2],format='(E11.4)'))+'	,'+$
              strip(string(entry[3],format='(E11.4)'))+') );'
        endif

        if not hit_ok and not jpl_ok then begin
            sr=1
            ;; quality not defined
            qual = '----'
            ;; entry string
            arts_ent = '"'+$
              strip(string(species_arr[loop_index[ii]].arts_tag))+'",	Qcoeff(	'+ $
              strip(string(hit_coeff[0],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[1],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[2],format='(E11.4)'))+'	,'+$
              strip(string(hit_coeff[3],format='(E11.4)'))+') );'
        endif



        arts_arr[ind1].sr_coeff[ind2] = sr
        arts_arr[ind1].quality[ind2] = qual
        arts_arr[ind1].iso[ind2] = arts_ent

    endelse

endfor


if keyword_set(make_arts_entry) then begin

    ;; write arts_arr info into file

    ;; resize array
    ind=where(arts_arr[*].name ne '',count)
    arts_arr = arts_arr[0:count-1]

    for i=0,n_elements(arts_arr)-1 do begin
        
        ;; how many entries are there
        ind=where(arts_arr[i].sr_coeff[*] ne 0,count)
        iso=count

        ;; print the header 
        printf,10,'  // '+arts_arr[i].name
        printf,10,'  // Coeff: '+string(arts_arr[i].sr_coeff[0:iso-1],$
                                        format='('+strip(string(iso))+'(I7))')
        
        printf,10,'  // Quality: '+string(arts_arr[i].quality[0:iso-1],$
                                        format='('+strip(string(iso))+'(A7))')

        printf,10,'  spec(it_species, it_isotope, "'+strip(arts_arr[i].name)+'");'
        printf,10,'  //			Name		c0		c1		c2		c3'
        printf,10,'  //			|		|		|		|		|'
        for j=0,count-1 do printf,10,'  iso(it_isotope,	'+arts_arr[i].iso[j]
        printf,10,'' & printf,10,'' & printf,10,''

    endfor

    close,10

endif

end 
