;;#################################################################################
;; 
FUNCTION WVSatPressureLiquidWater, TK, punit
;; 
;; The recommended formulas for the calculation of equilibrium water vapor pressure 
;; in the pure phase (no other gases present) over a plane surface of liquid water 
;; or ice are those of Sonntag, 1994. (References below.) The formulas are quoted 
;; also in the papers of Leiterer et al., 1997 and Helten et al., 1999. I have 
;; verified that they are exactly the same in all three papers. (Sonntag gives 
;; the formulas for partial pressure in Pa and hPa, Leiterer et al. for hPa, 
;; Helten et al. for Pa.)
;; 
;; The equation and coefficients below are for the equilibrium water vapor pressure 
;; in Pa/hPa and the temperature in K.
;; 
;; The temperature of 0C corresponds to 273.15K. (Not 273.16K, as stated in the 
;; Leiterer paper.)
;;
;; References:
;; Sonntag, D., Advancements in the field of hygrometry, 
;; Meteorologische Zeitschrift, 3, 51-66, 1994.
;;
;; Helten, M. et al, In-flight comparison of MOZAIC and POLINAT water
;; vapor measurements, JGR, 104, 26.087-26.096, 1999.
;;
;; Leiterer U. et al, Improvements in Radiosonde Humidity Profiles 
;; Using RS80/RS90 Radiosondes of Vaisala, 
;; Beitr. Phys. Atmosph., 70(4), 319-336, 1997.

;; ---------------------------------------------------------------------
;; INPUT  :  TK      DOUBLE  temperature [K]
;;           punit   STRING  pressure unit 'Pa' or 'hPa'     
;;
;; OUTPUT :  equilibrium water vapor pressure in the pure phase 
;;           (no other gases present) over a plane surface of 
;;           liquid water                                          [Pa]
;;
;; HISTORY:  alpha version,  2003-04-02, TKS
;; ---------------------------------------------------------------------

a  = DOUBLE(-6096.9385)
if (STRUPCASE(punit) EQ 'PA') then begin
    b  = DOUBLE(21.2409642)
endif else begin
    b  = DOUBLE(16.635794)
endelse
c  = DOUBLE(-2.711193e-2)
d  = DOUBLE(1.673952e-5)
e  = DOUBLE(2.433502)

ew = EXP( DOUBLE(a/TK)    + $
          b               + $ 
          DOUBLE(c*TK)    + $
          DOUBLE(d*TK*TK) + $
          DOUBLE(e*ALOG(TK)) )


return, eW ;; saturation pressure in Pa or hPa
end
;;
;; #################################################################################
;;
FUNCTION WVSatPressureIce, TK, punit
;;
;; The recommended formulas for the calculation of equilibrium water vapor pressure 
;; in the pure phase (no other gases present) over a plane surface of liquid water 
;; or ice are those of Sonntag, 1994. (References below.) The formulas are quoted 
;; also in the papers of Leiterer et al., 1997 and Helten et al., 1999. I have 
;; verified that they are exactly the same in all three papers. (Sonntag gives 
;; the formulas for partial pressure in Pa and hPa, Leiterer et al. for hPa, 
;; Helten et al. for Pa.)
;; 
;; The equation and coefficients below are for the equilibrium water vapor pressure 
;; in Pa/hPa and the temperature in K.
;; 
;; The temperature of 0C corresponds to 273.15K. (Not 273.16K, as stated in the 
;; Leiterer paper.)
;;
;; References:
;; Sonntag, D., Advancements in the field of hygrometry, 
;; Meteorologische Zeitschrift, 3, 51-66, 1994.
;;
;; Helten, M. et al, In-flight comparison of MOZAIC and POLINAT water
;; vapor measurements, JGR, 104, 26.087-26.096, 1999.
;;
;; Leiterer U. et al, Improvements in Radiosonde Humidity Profiles 
;; Using RS80/RS90 Radiosondes of Vaisala, 
;; Beitr. Phys. Atmosph., 70(4), 319-336, 1997.

;; ---------------------------------------------------------------------
;; INPUT  :  TK      DOUBLE  temperature [K]
;;           punit   STRING  pressure unit 'Pa' or 'hPa'     
;;
;; OUTPUT :  equilibrium water vapor pressure in the pure phase 
;;           (no other gases present) over ice                     [Pa]
;;
;; HISTORY:  alpha version,  2003-04-02, TKS
;; ---------------------------------------------------------------------

a  = DOUBLE(-6024.5282)
if (STRUPCASE(punit) EQ 'PA') then begin
    b  = DOUBLE(29.32707)
endif else begin
    b  = DOUBLE(24.7219)
endelse
c  = DOUBLE(1.0613868e-2)
d  = DOUBLE(-1.3198825e-5)
e  = DOUBLE(-0.49382577)

ei = EXP( DOUBLE(a/TK)    + $
          b               + $ 
          DOUBLE(c*TK)    + $
          DOUBLE(d*TK*TK) + $
          DOUBLE(e*ALOG(TK)) )


return, ei ;; saturation pressure in Pa or hPa
end
;;
;; #################################################################################
;;
FUNCTION AirCorrFunWater, P, TK
;;
;; ----------------------------------------------------------------
;; INPUT
;;       P    DOUBLE   total pressure  [Pa]
;;       TK   DOUBLE   temperature  [K]
;;
;; OUTPUT
;;       fw   DOUBLE   correction function for moist air over water
;; ----------------------------------------------------------------

P  = (P * 1.00000E-2)     ;; pressure in hPa 
TC = DOUBLE(TK - 273.15)  ;; temperature in Celsius
ew = WVSatPressureLiquidWater( DOUBLE(TK), 'hPa' )

a  = ew * 1.000E-4 / (273.000E0+TC)

b  = 38.000E0 +                                             $
     173.000E0*EXP(-TC/43.000E0)*(1.000E0-(ew/p)) +         $
     (6.390E0+4.280E0*EXP(-TC/107.000E0))*((p/ew)-1.000E0)

fw = 1.000E0 + (a * b)

return, fw
END
;;
;;#################################################################################
;;
FUNCTION AirCorrFunIce, P, TK
;;
;; --------------------------------------------------------------
;; INPUT
;;       P    DOUBLE   total pressure  [Pa]
;;       TK   DOUBLE   temperature  [K]
;;
;; OUTPUT
;;       fi   DOUBLE   correction function for moist air over ice
;; --------------------------------------------------------------

P  = (P * 1.00000E-2)     ;; pressure in hPa 
TC = DOUBLE(TK - 273.15)  ;; temperature in Celsius
ei = WVSatPressureLiquidIce( DOUBLE(TK), 'hPa' )

a  = ei * 1.000E-5 / (273.000E0+TC)

b  = (2100.000E0-65.000E0*TC)*(1.000E0-(ei/p))  +          $
     (109.000E0-0.350E0*TC+(TC*TC/338.000E0))*((p/ei)-1.000E0)

fi = 1.000E0 + (a * b)

return, fi
END
;;
;; #################################################################################
;;
FUNCTION WaterVaporSatPressure, TK, P, phase=phase, punit=punit, corr=corr
;;
;; ---------------------------------------------------------------------
;; INPUT:
;;        TK     DBLE    temperature [T]
;;        P      DBLE    total pressure [Pa]
;;        phase  STRING  water phase, possible phases are 'water' and 'ice'
;;        punit  STRING  pressure unit of the saturation pressure
;;                       possible units are 'Pa' and 'hPa' and 'mbar'
;;        corr   STRING  consider/not consider the correction term for
;;                       moist air instead of pure water vapor over waterr/ice
;;                       possible values are 'yes' and 'no'
;;
;; OUTPUT:
;;        es     DBLE    water vapor saturtation pressure [Pa, hPa, or mbar]
;;                       if the calculation had some error detected, 
;;                       the return value is es = -9.99999e9
;;
;; COMMENT:
;;        formula taken from reference:
;;        Sonntag, D., Advancements in the field of hygrometry, 
;;        Meteorologische Zeitschrift, 3, 51-66, 1994.
;;
;; HISTORY:  
;;        alpha version,  2003-04-02, TKS
;; ---------------------------------------------------------------------

;; default value (return value for erroneous calculation)
es = -9.99999e9

;; check the input variable 'phase'
if not KEYWORD_SET(phase) then phase='WATER'


;; check the input variable 'pressure unit'
if not KEYWORD_SET(punit) then phase='PA'
if (STRUPCASE(punit) EQ 'MBAR') then punit='HPA'


;; check the input variable for 'correction for moist air'
if not KEYWORD_SET(corr) then corr='no'


if ( (STRUPCASE(phase) NE 'WATER') AND (STRUPCASE(phase) NE 'ICE') ) then begin
    print,'WaterVaporSatPressure> !!! ERROR!'
    print,'WaterVaporSatPressure> !!! wrong water phase unit selected,'
    print,'WaterVaporSatPressure> !!! possible units are "water" and "ice"'
    return, es
endif 


if ( (STRUPCASE(punit) NE 'PA') AND (STRUPCASE(punit) NE 'HPA') ) then begin
    print,'WaterVaporSatPressure> !!! ERROR!'
    print,'WaterVaporSatPressure> !!! wrong pressure unit selected,'
    print,'WaterVaporSatPressure> !!! possible units are "Pa" and "hPa"'
    return, es
endif 


;; check the input variable temperature
if ( STRUPCASE(phase) EQ 'WATER' ) then begin
    if ((TK LT 173.15) OR (TK GT 373.15)) then begin
        print,'WaterVaporSatPressure> !!! ERROR!'
        print,'WaterVaporSatPressure> !!! temperature is out of the validity range'
        print,'WaterVaporSatPressure> !!! for the saturation over water: T=',TK,'K'
        print,'WaterVaporSatPressure> !!! valid range: 173.15K < T < 373.15K'
        return, es
    endif
endif
if ( STRUPCASE(phase) EQ 'ICE' ) then begin
    if ((TK LT 173.15) OR (TK GT 273.16)) then begin
        print,'WaterVaporSatPressure> !!! ERROR!'
        print,'WaterVaporSatPressure> !!! temperature is out of the validity range'
        print,'WaterVaporSatPressure> !!! for the saturation over ice: T=',TK,'K'
        print,'WaterVaporSatPressure> !!! valid range: 173.15K < T < 273.16K'
        return, es
    endif
endif


;; calculation for pure water vapor over liquid water [Pa or hPa]
if ( STRUPCASE(phase) EQ 'WATER' ) then begin
    es = WVSatPressureLiquidWater( DOUBLE(TK), punit )
;;  apply the correction term for moist air if necessary
    if (STRUPCASE(corr) EQ 'YES') then begin
        es = es * AirCorrFunWater(P, TK)
    endif
    return, es;; [Pa or hPa or mbar]
endif


;; calculation for pure water vapor over ice [Pa or hPa]
if ( STRUPCASE(phase) EQ 'ICE' ) then begin
    es = WVSatPressureIce( DOUBLE(TK), punit)
;;  apply the correction term for moist air if necessary
    if (STRUPCASE(corr) EQ 'YES') then begin
        es = es * AirCorrFunIce(P, TK)
    endif
    return, es;; [Pa or hPa or mbar]
endif


return, es; water vapor saturation pressure [Pa or hPa or mbar]
end


;;#################################################################################
