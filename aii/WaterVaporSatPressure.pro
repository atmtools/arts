;;#################################################################################
;;
FUNCTION WVSatPressureLiquidWater, t
; ---------------------------------------------------------------------
;; functions to calculate saturation water pressure over liquid water
;; and ice
;; history:  alpha version,  2001-06-05, TKS
; ---------------------------------------------------------------------
  ;  Rosenkranz approximation
  ;  COMPUTES SATURATION H2O VAPOR PRESSURE (OVER LIQUID)
  ;  USING LIEBE'S APPROXIMATION (CORRECTED)
  ;  input : T  in Kelvin
  ;  ouput : es in Pa
  ;  PWR 4/8/92
  ;
   TH       = 300.0 / t;
   es_PWR98 = 100.00 * 35.3 * EXP(22.64*(1.-TH)) * TH^5; saturation pressure in Pa
  ;
; ---------------------------------------------------------------------

  ; MPM93 calculation
   theta    = 373.16 / t;

   a = (11.344*(1.00-(1.00/theta)))
   b = a * ALOG(10.00)
   c = (-3.49149*(theta-1.00))
   d = c * ALOG(10.00)

   exponent = ( -7.90298 * (theta-1.000) +       $
		5.02808 * ALOG10(theta) -        $
	        1.3816e-7 * ( EXP(b) - 1.000 ) + $
	        8.1328e-3 * ( EXP(d) - 1.000) +  $
	        ALOG10(1013.246) );

   es_MPM93 = 100.000 * EXP(exponent * ALOG(10.00))

   return, es_MPM93 ; saturation pressure in Pa
end
;
; #################################################################################
;
FUNCTION WVSatPressureIce, t
; ---------------------------------------------------------------------
; from MPM93 model
; saturation water vapor pressure over ice,
; calculated according to Goff and Gratch formula.
; The saturation water vapor pressure is in units of Pa
; input is the temperature in Kelvin
; ---------------------------------------------------------------------
  ; MPM93 calculation
   theta    = 273.16 / t

   exponent = (-9.09718  * (theta-1.000) -          $
	        3.56654  * ALOG10(theta)  +         $
		0.876793 * (1.000-(1.000/theta)) +  $
	        ALOG10(6.1071) )

   es_MPM93 = 100.000 * EXP(exponent * ALOG(10.00))

  return, es_MPM93 ; saturation pressure in Pa
  end
;
; #################################################################################
;;
FUNCTION WaterVaporSatPressure, LWC, IWC, T
;; ---------------------------------------------------------------------
;; input:
;;   LWC : liquid water content  [g/m3]
;;   IWC : ice  water content    [g/m3]
;;     T : temperature           [T]
;;
;; output:
;;  es   : water vapor saturtation pressure [Pa]
;;
;; history:  alpha version,  2001-06-05, TKS
;; ---------------------------------------------------------------------
;;
es = -1.000e6 ; control value

;; liquid water cloud present:
;; calculation is for water vapor over liquid water
if ( (LWC GT 0.00) AND (IWC LE 0.00) ) then begin
    es = WVSatPressureLiquidWater( T )
endif


; ice cloud  present:
;; calculation is for water vapor over ice
if ( (LWC LE 0.00) AND (IWC GT 0.00) ) then begin
    es = WVSatPressureIce( T )
endif


; both cloud types are present, take then the lower saturation pressure
if ( (LWC GT 0.00) AND (IWC GT 0.00) ) then begin
    es1 = WVSatPressureLiquidWater( T )
    es2 = WVSatPressureIce( T )
    es = min([es1, es2])
endif

;; no cloud present:
;; below freezing point (0 deg Celcius) take saturation pressure
;; 'over ice', above take it from 'over liquid water'
if ( (LWC LE 0.00) AND (IWC LE 0.00) ) then begin
    if (T GE 273.15) then es = WVSatPressureLiquidWater( T )
    if (T LT 273.15) then es = WVSatPressureIce( T )
endif

return, es; water vapor saturation pressure [Pa]
end
;;#################################################################################
