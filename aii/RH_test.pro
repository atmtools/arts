;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
; H2O vapor saturation pressure over liquid water
;
FUNCTION  ESAT, T ;; = temperature in K

ON_ERROR,2
IF max(T) LT 105. THEN T0=273.16 ELSE T0=0.0

; Formula with T = temperature in K
;    esat = exp( -6763.6/(T+T0) - 4.9283*alog((T+T0)) + 54.2190 )

; Formula close to that of Magnus, 1844 with temperature TC in Celsius
;    ESAT = 6.1078 * EXP( 17.2693882 * TC / (TC + 237.3) ) ; TC in Celsius

; or Emanuel's formula (also approximation in form of Magnus formula, 1844)
;    esat=6.112*EXP(17.67*TC/(243.5+TC))

; WMO reference formula is that of Goff and Gratch (1946), slightly
; modified by Goff in 1965:

e1=1013.250
TK=273.16
esat=e1*10^(10.79586*(1-TK/(T+T0))-5.02808*alog10((T+T0)/TK)+$
            1.50474*1e-4*(1-10^(-8.29692*((T+T0)/TK-1)))+$
            0.42873*1e-3*(10^(4.76955*(1-TK/(T+T0)))-1)-2.2195983)

RETURN, (100.00*ESAT) ; [Pa]

END
;;
;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;;
FUNCTION PlotSetupA, text, n
;
;; flag: 0=ok, 1=error occured
ok = 0
;
;; use aii_plot_file for writing into plot output file
aii_plot_file, action='begin', fname=text, fformat=n
;
!P.MULTI     = 0
!P.FONT      = 1
!P.CHARSIZE  = 1.5
!X.CHARSIZE  = 1
!Y.CHARSIZE  = 1
!P.THICK     = 8
!X.THICK     = 5
!Y.THICK     = 5
!X.CHARSIZE  = 1.5
!Y.CHARSIZE  = 1.5
!P.CHARSIZE  = 1.5
!P.CHARTHICK = 4
;!X.MARGIN    = [0,0]
;!Y.MARGIN    = [0,0]
;
DEVICE, /times, font_size=12
;
!P.MULTI     = [0,1,1]
!P.POSITION  = [0.2, 0.2, 0.9, 0.8]
;
RETURN, ok
END
;
; ============================================================================
;
FUNCTION PlotSetupB, a, b, c
;
ok = 1  ;; flag: 0=ok, 1=error occured
aii_plot_file, action='end', show=a, print='no', $
               outdir=b, writedate=c
ok = 0
;
RETURN, ok
END



;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



PRO RH_test

close, /all

;; save settings
P_ini = !P

;; create color and lineshape arrays
ncolors = 25
colors = intarr(ncolors)
lsy    = intarr(ncolors)
for i=0, ncolors-1 do begin
    colors[i] = i mod 29    ;; we have 29 different colors in mycolor
    lsy[i]    = i mod  5    ;; we have  5 different line styles
endfor

;;  set general size aspects of characters and symbols
myxcharsize = 1.5
myycharsize = 1.5
mylcharsize = 1.5
mysymsize = 0.75

PSoutdir   = '.'
showPSplot = 'no'
WriteDate  = 'yes'



atmprof= ['tropical', $
          'midlatitude-summer', 'midlatitude-winter', $
          'subarctic-summer',   'subarctic-winter']
natm = N_ELEMENTS(atmprof)

for iatm = 0,natm-1 do begin

    FASCOD_ptz   = aa_read_general( '/pool/lookup2/arts-data/atmosphere/fascod/'+atmprof[iatm]+'.tz.aa' )
    N_FASCOD_ptz = N_ELEMENTS(FASCOD_ptz[0,*,0])
    
    FASCOD_H2O   = aa_read_general( '/pool/lookup2/arts-data/atmosphere/fascod/'+atmprof[iatm]+'.H2O.aa' )
    N_FASCOD_H2O = N_ELEMENTS(FASCOD_ptz[0,*,0])
    
    RH = dblarr(8, N_FASCOD_ptz)
    
    j = 0
    for i = 0,N_FASCOD_ptz-1 do begin
        
        p_Pa = FASCOD_ptz[0,i,0]
        TK   = FASCOD_ptz[0,i,1]
        z_km = (FASCOD_ptz[0,i,2]*1.00000e-3)
        
        if ( (ABS(FASCOD_H2O[0,i,0]-p_Pa)/p_Pa) LT 0.01 ) then begin
            
            ph2o = FASCOD_H2O[0,i,0] * FASCOD_H2O[0,i,1] ;; [Pa]
            
            es1 = WaterVaporSatPressure( TK, phase='WATER', corr='no'  )
;            es1 = WaterVaporSatPressure( TK, phase='WATER', punit='Pa', corr='no'  )
            es2 = WaterVaporSatPressure( TK, PTOT=p_Pa, phase='WATER', punit='Pa', corr='yes' )
            
            RH[0,j] = (FASCOD_H2O[0,i,0]*1.00000e-2)            ;; pressure [hPa]
            RH[1,j] = TK                                        ;; temperature [K]
            RH[2,j] = z_km                                      ;; altitude [km]
            RH[3,j] = 1.000E+2 * (ph2o / es1)                   ;; rel. humidity [%]
            RH[4,j] = 1.000E+2 * (ph2o / es2)                   ;; rel. humidity [%]
            RH[5,j] = 1.000E+2 * (ph2o / ESAT(TK))              ;; rel. humidity [%]
            RH[6,j] = es2/es1                                   ;; sat pressure ratio
            RH[7,j] = ESAT(TK) / es2                            ;; sat pressure ratio
            
            j = j + 1
            
        endif
        
    endfor
    
    
    plotok   = PlotSetupA( 'RH_test_A_'+atmprof[iatm], 1 )
    if (plotok NE 0) then begin
        print,'!!! ERROR occured while opening PS file'
        print,'!!! terminate without action'
        goto, endloop
    endif
    
;; plot Sonntag without correction
    plot,                                  $
      RH[3,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      /NORMAL,                             $
      title='Sonntag formula without moist air correction ('+atmprof[iatm]+')', $
      xcharsize=myxcharsize,               $
      ycharsize=myycharsize,               $
      xrange=[0.00, 120.0],                $
      yrange=[0.00, 15.0],                 $
      xtitle='relative humidity  (no corr.) [%]',          $
      ytitle='altitude [km]',              $
      color=colors[0],                     $
      thick=5,                             $
      linestyle=0,                         $
      xstyle=1,                            $
      ystyle=1,                            $
      /NODATA
    
    oplot,                                 $
      RH[3,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      thick=5,                             $
      color=colors[1],                     $
      linestyle=1
    
    
    
;; plot Sonntag with correction
    plot,                                  $
      RH[4,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      /NORMAL,                             $
      title='Sonntag formula with moist air correction ('+atmprof[iatm]+')', $
      xcharsize=myxcharsize,               $
      ycharsize=myycharsize,               $
      xrange=[0.00, 120.0],                $
      yrange=[0.00, 15.0],                 $
      xtitle='relative humidity (with corr.) [%]',          $
      ytitle='altitude [km]',              $
      color=colors[0],                     $
      thick=5,                             $
      linestyle=0,                         $
      xstyle=1,                            $
      ystyle=1,                            $
      /NODATA
    
    oplot,                                 $
      RH[4,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      thick=5,                             $
      color=colors[2],                     $
      linestyle=2
    
    
;; plot ESAT (NASA) function
    plot,                                  $
      RH[5,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      /NORMAL,                             $
      title='NASA (ESAT) formula ('+atmprof[iatm]+')', $
      xcharsize=myxcharsize,               $
      ycharsize=myycharsize,               $
      xrange=[0.00, 120.0],                $
      yrange=[0.00, 15.0],                 $
      xtitle='relative humidity  (ESAT) [%]',  $
      ytitle='altitude [km]',              $
      color=colors[0],                     $
      thick=5,                             $
      linestyle=0,                         $
      xstyle=1,                            $
      ystyle=1,                            $
      /NODATA
    
    oplot,                                 $
      RH[5,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      thick=5,                             $
      color=colors[3],                     $
      linestyle=3
    
    
    
;;  plot sat. pressure ratios 
    plot,                                  $
      RH[6,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      /NORMAL,                             $
      title='sat. pressure ratios ('+atmprof[iatm]+')', $
      xcharsize=myxcharsize,               $
      ycharsize=myycharsize,               $
      xrange=[0.95, 1.05],                 $
      yrange=[0.0, 15.0],                  $
      xtitle='sat. pressure ratios [1]',   $
      ytitle='altitude [km]',              $
      color=colors[0],                     $
      thick=5,                             $
      linestyle=1,                         $
      xstyle=1,                            $
      ystyle=1,                            $
      /NODATA
    
;; Sonntag with corr. / Sonntag without corr.
    oplot,                                 $
      RH[6,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      color=colors[4],                     $
      thick=5,                             $
      linestyle=4
    
;; NASA / Sonntag with corr.
    oplot,                                 $
      RH[7,0:N_FASCOD_ptz-1],              $
      RH[2,0:N_FASCOD_ptz-1],              $
      color=colors[5],                     $
      thick=5,                             $
      linestyle=5
    
    xyouts, 0.955, 1.50, $
      'Sonntag with/without corr', charsize=1.50, color=colors[4], /DATA

    xyouts, 0.955, 3.00, $
      'NASA / Sonntag with corr', charsize=1.50, color=colors[5], /DATA
    
    plots, [1.0, 1.0], [0.0, 15.0],  color=colors[0], linestyle=0
    
    
    
endloop:
;; CLOSE PS FILE:
    ok  = PlotSetupB( showPSplot, PSoutdir, WriteDate)
    
    
endfor




ende:
close, /all
!P = P_ini

END
