;; #################################################################################
;; functions to calculate saturation water pressure over liquid water
;; and ice
;; history:  alpha version,  2001-06-05, TKS
FUNCTION WVSatPressureLiquidWater, t
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
; from MPM93 model
; saturation water vapor pressure over ice,
; calculated according to Goff and Gratch formula.
; The saturation water vapor pressure is in units of Pa
; input is the temperature in Kelvin
FUNCTION WVSatPressureIce, t
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
;
function aii_file_exists, filename
;; checks whether the file filename exists and returns 1 for yes and 0
;; for no

YES = 1
NO  = 0

get_lun,unit
openr,unit,filename,error=err
free_lun,unit

if (err ne 0) then return,NO   $
else return,YES
end


PRO aii_readfile,filename,output
;; checks whether the is file present by using aii_file_exists,
;; compressed, or gzipped, and reads it if possible by using
;; read_datafile into output variable

;; do we have to uncompress?
com=-1

;; uncompressed
if aii_file_exists(filename) then begin
    filename = filename
    com = 0
endif 

;; compressed
if aii_file_exists(filename+'.Z') then begin
    filename = filename+'.Z'
    com = 1
endif 

;; gzipped
if aii_file_exists(filename+'.gz') then begin
    filename = filename+'.gz'
    com = 1
endif 

;; file not found
if com eq -1 then begin
    print,'Error: File not found: '+filename
    stop
endif

; read the file
if ( com ) then begin
	dummyname = './'+filename+'.arts.dummy.this.should.not.be.here'	
	print,'  decompressing...'
	spawn,'zcat '+filename+' > '+dummyname

	; read in from dummyname
        output=read_datafile(dummyname, /check)

	;remove the dummy matrix file
	spawn,'rm '+dummyname
endif else begin
	; read in from filename
        output=read_datafile(filename, /check)
endelse

end


Function aii_where_str,str1,str2,dim
;; similar to where, but works for string array str1 and str2
;;
;; returns -1 if not found, otherwise the index of str1 where entries
;; of str2 were found

;; size of str1, str2
ns1=(size(str1,/dime))[0]
ns2=(size(str2,/dime))[0]

;; allocate an int array, resize it later to dim
x=intarr(ns1)
dim=0

;; now loop the 2 arrays
for i=0,ns1-1 do begin
    for j=0,ns2-1 do begin
        if str1[i] eq str2[j] then begin
            x[dim]=i
            dim=dim+1
        endif
    endfor
endfor

;; now resize and return
if dim eq 0 then x=-1 else x=x[0:dim-1]
return, x
end


PRO plot_vmr_per_tg, jobname=jobname, $
                     add_to_title=add_to_title,$
                     avoid_tg=avoid_tg, $
                     color=color, $
                     cm=cm,$
                     xrange=xrange, $
                     yrange=yrange,$
                     nostamp=nostamp, $
                     jobdir=jobdir, $
                     punit=punit, $
                     tunit=tunit, $
                     vmrunit=vmrunit, $
                     altunit=altunit, $
                     plotfilename=plotfilename,$
                     plotfileformat=plotfileformat

;; reads an arts abs_per_tg file, frequency file, and altitude file ,
;; and plots the different tag group absorption, the tag groups are
;; identified by reading the arts controlfile. procedure does not plot
;; zero absorption of tag groups. Procedure can handle compressed and
;; gzipped files without uncompressing them.
;;
;; INPUT KEYWORDS:
;;
;;     jobname         : jobname of calculation
;;     jobdir          : sting containing the directory where the
;;                       jobfiles are located 
;;     add_to_title    : string-is added to default title
;;     avoid_tg        : array of strings-which tag groups to plot
;;     color           : int-make color plot (actually always color)
;;     nostamp         : int-produce no stamp info on output
;;     pressure        : selects pressure information
;;                       Possible values are: hPa, mbar, bar, Pa
;;     altunit        :  selects altitude information
;;                       Possible values are: km, m
;;     temperature     : selects temperature information
;;                       Possible values are: K, C
;;     vmrunit         : selection of the absorption units
;;                       Possible units are:
;;                       VMR, kg/m3
;;     plotfilename    : string containing output file name
;;     plotfileformat  : integer variable for the output file format.
;;                       Possible output files formats are
;;                          1: Postscript portrait mode, 
;;                          2: Postscript landscape mode, 
;;                          3: encapsulated Postscript portrait mode,
;;                          4: encapsulated Postscript landscape mode,
;;                          5: window
;; HISTORY:
;;     2001-06-04 : 1. version TKS
;;
;;
;; read the tag groups from controlfile
read_tag_groups,jobdir+jobname+'.arts','tgsDefine',tg
print, 'tag groups =', tg
;; total number of tag groups
tot_tg= (size(tg,/DIMENSIONS))[0]

;; set directory name apropriate
if keyword_set(jobdir) then begin 
    jobdirname = jobdir+jobname
endif else begin
    jobdirname = jobname
endelse

;; now reads the absorption per tag, and the altitudes and
;; frequencies, only ascii arrays are allowed currently
if not keyword_set(read) then begin
;    aii_readfile,jobdirname+'.f_mono.aa',f
;    aii_readfile,jobdirname+'.abs_per_tg.aa',abs
    aii_readfile,jobdirname+'.t_abs.aa',temperature
    aii_readfile,jobdirname+'.z_abs.aa',alt
    aii_readfile,jobdirname+'.p_abs.aa',pre
    aii_readfile,jobdirname+'.vmrs.aa',vmr
endif

;; are certain tag groups selected or should all be plotted?
;; we use min_arr values of zero for avoiding
if keyword_set(avoid_tg) then begin
    avoid_tg_in = aii_where_str(tg1,avoid_tg,dim)
    if dim eq 0 then begin
        print,'Error: No matching avoid_tg found.'
    endif else begin
        min_arr[avoid_tg_in] = 0
    endelse
endif


;; reading pressure information
string_title_p = ' '
if keyword_set(punit) then begin
    print,'Reading pressure grid data...'
    aii_readfile,jobdirname+'.p_abs.aa',p
    case punit of
        'Pa'   : begin
                 ppscale = 1.000
                 ppunit = '[Pa]'
                 end
        'hPa'  : begin
                 ppscale = 1.000/100.000
                 ppunit = '[hPa]'
                 end
        'kPa'  : begin
                 ppscale = 1.000/1000.000
                 ppunit = '[kPa]'
                 end
        else:    begin
                 ppscale = 1.000
                 ppunit = '[Pa]'
                 end
   endcase
endif


;; altitude unit selection
if keyword_set(altunit) then begin
    case altunit of
        'km'  : begin
                 zscale = 1.000/1000.0
                 zunit = '[km]'
                 end
        'm'   : begin
                 zscale = 1.000
                 zunit = '[m]'
                 end
        else:    begin
                 zscale = 1.000/1000.0
                 zunit = '[km]'
                 end
    endcase
endif else begin
    zscale = 1.000000
    zunit = 'VMR'
endelse
for i = 0, N_ELEMENTS(alt)-1 do begin
    alt[i] = alt[i] * zscale
endfor

;; vmr unit selection
if keyword_set(vmrunit) then begin
    case vmrunit of
        'VMR'  : begin
                 vmrscale = 1.000
                 vmrunit = '[1]'
                 end
        'kg/m3': begin
                 vmrscale = 1.000
                 vmrunit = '[kg/m!U3!N]'
                 end
        else:    begin
                 vmrscale = 1.000000
                 vmrunit = '[1]'
                 end
    endcase
endif else begin
    vmrscale = 1.000000
    vmrunit = 'VMR'
endelse


;; reading temperature information
string_title_T = ' '
if keyword_set(tunit) then begin
    print,'Reading temperature grid data...'
    aii_readfile,jobdirname+'.t_abs.aa',T
    case tunit of
        'K' : begin
              ttsub = 0.000
              ttunit = '[K]'
              end
        'C' : begin
              ttsub = 273.15
              ttunit = '[C]'
              end
        else: begin
              ttsub = 0.000
              ttunit = '[K]'
              end
    endcase
endif


;; check for color output? -> always use colors
if not keyword_set(color) then color=1

;; save settings
P_ini = !P

;; use aii_plot_file for writing into plot output file
if not keyword_set(plotfilename)   then plotfilename=jobname+'_plot'
if not keyword_set(plotfileformat) then plotfileformat=2
aii_plot_file, action='begin', fname=plotfilename, fformat=plotfileformat

;; create color and lineshape arrays
colors=intarr(tot_tg)
ls=intarr(tot_tg)
for i=0,tot_tg-1 do begin
    colors[i] = i mod 29    ;; we have 29 different colors in mycolor
    ls[i]     = i mod  5    ;; we have  5 different line styles
endfor

;; set line thickness
thick      = 3.5                ; standard value of thick
slidethick = 5.0                 ; value of thick to use for slides

;; 6 plots on a single page
!P.multi  = 0
posvec = [[0.2, 0.7, 0.45, 0.9], $
          [0.6, 0.7, 0.85, 0.9], $
          [0.2, 0.4, 0.45, 0.6], $
          [0.6, 0.4, 0.85, 0.6], $
          [0.2, 0.1, 0.45, 0.3], $
          [0.6, 0.1, 0.85, 0.3]]
!P.multi  = [0,2,3]


;; loop over all tags
index = 0 & posi = 0
for k = 0,N_ELEMENTS(tg)-1 do begin ; loop over all tags
;;  construct title of the plot:
    total_title= 'species="'+tg[k]+'"'

    s1='xmin = min(vmr.mat'+string(k,format='(I0)')+'[0,0:N_ELEMENTS(alt)-1])' 
    r1 = execute(s1)
    s1='xmax = max(vmr.mat'+string(k,format='(I0)')+'[0,0:N_ELEMENTS(alt)-1])' 
    r1 = execute(s1)
    if ((xmin/xmax) GT 0.10) then begin
        xmin = xmin *  0.100
        xmax = xmax * 10.000
    endif

    ax = 'volume mixing ratio '+vmrunit
    al = 'linestyle = ls[index], '
    if ( STRPOS(tg[k],'liquidcloud') GE 0) then begin
        xmin = xmax * 0.0100
        xmax = xmax * 10.00000
        ax = 'liquid water density '+'[kg/m!U3!N]'
    endif
    if ( STRPOS(tg[k],'icecloud') GE 0) then begin
        xmin = xmax * 0.00100
        xmax = xmax * 10.00000
        ax = 'ice particle density '+'[kg/m!U3!N]'
    endif
    ay = 'altitude '+zunit
    if (posi GE 6) then begin
        posi = 0
        !P.multi  = [0,2,3]
    endif

;;  plot VMR vs. altitude
    s1='plot, /XLOG, '+$
          'vmr.mat'+string(k,format='(I0)')+'[0,0:N_ELEMENTS(alt)-1], '+$
          'alt[0:N_ELEMENTS(alt)-1], '+$
          'title=total_title, '+$
          'yrange=[min(alt[0:N_ELEMENTS(alt)-1]), max(alt[0:N_ELEMENTS(alt)-1])], '+$
          'xrange=[xmin, xmax], '+$
          'xtitle=ax, '+$
          'ytitle=ay, '+$
          'color=colors[index], '+$
          al+$
          'thick=thick, '+$
          'xstyle=1, '+$
          'ystyle=1, '+$
          'POSITION=[posvec[0:3,'+string(posi,format='(I0)')+']]'
    r1 = execute(s1)
    posi = posi + 1
endfor

;; plot temperature
if (posi GE 6) then begin
    posi = 0
    !P.multi  = [0,2,3]
endif
xmin = min(temperature[0,0:N_ELEMENTS(alt)-1])-10.0
xmax = max(temperature[0,0:N_ELEMENTS(alt)-1])+10.0
plot, $
  temperature[0,0:N_ELEMENTS(alt)-1],$
  alt[0:N_ELEMENTS(alt)-1],$
  title='Temperature profile',$
  yrange=[min(alt[0:N_ELEMENTS(alt)-1]), max(alt[0:N_ELEMENTS(alt)-1])],$
  xrange=[xmin, xmax],$
  xtitle='temperaure [K]',$
  ytitle='altitude '+zunit,$
  color=colors[0],$
  linestyle=ls[0],$
  thick=thick,$
  xstyle=1,$
  ystyle=1, $
  POSITION=[posvec[0:3,posi]]
posi = posi+ 1
if (posi GE 6) then begin
    posi = 0
    !P.multi  = [0,2,3]
endif

;; relative humidity calculation if H2O is calculated in RT.
rhfound = 1
for k = 0,N_ELEMENTS(tg)-1 do begin
    if ( STRPOS(tg[k],'H2O') GE 0) then begin
        print, 'H2O tag found for RH calculation: ', tg[k]
        rhfound = 0
        tt = dblarr(N_ELEMENTS(alt))
        pp = dblarr(N_ELEMENTS(alt))
        zz = dblarr(N_ELEMENTS(alt))
        vv = dblarr(N_ELEMENTS(alt))
        for z = 0, N_ELEMENTS(alt)-1 do begin
            s1='vv['+string(z,format='(I)')+'] = vmr.mat'+$
               string(k,format='(I0)')+'[0,'+string(z,format='(I)')+']' 
            r1 = execute(s1)
        endfor
        goto, firh
    endif
endfor
firh:
if (rhfound EQ 0) then begin
    RH = dblarr(N_ELEMENTS(alt))
    for z = 0, N_ELEMENTS(alt)-1 do begin
        es = 0.000
        if (tt[z] LT 273.16) then es = WVSatPressureLiquidWater(temperature[0,z])
        if (tt[z] LT 273.16) then es = WVSatPressureIce(temperature[0,z])
        RH[z] = 100.00 * pre[z] * vv[z] / es
        ;print, z,' T=', temperature[0,z], ', vmr= ', vv[z],', es= ',  es,', RH=', RH[z]
    endfor
endif
; plot RH:
plot, $
  RH[0:N_ELEMENTS(alt)-1],$
  alt[0:N_ELEMENTS(alt)-1],$
  title='relative humidity profile',$
  yrange=[min(alt[0:N_ELEMENTS(alt)-1]), max(alt[0:N_ELEMENTS(alt)-1])],$
  xrange=[0.0, 110.0],$
  xtitle='relative humidity [%]',$
  ytitle='altitude '+zunit,$
  color=colors[0],$
  linestyle=ls[0],$
  thick=thick,$
  xstyle=1,$
  ystyle=1, $
  POSITION=[posvec[0:3,posi]]

;; close plot output file
aii_plot_file, action='end', show='yes', print='no'

;; cp plot output file to job directory
FINDFILENAMEVEC = FINDFILE(plotfilename+'.*', count=numfilefound)
if (numfilefound EQ 1) then begin
    spawn, 'cp '+FINDFILENAMEVEC[0]+' '+jobdir+'/'+FINDFILENAMEVEC[0]
endif

;; restore settings
!P = P_ini

ende:

end
