; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
PRO plot_abs_per_tg_2, jobname,                          $
                       abs=abs,                          $
                       f=f,                              $
                       alt=alt,                          $
                       altitude=altitude,                $
                       read=read,                        $
                       add_to_title=add_to_title,        $
                       avoid_tg=avoid_tg,                $
                       color=color,                      $
                       cm=cm,                            $
                       xrange=xrange,                    $
                       yrange=yrange,                    $
                       nostamp=nostamp,                  $
                       jobdir=jobdir,                    $
                       pressure=pressure,                $
                       temperature=temperature,          $
                       absunit=absunit,                  $
                       plotfilename=plotfilename,        $
                       plotfileformat=plotfileformat,    $
                       plotsum=plotsum,                  $
                       plotyaxis=plotyaxis
;
;***************************************************************************
;; reads an arts abs_per_tg file, frequency file, and altitude file ,
;; and plots the different tag group absorption, the tag groups are
;; identified by reading the arts controlfile. procedure does not plot
;; zero absorption of tag groups. Procedure can handle compressed and
;; gzipped files without uncompressing them.
;;
;; INPUT:
;;     jobname         : jobname of calculation
;; 
;; OUTPUT: (only needed to use the read keyword)
;;     f               : frequency grid
;;     abs             : altitude per tag group
;;     alt             : altitude grid
;;
;; KEYWORDS:
;;     altitude        : double-make abs plot near that altitude, 
;;                       default: lowest altitude
;;     read            : int-do not read the f, abs, alt again
;;     add_to_title    : string-is added to default title
;;     avoid_tg        : array of strings-which tag groups to plot
;;     color           : int-make color plot (actually always color)
;;     cm              : int-use cm^-1 instead of GHz
;;     xrange          : array of double-plot range, considers cm^-1 setting
;;     yrange          : array of double-plot range
;;     nostamp         : int-produce no stamp info on output
;;     jobdir          : sting containing the directory where the
;;                       jobfiles are located 
;;     pressure        : selects pressure information
;;                       Possible values are: hPa, mbar, bar, Pa
;;     temperature     : selects temperature information
;;                       Possible values are: K, C
;;     absunit         : selection of the absorption units
;;                       Possible units are:
;;                       1/m, 1/cm, 1/km, dB/km, Np/km
;;     plotfilename    : string containing output file name
;;     plotfileformat  : integer variable for the output file format.
;;                       Possible output files formats are
;;                          1: Postscript portrait mode, 
;;                          2: Postscript landscape mode, 
;;                          3: encapsulated Postscript portrait mode,
;;                          4: encapsulated Postscript landscape mode,
;;                          5: window
;;     plotsum         : flag for plotting in addition the total
;;                       absorption:
;;                       Possible values are:
;,                          1: plot additionally the total absorption
;;                          2: plot only tag absorption
;; HISTORY:
;;     2001-12-04 TKS  alpha version created 
;;
;***************************************************************************


;; set directory name apropriate
if keyword_set(jobdir) then begin 
    jobdirname = jobdir+jobname
    outdir     = jobdir
endif else begin
    jobdirname = jobname
    outdir     = '/home/home01/tkuhn/ARTS'
endelse


;; additional plotting of total absorption? 
if keyword_set(plotsum) then begin
    if ((plotsum NE 1) AND (plotsum NE 2)) then plotsum = 2
endif else begin
    plotsum = 2
endelse


;; read the tag groups from controlfile
read_tag_groups,jobdir+jobname+'.arts','tgsDefine',tg


;; total number of tag groups
tot_tg= (size(tg,/DIMENSIONS))[0]


;; now reads the absorption per tag, and the altitudes and
;; frequencies, only ascii arrays are allowed currently
if not keyword_set(abs) then begin
    print,'Reading absorption data...'
    aa = jobdirname+'.abs_per_tg.aa'
    abs = aa_read_general(aa)
endif
nmats = N_ELEMENTS(abs[*,0,0]) ; # of tags
nrows = N_ELEMENTS(abs[0,*,0]) ; # of frequencies
ncols = N_ELEMENTS(abs[0,0,*]) ; # of altitudes
print, 'abs: nmats=',nmats,', rows=',nrows,', colums=',ncols 
;; check that the tag groups and the absoprtion per tag group have the
;; same dimension
if (nmats ne tot_tg) then begin
    print,'Error: absorption per tag group and tag goups have not'
    print,'       the same dimension'
    stop
endif


;; 6) now read frequency data 
if not keyword_set(f) then begin
    print,'Reading frequency data...'
    f   = aa_read_general(jobdirname+'.f_mono.aa')
endif
if (nrows ne N_ELEMENTS(f[*])) then begin
    print,'Error: absorption per tag has not the same numbers of'
    print,'       frequencies like the frequency vector'
    stop
endif


;; now read altitude data 
if not keyword_set(alt) then begin
    print,'Reading altitude data...'
    alt = aa_read_general(jobdirname+'.z_abs.aa')
endif
if (ncols ne N_ELEMENTS(alt[*])) then begin
    print,'Error: absorption per tag has not the same numbers of'
    print,'       altitudes like the altitude vector'
    stop
endif


;; reading pressure information
print,'Reading pressure grid data...'
p = aa_read_general(jobdirname+'.p_abs.aa')


;; reading temperature information
print,'Reading temperature grid data...'
T = aa_read_general(jobdirname+'.t_abs.aa')


;; select altitude, if not given lowest altitude is taken
string_title_alt = ' '
if not keyword_set(altitude) then begin
    ialtitude = 0
endif else begin
    print,'preferred altitude [km]:',altitude
    ialtitude = 0
    delta_h = 1.000E5
    for i = 0, N_ELEMENTS(alt)-1 do begin 
        hz = alt[i]
        if ( ABS(hz-(altitude*1000.0)) LT delta_h) then begin 
            delta_h = ABS(hz-(altitude*1000.0))
            ialtitude = i
        endif
    endfor
;;    ialtitude = (where(alt gt (altitude*1000.0)))[0]
    if (ialtitude LT 0) then begin
        ialtitude = 0
    endif
endelse
print, 'selected altitude: ',(1.0E-3*alt[ialtitude]),'km'
string_title_alt = string(alt[ialtitude]/1000.0,format='(F5.1)')+' km'


;; frequency unit selection
if keyword_set(cm) then begin
    unit=2.9979E10
    unitstr='[cm!u-1!n]'
endif else begin
    unit=1E9
    unitstr='[GHz]'
endelse
for i = 0,N_ELEMENTS(f)-1 do f[i] = f[i] / unit

;; absorption unit selection
if keyword_set(absunit) then begin
    case absunit of
        '1/km' :  begin
                  absscale = 1000.000000
                  absunit = 'km!U-1!N'
                  end
        '1/cm' :  begin
                  absscale =    0.010000
                  absunit = 'cm!U-1!N'
                  end
        'dB/km' : begin
                  absscale = 4342.944819
                  absunit = 'dB/km'
                  end
        'Np/km' : begin
                  absscale = 1000.000000
                  absunit = 'dB/km'
                  end
        else:     begin
                  absscale =    1.000000
                  absunit = 'm!U-1!N'
                  end
    endcase
endif else begin
    absscale = 1.000000
    absunit = 'm!U-1!N'
endelse
for i = 0,nmats-1 do begin 
    for j = 0,nrows-1 do begin 
        for k = 0,ncols-1 do begin 
            abs[i,j,k] = abs[i,j,k] / absscale
        endfor
    endfor
endfor


;; pressure unit selection
string_title_p = ' '
if keyword_set(pressure) then begin
    case pressure of
        'bar'  : begin
                 pscale = 1.000/100000.000
                 punit = 'bar'
                 end
        'mbar' : begin
                 pscale = 1.000/100.000
                 punit = 'mbar'
                 end
        'Pa'   : begin
                 pscale = 1.000
                 punit = 'Pa'
                 end
        'hPa'  : begin
                 pscale = 1.000/100.000
                 punit = 'hPa'
                 end
        'kPa'  : begin
                 pscale = 1.000/1000.000
                 punit = 'kPa'
                 end
        else:    begin
                 pscale = 1.000
                 punit = 'Pa'
                 end
    endcase
    for i = 0,N_ELEMENTS(p)-1 do p[i] = p[i] * pscale
    string_title_p = string(p[ialtitude],format='(F7.2)')+punit
endif


;; temperature unit selection
if keyword_set(temperature) then begin
    case temperature of
        'K' : begin
              Tsub = 0.000
              Tunit = 'K'
              end
        'C' : begin
              Tsub = 273.15
              Tunit = 'C'
              end
        else: begin
              Tsub = 0.000
              Tunit = 'K'
              end
    endcase
    for i = 0,N_ELEMENTS(T)-1 do T[i] = T[i]-Tsub
    string_title_T = string(T[ialtitude],format='(F5.1)')+Tunit
endif


;; are certain tag groups selected or should all be plotted?
tag_index = intarr(N_ELEMENTS(tg))
tag_index_max = 0
if keyword_set(avoid_tg) then begin
    for j = 0,N_ELEMENTS(tg)-1 do begin
        flag = 0
        for i = 0,N_ELEMENTS(avoid_tg)-1 do begin
            if (avoid_tg[i] EQ tg[j]) then flag = 1
        endfor
        if (flag EQ 0) then begin
            tag_index[tag_index_max] = j
            tag_index_max = tag_index_max + 1
        endif
    endfor
    tag_index = tag_index[0:tag_index_max-1]
endif else begin
    tag_index_max = N_ELEMENTS(tg)
    tag_index     = indgen(tag_index_max)
endelse


;; sort the tags after max. absorption at the plotting altitude
absmax = dblarr(tag_index_max)
for j = 0,tag_index_max-1 do begin
    i = tag_index[j]
    absmax[j] = max(abs[i,*,ialtitude])
endfor
absmax_index = SORT(absmax)
index = intarr(tag_index_max)
for j = 0,tag_index_max-1 do begin
    index[tag_index_max-1-j] = tag_index[absmax_index[j]]
endfor
tag_index = index[0:tag_index_max-1]
print, 'selected tags for printing: '
for i = 0,tag_index_max-1 do print, ' ',i,': ',tg[tag_index[i]]


;                       ---------------------------------
; =================== plot the absorption vs. frequency ===================
;                       ---------------------------------

;; make 1 plot per page
!P.multi    = [0,1,1]


;; settings for the plot
!X.MARGIN   = [0,0]
!Y.MARGIN   = [0,0]
!X.OMARGIN  = [0,0]
!Y.OMARGIN  = [0,0]
!P.REGION   = [0.0, 0.0, 0.0, 0.0]
!P.POSITION = [0.1, 0.1, 0.7, 0.9]
plotpos = !P.POSITION
!P.CHARTHICK = 5.0
!P.FONT      = 1
!P.CHARSIZE  = 1.5
!X.CHARSIZE  = 1
!Y.CHARSIZE  = 1
!P.THICK = 5
thick = !P.THICK
!X.THICK = 5
!Y.THICK = 5

 
;; get datum to write it on the top of the plot:
spawn,'date +"%y"',year
spawn,'date +"%m"',month
spawn,'date +"%d"',day
spawn,'date +"%H"',hour
spawn,'date +"%M"',minute
spawn,'whoami',who
datum = who+':20'+string(year, FORMAT='(A2)')+'-'+$
string(month, FORMAT='(A2)')+'-'+string(day, FORMAT='(A2)')


;; check for color output? -> always use colors
if not keyword_set(color) then color=1

;; save settings
P_ini = !P

;; use aii_plot_file for writing into plot output file
if not keyword_set(plotfilename)   then plotfilename=jobname+'_plot'
if not keyword_set(plotfileformat) then plotfileformat=2
aii_plot_file, action='begin', fname=plotfilename, fformat=plotfileformat


;; set line thickness
;thick      = 4.                ; standard value of thick
;slidethick = 5.                 ; value of thick to use for slides


;; create color and lineshape arrays
if (plotsum EQ 1) then nn=1 else nn=0
colors=intarr(tot_tg+nn)
ls=intarr(tot_tg+nn)
for i=0,tot_tg+nn-1 do begin
    colors[i] = i mod 29    ;; we have 29 different colors in mycolor
    ls[i]     = i mod  5    ;; we have  5 different line styles
endfor


;; add to title string
if not keyword_set(add_to_title) then add_to_title=''


;; construct title of the plot:
total_title= 'z='+string_title_alt+' / '+ $
             'p='+string_title_p+' / '+$
             'T='+string_title_T+'  '+$
             add_to_title


;; calculate the total absorption in necessary
if (plotsum EQ 1) then begin
    abstotal = dblarr(nrows)
    for i = 0,nmats-1 do begin
        for j = 0,nrows-1 do begin
            abstotal[j] = abs[i,j,ialtitude]
        endfor
    endfor
endif

;; apropriate x range of plot
if not keyword_set(xrange) then begin
  xmin = min(f)
  xmax = max(f)
endif else begin
    xmin = xrange[0]
    xmax = xrange[1]
endelse


;; apropriate y range of plot
if not keyword_set(yrange) then begin
    abs_min = fltarr(tag_index_max)
    abs_max = fltarr(tag_index_max)
    for j = 0,tag_index_max-1 do begin
        i = tag_index[j]
        abs_min[j] = min(abs[i,*,ialtitude])
        abs_max[j] = max(abs[i,*,ialtitude])
    endfor
    if (plotsum EQ 1) then begin
        ymin = min([abstotal[*], abs_min[*]])
        ymax = max([abstotal[*], abs_max[*]])
    endif else begin 
        ymin = min(abs_min[*])
        ymax = max(abs_max[*])
    endelse
    if (ymin GT ymax) then begin
        ya   = ymin
        ymin = ymax
        ymax = ya
    endif
    if (ymin LT 0.0) then begin
        ymin = ymin * 1.05
    endif else begin
        ymin = ymin / 1.05
    endelse
    if (ymax LT 0.0) then begin
        ymax = ymax / 1.05
    endif else begin
        ymax = ymax * 1.05
    endelse
endif else begin
    ymin = yrange[0]
    ymax = yrange[1]
endelse

print,'x/y-axis range:',xmin,',',xmax,'/',ymin,',',ymax

;; set plot style 'log' or 'lin' for y-axis
if NOT keyword_set(plotyaxis) then plotyaxis='lin'
if ((plotyaxis NE 'lin') AND (plotyaxis NE 'log')) then plotyaxis='lin'

;; set frame of the plot
IF (plotyaxis EQ 'lin') then begin
    plot, f[0:nrows-1] , abs[0,0:nrows-1,ialtitude], $
          /NORMAL, $
          title=total_title,                         $
;          xcharsize=1.25,                            $
;          ycharsize=1.25,                            $
          xrange=[xmin, xmax],                       $
          yrange=[ymin, ymax],                       $
          xtitle='frequency '+unitstr,               $
          ytitle='absorption ['+absunit+']',         $
          color=colors[0],                           $
          linestyle=ls[0],                           $
          xstyle=1,                                  $
          ystyle=1,                                  $
          /nodata
endif
IF (plotyaxis EQ 'log') then begin
;; set frame of the plot
    plot, /YLOG, f[0:nrows-1] , abs[0,0:nrows-1,ialtitude], $
          /NORMAL, $
          title=total_title,                         $
;          xcharsize=1.25,                            $
;          ycharsize=1.25,                            $
          xrange=[xmin, xmax],                       $
          yrange=[ymin, ymax],                       $
          xtitle='frequency '+unitstr,               $
          ytitle='absorption ['+absunit+']',         $
          color=colors[0],                           $
          linestyle=ls[0],                           $
          xstyle=1,                                  $
          ystyle=1,                                  $
          /nodata
endif

;; scale absorption with apropriate unit factor and summ up all
;; tag absorption to total absorption if plotsum=0
col_start = 0
if (plotsum EQ 1) then begin
    oplot, f[0:nrows-1] , abstotal[0:nrows-1], $
           color=colors[col_start], $
           linestyle=ls[col_start]
    col_start = col_start + 1
endif


FOR i = 0,tag_index_max-1 DO BEGIN
    j = tag_index[i]
    oplot, f[0:nrows-1] , abs[j,0:nrows-1,ialtitude], $
           color=colors[col_start+i], $
           linestyle=ls[col_start+i]
ENDFOR


; print datum and user name:
xyouts, plotpos[0], plotpos[1]-0.10, datum, CHARSIZE=0.75, CHARTHICK=1.0, /NORMAL


;; legend of the plot
if (plotsum EQ 1) then begin
    charlegendarray = strarr(tag_index_max+1)
    charlegendarray[0] = 'total absorption'
    for i = 0,tag_index_max-1 do charlegendarray[i+1] = tg[tag_index[i]]
endif else begin
    charlegendarray = strarr(tag_index_max)
    for i = 0,tag_index_max-1 do charlegendarray[i] = tg[tag_index[i]]
endelse
ch_length = N_ELEMENTS(charlegendarray)
aii_plot_legend, charlegendarray[0:tag_index_max-1], $
               box=0,$
               usersym=usersym, $
               spacing=2.0,$
               pspacing=2.0, $
               psym=intarr(ch_length), $
               colors=colors[0:ch_length-1],$
               line=ls[0:ch_length-1], $
               position=[plotpos[2]-0.01,plotpos[3]+0.02], thick=thick*1.4


;; restore settings
P_ini = !P    

;
; --- end of procedure -----------------------------------------------------
ende:
;; close plot output file
aii_plot_file, action='end', show='yes', print='no', $
               outdir=outdir
;
END
; ==========================================================================
; ##########################################################################
; ==========================================================================
