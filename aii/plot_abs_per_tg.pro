; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
PRO plot_abs_per_tg, jobname, f, abs, alt, $
                     altitude=altitude, $
                     read=read,$
                     add_to_title=add_to_title,$
                     avoid_tg=avoid_tg, $
                     color=color, $
                     cm=cm,$
                     xrange=xrange, $
                     yrange=yrange,$
                     nostamp=nostamp, $
                     jobdir=jobdir, $
                     pressure=pressure, $
                     temperature=temperature, $
                     absunit=absunit, $
                     plotfilename=plotfilename,$
                     plotfileformat=plotfileformat, $
                     plotsum=plotsum
;;;;;;;;;;;;;;;
;+
;NAME:
;      plot_abs_per_tg
;
;PURPOSE:
; reads an arts abs_per_tg file, frequency file, and altitude file ,
; and plots the different tag group absorption, the tag groups are
; identified by reading the arts controlfile. procedure does not plot
; zero absorption of tag groups. Procedure can handle compressed and
; gzipped files without uncompressing them.
;
; INPUT:
;     jobname         : jobname of calculation
; 
; OUTPUT: (only needed to use the read keyword)
;     f               : frequency grid
;     abs             : altitude per tag group
;     alt             : altitude grid
;
; KEYWORDS:
;     altitude        : double-make abs plot near that altitude, 
;                       default: lowest altitude
;     read            : int-do not read the f, abs, alt again
;     add_to_title    : string-is added to default title
;     avoid_tg        : array of strings-which tag groups to plot
;     color           : int-make color plot (actually always color)
;     cm              : int-use cm^-1 instead of GHz
;     xrange          : array of double-plot range, considers cm^-1 setting
;     yrange          : array of double-plot range
;     nostamp         : int-produce no stamp info on output
;     jobdir          : sting containing the directory where the
;                       jobfiles are located 
;     pressure        : selects pressure information
;                       Possible values are: hPa, mbar, bar, Pa
;     temperature     : selects temperature information
;                       Possible values are: K, C
;     absunit         : selection of the absorption units
;                       Possible units are:
;                       1/m, 1/cm, 1/km, dB/km, Np/km
;     plotfilename    : string containing output file name
;     plotfileformat  : integer variable for the output file format.
;                       Possible output files formats are
;                          1: Postscript portrait mode, 
;                          2: Postscript landscape mode, 
;                          3: encapsulated Postscript portrait mode,
;                          4: encapsulated Postscript landscape mode,
;                          5: window
;     plotsum         : flag for plotting in addition the total
;                       absorption:
;                       Possible values are:
;,                          1: plot additionally the total absorption
;                          2: plot only tag absorption
; HISTORY:
;     2001-01-25 : Created AVE
;     2001-04-05 : additional functionalities TKS
;-
;;;;;;;;;;;;;;;;;;
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

;; set directory name apropriate
if keyword_set(jobdir) then begin 
    jobdirname = jobdir+jobname
endif else begin
    jobdirname = jobname
endelse

;; now reads the absorption per tag, and the altitudes and
;; frequencies, only ascii arrays are allowed currently
if not keyword_set(read) then begin
    print,'Reading data...'
    aii_readfile,jobdirname+'.f_mono.aa',f
    aii_readfile,jobdirname+'.abs_per_tg.aa',abs
    aii_readfile,jobdirname+'.z_abs.aa',alt
endif

;; check that the tag groups and the absoprtion per tag group have the
;; same dimension
if n_tags(abs) ne tot_tg then begin
    print,'Error: absorption per tag group and tag goups have not'
    print,'       the same dimension'
    stop
endif

;; altitude, if not given lowest altitude is taken
string_title_alt = ' '
if not keyword_set(altitude) then begin
    ialtitude = 0
endif else begin
    ialtitude = 0
    delta_h = 1000.0
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
        print, 'selected altitude: ',alt[ialtitude],'m'
    endif
endelse
string_title_alt = string(alt[ialtitude]/1000.0,format='(F6.2)')+' km'

;; rearrange the absorption and tag groups according to their
;; absorption magnitude, and get max and min values for plot range
sort_abs,abs,tg,ialtitude,abs1,tg1,max_arr,min_arr

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

;; frequency unit selection
if keyword_set(cm) then begin
    unit=30.0
    unitstr='[cm!u-1!n]'
endif else begin
    unit=1.0
    unitstr='[GHz]'
endelse

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

;; reading pressure information
string_title_p = ' '
if keyword_set(pressure) then begin
    print,'Reading pressure grid data...'
    aii_readfile,jobdirname+'.p_abs.aa',p
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
    string_title_p = string(p[ialtitude],format='(F8.3)')+punit
endif

;; reading temperature information
string_title_T = ' '
if keyword_set(temperature) then begin
    print,'Reading temperature grid data...'
    aii_readfile,jobdirname+'.t_abs.aa',T
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
    string_title_T = string(T[ialtitude],format='(F7.3)')+Tunit
endif

;; apropriate x and y range of plot
if not keyword_set(yrange) then begin
    ymin = min( min_arr(where(min_arr ne 0)) ) / 10.0
    ymax = max(max_arr)*10.0
endif
if not keyword_set(xrange) then begin
  xmin = min(f/1E9/unit)
  xmax = max(f/1E9/unit)
endif

;; check for color output? -> always use colors
if not keyword_set(color) then color=1

;; save settings
P_ini = !P

;; use aii_plot_file for writing into plot output file
if not keyword_set(plotfilename)   then plotfilename=jobname+'_plot'
if not keyword_set(plotfileformat) then plotfileformat=2
aii_plot_file, action='begin', fname=plotfilename, fformat=plotfileformat

;; set line thickness
thick      = 2.5                ; standard value of thick
slidethick = 5.                 ; value of thick to use for slides

;; create color and lineshape arrays

if (plotsum EQ 1) then nn=1 else nn=0
colors=intarr(tot_tg+nn)
ls=intarr(tot_tg+nn)
for i=0,tot_tg+nn-1 do begin
    colors[i] = i mod 29    ;; we have 29 different colors in mycolor
    ls[i]     = i mod  5    ;; we have  5 different line styles
endfor

;; make 2 plots, one with the legend, the other one with the data
;; data plot
!P.multi  = [0,2,1]
!X.margin = [10,-10]
!Y.margin = [4,4]

;; add to title string
if not keyword_set(add_to_title) then add_to_title=''

;; construct title of the plot:
total_title= 'z='+string_title_alt+'/'+ $
             'p='+string_title_p+'/'+$
             'T='+string_title_T+'  '+$
             add_to_title

;; scale absorption with apropriate unit factor and summ up all
;; tag absorption to total absorption if plotsum=0
if (plotsum EQ 1) then begin
    abstotal = dblarr(N_ELEMENTS(abs1.mat0[ialtitude,*]))
    for i = 0,N_ELEMENTS(abs1.mat0[ialtitude,*])-1 do begin
        abstotal[i] = 0.0
    endfor
endif
nameoftags = TAG_NAMES(abs1)
for k = 0,N_TAGS(abs1)-1 do begin ; loop over all tags
    ;; loop over all frequencies
    for i = 0,N_ELEMENTS(abs1.mat0[ialtitude,*])-1 do begin
        ;; single tag group manipulation
        s1 = 'abs1.'+nameoftags[k]+'[ialtitude,'+string(i,format='(I)')+'] = '+$
             'abs1.'+nameoftags[k]+'[ialtitude,'+string(i,format='(I)')+'] * '+$
             string(absscale,format='(F12.6)')
        r1 = execute(s1)
        if (plotsum EQ 1) then begin
             s2 = 'abstotal['+string(i,format='(I)')+'] = '+$
                  'abstotal['+string(i,format='(I)')+'] + '+$
                  'abs1.'+nameoftags[k]+'[ialtitude,'+string(i,format='(I)')+']'
             r2 = execute(s2)
        endif
    endfor
endfor
;; plot absorption vs. frequency
plot, f/1E9/unit, $
      abs1.mat0[ialtitude,*], $
      title=total_title, $
      yrange=[ymin*absscale, ymax*absscale], $
      ytype=1,$
      xrange=[xmin, xmax], $
      xtitle='Frequency '+unitstr,$
      ytitle='Absorption ['+absunit+']',$
      color=colors[0], $
      linestyle=ls[0], $
      thick=thick,$
      xstyle=1, $
      ystyle=1, $
      /nodata;;,charsize=0.7*!P.charsize

index=0 ;; count only the ones we plot
nozeros=where(min_arr ne 0)
for j=0,tot_tg-1 do begin
    ;; check that this absorption tag group is not zero
    dum=where(nozeros eq j,count)
    if count eq 1 then begin
        s1='oplot, f/1E9/unit, abs1.mat'+string(j,format='(I0)')+$
          '[ialtitude,*],linestyle='+string(ls[index],format='(I0)')+$
          ', color='+string(colors[index],format='(I0)')+', thick=2*thick'
        index=index+1
        r1 = execute(s1)
    endif else begin
        print,'plot_abs_per_tg> Avoided or found absorption is zero, no plotting: ',tg1[j]
    endelse
endfor

;; plot total absorption if necessary
if (plotsum EQ 1) then begin
    oplot, f/1E9/unit, $
           abstotal[*], $
           linestyle=ls[index], $
           color=colors[index], $
           thick=2*thick
    index = index + 1
    abstotchar = 'total absorption'
endif

;; legend plot
!X.margin = [-130,5]
!Y.margin = [4,4]
plot,[0,0],xstyle=4,ystyle=4,/nodata,/noerase,/noclip
if (plotsum EQ 1) then begin
    charlegendarray = strarr(index)
    charlegendarray[0:index-2] = tg1[where(min_arr ne 0)]
    charlegendarray[index-1]   = abstotchar
endif else begin
    charlegendarray = strarr(index)
    charlegendarray[0:index-1] = tg1[where(min_arr ne 0)]
endelse
;;  aii_klegend_d,tg1[where(min_arr ne 0)],/fill,box=0,$
aii_klegend_d, charlegendarray[0:index-1], $
               /fill, $
               box=0,$
               usersym=usersym, $
               charsize=1.1*!P.charsize, $
               spacing=3.2,$
               psym=intarr(index), $
               colors=colors[0:index-1],$
               line=ls[0:index-1], $
               position=[0.6,0.5+index*0.03], $
               thick=2*thick

!P.multi=0

;; close plot output file
aii_plot_file, action='end', show='yes', print='no', $
               outdir='/smiles_local/continuum/borysow'

;; restore settings
P_ini = !P    

end
