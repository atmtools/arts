; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
;
;
; ==========================================================================
;
FUNCTION datinfo, plotfilename
;
;; get datum to write it on the top of the plot:
spawn,'date +"%y"',year
spawn,'date +"%m"',month
spawn,'date +"%d"',day
spawn,'date +"%H"',hour
spawn,'date +"%M"',minute
spawn,'whoami',who
datum = who+':20'+string(year, FORMAT='(A2)')+'-'+$
string(month, FORMAT='(A2)')+'-'+string(day, FORMAT='(A2)')
;
if (STRLEN(plotfilename) GE 1) then datum = datum+'/'+plotfilename
;
datum = STRCOMPRESS(datum, /REMOVE_ALL)
datum = STRTRIM(datum, 2)
;
RETURN, datum
END
;
; ==========================================================================
;
FUNCTION PlotSetupA, plotfilename, plotfileformat
;
;; use aii_plot_file for writing into plot output file
aii_plot_file, action='begin', fname=plotfilename, fformat=plotfileformat
;
!P.MULTI     = 0
!P.MULTI     = [0,1,1]
!P.FONT      = 1
!P.CHARSIZE  = 1.5
!X.CHARSIZE  = 1.25
!Y.CHARSIZE  = 1.25
!P.THICK     = 8
!X.THICK     = 5
!Y.THICK     = 5
!P.CHARTHICK = 4
!X.MARGIN    = [0, 0]
!Y.MARGIN    = [0, 0]
!P.POSITION  = [0.2, 0.2, 0.8, 0.9]
print,'PlotSetupA: plot range: ',!P.POSITION
;
;
RETURN, !P.POSITION
END
;
; ==========================================================================
;
FUNCTION PlotSetupB, a, b
;
ok = 1
aii_plot_file, action='end', show=a, print='no', $
               outdir=b
ok = 0
;
RETURN, ok
END
;
;
; ==========================================================================
;
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
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
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
;
PRO plot_trans_1, jobname, f, abs, alt, $
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
                       altunit=altunit, $
                       plotfilename=plotfilename,$
                       plotfileformat=plotfileformat, $
                       plotyaxis=plotyaxis
;
;======================================================
;
;; save settings
P_ini = !P


;; set directory name apropriate
if keyword_set(jobdir) then begin 
    jobdirname = jobdir+jobname
    outdir     = jobdir
endif else begin
    jobdirname = jobname
    spawn, 'pwd', outdir
endelse


;; additional plotting of total absorption? 
YSCALE = ' '
if keyword_set(plotyaxis) then begin
    if (STRUPCASE(plotyaxis) EQ 'LOG') THEN YSCALE='/YLOG, '
endif

;; read the tag groups from controlfile
read_tag_groups,jobdir+jobname+'.arts','tgsDefine',tg


;; total number of tag groups
tot_tg= (size(tg,/DIMENSIONS))[0]


;; read the transmission per line of sight step width
if not keyword_set(read) then begin
    print,'Reading transmission data...'
    aa = jobdirname+'.trans.aa'
    trans = aa_read_general(aa)
endif
nmats = N_ELEMENTS(trans[*,0,0]) ; # of angles of line of sight (los)
nrows = N_ELEMENTS(trans[0,*,0]) ; # of frequencies
ncols = N_ELEMENTS(trans[0,0,*]) ; # of steps along los
print, 'trans: nmats=',nmats,', rows=',nrows,', colums=',ncols 
;; check that the tag groups and the absoprtion per tag group have the
;; same dimension

;; at the moment only one angle can be plotted due to the fact that
;; 'aa_read_general' can only handle several matrices if the single
;; matrices have the same size
IF (nmats GT 1) THEN BEGIN
    PRINT,'PLOT_TRANS_1> SORRY ONLY TRANSMISSION MATRICES WITH SAME SIZES CAN BE READ'
    PRINT,'              THEREFORE ONLY CALCULATIONS FOR A SINGLE ANGLE CAN BE PLOTTED'
    PRINT,'              SORRY - BUT WE HAVE TO TERMINATE NOW SINCE  # of angles=',nmats
    GOTO, ENDE
ENDIF

;; calculate total transmission
tottrans = dblarr(nmats,nrows)
tottrans[*,*] = 1.000e0
for i = 0,nmats-1 do begin
    for j = 0,nrows-1 do begin
        for k = 0,ncols-1 do begin
            tottrans[i,j] = tottrans[i,j] * trans[i,j,k]
        endfor
    endfor
endfor


;; 6a) now read viewing angle data 
print,'Reading viewing angle data...'
vangles = aa_read_general(jobdirname+'.za_pencil.aa')
if (nmats ne N_ELEMENTS(vangles)) then begin
    print,'Error: viewing angle vector and number of transmission matrices'
    print,'       have not the same dimension'
    stop
endif


;; 6b) now read frequency data 
if not keyword_set(read) then begin
    print,'Reading frequency data...'
    f   = aa_read_general(jobdirname+'.f_mono.aa')
endif
if (nrows ne N_ELEMENTS(f[*])) then begin
    print,'Error: absorption per tag has not the same numbers of'
    print,'       frequencies like the frequency vector'
    stop
endif


;; now read altitude data 
if not keyword_set(read) then begin
print,'Reading altitude data...'
alt = aa_read_general(jobdirname+'.z_abs.aa')
endif


;; reading pressure information
print,'Reading pressure grid data...'
p = aa_read_general(jobdirname+'.p_abs.aa')


;; reading temperature information
print,'Reading temperature grid data...'
T = aa_read_general(jobdirname+'.t_abs.aa')


;; frequency unit selection
if keyword_set(cm) then begin
    unit=2.9979E10
    unitstr='[cm!u-1!n]'
endif else begin
    unit=1E9
    unitstr='[GHz]'
endelse
for i = 0,N_ELEMENTS(f)-1 do f[i] = f[i] / unit



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
endif else begin
    pscale = 1.000/100.000
    punit = 'hPa'
endelse

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



;; altitude unit setting
if keyword_set(altunit) then begin
    case altunit of
        'km'  : begin
                 zscale = 0.001
                 zunit = 'km'
                 end
        'm'   : begin
                 zscale = 1.000
                 zunit = 'm'
                 end
        else:    begin
                 zscale = 1.000
                 zunit = 'm'
                 end
    endcase
endif else begin
    zscale = 0.001000
    zunit = 'km'
endelse
for i = 0, N_ELEMENTS(alt)-1 do begin
    alt[i] = alt[i] * zscale
endfor



;; altitude range selection
if keyword_set(yrange) then begin
    yunit = 'ALTITUDE'
    print,'yrange=',yrange
    if (N_ELEMENTS(yrange) LT 1) then begin
        print,'ERROR in finding y-axis range!!!'
        print,'too few elements of vector: yrange=',yrange
        print,'terminate now'
        goto, ende
    endif
    if (yunit EQ 'ALTITUDE') then begin
        if (N_ELEMENTS(yrange) EQ 1) then begin
            yrange0 = 0.00
            yrange1 = zscale*yrange*1.000e3
        endif
        if (N_ELEMENTS(yrange) EQ 2) then begin
            yrange0 = zscale*yrange[0]*1.000e3
            yrange1 = zscale*yrange[1]*1.000e3
        endif
        delta0 = 10000.00
        delta1 = 10000.00
        iyrange0 = 0
        iyrange1 = N_ELEMENTS(alt)-1
        for i = 0,N_ELEMENTS(alt)-1 do begin
            if (ABS(alt[i]-yrange0) LT delta0) then begin
                delta0   = ABS(alt[i]-yrange0)
                iyrange0 = i
            endif
            if (ABS(alt[i]-yrange1) LT delta1) then begin
                delta1   = ABS(alt[i]-yrange1)
                iyrange1 = i
            endif
        endfor
    endif
    if (yunit EQ 'PRESSURE') then begin
        if (N_ELEMENTS(yrange) EQ 1) then begin
            yrange0 = pscale*0.00
            yrange1 = pscale*yrange
        endif
        if (N_ELEMENTS(yrange) EQ 2) then begin
            yrange0 = pscale*yrange[0]
            yrange1 = pscale*yrange[1]
        endif
        delta0 = 10000.00
        delta1 = 10000.00
        iyrange0 = 0
        iyrange1 = N_ELEMENTS(p)-1
        for i = 0,N_ELEMENTS(p)-1 do begin
            if (ABS(p[i]-yrange0) LT delta0) then begin
                delta0   = ABS(p[i]-yrange0)
                iyrange0 = i
            endif
            if (ABS(p[i]-yrange1) LT delta1) then begin
                delta1   = ABS(p[i]-yrange1)
                iyrange1 = i
            endif
        endfor
    endif
endif else begin
    yunit = 'ALTITUDE'
    iyrange0 = 0
    iyrange1 = N_ELEMENTS(p)-1
    if (yunit EQ 'ALTITUDE') then begin
        iyrange0 = 0
        iyrange1 = N_ELEMENTS(alt)-1
    endif
    if (yunit EQ 'PRESSURE') then begin
        iyrange0 = 0
        iyrange1 = N_ELEMENTS(p)-1
    endif
endelse
print,' y range selected:'
print,'   pressure:  p_min=',p[iyrange0],punit,' p_max=',p[iyrange1],punit
print,'   altitude:  z_min=',alt[iyrange0],zunit,' z_max=',alt[iyrange1],zunit
;
;
;
;                     -----------------------------------
; =================== plot the transmission vs. frequency ===================
;                     -----------------------------------


;; check for color output? -> always use colors
if not keyword_set(color) then color=1


;; create color and lineshape arrays
colors=intarr(nrows+1)
ls=intarr(nrows+1)
for i=0,nrows do begin
    colors[i] = i mod 29    ;; we have 29 different colors in mycolor
    ls[i]     = i mod  5    ;; we have  5 different line styles
endfor


;; add to title string
total_title='Transmission'
if keyword_set(add_to_title) then total_title=', '+add_to_title



;; apropriate x range of plot
if not keyword_set(xrange) then begin
  xmin = min(f)
  xmax = max(f)
endif else begin
    xmin = (xrange[0]/ unit)
    xmax = (xrange[1]/ unit)
endelse


;; run over  all viewing angles
for iangle = 0,nmats-1 do begin


; open PS file
;; use aii_plot_file for writing into plot output file
    if not keyword_set(plotfilename)   then $
       plotfilename=jobname+'_plot_'+STRCOMPRESS(string(iangle, FORMAT='(I3)'))
    if not keyword_set(plotfileformat) then plotfileformat=4
    pos = PlotSetupA(plotfilename, plotfileformat)
    
    
;;  transmission axis range    
    ymin = MIN(tottrans[iangle,0:nrows-1], MAX=ymax)
    if ( (STRPOS(YSCALE, 'LOG') GE 0) AND (ymin LE 0.000e0) )then ymin=1.00e-20
    print,'[x],[y] plot range: [',xmin,',',xmax,'], [',ymin,',',ymax,']'


;; set frame of the plot
    xtitle='frequency '+unitstr
    ytitle='transmission [1]'
    RESULT= EXECUTE('plot, '+YSCALE+'f[0:nrows-1], '+             $
                    'tottrans[iangle,0:nrows-1], '+               $
                    '/NORMAL, '+                                  $
                    'title='+total_title+', '+                    $
                    'xrange=[xmin, xmax], '+                      $
                    'yrange=[ymin, ymax], '+                      $
                    'xtitle=xtitle, '+                            $
                    'ytitle=ytitle, '+                            $
                    'color=colors[0], '+                          $
                    'linestyle=ls[0], '+                          $
                    'xstyle=2, '+                                 $
                    'ystyle=1, '+                                 $
                    '/NODATA')
    print,'EXECUTE RESULT= ',RESULT
    col_start = 2
    oplot,                            $
      f[0:nrows-1],                   $
      tottrans[iangle,0:nrows-1],     $
      color=colors[col_start+iangle], $
      linestyle=ls[iangle]
;
;;  close PS file    
    ok = PlotSetupB('yes', outdir)
endfor
;    
;; restore settings
P_ini = !P
;
; --- end of procedure -----------------------------------------------------
ende:
;
END
; ==========================================================================
; ##########################################################################
; ==========================================================================
