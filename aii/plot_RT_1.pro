; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
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
;
; ==========================================================================
;
FUNCTION plotsetup0, text
;
;; use aii_plot_file for writing into plot output file
aii_plot_file, action='begin', fname=text, fformat=4
;
;; settings for the plot
!P.MULTI     = 0
!P.MULTI     = [0,1,1]
!P.FONT      = 1
!P.CHARSIZE  = 1.5
!X.CHARSIZE  = 1
!Y.CHARSIZE  = 1
!P.THICK     = 5
!X.THICK     = 5
!Y.THICK     = 5
!X.CHARSIZE  = 1.5
!Y.CHARSIZE  = 1.5
!P.CHARSIZE  = 1.5
!P.CHARTHICK = 4
!P.POSITION  = [0.2, 0.2, 0.8, 0.9]
;
;
RETURN, !P.POSITION
END
;
; ==========================================================================
;
FUNCTION plotsetup1, a, b
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
PRO plot_RT_1, jobname, $
               add_to_title=add_to_title,$
               color=color, $
               xrange=xrange, $
               yrange=yrange,$
               nostamp=nostamp, $
               jobdir=jobdir, $
               plotfilename=plotfilename,$
               plotfileformat=plotfileformat, $
               plotyaxis=plotyaxis
;======================================================
;
IF (N_ELEMENTS(jobname[*]) LT 1) THEN BEGIN
    print,'plot_RT_1: no input arts job names, terminate!'
    goto, ende
ENDIF


jobdirname = strarr(N_ELEMENTS(jobname[*])) 
;; set core file names apropriate
FOR i = 0, N_ELEMENTS(jobname[*])-1 DO BEGIN
    if keyword_set(jobdir) then begin 
        jobdirname[i] = jobdir+jobname[i]
        outdir     = jobdir
    endif else begin
        jobdirname[i] = jobname[i]
        outdir     = '/home/home01/tkuhn/ARTS'
    endelse
ENDFOR


;; read frequency data 
if not keyword_set(read) then begin
    print,'Reading frequency data...'
    f   = aa_read_general(jobdirname[0]+'.f_mono.aa')
endif


;; now read pencil beam zenit angles 
if not keyword_set(read) then begin
print,'Reading pencil beam angle data...'
pencils = aa_read_general(jobdirname[0]+'.za_pencil.aa')
endif


;; now reads the RT temperatures
RTT = dblarr(N_ELEMENTS(f[*]),       $
             N_ELEMENTS(pencils[*]), $
             N_ELEMENTS(jobname[*]))
FOR i = 0, N_ELEMENTS(jobname[*])-1 DO BEGIN
    aa = jobdirname[i]+'.y.aa'
    print,'Reading absorption data file:',aa
    RT = aa_read_general(aa)
    IF (N_ELEMENTS(RT[*]) NE $
        (N_ELEMENTS(f[*])*N_ELEMENTS(pencils[*])) ) THEN BEGIN
        print, 'plot_RT_1: inconsistent size of RT matrix with f and pencil matrices!'
        print, ' N_ELEMENTS(RT[*])     =',N_ELEMENTS(RT[*])
        print, ' N_ELEMENTS(f[*])      =',N_ELEMENTS(f[*])
        print, ' N_ELEMENTS(pencils[*])=',N_ELEMENTS(pencils[*])
        print, ' ==>terminate now!'
        goto, ende
    ENDIF
    FOR l = 0,N_ELEMENTS(pencils[*])-1 DO BEGIN
        RTT[0:N_ELEMENTS(f[*])-1, l, i] = RT[(l*(N_ELEMENTS(f[*])-1)):((l+1)*(N_ELEMENTS(f[*])-1))]
;        print,l,': RTT=',RTT[0:N_ELEMENTS(f[*])-1, l, i]
    ENDFOR
ENDFOR


;; frequency unit is [GHz]
funit=1.0e-9
funitstr='[GHz]'
for i = 0,N_ELEMENTS(f)-1 do f[i] = f[i] * funit



;                       ---------------------------------
; =================== plot the absorption vs. frequency ===================
;                       ---------------------------------

;; save settings
P_ini = !P


;; make 1 plot per page
!P.multi    = [0,1,1]


;; settings for the plot
!X.MARGIN   = [0,0]
!Y.MARGIN   = [0,0]
!X.OMARGIN  = [0,0]
!Y.OMARGIN  = [0,0]
!P.REGION   = [0.0, 0.0, 0.0, 0.0]
!P.POSITION = [0.025, 0.1, 0.6, 0.9]
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


;; create color and lineshape arrays
colors = intarr(N_ELEMENTS(jobname[*]))
ls     = intarr(N_ELEMENTS(jobname[*]))
for i=0, N_ELEMENTS(jobname[*])-1 do begin
    colors[i] = i mod 29    ;; we have 29 different colors in mycolor
    ls[i]     = i mod  5    ;; we have  5 different line styles
endfor


;; add to title string
if keyword_set(add_to_title) then total_title = total_title+add_to_title


;; apropriate x range of plot
if not keyword_set(xrange) then begin
  xmin = min(f)
  xmax = max(f)
endif else begin
    xmin = xrange[0]
    xmax = xrange[1]
endelse


imax = MIN([N_ELEMENTS(pencils[*])-1, 5])
FOR i = 0, imax DO BEGIN ; -----------------------------loop over all angles
;; use aii_plot_file for writing into plot output file
    ok = plotsetup0('RTvsf_plot_'+string(i,FORMAT='(I1)'))

    print, i,': pencils[i]=',pencils[i]
    pz = string(pencils[i],FORMAT='(F6.1)')
    total_title = TeXtoIDL(' pencil angle='+pz+'^{o}', font=0)


;; set plot style 'log' or 'lin' for y-axis
    if NOT keyword_set(plotyaxis) then plotyaxis='lin'
    if ((plotyaxis NE 'lin') AND (plotyaxis NE 'log')) then plotyaxis='lin'
    
    
    ymin = MIN([RTT[0:N_ELEMENTS(f[*])-1,i,0], RTT[0:N_ELEMENTS(f[*])-1,i,1]])
    ymax = MAX([RTT[0:N_ELEMENTS(f[*])-1,i,0], RTT[0:N_ELEMENTS(f[*])-1,i,1]])
    
;; set frame of the plot
    IF (plotyaxis EQ 'lin') then begin
        plot, f[0:N_ELEMENTS(f[*])-1] ,        $
          RTT[0:N_ELEMENTS(f[*])-1,i,0],       $
          /NORMAL,                             $
          title=total_title,                   $
;          xcharsize=1.25,                      $
;          ycharsize=1.25,                      $
          xrange=[xmin, xmax],                 $
          yrange=[ymin, ymax],                 $
          xtitle=TeXtoIDL('frequency '+funitstr, font=0), $
          ytitle=TeXtoIDL('T_B [K]', font=0),  $
          color=colors[0],                     $
          linestyle=ls[0],                     $
          xstyle=2,                            $
          ystyle=2,                            $
          /nodata
    endif
    IF (plotyaxis EQ 'log') then begin
;; set frame of the plot
        plot, /YLOG, $ 
          f[0:N_ELEMENTS(f[*])-1] ,            $
          RTT[0:N_ELEMENTS(f[*])-1,i,0],       $
          /NORMAL,                             $
          title=total_title,                   $
;          xcharsize=1.25,                      $
;          ycharsize=1.25,                      $
          xrange=[xmin, xmax],                 $
          yrange=[ymin, ymax],                 $
          xtitle=TeXtoIDL('frequency '+funitstr, font=0), $
          ytitle=TeXtoIDL('T_B [K]', font=0),  $
          color=colors[0],                     $
          linestyle=ls[0],                     $
          xstyle=2,                            $
          ystyle=1,                            $
          /nodata
    endif
    
    
    FOR j = 0,N_ELEMENTS(jobname[*])-1 DO BEGIN
        oplot, f[0:N_ELEMENTS(f[*])-1] ,      $
          RTT[0:N_ELEMENTS(f[*])-1,i,j],      $
          color=colors[j],                    $
          linestyle=ls[j]
    ENDFOR
    
    
; print datum and user name:
    xyouts, ok[0], ok[3]+0.005, datum, CHARSIZE=0.75, CHARTHICK=1.0, /NORMAL
    
    
;; legend of the plot
    charlegendarray = strarr(N_ELEMENTS(jobname[*]))
    ch_length       = N_ELEMENTS(charlegendarray)
    for k = 0,N_ELEMENTS(jobname[*])-1 do charlegendarray[k] = jobname[k]
    aii_plot_legend, charlegendarray[0:N_ELEMENTS(jobname[*])-1], $
      box=0,$
      usersym=usersym, $
      spacing=2.0,$
      pspacing=2.0, $
      psym=intarr(ch_length), $
      colors=colors[0:ch_length-1],$
      line=ls[0:ch_length-1], $
      position=[ok[2]+0.005, ok[3]], thick=thick*1.4
    
    
    ok = plotsetup1('yes', outdir)
    
ENDFOR


;; restore settings
P_ini = !P    
;
; --- end of procedure -----------------------------------------------------
ende:
close, /all
;
END
; ==========================================================================
; ##########################################################################
; ==========================================================================
