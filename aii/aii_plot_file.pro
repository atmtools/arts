; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;;
FUNCTION PSPlotOpen, text, n
;;
;; ---------------------------------------------------------------------
;; NAME    : PSPlotOpen
;; PURPOSE : opens a postscript file for IDL plots
;; EXTERNAL: calls procedure aii_plot_file
;; INPUT   : text  STRING   postscript file name (without extension) 
;;           n     INTEGER  postscript format 
;;                          1: postscript portrait
;;                          2: postscript landscape
;;                          3: encapsulated postscript portrait
;;                          4: encapsulated postscript landscape
;; OUTPUT  : ok    INTEGER  flag if everything went well (0=ok, 1=error)
;; ---------------------------------------------------------------------
;;
;; flag: 0=ok, 1=error occured
ok = 1
;
;; use aii_plot_file for writing into plot output file
aii_plot_file, action='begin', fname=text, fformat=n
;
ok = 0
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
;;
;; ============================================================================
;;
FUNCTION PSPlotClose, a, b, c
;;
;; ------------------------------------------------------------------------
;; NAME    : PSPlotClose
;; PURPOSE : close a postscript file for IDL plots
;; EXTERNAL: calls procedure aii_plot_file
;; INPUT   : a  STRING  open postscript file with ghostview a='yes' or 'no'
;;           b  STRING  send postscript file to printer     b='yes' or 'no'
;;           c  STRING  move postscript file to directory   c="DIRNAME"
;;           d  STRING  write datum lower left corner       d='yes' or 'no'
;; OUTPUT  : ok INTEGER flag if everything went well (0=ok, 1=error)
;; ------------------------------------------------------------------------
;;
ok = 1  ;; flag: 0=ok, 1=error occured
aii_plot_file, action='end', $
               show=a,       $
               print=b,      $
               outdir=c,     $
               writedate=d
ok = 0
;
RETURN, ok
END

;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


; ==========================================================================
;
pro aii_plot_file, action=action,      $
                   fname=fname,        $
                   fformat=fformat,    $
                   show=show,          $
                   print=print,        $
                   outdir=outdir,      $
                   writedate=writedate
;
;==========================================================================
;
; NAME:
;       aii_plot_file
; PURPOSE:
;       provide user specified color output plot file for ARTS jobs
; EXPLANATION:
;       aii_plot_file has tow parst. The first part is for opening an output
;       device. The user interactively decide which device should be
;       opened. The second part close the output device and send the 
;       output file to the printer if wanted.
;       Therefore the user has to call this procedure before the first
;       plot command is performed in the user procedure with
;       action='begin' and the at the very end of the user procedure
;       again with action='end'.
;
; CALLING EXAMPLES:
;         first call of aii_plot_file to open a Postscript output file
;         aii_plot_file, action='begin', fname='myfile', fformat=1
;
;         At the end of the user procedure call aii_plot_file again
;         to close the previously opend output file
;         aii_plot_file, action='end'
;
; INPUTS:
;       action  (string)   the content of the string decides the action, 
;                          if an output device should be opened or
;                          closed. Possible actions are 'begin' and 'end'.
;
;       fname   (string)   output file name without extension
;       fformat   (integer)  possible output files formats are
;                          1: Postscript portrait mode, 
;                          2: Postscript landscape mode, 
;                          3: encapsulated Postscript portrait mode,
;                          4: encapsulated Postscript landscape mode,
;                          5: window
;       show    (string)   parameter which handles the visualization
;                          of the output file. For Postscript or 
;                          encapsulated Postscript files it opens a
;                          gostview.
;                          Possible values are 'yes' and 'no'.
;       print   (string)   variable which handles the printing on the
;                          standard printer of the output file.
;                          Possible values are 'yes' and 'no'.
;       writedate (string) state if the user name and date should be
;                          written into the plot.  
;                          Possible values are 'yes' and 'no'.
;
; OUTPUTS:
;       no own output, but handling of graphical device output.
;
;
; MODIFICATION HISTORY:
;       03/04/01  TKS  alpha version created 
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
; define common block:
COMMON AII_PLOT_OUTPUT_CONTROL, ANTWORT, USERFILENAME, DEVICENAME, PUSERINITIAL
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
;
;
;                             ----------------
; ===========================< check keywords >=============================
;                             ----------------
;
; default values for input parameters:
IF NOT KEYWORD_SET(outdir) THEN BEGIN
    outdir = ''
ENDIF ELSE BEGIN
    outdir = check_backslash(outdir)
ENDELSE
;
IF NOT KEYWORD_SET(print) THEN BEGIN
    print = 'no' ; default value
ENDIF
;
IF NOT KEYWORD_SET(show) THEN BEGIN
    show = 'no'  ; default value
ENDIF
;
IF NOT KEYWORD_SET(writedate) THEN BEGIN
    writedate = 'NO' ; default value
ENDIF
;
; check input parameter action' of correctness:
IF NOT KEYWORD_SET(action) THEN BEGIN
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print,' PROCEDURE: AII_PLOT_FILE.pro                         '
    print,' ERROR:     no value for input variable action set    '
    print,'            RETURN without action                     '
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    RETURN
ENDIF
;
IF (STRLOWCASE(action) NE 'begin') AND (STRLOWCASE(action) NE 'end') THEN BEGIN
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print,' PROCEDURE  : AII_PLOT_FILE.pro                       '
    print,' ERROR      : wrong value for variable action         '
    print,' ALLOWED    : begin , end                             '
    print,' YOUR CHOICE:',action
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    RETURN
ENDIF
;
;
;                             ---------------
; ===========================< set constants >==============================
;                             ---------------
;
;
; possible extensions of the output file:
EXTENSION_VEC = ['.XX', $
                 '.ps', $
                 '.ps', $
                 '.eps',$
                 '.eps',$
                 '.XX']
;
FORMAT_VEC    = [' NOT DEFINED ',                          $
                 'POSTSCRIPT PORTRAIT MODE',               $
                 'POSTSCRIPT LANDSCAPE MODE',              $
                 'ENCAPSULATED POSTSCRIPT PORTRAIT MODE',  $
                 'ENCAPSULATED POSTSCRIPT LANDSCAPE MODE', $
                 'WINDOW']
;
GVOPT         = ['',           $
                 '-portrait',  $
                 '-landscape', $
                 '-portrait',  $
                 '-landscape', $
                 ''] 
;
;
;
;                           ------------------
; =========================< open output file >=============================
;                           ------------------
;
;
if (STRLOWCASE(action) EQ 'begin') then begin
;
;   a) save original user settings
;   ------------------------------
    PUSERINITIAL = !P
;
;   b) check input variable for output file name
;   --------------------------------------------
    IF (keyword_set(fname)) THEN BEGIN
        USERFILENAME = fname ; store file name in common block variable for later use
    ENDIF ELSE BEGIN
        USERFILENAME = 'aii_plot_file' ; default = name of this procedure
        print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print,' PROCEDURE: AII_PLOT_FILE.pro  '
        print,' ATTENTION: output file name is'
        print,' >> ',USERFILENAME,' <<'
        print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    ENDELSE
;
;   c) ask for input variable for output file format
;   ------------------------------------------------
    IF (keyword_set(fformat)) THEN BEGIN
        ANTWORT = fformat
    ENDIF ELSE BEGIN
        ANTWORT = 1  ; default = Postscript portrait mode
        print,'_________( aii_plot_file interactive mode )_________'
        print,'| select output file format  (default=1):          |'
        print,'| ==> (1) Postscript  portrait mode                |'
        print,'|     (2) Postscript  landscape mode               |'
        print,'|     (3) Encapsulated Postscript portrait mode    |'
        print,'|     (4) Encapsulated Postscript landscape mode   |'
        print,'|     (5) window                                   |'
        read, ANTWORT
        print,'|________( aii_plot_file interactive mode )________|'
        IF ((ANTWORT LT 1) OR (ANTWORT GT 5)) THEN BEGIN
            ANTWORT=1
        ENDIF
    ENDELSE
;
;   d) check input variable for output file format
;   ----------------------------------------------
    IF ( (ANTWORT NE 1) AND (ANTWORT NE 2) AND $
         (ANTWORT NE 3) AND (ANTWORT NE 4) AND $
         (ANTWORT NE 5) ) THEN BEGIN
       print, FORMAT='(A52)', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print, FORMAT='(A52)', '!! Procedure: AII_PLOT_FILE.pro                  !!'
       print, FORMAT='(A52)', '!! ATTENTION: wrong value for output file format !!'
       print, FORMAT='(A13,I1,A38)', '!! FFORMAT  :',ANTWORT,'                                    !!'
       print, FORMAT='(A52)', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       RETURN
   endif
;
;  e) save name of current graphics device
;  ---------------------------------------
   DEVICENAME = !D.name
;
;  f) remove files with the same output name
;  -----------------------------------------
   USERFILENAME = USERFILENAME+EXTENSION_VEC[ANTWORT]
   xf = check_input_file(USERFILENAME, outdir)
   IF (N_ELEMENTS(xf) GT 0) THEN $
     FOR i = 0, N_ELEMENTS(xf)-1 DO spawn,'rm -f '+xf[i]
;
;  g) open the output file according to input specifications
;  ---------------------------------------------------------
   case ANTWORT of
        1 : begin
              ; --- Postscript portrait mode ------------------------
;              dummy         = !P.color ; make background white
;              !P.color      = !P.background
;              !P.background = dummy
              SET_PLOT, 'PS'
              DEVICE, FILENAME=USERFILENAME, /times, /PORTRAIT, $
                      XSIZE=21.0, YSIZE=29.6, XOFFS=0.0,  YOFFS=0.0, /COLOR
;              DEVICE,/PORTRAIT
;              DEVICE, /COLOR
;              DEVICE, SET_FONT='Times'
              aii_color_table
              end
        2 : begin
              ; --- Postscript landscape mode -----------------------
;              dummy         = !P.color ; make background white
;              !P.color      = !P.background
;              !P.background = dummy
              SET_PLOT, 'PS'
              DEVICE, FILENAME=USERFILENAME, /times, /LANDSCAPE, $
                      YSIZE=21.0, XSIZE=29.6, XOFFS=0.0,  YOFFS=0.0, /COLOR  
;              DEVICE, /LANDSCAPE
;              DEVICE, /COLOR
;              DEVICE, SET_FONT='Times'
              aii_color_table
              end
        3 : begin
              ; --- encapsulated Postscript portrait mode -----------
;              dummy         = !P.color ; make background white
;              !P.color      = !P.background
;              !P.background = dummy
              SET_PLOT, 'PS'
              DEVICE, FILENAME=USERFILENAME, /times, /PORTRAIT, /ENCAPSULATED, /COLOR
;              DEVICE, /ENCAPSULATED
;              DEVICE, /COLOR
;              DEVICE, /PORTRAIT
;              DEVICE, SET_FONT='Times'
              aii_color_table
              end
        4 : begin
              ; --- encapsulated Postscript landscape mode ----------
;              dummy         = !P.color ; make background white
;              !P.color      = !P.background
;              !P.background = dummy
              SET_PLOT, 'PS'
              DEVICE, FILENAME=USERFILENAME, /times, /LANDSCAPE, /ENCAPSULATED, /COLOR
;              DEVICE, /ENCAPSULATED
;              DEVICE, /COLOR
;              DEVICE, /LANDSCAPE
;              DEVICE, SET_FONT='Times'
              aii_color_table
              end
        5 : begin
              ; --- window ------------------------------------------
              dummy         = !P.color ; make background white
              !P.color      = !P.background
              !P.background = dummy
              SET_PLOT, 'X'
              !P.FONT = 0
              device, retain=2,     $   ; manage window backing store
                      decomposed=1, $   ; only affects 24-bit mode sets decomposed color off
                      set_character_size=[10,12],$ ; vector font size
                      true_color=24     ; set true color mode for display mode
;                      pseudo_color=8
;              XPAGE=20.9  &  YPAGE=29.7  &  XOFFS=0.0  &  YOFFS=0.0
;              X0=1.374  &  Y0=1.283  &  XLEN=3.622  &  YLEN=6.157 ;Inches
;              !P.POSITION=[X0/XPAGE,Y0/YPAGE,(X0+XLEN)/XPAGE,(Y0+YLEN)/YPAGE]
              aii_color_table
              end
        else: begin
              ; --- default -----------------------------------------
              dummy         = !P.color ; make background white
              !P.color      = !P.background
              !P.background = dummy
              SET_PLOT, 'X'
              XPAGE=20.9  &  YPAGE=29.7  &  XOFFS=0.0  &  YOFFS=0.0
              X0=1.374  &  Y0=1.283  &  XLEN=3.622  &  YLEN=6.157  ;Inches
              !P.POSITION=[X0/XPAGE,Y0/YPAGE,(X0+XLEN)/XPAGE,(Y0+YLEN)/YPAGE]
              end
   endcase
endif
;
;
;
;                      -----------------------------
; ====================< close and print output file >=======================
;                      -----------------------------
;
;
if (action EQ 'end') then begin
;
; a) write date and user name into the lower left corner of the plot
; ------------------------------------------------------------------
    IF (STRUPCASE(writedate) EQ 'YES') THEN BEGIN
        XYOUTS, 0.1, 0.05, $
                aii_writedatum(), $
                CHARSIZE=0.45, $
                ALIGNMENT=0.0, $
                /NORMAL
    ENDIF

; b) close output file
; --------------------
    IF ((ANTWORT GE 1) AND (ANTWORT LE 4)) THEN DEVICE, /close_file
    SET_PLOT, DEVICENAME
;
; c) print info
; -------------
    print, ' * aii_plot_file> dir   : ','>>'+outdir+'<<'
    print, ' * aii_plot_file> file  : ','>>'+USERFILENAME+'<<'
    print, ' * aii_plot_file> format: ',FORMAT_VEC[ANTWORT]
;
; d) move output file into specified directory
; --------------------------------------------
    IF ((ANTWORT GE 1) AND (ANTWORT LE 4)) THEN BEGIN
        spawn, 'mv -u '+USERFILENAME+' '+outdir
    ENDIF
;
; e) close device and open ghostview
; ----------------------------------
    IF ((ANTWORT GE 1) AND (ANTWORT LE 4)) THEN BEGIN
        if (show EQ 'yes') then begin
            print, ' * show with ghostview: yes'
            spawn,'gv -swap -a4 -bg white -fg black '+GVOPT[ANTWORT]+' '+$
                  outdir+USERFILENAME+' &'
        endif else begin
            print, ' * show with ghostview: no'
        endelse
    ENDIF
;
; f) printing
; -----------
    IF ((ANTWORT GE 1) AND (ANTWORT LE 4)) THEN BEGIN
        if (print EQ 'yes') then begin
            print, ' * print with lpr: yes'
            spawn,'lpr '+outdir+USERFILENAME
        endif else begin
            print, ' * print with lpr: no'
        endelse
    ENDIF
;
; g) set back to saved original user settings
; -------------------------------------------
    !P = PUSERINITIAL
;
endif
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
; end of procedure aii_plot_file
end
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
