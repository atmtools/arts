; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
pro aii_plot_file, action=action, $
                   fname=fname, $
                   fformat=fformat, $
                   show=show, $
                   print=print, $
                   outdir=outdir
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
COMMON PLOT_OUTPUT_FUN_ANTWORT, ANTWORT, USERFILENAME
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
; default values for input parameters:
if not keyword_set(outdir) then begin
    outdir = '.'
endif
;
if not keyword_set(print) then begin
    print = 'no'
endif
;
if not keyword_set(show) then begin
    show = 'no'
endif
;
; check input parameter action' of correctness:
if not keyword_set(action) then begin
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print,' PROCEDURE: AII_PLOT_FILE.pro                         '
    print,' ERROR:     no value for input variable action        '
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    retall
endif
;
if (action NE 'begin') and (action NE 'end') then begin
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print,' PROCEDURE  : AII_PLOT_FILE.pro                       '
    print,' ERROR      : wrong value for variable action         '
    print,' ALLOWED    : begin , end                             '
    print,' YOUR CHOICE:',action
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    retall
endif
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
;
; set constants
; *************
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
; ==========================================================================
;
; open output file
; ****************
;
if (action EQ 'begin') then begin
;
;   check input variable for output file name:
    if (keyword_set(fname)) then begin
        USERFILENAME = fname ; store file name in common block variable for later use
    endif else begin
        USERFILENAME = 'aii_plot_file' ; default = name of this procedure
        print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print,' PROCEDURE: AII_PLOT_FILE.pro  '
        print,' ATTENTION: output file name is'
        print,' >> ',USERFILENAME,' <<'
        print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    endelse
;
;   check input variable for output file format:
    if (keyword_set(fformat)) then begin
        ANTWORT = fformat
    endif else begin
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
        if ((ANTWORT LT 1) or (ANTWORT GT 5)) then begin
            ANTWORT=1
        endif
    endelse
;
   if (ANTWORT NE 1) AND (ANTWORT NE 2) AND $
      (ANTWORT NE 3) AND (ANTWORT NE 4) AND $
      (ANTWORT NE 5) then begin
       print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print,'!! Procedure: AII_PLOT_FILE.pro                  !!'
       print,'!! ATTENTION: wrong value for output file format !!'
       print,'!! FFORMAT  :',ANTWORT,'                                 !!',$
             ANTWORT
       print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       retall
   endif
;
   case ANTWORT of
        1 : begin
              ; Postscript portrait mode
              dummy         = !P.color ; make background white
              !P.color      = !P.background
              !P.background = dummy
              aii_color_table
              spawn,'rm -f '+outdir+'/'+USERFILENAME+EXTENSION_VEC[ANTWORT]
              SET_PLOT, 'PS'
              DEVICE, FILENAME=USERFILENAME+EXTENSION_VEC[ANTWORT]
              DEVICE,/PORTRAIT
              DEVICE, /COLOR
              end
        2 : begin
              ; Postscript landscape mode
              dummy         = !P.color ; make background white
              !P.color      = !P.background
              !P.background = dummy
              spawn,'rm -f '+outdir+'/'+USERFILENAME+EXTENSION_VEC[ANTWORT]
              SET_PLOT, 'PS'
              DEVICE, FILENAME=USERFILENAME+EXTENSION_VEC[ANTWORT]
              DEVICE, /LANDSCAPE
              DEVICE, /COLOR
              aii_color_table
              end
        3 : begin
              ; encapsulated Postscript portrait mode
              dummy         = !P.color ; make background white
              !P.color      = !P.background
              !P.background = dummy
              color=2
              spawn,'rm -f '+outdir+'/'+USERFILENAME+EXTENSION_VEC[ANTWORT]
              SET_PLOT, 'PS'
              DEVICE, FILENAME=USERFILENAME+EXTENSION_VEC[ANTWORT]
              DEVICE, /ENCAPSULATED
              DEVICE, /COLOR
              DEVICE, /PORTRAIT
              aii_color_table
              end
        4 : begin
              ; encapsulated Postscript landscape mode
              dummy         = !P.color ; make background white
              !P.color      = !P.background
              !P.background = dummy
              color=2
              spawn,'rm -f '+outdir+'/'+USERFILENAME+EXTENSION_VEC[ANTWORT]
              SET_PLOT, 'PS'
              DEVICE, FILENAME=USERFILENAME+EXTENSION_VEC[ANTWORT]
              DEVICE, /ENCAPSULATED
              DEVICE, /COLOR
              DEVICE, /LANDSCAPE
              aii_color_table
              end
        5 : begin
              ; window
              SET_PLOT, 'X'
              end
        else: begin
              dummy         = !P.color ; make background white
              !P.color      = !P.background
              !P.background = dummy
              SET_PLOT, 'X'
              XPAGE=20.9  &  YPAGE=29.7  &  XOFFS=0.0  &  YOFFS=0.0
              X0=1.374  &  Y0=1.283  &  XLEN=3.622  &  YLEN=6.157  ;Inches
              !P.POSITION=[X0/XPAGE,Y0/YPAGE,(X0+XLEN)/XPAGE,(Y0+YLEN)/YPAGE]
              aii_color_table
              end
   endcase
endif
;
;
; ==========================================================================
;
;
; close and print output file
; ***************************
;
if (action EQ 'end') then begin
;
; a) close output file:
; ---------------------
    DEVICE,/CLOSE
;
; b) print info:
; --------------
    print, '-------------------------- aii_plot_file --------------------------'
    print, ' * dir   : ','"'+outdir+'"'
    print, ' * file  : ','"'+USERFILENAME+EXTENSION_VEC[ANTWORT]+'"'
    print, ' * format: ',FORMAT_VEC[ANTWORT]
;
; c) move output file into specified directory:
; ---------------------------------------------
    spawn, 'mv '+USERFILENAME+EXTENSION_VEC[ANTWORT]+' '+outdir+'/.'
;
; d) close device and open gostview:
; ----------------------------------
    if (show EQ 'yes') then begin
        print, ' * show with ghostview: yes'
        spawn,'gv -swap -a4 -bg white -fg black '+GVOPT[ANTWORT]+' '+$
              outdir+'/'+USERFILENAME+EXTENSION_VEC[ANTWORT]+' &'
    endif else begin
        print, ' * show with ghostview: no'
    endelse
;
; e) printing:
; ------------
    if (print EQ 'yes') then begin
        print, ' * print with lpr: yes'
        spawn,'lpr '+outdir+'/'+USERFILENAME+EXTENSION_VEC[ANTWORT]
    endif else begin
        print, ' * print with lpr: no'
    endelse
;
; f) print info:
; --------------
    print, '-------------------------- aii_plot_file --------------------------'
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
