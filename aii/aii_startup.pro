;; ==========================================================================
;; ####################### ARTS IDL INTERFACE PROCEDURE #####################
;; ==========================================================================
;+
;NAME:
;           aii_startup
;
;PURPOSE : this batch file sets the frame for an ARTS IDL session
;
;CALLING : in the IDL session type '@aii_startup'
;
;HISTORY : alpha version 2003-04-03, Thomas Kuhn, iup Bremen

;          added automatic creation of HTML documentation of aii
;          routines,  2003-09-08, Christian Melsheimer, iup Bremen     
;           
;-
;;
;; #########################################################################

print, ' '
print,'====================< begin of AII_STARTUP >===================='
print, ' '


; OS system specification to 8-bit mode for unix systems
;if ( !version.os_family EQ 'unix' ) then device, pseudo_color=8
;if ( !version.os_family EQ 'unix' ) then device, true_color=24

; set display mode of window
;window, /free, /pixmap, colors=-10 ;
;wdelete, !D.window


; --- device = X window ---------------------------------------------------------
;SET_PLOT, 'X'
;device, retain=2, $                  ; manage window backing store
;        decomposed=1, $              ; only affects 24-bit mode sets decomposed color off
;        set_character_size=[10,12],$ ; vector font size
;        true_color=24
;XPAGE=20.9  &  YPAGE=29.7  &  XOFFS=0.0  &  YOFFS=0.0,
;X0=1.374  &  Y0=1.283  &  XLEN=3.622  &  YLEN=6.157 ;Inches
;!P.POSITION=[X0/XPAGE,Y0/YPAGE,(X0+XLEN)/XPAGE,(Y0+YLEN)/YPAGE]


; --- device = Postscript file --------------------------------------------------
set_plot, 'PS'  ; setting graphics device (= WIN,MAC,X,PS,PRINTER,METAFILE,Z,...)
DEVICE,/PORTRAIT
DEVICE, /COLOR
DEVICE, /helvetica


; --- path ----------------------------------------------------------------------
spawn ,'echo $HOME',userhome
!PATH=!PATH+':'+userhome                                $
           +':'+userhome+'/pro'                         $
           +':'+'/pool/lookup/idl'                      $
           +':'+'/pool/lookup2/idl/TeXtoIDL'            $
           +':'+'/pool/lookup2/idl/astron/pro'          $
           +':'+'/pool/lookup2/idl/idlps'               $
           +':'+'/pool/lookup2/idl/idlps/support'

; increase size of history buffer
!EDIT_INPUT=500                    ; increase size of history buffer


; --- plot ----------------------------------------------------------------------
; set plotting system variables:
;!P.noerase                      ; screen/page will not be erased
!P.font       = 0                ; 0:device fonts, -1:vector fonts
!P.thick      = 2.5
!P.color      = 0                ; generel black lines/points
!P.background = 255              ; generel white background
!P.charthick  = 1.5              ;
!P.charsize   = 1.0
!X.charsize   = 1.0
!Y.charsize   = 1.0
!X.omargin    = [0,0]
!Y.omargin    = [0,0]
xstyle        = 2
ystyle        = 2

print, ' system options/parameters/variables setting:'
print, ' !D.name                   =',!D.name
print, ' !D.x_size,!D.y_size       =',!D.x_size,' ',!D.y_size
print, ' !D.x_vsize,!D.y_vsize     =',!D.x_vsize,!D.y_vsize
print, ' !D.x_ch_size,!D.y_ch_size =',!D.x_ch_size,' ',!D.y_ch_size
;
print, ' !P.font                   =',!P.font
print, ' !P.thick                  =',!P.thick
print, ' !P.color                  =',!P.color
print, ' !P.charsize               =',!P.charsize
print, ' !P.background             =',!P.background
print, ' !P.thick                  =',!P.thick
print, ' !P.charthick              =',!P.charthick
print, ' !P.charsize               =',!P.charsize
print, ' !X.charsize               =',!X.charsize
print, ' !Y.charsize               =',!Y.charsize

print, ' !PATH                     =',!PATH
print, ' '

; generating HTML help file with documentation of all aii files
helpfile='aii_help.html'
print,'Generating HTML documentation of AII procedures/functions'
mk_html_help, '.', helpfile, TITLE = "Arts IDL Interface Documentation"
print,'Documentation written to ' + helpfile
;
print,''
print,'======================< end of AII_STARTUP >======================'
print, ' '

;; #########################################################################

;; compiling all aii procedures/function:
@aii_compile

;; #########################################################################

