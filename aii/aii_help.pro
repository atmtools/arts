;; ==========================================================================
;; ####################### ARTS IDL INTERFACE PROCEDURE #####################
;; ==========================================================================
;;
;+
;NAME:
;          aii_help
;
;KEYWORDS: AIIDir         arts/aii/ directory path
;          HelpFileName   file name of the html help page
;          TITLE          title of mthe html page
;          
;
;
;PURPOSE : this procedure sets the frame for a help page in HTML format
;          for all the aii functions/procedures
;
;CALLING : aii_help, AIIDir='~/username/arts/aii/', $
;                    HelpFileName='aii_help_file.html'
;                    TITLE = 'Arts IDL Interface Documentation (aii doc)'
;
;HISTORY : 2003-09-08, Christian Melsheimer, iup Bremen     
;          added automatic creation of HTML documentation of aii routines. 
;
;          2003-11-14, Thomas Kuhn, iup Bremen     
;          put into a separate procedure with additional keywords
;           
;-
;;
;; #########################################################################


PRO aii_help, AIIDir=AIIDir, HelpFileName=HelpFileName, TITLE=TITLE

; generating HTML help file with documentation of all aii files

if not keyword_set(AIIDir) then AIIDir='~/arts/aii/'
AIIDir = STRCOMPRESS(STRTRIM(AIIDir, 2))

;; check if arts/aii directory is correctly given
Result = FINDFILE(AIIDir+'*.pro', count=k)
if (k lt 1) then begin
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print,'!! in procedure aii_help:                                  '
    print,'!! ERROR, input keyword AIIDir not appropriate!            '
    print,'!! input arts/aii directory path name: >>'+AIIDir+'<<' 
    print,'!! please check the path of your arts/aii directory again. '
    print,'!! termiate without creating the html file!                '
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    goto, ende
endif


if not keyword_set(HelpFileName) then HelpFileName='aii_help.html'
HelpFileName = STRCOMPRESS(STRTRIM(HelpFileName, 2))


if not keyword_set(TITLE) then $
  TITLE = 'Arts IDL Interface Documentation (aii doc)'


print,'aii_help> Generating HTML documentation of AII procedures/functions'

mk_html_help, AIIDir, AIIDir+HelpFileName, TITLE=TITLE

print,'aii_help> Documentation written to file '+AIIDir+HelpFileName


ende:
close, /all
END
