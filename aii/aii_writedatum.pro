; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================^X
;
FUNCTION aii_writedatum, TEXT=TEXT
;+
;NAME: 
;     aii_writedatum
;Purpose:
;     get date to write it on the top of the plot:
;     EXAMPLE OUTPUT OF SYSTIME: Wed Feb 19 13:39:49 2003
;Usage:
;     result = aii_writedatum()
;Arguments:
;     None
;Keywords:
;     TEXT : text to put before the date
;-
DATUM = SYSTIME()
STRVEC = STR_SEP(strcompress(DATUM), ' ')

W=GETENV('LOGNAME')

A = W+'/'+STRVEC[4]+'-'+STRVEC[1]+'-'+STRVEC[2]+'/'+STRVEC[3]

IF KEYWORD_SET(TEXT) THEN BEGIN 
 A = TEXT+A
ENDIF

RETURN, A
END
;
; ==========================================================================

