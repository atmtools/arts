FUNCTION AII_GET_HOSTNAME
;+
;NAME:
;      aii_get_hostname
;PURPOSE:
;      returns the hostname (without the domain part)
;CALLING:
;      result = aii_get_hostname()
;-
SPAWN,'echo $HOSTNAME',HOSTNAME
RESULT = STR_SEP( HOSTNAME[0], '.' )

RETURN, RESULT[0]
END
