;******************************************************************************
;+
;NAME:
;           read_artsvar
;PURPOSE:
;           Reads a ARTS variable.
;
;           The data is read from the file
;              basename.varname.am
;
;           See further read_datafile.
;
; Format:   x = read_artsvar(basename, varname [, /check])
;
;
; Inputs:   basename    the ARTS basename 
;           varname     variable name
; Optional: check       Keyword to check the data
;
; Output:   x           the data
;
; History:  28.02.01  Wolfram Haas
;-
;******************************************************************************

FUNCTION read_artsvar, basename, varname, check = check

; Create full file name
name = basename + '.' + varname + '.aa'

; Read the data by using read_datafile
IF keyword_set(check) THEN RETURN, read_datafile(name, /check) $
                      ELSE RETURN, read_datafile(name)

END
