;******************************************************************************
; Name:     read_artsvar
;
;           Reads a ARTS variable.
;
;           The data is read from the file
;              basename.varname.am
;
;           See further read_datafile.
;
; Format:   x = read_artsvar(basename, varname [, /optimize])
;
;
; Inputs:   basename    the ARTS basename 
;           varname     variable name
; Optional: optimize    Keyword to accelerate reading in
;
; Output:   x           the data
;
; History:  13.11.00  Wolfram Haas
;******************************************************************************

FUNCTION read_artsvar, basename, varname, optimize = optimize

; Create full file name
name = basename + '.' + varname + '.am'

; Read the data by using read_datafile
RETURN, read_datafile(name, /optimize)

END
