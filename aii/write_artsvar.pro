;******************************************************************************
; Name:     write_artsvar
;
;           Writes a ARTS variable to a file in ARTS format.
;           The data is written to a file called
;              basename.varname.am
;           See further write_datafile.
;
; Format:   write_artsvar, basename, varname, x [, prec]
;
; Inputs:   basename    the ARTS basename
;           varname     variable name
;           x           the data to store
; Optional: prec        number of digits to use, default 6
;                       If prec = 0, integer values are assumed.
;
; Output:   -
;
; History:  13.11.00  Wolfram Haas
;******************************************************************************

PRO write_artsvar, basename, varname, x, prec

; Create full file name
name = basename + '.' + varname + '.am'

; Create heading text
heading = 'This file contains the ARTS variable ' + varname + '.'

IF n_params() EQ 3 THEN write_datafile, name, x, heading       $
                   ELSE write_datafile, name, x, heading, prec

END
