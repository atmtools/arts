; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
;
; **************************************************************************
;+
; NAME:
;        check_backslash
;
; Purpose:  checks if the specified path name has a backslash at the
;           end or not. If  not the backslash will be added.
;
; Inputs:   pathname      path name
;
; Output:   pathname+'/'  path name with backslash at the end
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;-
; **************************************************************************
; 
FUNCTION check_backslash, pathname
;
; ---- SET PATH NAME-------------------------------------------
pathname = STRCOMPRESS(pathname, /REMOVE_ALL) ; remove whitespace
if (STRMID(pathname, STRLEN(pathname)-1, 1) NE '/') THEN pathname=pathname+'/'
;
RETURN, pathname
;
END
;
; **************************************************************************
; Name:     aii_check_unit
;
; Purpose:  checks if a unit is open for writing
;
; Inputs:   unit          unit
;
; Output:   flag          0=ok, 1=error occured
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
PRO aii_check_unit, unit, flag
;
; ---- GET INFORMATION ABOUT UNIT ------------------
A = FSTAT(unit)
;
; ---- CHECK UNIT ----------------------------------
flag = 1
IF (A.UNIT GT 0) THEN flag=0
;
END
;
; **************************************************************************
; Name:     check_output_unit
;
; Purpose:  checks if a unit is open for writing
;
; Inputs:   unit          unit
;
; Output:   flag          0=ok, 1=error occured
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION check_output_unit, unit
;
; ---- GET INFORMATION ABOUT UNIT ------------------
A = FSTAT(unit)
;
; ---- CHECK UNIT ----------------------------------
IF (A.UNIT GT 0) THEN RETURN, 0
RETURN, 1
;
END
;
; **************************************************************************
; Name:     check_input_file
;
; Purpose:  checks if the specified file exists in a specified subdirectory
;
; Inputs:   filename      name of the file to be read
;           pathname      path under which the file exists
;
; Output:   filename      vector of all the files detected
;
; History:  2001-01-04    Thomas Kuhn, iup Bremen
;
; **************************************************************************
; 
FUNCTION check_input_file, filename, pathname
;
; ---- SET PATH NAME-------------------------------------------
pathname = STRCOMPRESS(pathname, /REMOVE_ALL) ; remove whitespace
pathname = check_backslash(pathname)
;
; ---- CHECK FILE ---------------------------------------------
filename = STRCOMPRESS(filename, /REMOVE_ALL) ; remove whitespace
f = pathname+filename
filevec = FINDFILE(f,COUNT=n)
;
RETURN, filevec
;
END
;
; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
