; ==========================================================================
; ####################### ARTS IDL INTERFACE PROCEDURE #####################
; ==========================================================================
;
;###########################################################################
; Name:     aa_read_general
;
; Purpose:  Reads data from a file ARTS data format file into a matrix.
;
; Inputs:   filename      full file name
;
; Output:   matrix        the data matrix
;
; History:  2001-12-04    Thomas Kuhn, iup Bremen
;
;***************************************************************************
;
; 
FUNCTION aa_read_general, filename
;=================================
;
; ---- CHECKS ---------------------------------------------
filevec = FINDFILE(filename)
IF (N_ELEMENTS(filevec) GT 1) THEN BEGIN
    print,'aa_read_general> !!! ERROR: can not find data file!'
    print,'aa_read_general> !!! Multiple files found with name: ',filename,':',filevec,'! STOP!'
    dummy = 1
    return, dummy
ENDIF
IF (N_ELEMENTS(filevec) LT 1) THEN BEGIN
    print,'aa_read_general!!! ERROR: can not find data file!'
    print,'aa_read_general!!! No files found with name: ',filename,'! STOP!'
    dummy = 1
    return, dummy
ENDIF
;
; ---- OPEN FILE FOR READING ------------------------------
openr, unit, filename, error = err, /GET_LUN
;
IF err NE 0 THEN BEGIN
    print,' aa_read_general> !!! ERROR in procedure aa_read_general occured,'
    print,' aa_read_general> !!! could not open arts file >>'+filename+'<<'
    print,' aa_read_general> !!! please check the directory and file names.'
    print,' aa_read_general> !!! error message:', !ERR_STRING
    print,' aa_read_general> !!! STOP here!'
    dummy = 1
    return, dummy
ENDIF
;
; ---- CREATE MATRIX --------------------------------------
;matrix = dblarr(2, 2, 2)
g_nrows = ULONG(0)
g_ncols = ULONG(0)
;
; ---- READ NUMBER OF MATRICES ----------------------------
; Read until line does not begin with #
s = '#'
WHILE strmid(s, 0, 1) EQ '#' DO BEGIN 
    readf, unit, s
ENDWHILE
s = STRCOMPRESS(STRTRIM(s, 2))
;print,'s=',s

IF s EQ '' THEN BEGIN
    print, 'aa_read_general> Input line:',s
    print, 'aa_read_general> Blank lines are not allowed! STOP!'
    close, unit
    FREE_LUN, unit
    dummy = 1
    return, dummy
ENDIF

p = STR_SEP(s, ' ')

IF n_elements(p) GT 1 THEN BEGIN
    print,'aa_read_general> String array:',p
    print,'aa_read_general> Missing number of matrices! STOP! file: '+filename
    close, unit
    FREE_LUN, unit
    dummy = 1
    return, dummy
ENDIF

nmats = ULONG(0)
reads, s, nmats
IF nmats LT 1 THEN BEGIN
    print,'aa_read_general> no. of matrices:',nmats
    print, '-> Could not read number of matrices! STOP! file: '+filename
    close, unit
    FREE_LUN, unit
    dummy = 1
    return, dummy
ENDIF
;
; ---- READ EACH MATRIX -----------------------------------
FOR i = 0L, nmats - 1L DO BEGIN ; loop over each matrix
    ; a) ignore comment lines
    s = '#'
    WHILE strmid(s, 0, 1) EQ '#' DO $
      IF NOT eof(unit) THEN readf, unit, s ELSE BEGIN
        print, 'aa_read_general> One or more rows are missing! STOP! file: '+filename
        close, unit
        FREE_LUN, unit
        dummy = 1
        return, dummy
    ENDELSE
    
    ; b) read matrix size
    s = strcompress(strtrim(s, 2))
    IF s EQ '' THEN BEGIN
        print, 'aa_read_general> Input line:',s
        print, 'aa_read_general> Blank lines in matrices are not allowed! STOP! file: '+filename
        close, unit
        FREE_LUN, unit
        dummy = 1
        return, dummy
    ENDIF
;        
    p = str_sep(s, ' ')
    IF n_elements(p) NE 2 THEN BEGIN
        print,'aa_read_general> String array: rows=',p(0),'colums=',p(1)
        print,'aa_read_general> Could not read matrix size in loop! STOP!. file: '+filename
        close, unit
        FREE_LUN, unit
        dummy = 1
        return, dummy
    ENDIF
    ; b.1) read rows and colums info 
    nrows = ULONG(0)
    ncols = ULONG(0)
    reads, p[0], nrows
    reads, p[1], ncols

    ; c) define output matrix once
    IF (i EQ 0) THEN BEGIN
        matrix = dblarr(nmats, nrows, ncols)
        g_nrows = ULONG(nrows)
        g_ncols = ULONG(ncols)
    ENDIF

    ; d) check matric size each new matrix starts
    IF ((nrows NE g_nrows) OR (ncols NE g_ncols)) THEN BEGIN
        print,'aa_read_general> The number of rows/colums has changed from one matrix to'
        print,'aa_read_general> the next in one file! STOP! file: '+filename
        close, unit
        FREE_LUN, unit
        dummy = 1
        return, dummy
    ENDIF
;
    ; e) Read the matrix row by row
    FOR r = 0L, nrows - 1L DO BEGIN
        s = '#'
        WHILE strmid(s, 0, 1) EQ '#' DO $
          IF NOT eof(unit) THEN readf, unit, s ELSE BEGIN
            print, 'aa_read_general> One or more rows are missing! STOP! file: '+filename
            close, unit
            FREE_LUN, unit
            dummy = 1
            return, dummy
        ENDELSE

        s = strcompress(strtrim(s, 2))
        IF s EQ '' THEN BEGIN
            print, 'aa_read_general> Blank lines are not allowed within a matrix! STOP! file: '+filename
            close, unit
            FREE_LUN, unit
            dummy = 1
            return, dummy
        ENDIF

        p = str_sep(s, ' ')
        ;print,'p=',p
        IF N_ELEMENTS(p) NE ncols THEN BEGIN
            print, $
              'aa_read_general> Wrong number of column elements in matrix! STOP! file: '+filename
            print, 'aa_read_general> ',N_ELEMENTS(p), ncols
            close, unit
            FREE_LUN, unit
            dummy = 1
            return, dummy
        ENDIF
        FOR k = 0L, ncols-1L DO BEGIN
            dummy = 0.0e0
            reads, p[k], dummy
            matrix[i, r, k] = dummy
        ENDFOR
    ENDFOR 
;
ENDFOR
;
; ---- CHECK IF THE FILE IS NOW FINISHED ------------------
s = ''
WHILE NOT eof(unit) DO BEGIN
  readf, unit, s

  s = strcompress(s, /remove_all)

  IF s NE '' THEN BEGIN
    print, 'aa_read_general> WARNING: There is some additional lines at the end of the file: '+filename
    print, '"',s,'"'
    print, 'aa_read_general> WARNING: this part is not in the matrix'
;    close, unit
;    FREE_LUN, unit
;    return, matrix
  ENDIF

ENDWHILE
;
; ---- FREE_LUN THE FILE ----------------------------------
close, unit
FREE_LUN, unit  
;
; ---- RETURN MATRIX --------------------------------------
RETURN, matrix
;
END
;
; ==========================================================================
; ##########################################################################
; ==========================================================================
