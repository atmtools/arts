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
    stop
ENDIF
IF (N_ELEMENTS(filevec) LT 1) THEN BEGIN
    print,'aa_read_general!!! ERROR: can not find data file!'
    print,'aa_read_general!!! No files found with name: ',filename,'! STOP!'
    stop
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
    FREE_LUN, unit & STOP
ENDIF
;
; ---- CREATE MATRIX --------------------------------------
;matrix = dblarr(2, 2, 2)
g_nrows = 0
g_ncols = 0
;
; ---- READ NUMBER OF MATRICES ----------------------------
; Read until line does not begin with #
s = '#'
WHILE strmid(s, 0, 1) EQ '#' DO BEGIN 
    readf, unit, s
ENDWHILE
s = strcompress(strtrim(s, 2)) & p = str_sep(s, ' ')

IF s EQ '' THEN BEGIN
    print, 'aa_read_general> Input line:',s
    print, 'aa_read_general> Blank lines are not allowed! STOP!'
    FREE_LUN, unit & stop
ENDIF

IF n_elements(p) GT 1 THEN BEGIN
    print,'aa_read_general> String array:',p
    print,'aa_read_general> Missing number of matrices! STOP!'
    FREE_LUN, unit & stop
ENDIF

reads, s, nmats
IF nmats LT 1 THEN BEGIN
    print,'aa_read_general> no. of matrices:',nmats
    print, '-> Could not read number of matrices! STOP!'
    FREE_LUN, unit & stop
ENDIF
;
; ---- READ EACH MATRIX -----------------------------------
FOR i = 0, nmats - 1 DO BEGIN ; loop over each matrix
    ; a) ignore comment lines
    s = '#'
    WHILE strmid(s, 0, 1) EQ '#' DO $
      IF NOT eof(unit) THEN readf, unit, s ELSE BEGIN
        print, 'aa_read_general> One or more rows are missing! STOP!'
        FREE_LUN, unit & stop
    ENDELSE
    
    ; b) read matrix size
    s = strcompress(strtrim(s, 2))
    IF s EQ '' THEN BEGIN
        print, 'aa_read_general> Input line:',s
        print, 'aa_read_general> Blank lines in matrices are not allowed! STOP!'
        FREE_LUN, unit & stop
    ENDIF
;        
    p = str_sep(s, ' ')
    IF n_elements(p) NE 2 THEN BEGIN
        print,'aa_read_general> String array: rows=',p(0),'colums=',p(1)
        print,'aa_read_general> Could not read matrix size in loop! STOP!.'
        FREE_LUN, unit & stop
    ENDIF
    ; b.1) read rows and colums info 
    reads, p[0], nrows
    reads, p[1], ncols

    ; c) define output matrix once
    IF (i EQ 0) THEN BEGIN
        matrix = dblarr(nmats, nrows, ncols)
        g_nrows = nrows
        g_ncols = ncols
    ENDIF

    ; d) check matric size each new matrix starts
    IF ((nrows NE g_nrows) OR (ncols NE g_ncols)) THEN BEGIN
        print,'aa_read_general> The number of rows/colums has changed from one matrix to'
        print,'aa_read_general> the next in one file! STOP!'
        FREE_LUN, unit & STOP
    ENDIF
;
    ; e) Read the matrix row by row
    FOR r = 0, nrows - 1 DO BEGIN
        s = '#'
        WHILE strmid(s, 0, 1) EQ '#' DO $
          IF NOT eof(unit) THEN readf, unit, s ELSE BEGIN
            print, 'aa_read_general> One or more rows are missing! STOP!'
            FREE_LUN, unit & stop
        ENDELSE

        s = strcompress(strtrim(s, 2))
        IF s EQ '' THEN BEGIN
            print, 'aa_read_general> Blank lines are not allowed within a matrix! STOP!'
            FREE_LUN, unit & stop
        ENDIF

        p = str_sep(s, ' ')
        ;print,'p=',p
        IF N_ELEMENTS(p) NE ncols THEN BEGIN
            print, 'aa_read_general> Wrong number of column elements in matrix! STOP!'
            FREE_LUN, unit & stop
        ENDIF
        FOR k = 0, ncols-1 DO BEGIN
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
    print, 'aa_read_general> There is some garbage at the end of the file:'
    print, '"',s,'"'
    FREE_LUN, unit & stop
  ENDIF
ENDWHILE
;
; ---- FREE_LUN THE FILE ----------------------------------
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
