;; ==========================================================================
;; ####################### ARTS IDL INTERFACE PROCEDURE #####################
;; ==========================================================================
;;
;###########################################################################
;; Name:     aa_read_general
;;
;; Purpose:  Reads data from a file ARTS data format file into a matrix.
;;
;; Inputs:   filename      full file name
;;
;; Output:   matrix        the data matrix
;;
;; History:  2001-12-04    Thomas Kuhn, iup Bremen
;;           2003-08-04    AvE added verbose option
;;
;***************************************************************************

 
;;=============================================
FUNCTION aa_read_general, filename, $
                          COMMENTSYM=COMMENTSYM,$
                          verbose=verbose
;;=============================================


;; ---- CHECKS ---------------------------------------------
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


;; ---- OPEN FILE FOR READING ------------------------------

OPENR, unit, filename, error=err, /GET_LUN
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


;; ---- READ COMMENTS ----------------------------------------

;; DEFINE THE COMMENT LINE SYMBOL
IF NOT KEYWORD_SET(COMMENTSYM) THEN COMMENTSYM = '#'
COMMENTSYM = STRCOMPRESS( COMMENTSYM, /REMOVE_ALL )

s = '#'
WHILE (strmid(s, 0, 1) EQ COMMENTSYM) DO BEGIN 
    readf, unit, s
ENDWHILE

s = STRCOMPRESS(STRTRIM(s, 2)) ;; REMOVE BLANKS FROM THE STRING

IF s EQ '' THEN BEGIN
    print, 'aa_read_general> Blank lines are not allowed! STOP!'
    close, unit
    FREE_LUN, unit
    dummy = 1
    return, dummy
ENDIF


;; ---- CREATE MATRIX BY READING THE DIMENSIONS -------------

;; A) NUMBER OF MATRICES
nmats = ULONG(0)
reads, s, nmats
IF nmats LT 1 THEN BEGIN
    PRINT,'aa_read_general> no. of matrices:',nmats
    PRINT, '-> Could not read number of matrices! STOP! file: '+filename
    CLOSE, unit
    FREE_LUN, unit
    dummy = 1
    return, dummy
ENDIF


;; B) NUMBER OF ROWS AND COLUMNS
readf, unit, s
s = STRCOMPRESS(STRTRIM(s, 2))
IF s EQ '' THEN BEGIN
    PRINT, 'aa_read_general> Blank lines in matrices are not allowed! STOP! file: '+filename
    CLOSE, unit
    FREE_LUN, unit
    dummy = 1
    return, dummy
ENDIF

p = STRSPLIT(s, ' ', /EXTRACT)
IF N_ELEMENTS(p) NE 2 THEN BEGIN
    print,'aa_read_general> Row/column string array:',p
    print,'aa_read_general> Could not read matrix size in loop! STOP!. file: '+filename
    close, unit
    FREE_LUN, unit
    dummy = 1
    return, dummy
ENDIF
nrows = ULONG(0)
ncols = ULONG(0)
reads, p[0], nrows
reads, p[1], ncols

;; C) DEFINE OUTPUT MATRIX ONCE
matrix = DBLARR(nmats, nrows, ncols)
if keyword_set(verbose) then print,'defined matrix: nmats=',nmats,', nrows=',nrows,', ncols=',ncols


;; ---- READ EACH MATRIX -----------------------------------
FOR i = 0L, nmats - 1L DO BEGIN ;; loop over each matrix

    ;; B) READ THE MATRIX ROW BY ROW
    FOR r = 0L, nrows-1L DO BEGIN

        readf, unit, s ;; READ INPUT LINE

        s = strcompress(strtrim(s, 2)) ;; REMOVE BLANKS

        p = STRSPLIT(s, ' ', /EXTRACT)
        IF N_ELEMENTS(p) NE ncols THEN BEGIN
            print, $
              'aa_read_general> Wrong number of column elements in matrix! STOP! file: '+filename
            print, 'aa_read_general> N_ELEMENTS(p)=',N_ELEMENTS(p),' <-> ncols=',ncols
            print, 'aa_read_general> input string=', p
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

    
;;  C) CHECK EVERY MATRIX IF THE SIZE IS CONSISTENT
    IF (i LT (nmats-1L)) THEN BEGIN

        readf, unit, s ;; READ INPUT LINE

        ;; NUMBER OF ROWS AND COLUMNS:
        s = strcompress(strtrim(s, 2))
        p = STRSPLIT(s, ' ', /EXTRACT)
        IF n_elements(p) NE 2 THEN BEGIN
            print,'aa_read_general> !!ERROR!! in file: '+filename
            print,'aa_read_general> Matrix',i,'has not the appropriate rows/columns!'
            print,'aa_read_general> String array: rows=',p(0),'colums=',p(1)
            print,'aa_read_general> TERMINATE NOW'
            close, unit
            FREE_LUN, unit
            dummy = 1
            return, dummy
        ENDIF
        orows = ULONG(0)
        ocols = ULONG(0)
        reads, p[0], orows
        reads, p[1], ocols
;;      check matric size each new matrix starts
        IF ((nrows NE orows) OR (ncols NE ocols)) THEN BEGIN
            print,'aa_read_general> !!ERROR!! in file: '+filename
            print,'aa_read_general> Matrix',i,'has not the appropriate size!'
            print,'aa_read_general> The matrx is defined as:'
            print,'aa_read_general>  -> no. of matrices:',nmats
            print,'aa_read_general>  -> no. of rows    :',nrows
            print,'aa_read_general>  -> no. of columns :',ncols
            print,'aa_read_general> TERMINATE NOW'
            close, unit
            FREE_LUN, unit
            dummy = 1
            return, dummy
        ENDIF
    ENDIF

ENDFOR




;; ---- CHECK IF THE FILE IS NOW FINISHED ------------------
s = ''
WHILE NOT eof(unit) DO BEGIN
  readf, unit, s

  s = strcompress(s, /remove_all)

  IF s NE '' THEN BEGIN
    print, 'aa_read_general> WARNING: There is some additional lines at the end of the file: '+filename
    print, '"',s,'"'
    print, 'aa_read_general> WARNING: this part is not in the matrix'
  ENDIF

ENDWHILE


;; ---- FREE_LUN THE FILE ----------------------------------
close, unit
FREE_LUN, unit  


;; ---- RETURN MATRIX --------------------------------------
RETURN, matrix

END

;; ==========================================================================
;; ##########################################################################
;; ==========================================================================
