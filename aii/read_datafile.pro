;******************************************************************************
; Name:     read_datafile
;
;           Reads data from a file in ARTS data format.
;
;           The data is returned as an array or a structure of arrays.
;           For example, to get matrix 2, type
;              m = x.mat1
;
;
; Format:   x = read_datafile(filename [, /check])
;
; Inputs:   filename      full file name
; Optional: check         Keyword to check the data
;
; Output:   x             the data  
;
; History:  28.02.00  Wolfram Haas
;******************************************************************************

FUNCTION read_datafile, filename, check = check

; Open file for reading
openr, unit, filename, error = err, /get_lun

IF err NE 0 THEN printf, -2, !err_string

; Read until line does not begin with #
s = '#'
WHILE strmid(s, 0, 1) EQ '#' DO readf, unit, s

s = strcompress(strtrim(s, 2)) & p = str_sep(s, ' ')

IF s EQ '' THEN BEGIN
  print, 'Blank lines are not allowed.'
  stop
ENDIF

IF n_elements(p) GT 1 THEN BEGIN
  print, 'Missing number of matrices.'
  stop
ENDIF

; Read number of matrices
reads, s, nmat
IF nmat LT 1 THEN BEGIN
  print, 'Could not read number of matrices.'
  stop
ENDIF

pa = ptrarr(nmat)  ; returns a pointer array

nrow = intarr(nmat) & ncol = intarr(nmat)

;; loop over all matrices
FOR i = 0, nmat - 1 DO BEGIN
  ; Read size
  s = '#'
  WHILE strmid(s, 0, 1) EQ '#' DO $
    IF NOT eof(unit) THEN readf, unit, s ELSE BEGIN
      print, 'Wrong number of matrices.'
      stop
    ENDELSE

  ; Check size
  s = strcompress(strtrim(s, 2)) & p = str_sep(s, ' ')

  IF s EQ '' THEN BEGIN
    print, 'Blank lines are not allowed.'
    stop
  ENDIF

  IF n_elements(p) NE 2 THEN BEGIN
    print, 'Could not read matrix size.'
    stop
  ENDIF

  reads, s, j, k
  nrow(i) = j & ncol(i) = k

  IF keyword_set(check) THEN BEGIN
    sa = strarr(k, j)

    ; Read the matrix row by row
    FOR r = 0, j - 1 DO BEGIN
      s = '#'
      WHILE strmid(s, 0, 1) EQ '#' DO $
        IF NOT eof(unit) THEN readf, unit, s ELSE BEGIN
          print, 'One or more rows are missing.'
          stop
        ENDELSE

      s = strcompress(strtrim(s, 2)) & p = str_sep(s, ' ')

      IF s EQ '' THEN BEGIN
        print, 'Blank lines are not allowed within a matrix.'
        stop
      ENDIF

      IF n_elements(p) NE k THEN BEGIN
        print, 'Wrong number of column elements.'
        stop
      ENDIF

      sa(*, r) = p
    ENDFOR

    mat = dblarr(k, j)

    mat(*, *) = sa(*, *)
  ENDIF ELSE BEGIN
    mat = dblarr(k, j)
    readf, unit, mat
  ENDELSE
 
  pa(i) = ptr_new(mat)                     ; returns a pointer to pa(i)
ENDFOR

IF nmat EQ 1 THEN mat = *pa(0) ELSE BEGIN  ; create structure
  s1 = 'mat = {' & s2 = ''
  FOR i = 0, nmat - 1 DO BEGIN
    s1 = s1 + 'mat' + string(i, format = '(I0)') $
            + ': dblarr(' + string(ncol(i), format = '(I0)') + ', ' $
                          + string(nrow(i), format = '(I0)') + '), '
    s2 = s2 + 'mat.mat' + string(i, format = '(I0)') $
            + ' = *pa(' + string(i, format = '(I0)') + ') & '
  ENDFOR
  s1 = strmid(s1, 0, strlen(s1) - 2) + '}'
  s2 = strmid(s2, 0, strlen(s2) - 3)

  r1 = execute(s1)
  r2 = execute(s2)
ENDELSE

; Release memory used by the heap variable
ptr_free, pa

; Check if the file is now finished
s = ''
WHILE NOT eof(unit) DO BEGIN
  readf, unit, s

  s = strcompress(s, /remove_all)

  IF s NE '' THEN BEGIN
    print, 'There is some garbage at the end of the file.'
    stop
  ENDIF
ENDWHILE

; Close the file
free_lun, unit  
RETURN, mat

END
