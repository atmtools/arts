;******************************************************************************
;+
;NAME:
;           write_datafile
;Purpose:
;           Writes data to a file in ARTS format.
;
;           See also read_datafile.
;
; Format:   write_datafile(filename, x, heading [, prec])
;
; Inputs:   filename    full file name
;           x           the data to store
;           heading     heading text
;                       The function puts in '# ' at the start of each line. 
;                       Heading can be empty ('')
; Optional: prec        number of decimals to use, default 6
;                       If prec = 0, integer values are assumed.
;                             
; Output:   -
;
; History:  28.02.01  Wolfram Haas
;- 
;******************************************************************************

PRO write_mat, unit, x, prec

v = size(x) & ndim = v(0) & ncol = v(1) & nrow = v(2)

IF ndim EQ 0 OR ndim EQ 1 THEN BEGIN      ; x is a number or a line vector
  ncol = n_elements(x) & nrow = 1
ENDIF

sa = strarr(ncol, nrow)

IF prec GT 0 THEN $                       ; floating-point values
  FOR i = 0, nrow - 1 DO $
    FOR j = 0, ncol - 1 DO BEGIN
      s = 'sa(' + string(j, format = '(I0)') + ', ' $
                + string(i, format = '(I0)') + ') ' $
          + '= strtrim(string(x(' + string(j, format = '(I0)') + ', ' $
                                  + string(i, format = '(I0)') + '), ' $
          + "format = '(" + string(ncol, format = '(I0)') + 'e23.' $
                          + string(prec, format = '(I0)') + ")'), 1)"

      r = execute(s)
    ENDFOR $
ELSE $                                    ; integer values
  FOR i = 0, nrow - 1 DO $
    FOR j = 0, ncol - 1 DO BEGIN
      s = 'sa(' + string(j, format = '(I0)') + ', ' $
                + string(i, format = '(I0)') + ') ' $
          + '= string(x(' + string(j, format = '(I0)') + ', ' $
                          + string(i, format = '(I0)') + '), ' $
          + "format = '(" + string(ncol, format = '(I0)') + "I0)')"

      r = execute(s)
    ENDFOR

printf, unit, sa

END

PRO write_datafile, filename, x, heading, prec

; Check input
IF n_params() EQ 3 THEN prec = 6

IF prec LT 0 THEN print, 'The precision must be greater than or equal to 0.'

; Open file for writing
openw, unit, filename, error = err, /get_lun

IF err NE 0 THEN printf, -2, !err_string

; Print heading
IF keyword_set(heading) THEN $
  printf, unit, format = '("# ", A, /, "#")', heading

printf, unit, '# This file is created by IDL.'

; Write the data

v = size(x) & ndim = v(0) & type = v(2)

IF ndim EQ 1 AND type EQ 8 THEN BEGIN     ; x is a structure of arrays
  nmat = n_tags(x)
  name = tag_names(x)
  printf, unit, nmat, format = '(I0)'
 
  FOR i = 0, nmat - 1 DO BEGIN
    s = 'v' + string(i, format = '(I0)') + ' = size(x.' + name(i) + ') & ' $
        + 'ndim = v' + string(i, format = '(I0)') + '(0) & ' $
        + 'ncol = v' + string(i, format = '(I0)') + '(1) & ' $
        + 'nrow = v' + string(i, format = '(I0)') + '(2) & ' $
        + 'mat = x.' + name(i)

    r = execute(s)

    IF ndim EQ 0 OR ndim EQ 1 THEN BEGIN  ; mat is a number or a line vector
      ncol = n_elements(mat) & nrow = 1
    ENDIF
 
    printf, unit, nrow, ncol, format = '(I0, 2X, I0)'
    write_mat, unit, mat, prec
  ENDFOR
ENDIF ELSE BEGIN                          ; x is a matrix or a column vector
  ncol = v(1) & nrow = v(2)
  
  IF ndim EQ 0 OR ndim EQ 1 THEN BEGIN    ; x is a number or a line vector
    ncol = n_elements(x) & nrow = 1
  ENDIF

  printf, unit, '1'
  printf, unit, nrow, ncol, format = '(I0, 2X, I0)'
  write_mat, unit, x, prec
ENDELSE

; Close the file
free_lun, unit
END
