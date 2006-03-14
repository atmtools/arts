;******************************************************************************
;+
;NAME:
;           write_datafile
;Purpose:
;           Writes data to a file in ARTS format.
;
;           See also read_datafile.
;
;Usage:   write_datafile, filename, x, heading [, prec] [, /binary ][, /atype ]
;
;Inputs:   filename    full file name
;           x           the data to store
;           heading     heading text
;                       The function puts in '# ' at the start of each line. 
;                       Heading can be empty ('')
;Optional: prec        number of decimals to use, default 6
;                       If prec = 0, integer values are assumed.
;                             
;Keywords: binary      To write the data in ARTS binary format
;          atype       ARTS data type. The following types are allowed:
;                      NUMERIC, VECTOR, MATRIX, AOVECTOR, AOMATRIX.
;Output:   -
;
; History:  28.02.01  Wolfram Haas
;           22.02.06  Viju O. John : Included binary option
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

PRO write_binary, fid, data, dname, dtype, size, index = index
;; This procedure writes data into an already opened HDF file

IF KEYWORD_SET( index ) THEN BEGIN
    vdat = HDF_VD_ATTACH(fid, -1, /WRITE)
    HDF_VD_FDEFINE, vdat, dtype, /LONG
    HDF_VD_SETINFO, vdat, name = dname, class = 'UINT'
    HDF_VD_ATTRSET, vdat, -1, 'SIZE', size, /LONG
    HDF_VD_WRITE, vdat, dtype, data
    hdf_vd_detach, vdat
ENDIF ELSE BEGIN
    vdat = HDF_VD_ATTACH(fid, -1, /WRITE)
    HDF_VD_FDEFINE, vdat, dtype, /DOUBLE
    HDF_VD_SETINFO, vdat, name = dname, class = 'DOUBLE'
    HDF_VD_ATTRSET, vdat, -1, 'SIZE', size, /LONG
    HDF_VD_WRITE, vdat, dtype, data
    hdf_vd_detach, vdat
ENDELSE

END 

PRO write_datafile, filename, x, heading, prec, binary = binary, atype = atype

IF( N_ELEMENTS( x ) LT 1 ) THEN MESSAGE, 'x is undefined'

IF KEYWORD_SET( binary ) THEN BEGIN 
    ;; I am sure that this part can be written in a better way, for
    ;; example, a lot of duplication can be omitted. This is a quick
    ;; and dirty attemp.
    IF NOT KEYWORD_SET( atype ) THEN $
      MESSAGE, 'To write a binary file you need to specify ARTS data type'
    
    fid  =  HDF_OPEN( filename, /CREATE, /WRITE )

    atype  =  STRUPCASE( atype )
    
    type_check  =  STRCMP( atype, ['NUMERIC', 'VECTOR', 'MATRIX', 'AOVECTOR', 'AOMATRIX'] )
    IF( TOTAL( type_check ) NE 1.0 ) THEN MESSAGE, 'ARTS data type given is ' + atype + $
      ', but should be one of [NUMERIC, VECTOR, MATRIX, AOVECTOR, AOMATRIX]'

    IF( STRCMP( atype, 'NUMERIC' ) ) THEN BEGIN
        ;; The input data can be given as a pointer
        type  =  SIZE( x, /TYPE )
        IF type EQ 10 THEN x  =  *x 
      
        IF( N_ELEMENTS( x ) NE 1 ) THEN MESSAGE, 'Data should contain only 1 element.'
        dname  =  'NUMERIC'
        dtype  =  'SCALAR'
        nrow   =  1
        ncol   =  1
        size   =  LONG( [nrow, ncol] )
        x  =  REFORM( [x], nrow * ncol )
        write_binary, fid, x, dname, dtype, size
    ENDIF

    IF( STRCMP( atype, 'VECTOR' ) ) THEN BEGIN
        ;; The input data can be given as a pointer
        type  =  SIZE( x, /TYPE )
        IF type EQ 10 THEN x  =  *x 
        
        dname  =  'VECTOR'
        dtype  =  'VECTOR'
        dims   =  SIZE( x, /DIM )
        IF( SIZE( x, /N_DIM ) EQ 1 ) THEN BEGIN
            nrow  =  1 & ncol  =  dims 
        ENDIF ELSE BEGIN
            nrow  =  dims[1] & ncol  =  dims[0]
        ENDELSE
        size   =  LONG( [nrow, ncol] )
        x  =  REFORM( x, nrow * ncol )
        write_binary, fid, x, dname, dtype, size
    ENDIF

    IF( STRCMP( atype, 'MATRIX' ) ) THEN BEGIN
        ;; The input data can be given as a pointer
        type  =  SIZE( x, /TYPE )
        IF type EQ 10 THEN x  =  *x 
      
        dname  =  'MATRIX'
        dtype  =  'MATRIX'
        dims   =  SIZE( x, /DIM )
        IF( SIZE( x, /N_DIM ) NE 2 ) THEN $
          MESSAGE, 'A matrix should have 2 dimensions'
        nrow  =  dims[1] & ncol  =  dims[0]
        size  =  LONG( [nrow, ncol] )
        x  =  REFORM( x, nrow * ncol )
        write_binary, fid, x, dname, dtype, size
    ENDIF

    IF( STRCMP( atype, 'AOVECTOR' ) ) THEN BEGIN
        ;; The input data can be given as an IDL structure or pointer,
        ;; or sometimes it can be just a vector in that cane n_tag
        ;; will be 1 
        type  =  SIZE( x, /TYPE )
        
        IF type EQ 8 THEN ntag = N_TAGS( x ) ELSE $
          IF type EQ 10 THEN ntag = N_ELEMENTS( x ) ELSE ntag = 1
          
        
        write_binary, fid, ntag, 'N_VECTOR', 'SCALAR', LONG( [1, 1] ), /index
        dtype  =  'VECTOR'

        FOR itag = 0, ntag - 1 DO BEGIN
            dname  =  'VECTOR' + STRTRIM( STRING( itag ), 2 )
            
            IF type EQ 8 THEN data_tag = x.( itag ) ELSE $
              IF type EQ 10 THEN data_tag = *x[itag] ELSE data_tag = x

            dims  =  SIZE( data_tag, /DIM )
            IF( SIZE( data_tag, /N_DIM ) EQ 1 ) THEN BEGIN
                nrow  =  1 & ncol  =  dims 
            ENDIF ELSE BEGIN
                nrow  =  dims[1] & ncol  =  dims[0]
            ENDELSE
            size  =  LONG( [nrow, ncol] )
            data_tag  =  REFORM( data_tag, nrow * ncol )
            write_binary, fid, data_tag, dname, dtype, size
        ENDFOR
    ENDIF
    
    IF( STRCMP( atype, 'AOMATRIX' ) ) THEN BEGIN
        ;; The input data can be given as an IDL structure or pointer,
        ;; or sometimes it can be just a vector in that cane n_tag
        ;; will be 1 
        type  =  SIZE( x, /TYPE )
        
        IF type EQ 8 THEN ntag = N_TAGS( x ) ELSE $
          IF type EQ 10 THEN ntag = N_ELEMENTS( x ) ELSE ntag = 1
        
        write_binary, fid, ntag, 'N_MATRIX', 'SCALAR', LONG( [1, 1] ), /index
        dtype  =  'MATRIX'

        FOR itag = 0, ntag - 1 DO BEGIN
            dname  =  'MATRIX' + STRTRIM( STRING( itag ), 2 )
            
            IF type EQ 8 THEN data_tag = x.( itag ) ELSE $
              IF type EQ 10 THEN data_tag = *x[itag] ELSE data_tag = x

            dims   =  SIZE( data_tag, /DIM )
            IF( SIZE( data_tag, /N_DIM ) NE 2 ) THEN $
              MESSAGE, 'A matrix should have 2 dimensions'
            nrow  =  dims[1] & ncol  =  dims[0]
            size  =  LONG( [nrow, ncol] )
            data_tag  =  REFORM( data_tag, nrow * ncol )
            write_binary, fid, data_tag, dname, dtype, size
        ENDFOR
    ENDIF

    HDF_CLOSE, fid
    
ENDIF ELSE BEGIN ;; Code written by Haas

    ; Check input
    IF n_params() EQ 3 THEN prec = 6

    IF prec LT 0 THEN print, 'The precision must be greater than or equal to 0.'
    
    ;; Open file for writing
    openw, unit, filename, error = err, /get_lun

    IF err NE 0 THEN printf, -2, !err_string

    ;; Print heading
    IF keyword_set(heading) THEN $
      printf, unit, format = '("# ", A, /, "#")', heading

    printf, unit, '# This file is created by IDL.'

    ;; Write the data

    v = size(x) & ndim = v(0) & type = v(2)

    IF ndim EQ 1 AND type EQ 8 THEN BEGIN ; x is a structure of arrays
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
            
            IF ndim EQ 0 OR ndim EQ 1 THEN BEGIN ; mat is a number or a line vector
                ncol = n_elements(mat) & nrow = 1
            ENDIF
            
            printf, unit, nrow, ncol, format = '(I0, 2X, I0)'
            write_mat, unit, mat, prec
        ENDFOR
    ENDIF ELSE BEGIN            ; x is a matrix or a column vector
        ncol = v(1) & nrow = v(2)
        
        IF ndim EQ 0 OR ndim EQ 1 THEN BEGIN ; x is a number or a line vector
            ncol = n_elements(x) & nrow = 1
        ENDIF
        
        printf, unit, '1'
        printf, unit, nrow, ncol, format = '(I0, 2X, I0)'
        write_mat, unit, x, prec
    ENDELSE

    ;; Close the file
    free_lun, unit
ENDELSE

END
