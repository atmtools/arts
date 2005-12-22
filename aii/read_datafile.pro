;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Name:     read_datafile
;
;           Reads data from a file in ARTS data format (ASCII or BINARY).
;
;           The data is returned as an array or a structure or a
;           pointer array.
;           For example, to get the first matrix type 
;             m = x.(0) ( for ASCII )
;             m = *x[0] ( for Binary, the output is a pointer array )
;
; Usage:   x = read_datafile( filename [, /check, /binary ] )
;
; Inputs:   filename      full file name
; Optional: check         Keyword to check the data
;            
;           binary        Keyword to read ARTS binary file. The output
;                         is a pointer array. Be sure to free the
;                         memory after the use.  
;
; Output:   x             the data  
;
; History:  28.02.00  Wolfram Haas  Created    
;           22.12.00  Viju O. John  Included the option for reading
;                                   binary data
;
;------------------------------------------------------------------------------

FUNCTION read_datafile, filename, check = check, binary = binary

IF KEYWORD_SET( binary ) THEN BEGIN
    
    fid  =  HDF_OPEN( filename ) 
    
    IF( fid EQ -1 ) THEN $
      MESSAGE, filename + ' can not be opened, may be it does not exist ' + $
               'or it may not be a HDF file'

    lons  =  HDF_VD_LONE( fid )
    nlon  =  N_ELEMENTS( lons )

    IF( nlon LT 2 ) THEN $
      MESSAGE, 'The file may be emply'
    
    ;; If there are only two elements in lons the data is a scalar
    ;; vector or matrix. 
    IF( nlon EQ 2 ) THEN BEGIN 
        vid   =  HDF_VD_ATTACH( fid, lons[0] )
        nrec  =  HDF_VD_READ( vid, data )
        HDF_VD_DETACH, vid
        
        vid   =  HDF_VD_ATTACH( fid, lons[1] )
        nrec  =  HDF_VD_READ( vid, size )
        HDF_VD_DETACH, vid

        ;; Check whether it is a scalar
        IF TOTAL( size ) EQ 2.0 THEN $
          data  =  data ELSE $
          data  =  REFORM( data, REVERSE( size ) ) 
          
        pa  =  PTR_NEW( data )
    ENDIF ELSE BEGIN ;; It is an array of matrix
        ;; First two elements are for the number of matrices.
        vid   =  HDF_VD_ATTACH( fid, lons[0] )
        nrec  =  HDF_VD_READ( vid, nmat )
        HDF_VD_DETACH, vid
        
        pa = PTRARR( nmat ) ;; returns a pointer array

        vid   =  HDF_VD_ATTACH( fid, lons[1] )
        nrec  =  HDF_VD_READ( vid, size )
        HDF_VD_DETACH, vid        

        counter  =  0
        FOR ilon = 2L, N_ELEMENTS( lons ) - 2, 2 DO BEGIN  
            vid   =  HDF_VD_ATTACH( fid, lons[ilon + 1] )
            nrec  =  HDF_VD_READ( vid, size )
            HDF_VD_DETACH, vid
            
            vid   =  HDF_VD_ATTACH( fid, lons[ilon] )
            nrec  =  HDF_VD_READ( vid, data )
            data  =  REFORM( data, REVERSE( size ) )
            HDF_VD_DETACH, vid
            
            pa[counter] = PTR_NEW( data )
            
            counter  =  counter + 1
        ENDFOR
    ENDELSE
    HDF_CLOSE, fid 
    RETURN, pa
ENDIF ELSE BEGIN
    
    ;; Open file for reading
    openr, unit, filename, error = err, /get_lun
    
    IF err NE 0 THEN printf, -2, !err_string

    ;; Read until line does not begin with #
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
    
    ;; Read number of matrices
    reads, s, nmat
    IF nmat LT 1 THEN BEGIN
        print, 'Could not read number of matrices.'
        stop
    ENDIF
    

    pa = ptrarr(nmat) ;; returns a pointer array
    
    nrow = intarr(nmat) & ncol = intarr(nmat)
    
    ;; loop over all matrices
    FOR i = 0, nmat - 1 DO BEGIN
        ;; Read size
        s = '#'
        WHILE strmid(s, 0, 1) EQ '#' DO $
          IF NOT eof(unit) THEN readf, unit, s ELSE BEGIN
            print, 'Wrong number of matrices.'
            stop
        ENDELSE
        
        ;; Check size
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
            
            ;; Read the matrix row by row
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
        
        pa(i) = ptr_new(mat) ;; returns a pointer to pa(i)
    ENDFOR

    IF nmat EQ 1 THEN mat = *pa[0] ELSE BEGIN ;; create structure
        
        ;; I (VOJ) modified this part by making use of the IDL
        ;; function CREATE_STRUCT. Before It was done by creating a
        ;; string and then executing it. It had a limitation that
        ;; the string lenth could not exceed a certain limit.
        mat  =  CREATE_STRUCT( 'tag0',  *pa[0] ) 
        FOR imat = 1, nmat - 1 DO BEGIN
            tag  =  'tag' + STRTRIM( STRING( imat ), 2 )
            mat  =  CREATE_STRUCT( mat, tag, *pa[imat] )
        ENDFOR
        
    ENDELSE

    ;; Release memory used by the heap variable
    ptr_free, pa
    
    ;; Check if the file is now finished
    s = ''
    WHILE NOT eof(unit) DO BEGIN
        readf, unit, s
        
        s = strcompress(s, /remove_all)
        
        IF s NE '' THEN BEGIN
            print, 'There is some garbage at the end of the file.'
            stop
        ENDIF
    ENDWHILE
    
    ;; Close the file
    free_lun, unit  
ENDELSE


RETURN, mat

END
