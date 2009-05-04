PRO find_freq_H, doplot = doplot, xml=xml, aa = aa

SATID  =  ['TIROSN', 'NOAA6', 'NOAA7', 'NOAA8', 'NOAA9', 'NOAA10', 'NOAA11', $
           'NOAA12', 'NOAA14', 'NOAA15', 'NOAA16', 'NOAA17', 'NOAA18', $
           'NOAA19', 'METOPA']

;;centre_freq  =  []

nsat  =  n_elements( satid )

FOR isat = 0,  nsat - 1 DO BEGIN
    
    print, satid[isat]

    if keyword_set( xml ) then begin
        openw, srflun, 'xml/' + SATID[isat] + '_HIRS' + '.sideband_response.xml', /get_lun
        printf, srflun, '<?xml version="1.0"?>'
        printf, srflun, '<arts format="ascii" version="1">'
        printf, srflun, '<Array type="GriddedField1" nelem="19">'
        
        openw,  lolun, 'xml/' + SATID[isat] + '_HIRS' + '.f_backend.xml', /get_lun
        printf, lolun, '<?xml version="1.0"?>'
        printf, lolun, '<arts format="ascii" version="1">'
    endif

    if keyword_set( aa ) then begin
        openw, aasrflun, 'xml/' + SATID[isat] + '_HIRS' + '.srf.aa', /get_lun
        printf, aasrflun, '19'
    endif


    flt_file  =  'NOAA/' + SATID[isat] + '.FLT'
    all_data  =  read_flt_files( flt_file )
    
    nch  =  n_tags( all_data )
    
    str_nch  =  strtrim( string( nch ), 2 )
    if keyword_set( xml ) then printf, lolun, '<Vector nelem="' + str_nch + '">'

    for ich = 0, nch - 1 do begin
        print, 'Channel: ', ich + 1

        flt_data  =  all_data.( ich )

        wavnum  =  REFORM( flt_data[0, *] )
        filter  =  REFORM( flt_data[4, *] )
        
        freqnc  =  wavnum2freq( wavnum )
        sortin  =  SORT( freqnc )
        f_mono  =  TRANSPOSE( freqnc[sortin] )
        nfmono  =  n_elements( f_mono )
        
        hres_f_mono  =  transpose( vectornlinspace( f_mono[0, 0], f_mono[0, nfmono-1], 1001 ) )
        hres_filter  =  interpol( filter[sortin], f_mono, hres_f_mono )
        
        hmatrx  =  (1.0 / TOTAL( filter[sortin] )) * filter[sortin] 
        hres_hmatrx  =  transpose( (1.0 / TOTAL( hres_filter )) * hres_filter )
        
        print, hmatrx ## f_mono 
        
        centre_freq  =  hmatrx ## f_mono 
        
        print, 3e14 / ( hmatrx ## f_mono )
        print, 3e14 / ( hres_hmatrx ## hres_f_mono ) 
        
        ;; write_datafile, SATID[isat] + 'HIR12.f_mono.aa', f_mono, SATID[isat] + 'hires  f_mono '
        ;; write_datafile, SATID[isat] + 'HIR12.H.aa', hmatrx, SATID[isat] + ' H Matrix '
        ;; write_datafile, SATID[isat] + 'HIR12.hres_f_mono.aa', hres_f_mono, SATID[isat] + 'hres  f_mono '
        ;; write_datafile, SATID[isat] + 'HIR12.hres_H.aa', hres_hmatrx, SATID[isat] + ' hres H Matrix '
        
        if keyword_set( xml ) then begin
            
            ;;printf, lolun, centre_freq 
            
            ;;ch_no  =  strtrim( string( ich + 1 ), 2 )
            ;;if ich lt 9 then ch_no = '0' + ch_no
            
            ;;openw, lun, 'xml/' + SATID[isat] + '_HIRS' + ch_no + '.sideband_response.xml', /get_lun
            
            printf, srflun, '<GriddedField1 name="Backend channel response function">'
            printf, srflun, '<Vector name="Frequency" nelem="' + strtrim(string(n_elements(f_mono)),2) + '">'
            
            printf, srflun, f_mono - centre_freq[0]
            
            printf, srflun, '</Vector>'
            printf, srflun, '<Vector nelem="' + strtrim(string(n_elements(f_mono)),2) + '">'
            
            printf, srflun, transpose(hmatrx)
            
            printf, srflun, '</Vector>'
            printf, srflun, '</GriddedField1>'


            ;; writing centre frequency
            ;;printf, lolun, '<Vector nelem="1">'
            printf, lolun, centre_freq
            ;;printf, lolun, '</Vector>'
            
            
        endif


        if keyword_set( aa ) then begin
            printf, aasrflun, strtrim(string(n_elements(f_mono)),2), '  2'
            printf, aasrflun, [f_mono, transpose(hmatrx)]
        endif



        
        if keyword_set( doplot ) then begin
            set_my_eps_plot, 'n14_h12_srf.eps'
            !p.multi = [0, 3, 2]
            
            plot, hres_f_mono, hres_filter, pos = [0.2, 0.2, 0.8, 0.5], $
                  xtitle = 'Frequency [ GHz ]', ytitle = 'SRF [ ]'
            
            close_my_eps_plot, 'n14_h12_srf.eps';;, /show
        endif
    endfor

    if keyword_set( xml ) then begin
        
        printf, srflun, '</Array>'
        printf, srflun, '</arts>'
        
        free_lun, srflun


        printf, lolun, '</Vector>'
        printf, lolun, '</arts>'
        
        free_lun, lolun
    endif

    if keyword_set( aa ) then free_lun, aasrflun

ENDFOR

RETURN
END
