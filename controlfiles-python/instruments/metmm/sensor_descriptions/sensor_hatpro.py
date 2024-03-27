

# sensor_descriptions/sensor_hatpro.py
def sensor_hatpro(ws):
    """
    Sensor characteristics based on the pdf-presentation
    ftp://ftp.etl.noaa.gov/psd3/arctic/summit/mwr/0_docs/Summit_Datagrams_MicroWaveRadiometer.pdf
    
    Viewing angles( !Caution ground-based instrument! )
    This instrument has only 1 viewing angle
    """
    import numpy as np

    dlos = [180.0]  #0 (zenith)
        
    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)
    
    mm = [
        [22.240e9,   0., 0.,  230e6],   #0  (HATPRO-01)
        [23.040e9,   0., 0.,  230e6],   #1  (HATPRO-02)
        [23.840e9,   0., 0.,  230e6],   #2  (HATPRO-03)
        [25.440e9,   0., 0.,  230e6],   #3  (HATPRO-04)
        [26.240e9,   0., 0.,  230e6],   #4  (HATPRO-05)
        [27.840e9,   0., 0.,  230e6],   #5  (HATPRO-06)
        [31.400e9,   0., 0.,  230e6],   #6  (HATPRO-07)
        [51.260e9,   0., 0.,  182e6],   #7  (HATPRO-08)
        [52.280e9,   0., 0.,  179e6],   #8  (HATPRO-09)
        [53.860e9,   0., 0.,  188e6],   #9  (HATPRO-10)
        [54.940e9,   0., 0.,  170e6],   #10 (HATPRO-11)
        [56.660e9,   0., 0.,  704e6],   #11 (HATPRO-12)
        [57.300e9,   0., 0.,  927e6],   #12 (HATPRO-13)
        [58.000e9,   0., 0., 1854e6],   #13 (HATPRO-14)
        [90.000e9,   0., 0., 2000e6],   #14 (HATPRO-15)
        [150.000e9,  0., 0., 2000e6],   #15 (HATPRO-16)
    ]
    
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "AMSU-V", #0  (H-01)
            "AMSU-V", #1  (H-02)
            "AMSU-V", #2  (H-03)
            "AMSU-V", #3  (H-04)
            "AMSU-V", #4  (H-05)
            "AMSU-V", #5  (H-06)
            "AMSU-V", #6  (H-07)
            "AMSU-V", #7  (H-08)
            "AMSU-V", #8  (H-09)
            "AMSU-V", #9  (H-10)
            "AMSU-V", #10 (H-11)
            "AMSU-V", #11 (H-12)
            "AMSU-V", #12 (H-13)
            "AMSU-V", #13 (H-14)
            "AMSU-V", #14 (H-15)
            "AMSU-V", #15 (H-16)
        ])
        
    #ws.VectorSet( ws.met_mm_antenna, "" )
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(
        ws.met_mm_freq_number,
        [
            12, #0  (H-01)
            12, #1  (H-02)
            12, #2  (H-03)
            12, #3  (H-04)
            12, #4  (H-05)
            12, #5  (H-06)
            12, #6  (H-07)
            12, #7  (H-08)
            12, #8  (H-09)
            12, #9  (H-10)
            12, #10 (H-11)
            12, #11 (H-12)
            12, #12 (H-13)
            12, #13 (H-14)
            12, #14 (H-15)
            12, #15 (H-16)
         ])
    
    ws.VectorSet(ws.freq_spacing_tmp, np.array([1e9] * 16))
    
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)

    ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
    
    ws.VectorSetConstant(
        ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
    ws.Delete(ws.current_spacing)


if __name__ == "__main__":
    pass
