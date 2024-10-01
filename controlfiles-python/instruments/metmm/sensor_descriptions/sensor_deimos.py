
# sensor_descriptions/sensor_deimos.arts   
def sensor_deimos(ws):
    ws.Error("The MetMM file for DEIMOS is not yet finalised!!!")
    import numpy as np

    dlos = []
        
    ws.MatrixSet(ws.antenna_dlos, np.array( [dlos] ).T )
    mm = [ 
        [23.8e9,  0.07e9,  0.0e6,   127e6],  #0 (1, V only)
        [23.8e9,  0.07e9,  0.0e6,   127e6],  #1 (1, H only)
        [50.1e9,  0.08e9,  0.0e6,    82e6],  #2 (2, V only)
        [50.1e9,  0.08e9,  0.0e6,    82e6],  #3 (2, H only)
    ]
    
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "ISMAR-V", #0 (1, V only)
            "ISMAR-H", #1 (1, H only)
            "ISMAR-V", #2 (2, V only)
            "ISMAR-H"  #3 (2, H only)
        ])
        
    #ws.VectorSet( ws.met_mm_antenna, "" )
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(
        ws.met_mm_freq_number,
        [12, 12, 12, 12]   #0 (1, V only) #1 (1, H only) #2 (2, V only) 
#3 (2, H only)
    )
    
    
    ws.VectorSet( 
        ws.freq_spacing_tmp, 
        np.array( [10e9, 1e9, 1e9, 1e9] )
    )
    
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
    
    ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
    ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, 
                             ws.current_spacing)
    ws.Delete(ws.current_spacing)   


if __name__ == "__main__":
    pass
