#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: obobryshev

def sensor_amsua(ws):
    """
    # sensor_descriptions/sensor_amsua.py
    """
    import numpy as np

    dlos = [
        -48.330000, #0
        -44.996897, #1
        -41.663793, #2
        -38.330690, #3
        -34.997586, #4
        -31.664483, #5
        -28.331379, #6
        -24.998276, #7
        -21.665172, #8
        -18.332069, #9
        -14.998966, #10
        -11.665862, #11
        -8.332759,  #12
        -4.999655,  #13
        -1.666552,  #14 (nadir)
    ]
    
    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)
        
    mm = [
        [23.8e9,         0.0e6,  0.0e6,  270e6],    #0  (1)
        [31.4e9,         0.0e6,  0.0e6,  180e6],    #1  (2)
        [50.3e9,         0.0e6,  0.0e6,  180e6],    #2  (3)
        [52.8e9,         0.0e6,  0.0e6,  400e6],    #3  (4)
        [53.596115e9,  115.0e6,  0.0e6,  170e6],    #4  (5)
        [54.400e9,       0.0e6,  0.0e6,  400e6],    #5  (6)
        [54.940e9,       0.0e6,  0.0e6,  400e6],    #6  (7)
        [55.50e9,        0.0e6,  0.0e6,  330e6],    #7  (8)
        [57.290344e9,    0.0e6,  0.0e6,  330e6],    #8  (9)
        [57.290344e9,  217.0e6,  0.0e6,   78e6],    #9  (10)
        [57.290344e9,  322.2e6, 48.0e6,   36e6],    #10 (11)
        [57.290344e9,  322.2e6, 22.0e6,   16e6],    #11 (12)
        [57.290344e9,  322.2e6, 10.0e6,    8e6],    #12 (13)
        [57.290344e9,  322.2e6,  4.5e6,    3e6],    #13 (14)
        [89.0e9,         0.0e9,  0.0e6, 2000e6],    #14 (15)    
    ]
    
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "AMSU-V", #0  (1)
            "AMSU-V", #1  (2)
            "AMSU-V", #2  (3)
            "AMSU-V", #3  (4)
            "AMSU-H", #4  (5)
            "AMSU-H", #5  (6)
            "AMSU-V", #6  (7)
            "AMSU-H", #7  (8)
            "AMSU-H", #8  (9)
            "AMSU-H", #9  (10)
            "AMSU-H", #10 (11)
            "AMSU-H", #11 (12)
            "AMSU-H", #12 (13)
            "AMSU-H", #13 (14)
            "AMSU-V"  #14 (15)
        ])
        
    #ws.VectorSet( ws.met_mm_antenna, "" )
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    #0 (1)
    #1 (2)
    #2 (3)
    #3 (4)
    #4 (5)
    #5 (6)
    #6 (7)
    #7 (8)
    #8 (9)
    #9 (10)
    #10 (11)
    #11 (12)
    #12 (13)
    #13 (14)     
    #14 (15) 
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    )
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for second accuracy (normal)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [1, 1, 1, 3, 3, 5, 5, 4, 4, 3, 3, 3, 4, 2, 1] 
    )
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for third accuracy (high)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [1, 1, 1, 8, 8, 16, 15, 11, 13, 9, 9, 9, 11, 6, 1] 
    )
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for fourth accuracy (reference)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [6, 1, 3, 23, 24, 44, 43, 34, 38, 26, 26, 27, 31, 17, 4] 
    )
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    
    ws.VectorSet(
        ws.freq_spacing_tmp, 
        np.array( [10e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9])
    )
    
    ws.Delete(ws.freq_number_tmp)

    ws.Extract(
        ws.met_mm_freq_number, ws.met_mm_available_accuracies, ws.met_mm_accuracy)
    
    ws.Extract(
        ws.current_spacing, ws.freq_spacing_tmp,            ws.met_mm_accuracy)

    ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
    
    ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
    ws.Delete(ws.current_spacing)
   

if __name__ == "__main__":
    pass   
