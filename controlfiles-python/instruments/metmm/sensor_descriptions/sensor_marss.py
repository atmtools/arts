#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: obobryshev

# sensor_descriptions/sensor_marss.py
def sensor_marss(ws):
    """
    ARTS sensor description for MARSS simulations
    
    This requires to run prepare_metmm.arts beforehand.
    
    This expects the following workspace variables to exist and to be set:
       met_mm_accuracy (Index)    Selection of accuracy level.
    
    The following variables are set:
       antenna_dlos
       met_mm_backend
       met_mm_polarisation
       met_mm_freq_number
       met_mm_freq_spacing
       met_mm_antenna
    """
    import numpy as np

    dlos = [
        -52.777777, #0 (off-nadir)
        -51.666666, #1
        -50.555555, #2
        -49.444444, #3
        -48.333333, #4
        -47.222222, #5
        -46.111111, #6
        -45.000000, #7
        -43.888888, #8
        -42.777777, #9
        -41.666666, #10
        -40.555555, #11
        -39.444444, #12
        -38.333333, #13
        -37.222222, #14
        -36.111111, #15
        -35.000000, #16
        -33.888889, #17
        -32.777777, #18
        -31.666666, #19
        -30.555555, #20
        -29.444444, #21
        -28.333333, #22
        -27.222222, #23
        -26.111111, #24
        -25.000000, #25
        -23.888889, #26
        -22.777778, #27
        -21.666666, #28
        -20.555555, #29
        -19.444444, #30
        -18.333333, #31
        -17.222222, #32
        -16.111111, #33
        -15.000000, #34
        -13.888889, #35
        -12.777778, #36
        -11.666667, #37
        -10.555555, #38
         -9.444444, #39
         -8.333333, #40
         -7.222222, #41
         -6.111111, #42
         -5.000000, #43
         -3.888889, #44 
         -2.777778, #45 
         -1.666667, #46
         -0.555556, #47 (nadir)
    ]
    
    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)
        
    mm = [ 
        [88.992e9,   1.075e9, 0.0,   650e6],   #0 (1)
        [157.075e9,  2.600e9, 0.0,  2600e6],   #1 (2)
        [183.248e9,  0.975e9, 0.0,   450e6],   #2 (3)
        [183.248e9,  3.000e9, 0.0,  1000e6],   #3 (4)
        [183.248e9,  7.000e9, 0.0,  2000e6],   #4 (5)
    ]
    
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "MARSS-V", #0 (1)
            "MARSS-H", #1 (2)
            "MARSS-H", #2 (3)
            "MARSS-H", #3 (4)
            "MARSS-H", #4 (5)
        ])
        
    #ws.VectorSet( ws.met_mm_antenna, "" )
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(
        ws.met_mm_freq_number,
        [
            12, #0  (1)
            12, #1  (2)
            12, #2  (3)
            12, #3  (4)
            12, #4  (5)
         ])
    
    
    ws.VectorSet(ws.freq_spacing_tmp, np.array([10e9, 1e9, 1e9, 1e9, 1e9]))
    
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
    
    ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
    ws.VectorSetConstant(
        ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
    ws.Delete(ws.current_spacing)   


if __name__ == "__main__":
    pass
