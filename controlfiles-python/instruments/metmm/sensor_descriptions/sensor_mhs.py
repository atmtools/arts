#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: obobryshev

# sensor_descriptions/sensor_mhs.py
def sensor_mhs(ws):
    """
    ARTS sensor description for MHS simulations
    
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
    """
    Sensor characteristics based on KLM User's Guide at
    http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c3/sec3-9.htm
    
    Viewing angles
    There are 45 different angles, corresponding to one side of the MHS scan.
    """
    import numpy as np

    dlos = [
        -49.444444, #0 (off-nadir)
        -48.333333, #1
        -47.222222, #2
        -46.111111, #3
        -45.000000, #4
        -43.888888, #5
        -42.777777, #6
        -41.666666, #7
        -40.555555, #8
        -39.444444, #9
        -38.333333, #10
        -37.222222, #11
        -36.111111, #12
        -35.000000, #13
        -33.888889, #14
        -32.777777, #15
        -31.666666, #16
        -30.555555, #17
        -29.444444, #18
        -28.333333, #19
        -27.222222, #20
        -26.111111, #21
        -25.000000, #22
        -23.888889, #23
        -22.777778, #24
        -21.666666, #25
        -20.555555, #26
        -19.444444, #27
        -18.333333, #28
        -17.222222, #29
        -16.111111, #30
        -15.000000, #31
        -13.888889, #32
        -12.777778, #33
        -11.666667, #34
        -10.555555, #35
         -9.444444, #36
         -8.333333, #37
         -7.222222, #38
         -6.111111, #39
         -5.000000, #40
         -3.888889, #41
         -2.777778, #42
         -1.666667, #43
         -0.555556, #44 (nadir)
        ]
    
    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)

    # all frequencies are in Hz 
    # CenterFreq, Offset1, Offset2, Bandwidth; #ARTS channel index
    #                                          (Instrument channel) 
    mm = [  
        [89.000e9,   0.,      0.,  2800e6],   #0 (H1)
        [157.000e9,  0.,      0.,  2800e6],   #1 (H2)
        [183.311e9,  1.00e9,  0.,   500e6],   #2 (H3)
        [183.311e9,  3.00e9,  0.,  1000e6],   #3 (H4)
        [190.311e9,  0.,      0.,  2200e6],   #4 (H5)
    ]   
    
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "AMSU-V", #0 (H1)
            "AMSU-V", #1 (H2)
            "AMSU-H", #2 (H3)
            "AMSU-H", #3 (H4)
            "AMSU-V"  #4 (H5)
        ])
    
    # Antenna is not supported for now 
    #ws.VectorSet( ws.met_mm_antenna, "" )
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(ws.met_mm_freq_number,
        [12, 12, 12, 12, 12]) #0 (H1) #1 (H2) #2 (H3) #3 (H4) #4 (H5)

    ws.VectorSet(ws.freq_spacing_tmp, np.array([10e9] + [1e9] * 4))
    
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
    
    ws.nrowsGet(ws.nrows, ws.met_mm_backend) 

    ws.VectorSetConstant(
        ws.met_mm_freq_spacing, ws.nrows, ws.current_spacing)
    
    ws.Delete(ws.current_spacing)  




if __name__ == "__main__":
    pass
