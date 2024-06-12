#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: obobryshev

# sensor_descriptions/sensor_saphir.py   
def sensor_saphir(ws):
    """
    ARTS sensor description for SAPHIR L1A2 simulations
    
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
    Sensor characteristics based 
    
    Viewing angles
    There are 65 different angles, corresponding to one side of the SAPHIR Scan.
    Viewing angles definition from Table 3.2-2 Scan angle coverage. Distance
    between pixels is always 0.66:
    https://cnes.fr/fr/media/20130117level-1productdefed3rev4pdf
    """
    import numpy as np

    dlos = [
        -42.96, #0 (off-nadir)
        -42.30, #1
        -41.64, #2
        -40.98, #3
        -40.31, #4
        -39.65, #5
        -38.99, #6
        -38.33, #7
        -37.67, #8
        -37.01, #9
        -36.34, #10
        -35.68, #11
        -35.02, #12
        -34.36, #13
        -33.70, #14
        -33.04, #15
        -32.37, #16
        -31.71, #17
        -31.05, #18
        -30.39, #19
        -29.73, #20
        -29.07, #21
        -28.41, #22
        -27.74, #23
        -27.08, #24
        -26.42, #25
        -25.76, #26
        -25.10, #27
        -24.44, #28
        -23.77, #29
        -23.11, #30
        -22.45, #31
        -27.79, #32
        -21.13, #33
        -20.47, #34
        -19.80, #35
        -17.14, #36
        -18.48, #37
        -17.82, #38
        -17.16, #39
        -16.50, #40
        -15.83, #41 
        -15.17, #42
        -14.51, #43
        -13.85, #44
        -13.19, #45
        -12.53, #46
        -11.87, #47
        -11.20, #48
        -10.54, #49
        -9.88,  #50
        -9.22,  #51
        -8.56,  #52
        -7.90,  #53
        -7.23,  #54
        -6.57,  #55
        -5.91,  #56
        -5.25,  #57
        -4.59,  #58
        -3.93,  #59
        -3.26,  #60
        -2.60,  #61
        -1.94,  #62
        -1.27,  #63
        -0.61,  #64 (nadir)
    ]
    
    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)

    # all frequencies are in Hz 
    # CenterFreq, Offset1, Offset2, Bandwidth; #ARTS channel index
    #                                         (Instrument channel)
    mm = [  
        [183.31e9,  0.20e9,  0.,   200e6],   #0 (1)
        [183.31e9,  1.10e9,  0.,   350e6],   #1 (2)
        [183.31e9,  2.80e9,  0.,   500e6],   #2 (3)
        [183.31e9,  4.20e9,  0.,   700e6],   #3 (4)
        [183.31e9,  6.80e9,  0.,  1200e6],   #4 (5)
        [183.31e9,  11.00e9, 0.,  2000e6],   #5 (6)
    ]   
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "AMSU-V", #0 (1)
            "AMSU-V", #1 (2)
            "AMSU-V", #2 (3)
            "AMSU-V", #3 (4)
            "AMSU-V", #4 (5)
            "AMSU-V", #5 (6)
        ])
    
    # Antenna is not supported for now 
    #ws.VectorSet( ws.met_mm_antenna, "" )
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(
        ws.met_mm_freq_number,
        [    
            12, #0 (1)
            12, #1 (2)
            12, #2 (3)
            12, #3 (4)
            12, #4 (5)
            12, #5 (6)
        ])
    
    ws.VectorSet(ws.freq_spacing_tmp, np.array([1e9] * 6))
    
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
    
    ws.nrowsGet(ws.nrows, ws.met_mm_backend)
    
    ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.nrows, 
                             ws.current_spacing)
    ws.Delete(ws.current_spacing)   


if __name__ == "__main__":
    pass
