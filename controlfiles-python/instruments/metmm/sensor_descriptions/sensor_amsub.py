# -*- coding: utf-8 -*-
# author: obobryshev

def sensor_amsub(ws):
    """
    ARTS sensor description for AMSU-B simulations
    
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
    http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c3/sec3-4.htm
    
    Viewing angles
    There are 45 different angles, corresponding to one side of the AMSU-B Scan.
    """
    import numpy as np

    dlos = [
        -48.95,  #0 (off-nadir)
        -47.85,  #1
        -46.75,  #2
        -45.65,  #3
        -44.55,  #4
        -43.45,  #5
        -42.35,  #6
        -41.25,  #7
        -40.15,  #8
        -39.05,  #9
        -37.95,  #10
        -36.85,  #11
        -35.75,  #12
        -34.65,  #13
        -33.55,  #14
        -32.45,  #15
        -31.35,  #16
        -30.25,  #17
        -29.15,  #18
        -28.05,  #19
        -26.95,  #20
        -25.85,  #21
        -24.75,  #22
        -23.65,  #23
        -22.55,  #24
        -21.45,  #25
        -20.35,  #26
        -19.25,  #27
        -18.15,  #28
        -17.05,  #29
        -15.95,  #30
        -14.85,  #31
        -13.75,  #32
        -12.65,  #33
        -11.55,  #34
        -10.45,  #35
        -9.35,   #36
        -8.25,   #37
        -7.15,   #38
        -6.05,   #39
        -4.95,   #40
        -3.85,   #41
        -2.75,   #42
        -1.65,   #43
        -0.55,   #44 (nadir) 
    ]

    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)
    
    # all frequencies are in Hz 
    # CenterFreq, Offset1, Offset2, Bandwidth; #ARTS channel index
    #                                         (Instrument channel)
    mm = [
        [89.00e9,   0.90e9,  0.,  1000e6],
        [150.00e9,  0.90e9,  0.,  1000e6],
        [183.31e9,  1.00e9,  0.,  500e6],
        [183.31e9,  3.00e9,  0.,  1000e6],
        [183.31e9,  7.00e9,  0.,  2000e6],
    ]   
    
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "AMSU-V",  #0 (16)
            "AMSU-V",  #1 (17)
            "AMSU-V",  #2 (18)
            "AMSU-V",  #3 (19)
            "AMSU-V",  #4 (20)
        ])
    # Antenna is not supported for now     
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [1, 1, 1, 1, 1], #0 (16) #1 (17) #2 (18) #3 (19) #4 (20)
    )
    
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for second accuracy (normal)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [1, 2, 2, 2, 3], #0 (16) #1 (17) #2 (18) #3 (19) #4 (20)
    )
       
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for third accuracy (high)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [1, 18, 20, 7, 10], #0 (16) #1 (17) #2 (18) #3 (19) #4 (20)
    )
       
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for fourth accuracy (reference)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [2, 23, 67, 19, 25], #0 (16) #1 (17) #2 (18) #3 (19) #4 (20)
    )
       
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    
    ws.VectorSet(ws.freq_spacing_tmp, np.array([10e9, 1e9, 1e9, 1e9, 1e9]))
    
    ws.Extract(ws.met_mm_freq_number, ws.met_mm_available_accuracies, ws.met_mm_accuracy)
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
    
    ws.nrowsGet(ws.nrows, ws.met_mm_backend)
    ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.nrows, ws.current_spacing)
    ws.Delete(ws.current_spacing)


if __name__ == "__main__":
    pass
