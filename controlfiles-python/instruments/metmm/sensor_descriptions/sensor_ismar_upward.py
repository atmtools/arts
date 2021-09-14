#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: obobryshev

# sensor_descriptions/sensor_ismar_upward.py   
def sensor_ismar_upward(ws):
    """
    ARTS sensor description for ISMAR simulations
    
    This requires to run prepare_metmm.arts beforehand.
    
    This expects the following workspace variables to exist and to be set:
       met_mm_accuracy (Index)    Selection of accuracy level.
    
    The following variables are set:
       antenna_dlos
       met_mm_backend
       met_mm_polarisation
       met_mm_freq_number
       met_mm_freq_spacing
       met_mm_antenna
    """
    import numpy as np

    dlos = [
        -180.00, #0  
        -170.00, #1 
        -110.00, #2   
         -50.00, #3   
         -40.00, #4   
         -30.00, #5   
         -20.00, #6   
         -10.00, #7     
           0.00, #8 (nadir)    
          10.00, #9    
          70.00, #10   
         140.00, #11   
         150.00, #12   
         160.00, #13   
         170.00, #14 
    ]
    
    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)
        
    mm = [ 
        [118.7503e9,  1.10e9,  0.0e6,   400e6],  #0 
        [118.7503e9,  1.50e9,  0.0e6,   400e6],  #1 
        [118.7503e9,  2.10e9,  0.0e6,   800e6],  #2 
        [118.7503e9,  3.00e9,  0.0e6,  1000e6],  #3 
        [118.7503e9,  5.00e9,  0.0e6,  2000e6],  #4 
        [243.2000e9,  2.50e9,  0.0e6,  3000e6],  #5 
        [243.2000e9,  2.50e9,  0.0e6,  3000e6],  #6 
        [325.1500e9,  1.50e9,  0.0e6,  1600e6],  #7 
        [325.1500e9,  3.50e9,  0.0e6,  2400e6],  #8 
        [325.1500e9,  9.50e9,  0.0e6,  3000e6],  #9 
        [424.0000e9,      -1,     -1,      -1],  #10 Not yet implemented
        [424.0000e9,      -1,     -1,      -1],  #11 Not yet implemented
        [424.0000e9,      -1,     -1,      -1],  #12 Not yet implemented
        [424.0000e9,      -1,     -1,      -1],  #13 Not yet implemented
        [448.0000e9,  1.40e9,  0.0e6,  1200e6],  #14
        [448.0000e9,  3.00e9,  0.0e6,  2000e6],  #15
        [448.0000e9,  7.20e9,  0.0e6,  3000e6],  #16
        [664.0000e9,  4.20e9,  0.0e6,  5000e6],  #17
        [664.0000e9,  4.20e9,  0.0e6,  5000e6],  #18
        [874.0000e9,      -1,     -1,      -1],  #19 Not yet implemented
        [874.0000e9,      -1,     -1,      -1],  #20 Not yet implemented
    ]        
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "ISMAR-V", #0 
            "ISMAR-V", #1 
            "ISMAR-V", #2 
            "ISMAR-V", #3 
            "ISMAR-V", #4 
            "ISMAR-H", #5 
            "ISMAR-V", #6 
            "ISMAR-V", #7 
            "ISMAR-V", #8 
            "ISMAR-V", #9 
            "?",       #10
            "?",       #11
            "?",       #12
            "?",       #13
            "ISMAR-V", #14
            "ISMAR-V", #15
            "ISMAR-V", #16
            "ISMAR-H", #17
            "ISMAR-V", #18
            "?",       #19
            "?"        #20
        ])
        
    #ws.VectorSet( ws.met_mm_antenna, "" )
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [
            1,     #0 
             1,    #1 
             1,    #2 
             1,    #3 
             1,    #4 
             1,    #5 
             1,    #6 
             1,    #7 
             1,    #8 
             1,    #9 
            -1,    #10
            -1,    #11
            -1,    #12
            -1,    #13
             1,    #14
             1,    #15
             1,    #16
             1,    #17
             1,    #18
            -1,    #19
            -1,    #20
         ])
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for second accuracy (normal)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [
            2,     #0 
             1,    #1 
             3,    #2 
             4,    #3 
             5,    #4 
             6,    #5 
             6,    #6 
             11,   #7 
             7,    #8 
             5,    #9 
            -1,    #10
            -1,    #11
            -1,    #12
            -1,    #13
             2,    #14
             3,    #15
             5,    #16
             14,   #17
             14,   #18
            -1,    #19
            -1     #20
        ])
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for third accuracy (high)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [
            6,     #0 
            3,     #1 
            9,     #2 
            10,    #3 
            11,    #4 
            78,    #5 
            78,    #6 
            27,    #7 
            10,    #8 
            55,    #9 
            -1,    #10
            -1,    #11
            -1,    #12
            -1,    #13
             3,    #14
             5,    #15
             33,   #16
             58,   #17
             58,   #18
            -1,    #19
            -1     #20
        ])
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    # Number of frequencies for fourth accuracy (reference)
    ws.ArrayOfIndexSet(
        ws.freq_number_tmp,
        [
            18,    #0 
             9,    #1 
             28,   #2 
             36,   #3 
            -1,    #4 
            -1,    #5 
            -1,    #6 
             66,   #7 
             84,   #8 
            -1,    #9 
            -1,    #10
            -1,    #11
            -1,    #12
            -1,    #13
             6,    #14
             98,   #15
            -1,    #16
            -1,    #17
            -1,    #18
            -1,    #19
            -1     #20
        ])
    ws.Append(ws.met_mm_available_accuracies, ws.freq_number_tmp)
    
    ws.VectorSet(ws.freq_spacing_tmp, np.array([10e9, 1e9, 1e9, 1e9]))
    
    ws.Extract(ws.met_mm_freq_number, ws.met_mm_available_accuracies, ws.met_mm_accuracy)
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
    
    ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
    ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, 
                             ws.current_spacing)
    ws.Delete(ws.current_spacing)   


if __name__ == "__main__":
    pass
