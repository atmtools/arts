
# sensor_descriptions/sensor_atms.py   
def sensor_atms(ws):
    """
    ARTS sensor description for ATMS simulations
    
    The following variables are set:
    
       antenna_dlos
       met_mm_backend
       met_mm_polarisation
       met_mm_freq_number
       met_mm_antenna
    
    """
    """
    Sensor characteristics based on the WMO-OSCAR at
    https://www.wmo-sat.info/oscar/instruments/view/53
    OSCAR stands for Observing Systems Capability Analysis and Review Tool
    Central frequency of Channel 16 (88.2GHz) is unclear (different sources
    give different numbers: either 88.2 or 89.5.
    
    Viewing angles
    There are 48 different angles, corresponding to one side of the ATMS scan.
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
         -0.555556  #47 (nadir)
    ]
        
    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)

    # Sensor response setup
    # ---
    # all frequencies are in Hz 
    # CenterFreq, Offset1, Offset2, Bandwidth; #ARTS channel index
    #                                         (Instrument channel)
    mm = [ 
        [23.800e9,     0.,       0.,       270e6 ],  #0  (1)
        [31.400e9,     0.,       0.,       180e6 ],  #1  (2)
        [50.300e9,     0.,       0.,       180e6 ],  #2  (3)
        [51.760e9,     0.,       0.,       400e6 ],  #3  (4)
        [52.800e9,     0.,       0.,       400e6 ],  #4  (5)
        [53.596e9,     0.115e9,  0.,       170e6 ],  #5  (6)
        [54.400e9,     0.,       0.,       400e6 ],  #6  (7)
        [54.940e9,     0.,       0.,       400e6 ],  #7  (8)
        [55.500e9,     0.,       0.,       330e6 ],  #8  (9)
        [57.290344e9,  0.,       0.,       330e6 ],  #9  (10)
        [57.290344e9,  0.217e9,  0.,       78e6  ],  #10 (11)
        [57.290344e9,  0.3222e9, 0.048e9,  36e6  ],  #11 (12)
        [57.290344e9,  0.3222e9, 0.022e9,  16e6  ],  #12 (13)
        [57.290344e9,  0.3222e9, 0.010e9,  8e6   ],  #13 (14)
        [57.290344e9,  0.3222e9, 0.0045e9, 3e6   ],  #14 (15)
        [88.200e9,     0.,       0.,       5000e6],  #15 (16)
        [165.500e9,    0.,       0.,       3000e6],  #16 (17)
        [183.310e9,    7.00e9,   0.,       2000e6],  #17 (18)
        [183.310e9,    4.50e9,   0.,       2000e6],  #18 (19)
        [183.310e9,    3.00e9,   0.,       1000e6],  #19 (20)
        [183.310e9,    1.80e9,   0.,       1000e6],  #20 (21)
        [183.310e9,    1.00e9,   0.,       500e6 ],  #21 (22)
    ]
    
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "AMSU-V", #0  (1)
            "AMSU-V", #1  (2)
            "AMSU-H", #2  (3)
            "AMSU-H", #3  (4)
            "AMSU-H", #4  (5)
            "AMSU-H", #5  (6)
            "AMSU-H", #6  (7)
            "AMSU-H", #7  (8)
            "AMSU-H", #8  (9)
            "AMSU-H", #9  (10)
            "AMSU-H", #10 (11)
            "AMSU-H", #11 (12)
            "AMSU-H", #12 (13)
            "AMSU-H", #13 (14)
            "AMSU-H", #14 (15)
            "AMSU-V", #15 (16)
            "AMSU-H", #16 (17)
            "AMSU-H", #17 (18)
            "AMSU-H", #18 (19)
            "AMSU-H", #19 (20)
            "AMSU-H", #20 (21)
            "AMSU-H"  #21 (22)
        ])
        
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(
        ws.met_mm_freq_number,
        [12] * 22, 
    )
    
    
    ws.VectorSet( 
        ws.freq_spacing_tmp, 
        np.array([1e9] *21 + [10e9])
    )
                
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)
    
    ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
    ws.VectorSetConstant(ws.met_mm_freq_spacing, ws.met_mm_nchannels, 
        ws.current_spacing)
    ws.Delete(ws.current_spacing)   


if __name__ == "__main__":
    pass
