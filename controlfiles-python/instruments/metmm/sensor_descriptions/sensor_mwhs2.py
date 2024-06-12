#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: obobryshev

# sensor_descriptions/sensor_mwhs2.py   
def sensor_mwhs2(ws):
    """
    ARTS sensor description for MWS-2 simulations onboard FY-3C satellite.
    Lanuched in Sept. 2013
    
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
    # Viewing angle information is taken from Table 15:
    # https://directory.eoportal.org/web/eoportal/satellite-missions/f/fy-3#3Y4Y413eKram
    # There are 48 different angles, corresponding to one side of the MWHS-2 scan.
    """
    import numpy as np

    dlos = [
        -53.3500,  #0 (off-nadir)
        -52.2500,  #1
        -51.1500,  #2
        -50.0500,  #3
        -48.9500,  #4
        -47.8500,  #5
        -46.7500,  #6
        -45.6500,  #7
        -44.5500,  #8
        -43.4500,  #9
        -42.3500,  #10
        -41.2500,  #11
        -40.1500,  #12
        -39.0500,  #13
        -37.9500,  #14
        -36.8500,  #15
        -35.7500,  #16
        -34.6500,  #17
        -33.5500,  #18
        -32.4500,  #19
        -31.3500,  #20
        -30.2500,  #21
        -29.1500,  #22
        -28.0500,  #23
        -26.9500,  #24
        -25.8500,  #25
        -24.7500,  #26
        -23.6500,  #27
        -22.5500,  #28
        -21.4500,  #29
        -20.3500,  #30
        -19.2500,  #31
        -18.1500,  #32
        -17.0500,  #33
        -15.9500,  #34
        -14.8500,  #35
        -13.7500,  #36
        -12.6500,  #37
        -11.5500,  #38
        -10.4500,  #39
         -9.3500,  #40
         -8.2500,  #41
         -7.1500,  #42
         -6.0500,  #43
         -4.9500,  #44
         -3.8500,  #45
         -2.7500,  #46
         -1.6500,  #47
         -0.5500,  #48 (nadir)
    ]
   
    ws.MatrixSet(ws.antenna_dlos, np.array([dlos]).T)

    # Sensor response setup
    # ---
    # all frequencies are in Hz 
    # CenterFreq, Offset1, Offset2, Bandwidth; #ARTS channel index
    #                                         (Instrument channel)
    mm = [ 
        [89.000e9,   0.,      0.,  1500e6],  #0  (H1 )
        [118.750e9,  5.00e9,  0.,  2000e6],  #1  (H2 )
        [118.750e9,  3.00e9,  0.,  1000e6],  #2  (H3 )
        [118.750e9,  2.50e9,  0.,   200e6],  #3  (H4 )
        [118.750e9,  1.10e9,  0.,   200e6],  #4  (H5 )
        [118.750e9,  0.80e9,  0.,   200e6],  #5  (H6 )
        [118.750e9,  0.30e9,  0.,   165e6],  #6  (H7 )
        [118.750e9,  0.20e9,  0.,   100e6],  #7  (H8 )
        [118.750e9,  0.08e9,  0.,    20e6],  #8  (H9 )
        [150.000e9,  0.,      0.,  1500e6],  #9  (H10)
        [183.311e9,  7.00e9,  0.,  2000e6],  #10 (H11)
        [183.311e9,  4.50e9,  0.,  2000e6],  #11 (H12)
        [183.311e9,  3.00e9,  0.,  1000e6],  #12 (H13)
        [183.311e9,  1.80e9,  0.,   700e6],  #13 (H14)
        [183.311e9,  1.00e9,  0.,   500e6],  #14 (H15)
    ]
    
    ws.MatrixSet(ws.met_mm_backend, np.array(mm))
    
    ws.ArrayOfStringSet(
        ws.met_mm_polarisation, 
        [
            "AMSU-V", #0  (H1 )
            "AMSU-H", #1  (H2 )
            "AMSU-H", #2  (H3 )
            "AMSU-H", #3  (H4 )
            "AMSU-H", #4  (H5 )
            "AMSU-H", #5  (H6 )
            "AMSU-H", #6  (H7 )
            "AMSU-H", #7  (H8 )
            "AMSU-H", #8  (H9 )
            "AMSU-V", #9  (H10)
            "AMSU-H", #10 (H11)
            "AMSU-H", #11 (H12)
            "AMSU-H", #12 (H13)
            "AMSU-H", #13 (H14)
            "AMSU-H", #14 (H15)
        ])

    # Antenna is not supported for now
    #ws.VectorSet( ws.met_mm_antenna, "" )
    ws.Touch(ws.met_mm_antenna)
    
    # Number of frequencies for first accuracy (fast)
    ws.ArrayOfIndexSet(
        ws.met_mm_freq_number,
        [12] * 15,
    )
    
    ws.VectorSet( 
        ws.freq_spacing_tmp, 
        np.array([10e9] + [1e9] * 14))
    
    ws.Extract(ws.current_spacing, ws.freq_spacing_tmp, ws.met_mm_accuracy)

    ws.nrowsGet(ws.met_mm_nchannels, ws.met_mm_backend)
    
    ws.VectorSetConstant(
        ws.met_mm_freq_spacing, ws.met_mm_nchannels, ws.current_spacing)
    ws.Delete(ws.current_spacing)   


if __name__ == "__main__":
    pass
