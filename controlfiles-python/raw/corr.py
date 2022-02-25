import os
import pyarts
import numpy as np

# Size of problem
N = 500
n = 51
K = 3

# Workspace
ws = pyarts.workspace.Workspace()

# A simple noisy "signal" generator
def f(sec: float, n: int):
    eps = np.random.normal(size=(n))
    off = 20 + 5*abs(sec)/86000
    
    return 30*np.exp(-(np.linspace(-3, 3, n))**2) + off + eps

# Help to access
ws.sensor_time = []
ws.time_stamps = []
ws.ybatch = []

# Current time
some_time = pyarts.classes.Time()
some_time.sec = 1589895540.0

# Set up the signal every 20 minutes from now
for i in range(N):
    sec = i*60*20
    ws.time_stamps.value.append(some_time + sec)
    np.random.seed(i)
    ws.ybatch.value.append(f(sec, n))

# Sort the time and signal
ws.time_stampsSort(ws.sensor_time, ws.time_stamps, ws.time_stamps)
ws.time_stampsSort(ws.ybatch, ws.time_stamps, ws.ybatch)
ws.time_stampsSort(ws.time_stamps, ws.time_stamps, ws.time_stamps)  # Must be last...

# Apply a simplified troposåheric correction
ws.range = pyarts.classes.ArrayOfIndex()
ws.trop_temp = pyarts.classes.Vector()
ws.targ_temp = pyarts.classes.Numeric()
ws.targ_temp = 2.73
ws.trop_temp = [273.15 for i in range(N)]
trp_range = list(np.append(np.array([i for i in range(K)]), np.array([n-i-1 for i in range(K)])))
trp_range.sort()
ws.range = trp_range

# Apply a simple tropospheric correction and then average the
ws.ybatchTroposphericCorrectionNaiveMedianForward(ws.ybatch_corr, ws.ybatch, ws.range, ws.trop_temp, ws.targ_temp)
ws.ybatchTimeAveraging(time_step="24 h", disregard_first=1, disregard_last=1)

