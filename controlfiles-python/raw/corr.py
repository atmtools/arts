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
tg = pyarts.classes.from_workspace(ws.time_grid)
ts = pyarts.classes.from_workspace(ws.time_stamps)
ys = pyarts.classes.from_workspace(ws.ybatch)
ceps = pyarts.classes.from_workspace(ws.covmat_sepsbatch)

# Current time
some_time = pyarts.classes.Time()
some_time.sec = 1589895540.0

# Set up the signal every 20 minutes from now
for i in range(N):
    sec = i*60*20
    ts.append(some_time + sec)
    np.random.seed(i)
    ys.append(f(sec, n))

# Sort the time and signal
ws.time_stampsSort(ws.time_grid, ws.time_stamps, ws.time_stamps)
ws.time_stampsSort(ws.ybatch, ws.time_stamps, ws.ybatch)
ws.time_stampsSort(ws.time_stamps, ws.time_stamps, ws.time_stamps)  # Must be last...

# Apply a simplified troposåheric correction
ws.create_variable("ArrayOfIndex", "range")
ws.create_variable("Vector", "trop_temp")
ws.create_variable("Numeric", "targ_temp")
rang = pyarts.classes.from_workspace(ws.range)
trop_temp = pyarts.classes.from_workspace(ws.trop_temp)
targ_temp = pyarts.classes.from_workspace(ws.targ_temp)
targ_temp.set(2.73)
trop_temp.data = [273.15 for i in range(N)]
trp_range = list(np.append(np.array([i for i in range(K)]), np.array([n-i-1 for i in range(K)])))
trp_range.sort()
rang.set(trp_range)

# Apply a simple tropospheric correction and then average the
ws.ybatchTroposphericCorrectionNaiveMedianForward(ws.ybatch_corr, ws.ybatch, ws.range, ws.trop_temp, ws.targ_temp)
ws.ybatchTimeAveraging(time_step="24 h", disregard_first=1, disregard_last=1)

# A simple
ncfn = os.path.splitext(__file__)[0]+'_ref.nc'
# pyarts.classes.netcdf.save(ncfn, [ws.covmat_sepsbatch, ws.ybatch])
refs = pyarts.classes.netcdf.load(ncfn)
