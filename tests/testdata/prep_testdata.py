import numpy as np
import pyarts
import sys
import os

nargs = len(sys.argv)

if nargs < 2:
    sys.exit(0)

dir = sys.argv[1]
file = f"{dir}/abs_lines.xml"
mode = sys.argv[2]

moldata = pyarts.classes.ArrayOfAbsorptionLines()
molpath = sys.argv[3] if os.path.isfile(sys.argv[3]) else sys.argv[4]
moldata.readxml(molpath)

if mode == "full":
    fl = float(sys.argv[5])
    fu = float(sys.argv[6])
    sl = float(sys.argv[7])
elif mode == "freq":
    fl = float(sys.argv[5])
    fu = float(sys.argv[6])
    sl = 0.0
elif mode == "str":
    fl = -np.infty
    fu = np.infty
    sl = float(sys.argv[5])
else:
    fl = -np.infty
    fu = np.infty
    sl = 0.0

for i in range(len(moldata)):
    rem = []
    for j in range(len(moldata[i].lines)):
        if moldata[i].lines[j].f0 < fl or moldata[i].lines[j].f0 > fu or moldata[i].lines[j].i0 < sl:
            rem.append(j)
    print(rem)

linedata = pyarts.classes.ArrayOfAbsorptionLines()
if sys.argv[-1] == "load":
    linedata.readxml(file)

for x in moldata:
    linedata.append(x)

linedata.savexml(file)
