#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

from SphereDataPhysical import SphereDataPhysical


if len(sys.argv) <= 1:
    print("Usage:")
    print(" "+sys.argv[0]+" [file.csv] [outputfile.ext] [min] [max]")
    sys.exit(1)


print("Loading data from '"+sys.argv[1]+"'")

spheredata = SphereDataPhysical(sys.argv[1])
data = np.flip(spheredata.data, axis=0)


print("min: "+str(np.min(data)))
print("max: "+str(np.max(data)))

vmin = None
vmax = None


if len(sys.argv) > 4:
    vmin = float(sys.argv[3])
    vmax = float(sys.argv[4])


fig, ax = plt.subplots()

im = plt.imshow(data, vmin=vmin, vmax=vmax)

cbar = ax.figure.colorbar(im, ax=ax)

#filename = sys.argv[1].split('/')[0]
#plt.title(filename, fontsize=8)

plt.tight_layout()

if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()
