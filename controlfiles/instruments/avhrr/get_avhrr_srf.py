#!/usr/bin/env python

"""Downloads AVHRR SRF from NOAA KLM User's Guide
"""

# Created by Gerrit Holl, July 2011
#
# Last updated: $Date$

import commands
import datetime

now = datetime.datetime.now

# import scipy.integrate

from string import printable, digits
from numpy import array

import PyARTS.artsXML

from PyARTS.arts_types import GriddedField1 as GF
from PyARTS.physics import wavelength2frequency, wavenumber2frequency

baseurl = "http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/d/app-d%d.htm"
urls = dict(
    NOAA15=baseurl % 1,
    NOAA16=baseurl % 2,
    NOAA17=baseurl % 3,
    NOAA18=baseurl % 4,
    METOPA=baseurl % 5,
    NOAA19=baseurl % 6,
)

start_marker = "Channel 3B Channel 4 Channel 5"
channels = "3B 4 5".split()

cmd = "lynx -dump -nolist %s"
out = "AVHRR_SRF_%(sat)s_%(channel)s.xml"

# HTML bugs...
repl_table = {
    "< /TR>": "",
    " e": " ",
    "R-": "E-",
    "4X": "4",
    "9090E-01": "9.90E-1",
    "1.00E+01": "1.00E-00",
    "3025E-04": "3.25E-04",
}


def html2srf(s):
    """Gets AVHRR 3B, 4, 5 SRF from KLM User's Guide HTML

    Returns 3-tuple with SRF as GF1 (3B, 4, 5)
    """
    start = s.find(start_marker)
    started = False
    ch = [[[], []], [[], []], [[], []]]
    for line in s[start:].split("\n"):
        if line.strip() and line.strip()[0].isdigit():
            # bugs in html-code necessisates:
            for k, v in repl_table.items():
                line = line.replace(k, v)
            started = True
            L = line.strip().split()
            for i in range(3):
                try:
                    wl = float(L[2 * i])
                    freq = wavelength2frequency(wl / 1e6)
                    wgt = float(L[2 * i + 1])
                except (ValueError, IndexError):  # probably only 3B has values here
                    # print "only partly understood:", line
                    pass
                else:
                    # some double entries appear
                    if len(ch[i][0]) == 0 or not freq == ch[i][0][-1]:
                        ch[i][0].append(freq)
                        ch[i][1].append(wgt)
        elif started:  # finished
            break

    return tuple(GF(array(ch[i][0]), array(ch[i][1])) for i in xrange(3))


def cleanup(srf):
    """Cleans up srf (GF1), removing everything beyond the negatives
    """
    data = srf.data
    middle = data.size // 2
    if (data < 0).any():
        first = (data[:middle] < 0).nonzero()[0][-1]
        last = (data[middle:] < 0).nonzero()[0][0] + middle + 1
        return GF(srf.axes[0][first:last], srf.data[first:last])
    else:
        print data
        return srf


def centres(s):
    """Gets centre frequencies for channels 3B, 4, 5 from HTML
    """

    L = []
    tostrip = str(set(printable) - set(digits))
    for line in s.split("\n"):
        if "Wavenumber =" in line:
            wn_cm1 = float(line.strip(tostrip))
            wn_m1 = wn_cm1 * 1e2
            freq = wavenumber2frequency(wn_m1)
            L.append(freq)
    return L[3:]


def ch2gf(sat):
    (status, output) = commands.getstatusoutput(cmd % urls[sat])
    c = centres(output)
    allnewgf = PyARTS.arts_types.ArrayOfGriddedField1()
    for (i, gf) in enumerate(html2srf(output)):
        # sort according to relative frequency
        relfreq = gf.axes[0] - c[i]
        ii = relfreq.argsort()
        newgf = cleanup(GF(relfreq[ii], gf.data[ii]))
        # integral = scipy.integrate.trapz(newgf.axes[0], newgf.data)
        # newgf.data /= integral
        # newgf.data[newgf.data<0]=0
        allnewgf.append(newgf)
        # print c[i], newgf.axes[0]
        # gf.save("AVHRR_SRF_%s_%s" % (sat, channels[i]))
    PyARTS.artsXML.save(array(c), "%s_AVHRR.f_backend.xml" % sat)
    allnewgf.save("%s_AVHRR.backend_channel_response.xml" % sat)


for sat in urls.keys():
    print now(), "doing", sat
    ch2gf(sat)
# print ch
# end = output.find(end_marker)
# print start, end
# print output[start:end]
