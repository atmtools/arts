import sys
import numpy as np

files = sys.argv[1:]

def treat_file(f):
    text = open(f).read().split("\n")
    n, title = text.pop(0).split(" ")
    n = int(n)
    out = {}

    new = False
    for line in text:
        s = line.split(' ')
        if len(s) != 2:
            new = True
            continue

        if new:
            new = False
            m, test = line.split(" ")
            if test not in out:
                out[test] = {}
            cur = out[test]
            continue

        name, time = line.split(" ")
        if name not in cur:
            cur[name] = []
        cur[name].append(float(time))

    for test in out:
        for name in out[test]:
            if len (out[test][name]) != n:
                raise ValueError("Number of runs is not consistent")
            out[test][name] = np.min(out[test][name]), np.max(out[test][name])
    
    return title, out


out = {}
for f in files:
    out[f] = treat_file(f)

for f in out:
    title = out[f][0]
    data = out[f][1]
    print(f".. list-table:: {title}")
    print( "    :widths: 50 10 10")
    print( "    :header-rows: 1")
    print()
    print( "    * - Test")
    print( "      - Min (μs)")
    print( "      - Max (μs)")
    for test in data:
        for name in data[test]:
            min_us = round(1e6 * data[test][name][0])
            max_us = round(1e6 * data[test][name][1])
            print(f"    * - ``{test}``")
            print(f"      - {min_us}")
            print(f"      - {max_us}")
    print()
