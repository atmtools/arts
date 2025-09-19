import pyarts3 as pyarts
import os
import sys


def recurse(path):
    paths = []
    if os.path.isdir(path):
        items = os.listdir(path)
        items.sort()
        for item in items:
            paths.extend(recurse(os.path.join(path, item)))
    elif path.endswith(".xml"):
        paths.append(path)
    return paths


resave = bool(sys.argv[1] == "save") if len(sys.argv) > 1 else False


test = False
for path in pyarts.arts.globals.parameters.datapath:
    if path.endswith("arts-xml-data"):
        test = True
        break

if test:
    print("arts-xml-data found in datapath - commenceing test run")
    x = recurse(path)
    v = pyarts.arts.WsvMap(x)
    if resave:
        v.write_split()
    print("All XML files read successfully!")
else:
    print("arts-xml-data not found in datapath - no test run")
