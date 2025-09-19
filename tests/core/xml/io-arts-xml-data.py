import pyarts3 as pyarts
import os
import sys

resave = bool(sys.argv[1] == "save") if len(sys.argv) > 1 else False


def recursexml(file):
    for x in dir(pyarts.arts):
        attr = getattr(pyarts.arts, x)
        if hasattr(attr, "fromxml"):
            try:
                obj = attr.fromxml(file)
                global resave
                if resave:
                    obj.savexml(file)
                print(f"Read {file} as {x}")
                return True
            except Exception as e:
                pass
    return False


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


test = False
for path in pyarts.arts.globals.parameters.datapath:
    if path.endswith("arts-xml-data"):
        test = True
        break

if test:
    print("arts-xml-data found in datapath - commenceing test run")
    x = recurse(path)
    v = pyarts.arts.WsvMap(x)
    print("All XML files read successfully!")
else:
    print("arts-xml-data not found in datapath - no test run")
