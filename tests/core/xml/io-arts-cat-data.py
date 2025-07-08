import pyarts
import os
import sys

resave = bool(sys.argv[1] == "save") if len(sys.argv) > 1 else False

test = False
for path in pyarts.arts.globals.parameters.datapath:
    if path.endswith("arts-cat-data"):
        test = True
        break

if test:
    print("arts-cat-data found in datapath - commenceing test run")

    cia = os.path.join(path, "cia")
    for file in os.listdir(cia):
        if file.endswith(".xml"):
            filepath = os.path.join(cia, file)
            print(f"Reading {filepath}")
            x = pyarts.arts.CIARecord.fromxml(filepath)
            if resave:
                x.savexml(filepath)
    print("All CIA files read successfully.")
    print()

    xsec = os.path.join(path, "xsec")
    for file in os.listdir(xsec):
        if file.endswith(".xml"):
            filepath = os.path.join(xsec, file)
            print(f"Reading {filepath}")
            x = pyarts.arts.XsecRecord.fromxml(filepath)
            if resave:
                x.savexml(filepath)
    print("All xsec files read successfully.")
    print()

    predef = os.path.join(path, "predef")
    for file in os.listdir(predef):
        if file.endswith(".xml"):
            filepath = os.path.join(predef, file)
            print(f"Reading {filepath}")
            x = pyarts.arts.PredefinedModelData.fromxml(filepath)
            if resave:
                x.savexml(filepath)
    print("All predef files read successfully.")
    print()

    lines = os.path.join(path, "lines")
    for file in os.listdir(lines):
        if file.endswith(".xml"):
            filepath = os.path.join(lines, file)
            print(f"Reading {filepath}")
            x = pyarts.arts.AbsorptionBands.fromxml(filepath)
            if resave:
                x.savexml(filepath)
    print("All line files read successfully.")
    print()

else:
    print("arts-cat-data not found in datapath - no test run")
