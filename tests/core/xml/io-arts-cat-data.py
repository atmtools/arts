import pyarts3 as pyarts
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

    cia_files = []
    cia = os.path.join(path, "cia")
    for file in os.listdir(cia):
        if file.endswith(".xml"):
            filepath = os.path.join(cia, file)
            print(f"Found {filepath}")
            cia_files.append(filepath)
    cia = pyarts.arts.ArrayOfCIARecord.fromxmls(cia_files)
    if resave:
        for i, file in enumerate(cia_files):
            cia[i].savexml(file, type='binary')
    print("All CIA files read successfully.")
    print()

    xsec_files = []
    xsec = os.path.join(path, "xsec")
    for file in os.listdir(xsec):
        if file.endswith(".xml"):
            filepath = os.path.join(xsec, file)
            print(f"Found {filepath}")
            xsec_files.append(filepath)
    for file in xsec_files:
        xsec = pyarts.arts.XsecRecord.fromxml(file)
        if resave:
            xsec.savexml(file, type='binary')
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
    ws = pyarts.Workspace()
    ws.abs_bandsReadSplit(dir=lines)
    print("All line files read successfully.")
    print()

else:
    print("arts-cat-data not found in datapath - no test run")
