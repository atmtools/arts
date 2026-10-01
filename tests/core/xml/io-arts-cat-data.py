import pyarts3 as pyarts
import os
import sys

resave = bool(sys.argv[1] == "save") if len(sys.argv) > 1 else False

test = False
for path in pyarts.arts.globals.parameters.datapath:
    if "arts-cat-data" in path:
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
    for file in cia_files:
        cia = pyarts.arts.CIARecord.fromxml(file)
        if resave:
            cia.savexml(file, type='binary')
        del cia
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
        del xsec
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
            del x
    print("All predef files read successfully.")
    print()

    lines = os.path.join(path, "lines")
    # Validate every isotope file without retaining the full line catalogue.
    # The complete catalogue can expand to many gigabytes in memory.
    for file in sorted(os.listdir(lines)):
        filepath = os.path.join(lines, file)
        if not file.endswith(".xml") or not os.path.isfile(filepath):
            continue
        print(f"Reading {filepath}", flush=True)
        bands = pyarts.arts.AbsorptionBands.fromxml(filepath)
        if resave:
            bands.savexml(filepath)
        del bands
    print("All line files read successfully.")
    print()

else:
    print("arts-cat-data not found in datapath - no test run")
    sys.exit(1)
