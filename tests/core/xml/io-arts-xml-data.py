import sys

from xmlload_helper import pyarts, xml_files

resave = bool(sys.argv[1] == "save") if len(sys.argv) > 1 else False


test = False
for path in pyarts.arts.globals.parameters.datapath:
    if path.endswith("arts-xml-data"):
        test = True
        break

if test:
    print("arts-xml-data found in datapath - commenceing test run")
    # Preserve coverage and type detection while keeping only one decoded file.
    for file in xml_files(path):
        print(f"Reading {file}", flush=True)
        v = pyarts.arts.WsvMap([file])
        if resave:
            v.write_split()
        del v
    print("All XML files read successfully!")
else:
    print("arts-xml-data not found in datapath - no test run")
