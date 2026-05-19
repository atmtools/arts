import sys

from xmlload_helper import pyarts, recurse

testdata_found = False
for path in pyarts.arts.globals.parameters.datapath:
    if "arts-cat-data" in path:
        print(f"arts-cat-data found at {path} - starting test run")
        recurse(path)
        print("All XML files read successfully!")
        testdata_found = True

if not testdata_found:
    print("arts-cat-data not found")
    sys.exit(1)
