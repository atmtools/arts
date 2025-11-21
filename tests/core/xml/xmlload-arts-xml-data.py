from xmlload_helper import pyarts, recurse

for path in pyarts.arts.globals.parameters.datapath:
    if "arts-xml-data" in path:
        print(f"arts-xml-data found at {path} - starting test run")
        recurse(path)
        print("All XML files read successfully!")
    else:
        print("arts-xml-data not found - skipping test run")
