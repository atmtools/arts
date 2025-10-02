from xmlload_helper import pyarts, recurse

for path in pyarts.arts.globals.parameters.datapath:
    if "arts-cat-data" in path:
        print(f"arts-cat-data found at {path} - starting test run")
        recurse(path)
        print("All XML files read successfully!")
else:
    print("arts-cat-data not found - skipping test run")
