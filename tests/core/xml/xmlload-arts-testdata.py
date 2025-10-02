from xmlload_helper import pyarts, recurse

path = pyarts.arts.globals.data.arts_source_dir
if len(path):
    print(f"arts source found at {path} - starting test run")
    recurse(path)
    print("All XML files read successfully!")
else:
    print("arts source dir not found - skipping test run")
