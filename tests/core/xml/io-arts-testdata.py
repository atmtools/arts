import pyarts3 as pyarts
import os
import sys

from xmlload_helper import xml_files

resave = bool(sys.argv[1] == "save") if len(sys.argv) > 1 else False

module = type(pyarts.arts.path)

def recursexml(file, X=pyarts.arts):
    for x in dir(X):
        attr = getattr(X, x)
        if hasattr(attr, "fromxml"):
            try:
                obj = attr.fromxml(file)
                if resave:
                    obj.savexml(file)
                print(f"Read {file} as {x}")
                return True
            except Exception as e:
                pass
        elif isinstance(attr, module):
            if recursexml(file, attr):
                return True
    return False


path = pyarts.arts.globals.data.arts_source_dir
if len(path) != 0:
    print(f"arts source found at {path} - commenceing test run")
    subs = ["tests", "examples", "doc", "python", 'src', 'testdata']
    paths = [os.path.join(path, sub) for sub in subs]
    for p in paths:
        for file in xml_files(p):
            if not recursexml(file):
                print(f"Failed to read {file}")
                assert False
    print("All XML files read successfully!")
