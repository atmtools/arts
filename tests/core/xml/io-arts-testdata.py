import pyarts3 as pyarts
import os
import sys

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


def recurse(path):
    if os.path.isdir(path) and not os.path.basename(path).startswith("."):
        items = os.listdir(path)
        items.sort()
        for item in items:
            if item == "3rdparty":
                continue
            if "build" in item:
                continue
            recurse(os.path.join(path, item))
    elif path.endswith(".xml"):
        if not recursexml(path):
            print(f"Failed to read {path}")
            assert False


path = pyarts.arts.globals.data.arts_source_dir
if len(path) != 0:
    print(f"arts source found at {path} - commenceing test run")
    recurse(path)
    print("All XML files read successfully!")
