import os

import pyarts3 as pyarts


def recursexml(file, X=pyarts.arts):
    print(f"Loading {file}")
    try:
        pyarts.xml.load(file)
    except TypeError:
        print(f"Ignoring Map type in {file}")
    except Exception as e:
        raise RuntimeError(f"Failed to xml.load {file}:\n {e}")


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
        recursexml(path)
