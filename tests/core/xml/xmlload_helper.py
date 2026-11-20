import os

import pyarts3 as pyarts


def xml_files(path):
    """All XML files below path, sorted, skipping hidden, 3rdparty and build items"""
    if os.path.isdir(path) and not os.path.basename(path).startswith("."):
        for item in sorted(os.listdir(path)):
            if item == "3rdparty" or "build" in item:
                continue
            yield from xml_files(os.path.join(path, item))
    elif path.endswith(".xml"):
        yield path


def recursexml(file):
    print(f"Loading {file}")
    try:
        pyarts.xml.load(file)
    except TypeError:
        print(f"Ignoring Map type in {file}")
    except Exception as e:
        raise RuntimeError(f"Failed to xml.load {file}:\n {e}")


def recurse(path):
    for file in xml_files(path):
        recursexml(file)
