#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 10:59:55 2023

@author: richard
"""

import os
import sys


def recursive_python_files(path):
    out = {}
    files_and_dirs = os.listdir(path)
    for part_path in files_and_dirs:
        part = os.path.abspath(os.path.join(path, part_path))

        if os.path.isdir(part):
            out[part] = recursive_python_files(part)
        elif part.endswith(".py"):
            out[part] = part

    return out


def title_from_path(path, origpath, extra):
    title = path
    for c in origpath:
        if len(title):
            if title[0] == c:
                title = title[1:]
    for c in "yp.":
        if len(title):
            if title[-1] == c:
                title = title[:-1]
    title = title.replace("/", ".").replace("\\", ".")
    if len(title):
        while title[0] == '.':
            title = title[1:]
        return extra + "." + title
    else:
        return extra  # For empty, i.e., main folder


def rstfn_from_title(title, outpath):
    return os.path.join(outpath, title + ".rst")


def readme_path(path):
    if os.path.isdir(path):
        readme = os.path.join(path, "README.txt")
    else:
        d, f = os.path.split(path)
        readme = os.path.join(d, "README." + f.rstrip(".py") + ".txt")

    if os.path.exists(readme):
        print(f"Found readme: {readme}")
        return readme
    print(f"Did not find readme: {readme}")
    return None


def print_pyfile(path, origpath, extra, outpath):
    title = title_from_path(path, origpath, extra)
    rstfn = rstfn_from_title(title, outpath)
    readme = readme_path(path)

    print(f"Creating {rstfn} from {path}")
    with open(rstfn, "w") as rstfile:
        rstfile.write(f"{title}\n{'=' * len(title)}\n")

        if readme:
            with open(readme, "r") as r:
                rstfile.write(r.read())
                rstfile.write("\n")

        rstfile.write(".. code-block:: python\n")
        rstfile.write("    :linenos:\n\n")
        with open(path, "r") as pyfile:
            for line in pyfile.read().split('\n'):
                rstfile.write(f"    {line}\n")


def print_folder(path, paths, origpath, extra, outpath):
    keys = paths.keys()
    title = title_from_path(path, origpath, extra)
    rstfn = rstfn_from_title(title, outpath)
    readme = readme_path(path)

    print(f"Creating {rstfn} from {path}")
    with open(rstfn, "w") as rstfile:
        rstfile.write(f"{title}\n{'=' * len(title)}\n")

        if readme:
            with open(readme, "r") as r:
                rstfile.write(r.read())
                rstfile.write("\n")

        rstfile.write(".. toctree::\n")
        for key in keys:
            rstfile.write(f"    {title_from_path(key, origpath, extra)}\n")


def folders(paths, origpath, extra, outpath):
    keys = paths.keys()
    for key in keys:
        if os.path.isdir(key):
            folders(paths[key], origpath, extra, outpath)
            print_folder(key, paths[key], origpath, extra, outpath)
        else:
            print_pyfile(paths[key], origpath, extra, outpath)


if __name__ == "__main__":
    extra, outpat = sys.argv[2], sys.argv[3]
    origpath = os.path.abspath(sys.argv[1])
    print(f"Scanning and printing {origpath}")
    paths = recursive_python_files(origpath)
    folders(paths, origpath, extra, outpat)
    print_folder(origpath, paths, origpath, extra, outpat)
