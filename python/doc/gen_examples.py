#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 10:59:55 2023

@author: richard
"""

import os
import re
import sys


MAXDEPTH = 2

PRIORITIZE = ["getting-started"]


def priority_sort(paths):
    paths.sort()
    basenames = [os.path.basename(p) for p in paths]
    for prio in PRIORITIZE[::-1]:
        try:
            i = basenames.index(prio)
        except ValueError:
            pass
        else:
            basenames.insert(0, basenames.pop(i))
            paths.insert(0, paths.pop(i))


def recursive_python_files(path):
    out = {}
    files_and_dirs = os.listdir(path)
    priority_sort(files_and_dirs)
    for part_path in files_and_dirs:
        part = os.path.abspath(os.path.join(path, part_path))

        if os.path.isdir(part):
            out[part] = recursive_python_files(part)
        elif part.endswith(".py"):
            out[part] = part

    return out


def title_to_heading(title):
    pos = title.rfind(".")
    if pos == -1:
        ret = title[0].upper() + title[1:]
    else:
        ret = title[pos + 1].upper() + title[pos + 2:]
    ret = re.sub(r"^(\d+)-", r"\1. ", ret)
    m = re.match(r"^\d+\. ", ret)
    if m:
        pos = len(m.group(0))
        ret = ret[:pos] + ret[pos].upper() + ret[pos + 1:]
    return ret


def title_from_path(path, origpath, extra=None):
    title = path.removeprefix(origpath).removesuffix(".py")
    title = title.replace("/", ".").replace("\\", ".")
    if len(title):
        title = title.lstrip(".")
        return extra + "." + title
    else:
        return extra  # For empty, i.e., main folder


def rstfn_from_title(title, outpath):
    return os.path.join(outpath, title + ".rst")


def readme_path(path):
    if os.path.isdir(path):
        readme = os.path.join(path, "README.rst")
    else:
        d, f = os.path.split(path)
        readme = os.path.join(d, "README." + f.removesuffix(".py") + ".rst")

    if os.path.exists(readme):
        print(f"Found readme: {readme}")
        return readme
    print(f"Did not find readme: {readme}")
    return None


def output_readme(rstfile, path):
    readme = readme_path(path)
    if readme:
        with open(readme, "r") as r:
            rstfile.write(r.read())
            rstfile.write("\n")


def combine_pyfiles(paths, origpath, extra, outpath):
    parentpath = os.path.dirname(paths[0])
    title = title_from_path(parentpath, origpath, extra)
    rstfn = rstfn_from_title(title, outpath)
    print(f"Creating {rstfn} from {parentpath}")
    with open(rstfn, "w") as rstfile:
        heading = title_to_heading(title)
        rstfile.write(f"{heading}\n{'=' * len(heading)}\n")
        output_readme(rstfile, parentpath)
        for path in paths:
            title = title_from_path(path, origpath, extra)

            if len(paths) > 1:
                heading = title_to_heading(title)
                rstfile.write(f"{heading}\n{'-' * len(heading)}\n")

            output_readme(rstfile, path)

            rstfile.write(".. code-block:: python\n")
            rstfile.write(f"    :name: {heading}\n")
            rstfile.write("    :linenos:\n\n")
            with open(path, "r") as pyfile:
                for line in pyfile.read().split('\n'):
                    rstfile.write(f"    {line}\n")


def print_pyfile(path, origpath, extra, outpath):
    title = title_from_path(path, origpath, extra)
    rstfn = rstfn_from_title(title, outpath)

    print(f"Creating {rstfn} from {path}")
    with open(rstfn, "w") as rstfile:
        heading = title_to_heading(title)
        rstfile.write(f"{heading}\n{'=' * len(heading)}\n")

        output_readme(rstfile, path)

        rstfile.write(".. code-block:: python\n")
        rstfile.write("    :linenos:\n\n")
        with open(path, "r") as pyfile:
            for line in pyfile.read().split('\n'):
                rstfile.write(f"    {line}\n")


def print_folder(path, paths, origpath, extra, outpath):
    keys = paths.keys()
    title = title_from_path(path, origpath, extra)
    rstfn = rstfn_from_title(title, outpath)

    print(f"Creating {rstfn} from {path}")
    with open(rstfn, "w") as rstfile:
        heading = title_to_heading(title)
        rstfile.write(f"{heading}\n{'=' * len(heading)}\n")

        output_readme(rstfile, path)

        rstfile.write(".. toctree::\n")
        for key in keys:
            rstfile.write(f"    {title_from_path(key, origpath, extra)}\n")


def folders(paths, origpath, extra, outpath, depth=0):
    keys = paths.keys()
    if depth < MAXDEPTH:
        for key in keys:
            if os.path.isdir(key):
                print(f"Processing {key}")
                print(f"FolderDepth {depth}")
                folders(paths[key], origpath, extra, outpath, depth+1)
                if depth+1 < MAXDEPTH:
                    print_folder(key, paths[key], origpath, extra, outpath)
            else:
                print_pyfile(paths[key], origpath, extra, outpath)
    else:
        combine_files = []
        for key in keys:
            if os.path.isdir(key):
                print(f"Ignoring folder {key}")
            else:
                combine_files.append(paths[key])
        print(f"Combining {combine_files}")
        combine_pyfiles(combine_files, origpath, extra, outpath)


if __name__ == "__main__":
    extra, outpat = sys.argv[2], sys.argv[3]
    origpath = os.path.abspath(sys.argv[1])
    print(f"Scanning and printing {origpath}")
    paths = recursive_python_files(origpath)
    folders(paths, origpath, extra, outpat)
    print_folder(origpath, paths, origpath, extra, outpat)
