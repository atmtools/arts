#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 10:59:55 2023

@author: richard
"""

import os
import re
import sys


INTROFILE = "README.rst"


def title_to_heading(title):
    pos = title.rfind(".")
    if pos == -1:
        ret = title[0].upper() + title[1:]
    else:
        ret = title[pos + 1].upper() + title[pos + 2 :]
    ret = re.sub(r"^(\d+)-", r"", ret)
    m = re.match(r"^\d+\. ", ret)
    if m:
        pos = len(m.group(0))
        ret = ret[:pos] + ret[pos].upper() + ret[pos + 1 :]
    else:
        ret = ret[0].upper() + ret[1 :]
    return ret.replace("-", " ")


def all_files(path):
    """Returns a dict of all files in the path as a list and all folders as a dict"""
    out = {path: []}

    for file in os.listdir(path):
        full_path = os.path.join(path, file)
        if os.path.isdir(full_path):
            out[file] = all_files(full_path)
        else:
            out[path].append(file)

    return out


def filename_from_path(path, prefix, save_path, fileending, suffix):
    path = path.removeprefix(prefix).removesuffix(suffix) + fileending
    path = path.replace(os.path.sep, ".")
    return os.path.join(save_path, path)


class notebook:
    def __init__(self, str):
        self.str = str


class rst:
    def __init__(self, str):
        self.str = str


def pyrst(str):
    if str.endswith(".py"):
        stem = str.rstrip(".py")
    elif str.endswith(".rst"):
        stem = str.rstrip(".rst")
    return stem + ".py", stem + ".rst"


def combine_rstpy(rsts, pyfiles, out=None):
    out = []

    if INTROFILE in rsts:
        rsts.remove(INTROFILE)
        out.append([INTROFILE, None])

    while len(rsts) != 0 or len(pyfiles) != 0:
        if len(rsts) == 0:
            file = pyfiles[0]
        elif len(pyfiles) == 0:
            file = rsts[0]
        else:
            file = rsts[0] if rsts[0] < pyfiles[0] else pyfiles[0]

        py, rst = pyrst(file)

        haspy = pyfiles.count(py) == 1
        hasrst = rsts.count(rst) == 1
        out.append([rst if hasrst else None, py if haspy else None])

        if haspy:
            pyfiles.remove(py)
        if hasrst:
            rsts.remove(rst)

    return out


def from_list(lst, path):
    """Converts a list of paths to rst string"""

    assert isinstance(lst, list), "Expected a list of paths"

    # Gather all RSTs
    rsts = [item for item in lst if item.endswith(".rst")]
    rsts.sort()

    # Gather all python files
    pyfiles = [item for item in lst if item.endswith(".py")]
    pyfiles.sort()

    # Gather all notebooks
    notebooks = [item for item in lst if item.endswith(".ipynb")]
    notebooks.sort()

    # E.g., if it is a test-data folder, we can have no files
    if len(rsts) == 0 and len(pyfiles) == 0 and len(notebooks) == 0:
        return None

    assert len(notebooks) <= 1, "Max 1 ipynb file per folder"
    assert len(rsts) * len(notebooks) == 0, "Choose rst or ipynb"
    assert len(pyfiles) * len(notebooks) == 0, "Choose py or ipynb"
    if len(notebooks):
        with open(os.path.join(path, notebooks[0]), "r") as nbfile:
            return notebook(nbfile.read())

    # Combine RST and Python files
    files = combine_rstpy(rsts, pyfiles)

    heading = title_to_heading(os.path.basename(os.path.normpath(path)))
    out = f"{heading}\n{'=' * len(heading)}\n\n"

    for item in files:
        rstf, py = item

        if len(files) > 1 and rstf != INTROFILE:
            heading = title_to_heading(
                rstf.rstrip(".rst") if rstf else py.rstrip(".py")
            )
            out += f"{heading}\n{'-' * len(heading)}\n\n"

        if rstf:
            with open(os.path.join(path, rstf), "r") as f:
                out += f.read() + "\n"
            if rstf == INTROFILE:
                out += """.. contents::
   :depth: 1
   :local:
   :backlinks: none

"""

        if py:
            out += f".. code-block:: python\n"
            out += "    :linenos:\n\n"
            with open(os.path.join(path, py), "r") as pyfile:
                for line in pyfile.read().split("\n"):
                    out += f"    {line}\n"

    return rst(out)


def generate_texts(paths):
    out = {}

    for path in paths:
        assert path not in out, f"{path} in {list(out.keys())}"

        data = paths[path]

        if isinstance(data, dict):
            x = generate_texts(data)
            if len(x) != 0:
                out[path] = x
        elif isinstance(data, list):
            x = from_list(data, path)

            if x is not None:
                out[path] = x
        else:
            assert False, "Unexpected data type in paths: " + str(type(data))

    return out


def generate_filetrees(paths, corepath):
    out = {corepath: []}
    for path in paths:
        data = paths[path]

        if isinstance(data, dict):
            out[corepath].append(path)
            x = generate_filetrees(data, os.path.join(corepath, path))
            if x is not None:
                out[path] = x

    if len(out) == 1 and len(out[corepath]) == 0:
        return None
    out[corepath].sort()
    return out


def filetrees_to_toctrees(filetree, arts_path):
    out = {}

    for path in filetree:
        assert path not in out, f"{path} in {list(out.keys())}"
        data = filetree[path]

        if isinstance(data, dict):
            out[path] = filetrees_to_toctrees(data, arts_path)
        elif isinstance(data, list):
            text =  ".. toctree::\n"
            text += "   :maxdepth: 2\n\n"
            for item in data:
                text += (
                    f"   {filename_from_path(path+"."+item, arts_path, "", "", "")}\n"
                )
            out[path] = rst(text)

    return out


def flatten_map(d):
    out = {}
    for key in d:
        assert key not in out, f"Duplicate key {key} in flatten_map"

        data = d[key]
        if isinstance(data, dict):
            v = flatten_map(data)
            for k in v:
                assert k not in out, f"Duplicate key {k} in flatten_map"
                out[k] = v[k]
        else:
            out[key] = data

    return out


def generate_docfiles(flat_toc, flat_txt, arts_path, save_path):
    for key in flat_txt:
        data = flat_txt[key]

        if isinstance(data, notebook):
            with open(
                filename_from_path(key, arts_path, save_path, ".ipynb", ".py"), "w"
            ) as f:
                f.write(data.str)
        elif isinstance(data, rst):
            extr = rst("")

            if key in flat_toc:
                extr = flat_toc[key]
                del flat_toc[key]

            test = data.str + "\n" + extr.str + "\n"
            with open(
                filename_from_path(key, arts_path, save_path, ".rst", ".py"), "w"
            ) as f:
                f.write(test)
        else:
            assert False, "Unexpected data type in flat_text: " + str(type(data))

    for key in flat_toc:
        data = flat_toc[key]

        assert isinstance(data, rst), "Expected rst type in flat_toc, got " + str(
            type(data)
        )

        with open(
            filename_from_path(key, arts_path, save_path, ".rst", ".py"), "w"
        ) as f:
            f.write(data.str)


if __name__ == "__main__":
    arts_path = sys.argv[1]
    save_path = sys.argv[2]
    examplepath = os.path.join(arts_path, "examples")

    files = all_files(examplepath)
    text = generate_texts(files)
    filetrees = generate_filetrees(files, examplepath)
    toctrees = filetrees_to_toctrees(filetrees, arts_path)

    generate_docfiles(flatten_map(toctrees), flatten_map(text), arts_path, save_path)
