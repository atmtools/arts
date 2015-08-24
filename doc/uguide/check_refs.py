#!/usr/bin/env python
r"""
Checks that WSVs and WSMs used in tex macros \wsmindex{}, \wsmindex{},
\wsaindex{} and \builtindoc{} really exist.

__author__ = "Oliver Lemke"
"""

import os
import sys
import glob
import re
import subprocess as sp

if sys.version_info < (3,):
    from io import open


class ArtsWorkspace:
    def __init__(self, bindir):
        self.arts_exe = bindir + "/src/arts"
        self.wsms = self.parse_arts_list(["-m", "all"])
        self.wsvs = self.parse_arts_list(["-w", "all"])
        self.groups = self.parse_arts_list(["-g", "all"])

    def parse_arts_list(self, arts_arg):
        return [l[2:].decode()
                for l in
                sp.Popen([self.arts_exe] + arts_arg,
                         stdout=sp.PIPE).communicate()[0].splitlines()
                if l[0:2] == b'- ']


def check_tex_file(ws, filename):
    l = 0
    rebuiltin = re.compile(r"\\builtindoc{(.*?)}")
    rewsm = re.compile(r"\\wsmindex{(.*?)}")
    rewsv = re.compile(r"\\wsvindex{(.*?)}")
    rewag = re.compile(r"\\wsaindex{(.*?)}")
    wsms = []
    wsvs = []
    wags = []
    builtins = []

    f = open(filename, 'r', encoding='utf-8')
    for line in f:
        l += 1
        binmatches = rebuiltin.findall(line)
        for m in binmatches:
            builtins.append((l, m.replace("\\", "")))

        wsmmatches = rewsm.findall(line)
        for m in wsmmatches:
            wsms.append((l, m.replace("\\", "")))

        wsvmatches = rewsv.findall(line)
        for m in wsvmatches:
            wsvs.append((l, m.replace("\\", "")))

        wagmatches = rewag.findall(line)
        for m in wagmatches:
            wags.append((l, m.replace("\\", "")))

    f.close()

    err = []
    err.extend([(l, "WSM", m) for l, m in wsms if m not in ws.wsms])
    err.extend([(l, "WSV", m) for l, m in wsvs if m not in ws.wsvs])
    err.extend([(l, "agenda", m) for l, m in wags if m not in ws.wsvs])
    err.extend([(l, "workspace entity", m)
                for l, m in builtins
                if m not in ws.wsms
                and m not in ws.wsvs and m not in ws.groups])

    if err:
        sys.stderr.write("%s:\n" % (filename))
        for e in err:
            sys.stderr.write("l. %d: Unknown %s %s\n" % e)
        sys.stderr.write("\n")

    return len(err) == 0


def check_arts_tex_files():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s CMAKE_BINARY_DIR CMAKE_CURRENT_SOURCE_DIR"
                         % (sys.argv[0]))
        sys.exit(1)

    bindir = sys.argv[1]
    srcdir = sys.argv[2]

    ws = ArtsWorkspace(bindir)

    fail = False
    texdir = srcdir + "/*.tex"
    texfiles = glob.glob(texdir)
    if texfiles:
        for filename in texfiles:
            if os.path.basename(filename) != "common.tex":
                if not check_tex_file(ws, filename):
                    fail = True
    else:
        sys.stderr.write("Error: No tex files found in %s!\n" % (texdir))

    if fail:
        sys.exit(1)


if __name__ == "__main__":
    check_arts_tex_files()
