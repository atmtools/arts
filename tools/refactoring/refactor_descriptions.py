#!/usr/bin/env python3

"""
Converts old string descriptions to string literals.
Currently works for methods.cc and agendas.cc.
"""

import re
import sys

with open(sys.argv[1], "r") as f:
    lines = f.readlines()

with open(sys.argv[1] + ".new", "w") as f:
    inside = False  # inside DESCRIPTION block
    inside_r = False  # inside R"() block inside DESCRIPTION
    i = 0
    numlines = len(lines)
    while i < numlines:
        line = lines[i]
        if re.match(r' *//', line):
            f.write(line)
            i += 1
            continue
        if i < numlines - 1:
            nextline = lines[i + 1]
            while re.match(r' *//', line):
                f.write(line)
                i += 1
                nextline = lines[i + 1]
            end_description = re.search(r"(GROUP|OUTPUT|AUTHORS)\(", nextline)
        else:
            end_description = False
        if inside:
            s = line[:-1]
            if re.match(r' *R"\(', s) or re.match(r' *R"--\(', s):
                inside_r = True
                f.write('\n')
                i += 1
                continue
            if inside_r:
                if s == ')"),' or s == ')--"),':
                    inside = False
                    inside_r = False
                    f.write(')--"),\n')
                elif s == ')"':
                    inside_r = False
                    f.write('\n')
                else:
                    f.write(s + '\n')
                i += 1
                continue
            # Ignore empty lines inside descriptions
            if s == "":
                i += 1
                continue
            s = re.sub('//.*$', '', s)
            s = re.sub('^[^"]*"', '', s)
            s = re.sub('" *$', '', s)
            s = s.replace('\\n\\n', '\n')
            s = s.replace('\\n', '')
            if end_description:
                if re.match(r' *\), *', s):
                    s = re.sub(r'[ "]*\) *, *$', ')--"),', s)
                else:
                    s = re.sub(r'[ "]*\) *, *$', '\n)--"),', s)
                inside = False
            f.write(s + '\n')
        else:
            m = re.search(r"(^ *)(DESCRIPTION\()(.*)", line)
            # Match DESCRIPTION if it doesn't already contain a literal string
            if ((m and len(m.group(3)) == 0)
                    or (m and len(m.group(3)) and m.group(3)[0] != 'R')):
                # Ignore descriptions that are not plain strings or
                # already literal strings
                if (re.search(r"^ *String\(", lines[i + 1])
                        or re.search(r'^ *R"(--)?\(', lines[i + 1])):
                    f.write(line)
                    i += 1
                    continue
                inside = True
                # Does the description start on the same line?
                if len(m.group(3)) and m.group(3)[0] == '"':
                    s = m.group(3)
                else:
                    i += 1
                    s = lines[i][:-1]
                s = re.sub('^ *"', 'R"--(', s)
                s = re.sub('" *$', '', s)
                s = s.replace('\\n\\n', '\n')
                s = s.replace('\\n', '')
                # Is this the last line of the description?
                if i < numlines - 1:
                    nextline = lines[i + 1]
                    end_description = re.search(r"(GROUP|OUTPUT|AUTHORS)\(", nextline)
                else:
                    end_description = False
                if end_description:
                    s = re.sub(r'[ "]*\) *, *$', '\n)--"),', s)
                    inside = False
                    end_description = False
                f.write(m.group(1) + "DESCRIPTION(" + s + '\n')
            else:
                f.write(line)
        i += 1
