import pyarts3 as pyarts
from warnings import warn


def matches(x):
    wsvs = list(pyarts.arts.globals.workspace_variables().keys())

    ret = []
    for v in wsvs:
        if x == v:
            ret.append(v)
        elif v.startswith(x + "_"):
            ret.append(v)
        elif v.endswith("_" + x):
            ret.append(v)
        elif ("_" + x + "_") in v:
            ret.append(v)
    return sorted(ret)


def collect():
    shorts = pyarts.arts.globals.workspace_variables_shortnames()
    res = {}

    for name in shorts:
        shortform = shorts[name]
        title = shortform.title
        if title not in res:
            res[title] = {}

        res[title][name] = (shortform.desc, matches(name))

    return res


def create_rst():
    shortforms = collect()
    rst = """
Short-forms for naming workspace variables
##########################################

This is a mostly auto-generated list of short-forms for naming workspace
variables in ARTS.

All variables in ARTS use short-hand names to identify them and
simplify writing and using them. The keys to read these short-hand
names are listed below.

All complete short-hand names are separated by an underscore (``_``), or
equals the name itself.
"""

    all_wsvs = []

    main_keys = sorted(list(shortforms.keys()))
    for title in main_keys:
        rst += f"\n{title}\n{'=' * len(title)}\n\n"

        secondary_keys = sorted(list(shortforms[title].keys()))
        for name in secondary_keys:
            desc, wsvs = shortforms[title][name]

            rst += f"\n``{name}``\n{'-' * (len(name) + 4)}\n\n{desc}\n\n"
            all_wsvs.extend(wsvs)

            if len(wsvs) > 0:
                rst += """.. hlist::
    :columns: 2

"""
                for w in wsvs:
                    rst += f"    * :attr:`~pyarts3.workspace.Workspace.{w}`\n"
                rst += "\n"
            else:
                msg = Warning(
                    f"No workspace variables found documenting use of shortform '{name}'.")
                warn(msg, stacklevel=3)

    rst += ".. rubric:: Unshortened Variables\n\n"

    all_wsvs = set(list(pyarts.arts.globals.workspace_variables().keys())
                   ).symmetric_difference(set(all_wsvs))
    rst += """.. hlist::
    :columns: 2

"""
    for w in all_wsvs:
        rst += f"    * :attr:`~pyarts3.workspace.Workspace.{w}`\n"

    return rst


print(create_rst())
