reason = """
We do not want any workspace groups that are not used in any workspace method or variable.

You have three options to fix this issue:
1.  Remove the unused workspace groups.
2.  Add a dependency of the workspace group to a workspace method or agenda.
3.  Move the declaration of the workspace group to the group friend list.
    This is required for one of two reasons:
    1.  The *Group* syntax is used to reference the group in some other
        documentation, such as other workspace groups, variables, methods,
        or agendas.
    2.  The group is used in a way that requires it to be opaque in the
        Python bindings.  This is required for all enum-types, all std::vector,
        and all std::*map types.
"""

import pyarts3 as pyarts

wsvs = pyarts.arts.globals.workspace_variables()
wsgs = pyarts.arts.globals.workspace_groups()
wsms = pyarts.arts.globals.workspace_methods()

wsgs = {}.fromkeys(wsgs, 0)
wsgs["CallbackOperator"] = 1  # Exception since pure user-agenda-method

for wsv in wsvs:
    wsgs[wsvs[wsv].type] += 1

for wsm in wsms:
    for nwsv in wsms[wsm].gin_type:
        if "," in nwsv:
            for sub_nwsv in nwsv.split(","):
                wsgs[sub_nwsv.strip()] += 1
        else:
            wsgs[nwsv] += 1

    for nwsv in wsms[wsm].gout_type:
        if "," in nwsv:
            for sub_nwsv in nwsv.split(","):
                wsgs[sub_nwsv.strip()] += 1
        else:
            wsgs[nwsv] += 1

keys = list(wsgs.keys())
keys.sort()

errors = [f"\n{wsg} should not be a workspace group"
          for wsg in keys if wsgs[wsg] == 0]

assert len(errors) == 0, ", and".join(errors) + '\n' + reason

print("All workspace groups are used in workspace variables or methods")
