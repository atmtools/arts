import pyarts

wsvs = pyarts.arts.globals.workspace_variables()
wsgs = pyarts.arts.globals.workspace_groups()
wsms = pyarts.arts.globals.workspace_methods()

wsgs = {}.fromkeys(wsgs, 0)
wsgs["CallbackOperator"] = 1

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

for wsg in keys:
    if wsgs[wsg] == 0:
        print(f"{wsg} is not used in any workspace method and is not a variable")
