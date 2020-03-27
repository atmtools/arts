import os


def dict2lfd(d):
    if "HTP" in d:
        key = "HTP"
    elif "LP" in d:
        key = "LP"
    elif "VP" in d:
        key = "VP"
    else:
        key = "DP"

    if "LM2" in d:
        lm = "LM2"
    else:
        lm = "#"

    start = "LF {} {} {} ".format(key, lm, "0" if key == "DP" else "1 AIR")

    if key == "HTP":
        values = [
            d[key].get("G0", None),
            d[key].get("D0", None),
            d[key].get("G2", None),
            d[key].get("D2", None),
            d[key].get("FVC", None),
            d[key].get("ETA", None),
        ]
    elif key == "VP" or key == "LP":
        values = [d[key].get("G0", None), d[key].get("D0", None)]
    else:
        values = []

    if lm == "LM2":
        values.append(d[lm].get("Y", None))
        values.append(d[lm].get("G", None))
        values.append(d[lm].get("DV", None))

    end = ""
    for x in values:
        if x is None:
            end += "T0 0 "
        else:
            for y in x:
                end += "{} ".format(y)
    return start + end


x1 = '<?xml version="1.0"?>'
x2 = '<arts format="ascii" version="1">'
x3 = '<ArrayOfLineRecord version="ARTSCAT-5" nelem="1">'
x4 = "@ O2-66 100000000000 1e-15 296 3e-20 0 3 1 "
x6 = " QN UP J 1 LO J 0"
x7 = "</ArrayOfLineRecord>"
x8 = "</arts>"

LM_HTP_params = {
    "HTP": {
        "G0": ("T1", 10000, 0.7),
        "D0": ("T5", 1000, 0.7),
        "G2": ("T3", 150, 0.1),
        "D2": ("T2", 100, 1, 0.1),
        "FVC": ("T3", 800, 250),
        "ETA": ("T3", 0.5, 0.001),
    },
    "LM2": {
        "Y": ("T4", 1e-7, 1e-9, 0.8),
        "G": ("T4", 1e-11, 1e-13, 1.6),
        "DV": ("T4", 1e-4, 1e-6, 1.6),
    },
}
HTP_params = {"HTP": LM_HTP_params["HTP"]}
HTP_SDVP_params = {
    "HTP": {
        "G0": HTP_params["HTP"]["G0"],
        "D0": HTP_params["HTP"]["D0"],
        "G2": HTP_params["HTP"]["G2"],
        "D2": HTP_params["HTP"]["D2"],
    }
}
HTP_VP_params = {"HTP": {"G0": HTP_params["HTP"]["G0"], "D0": HTP_params["HTP"]["D0"]}}
VP_params = {"VP": {"G0": HTP_params["HTP"]["G0"], "D0": HTP_params["HTP"]["D0"]}}
LP_params = {"LP": {"G0": HTP_params["HTP"]["G0"], "D0": HTP_params["HTP"]["D0"]}}
LM_VP_params = {
    "VP": {"G0": HTP_params["HTP"]["G0"], "D0": HTP_params["HTP"]["D0"]},
    "LM2": LM_HTP_params["LM2"],
}

LM_LP_params = {
    "LP": {"G0": HTP_params["HTP"]["G0"], "D0": HTP_params["HTP"]["D0"]},
    "LM2": LM_HTP_params["LM2"],
}
DP_params = {}

lsd = [
    LM_HTP_params,
    HTP_params,
    HTP_SDVP_params,
    HTP_VP_params,
    LM_VP_params,
    VP_params,
    LM_LP_params,
    LP_params,
    DP_params,
]
lsn = ["lm-htp", "htp", "htp-sdvp", "htp-vp", "lm-vp", "vp", "lm-lp", "lp", "dp"]
for i in range(len(lsd)):
    f = open("{}-line.xml".format(lsn[i]), "w")
    f.write(x1)
    f.write("\n")
    f.write(x2)
    f.write("\n")
    f.write(x3)
    f.write("\n")
    f.write(x4)
    f.write(dict2lfd(lsd[i]))
    f.write(x6)
    f.write("\n")
    f.write(x7)
    f.write("\n")
    f.write(x8)
    f.write("\n")
    f.close()
    try:
        os.mkdir("test-{}".format(lsn[i]))
    except FileExistsError:
        print("test-{}".format(lsn[i]), "already exists, will ignore mkdir")
