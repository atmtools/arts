from urllib.request import urlopen

# Map Hitran to ARTS species names
_HITRAN_TO_ARTS_NAMES = {
    "CH3CN-2124": "CH3CN-211124",
    "CO2-827": "CO2-728",
    "H2CO-126": "H2CO-1126",
    "H2CO-128": "H2CO-1128",
    "H2CO-136": "H2CO-1136",
    "HC3N-1224": "HC3N-12224",
    "HCOOH-126": "HCOOH-1261",
}


def gen_latest_molparam_map(molparam_txt_file=None):
    """ Generates a version of latest_molparam_map used in hitran_species.cc

    The variable is simply printed to stream with print() as the intent is
    to use this output in ARTS directly.  ARTS does not use AFGL notation
    internally, but species names that are different from Hitran are mapped to
    ARTS names in the output.
    """
    def pos2char(ind):
        """ Convert an isotoplogue index to a HITRAN char for that index """
        if ind < 10:
            return "'{}'".format(ind)
        elif ind == 10:
            return "'0'"
        elif ind == 11:
            return "'A'"
        elif ind == 12:
            return "'B'"

    if molparam_txt_file:
        molparam_txt = open(molparam_txt_file, 'r').read().split('\n')
    else:
        molparam_txt = urlopen("https://hitran.org/media/molparam.txt").read(
        ).decode("utf-8").split('\n')
    molparam_txt.pop(0)  # erase header

    out = {}
    for i in range(len(molparam_txt)):
        line = molparam_txt[i]
        vec = line.split()
        if len(vec) == 0:
            continue
        elif len(vec) == 2:
            spec = vec[0]
            specnum = int(eval(vec[1]))
            pos = 1
            out [spec] = []
        else:
            out[spec].append([specnum, pos, vec[0]])
            pos += 1

    print ('static const std::map<Index, std::map<char, const char * const>> latest_molparam_map{')
    for spec in out:
        print ('{',out[spec][0][0], ', {  // ', spec, sep='')
        for isot in out[spec]:
            isoname = f"{spec}-{isot[2]}"
            if isoname in _HITRAN_TO_ARTS_NAMES:
                isoname = _HITRAN_TO_ARTS_NAMES[isoname]
            print ('{', pos2char(isot[1]), ', "', isoname, '"},', sep='')
        print ('}},')
    print ('};')
