def gen_latest_molparam_map(molparam_txt_file):
    """ Generates a version of latest_molparam_map used in hitran_species.cc

    The variable is simply printed to stream with print() as the intent is
    to use this output in ARTS directly.  Note, to keep this simple, however,
    that since HITRAN does not know about ARTS tags, the species name is
    given as Species-AFGL.  ARTS does not use AFGL notation internally so
    some of these have to be manually changed.
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
    
    molparam_txt = open(molparam_txt_file, 'r').read().split('\n')
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
            print ('{', pos2char(isot[1]), ', "{}-{}"'.format(spec, isot[2]), '},', sep='')
        print ('}},')
    print ('};')
