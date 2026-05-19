"""This test has 2 parts:
1. Check that the naming of species data is consistent.
2. Check that the convention of the species data is correct.  The latter should only fail if you have an out-of-date version of the data.
"""

import os
import pyarts3 as pyarts


def get_arts_data_path():
    path = None
    test = False
    for p in pyarts.arts.globals.parameters.datapath:
        if "arts-cat-data" in p:
            path = p
            test = True
            break

    if not test:
        raise RuntimeError(
            "Data path does not include arts-cat-data. Please update your data path to include arts-cat-data."
        )

    return path


def test_iso(fn: str):
    data = pyarts.xml.load(fn)

    assert isinstance(data, pyarts.arts.SpeciesIsotopologueInfo), (
        f"Wrong type for {fn}: {type(data)}"
    )

    iso = os.path.basename(fn).split(".")[0]
    [spec, code] = iso.split("-")

    assert data.species == spec or data.code == code, (
        f"[Unmaintainable isotopologues data]: {fn} has species {data.species} and code {data.code}, but expected {spec} and {code} from filename"
    )

    assert data.mass > 0, (
        f"[Unmaintainable isotopologues data]: {fn} has bad mass {data.mass}.  Set to round atomic weights if you do not know."
    )

    assert data.default_ratio >= 0, (
        f"[Unmaintainable isotopologues data]: {fn} has bad default_ratio {data.default_ratio}.  Set to zero to ignore by default."
    )

    assert data.degeneracy > 0 or data.degeneracy == -1, (
        f"[Unmaintainable isotopologues data]: {fn} has invalid degeneracy {data.degeneracy}.  Set to -1 if unknown."
    )

    return iso, data


def test_hitran(fn: str):
    data = pyarts.xml.load(fn)

    assert isinstance(data, pyarts.arts.HitranSpeciesInfo), (
        f"Wrong type for {fn}: {type(data)}"
    )

    try:
        iso = os.path.basename(fn).split(".")[0]
        [spec, code] = iso.split("-")
    except (IndexError, ValueError):
        raise ValueError(
            f"Invalid filename format for {fn}.  Should be <spec>-<code>.xml, where spec is a string and code is the isotopologue code."
        )

    assert data.spec.spec == spec or data.spec.isotname == code, (
        f"[Unmaintainable hitran data]: {fn} has species {data.species} and code {data.code}, but expected {spec} and {code} from filename"
    )

    assert data.ratio >= 0, (
        f"[Unmaintainable hitran data]: {fn} has bad ratio {data.ratio}.  Must be positive."
    )

    return iso, data


def test_jpl(fn: str):
    data = pyarts.xml.load(fn)

    assert isinstance(data, pyarts.arts.JplSpeciesInfo), (
        f"Wrong type for {fn}: {type(data)}"
    )

    try:
        iso = os.path.basename(fn).split(".")
        code = int(iso[0])
    except (IndexError, ValueError):
        raise ValueError(
            f"Invalid filename format for {fn}.  Should be <code>.xml, where code is an integer from JPL."
        )

    assert data.id == code, (
        f"[Unmaintainable jpl data]: {fn} has id {data.id}, but expected {code} from filename"
    )

    assert data.T0 > 0, (
        f"[Unmaintainable jpl data]: {fn} has bad T0 {data.T0}.  Must be positive."
    )

    assert data.QT0 > 0, (
        f"[Unmaintainable jpl data]: {fn} has bad QT0 {data.QT0}.  Must be positive."
    )

    return str(data.spec), data


def test_species(fn: str):
    data = pyarts.xml.load(fn)

    assert isinstance(data, pyarts.arts.SpeciesEnumInfo), (
        f"Wrong type for {fn}: {type(data)}"
    )

    try:
        iso = os.path.basename(fn).split(".")
        code = int(iso[0])
        spec = iso[1]
    except (IndexError, ValueError):
        raise ValueError(
            f"Invalid filename format for {fn}.  Should be <code>.<spec>.xml, where code is an integer and spec is a string."
        )

    assert data.enum_value == code, (
        f"[Unmaintainable species data]: {fn} has enum_value {data.enum_value}, but expected {code} from filename"
    )

    assert data.shortname == spec, (
        f"[Unmaintainable species data]: {fn} has shortname {data.shortname}, but expected {spec} from filename"
    )

    return iso[1], data


def test_partfun(fn: str):
    data = pyarts.xml.load(fn)

    assert isinstance(data, pyarts.arts.PartitionFunctionsData), (
        f"Wrong type for {fn}: {type(data)}"
    )

    try:
        iso = os.path.basename(fn).split(".")
        spec = iso[0]
    except (IndexError, ValueError):
        raise ValueError(
            f"Invalid filename format for {fn}.  Should be <spec>.xml, where spec is a string."
        )

    return spec, data


def all_isotopologues():
    path = get_arts_data_path()

    out = {}
    isot_path = os.path.join(path, "isotopologues")
    for file in os.listdir(isot_path):
        if not file.endswith(".xml"):
            continue

        key, ext = test_iso(os.path.join(isot_path, file))
        out[key] = ext

    return out


def all_hitran():
    path = get_arts_data_path()

    out = {}
    hitran_path = os.path.join(path, "hitran")
    for file in os.listdir(hitran_path):
        if not file.endswith(".xml"):
            continue

        key, ext = test_hitran(os.path.join(hitran_path, file))
        out[key] = ext

    testspec = {}
    testind = {}
    testchar = {}

    for spec in out:
        data = out[spec]

        if data.spec.spec not in testspec:
            testspec[data.spec.spec] = []

        if data.hitind not in testind:
            testind[data.hitind] = []

        if data.hitchar not in testchar:
            testchar[data.hitchar] = []

        # unique characters for each species
        testspec[data.spec.spec].append(data.hitchar)
        testind[data.hitind].append(data.spec.spec)  # same species for same
        # different species for same character
        testchar[data.hitchar].append(data.hitind)

    for spec in testspec:
        n = len(testspec[spec])
        x = set(testspec[spec])
        assert len(x) == n, (
            f"[Unmaintainable isotopologues data]: Species {spec} has duplicate characters: {testspec[spec]}"
        )

    for ind in testind:
        n = len(testind[ind])
        x = set(testind[ind])
        assert len(x) == 1, (
            f"[Unmaintainable isotopologues data]: Hitran index {ind} corresponds to multiple species: {testind[ind]}"
        )

    for char in testchar:
        n = len(testchar[char])
        x = set(testchar[char])
        assert len(x) == n, (
            f"[Unmaintainable isotopologues data]: Character {char} corresponds to multiple hitran indices: {testchar[char]}"
        )

    return out


def all_jpl():
    path = get_arts_data_path()

    out = {}
    jpl_path = os.path.join(path, "jpl")
    for file in os.listdir(jpl_path):
        if not file.endswith(".xml"):
            continue

        key, ext = test_jpl(os.path.join(jpl_path, file))
        out[key] = ext

    testcode = [spec for spec in out]
    assert len(testcode) == len(set(testcode)), (
        f"[Unmaintainable jpl data]: Duplicate codes found: {testcode}"
    )

    return out


def all_species():
    path = get_arts_data_path()

    out = {}
    species_path = os.path.join(path, "species")
    for file in os.listdir(species_path):
        if not file.endswith(".xml"):
            continue

        key, ext = test_species(os.path.join(species_path, file))
        out[key] = ext

    assert len(out) == len(set(out.keys())), (
        f"[Unmaintainable species data]: Duplicate species found: {list(out.keys())}"
    )

    return out


def all_partfun():
    path = get_arts_data_path()

    out = {}
    partfun_path = os.path.join(path, "partition-functions")
    for file in os.listdir(partfun_path):
        if not file.endswith(".xml"):
            continue

        key, ext = test_partfun(os.path.join(partfun_path, file))
        out[key] = ext

    return out


iso = all_isotopologues()
hit = all_hitran()
jpl = all_jpl()
specs = all_species()
pfun = all_partfun()

iso_keys = set(iso.keys())
hit_keys = set(hit.keys())
jpl_keys = set(jpl.keys())
specs_keys = set(specs.keys())
pfun_keys = set(pfun.keys())

specs_from_iso = set([s.split("-")[0] for s in iso_keys])
specs_from_hit = set([s.split("-")[0] for s in hit_keys])
specs_from_jpl = set([s.split("-")[0] for s in jpl_keys])

assert iso_keys.issubset(pfun_keys), (
    f"[Unmaintainable data]: Partition function isotopologues do not match isotopologues.  Missing: {iso_keys - pfun_keys}"
)

assert hit_keys.issubset(iso_keys), (
    f"[Unmaintainable data]: HITRAN isotopologues do not match isotopologues.  Missing: {hit_keys - iso_keys}"
)

assert jpl_keys.issubset(iso_keys), (
    f"[Unmaintainable data]: JPL isotopologues do not match isotopologues.  Missing: {jpl_keys - iso_keys}"
)

assert specs_from_iso.issubset(specs_keys), (
    f"[Unmaintainable data]: Isotopologue species do not match species.  Missing: {specs_from_iso - specs_keys}"
)

assert specs_from_hit.issubset(specs_keys), (
    f"[Unmaintainable data]: HITRAN species do not match species.  Missing: {specs_from_hit - specs_keys}"
)

assert specs_from_iso.issubset(specs_keys), (
    f"[Unmaintainable data]: JPL species do not match species.  Missing: {specs_from_iso - specs_keys}"
)
