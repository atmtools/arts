"""The licences of the bundled code are compiled into ARTS.

Every component registered with arts_add_bundled_license
(cmake/modules/ArtsBundledLicenses.cmake) must be reported by
pyarts3.arts.globals.bundled_components() with its complete licence files,
exactly when this build compiles it in, and the licence expression must
cover all of them.  In a build tree, the licence files of the Python package
(build/python/licenses) must match the compiled-in ones.
"""

from pathlib import Path

import pyarts3 as pyarts

g = pyarts.arts.globals
components = {c.name: c for c in g.bundled_components()}
expression = g.license_expression()
print(expression)
for c in components.values():
    print(f"  {c!r}: {sorted(c.files)}")

assert len(components) == len(g.bundled_components()), "component names must be unique"
assert expression.startswith("(LGPL-3.0-or-later OR GPL-3.0-or-later)"), expression

always = {"invlib", "mdspan", "Faddeeva", "wigxjpf"}
optional = {
    "cdisort": g.data.has_cdisort,
    "tmatrix": pyarts.arts.tmatrix.available(),
    "polradtran": pyarts.arts.rt3.available() or pyarts.arts.rt4.available(),
}
assert always <= set(components), f"missing always-bundled components: {always - set(components)}"
for name, built in optional.items():
    assert (
        name in components
    ) == built, f"{name} is compiled in: {built}, but registered: {name in components}"

for c in components.values():
    assert c.files, f"{c.name} has no licence files"
    for file, text in c.files.items():
        assert text.strip(), f"{c.name}/{file} is empty"
    assert c.spdx in expression, f"the licence expression must include {c.name}'s {c.spdx}"

if "polradtran" in components:
    assert "K. Franklin Evans" in components["polradtran"].files["LICENSE"]

# In a build tree, the Python package's licence files are the compiled-in ones, byte for byte
# (read_text would decode with the locale's encoding and translate CRLF line endings, as a
# Windows checkout has them, while the compiled-in text keeps the file's bytes)
packaged = Path(pyarts.__file__).resolve().parents[2] / "licenses"
if packaged.is_dir():
    assert {p.name for p in packaged.iterdir() if p.is_dir()} == set(components)
    for c in components.values():
        for file, text in c.files.items():
            assert (packaged / c.name / file).read_bytes() == text.encode("utf-8"), f"{c.name}/{file} differs"
