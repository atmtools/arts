from inspect import getmembers, ismodule
from pathlib import Path

import pyarts
from pyarts import arts


class TestStubFiles:
    def walk_modules(self, m):
        submodules = getmembers(m, ismodule)
        for submodule in submodules:
            if submodule[1].__name__.startswith("pyarts.arts."):
                yield submodule[1].__name__
                yield from self.walk_modules(submodule[1])

    def test_pyi(self):
        missing_stubs = []
        for module in self.walk_modules(arts):
            if (
                not Path(pyarts.__path__[0])
                .joinpath(module.replace("pyarts.", "").replace(".", "/") + ".pyi")
                .is_file()
            ):
                missing_stubs.append(module.replace("pyarts.arts.", ""))

        assert len(missing_stubs) == 0, (
            f"Missing stub file(s) for arts submodule(s): \n"
            f"  {"\n  ".join(missing_stubs)}\n"
            f"Add them to the `arts_add_cpp_stubs` call in "
            f"src/python_interface/CMakeLists.txt"
        )


if __name__ == "__main__":
    a = TestStubFiles()
    a.test_pyi()
