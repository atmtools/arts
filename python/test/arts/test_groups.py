from pyarts.arts import *  # noqa
from copy import copy, deepcopy


class TestGroups:

    def test_xml(self):
        ignore_groups = ["CallbackFunction", "NumericUnaryOperator",
                         "SpectralRadianceProfileOperator", "SingleScatteringData"]

        groups = list(globals.workspace_groups().keys())
        groups.sort()

        for group in ignore_groups:
           if group not in groups:
               raise Exception(f"Ignored group {group} is not in workspace")

        fail = []
        for group in groups:
            if group in ignore_groups:
                print(f"Skipping group {group}")
                continue

            try:
              t = eval(group)
              print(f"Testing group {group}")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, does not exist:\n{e}\n\n")
              continue

            try:
              x = t()
              print(f"  Can init as {group}()")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot init:\n{e}\n\n")
              continue

            try:
              x.savexml("test.xml")
              print(f"  Can save as xml ({group}.savexml())")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot save:\n{e}\n\n")
              continue

            try:
              x.readxml("test.xml")
              print(f"  Can read as xml ({group}.readxml())")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot read:\n{e}\n\n")
              continue

            try:
              t.fromxml("test.xml")
              print(f"  Can init ax xml ({group}.fromxml())")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot init from file:\n{e}\n\n")
              continue

        if len(fail) > 0:
          for g in fail:
            print(g)
          raise Exception(f"There are {len(fail)} tests failing")

    def test_print(self):
        ignore_groups = []

        groups = list(globals.workspace_groups().keys())
        groups.sort()

        for group in ignore_groups:
           if group not in groups:
               raise Exception(f"Ignored group {group} is not in workspace")

        fail = []
        for group in groups:
            if group in ignore_groups:
                print(f"Skipping group {group}")
                continue

            try:
              t = eval(group)
              print(f"Testing group {group}")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, does not exist:\n{e}\n\n")
              continue

            try:
              x = t()
              print(f"  Can init as {group}()")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot init:\n{e}\n\n")
              continue

            try:
              print(x)
              print(f"  Can print (print({group}()))")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot print:\n{e}\n\n")
              continue

            try:
              assert isinstance(repr(x), str), f"repr not a str but {type(repr(x))}"
              print(f"  Can repr (isinstance(repr({group}()), str))")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot repr:\n{e}\n\n")
              continue

        if len(fail) > 0:
          for g in fail:
            print(g)
          raise Exception(f"There are {len(fail)} tests failing")

    def test_copy(self):
        ignore_groups = []

        groups = list(globals.workspace_groups().keys())
        groups.sort()

        for group in ignore_groups:
           if group not in groups:
               raise Exception(f"Ignored group {group} is not in workspace")

        fail = []
        for group in groups:
            if group in ignore_groups:
                print(f"Skipping group {group}")
                continue

            try:
              t = eval(group)
              print(f"Testing group {group}")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, does not exist:\n{e}\n\n")
              continue

            try:
              x = t()
              print(f"  Can init as {group}()")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot init:\n{e}\n\n")
              continue

            try:
              y = t(x)
              print(f"  Can copy init (x = {group}({group}()))")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot copy init:\n{e}\n\n")
              continue

            try:
              y = copy(x)
              print(f"  Can copy init (x = copy({group}()))")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot copy:\n{e}\n\n")
              continue

            try:
              y = deepcopy(x)
              print(f"  Can copy init (x = deepcopy({group}()))")
            except Exception as e:
              fail.append(f"\n\nFailed {group}, cannot deepcopy:\n{e}\n\n")
              continue

        if len(fail) > 0:
          for g in fail:
            print(g)
          raise Exception(f"There are {len(fail)} tests failing")


if __name__ == "__main__":
    x = TestGroups()
    x.test_xml()
    x.test_copy()
