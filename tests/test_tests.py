import pyarts
import os

def cmake_tests(cmake_list):
    tests = []
    
    for x in cmake_list:
        if x.startswith('pyfile') or x.startswith('ctlfile'):
            x = x.removeprefix('pyfile').removeprefix('ctlfile').strip()
            x = x.removeprefix('(').strip().split()
            x[1] = x[1].split(')')[0]
            assert x[0] == "fast", f"We only accept fast tests, not {x[0]}"
            tests.append(x[1])
            
    return tests

def arts_tests(dir):
    tests = []
    if os.path.isdir(dir):
        for d in os.listdir(dir):
            if d.endswith(".arts") or d.endswith(".py"):
                tests.append(os.path.join(dir, d))
            else:
                for x in arts_tests(os.path.join(dir, d)):
                    tests.append(x)
    elif dir.endswith(".arts") or dir.endswith(".py"):
        tests.append(dir)
    return tests

def curdir_has_all_tests():
    # Get all the tests run by cmake
    cmake = "CMakeLists.txt"
    assert os.path.isfile(cmake)
    cmake_list = open(cmake, 'r').read().split("arts_test_run_")
    run_tests = cmake_tests(cmake_list)
    assert len(run_tests) > 0, "Error in test, cannot find any cmake tests"
    
    # Find all the python and arts file in this folder and those below
    all_tests = []
    for dir in os.listdir():
        for x in arts_tests(dir):
            all_tests.append(x)
    assert len(all_tests) > 0, "Error in test, cannot find any test files"
    
    # All tests must be run
    for test in all_tests:
        assert test in run_tests, f"Cannot find {test} in cmake list of tests"
    
    # As a catch all, these lists should contain the same items but we run it
    # last for a better error message on failure
    all_tests.sort()
    run_tests.sort()
    assert all_tests == run_tests, f"{all_tests} != {run_tests}"


curdir = os.getcwd()
main_dirs = pyarts.environment.environ.get("ARTS_INCLUDE_PATH", '.').split(':')

test_dir = '.'
for dir in main_dirs:
    content = os.listdir(dir)
    if "this_is_the_main_test_directory" in content:
        test_dir = dir
        break

content = os.listdir(test_dir)
assert "this_is_the_main_test_directory" in content, \
    f"Cannot find a valid scenario in test directory: {os.path.abspath(test_dir)}"

# Running the tests will have you change directories, but we want to be back
# here at the end regardless
try:
    os.chdir(test_dir)
    curdir_has_all_tests()
    os.chdir("../examples")
    curdir_has_all_tests()
finally:
    os.chdir(curdir)
