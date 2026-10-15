"""
Two workspaces must not share workspace-variable payloads.

Wsv holds its value through a std::shared_ptr, so a plain copy of the global
defaults map shares the underlying value objects between workspaces.  A
workspace method writing its `out` variable in place then silently changes
that variable in every other workspace as well.  The Workspace constructor
must deep-copy the defaults (Wsv::copied) instead.
"""

import pyarts3 as pyarts

a = pyarts.Workspace()
b = pyarts.Workspace()

assert b.jac_targets.target_count() == 0, "the default jac_targets is not empty"

a.jac_targetsAddTemperature()  # writes the `out` variable `jac_targets` of `a`

assert a.jac_targets.target_count() == 1, "the method did not write its out variable"
assert b.jac_targets.target_count() == 0, "workspaces share their `jac_targets`"

# The global defaults must stay untouched too, so a workspace created after
# the method call is still pristine.
c = pyarts.Workspace()
assert c.jac_targets.target_count() == 0, "the global defaults were modified"

print("Workspace defaults are deep-copied; workspaces are independent")
