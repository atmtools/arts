Zeeman radiative-transfer examples
==================================

These examples exercise polarized clear-sky radiative transfer around the
118.75 GHz oxygen line.  In addition to geometry, sensors, refraction, solar
backgrounds, and adaptive paths, the directory compares the available
within-layer propagation and source approximations selected by
``ws.rte_option``.

The relevant examples are:

* ``2-zeeman.py`` uses ``constant``;
* ``9-zeeman-linear-src.py`` uses ``lintau``;
* ``10-zeeman-linear-prop.py`` uses ``linprop``;
* ``11-zeeman-magnus.py`` uses ``magop``; and
* ``12-zeeman-magnus-linear-src.py`` uses ``magop_linsrc``.

``magop`` represents the propagation matrix as linear across each layer and
includes the first commutator correction of the Magnus expansion, while using
the endpoint-average source.  ``magop_linsrc`` additionally represents the
source as linear across the layer using the augmented Magnus source operator.
Both options support polarized analytical Jacobians.

The regression arrays sampled in these examples are intentional: the Magnus
results need not equal ``constant``, ``lintau``, or ``linprop`` when adjacent
polarized propagation matrices do not commute.  For a constant or commuting
propagation matrix, ``magop`` reduces to the ordinary endpoint-average matrix
exponential.
