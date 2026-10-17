"""Generate (without replacing) the ARTS2 MC oblate ice particle.

Run with a T-matrix-enabled pyarts3, e.g.:
  PYTHONPATH=build/python/src python tools/compare_mc_tmatrix.py --reference /path/to/original.xml

The output is a comparison artifact in tmp/, not a new test fixture. No optical
values are read from the reference until after generation. The small REFICE
material table in the test is evaluated independently of the scattering reference.
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pyarts3 as pa

A = pa.arts
ROOT = Path(__file__).resolve().parents[1]
# Load definitions without executing the MC acceptance cases.
import runpy

generate = runpy.run_path(str(ROOT / "tests/core/scat/mc_general_arts2.tmatrix.py"))["generate"]

def compare(ssd, reference):
    report = {}
    for name in ['pha_mat_data', 'ext_mat_data', 'abs_vec_data']:
        actual, expected = np.asarray(getattr(ssd, name)), np.asarray(getattr(reference, name))
        assert actual.shape == expected.shape
        delta = actual - expected
        scale = np.max(np.abs(expected))
        report[name] = dict(
            max_absolute=float(np.max(np.abs(delta))),
            max_error_over_peak=float(np.max(np.abs(delta)) / scale),
            relative_l2=float(np.linalg.norm(delta.ravel()) / np.linalg.norm(expected.ravel())),
            worst_index=list(map(int, np.unravel_index(np.argmax(np.abs(delta)), delta.shape))),
        )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True, help='Original ARTS2 particle XML')
    parser.add_argument('--output', type=Path, default=ROOT / 'tmp/mc-tmatrix-comparison/azi-random.arts3.xml')
    parser.add_argument('--legacy-numerics', action='store_true')
    parser.add_argument('--quadrature', type=int, default=10)
    args = parser.parse_args()
    if args.quadrature < 2:
        parser.error('--quadrature must be at least 2')
    if args.output.resolve() == args.reference.resolve():
        parser.error('the original MC input must not be overwritten')
    start = perf_counter()
    ssd, orders = generate(args.legacy_numerics, args.quadrature)
    seconds = perf_counter() - start
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pa.xml.save(ssd, str(args.output), precision='.17g')
    # Check serialization independently of comparisons to old floating-point data.
    saved = pa.xml.load(str(args.output))
    for name in ['pha_mat_data', 'ext_mat_data', 'abs_vec_data']:
        np.testing.assert_array_equal(np.asarray(getattr(ssd, name)), np.asarray(getattr(saved, name)))
    report = dict(seconds=seconds, orders=orders,
                  legacy_numerics=args.legacy_numerics, quadrature=args.quadrature,
                  comparison=compare(saved, pa.xml.load(str(args.reference))))
    args.output.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
