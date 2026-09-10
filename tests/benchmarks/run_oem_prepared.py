"""Benchmark the current prepared implementation against the preserved results; restore the temporary test harness.
Run from repository root. Do not edit repository files while this runs.
Untracked files, the index, and branches are never changed. Backups persist in /tmp.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

root=Path.cwd()
backup=Path(tempfile.mkdtemp(prefix='arts-oem-comparison-'))
print(f'BACKUP={backup}',flush=True)
out=root/'doc/arts/benchmarks/oem_48'
out.mkdir(parents=True,exist_ok=True)
paths=subprocess.check_output(['git','diff','HEAD','--name-only','-z']).decode().split('\0')
paths=[p for p in paths if p]
target='src/tests/test_oem_methods.cc'
if target not in paths: paths.append(target)
original={p:(root/p).read_bytes() if (root/p).exists() else None for p in paths}
for p,data in original.items():
 if data is not None:
  dest=backup/p;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
head=subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()
metadata={'head':head,'backup':str(backup),'sha256':{p:hashlib.sha256(d).hexdigest() for p,d in original.items() if d is not None}}
(out/'prepared_sources.json').write_text(json.dumps(metadata,indent=2)+'\n')
base_fixture=subprocess.check_output(['git','show',f'HEAD:{target}']).decode().split('void check_history(')[0]
harness='#include <chrono>\n#include <iomanip>\n#include <iostream>\n'+base_fixture+(root/'tests/benchmarks/oem_scenarios.inc').read_text()+'\n} // namespace\nint main() try { scenario_benchmark(); } catch(const std::exception& e) {std::cerr<<e.what()<<"\\n";return 1;}\n'
env=dict(os.environ,CCACHE_READONLY='1',CCACHE_TEMPDIR='/tmp/arts-oem-ccache',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
def restore():
 for p,data in [(target,original[target])]:
  dest=root/p
  if data is None: dest.unlink(missing_ok=True)
  elif not dest.exists() or dest.read_bytes()!=data: dest.write_bytes(data)
def run(cmd,log):
 with (out/log).open('w') as f: subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
try:
 (root/target).write_text(harness)
 print('BUILD prepared',flush=True)
 run(['ninja','-C','build','-j6','test_oem_methods'],'prepared_build.log')
 exe=backup/'oem_prepared';shutil.copy2(root/'build/src/tests/test_oem_methods',exe)
 print('RUN prepared',flush=True)
 run([str(exe)],'prepared.csv')
finally:
 restore()
 print('RESTORED regression test source',flush=True)
print('REBUILD regression executable',flush=True)
run(['ninja','-C','build','-j6','test_oem_methods'],'restored_build.log')
print('DONE',flush=True)
