"""Compare complete OEM timings; emit all rows, including regressions."""
import csv
import sys
from pathlib import Path
import statistics
root=Path('doc/arts/benchmarks/oem_48')
keys=['scenario','measurements','states','cache_case','mixed','method','gain','phase']
def load(name):
    rows=list(csv.DictReader((root/name).open()))
    assert len(rows)==1920,(name,len(rows))
    return {tuple(r[k] for k in keys):r for r in rows}
before=load(sys.argv[3] if len(sys.argv)>3 else 'head.csv');after=load(sys.argv[1] if len(sys.argv)>1 else 'working.csv')
prefix=sys.argv[2] if len(sys.argv)>2 else ''
assert before.keys()==after.keys()
output=[]
for k,a in before.items():
    b=after[k]
    xa=[float(x) for x in a['state'].split(';') if x]
    xb=[float(x) for x in b['state'].split(';') if x]
    assert len(xa)==len(xb)
    error=max(abs(x-y) for x,y in zip(xa,xb))
    ratio=float(b['median_ms'])/float(a['median_ms'])
    numerical=(error<=1e-6 and abs(float(a['cost'])-float(b['cost']))<=1e-7*(1+abs(float(a['cost']))) and a['status']==b['status'])
    flag='numerical_change' if not numerical else ('regression' if ratio>1.10 else 'improvement' if ratio<0.90 else 'within_10_percent')
    separated=float(b['min_ms'])>float(a['max_ms'])
    output.append(dict(zip(keys,k),before_ms=a['median_ms'],after_ms=b['median_ms'],after_over_before=ratio,max_state_difference=error,status_match=a['status']==b['status'],before_iterations=a['iterations'],after_iterations=b['iterations'],se_inverse_before=a['se_inverse_blocks'],se_inverse_after=b['se_inverse_blocks'],flag=flag,slower_ranges_disjoint=separated))
with (root/(prefix+'comparison.csv')).open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
lines=['# Complete OEM comparison','', '48 scenarios; 10 methods × gain on/off × cold/warm = 1,920 matched rows.', '', 'Ratio is current / HEAD: greater than one is slower. Flags use a 10% descriptive threshold, not a statistical significance test. Disjoint sample ranges are recorded separately.','']
for label in ['regression','improvement','within_10_percent','numerical_change']:
    lines.append(f'- {label}: {sum(r["flag"]==label for r in output)}')
lines+=['','| Cache case | Phase | Median ratio across rows | >10% regressions / rows |','|---|---|---:|---:|']
for cache in ['dense_inverse','sparse_inverse','diagonal','cholesky']:
    for phase in ['cold','warm']:
        group=[r for r in output if r['cache_case']==cache and r['phase']==phase]
        lines.append(f'| {cache} | {phase} | {statistics.median(r["after_over_before"] for r in group):.3f} | {sum(r["flag"]=="regression" for r in group)}/{len(group)} |')
lines+=['','## Largest absolute slowdowns','', '| N | M | Cache | Mixed | Method | Gain | Phase | HEAD ms | Current ms | Ratio |', '|---:|---:|---|---|---|---|---|---:|---:|---:|']
for r in sorted(output,key=lambda r:float(r['after_ms'])-float(r['before_ms']),reverse=True)[:20]:
    lines.append('| '+' | '.join(str(r[k]) for k in ['measurements','states','cache_case','mixed','method','gain','phase','before_ms','after_ms'])+f' | {r["after_over_before"]:.2f} |')
(root/(prefix+'summary.md')).write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[:20]))
with (root/(prefix+'scenario_summary.csv')).open('w') as f:
    w=csv.writer(f)
    w.writerow(['scenario','measurements','states','cache_case','mixed','rows','regressions','improvements','median_current_over_head','worst_ratio','worst_method','worst_gain','worst_phase'])
    for scenario in range(1,49):
        g=[r for r in output if r['scenario']==str(scenario)]
        a=g[0]; worst=max(g,key=lambda r:r['after_over_before'])
        w.writerow([scenario,a['measurements'],a['states'],a['cache_case'],a['mixed'],len(g),sum(r['flag']=='regression' for r in g),sum(r['flag']=='improvement' for r in g),statistics.median(r['after_over_before'] for r in g),worst['after_over_before'],worst['method'],worst['gain'],worst['phase']])
