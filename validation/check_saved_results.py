#!/usr/bin/env python3
"""Read-only checks of saved CSV outputs. No simulation, no dependencies."""
import csv, math, pathlib, statistics, collections, sys
ROOT=pathlib.Path(__file__).resolve().parents[1]
TABLE=ROOT/'results'/'tables'

def rows(name):
    with (TABLE/name).open(newline='') as f:return list(csv.DictReader(f))

def near(a,b,tol=1e-9):return math.isclose(float(a),float(b),rel_tol=tol,abs_tol=tol)

def check_e6(w):
    raw=rows(f'E6_replications_analytic_W{w}.csv')
    summary=rows(f'E6_summary_analytic_W{w}.csv')
    paired=rows(f'E6_paired_differences_analytic_W{w}.csv')
    assert len(raw)==13500,(w,len(raw))
    groups=collections.defaultdict(list)
    for r in raw:groups[(r['scenario'],r['depth'],r['method'])].append(r)
    assert len(groups)==27
    for g,arr in groups.items():
        assert len(arr)==500,(w,g,len(arr))
    for r in summary:
        g=groups[(r['scenario'],r['depth'],r['method'])]
        x=[float(z[r['metric']]) for z in g if z[r['metric']] not in ('','NA','NaN')]
        assert len(x)==int(r['n_valid'])
        assert near(statistics.mean(x),r['mean'])
        assert near(statistics.stdev(x),r['sd'])
        assert near(statistics.stdev(x)/math.sqrt(len(x)),r['mc_se'])
    assert len(paired)==54,(w,len(paired))
    for r in paired:
        g=(r['scenario'],r['depth']);metric=r['metric']
        aa=sorted(groups[(*g,'SMPM')],key=lambda z:int(z['replication']))
        bb=sorted(groups[(*g,r['competitor'])],key=lambda z:int(z['replication']))
        assert [z['replication'] for z in aa]==[z['replication'] for z in bb]
        d=[float(a[metric])-float(b[metric]) for a,b in zip(aa,bb)]
        assert near(statistics.mean(d),r['mean_difference'])
        assert near(statistics.stdev(d)/math.sqrt(len(d)),r['mc_se_difference'])
    print(f'E6 W={w}: {len(raw)} replicate-method rows; {len(summary)} summary rows; {len(paired)} paired rows — PASS')

def quantile_type7(values,p):
    x=sorted(values);h=(len(x)-1)*p;lo=math.floor(h);hi=math.ceil(h)
    return x[lo]+(h-lo)*(x[hi]-x[lo])

def check_e4():
    raw=rows('E4_replications.csv');summary=rows('E4_summary.csv')
    assert len(raw)==2000
    groups=collections.defaultdict(list)
    for r in raw:groups[r['rho']].append(r)
    assert len(groups)==4 and all(len(v)==500 for v in groups.values())
    for r in summary:
        x=[float(z[r['metric']]) for z in groups[r['rho']]]
        assert near(statistics.mean(x),r['mean'])
        assert near(statistics.stdev(x),r['sd'])
        assert near(quantile_type7(x,.025),r['q025'])
        assert near(quantile_type7(x,.975),r['q975'])
    print(f'E4: {len(raw)} replicate rows; {len(summary)} summary rows recalculated — PASS')

def check_e5():
    raw=rows('E5_replications.csv');summary=rows('E5_summary.csv')
    assert len(raw)==4500
    groups=collections.defaultdict(list)
    for r in raw:groups[(r['method'],r['depth'])].append(r)
    assert len(groups)==9 and all(len(v)==500 for v in groups.values())
    for r in summary:
        x=[float(z[r['metric']]) for z in groups[(r['method'],r['depth'])]]
        assert near(statistics.mean(x),r['mean'])
        assert near(statistics.stdev(x),r['sd'])
        assert near(statistics.stdev(x)/math.sqrt(len(x)),r['mc_se'])
    print(f'E5: {len(raw)} replicate rows; {len(summary)} summary rows recalculated — PASS')

if __name__=='__main__':
    check_e4();check_e5();check_e6(20);check_e6(40)
    print('PASS: archived CSV consistency checks completed (not an RDS or manuscript audit).')
