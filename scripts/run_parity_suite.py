#!/usr/bin/env python3
"""Run immutable upstream R source against Python on nine shared-input cases."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from datetime import datetime, timezone

from compare_original_r_python import compare, TOLERANCES

ROOT=Path(__file__).resolve().parents[1]
UPSTREAM_SHA='1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66'
SOURCES=('SpatialPCA.R','SpatialPCA_EstimateLoading.R','SpatialPCA_SpatialPCs.R')


def run(command, log):
    with log.open('w') as stream:
        subprocess.run([str(x) for x in command],cwd=ROOT,stdout=stream,
                       stderr=subprocess.STDOUT,check=True)


def suite(args):
    upstream=args.upstream.resolve()
    revision=subprocess.check_output(['git','-C',str(upstream),'rev-parse','HEAD'],text=True).strip()
    if revision!=UPSTREAM_SHA:
        raise ValueError(f'Expected upstream {UPSTREAM_SHA}, got {revision}')
    if subprocess.check_output(['git','-C',str(upstream),'status','--porcelain','--','R'],text=True).strip():
        raise ValueError('Upstream R source was modified')
    hashes={name:hashlib.sha256((upstream/'R'/name).read_bytes()).hexdigest() for name in SOURCES}
    out=args.output.resolve()
    # Avoid accidentally comparing against outputs from a previous invocation.
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Output must be new or empty: {out}')
    out.mkdir(parents=True,exist_ok=True)
    manifest={'started_utc':datetime.now(timezone.utc).isoformat(),'upstream_revision':revision,
              'source_sha256':hashes,'python':sys.version,'platform':platform.platform(),
              'github_actions':os.environ.get('GITHUB_ACTIONS')=='true',
              'github_run_id':os.environ.get('GITHUB_RUN_ID'),
              'github_repository':os.environ.get('GITHUB_REPOSITORY'),
              'github_sha':os.environ.get('GITHUB_SHA'),
              'tolerances':TOLERANCES,'scope':'synthetic full-rank, intercept-only model core'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    run([args.rscript,'-e',"parse(file='scripts/run_author_original_r.R'); sessionInfo()"],
        out/'r_preflight.log')
    run([sys.executable,'-m','pip','freeze'],out/'python_dependencies.txt')
    rows=[]
    failures=[]
    for seed in args.seeds:
        for scale in args.scales:
            case=out/f'seed_{seed}_scale_{scale:g}'
            case.mkdir()
            print(f'Running seed={seed}, scale={scale:g}',flush=True)
            try:
                run([sys.executable,'scripts/prepare_original_r_input.py','--seed',seed,
                     '--scale',scale,'--output',case],case/'python.log')
                run([args.rscript,'scripts/run_author_original_r.R',upstream,case,2],case/'original_r.log')
                metadata=json.loads((case/'input_metadata.json').read_text())
                for name,digest in metadata['sha256'].items():
                    if hashlib.sha256((case/f'{name}.csv').read_bytes()).hexdigest()!=digest:
                        raise ValueError('Shared inputs changed during execution')
                report=compare(case)
                rows.extend(dict(seed=seed,scale=scale,**row) for row in report['rows'])
                if not report['profile_parity_passed']:
                    failures.append({'seed':seed,'scale':scale,'error':'profile parity tolerance exceeded'})
            except Exception as exc:
                failures.append({'seed':seed,'scale':scale,'error':str(exc)})
    manifest['finished_utc']=datetime.now(timezone.utc).isoformat()
    manifest['upstream_source_unchanged']=all(
        hashlib.sha256((upstream/'R'/name).read_bytes()).hexdigest()==digest
        for name,digest in hashes.items())
    if not manifest['upstream_source_unchanged']:
        failures.append({'error':'Upstream source changed'})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    (out/'suite_results.json').write_text(json.dumps({'manifest':manifest,'rows':rows,
        'failures':failures,'passed':not failures},indent=2,allow_nan=False))
    if rows:
        with (out/'suite_results.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    lines=['# Original R and Python comparison','',
           f"Execution: {'GitHub Actions' if manifest['github_actions'] else 'local runtime'}",
           f'Upstream commit: {revision}','','Scope: synthetic full-rank model core with identical Y and K; no Seurat/SPARK preprocessing.','',
           'Correction: the official upstream objective includes log(S). The earlier missing-log claim is withdrawn. no_log_ablation is deliberately altered Python.','',
           '| Seed | Scale | R tau | Python tau | Relative tau gap | W angle (deg) | Relative WZ gap | Profile parity |',
           '|---:|---:|---:|---:|---:|---:|---:|:---|']
    for r in rows:
        if r['python_mode']=='profile':
            lines.append(f"| {r['seed']} | {r['scale']:g} | {r['R_tau']:.9g} | {r['Python_tau']:.9g} | {r['tau_relative_gap']:.3g} | {r['max_subspace_angle_deg']:.3g} | {r['reconstruction_relative_gap']:.3g} | {'PASS' if r['within_profile_tolerances'] else 'FAIL'} |")
    lines.extend(['','## Scale invariance',''])
    for seed in args.seeds:
        selected=[r for r in rows if r['seed']==seed and r['python_mode']=='profile']
        if len(selected)==len(args.scales):
            for label,key in [('R','R_tau'),('Python','Python_tau')]:
                vals=[r[key] for r in selected]
                lines.append(f'- Seed {seed}, {label}: max/min tau = {max(vals)/min(vals):.9g}.')
    lines.extend(['','## Failures','',json.dumps(failures,indent=2) if failures else 'None.',
                  '','Boundary solutions are explicitly flagged in the CSV/JSON and do not establish an unconstrained optimum.',
                  'No claim is made about real biological datasets or the full published pipeline.'])
    (out/'experiment_results.md').write_text('\n'.join(lines)+'\n')
    print(f'Completed {len(rows)//2}/{len(args.seeds)*len(args.scales)} cases; failures={len(failures)}',flush=True)
    if failures:raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--upstream',type=Path,default=ROOT/'upstream/SpatialPCA')
    p.add_argument('--output',type=Path,default=ROOT/'results/original_r_suite')
    p.add_argument('--rscript',default='Rscript')
    p.add_argument('--seeds',type=int,nargs='+',default=[23,41,2026])
    p.add_argument('--scales',type=float,nargs='+',default=[.5,1.,2.])
    suite(p.parse_args())
