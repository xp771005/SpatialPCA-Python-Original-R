#!/usr/bin/env python3
"""Run seeded full-rank comparison and save CSV, JSON, and a diagnostic plot."""
import argparse
import csv
import json
import pathlib
import sys

root = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / 'src'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from spatialpca_verify import make_synthetic, optimize_tau, compare_fits


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=pathlib.Path, default=root / 'results')
    parser.add_argument('--seeds', type=int, nargs='+', default=[23, 41, 2026])
    parser.add_argument('--scales', type=float, nargs='+', default=[0.5, 1.0, 2.0])
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for seed in args.seeds:
        for scale in args.scales:
            X,Y,K,_,_,_ = make_synthetic(seed=seed, scale=scale)
            a = optimize_tau(X,Y,K,2,objective='profile')
            b = optimize_tau(X,Y,K,2,objective='no_log_ablation')
            metrics = compare_fits(a,b)
            rows.append({
                'seed':seed, 'scale':scale, 'profile_tau':a.tau,
                'no_log_ablation_tau':b.tau,
                'profile_sigma2':a.sigma2,
                'no_log_ablation_sigma2':b.sigma2,
                'profile_S':a.S, 'no_log_ablation_S':b.S,
                'ablation_bound_hit':bool(abs(b.log_tau) >= 9.999),
                'profile_objective_at_profile_tau':a.objective_value,
                # True profiled objective at the alternative tau:
                'profile_objective_at_ablation_tau':
                    __import__('spatialpca_verify').fit_at_log_tau(
                        X,Y,K,2,b.log_tau,'profile').objective_value,
                **metrics})
            print(f'seed={seed:4d} scale={scale:3g} | '
                  f'profile tau={a.tau:10.3f} | Ablation tau={b.tau:10.3f} | '
                  f'W angle={metrics["largest_subspace_angle_degrees"]:7.2f} deg | '
                  f'reconstruction gap={metrics["relative_reconstruction_difference"]:7.3f}')

    with (args.output_dir/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    with (args.output_dir/'summary.json').open('w',encoding='utf-8') as f:
        json.dump({'experiment': 'hypothetical no-log ablation; not original R output',
                   'formula_source':'profile likelihood and deliberately altered no-log objective',
                   'seeds':args.seeds,'scales':args.scales,'rows':rows},f,indent=2)

    # One standalone plot, with all cases from one representative seed.
    seed=args.seeds[-1]
    subset=sorted((r for r in rows if r['seed']==seed),key=lambda x:x['scale'])
    x=[r['scale'] for r in subset]
    plt.figure(figsize=(7.2,4.4))
    plt.plot(x,[r['profile_tau'] for r in subset],marker='o',label='Profile likelihood')
    plt.plot(x,[r['no_log_ablation_tau'] for r in subset],marker='s',label='No-log ablation (not R)')
    plt.yscale('log')
    plt.xlabel('Multiplicative scaling of the same Y')
    plt.ylabel('Estimated tau (log scale)')
    plt.title(f'SpatialPCA objective comparison: fixed seed {seed}')
    plt.grid(True,which='both',alpha=0.22)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_dir/'scale_sensitivity.png',dpi=170)
    plt.close()
    print('Saved to',args.output_dir)

if __name__=='__main__':
    main()
