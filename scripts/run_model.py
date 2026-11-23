#!/usr/bin/env python3
"""Fit SpatialPCA Python model-core to NPZ data (Y, coords, optional covariates)."""
import argparse
import pathlib
import sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from spatialpca_verify.model import fit_spatialpca, standardize_genes, log_normalize_counts


def main():
    p = argparse.ArgumentParser(description='Independent SpatialPCA model-core reproduction')
    p.add_argument('--input', required=True, type=pathlib.Path,
                   help='NPZ with Y (genes x spots) and coords (spots x spatial dimensions)')
    p.add_argument('--output', required=True, type=pathlib.Path)
    p.add_argument('--components', type=int, default=2)
    p.add_argument('--bandwidth', type=float, default=.3)
    p.add_argument('--kernel-rank', type=int, default=None)
    p.add_argument('--objective', choices=['profile','no_log_ablation'], default='profile')
    p.add_argument('--counts-log-normalize', action='store_true',
                   help='Demonstration only; NOT SCTransform')
    p.add_argument('--standardize', action='store_true', help='Per-gene z-score')
    a = p.parse_args()
    data = np.load(a.input)
    Y, coords = data['Y'], data['coords']
    if a.counts_log_normalize:
        Y = log_normalize_counts(Y)
    cov = data['covariates'] if 'covariates' in data else None
    result = fit_spatialpca(Y, coords, n_components=a.components,
                            bandwidth=a.bandwidth, covariates=cov,
                            kernel_rank=a.kernel_rank, objective=a.objective,
                            standardize=a.standardize)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(a.output, W=result.W, Z=result.Z, beta=result.beta,
                        tau=result.tau, sigma2=result.sigma2,
                        Y_fitted=result.fitted_expression, objective=result.objective,
                        residual_sum_S=result.residual_sum_S,
                        kernel_rank=result.rank)
    print(f'W {result.W.shape}, Z {result.Z.shape}; tau {result.tau:.6g}; '
          f'sigma2 {result.sigma2:.6g}; saved {a.output}')

if __name__ == '__main__':
    main()
