#!/usr/bin/env python3
"""Export precisely the same normalized expression and Gaussian kernel for upstream R and Python.

This bypasses Seurat SCTransform/SPARK deliberately: compare the original model-core R
implementation (sourced from the author's repository) to our Python model-core first.
"""
from pathlib import Path
import argparse
import json
import sys
import hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from spatialpca_verify.model import gaussian_kernel, fit_spatialpca
from spatialpca_verify.core import make_synthetic


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,default=ROOT/'data'/'demo_synthetic.npz')
    parser.add_argument('--output',type=Path,default=ROOT/'results'/'original_r_parity')
    parser.add_argument('--bandwidth',type=float,default=0.19)
    parser.add_argument('--components',type=int,default=2)
    parser.add_argument('--seed',type=int,default=None,
                        help='Generate 16 genes x 25 spots instead of reading --input')
    parser.add_argument('--scale',type=float,default=1.)
    args=parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if not np.isfinite(args.scale) or args.scale <= 0:
        raise ValueError('Scale must be positive and finite')
    if args.seed is None:
        raw=np.load(args.input)
    else:
        _, y, _, _, _, coords = make_synthetic(seed=args.seed)
        raw={'Y':y,'coords':coords}
    Y=np.asarray(raw['Y'],dtype=float)
    Y=Y*args.scale
    coords=np.asarray(raw['coords'],dtype=float)
    if 'covariates' in raw:
        raise ValueError('This initial upstream-R parity harness uses intercept-only data; extra covariates need explicit parity handling.')
    K=gaussian_kernel(coords,bandwidth=args.bandwidth)+1e-8*np.eye(coords.shape[0])
    for name,val in [('Y',Y),('coords',coords),('K',K)]:
        np.savetxt(args.output/f'{name}.csv',val,delimiter=',',fmt='%.17g')
    for objective in ('no_log_ablation','profile'):
        result=fit_spatialpca(Y,coords,n_components=args.components,bandwidth=args.bandwidth,
                              objective=objective)
        p=args.output/f'python_{objective}'
        p.mkdir(exist_ok=True)
        for name,val in [('W',result.W),('Z',result.Z)]:
            np.savetxt(p/f'{name}.csv',val,delimiter=',',fmt='%.17g')
        with open(p/'parameters.json','w') as f:
            json.dump({'tau':result.tau,'sigma2':result.sigma2,'objective':objective,
                       'log_tau_bounds':list(result.log_tau_bounds),
                       'bound_hit':bool(min(abs(np.log(result.tau)-b)
                                           for b in result.log_tau_bounds)<1e-3)},f,indent=2)
    with open(args.output/'input_metadata.json','w') as f:
        json.dump({'source_input':args.input.name,'m':int(Y.shape[0]),'n':int(Y.shape[1]),
                   'components':args.components,'bandwidth':args.bandwidth,'intercept_only':True,
                   'seed':args.seed,'scale':args.scale,'jitter':1e-8,
                   'data_status':'synthetic' if args.seed is not None or args.input.name=='demo_synthetic.npz' else 'user_supplied',
                   'sha256':{name:hashlib.sha256((args.output/f'{name}.csv').read_bytes()).hexdigest()
                             for name in ('Y','coords','K')}},f,indent=2)
    print('Wrote shared Y, coords, K and independently fitted Python outputs to',args.output)

if __name__=='__main__': main()
