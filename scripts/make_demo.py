#!/usr/bin/env python3
"""Create a deterministic synthetic SpatialPCA input file."""
import pathlib, sys, numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from spatialpca_verify.core import make_synthetic
x,y,k,w,z,coords=make_synthetic(seed=2026,n_side=4,m=12,d=2)
out=ROOT/'data'/'demo_synthetic.npz'
out.parent.mkdir(parents=True,exist_ok=True)
np.savez_compressed(out,Y=y,coords=coords)
print('Saved',out,'Y',y.shape,'coords',coords.shape)
