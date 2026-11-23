#!/usr/bin/env python3
"""Measure actual R/Python agreement, with a required profile-parity gate."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import subspace_angles

TOLERANCES = {'tau_relative_gap':1e-3, 'sigma2_relative_gap':1e-5,
              'max_subspace_angle_deg':1e-2, 'reconstruction_relative_gap':1e-4}


def load_result(directory):
    params=json.loads((directory/'parameters.json').read_text())
    W=np.loadtxt(directory/'W.csv',delimiter=',',ndmin=2)
    Z=np.loadtxt(directory/'Z.csv',delimiter=',',ndmin=2)
    if (not np.isfinite(W).all() or not np.isfinite(Z).all()
            or not np.isfinite(params['tau']) or params['tau']<=0
            or not np.isfinite(params['sigma2']) or params['sigma2']<=0):
        raise ValueError(f'Invalid fit in {directory}')
    if not np.allclose(W.T@W,np.eye(W.shape[1]),atol=1e-7):
        raise ValueError(f'Non-orthogonal loadings in {directory}')
    if W.shape[1]!=Z.shape[0]:
        raise ValueError(f'Incompatible factor dimensions in {directory}')
    return params,W,Z


def compare(base):
    R,Wr,Zr=load_result(base/'author_R_original')
    if R.get('original_R_executed') is not True:
        raise ValueError('Original R execution receipt missing')
    rows=[]
    for mode in ('profile','no_log_ablation'):
        P,Wp,Zp=load_result(base/f'python_{mode}')
        if Wr.shape!=Wp.shape or Zr.shape!=Zp.shape:
            raise ValueError('R/Python shape mismatch')
        Qr=np.linalg.qr(Wr)[0]; Qp=np.linalg.qr(Wp)[0]
        row={'python_mode':mode,'R_tau':float(R['tau']),'Python_tau':float(P['tau']),
             'R_sigma2':float(R['sigma2']),'Python_sigma2':float(P['sigma2']),
             'tau_relative_gap':float(abs(R['tau']-P['tau'])/R['tau']),
             'sigma2_relative_gap':float(abs(R['sigma2']-P['sigma2'])/R['sigma2']),
             'max_subspace_angle_deg':float(np.max(subspace_angles(Qr,Qp))*180/np.pi),
             'subspace_projector_gap_fro':float(np.linalg.norm(Qr@Qr.T-Qp@Qp.T,'fro')),
             'reconstruction_relative_gap':float(np.linalg.norm(Wr@Zr-Wp@Zp,'fro')/
                                                  max(np.linalg.norm(Wr@Zr,'fro'),1e-12)),
             'R_bound_hit':bool(abs(np.log(R['tau']))>9.999),
             'Python_bound_hit':P['bound_hit']}
        row['within_profile_tolerances']=all(row[k]<=v for k,v in TOLERANCES.items())
        rows.append(row)
    report={'original_R_executed':True,'tolerances':TOLERANCES,
            'profile_parity_passed':rows[0]['within_profile_tolerances'],'rows':rows,
            'note':'no_log_ablation is deliberately altered Python, not the original R objective.'}
    (base/'original_r_comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--dir',type=Path,default=Path('results/original_r_parity'))
    a=p.parse_args()
    report=compare(a.dir)
    print(json.dumps(report,indent=2))
    if not report['profile_parity_passed']:
        raise SystemExit('Original R/Python profile parity outside tolerances; inspect the report.')


if __name__=='__main__': main()
