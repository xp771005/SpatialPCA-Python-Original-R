"""Test the numerical gate using explicit fixtures, not original R execution."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_original_r_python import compare


def fixture(base, perturb=False):
    rng=np.random.default_rng(19)
    W=np.linalg.qr(rng.normal(size=(8,2)))[0]
    Z=rng.normal(size=(2,7))
    rotation=np.array([[0.,-1.],[1.,0.]])
    for name,w,z in [('author_R_original',W,Z),
                      ('python_profile',W@rotation,rotation.T@Z),
                      ('python_no_log_ablation',W@rotation,rotation.T@Z)]:
        path=base/name;path.mkdir()
        if perturb and name=='python_profile':z=z+0.1
        np.savetxt(path/'W.csv',w,delimiter=',')
        np.savetxt(path/'Z.csv',z,delimiter=',')
        (path/'parameters.json').write_text(json.dumps({
            'tau':3.,'sigma2':.4,'bound_hit':False,'original_R_executed':True}))


def test_joint_rotation_preserves_comparison(tmp_path):
    fixture(tmp_path)
    assert compare(tmp_path)['profile_parity_passed']


def test_reconstruction_mismatch_fails_gate(tmp_path):
    fixture(tmp_path,perturb=True)
    assert not compare(tmp_path)['profile_parity_passed']


def test_missing_r_execution_receipt_rejected(tmp_path):
    fixture(tmp_path)
    p=tmp_path/'author_R_original/parameters.json'
    params=json.loads(p.read_text());params.pop('original_R_executed')
    p.write_text(json.dumps(params))
    with pytest.raises(ValueError,match='execution receipt'):compare(tmp_path)


def test_nonfinite_output_rejected(tmp_path):
    fixture(tmp_path)
    p=tmp_path/'author_R_original/Z.csv'
    Z=np.loadtxt(p,delimiter=',');Z[0,0]=np.nan
    np.savetxt(p,Z,delimiter=',')
    with pytest.raises(ValueError,match='Invalid fit'):compare(tmp_path)
