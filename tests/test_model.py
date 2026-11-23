import numpy as np
import pytest
from numpy.testing import assert_allclose
from spatialpca_verify.model import (gaussian_kernel, fit_spatialpca,
                                      standardize_genes, log_normalize_counts)
from spatialpca_verify.core import make_synthetic, fit_at_log_tau


def example():
    X,Y,K,W,Z,coords=make_synthetic(seed=2026, n_side=4, m=12,d=2)
    return Y, coords, K


def test_kernel_dimensions_and_symmetry():
    y,c,k=example()
    kk=gaussian_kernel(c,bandwidth=.19)
    assert kk.shape == (16,16)
    assert_allclose(kk, kk.T, atol=1e-13)
    assert_allclose(np.diag(kk),1.)
    assert gaussian_kernel(c[:3],c[:5],bandwidth=.19).shape == (3,5)


def test_full_rank_end_to_end_and_predict():
    y,c,k=example()
    fit=fit_spatialpca(y,c,n_components=2,bandwidth=.19)
    assert fit.W.shape==(12,2)
    assert fit.Z.shape==(2,16)
    assert fit.beta.shape==(1,12)
    assert_allclose(fit.W.T@fit.W,np.eye(2),atol=1e-10)
    assert fit.sigma2>0 and fit.tau>0 and np.isfinite(fit.objective_value)
    znew,ynew=fit.predict_locations(c[:3])
    assert znew.shape==(2,3) and ynew.shape==(12,3)
    assert_allclose(znew,fit.Z[:,:3],atol=1e-5)
    assert np.isfinite(fit.fitted_expression).all()


def test_multiple_covariates():
    y,c,k=example()
    xextra=c[:,0]
    fit=fit_spatialpca(y,c,n_components=2,bandwidth=.19,covariates=xextra)
    assert fit.X.shape==(16,2) and fit.beta.shape==(2,12)
    with pytest.raises(ValueError): fit.predict_locations(c[:2])
    z,yp=fit.predict_locations(c[:2],new_covariates=xextra[:2])
    assert_allclose(z,fit.Z[:,:2],atol=1e-5)


def test_lowrank_predict_and_covariance():
    y,c,k=example()
    fit=fit_spatialpca(y,c,n_components=2,bandwidth=.19,kernel_rank=8)
    assert fit.rank==8
    assert_allclose(fit.W.T@fit.W,np.eye(2),atol=1e-10)
    z, yp=fit.predict_locations(c[:4])
    assert_allclose(z,fit.Z[:,:4],atol=1e-5)


def test_objective_modes_and_preprocessing():
    y,c,k=example()
    f1=fit_spatialpca(y,c,n_components=2,bandwidth=.19,objective='profile')
    f2=fit_spatialpca(y,c,n_components=2,bandwidth=.19,objective='no_log_ablation')
    assert np.isfinite(f1.objective_value) and np.isfinite(f2.objective_value)
    z=standardize_genes(y)
    assert_allclose(z.mean(axis=1),0.,atol=1e-13)
    counts=np.array([[1.,2.,3.],[3.,4.,5.]])
    assert log_normalize_counts(counts).shape==(2,3)
    with pytest.raises(ValueError): fit_spatialpca(y,c,covariates=np.ones((16,1)))


def test_match_existing_dense_formula_at_same_tau():
    y,c,k=example()
    f=fit_spatialpca(y,c,n_components=2,bandwidth=.19)
    # Our model includes an added diagonal jitter of 1e-8; recompute the
    # reference evaluation with the identical kernel.
    x=np.ones((c.shape[0],1))
    ref=fit_at_log_tau(x,y,gaussian_kernel(c,bandwidth=.19)+1e-8*np.eye(len(c)),
                       2,np.log(f.tau),'profile')
    assert_allclose(f.sigma2,ref.sigma2,atol=1e-8)
    assert_allclose(f.objective_value,ref.objective_value,atol=1e-7)
    assert_allclose(f.W@f.Z,ref.W@ref.Z,atol=1e-7)


def test_scaling_invariance_profile_tau():
    y,c,k=example()
    a=fit_spatialpca(y,c,n_components=2,bandwidth=.19)
    b=fit_spatialpca(y*2,c,n_components=2,bandwidth=.19)
    assert_allclose(a.tau,b.tau,rtol=1e-5)
    assert_allclose(b.sigma2,4*a.sigma2,rtol=1e-5)
