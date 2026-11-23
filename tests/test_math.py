import sys
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from scipy.integrate import quad

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from spatialpca_verify import make_synthetic, optimize_tau, fit_at_log_tau


def test_matrix_dimensions_and_M_projection():
    X,Y,K,*_=make_synthetic(seed=7)
    n=X.shape[0]
    M=np.eye(n)-X@solve(X.T@X,X.T,assume_a='pos')
    assert np.allclose(M@X,0,atol=1e-11)
    assert np.allclose(M@M,M,atol=1e-11)


def test_integral_over_B_has_determinant_power_m():
    X=np.ones((3,1));Ystar=np.array([[4.,2.],[6.,3.],[8.,7.]])
    sig2=1.7;n,m=Ystar.shape;q=X.shape[1]
    M=np.eye(n)-X@solve(X.T@X,X.T)
    sse=np.linalg.norm(M@Ystar,'fro')**2
    analytic=(2*np.pi*sig2)**(-m*(n-q)/2)*np.linalg.det(X.T@X)**(-m/2)*np.exp(-sse/(2*sig2))
    integrals=[]
    for j in range(m):
        f=lambda b:(2*np.pi*sig2)**(-n/2)*np.exp(-np.sum((Ystar[:,j]-b)**2)/(2*sig2))
        integrals.append(quad(f,-np.inf,np.inf)[0])
    assert np.isclose(np.prod(integrals),analytic,rtol=1e-8)


def test_f_identity_full_rank():
    X,Y,K,*_=make_synthetic(seed=8)
    n=X.shape[0]
    M=np.eye(n)-X@solve(X.T@X,X.T)
    tau=1.3
    lhs=solve(M+solve(K,np.eye(n))/tau,np.eye(n))
    rhs=solve(np.eye(n)+tau*K@M,tau*K)
    assert np.allclose(lhs,rhs,atol=1e-7,rtol=1e-7)


def test_profile_objective_invariant_to_scale_at_fixed_tau():
    X,Y,K,*_=make_synthetic(seed=11,scale=1.0)
    _,Ys,_,*_=make_synthetic(seed=11,scale=2.0)
    a=fit_at_log_tau(X,Y,K,2,np.log(5),'profile')
    b=fit_at_log_tau(X,Ys,K,2,np.log(5),'profile')
    n=Y.shape[1];m=Y.shape[0];q=X.shape[1]
    # Difference N/2*log(4), independent of tau.
    assert np.isclose(b.objective_value-a.objective_value,m*(n-q)*np.log(2),rtol=1e-8)


def test_different_objectives_and_scale_behavior():
    X,Y,K,*_=make_synthetic(seed=2026)
    a=optimize_tau(X,Y,K,2,'profile')
    b=optimize_tau(X,Y,K,2,'no_log_ablation')
    assert np.isclose(a.W.T@a.W,np.eye(2)).all()
    assert abs(np.log(a.tau/b.tau))>1e-2
    _,Y2,_,*_=make_synthetic(seed=2026,scale=0.5)
    a2=optimize_tau(X,Y2,K,2,'profile')
    b2=optimize_tau(X,Y2,K,2,'no_log_ablation')
    assert np.isclose(a.tau,a2.tau,rtol=3e-4)
    assert abs(np.log(b.tau/b2.tau))>0.1


def test_lowrank_prior_vs_truncated_expression_not_identical():
    U=np.array([[1.],[0.],[0.]])
    D=np.array([[2.]])
    tau=1.
    M=np.eye(3)
    exact=U@solve(np.linalg.inv(tau*D)+U.T@M@U,U.T)
    truncated=np.eye(3)-U@solve(tau*D+U.T@U,U.T)
    assert not np.allclose(exact,truncated)


def test_A7_determinant_identity():
    X,Y,K,*_=make_synthetic(seed=4)
    n=X.shape[0];tau=1.7
    M=np.eye(n)-X@solve(X.T@X,X.T)
    a=np.linalg.slogdet(np.eye(n)+tau*M@K)[1]
    B=np.eye(n)+tau*K
    H=X.T@solve(B,X)
    b=(np.linalg.slogdet(B)[1]+np.linalg.slogdet(H)[1]
       -np.linalg.slogdet(X.T@X)[1])
    assert np.isclose(a,b,atol=2e-9)
