"""Small dense-matrix formula audit for SpatialPCA.

The implementation below was written independently from the model equations.
It does NOT implement or run the authors' full R pipeline. In particular,
we compare the profiled negative log likelihood with a deliberately altered
no-log objective. The altered objective is an ablation, NOT the author's R
objective: the pinned official source includes log(S).

Dimensions: Y (m x n) genes x positions; X (n x q) covariates;
K (n x n) positive-definite spatial kernel; W (m x d); Z (d x n).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh, solve, subspace_angles
from scipy.optimize import minimize_scalar


@dataclass(frozen=True)
class Fit:
    objective_name: str
    log_tau: float
    tau: float
    sigma2: float
    S: float
    logdet: float
    objective_value: float
    W: np.ndarray
    Z: np.ndarray
    lowrank_reconstruction: np.ndarray


def make_synthetic(seed: int = 2026, scale: float = 1.0,
                   n_side: int = 5, m: int = 16, d: int = 2,
                   noise_variance: float = 0.35, tau_true: float = 3.5):
    """Sample Y from the intercept-only, full-kernel SpatialPCA model.

    The same base random realization is used for any `scale` with a fixed seed:
    scale=2 doubles every entry of Y (both signal and noise).
    """
    if n_side < 2 or m <= d or scale <= 0 or noise_variance <= 0 or tau_true <= 0:
        raise ValueError('Invalid simulation parameters')
    rng = np.random.default_rng(seed)
    grid = np.linspace(0.0, 1.0, n_side)
    coords = np.array([(x, y) for x in grid for y in grid], dtype=float)
    n = coords.shape[0]
    dist2 = np.sum((coords[:, None, :] - coords[None, :, :])**2, axis=2)
    K = np.exp(-dist2 / 0.19) + 1e-6 * np.eye(n)
    X = np.ones((n, 1), dtype=float)
    W_true, _ = np.linalg.qr(rng.standard_normal((m, d)))
    # Each latent factor has covariance sigma2 * tau * K.
    Z_true = (rng.standard_normal((d, n)) @ np.linalg.cholesky(K).T
              * np.sqrt(noise_variance * tau_true))
    gene_intercepts = rng.normal(0, 0.7, size=(m, 1))
    E = np.sqrt(noise_variance) * rng.standard_normal((m, n))
    Y = (gene_intercepts + W_true @ Z_true + E) * scale
    return X, Y, K, W_true, Z_true, coords


def _validate(X, Y, K, d):
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    K = np.asarray(K, float)
    if X.ndim != 2 or Y.ndim != 2 or K.ndim != 2:
        raise ValueError('X, Y, K must be two-dimensional')
    m, n = Y.shape
    if X.shape[0] != n or K.shape != (n, n):
        raise ValueError('Incompatible dimensions: Y (m,n), X (n,q), K (n,n)')
    q = X.shape[1]
    if not (0 < q < n and 1 <= d <= m):
        raise ValueError('Requires 0<q<n and 1<=d<=m')
    if np.linalg.matrix_rank(X) != q:
        raise ValueError('X must have full column rank')
    if not np.allclose(K, K.T, atol=1e-10):
        raise ValueError('K must be symmetric')
    if np.linalg.eigvalsh(K)[0] <= 0:
        raise ValueError('K must be positive definite for this full-rank audit')
    return X, Y, K


def _terms(X: np.ndarray, Y: np.ndarray, K: np.ndarray,
           d: int, log_tau: float):
    tau = float(np.exp(log_tau))
    m, n = Y.shape
    q = X.shape[1]
    P = X @ solve(X.T @ X, X.T, assume_a='pos')
    M = np.eye(n) - P
    M = (M + M.T) / 2
    # Full-rank identity: F = (M + tau^{-1} K^{-1})^{-1}
    # = (I + tau K M)^{-1} tau K, avoiding K^{-1} explicitly.
    F = solve(np.eye(n) + tau * K @ M, tau * K, assume_a='gen')
    F = (F + F.T) / 2  # round-off symmetry only
    G = Y @ M @ F @ M @ Y.T
    G = (G + G.T) / 2
    eigvals, eigvecs = eigh(G)
    W = eigvecs[:, -d:]
    S = float(np.trace(Y @ M @ Y.T) - np.sum(eigvals[-d:]))
    if not (S > 0):
        raise FloatingPointError(f'Residual sum term must be positive, got {S}')
    # determinant(I+tau M K) is positive; compute via its SPD
    # restriction to the orthogonal complement of X.
    # This also avoids the nonsymmetric MK numerical determinant.
    Q_full, _ = np.linalg.qr(X, mode='complete')
    Q_perp = Q_full[:, q:]
    red = np.eye(n - q) + tau * (Q_perp.T @ K @ Q_perp)
    sign, logdet = np.linalg.slogdet(red)
    if sign <= 0:
        raise FloatingPointError('Expected a positive determinant')
    Z = W.T @ Y @ M @ F
    return dict(tau=tau, S=S, logdet=float(logdet), W=W, Z=Z,
                reconstruction=W @ Z, N=m*(n-q))


def fit_at_log_tau(X: np.ndarray, Y: np.ndarray, K: np.ndarray,
                   d: int, log_tau: float, objective: str = 'profile') -> Fit:
    """Evaluate parameters at tau for profiled likelihood or a deliberate no-log ablation.

    profile: N/2 log S(tau) + d/2 log det(I+tau M K)
    no_log_ablation: N/2 S(tau) + d/2 log det(I+tau M K)

    The no-log ablation is not the original R implementation.
    """
    X, Y, K = _validate(X, Y, K, d)
    if objective not in ('profile', 'no_log_ablation'):
        raise ValueError('Unknown objective')
    t = _terms(X, Y, K, d, log_tau)
    residual = np.log(t['S']) if objective == 'profile' else t['S']
    val = t['N'] * residual / 2 + d * t['logdet'] / 2
    return Fit(objective_name=objective, log_tau=log_tau, tau=t['tau'],
               sigma2=t['S'] / t['N'], S=t['S'], logdet=t['logdet'],
               objective_value=float(val), W=t['W'], Z=t['Z'],
               lowrank_reconstruction=t['reconstruction'])


def optimize_tau(X: np.ndarray, Y: np.ndarray, K: np.ndarray, d: int = 2,
                 objective: str = 'profile', bounds=(-10.0, 10.0)) -> Fit:
    """One-dimensional bounded optimization in log(tau), as in the R interface.

    Uses a grid of starts plus local bounded refinement to reduce the risk of
    mistaking a local minimum for the best 1D candidate.
    """
    X, Y, K = _validate(X, Y, K, d)
    grid = np.linspace(bounds[0], bounds[1], 81)
    vals = np.array([fit_at_log_tau(X,Y,K,d,x,objective).objective_value for x in grid])
    candidates = [(float(vals[i]), float(grid[i])) for i in [0, len(grid)-1]]
    for i in range(1, len(grid)-1):
        if vals[i] <= vals[i-1] and vals[i] <= vals[i+1]:
            result = minimize_scalar(
                lambda z: fit_at_log_tau(X,Y,K,d,z,objective).objective_value,
                bounds=(grid[i-1], grid[i+1]), method='bounded',
                options={'xatol': 1e-9})
            candidates.append((result.fun, float(result.x)))
    _, best = min(candidates)
    return fit_at_log_tau(X,Y,K,d,best,objective)


def compare_fits(profile: Fit, no_log_ablation: Fit):
    """Metrics invariant to sign changes and joint orthogonal rotations of W,Z."""
    angles = np.rad2deg(subspace_angles(profile.W, no_log_ablation.W))
    rel_recon = float(np.linalg.norm(profile.lowrank_reconstruction-no_log_ablation.lowrank_reconstruction,
                                    'fro') / np.linalg.norm(profile.lowrank_reconstruction, 'fro'))
    return {'largest_subspace_angle_degrees': float(np.max(angles)),
            'relative_reconstruction_difference': rel_recon,
            'tau_ratio_ablation_over_profile': float(no_log_ablation.tau/profile.tau)}
