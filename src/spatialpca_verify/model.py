"""Independent, usable SpatialPCA *model-core* implementation.

Reference: Shang & Zhou, Nature Communications (2022),
https://doi.org/10.1038/s41467-022-34879-1

This is an original mathematical reimplementation, not a line-by-line R port.
Y: genes x observed positions; coordinates: observed positions x spatial dimensions.
X: observed positions x covariates (intercept added by default).

Two objectives:
    profile: mathematically profiled negative log likelihood;
    no_log_ablation: deliberately omit log(S) as a sensitivity ablation.
                     This is NOT the author's R objective: the pinned upstream
                     source includes log(S).

Input is already normalized expression, optionally standardized by gene.
The authors' full preprocessing/feature selection uses Seurat SCTransform/SPARK
and is NOT reproduced here.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.linalg import eigh, solve
from scipy.optimize import minimize_scalar
from scipy.spatial.distance import cdist


def gaussian_kernel(coords_a: np.ndarray, coords_b: np.ndarray | None = None,
                    bandwidth: float = 1.0) -> np.ndarray:
    """exp(-squared Euclidean distance/bandwidth), matching article convention."""
    a = np.asarray(coords_a, dtype=float)
    b = a if coords_b is None else np.asarray(coords_b, dtype=float)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1]:
        raise ValueError('Coordinates must be 2D with matching coordinate dimensions')
    if a.shape[0] == 0 or b.shape[0] == 0 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Coordinates must be nonempty and finite')
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError('bandwidth must be positive')
    return np.exp(-cdist(a, b, metric='sqeuclidean') / bandwidth)


def standardize_genes(Y: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Center and scale each gene, after external library normalization if needed."""
    Y = np.asarray(Y, float)
    if Y.ndim != 2 or not np.isfinite(Y).all():
        raise ValueError('Y must be a finite 2D array')
    means = Y.mean(axis=1, keepdims=True)
    std = Y.std(axis=1, keepdims=True, ddof=1)
    if np.any(std < eps):
        raise ValueError('Constant genes must be removed before standardization')
    return (Y - means) / std


def log_normalize_counts(counts: np.ndarray, scale_factor: float = 1e4) -> np.ndarray:
    """Convenience library-size/log1p normalization; NOT Seurat SCTransform.

    Use for a first demonstration only. For R parity, import the same
    normalized expression from R/Seurat instead.
    """
    counts = np.asarray(counts, float)
    if counts.ndim != 2 or np.any(counts < 0) or not np.isfinite(counts).all():
        raise ValueError('counts must be a 2D nonnegative finite matrix')
    library_size = counts.sum(axis=0, keepdims=True)
    if np.any(library_size <= 0):
        raise ValueError('Remove empty spots before library-size normalization')
    return np.log1p((counts / library_size) * scale_factor)


@dataclass
class SpatialPCAFit:
    W: np.ndarray
    Z: np.ndarray
    beta: np.ndarray
    tau: float
    sigma2: float
    objective_value: float
    residual_sum_S: float
    objective: str
    rank: int
    bandwidth: float
    coords: np.ndarray
    X: np.ndarray
    gene_expression: np.ndarray
    basis: np.ndarray
    eigenvectors: np.ndarray
    eigenvalues: np.ndarray
    posterior_coeff: np.ndarray
    n_components: int
    log_tau_bounds: tuple[float, float]

    @property
    def fitted_expression(self) -> np.ndarray:
        return (self.X @ self.beta).T + self.W @ self.Z

    def predict_locations(self, new_coords: np.ndarray, new_covariates: Optional[np.ndarray] = None):
        """Posterior mean spatial PCs and expression at new coordinates.

        For truncated kernel, use Nystrom cross-kernel extension. Includes
        covariate contribution if new covariates are given and were used in fit.
        Returns (Z_new: d x n_new, Y_new: m x n_new).
        """
        new_coords = np.asarray(new_coords, float)
        if new_coords.ndim != 2 or new_coords.shape[1] != self.coords.shape[1]:
            raise ValueError('new_coords must be n_new x coordinate_dimension')
        cross = gaussian_kernel(new_coords, self.coords, self.bandwidth)
        # If K ~= U diag(eigvals) U^T, the Nystrom basis at a new spot is
        # K(new, train) U diag(1/sqrt(eigvals)).
        basis_new = (cross @ self.eigenvectors) / np.sqrt(self.eigenvalues)[None, :]
        Znew = self.posterior_coeff @ basis_new.T
        nc = new_coords.shape[0]
        q = self.X.shape[1]
        if new_covariates is None:
            if q != 1:
                raise ValueError('Training included additional covariates; supply new_covariates')
            Xnew = np.ones((nc, 1))
        else:
            cov = np.asarray(new_covariates, float)
            if cov.ndim == 1:
                cov = cov[:, None]
            if cov.shape != (nc, q - 1) or not np.isfinite(cov).all():
                raise ValueError(f'new_covariates must have shape {(nc, q - 1)}')
            Xnew = np.column_stack([np.ones(nc), cov])
        Ynew = (Xnew @ self.beta).T + self.W @ Znew
        return Znew, Ynew


def fit_spatialpca(Y: np.ndarray, coords: np.ndarray, n_components: int = 2,
                   bandwidth: float = 0.3, covariates: np.ndarray | None = None,
                   kernel_rank: int | None = None, objective: str = 'profile',
                   log_tau_bounds: tuple[float, float] = (-10., 10.),
                   standardize: bool = False, jitter: float = 1e-8) -> SpatialPCAFit:
    """Fit the full-rank or truncated-kernel SpatialPCA model.

    Dense matrices are used for pedagogical traceability. Full rank is
    intended for small/moderate n; kernel_rank=r is a reduced-rank GP model.
    The optimization searches a grid followed by bounded scalar refinements.
    """
    Y = np.array(Y, dtype=float, copy=True)
    coords = np.asarray(coords, float)
    if Y.ndim != 2 or coords.ndim != 2 or Y.shape[1] != coords.shape[0]:
        raise ValueError('Y: genes x positions; coords: positions x dimensions')
    if not np.isfinite(Y).all():
        raise ValueError('Expression must be finite')
    m, n = Y.shape
    if n < 3 or m < 2 or n_components < 1 or n_components > m:
        raise ValueError('Require at least 2 genes, 3 spots, and 1 <= components <= genes')
    if objective not in ('profile', 'no_log_ablation'):
        raise ValueError('objective must be profile or no_log_ablation')
    if standardize:
        Y = standardize_genes(Y)
    if covariates is None:
        X = np.ones((n, 1))
    else:
        cov = np.asarray(covariates, float)
        if cov.ndim == 1:
            cov = cov[:, None]
        if cov.ndim != 2 or cov.shape[0] != n or not np.isfinite(cov).all():
            raise ValueError('covariates must be n_spots x q_extra')
        X = np.column_stack([np.ones(n), cov])
    q = X.shape[1]
    if q >= n or np.linalg.matrix_rank(X) != q:
        raise ValueError('Design matrix with intercept must be full column rank and q<n')
    if kernel_rank is not None and (not isinstance(kernel_rank, (int, np.integer))
                                    or not 1 <= kernel_rank <= n):
        raise ValueError('kernel_rank must be integer in [1, n]')
    if not (0 <= jitter < 1):
        raise ValueError('jitter must be in [0,1)')
    if len(log_tau_bounds) != 2 or not log_tau_bounds[0] < log_tau_bounds[1]:
        raise ValueError('Invalid bounds')

    kernel = gaussian_kernel(coords, bandwidth=bandwidth)
    kernel = kernel + jitter * np.eye(n)
    eigval, eigvec = eigh(kernel)
    ix = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[ix], eigvec[:, ix]
    r = n if kernel_rank is None else kernel_rank
    eigval, eigvec = eigval[:r], eigvec[:, :r]
    if eigval[-1] <= 0:
        raise ValueError('Selected kernel eigenvalues must be positive. Increase jitter or reduce rank')
    basis = eigvec * np.sqrt(eigval)[None, :]
    # QR is numerically stable for projectors even with non-orthogonal columns X.
    Q, _ = np.linalg.qr(X, mode='reduced')
    M = np.eye(n) - Q @ Q.T
    M = (M + M.T) / 2
    R = Y @ M
    RB = R @ basis
    BtMB = basis.T @ M @ basis
    BtMB = (BtMB + BtMB.T) / 2
    total_ss = float(np.sum(R * R))
    N = m * (n - q)
    eye = np.eye(r)

    def evaluate(log_tau: float, all_terms=False):
        tau = float(np.exp(log_tau))
        H = eye + tau * BtMB
        # H is SPD; Cholesky solve (no explicit inverse).
        RBS = solve(H, RB.T, assume_a='pos')
        G = tau * RB @ RBS
        G = (G + G.T) / 2
        ev, U = eigh(G)
        W = U[:, -n_components:]
        S = float(total_ss - np.sum(ev[-n_components:]))
        if S <= 0:
            raise FloatingPointError(f'Non-positive S={S}; check input/model dimension')
        logdet = float(np.linalg.slogdet(H)[1])
        obj = (N/2) * (np.log(S) if objective == 'profile' else S) + n_components/2 * logdet
        if not all_terms:
            return float(obj)
        # Posterior mean for latent coefficients v in Z = v @ basis.T
        # E[v|Y] = tau * W.T @ (Y M basis) @ H^-1.
        posterior_coeff = tau * solve(H, RB.T @ W, assume_a='pos').T
        Z = posterior_coeff @ basis.T
        beta = np.linalg.lstsq(X, (Y - W @ Z).T, rcond=None)[0]
        return float(obj), S, tau, S/N, W, Z, beta, posterior_coeff

    grid = np.linspace(log_tau_bounds[0], log_tau_bounds[1], 61)
    vals = np.array([evaluate(x) for x in grid])
    candidates = [(float(vals[i]), float(grid[i])) for i in (0, len(grid)-1)]
    for i in range(1, len(grid) - 1):
        if vals[i] <= vals[i - 1] and vals[i] <= vals[i + 1]:
            fit = minimize_scalar(evaluate, bounds=(grid[i-1],grid[i+1]), method='bounded',
                                  options={'xatol': 1e-9})
            candidates.append((float(fit.fun), float(fit.x)))
    _, logtau = min(candidates)
    obj, S, tau, sigma2, W, Z, beta, postcoef = evaluate(logtau, True)
    return SpatialPCAFit(W=W, Z=Z, beta=beta, tau=tau, sigma2=sigma2,
                         objective_value=obj, residual_sum_S=S, objective=objective,
                         rank=r, bandwidth=bandwidth, coords=np.array(coords,copy=True),
                         X=X, gene_expression=Y, basis=basis,
                         eigenvectors=eigvec, eigenvalues=eigval,
                         posterior_coeff=postcoef, n_components=n_components,
                         log_tau_bounds=log_tau_bounds)
