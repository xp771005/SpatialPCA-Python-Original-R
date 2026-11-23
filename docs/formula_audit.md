# Mathematical checks (working notes)

## Notation

- `Y` is m genes by n spatial positions.
- `X` is n by q, full column rank, with q < n.
- `W` is m by d with `W.T @ W = I_d`.
- `Z` is d by n, with independent Gaussian spatial prior rows and a shared kernel `K`.
- `M = I_n - X (X.T X)^-1 X.T`, symmetric and idempotent.
- `F_tau = (M + tau^-1 K^-1)^-1` when `K` is positive definite.

## 1. Integral over B

Given `Y* = (Y-WZ).T`, the Gaussian residual sum of squares admits the orthogonal decomposition:

`||Y* - X B||_F² = ||M Y*||_F² + tr[(B - Bhat).T (X.T X) (B - Bhat)]`,

where `Bhat = (X.T X)^-1 X.T Y*`.

Integrating each of the **m independent column-coefficient vectors** `b_j ∈ R^q` yields

`∫ exp(-tr[(B-Bhat).T (X.T X) (B-Bhat)]/(2 sigma²)) dB`

`= (2 pi sigma²)^(qm/2) |X.T X|^(-m/2)`.

The exponent `-m/2` is required when integrating the complete coefficient matrix. A prior audit noted that a thesis appendix version appears to give `-1/2` in one displayed line. This matters to the *full expression*; for fixed X it is constant in W, sigma² and tau.

## 2. Integrating out Z and profiling sigma²

Assume `K` is SPD, tau > 0. Integrating Z leads to

`L ∝ (sigma²)^(-N/2) |I + tau M K|^(-d/2) exp[-S(tau)/(2 sigma²)]`,

where `N=m(n-q)`, `G_tau=Y M F_tau M Y.T`, and

`S(tau) = tr(Y M Y.T) - sum(top d eigenvalues of G_tau)`.

Set the derivative of `-N/2 log(sigma²) - S/(2 sigma²)` to zero:

`-N/(2 sigma²) + S/(2 sigma^4) = 0 => sigma²_hat = S/N`.

Substituting back (including constants) gives the tau-dependent objective

`J_profile(tau) = N/2 log S(tau) + d/2 log|I+tau M K| + const`.

## 3. Corrected comparison with the actual R source

The earlier missing-log claim was incorrect and is withdrawn. At pinned upstream commit
1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66, the final expression in
SpatialPCA_estimate_parameter() explicitly contains log(S).

For a full-rank kernel the R determinant expression differs from the reduced
profile determinant only by the tau-independent constant log|X.T X|.
Its G construction uses I minus a resolvent with a regression correction;
this is algebraically compatible with the full-rank projected dense formula.
The actual R/Python numerical comparison is in results/local_reference.

The optional Python no_log_ablation deliberately replaces log(S) by S.
This is a hypothetical sensitivity experiment, not evidence about the original R code.
Under Y -> cY, S -> c^2 S, so the true profile objective changes by N log(c),
leaving tau invariant. The no-log ablation generally changes shape.

## 4. Subspace comparison

Eigenvectors can differ by sign or by a joint rotation within a loading subspace. Use principal angles between `span(W_profile)` and `span(W_Rexpr)` as well as `W @ Z` reconstruction differences, rather than comparing W elements directly.

## 5. Kernel low rank caveat

If `K_r = U_r D_r U_r.T` with `r<n`, then `K_r^-1` does not exist. Replacing the prior covariance with this rank-deficient kernel gives an extension expressible as

`F_r = U_r (tau^-1 D_r^-1 + U_r.T M U_r)^-1 U_r.T`,

where only the strictly positive retained eigenvalues are inverted. This is not necessarily identical to approximating an intermediate inverse from a *full-rank* derivation by truncating its eigendecomposition. These are different approximations and must be labeled accordingly; a discrepancy alone is not proof of a software bug.

## Source links

- [Journal article](https://www.nature.com/articles/s41467-022-34879-1)
- [Shang thesis appendix](https://deepblue.lib.umich.edu/bitstream/2027.42/177856/1/shanglu_1.pdf)
- [Public R objective source](https://rdrr.io/github/shangll123/SpatialPCA/src/R/SpatialPCA_EstimateLoading.R)

These notes are an independent verification work-in-progress, not peer-reviewed findings.

## Verification boundary

This round checked the fixed upstream source and actual runtime, not every page of the journal supplement. Earlier suspected paper/thesis typographical errors still need version-specific page references before presentation as confirmed corrections.
