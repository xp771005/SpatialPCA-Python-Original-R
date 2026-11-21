# SpatialPCA Python + Original R Runner

Independent Python model-core reproduction of SpatialPCA, with a runner for the authors' unmodified R functions.

**Verified locally:** nine synthetic R/Python comparisons passed. **GitHub upload and GitHub Actions execution are pending GitHub authentication.** Local execution must not be reported as a GitHub Actions run.

## Correction to the earlier package

The earlier report incorrectly described the public R objective as omitting log(S). The immutable upstream source at commit 1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66 explicitly includes log(S). That claim is withdrawn. The former Python mode r_expression has been renamed **no_log_ablation** because it deliberately alters the likelihood and does not represent the original R algorithm.

The original R runner also had an unmatched parenthesis and incorrectly transposed SpatialPCs. Both are fixed. The upstream source is not patched.

## What has been verified

- Original upstream R model-core functions execute on the same Y and K as Python.
- Seeds: 23, 41, 2026. Scales of Y: 0.5, 1, 2. Each case has 16 genes, 25 spots and 2 components.
- Full-rank Gaussian kernel; bandwidth 0.19, diagonal jitter 1e-8; intercept only; log(tau) bounds [-10, 10].
- Compare tau, sigma2, loading subspaces, and W @ Z, which avoids arbitrary signs/rotations of latent factors.
- Numerical tolerances: relative tau gap <= 1e-3; relative sigma2 gap <= 1e-5; maximum loading-subspace angle <= 0.01 degrees; relative WZ gap <= 1e-4.
- Inputs, outputs, SHA-256 hashes, upstream revision, dependency versions, sessionInfo, and per-case logs are retained in results/local_reference.
- The mathematical and model tests pass. See results/verification_status.md for final status.

This is a **synthetic model-core comparison**. It does not reproduce SCTransform/SPARK preprocessing, real biological data analyses, every paper figure, or the upstream accelerated approximation.

## Install and use Python

Python 3.10+:

~~~bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
PYTHONPATH=src python -m pytest -q
python scripts/run_model.py --input data/demo_synthetic.npz --output results/demo_fit_profile.npz --components 2 --bandwidth 0.19
~~~

The default **profile** objective follows the likelihood. **no_log_ablation** is a hypothetical sensitivity control.

~~~python
from spatialpca_verify.model import fit_spatialpca
fit = fit_spatialpca(Y, coords, n_components=2, bandwidth=0.19)
# Y: genes x spots; coords: spots x dimensions
print(fit.tau, fit.sigma2, fit.W.shape, fit.Z.shape)
Z_new, Y_new = fit.predict_locations(coords_new)
~~~

Optional extra covariates, truncated-kernel fitting, and GP prediction remain available. Truncated-kernel Python support is not a claim of parity with the authors' fast R variant.

## Run the original R comparison locally

Install R 4.4.3 and RSpectra. Then:

~~~bash
git clone https://github.com/shangll123/SpatialPCA.git upstream/SpatialPCA
git -C upstream/SpatialPCA checkout 1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66
python scripts/run_parity_suite.py --output results/new_parity_run
~~~

Use a new or empty output directory. This prevents old outputs being mistaken for new R execution. The suite checks the upstream commit and source hashes, parses the R runner, executes all nine cases, retains failure logs, and exits unsuccessfully if profile parity fails. It never modifies the original R functions.

## GitHub Actions

Upload to a **private** repository, then open Actions -> Original R vs Python model-core comparison -> Run workflow.

The workflow refuses to run in a public repository, has read-only repository permissions, pins the upstream revision, runs Python tests, and executes the nine comparisons. It uploads matrices, metrics and logs even after a failure, and writes a Markdown result table to the Actions summary. Artifacts are retained for 30 days.

No paid plan or billing change is configured. No GitHub Pages deployment, release or repository visibility change is configured.

## Files

| Path | Purpose |
|---|---|
| src/spatialpca_verify/model.py | Python model, posterior PCs, covariates, prediction |
| src/spatialpca_verify/core.py | Independent dense-form mathematical baseline |
| scripts/run_author_original_r.R | Loads and executes unmodified upstream R model core |
| scripts/prepare_original_r_input.py | Shared input/kernel and Python fits |
| scripts/compare_original_r_python.py | Numerical comparison and pass/fail gate |
| scripts/run_parity_suite.py | Nine-case suite, provenance, logs and report |
| results/local_reference/experiment_results.md | Actual local R/Python results |
| docs/formula_audit.md | Derivation and corrected R-source audit |
| docs/chinese_notes.md | Chinese explanation and remaining scope |

## Attribution

SpatialPCA was developed by Lulu Shang and Xiang Zhou. [Paper](https://www.nature.com/articles/s41467-022-34879-1); [official R source](https://github.com/shangll123/SpatialPCA). The upstream repository is GPL-3.0. No upstream R source is bundled in this deliverable. This repository's Python code was independently written from the model equations. Repository publication has not been authorized.
