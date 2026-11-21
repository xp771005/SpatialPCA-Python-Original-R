# SpatialPCA Python + Original R Runner

Private NJIT research project: an independent Python model-core reproduction with a runner for the authors' unmodified R functions. Repository visibility must remain private unless the owner explicitly approves publication.

## Project files

The complete source, derivation notes, tests, demo data and previous local results are in **SpatialPCA-Python-with-Original-R-Runner.zip**. Extract it before using the source paths below. The root Actions workflow extracts the same archive automatically.

~~~bash
unzip SpatialPCA-Python-with-Original-R-Runner.zip
cd SpatialPCA-Mathematical-Verification
python -m venv .venv
# Activate .venv for your operating system, then:
python -m pip install -r requirements.txt
python -m pytest -q
~~~

## GitHub Actions evidence

[Actions run #1](https://github.com/xp771005/SpatialPCA-Python-Original-R/actions/runs/44842777863) is the actual GitHub execution record. Use its final status, summary and downloadable artifact as evidence. Historical reports inside the archive describe the earlier local execution; they do not establish the outcome of this Actions run.

The workflow checks that the repository is private, checks out immutable upstream commit 1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66, runs 18 mathematical/model tests and nine matched-input R/Python experiments. It retains logs, matrices, input hashes, dependency versions and metrics in the original-r-python-comparison artifact for 30 days. Workflow permissions are contents: read.

To rerun, open Actions -> Original R vs Python model-core comparison -> Run workflow.

## Experimental design

Seeds 23, 41 and 2026, crossed with expression scales 0.5, 1 and 2. Each case uses 16 genes, 25 spatial positions, 2 components, an intercept-only design, a full-rank Gaussian kernel, bandwidth 0.19 and jitter 1e-8. R and Python receive identical Y, coordinates and K. Optimizer bounds are log(tau) in [-10,10].

| Quantity | Pass threshold |
|---|---:|
| tau relative gap | 1e-3 |
| sigma2 relative gap | 1e-5 |
| Maximum W subspace angle | 0.01 degrees |
| WZ relative gap | 1e-4 |

Subspaces and WZ are compared because latent factors may differ by signs or rotations.

## Corrections and scope

The earlier claim that the official R objective omits log(S) was incorrect and is withdrawn. The pinned upstream code explicitly includes log(S). Python no_log_ablation deliberately changes the objective and is only a hypothetical sensitivity control. The runner's unmatched parenthesis and incorrect SpatialPCs transpose were fixed before upload; these were local findings, not Actions failures. Upstream R source is not patched.

This verifies the synthetic full-rank, intercept-only model core. It does not establish replication of SCTransform/SPARK preprocessing, real tissue analyses, the upstream fast approximation, every published figure, or a complete version-specific audit of all paper appendices.

## Source navigation after extraction

| Path | Purpose |
|---|---|
| src/spatialpca_verify/model.py | Python model, posterior PCs and prediction |
| src/spatialpca_verify/core.py | Independent dense mathematical baseline |
| scripts/run_author_original_r.R | Executes unmodified original R functions |
| scripts/run_parity_suite.py | Nine cases, provenance, logs and report |
| scripts/compare_original_r_python.py | Numerical parity gate |
| docs/formula_audit.md | Derivation and corrected source audit |
| docs/chinese_notes.md | Chinese explanation and remaining scope |
| results/local_reference | Previous local execution evidence |

For a local original-R comparison, install R 4.4.3 and RSpectra, then run:

~~~bash
git clone https://github.com/shangll123/SpatialPCA.git upstream/SpatialPCA
git -C upstream/SpatialPCA checkout 1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66
python scripts/run_parity_suite.py --output results/new_parity_run
~~~

Use a new or empty output directory. The suite rejects stale output reuse and modified upstream R files.

SpatialPCA was developed by Lulu Shang and Xiang Zhou. Paper: https://www.nature.com/articles/s41467-022-34879-1 . Official GPL-3.0 source: https://github.com/shangll123/SpatialPCA . The archive does not bundle the upstream R source.
