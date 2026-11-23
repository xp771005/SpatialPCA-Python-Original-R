# SpatialPCA Python + Original R

Python reproduction of the SpatialPCA model with mathematical derivations and numerical validation against the original R implementation.

## Code and mathematical note

| Path | Contents |
|---|---|
| [src/spatialpca_verify](src/spatialpca_verify) | Python model, spatial PCs, prediction, and dense mathematical reference |
| [scripts](scripts) | Original R runner, shared-input generation, numerical comparison, and nine-case suite |
| [tests](tests) | 18 mathematical, model, and comparison-gate tests |
| [docs/SpatialPCA_Short_Note.pdf](docs/SpatialPCA_Short_Note.pdf) | Mathematical note: notation issues, checks, and full derivation |
| [docs/SpatialPCA_Short_Note.tex](docs/SpatialPCA_Short_Note.tex) | Editable LaTeX source for the note |
| [docs/formula_audit.md](docs/formula_audit.md) | Working derivation and source-code audit |
| [data](data) | Synthetic demonstration input |
| [experiment_results.md](experiment_results.md) | Numerical results and verification scope |

The source is directly browsable in this repository. The ZIP is retained as a historical snapshot of the earlier package; its bundled README, workflow, and status reports are superseded by the expanded files on `main`.

## Verified execution

**[Actions #2 — successful R/Python validation](https://github.com/xp771005/SpatialPCA-Python-Original-R/actions/runs/44843893756)**: 18 tests and all nine shared-input comparisons passed. The runtime deprecation warning from the first run was fixed before this run.

| Quantity | Maximum error in Actions #2 | Pass threshold |
|---|---:|---:|
| tau relative gap | 2.35e-7 | 1e-3 |
| sigma2 relative gap | 1.51e-8 | 1e-5 |
| W subspace angle | 1.61e-6 degrees | 0.01 degrees |
| WZ relative gap | 3.17e-8 | 1e-4 |

**[Actions #3 — expanded-source public rerun](https://github.com/xp771005/SpatialPCA-Python-Original-R/actions/runs/45236586819)** also passed all 18 tests and all nine comparisons on November 23, 2026 (attempt 2, 7m 31s). Its maximum errors match Actions #2. Attempt 1 was manually canceled during R environment setup; it did not reach the parity suite.

The public-repository workflow now runs the expanded source directly, with read-only repository permissions. It does not extract the historical ZIP. To run it, open **Actions → Original R vs Python model-core comparison → Run workflow**. Matrices, metrics, per-case logs, dependency versions, input hashes, and R session information are uploaded as `original-r-python-comparison` for 30 days.

## Run Python locally

From the repository root, with Python 3.10 or newer:

```bash
python -m venv .venv
# Activate .venv for your operating system.
python -m pip install -r requirements.txt
python -m pytest -q
python scripts/run_model.py --input data/demo_synthetic.npz --output results/demo_fit_profile.npz --components 2 --bandwidth 0.19
```

For direct imports, add `src` to the Python module search path, for example `PYTHONPATH=src` on Linux/macOS. Inputs are normalized expression `Y` (genes × locations) and coordinates (locations × spatial dimensions).

## Compare against the original R implementation

Install R 4.4.3 and RSpectra, then run from the repository root:

```bash
Rscript -e 'install.packages("RSpectra", repos="https://cloud.r-project.org")'
git clone https://github.com/shangll123/SpatialPCA.git upstream/SpatialPCA
git -C upstream/SpatialPCA checkout 1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66
python scripts/run_parity_suite.py --output results/new_parity_run
```

Use a new or empty output directory. The suite rejects stale output reuse and modified upstream R files. The authors' original model-core functions are loaded without patching their source.

## Experimental design and scope

Seeds 23, 41, and 2026 are crossed with expression scales 0.5, 1, and 2. Each case uses 16 genes, 25 locations, 2 factors, an intercept-only design, and a full-rank Gaussian kernel (bandwidth 0.19, jitter 1e-8). R and Python receive identical Y, coordinates, and K. Optimization uses log(tau) in [-10, 10]. W is compared through principal angles and WZ because latent factors may differ by signs or rotations.

This validates the **synthetic, intercept-only, full-kernel model core**. It does not reproduce SCTransform/SPARK preprocessing, real tissue analyses, the authors' fast approximation, or every published biological result. The LaTeX note distinguishes dimension-based notation checks from the dissertation determinant issue; it does not claim that every journal-supplement equation has been audited.

## Correction and attribution

An earlier claim that the original R objective omitted log(S) was incorrect and is withdrawn. The pinned R source contains log(S). Python `no_log_ablation` deliberately changes that objective and is only a sensitivity control. Our runner's parenthesis and SpatialPCs-transpose errors were fixed before the recorded Actions runs; the authors' R source was not modified.

SpatialPCA was developed by Lulu Shang and Xiang Zhou. [Paper](https://www.nature.com/articles/s41467-022-34879-1) · [Original R repository](https://github.com/shangll123/SpatialPCA). The upstream repository uses GPL-3.0; upstream R source is fetched at runtime and is not bundled here.
