# Verification status

- Repository: public; source expanded under src/, scripts/, tests/, docs/ and data/.
- Original R source: unchanged at 1d4dbb1fe6f23b77bb7ed805ce2de8bbc017ef66.
- Recorded GitHub Actions #2: [44843893756](https://github.com/xp771005/SpatialPCA-Python-Original-R/actions/runs/44843893756), success; 18 tests and 9/9 parity cases passed.
- Node runtime deprecation warning from #1: fixed before #2.
- Current workflow: uses expanded root source, read-only permissions, no ZIP extraction or private-repository gate.
- Mathematical note: [PDF](../docs/SpatialPCA_Short_Note.pdf) and [LaTeX](../docs/SpatialPCA_Short_Note.tex).
- Scope: synthetic full-kernel, intercept-only model core; fast R mode and real-data preprocessing are not covered.
- Missing-log claim: withdrawn; no_log_ablation is hypothetical.
- Historical local evidence: local_reference/. Current cross-platform report: ../experiment_results.md.
