"""Independent, small-scale verification of SpatialPCA mathematical objectives."""
from .core import fit_at_log_tau, optimize_tau, make_synthetic, compare_fits
__all__ = ['fit_at_log_tau', 'optimize_tau', 'make_synthetic', 'compare_fits']
from .model import SpatialPCAFit, fit_spatialpca, gaussian_kernel, standardize_genes, log_normalize_counts
__all__ += ['SpatialPCAFit', 'fit_spatialpca', 'gaussian_kernel',
            'standardize_genes', 'log_normalize_counts']
