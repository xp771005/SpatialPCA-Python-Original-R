#!/usr/bin/env Rscript
# Execute ORIGINAL R MODEL-CORE SOURCE as written by Shang & Zhou.
# Avoids installation of the full Seurat/SPARK preprocessing stack by providing the
# identical pre-normalized Y and kernel K that were exported for Python.
# The official .R source files are not modified, copied into, or bundled with this repo.
# Usage: Rscript scripts/run_author_original_r.R <upstream_directory> <parity_directory> [components]
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<2L) stop('Usage: Rscript scripts/run_author_original_r.R upstream_repo parity_directory [components]')
upstream <- normalizePath(args[[1]],mustWork=TRUE)
workdir <- normalizePath(args[[2]],mustWork=TRUE)
components <- if(length(args)>=3L) as.integer(args[[3]]) else 2L

required <- c('SpatialPCA.R','SpatialPCA_EstimateLoading.R','SpatialPCA_SpatialPCs.R')
for (f in required) {
  path <- file.path(upstream,'R',f)
  if (!file.exists(path)) stop(paste('Missing official upstream R source:',path))
}
if (!requireNamespace('RSpectra',quietly=TRUE)) stop('Please install the RSpectra R package')
# Keep a record of the source commit separately in the GitHub Actions logs.
# Source the actual original code, including original class declaration, unmodified.
library(methods)
library(RSpectra)
for (f in required) source(file.path(upstream,'R',f),local=FALSE,echo=FALSE)

read_numeric <- function(name) {
  out <- as.matrix(read.csv(file.path(workdir,paste0(name,'.csv')),
                            header=FALSE,check.names=FALSE))
  storage.mode(out) <- 'double'
  out
}
Y <- read_numeric('Y')
coords <- read_numeric('coords')
K <- read_numeric('K')
n <- ncol(Y)
if (nrow(coords)!=n || !all(dim(K)==c(n,n))) stop('Shared input dimensions inconsistent')
if (any(!is.finite(Y)) || any(!is.finite(K))) stop('Nonfinite inputs')
if (components<1L || components>=nrow(Y)) stop('Invalid component count')

# Construct original S4 object without calling SCTransform/SPARK: the objective
# is R-vs-Python *model-core* parity given identical already normalized Y and K.
obj <- methods::new('SpatialPCA',
  project='R-core-parity', counts=Y, normalized_expr=Y,
  location=coords, kernelmat=K, covariate=NULL,
  sparseKernel=FALSE, fast=FALSE,
  params=list(expr=Y))
cat('Original package version/source is from:',upstream,'\n')
cat('Original R input matrix:',nrow(Y),'genes x',ncol(Y),'spots\n')

obj <- SpatialPCA_EstimateLoading(obj, maxiter=300,
      initial_tau=1, fast=FALSE, SpatialPCnum=components)
obj <- SpatialPCA_SpatialPCs(obj,fast=FALSE)
W <- as.matrix(obj@W)
# At pinned upstream revision, SpatialPCs is already PCs x spots.
Z <- as.matrix(obj@SpatialPCs)
if (!identical(dim(W),c(as.integer(nrow(Y)),as.integer(components))))
    stop(paste('Unexpected W shape:',paste(dim(W),collapse='x')))
if (!identical(dim(Z),c(as.integer(components),as.integer(ncol(Y)))))
    stop(paste('Unexpected Z shape:',paste(dim(Z),collapse='x')))

outdir <- file.path(workdir,'author_R_original')
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
if (any(!is.finite(W)) || any(!is.finite(Z)) ||
    !is.finite(obj@tau) || !is.finite(obj@sigma2_0) || obj@sigma2_0 <= 0)
    stop('Invalid original R outputs')
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))
write.table(W,file.path(outdir,'W.csv'),sep=',',row.names=FALSE,col.names=FALSE)
write.table(Z,file.path(outdir,'Z.csv'),sep=',',row.names=FALSE,col.names=FALSE)
writeLines(sprintf('{"tau":%.17g,"sigma2":%.17g,"original_R_executed":true}',
                   obj@tau,obj@sigma2_0),file.path(outdir,'parameters.json'))
cat('ORIGINAL R OUTPUT: tau=',obj@tau,' sigma2=',obj@sigma2_0,' W=',
    paste(dim(W),collapse='x'),' Z=',paste(dim(Z),collapse='x'),'\n')
