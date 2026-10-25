#!/bin/sh
# Run gene_dm.R in the pinned DSS container (Bioconductor DSS 2.58.0, R 4.5)
# from the repository root; results land in analysis/ROBERTS2022_PH_PLOIDY_WGBS/results/.
set -eu
repo=$(cd "$(dirname "$0")/../.." && pwd)
docker run --rm --platform linux/amd64 -v "$repo":/repo -w /repo \
  quay.io/biocontainers/bioconductor-dss:2.58.0--r45h01b2380_0 \
  Rscript analysis/ROBERTS2022_PH_PLOIDY_WGBS/gene_dm.R /repo
