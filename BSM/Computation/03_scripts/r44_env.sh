#!/usr/bin/env bash
set -euo pipefail

# Bypass the shared R wrapper because its etc/ldpaths contains a stray final
# `_LIBRARY_PATH` command.  The underlying R executable is used read-only with
# all relevant paths declared here.
export R_HOME=/home/data/t200301/R-4.4.0-install/lib/R
export R_LIBS_USER=/home/data/t200301/R/x86_64-pc-linux-gnu-library/4.4
export LD_LIBRARY_PATH=/home/data/t200301/R-4.4.0-install/lib/R/lib:/home/data/t200301/icu/lib:/home/data/t200301/libiconv/lib:/usr/lib/x86_64-linux-gnu:/home/data/t200301/miniconda3/envs/phylo-core/lib/jvm/lib/server

# Avoid nested BLAS/OpenMP oversubscription; BioGeoBEARS controls its own worker
# count through num_cores_to_use.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

exec /home/data/t200301/R-4.4.0-install/lib/R/bin/exec/R --vanilla --slave "$@"
