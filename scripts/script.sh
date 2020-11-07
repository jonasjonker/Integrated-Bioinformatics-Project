#!/usr/bin/env bash

# directory used to store data
DATA="$1"

# only execute scripts if previous script returned 0
Rscript arrowmaker.R $DATA        && \
Rscript doublet_inference.R $DATA && \
Rscript createproject.R $DATA     && exit 0

exit 1
