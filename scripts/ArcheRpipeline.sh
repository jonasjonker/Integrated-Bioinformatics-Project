#!/usr/bin/env bash

# directory used to store data
DATADIR="$1"

# only execute scripts if previous script returned 0
Rscript makeArrowFiles.R     $DATADIR && \
Rscript inferDoublets.R      $DATADIR && \
Rscript createArchRProject.R $DATADIR && exit 0

exit 1
