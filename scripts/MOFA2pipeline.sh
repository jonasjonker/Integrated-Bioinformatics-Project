#!/usr/bin/env bash

usage () {
    echo "Usage: `basename $0` [DATADIR]"
    echo "PreProcess and train MOFA2 object with files in [DATADIR]"
    echo
    echo "Options:"
    echo "  -h        Show this message."
    exit 1
}

[[ "$#" != 1 ]] && usage
[[ "$1" == "-h" ]] && usage

# only execute scripts if previous script returned 0
Rscript createSeuratObject.R    $DATADIR && \
Rscript preProcessMOFA2Object.R $DATADIR && \
Rscript trainMODA2Object.R      $DATADIR && exit 0

exit 1
