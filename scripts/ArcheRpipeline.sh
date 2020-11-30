#!/usr/bin/env bash

usage () {
    echo "Usage: `basename $0` <DATADIR>"
    echo "Create and process ArrowFiles from data in <DATADIR>."
    echo
    echo "Options:"
    echo "  -h, --help                 Show this message."
    exit 0
}

[[ "$#" != 1 ]] && usage
[[ "$1" == "-h" ]] && usage
[[ "$1" == "--help" ]] && usage
[[ ! -d "$1" ]] && echo "\"$1\" is not a valid directory." && exit 1

DATADIR="$1"

Rscript makeArrowFiles.R     $DATADIR && \
Rscript inferDoublets.R      $DATADIR && \
Rscript createArchRProject.R $DATADIR && exit 0

exit 1
