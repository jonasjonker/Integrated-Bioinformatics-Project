#!/usr/bin/env bash

usage () {
    echo "Usage: `basename $0` <10XDIR> <CELLANNOTATION> <PEAKDATA> <FACTORS> <OUTFILE>"
    echo "Pre-process and train MOFA2 object."
    echo
    echo "Options:"
    echo "  -h, --help                 Show this message."
    exit 0
}

[[ "$#" != 5 ]] && usage
[[ "$1" == "-h" ]] && usage
[[ "$1" == "--help" ]] && usage
[[ ! -d "$1" ]] && echo "\"$1\" is not a valid directory." && exit 1
[[ ! -f "$2" ]] && echo "\"$2\" is not a valid file." && exit 1
[[ ! -f "$3" ]] && echo "\"$3\" is not a valid file." && exit 1
[[ $4 =~ ^[0-9]+$ ]] || (echo "\"$4\" is not an integer." && exit 1)

XXDIR="$1"
CELLANNOTATION="$2"
PEAKDATA="$3"
FACTORS="$4"
OUTFILE="$5"

echo "start analysis"

Rscript createSeuratObject.R    $XXDIR $CELLANNOTATION $PEAKDATA $OUTFILE && \
Rscript createMOFA2Object.R     $FACTORS $OUTFILE && exit 0

exit 1
