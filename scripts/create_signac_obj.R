#!/usr/bin/Rscript
#
# Command line interface 
###############################################################################
library(docopt)

'Create Signac object

Usage:
  create_signac_obj.R <feature_matrix> <metrics>

Options:
  -h --help     Show this screen.
  --version     Show version.
' -> doc

args <- docopt(doc, version = '0.1')

###############################################################################
# 
###############################################################################
