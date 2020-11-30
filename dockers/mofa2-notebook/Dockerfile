FROM jonasjonker/r-bioc:0.0.1

RUN     R -e "devtools::install_github('bioFAM/MOFA2', build_opts = c('--no-resave-data --no-build-vignettes'))"

CMD jupyter notebook --ip 0.0.0.0
