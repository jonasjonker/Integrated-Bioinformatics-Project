FROM jonasjonker/signac-notebook:0.0.2

RUN     R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())" \
    &&  R -e "ArchR::installExtraPackages()"

RUN     R -e "devtools::install_github('bioFAM/MOFA2', build_opts = c('--no-resave-data --no-build-vignettes'))"

CMD jupyter notebook --ip 0.0.0.0
